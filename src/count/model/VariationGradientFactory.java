package count.model;
/*
 * Copyright 2024 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_MODEL_ROOT_PRIOR;
import static count.io.CommandLine.OPT_TRUNCATE;
import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_LOSS;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;
import java.util.function.IntFunction;

import count.Count;
import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.ds.UniqueProfileTable;
import count.io.CommandLine;
import count.io.RateVariationParser;
import count.matek.Logarithms;
import count.model.LogGradient.PosteriorStatistics;
import count.model.RateVariationLogGradient.SampleLogGradient;

/**
 * Computations for the gradient in rate-variation model.  
 * Used in main BFGS implementation {@link MLRateVariation}.
 * 
 * @author csuros
 *
 */
public class VariationGradientFactory implements Count.UsesThreadpool
{
	public static final int PARAMETER_MOD_LENGTH = 0;
	public static final int PARAMETER_MOD_DUPLICATION = 1;
	/**
	 * Thread pool used across different calls 
	 * 
	 */
	private static ForkJoinPool thread_pool=null;
	/**
	 * Initialized only once, if {@link Count#THREAD_PARALLELISM} is greater than 1.
	 * 
	 * @return
	 */
	protected synchronized static ForkJoinPool threadPool()
	{
		if (thread_pool == null && 1<Count.THREAD_PARALLELISM) 
		{
			thread_pool = Count.threadPool(); 
		}
		return thread_pool;
	}
	
	/**
	 * The model we work with
	 */
	private final RateVariationModel variation_model;
	private final List<CategoryLogGradient> category_gradients;
	
	private final ProfileTable orig_table;
	private final UniqueProfileTable utable;

	/**
	 * Stored here for simpler access
	 */
	private final int num_nodes;
	private int min_copies;
	
	
	public VariationGradientFactory(RateVariationModel variation_model, ProfileTable table)
	{
		this.orig_table = table;
		if (table instanceof UniqueProfileTable)
		{
			this.utable = (UniqueProfileTable) table;
		} else
		{
			this.utable = new UniqueProfileTable(table);
		}
	
		this.variation_model = variation_model;
		this.min_copies = Integer.min(2,utable.minCopies());
		this.num_nodes = variation_model.getBaseModel().getTree().getNumNodes();
		this.category_gradients = new ArrayList<>();
		this.computeClasses();
	}
	
	/**
	 * Convenience method for no rate variation (i.e., 1 class with base model's rates)
	 * @param rates
	 * @param table
	 * @return
	 */
	public static VariationGradientFactory newFactory(TreeWithLogisticParameters rates, ProfileTable table)
	{
		final RateVariationModel variation_model = new RateVariationModel(rates);
		variation_model.initConstantRates();
		VariationGradientFactory newGradientFactory = new VariationGradientFactory(variation_model, table);
		return newGradientFactory;
	}
	
	private int threshold_width_absolute = Integer.MAX_VALUE;
	private double threshold_width_relative = Double.POSITIVE_INFINITY;
	
	private double ancestor_deviation = 0.0; // set to 0 for don't care
	
	
	
	
	

	public void setAncestorWidthThreshold(double deviation)
	{
		this.ancestor_deviation = deviation;
	}
	

	/**
	 * Sets the observation bias : minimum number of observed copies in the families. 
	 * @param min_copies 0, 1, or 2
	 */
	public void setMinimumObservedCopies(int min_copies)
	{
		if (min_copies<0 || min_copies>2)
			throw new UnsupportedOperationException("Minimum copies must be 0,1 or 2.");
		this.min_copies = min_copies;
	}
	
	public int getMinimumObservedCopies()
	{
		return min_copies;
	}
	
	/**
	 * Sets up truncated likelihood calculations: 
	 * at ancestors, the maximum copy number assumed
	 * is given by the minimum between the absolute
	 * threshold and the product of the relative threshold 
	 * times maximum copy number observed. 
	 * 
	 * @param absolute
	 * @param relative
	 */
	public void setCalculationWidthThresholds(int absolute, double relative)
	{
		this.threshold_width_absolute = absolute;
		this.threshold_width_relative = relative;
	}
	
	public int getCalculationWidthAbsolute()
	{
		return threshold_width_absolute; 
	}
	
	public double getCalculationWidthRelative()
	{
		return threshold_width_relative;
	}
	
	/**
	 * Direct access to underlying model; call {@link #computeClasses()} if changes
	 * @return
	 */
	public RateVariationModel getVariationModel()
	{
		return this.variation_model;
	}
	

	private int getNumCategories()
	{
		return variation_model.getNumClasses();
	}

	/**
	 * Called at instantiation, and whenever base model or the category modifiers 
	 * change.
	 */
	public void computeClasses()
	{
		category_gradients.clear();
		for (int k=0; k<getNumCategories(); k++)
		{
			CategoryLogGradient CG = new CategoryLogGradient(k);
			category_gradients.add(CG);
		}
	}
	
	
	/*
	 * 
	 * Computations
	 * 
	 */
	
	
	/**
	 * Corrected log-likelihood, by minimum unobserved copies.
	 * Uses multi-threading, if enabled. 
	 * 
	 * @return
	 */
	public double getCorrectedLL()
	{
		double LLuncorr = computeLogLikelihood();

		MixedProfile empty = new MixedProfile(-1);
		double LL0 = empty.getLL();
		double logitL0 = Logarithms.logToLogit(LL0);
		
		double F = utable.getTotalFamilyCount();
		
		double getCorrectedLL = LLuncorr - F*Logarithms.logitToLogComplement(logitL0);
		
		return getCorrectedLL;
	}

	public double[] getCorrectedLogPosteriorCatProbabilities()
	{
		double[] log_c = computeLogPosteriorCatProbabilities();
		MixedProfile empty = new MixedProfile(-1);
		double LL0 = empty.getLL();
		double logitL0 = Logarithms.logToLogit(LL0);
		
		double F = utable.getTotalFamilyCount();
		double corr0 = Math.log(F)+logitL0; // this many unobserved families added
		int ncat = getNumCategories();
		for (int k=0; k<ncat; k++)
		{
			log_c[k] = Logarithms.add(log_c[k], corr0+empty.getLogPosteriorCatProbability(k));
		}
		return log_c;
	}
	
	public double[] getCorrectedGradient(boolean want_duplication)
	{
		SampleGradient SG = computeSampleGradient(want_duplication);
		SG.correctForUnobserved(want_duplication);
		return SG.get();
	}
	
	private double[] computeLogPosteriorCatProbabilities()
	{
		int nF = utable.getFamilyCount();
		
		final int unit_task = Count.unitTask(nF);
		ForkJoinPool thread_pool = threadPool();

		class PartialP extends RecursiveTask<double[]>
		{
			/**
			 * First family index.
			 */
			private final int minF;
			/**
			 * Last family index, exclusive.
			 */
			private final int maxF;

			PartialP(int min, int max)
			{
				this.minF = min;
				this.maxF = max;
			}
			
			@Override
			public double[] compute()
			{
				try
				{
					if (maxF-minF>unit_task)
					{
						int medF = (minF+maxF)/2;
						PartialP left = new PartialP(minF, medF);
						PartialP right = new PartialP(medF, maxF);
						right.fork();
						double[] log_pL = left.compute();
						double[] log_pR = right.join();
						int ncat = log_pL.length;
						for (int k=0; k<ncat; k++)
							log_pL[k] = Logarithms.add(log_pL[k], log_pR[k]);
						return log_pL;
					} else
					{
						int ncat = getNumCategories(); //rate_categories.size();
						double[] log_p = new double[ncat];
						for (int k=0; k<ncat; k++)
							log_p[k] = Double.NEGATIVE_INFINITY;
						for (int f=minF; f<maxF; f++)
						{
							MixedProfile P = new MixedProfile(f);
							double mul = utable.getMultiplicity(f);
							double log_mul = Math.log(mul);
							for (int k=0; k<ncat; k++)
							{
								log_p[k] = Logarithms.add(log_p[k], P.getLogPosteriorCatProbability(k)+log_mul);
							}
						}
						return log_p;
					}
				} catch (Throwable t)
				{
					throw t;
				}
			} // compute
			
		} // task PartialP
		PartialP bigjob = new PartialP(0,nF);
		double[] log_cnt;
		try 
		{
			if (nF>unit_task)
			{
				log_cnt = thread_pool.invoke(bigjob);
			} else
			{
				log_cnt = bigjob.compute();
			}
		} catch (Throwable t)
		{
			// could be out of memory, or out of threads, or numerical error, or whatever else			
			throw new RuntimeException(t);  
		}
		return log_cnt;
	}
	
	/**
	 * Computes the uncorrected log-likelihood across the entire table.
	 * Uses multi-threading, if enabled. 
	 * 
	 * @return
	 */
	private double computeLogLikelihood()
	{
		int nF = utable.getFamilyCount();
		
		final int unit_task = Count.unitTask(nF);
		ForkJoinPool thread_pool = threadPool();
		
		class PartialL extends RecursiveTask<Double>
		{
			/**
			 * First family index.
			 */
			private final int minF;
			/**
			 * Last family index, exclusive.
			 */
			private final int maxF;

			PartialL(int min, int max)
			{
				this.minF = min;
				this.maxF = max;
			}
			
			@Override
			protected Double compute() 
			{
				try
				{
					if (maxF-minF>unit_task)
					{
						int medF = (minF+maxF)/2;
						PartialL left = new PartialL(minF, medF);
						PartialL right = new PartialL(medF, maxF);
						right.fork();
						return left.compute()+right.join();
					} else
					{
						double LL = 0.0;
						for (int f =minF; f<maxF; f++)
						{
							MixedProfile P = new MixedProfile(f);
							double mul = utable.getMultiplicity(f);
							LL += mul*P.getLL();
						}
						return LL;
					}
				} catch (Throwable t)
				{
					throw t;
				}
			}
		} // task PartialL
		PartialL bigjob = new PartialL(0,nF);
		double LL;
		try
		{
			if (nF>unit_task)
			{
				LL = thread_pool.invoke(bigjob);
			} else
			{
				LL = bigjob.compute();
			}
		} catch (Throwable t)
		{
			// could be out of memory, or out of threads, or numerical error, or whatever else			
			throw new RuntimeException(t);  
		}
		return LL;
	}
	
	/** 
	 * Computes the gradient across the input table
	 * 
	 * @return
	 */
	public SampleGradient computeSampleGradient(boolean want_duplication)
	{
		int nF = utable.getFamilyCount();
		
		final int unit_task = Count.unitTask(nF);
		ForkJoinPool thread_pool = threadPool();
		
		class PartialG extends RecursiveTask<SampleGradient>
		{
			/**
			 * First family index.
			 */
			private final int minF;
			/**
			 * Last family index, exclusive.
			 */
			private final int maxF;
			
			PartialG(int min, int max)
			{
				this.minF = min;
				this.maxF = max;
			}
			
			@Override 
			protected SampleGradient compute()
			{
				try
				{
					SampleGradient G;
					if (maxF-minF > unit_task)
					{
						int medF = (minF+maxF)/2;
						PartialG left = new PartialG(minF, medF);
						PartialG right = new PartialG(medF, maxF);
						left.fork();
						G = right.compute();
						SampleGradient G2 = left.join();
						G.add(G2, 0.0);
					} else
					{
						G = new SampleGradient();
						for (int f=minF; f<maxF; f++)
						{
							double log_mul = Math.log(utable.getMultiplicity(f));
							G.addUniq(f, log_mul, want_duplication);
						}
					}
					return G;
				} catch (Throwable t)
				{
					throw t;
				}
			} // compute()
		} // end of class def
		
		PartialG bigjob = new PartialG(0, nF);
		try
		{
			SampleGradient G;
			if (nF>unit_task)
			{
				G = thread_pool.invoke(bigjob);
			} else
			{
				G = bigjob.compute();
			}
			return G;
		} catch (Throwable t)
		{
			throw new RuntimeException(t);
		}
	}
	/*
	 * Helper classes
	 */
	
	/**
	 * 
	 * Gradient computations in the variation categories and conversion 
	 * to gradient by base model parameters  
	 */
	private class CategoryLogGradient
	{
		private CategoryLogGradient(int cat)
		{
			this.cat_idx = cat;
			RateVariationModel.Category C = variation_model.getCategory(cat);
			this.lrates = C.getRates();
			LikelihoodParametrized cfactory = new LikelihoodParametrized(lrates, utable);
			this.gradient = new LogGradient(cfactory);
			gradient.setMinimumObservedCopies(VariationGradientFactory.this.min_copies);
			gradient.setCalculationWidthThresholds(VariationGradientFactory.this.threshold_width_absolute, VariationGradientFactory.this.threshold_width_relative);
			gradient.setAncestorWidthThreshold(VariationGradientFactory.this.ancestor_deviation);
			gradient.computeParameters(); // needed?
		}
		final int cat_idx;
		final LogGradient gradient;
		final TreeWithLogisticParameters lrates;
		
		protected double[][] getLogMainGradientCommonGain(double[][] logSD)
		{
			double[][] log_Dp = gradient.convertToLogDistributionGradient(logSD);
			
			double[][] log_Dpλ = LogGradient.convertToLogRelativeRateGradient(lrates, log_Dp);
			log_Dpλ = LogGradient.convertRelativeRateToLogGainGradient(lrates, log_Dpλ, variation_model.getCommonGainType(), variation_model.isUniversalGain());
			
			
			// HERE: check this.convertToLogMainGradient
			double[][] logD = this.convertToLogMainGradient(log_Dpλ);
			
			return logD;
		}
		
		protected double[][] getLogMainDuplicationGradient(double[][] logS_D)
		{
			double[][] log_Dp = gradient.convertToLogDistributionGradient(logS_D);
			double[][] log_Dpq = LogGradient.convertDistributionToLogGainGradient(lrates, log_Dp, variation_model.getCommonGainType(), variation_model.isUniversalGain());
			double[][] logD = this.convertToLogDuplicationGradient(log_Dpq);
			
			return logD;
		}
	
		/**
		 * Calculates gradient by base model's logit-p, logit-lambda, log-gain; 
		 * and rate modifiers. The last two entries of the returned array
		 * are the derivatives by the length-modifier and the duplication-modifier
		 * with indexes 3<var>n</var>+{@link RateVariationLogGradient#PARAMETER_MOD_LENGTH}
		 * and 3<var>n</var>+ {@link RateVariationLogGradient#PARAMETER_MOD_DUPLICATION}.
		 * 
		 * Default implementation just adds then last two cells with 0 gradient; useful for 
		 * initialization.
		 * 
		 * @param logCatGradient category gradient ; unaffected
		 * 	
		 * @return array of length 3<var>n</var>+2 for <var>n</var> nodes in the tree
		 */
		protected double[][] convertToLogMainGradient(double[][] logCatGradient)
		{
			int num_nodes = lrates.getTree().getNumNodes();
			double[][] logMD = new double[3*num_nodes+2][];
			int j=0; 
			while (j<3*num_nodes)
			{
				logMD[j]=logCatGradient[j].clone();
				j++;
			}
			logMD[j++] = Logarithms.ldiff();
			logMD[j++] = Logarithms.ldiff();		
			
			
			TreeWithLogisticParameters main_rates = variation_model.getBaseModel();
			IndexedTree tree = main_rates.getTree();
			
			RateVariationModel.Category C = variation_model.getCategory(cat_idx);
			double mod_duplication = C.getModDuplication();
			double mod_length = C.getModLength();
			if (mod_duplication!=0.0 || mod_length != 0.0)
			{
				// mod indices
				final int j_mod_len = 3*num_nodes + PARAMETER_MOD_LENGTH;
				final int j_mod_dup = 3*num_nodes + PARAMETER_MOD_DUPLICATION;
				// mod gradients
				double[] dmod_len = logMD[j_mod_len];
				double[] dmod_dup = logMD[j_mod_dup];
				for (int v=0; v<num_nodes; v++)
				{
					final int j_loss = 3*v + PARAMETER_LOSS;
					final int j_dup  = 3*v + PARAMETER_DUPLICATION;
					// class log-gradient by logit(p_c), logit(lambda_c)
					double[] dlp = logMD[j_loss];
					double[] dllm = logMD[j_dup];
					
					final double log1_p = main_rates.getLogLossComplement(v);
					final double log_λ = main_rates.getLogRelativeRate(v);
					final double logit_λ = main_rates.getLogitRelativeRate(v);
					if (RateVariationModel.Multiplier.class.equals(C.getClass()))
					{
						if (true) throw new UnsupportedOperationException("Not tested yet");
						if (logit_λ==Double.POSITIVE_INFINITY) // lambda == 1 
						{
							// we go directly 
							if (mod_length != 0.0)
								dmod_len = Logarithms.ldiffAddMultiply(dmod_len, 0.0, dlp, dmod_len);
						} else if (log1_p == Double.NEGATIVE_INFINITY) // p==1.0
						{
							if (mod_duplication != 0.0)
								dmod_dup = Logarithms.ldiffAddMultiply(dmod_dup, 0.0, dllm, dmod_dup);
						} else
						{
							double log_λcat = lrates.getLogRelativeRate(v);
							// convert to transient-gradient
							if (log_λcat != Double.NEGATIVE_INFINITY)
								dllm = Logarithms.ldiffAddMultiply(dllm, log_λcat, dlp, dllm);
							dlp = Logarithms.ldiffInverse(dlp, dlp);

							// multipliers in the Jacobian 
							double logit_σ = main_rates.getLogitTransientParameter(v);
							double loglog_σ = Math.log(-Logarithms.logitToLogValue(logit_σ));
							double log1_σ = Logarithms.logitToLogComplement(logit_σ);
							double log1_σcat = Logarithms.logitToLogComplement(lrates.getLogitTransientParameter(v));
							
							double log1_λ = main_rates.getLogRelativeComplement(v);
							double log1_λcat = lrates.getLogRelativeComplement(v);
							double d1_λ = log1_λcat-log1_λ;
							double d1_σ = log1_σ-log1_σcat;
							double d_lnσ = loglog_σ-log1_σcat;
							
							double ds_ds = mod_length + d1_σ + d1_λ;
							double minus_ds_dmlen = mod_length + d_lnσ + d1_λ;
	
							
							// convert to base-parameter gradients 
							if (log_λ != Double.NEGATIVE_INFINITY)
							{
								if (mod_duplication != 0.0)
								{
									dmod_dup = Logarithms.ldiffAddMultiply(dmod_dup, 0.0, dllm, dmod_dup);
									
									double ds_dmdup = minus_ds_dmlen + d1_λ + log_λ + mod_duplication;
									dmod_dup = Logarithms.ldiffAddMultiply(dmod_dup, ds_dmdup, dlp, dmod_dup);
								}
								
								if (0.0 < mod_duplication)
								{
									assert (log_λ < log_λcat);
									double log_diff =  log_λcat + Logarithms.logToLogComplement(log_λ-log_λcat);
									dllm = Logarithms.ldiffSubtractMultiply(dllm, minus_ds_dmlen+log_diff, dlp, dllm);
								} else if (mod_duplication != 0.0)
								{
									double minus_log_diff = log_λ + Logarithms.logToLogComplement(log_λcat-log_λ);
									dllm = Logarithms.ldiffAddMultiply(dllm, minus_ds_dmlen+minus_log_diff, dlp, dllm);
								}
							}
							
							if (mod_length != 0.0)
							{
								dmod_len = Logarithms.ldiffSubtractMultiply(dmod_len, minus_ds_dmlen, dlp, dmod_len);
							}
							dlp = Logarithms.ldiffMultiply(dlp, ds_ds, dlp);
							
							// convert to loss-gradient
							if (log_λ != Double.NEGATIVE_INFINITY)
								logMD[j_dup] = Logarithms.ldiffAddMultiply(dllm, log_λ, dlp, dllm);
							logMD[j_loss] = Logarithms.ldiffInverse(dlp, dlp);
						} // lambda,p<1 		
					} else
				 	{
						assert (RateVariationModel.LogisticShift.class.equals(C.getClass()));
						
						if (mod_duplication!=0.0 )
						{
							if (log_λ != Double.NEGATIVE_INFINITY)
							{
								if (logit_λ<Double.POSITIVE_INFINITY) // (log_λ != 0.0)  // lambda==1 has no component here)
									dmod_dup = Logarithms.ldiffAddMultiply(dmod_dup, 0.0, dllm, dmod_dup);
								double log_λcat = lrates.getLogRelativeRate(v);
								
//								{ // DEBUG
//									if (!Double.isFinite(log_λcat) || (!(0.0 <mod_duplication) && !(log_λcat-log_λ<=0.0)))
//									{
//										System.out.println("#**VGF.CLG.cTLMG bad cat lambda "+log_λcat+" at node "+v+"\tlm "+log_λ+"\tmoddup "+mod_duplication
//												+"\tlogitlm "+logit_λ
//												+"\tlogitclm "+lrates.getLogitRelativeRate(v)
//												+"\t// "+lrates.toString(v));
//									}
//									
//									assert Double.isFinite(log_λ);
//									
//									
//									assert Double.isFinite(log_λcat);
//								}
								
								if (!tree.isRoot(v) && log1_p != Double.NEGATIVE_INFINITY)
								{
									dmod_dup = Logarithms.ldiffAddMultiply(dmod_dup, log_λcat, dlp, dmod_dup);
									
									if (logit_λ<Double.POSITIVE_INFINITY)
									//if (log_λ != 0.0)
									{
										if (0.0< mod_duplication)
										{
											assert (log_λ < log_λcat);
											double log_diff =  log_λcat + Logarithms.logToLogComplement(log_λ-log_λcat);
											logMD[j_dup] = Logarithms.ldiffAddMultiply(dllm, log_diff, dlp, dllm);
										} else
										{
											double minus_log_diff = log_λ + Logarithms.logToLogComplement(log_λcat-log_λ);
											logMD[j_dup] = Logarithms.ldiffSubtractMultiply(dllm, minus_log_diff, dlp, dllm);
										}
									}
								}
							}
						}
						if (mod_length != 0.0)
						{
							if (!tree.isRoot(v) && log1_p != Double.NEGATIVE_INFINITY)
							{
								dmod_len = Logarithms.ldiffAddMultiply(dmod_len, 0.0, dlp, dmod_len);
							}
						}
				 	}
					
				} // for v 	
				logMD[j_mod_len] = dmod_len;
				logMD[j_mod_dup] = dmod_dup;					
			} // if any mods 

			return logMD;
		}
		
		/**
		 * Calculates gradient by base model's logit(q); other derivatives are 
		 * unchanged (and possibly incorrect).
		 * 
		 * Default implementation just adds the last two cells with 0 gradient.
		 * 
		 * @param logCatGradient category gradient from {@link #convertToLogCatGradient(double[][])}; unaffected
		 * 	
		 * @return array of length 3<var>n</var>+2 for <var>n</var> nodes in the tree		 
		 */
		protected double[][] convertToLogDuplicationGradient(double[][] logCatGradient)
		{
			int num_nodes = lrates.getTree().getNumNodes();
			double[][] logMD = new double[3*num_nodes+2][];
			int j=0; 
			while (j<3*num_nodes)
			{
				logMD[j]=logCatGradient[j].clone();
				j++;
			}
			logMD[j++] = Logarithms.ldiff();
			logMD[j++] = Logarithms.ldiff();		

			TreeWithLogisticParameters main_rates = variation_model.getBaseModel();
			
			RateVariationModel.Category C = variation_model.getCategory(cat_idx);
			double mod_duplication = C.getModDuplication();
			if (mod_duplication!=0.0 )
			{
				for (int v=0; v<num_nodes; v++)
				{
					final int j_loss = 3*v + PARAMETER_LOSS;
					final int j_dup  = 3*v + PARAMETER_DUPLICATION;
					// class log-gradient by logit(p_c), logit(lambda_c)
					double[] dlp = logMD[j_loss];
					double[] dlq = logMD[j_dup];
					
					final double log1_q = main_rates.getLogDuplicationComplement(v);
					if (main_rates.getLogitDuplicationParameter(v)<Double.POSITIVE_INFINITY)   // (log1_q != 0.0)
					{
						// Polya
						double log1_qcat = lrates.getLogDuplicationComplement(v);
						logMD[j_dup] = Logarithms.ldiffMultiply(dlq, log1_q-log1_qcat, dlq);
					}
				}
			}
			return logMD;
		}
	}

	/**
	 * Class for computing likelihoods, posteriors and gradients for a family  
	 * 
	 */
	public class MixedProfile
	{
		private final Posteriors.Profile[] cat_profiles;
		private final PosteriorStatistics[] cat_stats;
		private final double[] cat_LL;
		private final double LL;
		private final int uniq_fam;
		private final double[] cat_log_post;
		
		private double[][] log_events = null; // lazy initialization
		
		/**
		 * 
		 * @param uniq_f -1 for unobserved profiles
		 */
		private MixedProfile(int uniq_f)
		{
			int ncat = getNumCategories();
			this.cat_LL = new double[ncat];
			this.cat_stats = new PosteriorStatistics[ncat];
			this.cat_log_post = new double[ncat];
			
			if (0<=uniq_f)
			{
				this.uniq_fam=uniq_f;
				this.cat_profiles = new Posteriors.Profile[ncat];
				
				
				double LLf = Double.NEGATIVE_INFINITY;
				for (int k=0; k<ncat; k ++)
				{
					RateVariationModel.Category C = variation_model.getCategory(k);
					double log_pclass = C.getLogCatProbability();
					LogGradient G = category_gradients.get(k).gradient;
					
					Posteriors.Profile P = G.getPosteriors(uniq_f);
					PosteriorStatistics S = G.getProfileStatistics(P);
					
					double cLL = P.inside.getLogLikelihood();
					if (log_pclass != Double.NEGATIVE_INFINITY)
						LLf = Logarithms.add(LLf, cLL+log_pclass);
	
					this.cat_profiles[k] = P;
					this.cat_LL[k] = cLL;
					this.cat_stats[k] = S;
				}
				this.LL = LLf;
			} else
			{ // unobserved
				this.uniq_fam=-1;
				this.cat_profiles = null;
				double LL0 = Double.NEGATIVE_INFINITY;
				for (int k=0; k<ncat; k ++)
				{
					RateVariationModel.Category C = variation_model.getCategory(k);
					double log_pclass = C.getLogCatProbability();
					LogGradient G = category_gradients.get(k).gradient;
					PosteriorStatistics S = G.getUnobservedStatistics();
					double cLL = S.getLogLikelihood();
					if (log_pclass != Double.NEGATIVE_INFINITY)
						LL0 = Logarithms.add(LL0, cLL+log_pclass);
					this.cat_LL[k] = cLL;
					this.cat_stats[k] = S;
				}			
				this.LL = LL0;
			}
			for (int k=0; k<ncat; k ++)
			{
				cat_log_post[k] = variation_model.getCategory(k).getLogCatProbability()+ this.cat_LL[k]-this.LL;
			}
		}
		
		public boolean isUnobserved()
		{
			return this.uniq_fam == -1;
		}
		
		public double getLogPosteriorCatProbability(int cat)
		{
			return cat_log_post[cat]; // variation_model.getCategory(cat).getLogCatProbability()+ this.cat_LL[cat]-this.LL;
		}
		
		/**
		 * Posterior category probability
		 * 
		 * @param cat
		 * @return
		 */
		public double getClassProbability(int cat)
		{
			return Math.exp(getLogPosteriorCatProbability(cat));
		}		

		/**
		 * Expected ancestral copy number at the node
		 * 
		 * @param node
		 * @return
		 */
		public double getNodeMean(int node)
		{
			assert !isUnobserved();
			
			double logN = sumPostLog(k->cat_profiles[k].getLogNodeMean(node));
			
			return Math.exp(logN);
		}		
		
		/**
		 * Expected number of ancestral copies inherited from parent
		 * 
		 * @param node
		 * @return
		 */
		public double getEdgeMean(int node)
		{
			assert !isUnobserved();

			double logS = sumPostLog(k->cat_profiles[k].getLogEdgeMean(node));
			return Math.exp(logS);
		}
		
		/**
		 * Expected number of new ancestral copies at node (by duplication and gain)
		 * 
		 * @param node
		 * @return
		 */
		public double getBirthMean(int node)
		{
			double logN_S = this.sumPostLog(k->cat_profiles[k].getLogNodeIncrease(node));
			return Math.exp(logN_S);
		}
		
		/**
		 * Expected number of non-inherited ancestral copies at parent.
		 * 
		 * @param node
		 * @return
		 */
		public double getDeathMean(int node)
		{
			double logN_S = this.sumPostLog(k->cat_profiles[k].getLogEdgeDecrease(node));
			return Math.exp(logN_S);
		}

		/**
		 * Posterior probabilities for copy numbers 0 and 1
		 * 
		 * @param node
		 * @return
		 */
		public double[] getNodeAncestor(int node)
		{
			
			double[] pN=new double[2];
			for (int c=0; c<cat_profiles.length; c++)
			{
				Posteriors.Profile P = cat_profiles[c];
				double[] cN = P.getNodeAncestorPosteriors(node);
				if (pN.length<cN.length)
					pN = Arrays.copyOf(pN, cN.length);
				for (int n=0; n<cN.length; n++)
				{
					pN[n] += Math.exp(cat_log_post[c])*cN[n];
				}
			}
			return pN;
		}
				
		
		/** 
		 * Posterior probability for family-level events.
		 * 
		 * @param node
		 * @param event_type gain, loss, expansion or contraction
		 * @return
		 */
		public double getFamilyEvent(int node, Posteriors.FamilyEvent event_type)
		{
			if (log_events == null)
				initEvents();
			final int event_idx = event_type.ordinal();
			double log_E = log_events[node][event_idx];
			return Math.exp(log_E);
		}		
				
		private void initEvents()
		{
			log_events = new double[num_nodes][];
			int ncat = getNumCategories();
			
			for (int v=0; v<num_nodes; v++)
			{
				double[] node_ev = new double[Posteriors.FamilyEvent.values().length];
				Arrays.fill(node_ev, Double.NEGATIVE_INFINITY);

				for (int k=0; k<ncat; k++)
				{
					Posteriors.Profile P = cat_profiles[k];
					
					double[] ev = P.getFamilyEventPosteriors(v);
					for (int j=0; j<ev.length; j++)
					{
						node_ev[j] = Logarithms.add(node_ev[j], cat_log_post[k]+Math.log(ev[j]));
					}
				}
				log_events[v] = node_ev;
			}
		}
		
		private double sumPostLog(IntFunction<Double> cat_log_stats)
		{
			double sumLog = Double.NEGATIVE_INFINITY;
			int ncat = cat_log_post.length;
			for (int k=0; k<ncat; k++)
			{
				sumLog = Logarithms.add(sumLog, cat_log_post[k]+cat_log_stats.apply(k));
			}
			return sumLog;
		}
		
		/**
		 * Log-gradient by base model's logit(p), logit(lambda), log(gain)
		 * and rate modifiers. 
		 * For <var>n</var> tree nodes and <var>K</var> categories, 
		 * the returned array has 3<var>n</var>+2<var>K</var> entries.  
		 * The last 
		 * 2<var>K</var> entries of the returned array
		 * are the derivatives by the modifiers
		 * with indexes 3<var>n</var>+2<var>k</var> + {@link RateVariationLogGradient#PARAMETER_MOD_LENGTH}
		 * and 3<var>n</var>+2<var>k</var> + {@link RateVariationLogGradient#PARAMETER_MOD_DUPLICATION}
		 * for class <var>k</var>=0,...,<var>K</var>-1.
		 * 
		 * @return array of partial derivatives on logdiff scale 
		 */
		protected double[][] getLogGradient()
		{
			int ncat = getNumCategories();
			
			double[][] log_D = new double[3*num_nodes+2*ncat][];
			for (int j=0; j<3*num_nodes; j++)
			{
				log_D[j] = Logarithms.ldiff();
			}
			for (int k=0; k<ncat; k ++)
			{
				RateVariationModel.Category C =  variation_model.getCategory(k);
				PosteriorStatistics S = this.cat_stats[k];
				
				double[][] log_SD = S.getLogSurvivalGradient();
				CategoryLogGradient cfactory = category_gradients.get(k);
				double[][] log_Dcat = cfactory.getLogMainGradientCommonGain(log_SD);
				double log_pc = this.getLogPosteriorCatProbability(k);
				for (int j=0; j<3*num_nodes; j++)
				{
					log_D[j] = Logarithms.ldiffAddMultiply(log_D[j], log_pc, log_Dcat[j], log_D[j]);
				}
				int j_modlen = 3*num_nodes + PARAMETER_MOD_LENGTH;
				int j_moddup = 3*num_nodes + PARAMETER_MOD_DUPLICATION;
				
				// shift to position
				log_D[j_modlen+2*k] = log_Dcat[j_modlen];
				log_D[j_moddup+2*k] = log_Dcat[j_moddup];
			}	
//			// convert to gain
//			log_D = LogGradient.convertFromLogGainGradient(main_rates, log_D, common_gain_by, use_universal_gain);
			
			return log_D;
		}
		
		/**
		 * Log-gradient by logit(p), logit(q), log(commongain).
		 * 
		 * @return
		 */
		protected double[][] getLogDuplicationGradient()
		{
			int ncat = getNumCategories();
			
			double[][] log_D = new double[3*num_nodes+2*ncat][];
			for (int j=0; j<3*num_nodes; j++)
			{
				log_D[j] = Logarithms.ldiff();
			}
			for (int k=0; k<ncat; k ++)
			{
				RateVariationModel.Category C =  variation_model.getCategory(k);
				PosteriorStatistics S = this.cat_stats[k];
				
				double[][] log_SD = S.getLogSurvivalGradient();
				CategoryLogGradient cfactory = category_gradients.get(k);
				double[][] log_Dcat = cfactory.getLogMainDuplicationGradient(log_SD);
				double log_pc = this.getLogPosteriorCatProbability(k);
				for (int j=0; j<3*num_nodes; j++)
				{
					log_D[j] = Logarithms.ldiffAddMultiply(log_D[j], log_pc, log_Dcat[j], log_D[j]);
				}
				int j_modlen = 3*num_nodes + PARAMETER_MOD_LENGTH;
				int j_moddup = 3*num_nodes + PARAMETER_MOD_DUPLICATION;
				
				// shift to position
				log_D[j_modlen+2*k] = log_Dcat[j_modlen];
				log_D[j_moddup+2*k] = log_Dcat[j_moddup];
			}	
			
			return log_D;
		}
		
		protected double getLL()
		{
			return this.LL;
		}
	}
	
	/**
	 * Class for summing the log-likelihood gradients across the sample
	 * and for storing intermediate results.  
	 */
	public class SampleGradient
	{
		private double profile_count;
		private double LL;
		private final double[] G;
		private final double[] log_post;
		
		private SampleGradient()
		{
			int ncat = getNumCategories();
			this.G = new double[3*num_nodes+2*ncat];
			this.log_post = new double[ncat];
			
			this.init();
		}
		
		private void init()
		{
//			for (int j=0; j<log_D.length; j++)
//				log_D[j]=Logarithms.ldiff();

			this.profile_count = 0.0;
			this.LL = 0.0;
			for (int k=0; k<log_post.length; k++)
				log_post[k] = Double.NEGATIVE_INFINITY;
			Arrays.fill(this.G, 0.0);
		}
		
		public double getLogLikelihood()
		{
			return this.LL;
		}
		
		
		public double getLogPosteriorCount(int cat)
		{
			return this.log_post[cat];
		}
		
		public double getLogPosteriorFrequency(int cat)
		{
			return (profile_count==0.0
					?Double.NEGATIVE_INFINITY
					:this.log_post[cat]-Math.log(profile_count));
		}
		
		public double getProfileCount()
		{
			return this.profile_count;
		}
		
		/**
		 * Adds a family with given multiplicity to the sample; log-likelihood 
		 * is updated to log-product of family likelihoods
		 * 
		 * @param uniq_f -1 if unobserved
		 * @param log_multiplier
		 */
		private void addUniq(int uniq_f, double log_multiplier, boolean want_duplication)
		{
			double mul = Math.exp(log_multiplier);

			MixedProfile MP = new MixedProfile(uniq_f);
			double[][] log_D;
			if (want_duplication)
				log_D = MP.getLogDuplicationGradient();
			else
				log_D = MP.getLogGradient();
			double LLf = MP.getLL();
			
			this.profile_count += mul;
			this.LL += mul*LLf; 
			
			for (int k=0; k<log_post.length; k++)
			{
				double log_pcat = MP.getLogPosteriorCatProbability(k);
				this.log_post[k] = Logarithms.add(log_post[k],log_pcat+log_multiplier);
			}
			this.addLogGradient(log_D, log_multiplier);
		}
		
		protected void add(SampleGradient that, double log_multiplier)
		{
			double mul = Math.exp(log_multiplier);
			this.profile_count += mul * that.profile_count ;
			this.LL += mul*that.LL;
			for (int k=0; k<log_post.length; k++)
			{
				this.log_post[k] = Logarithms.add(this.log_post[k], that.log_post[k]+log_multiplier);
			}		
			
			assert (this.G.length == that.G.length);
			for (int j=0; j<G.length; j++)
			{
				this.G[j] += mul*that.G[j];
			}
		}
//		/**
//		 * Sets the corrected gradient with unobserved profiles
//		 */
//		public void correctForUnobserved()
//		{
//			this.correctForUnobserved(false);
//		}
		
		/**
		 * Sets the corrected gradient with unobserved profiles
		 */
		public void correctForUnobserved(boolean want_duplication)
		{
			if (0<min_copies) // otherwise no correction
			{
				SampleGradient G0 = new SampleGradient();
				G0.addUniq(-1, 0.0, want_duplication); // unobserved
				
				double LL0 = G0.LL;
				double logitL0 = Logarithms.logToLogit(LL0);
				
				double LLuncorr = this.LL;
				double F = this.profile_count;
				double logF = Math.log(F);
				
				
				
				// calculate corrected gradient and
				// posterior counts
				this.add(G0, logF+logitL0);
				
				this.LL = LLuncorr - F*Logarithms.logitToLogComplement(logitL0);
			}
		}	
		/**
		 * The gradient of the log-likelihood
		 * @return
		 */
		public double[] get()
		{
			return this.G;
		}
			
		/**
		 * Auxiliary method for aggregating log-gradients
		 * 
		 * @param log_D
		 * @param log_multiplier
		 */
		private void addLogGradient(double[][] log_D, double log_multiplier)
		{
			assert (G.length == 3*num_nodes+2*getNumCategories());
			assert (G.length == log_D.length);

			for (int j=0; j<G.length; j++)
			{
				double[] log_dD = Logarithms.ldiffMultiply(log_D[j], log_multiplier, null);
				double dD =  Logarithms.ldiffValue(log_dD);
				this.G[j] += dD;
			}
		}
		
	}
	
	
	/*
	 * 
	 * Test code
	 * 
	 */
	private void testConstantRates()
	{
		variation_model.initConstantRates();
		this.computeClasses();
		TreeWithLogisticParameters main_rates = variation_model.getBaseModel();
		
		java.io.PrintStream out = System.out;
		int num_nodes = main_rates.getTree().getNumNodes();
		
		SampleGradient SG = computeSampleGradient(false);
		double uncorrLL = SG.getLogLikelihood();
		SG.correctForUnobserved(false);
		double corrLL = SG.getLogLikelihood();
		
		double[] dLL =  SG.get(); 
		
		LogGradient altG = new LogGradient(main_rates, this.orig_table);
		altG.setMinimumObservedCopies(this.getMinimumObservedCopies());
		altG.setCalculationWidthThresholds(this.getCalculationWidthAbsolute(), this.getCalculationWidthRelative());
		PosteriorStatistics S = altG.getSampleStatistics();
		
		out.println("# Uncorrected likelihood "+uncorrLL+"\t// "+S.getLogLikelihood());
		double[][] alt_logD = altG.correctedLogSurvivalGradient(S);
		out.println("# Corrected likelihood "+corrLL+"\t// "+S.getLogLikelihood()+"\t// recomputed "+getCorrectedLL());
		
		alt_logD = altG.convertToLogDistributionGradient(alt_logD);		
		double[][] alt_log_Dpλ = LogGradient.convertToLogRelativeRateGradient(main_rates, alt_logD);
		alt_log_Dpλ = LogGradient.convertRelativeRateToLogGainGradient(main_rates, alt_log_Dpλ, variation_model.getCommonGainType(), variation_model.isUniversalGain());
		double[] A = LogGradient.toGradient(alt_log_Dpλ);
		
		out.println("# Model logistic/log gradient (GLD order)"
				+"; gain-by "+variation_model.getCommonGainType()+"/"+(variation_model.isUniversalGain()?"univ":"lin"));
		Gradient.printGradient(out, dLL, num_nodes);
		
		// DEBUG: via LogGradient 		
		out.println("# Alt logistic/log gradient gradient (GLD order) via "+altG.getClass().getName());
		Gradient.printGradient(out, A, num_nodes);

		for (int j=0; j<3*num_nodes; j++)
		{
			double delta = A[j]-dLL[j];
			double rd = delta/Math.abs(dLL[j]);
			out.println("#**LG.mm par "+(j/3)+"/"+(j%3)+"\tdelta "+delta+"\trd "+rd);
		}
	}


	/**
	 * Command-line execution: for testing during development
	 * (calculates gradients in two ways). 
	 *  
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception
	{
		Class<?> us = java.lang.invoke.MethodHandles.lookup().lookupClass();
		count.io.CommandLine cli = new count.io.CommandLine(args,us);
		
		java.io.PrintStream out = System.out;
		
		
		RateVariationModel model = null;
		
    	TreeWithRates starting_rates;
		count.matek.DiscreteDistribution root_prior = cli.getOptionDistribution(OPT_MODEL_ROOT_PRIOR, null);
    	if (cli.getVariationModel()==null)
    	{
    		Random RND = cli.getOptionRND(out);
    		starting_rates = new TreeWithRates(cli.getTree(), RND);
    		if (root_prior != null)
    		{
    			starting_rates.setRootDistribution(root_prior);
    			out.println(CommandLine.getStandardHeader("Root prior preset: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
    		}
    		starting_rates.initNodeParameters(starting_rates.getTree().getRoot());
			out.println(CommandLine.getStandardHeader("(Root prior random: "+starting_rates.getRootDistribution()+")"));
    		model = new RateVariationModel(starting_rates); // 
    		model.initConstantRates();
    		starting_rates = model.getBaseModel();
    	} else
    	{
    		starting_rates = cli.getRates(); 
    		model = cli.getVariationModel();
    		assert starting_rates == model.getBaseModel();
    		Random RND = cli.getOptionRND(out);
    		if (RND!=null)
    			starting_rates.setRandom(RND);
    		if (root_prior != null)
    		{
    			starting_rates.setRootDistribution(root_prior);
    			out.println(CommandLine.getStandardHeader("Root prior preset: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
    		}
    	}
		
		
		VariationGradientFactory G = new VariationGradientFactory(model, cli.getTable());
        if (cli.getOptionValue(OPT_TRUNCATE)!=null)
        {
        	int absolute = cli.getOptionTruncateAbsolute();
        	double relative = cli.getOptionTruncateRelative();
        	G.setCalculationWidthThresholds(absolute, relative);
    		System.out.println(CommandLine.getStandardHeader("Truncated computation (absolute,relative)="
    				+absolute+","+relative));
        } 
		int min_copies = cli.getOptionInt(OPT_MINCOPY, G.min_copies);
		G.setMinimumObservedCopies(min_copies);
		
		
		
		if (model.getNumClasses()==1)
		{
			G.testConstantRates();
		} 
//		else
//		{
//			G.testCenter();
//		}
		//G.testVariableRates();
		
		out.println(RateVariationParser.printRates(G.getVariationModel()));
	}		
	

}
