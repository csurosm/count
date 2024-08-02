
 /* Copyright 2021 Mikl&oacute;s Cs&#369;r&ouml;s.
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
package count.model;


import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;
import java.util.Arrays;

import count.Count;
import count.ds.ProfileTable;
import count.ds.UniqueProfileTable;
import count.io.CommandLine;
import count.io.NewickParser;
import count.matek.ArraySum;
import count.matek.FunctionMinimization;
import count.matek.Logarithms;

import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_TRUNCATE;
import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LENGTH;
import static count.model.GLDParameters.PARAMETER_LOSS;


/**
 * Computation of likelihood derivatives.
 * 
 * @author Miklós Csűrös
 *
 */
public class Gradient extends Posteriors implements Count.UsesThreadpool
{
	public Gradient(TreeWithRates rates, ProfileTable table)
	{
		super(rates, table);
		if (table instanceof UniqueProfileTable)
		{
			this.utable = (UniqueProfileTable) table; 
		}else
		{
			this.utable = null;
		}
		this.min_copies = Integer.min(2,table.minCopies());
	}
	
	public Gradient(Likelihood factory)
	{
		super(factory);
		if (factory.table instanceof UniqueProfileTable)
		{
			this.utable = (UniqueProfileTable) factory.table; 
		}else
		{
			this.utable = null;
		}
		this.min_copies = Integer.min(2,factory.table.minCopies());
	}
	
	/**
	 * Same as the instantiating table, if it 
	 * was a UniqueProfileTable; or else null. 
	 */
	private UniqueProfileTable utable;

	/**
	 * Observation bias : minimum number of observed copies in the families. 
	 */
	private int min_copies = 1;

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
	
	protected int getMinimumObservedCopies()
	{
		return min_copies;
	}
	
	public void computeParameters()
	{
		factory.computeParameters();
		this.cachedLL = 0.0;
	}

	@Override
	public void setCalculationWidthThresholds(int absolute, double relative)
	{
		super.setCalculationWidthThresholds(absolute, relative);
		this.cachedLL = 0.0;
	}
	
	
	private static int DEFAULT_TRUNCATE_ABSOLUTE = 6;
	private static double DEFAULT_TRUNCATE_RELATIVE = 1.0;
	/**
	 * Tries to increase/decrease absolute an relative calculation width parameters 
	 * for better estimation of true likelihood and gradient (with max 
	 * calculation width).
	 * 
	 * @param tol
	 * @return log-likelihood
	 */
	public double adjustCalculationWidth(double tol)
	{
		int absolute = getCalculationWidthAbsolute();
		double relative = getCalculationWidthRelative();
		
		double current_LL = getCorrectedLL();
		// need to increase?
		double step_size = Math.log(2.0)/3.0;
		
		
		boolean adjustCalculationWidth=false;
		
		int num_adjustments = 0;
		
		int dabs =0, drel = 0;
		
		final int maxiter = 8;
		boolean adjusted_in_iteration=true;
		while (adjusted_in_iteration && num_adjustments < maxiter)
		{
			adjusted_in_iteration = false;
			
			{ // try changing absolute
				
				int next_absolute = Integer.max(absolute+1,(int)Math.ceil(Math.exp(Math.log(absolute)+step_size)));
				this.setCalculationWidthThresholds(next_absolute, relative);
				
				double next_LL = getCorrectedLL();
				
				double next_delta = next_LL-current_LL;
				double next_rdiff = Math.abs(next_delta/current_LL); 
				
				if (tol < next_rdiff )
				{
					absolute = next_absolute;
					current_LL = next_LL;
					adjusted_in_iteration = adjustCalculationWidth = true;
					++num_adjustments;
					++dabs;
				} else if (DEFAULT_TRUNCATE_ABSOLUTE < absolute)
				{
					int prev_absolute = Integer.max(Integer.min((int)Math.ceil(Math.exp(Math.log(absolute)-step_size)), absolute-1), DEFAULT_TRUNCATE_ABSOLUTE);
					this.setCalculationWidthThresholds(prev_absolute, relative);
					
					double prev_LL = getCorrectedLL();
					double prev_delta = prev_LL-current_LL;
					double prev_rdiff = Math.abs(prev_delta/current_LL);
					
					if (prev_rdiff < tol)
					{
						absolute = prev_absolute;
						current_LL = prev_LL;
						adjusted_in_iteration = adjustCalculationWidth = true;
						++num_adjustments;
						--dabs;
					}
				}
				
			}
			{
				// try changing relative
				double next_rel = Math.exp(Math.log(relative)+step_size);
				this.setCalculationWidthThresholds(absolute, next_rel);
				
				double next_LL = getCorrectedLL();
				
				double next_delta = next_LL-current_LL;
				double next_rdiff = Math.abs(next_delta/current_LL); 
				
				if (tol < next_rdiff)
				{
					relative = next_rel;
					current_LL = next_LL;
					adjusted_in_iteration = adjustCalculationWidth = true;
					++num_adjustments;
					++drel;
				} else if (DEFAULT_TRUNCATE_RELATIVE < relative)
				{
					double prev_rel = Double.max(DEFAULT_TRUNCATE_RELATIVE, Math.exp(Math.log(relative)-step_size));
					this.setCalculationWidthThresholds(absolute, prev_rel);
					
					double prev_LL = getCorrectedLL();
					double prev_delta = prev_LL-current_LL;
					double prev_rdiff = Math.abs(prev_delta/current_LL);
					if (prev_rdiff < tol)
					{
						relative = prev_rel;
						current_LL = prev_LL;
						adjusted_in_iteration = adjustCalculationWidth = true;
						++num_adjustments;
						--drel;
					}
				}
			} 			
		} // while changes		
		this.setCalculationWidthThresholds(absolute, relative);
		System.out.println("#**G.aCW setting "+absolute+","+relative);		
		
		return current_LL;
		
	}
	
	
	/**
	 * 
	 * Multiplicity of a family. 
	 * 
	 * @param family
	 * @return 1 if instantiated with ProfileTable, or the multiplicity from UniqueProfileTable
	 */
	protected int getMultiplicity(int family)
	{
		return utable==null?1:utable.getMultiplicity(family);
	}
	
	/**
	 * Sum of family multiplicities
	 * 
	 * @return
	 */
	public int getTotalFamilyCount()
	{
		return utable==null?factory.table.getFamilyCount():utable.getTotalFamilyCount();
	}
	private static ForkJoinPool thread_pool=null;
	/**
	 * Initializad only once, if {@link Count#THREAD_PARALLELISM} is greater than 1.
	 * 
	 * @return
	 */
	protected synchronized static ForkJoinPool threadPool()
	{
		if (thread_pool == null && 1<Count.THREAD_PARALLELISM) // && Count.THREAD_UNIT_TASK<Integer.MAX_VALUE)
		{
//			System.out.println("#**G.threadPool init: "+Count.THREAD_PARALLELISM+" threads on "+Thread.currentThread());
			thread_pool = Count.threadPool(); // new ForkJoinPool(Count.THREAD_PARALLELISM);	
		}
		return thread_pool;
	}
	
//	/**
//	 * Same thread pool across different instantiations.
//	 * (Not likely to be an issue through GUI but 
//	 * in a big model-selection procedure, this class might be instantiated many times.)  
//	 */
//	private static final ForkJoinPool thread_pool;
//	static 
//	{
//		if (THREAD_PARALLELISM>1)
//			thread_pool = new ForkJoinPool(THREAD_PARALLELISM);			
//		else
//			thread_pool = null;
//	}
//	
//	protected ForkJoinPool threadPool()
//	{
//		return thread_pool;
//	}
	
	private double cachedLL=0.0;
	/**
	 * Calculation of the raw likelihood, 
	 * with multi-threading, if enabled.
	 * 
	 * @return log-likelihood (recognizing multplicities in the uniqueProfileTable if initialize so)
	 */
	public double getLL()
	{
		if (cachedLL != 0.0) return cachedLL;
		
		int nF = factory.table.getFamilyCount();
		final ForkJoinPool thread_pool = threadPool();
		final int unit_task = Count.unitTask(nF);
		
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
							Likelihood.Profile P = factory.getProfileLikelihood(f);
							double fLL = P.getLogLikelihood();
							LL += fLL * getMultiplicity(P.family_idx);
						}
						return LL;
					}
				} catch (Throwable t)
				{
					throw t;
				}
			}
		}
		
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
		
		
		if (Double.isNaN(LL)) // somebody will complain with an assertion
		{
			System.out.println("#*G.gL badrates "+NewickParser.printTree(factory.tree));
			for (int node=0; node<factory.tree.getNumNodes(); node++)
			{
				System.out.println("#*MLD.oDF badrates "+node+"\t"+factory.rates.toString(node)+"\t"+factory.toString(node));
			}
			
			LL = Double.MAX_VALUE/-8.0;
		}
		
		
		cachedLL = LL;
		
		
		
		return LL;
	}
	
	public double getUnobservedLL()
	{
		double unobservedLL = Double.NEGATIVE_INFINITY;
		if (min_copies>0)
		{
			unobservedLL = factory.getEmptyLL();
			if (min_copies>1)
			{
				assert (min_copies==2);
				unobservedLL = Logarithms.add(unobservedLL, factory.getSingletonLL());
			}
		}
		return unobservedLL;
	}

	/**
	 * Corrected log-likelihood, by minimum unobserved copies.
	 * Uses multi-threading, if enabled. 
	 * 
	 * @return
	 */
	public double getCorrectedLL()
	{
		double LL = getLL(); 
		// true number of families 
		int nF = utable==null
				?factory.table.getFamilyCount()
				:utable.getTotalFamilyCount();
		double L0 = getUnobservedLL();
		// LL-F*log(1-exp(L0))
		double p_not0  = -Math.expm1(L0); // ok with L0== -infinity
		
		
		double LLcorr = LL-nF*Math.log(p_not0);
//		System.out.println("#**G.gCLL uncorr "+LL+"\tunobs "+L0+"\tp "+p_not0+"\tLLcorr "+LLcorr);
		
		assert (!Double.isNaN(LLcorr));
		return LLcorr;
	}
	
	/**
	 * Variance of the LL (by the same model)
	 * on bootstrap samples (profiles drawn with replacement).  
	 *  
	 * @return
	 */
	public double bootstrapLLVariance()
	{
		int nF = factory.table.getFamilyCount();
		final int unit_task = Count.unitTask(nF);
		final ForkJoinPool thread_pool = threadPool();
//		if (thread_pool != null)
//		{
//			unit_task = THREAD_UNIT_TASK;
//		}
//		else
//		{
//			unit_task = nF; // do not fork
//		}
		
		class PartialV extends RecursiveTask<double[]>
		{
			/**
			 * First family index.
			 */
			private final int minF;
			/**
			 * Last family index, exclusive.
			 */
			private final int maxF;

			PartialV(int min, int max)
			{
				this.minF = min;
				this.maxF = max;
			}
			@Override
			protected double[] compute() 
			{
				try
				{
					if (maxF-minF>unit_task)
					{
						int medF = (minF+maxF)/2;
						PartialV left = new PartialV(minF, medF);
						PartialV right = new PartialV(medF, maxF);
						right.fork();
						double[] sumsL = left.compute();
						double[] sumsR = right.join();
						sumsL[0]+=sumsR[0];
						sumsL[1]+=sumsR[1];
						return sumsL;
					} else
					{
						double LL = 0.0;
						double L2 = 0.0;
						for (int f =minF; f<maxF; f++)
						{
							Likelihood.Profile P = factory.getProfileLikelihood(f);
							double fLL = P.getLogLikelihood();
							int m = getMultiplicity(P.family_idx);
							LL += fLL * m;
							L2 += fLL*fLL*m;
						}
						double[] sums = new double[2];
						sums[0]=LL;
						sums[1]=L2;
						return sums;
					}
				} catch (Throwable t)
				{
					throw t;
				}
			}
		}
		
		PartialV bigjob = new PartialV(0,nF);
		double[] sums;
		try
		{
			if (nF>unit_task)
			{
				sums = thread_pool.invoke(bigjob);
			} else
			{
				sums = bigjob.compute();
			}
		} catch (Throwable t)
		{
			// could be out of memory, or out of threads, or numerical error, or whatever else			
			throw new RuntimeException(t);  
		}
		int nT = getTotalFamilyCount();
		double avgL = sums[0]/nT;
		double avgL2 = sums[1]/nT;
		
		
		double var = avgL2 - sums[0]*sums[0]/(nT*(nT-1.0));
		
//		System.out.println("#**G.bLLV avgL "+avgL+"\tavgL2 "+avgL2+"\tsqrt "+Math.sqrt(avgL2)+"\tvar "+var+"\tnT "+nT);
		return nT*var;
	}
	
	
	
	
	
	Profile getGradient(int family_idx)
	{
		return new Profile(family_idx);
	}
	
	
	class Profile 
	{
		private Profile(int family_idx)
		{
			this.post = getPosteriors(family_idx);
		}
		final Posteriors.Profile post;
		
		
		
		/**
		 * Gradient of the uncorrected log-likelihood.
		 * Index is 3*<var>node</var>+<var>i</var> where <var>i</var>={@link GLDParameters#PARAMETER_DUPLICATION},
		 * {@link GLDParameters#PARAMETER_GAIN}, {@link GLDParameters#PARAMETER_LOSS}.
		 * 
		 * Numerically stable version that works with very small 
		 * parameter values too. 
		 * 
		 * @return gradient array
		 */
		protected double[] getSurvivalGradient()
		{
			int num_nodes = factory.tree.getNumNodes();
			double[] dL = new double[num_nodes*3];// return value

			for (int node=0; node<num_nodes; node++)
			{
				double log1_q = factory.getLogDuplicationComplement(node); // Math.log(q1); // Math.log1p(-q); //Math.log(1.0-q); // 
				double log_q = factory.getLogDuplicationParameter(node);
				
				double log_mean_Sv = post.getLogEdgeMean(node);
				double log_mean_Nv = post.getLogNodeMean(node);

				// II.a gain
				if (log_q==Double.NEGATIVE_INFINITY) 
				{
					// Poisson
					double r = factory.getGainParameter(node);
					if (r!=0.0)
					{
						double log_mean_Nv_Sv = post.getLogNodeIncrease(node);
						double log_r = Math.log(r);
						//dL[3*node+PARAMETER_GAIN] = Math.expm1(log_mean_Nv_Sv-log_r);
						
						double dLdg = Logarithms.ldiffValue(Logarithms.ldiff(log_mean_Nv_Sv-log_r, 0.0));
						dL[3*node+PARAMETER_GAIN] = dLdg;
					}
				} else
				{
					// Pólya
					double κ = factory.getGainParameter(node);
					double[] log_Nv_Sv = post.getLogNodePosteriorIncrease(node);
					double log_kappa = Math.log(κ);
					if (κ!=0.0)
					{
						double log_dLdk = Double.NEGATIVE_INFINITY;
						for (int i=0; i<log_Nv_Sv.length; i++)
						{
							double log_denom = (κ<i)
									?Math.log(i)+Math.log1p(κ/i)
									:log_kappa+Math.log1p(i/κ);
							
							log_dLdk = Logarithms.add(log_dLdk, log_Nv_Sv[i]-log_denom);
						}
						double loglog1_q = Logarithms.logitToLogLogComplement(log_q-log1_q);
						dL[3*node+PARAMETER_GAIN] 
								= Math.exp(log_dLdk) + log1_q;
						double dLdk = Logarithms.ldiffValue(Logarithms.ldiff(log_dLdk, loglog1_q));
						dL[3*node+PARAMETER_GAIN] = dLdk;
					}
					double log_mean_Nv_Sv = Double.NEGATIVE_INFINITY;
					for (int i=0; i<log_Nv_Sv.length; i++)
					{
						log_mean_Nv_Sv = Logarithms.add(log_mean_Nv_Sv, log_Nv_Sv[i]);
					}
					// DEBUG
					
//					dL[3*node+PARAMETER_DUPLICATION]
//							= Math.exp(log_mean_Nv_Sv-log_q)
//								-Math.exp(Logarithms.add(log_mean_Sv, log_kappa)-log1_q);
					
					//log_mean_Nv_Sv = Math.log(Double.max(0.0,Math.exp(log_mean_Nv)-Math.exp(log_mean_Sv)));
					
					double[] log_dLdq = Logarithms.ldiff(log_mean_Nv_Sv-log_q, Logarithms.add(log_mean_Sv, log_kappa)-log1_q);
					double dLdq = Logarithms.ldiffValue(log_dLdq);
					dL[3*node+PARAMETER_DUPLICATION] = dLdq;
					
//					double diff_means = Math.exp(log_mean_Nv_Sv)-((Math.exp(log_mean_Nv)-Math.exp(log_mean_Sv)));
//					double err_means = Math.abs(diff_means/Math.exp(log_mean_Nv_Sv));
//					
//					if (dLdq == 0.0 && !Logarithms.ldiffIsZero(log_dLdq))
////							|| err_means>1.0)
//					{
//						int uniqf = post.inside.family_idx;
//						System.out.println("#**G.P.gSG "
//								+uniqf+"("+utable.getLineageCount(uniqf)+"/"+utable.getMemberCount(uniqf)+")"
//								+"\t"+node
//								+"\t"+err_means+"/"+diff_means
//								+"\tNvSv "+log_mean_Nv_Sv+"/"+Math.exp(log_mean_Nv_Sv)
//									+"\tSv "+log_mean_Sv+"/"+Math.exp(log_mean_Sv)
//									+"\tNv "+log_mean_Nv+"/"+Math.exp(log_mean_Nv)
//									+"\tNv-Sv "+(Math.exp(log_mean_Nv)-Math.exp(log_mean_Sv))
//									+"\tNvSv/q "+Math.exp(log_mean_Nv_Sv-log_q)
//									+"\tSv/1_q "+Math.exp(log_mean_Sv-log1_q)
//									+"\tkp/1_q "+Math.exp(log_kappa-log1_q)
//								+"\tlq "+log_q+"/"+Math.exp(log_q)
//								+"\tl1q "+log1_q+"/"+Math.exp(log1_q)
//								+"\tlkappa "+log_kappa
//								+"\tdL "+dL[3*node+PARAMETER_DUPLICATION]
//								+"\tdLlogit "+dL[3*node+PARAMETER_DUPLICATION]
//										*Math.exp(log_q+log1_q));
//					}
				}
				
				double logp = factory.getLogLossParameter(node);
				double log1_p = factory.getLogLossComplement(node);
				if (log1_p!=Double.NEGATIVE_INFINITY)   // (logp!=0.0)
				{
					double log_mean_Nu_Sv = post.getLogEdgeDecrease(node); 
					
					
					if (factory.tree.isRoot(node)) // untested model setting
					{
//						dL[3*node+PARAMETER_LOSS] =
//								Math.exp(log_mean_Nu_Sv-logp)-Math.exp(log_mean_Sv-log1_p);
						double dLdp =  Logarithms.ldiffValue(Logarithms.ldiff(log_mean_Nu_Sv-logp, log_mean_Sv-log1_p));
						dL[3*node+PARAMETER_LOSS] = dLdp;
						throw new UnsupportedOperationException("Root with loss<1 was never tested.");
					} else
					{
						int parent = factory.tree.getParent(node);
						
						// do not enforce high precision
						// log_mean_Nu_Sv = Math.log(Double.max(0.0,Math.exp(post.getLogNodeMean(parent)-Math.exp(log_mean_Sv))));
//						double log_szer2; 
//						if (factory.tree.getNumChildren(parent)==2)
//						{
//							log_szer2 = factory.getLogLossComplement(factory.tree.getSibling(node));
//						} else
//						{
//							double epsi = factory.getExtinction(parent);
//							log_szer2 = Math.log1p(-Math.exp(factory.getLogExtinction(parent)-logp));
//						}			
						double log1_e;
						if (factory.tree.getNumChildren(parent)==2)
						{
							log1_e = factory.getLogLossComplement(factory.tree.getSibling(node));
						} else
						{
							log1_e = Logarithms.logToLogComplement(factory.getLogExtinction(node)-logp);
						}			
						
						
						double lt1 = log_mean_Nu_Sv-logp;
						//double lt2 = log_mean_Sv+log_szer2-log1_p;
						double lt2 = log_mean_Sv+log1_e-log1_p;
						
						// fixing a math error in the formulas: a mising 1-p*epsilon denominator
						// 1-pe = 1-p + p*(1-e)
						double log1_pe = Logarithms.add(log1_p, logp+log1_e);
						lt1-=log1_pe;
						lt2-=log1_pe;
						
//						dL[3*node+PARAMETER_LOSS] =
//								Math.exp(lt1)
//								-Math.exp(lt2);	
						
						double[] log_dLdp = Logarithms.ldiff(lt1, lt2);
						double dLdp = Logarithms.ldiffValue(log_dLdp);
						dL[3*node+PARAMETER_LOSS] = dLdp;
						
//						if (dLdp == 0.0 && !Logarithms.ldiffIsZero(log_dLdp))
//						{
//							System.out.println("#**G.P.gSG "
//							+post.inside.family_idx
//							+"\t"+node
//							+"\tNuSv "+log_mean_Nu_Sv+"\tSv "+log_mean_Sv
//							+"\tlp "+logp+"\tl1p "+log1_p+"\tl1e "+log1_e
//							+"\t"+Arrays.toString(log_dLdp)
//							+"\tdL "+dL[3*node+PARAMETER_LOSS]
//							+"\tdLlogit "+dL[3*node+PARAMETER_LOSS]
//									*Math.exp(logp+log1_p));
//							
//						}
						
					}
				}
				
			} // for node
			return dL;
		}

		/**
		 * Gradient by logit p~, logit q~, log kappa/log r~
		 * 
		 * @return
		 */
		protected double[] getLogitSurvivalGradient()
		{
			int num_nodes = factory.tree.getNumNodes();
			double[] dL = new double[num_nodes*3];// return value

			for (int v=0; v<num_nodes; v++) // for all nodes v 
			{		
				double log_Sv = post.getLogEdgeMean(v);
				/*
				 * derivative by logit p
				 */
				int jloss = 3*v+PARAMETER_LOSS; // parameter index
				double log_Nu_Sv = post.getLogEdgeDecrease(v); 
				double log_p = factory.getLogLossParameter(v);
				double log1_p = factory.getLogLossComplement(v);
				if (factory.tree.isRoot(v) || log1_p == Double.NEGATIVE_INFINITY) // p==1.0 
				{
					dL[jloss] = 0.0;
				} else
				{
					int u = factory.tree.getParent(v);
					double log1_e; // non-extinction at v's siblings
					if (factory.tree.getNumChildren(u)==2)
					{
						log1_e = factory.getLogLossComplement(factory.tree.getSibling(v));
					} else
					{
						log1_e = Logarithms.logToLogComplement(factory.getLogExtinction(u)-log_p);
					}		
					double log1_eu = factory.getLogExtinctionComplement(u);
					double dpos = log_Nu_Sv+log1_p-log1_eu;
					double dneg = log_Sv+log_p+log1_e-log1_eu;
					dL[jloss]  = Logarithms.ldiffValue(Logarithms.ldiff(dpos, dneg));
				}
				/*
				 * derivative by log-gain, and logit-duplication
				 */
				int jgain = 3*v+PARAMETER_GAIN;
				int jdup = 3*v+PARAMETER_DUPLICATION;
				double log1_q = factory.getLogDuplicationComplement(v);  
				double log_q = factory.getLogDuplicationParameter(v);
				double log_Nv_Sv = post.getLogNodeIncrease(v);
				if (log_q==Double.NEGATIVE_INFINITY)
				{
					// Poisson model
					double log_r = factory.getLogGainParameter(v);
					
					dL[jgain] = Logarithms.ldiffValue(Logarithms.ldiff(log_Nv_Sv, log_r));
					dL[jdup] = 0.0;
				} else
				{
					assert (Double.isFinite(log_q));
					// Polya model
					double log_kappa = factory.getLogGainParameter(v);
					double κ = Math.exp(log_kappa);
					double[] tNv_Sv = post.getLogNodePosteriorIncrease(v); // difference of tail probabilities					
					double dpos = tNv_Sv[0]; // first term
					for (int i=1; i<tNv_Sv.length; i++)
					{
						// ln(k/(k+i)) = ln(1/(1+i/k)) =-ln(1+i/k) for i<k
						// = ln (k/i)-ln(1+k/i) for k<i
						double log_k_ki;
						if (i<κ)
						{
							log_k_ki = -Math.log1p(i/κ);
						} else
						{
							log_k_ki = log_kappa - Math.log(i) -Math.log1p(κ/i);
						}
						dpos = Logarithms.add(dpos, tNv_Sv[i]+log_k_ki);
					}
					// gain
					double loglog1_q = Logarithms.logitToLogLogComplement(log_q-log1_q);	
					double dneg = log_kappa + loglog1_q;
					dL[jgain] = Logarithms.ldiffValue(Logarithms.ldiff(dpos,dneg));
					// duplication
					dpos = log_Nv_Sv+log1_q;
					dneg = Logarithms.add(log_Sv, log_kappa) + log_q;			
					dL[jdup] = Logarithms.ldiffValue(Logarithms.ldiff(dpos, dneg));
				}
			} // for all nodes v
			return dL;
			
		}
		
		
		/**
		 * Gradient of the uncorrected log-likelihood.
		 * Index is 3*<var>node</var>+<var>i</var> where <var>i</var>={@link GLDParameters#PARAMETER_DUPLICATION},
		 * {@link GLDParameters#PARAMETER_GAIN}, {@link GLDParameters#PARAMETER_LOSS}.
		 * 
		 * Derivative by p is not correct. 
		 * 
		 * @return gradient array
		 */
		double[] oldSurvivalGradient()
		{
			// I. compute the posterior means and tails
			double[] node_means = post.getNodeMeans(); // array of expected values across the nodes
			double[] edge_means = post.getEdgeMeans(); 
			int num_nodes = factory.tree.getNumNodes();
			double[][] node_tails = new double[num_nodes][];
			double[][] edge_tails = new double[num_nodes][];		
			
			for (int node=0; node<num_nodes; node++)
			{
				node_tails[node] = post.getNodeTail(node);
				edge_tails[node] = post.getEdgeTail(node);
//				double[] Ncdf = post.getNodeCDF(node);
//				double[] Ntail 
//					= node_tails[node] = new double[Ncdf.length];
//				for (int ell=0; ell<Ncdf.length-1; ell++) // last entry is 1.0
//					Ntail[ell]=1.0-Ncdf[ell];
//				double[] Scdf = post.getEdgeCDF(node);
//				double[] Stail 
//					= edge_tails[node] = new double[Ntail.length];
//				
//				assert (Scdf.length <= Stail.length); 
//				
//				for (int s=0; s<Scdf.length-1; s++)
//					Stail[s]=1.0-Scdf[s];
			}
			
			
			// II. compute the gradients 
			double[] dL = new double[num_nodes*3];// return value
			
			for (int node=0; node<num_nodes; node++)
			{
				double q = factory.getDuplicationParameter(node);
				double q1 = factory.getDuplicationParameterComplement(node);
				if (q==0.0)
				{
					// Poisson
					double r = factory.getGainParameter(node);
					if (r!=0.0)
					{
						dL[3*node+PARAMETER_GAIN] 
								= (node_means[node] - edge_means[node])/r - 1.0; // d/dr
					}
				} else
				{
					// Pólya
					double κ = factory.getGainParameter(node);
					if (κ!=0.0)
					{
						double dLdk = 0.0;
						double[] Ntail = node_tails[node];
						double[] Stail = edge_tails[node];
						assert (Stail.length<=Ntail.length);
						{
							int i; 
							for (i=0; i<Stail.length; i++)
							{
								double Nu_Su = Double.max(0.0, Ntail[i]-Stail[i]);
								dLdk += Nu_Su/(κ + i);
							}
							for (; i<Ntail.length; i++)
							{
								double Nu_Su = Ntail[i];
								dLdk += Nu_Su/(κ + i);
							}
						}
						double log1_q = Math.log(q1); // Math.log1p(-q); //Math.log(1.0-q); // 
						dL[3*node+PARAMETER_GAIN] 
								= dLdk + log1_q;// d/dκ
						
						assert !Double.isNaN(dL[3*node+PARAMETER_GAIN] );
					}
					
					
					dL[3*node+PARAMETER_DUPLICATION]
							= (node_means[node]-edge_means[node])/q
							- (edge_means[node]+ κ)/q1; // d/dq
				}
				double p = factory.getLossParameter(node);
				double p1 = factory.getLossParameterComplement(node);
				if (p!=1.0)
				{
					if (factory.tree.isRoot(node)) // untested model setting
					{
						dL[3*node+PARAMETER_LOSS]
								= (1.0-edge_means[node])/p - edge_means[node]/p1;
//						if (! Double.isFinite(dL[3*node+PARAMETER_LOSS] ))
//						{
//							System.out.println("#**G.P.gSG dloss "+dL[3*node+PARAMETER_LOSS]+"\tSn "+edge_means[node]+"\tp "+p+"\t// "+factory.rates.toString(node));
//						}
					} else
					{
						int parent = factory.tree.getParent(node);
						double omeu =  factory.getExtinctionComplement(parent); // 1.0 - factory.extinction[parent];
						// double epsi = factory.extinction[parent]/p;
						dL[3*node+PARAMETER_LOSS]
								= (node_means[parent]-omeu*edge_means[node])/p - omeu*edge_means[node]/p1;
//						if (! Double.isFinite(dL[3*node+PARAMETER_LOSS] ))
//						{
//							System.out.println("#**G.P.gSG dloss "+dL[3*node+PARAMETER_LOSS]+"\tNp "+node_means[parent]+"\tSn "+edge_means[node]+"\tome "+omeu+"\tp "+p+"\t// "+factory.rates.toString(node)
//							+"\t"+factory.tree.toString(node));
//						}
					}
				} else // with p==1.0, keep dLdp=0
				{
				}
				assert Double.isFinite(dL[3*node+PARAMETER_DUPLICATION] );
				assert Double.isFinite(dL[3*node+PARAMETER_LOSS] );
				assert Double.isFinite(dL[3*node+PARAMETER_GAIN] );
			}
			
//			System.out.println("#**G.P.gSG "+Arrays.toString(dL));
			
			return dL;
		}		
	}
	

	public double[] getEmptySurvivalGradient()
	{
		int num_nodes = factory.tree.getNumNodes();
		double[] dL = new double[num_nodes*3];// return value
		
		for (int node=0; node<num_nodes; node++)
		{
			//double q = factory.getDuplicationParameter(node);
			double log1_q = factory.getLogDuplicationComplement(node) ;
			double log_q = factory.getLogDuplicationParameter(node);
			
			if  (log_q == Double.NEGATIVE_INFINITY) // (log1_q==0.0)
			{
				// Poisson
//				double r = factory.getGainParameter(node);
//				if (r!=0.0)
//				{
					dL[3*node+PARAMETER_GAIN] 
							=  - 1.0; // d/dr
//				}
			} else
			{
				// Pólya
				double κ = factory.getGainParameter(node);
				//double log1_q = factory.getLogDuplicationComplement(node) ; // Math.log1p(-q); // Math.log(1.0-q); // Mat.log1p(-q)
//				if (κ!=0.0)
//				{
					dL[3*node+PARAMETER_GAIN] 
							= log1_q;// d/dκ
//				}
				double log_kappa = factory.getLogGainParameter(node);
				
				dL[3*node+PARAMETER_DUPLICATION]
						= - Math.exp(log_kappa-log1_q);//   κ/(1.0-q); // d/dq
			}
		}	
		int root = num_nodes-1; assert (root == factory.tree.getRoot());
		double proot = factory.getLossParameter(root);
		if (proot != 1.0) // proot==1.0 does not change
		{
			dL[3*root+PARAMETER_LOSS] = 1.0/proot;
		}
		return dL;
	}
	
	private ProfileTable singletons = null;
	
	/**
	 * A table of singleton profiles; instantiated here if not before.  
	 * 
	 * @return always non-null
	 */
	protected ProfileTable singletons() 
	{ 
		if (singletons == null)
			singletons = ProfileTable.singletonTable(factory.tree);
		return singletons;
	}
	
	
	public double[] getSingletonSurvivalGradient()
	{	/*
		 d log(p1+p2+...) / dx 
		 = (d p1/dx + d p2/dx +...)/(p1+p2+...)
		 
		 and d pj / dx = pj * d log pj/dx 
		*/
		if (singletons == null)
			singletons = ProfileTable.singletonTable(factory.tree);
		Gradient S = new Gradient(factory.rates, singletons);
		int num_nodes = factory.tree.getNumNodes();
		double[] dLs = new double[num_nodes*3];
		double log_ps = factory.getSingletonLL();
		for (int f=0; f<singletons.getFamilyCount(); f++)
		{
			Profile P = S.getGradient(f);
			Likelihood.Profile LP = S.factory.getProfileLikelihood(f);
			double log_pf = LP.getLogLikelihood();
			double wf = Math.exp(log_pf-log_ps);
			
			double[] dLf = P.getSurvivalGradient();
			
			assert (dLf.length == dLs.length);
			for (int i=0; i<dLs.length; i++)
			{
				dLs[i] += dLf[i] * wf;
			}
		}
		return dLs;
	}
	
	public double[] getUnobservedSurvivalGradient()
	{
		double[] survival_gradient ;
		if (min_copies == 0)
		{
			int num_nodes = factory.tree.getNumNodes();
			survival_gradient = new double[3*num_nodes];
		} else
		{
			survival_gradient = getEmptySurvivalGradient();
			if (min_copies>1)
			{
				double[] singleton_gradient = getSingletonSurvivalGradient();
				assert (min_copies==2);
				double L0 = factory.getEmptyLL();
				double L1 = factory.getSingletonLL();
				double unobsL = Logarithms.add(L0, L1);
				double w0 = Math.exp(L0-unobsL); 
				double w1 = Math.exp(L1-unobsL);

				for (int i=0; i<survival_gradient.length; i++)
				{
					survival_gradient[i]
							= survival_gradient[i]*w0
							+ singleton_gradient[i]*w1;
				}
			}
		}
		return survival_gradient;
	}
	
	/**
	 * Calculates gradient by logit-p, logit-q, log-kappa/log-r; 
	 * replaces entries of the input array.
	 * 
	 * @param survival_gradient from {@link Profile#getLogitSurvivalGradient()} 
	 * @return same array, with updated entries
	 */
	public double[] convertToLogitDistributionGradient(double[] survival_gradient)
	{
		double[] dL = survival_gradient; // return value
		int num_nodes = factory.tree.getNumNodes();
		double[] dde = new double[num_nodes]; // d/dε
		
		int v=factory.tree.getRoot();
		while (0<=v)
		{
			// parameter indices
			int jloss = 3*v+PARAMETER_LOSS;
			int jdup  = 3*v+PARAMETER_DUPLICATION;
			int jgain = 3*v+PARAMETER_GAIN;
			
			double log_e = factory.getLogExtinction(v);
			if (factory.tree.isRoot(v) || factory.rates.getLogLossComplement(v) ==Double.NEGATIVE_INFINITY)
			{
				if (factory.rates.getLogDuplicationParameter(v)==Double.NEGATIVE_INFINITY) // (q==0.0)			
				{
					// Poisson
					dde[v] = -Math.exp(log_e)*dL[jgain];
					assert (dL[jloss]==0.0);
					assert (dL[jdup] ==0.0);
				} else
				{
					// Polya
					dde[v] = -Math.exp(log_e)*dL[jdup];
				}
			} else
			{
				int u=factory.tree.getParent(v);

				double log_p = factory.getLogLossParameter(v);
				double log1_p = factory.getLogLossComplement(v);
				
				double log_dedp = log1_p-factory.getLogExtinctionComplement(u); 
				double dLdpe = dL[jloss] + Math.exp(log_dedp)*dde[u];
				
				double log_dpdp = factory.rates.getLogLossParameter(v)-log_p; 
				dL[jloss] = dL[jloss] * Math.exp(log_dpdp);
				
				if (factory.rates.getLogDuplicationParameter(v)==Double.NEGATIVE_INFINITY) // (q==0.0)
				{
					// Poisson
					double log_dpde = log_e-log_p;
					dde[v] = Math.exp(log_dpde) * dLdpe - Math.exp(log_e)*dL[jgain];
					assert (dL[jdup]==0.0);
				} else
				{
					// Polya
					double log_dpde = factory.getLogDuplicationComplement(v)+log_e-log_p;
					double log_dpdq = log_dpde + factory.rates.getLogDuplicationParameter(v);
					
					double dLdq = dL[jdup];
					dL[jdup] = dLdq - Math.exp(log_dpdq)*dLdpe;
					dde[v] = Math.exp(log_dpde)*dLdpe - Math.exp(log_e)*dLdq;
				}
			}
			
			--v;
		} // for all nodes in preorder	
		return dL;
	}
	
	/**
	 * Calculates the gradient by distribution parameters from the gradient by survival parameters. 
	 * 
	 * @param survival_gradient content is rewritten and returned 
	 * @return same as input array, entries updated 
	 */
	public double[] getDistributionGradient(double[] survival_gradient)
	{
		final double[] dL = survival_gradient; 
		int num_nodes = factory.tree.getNumNodes();
		double[] de = new double[num_nodes]; // d/dε
		
		int node = factory.tree.getRoot();
		// init at the root 
		double q = factory.rates.getDuplicationParameter(node);
		double log_q = factory.rates.getLogDuplicationParameter(node);
		if (log_q==Double.NEGATIVE_INFINITY) // (q==0.0)
		{
			// Poisson
			double r = factory.rates.getGainParameter(node);
			if (r!=0.0)
			{
				double dLdr =  dL[3*node+PARAMETER_GAIN];
				dL[3*node+PARAMETER_GAIN] = dLdr*factory.getExtinctionComplement(node); // (1.0-factory.extinction[node]); 
				de[node] = dLdr*(-r);
			}
		} else
		{
			// Pólya or geometric
			
			double q_1 = factory.rates.getDuplicationParameterComplement(node);
			double dLdq = dL[3*node+PARAMETER_DUPLICATION];
			double e_1 = factory.getExtinctionComplement(node);
			double a = q_1 + q*e_1; // 1-qe = 1-q + q(1-e)
					   //1.0-q*factory.getExtinction(node);
			double a2 = a*a;
			dL[3*node+PARAMETER_DUPLICATION] = dLdq * e_1 //  (1.0-factory.extinction[node])
												/a2;
			// dL[3*node+PARAMETER_GAIN] does not change
			de[node] = -dLdq*q_1*q/a2;
		}
		assert Double.isFinite(dL[3*node+PARAMETER_DUPLICATION] );
		assert Double.isFinite(dL[3*node+PARAMETER_LOSS] );
		assert Double.isFinite(dL[3*node+PARAMETER_GAIN] );
//		{ // DEBUG
//			System.out.println("#**G.gDG "+node
//					+"\tdLdlp "+dL[3*node+PARAMETER_LOSS]*factory.rates.getLossParameter(node)*factory.rates.getLossParameterComplement(node)
//					+"\tdLdlq "+dL[3*node+PARAMETER_DUPLICATION]*factory.rates.getDuplicationParameter(node)*factory.rates.getDuplicationParameterComplement(node)
//					+"\tdLdlg "+dL[3*node+PARAMETER_GAIN]*factory.rates.getGainParameter(node)
//					+"\tdLdle "+de[node]*factory.getExtinction(node)*factory.getExtinctionComplement(node)
//				);
//		}
		while (node>0)
		{
			--node;
			int parent = factory.tree.getParent(node);
			q = factory.rates.getDuplicationParameter(node);
			double p = factory.rates.getLossParameter(node);
			double ε = 
				    factory.tree.getNumChildren(parent)==2
				    ?factory.getLossParameter(factory.tree.getSibling(node))
				    :factory.getExtinction(parent)/factory.getLossParameter(node);
			double dLdpe = dL[3*node+PARAMETER_LOSS]+ε*de[parent];
			double epsi1 = factory.getExtinctionComplement(node);
			if (q==0.0)
			{
				// Poisson
				double r = factory.rates.getGainParameter(node);
				if (r!=0.0)
				{
					double dLdr =  dL[3*node+PARAMETER_GAIN];
					dL[3*node+PARAMETER_GAIN] = dLdr *  epsi1; // (1.0-factory.extinction[node]); 
					dL[3*node+PARAMETER_LOSS] = dLdpe * epsi1; //(1.0-factory.extinction[node]);
					de[node] = dLdpe * factory.rates.getLossParameterComplement(node) - dLdr * r;
				}
				
			} else
			{
				// Pólya or geometric
				double q_1 = factory.rates.getDuplicationParameterComplement(node);
				double p_1 = factory.rates.getLossParameterComplement(node); 
				//double a = 1.0-q*factory.getExtinction(node);
				double dp_dp = Math.exp(factory.getLogLossComplement(node)-factory.rates.getLogLossComplement(node));
				
				double dLdp = dL[3*node+PARAMETER_LOSS];
				dL[3*node+PARAMETER_LOSS] = dLdpe *dp_dp;  // epsi1/a; // (1.0-factory.extinction[node]) 
//				System.out.println("#**G.gDG   "+node+"\tdpdp "+dp_dp+"\tdLdpe "+dLdpe
//						+"\tdLdp~ "+dLdp
//						+"\tdeu "+de[parent]+"\tepsi "+ε
//						+"\tdLdp "+dL[3*node+PARAMETER_LOSS]
//					);
											
				double dLdq = dL[3*node+PARAMETER_DUPLICATION];
				//double a2 = a*a;
				
				double dq_dq = Math.exp(factory.getLogDuplicationParameter(node)+factory.getLogDuplicationComplement(node)
						-factory.rates.getLogDuplicationParameter(node)-factory.rates.getLogDuplicationComplement(node));


				
				dL[3*node+PARAMETER_DUPLICATION] = (dLdq - dLdpe * p_1 * factory.getExtinction(node)) 
							*dq_dq;
							// * epsi1 //(1.0-factory.extinction[node])
							// /a2;
				
//				System.out.println("#**G.gDG   "+node+"\tdqdq "+dq_dq+"\tdLdq~ "+dLdq+"\tdLdpe "+dLdpe
//							+"\tp1e "+ (p_1 * factory.getExtinction(node))
//							+"\tdLdq "+dL[3*node+PARAMETER_DUPLICATION]);
				
				double d_de = Math.exp(2.0*factory.getLogDuplicationComplement(node)-factory.rates.getLogDuplicationComplement(node));
				//System.out.println("#***G.gDG "+node+"\td_de "+d_de+"\t/"+(q_1/a2));
				de[node] = (dLdpe * p_1 - dLdq * q) * d_de ; // q_1/a2;

//				System.out.println("#**G.gDG   "+node+"\td_de "+d_de+"\tdLdq~ "+dLdq+"\tdLdpe "+dLdpe
//						+"\tdLde "+de[node]);
			}
			assert Double.isFinite(dL[3*node+PARAMETER_DUPLICATION] );
			assert Double.isFinite(dL[3*node+PARAMETER_LOSS] );
			assert Double.isFinite(dL[3*node+PARAMETER_GAIN] );
			
//			{ // DEBUG
//				System.out.println("#**G.gDG   "+node
//						+"\tdLdlp "+dL[3*node+PARAMETER_LOSS]*factory.rates.getLossParameter(node)*factory.rates.getLossParameterComplement(node)
//						+"\tdLdlq "+dL[3*node+PARAMETER_DUPLICATION]*factory.rates.getDuplicationParameter(node)*factory.rates.getDuplicationParameterComplement(node)
//						+"\tdLdlg "+dL[3*node+PARAMETER_GAIN]*factory.rates.getGainParameter(node)
//						+"\tdLdle "+de[node]*factory.getExtinction(node)*factory.getExtinctionComplement(node)
//					);
//			}
		}
		return dL;
	}
	
	/**
	 * 
	 * @param node
	 * @param dL distribution-parameter gradient
	 * @return
	 */
	public double inferDuplicationRateGradient(int node, double[] dL)
	{
		// TODO use base parameter's values
		
		final TreeWithRates rates = factory.rates;

		double μ = rates.getLossRate(node);
		double λ = rates.getDuplicationRate(node);
		double t = rates.getEdgeLength(node);	
		
		double dLdp = dL[3*node+PARAMETER_LOSS] ;
		double dLdq = dL[3*node+PARAMETER_DUPLICATION] ;
		
		double dLdλ; // return value
		
		if (Double.isInfinite(t))
		{
			// q = λ/μ 
			double dqdλ = 1.0/μ;
			dLdλ = dLdq * dqdλ;
		} else
		{
			double μt = μ*t;
			double λt = λ*t;
			if (μ == λ) // symmetric rates p=q=mu *t/(1+mu*t); should not be needed since p can be optimized directly with fixed λ=1
			{
				double denom = 1.0+μt;
				double divby = denom*denom;
				
				double dpdλ = 1.0/divby;
				double dqdλ = dpdλ;
				
				dLdλ = //dLdp*dpdλ + 
							dLdq*dqdλ;  // derivative by (λt) ??
			} else 
			{ 
				double gap = rates.getRateGap(node);
				double d = gap*t; //(μ-λ)*t;
				double E = Math.exp(-d);
				double E1 = -Math.expm1(-d); // 1.0-E;
				double delta = gap/μ;
				double denom;
				if (delta<0.5)
				{
					denom = μt * (-Math.expm1(-d+Math.log1p(-delta)));
				} else
				{
					denom = μt-λt*E; // = μt-λt*E;      // and not d*E = (μt-λt)*E;
				}
				double divby = denom*denom;
				
				double dpdλ = E*μt*(E1-d)/divby;
				double dqdλ = (E1*μt-λt*d*E)/divby;
				dLdλ = dLdp*dpdλ 
						 + dLdq*dqdλ
						;	// derivative by (λt)
			}	
			dLdλ *= t;
		}
		
		return dLdλ; 
	}
	
//	public double inferRateGradient(int  node, int param_type, double[] distribution_gradient)
//	{
//		final double[] dL = distribution_gradient; // shorter name
//		TreeWithRates rates = factory.rates;
//		double dLdr = dL[3*node+PARAMETER_GAIN] ;
//		double dLdp = dL[3*node+PARAMETER_LOSS] ;
//		double dLdq = dL[3*node+PARAMETER_DUPLICATION] ;
//		
//		
//		double dLdθ ; // return value
//
//		double μ = rates.getLossRate(node);
//		double λ = rates.getDuplicationRate(node);
//		double t = rates.getEdgeLength(node);	
//		
//		if (param_type == PARAMETER_GAIN)
//		{
//			if (Double.isInfinite(t) || factory.tree.isRoot(node))
//			{
//				// Pólya or Poisson
//				// nothing to do : r=γ for Poisson, or κ for Pólya
//				dLdθ = dLdr;
//			} else
//			{
//				if (λ==0.0 && rates.getGainParameter(node)!=0.0) // Poisson
//				{
//					double mu_t = μ*t;
//					double drdγ = -Math.expm1(-mu_t); // 1.0-dpdm;
//					dLdθ = dLdr*drdγ;
//				} else // Polya 
//				{
//					dLdθ = dLdr;
//				}
//			}
//		} else if (Double.isInfinite(t))
//		{
//			assert (dLdp==0.0);
//			// but q = λ/μ
//			if (param_type == PARAMETER_LOSS)
//			{
//				double dqdμ = -λ/(μ*μ);
//				dLdθ =  dLdq * dqdμ;
//			}
//			else
//			{
//				assert (param_type == PARAMETER_DUPLICATION);
//				double dqdλ = 1.0/μ;
//				dLdθ = dLdq * dqdλ;
//			}
//		} else if (factory.tree.isRoot(node)) // and finite t 
//		{
//			// shifted geometric
//
//			double μt = μ*t;
//			double λt = λ*t;
//			double d = (μ-λ)*t;
//			double E = Math.exp(-d);
//			double E1 = -Math.expm1(-d); // 1.0-E;
//			
//			double denom = μt-λt*E;
//			double divby = denom*denom;
//			
//			if (param_type == PARAMETER_LOSS)
//			{
//				double dpdμ = E*(μt*d-λt*E1)/divby;
//				double dqdμ = λt*(E*d-E1)/divby;
//				dLdθ = dLdp*dpdμ + dLdq*dqdμ;
//			} else
//			{
//				assert (param_type == PARAMETER_DUPLICATION);
//				double dpdλ = E*μt*(E1-d)/divby;
//				double dqdλ = (E1*μt-λt*d*E)/divby;
//				dLdθ = dLdp*dpdλ + dLdq*dqdλ;
//			}
//		} else // inner node and finite t 
//		{ 
//			if (λ==0.0 && rates.getGainParameter(node)!=0.0) // Poisson
//			{ // Poisson
//				if (param_type == PARAMETER_LOSS)
//				{
//					double mu_t = μ*t;
//					double dpdμ = Math.exp(-mu_t);
//					double drdμ = rates.getGainRate(node)*dpdμ;
//					dLdθ = dLdp * dpdμ + dLdr * drdμ;
//				} else
//				{
//					assert (param_type == PARAMETER_DUPLICATION);
//					dLdθ = 0.0;
//				}
//			} else
//			{
//				// Pólya
//				double μt = μ*t;
//				double λt = λ*t;
//				
//
//				if (μ == λ) // symmetric rates p=q=mu *t/(1+mu*t)
//				{
//					double denom = 1.0+μt;
//					double divby = denom*denom;
//					
//					double dpdμ = 1.0/divby;
//					double dqdμ = dpdμ;
//					
//					dLdθ = dLdp*dpdμ + dLdq*dqdμ; // same for both
//				} else 
//				{ 
//					double d = μt-λt;
//					double E = Math.exp(-d);
//					double E1 = -Math.expm1(-d); // 1.0-E;
//					
//					double denom = μt-λt*E;      // and not d*E = (μt-λt)*E; // BUG 4/13/2023
//					double divby = denom*denom;
//					
//					if (param_type == PARAMETER_LOSS)
//					{
//						double dpdμ = E*(μt*d-λt*E1)/divby;
//						double dqdμ = λt*(E*d-E1)/divby;
//						dLdθ = dLdp*dpdμ + dLdq*dqdμ;
//					} else
//					{
//						assert (param_type == PARAMETER_DUPLICATION);
//						double dpdλ = E*μt*(E1-d)/divby;
//						double dqdλ = (E1*μt-λt*d*E)/divby;
//						dLdθ = dLdp*dpdλ + dLdq*dqdλ;	
//					}
//				}
//			}
//		}
//		return dLdθ;
//	}
	
	/**
	 * Transforms distribution parameter gradient to (total) rate gradients.
	 * 
	 * @param distribution_gradient 3 cells per node 
	 * @return 3 cells per node 
	 */
	public double[] getTotalRateGradient(double[] distribution_gradient)
	{
		final double[] dL = distribution_gradient;
		int num_nodes = factory.tree.getNumNodes();
		
		TreeWithRates rates = factory.rates;
		int node = 0;
		while (node < num_nodes-1) // but not for the root
		{
			double t = rates.getEdgeLength(node);	
			double μ = rates.getLossRate(node);
			double λ = rates.getDuplicationRate(node);

			double dLdp = dL[3*node+PARAMETER_LOSS] ;
			double dLdq = dL[3*node+PARAMETER_DUPLICATION] ;

			if (Double.isInfinite(t))
			{
				assert (dLdp==0.0);
				// Pólya or Poisson
				// nothing to do : r=γ for Poisson, or κ for Pólya
				// but q = λ/μ
				double dqdλ = 1.0/μ;
				double dqdμ = -λ/(μ*μ);
				dL[3*node+PARAMETER_DUPLICATION] = dLdq * dqdλ;
				dL[3*node+PARAMETER_LOSS] = dLdq * dqdμ;
			}  else
			{
				double q = rates.getDuplicationParameter(node);
				if (q==0.0)
				{
					// Poisson
					double r = rates.getGainParameter(node);
	
					double mu_t = μ*t;
					double dpdμ = Math.exp(-mu_t);
					double drdμ = rates.getGainRate(node)*dpdμ;
					double drdγ = -Math.expm1(-mu_t); // 1.0-dpdm;
	
					if (r!=0.0)
					{
						double dLdr = dL[3*node+PARAMETER_GAIN] ;
						
						dL[3*node+PARAMETER_GAIN] = dLdr*drdγ;
						dL[3*node+PARAMETER_LOSS] = dLdp * dpdμ + dLdr * drdμ;
					} else
					{
						dL[3*node+PARAMETER_LOSS] = dLdp * dpdμ;
					}
				} else
				{
					// Pólya
					double μt = μ*t;
					double λt = λ*t;
					
	
					if (μ == λ) // symmetric rates p=q=mu *t/(1+mu*t)
					{
						double denom = 1.0+μt;
						double divby = denom*denom;
						
						double dpdμ = 1.0/divby;
						double dqdμ = dpdμ;
						
						dL[3*node+PARAMETER_LOSS] 
						= dL[3*node+PARAMETER_DUPLICATION]
						= dLdp*dpdμ + dLdq*dqdμ;
					} else 
					{ 
						double d = (μ-λ)*t;
						double E = Math.exp(-d);
						double E1 = -Math.expm1(-d); // 1.0-E;
						
						double denom = μt-λt*E;      // and not d*E = (μt-λt)*E; // BUG 4/13/2023
						double divby = denom*denom;
						
						
						double dpdμ = E*(μt*d-λt*E1)/divby;
						double dpdλ = E*μt*(E1-d)/divby;
						double dqdμ = λt*(E*d-E1)/divby;
						double dqdλ = (E1*μt-λt*d*E)/divby;
						
						dL[3*node+PARAMETER_LOSS] = dLdp*dpdμ + dLdq*dqdμ;
						dL[3*node+PARAMETER_DUPLICATION] = dLdp*dpdλ + dLdq*dqdλ;
	
					} // mu/lambda
				} // q>0.0
			} // finite edge length
			assert Double.isFinite(dL[3*node+PARAMETER_DUPLICATION] );
			assert Double.isFinite(dL[3*node+PARAMETER_LOSS] );
			assert Double.isFinite(dL[3*node+PARAMETER_GAIN] );
			
			++node;
		}			
//		if (!factory.tree.isRoot(node))
//		{
//			System.out.println("#**G.gRG expecting root @ "+node+"\trates "+rates.getTree()+"\tfactory "+factory.tree);
//			Phylogeny phylo = (Phylogeny) factory.tree;
//			Phylogeny.Node phylo_root = phylo.getRootNode();
//			System.out.println("#**G.gRG rindex "+factory.tree.getRoot()+"\trnode "+phylo_root);
//			System.out.print(phylo_root.prettySubtree());
//			for (int pnode=0; pnode<factory.tree.getNumNodes(); pnode++)
//			{
//				System.out.println("#**G.gRG node "+pnode+"\t"+factory.tree.toString(pnode));
//			}
//			
//		}
		assert (factory.tree.isRoot(node)); // no rate transformation for root
		{
			double dLdq = dL[3*node+PARAMETER_DUPLICATION] ;
			double dLdp = dL[3*node+PARAMETER_LOSS] ;
			double μ = rates.getLossRate(node);
			double λ = rates.getDuplicationRate(node);
			double t = rates.getEdgeLength(node);	
			if (Double.isInfinite(t))
			{
				assert (dLdp==0.0);
				// Pólya or Poisson
				// nothing to do : r=γ for Poisson, or κ for Pólya
				// but q = λ/μ
				double dqdλ = 1.0/μ;
				double dqdμ = -λ/(μ*μ);
				dL[3*node+PARAMETER_DUPLICATION] = dLdq * dqdλ;
				dL[3*node+PARAMETER_LOSS] = dLdq * dqdμ;
			} else
			{
				// shifted geometric

				double μt = μ*t;
				double λt = λ*t;
				double d = (μ-λ)*t;
				double E = Math.exp(-d);
				double E1 = -Math.expm1(-d); // 1.0-E;
				
				double denom = μt-λt*E;
				double divby = denom*denom;
				
				
				double dpdμ = E*(μt*d-λt*E1)/divby;
				double dpdλ = E*μt*(E1-d)/divby;
				double dqdμ = λt*(E*d-E1)/divby;
				double dqdλ = (E1*μt-λt*d*E)/divby;
				
				dL[3*node+PARAMETER_LOSS] = dLdp*dpdμ + dLdq*dqdμ;
				dL[3*node+PARAMETER_DUPLICATION] = dLdp*dpdλ + dLdq*dqdλ;
			}
			assert Double.isFinite(dL[3*node+PARAMETER_DUPLICATION] );
			assert Double.isFinite(dL[3*node+PARAMETER_LOSS] );
			assert Double.isFinite(dL[3*node+PARAMETER_GAIN] );
		}
		
		return dL;
	}
	
	
	
	public double[] getRateComponentGradient(double[] distribution_gradient)
	{
		int num_nodes = factory.tree.getNumNodes();
		final double[] dL = new double[4*num_nodes];
		
		TreeWithRates rates = factory.rates;
		int node = 0;
		while (node < num_nodes-1) // but not for the root
		{
			double t = rates.getEdgeLength(node);	
			double μ = rates.getLossRate(node);
			double λ = rates.getDuplicationRate(node);

			double dLdp = distribution_gradient[3*node+PARAMETER_LOSS] ;
			double dLdq = distribution_gradient[3*node+PARAMETER_DUPLICATION] ;
			double dLdr = distribution_gradient[3*node+PARAMETER_GAIN]; 
			dL[4*node+PARAMETER_GAIN] = dLdr; // same, except for Poisson 
			
			if (Double.isInfinite(t))
			{
				assert (dLdp==0.0);
				// Pólya or Poisson
				// nothing to do : r=γ for Poisson, or κ for Pólya
				// but q = λ/μ
				double dqdλ = 1.0/μ;
				double dqdμ = -λ/(μ*μ);
				dL[4*node+PARAMETER_DUPLICATION] = dLdq * dqdλ;
				dL[4*node+PARAMETER_LOSS] = dLdq * dqdμ;
			}  else
			{
				double q = rates.getDuplicationParameter(node);
				if (q==0.0)
				{
					// Poisson
					double r = rates.getGainParameter(node);
	
					double mu_t = μ*t;
					double dpdμ = Math.exp(-mu_t);
					double drdμ = rates.getGainRate(node)*dpdμ;
					double drdγ = -Math.expm1(-mu_t); // 1.0-dpdm;
	
					double dLdμ; 
					if (r!=0.0)
					{
						;
						
						dL[4*node+PARAMETER_GAIN] = dLdr*drdγ;
						//dL[4*node+PARAMETER_LOSS] 
						dLdμ = dLdp * dpdμ + dLdr * drdμ;
						
					} else
					{
						dLdμ =  dLdp * dpdμ;
					}
					dL[4*node+PARAMETER_LOSS] = dLdμ * t;
					dL[4*node+PARAMETER_LENGTH] = dLdμ * μ;
				} else
				{
					// Pólya
					double μt = μ*t;
					double λt = λ*t;
					
	
					if (μ == λ) // symmetric rates p=q=mu *t/(1+mu*t)
					{
						double denom = 1.0+μt;
						double divby = denom*denom;
						
						double dpdμ = 1.0/divby;
						double dqdμ = dpdμ;
						
						double dLdμ = dLdp*dpdμ + dLdq*dqdμ; 
						dL[4*node+PARAMETER_LOSS] 
						= dL[4*node+PARAMETER_DUPLICATION]
						= dLdμ * t;
						dL[4*node+PARAMETER_LENGTH] = dLdμ * μ;
					} else 
					{ 
						double gap = rates.getRateGap(node);
						double d = t*gap; //    (μ-λ)*t; // μt * gap; // (μ-λ)*t;
						double E = Math.exp(-d);
						double E1 = -Math.expm1(-d); // 1.0-E;

						double denom;
						double delta = gap/μ; //positive or not 
						
						if (delta<0.5 && 0.0<delta)
						{
							denom = μt * (-Math.expm1(-d+Math.log1p(-delta)));
						} else
						{
							denom = μt-λt*E; // and not d*E = (μt-λt)*E; // BUG 4/13/2023
						}
						double divby = denom*denom;
						
						double dpdμ = E*(μt*d-λt*E1)/divby;
						double dpdλ = E*μt*(E1-d)/divby;
						double dqdμ = λt*(E*d-E1)/divby;
						double dqdλ = (E1*μt-λt*d*E)/divby;
						
						double zdt = E*d*d/divby;
						
						double dLdμ = dLdp*dpdμ + dLdq*dqdμ;
						double dLdλ = dLdp*dpdλ + dLdq*dqdλ;
						
						dL[4*node+PARAMETER_LOSS] =  t*dLdμ;
						dL[4*node+PARAMETER_DUPLICATION] = t*dLdλ;
						dL[4*node+PARAMETER_LENGTH] = (dLdp * μ + dLdq * λ) * zdt; // dLdμ*μ+ dLdλ*λ; //  	
						
					} // mu/lambda
				} // q>0.0
			} // finite edge length
			assert Double.isFinite(dL[4*node+PARAMETER_DUPLICATION] );
			assert Double.isFinite(dL[4*node+PARAMETER_LOSS] );
			assert Double.isFinite(dL[4*node+PARAMETER_GAIN] );
			assert Double.isFinite(dL[4*node+PARAMETER_LENGTH] );
			
			++node;
		}			
		assert (factory.tree.isRoot(node)); // no rate transformation for root
		{
			double dLdp = distribution_gradient[3*node+PARAMETER_LOSS] ;
			double dLdq = distribution_gradient[3*node+PARAMETER_DUPLICATION] ;
			double dLdr = distribution_gradient[3*node+PARAMETER_GAIN]; 
			dL[4*node+PARAMETER_GAIN] = dLdr; // same

			double μ = rates.getLossRate(node);
			double λ = rates.getDuplicationRate(node);
			double t = rates.getEdgeLength(node);	
			if (Double.isInfinite(t))
			{
				assert (dLdp==0.0);
				// Pólya or Poisson
				// nothing to do : r=γ for Poisson, or κ for Pólya
				// but q = λ/μ
				double dqdλ = 1.0/μ;
				double dqdμ = -λ/(μ*μ);
				dL[4*node+PARAMETER_DUPLICATION] = dLdq * dqdλ;
				dL[4*node+PARAMETER_LOSS] = dLdq * dqdμ;
			} else
			{
				// shifted geometric

				double μt = μ*t;
				double λt = λ*t;
				double d = (μ-λ)*t;
				double E = Math.exp(-d);
				double E1 = -Math.expm1(-d); // 1.0-E;
				
				double denom = μt-λt*E;
				double divby = denom*denom;
				
				
				double dpdμ = E*(μt*d-λt*E1)/divby;
				double dpdλ = E*μt*(E1-d)/divby;
				double dqdμ = λt*(E*d-E1)/divby;
				double dqdλ = (E1*μt-λt*d*E)/divby;
				
				double zdt = E*d*d/divby;
				
				dL[4*node+PARAMETER_LOSS] = dLdp*dpdμ + dLdq*dqdμ;
				dL[4*node+PARAMETER_DUPLICATION] = dLdp*dpdλ + dLdq*dqdλ;
				dL[4*node+PARAMETER_LENGTH] = (dLdp * μt + dLdq * λt) * zdt;
				
			}
			assert Double.isFinite(dL[4*node+PARAMETER_DUPLICATION] );
			assert Double.isFinite(dL[4*node+PARAMETER_LOSS] );
			assert Double.isFinite(dL[4*node+PARAMETER_GAIN] );
			assert Double.isFinite(dL[4*node+PARAMETER_LENGTH] );
		}
		
		return dL;
		
		
	}
	
	
//	double[][][] getPosteriorSums()
//	{
//		class PartialD extends RecursiveTask<double[][][]>
//		{
//			double[] node_means = post.getNodeMeans(); // array of expected values across the nodes
//			double[] edge_means = post.getEdgeMeans(); 
//			
//		}
//	}
//	

	private class PosteriorStatistics 
	{
		PosteriorStatistics()
		{
			int num_nodes = factory.tree.getNumNodes();
			this.node_posteriors = new double[num_nodes][];
			this.edge_posteriors = new double[num_nodes][];
			this.birth_tails = new double[num_nodes][];
			this.death_tails = new double[num_nodes][];
			this.profile_count = 0.0;
			this.LL = 0.0;
		}
		
		// TODO: use birth and death tails
		
		private final double[][] node_posteriors;
		private final double[][] edge_posteriors;
		
		private final double[][] birth_tails;
		private final double[][] death_tails;
		
		
		private double profile_count;
		private double LL;
		private double LLunobs = Double.NEGATIVE_INFINITY;
		
		void add(Profile P, double multiplier)
		{
			// TODO: use birth and death tails
			
			for (int node=0; node<node_posteriors.length; node++)
			{
				double[] pN = P.post.getNodePosteriors(node);
				double[] pS = P.post.getEdgePosteriors(node);
				
				node_posteriors[node] = ArraySum.addCells(node_posteriors[node], pN, multiplier);
				edge_posteriors[node] = ArraySum.addCells(edge_posteriors[node], pS, multiplier);
				
//				{ // DEBUG
//					double S = P.post.getEdgeMean(node);
//					double logS = P.post.getLogEdgeMean(node);
//					double N = P.post.getNodeMean(node);
//					double logN_S = P.post.getLogNodeIncrease(node);
//					
//					StringBuilder sb = new StringBuilder();
//					if (!factory.tree.isRoot(node))
//					{
//						int parent = factory.tree.getParent(node);
//						double Nu = P.post.getNodeMean(parent);
//						double logNu_S = P.post.getLogEdgeDecrease(node);
//						sb.append("\tNu-S "+(Nu-S)+"/ "+Math.exp(logNu_S));
//					}
//					
//					System.out.println("#**G.PS.a "+P.post.inside.family_idx+"\tnode "+node+"\tS "+S+"\t/ "+Math.exp(logS)
//							+"\tN-S "+(N-S)+"/ "+Math.exp(logN_S)+sb.toString());
//				}
				
				double[] pNv_Sv = P.post.getNodeBirthTails(node);
				double[] pNu_Sv = P.post.getEdgeDeathTails(node);
				birth_tails[node] = ArraySum.addCells(birth_tails[node], pNv_Sv, multiplier);
				death_tails[node] = ArraySum.addCells(death_tails[node], pNu_Sv, multiplier);
			}
			profile_count += multiplier;
			LL = Math.fma(multiplier, P.post.inside.getLogLikelihood(), LL);
		}
		
		void add(PosteriorStatistics that)
		{
			// TODO: use birth and death tails
			for (int node=0; node<node_posteriors.length; node++)
			{
				node_posteriors[node] = ArraySum.addCells(this.node_posteriors[node], that.node_posteriors[node]);
				edge_posteriors[node] = ArraySum.addCells(this.edge_posteriors[node], that.edge_posteriors[node]);
				birth_tails[node] = ArraySum.addCells(this.birth_tails[node], that.birth_tails[node]);
				death_tails[node] = ArraySum.addCells(this.death_tails[node], that.death_tails[node]);
			}
			this.profile_count += that.profile_count;
			this.LL += that.LL;
		}
		
		/**
		 * Gradient of the uncorrected log-likelihood.
		 * Index is 3*<var>node</var>+<var>i</var> where <var>i</var>={@link GLDParameters#PARAMETER_DUPLICATION},
		 * {@link GLDParameters#PARAMETER_GAIN}, {@link GLDParameters#PARAMETER_LOSS}.
		 * 
		 * Problematic when p/q are too close to 0, because of (a-b)/p style terms 
		 * where a-b is not calculated precisely.  
		 * 
		 * @deprecated
		 * @return gradient array
		 */
		double[] getSurvivalGradient()
		{
			// TODO: use birth and death tails

			int num_nodes = factory.tree.getNumNodes();

			// I. compute the posterior tails and means
			double[][] node_tails = new double[num_nodes][];
			double[][] edge_tails = new double[num_nodes][];	
			double[] node_means = new double[num_nodes];
			double[] edge_means = new double[num_nodes];
			
			
			for (int node=0; node<num_nodes; node++)
			{
				node_tails[node] = ArraySum.tail(node_posteriors[node]);
				node_means[node] = ArraySum.sum(node_tails[node]);
				edge_tails[node] = ArraySum.tail(edge_posteriors[node]);
				edge_means[node] = ArraySum.sum(edge_tails[node]);
			}			
			
			// II. compute the gradients 
			double[] dL = new double[num_nodes*3];// return value
			
			for (int node=0; node<num_nodes; node++)
			{
				double q = factory.getDuplicationParameter(node);
				double q1 = factory.getDuplicationParameterComplement(node);
				double log1_q = factory.getLogDuplicationComplement(node); // Math.log(q1); // Math.log1p(-q); //Math.log(1.0-q); // 
				double log_q = factory.getLogDuplicationParameter(node);
				
				// II.a gain
				if (log_q==Double.NEGATIVE_INFINITY) 
				{
					// Poisson // untested
					double r = factory.getGainParameter(node);
					if (r!=0.0)
					{
						double Nu_Su = 0.0; // difference of means 
						double[] tNu = node_tails[node];
						double[] tSu = edge_tails[node];
						assert (tSu.length<=tNu.length);
						int i; 
						for (i=0; i<tSu.length; i++)
						{
							double diff = Double.max(0.0, tNu[i]-tSu[i]);
							Nu_Su += diff;
						}
						for (; i<tNu.length; i++)
						{
							Nu_Su += tNu[i];
						}						
						dL[3*node+PARAMETER_GAIN] 
//								= (node_means[node] - edge_means[node])/r - profile_count; // d/dr
								= Nu_Su/r - profile_count; // d/dr
					}
					dL[3*node+PARAMETER_DUPLICATION] = 0.0; // nothing to do 
				} else
				{
					// Pólya
					double κ = factory.getGainParameter(node);
					double Nu_Su;
					double[] tNu = node_tails[node];
					double[] tSu = edge_tails[node];
					assert (tSu.length<=tNu.length);

					if (κ!=0.0)
					{
						double dLdk = 0.0;
						Nu_Su = 0.0; // difference of means
						{
							int i; 
							for (i=0; i<tSu.length; i++)
							{
								double diff = Double.max(0.0, tNu[i]-tSu[i]);
								dLdk += diff/(κ + i);
								Nu_Su += diff;
							}
							for (; i<tNu.length; i++)
							{
								dLdk += tNu[i]/(κ + i);
								Nu_Su += tNu[i];
							}
						}
//						double log1_q = factory.getLogDuplicationComplement(node); // Math.log(q1); // Math.log1p(-q); //Math.log(1.0-q); // 
						dL[3*node+PARAMETER_GAIN] 
								= dLdk + profile_count*log1_q;// d/dκ
						
						assert !Double.isNaN(dL[3*node+PARAMETER_GAIN] );
					} else // kappa = 0?
					{
						Nu_Su = 0.0;
						{
							int i; 
							for (i=0; i<tSu.length; i++)
							{
								double diff = Double.max(0.0, tNu[i]-tSu[i]);
								Nu_Su += diff;
							}
							for (; i<tNu.length; i++)
							{
								Nu_Su += tNu[i];
							}
						}
					}
					// if (Nu_Su<TOO_SMALL) Nu_Su = 0.0;
					// II.b duplication

					dL[3*node+PARAMETER_DUPLICATION]
							= Math.exp(Math.log(Nu_Su)-log_q)
								- Math.exp(Math.log(edge_means[node]+ profile_count*κ)-log1_q);
//							(Nu_Su==0.0?0.0:Nu_Su/q)
//							- (edge_means[node]+ profile_count*κ)/q1; // d/dq
				}
				double logp = factory.getLogLossParameter(node);
				double log1_p = factory.getLogLossComplement(node);
				if (logp!=0.0)
				{
					if (factory.tree.isRoot(node)) // untested model setting
					{
						double p = factory.getLossParameter(node);
						double p1 = factory.getLossParameterComplement(node);
						dL[3*node+PARAMETER_LOSS]
								= (profile_count-edge_means[node])/p - edge_means[node]/p1;
					} else
					{
						int parent = factory.tree.getParent(node);
						
						double[] tNu = node_tails[parent];
						double[] tSv = edge_tails[node];
						double Nu_Sv = 0.0;
						double Sv = 0.0;
						assert (tSv.length<=tNu.length);
						int i; 
						for (i=0; i<tSv.length; i++)
						{
							double diff = Double.max(0.0, tNu[i]-tSv[i]);
							Nu_Sv += diff;
							Sv += tSv[i];
						}
						for (; i<tNu.length; i++)
						{
							Nu_Sv += tNu[i];
						}
						//if (Nu_Sv<TOO_SMALL) Nu_Sv = 0.0;
						

						double log_szer;
//						double p = factory.getLossParameter(node);
//						double p1 = factory.getLossParameterComplement(node);
//						double szer;
						if (factory.tree.getNumChildren(parent)==2)
						{
							// epsi = pv * pw if 2 children, so pv-epsi = pv*(1-pw)
//							szer = p * factory.getLossParameterComplement(factory.tree.getSibling(node))/p1;
							log_szer = logp-log1_p +factory.getLogLossComplement(factory.tree.getSibling(node));
						} else
						{
							double epsi = factory.getExtinction(parent);
//							szer =   (p - epsi)/p1; // = p*(1-epsi/p)/p1 // 
							log_szer = logp-log1_p+Math.log1p(-Math.exp(factory.getLogExtinction(parent)-logp));
						}
//						double szd = Math.exp(log_szer)-szer;
						double ldiff = Nu_Sv - Math.exp(log_szer + Math.log(Sv)); // (Sv==0.0?0.0:szer*Sv);
						
//						dL[3*node+PARAMETER_LOSS]
//						= ldiff==0.0?0.0:ldiff/p;	

						if (ldiff<0.0)
						{
							dL[3*node+PARAMETER_LOSS] = 
									-Math.exp(Math.log(-ldiff)-logp);
						} else if (0.0<ldiff)
						{
							dL[3*node+PARAMETER_LOSS] = 
									Math.exp(Math.log(ldiff)-logp);
						} else // ldiff==0.0
						{
							assert (ldiff==0.0);
							dL[3*node+PARAMETER_LOSS] = 0.0;
						}
						
						// Nu_Sv/p - szer*Sv/p
						// = Nu_Sv/p - Sv*(1-epsi/p)/(1-p)
						
						double log_szer2; 
						if (factory.tree.getNumChildren(parent)==2)
						{
							// epsi = pv * pw if 2 children, so pv-epsi = pv*(1-pw)
//							szer = p * factory.getLossParameterComplement(factory.tree.getSibling(node))/p1;
							log_szer2 = factory.getLogLossComplement(factory.tree.getSibling(node));
						} else
						{
							double epsi = factory.getExtinction(parent);
//							szer =   (p - epsi)/p1; // = p*(1-epsi/p)/p1 // 
							log_szer2 = Math.log1p(-Math.exp(factory.getLogExtinction(parent)-logp));
						}			
						double dl2
							=
								Math.exp(Math.log(Nu_Sv)-logp)
								-Math.exp(Math.log(Sv)+log_szer2-log1_p);
//						double logt1=Math.log(Nu_Sv)-logp;
//						double szd = dl2-dL[3*node+PARAMETER_LOSS];
//						System.out.println("#**G.gSG "+node+"\tszd "+szd+"\tszer2 "+log_szer2+"\tdl2 "+dl2+"\tdlL "+-dL[3*node+PARAMETER_LOSS]);
//						if (Math.abs(dl2)>1e10 || (Nu_Sv<1e-12 && Nu_Sv!=0.0) || logp<Math.log(0.001))
//						{
//							System.out.println("#**G.PS.gSG "+node+"\tszer2 "+log_szer2+"\tdl2 "+dl2+"\tdlL "+dL[3*node+PARAMETER_LOSS]);
//							System.out.println("#**G.PS.gSG "+node+"\tldiff "+ldiff+"\tlog_szer "+log_szer+"\tSv "+Sv+"\tNu-Sv "+Nu_Sv+"\tlogp "+logp+"\tlog1p "+log1_p+"\tdl2 "+dl2
//									+"\tt+ "+Math.exp(Math.log(Nu_Sv)-logp)
//									+"\tt- "+Math.exp(Math.log(Sv)+log_szer2-log1_p));
//							
//							// problem p too small:
//							// Nu_Sv is very small too; their ratio 
//							// Nu_Sv/p is hard to estimate
//						}
						dL[3*node+PARAMETER_LOSS] = dl2;
//						if (factory.rates.getLossParameter(node)<1e-4 )
//						{
//							System.out.println("#**G.PS.gSG node "+node+"\tSv "+Sv+"\tNu_Sv "+Nu_Sv+"\tszer "+szer+"\tlik "+factory.toString(node)+"\trates "+factory.rates.toString(node));
//						}
					}
				} else // with p==1.0, keep dLdp=0
				{
				}
				
//				if (p<1e-3 )
//				{
//					System.out.println("#**G.PS.gSG node "+node+"\tdL "+dL[3*node+PARAMETER_LOSS]+"\t"+factory.toString(node));
//				}				
				assert Double.isFinite(dL[3*node+PARAMETER_DUPLICATION] );
				assert Double.isFinite(dL[3*node+PARAMETER_LOSS] );
				assert Double.isFinite(dL[3*node+PARAMETER_GAIN] );
			}
			
//			System.out.println("#**G.P.gSG "+Arrays.toString(dL));
			
//			if (FunctionMinimization.euclideanNorm(dL)>1e10)
//			{ // DEBUG
//				printGradient(System.out, dL);
//			}
			
			return dL;
		}
		
	}
	
	private PosteriorStatistics computePosteriorStatistics()
	{
		int nF = factory.table.getFamilyCount();
		final int unit_task = Count.unitTask(nF);
		final ForkJoinPool thread_pool = threadPool();
		
		class PartialS extends RecursiveTask<PosteriorStatistics>
		{
			/**
			 * First family index.
			 */
			private final int minF;
			/**
			 * Last family index, exclusive.
			 */
			private final int maxF;

			PartialS(int min, int max)
			{
				this.minF = min;
				this.maxF = max;
			}
			@Override
			protected PosteriorStatistics compute() 
			{
				try
				{
					PosteriorStatistics S;
					if (maxF-minF > unit_task)
					{
						int medF = (minF+maxF)/2;
						PartialS left = new PartialS(minF, medF);
						PartialS right = new PartialS(medF, maxF);
						left.fork();
						S = right.compute();
						PosteriorStatistics S2 = left.join();
						S.add(S2);
					} else
					{
						S = new PosteriorStatistics();
						for (int f=minF; f<maxF; f++)
						{
							Profile G = getGradient(f);
							S.add(G, getMultiplicity(f));
						}
					}
					return S;
				} catch (Throwable t)
				{
					throw t;
				}
			}
		}
			
		PartialS bigjob = new PartialS(0, nF);
		try
		{
			PosteriorStatistics S;
			if (nF>unit_task)
			{
				S = thread_pool.invoke(bigjob);
			} else
			{
				S = bigjob.compute();
			}
			return S;
		} catch (Throwable t)
		{
			throw new RuntimeException(t);
		}
	}
	
	public double[] altCorrectedGradient()
	{
		PosteriorStatistics stats = computePosteriorStatistics();
		double[] dL = stats.getSurvivalGradient();
		
		double L0 = getUnobservedLL();
		double corr0 = -Math.exp(L0)/Math.expm1(L0); // ==0 if L0==Double.NEGATIVE_INFINITY;
		double[] emptyD = getUnobservedSurvivalGradient();
		
//		System.out.println("#**G.gCG dL "+FunctionMinimization.euclideanNorm(dL)
//			+"\tdL0 "+FunctionMinimization.euclideanNorm(emptyD)
//			+"\tcorr0 "+corr0);
//				
//				nF "+stats.profile_count+"\tcorr0 "+corr0+"\tL0 "+L0);

		
//		double nF = getTotalFamilyCount();
		corr0 *= stats.profile_count;
		
		for (int pidx=0; pidx<dL.length; pidx++)
		{
			dL[pidx]+=corr0*emptyD[pidx];
		}

		return dL;
	}
	
	/**
	 * Computes the array of partial derivatives by survival parameters.
	 * Index is 3*<var>node</var>+<var>i</var> where <var>i</var>={@link GLDParameters#PARAMETER_DUPLICATION},
	 * {@link GLDParameters#PARAMETER_GAIN}, {@link GLDParameters#PARAMETER_LOSS}.
	 * 
	 * @return gradient array
	 */
	public double[] getCorrectedGradient()
	{
		int nF = factory.table.getFamilyCount();
		final int unit_task = Count.unitTask(nF);
		final ForkJoinPool thread_pool = threadPool();
//		if (thread_pool != null)
//		{
//			unit_task = THREAD_UNIT_TASK;
//		}
//		else
//		{
//			unit_task = nF; // do not fork
//		}

		double L0 = getUnobservedLL();
		double corr0 = -Math.exp(L0)/Math.expm1(L0); // ==0 if L0==Double.NEGATIVE_INFINITY;
		double[] emptyD = getUnobservedSurvivalGradient();
//		{ // DEBUG
//			for (int j=0; j<emptyD.length; j++)
//			{
//				System.out.printf("#**G.gCG %d/%d\temptyD %.12f\n", j/3, j%3, emptyD[j]);
//			}
//		}
		
		class PartialD extends RecursiveTask<double[]>
		{
			/**
			 * First family index.
			 */
			private final int minF;
			/**
			 * Last family index, exclusive.
			 */
			private final int maxF;

			PartialD(int min, int max)
			{
				this.minF = min;
				this.maxF = max;
			}
			
			@Override
			protected double[] compute() 
			{
				try
				{
					double[] dL;
					if (maxF-minF>unit_task)
					{
						int medF = (minF+maxF)/2;
						PartialD left = new PartialD(minF, medF);
						PartialD right = new PartialD(medF, maxF);
						left.fork();
						dL = right.compute();
						double[] dL2 = left.join();
						assert (dL.length == dL2.length);
						for (int pidx=0; pidx<dL.length; pidx++)
						{
							dL[pidx]+=dL2[pidx];
						}
					} else
					{
						int num_nodes= factory.tree.getNumNodes();
						dL = new double[3*num_nodes];
						assert (dL.length == emptyD.length);
						for (int f = minF; f<maxF; f++)
						{
							Profile G = getGradient(f);
							double[] dF = G.getSurvivalGradient();
							assert (dF.length == dL.length);
							for (int pidx=0; pidx<dL.length; pidx++)
							{
								dL[pidx]+=(dF[pidx]+corr0*emptyD[pidx])*getMultiplicity(f);
							}
							
						}
					}
					return dL;
				} catch (Throwable t)
				{
					throw t;
				}
			}
		}
		
		PartialD bigjob = new PartialD(0, nF);
		double[] dL;
		try
		{
			if (nF>unit_task)
			{
				dL = thread_pool.invoke(bigjob);
			} else
			{
				dL = bigjob.compute();
			}
		} catch (Throwable t)
		{
			throw new RuntimeException(t);
		}
//		double[] dL = new double[factory.tree.getNumNodes()*3];// return value
//		
//		
//		int nF = factory.table.getFamilyCount();
//		for (int family_idx=0; family_idx<nF; family_idx++)
//		{
//			Profile G = new Profile(family_idx);
//			double[] D = G.getSurvivalGradient();
//			assert (D.length == dL.length);
//			for (int i=0; i<D.length; i++)
//				dL[i] += D[i];
//		}		
//		double logL0 = factory.getEmptyLL();
//		double p0 = Math.exp(logL0);
//		double p_not0 =-Math.expm1(logL0); // 1.0-Math.exp(logL0); // 1.0-Math.exp(logL0)
//		
//		
//		double[] emptyD = getEmptySurvivalGradient();
//		assert (emptyD.length == dL.length);
//		
//		double corr0 = nF * p0/p_not0;
//		for (int i=0; i<dL.length; i++)
//			dL[i] += corr0 * emptyD[i];
//
////		System.out.println("#**G.gCG corr "+corr0+"\tdL "+Arrays.toString(dL));
//		
		return dL;
	}
	
	protected static void printGradient(java.io.PrintStream out, double[] dL, int num_nodes)
	{
		for (int node =0; node<num_nodes; node++)
		{
			double dLdg = dL[3*node + PARAMETER_GAIN];
			double dLdp = dL[3*node + PARAMETER_LOSS];
			double dLdq = dL[3*node + PARAMETER_DUPLICATION];
			out.println(node+"\t"+dLdg+"\t"+dLdp+"\t"+dLdq);
		}
		for (int j=3*num_nodes; j<dL.length; j++)
		{
			int jd = j-3*num_nodes;
			out.println("p+"+jd+"\t"+dL[j]);
		}
		out.println("LENGTH\t"+FunctionMinimization.euclideanNorm(dL));
	}
	
	protected void printGradient(java.io.PrintStream out, double[] dL)
	{
		for (int node =0; node<factory.tree.getNumNodes(); node++)
		{
			double dLdg = dL[3*node + PARAMETER_GAIN];
			double dLdp = dL[3*node + PARAMETER_LOSS];
			double dLdq = dL[3*node + PARAMETER_DUPLICATION];
			out.println(node+"\t"+dLdg+"\t"+dLdp+"\t"+dLdq);
		}
		out.println("LENGTH\t"+FunctionMinimization.euclideanNorm(dL));
	}
	
	private void mainmain()
	{
		java.io.PrintStream out = System.out;
		
		double LL = getLL();
		out.println("# Uncorrected likelihood "+LL); //+"\t// "+factory.getLL());
		out.println("# Corrected likelihood "+getCorrectedLL());
		double varLL = bootstrapLLVariance();
		double sdLL = Math.sqrt(varLL);
		out.println("# Bootstrap variance "+varLL+"\tsd "+sdLL);
		out.println("# Survival parameter gradient (gain, loss, duplication order)");
		double[] D = getCorrectedGradient();
		printGradient(out, D);
		{ // DEBUG: via PosteriorStatistics
			out.println("# Alt survival gradient (GLD order)");
			double[] A = altCorrectedGradient();
			printGradient(out, A);
			for (int j=0; j<D.length; j++)
			{
				double delta = A[j]-D[j];
				double rd = delta/Math.abs(D[j]);
				out.println("#**G.mm par "+(j/3)+"/"+(j%3)+"\tdelta "+delta+"\trd "+rd);
			}
		}

		out.println("# Model parameter gradient (GLD order)");
		D = getDistributionGradient(D);
		printGradient(out, D);
		
		out.println("# Logistic parameter gradient (LD order)");
		for (int node =0; node<factory.tree.getNumNodes(); node++)
		{
			double dLdp = D[3*node + PARAMETER_LOSS];
			double dLdq = D[3*node + PARAMETER_DUPLICATION];
			double dLdx = dLdp * factory.rates.getLossParameter(node)*factory.rates.getLossParameterComplement(node); //  Math.exp(factory.getLogLossParameter(node)+factory.getLogLossComplement(node));
			double dLdy = dLdq * factory.rates.getDuplicationParameter(node)*factory.rates.getDuplicationParameterComplement(node); // Math.exp(factory.getLogDuplicationParameter(node)+factory.getLogDuplicationComplement(node));
			out.println(node+"\t"+dLdx+"\t"+dLdy);
		}
		
//		D = getTotalRateGradient(D);
//		out.println("# Rate parameter gradient (GLD order)");
//		printGradient(out, D);
	}
	
	public static void main(String[] args) throws Exception
	{
		Class<?> us = java.lang.invoke.MethodHandles.lookup().lookupClass();
		CommandLine cli = new CommandLine(args,  us);
		
//		boolean use_logit = cli.getOptionBoolean("logit", false);
		
		Gradient G;
//		if (use_logit)
//		{
//    		System.out.println(CommandLine.getStandardHeader("Using logit scale: -logit "+use_logit));
//			LikelihoodParametrized factory = new LikelihoodParametrized(cli.getRates(),cli.getTable());
//			G = new Gradient(factory);
//		} else
		{
			G = new Gradient(cli.getRates(), cli.getTable());
		}
		
		
		int absolute = G.getCalculationWidthAbsolute();
		double relative = G.getCalculationWidthRelative();
		
        if (cli.getOptionValue(OPT_TRUNCATE)!=null)
        {
        	absolute = cli.getOptionTruncateAbsolute();
        	relative = cli.getOptionTruncateRelative();
        	G.setCalculationWidthThresholds(absolute, relative);
        } 
		System.out.println(CommandLine.getStandardHeader("Truncated computation (absolute,relative)="
				+absolute+","+relative));
		int min_copies = cli.getOptionInt(OPT_MINCOPY, G.min_copies);
		G.setMinimumObservedCopies(min_copies);

		G.mainmain();
	}
	
}
