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


import java.util.Arrays;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;

import count.Count;
import count.ds.ProfileTable;
import count.ds.UniqueProfileTable;
import count.matek.Logarithms;
import count.model.StraightLikelihood.Profile;

import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LOSS;

/**
 * Unused yet/
 * 
 *
 */
public class StraightGradient implements Count.UsesThreadpool
{
//	private static final double MAX_RATE = 33.0; 
//	private static final double MIN_GAIN_PARAMETER = 1.0/(1L<<40);
		
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
	
	public StraightGradient(TreeWithRates rates, ProfileTable table)
	{
		if (table instanceof UniqueProfileTable)
		{
			this.utable = (UniqueProfileTable) table; 
		}else
		{
			this.utable = new UniqueProfileTable(table);
		}
		this.min_copies = Integer.min(2,table.minCopies());
		
		this.factory = new StraightLikelihood(rates, utable);
	}
	
	private final StraightLikelihood factory;

	/**
	 * Same as the instantiating table, if it 
	 * was a UniqueProfileTable; or else its unique-profile version. 
	 */
	private final UniqueProfileTable utable;
	
	private int min_copies;
	/**
	 * A table of singleton profiles  
	 */
	private ProfileTable singletons = null;
		
	/* ========================= */
	/* Setter and getter methods */ 
	/* ========================= */
	
	public void setMinimumObservedCopies(int min_copies)
	{
		if (min_copies<0 || min_copies>2)
			throw new UnsupportedOperationException("Minimum copies must be 0,1 or 2.");
		this.min_copies = min_copies;
	}

	public void setCalculationWidth(int absolute, double relative)
	{
		factory.setCalculationWidthThresholds(absolute, relative);
	}
	
	public void computeParameters()
	{
		factory.computeParameters();
	}
	
	/* ================================ */
	/* Calculating the gradient         */ 
	/* ================================ */

	/**
	 * Log-gradient 
	 * by logit(p), logit(q), log(kappa)/log(r)
	 * 
	 * @param S
	 * @return logD[][] array of derivatives as log-differences  
	 */
	public double[][] getLogDistributionGradient(PosteriorStatistics S)
	{
		int num_nodes = factory.tree.getNumNodes();
		double[][] logD = new double[3*num_nodes][]; // return value
	
		double F = S.profile_count;
		double logF = Math.log(F);
		for (int v=0; v<num_nodes; v++)
		{
			double[] pSv = S.log_edge_posteriors[v].clone();
			{ // calculating tails and then summing to get the mean 
				double log_tail = Double.NEGATIVE_INFINITY;
				int j=pSv.length;
				while (0<j)
				{
					--j;
					double x = pSv[j];
					pSv[j] = log_tail;
					log_tail = Logarithms.add(log_tail, x);
				}
			}
			double logSv = Logarithms.sum(pSv, pSv.length); 
			
			// 
			// 1. set loss-derivative
			// 
			
			//(1-p)*Nu_Sv - p*S
			double logp = factory.getLogLossParameter(v);
			double log1_p = factory.getLogLossComplement(v);
			if (factory.tree.isRoot(v) || logp==0.0)
			{
				logD[3*v + PARAMETER_LOSS] = Logarithms.ldiff(); 
			} else
			{
				double[] tNu_Sv = S.log_death_tails[v];
				double logNu_Sv = Logarithms.sum(tNu_Sv, tNu_Sv.length);
				double dpos = logNu_Sv+log1_p;
				double dneg = logSv+logp;
				logD[3*v + PARAMETER_LOSS] = Logarithms.ldiff(dpos,dneg); 
			}
			
			// 
			// 2. set gain-derivative
			// 
			double log1_q = factory.getLogDuplicationComplement(v);  
			double log_q = factory.getLogDuplicationParameter(v);
			
			double[] tNv_Sv = S.log_birth_tails[v];
			double logNv_Sv = Logarithms.sum(tNv_Sv, tNv_Sv.length);
			
			// TODO log gain params ?
			
			if (log_q==Double.NEGATIVE_INFINITY) // q==0.0
			{
				// Poisson
				double r = factory.getGainParameter(v);
				if (r==0.0)
				{
					logD[3*v+PARAMETER_GAIN] = Logarithms.ldiff();
				} else
				{
					double log_r = Math.log(r);
					logD[3*v+PARAMETER_GAIN] = Logarithms.ldiff(logNv_Sv, log_r+logF);
					logD[3*v+PARAMETER_DUPLICATION] = Logarithms.ldiff();
				}
			} else
			{
				// Polya
				double κ = factory.getGainParameter(v);
				double log_kappa = Math.log(κ);
				if (κ==0.0)
				{
					logD[3*v+PARAMETER_GAIN] = Logarithms.ldiff();
				} else
				{
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
					double dneg = (Math.log(-log1_q)+log_kappa)+logF; // +F*kappa*ln(1-q)
					logD[3*v+PARAMETER_GAIN] = Logarithms.ldiff(dpos,dneg);
				}
				// 
				// 3. set duplication- derivative
				//				
				double dpos = logNv_Sv+log1_q;
				double dneg = Logarithms.add(logSv, log_kappa+logF) + log_q;
				logD[3*v+PARAMETER_DUPLICATION] = Logarithms.ldiff(dpos,dneg);
			}
		}		
		return logD;
	}	

	/* ================================ */
	/* Collecting posterior statistics  */ 
	/* ================================ */
	
	
	/**
	 * Collects the posterior statistics across the sample.
	 * 
	 * 
	 * @return
	 */
	public PosteriorStatistics getSampleStatistics()
	{
		int nF = utable.getFamilyCount();
				
		final int unit_task = Count.unitTask(nF);
		ForkJoinPool thread_pool = threadPool();
		
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
						S.add(S2, 1.0);
					} else
					{
						S = new PosteriorStatistics();
						for (int f=minF; f<maxF; f++)
						{
//							System.out.println("#**DEM.gSS "+f+"\t"+Thread.currentThread().toString());
							Profile P = factory.getProfile(f);
							S.add(P, Math.log(utable.getMultiplicity(f)));
						}
					}
					
//					{ // DEBUG
//						for (int node=0; node<factory.tree.getNumNodes(); node++)
//						{
//							double[] pS = S.log_edge_posteriors[node];
//							double[] pN = S.log_node_posteriors[node];
//							double mS=0.0;
//							double tS=0.0; 
//							for (int s=0; s<pS.length; s++) 
//							{
//								pS[s]=Math.exp(pS[s]);
//								mS += s*pS[s];
//								tS += pS[s];
//							}
//							double mN=0.0;
//							double tN=0.0;
//							for (int n=0; n<pN.length; n++) 
//							{
//								pN[n]=Math.exp(pN[n]);
//								mN += n*pN[n];
//								tN += pN[n];
//							}
//							System.out.println("#**SEM.gSS "+node+"\tS "+mS+"("+tS+")"+"\tN "+mN+"("+tN+")"+"\t// "+Arrays.toString(pS)+"\t"+Arrays.toString(pN));
//						}
//					}
					
					
					
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
	
	/**
	 * Produces the survival posterior statistics for one unobserved 
	 * family. 
	 * 
	 * @return {@link PosteriorStatistics#LL} is set to the log-likelihood of unobserved profile
	 */
	public PosteriorStatistics getUnobservedStatistics()
	{
		
		PosteriorStatistics Sunobs = new PosteriorStatistics();
		Sunobs.LL = Double.NEGATIVE_INFINITY;
		
		if (0<this.min_copies)
		{
			int miss_calc_absolute = Integer.max(5,factory.getCalculationWidthAbsolute());
			double miss_calc_relative = Double.max(1.0, factory.getCalculationWidthRelative());
			
			int num_unobserved_profiles; 
			if (this.min_copies==1)
			{
				num_unobserved_profiles = 1;
			} else
			{
				assert (min_copies==2);
				if (this.singletons == null)	
					this.singletons = ProfileTable.singletonTable(factory.tree);
				int n1 = singletons.getFamilyCount(); 
				num_unobserved_profiles = n1+1;
			}
			Profile[] SP = new Profile[num_unobserved_profiles];
			
			int ui=0;
			StraightLikelihood SL0 = new StraightLikelihood(factory, ProfileTable.emptyProfile(factory.tree));
			SL0.setCalculationWidthThresholds(miss_calc_absolute, miss_calc_relative);
			Profile P0 = SL0.getProfile(0);
			SP[ui++] = P0;
			
			if (this.min_copies==2)
			{
				StraightLikelihood SL1 = new StraightLikelihood(factory, singletons);
				int n1 = singletons.getFamilyCount(); 
				for (int f=0; f<n1; f++)
				{
					Profile P1 = SL1.getProfile(f);
					SP[ui++] = P1;
				}
			}
			double LLunobs = Double.NEGATIVE_INFINITY;
			double[] unobsLL = new double[num_unobserved_profiles];
			while (0<ui)
			{
				--ui;
				double SPLL = 
				unobsLL[ui] = SP[ui].getLogLikelihood();
				LLunobs = Logarithms.add(LLunobs,SPLL );
			}
			while (ui<SP.length)
			{
				Sunobs.add(SP[ui], unobsLL[ui] - LLunobs);
				ui++;
			}
			Sunobs.LL = LLunobs;
		}
		return Sunobs;
	}
	
	/**
	 * Inner class for collecting posterior statistics 
	 * across family profiles.
	 * 
	 * @author csuros
	 *
	 */
	public class PosteriorStatistics
	{
		// private final double[][] log_node_posteriors; // not needed
		private final double[][] log_edge_posteriors;
		private final double[][] log_birth_tails;
		private final double[][] log_death_tails;
		private double profile_count;
		private double LL;
		
		private PosteriorStatistics()
		{
			int num_nodes = factory.tree.getNumNodes();
			// this.log_node_posteriors = new double[num_nodes][];
			this.log_edge_posteriors = new double[num_nodes][];
			this.log_birth_tails = new double[num_nodes][];
			this.log_death_tails = new double[num_nodes][];
			this.profile_count = 0.0;
			this.LL = 0.0;
		}
		
		public void add(Profile P, double log_multiplier)
		{
			for (int node=0; node<log_edge_posteriors.length; node++)
			{
				//double[] pN = P.getLogNodePosteriors(node);
				double[] pS = P.getLogEdgePosteriors(node);
				double[] pNv_Sv = P.getLogBirthDifferenceTails(node);
				double[] pNu_Sv = P.getLogDeathDifferenceTails(node);
				//log_node_posteriors[node] = addCells(log_node_posteriors[node], pN, log_multiplier);
				log_edge_posteriors[node] = addCells(log_edge_posteriors[node], pS, log_multiplier);
				log_birth_tails[node] = addCells(log_birth_tails[node], pNv_Sv, log_multiplier);
				log_death_tails[node] = addCells(log_death_tails[node], pNu_Sv, log_multiplier);
			}
			
			double mul = Math.exp(log_multiplier);
			
			profile_count += mul;
			LL = Math.fma(mul, P.getLogLikelihood(), LL); // fma(a,b,c)=a*b+c
		}
		
		public void add(PosteriorStatistics that, double multiplier)
		{
			double logm = Math.log(multiplier);
			
			for (int node=0; node<log_edge_posteriors.length; node++)
			{
				//this.log_node_posteriors[node] = addCells(this.log_node_posteriors[node],that.log_node_posteriors[node],logm);
				this.log_edge_posteriors[node] = addCells(this.log_edge_posteriors[node],that.log_edge_posteriors[node],logm);
				this.log_birth_tails[node] = addCells(this.log_birth_tails[node],that.log_birth_tails[node],logm);
				this.log_death_tails[node] = addCells(this.log_death_tails[node],that.log_death_tails[node],logm);
			}
			this.profile_count = Math.fma(multiplier, that.profile_count, this.profile_count); // += that.profile_count * multiplier;
			this.LL = Math.fma(multiplier, that.LL, this.LL); //  += that.LL*multiplier;
		}
		
		/**
		 * Adds b to a, cell by cell, reusing a and/or b if possible. 
		 * If a is null, then b is returned.   
		 * 
		 * @param a expanded if necessary to match length of b
		 * @param b untouched if a is  not null; 
		 * @return cell-by-cell sums in array a (or b if a is null), or the resized version  of a
		 */
		private double[] addCells(double[] a, double[] b, double log_mul)
		{
			if (b==null) return a;
			if (a==null)
			{
				a = b.clone();
				for (int i=0; i<b.length; i++)
					b[i]+=log_mul;
				return b;
			}
			if (a.length<b.length)
			{
				int alen = a.length;
				a = Arrays.copyOf(a, b.length);
				Arrays.fill(a, alen, a.length, Double.NEGATIVE_INFINITY);
			}
			for (int i=0; i<b.length; i++) 
				a[i] = Logarithms.add(a[i], b[i]+log_mul);
			return a;
		}
	}


}
