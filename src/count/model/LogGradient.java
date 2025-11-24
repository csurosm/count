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
import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.ds.UniqueProfileTable;
import count.io.CommandLine;
import count.matek.Logarithms;

import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_TRUNCATE;
import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LOSS;

/**
 * High-precision gradient computations.
 * <ol>
 * <li> collect the posteriors statistics with {@link #getSampleStatistics()}</li>
 * <li> calculate the survival gradient with {@link PosteriorStatistics#getLogSurvivalGradient()}</li>
 * <li> convert to distribution gradient with {@link #convertToLogDistributionGradient(double[][])}</li>
 * </ol>
 *
 * 
 *
 */
public class LogGradient extends Posteriors implements Count.UsesThreadpool
{
	private static final boolean LOGIT_TO_LOGIT = true;
	/**
	 * Minimum observation for generic correction 
	 */
	private static final int USE_FAMILY_SIZE_LIKELIHOOD = 3;
	
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

	
	public LogGradient(TreeWithRates rates, ProfileTable table)
	{
		this(new LikelihoodParametrized(rates, table));
	}
	
	public LogGradient(LikelihoodParametrized factory)
	{
		super(factory);
		this.log_factory = factory;
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
	
	private final LikelihoodParametrized log_factory;
	
	/**
	 * Observation bias : minimum number of observed copies in the families. 
	 */
	private int min_copies = 1;
	
	/**
	 * A table of singleton profiles  
	 */
	private ProfileTable singletons = null;
	
	
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
	 * Sets the observation bias : minimum number of observed copies in the families. 
	 * @param min_copies 0, 1, or 2
	 */
	public void setMinimumObservedCopies(int min_copies)
	{
//		if (min_copies<0 || min_copies>2)
//			throw new UnsupportedOperationException("Minimum copies must be 0,1 or 2.");
		this.min_copies = min_copies;
	}
	
	public int getMinimumObservedCopies()
	{
		return min_copies;
	}
	
	public void computeParameters()
	{
		log_factory.computeParameters();
	}
	
	/* ================================ */
	/* Collecting survival gradient  */ 
	/* ================================ */
	
//	public double[][] getLogSurvivalGradient(Profile P, double[][] logSD)
//	{
//		int num_nodes = factory.tree.getNumNodes();
//		if (logSD == null)
//		{
//			logSD = new double[3*num_nodes][];
//		} else
//		{
//			assert logSD.length == 3*num_nodes;
//		}
//		
//		for (int v=0; v<num_nodes; v++)
//		{
//			double[] pSv = P.getLogEdgePosteriors(v); // clone bc we destroy it 
//			{ // calculating tails and then summing to get the mean 
//				double log_tail = Double.NEGATIVE_INFINITY;
//				int j=pSv.length;
//				while (0<j)
//				{
//					--j;
//					double x = pSv[j];
//					pSv[j] = log_tail;
//					log_tail = Logarithms.add(log_tail, x);
//				}
//			}
//			double logSv = Logarithms.sum(pSv, pSv.length); 
//			
//			// set loss-derivative
//			
//			//(1-p)*Nu_Sv - p*S
//			
//			double logp = log_factory.getLogLossParameter(v);
//			double log1_p = log_factory.getLogLossComplement(v);
//			
//			if (factory.tree.isRoot(v) || logp==0.0)
//			{
//				logSD[3*v + PARAMETER_LOSS] = Logarithms.ldiff(); 
//			} else
//			{
//				double[] tNu_Sv = P.getLogEdgePosteriorDecrease(v);
//				double logNu_Sv = Logarithms.sum(tNu_Sv, tNu_Sv.length);
//				
//				int u = factory.tree.getParent(v);
//				double log1_e;
//				if (factory.tree.getNumChildren(u)==2)
//				{
//					log1_e = factory.getLogLossComplement(factory.tree.getSibling(v));
//				} else
//				{
//					log1_e = Logarithms.logToLogComplement(log_factory.getLogExtinction(u)-logp);
//				}
//				double dpos = logNu_Sv+log1_p;
//				double dneg = logSv+logp+log1_e;
//				logSD[3*v + PARAMETER_LOSS] = Logarithms.ldiff(dpos,dneg); 
//			}
//			
//			
//			double log1_q = log_factory.getLogDuplicationComplement(v);  
//			double log_q = log_factory.getLogDuplicationParameter(v);
//			
//			double[] tNv_Sv = P.getLogNodePosteriorIncrease(v); 
//			double logNv_Sv = Logarithms.sum(tNv_Sv, tNv_Sv.length);
//			
//			// gain
//			if (log_q==Double.NEGATIVE_INFINITY) 
//			{
//				// Poisson
//				double r = factory.getGainParameter(v);
//				if (r==0.0)
//				{
//					logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff();
//				} else
//				{
//					double log_r = Math.log(r);
//					logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff(logNv_Sv, log_r);
//					logSD[3*v+PARAMETER_DUPLICATION] = Logarithms.ldiff();
//				}
//			} else
//			{
//				// Polya
//				double κ = factory.getGainParameter(v);
//				double log_kappa = Math.log(κ);
//				if (κ==0.0)
//				{
//					logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff();
//				} else
//				{
//					double dpos = tNv_Sv[0]; // first term
//					for (int i=1; i<tNv_Sv.length; i++)
//					{
//						// ln(k/(k+i)) = ln(1/(1+i/k)) =-ln(1+i/k) for i<k
//						// = ln (k/i)-ln(1+k/i) for k<i
//						
//						double log_k_ki;
//						if (i<κ)
//						{
//							log_k_ki = -Math.log1p(i/κ);
//						} else
//						{
//							log_k_ki = log_kappa - Math.log(i) -Math.log1p(κ/i);
//						}
//						dpos = Logarithms.add(dpos, tNv_Sv[i]+log_k_ki);
//					}
//					double dneg = (Math.log(-log1_q)+log_kappa); // +F*kappa*ln(1-q)
//					logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff(dpos,dneg);
//					
//					// now for duplication
//					dpos = logNv_Sv+log1_q;
//					dneg = Logarithms.add(logSv, log_kappa) + log_q;
//					logSD[3*v+PARAMETER_DUPLICATION] = Logarithms.ldiff(dpos,dneg);
//				}
//			}
//		}		
//		return logSD;		
//	}

	
	
	/* ================================ */
	/* Collecting posterior statistics  */ 
	/* ================================ */
	
	/**
	 * For a single profile 
	 * 
	 * @param P
	 * @return
	 */
	public PosteriorStatistics getProfileStatistics(Profile P)
	{
		
		PosteriorStatistics S = new PosteriorStatistics(P);
		//S.add(P, 0.0);
		return S;
	}
	
	/**
	 * Collects the posterior statistics across the entire sample;
	 * used in testing.
	 * 
	 * 
	 * @return
	 */
	public PosteriorStatistics getSampleStatistics()
	{
		int nF = factory.table.getFamilyCount();
				
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
							Profile P = getPosteriors(f);
							S.add(P, Math.log(getMultiplicity(f)));
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
		
	/**
	 * Calculates statistics across unobserved profiles 
	 * 
	 * @return {@link PosteriorStatistics#LL} is set to the log-likelihood of unobserved profile
	 */
	public PosteriorStatistics getUnobservedStatistics()
	{
		PosteriorStatistics Sunobs;
		
		
		int min_copies = this.getMinimumObservedCopies();
		if (Integer.min(3,USE_FAMILY_SIZE_LIKELIHOOD)<=min_copies) {
			Posteriors.Profile unobsProfile = FamilySizeLikelihood.getUnobservedProfile(factory.rates, min_copies-1);
			Sunobs = new PosteriorStatistics(unobsProfile);
		} else {
			assert (min_copies <= 2);
			// constructing by collecting the possible profiles
			Sunobs = new PosteriorStatistics();
			Sunobs.LL = Double.NEGATIVE_INFINITY;
			if (0<min_copies)
			{
				int num_unobs ; // unobserved profiles
				if (min_copies == 1)
				{
					num_unobs  = 1; // empty profile
				} else
				{
					assert (min_copies == 2);
					if (this.singletons == null)	
						this.singletons = ProfileTable.singletonTable(factory.tree);
					num_unobs = 1 + singletons.getFamilyCount();
				}
				Profile[] SP = new Profile[num_unobs];
				LogGradient SG0 = new LogGradient(new LikelihoodParametrized(log_factory.rates, ProfileTable.emptyProfile(factory.tree)));
				// SG0.setCalculationWidthThresholds(this.getCalculationWidthAbsolute(), this.getCalculationWidthAbsolute());
				int ui = 0;
				SP[ui++] = SG0.getPosteriors(0);
				if (min_copies==2)
				{
					LogGradient SG1 = new LogGradient(new LikelihoodParametrized(log_factory.rates, singletons));
					// SG1.setCalculationWidthThresholds(this.getCalculationWidthAbsolute(), this.getCalculationWidthAbsolute());
					for (int f=0; f<singletons.getFamilyCount(); f++)
						SP[ui++] = SG1.getPosteriors(f);
				}
				double Lunobs = Double.NEGATIVE_INFINITY;
				double[] Pll = new double[SP.length];
				while (0<ui)
				{
					--ui;
					Profile P = SP[ui];
					Pll[ui]  = P.getLogLikelihood();
					Lunobs = Logarithms.add(Lunobs, Pll[ui]);
				}
				while (ui<SP.length)
				{
					Sunobs.add(SP[ui], Pll[ui]-Lunobs);
					++ui;
				}
				Sunobs.LL = Lunobs;
			} 
		}
		return Sunobs;
	}
	
	
	
	
	
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
		
		/**
		 * Instantiation for a single Profile 
		 * 
		 * @param P
		 */
		private PosteriorStatistics(Profile P)
		{
			int num_nodes = factory.tree.getNumNodes();
			// this.log_node_posteriors = new double[num_nodes][];
			this.log_edge_posteriors = new double[num_nodes][];
			this.log_birth_tails = new double[num_nodes][];
			this.log_death_tails = new double[num_nodes][];
			this.profile_count = 1.0;
			this.LL = P.getLogLikelihood();
			for (int node=0; node<log_edge_posteriors.length; node++)
			{
				log_edge_posteriors[node] = P.getLogEdgePosteriors(node);
				log_birth_tails[node] = P.getLogNodePosteriorIncrease(node);
				log_death_tails[node] = P.getLogEdgePosteriorDecrease(node);
				
//				// DEBUG
//				if (hasNaN(log_death_tails[node])) {
//					System.out.println("#**LG.PS(Profile#"+P.inside.family_idx+")"
//							+"\tnode "+node
//							+"\tld "+Arrays.toString(log_death_tails[node])
//							+"\t// "+P.inside.getOwner().toString(node)
//							+"\t// "+P.inside.getOwner().getRates().toString(node)
//							);
//				}
				
				
			}			
		}
		
		/*
		 * Access methods 
		 */
		public double getLogLikelihood(){
			return this.LL;
		}
		
		public double[] getLogEdgePosteriors(int node) {
			return this.log_edge_posteriors[node];
		}
		
		public double[] getLogBirthTails(int node) {
			return this.log_birth_tails[node];
		}
		
		public double[] getLogDeathTails(int node) {
			return this.log_death_tails[node];
		}
		
		/*
		 * Adding data
		 */
		
		public void add(Profile P, double log_multiplier)
		{
//			int Psize =0;
//			for (int n: P.inside.get()) Psize += n;
			for (int node=0; node<log_edge_posteriors.length; node++)
			{
				//double[] pN = P.getLogNodePosteriors(node);
				double[] pS = P.getLogEdgePosteriors(node);
				double[] pNv_Sv = P.getLogNodePosteriorIncrease(node);
				double[] pNu_Sv = P.getLogEdgePosteriorDecrease(node);
				
//				// DEBUG
//				if (hasNaN(pNu_Sv)) {
//					System.out.println("#**LG.PS.add "+P.inside.family_idx
//							+"*"+Math.exp(log_multiplier)+"\tnode "+node+"\tpNuSv "+Arrays.toString(pNu_Sv)
//							+"\t// ");
//				}
				
				
				
				//log_node_posteriors[node] = addCells(log_node_posteriors[node], pN, log_multiplier);
				log_edge_posteriors[node] = addCells(log_edge_posteriors[node], pS, log_multiplier);
				log_birth_tails[node] = addCells(log_birth_tails[node], pNv_Sv, log_multiplier);
				log_death_tails[node] = addCells(log_death_tails[node], pNu_Sv, log_multiplier);

//				// DEBUG
//				if (hasNaN(log_death_tails[node])) {
//					System.out.println("#**LG.PS.add "+P.inside.family_idx
//							+"*"+Math.exp(log_multiplier)+"\tnode "+node+"\tpNuSv "+Arrays.toString(pNu_Sv)
//							+"\tld "+Arrays.toString(log_death_tails[node])
//							);
//				}
				
				
//				if (Psize==0)// DEBUG
//				{
//					System.out.println("#**LG.PS.a empty\tpS "+Arrays.toString(pS)
//							+"\tNvSv "+Arrays.toString(pNv_Sv)
//							+"\tNu_Sv "+Arrays.toString(pNu_Sv));
//				}
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

//				// DEBUG
//				if (hasNaN(log_death_tails[node])) {
//					System.out.println("#**LG.PS.add this("+this.profile_count+")+that("+that.profile_count+")*"+multiplier 
//							+"\tnode "+node
//							+"\tld "+Arrays.toString(log_death_tails[node])
//							);
//				}
			
			
			}
			this.profile_count = Math.fma(multiplier, that.profile_count, this.profile_count); // += that.profile_count * multiplier;
			this.LL = Math.fma(multiplier, that.LL, this.LL); //  += that.LL*multiplier;
		}
		
		
		
		/*
		 * Gradient computation
		 */
		
		/**
		 * Log-gradient 
		 * by logit(p~), logit(q~), log(kappa)/log(r~)
		 * 
		 * @return
		 */
		public double[][] getLogSurvivalGradient()
		{
			int num_nodes = factory.tree.getNumNodes();
			double[][] logSD = new double[3*num_nodes][]; // return value
			
			double logF = Math.log(this.profile_count);
			for (int v=0; v<num_nodes; v++)
			{
				double[] pSv = log_edge_posteriors[v].clone(); // clone bc we destroy it 
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
				
				// set loss-derivative
				
				//(1-p)*Nu_Sv - p*S
				
				double logp = log_factory.getLogLossParameter(v);
				double log1_p = log_factory.getLogLossComplement(v);
				
				if (log1_p==Double.NEGATIVE_INFINITY) //  logp==0.0) // || factory.tree.isRoot(v) || 
				{
					logSD[3*v + PARAMETER_LOSS] = Logarithms.ldiff(); 
				} else
				{
					double[] tNu_Sv = log_death_tails[v];
					double logNu_Sv = Logarithms.sum(tNu_Sv, tNu_Sv.length);
					
					double dpos, dneg;
					double log1_e;
					
					if (factory.tree.isRoot(v) ) {
						// experimental ROOTLOSS  
						// root with loss from N=1
						dpos = logNu_Sv + log1_p;
						dneg = logSv + logp;
						log1_e = 0.0; // e==0.0. log(1-pe) = 0
						
//						double[] ldiff = Logarithms.ldiff(dpos, dneg);
//						
//						// DEBUG // prints for each profile ...
//						System.out.println("#**LG.PS.gLSG root "+v
//								+"\tnfam "+profile_count
//								+"\tNuSv "+Math.exp(logNu_Sv)
//								+"\tSv "+Math.exp(logSv)
//								+"\t1_p "+Math.exp(log1_p)
//								+"\tdLLdp  "+Logarithms.ldiffValue(ldiff)
//								);
						
					} else {
						int u = factory.tree.getParent(v);
						if (factory.tree.getNumChildren(u)==2)
						{
							log1_e = factory.getLogLossComplement(factory.tree.getSibling(v));
						} else
						{
							log1_e = Logarithms.logToLogComplement(log_factory.getLogExtinction(u)-logp);
						}
						double log1_pe = Logarithms.add(log1_p, logp+log1_e);
						dpos = logNu_Sv+log1_p-log1_pe;
						dneg = logSv+logp+log1_e-log1_pe;
					}
					
//					if (Double.isNaN(dpos)||Double.isNaN(dneg)) {
//						// DEBUG
//						System.out.println("#**LG.gLSG\t"+v+"\tSv "+logSv+"\tNuSv "+logNu_Sv+"\tl1e "+log1_e
//							+"\ttNuSv "+Arrays.toString(tNu_Sv)
//								+"\t// "+log_factory.toString(v)+"\t// "+factory.rates.toString(v)
//						);
//					}
					
					logSD[3*v + PARAMETER_LOSS] = Logarithms.ldiff(dpos,dneg); 
					if (!Logarithms.ldiffIsFinite(logSD[3*v+PARAMETER_LOSS])) {
						// DEBUG
						System.out.println("#**LG.gLSG\t"+v+"\tSv "+logSv+"\tNuSv "+logNu_Sv+"\tl1e "+log1_e+"\t// "+log_factory.toString(v)+"\t// "+factory.rates.toString(v));
						// code dies here with 	Sv -Infinity	NuSv NaN
					}
					
					assert (Logarithms.ldiffIsFinite(logSD[3*v+PARAMETER_LOSS]));
					
				}
				
				
				double log1_q = log_factory.getLogDuplicationComplement(v);  
				double log_q = log_factory.getLogDuplicationParameter(v);
				
				double[] tNv_Sv = log_birth_tails[v];
				double logNv_Sv = Logarithms.sum(tNv_Sv, tNv_Sv.length);
//				System.out.println("#**LG.gLASG\t"+v+"\tSv "+logSv+"\tNvSv "+logNv_Sv+"\tF "+F);
				
				// gain
				if (log_q==Double.NEGATIVE_INFINITY) 
				{
					// Poisson
					double r = factory.getGainParameter(v);
					double log_r = factory.getLogGainParameter(v);
					if (log_r==Double.NEGATIVE_INFINITY)
					{
						logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff();
					} else
					{
						logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff(logNv_Sv, log_r+logF);
						logSD[3*v+PARAMETER_DUPLICATION] = Logarithms.ldiff();
						
						assert (Logarithms.ldiffIsFinite(logSD[3*v+PARAMETER_GAIN]));
					}
				} else
				{
					// Polya
					double κ = factory.getGainParameter(v);
					double log_kappa = factory.getLogGainParameter(v);
					if (log_kappa==Double.NEGATIVE_INFINITY)
					{
						logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff();
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
						double loglog1_q = Math.log(-log1_q);
						if (loglog1_q==Double.NEGATIVE_INFINITY) 
							loglog1_q = log_q;
						double log_kappa_log1_q = log_kappa+loglog1_q;
						
						double dneg = log_kappa_log1_q+logF;
						
						//dneg = (Math.log(-log1_q)+log_kappa)+logF; // +F*kappa*ln(1-q)

						logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff(dpos,dneg);
						assert (Logarithms.ldiffIsFinite(logSD[3*v+PARAMETER_GAIN]));
						
						
//						// DEBUG
//						System.out.println("#**LG.gLASG "+v+"/"+PARAMETER_DUPLICATION
////								+"\tNvSv "+logNv_Sv+"\tSv "+logSv
////								+"\tlq "+log_q+"\tl1q "+log1_q
//								+"\tdL "+Logarithms.ldiffValue(
//										Logarithms.ldiffMultiply(
//										logSD[3*v+PARAMETER_DUPLICATION],
//										-log_q-log1_q, null))
//								+"\tdLlogit "+
//								Logarithms.ldiffValue(logSD[3*v+PARAMETER_DUPLICATION]));
						
					}
					// now for duplication
					double dpos = logNv_Sv+log1_q;
					double dneg = Logarithms.add(logSv, log_kappa+logF) + log_q;

					logSD[3*v+PARAMETER_DUPLICATION] = Logarithms.ldiff(dpos,dneg);
					assert (Logarithms.ldiffIsFinite(logSD[3*v+PARAMETER_DUPLICATION]));
				}
			}		
			return logSD;
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

	
	/* ================================ */
	/* Calculating the gradient         */ 
	/* ================================ */

//	/**
//	 * Log-gradient 
//	 * by logit(p~), logit(q~), log(kappa)/log(r~)
//	 * 
//	 * @param S
//	 * @return
//	 */
//	public double[][] getLogSurvivalGradient(PosteriorStatistics S)
//	{
////		int num_nodes = factory.tree.getNumNodes();
////		double[][] logSD = new double[3*num_nodes][]; // return value
////		
////		double F = S.profile_count;
////		double logF = Math.log(F);
////		for (int v=0; v<num_nodes; v++)
////		{
////			double[] pSv = S.log_edge_posteriors[v].clone(); // clone bc we destroy it 
////			{ // calculating tails and then summing to get the mean 
////				double log_tail = Double.NEGATIVE_INFINITY;
////				int j=pSv.length;
////				while (0<j)
////				{
////					--j;
////					double x = pSv[j];
////					pSv[j] = log_tail;
////					log_tail = Logarithms.add(log_tail, x);
////				}
////			}
////			double logSv = Logarithms.sum(pSv, pSv.length); 
////			
////			// set loss-derivative
////			
////			//(1-p)*Nu_Sv - p*S
////			
////			double logp = log_factory.getLogLossParameter(v);
////			double log1_p = log_factory.getLogLossComplement(v);
////			
////			if (factory.tree.isRoot(v) || logp==0.0)
////			{
////				logSD[3*v + PARAMETER_LOSS] = Logarithms.ldiff(); 
////			} else
////			{
////				double[] tNu_Sv = S.log_death_tails[v];
////				double logNu_Sv = Logarithms.sum(tNu_Sv, tNu_Sv.length);
////				
////				int u = factory.tree.getParent(v);
////				double log1_e;
////				if (factory.tree.getNumChildren(u)==2)
////				{
////					log1_e = factory.getLogLossComplement(factory.tree.getSibling(v));
////				} else
////				{
////					log1_e = Logarithms.logToLogComplement(log_factory.getLogExtinction(u)-logp);
////				}
////				double dpos = logNu_Sv+log1_p;
////				double dneg = logSv+logp+log1_e;
////				logSD[3*v + PARAMETER_LOSS] = Logarithms.ldiff(dpos,dneg); 
////			}
////			
////			
////			double log1_q = log_factory.getLogDuplicationComplement(v);  
////			double log_q = log_factory.getLogDuplicationParameter(v);
////			
////			double[] tNv_Sv = S.log_birth_tails[v];
////			double logNv_Sv = Logarithms.sum(tNv_Sv, tNv_Sv.length);
//////			System.out.println("#**LG.gLASG\t"+v+"\tSv "+logSv+"\tNvSv "+logNv_Sv+"\tF "+F);
////			
////			// gain
////			if (log_q==Double.NEGATIVE_INFINITY) 
////			{
////				// Poisson
////				double r = factory.getGainParameter(v);
////				if (r==0.0)
////				{
////					logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff();
////				} else
////				{
////					double log_r = Math.log(r);
////					logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff(logNv_Sv, log_r+logF);
////					logSD[3*v+PARAMETER_DUPLICATION] = Logarithms.ldiff();
////				}
////			} else
////			{
////				// Polya
////				double κ = factory.getGainParameter(v);
////				double log_kappa = Math.log(κ);
////				if (κ==0.0)
////				{
////					logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff();
////				} else
////				{
////					double dpos = tNv_Sv[0]; // first term
////					for (int i=1; i<tNv_Sv.length; i++)
////					{
////						// ln(k/(k+i)) = ln(1/(1+i/k)) =-ln(1+i/k) for i<k
////						// = ln (k/i)-ln(1+k/i) for k<i
////						
////						double log_k_ki;
////						if (i<κ)
////						{
////							log_k_ki = -Math.log1p(i/κ);
////						} else
////						{
////							log_k_ki = log_kappa - Math.log(i) -Math.log1p(κ/i);
////						}
////						dpos = Logarithms.add(dpos, tNv_Sv[i]+log_k_ki);
////					}
////					double dneg = (Math.log(-log1_q)+log_kappa)+logF; // +F*kappa*ln(1-q)
////					logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff(dpos,dneg);
////					
////					// now for duplication
////					dpos = logNv_Sv+log1_q;
////					dneg = Logarithms.add(logSv, log_kappa+logF) + log_q;
////					logSD[3*v+PARAMETER_DUPLICATION] = Logarithms.ldiff(dpos,dneg);
////					
//////					// DEBUG
//////					System.out.println("#**LG.gLASG "+v+"/"+PARAMETER_DUPLICATION
////////							+"\tNvSv "+logNv_Sv+"\tSv "+logSv
////////							+"\tlq "+log_q+"\tl1q "+log1_q
//////							+"\tdL "+Logarithms.ldiffValue(
//////									Logarithms.ldiffMultiply(
//////									logSD[3*v+PARAMETER_DUPLICATION],
//////									-log_q-log1_q, null))
//////							+"\tdLlogit "+
//////							Logarithms.ldiffValue(logSD[3*v+PARAMETER_DUPLICATION]));
////					
////				}
////			}
////		}		
//		
//		double[][] logSD = S.getLogSurvivalGradient();
//		return logSD;
//	}
//
	/**
	 * Log-gradient by logit(p), logit(q), log(kappa)/log(r); 
	 * replaces entries of the survival gradient.
	 * 
	 * @param logSD survival gradient from {@link #correctedLogSurvivalGradient(PosteriorStatistics)}
	 * @return same logSD array, with updated entries
	 */ 
	public double[][] convertToLogDistributionGradient(double[][] logSD)
	{
		int num_nodes = factory.tree.getNumNodes();
		double[][] dde = new double[num_nodes][]; // d/dε
		
		int v=factory.tree.getRoot();
		while (0<=v)
		{
			double log_p = log_factory.getLogLossParameter(v);
			double log_e = log_factory.getLogExtinction(v);
			double log_q = log_factory.getLogDuplicationParameter(v);
			if (log_factory.rates.getLogLossComplement(v) ==Double.NEGATIVE_INFINITY) // || factory.tree.isRoot(v)  
			{
				if (log_factory.rates.getLogDuplicationParameter(v)==Double.NEGATIVE_INFINITY) // (q==0.0)
				{
					// Poisson
					double logr = factory.rates.getLogGainParameter(v); // .getGainParameter(v);
					if (logr==Double.NEGATIVE_INFINITY)
					{
						dde[v] = Logarithms.ldiff();
						logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiff();
					} else
					{
						double[] dLdr =  logSD[3*v+PARAMETER_GAIN]; 
						// logSD[3*v+PARAMETER_GAIN] = dLdr;   // stays the same
						dde[v] =  Logarithms.ldiffMultiply(dLdr, log_e, null);
						dde[v] = Logarithms.ldiffInverse(dde[v], dde[v]);
					}
					logSD[3*v+PARAMETER_DUPLICATION] = Logarithms.ldiff();
				} else
				{
					// Polya
					double[] dLdq = logSD[3*v+PARAMETER_DUPLICATION];
					// stays same 
					dde[v] = Logarithms.ldiffMultiply(dLdq, log_e, null);
					dde[v] = Logarithms.ldiffInverse(dde[v], dde[v]);
				}
				logSD[3*v+PARAMETER_LOSS]=Logarithms.ldiff();
			} else
			{
				int u=factory.tree.getParent(v);
				
				// TODO: how to deal with the root? 
				
				if (LOGIT_TO_LOGIT)
				{ // logit-to-logit: faster calculations
					double log1_p = log_factory.getLogLossComplement(v);
					
					
					double[] dLdp = logSD[3*v+PARAMETER_LOSS];

					
					double[] dLdpe;
					if (factory.tree.isRoot(v)) { // ROOTLOSS 
						dLdpe = dLdp;
					} else {
						double dedp = log1_p-log_factory.getLogExtinctionComplement(u); 
						dLdpe = Logarithms.ldiffAddMultiply(dLdp, dedp, dde[u], null);
					}
					double dpdp = log_factory.rates.getLogLossParameter(v)-log_p; 
					
					logSD[3*v+PARAMETER_LOSS] = Logarithms.ldiffMultiply(dLdpe, dpdp, null);
					
//					System.out.println("#**LG.gLDG "+v
//						+"\tdLdpe "+Logarithms.ldiffValue(Logarithms.ldiffMultiply(dLdpe, -log_p-log1_p, null))
//						+"\tdLdp~ "+Logarithms.ldiffValue(Logarithms.ldiffMultiply(dLdp, -log_p-log1_p, null))
//						+"\tdeu "+Logarithms.ldiffValue(Logarithms.ldiffMultiply(dde[u], -log_factory.getLogExtinction(u)-log_factory.getLogExtinctionComplement(u), null))
//						+"\tdLdp "+Logarithms.ldiffValue(Logarithms.ldiffMultiply(logSD[3*v+PARAMETER_LOSS], -log_factory.rates.getLogLossParameter(v)-log_factory.rates.getLogLossComplement(v), null))
//						);
					
					
					if (log_factory.rates.getLogDuplicationParameter(v)==Double.NEGATIVE_INFINITY) // (q==0.0)
					{
						// Poisson
						// double r = factory.rates.getGainParameter(v);
						double[] dLdr =  logSD[3*v+PARAMETER_GAIN];
						double dpde = log_e-log_p;
						
						dde[v] = Logarithms.ldiffSubtractMultiply(dpde, dLdpe, log_e, dLdr, null);
						
						// logSD[3*v+PARAMETER_GAIN] = dLdr; // stays the same
						logSD[3*v+PARAMETER_DUPLICATION] = Logarithms.ldiff();
					} else
					{
						double dpde = log_factory.getLogDuplicationComplement(v)+log_e-log_p;
						
						double dpdq = dpde + factory.rates.getLogDuplicationParameter(v);
						
						double[] dLdq = logSD[3*v+PARAMETER_DUPLICATION] ;
						logSD[3*v+PARAMETER_DUPLICATION] = Logarithms.ldiffSubtractMultiply(dLdq, dpdq, dLdpe, null);
						
						dde[v] = Logarithms.ldiffSubtractMultiply(dpde, dLdpe, log_e, dLdq, null);
						
					}
				} else
				{ // development code for debugging 
					// logit-to-parameter-to-parameter-to-logit
					double dedp = factory.tree.getNumChildren(u)==2
									?log_factory.getLogLossParameter(factory.tree.getSibling(v))
									:log_factory.getLogExtinction(u)-log_p;
					
					double log1_p = log_factory.getLogLossComplement(v);
					double log1_e = log_factory.getLogExtinctionComplement(v);
					double[] dLdp = Logarithms.ldiffMultiply(logSD[3*v+PARAMETER_LOSS],-log_p-log1_p, null);
					double[] dLdeu = Logarithms.ldiffMultiply(dde[u],-log_factory.getLogExtinction(u)-log_factory.getLogExtinctionComplement(u), null);
					double[] dLdpe = Logarithms.ldiffAddMultiply(dLdp, dedp, dLdeu, null); 
					double[] dLdq, dLdg, dLde;
					
					if (log_factory.rates.getLogDuplicationParameter(v)==Double.NEGATIVE_INFINITY) // (q==0.0)
					{
						// Poisson
						double logr = log_factory.rates.getLogGainParameter(v);// .getGainParameter(v);
						double[] dLdr = Logarithms.ldiffMultiply(logSD[3*v+PARAMETER_GAIN], -log_factory.getLogGainParameter(v), null);

//						double[] dLdr = Logarithms.ldiffMultiply(logSD[3*v+PARAMETER_GAIN], -Math.log(log_factory.getGainParameter(v)), null);
						
						dLdg = Logarithms.ldiffMultiply(dLdr, log1_e, null);
						dLdp = Logarithms.ldiffMultiply(dLdpe, log1_e, null);
						dLdq = logSD[3*v+PARAMETER_DUPLICATION];
						assert (Logarithms.ldiffIsZero(dLdq));
						dLde = Logarithms.ldiffSubtractMultiply(log_factory.rates.getLogLossComplement(v), dLdpe, logr, dLdr, null);
					} else
					{
						double log1_q = log_factory.getLogDuplicationComplement(v);

						dLdg = Logarithms.ldiffMultiply(logSD[3*v+PARAMETER_GAIN], -log_factory.getLogGainParameter(v), null);
//						dLdg = Logarithms.ldiffMultiply(logSD[3*v+PARAMETER_GAIN], -Math.log(log_factory.getGainParameter(v)), null);

						double dp_dp = log_factory.getLogLossComplement(v)-log_factory.rates.getLogLossComplement(v);						
						dLdp = Logarithms.ldiffMultiply(dLdpe, dp_dp, null);
						
//						System.out.println("#**LG.gLDG "+v+"\tdpdp "+Math.exp(dp_dp)+"\tdLdpe "+Logarithms.ldiffValue(dLdpe)
//								+"\tdeu "+Logarithms.ldiffValue(dLdeu)+"\tepsi "+Math.exp(dedp)
//								+"\tdLdp "+Logarithms.ldiffValue(dLdp));
						

						dLdq =  Logarithms.ldiffMultiply(logSD[3*v+PARAMETER_DUPLICATION], -log_q-log1_q, null);
						
						double d_de = log1_q + log1_p - log_factory.rates.getLogLossComplement(v)-log1_e;
						dLde = Logarithms.ldiffMultiply(Logarithms.ldiffSubtractMultiply(log_factory.rates.getLogLossComplement(v), dLdpe, log_factory.rates.getLogDuplicationParameter(v), dLdq, null), d_de, null);
						
//						System.out.println("#**LG.gLDG "+v+"\td_de "+d_de+"\tdLdq~ "+Logarithms.ldiffValue(dLdq)+"\tdLdpe "+Logarithms.ldiffValue(dLdpe)
//								+"\tdLde "+Logarithms.ldiffValue(dLde));
						
						
						
						double dq_dq = log_q+log1_q
							-log_factory.rates.getLogDuplicationParameter(v)-log_factory.rates.getLogDuplicationComplement(v);
						double log1pe = log_factory.rates.getLogLossComplement(v)+log_e;
						double[] dLdqd = Logarithms.ldiffSubtractMultiply(dq_dq, dLdq, dq_dq+log1pe, dLdpe, null);
						
//						System.out.println("#**LG.gLDG "+v+"\tdqdq "+Math.exp(dq_dq)
//								+"\tdLdq~ "+Logarithms.ldiffValue(dLdq)+"\tdLdpe "+Logarithms.ldiffValue(dLdpe)
//								+"\tp1e "+ Math.exp(log1pe)
//								+"\tdLdq "+Logarithms.ldiffValue(dLdqd));
						dLdq = dLdqd;
					}
					logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiffMultiply(dLdg, log_factory.rates.getLogGainParameter(v), null); 
//					logSD[3*v+PARAMETER_GAIN] = Logarithms.ldiffMultiply(dLdg, Math.log(log_factory.rates.getGainParameter(v)), null); 
					logSD[3*v+PARAMETER_LOSS] = Logarithms.ldiffMultiply(dLdp, log_factory.rates.getLogLossParameter(v)+log_factory.rates.getLogLossComplement(v), null); 
					if (log_factory.rates.getLogDuplicationParameter(v)==Double.NEGATIVE_INFINITY)
					{
						
					} else
					{
						logSD[3*v+PARAMETER_DUPLICATION] = Logarithms.ldiffMultiply(dLdq, log_factory.rates.getLogDuplicationParameter(v)+log_factory.rates.getLogDuplicationComplement(v), null); 
					}
					dde[v] = Logarithms.ldiffMultiply(dLde, log_e+log1_e, null);
				}
			}
//			{ // DEBUG
//				System.out.println("#**LG.gLDG "+v
//						+"\tdLdlp "+Logarithms.ldiffValue(logSD[3*v+PARAMETER_LOSS])
//						+"\tdLdlq "+Logarithms.ldiffValue(logSD[3*v+PARAMETER_DUPLICATION])
//						+"\tdLdlg "+Logarithms.ldiffValue(logSD[3*v+PARAMETER_GAIN])
//						+"\tdLdle "+Logarithms.ldiffValue(dde[v])
//					);				
//			}
			v--;
		}
		
		return logSD;
	}
	
	
	/**
	 * Conversion between log-gradient by logit(p), 
	 * and by logit(transient)
	 * 
	 * @param lrates model parameters
	 * @param logDλ log-gradient by logit(p)/logit(transient), logit(lambda), log(gain)
	 * @return logDλ updated: log-gradient by logit(transient)/logit(p), logit(lambda), log(gain)
	 */
	public static double[][] convertLogTransientGradient(TreeWithLogisticParameters lrates, double[][] logDλ)
	{
		IndexedTree tree = lrates.getTree();
		int num_nodes = tree.getNumNodes();
		
		for (int v=0; v<num_nodes; v++)
		{
			// indexes for parameters
			final int j_loss = 3*v + PARAMETER_LOSS;
			final int j_dup  = 3*v + PARAMETER_DUPLICATION;
			// input gradient 
			double[] dlp = logDλ[j_loss];
			double[] dlq = logDλ[j_dup];
			
			final double log_λ = lrates.getLogRelativeRate(v);
			if (log_λ != Double.NEGATIVE_INFINITY)
			{
				logDλ[j_dup] = Logarithms.ldiffAddMultiply(dlq, log_λ, dlp, dlq);
			}
			logDλ[j_loss] = Logarithms.ldiffInverse(dlp, dlp);
		}
		return logDλ;
	}
	
	
	/**
	 * Log-gradient by logit(p), logit(lambda), log(kappa)/log(r)
	 * 
	 * @param lrates model parameters
	 * @param log_Dpq log-gradient by logit(p), logit(q), log(kappa)/log(r) from {@link LogGradient#convertToLogDistributionGradient(double[][])}
	 * 	
	 * @return updated log_Dpq (at loss- and duplication-entries)
	 */
	public static double[][] convertToLogRelativeRateGradient(TreeWithLogisticParameters lrates,  double[][] log_Dpq)
	{
		IndexedTree tree = lrates.getTree();
		int num_nodes = tree.getNumNodes();
		
		for (int v=0; v<num_nodes; v++)
		{
			final int j_loss = 3*v + PARAMETER_LOSS;
			final int j_dup  = 3*v + PARAMETER_DUPLICATION;
			final double log1_λ = lrates.getLogRelativeComplement(v);
			
			if (lrates.getLogDuplicationParameter(v)==Double.NEGATIVE_INFINITY) //  (log1_λ==0.0) // Poisson
			{
				// nothing to do
			} else 
			{
				// indexes for parameters
				
				// distribution gradient 
				double[] dlp = log_Dpq[j_loss];
				double[] dlq = log_Dpq[j_dup];

				final double log1_p = lrates.getLogLossComplement(v);
				final double log1_q = lrates.getLogDuplicationComplement(v);

				
				if (log1_p != Double.NEGATIVE_INFINITY) // && !tree.isRoot(v) 
				{
					log_Dpq[j_loss] = Logarithms.ldiffAddMultiply(dlp, log1_p-log1_q, dlq, dlp);
				}

				if (log1_λ == Double.NEGATIVE_INFINITY) // lambda==1 gives log1_λ==-infty
				{
					 log_Dpq[j_dup] = Logarithms.ldiff();
				} else
				{
					log_Dpq[j_dup] = Logarithms.ldiffMultiply(dlq, log1_λ-log1_q, dlq);
				}
			}
			
			if (!Logarithms.ldiffIsFinite(log_Dpq[j_dup]) || !Logarithms.ldiffIsFinite(log_Dpq[j_loss])) {
				// DEBUG
				System.out.println("#**LG.cTLRRG node "+v+"\tloss "+Arrays.toString(log_Dpq[j_dup])+"\tdup "+Arrays.toString(log_Dpq[j_dup])
						+"\t// "+lrates.toString(v));
			}
			
			assert (!Logarithms.ldiffIsNaN(log_Dpq[j_dup]) && !Logarithms.ldiffIsNaN(log_Dpq[j_loss]));
			
		}
		return log_Dpq;
	}
	
	
	/**
	 * Log-gradient by logit(p), logit(lambda), log(gain).
	 * Universal gain: r=kappa*(-ln(1-q)); linear gain: r=kappa * q.
	 * 
	 * @param lrates model parameters
	 * @param log_Dpλ log-gradient by logit(p), logit(lambda), log(kappa)/log(r) from {@link #convertToLogRelativeRateGradient(TreeWithLogisticParameters, double[][])}
	 * @param gain_parameter_by one of PARAMETER_GAIN (r), PARAMETER_DUPLICATION (kappa/r), PARAMETER_LOSS (gamma) 
	 * @param universal_gain whether juniversal or linear gain
	 * @return updated log_Dpλ (at loss- and duplication-entries)
	 */
	public static double[][] convertRelativeRateToLogGainGradient(TreeWithLogisticParameters lrates,  double[][] log_Dpλ, int gain_parameter_by, boolean universal_gain)
	{
		IndexedTree tree = lrates.getTree();
		int num_nodes = tree.getNumNodes();
		
		if (gain_parameter_by == PARAMETER_DUPLICATION)
		{
			// keep kappa/r : nothing to do 
		} else // if 
		{
			// gradient by r or gamma 
			assert (gain_parameter_by == PARAMETER_GAIN)  || (gain_parameter_by == PARAMETER_LOSS);
			
			for (int v=0; v<num_nodes; v++)
			{
				// indexes for parameters
				final int j_gain = 3*v + PARAMETER_GAIN;
				final int j_loss = 3*v + PARAMETER_LOSS;
				final int j_dup  = 3*v + PARAMETER_DUPLICATION;
				
				// distribution gradient 
				double[] dlgn = log_Dpλ[j_gain];
				double[] dlp = log_Dpλ[j_loss];
				double[] dllm = log_Dpλ[j_dup];
				
				final double log1_p = lrates.getLogLossComplement(v);
				final double log1_λ = lrates.getLogRelativeComplement(v);

				if (lrates.getDuplicationParameter(v)==Double.NEGATIVE_INFINITY) //  (log1_λ==0.0) // Poisson
				{
					if (gain_parameter_by == PARAMETER_LOSS)
					{
						if (log1_p!=Double.NEGATIVE_INFINITY) // && !tree.isRoot(v)  
						{
							log_Dpλ[j_loss]  = Logarithms.ldiffAddMultiply(dlp, log1_p, dlgn, dlp);
						}
					}
				} else
				{
					final double log1_q = lrates.getLogDuplicationComplement(v);
					// Polya
					if (universal_gain)
					{
						final double log_q = lrates.getLogDuplicationParameter(v);
						final double log_b = log_q-Math.log(-log1_q);							

						// adjust derivative by logit(p)
						if (log1_p == Double.NEGATIVE_INFINITY) // || (tree.isRoot(v) 
						{
							// nothing to do 
						} else
						{						
							final double minus_dgn_dp;
							if (gain_parameter_by == PARAMETER_GAIN) 
							{
								// by r
								minus_dgn_dp = log1_p-log1_q + log_b;
							} else
							{
								// by gamma
								assert (gain_parameter_by == PARAMETER_LOSS);
								double log_b_1_q = log_b + Logarithms.logToLogComplement(log1_q-log_b);
								minus_dgn_dp = log1_p-log1_q + log_b_1_q;
							}
							log_Dpλ[j_loss] = Logarithms.ldiffSubtractMultiply(dlp, minus_dgn_dp, dlgn, dlp);
						}
						
						// adjust derivative by logit(lambda)
						if (log1_λ==Double.NEGATIVE_INFINITY) //lambda==1
						{
							// nothing to do
						} else
						{
							double minus_dgn_dλ = log1_λ-log1_q+log_b;

							log_Dpλ[j_dup] = Logarithms.ldiffSubtractMultiply(dllm, minus_dgn_dλ, dlgn, dllm);
						}
						
					} else // linear gain
					{
						// adjust derivative by logit(p)
						if (gain_parameter_by == PARAMETER_LOSS 
								// || tree.isRoot(v) 
								|| log1_p == Double.NEGATIVE_INFINITY )
						{					
							// nothing to do 
						} else
						{
							log_Dpλ[j_loss] = Logarithms.ldiffSubtractMultiply(dlp, log1_p, dlgn, dlp);
						}
						// adjust derivative by logit(lambda)
						if (log1_λ==Double.NEGATIVE_INFINITY) //lambda==1
						{
							// nothing to do
						} else
						{
							log_Dpλ[j_dup] = Logarithms.ldiffSubtractMultiply(dllm, log1_λ, dlgn, dllm);
						}
					} // linear gain
				} // Polya
			} // for v
		} // not duplication-bound
		return log_Dpλ;
	}
	
	
	/**
	 * Log-gradient by logit(p), logit(q), log(gain).
	 * Universal gain: r=kappa*(-ln(1-q)); linear gain: r=kappa * q.
	 * 
	 * @param lrates model parameters
	 * @param log_Dpq log-gradient by logit(p), logit(q), log(kappa)/log(r) from {@link #convertToLogDistributionGradient(double[][])}
	 * @param gain_parameter_by one of PARAMETER_GAIN (r), PARAMETER_DUPLICATION (kappa/r), PARAMETER_LOSS (gamma) 
	 * @param universal_gain whether juniversal or linear gain
	 * @return updated log_Dpq (at loss- and duplication-entries)
	 */
	public static double[][] convertDistributionToLogGainGradient(TreeWithLogisticParameters lrates,  double[][] log_Dpq, int gain_parameter_by, boolean universal_gain)
	{
		IndexedTree tree = lrates.getTree();
		int num_nodes = tree.getNumNodes();
		
		if (gain_parameter_by == PARAMETER_DUPLICATION)
		{
			// keep kappa/r : nothing to do 
		} else // if 
		{
			// gradient by r or gamma 
			assert (gain_parameter_by == PARAMETER_GAIN)  || (gain_parameter_by == PARAMETER_LOSS);
			
			for (int v=0; v<num_nodes; v++)
			{
				// indexes for parameters
				final int j_gain = 3*v + PARAMETER_GAIN;
				final int j_loss = 3*v + PARAMETER_LOSS;
				final int j_dup  = 3*v + PARAMETER_DUPLICATION;
				
				// distribution gradient 
				double[] dlgn = log_Dpq[j_gain];
				double[] dlp = log_Dpq[j_loss];
				double[] dlq = log_Dpq[j_dup];
				
				final double log1_p = lrates.getLogLossComplement(v);
				final double log1_q = lrates.getLogDuplicationComplement(v);
				
				
				// adjust derivative by logit(p)
				if (log1_p == Double.NEGATIVE_INFINITY) // ||tree.isRoot(v) 
				{
					// nothing to do 
				} else
				{						
					if (gain_parameter_by== PARAMETER_LOSS)
					{
						log_Dpq[j_loss]  = Logarithms.ldiffAddMultiply(dlp, log1_p, dlgn, dlp);
					}
				}
				// adjust derivative by logit(q)
				if (lrates.getLogDuplicationParameter(v)==Double.NEGATIVE_INFINITY) //(log1_q==0.0) // Poisson
				{
					// nothing to do
				} else
				{
					// Polya
					final double log_minus_dgn_dq;
					if (universal_gain)
					{
						double log_q = lrates.getLogDuplicationParameter(v);
						double log_beta = log_q-Math.log(-log1_q);	
						log_minus_dgn_dq = log_beta;
					} else
					{
						log_minus_dgn_dq = log1_q;
					}
					log_Dpq[j_dup] = Logarithms.ldiffSubtractMultiply(dlq, log_minus_dgn_dq, dlgn, dlq);
				}
			} // for v			
		}
		return log_Dpq;
	}
	
	/**
	 * Log-gradient by logit(p), logit(lambda), log(gamma).
	 * Universal gain: r=kappa*(-ln(1-q)); linear gain: r=kappa * q.
	 */
	public static double[][] convertDistributionFromLogGainGradient(TreeWithLogisticParameters lrates,  double[][] log_Dpq, int gain_parameter_by, boolean universal_gain)
	{
		if (gain_parameter_by == PARAMETER_LOSS)
		{
			// keep gamma : nothing to do 
		} else if (gain_parameter_by == PARAMETER_DUPLICATION)
		{
			log_Dpq = convertDistributionToLogGainGradient(lrates, log_Dpq, PARAMETER_LOSS, universal_gain);
		} else
		{
			// input by r , to gamma: only dldp is affected 
			assert (gain_parameter_by == PARAMETER_GAIN); 
			// same Jacobian as with lambda 
			
			
			
			
			log_Dpq = convertRelativeRateFromLogGainGradient(lrates, log_Dpq, gain_parameter_by, universal_gain);
		}
		return log_Dpq;
	}
	
	/**
	 * Log-gradient by logit(p), logit(lambda), log(gamma).
	 * Universal gain: r=kappa*(-ln(1-q)); linear gain: r=kappa * q.
	 * 
	 * @param lrates model parameters
	 * @param log_Dpλ log-gradient by logit(p), logit(lambda), log(gain) from {@link #convertRelativeRateToLogGainGradient(TreeWithLogisticParameters, double[][], int, boolean)}
	 * @param gain_parameter_by one of PARAMETER_GAIN (r), PARAMETER_DUPLICATION (kappa/r), PARAMETER_LOSS (gamma) 
	 * @param universal_gain whether universal or linear gain
	 * @return updated log_Dpλ (at loss- and duplication-entries)
	 */
	public static double[][] convertRelativeRateFromLogGainGradient(TreeWithLogisticParameters lrates,  double[][] log_Dpλ, int gain_parameter_by, boolean universal_gain)
	{
		IndexedTree tree = lrates.getTree();
		int num_nodes = tree.getNumNodes();
		
		
		
		if (gain_parameter_by == PARAMETER_LOSS)
		{
			// keep gamma : nothing to do 
		} else if (gain_parameter_by == PARAMETER_DUPLICATION)
		{
			log_Dpλ = convertRelativeRateToLogGainGradient(lrates, log_Dpλ, PARAMETER_LOSS, universal_gain);
		} else
		{
			// input by r 
			assert (gain_parameter_by == PARAMETER_GAIN); 
			for (int v=0; v<num_nodes; v++)
			{
				// indexes for parameters
				final int j_gain = 3*v + PARAMETER_GAIN;
				final int j_loss = 3*v + PARAMETER_LOSS;
				
				// distribution gradient 
				double[] dlgn = log_Dpλ[j_gain];
				double[] dlp = log_Dpλ[j_loss];
				
				final double log1_p = lrates.getLogLossComplement(v);
				if (log1_p!=Double.NEGATIVE_INFINITY) // && !tree.isRoot(v) 
				{
					log_Dpλ[j_loss]  = Logarithms.ldiffAddMultiply(dlp, log1_p, dlgn, dlp);
				}
			}
		} 

		return log_Dpλ;
	}
	
	
//	public static double[][] getLogParameterGradient(TreeWithLogisticParameters lrates,  double[][] log_Dpλ, int gain_parameter_by, boolean universal_gain)
//	{
//		IndexedTree tree = lrates.getTree();
//		for (int node=0; node<tree.getNumNodes(); node++)
//		{
//			double xp, xlm, xr;
//			xp = -lrates.getLogLossParameter(node)-lrates.getLogLossComplement(node);
//			xlm = -lrates.getLogRelativeRate(node)-lrates.getLogRelativeComplement(node);				
//			xr = -lrates.getLogGainParameter(node, gain_parameter_by, !universal_gain);
//			int j_loss = 3*node+PARAMETER_LOSS;
//			int j_dup = 3*node+PARAMETER_DUPLICATION;
//			int j_gain = 3*node + PARAMETER_GAIN;
//			
//			double[] dlgn = log_Dpλ[j_gain];
//			double[] dlp = log_Dpλ[j_loss];
//			double[] dllm = log_Dpλ[j_dup];
//			
//			log_Dpλ[j_gain] = Logarithms.ldiffMultiply(dlgn, xr, null);
//			log_Dpλ[j_loss] = Logarithms.ldiffMultiply(dlp, xp, null);
//			log_Dpλ[j_dup] = Logarithms.ldiffMultiply(dllm, xlm, null);
//		}
//		return log_Dpλ;
//	}
	
	/**
	 * Conversion from log-gradient to gradient 
	 * 
	 * @param log_G arbitrary size 
	 * @return
	 */
	public static double[] toGradient(double[][] log_G)
	{
		double[] dLL = new double[log_G.length];
		for (int j=0; j<dLL.length; j++)
		{
			double[] logD = log_G[j];
			
			dLL[j] = logD==null?0.0:Logarithms.ldiffValue(log_G[j]);
		}
		return dLL;
	}
	
	/**
	 * Corrected log-gradient 
	 * by logit(p), logit(q) and log(kappa)/log(q); updates 
	 * the log-likelihood for corrected log-likelihood
	 * 
	 * @param S sample statistics from {@link #getSampleStatistics()}; its .LL is updated
	 * @return
	 */
	protected double[][] correctedLogSurvivalGradient(PosteriorStatistics S )
	{
		PosteriorStatistics S0 = getUnobservedStatistics();
		double L0 = S0.LL;
		double F = S.profile_count; // sample size
		
		double log_corr0 = Logarithms.logToLogit(L0);
		
		
		double[][] logD = S.getLogSurvivalGradient();
		double[][] logD0 = S0.getLogSurvivalGradient();
		
		double[][] logC = new double[logD.length][];
		
		double log_Fcorr0 =log_corr0+Math.log(F);
		
		for (int j=0; j<logC.length; j++)
		{
			logC[j] = Logarithms.ldiffAddMultiply(logD[j], log_Fcorr0, logD0[j], null);
		}
		
		S.LL = S.LL-F*Logarithms.logitToLogComplement(log_corr0);
		return logC;
	}
	
	/**
	 * Gradient of the log-likelihood.
	 * 
	 * @param log_gradient log-gradient by logit(p), logit(q), log(kappa)/log(r)
	 * @return gradient by p, q, kappa/r
	 */
	private double[] getParameterGradient(double[][] log_gradient, boolean survival)
	{
//		double F = factory.table.getFamilyCount();
//		double logF = Math.log(F);
		
		int num_nodes = factory.tree.getNumNodes();
		double[] dLL = new double[3*num_nodes];
		
		for (int node=0; node<num_nodes; node++)
		{
			double xp, xq, xr;
			if (survival)
			{
				xp = -log_factory.getLogLossParameter(node)-log_factory.getLogLossComplement(node);
				xq = -log_factory.getLogDuplicationParameter(node)-log_factory.getLogDuplicationComplement(node);				
				xr = log_factory.getLogGainParameter(node);
//				xr = -Math.log(log_factory.getGainParameter(node));
			} else
			{
				xp = -factory.rates.getLogLossParameter(node)-factory.rates.getLogLossComplement(node);
				xq = -factory.rates.getLogDuplicationParameter(node)-factory.rates.getLogDuplicationComplement(node);				
				xr = -factory.rates.getLogGainParameter(node);
//				xr = -Math.log(factory.rates.getGainParameter(node));
			}
			int j = 3*node+PARAMETER_LOSS;
			
//			if (Double.isInfinite(xp)) // DEBUG
//			{
//				System.out.println("#**LG.cTG "+node+"\txp "+xp+"\tld "+Arrays.toString(average_log_gradient[j])+"\t// "+log_factory.rates.toString(node));
//			}
			dLL[j] = Logarithms.ldiffValue(Logarithms.ldiffMultiply(log_gradient[j], xp, null));
			j = 3*node+PARAMETER_DUPLICATION;
			dLL[j] = Logarithms.ldiffValue(Logarithms.ldiffMultiply(log_gradient[j], xq, null));
			j = 3*node+PARAMETER_GAIN;
			dLL[j] = Logarithms.ldiffValue(Logarithms.ldiffMultiply(log_gradient[j], xr, null));
		}
		return dLL;
	}
		

	
	/**
	 * Test code: computing corrected survival and distribution gradients 
	 * by this class and by {@link count.model.Gradient}.  
	 */
	private void mainmain()
	{
		java.io.PrintStream out = System.out;
		int num_nodes = factory.tree.getNumNodes();
		
		Gradient G = new Gradient(log_factory);
		G.setMinimumObservedCopies(this.getMinimumObservedCopies());
		G.setCalculationWidthThresholds(this.getCalculationWidthAbsolute(), this.getCalculationWidthRelative());
		
		PosteriorStatistics S = this.getSampleStatistics();
		double F = S.profile_count;
		double uncorrLL = S.LL;
		double[][] logD = correctedLogSurvivalGradient(S);
		double corrLL = S.LL;
		out.println("# Uncorrected likelihood "+uncorrLL+"\t// "+log_factory.getLL());
		out.println("# Corrected likelihood "+corrLL);
		
		out.println("# Survival parameter gradient (gain, loss, duplication order)");
		double[] dLL = getParameterGradient(logD, true);
		G.printGradient(out, dLL);
		// DEBUG: via Gradient
		out.println("# Alt survival gradient (GLD order) via "+G.getClass().getName());
		double[] A = G.getCorrectedGradient();
		G.printGradient(out, A);
		for (int j=0; j<3*num_nodes; j++)
		{
			double delta = A[j]-dLL[j];
			double rd = delta/Math.abs(dLL[j]);
			out.println("#**LG.mm par "+(j/3)+"/"+(j%3)+"\tdelta "+delta+"\trd "+rd);
		}
		out.println("# Model parameter gradient (GLD order)");
		logD = convertToLogDistributionGradient(logD);
		dLL = getParameterGradient(logD, false);
		G.printGradient(out, dLL);
		// DEBUG: via Gradient
		out.println("# Alt parameter gradient (GLD order) via "+G.getClass().getName());
		A = G.getDistributionGradient(A);
		G.printGradient(out, A);
		for (int j=0; j<3*num_nodes; j++)
		{
			double delta = A[j]-dLL[j];
			double rd = delta/Math.abs(dLL[j]);
			out.println("#**LG.mm par "+(j/3)+"/"+(j%3)+"\tdelta "+delta+"\trd "+rd);
		}
		
		out.println("# Model logistic/log gradient (GLD order)");
		for (int j=0; j<3*num_nodes; j++)
		{
			dLL[j] = Logarithms.ldiffValue(logD[j]);
		}
		G.printGradient(out, dLL);
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
		CommandLine cli = new CommandLine(args,  us);
		
		LogGradient G = new LogGradient(cli.getRates(), cli.getTable());
		
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
