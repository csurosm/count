
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

import count.Count;
import count.ds.ProfileTable;
import count.ds.UniqueProfileTable;
import count.io.CommandLine;
import count.io.NewickParser;
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
		
//		System.out.println("#**MRG.gCLL uncorr "+LL+"\tunobs "+L0+"\tp "+p_not0);
		LL -= nF*Math.log(p_not0);
		
		assert (!Double.isNaN(LL));
		return LL;
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
		
//		/**
//		 * Multiplicity of a 
//		 * family: 
//		 * if owner Gradient was instantiated
//		 * with a UniqueProfileTable, 
//		 * then the family's multiplicity, or else 1. 
//		 * 
//		 * @return
//		 */
//		int getMultiplicity()
//		{
//			return utable == null?1:utable.getMultiplicity(post.inside.family_idx);
//		}
		
		/**
		 * Gradient of the uncorrected log-likelihood.
		 * Index is 3*<var>node</var>+<var>i</var> where <var>i</var>={@link GLDParameters#PARAMETER_DUPLICATION},
		 * {@link GLDParameters#PARAMETER_GAIN}, {@link GLDParameters#PARAMETER_LOSS}.
		 * 
		 * @return gradient array
		 */
		double[] getSurvivalGradient()
		{
			// I. compute the posterior means and tails
			double[] node_means = post.getNodeMeans(); // array of expected values across the nodes
			double[] edge_means = post.getEdgeMeans(); 
			int num_nodes = factory.tree.getNumNodes();
			double[][] node_tails = new double[num_nodes][];
			double[][] edge_tails = new double[num_nodes][];		
			
			for (int node=0; node<num_nodes; node++)
			{
				double[] Ncdf = post.getNodeCDF(node);
				double[] Ntail 
					= node_tails[node] = new double[Ncdf.length];
				for (int ell=0; ell<Ncdf.length-1; ell++) // last entry is 1.0
					Ntail[ell]=1.0-Ncdf[ell];
				double[] Scdf = post.getEdgeCDF(node);
				double[] Stail 
					= edge_tails[node] = new double[Ntail.length];
				
				assert (Scdf.length <= Stail.length); 
				
				for (int s=0; s<Scdf.length-1; s++)
					Stail[s]=1.0-Scdf[s];
			}
			
			
			// II. compute the gradients 
			double[] dL = new double[num_nodes*3];// return value
			
			for (int node=0; node<num_nodes; node++)
			{
				double q = factory.getDuplicationParameter(node);
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
						
						for (int i=0; i<Ntail.length; i++)
						{
							dLdk += (Ntail[i]-Stail[i])/(κ + i);
						}
						double log1_q = Math.log1p(-q); //Math.log(1.0-q); // 
						dL[3*node+PARAMETER_GAIN] 
								= dLdk + log1_q;// d/dκ
						
						assert !Double.isNaN(dL[3*node+PARAMETER_GAIN] );
					}
					dL[3*node+PARAMETER_DUPLICATION]
							= (node_means[node]-edge_means[node])/q
							- (edge_means[node]+ κ)/(1.0-q); // d/dq
				}
				double p = factory.getLossParameter(node);
				if (p!=1.0)
				{
					if (factory.tree.isRoot(node))
					{
						dL[3*node+PARAMETER_LOSS]
								= (1.0-edge_means[node])/p - edge_means[node]/(1.0-p);
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
								= (node_means[parent]-omeu*edge_means[node])/p - omeu*edge_means[node]/(1.0-p);
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
			double q = factory.getDuplicationParameter(node);
			if (q==0.0)
			{
				// Poisson
				double r = factory.getGainParameter(node);
				if (r!=0.0)
				{
					dL[3*node+PARAMETER_GAIN] 
							=  - 1.0; // d/dr
				}
			} else
			{
				// Pólya
				double κ = factory.getGainParameter(node);
				if (κ!=0.0)
				{
					double log1_q = Math.log1p(-q); // Math.log(1.0-q); // Mat.log1p(-q)
					dL[3*node+PARAMETER_GAIN] 
							= log1_q;// d/dκ
				}
				dL[3*node+PARAMETER_DUPLICATION]
						= -κ/(1.0-q); // d/dq
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
		if (q==0.0)
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
			double a = 1.0-q*factory.extinction[node];
			double a2 = a*a;
			dL[3*node+PARAMETER_DUPLICATION] = dLdq * factory.getExtinctionComplement(node) //  (1.0-factory.extinction[node])
												/a2;
			// dL[3*node+PARAMETER_GAIN] does not change
			de[node] = -dLdq*q_1*q/a2;
		}
		assert Double.isFinite(dL[3*node+PARAMETER_DUPLICATION] );
		assert Double.isFinite(dL[3*node+PARAMETER_LOSS] );
		assert Double.isFinite(dL[3*node+PARAMETER_GAIN] );
		while (node>0)
		{
			--node;
			int parent = factory.tree.getParent(node);
			q = factory.rates.getDuplicationParameter(node);
			double p = factory.rates.getLossParameter(node);
			double ε = factory.extinction[parent]/factory.getLossParameter(node);
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
					de[node] = dLdpe * (1.0-p) - dLdr * r;
				}
			} else
			{
				// Pólya or geometric
				double q_1 = factory.rates.getDuplicationParameterComplement(node);
				double p_1 = factory.rates.getLossParameterComplement(node); 
				double a = 1.0-q*factory.extinction[node];
				dL[3*node+PARAMETER_LOSS] = dLdpe * epsi1 // (1.0-factory.extinction[node]) 
											/ a;
				double dLdq = dL[3*node+PARAMETER_DUPLICATION];
				double a2 = a*a;
				dL[3*node+PARAMETER_DUPLICATION] = (dLdq - dLdpe * p_1 * factory.extinction[node]) 
							* epsi1 //(1.0-factory.extinction[node])
							/a2;
				de[node] = (dLdpe * p_1 - dLdq * q) * q_1/a2;
			}
			assert Double.isFinite(dL[3*node+PARAMETER_DUPLICATION] );
			assert Double.isFinite(dL[3*node+PARAMETER_LOSS] );
			assert Double.isFinite(dL[3*node+PARAMETER_GAIN] );
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
	
	
	private void printGradient(java.io.PrintStream out, double[] dL)
	{
		for (int node =0; node<factory.tree.getNumNodes(); node++)
		{
			double dLdg = dL[3*node + PARAMETER_GAIN];
			double dLdp = dL[3*node + PARAMETER_LOSS];
			double dLdq = dL[3*node + PARAMETER_DUPLICATION];
			out.println(node+"\t"+dLdg+"\t"+dLdp+"\t"+dLdq);
		}
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
		out.println("# Survival parameter gradient (gain. loss, duplication order)");
		double[] D = getCorrectedGradient();
		printGradient(out, D);
		out.println("# Model parameter gradient (GLD order)");
		D = getDistributionGradient(D);
		printGradient(out, D);
		D = getTotalRateGradient(D);
		out.println("# Rate parameter gradient (GLD order)");
		printGradient(out, D);
	}
	
	public static void main(String[] args) throws Exception
	{
		count.io.CommandLine cli = new count.io.CommandLine(args, Gradient.class);
		
		
		Gradient G = new Gradient(cli.getRates(), cli.getTable());
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

		G.mainmain();
	}
	
}
