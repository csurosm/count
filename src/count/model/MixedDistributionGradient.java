package count.model;
/*
 * Copyright 2023 Mikl&oacute;s Cs&#369;r&ouml;s.
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

import static count.Count.THREAD_PARALLELISM;
import static count.Count.THREAD_UNIT_TASK;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;

import count.ds.AnnotatedTable;
import count.ds.Phylogeny;
import count.ds.ProfileTable;
import count.ds.UniqueProfileTable;
import count.io.CommandLine;
import count.matek.Logarithms;
import static count.io.CommandLine.OPT_TRUNCATE;
import static count.io.CommandLine.OPT_MINCOPY;

/**
 * Experimental code. 
 * 
 * @author csuros
 *
 */
public class MixedDistributionGradient 
{
	public MixedDistributionGradient(MixedRateModel mixed_model, ProfileTable table)	
	{
		this.model = mixed_model;
		if (table instanceof UniqueProfileTable)
		{
			this.utable = (UniqueProfileTable)table;
		} else
		{
			this.utable = new UniqueProfileTable(table);
		}

		this.min_copies = Integer.min(2, table.minCopies());
		this.num_classes = 0;
	}
	
	private final MixedRateModel model;
	private final UniqueProfileTable utable;
	
	/**
	 * Same thread pool across different instantiations.
	 * (Not likely to be an issue through GUI but 
	 * in a big model-selection procedure, this class might be instantiated many times.)  
	 */
	private static final ForkJoinPool thread_pool;
	static 
	{
		if (THREAD_PARALLELISM>1)
			thread_pool = new ForkJoinPool(THREAD_PARALLELISM);			
		else
			thread_pool = null;
	}
	
	/**
	 * Minimum total number of copies in observed profiles. 
	 */
	private int min_copies;
	
	public void setMinimumObservedCopies(int min_copies)
	{
		if (min_copies<0 || min_copies>2)
			throw new UnsupportedOperationException("Minimum copies must be 0,1 or 2.");
		this.min_copies = min_copies;
	}
	private int threshold_width_absolute = Integer.MAX_VALUE;
	private double threshold_width_relative = Double.POSITIVE_INFINITY;
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

	/**
	 * Number of non-zero probability classes.
	 */
	private int num_classes;
	/*
	 * Arrays filled in for non-zero classes/categories.
	 */
	/**
	 * Original index of this non-zero category.
	 */
	private int[] cat_og_idx;
	/**
	 * Log-probability of category.
	 */
	private double[] cat_log_probs;
	/**
	 * Gradient calculations for category
	 */
	private Gradient[] cat_gradients;
	
	
	public int getClassIndex(int cat)
	{
		return cat_og_idx[cat];
	}
	/**
	 * Sets up the data structures for the given model. 
	 */
	public void computeClasses()
	{
		if (num_classes != 0)
		{
			int nc = 0;
			for (int c=0; c<model.getNumClasses(); c++)
			{
				double pc = model.getClassProbability(c);
				if (pc!=0.0)
				{
					if (cat_og_idx[nc]!=c) // need to recompute 
					{
						num_classes = 0;
						break; 
					}
					nc++;
				}
			}
			if (nc!=num_classes)
				num_classes = 0; // need to recompute
		}
		
		if (num_classes==0)
		{
			int nc = model.getNumClasses();

			cat_log_probs= new double[nc];
			cat_gradients = new Gradient[nc];
			cat_og_idx = new int[nc];
	
			num_classes = 0;
			for (int c=0; c<model.getNumClasses(); c++)
			{
				double pc = model.getClassProbability(c);
				if (pc!=0.0)
				{
					TreeWithRates rates = model.getClassModel(c);
					cat_og_idx   [num_classes] = c;
					cat_log_probs[num_classes] = Math.log(pc);
					Gradient g = 
					cat_gradients[num_classes] = new Gradient(rates, utable);
					g.setCalculationWidthThresholds(threshold_width_absolute, threshold_width_relative);
					g.setMinimumObservedCopies(min_copies);
					
					num_classes++;
				}
			}
			cat_og_idx = Arrays.copyOf(cat_og_idx, num_classes);
			cat_log_probs = Arrays.copyOf(cat_log_probs, num_classes);
			cat_gradients = Arrays.copyOf(cat_gradients, num_classes);
		} else
		{
			// update class probabilities only; refresh gradient parameters 
			for (int c=0; c<num_classes; c++)
			{
				int og_c = cat_og_idx[c];
				double pc = model.getClassProbability(og_c);
				cat_log_probs[c]= Math.log(pc);
				cat_gradients[c].computeParameters();
			}
		}
		this.cachedLL = 0.0;
	}

	private double sumClassLL(double[] classLL)
	{
		assert (classLL.length == num_classes);
		double[] terms = new double[num_classes];
		for (int c=0; c<num_classes; c++)
		{
			terms[c] = cat_log_probs[c] + classLL[c];
		}
		double sum = Logarithms.sum(terms, num_classes);
		return sum;
	}
	
	/**
	 * Posterior probability for each nonzero class/category.
	 * 
	 * @param classLL
	 * @return
	 */
	private double[] posteriors(double[] classLL)
	{
		double[] posteriors = new double[num_classes];
		for (int c=0; c<num_classes; c++)
		{
			posteriors[c] = cat_log_probs[c] + classLL[c];
		}
		double LL = Logarithms.sum(posteriors, num_classes);
		for (int c=0; c<num_classes; c++)
		{
			posteriors[c]=Math.exp(posteriors[c]-LL);
		}
		return posteriors;
	}
	
	private double[] getEmptyClassLikelihoods()
	{
		if (num_classes == 0) computeClasses();
		double[] getEmptyClassLikelihoods = new double[num_classes];
		for (int c=0; c<num_classes; c++)
			getEmptyClassLikelihoods[c] = cat_gradients[c].factory.getEmptyLL();
		return getEmptyClassLikelihoods;
	}
	
	private double[] getSingletonClassLikelihoods()
	{
		if (num_classes == 0) computeClasses();
		double[] getSingletonClassLikelihoods = new double[num_classes];
		for (int c=0; c<num_classes; c++)
			getSingletonClassLikelihoods[c] = cat_gradients[c].factory.getSingletonLL();
		return getSingletonClassLikelihoods;
	}
	
	private double[] getUnobservedClassLikelihoods()
	{
		double[] unobserved;
		if (min_copies==0)
		{
			unobserved = new double[num_classes];
			Arrays.fill(unobserved, Double.NEGATIVE_INFINITY);
		} else 
		{
			unobserved = getEmptyClassLikelihoods();
			if (min_copies>1)
			{
				assert (min_copies==2);

				double[] singleL = getSingletonClassLikelihoods();
				for (int c=0; c<num_classes; c++)
					unobserved[c] = Logarithms.add(unobserved[c], singleL[c]);
			}
		}
		return unobserved;
	}
	
	private double[][] getUnobservedDistributionGradients()
	{
		if (num_classes == 0) computeClasses();
		
		double[][] distr_gradients = new double[num_classes][];

		for (int c=0; c<num_classes; c++)
		{
			Gradient G = cat_gradients[c];
			if (min_copies==0)
			{
				int num_nodes = G.factory.tree.getNumNodes();
				distr_gradients[c] = new double[3*num_nodes];
			} else 
			{
				assert (min_copies>=1);
				double[] survival_gradient =G.getEmptySurvivalGradient();
				if (min_copies>1)
				{
					double[] singleton_gradient = G.getSingletonSurvivalGradient();
					assert (min_copies==2);
					double L0 = G.factory.getEmptyLL();
					double L1 = G.factory.getSingletonLL();
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
				distr_gradients[c] = G.getDistributionGradient(survival_gradient);
			}
		}
		return distr_gradients;
	}
	
	
	class Profile
	{
		private Profile(int family_idx)
		{
			dLL = new Gradient.Profile[num_classes];
			for (int c=0; c<num_classes; c++)
			{
				dLL[c] = cat_gradients[c].getGradient(family_idx);
			}
		}
		
		private final Gradient.Profile[] dLL;
		
		double[] getClassLL()
		{
			double[] getClassLL = new double[num_classes];
			for (int c=0; c<num_classes; c++)
			{
				getClassLL[c]=dLL[c].post.inside.getLogLikelihood();
			}
			return getClassLL;
		}
		
		double getLL()
		{
			return sumClassLL(getClassLL());
		}

		double[][] getDistributionGradients()
		{
			double[][] distr_gradients = new double[num_classes][];
			for (int c=0; c<num_classes; c++)
			{
				Gradient G = cat_gradients[c];
				double[] survival_gradient = dLL[c].getSurvivalGradient();
				distr_gradients[c] = G.getDistributionGradient(survival_gradient);
			}
			return distr_gradients;
		}
		
		double[] getClassPosteriors()
		{
			return posteriors(getClassLL());
		}
	}


	private double cachedLL=0.0;
	
	/**
	 * Calculates the uncorrected log-likelihood, 
	 * with multi-threading if enabled. 
	 *  
	 * @return
	 */
	public double getLL()
	{
		if (num_classes == 0) computeClasses();
		if (cachedLL != 0.0) return cachedLL;

		int nF = utable.getFamilyCount();
		final int unit_task;
		if (thread_pool != null)
		{
			unit_task = THREAD_UNIT_TASK;
		}
		else
		{
			unit_task = nF; // do not fork
		}
		
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
							Profile P = new Profile(f);
							double fLL = P.getLL();
							LL += fLL * utable.getMultiplicity(f);
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
					
		this.cachedLL = LL;
		return LL;
	}

	public double getCorrectedLL()
	{
		double LL = getLL(); 
		// true number of families 
		int nT = utable.getTotalFamilyCount();
		double[] unobsLL = getUnobservedClassLikelihoods();
		double L0 = sumClassLL(unobsLL);
		// LL-F*log(1-exp(L0))
		double p_obs  = -Math.expm1(L0); // 1-exp(L0); ok with L0== -infinity
		
		LL -= nT*Math.log(p_obs);
		
		assert (!Double.isNaN(LL));
		return LL;		
	}
	
	public double[][] getCorrectedGradient()
	{
		if (num_classes == 0) computeClasses();
		double[] unobsLL = getUnobservedClassLikelihoods();
		double[] unobsP = posteriors(unobsLL);
		
		
		double L0 = sumClassLL(unobsLL);
		
		double //corr0 = 0.5*Math.exp(L0/2.0)/Math.sinh(-L0/2.0); // Math.exp(L0)/(1-Math.exp(L0))
		corr0 = -Math.exp(L0)/Math.expm1(L0); // ==0 if L0==Double.NEGATIVE_INFINITY;
		
		double[][] unobsD = getUnobservedDistributionGradients();
		
		int nF = utable.getFamilyCount();
		final int unit_task;
		if (thread_pool != null)
		{
			unit_task = THREAD_UNIT_TASK;
		}
		else
		{
			unit_task = nF; // do not fork
		}
		
		class PartialD extends RecursiveTask<double[][]>
		{
			PartialD(int min_idx, int max_idx)
			{
				this.min_idx = min_idx;
				this.max_idx = max_idx;
			}
			private final int min_idx;
			private final int max_idx;
			
			@Override
			protected double[][] compute()
			{
				try
				{
					double[][] D;
					if (max_idx-min_idx>unit_task)
					{
						int mid = (min_idx+max_idx)/2;
						PartialD left = new PartialD(min_idx, mid);
						left.fork();
						PartialD right = new PartialD(mid, max_idx);
						D = right.compute();
						double[][] Dl = left.join();
						assert (D.length == Dl.length); // == num_classes
						assert (D.length == num_classes);
						for (int c=0; c<num_classes; c++)
						{
							assert (D[c].length == Dl[c].length);
							for (int pidx=0; pidx<D[c].length; pidx++)
							{
								D[c][pidx]+=Dl[c][pidx];
							}
						}
					} else
					{
						D = new double[num_classes][];
						for (int c=0; c<num_classes; c++)
						{
							D[c] = new double[3*cat_gradients[c].factory.tree.getNumNodes()];
						}
						
						for (int f=min_idx; f<max_idx; f++)
						{
							Profile G = new Profile(f); // recalculate instead of saving them across all families
							double[][] famD = G.getDistributionGradients();
							
							double[] famP = posteriors(G.getClassLL());
							
							assert (famD.length == D.length);
							for (int c=0; c<num_classes; c++)
							{
								for (int pidx=0; pidx<D[c].length; pidx++)
								{
									D[c][pidx] += (famP[c]*famD[c][pidx]+corr0*unobsP[c]*unobsD[c][pidx])*utable.getMultiplicity(f);
								}							
							}
						}
					}
					return D;
				} catch (Throwable t)
				{
					throw t;
				}
			} // compute
		} // Task
		
		PartialD bigjob = new PartialD(0, nF);
		double[][] D;
		try
		{
			if (nF > unit_task)
			{
				D = thread_pool.invoke(bigjob);
			} else
			{
				D = bigjob.compute(); // no threads
			}
		} catch (Throwable t)
		{
			throw new RuntimeException(t);  			
		}
		return D;
	}
	
	public double[] getPosteriorFamilyCount()
	{
		if (num_classes == 0) computeClasses();

		double[] count = new double[num_classes];
		
		int nF = utable.getFamilyCount();
		final int unit_task;
		if (thread_pool != null)
		{
			unit_task = THREAD_UNIT_TASK;
		}
		else
		{
			unit_task = nF; // do not fork
		}
		
		class PartialP extends RecursiveTask<double[]>
		{
			PartialP(int min_idx, int max_idx)
			{
				this.min_idx = min_idx;
				this.max_idx = max_idx;
			}
			private final int min_idx;
			private final int max_idx;
			@Override
			protected double[] compute() 
			{
				try
				{
					double[] C;
					if (max_idx-min_idx>unit_task)
					{
						int mid = (min_idx+max_idx)/2;
						PartialP left = new PartialP(min_idx, mid);
						left.fork();
						PartialP right = new PartialP(mid, max_idx);
						C = right.compute();
						double[] Cl = left.join();
						assert (C.length == Cl.length);
						assert (C.length == num_classes);
						for (int c=0; c<num_classes; c++)
						{
							C[c] += Cl[c];
						}
					} else
					{
						C = new double[num_classes];
	
						for (int f=min_idx; f<max_idx; f++)
						{
							Profile G = new Profile(f); // recalculate instead of saving them across all families
							double[] classLL = G.getClassLL();
							double[] p = posteriors(classLL);
							int m = utable.getMultiplicity(f);
							for (int c=0; c<num_classes; c++)
							{
								C[c] += m*p[c];
							}
						}
					}
					return C;		
				} catch (Throwable t)
				{
					throw t;
				}
			}
		}	
		PartialP bigjob = new PartialP(0, nF);
		double[] C;
		try
		{
			if (nF > unit_task)
			{
				C = thread_pool.invoke(bigjob);
			} else
			{
				C = bigjob.compute(); // no threads
			}
		} catch (Throwable t)
		{
			throw new RuntimeException(t);  			
		}
		return C;		
	}
	
	public double[] getEmpiricalClassProbabilities()
	{
		double[] emptyLL = getUnobservedClassLikelihoods();
		double L0 = sumClassLL(emptyLL);
		double[] post0 = posteriors(emptyLL);
		double p_unobs = Math.exp(L0);
		double p_obs = -Math.expm1(L0); // 1-exp(L0)
		
		double[] C = getPosteriorFamilyCount();
		int nT = utable.getTotalFamilyCount();
		
		double[] p = new double[num_classes];
		for (int c=0; c<num_classes; c++)
		{
			double obs = C[c]/nT;
			double unobs = post0[c];
			p[c] = p_obs * obs + p_unobs * unobs;
			System.out.println("#**MDG.gELP "+c+"\tobs "+obs+"\tp "+p[c]+"\tpunobs "+p_unobs+"\tunobs "+unobs);
		}
		return p;
	}
	
	public static void main(String[] args) throws Exception
	{
		CommandLine cli = new CommandLine(args, MixedDistributionGradient.class);
		Phylogeny tree = cli.getTree();
		AnnotatedTable table = cli.getTable();
		GammaInvariant input_model = cli.getModel();
		
		MixedDistributionGradient mixedG = new MixedDistributionGradient(input_model, table);
        if (cli.getOptionValue(OPT_TRUNCATE)!=null)
        {
        	int absolute = cli.getOptionTruncateAbsolute();
        	double relative = cli.getOptionTruncateRelative();
        	mixedG.setCalculationWidthThresholds(absolute, relative);
        	System.out.println(CommandLine.getStandardHeader("Truncated computation: -"+OPT_TRUNCATE+" "
        				+absolute+","+relative));
        }
        if (cli.getOptionValue(OPT_MINCOPY)!=null)
        {
        	int min_copies = cli.getOptionInt(OPT_MINCOPY, table.minCopies());
        	mixedG.setMinimumObservedCopies(min_copies);
        	System.out.println(CommandLine.getStandardHeader("Minimum observed: -"+OPT_MINCOPY+" "
    				+min_copies));
        }
        double corrLL = mixedG.getCorrectedLL();
        System.out.println("Corrected log-likelihood "+corrLL);
        
        double[][] corrD = mixedG.getCorrectedGradient();
        
        double[] post = mixedG.getEmpiricalClassProbabilities();;
        
        for (int c=0; c<post.length; c++)
        {
        	int og_c = mixedG.getClassIndex(c);
        	System.out.println("Class "+c+"/"+og_c
        				+"\tprior "+input_model.getClassProbability(og_c)
        				+"\temp "+post[c]);
        }
		
	}
}
