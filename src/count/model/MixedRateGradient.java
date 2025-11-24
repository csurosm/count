package count.model;
/*
 * Copyright 2022 Mikl&oacute;s Cs&#369;r&ouml;s.
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
import count.io.RateVariationParser;
import count.matek.Logarithms;

import static count.Count.THREAD_UNIT_TASK;
import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LENGTH;
import static count.model.GLDParameters.PARAMETER_LOSS;


/**
 * Computing the likelihood and its base-rate gradient in a 
 * rate-multiplier model. 
 * 
 * @author csuros
 *
 */
public class MixedRateGradient implements Count.UsesThreadpool
{
	public MixedRateGradient(MixedRateModel.RateMultipliers mixed_model, ProfileTable table)
	{
		this.mixed_model = mixed_model;
		this.class_gradients = new Gradient[mixed_model.getNumClasses()];
		this.class_index = new int[class_gradients.length];
		Arrays.fill(class_index,-1);
		this.log_class_prob = new double[class_gradients.length]; // precomputed
		this.table = new UniqueProfileTable(table);
		this.num_classes = 0;
		this.cachedLL = 0.0;
//		pool = new ForkJoinPool(THREAD_PARALLELISM); // ForkJoinPool.commonPool();
//		System.out.println("#**MRG.gLL forkjoin par "+pool.getParallelism()+"\tnproc "+
//				Runtime.getRuntime().availableProcessors());
	}
	
		
	private final MixedRateModel.RateMultipliers mixed_model;
	private final UniqueProfileTable table;
	
//	/**
//	 * Same thread pool across different instantiations.
//	 * (Not likely to be an issue through GUI but 
//	 * in a big model-selection procedure, this class might be instantiated many times.)  
//	 */
//	private static final ForkJoinPool pool;
//	static 
//	{
//		if (THREAD_PARALLELISM>1)
//			pool = new ForkJoinPool(THREAD_PARALLELISM);			
//		else
//			pool = null;
//	}

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
	
	/**
	 * First cells for classes with non-0 class probability; empty cells up to @link #mixed_model} class count 
	 */
	private final Gradient[] class_gradients;
	/**
	 * Mapping from our class indices to {@link #mixed_model} class indices.
	 */
	private final int[] class_index;
	private final double[] log_class_prob;
	/**
	 * Non-0 classes. 
	 */
	private int num_classes;
	
	private int min_copies=1;
	
	/**
	 * Call if the underlying model changes.
	 */
	public void computeClasses()
	{
		Arrays.fill(class_index,-1);
		Arrays.fill(log_class_prob,Double.NEGATIVE_INFINITY);
		
		num_classes=0;
		for (int c = 0; c<class_gradients.length; c++)
		{
			double pc = mixed_model.getClassProbability(c);
			if (pc != 0.0)
			{
				TreeWithRates rates =mixed_model.getClassModel(c);
				
				Gradient grad = new Gradient(rates, table);
				grad.setCalculationWidthThresholds(threshold_width_absolute,  threshold_width_relative);
				
				class_gradients[num_classes] = grad;
				class_index[num_classes] = c;
				log_class_prob[num_classes] = Math.log(pc);
				num_classes++;
			}
		}
		cachedLL = 0.0;
		
	}
	
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
	
	
	private Profile getProfile(int family_idx) { return new Profile(family_idx);}
	
	class Profile
	{
		private Profile(int family_idx)
		{
			if (num_classes==0) computeClasses();
			dLL = new Gradient.Profile[num_classes];
			for (int c=0; c<num_classes; c++)
			{
				dLL[c] = class_gradients[c].getGradient(family_idx);
			}
		}
		
		private final Gradient.Profile[] dLL;
		
		double[][] getRateGradients()
		{
			double[][] rate_gradients = new double[num_classes][];
			for (int c=0; c<num_classes; c++)
			{
				Gradient G = class_gradients[c];
				double[] survival_gradient = dLL[c].getSurvivalGradient();
				rate_gradients[c] = G.getTotalRateGradient(G.getDistributionGradient(survival_gradient));
				//System.out.println("#**MRG.gRG "+c+"\trg "+Arrays.toString(rate_gradients[c]));
			}
			return rate_gradients;
		}
		
		double[] getClassLikelihoods()
		{
			double[] getClassLikelihoods = new double[num_classes];
			for (int c=0; c<num_classes; c++)
			{
				//Likelihood.Profile cLik = dLL[c].post.inside;
				
				
				double ll = dLL[c].post.getLogLikelihood();
				if (!Double.isFinite(ll) )
				{
//					synchronized(System.out)
//					{
//						Likelihood Lik = cLik.getOwner();
//						System.out.println("#**MRG.P.gCL family "+dLL[c].post.inside.family_idx+
//								"\tlik "+ll+"\t"+Arrays.toString(cLik.get()));
//	//					System.out.println(count.io.RateVariationParser.printRates(Lik.rates));
//						for (int node=0; node<Lik.tree.getNumNodes(); node++)
//						{
//							System.out.println("#**MRG.P.gCL "+dLL[c].post.inside.family_idx+"\tnlik "+node+"\t"+Arrays.toString(cLik.getNodeLikelihoods(node))
//									);
//							System.out.println("#**MRG.P.gCL "+dLL[c].post.inside.family_idx+"\telik "+node+"\t"+Arrays.toString(cLik.getEdgeLikelihoods(node))
//									+"\tgdist "+Lik.getGainDistribution(node)
//									+"\t"+Arrays.toString(Lik.getGainDistribution(node).getPointMassFunction(5)));
//							System.out.println("#**MRG.P.gCL rate "+node+"\t"+Lik.rates.toString(node));
//						}
//						Lik.printParameters(System.out);
//					}
					ll = -1e22; // instead of log(0)
				}
				getClassLikelihoods[c] = ll;
			}
			return getClassLikelihoods;
		}
		
		double getLogLikelihood()
		{
			return calculateLL(getClassLikelihoods());
		}
	}
	
	
	
	private double calculateLL(double[] classLL)
	{
		double[] sum=new double[classLL.length];
		
		for (int c=0; c<num_classes; c++)
		{
			sum[c] = classLL[c]+(log_class_prob[c]);//-log_class_prob[0]); // better if same-prob classes
		}
		double LL=Logarithms.sum(sum, num_classes); //log_class_prob[0]  
		return LL;
	}
	
	private double[] getEmptyClassLikelihoods()
	{
		if (num_classes==0) computeClasses();
		double[] getEmptyClassLikelihoods = new double[num_classes];
		for (int c=0; c<num_classes; c++)
			getEmptyClassLikelihoods[c] = class_gradients[c].factory.getEmptyLL();
		return getEmptyClassLikelihoods;
	}
	
	private double[] getSingletonClassLikelihoods()
	{
		if (num_classes==0) computeClasses();
		double[] getSingletonClassLikelihoods = new double[num_classes];
		for (int c=0; c<num_classes; c++)
			getSingletonClassLikelihoods[c] = class_gradients[c].factory.getSingletonLL();
		return getSingletonClassLikelihoods;
	}
	
	private double[] getUnobservedClassLikelihoods()
	{
		if (num_classes==0) computeClasses();
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
	
//	private double[][] getEmptyRateGradients()
//	{
//		if (num_classes==0) computeClasses();
//		double[][] rate_gradients = new double[num_classes][];
//		for (int c=0; c<num_classes; c++)
//		{
//			Gradient G = class_gradients[c];
//			double[] survival_gradient =G.getEmptySurvivalGradient();
//			rate_gradients[c] = G.getRateGradient(G.getDistributionGradient(survival_gradient));
//		}
//		return rate_gradients;
//	}
//	
//	private double[][] getSingletonRateGradients()
//	{
//		if (num_classes==0) computeClasses();
//		double[][] rate_gradients = new double[num_classes][];
//		for (int c=0; c<num_classes; c++)
//		{
//			Gradient G = class_gradients[c];
//			double[] survival_gradient =G.getSingletonSurvivalGradient();
//			rate_gradients[c] = G.getRateGradient(G.getDistributionGradient(survival_gradient));
//		}
//		return rate_gradients;
//	}
	
	private double[][] getUnobservedRateGradients()
	{
		if (num_classes==0) computeClasses();
		double[][] rate_gradients = new double[num_classes][];

		for (int c=0; c<num_classes; c++)
		{
			Gradient G = class_gradients[c];
			if (min_copies==0)
			{
				int num_nodes = G.factory.tree.getNumNodes();
				rate_gradients[c] = new double[4*num_nodes];
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
				rate_gradients[c] = G.getTotalRateGradient(G.getDistributionGradient(survival_gradient));
			}
		}
		return rate_gradients;
		
	}
	
	/**
	 * Likelihood for all-0 profile
	 * @return
	 */
	public double getEmptyLL()
	{
		double[] class_ll = getEmptyClassLikelihoods();
		double emptyLL = calculateLL(class_ll);
		
		if (Double.isNaN(emptyLL) || 0.<emptyLL)
		{
			System.out.println("#**MRG.gEL numerr emptyLL "+emptyLL+"\tnc "+num_classes+"\tclassLL "+Arrays.toString(class_ll));
		}
		return emptyLL;
	}
	
	/**
	 * Likelihood for single-gene profile
	 * @return
	 */
	public double getSingletonLL()
	{
		return calculateLL(getSingletonClassLikelihoods());	
	}
	
	/**
	 * Likelihood for unobserved (0 or  0 and 1) profile  
	 * @return
	 */
	public double getUnobservedLL()
	{
		double unobservedLL = Double.NEGATIVE_INFINITY;
		if (min_copies>0)
		{
			unobservedLL = getEmptyLL();
			if (min_copies>1)
			{
				assert (min_copies==2);
				unobservedLL = Logarithms.add(unobservedLL, getSingletonLL());
			}
		}
		return unobservedLL;
	}
	
	/**
	 * 
	 * 
	 * 
	 * @return
	 */
	public double getObservedLL()
	{
		double[] unobs = getUnobservedClassLikelihoods();
		double[] obs = new double[num_classes];
		for (int c=0; c<num_classes; c++)
		{
			// want log(1-exp(unobs))
			double p = -Math.expm1(unobs[c]);
			obs[c] = Math.log(p);		
//			
//			System.out.println("#**MRG.gOL c "+c+"\tunobs "+unobs[c]+"\tp0 "+Math.exp(unobs[c])+"\tobs "+obs[c]+"\tp "+p+"\t1-p "+(1.0-p));
		}
		return calculateLL(obs);
	}
	
	private double cachedLL=0.0;
	
	/**
	 * Calculates the log-likelihood, 
	 * with multi-threading if enabled. 
	 *  
	 * @return
	 */
	public double getLL()
	{
		if (cachedLL != 0.0) return cachedLL;
		int F = table.getFamilyCount();
		final ForkJoinPool pool = threadPool();
		final int unit_task = Count.unitTask(F);

		class PartialSum extends RecursiveTask<Double>
		{
			PartialSum(int min_idx, int max_idx)
			{
				this.max_idx = max_idx;
				this.min_idx = min_idx;
			}
			
			private final int min_idx;
			private final int max_idx;
			
			@Override
			public Double compute()
			{
				if (max_idx-min_idx>THREAD_UNIT_TASK)
				{
					int mid = (min_idx+max_idx)/2;
					PartialSum left = new PartialSum(min_idx, mid);
					left.fork();
					PartialSum right = new PartialSum(mid, max_idx);
					return right.compute()+left.join();
				} else
				{
					double LL=0.0;
					for (int f=min_idx; f<max_idx; f++)
					{
						Profile G = getProfile(f); // recalculate instead of saving them across all families
						double Gll = G.getLogLikelihood();
						LL += Gll*table.getMultiplicity(f);
					}
					return LL;
				}
			}
		}
		PartialSum bigjob = new PartialSum(0, F);
		
		double LL;
		try
		{
			if (num_classes==0) computeClasses(); // do not initialize from a forked thread
			if (F>unit_task)
			{
				LL = pool.invoke(bigjob);
			} else
			{
				LL = bigjob.compute();
			}
			cachedLL = LL;
		} catch (Throwable t)
		{
			// could be out of memory, or out of threads, or numerical error, or whatever else			
			throw new RuntimeException(t);  
		}
		
//		if (pool != null)
//		{
//			// thread task for a subset of families
//			
//			LL = pool.invoke(new PartialSum(0,F));
//		} else
//		{
//			for (int f=0; f<F; f++)
//			{
//				Profile G = getProfile(f); // recalculate instead of saving them across all families
//				double Gll = G.getLogLikelihood();
//				if (Double.isNaN(Gll))
//				{
//					System.out.println("#**MRG.gLL "+f+"\t"+Gll);
//				}
//				
//				LL += Gll*table.getMultiplicity(f);
//			}
//		}
//		cachedLL = LL;
		
		return LL;
	}
	
	/**
	 * Corrected log-likelihood, by minimum unobserved copies.
	 * @return
	 */
	public double getCorrectedLL()
	{
		double getCorrectedLL; // return value
		
		double LL = getLL(); 
		int nF = table.getTotalFamilyCount();
		double log_p_not0 = getObservedLL();
		getCorrectedLL = LL-nF*log_p_not0; //    Math.log(p_not0);

		{ // DEBUG code  
			double L0 = getUnobservedLL();
			// LL-F*log(1-exp(L0))
			double p_not0  = -Math.expm1(L0); // ok with L0== -infinity
//			System.out.println("#**MRG.cLL uncorr "+LL+"\tunobs "+L0+"\tp0 "+Math.exp(L0)+"\tp "+p_not0+"\tlogp "+log_p_not0+"\tlog(p) "+Math.log(p_not0));

			if (Double.isNaN(getCorrectedLL) || Double.isInfinite(getCorrectedLL))
			{
				System.out.println("#**MRG.gCLL numerr uncorr "+LL+"\tunobs "+L0+"\tp "+p_not0+"\tlogp "+log_p_not0+"\tlog(p) "+Math.log(p_not0)+"\tcorr "+getCorrectedLL+"\tnF "+nF+"\tclassprob "+Arrays.toString(log_class_prob)+"\temptyprob "+Arrays.toString(getEmptyClassLikelihoods()));
				
				double[] empty = getEmptyClassLikelihoods();
				assert (empty.length == num_classes);
				double[] terms = new double[num_classes];
				for (int c=0; c<num_classes; c++)
				{
					terms[c] = (log_class_prob[c]-log_class_prob[0])+empty[c];
				}
				
				double s = 0.0;
				
				for (int c=0; c<num_classes; c++)
				{
					double dt = terms[c];
					double dtz = Math.exp(dt);
					s += dtz;
					System.out.println("#**MRG.gCLL numerr c "+c+"\tterm "+terms[c]+"\tdt "+dt+"\tdtz "+dtz+"\ts "+s);			
				}
				
				double ss = Math.exp(log_class_prob[0])*s;
				
				double ls = Logarithms.sum(terms, num_classes);
				
				System.out.println("#**MRG.gCLL numerr ss "+s+"\tc0 "+Math.exp(log_class_prob[0])+"\tls "+ls+"\tels "+Math.exp(ls));
				
				
				
				double[] log_pnz = new double[num_classes];
				for (int c=0; c<num_classes; c++)
				{
					log_pnz[c] = -Math.expm1(empty[c]); // 1-e^L
					terms[c] = (log_class_prob[c]-log_class_prob[0])+log_pnz[c];
				}
				ls = Logarithms.sum(terms, num_classes);
				ss = Math.exp(ls + log_class_prob[0]);
	
				System.out.println("#**MRG.gCLL numerr ss.pnz "+ss+"\tls "+ls+"\tlpnz "+Arrays.toString(log_pnz)+"\tterms "+Arrays.toString(terms));
	
				
				System.out.println("#**MRG.cLL numerr obsLL "+getObservedLL());
			}
		}
		
		assert (!Double.isNaN(getCorrectedLL));
		
		return getCorrectedLL;
	}
	
	public double[][] getCorrectedClassGradients()
	{
		int num_nodes = mixed_model.getBaseModel().getTree().getNumNodes();
		
		double[] emptyLL = getUnobservedClassLikelihoods();
		double L0 = calculateLL(emptyLL);
		
		// exp(L0)/(1-exp(L0))
		double //corr0 = 0.5*Math.exp(L0/2.0)/Math.sinh(-L0/2.0); // Math.exp(L0)/(1-Math.exp(L0))
		corr0 = -Math.exp(L0)/Math.expm1(L0); // ==0 if L0==Double.NEGATIVE_INFINITY;
		
		double[][] emptyRD = getUnobservedRateGradients();
		double[] emptyD = getFullParameterGradient(emptyLL, emptyRD);
//		System.out.println("#*MRG.gCG\tL0 "+L0+"\tp0 "+Math.exp(L0)+"\tcorr0 "+corr0+"\temptyLL "+Arrays.toString(emptyLL)+"\temptyRD "+Arrays.toString(emptyRD[0])+"\temptyD "+Arrays.toString(emptyD));
		int nF = table.getFamilyCount();
		
		double classD[][] = new double[num_classes][3*num_nodes];
		
		for (int f=0; f<nF; f++)
		{
			Profile G = getProfile(f); // recalculate instead of saving them across all families
			double[] classLL = G.getClassLikelihoods();
			double[][] rate_gradients = G.getRateGradients();
			double LL = calculateLL(classLL);
			for (int c=0; c<num_classes; c++)
			{
				double wt = Math.exp(classLL[c]-LL)*table.getMultiplicity(f);
				for (int pi=0; pi<rate_gradients[c].length; pi++)
				{
					classD[c][pi]+=wt*rate_gradients[c][pi];
				}
			}
		}
		
//		for (int pi=0; pi<3*num_nodes; pi++)
//		{
//			int node = pi/3;
//			StringBuilder sb = new StringBuilder();
//			for (int c=0; c<num_classes; c++)
//			{
//				sb.append("\tclass"+c+"/"+node+"."+(pi%3)+"\t"+classD[c][pi]);
//			}
//			System.out.println("#**MRG.gCCG"+sb);
//		}
		
		return classD;
		
	}
	
		
	
	/**
	 * Gradient with 4 parameters per node (rates and length).
	 * @return
	 */
	public double[] getCorrectedGradient()
	{
		int num_nodes = mixed_model.getBaseModel().getTree().getNumNodes();
		
		double[] emptyLL = getUnobservedClassLikelihoods();
		double L0 = calculateLL(emptyLL);
		
		// exp(L0)/(1-exp(L0))
		double //corr0 = 0.5*Math.exp(L0/2.0)/Math.sinh(-L0/2.0); // Math.exp(L0)/(1-Math.exp(L0))
		corr0 = -Math.exp(L0)/Math.expm1(L0); // ==0 if L0==Double.NEGATIVE_INFINITY;
		
		double[][] emptyRD = getUnobservedRateGradients();
		double[] emptyD = getFullParameterGradient(emptyLL, emptyRD);
//		System.out.println("#*MRG.gCG\tL0 "+L0+"\tp0 "+Math.exp(L0)+"\tcorr0 "+corr0+"\temptyLL "+Arrays.toString(emptyLL)+"\temptyRD "+Arrays.toString(emptyRD[0])+"\temptyD "+Arrays.toString(emptyD));
		int nF = table.getFamilyCount();
		
		final ForkJoinPool pool = threadPool();
		final int unit_task = Count.unitTask(nF);
		
//		if (pool != null &&  nF>THREAD_UNIT_TASK)
//		{
			class PartialD extends RecursiveTask<double[]>
			{
				PartialD(int min_idx, int max_idx)
				{
					this.min_idx = min_idx;
					this.max_idx = max_idx;
				}
				private final int min_idx;
				private final int max_idx;
				
				@Override
				protected double[] compute()
				{
					double[] D;
					if (max_idx-min_idx>THREAD_UNIT_TASK)
					{
						int mid = (min_idx+max_idx)/2;
						PartialD left = new PartialD(min_idx, mid);
						left.fork();
						PartialD right = new PartialD(mid, max_idx);
						D = right.compute();
						double[] Dl = left.join();
						assert (D.length == Dl.length);
						for (int pidx=0; pidx<D.length; pidx++)
						{
							D[pidx]+=Dl[pidx];
						}
					} else
					{
						D = new double[4*num_nodes];
						assert (emptyD.length == D.length);
						
						for (int f=min_idx; f<max_idx; f++)
						{
							Profile G = getProfile(f); // recalculate instead of saving them across all families
							double[] classLL = G.getClassLikelihoods();
							double[][] rate_gradients = G.getRateGradients();
							double[] famD = getFullParameterGradient(classLL, rate_gradients);
							
							assert (famD.length == D.length);
							for (int pidx=0; pidx<D.length; pidx++)
							{
								D[pidx] += (famD[pidx]+corr0*emptyD[pidx])*table.getMultiplicity(f);
							}							
						}
					}
					return D;
				} // compute
			} // Task
			
			if (num_classes==0) computeClasses(); // do not initialize from a forked thread
			PartialD bigjob = new PartialD(0, nF);
			double[] getCorrectedGradient;
			try
			{
				if (nF>unit_task)
				{
					getCorrectedGradient = pool.invoke(bigjob);
				} else
				{
					getCorrectedGradient = bigjob.compute();
				}
			}catch (Throwable t)
			{
				// could be out of memory, or out of threads, or numerical error, or whatever else			
				throw new RuntimeException(t);  
			}
			return getCorrectedGradient;
//		} else
//		{
//			double[] D = new double[4*num_nodes];
//			
//			for (int f=0; f<nF; f++)
//			{
//				Profile G = getProfile(f); // recalculate instead of saving them across all families
//				double[] classLL = G.getClassLikelihoods();
//				double[][] rate_gradients = G.getRateGradients();
//				double[] famD = getFullParameterGradient(classLL, rate_gradients);
//				
//				assert (famD.length == D.length);
//				for (int pidx=0; pidx<D.length; pidx++)
//				{
//					D[pidx] += (famD[pidx]+corr0*emptyD[pidx])*table.getMultiplicity(f);
//				}
//			}
////			System.out.println("#**MRG.gCG "+corr0+"\t"+Arrays.toString(D));
//			
//			return D;
//		}
	}
	
		
	private double[] getFullParameterGradient(double[] classLL, double[][] rate_gradients)
	{
		
		double[] sum = new double[classLL.length];
		for (int c=0; c<sum.length; c++)
			sum[c]  = classLL[c]+log_class_prob[c];
		double LL = Logarithms.sum(sum,sum.length);
		
		IndexedTree tree = mixed_model.getBaseModel().getTree();
		int num_nodes = tree.getNumNodes();
		double[] D=new double[4*num_nodes]; // return value

		if (LL == Double.NEGATIVE_INFINITY)
		{
			return D; // gradients are 0.0 if impossible in all classes (unobserved gradient when min_copies==0)
		}

		for (int c=0; c<classLL.length; c++)
		{
			
			Gradient G = class_gradients[c];
			int model_class_idx = class_index[c];

			double wt = Math.exp(log_class_prob[c] + classLL[c]-LL);

			//System.out.println("#**MRG.gFPG class "+c+"("+model_class_idx+")"+"\twt "+wt);
			
			for (int node=0; node<num_nodes; node++)
			{
				
				double dLdκ = rate_gradients[c][3*node+PARAMETER_GAIN];
				double dLdμ = rate_gradients[c][3*node+PARAMETER_LOSS];
				double dLdλ = rate_gradients[c][3*node+PARAMETER_DUPLICATION];
				double t = G.factory.rates.getEdgeLength(node);
				
				if (!tree.isRoot(node))
				{
					if (t==Double.POSITIVE_INFINITY)
					{
						D[4*node+PARAMETER_GAIN] 
								+= wt * mixed_model.getGainRateMultiplier(model_class_idx) * dLdκ;
						D[4*node+PARAMETER_LOSS]
								+= wt * mixed_model.getLossRateMultiplier(model_class_idx) * dLdμ ;
						D[4*node+PARAMETER_DUPLICATION]
								+= wt * mixed_model.getDuplicationRateMultiplier(model_class_idx) * dLdλ ;
						D[4*node+PARAMETER_LENGTH] = 0.0;
					} else
					{
						D[4*node+PARAMETER_GAIN] 
								+= wt * mixed_model.getGainRateMultiplier(model_class_idx) * dLdκ;
						D[4*node+PARAMETER_LOSS]
								+= wt * mixed_model.getLossRateMultiplier(model_class_idx) * dLdμ * t;
						D[4*node+PARAMETER_DUPLICATION]
								+= wt * mixed_model.getDuplicationRateMultiplier(model_class_idx) * dLdλ * t;
						
						double μ = G.factory.rates.getLossRate(node);
						double λ = G.factory.rates.getDuplicationRate(node);
						
						D[4*node+PARAMETER_LENGTH]
								+= wt * mixed_model.getEdgeLengthMultiplier(model_class_idx) 
									* (dLdμ * μ + dLdλ * λ);
					}
				} else
				{
					// no multipliers at the root
					D[4*node+PARAMETER_GAIN] 
							+= wt * dLdκ;
					if (Double.isInfinite(t))
					{
						D[4*node+PARAMETER_LOSS]
								+= wt * dLdμ ;
						D[4*node+PARAMETER_DUPLICATION]
								+= wt * dLdλ ;		
						D[4*node+PARAMETER_LENGTH] = 0.0;
					} else
					{
						D[4*node+PARAMETER_LOSS]
								+= wt * dLdμ * t;						
						D[4*node+PARAMETER_DUPLICATION]
								+= wt * dLdλ * t;
						double μ = G.factory.rates.getLossRate(node);
						double λ = G.factory.rates.getDuplicationRate(node);
						
						D[4*node+PARAMETER_LENGTH]
								+= wt * (dLdμ * μ + dLdλ * λ);
					}
				}
			}
		}
		
		return D;
	}
	
	/**
	 * Inverts 4-parameter rate gradients into 3-parameter distribution gradients.
	 * Low-precision calculations prone to numerical error for extreme values of the parameters.
	 * 
	 * @param node
	 * @param param_type {@link GLDParameters#PARAMETER_GAIN}, {@link GLDParameters#PARAMETER_LOSS} or {@link GLDParameters#PARAMETER_DUPLICATION}
	 * @param gradient 4-parameter per node 
	 * @return partial derivative by the base model's distribution parameter for the given node
	 * @deprecated
	 */
	public double inferDistributionGradient(int node, int param_type, double[] gradient)
	{
		TreeWithRates base = mixed_model.getBaseModel();
		
		assert (gradient.length == 4*base.getTree().getNumNodes());
		// dLdμ dLdλ dLdt 
		
//		double dLdκ = gradient[4*node+PARAMETER_GAIN];
//		double dLdμ = gradient[4*node+PARAMETER_LOSS];
//		double dLdλ = gradient[4*node+PARAMETER_DUPLICATION];
//		double dLdt = gradient[4*node+PARAMETER_LENGTH];
		
		double t = base.getEdgeLength(node);	
		double μ = base.getLossRate(node);
		double λ = base.getDuplicationRate(node);
		
		double dLdθ;
		if (Double.isInfinite(t))
		{
//			assert (dLdμ==0.0);
			// Pólya or Poisson
			
			if (param_type == PARAMETER_GAIN)
				// nothing to do : r=γ for Poisson, or κ for Pólya
				dLdθ = gradient[4*node+PARAMETER_GAIN];
			else if (param_type == PARAMETER_LOSS) // p=1  
				dLdθ = 0.;
			else 
			{
				assert (param_type == PARAMETER_DUPLICATION); // q=λ/μ; p=1
				
				double q = base.getDuplicationParameter(node);
				double dLdλ = gradient[4*node+PARAMETER_DUPLICATION];
				double dLdμ = gradient[4*node+PARAMETER_LOSS];
				// λ = q*μ; μ = λ/q  
				double dμdq = -λ/q/q;
				double dλdq = μ;

				dLdθ = dLdλ * dλdq + dLdμ * dμdq;
			}
		}  else
		{
			double q = base.getDuplicationParameter(node);
			double p = base.getLossParameter(node);

			double dLdμ = gradient[4*node+PARAMETER_LOSS];
			double dLdt = gradient[4*node+PARAMETER_LENGTH];
			double dLdmut = dLdμ / t + dLdt / μ; 	

			if (q==0.0)
			{
				// Poisson
				double r = base.getGainParameter(node);
				double dmutdp = 1.0/base.getLossParameterComplement(node);
				if (r!=0.0)
				{
					double dLdγ = gradient[4*node+PARAMETER_GAIN];
					if (param_type == PARAMETER_GAIN)
						dLdθ = dLdγ / p;
					else if (param_type == PARAMETER_LOSS)
					{
						double dγdp = -r/p/p;
						dLdθ = dLdmut * dmutdp + dLdγ * dγdp;
					}
					else 
					{
						assert (param_type==PARAMETER_DUPLICATION);
						dLdθ = 0.0;
					}
				} else
				{
					if (param_type == PARAMETER_LOSS)
					{
						dLdθ = dLdmut * dmutdp;
					} else // DUP or GAIN 
					{
						dLdθ = 0.0;
					}
				}
			} else 				// Pólya
			{
				if (param_type == PARAMETER_GAIN)
				{
					dLdθ = gradient[4*node+PARAMETER_GAIN];
				} else 
				{
					if (μ == λ) // symmetric rates p=q=mu *t/(1+mu*t)
					{
						double p_1 = base.getLossParameterComplement(node);
						
						double dmutdp = 1.0/p_1/p_1;
						
						dLdθ = dLdmut * dmutdp; // same for lambda 
					} else 
					{
						double dLdλ = gradient[4*node+PARAMETER_DUPLICATION];
						double dLdlmt = dLdλ / t + dLdt / λ;
						
						double p_1 = base.getLossParameterComplement(node);
						double q_1 = base.getDuplicationParameterComplement(node);
						double p_q = p-q;
	
						double lnf = Math.log(q_1/p_1);
						
						if (param_type == PARAMETER_LOSS)
						{
							double z = q / p_q * lnf;
							double dmutdp = (p/p_1 -z)/p_q;
							double dlmtdp = (q/p_1 -z)/p_q;
	
							dLdθ = dLdmut * dmutdp + dLdlmt * dlmtdp; 	
						} else
						{
							double z = p / p_q * lnf;
							double dmutdq = (z-p/q_1)/p_q;
							double dlmtdq = (z-q/q_1)/p_q;
							
							dLdθ = dLdmut * dmutdq + dLdlmt * dlmtdq;
						}
					} // mu/lambda
				}
			} // q>0.0
		} // finite edge length

		return dLdθ;
	}
	
	protected void reportClassParameters()
	{
		int num_nodes = mixed_model.getBaseModel().getTree().getNumNodes();
		
		for (int c=0; c<num_classes; c++)
		{
			int class_idx = class_index[c];
			TreeWithRates class_model   = mixed_model.getClassModel(class_idx);
			for (int node=0; node<num_nodes; node++)
			{
				System.out.println("#**MRG.rCP class "+c+"("+class_idx+")\tnode "+node+"\t"+class_model.toString(node));
			}
		}
	}
	
	protected void reportRatesGradient(double[] dL)
	{
		int num_nodes = mixed_model.getBaseModel().getTree().getNumNodes();
		assert(dL.length == 4*num_nodes);
		
		for (int node =0; node<num_nodes; node++)
		{
			StringBuilder sb = new StringBuilder();
			double dLdμ = dL[4*node+PARAMETER_LOSS];
			double dLdt = dL[4*node+PARAMETER_LENGTH];
			double dLdλ = dL[4*node+PARAMETER_DUPLICATION];
			double dLdκ = dL[4*node+PARAMETER_GAIN];
			
			double glen = Math.sqrt(dLdμ*dLdμ+dLdt*dLdt+dLdλ*dLdλ+dLdκ*dLdκ);
			
			sb.append("dloss ").append(dLdμ)
			.append("\tdlen ").append(dLdt)
			.append("\tddup ").append(dLdλ)
			.append("\tdgain ").append(dLdκ);
			
			if (glen > 1e6) sb.append("\t***BIG***");
			
			System.out.println("#**MRG.rRG "+node+"\tglen "+glen+"\t"+sb.toString()+"\t// "+mixed_model.getBaseModel().toString(node));
		}
	}
	
	
	

	private void printSingletons(java.io.PrintStream out)
	{
        double L1 = getSingletonLL();
        double p1 = Math.exp(L1);
        int nF = table.getFamilyCount();
        double Exp_n1 = p1*table.getTotalFamilyCount();
        int n1=0;
        for (int f=0; f<nF; f++)
        {
        	int[] copies = table.getFamilyProfile(f);
        	int tot_copies=0;
        	for (int leaf=0; leaf<copies.length && tot_copies<2; leaf++)
        	{
        		tot_copies += copies[leaf];
        	}
        	if (tot_copies==1)
        		n1+=table.getMultiplicity(f);
        }
        out.println("Singletons L1 "+L1+" (p1 "+p1+")\tn1 "+n1+"\tExp "+Exp_n1);
	
	}
	
//	public static void main(String[] args) throws Exception
//	{
//		CommandLine cli = new CommandLine(args, MixedRateGradient.class);
////        if (args.length!=3)
////        {
////            System.err.println("Call as $0 phylogeny table rates");
////            System.exit(2008);
////        }
////        String tree_file = args[0];
////        String table_file = args[1];
////        String rates_file = args[2];
////        count.ds.Phylogeny tree = NewickParser.readTree(
////        		GeneralizedFileReader.guessReaderForInput(tree_file));
////        count.ds.AnnotatedTable table = TableParser.readTable(tree.getLeafNames(),
////        		GeneralizedFileReader.guessReaderForInput(table_file),true);
////        GammaInvariant input_model = RateVariationParser.readRates(
////        		GeneralizedFileReader.guessReaderForInput(rates_file)
////        		, tree);
//		count.ds.Phylogeny tree = cli.getTree();
//		count.ds.AnnotatedTable table = cli.getTable();
//		GammaInvariant input_model = cli.getGammaModel();
//        MixedRateGradient mixedG = new MixedRateGradient(input_model, table);
//        if (cli.getOptionValue(CommandLine.OPT_TRUNCATE)!=null)
//        {
//        	int absolute = cli.getOptionTruncateAbsolute();
//        	double relative = cli.getOptionTruncateRelative();
//        	mixedG.setCalculationWidthThresholds(absolute, relative);
//        	System.out.println(CommandLine.getStandardHeader("Truncated computation (absolute,relative)="
//        				+absolute+","+relative));
//        }
//        System.out.println(RateVariationParser.printRates(input_model));
//        
//        
//        
//        double LL = mixedG.getLL();
//        double L0 = mixedG.getEmptyLL();
//        double p0 = Math.exp(L0);
//        double p_not0 = -Math.expm1(L0);
//        double corrL = LL-table.getFamilyCount()*Math.log(p_not0);
//        System.out.println("log-likelihood "+LL+"\tL0 "+L0+"\tp0 "+p0+"\tcorr "+corrL);
//        mixedG.printSingletons(System.out);
//        
//        double[] D = mixedG.getCorrectedGradient();
//        
//        int num_nodes = tree.getNumNodes();
//    	System.out.println("#GRADIENT\tidx\tnode\tgain\tloss\tdup\tlength");
//        for (int node=0; node<num_nodes; node++)
//        {
//        	double dκ = D[4*node + PARAMETER_GAIN];
//        	double dμ = D[4*node + PARAMETER_LOSS];
//        	double dλ = D[4*node + PARAMETER_DUPLICATION];
//        	double dt = D[4*node + PARAMETER_LENGTH];
//			//sb.append(len+"\t"+drate+"\t"+lrate+"\t"+grate+"// "+tree.toString(node)+"\n");
//        	
//        	double dr = mixedG.inferDistributionGradient(node, PARAMETER_GAIN, D);
//        	double dp = mixedG.inferDistributionGradient(node, PARAMETER_LOSS, D);
//        	double dq = mixedG.inferDistributionGradient(node, PARAMETER_DUPLICATION, D);
//        	
//        	System.out.printf("#GRADIENT\t%d\t%s\t%g\t%g\t%g\t%g\t// dp %g\tdq %g\tdr %g\n", node, 
//        			tree.getNode(node).getFullName(), dκ, dμ, dλ, dt, dp, dq, dr);
//        }
//		
//	}
}
