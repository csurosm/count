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


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;
import java.util.function.DoubleFunction;

import count.Count;
import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.ds.UniqueProfileTable;
import count.io.CommandLine;
import count.io.RateVariationParser;
import count.matek.FunctionMinimization;
import count.matek.Functions;
import count.matek.Logarithms;
import count.model.StraightLikelihood.Profile;

import static count.io.CommandLine.OPT_EPS;
import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_MODEL_ROOT_PRIOR;
import static count.io.CommandLine.OPT_OUTPUT;
import static count.io.CommandLine.OPT_PVALUE;
import static count.io.CommandLine.OPT_ROUNDS;
import static count.io.CommandLine.OPT_TRUNCATE;


/**
 * Expectation-maximization algorithm with ancestor copies, 
 * using {@link StraightLikelihood}. 
 */
public class StraightEM extends ML implements GLDParameters, Count.UsesThreadpool
{
	private static boolean PRINT_OPTIMIZATION_MESSAGES = false;
	private static final double MAX_GAIN_RATE = 33.0; 
	private static final double MIN_GAIN_PARAMETER = 1.0/(1L<<40);
	
	/**
	 * Does not work well with universal gain param
	 */
	private static boolean OPT_UNIVERSAL_GAIN_PARAMETER = false;
	
	
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
	
	public StraightEM(TreeWithRates rates, ProfileTable table)
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
		this.optimize_node = new boolean[factory.tree.getNumNodes()];
		Arrays.fill(optimize_node, true);
		this.fixed_parameter=new boolean[3*optimize_node.length];
		final int root = factory.tree.getRoot();
		this.fixLoss(root, true);
		// keep dup=0.0 if Poisson
		this.fixDuplication(root, factory.getLogitDuplicationParameter(root)==Double.NEGATIVE_INFINITY);
			
	}
	
	private final StraightLikelihood factory;
	private final boolean[] optimize_node;	
	private final boolean[] fixed_parameter; 
	
	/**
	 * Same as the instantiating table, if it 
	 * was a UniqueProfileTable; or else its unique-profile version. 
	 */
	private final UniqueProfileTable utable;
	
	private int min_copies;
	 
	private boolean is_duprate_bounded = true; // if bounded, <1.0 is enforced
	
	/* ========================= */
	/* Setter and getter methods */ 
	/* ========================= */
	
	@Override
	public void setMinimumObservedCopies(int min_copies)
	{
		if (min_copies<0 || min_copies>2)
			throw new UnsupportedOperationException("Minimum copies must be 0,1 or 2.");
		this.min_copies = min_copies;
	}
	@Override
	public void fixNodeParameters(int node, boolean do_not_optimize)
	{
		this.optimize_node[node] = !do_not_optimize;
		fixNodeParameter(node, PARAMETER_GAIN, do_not_optimize);
		fixNodeParameter(node, PARAMETER_LOSS, do_not_optimize);
		fixNodeParameter(node, PARAMETER_DUPLICATION, do_not_optimize);
	}
	
	private void fixNodeParameter(int node, int param_type, boolean do_not_optimize)
	{
		this.fixed_parameter[3*node+param_type]=do_not_optimize;
		this.optimize_node[node] = 
				!this.isFixedParameter(node,PARAMETER_GAIN)
				|| !this.isFixedParameter(node,PARAMETER_LOSS)
				|| !this.isFixedParameter(node,PARAMETER_DUPLICATION);
	}
	private boolean isFixedParameter(int node, int param_type)
	{
		return this.fixed_parameter[3*node+param_type];
	}
	
	public void fixGain(int node, boolean not_optimized)
	{
		this.fixNodeParameter(node,PARAMETER_GAIN,not_optimized);
	}
	public void fixLoss(int node, boolean not_optimized)
	{
		this.fixNodeParameter(node,PARAMETER_LOSS,not_optimized);
	}
	public void fixDuplication(int node, boolean not_optimized)
	{
		this.fixNodeParameter(node,PARAMETER_DUPLICATION,not_optimized);
	}
	
	
	@Override
	public int getModelParameterCount()
	{
		IndexedTree tree = factory.tree;
		int np = 0;
		int u=tree.getNumNodes();
		while (u>0)
		{
			--u;
			if (optimize_node[u])
			{
				if (factory.getLogitDuplicationParameter(u)!=Double.NEGATIVE_INFINITY)//   rates.getDuplicationRate(u)!=0.0)
					np++; // dup
				if (!tree.isRoot(u)) np++; // loss
				if (factory.getGainParameter(u)!=0.0) //    rates.getGainRate(u)!=0.0)
					np++; // gain
			}
		}
		return np;
	}

	@Override
	public void setCalculationWidth(int absolute, double relative)
	{
		factory.setCalculationWidthThresholds(absolute, relative);
	}
	
	@Override
	public double getGainParameter(int node)
	{
		return factory.getGainParameter(node);
	}
	
	@Override
	public double getLossParameter(int node)
	{
		return factory.getLossParameter(node);
	}
	
	@Override
	public double getDuplicationParameter(int node)
	{
		return factory.getDuplicationParameter(node);
	}
	
	/* ========================= */
	/* E and M steps             */ 
	/* ========================= */
	
	
	/**
	 * Calculates sample posterior statistics and corrects them with 
	 * unobserved profiles. 
	 * 
	 * @return statistic for M step
	 */
	public PosteriorStatistics Estep()
	{
		long time0 = System.nanoTime(); // TIMING
		PosteriorStatistics S = getSampleStatistics();
		double LLuncorr = S.LL;
		double Fobs = S.profile_count;
		
		// add correction for unobserved
		if (0<min_copies)
		{
			int miss_calc_absolute = Integer.max(5,factory.getCalculationWidthAbsolute());
			double miss_calc_relative = Double.max(1.0, factory.getCalculationWidthRelative());
			
			int num_unobserved_profiles; 
			if (min_copies==1)
			{
				num_unobserved_profiles = 1;
			} else
			{
				assert (min_copies==2);
				int n1 = factory.tree.getNumLeaves(); 
				num_unobserved_profiles = n1+1;
			}
			double[] unobsLL = new double[num_unobserved_profiles];
			int ui=0;
			PosteriorStatistics Sunobs = new PosteriorStatistics();
			StraightLikelihood SL0 = new StraightLikelihood(factory, ProfileTable.emptyProfile(factory.tree));
			SL0.setCalculationWidthThresholds(miss_calc_absolute, miss_calc_relative);
			Profile P0 = SL0.getProfile(0);
			double L0 = 
			unobsLL[ui++] = P0.getLogLikelihood();
			Sunobs.add(P0, L0);
			
			double Lunobs;
			if (min_copies == 1)
			{
				Lunobs = L0;
			} else
			{
				ProfileTable singletons = ProfileTable.singletonTable(factory.tree);
				StraightLikelihood SL1 = new StraightLikelihood(factory, singletons);
				SL1.setCalculationWidthThresholds(factory.getCalculationWidthAbsolute(), factory.getCalculationWidthRelative());
				int ns = singletons.getFamilyCount(); 
				for (int f=0; f<ns; f++)
				{
					Profile P1 = SL1.getProfile(f);
					double L1 = 
					unobsLL[ui++] = P1.getLogLikelihood();
					Sunobs.add(P1, L1);
				}
				Lunobs = Logarithms.sum(unobsLL, ui);
			}			
			
			double pobs = -Math.expm1(Lunobs);
			int F = utable.getTotalFamilyCount();
			double LLcorr = LLuncorr - Fobs * Math.log(pobs);
			
			S.add(Sunobs, F/pobs);
			
//			System.out.println("#**SEM.E LLcorr "+LLcorr+"\tuncorr "+LLuncorr+"\tS.LL "+S.LL);
			S.LL = LLcorr; 
		}
		
		timeE += System.nanoTime()-time0;
		
		return S;
	}	
	
	private static final boolean MSTEP_POISSON_AT_SMALL_DUP = true;
	private static final boolean MSTEP_POISSON_AT_BIG_KAPPA = false;
	private static final boolean MSTEP_ENFORCE_MIN_GAIN = false;
	private static final boolean MSTEP_KEEP_POISSON = true;

	/**
	 * M-step: parameter estimation from posterior statistics. Resets {@link #factory} 
	 * parameters and copies them to the underlying rate model. 
	 * 
	 * @param E
	 */
	public void Mstep(PosteriorStatistics E)
	{
		long time0 = System.nanoTime(); // TIMING
		int num_nodes = factory.tree.getNumNodes();

		final double F = E.profile_count; // already corrected for unobserved profiles
		final double kappa_tol = 1.0/(1L<<30); // tolerance for line minimization
//		int maxiter = 60; 
		final double small_dL =  kappa_tol; //1.0/(1L<<26); //1.0/(1L<<48); //  1.0/(1L<<40); // 0.0; //1e-12;
		final double bracketing_factor = 2.0;
		
		for (int v = num_nodes-1; v>=0; v--) // in preorder
		{
			if (optimize_node[v])
			{
				if (OPT_UNIVERSAL_GAIN_PARAMETER)
					MStepWithUniversalGain(v, E);
				else
				{
					final boolean fixed_loss =  isFixedParameter(v, PARAMETER_LOSS);

					double[] pSv = E.log_edge_posteriors[v].clone();
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
					
					
					// set loss
					double logit_p;
					double logNu_Sv;
					if (factory.tree.isRoot(v) || factory.getLogitLossParameter(v)==Double.POSITIVE_INFINITY)
					{
						logit_p = Double.POSITIVE_INFINITY; // p=1.0
						logNu_Sv = Double.NEGATIVE_INFINITY; // not needed
					} else 
					{
	//					double[] pNu = E.log_node_posteriors[factory.tree.getParent(v)].clone();
	//					{
	//						double log_tail = Double.NEGATIVE_INFINITY;
	//						int j=pNu.length;
	//						while (0<j)
	//						{
	//							--j;
	//							double x = pNu[j];
	//							pNu[j] = log_tail;
	//							log_tail = Logarithms.add(log_tail, x);
	//						}
	//					}
	//					double logNu = Logarithms.sum(pNu, pNu.length);
						
						// set loss parameter
						double[] tNu_Sv = E.log_death_tails[v].clone();
						logNu_Sv = Logarithms.sum(tNu_Sv, tNu_Sv.length);
						logit_p = logNu_Sv-logSv;
						
	//					System.out.println("#**SEM.M "+v
	//							+"\tp "+Math.exp(Logarithms.logitToLogValue(logit_p))
	//							+"/1-"+Math.exp(Logarithms.logitToLogComplement(logit_p))
	//							+"\tSv "+Math.exp(logSv)
	//							+"\tdiff "+Math.exp(logNu_Sv)
	//							+"\tNu "+Math.exp(logNu));
					}
					if (fixed_loss)
					{
						logit_p = factory.getLogitLossParameter(v);
					}
					// set duplication and gain
					final boolean fixed_gain = this.isFixedParameter(v, PARAMETER_GAIN);
					final boolean fixed_dup = this.isFixedParameter(v, PARAMETER_DUPLICATION);
					
					double[] tNv_Sv = E.log_birth_tails[v];
					double logNv_Sv = Logarithms.sum(tNv_Sv, tNv_Sv.length);
	
					double rmax =  Math.exp(logNv_Sv)/F; // with q=0
					
					double logit_q, r;
	
					logit_q = factory.getLogitDuplicationParameter(v);
					r = factory.getUniversalGainParameter(v);

					final double logF = Math.log(F);

					boolean fixed_dup0 = fixed_dup && logit_q==Double.NEGATIVE_INFINITY;
					
					if (fixed_dup0 || (MSTEP_KEEP_POISSON && logit_q==Double.NEGATIVE_INFINITY )) // q==0.0
					{
						// Poisson 
						// logit_q = Double.NEGATIVE_INFINITY;
						
						if (!fixed_gain)
							r = Math.exp(logNv_Sv-logF);
						
	//					System.out.println("#**SEM.M Poisson "+v+"\ty "+logit_q+"\tr "+r+"\t"+Arrays.toString(tNv_Sv));
					} else
					{
						DoubleFunction<Double> optLogitq = κ-> logNv_Sv-Logarithms.add(logSv, logF+Math.log(κ));
						// to be used when p=q
						DoubleFunction<Double> optLogitp = κ-> Logarithms.add(logNu_Sv, logNv_Sv)-Logarithms.add(Math.log(2.0)+logSv, logF+Math.log(κ));
	

						if (fixed_gain)
						{
							if (!fixed_dup)
							{
								// we set q 
								//current rate is fixed 
								double kappa = -r/factory.getLogDuplicationComplement(v);
								assert (kappa != Double.POSITIVE_INFINITY); 
								logit_q = optLogitq.apply(kappa);
								if (is_duprate_bounded && logit_p < logit_q)
								{
									if (fixed_loss)
										logit_q = logit_p;
									else
									{
										logit_p = logit_q = optLogitp.apply(kappa);
									}
								}
								r = -kappa*Logarithms.logitToLogComplement(logit_q);								
							}
						} else // gain is not fixed
						{
							// fixed duplication is handled when best q is set for a given kappa
							
							
							// Polya
							// set q and κ with numerical root finding:
							// find the sign change (positive to negative) 
							// of dLdkappa along 
							// the curve for q(kappa) with dLdq=0.0, or q=p if bounded duplication rate 
							// (we are looking for the unique intersection of two convex decreasing 
							// functions for dLdq=0 and dLdkappa=0)
							
							DoubleFunction<Double> logSumTails
							= new DoubleFunction<>()
							{
								@Override
								public Double apply(double κ)
								{
									double log_sum=Double.NEGATIVE_INFINITY;
									double log_κ=Math.log(κ);
									if (1.0<κ)
									{
										
										for (int i=0; i<tNv_Sv.length; i++)
										{
											// ln (e^d/(κ+i)) = d-ln(κ)-log1p(i/κ)
											double t = tNv_Sv[i]-log_κ-Math.log1p(i/κ);
											log_sum = Logarithms.add(log_sum, t);
										}
									} else
									{
										// take first term
										int i=0; 
										log_sum = tNv_Sv[i]-log_κ;
										i++;
										while (i<tNv_Sv.length)
										{
											// ln(e^d/(κ+i)) = d-ln(i)-log1p(κ/i)
											double t = tNv_Sv[i]-Math.log(i)-Math.log1p(κ/i);
											log_sum = Logarithms.add(log_sum, t);
											i++;
										}
									}
									return log_sum;
								}
							};
							
							final double free_logit_p = logit_p;
							double logNu = Logarithms.add(logSv, logNu_Sv);
							final double log1_p = logSv-logNu;
							final double log1_p0 = factory.getLogLossComplement(v);
							final double log1_q0 = factory.getLogDuplicationComplement(v);
							
							final double logit_p0 = factory.getLogitLossParameter(v);
							final double logit_q0 = factory.getLogitDuplicationParameter(v);
							double kappa0 = -r/log1_q0;							
							
							final int node = v;
							DoubleFunction<Double> rel_dLdκ  
							= new DoubleFunction<>()
							{
								@Override
								public Double apply(double κ)
								{
									double log1_q;
									double logit_q;
									if (fixed_dup)
									{
										log1_q = log1_q0;
										logit_q = logit_q0;
									} else
									{
										double q1term = Logarithms.add(logSv, logF+Math.log(κ));
										double qdenom = Logarithms.add(q1term, logNv_Sv);
										log1_q = q1term-qdenom;
										logit_q = logNv_Sv-q1term;
										double max_log1_q = fixed_loss?log1_p0:log1_p;
										if (is_duprate_bounded && log1_q<max_log1_q)
										{
											if (fixed_loss)
											{
												log1_q = log1_p0;
												logit_q = logit_p0;
											}
											else
											{
												// work with q=p
												qdenom = Logarithms.add(logNu,qdenom);
												q1term = Logarithms.add(logSv,q1term);
												log1_q = q1term-qdenom;
												logit_q = Logarithms.add(logNv_Sv,logNu_Sv)-q1term;
											}
										}
									}
									double log_t = logSumTails.apply(κ);
									// divide by (-F*ln(1-q))
									double log_rel 
										= log_t - logF - Logarithms.logitToLogLogComplement(logit_q); //      Math.log(-log1_q);
									
									double retval = Math.expm1(log_rel);
									if (Double.isInfinite(retval)) // DEBUG
									{
										System.out.println("#**SEM.M.rddk node "+node+"\tkappa "+κ+"\tlogrel "+log_rel+"\tlogt "+log_t+"\tlogf "+logF+"\tllq "+Math.log(-log1_q));
									}
									
									return retval;
								}
							};
							
	//						DoubleFunction<Double> dLdκ  
	//						= new DoubleFunction<>()
	//						{
	//							@Override
	//							public Double apply(double κ)
	//							{
	//								double olq = optLogitq.apply(κ);
	//								double t = Math.exp(logSumTails.apply(κ));
	//								if (is_duprate_bounded && free_logit_p < olq ) // certainly not when p=1
	//								{
	//									// work with q=p instead
	//									// and so optimize p
	//									double olp = optLogitp.apply(κ);
	//									return t+F*Logarithms.logitToLogComplement(olp);
	//								} else
	//								{
	//									return t+F*Logarithms.logitToLogComplement(olq);
	//								}
	//							}
	//						};
		
							
	//						double kappa0 = factory.getGainParameter(v);
	//						double dL0 = rel_dLdκ.apply(kappa0);
	//						double dLmin = dL0;
	//						double kappamin = kappa0;
							//	bracketing
							double kappa1; // left endpoint with non-negative derivative
							
							kappa1 = 1.0/256.0; //   Math.min(1.0/128.0, kappa0);	
							
							double d1 = rel_dLdκ.apply(kappa1);
							double kappamin = kappa1;
							double dLmin = d1;
							
							
							
							double kappa2, d2; // right endpoint with non-positive derivative
							
							int iter=0;
							if (0.0<d1) // increase kappa until derivative is negative at the right bracket 
							{
								
								while (0.0<(d2 = rel_dLdκ.apply(kappa2=bracketing_factor*kappa1)) && 
										!(MSTEP_POISSON_AT_BIG_KAPPA && MAX_GAIN_RATE<=kappa2) 
										&& Math.abs(dLmin)>small_dL)
								{
									kappa1 = kappa2;
									d1 = d2;
									if (60< ++iter)
									{
										// throw new RuntimeException("Bracketing failed "+kappa1+"\td1 "+d1);
										
										System.out.println("#**SEM.M max-bracketing failed at node "+v
													+ "\tkappa "+kappa2+"\tdL "+d2
													+"\tbestk "+kappamin+"\tdlmin "+dLmin
													+"\t//" //+" kappa0 "+kappa0
													+"\ty "+logit_q
													+"\tq "+Math.exp(Logarithms.logitToLogValue(logit_q))+"/1-"+Math.exp(Logarithms.logitToLogComplement(logit_q))
													+"\tx "+logit_p
													);
										break;
									}
									if (Math.abs(d1)<Math.abs(dLmin))
									{
//										if (d1<1e-12)
//										{
//											double kappaj = 16.0*kappa1;
//											for (int j=0; j<60 && 1.0/1024<kappaj; j++)
//											{
//												System.out.println("#**SEM.M.dlk small dL "+kappaj+"/"+j+"\tdL "+rel_dLdκ.apply(kappaj)+"\tlq "+optLogitq.apply(kappaj)+"\td1 "+d1);
//												kappaj = kappaj/2.0; 
//											}
//										}
										
										kappamin = kappa1;
										dLmin = d1;
									}
									if (Double.isInfinite(d2) || Double.isInfinite(d1)) // DEBUG
									{
										System.out.println("#**SEM.M node "+v+"\tk1 "+kappa1+"\td1 "+d1+"\tkappa2 "+kappa2+"\td2 "+d2);
									}
	 							}
								if (0.0<d2 )
								{
									assert (0.0<d1);
									kappa1=kappa2 = kappamin;
									if (MSTEP_POISSON_AT_BIG_KAPPA && MAX_GAIN_RATE<kappamin)
										kappamin = MAX_GAIN_RATE;
//									System.out.println("#**SEM.M node "+v+"\tsetting max kappa "+kappa1+"\tdlmin "+dLmin+"\td1 "+d1+"\td2 "+d2+"\tkappa0 "+kappa0+"\td0 "+rel_dLdκ.apply(kappa0));
								}
							} else if (d1==0.0) // lucky lucky, got the minimum 
							{
								kappa2 = kappa1;
								d2 = d1;
							} else // d1<0.0, so kappa should be smaller 
							{
								kappa2 = kappa1;
								d2 = d1;
								double min_rate = MSTEP_ENFORCE_MIN_GAIN?1.0-MAX_PROB_NOT1:0.0;
								while ((d1 = rel_dLdκ.apply(kappa1=kappa2/bracketing_factor))<0.0 
										&& small_dL<Math.abs(dLmin) 
										&& !(MSTEP_ENFORCE_MIN_GAIN && kappa1<=min_rate))
								{
									kappa2 = kappa1;
									d2 = d1;
									if (Math.abs(d2)<Math.abs(dLmin))
									{
										kappamin = kappa2;
										dLmin = d2;
									}
									
									if (60< ++iter)
									{
										System.out.println("#**SEM.M min-bracketing failed at node "+v
												+"\tkappa "+kappa2+"\tdL "+d2
												+"\tbestk "+kappamin+"\tdlmin "+dLmin
												+"\t//" // kappa0 "+kappa0
												+"\ty "+logit_q
												+"\tq "+Math.exp(Logarithms.logitToLogValue(logit_q))+"/1-"+Math.exp(Logarithms.logitToLogComplement(logit_q))
												+"\tx "+logit_p
												);
										kappa1 = kappa2 = kappamin;
										break;
									}
								}
								if (d1<0.0)
								{
									assert (d2<0.0);
									kappa1=kappa2=Double.max(min_rate,kappamin);
//									System.out.println("#**SEM.M node "+v+"\tsetting min kappa "+kappa1);
								}
							}
							
							if (Math.abs(dLmin)<=small_dL)
							{
//								double q = Math.exp(Logarithms.logitToLogValue(logit_q));
//								System.out.println("#**SEM.M "+v+"\tsmalldL "+dLmin+"\tkappa "+kappamin
//										+"\t//" // kappa0 "+kappa0
//										+"\ty "+logit_q
//										+"\tq "+q+"/1-"+Math.exp(Logarithms.logitToLogComplement(logit_q))
//										+"\tx "+logit_p
//										+"\tr "+(-kappamin*Logarithms.logitToLogComplement(logit_q))+"/"+(kappamin*q)
//										);
								kappa1=kappa2 = kappamin;
							}
							
							double kappa;	
							if (kappa1==kappa2)
							{
								kappa = kappa1;
							} else
							{
								// risky, but we look for the minimum as 
								// the zero of the gradient
								kappa = FunctionMinimization.zbrent(
										rel_dLdκ, kappa1, kappa2, kappa_tol);
								
								
								
								// TODO
								// ? use dbrent but then we need to calculate function value : log-likelihood with given kappa and posteriors
							}
							if (!fixed_dup)
							{
								logit_q = optLogitq.apply(kappa);
								if (is_duprate_bounded && logit_p<logit_q)
								{
									if (fixed_loss)
									{
										logit_q = logit_p;
									} else
									{
										logit_p = logit_q = optLogitp.apply(kappa);
									}
								}
							}

							r = -kappa*Logarithms.logitToLogComplement(logit_q);
							
							if (!fixed_dup)
							{
								double rmin=Math.exp(tNv_Sv[0]-logF);
								if (rmax<=r 
										|| (MSTEP_POISSON_AT_SMALL_DUP && Math.exp(Logarithms.logitToLogComplement(logit_q))==1.0) 
										|| (MSTEP_POISSON_AT_BIG_KAPPA && MAX_GAIN_RATE<=kappa) // disallow very large kappa 
									)
								{
									System.out.println("#**SEM.M "+v+"\tbigkappa/smalldup "+kappa+"\ty "+logit_q+"\tr "+r+"\trmax "+rmax+"\tx "+logit_p+"\t(free "+free_logit_p+")"
											+"\topt_y "+logit_q+"\tq "+Math.exp(Logarithms.logitToLogValue(logit_q))+"/1-"+Math.exp(Logarithms.logitToLogComplement(logit_q)));
									
									// should be Poisson instead
									
									r = rmax;
									logit_p = free_logit_p;
									logit_q = Double.NEGATIVE_INFINITY;
								} else if (r<MIN_GAIN_PARAMETER)
								{
									if (MSTEP_ENFORCE_MIN_GAIN)
									{
										System.out.println("#**SEM.M "+v+"\tsmallr/mingain "+r+"\tkappa "+kappa+"\ty "+logit_q+"\trmax "+rmax+"\trmin "+rmin+"\tx "+logit_p+"\t(free "+free_logit_p+")"
												+"\topt_y "+logit_q+"\tq "+Math.exp(Logarithms.logitToLogValue(logit_q))+"/1-"+Math.exp(Logarithms.logitToLogComplement(logit_q))
												);
										r = MIN_GAIN_PARAMETER;
										logit_p = free_logit_p;
										logit_q = Double.NEGATIVE_INFINITY;
									}
								} else if (r<=rmin)
								{
									System.out.println("#**SEM.M "+v+"\tsmallr/rmin "+r+"\tkappa "+kappa+"\ty "+logit_q+"\trmax "+rmax+"\trmin "+rmin+"\tx "+logit_p+"\t(free "+free_logit_p+")"
											+"\topt_y "+logit_q+"\tq "+Math.exp(Logarithms.logitToLogValue(logit_q))+"/1-"+Math.exp(Logarithms.logitToLogComplement(logit_q)));
//									r = rmin;
//									logit_p = free_logit_p;
//									logit_q = Double.NEGATIVE_INFINITY;
								}
							}
						} // Polya
					} // not fixed gain
					factory.setUniversalNodeParameters(v, logit_p, logit_q, r);
				}
			}			
		}	// for nodes
		factory.copyParametersToModel();
		timeM += System.nanoTime()-time0; // TIMING
	}	
	
	/**
	 * Whether start with an adjustment of the input model in {@link #findLikelihoodInterval(int, int, double, double, int)}
	 */
	private static final boolean LIKELIHOOD_INTERVAL_OPTIMIZE_START = false;

	
	/**
	 * LRT test: set max or min adjustment of a parameter 
	 * for log-likelihood increase bounded by chi-square-test
	 * significance level. 
	 * 
	 * @param node lineage
	 * @param param_type loss,gain,or duplication from {@link GLDParameters}
	 * @param significance_alpha chi-square-test significance level: negative for lower tail; positive for upper tail 
	 * @param delta convergence for adjusted model optimization
	 * @param itmax max iteration for adjusted model optimization
	 */
	private void findLikelihoodInterval(int node, int param_type, double significance_alpha, double delta, int itmax)
	{
		PrintStream out = System.out;
		
		final double negLL0 = optimize(delta, LIKELIHOOD_INTERVAL_OPTIMIZE_START?itmax:0);
		double[] model0 = factory.getParameters();
		
		double logit_p0 = factory.getLogitLossParameter(node);
		double logit_q0 = factory.getLogitDuplicationParameter(node);
		double logit_lm0 = Logarithms.logitParameterRatio(logit_p0, logit_q0);
		
		double r0 = factory.getUniversalGainParameter(node);
		double kappa0 = -r0/Logarithms.logitToLogComplement(logit_q0);
		
		double x0;
		if (PARAMETER_LOSS == param_type)
		{
			x0 = logit_p0;
		} else if (PARAMETER_DUPLICATION == param_type)
		{
			if (is_duprate_bounded)
			{
				x0 = logit_lm0;
			} else
				x0 = logit_q0;
		} else
		{
			assert (PARAMETER_GAIN == param_type);
			if (factory.getLogitDuplicationParameter(node)==Double.NEGATIVE_INFINITY)
			{
				x0 = r0;
			}
			else
			{
				x0 = kappa0; // kappa
			}
		}
		if (PRINT_OPTIMIZATION_MESSAGES)
			out.println("#**SEM.fLI start\t"+node+"/"+param_type
					+"\tx0 "+x0
					+"\tlp "+logit_p0+"("+factory.getLogitLossParameter(node)+")"
					+"\tlq "+logit_q0+"("+factory.getLogitDuplicationParameter(node)+")"  
					+"\tr " +r0+"("+factory.getUniversalGainParameter(node)+")"
					+"\t// "+factory.toString(node));
		
		// diff_cache stores difference in log-likelihood for x settings  
		final Map<Double,Double> diff_cache = new HashMap<>();
		diff_cache.put(x0, 0.0);
		// param_cache stores best model parameters for given x settings
		final Map<Double, double[]> param_cache = new HashMap<>();
		param_cache.put(x0, model0);
		
		DoubleFunction<Double> diffLL  = new DoubleFunction<>()
		{
			@Override
			public Double apply(double x)
			{
				double diffLL;
				if (diff_cache.containsKey(x)) 
				{
					diffLL= diff_cache.get(x);
				} else
				{
					
					double logit_p, logit_q, r;
					if (PARAMETER_LOSS == param_type)
					{
						logit_p = x;
					} else
					{
						logit_p = logit_p0;
					}
					if (PARAMETER_DUPLICATION == param_type)
					{
						if (is_duprate_bounded)
						{
							logit_q = Logarithms.mulLogit(logit_p, x);
						} else
						{
							logit_q = x;
						}
					} else
					{
						if (is_duprate_bounded)
						{
							logit_q = Logarithms.mulLogit(logit_p, logit_lm0);
						} else
						{
							logit_q = logit_q0;
						}
					}
					if (PARAMETER_GAIN == param_type)
					{
						if (logit_q == Double.NEGATIVE_INFINITY)
						{
							r = x;
						} else
						{
							// x is kappa
							r = -x*Logarithms.logitToLogComplement(logit_q);
						}
					} else
					{
						if (logit_q == Double.NEGATIVE_INFINITY)
						{
							r = r0;
						} else
						{
							assert (kappa0 != 0.0);
							r = -kappa0*Logarithms.logitToLogComplement(logit_q);
						}
					}
					factory.setUniversalNodeParameters(node, logit_p, logit_q, r);
					factory.copyParametersToModel();
					factory.computeParameters();

					if (PRINT_OPTIMIZATION_MESSAGES)
						out.println("#**SEM.fLI.dL\tset "+node
								+"\tlp "+logit_p+"("+factory.getLogitLossParameter(node)+")"
								+"\tlq "+logit_q+"("+factory.getLogitDuplicationParameter(node)+")"  
								+"\tr " +r+"("+factory.getUniversalGainParameter(node)+")"
								+"\tx "+x+"\t// now "+factory.toString(node));
					
					
					if (PARAMETER_LOSS == param_type)
					{
						fixLoss(node,true);
					} else if (PARAMETER_DUPLICATION == param_type)
					{
						fixDuplication(node, true);
					} else
					{
						fixGain(node, true);
					}
					if (PRINT_OPTIMIZATION_MESSAGES)
						out.println("#**SEM.fLI.dL\teval "+x+"\tx0 "+x0+"\t// start "+factory.toString(node));
					double[] model_start = factory.getParameters();
					
					PosteriorStatistics S = Estep();
					double debugLL = -S.LL;
					if (PRINT_OPTIMIZATION_MESSAGES)
						out.println("#**SEM.fLI.dL\tdebugE negLL "+debugLL); 

					Mstep(S);
					S = Estep();
					debugLL = -S.LL;
					if (PRINT_OPTIMIZATION_MESSAGES)
						out.println("#**SEM.fLI.dL\tdebugM+E negLL "+debugLL); 					
					
					double negLL = optimize(delta, itmax);
					diffLL = negLL-negLL0;
					double[] model = factory.getParameters();
					diff_cache.put(x, diffLL);
					param_cache.put(x, model);
					
					if (PRINT_OPTIMIZATION_MESSAGES)
						out.println("#**SEM.fLI.dL put\t"+x+"\tdiff "+diffLL+"\tx0 "+x0+"\tLL "+negLL+"\tLL0 "+negLL0
									+"\tlp "+logit_p+"("+factory.getLogitLossParameter(node)+")"
									+"\tlq "+logit_q+"("+factory.getLogitDuplicationParameter(node)+")"  
									+"\tr " +r+"("+factory.getUniversalGainParameter(node)+")"
									+"\t// "+factory.toString(node));
					
					
					if (diffLL<0.0) // we found a better model than the input?? 
						// DEBUG
					{
						double[] model_end = factory.getParameters();
						StringBuilder sb = new StringBuilder();
						for (int j=0; j<model_start.length; j++)
						{
							if (model_start[j]!=model_end[j])
							{
								if (0<sb.length())
									sb.append("\n#**SEM.fLI.dL ; ");
								int node = j/3;
								int type = j%3;
								sb.append("par").append(j).append('(').append(node).append('/').append(type)
								.append("=").append(model_end[j]).append(" was ").append(model_start[j]);
							}
						}
						if (sb.length()==0)
						{
							sb.append("none");
						}
						out.println("#**SEM.fLI.dL decrease "+diffLL+"\tchanges "+sb.toString());
						
						
						
						
						// we might as well quit now
						RateVariationModel vmodel = new RateVariationModel(factory.rates);
						vmodel.initConstantRates();
						out.println(RateVariationParser.printRates(vmodel));
						out.println("#SCORE\t"+(negLL0+diffLL)+"\twas "+negLL0);
						System.exit(99);
					}
					
					// reset
					if (PARAMETER_LOSS == param_type)
					{
						fixLoss(node,false);
					} else if (PARAMETER_DUPLICATION == param_type)
					{
						fixDuplication(node, false);
					} else
					{
						fixGain(node, false);
					}
					
					factory.setParametersAndCopy(model0);
				}
				// LRT test
				double chi_square_p = Functions.Chi_square_tail(1, 2.0*Math.max(0.0, diffLL)); 
				if (PRINT_OPTIMIZATION_MESSAGES)
					out.println("#**SEM.fLI.dL\t"+x+"\tdiff "+diffLL+"\tx0 "+x0+"\tLL0 "+negLL0+"\tp "+chi_square_p);
				return Math.abs(significance_alpha) - chi_square_p;
			}
		};
		
		double step = (PARAMETER_GAIN == param_type?1.25:1.0);
		int bracket_iter = 30;
		double xtol = 1.0/(1<<20); 

		if (0.0<significance_alpha)
		{
			// bracket for xmax
			
			double x2,d2;
			if (x0==Double.NEGATIVE_INFINITY)
			{
				x2 = -Math.log(1L<<52);
				d2 = diffLL.apply(x2);
			} else
			{
				x2 = x0;
				d2 = significance_alpha-1.0;
			}
			assert (d2<0.0);

			double x1,d1; // left bracket 
			int iter=0;
			do
			{
				x1 = x2;
				d1 = d2;
				if (PARAMETER_GAIN == param_type)
					x2 = x1*step;
				else
					x2 = x1+step;
				d2 = diffLL.apply(x2);
				++iter;
			} while (d2<0.0 && iter<bracket_iter);
			
			double xmax;
			if (d2==0.0) // unlikely
			{
				// got rmax 
				xmax = x2;
			} else if (d2<0) // at maxiter 
			{
				xmax = x2;
			} else 
			{
				// d2>0
				if (PRINT_OPTIMIZATION_MESSAGES)
					out.println("#**SEM.fLI brackethi\t"+x1+"\t"+x2);
				xmax = FunctionMinimization.zbrent(diffLL, x1, x2, xtol);
			}
			
			if (PRINT_OPTIMIZATION_MESSAGES)
			{
				Double[] cached_x =  diff_cache.keySet().toArray(new Double[0]);
				Arrays.sort(cached_x);
				for (double x: cached_x)
				{
					out.println("#*SEM.fLI\t"+node+"/"+param_type+"\t"+x
							+"\t"+diff_cache.get(x)
							+"\t"+(Math.abs(significance_alpha)-diffLL.apply(x))
							);
				}
				
			}
			
			
			
			double dmax = diffLL.apply(xmax); // caches model for rmax
			double[] model = param_cache.get(xmax);
			
			factory.setParametersAndCopy(model);
			if (PRINT_OPTIMIZATION_MESSAGES)
				out.println("#**SEM.fLI\txmax "+xmax+"\tmax "+dmax+"\tx0 "+x0+"\t// "+factory.toString(node));
		} else
		{
			// bracket for xmin
			
			double x1,d1;
			if (x0==Double.POSITIVE_INFINITY)
			{
				x1 = Math.log(1L<<52);
				d1 = diffLL.apply(x1); 
			} else
			{
				x1 = x0;
				d1 = -significance_alpha-1.0;
			}
			double x2,d2;
			int iter = 0;
			do
			{
				x2 = x1;
				d2 = d1;
				if (PARAMETER_GAIN == param_type)
					x1 = x2/step;
				else
					x1 = x2-step;
				d1 = diffLL.apply(x1);
				
				if (diff_cache.get(x1)<0.0) break; // what?
				
				++iter;
			} while (d1<0.0 && iter<bracket_iter); 
			
			double xmin;
			if (d1<=0.0) xmin = x1;
			else
			{
				if (PRINT_OPTIMIZATION_MESSAGES)
					out.println("#**SEM.fLI bracketlow\t"+x1+"\t"+x2);
				xmin = FunctionMinimization.zbrent(diffLL, x1, x2, xtol);
			}

			
			if (PRINT_OPTIMIZATION_MESSAGES)
			{
				Double[] cached_x =  diff_cache.keySet().toArray(new Double[0]);
				Arrays.sort(cached_x);
				for (double x: cached_x)
				{
					out.println("#*SEM.fLI\t"+node+"/"+param_type+"\t"+x
							+"\t"+diff_cache.get(x)
							+"\t"+(Math.abs(significance_alpha)-diffLL.apply(x))
							);
				}
				
			}
			
			
			double dmin = diffLL.apply(xmin);
			double[] model = param_cache.get(xmin);
			factory.setParametersAndCopy(model);
			if (PRINT_OPTIMIZATION_MESSAGES)
				out.println("#**SEM.fLI\txmin "+xmin+"\tmin "+dmin+"\tx0 "+x0+"\t// "+factory.toString(node));
		}
		
		
		
		double xbest=x0;
		double dbest=0.0;
		
		for (double x: diff_cache.keySet())
		{
			if (diff_cache.get(x)<dbest)
			{
				dbest = diff_cache.get(x);
				xbest = x;
			}
		}
		
		if (dbest<0.0)
		{
			// a better model was found
			double[] model = param_cache.get(xbest);
			factory.setParametersAndCopy(model);
			out.println("#*SEM.fLI\txbest "+xbest+"\tbest "+dbest+"\tx0 "+x0+"\t// "+factory.toString(node));
			
			
		}
		
	}
	

	/**
	 * Does not converge as well, needs more tuning.
	 * 
	 * @param v
	 * @param E
	 * @deprecated
	 */
	private void MStepWithUniversalGain(int v, PosteriorStatistics E)
	{
		final int maxiter = 256;
		final double rtol = 1.0/(1<<40); // tolerance for zero-finding
		final double gtol = 0.0; //1e-6; // tolerance for gradient
		
		double[] pSv = E.log_edge_posteriors[v].clone();
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
		
		// set loss
		double logit_p;
		double logNu_Sv;
		if (factory.tree.isRoot(v))
		{
			logit_p = Double.POSITIVE_INFINITY; // p=1.0
			logNu_Sv = Double.NEGATIVE_INFINITY; // not needed
		} else
		{
//			double[] pNu = E.log_node_posteriors[factory.tree.getParent(v)].clone();
//			{
//				double log_tail = Double.NEGATIVE_INFINITY;
//				int j=pNu.length;
//				while (0<j)
//				{
//					--j;
//					double x = pNu[j];
//					pNu[j] = log_tail;
//					log_tail = Logarithms.add(log_tail, x);
//				}
//			}
//			double logNu = Logarithms.sum(pNu, pNu.length);
			
			// optimal loss parameter
			double[] tNu_Sv = E.log_death_tails[v].clone();
			logNu_Sv = Logarithms.sum(tNu_Sv, tNu_Sv.length);
			logit_p = logNu_Sv-logSv;
		}
		// save the best logit p
		double free_logit_p = logit_p; 
		
		
		double F = E.profile_count; // already corrected for unobserved profiles
		double logF = Math.log(F);
		double[] tNv_Sv = E.log_birth_tails[v];
		double logNv_Sv = Logarithms.sum(tNv_Sv, tNv_Sv.length);
		
		double rmax = Math.exp(logNv_Sv-logF); // with q=0
		double rmin = Math.exp(tNv_Sv[0]-logF); // with q->1 and -ln(1-q) going to infinitiy
		
//		boolean optimize_dup = !factory.tree.isRoot(v) || factory.getLogitDuplicationParameter(v)!=Double.NEGATIVE_INFINITY;
		
		DoubleFunction<Double> optLogitq = κ-> logNv_Sv-Logarithms.add(logSv, logF+Math.log(κ));
		// to be used when p=q
		DoubleFunction<Double> optLogitp = κ-> Logarithms.add(logNu_Sv, logNv_Sv)-Logarithms.add(Math.log(2.0)+logSv, logF+Math.log(κ));
		
		double logit_q =factory.getLogitDuplicationParameter(v); 

		boolean dup_fixed0 = factory.tree.isRoot(v) && logit_q==Double.NEGATIVE_INFINITY;
		
		
		double r;
		
		if (dup_fixed0)
		{
			r = rmax;
			assert (logit_q == Double.NEGATIVE_INFINITY);
		} else
		{
			// optimization
			int iter=0; 
			
			// initialize r 
			r = factory.getUniversalGainParameter(v);
//			if (logit_q==Double.NEGATIVE_INFINITY)
//			{
//			} else
//			{
//				double κ=factory.getGainParameter(v);
//				r = -κ*Logarithms.logitToLogComplement(logit_q);
//			}
			if (r<rmin || rmax<r) // must stay between brackets
			{
				r = (rmax+rmin)/2.0;
			}

			// initialize logit_q
			// do not keep q=0
			if (logit_q==Double.NEGATIVE_INFINITY || is_duprate_bounded && free_logit_p<logit_q)
			{
				if (logit_p == Double.POSITIVE_INFINITY)
				{
					logit_q = Math.log(TreeWithRates.DEFAULT_DUPLICATION_RATE)-Math.log(1.0-TreeWithRates.DEFAULT_DUPLICATION_RATE);
				} else
				{
					logit_q = logit_p-1.0;
				}
			}
			

			
			while (iter<maxiter)
			{
				// find best r for this q 
				final double log1_q = Logarithms.logitToLogComplement(logit_q);
				DoubleFunction<Double> dLdr = new DoubleFunction<>() // with fixed q
				{
					@Override
					public Double apply(double r)
					{
						// want sum_i T_i/(r-i*ln(1-q))
						// kappa < 1.0: r<-ln(1-q)
						double log_sum;
						double log_r = Math.log(r);
						if (-log1_q<=r) // 
						{
							// 1.0<=kappa
							log_sum=Double.NEGATIVE_INFINITY;
							for (int i=0; i<tNv_Sv.length; i++)
							{
								// ln (e^T/(r +i*lq))=T-ln(r)-log1p(i*lq/r)
								double h = tNv_Sv[i]-log_r-Math.log1p(-i*log1_q/r);
								log_sum = Logarithms.add(log_sum, h);
							}
						} else
							{
								int i=0; 
								log_sum = tNv_Sv[i]-log_r;
								i++;
								while (i<tNv_Sv.length)
								{
									// ln(e^T/(r+i*lq)) = d-ln(i*lq)-log1p(r/(i*lq))
									double t = tNv_Sv[i]-Math.log(-i*log1_q)-Math.log1p(-r/(i*log1_q));
									log_sum = Logarithms.add(log_sum, t);
									i++;
								}
							}
							double sum =Math.exp(log_sum);
							return sum-F;
						}
				};
				double prev_r = r;
				double grad_r = dLdr.apply(r);
//				{ // DEBUG
//					double drmin = dLdr.apply(rmin*0.95);
//					double drmax = dLdr.apply(rmax*1.05);
//					if (drmin<0.0 || drmax>0.0 || iter==0)
//					{
//						String badbracket = (drmin<0.0?".BADmin":"")
//								+(drmax>0.0?".BADmax":"");
//						System.out.println("#**SEM.M "
//								+iter
//								+"\tnode "
//								+v+"\tx "+logit_p+"\ty0 "+logit_q+"\tr "+r+"\tdr "+grad_r+"\trmin "+rmin+"\trmax "+rmax
//								+"\tdrmin "+drmin+"\tdrmax "+drmax
//								+"\tNvSv "+logNv_Sv+"\tT0 "+tNv_Sv[0]+"\tlogF "+logF
//								+"\tSv "+logSv
//								+"\tNu_Sv "+logNu_Sv
//								+"\t// rrng "+(rmax-rmin)+"\trtol "+rtol
//								+"\t"+badbracket);
//					}
//				}
				
				// brackets for zero-finding: make it a little bit bigger 
				// than [rmin,rmax], in order to avoid a sign mistake 
				// with dLdr ~ 0 close to the boundary
				// 
				
				if (gtol < Math.abs(grad_r))
				{
					if (grad_r<0.0)
					{
						r = FunctionMinimization.zbrent(dLdr, rmin*0.95, r, rtol);
					} else
					{
						r = FunctionMinimization.zbrent(dLdr, r, rmax*1.05, rtol);
					}
				}
				double delta_r = r-prev_r;
				
				final double κ = -r/log1_q;

				double prev_y = logit_q;

				logit_q = optLogitq.apply(κ);
				if (is_duprate_bounded && free_logit_p<logit_q)
				{
					logit_p=logit_q = optLogitp.apply(κ);
				} else
				{
					logit_p = free_logit_p;
				}

				double delta_y = logit_q-prev_y;
				
				// this changes r as well
				double set_r = r; // before changing q 
				r = -κ*Logarithms.logitToLogComplement(logit_q);
				delta_r += r-set_r;
				
				double rd_r = Math.abs(delta_r)/Math.max(prev_r,1.0);
				double rd_y = Math.abs(delta_y)/Math.max(Math.abs(prev_y), 1.0);
				
				boolean done_delta_ry = rd_r<rtol && rd_y<rtol;
				
				double dup_ratio = Math.exp(Logarithms.logitToLogValue(logit_q)-Logarithms.logitToLogValue(logit_p));
				
//					if (done_delta_ry || iter==maxiter-1)
//						System.out.println("#**SE.M "+v+"\t"+iter+"\tsety "+logit_q+"\tdy "+delta_y
//								+"\tsetr "+r+"\tdr "+delta_r+"\t(bykappa "+(r-set_r)+")"
//								+"\tdrate "+dup_ratio
//								+"\t"+(done_delta_ry?"DONE(delta)":"loop"));
				
				if (done_delta_ry)  break;
				
				iter++;
			} // while 
		} // unless fixed no-duplication

		// we have r and logitq and logitp
		
		if (logit_q != Double.NEGATIVE_INFINITY)
		{
			// calculate kappa
			double κ = -r/Logarithms.logitToLogComplement(logit_q);
			
			
			if (MAX_GAIN_RATE<=κ || rmax<=r)
			{
				System.out.println("#**SEM.M "+v+"\tbigkappa "+κ+"\ty "+logit_q+"\tr "+r+"\trmax "+rmax+"\tx "+logit_p+"\t(free "+free_logit_p+")"
						+"\topt_y "+logit_q+"\tq "+Math.exp(Logarithms.logitToLogValue(logit_q))+"/1-"+Math.exp(Logarithms.logitToLogComplement(logit_q)));
				
				// should be Poisson instead?
				
				r = rmax;
				logit_p = free_logit_p;
				logit_q = Double.NEGATIVE_INFINITY;
//				factory.setNodeParameters(v, logit_p, logit_q, r);
				
			} 
			factory.setUniversalNodeParameters(v, logit_p, logit_q, r);
		} else
		{
			factory.setUniversalNodeParameters(v, logit_p, logit_q, r);
		}
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
			LL = LL + mul*P.getLogLikelihood();// Math.fma(mul, P.getLogLikelihood(), LL); // fma(a,b,c)=a*b+c
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
	
	@Override
	public double optimize(double eps)
	{
		return this.optimize(eps, Integer.MAX_VALUE);
	}
	
	private long timeM; // TIMING
	private long timeE; // TIMING
	final static double nano = 1e-9; // TIMING

	public double optimize(double eps, int maxiter)
	{
		timeM = timeE = 0L; // TIMING
		
		List<Double> history = getOptimizationHistory();
		if (history == null)
		{
			history = new ArrayList<>();
			setOptimizationHistory(history);
		}
		
		//firePropertyChange​(PROPERTY_OPTIMIZATION_PHASE, "E-step");
		PosteriorStatistics E = Estep();
		double LLstart = E.LL;
		if (history!=null) history.add(-LLstart);
		if (PRINT_OPTIMIZATION_MESSAGES)
		{
			Count.out.println("#*SEM.o start "+LLstart);
		}
		double LLprev = LLstart;
		double[] xprev = factory.getParameters();
		int iter = 0;
		while (iter<maxiter)
		{
			//firePropertyChange​(PROPERTY_OPTIMIZATION_PHASE, "M-step");
			Mstep(E);
			//firePropertyChange​(PROPERTY_OPTIMIZATION_PHASE, "E-step");
			E = Estep();
			
			// check convergence
			double LLnow = E.LL;
			double diff = LLprev-LLnow; // should be negative
			double delta = -diff/LLprev; // should be negative 
			double[] x = factory.getParameters();
			// test convergence on x 
			double max_xd = 0.0;
			for (int i=0; i<x.length; i++)
			{
				double xd;
				if (Double.isInfinite(x[i])||Double.isInfinite(xprev[i]))
				{ // with p==1.0 at the root, or q==0.0 
					if (x[i]==xprev[i])
						xd = 0.0;
					else if (x[i]==-xprev[i])
							xd=2.0;
					else
					{
						assert (Double.isInfinite(x[i]) && Double.isFinite(xprev[i]))
							|| (Double.isFinite(x[i]) && Double.isInfinite(xprev[i]));
						
						xd = x[i]-xprev[i]; // Double.POSITIVE_INFINITY;
						{ // DEBUG
							int node = i/3;
							int partype = i%3;
							
//							System.out.println("#**SEM.o "+iter+"/par"+i+"("+node+"/"+partype+")"+"\twas "+xprev[i]+"\tnow "+x[i]
//									+"\t// "+factory.toString(node));

						}
					}
				} else
				{
					double xdiff =x[i]-xprev[i];  
					xd = xdiff/Double.max(Math.abs(xprev[i]),1.0);
				}
				max_xd = Double.max(max_xd, Math.abs(xd));
//				System.out.println("#**SEM.o xd "+i+"\t"+xd+"\tx "+x[i]+"\txp "+xprev[i]);
			}
			String timing_info = "\ttiming\tavgE "+(nano*timeE/(iter+1.0))+"\ttotE "+(nano*timeE)
					+"\tavgM "+(nano*timeM/(iter+0.0))+"\ttotM "+(nano*timeM);
			if (PRINT_OPTIMIZATION_MESSAGES)
			{
				Count.out.println("#*SEM.o "+iter+"\tLL "+LLnow+"\tincrease "+(-diff)+"\trdiff "+(-delta)
						+"\tmaxxd "+max_xd
						);				
			}
			if (LLnow<LLprev) // dubious M step
			{
				Count.out.println("#*SEM.o done/decrease ("+(-diff)+")"+"\t LL "+LLnow+"\ttotincrease "+(LLprev-LLstart)+timing_info);
				factory.setParametersAndCopy(xprev);
				iter++;
				break;
			}
			LLprev = LLnow;
			history.add(-LLprev);
			xprev = x;
			
			// check convergence
			boolean done_dx = max_xd<=FunctionMinimization.DFP_TOLX;
			boolean done_dL = (-delta<eps);
			
			if ( done_dx || done_dL)
			{
				String reason = (done_dx?"dx.":"") + (done_dL?"dL.":"");
				Count.out.println("#*SEM.o done/converged ("+reason+"@"+iter+")"+"\t LL "+LLprev+"\tdiff "+(LLstart-LLprev)+"\tdrop "+diff+"\trdiff "+delta+timing_info);
				iter++;
				break;
			}
			++iter;
			if (iter==maxiter)
			{
				Count.out.println("#*SEM.o done/iterations ("+iter+")"+"\t LL "+LLprev+"\tdiff "+(LLstart-LLprev)+"\tdrop "+diff+"\trdiff "+(-delta)+timing_info);
			}
		}
		return -LLprev;
	}	
	
	
	
	
	
	
	public static void main(String[] args) throws Exception
	{
		PRINT_OPTIMIZATION_MESSAGES = true;

		Class<?> our_class = java.lang.invoke.MethodHandles.lookup().lookupClass();
		CommandLine cli = new CommandLine(args, our_class);

        PrintStream out = System.out; 
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(our_class));
    	    out.println(CommandLine.getStandardRuntimeInfo(our_class, args));
    	}
        AnnotatedTable table = cli.getTable();
    	MixedRateModel model = null; 
    	TreeWithRates rates;
		count.matek.DiscreteDistribution root_prior = cli.getOptionDistribution(OPT_MODEL_ROOT_PRIOR, null);
    	if (cli.getMixedrateModel()==null)
    	{
    		int rnd_seed = cli.getOptionInt(CommandLine.OPT_RND, 0);
    		Random RND = cli.getOptionRND(out);
    		//Count.out.println("#**SEM.main random starting rates: -"+CommandLine.OPT_RND+" "+rnd_seed);
    		
    		rates = new TreeWithRates(cli.getTree(), RND);
    		if (root_prior != null)
    		{
    			rates.setRootDistribution(root_prior);
    			out.println(CommandLine.getStandardHeader("Root prior preset: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
    		}
    		rates.initNodeParameters(rates.getTree().getRoot());
			out.println(CommandLine.getStandardHeader("(Root prior random: "+rates.getRootDistribution()+")"));
    		
//    		model = GammaInvariant.nullModel(cli.getTree(), RND);
    		model = new RateVariationModel(rates);
    		((RateVariationModel)model).initConstantRates();
    		
    		rates = model.getBaseModel();
    		
    	} else
    	{
    		rates = cli.getRates(); 
    		RateVariationModel constant_rates = new RateVariationModel(rates);
    		constant_rates.initConstantRates();
    		model = constant_rates;
    		rates = model.getBaseModel();

    		Random RND = cli.getOptionRND(out);
    		if (RND!=null)
    			rates.setRandom(RND);
    		if (root_prior != null)
    		{
    			rates.setRootDistribution(root_prior);
    			out.println(CommandLine.getStandardHeader("Root prior set: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
    		}
    	}    	
    	
    	assert (rates == model.getBaseModel());
    	
    	StraightEM O= new StraightEM(rates, table);

        int min_copies = Integer.min(2, table.minCopies());
        min_copies = cli.getOptionInt(OPT_MINCOPY, min_copies);
		O.setMinimumObservedCopies(min_copies);
		out.println(CommandLine.getStandardHeader("Minimum observed copies: -"+OPT_MINCOPY+" "+min_copies));
		
		
		int absolute = 1;
		double relative = 1.0;
        if (cli.getOptionValue(OPT_TRUNCATE)!=null)
        {
        	absolute = cli.getOptionTruncateAbsolute();
        	relative = cli.getOptionTruncateRelative();
        } 
    	O.setCalculationWidth(absolute, relative);
        out.println(CommandLine.getStandardHeader("Truncated computation: -"
        		+OPT_TRUNCATE+" "+absolute+","+relative));
		
        O.is_duprate_bounded =  cli.getOptionBoolean("opt.dupbound", O.is_duprate_bounded);
        out.println(CommandLine.getStandardHeader("Bounded duplication rate: -opt.dupbound "+O.is_duprate_bounded));
        
		int maxiter = cli.getOptionInt(OPT_ROUNDS, Integer.MAX_VALUE);
		double eps = cli.getOptionDouble(OPT_EPS, 1.0/(1L<<30));
        out.println(CommandLine.getStandardHeader("Convergence: -"+OPT_ROUNDS+" "+maxiter+" -"+OPT_EPS+" "+eps));

        {
        	double startscore = -O.Estep().LL;
	        // test: compare to real likelihood
	        //Gradient G = new Gradient(rates,table);
        	
        	
	        VariationGradientFactory G = new VariationGradientFactory((RateVariationModel)model,table);
	        
	        
	        G.setMinimumObservedCopies(min_copies);
	        G.setCalculationWidthThresholds(Integer.MAX_VALUE, Double.POSITIVE_INFINITY);
	        double trueLL = G.getCorrectedLL();
	        double diff = (startscore - (-trueLL));
	        
	        
	        
	        double estimated256 = O.timeE*256*nano; 
	        
	        Count.out.println("#*SEM.main startscore "+startscore+"\ttruescore "+(-trueLL)+"\tdiff "+diff+"\trdiff "+diff/(-trueLL)
	        			+"\testimated: "+ ((int)estimated256)+" seconds/256 iterations");
	        O.timeE = 0L; // reset
        }

        int testpnode = cli.getOptionInt("testp", -1);
        int testrnode = cli.getOptionInt("testr", -1);
        int testqnode = cli.getOptionInt("testq", -1);
        
        double score;
        
        if (testpnode!=-1 || testrnode != -1 ||  testqnode != -1)
        {
        	double pvalue = cli.getOptionDouble(OPT_PVALUE, 0.05);
        	int node, param_type;
        	if (testpnode != -1)
        	{
        		node = testpnode;
        		param_type = PARAMETER_LOSS;
        	} else if (testqnode != -1)
        	{
        		node = testqnode;
        		param_type = PARAMETER_DUPLICATION;
        	} else
        	{
        		node = testrnode;
        		param_type = PARAMETER_GAIN;
        	}
        	O.findLikelihoodInterval(node, param_type, pvalue, eps, maxiter);
        	score = O.optimize(eps,0);
        } else
        {
	        // now we optimize the model parameters
	       score = O.optimize(eps, maxiter);
        }

        // report the optimal values
		out.println(count.io.RateVariationParser.printRates(model));

        {
	        // test: compare to real likelihood
	        //Gradient G = new Gradient(model.getBaseModel(),table);
	        VariationGradientFactory G = new VariationGradientFactory((RateVariationModel)model,table);

	        G.setMinimumObservedCopies(min_copies);
	        G.setCalculationWidthThresholds(Integer.MAX_VALUE, Double.POSITIVE_INFINITY);
	        double trueLL = G.getCorrectedLL();
	        double diff = (score - (-trueLL));
	        
	        Count.out.println("#*SEM.main truescore "+(-trueLL)+"\tdiff "+diff+"\trdiff "+diff/(-trueLL));
        }
        
		
		
		// scoring info
		double ascore = score/table.getFamilyCount();

        double bic_pty = 0.5*O.getModelParameterCount()*Math.log(table.getFamilyCount());
		
		out.println("#SCORE "+score+"\tBICpty "+bic_pty+"\tregularized "+(score+bic_pty));
        out.println("#AVGSCORE "+ascore);
        
		if (out != System.out)
		{
			if (out.checkError())
				throw new java.io.IOException("Write failed.");
			out.close();
		}
        
        
	}
    		
}
