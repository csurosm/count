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
import java.util.List;
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
import count.matek.FunctionMinimization;
import count.matek.Logarithms;
import count.model.DirectLikelihood.Profile;

import static count.io.CommandLine.OPT_EPS;
import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_OUTPUT;
import static count.io.CommandLine.OPT_ROUNDS;
import static count.io.CommandLine.OPT_TRUNCATE;



/**
 * Expectation-maximization algorithm using ancestor copies, 
 * using {@link DirectLikelihood}. 
 */
public class DirectEM extends ML implements GLDParameters, Count.UsesThreadpool
{
	private static boolean PRINT_OPTIMIZATION_MESSAGES = false;
	
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
	
	/**
	 * Thread pool used across different calls to {@link #getSampleStatistics()} and {@link #getTransitionCounts()}
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
		if (thread_pool == null && 1<Count.THREAD_PARALLELISM) // && Count.THREAD_UNIT_TASK<Integer.MAX_VALUE)
		{
//			System.out.println("#**DEM.threadPool init: "+Count.THREAD_PARALLELISM+" threads on "+Thread.currentThread());
			thread_pool = Count.threadPool(); // new ForkJoinPool(Count.THREAD_PARALLELISM);	
		}
		return thread_pool;
	}

	public DirectEM(TreeWithRates rates, ProfileTable table)
	{
		this(rates,table,false);
	}
	protected DirectEM(TreeWithRates rates, ProfileTable table, boolean want_parsimony_fitting)
	{
		if (table instanceof UniqueProfileTable)
		{
			this.utable = (UniqueProfileTable) table; 
		}else
		{
			this.utable = new UniqueProfileTable(table);
		}
		this.min_copies = Integer.min(2,table.minCopies());
		
		this.factory = new DirectLikelihood(rates, utable, want_parsimony_fitting);
		
		this.optimize_node = new boolean[factory.tree.getNumNodes()];
		Arrays.fill(optimize_node, true);
	}
	
	/**
	 * Same as the instantiating table, if it 
	 * was a UniqueProfileTable; or else null. 
	 */
	private UniqueProfileTable utable;
	
	private int min_copies;

	@Override
	public void setMinimumObservedCopies(int min_copies)
	{
		if (min_copies<0 || min_copies>2)
			throw new UnsupportedOperationException("Minimum copies must be 0,1 or 2.");
		this.min_copies = min_copies;
	}
	
	private final DirectLikelihood factory;
	private final boolean[] optimize_node;	
	
	@Override
	public void fixNodeParameters(int node, boolean do_not_optimize)
	{
		this.optimize_node[node] = !do_not_optimize;
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
				np += tree.isRoot(u)?2:3;
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
	
	
	private double[] getDistributionParameters()
	{
		int num_nodes = factory.tree.getNumNodes();
		double[] params = new double[3*num_nodes];
		
		for (int node = 0; node<num_nodes; node++)
		{
			double p = factory.getLossParameter(node);
			double q = factory.getDuplicationParameter(node);
			double r = factory.getGainParameter(node);
			
			params[3*node+PARAMETER_GAIN]=r;
			params[3*node+PARAMETER_LOSS]=p;
			params[3*node+PARAMETER_DUPLICATION]=q;
		}
		return params;
	}
	
	private void setDistributionParameters(double[] params)
	{
		int num_nodes = factory.tree.getNumNodes();
		for (int node=0; node<num_nodes; node++)
		{
			double r = params[3*node+PARAMETER_GAIN];
			double p = params[3*node+PARAMETER_LOSS];
			double q = params[3*node+PARAMETER_DUPLICATION];
			
			factory.rates.setParameters(node, r, p, q);
		}
		factory.computeParameters();
	}
	
	@Override
	public double optimize(double eps)
	{
		return this.optimize(eps, Integer.MAX_VALUE);
	}
	
	private long timeM; // TIMING
	private long timeE; // TIMING
	final double nano = 1e-9; // TIMING

	public double optimize(double eps, int maxiter)
	{
		timeM = timeE = 0L; // TIMING

		List<Double> history = getOptimizationHistory();
		if (history == null)
		{
			history = new ArrayList<>();
			setOptimizationHistory(history);
		}
		
		PosteriorStatistics E = Estep();
		
		double LLstart = E.LL;
		if (history!=null) history.add(-LLstart);
		
		if (PRINT_OPTIMIZATION_MESSAGES)
		{
			Count.out.println("#*DEM.o start "+LLstart);
		}
		
		double LLprev = LLstart;
		double[] xprev = getDistributionParameters();
		int iter = 0;
		
		while(iter<maxiter)
		{
			Mstep(E);
			E = Estep();

			double LLnow = E.LL;
			

			double diff = LLprev-LLnow;
			
			double delta = -diff/LLprev;
			
			double[] x = getDistributionParameters();
			// test convergence on x 
			double max_xd = 0.0;
			for (int i=0; i<x.length; i++)
			{
				double xdiff = Math.abs(x[i]-xprev[i]); // p,q less than 1; r is max O(1) 
				double xd = xdiff/Double.max(x[i],1.0);
				max_xd = Double.max(max_xd, xdiff);
			}
			
			String timing_info = "\ttiming\tavgE "+(nano*timeE/(iter+1.0))+"\ttotE "+(nano*timeE)
					+"\tavgM "+(nano*timeM/(iter+0.0))+"\ttotM "+(nano*timeM);
			
			if (PRINT_OPTIMIZATION_MESSAGES)
			{
				Count.out.println("#*DEM.o "+iter+"\tLL "+LLnow+"\tdiff "+diff+"\trdiff "+delta
						+"\tmaxxd "+max_xd
						);				
			}
			if (LLnow<LLprev) // dubious M step
			{
				Count.out.println("#*DEM.o done/decrease ("+iter+")"+"\t LL "+LLprev+"\tdiff "+(LLstart-LLprev)+timing_info);
				setDistributionParameters(xprev);
				iter++;
				break;
			}
			LLprev = LLnow;
			if (history!=null) history.add(-LLprev);
			xprev = x;
			if (max_xd<=FunctionMinimization.DFP_TOLX || -delta<eps)
			{
				Count.out.println("#*DEM.o done/converged ("+iter+")"+"\t LL "+LLprev+"\tdiff "+(LLstart-LLprev)+"\tdrop "+diff+"\trdiff "+delta+timing_info);
				iter++;
				break;
			}
			
			iter++;
			if (iter==maxiter)
			{
				Count.out.println("#*DEM.o done/iterations ("+iter+")"+"\t LL "+LLprev+"\tdiff "+(LLstart-LLprev)+"\tdrop "+diff+"\trdiff "+(-delta)+timing_info);
			}
			
		}
		
		if (!PRINT_OPTIMIZATION_MESSAGES) // DEBUG
		{
			Gradient G = new Gradient(factory.rates, utable);
			G.setCalculationWidthThresholds(12, 3.0);
			G.setMinimumObservedCopies(min_copies);
			double Lunobs2 = G.getUnobservedLL();
			double LLcorr2 = G.getCorrectedLL();
			double pobs2 = -Math.expm1(Lunobs2);
			
			double punobs2 = Math.exp(Lunobs2);
			double pobs = -Math.expm1(E.LLunobs);
			double punobs = Math.exp(E.LLunobs);
			
			double dLL = E.LL-LLcorr2; 
					
			Count.out.println("#**DEM.o/debugLL LLcorr "+E.LL+"/"+LLcorr2+"\tLunobs "+E.LLunobs+"/"+Lunobs2+"\tpunobs "+punobs+"/"+punobs2+"\tpobs "+pobs+"/"+pobs2
						+"\t// diff "+dLL+"\trdiff "+dLL/LLcorr2);
			
		}
		
		return -LLprev;
		
	}
	
	
	/**
	 * Adds b to a, cell by cell, reusing a and/or b if possible. 
	 * If a is null, then b is returned.   
	 * 
	 * @param a expanded if necessary to match length of b
	 * @param b untouched
	 * @return cell-by-cell sums in array a (or b if a is null), or the resized version  of a
	 */
	private static double[] addCells(double[] a, double[] b, double bmul)
	{
		if (b==null) return a;
		if (a==null)
		{
			a = b.clone();
			for (int i=0; i<b.length; i++)
				b[i]*=bmul;
			return b;
		}
		if (a.length<b.length)
			a = Arrays.copyOf(a, b.length);
		for (int i=0; i<b.length; i++) 
			a[i] = Math.fma(bmul, b[i], a[i]); // == bmul*b[i]+a[i]
		return a;
	}

	
	/**
	 * Adds b to a, cell by cell, reusing a and/or b if possible. 
	 * If a is null, then b is returned.   
	 * 
	 * @param a expanded if necessary to match length of b
	 * @param b untouched
	 * @return cell-by-cell sums in array a (or b if a is null), or the resized version  of a
	 */
	private static double[] addCells(double[] a, double[] b)
	{
		if (a==null)
			return b;
		if (a.length<b.length)
			a = Arrays.copyOf(a, b.length);
		for (int i=0; i<b.length; i++)
			a[i] += b[i];
		return a;
	}	
	
	/**
	 * Array of tail sums: t[i]=sum_{j&gt;i} p[j]   
	 * 
	 * @param p may be null
	 * @return array of tail sums 
	 */
	private static double[] tail(double[] p)
	{
		if (p==null) return new double[1];
		int j = p.length;
		double[] t = new double[j];
		double x = 0.0;
		while (j>0)
		{
			--j;
			t[j] = x;
			x+=p[j];
		}
		return t;
	}
	
	/**
	 * Sum across the array. 
	 * 
	 * @param x
	 * @return sum of x[i]
	 */
	private static double sum(double[] x)
	{
		double sum=0.0;
		if (x!=null)
		{
			int i=x.length;
			while (0<i)
			{
				--i;
				double y = x[i];
				sum += y;
			}
		}
		return sum;
	}
	
	
	
	public class PosteriorStatistics
	{
		private PosteriorStatistics()
		{
			int num_nodes = factory.tree.getNumNodes();
			this.node_posteriors = new double[num_nodes][];
			this.edge_posteriors = new double[num_nodes][];
			this.profile_count = 0.0;
			this.LL = 0.0;
		}
		
		private final double[][] node_posteriors;
		private final double[][] edge_posteriors;
		private double profile_count;
		private double LL;
		private double LLunobs = Double.NEGATIVE_INFINITY;
		
		public void add(Profile P, double multiplier)
		{
			for (int node=0; node<node_posteriors.length; node++)
			{
				double[] pN = P.getNodePosteriors(node);
				double[] pS = P.getEdgePosteriors(node);
				
				node_posteriors[node] = addCells(node_posteriors[node], pN, multiplier);
				edge_posteriors[node] = addCells(edge_posteriors[node], pS, multiplier);
			}
			profile_count += multiplier;
			LL = Math.fma(multiplier, P.getLogLikelihood(), LL);
		}
		
		public void add(PosteriorStatistics stats)
		{
			for (int node=0; node<node_posteriors.length; node++)
			{
				node_posteriors[node] = addCells(node_posteriors[node], stats.node_posteriors[node]);
				edge_posteriors[node] = addCells(edge_posteriors[node], stats.edge_posteriors[node]);
			}
			this.profile_count += stats.profile_count;
			this.LL += stats.LL;
		}
		public void add(PosteriorStatistics stats, double multiplier)
		{
			for (int node=0; node<node_posteriors.length; node++)
			{
				node_posteriors[node] = addCells(node_posteriors[node], stats.node_posteriors[node], multiplier);
				edge_posteriors[node] = addCells(edge_posteriors[node], stats.edge_posteriors[node], multiplier);
			}
			this.profile_count += multiplier * stats.profile_count;
			this.LL += multiplier * stats.LL;
		}
	}
	
	protected void reportTransitionCounts(PrintStream out, String prefix)
	{
		double[][] T = getTransitionCounts();
		double[] pty = fitParsimonyPenalty(T);
		int max_m = 0;
		for (double[] rows: T) max_m = Integer.max(rows.length, max_m);
		out.printf("%s\tfrom",prefix);
		for (int m=0; m<max_m; m++)
			out.printf("\t%d", m);
		out.println();
		for (int n=0; n<T.length; n++)
		{
			int m=0; 
			out.printf("%s\t%d", prefix, n);
			while (m<T[n].length)
			{
				out.printf("\t%g", T[n][m]);
				m++;
			}
			while (m<max_m)
			{
				out.print("\t0");
				m++;
			}
			out.println();
		}
	}
	
	/**
	 * Parsimony penalties for {@link Parsimony#setPenalties(double, double, double)}
	 * calculated by fitting to log-odds scores from {@link #getTransitionCounts()}.
	 * 
	 * @return [gain, loss/death, dup]
	 */
	protected double[] fitParsimonyPenalty()
	{
		double[][] T = getTransitionCounts();
		double[] pty = fitParsimonyPenalty(T);
		return pty;
	}
	
	/**
	 * Parsimony penalties for {@link Parsimony#setPenalties(double, double, double)}
	 * calculated by fitting to log-odds scores from the argument.
	 * 
	 * @param trans n-to-m transition counts (from {@link #getTransitionCounts()})
	 * @return [gain, loss/death, dup]
	 */
	private double[] fitParsimonyPenalty(double[][] trans)
	{
		double[][] w=new double[trans.length][];
		for (int n=0; n<w.length; n++)
		{
			w[n]=new double[trans[n].length];
			for (int m=0; m<w[n].length; m++)
			{
				w[n][m] = Math.log(trans[n][n]/trans[n][m]);
//				System.out.println("#**DEM.fPP n "+n+"\tm "+m+"\tw "+w[n][m]);
			}
		}
		
		double avg_m=0.0;
		double avg_w0m = 0.0;
		double tsum = 0.0;
		for (int m=1; m<w[0].length; m++)
		{
			double t = trans[0][m];
			if (t!=0.0)
			{
				tsum += t;
				avg_m += t*m;
				avg_w0m += t*w[0][m];
			}
		}
//		System.out.println("#**DEM.fPP summ "+avg_m+"\tsumw "+avg_w0m+"\ttsum "+tsum);
		
		avg_m /= tsum;
		avg_w0m /= tsum;
		
		double inc_num = 0.0;
		double inc_denom = 0.0;
		for (int m=1; m<w[0].length; m++)
		{
			double t = trans[0][m];
			if (t!=0.0)
			{
				double dm = (m-avg_m);
				inc_num += t*w[0][m]*dm;
				inc_denom += t*dm*dm;
			}
		}		
		for (int n=1; n<w.length; n++)
		{
			for (int d=1, m=n+d; m<w[n].length; m++, d++)
			{
				double t = trans[n][m];
				if (t!=0.0)
				{
					inc_num += t*w[n][m]*d;
					inc_denom += t*d*d;
				}
			}
		}
		double opt_inc = inc_num/inc_denom;
		double opt_gain = avg_w0m-opt_inc*avg_m;
//		System.out.println("#**DEM.fPP inc_num "+inc_num+"\tinc_den "+inc_denom
//				+"\topt_inc "+opt_inc+"\topt_gain "+opt_gain);
		
		tsum=0.0;
		double avg_n=0.0;
		double avg_wn0 = 0.0;
		for (int n=1; n<w.length; n++)
		{
			double t = trans[n][0];
			if (t!=0.0)
			{
				tsum += t;
				avg_n += t*n;
				avg_wn0 += t*w[n][0];
			}
		}
//		System.out.println("#**DEM.fPP sumn "+avg_n+"\tsumw "+avg_wn0+"\ttsum "+tsum);

		avg_n /= tsum;
		avg_wn0 /= tsum;
		
		double dec_num=0.0;
		double dec_denom=0.0;
		for (int n=1; n<w.length; n++)
		{
			double t = trans[n][0];
			if (t!=0.0)
			{
				double dn = n-avg_n;
				dec_num += t*w[n][0]*dn;
				dec_denom += t*dn*dn;
			}
		}
		for (int n=2; n<w.length; n++)
		{
			for (int d=1, m=n-d; m>0; d++, m--)
			{
				double t = trans[n][m];
				if (t!=0.0)
				{
					dec_num += t*w[n][m]*d;
					dec_denom += t*d*d;
				}
			}
		}

		
		double opt_dec = dec_num/dec_denom;
		double opt_loss = avg_wn0-opt_dec*avg_n;
//		System.out.println("#**DEM.fPP dec_num "+dec_num+"\tdec_den "+dec_denom
//				+"\topt_dec "+opt_dec+"\topt_loss "+opt_loss);
		
		double pty_dup = opt_inc/opt_dec;
		double pty_gain = (opt_gain+opt_inc)/opt_dec;
		double pty_loss = (opt_loss+opt_dec)/opt_dec;
		
		
		Count.out.println("#**DEM.fitP gain "+pty_gain+"\tloss "+pty_loss+"\tdup "+pty_dup);
		double[] pty = new double[3];
		pty[PARAMETER_GAIN] = pty_gain;
		pty[PARAMETER_LOSS] = pty_loss;
		pty[PARAMETER_DUPLICATION] = pty_dup;
		
//		for (int n=0; n<w.length; n++)
//		{
//			for (int m=0; m<w[n].length && m<n; m++)
//			{
//				double c =   (n-m)*opt_dec;
//				double v = (n-m)*1.0;
//				
//				if (m==0)
//				{
//					c += opt_loss;
//					v += pty_loss-1.0;
//				}
//				double e = (w[n][m]-c)/w[n][m];
//				System.out.println("#**DEM.fitP\t"+n+"\t"+m+"\t"+trans[n][m]+"\t"+w[n][m]+"\t"+c+"\t"+e+"\t"+(v*opt_dec)+"\t"+v);
//			}
//			for (int m=n+1; m<w[n].length; m++)
//			{
//				double c = (m-n)*opt_inc;
//				double v = (m-n)*pty_dup;
//				if (n==0)
//				{
//					c += opt_gain;
//					v += pty_gain-pty_dup;
//				}
//
//				double e = (w[n][m]-c)/w[n][m];
//				System.out.println("#**DEM.fitP\t"+n+"\t"+m+"\t"+trans[n][m]+"\t"+w[n][m]
//						+"\t"+c+"\t"+e+"\t"+(v*opt_dec)+"\t"+v);
//			}
//		}
		
		
		return pty;
	}
	

	/**
	 * n to m transitions summed across all edges and families
	 * 
	 * @return
	 */
	public double[][] getTransitionCounts()
	{
		int nF = utable.getFamilyCount();
		final int unit_task;
		ForkJoinPool thread_pool = threadPool();
		if (thread_pool != null)
		{
//			unit_task = THREAD_UNIT_TASK;
			unit_task = Count.unitTask(nF);
		}
		else
		{
			unit_task = nF; // do not fork
		}

		class PartialC extends RecursiveTask<double[][]>
		{
			/**
			 * First family index.
			 */
			private final int minF;
			/**
			 * Last family index, exclusive.
			 */
			private final int maxF;

			PartialC(int min, int max)
			{
				this.minF = min;
				this.maxF = max;
			}
			
			@Override
			protected double[][] compute() 
			{
				try
				{
					double[][] T;
					if (maxF-minF > unit_task)
					{
						int medF = (minF+maxF)/2;
						PartialC left = new PartialC(minF, medF);
						PartialC right = new PartialC(medF, maxF);
						left.fork();
						double[][] Tr = right.compute();
						double[][] Tl = left.join();
						
						T = add(Tl, Tr, 1);
					} else
					{
						T = null;
						for (int f=minF; f<maxF; f++)
						{
							double[][] fT = sumNodeTransitions(f); // without multiplicity included
							int mul = utable.getMultiplicity(f);
							T = add(T, fT, mul);
							
						}
					}
					return T;
				} catch (Throwable t)
				{
					throw t;
				}
			}
			
			
			private double[][] sumNodeTransitions(int family)
			{
				DirectLikelihood.Profile P = factory.getProfile(family);

				int num_nodes = factory.tree.getNumNodes();
				double[][] T = null;
				for (int node=num_nodes-1; node>=0; node--)
				{
					double[][] pT = P.getTransitionPosteriors(node);
					T = add(T, pT, 1);
				}
				return T;
			}
			
			private double[][] add(double[][] Tl, double[][] Tr, int rmul)
			{
				if (Tr==null) return Tl;
				double[][] T;
				if (Tl==null)
					T=new double[Tr.length][];
				else
					T = new double[Integer.max(Tl.length, Tr.length)][];
				for (int n=0; n<T.length; n++)
				{
					if (n<Tr.length)
					{
						if (Tl!=null && n<Tl.length)
						{
							int len = Integer.max(Tr[n].length,Tl[n].length);
							T[n] = new double[len];
							for (int m=0; m<len; m++)
							{
								if (m<Tr[n].length)
								{
									if (m<Tl[n].length)
									{
										T[n][m] = Tl[n][m]+Tr[n][m]*rmul;
									} else
									{
										T[n][m] = Tr[n][m]*rmul;
									}
								} else
								{
									T[n][m] = Tl[n][m];
								}
							}
						} else
						{
							if (rmul==1)
								T[n]=Tr[n];
							else
							{
								T[n] = new double[Tr[n].length];
								for (int m=0; m<T[n].length; m++)
									T[n][m] = rmul*Tr[n][m];
							}
						}
					} else
					{
						T[n] = Tl[n];
					}
				}
				return T;
			}
		
		}
		
		PartialC bigjob = new PartialC(0, nF);
		try
		{
			double[][] T;
			if (nF>unit_task)
			{
				T = thread_pool.invoke(bigjob);
			} else
			{
				T = bigjob.compute();
			}
			return T;
		} catch (Throwable t)
		{
			throw new RuntimeException(t);
		}
		
	
	}
	
	
	public void Mstep(PosteriorStatistics stats)
	{
		long time0 = System.nanoTime(); // TIMING
		int num_nodes = factory.tree.getNumNodes();
		double F = stats.profile_count; // already corrected for unobserved profiles
		double kappa_tol = 1.0/(1<<24);
		double[] Nmeans = new double[num_nodes]; // filled in during preorder traversal
		double[][] Ntails = new double[num_nodes][]; // filled in during traversal
		
		final double minp = 1.0-MAX_PROB_NOT1;
		final double minq = minp*TreeWithRates.DEFAULT_DUPLICATION_RATE;
		
		for (int v = num_nodes-1; v>=0; v--) // in preorder, so that parent's posteriors are already stored  
		{
			double[] pNv = stats.node_posteriors[v];
			double[] tNv = tail(pNv);
			double Nv = sum(tNv); // mean
			Nmeans[v] = Nv;
			Ntails[v] = tNv;
			
			if (optimize_node[v])
			{
				double[] pSv = stats.edge_posteriors[v];
				double[] tSv = tail(pSv);
				double Sv = sum(tSv); // mean
				//assert (Sv <= Nv);
				
				double p = factory.getLossParameter(v);
				double p_1 = factory.getLossParameterComplement(v);
				if (factory.tree.isRoot(v) || p_1==0.0)
				{
					p = 1.0;
					p_1 = 0.0;
				} else   // set loss parameter
				{				
					int u = factory.tree.getParent(v);
					double Nu = Nmeans[u]; 
					double[] tNu = Ntails[u];
					double Nu_Sv = 0.0;
					for (int j=0; j<tNu.length && j<tSv.length; j++)
					{
						Nu_Sv += Double.max(0.0, tNu[j]-tSv[j]);
					}
					
					double opt_p = Nu_Sv/Nu;
					
					
					if (opt_p < minp)
					{
//						System.out.println("#**DEM.M node "+v+"\tsmallp "+opt_p+"\tNu-Sv "+Nu_Sv+"\tNu "+Nu+"\tu "+u
//									+"\t// tNu "+Arrays.toString(tNu)+"\tSv "+Arrays.toString(tSv)
//									+"\t// rates "+factory.rates.toString(v)); // DEBUG
						p = minp;
						p_1 = 1.0-p;
					} else
					{
						p = opt_p;
						p_1 = Sv/Nu;
					}
				}
				factory.setLossParameter(v, p);
				factory.setLossComplement(v, p_1);
					
//					if (p==0.0)
//						System.out.println("#**DEM.M "+v+"\tp "+p+"/"+p_1+"\tNu "+Nu+"\tSv "+Sv+"\tsubtract "+(Nu-Sv)+"\tdiff "+Nu_Sv+"\t// "+factory.rates.toString(v));
				
				double[] t = tNv.clone(); // tail differences
				for (int i=0; i<tSv.length; i++)
				{
					t[i] = t[i]-tSv[i];
					t[i] = Double.max(t[i],0.0);
				}
				double Nv_Sv = sum(t);  // mean Nv - mean Sv

				if (factory.getDuplicationParameter(v)==0.0)
				{
					// Poisson
					double r = (Nv_Sv)/F;
					factory.setGainParameter(v, r);
				} else
				{
					// set q and κ with numerical root finding
					
					DoubleFunction<Double> q_1  = κ->(Sv+κ*F)/(Nv+κ*F);
					DoubleFunction<Double> dL  
						= new DoubleFunction<>()
						{
							@Override
							public Double apply(double κ)
							{
								double sum_tails = 0.0;
								for (int i=0; i<t.length; i++)
								{
									sum_tails += t[i]/(κ+i);
								}
								double qterm = F*Math.log(q_1.apply(κ));
								double dLdκ = sum_tails+qterm;
								return dLdκ;
							}
						};
						
					//	bracketing
					double kappa1; // left endpoint with non-positive derivative
					kappa1 = 1.0/128.0;				
					double d1 = dL.apply(kappa1);
					double kappa2, d2; // right endpoint with non-negative derivative
					
					int iter=0;
					if (0.0<d1) // increase kappa until derivative is negative at the right bracket 
					{
						while (0.0<(d2 = dL.apply(kappa2=2.0*kappa1)) && kappa2<MAX_RATE)
						{
							kappa1 = kappa2;
							d1 = d2;
							if (60< ++iter) throw new RuntimeException("Bracketing failed");
						}
						if (0.0<d2)
						{
							assert (0.0<d1);
							kappa1=kappa2 = MAX_RATE;
	//						System.out.println("#**DEM.M node "+v+"\tsetting max kappa "+kappa1);
						}
							
					} else if (d1==0.0) // boom, got the minimum 
					{
						// unlikely
						kappa2 = kappa1;
						d2 = d1;
					} else // d1<0.0, so this is right bracket
					{
						kappa2 = kappa1;
						d2 = d1;
						double min_rate = 1.0-MAX_PROB_NOT1;
						while ((d1 = dL.apply(kappa1=0.5*kappa2))<0.0 && min_rate<kappa1)
						{
							kappa2 = kappa1;
							d2 = d1;
							if (60< ++iter) throw new RuntimeException("Bracketing failed");
						}
						if (d1<0.0)
						{
							assert (d2<0.0);
							kappa1=kappa2=min_rate;
	//						System.out.println("#**DEM.M node "+v+"\tsetting min kappa "+kappa1);
						}
					}
					double κ;
						
					if (kappa1==kappa2)
					{
						κ = kappa1;
					} else
					{
						κ = FunctionMinimization.zbrent(dL, kappa1, kappa2, kappa_tol);
						
						if (κ<1.0-MAX_PROB_NOT1) // DEBUG
						{
							Count.out.println("#**DEM.M node "+v+"\tsmallkappa "+κ+"\tbrkt "+kappa1+","+kappa2);
						}
					}
					
					double opt_q = (Nv_Sv)/(Nv+κ*F);
					double q;
					double one_minus_q;
					
					
					if (opt_q<minq)
					{
//						System.out.println("#**DEM.M node "+v+"\tsmallq "+opt_q+"\tNv-Sv "+Nv_Sv+"\tNv "+Nv
//								+"\tkappa "+κ+"\tF "+F
//								+"\t// tNv "+Arrays.toString(tNv)+"\tSv "+Arrays.toString(tSv)
//								+"\tt "+Arrays.toString(t)
//								+"\t// rates "+factory.rates.toString(v)); // DEBUG
						q = minq;
						one_minus_q = 1.0-q;
					} else if (opt_q>MAX_PROB_NOT1)
					{
						Count.out.println("#**DEM.M node "+v+"\tsmallq1 "+q_1.apply(κ)
								+"\tNv-Sv "+Nv_Sv+"\tNv "+Nv
								+"\tkappa "+κ+"\tF "+F
								+"\t// tNv "+Arrays.toString(tNv)+"\tSv "+Arrays.toString(tSv)
								+"\tt "+Arrays.toString(t)); // DEBUG
						q = MAX_PROB_NOT1;
						one_minus_q = 1.0-q;
					} else
					{
						q = opt_q;
						one_minus_q = q_1.apply(κ);
					}
					
					
					factory.setDuplicationParameter(v, q);
					factory.setDuplicationComplement(v, one_minus_q);
					factory.setGainParameter(v, κ); // recomputes the associated factorials
				}
				factory.setRatesForParameters(v);
			}
		}
		timeM += System.nanoTime()-time0; // TIMING
	}
	
	
	
	/**
	 * Calculates sample posterior statistics and corrects them with 
	 * unobserved profiles. 
	 * 
	 * @return
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
			DirectLikelihood D0 = new DirectLikelihood(factory.rates, ProfileTable.emptyProfile(factory.tree));
			D0.setCalculationWidthThresholds(miss_calc_absolute, miss_calc_relative);
//			D0.setSameCalculationWidthThresholds(factory);
			DirectLikelihood.Profile P0 = D0.getProfile(0);
			double L0 = 
			unobsLL[ui++] = P0.getLogLikelihood();
			double p0 = Math.exp(L0);
			Sunobs.add(P0, p0);
			double Lunobs;
			if (min_copies == 1)
			{
				Lunobs = L0;
			} else
			{
				ProfileTable singletons = ProfileTable.singletonTable(factory.tree);
				DirectLikelihood D1 = new DirectLikelihood(factory.rates, singletons);
				D1.setCalculationWidthThresholds(miss_calc_absolute, miss_calc_relative);
//				D1.setSameCalculationWidthThresholds(factory);
				DirectLikelihood.Profile[] P1 = new DirectLikelihood.Profile[singletons.getFamilyCount()];
				for (int i=0; i<P1.length; i++, ui++)
				{
					P1[i] = D1.getProfile(i);
					double L1 = unobsLL[ui] = P1[i].getLogLikelihood();
					double p1 = Math.exp(L1);
					Sunobs.add(P1[i],p1);
				}
				Lunobs = Logarithms.sum(unobsLL, ui);
			}
			
			double pobs = -Math.expm1(Lunobs);
			int F = utable.getTotalFamilyCount();
			double LLcorr = LLuncorr - Fobs * Math.log(pobs);
			double punobs = Math.exp(Lunobs);
			S.add(Sunobs, F/pobs);
			double Fcorr = S.profile_count;
			if (PRINT_OPTIMIZATION_MESSAGES)
			{
				Gradient G = new Gradient(factory.rates, utable);
				G.setCalculationWidthThresholds(12, 3.0);
				G.setMinimumObservedCopies(min_copies);
				double Lunobs2 = G.getUnobservedLL();
				double pobs2 = -Math.expm1(Lunobs2);
				double LLcorr2 = LLuncorr - Fobs * Math.log(pobs2);
				double punobs2 = Math.exp(Lunobs2);
				Count.out.println("#**DEM.E LLcorr "+LLcorr+"/"+LLcorr2+"\tLLun "+LLuncorr+"\tLunobs "+Lunobs+"/"+Lunobs2+"\tpunobs "+punobs+"/"+punobs2+"\tpobs "+pobs+"/"+pobs2+"\tF* "+Fcorr);
				
			}
//			if (PRINT_OPTIMIZATION_MESSAGES)
//				System.out.println("#**DEM.E LLcorr "+LLcorr+"\tLLun "+LLuncorr+"\tLunobs "+Lunobs+"\tpunobs "+punobs+"\tpobs "+pobs+"\tF* "+Fcorr);
			
			S.LL = LLcorr;
			S.LLunobs = Lunobs;
		}
		timeE += System.nanoTime()-time0;
		
		return S;
	}
	
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
//		System.out.println("#**DEM.gSS start pool "+thread_pool+"\t"+Thread.currentThread().toString());

//		if (thread_pool != null)
//		{
////			unit_task = THREAD_UNIT_TASK;
//			unit_task ;
//		}
//		else
//		{
//			unit_task = nF; // do not fork
//		}
		
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
//							System.out.println("#**DEM.gSS "+f+"\t"+Thread.currentThread().toString());
							DirectLikelihood.Profile P = factory.getProfile(f);
							S.add(P, utable.getMultiplicity(f));
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
			
//			System.out.println("#**DEM.gSS launching 0.."+(nF-1)+"\tunit "+unit_task);
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
	
	
//	private void debugLL(PosteriorStatistics E, String reason)
//	{
//		Gradient G = new Gradient(factory.rates, utable.mappedToTree(factory.tree));
//		G.setMinimumObservedCopies(min_copies);
//		G.setCalculationWidthThresholds(12, 3.0);
//		double gscore = -G.getCorrectedLL();
//		double score = -E.LL;
//		double gLL = G.getLL();
//		double g0 = G.getUnobservedLL();
//		double L0 = E.LLunobs;
//		double dg = (score-gscore);
//		System.out.println("#**DEM.debugLL "+reason+"\tE "+score+"\tG "+gscore+"\tuncorr "+gLL+"\tgunobs "+g0
//						+"\tL0 "+L0
//						+"\t// diff "+dg+"\trdiff "+dg/gscore);
//	}
	
	public static void main(String[] args) throws Exception
	{
		PRINT_OPTIMIZATION_MESSAGES = true;

		CommandLine cli = new CommandLine(args, DirectEM.class);

        PrintStream out = System.out; 
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(DirectEM.class));
    	    out.println(CommandLine.getStandardRuntimeInfo(DirectEM.class, args));
    	}
        AnnotatedTable table = cli.getTable();
    	GammaInvariant model = null; 
    	TreeWithRates rates;
    	if (cli.getModel()==null)
    	{
    		Random RND = cli.getOptionRND(out);
    		System.out.println("#**DEM.main random starting rates: "+RND);
    		model = GammaInvariant.nullModel(cli.getTree(), RND);
    		rates = model.getBaseModel();
    	} else
    	{
    		rates = cli.getRates(); 
    		model = cli.getModel();
    	}
        
        DirectEM D = new DirectEM(rates, table);
        int min_copies = Integer.min(2, table.minCopies());
        min_copies = cli.getOptionInt(OPT_MINCOPY, min_copies);
		D.setMinimumObservedCopies(min_copies);
		
		int absolute = 1;
		double relative = 1.0;
        if (cli.getOptionValue(OPT_TRUNCATE)!=null)
        {
        	absolute = cli.getOptionTruncateAbsolute();
        	relative = cli.getOptionTruncateRelative();
        } 
    	D.setCalculationWidth(absolute, relative);
        out.println(CommandLine.getStandardHeader("Truncated computation: -"
        		+OPT_TRUNCATE+" "+absolute+","+relative));
		

		int maxiter = cli.getOptionInt(OPT_ROUNDS, Integer.MAX_VALUE);
		double eps = cli.getOptionDouble(OPT_EPS, 1.0/(1<<20));
        out.println(CommandLine.getStandardHeader("Convergence: -"+OPT_ROUNDS+" "+maxiter+" -"+OPT_EPS+" "+eps));

        // now we optimize the model parameters
        double score = D.optimize(eps, maxiter);
        
        // test: compare to real likelihood
        Gradient G = new Gradient(rates,table);
        G.setMinimumObservedCopies(min_copies);
        G.setCalculationWidthThresholds(12, 3.0);
        double trueLL = G.getCorrectedLL();
        double diff = (score - (-trueLL));
        System.out.println("#*DEM.main truescore "+(-trueLL)+"\tdiff "+diff+"\trdiff "+diff/(-trueLL));
        // report the optimal values
		out.println(count.io.RateVariationParser.printRates(model));
		
		// scoring info
		double ascore = score/table.getFamilyCount();

        double bic_pty = 0.5*D.getModelParameterCount()*Math.log(table.getFamilyCount());
		
//		out.println("#SCORE "+score);
        
        
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
