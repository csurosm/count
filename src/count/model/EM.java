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


import static count.io.CommandLine.OPT_EPS;
import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_MODEL_ROOT_PRIOR;
import static count.io.CommandLine.OPT_OUTPUT;
import static count.io.CommandLine.OPT_ROUNDS;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;
import java.util.function.DoubleFunction;

import count.Count;
import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.io.CommandLine;
//import count.ds.UniqueProfileTable;
//import count.model.Gradient.Profile;
import count.matek.FunctionMinimization;

import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LOSS;

/**
 * Experimental code for expectation-maximization with ancestral (surviving) copies. 
 * Does not work very robustly. 
 *
 * @deprecated
 */
public class EM extends ML
{
	
	private static boolean PRINT_OPTIMIZATION_MESSAGES= false;
	public EM(TreeWithRates rates, ProfileTable table)
	{
		this(new Gradient(rates, table));
	}
	protected EM(Gradient G)
	{
		this.gradient = G;
		this.optimize_node = new boolean[G.factory.tree.getNumNodes()];
		Arrays.fill(optimize_node, true);
	}
	private Gradient gradient;

	
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
	public void setCalculationWidth(int absolute, double relative)
	{
		gradient.setCalculationWidthThresholds(absolute, relative);
	}
	
	public void setMinimumObservedCopies(int m)
	{
		gradient.setMinimumObservedCopies(m);
	}

	private final boolean[] optimize_node;
	
	
	@Override
	public void fixNodeParameters(int node, boolean do_not_optimize)
	{
		this.optimize_node[node] = !do_not_optimize;
	}
	
	@Override
	public int getModelParameterCount()
	{
		IndexedTree tree = gradient.factory.tree;
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
	/**
	 * Adds b to a, cell by cell, reusing a and/or b if possible. 
	 * If a is null, then b is returned.   
	 * 
	 * @param a expanded if necessary to match length of b
	 * @param b untouched
	 * @return cell-by-cell sums in array a (or b if a is null), or the resized version  of a
	 */
	private static double[] add(double[] a, double[] b)
	{
		if (a==null)
			return b;
		if (a.length<b.length)
			a = Arrays.copyOf(a, b.length);
		for (int i=0; i<b.length; i++)
			a[i] += b[i];
		return a;
	}
	
	private static double mean(double[] p)
	{
		double m = 0.0;
		if (p==null) return m;
		int j = p.length;
		while (1<j)
		{
			--j;
			m = Math.fma(j, p[j], m);
		}
		return m;
	}
	
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
	
	private static double sum (double[] x)
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
	
	/**
	 * Adds b to a, cell by cell, reusing a and/or b if possible. 
	 * If a is null, then b is returned.   
	 * 
	 * @param a expanded if necessary to match length of b
	 * @param b untouched
	 * @return cell-by-cell sums in array a (or b if a is null), or the resized version  of a
	 */
	private static double[] add(double[] a, double[] b, double bmul)
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

	
	// sum posteriors: 
	// N[][...m]  = summing node_posteriors[]
	// S[][...]   = summing edge_posteriors[]
	// Nu-Su at every node (dup) 
	// Nu-Sv on every edge uv (loss)
	private class NodeStatistics
	{
		NodeStatistics()
		{
//			Nu = new double[1];
//			Su = new double[1];
//			Nu_Su = new double[1];
//			Nu_Sv = new double[1];
		}
		private double[] Nv;
		private double[] Sv;
//		private double[] Nv_Sv;
//		private double[] Nu_Sv;
		
		void addStats(NodeStatistics PS, double multiplier)
		{
			this.Nv = add(this.Nv, PS.Nv, multiplier);
			this.Sv = add(this.Sv, PS.Sv, multiplier);
//			this.Nv_Sv = add(this.Nv_Sv, PS.Nv_Sv, multiplier);
//			this.Nu_Sv = add(this.Nu_Sv, PS.Nu_Sv, multiplier); // ==null at root OK
			
		}
		void addStats(NodeStatistics PS)
		{
			this.Nv = add(this.Nv, PS.Nv);
			this.Sv = add(this.Sv, PS.Sv);
//			this.Nv_Sv = add(this.Nv_Sv, PS.Nv_Sv);
//			this.Nu_Sv = add(this.Nu_Sv, PS.Nu_Sv); // ==null at root OK
			
		}
		@Override
		public String toString()
		{
			StringBuilder sb=new StringBuilder();
			double Nmean = mean(Nv);
			double Smean = mean(Sv);
			sb.append("Nm ").append(Nmean);
			sb.append(", Sm ").append(Smean);
			sb.append("\t// N ").append(Arrays.toString(Nv));
			sb.append("\tS ").append(Arrays.toString(Sv));
			return sb.toString();
		}
	}
	
	private class EStatistics
	{
		EStatistics()
		{
			int num_nodes = gradient.factory.tree.getNumNodes();
	 		node_stats = new NodeStatistics[num_nodes];
	 		for (int node=0; node<num_nodes; node++)
	 			node_stats[node] = new NodeStatistics();
		}
		
		final NodeStatistics[] node_stats;
		
		void addProfile(Posteriors.Profile post, double multiplier)
		{
			IndexedTree tree = gradient.factory.tree;
			int num_nodes = tree.getNumNodes();
			NodeStatistics[] profile_stats = new NodeStatistics[num_nodes];
			int node=num_nodes;
			while (0<node) // in postorder
			{
				--node;
				NodeStatistics stats = profile_stats[node] = new NodeStatistics();
				stats.Nv = post.getNodePosteriors(node);
				stats.Sv = post.getEdgePosteriors(node);
//				stats.Nv_Sv = stats.Nv.clone();
//				// assert (stats.Nv.length<=stats.Sv.length);
//				
//				for (int i=0; i<stats.Sv.length && i<stats.Nv.length; i++)
//				{
//					stats.Nv_Sv[i] -= stats.Sv[i];
////					if (stats.Nv[i]<stats.Sv[i])
////					{
////						System.out.println("#**EM.ES.aP node "+node+"\ti "+i+"\tNvi "+stats.Nv[i]+"\tSvi "+stats.Sv[i]);
////					}
//				}
//				if (!tree.isRoot(node))
//				{
//					int parent = tree.getParent(node);
//					stats.Nu_Sv = profile_stats[parent].Nv.clone();
//					for (int i=0; i<stats.Sv.length; i++)
//					{
//						stats.Nu_Sv[i] -= stats.Sv[i];
//					}
//				}
				
			}
			// now add to node_stats
			while (node<num_nodes)
			{
				NodeStatistics NS = node_stats[node];
				NodeStatistics PS = profile_stats[node];
				NS.addStats(PS, multiplier);
//				NS.Nv = add(NS.Nv, PS.Nv, multiplier);
//				NS.Sv = add(NS.Nv, PS.Nv, multiplier);
//				NS.Nv_Sv = add(NS.Nv_Sv, PS.Nv_Sv, multiplier);
//				NS.Nu_Sv = add(NS.Nu_Sv, PS.Nu_Sv, multiplier); // ==null at root
				
				node++;
			}
			
			if (false) { // DEBUG
				int[] copy_numbers = post.inside.get();
				String cns = AnnotatedTable.getPatternString(copy_numbers);
				while (0<node)
				{
					--node;
					if (!tree.isLeaf(node))
						System.out.println("#**EM.ES.aP "+cns+"/"+node+"\t"+multiplier+"\t"+profile_stats[node]
								);
				}
			}
			
			
			
		}
		
		void addStats(EStatistics E)
		{
			int num_nodes = gradient.factory.tree.getNumNodes();
			for (int node=0; node<num_nodes; node++)
			{
				NodeStatistics NS = node_stats[node];
				NodeStatistics PS = E.node_stats[node];
				NS.addStats(PS);
			}
		}
		void addStats(EStatistics E, double multiplier)
		{
			int num_nodes = gradient.factory.tree.getNumNodes();
			for (int node=0; node<num_nodes; node++)
			{
				NodeStatistics NS = node_stats[node];
				NodeStatistics PS = E.node_stats[node];
				NS.addStats(PS, multiplier);
			}
		}
		void addEmpty(double multiplier)
		{
			int num_nodes = gradient.factory.tree.getNumNodes();
			double[] sure0 = new double[1];
			sure0[0]=multiplier;
			for (int node=0; node<num_nodes; node++)
			{
				NodeStatistics NS = node_stats[node];
				NS.Nv = sure0.clone();
				NS.Sv = sure0.clone();
//				NS.Nv_Sv = new double[1];
//				NS.Nu_Sv = new double[1];
			}
		}
		
		void reportStats()
		{
			int num_nodes = gradient.factory.tree.getNumNodes();
			for (int node=0; node<num_nodes; node++)
			{
				NodeStatistics NS = node_stats[node];
				System.out.println("#*EM.ES.rS "+node+"\t"+NS.toString());
			}
		}
				
	}
	
	
	private EStatistics EStep()
	{
		EStatistics oES = getSampleStats();
		//oES.reportStats();
		// now we add unobserved profiles 
		// first (M-)estimate 
		// unobserved profile counts 
		// as hidden variables
		
		int min_copies = gradient.getMinimumObservedCopies();
		int num_unobserved_profiles = 0; 
		if (min_copies == 0)
		{
			num_unobserved_profiles=0;
		} else if (min_copies==1)
		{
			num_unobserved_profiles = 1;
		} else
		{
			assert (min_copies==2);
			int n1 = gradient.singletons().getFamilyCount();
			num_unobserved_profiles = n1+1;
		}
		EStatistics[] uES = new EStatistics[num_unobserved_profiles];
		double Lunobs = gradient.getUnobservedLL();
		double denom =  -Math.expm1(Lunobs); // ==(1-Pr(unobs))
		if (min_copies > 0)
		{
			EStatistics E0 = uES[0] = new EStatistics();
			double L0 = gradient.factory.getEmptyLL();
			double m0 = Math.exp(L0)/denom;
			E0.addEmpty(m0);
			if (min_copies==2)
			{
				ProfileTable singletons = gradient.singletons();
				Posteriors spost = new Posteriors(gradient.factory.rates, singletons);
				for (int i=1; i<uES.length; i++)
				{
					int si = i-1;
					Posteriors.Profile P = spost.getPosteriors(si);
					Likelihood.Profile LP = spost.factory.getProfileLikelihood(si);
					double Ls = LP.getLogLikelihood();
					double m1 = Math.exp(Ls)/denom;
					EStatistics E1 = uES[i] = new EStatistics();
					E1.addProfile(P, m1);
				}
			}
		}
		int nF = gradient.getTotalFamilyCount();
		EStatistics corrES = oES;
		for (int u=0; u<uES.length; u++)
		{
			corrES.addStats(uES[u], nF);
		}
		return corrES;
	}
	private int calls_kappa = 0;
	
	@Override
	public double optimize(double eps)
	{
		return this.optimize(eps, Integer.MAX_VALUE);
	}
	
	private void setSurvivalParameters(double[] survival_params)
	{
		double[] distribution_params = invertSurvivalParameters(survival_params);
		int u=gradient.factory.tree.getNumNodes();
		while (0<u)
		{
			--u;
			double r = distribution_params[3*u+PARAMETER_GAIN]	;
			double p = distribution_params[3*u+PARAMETER_LOSS]	;
			double q = distribution_params[3*u+PARAMETER_DUPLICATION]	;
			
//			q=factory.rates.getDuplicationParameter(u);
//			r=factory.rates.getGainParameter(u);
			
			gradient.factory.rates.setParameters(u, r, p, q);
		}
		gradient.computeParameters();		
	}
	
	private double[] getSurvivalParameters()
	{
		int num_nodes = gradient.factory.tree.getNumNodes();
		//TreeWithRates rates = gradient.factory.rates;
		Likelihood L = gradient.factory;
		double[] parameters = new double[3*num_nodes];
		for (int node=0; node<num_nodes; node++)
		{
			double p = L.getLossParameter(node);
			double q = L.getDuplicationParameter(node);
			double r = L.getGainParameter(node);
			parameters[3*node+PARAMETER_GAIN]=r;
			parameters[3*node+PARAMETER_LOSS]=p;
			parameters[3*node+PARAMETER_DUPLICATION]=q;
		}
		return parameters;
	}
	
	public double optimize(double eps, int maxiter)
	{
		double LLstart = gradient.getCorrectedLL();
		if (PRINT_OPTIMIZATION_MESSAGES)
		{
			System.out.println("#*EM.o start "+LLstart);
		}
		
		double LLprev = LLstart;
		double[] xprev = getSurvivalParameters();
		int iter = 0;
		while(iter<maxiter)
		{
			EStatistics E = EStep();
			MStep(E);
			double LLnow = gradient.getCorrectedLL();
			
			double diff = LLprev-LLnow;
			
			double delta = -diff/LLprev;
			
			double[] x = getSurvivalParameters();
			
			// test convergence on x 
			double max_xd = 0.0;
			for (int i=0; i<x.length; i++)
			{
				double xdiff = Math.abs(x[i]-xprev[i]); // p,q less than 1; r is max O(1) 
				double xd = xdiff/Double.max(x[i],1.0);
				max_xd = Double.max(max_xd, xdiff);
			}
			// test convergence on gradient 
			double[] dLdx = gradient.getCorrectedGradient();
			double max_dL=0.0;
			double score = -LLnow; // 0<score
			for (int i=0; i<x.length; i++)
			{
				double delta_dL = Math.abs(dLdx[i])*Double.max(x[i],1.0); // all x[i] >= 0.0 
				double dL = delta_dL/score;
				max_dL = Double.max(max_dL, dL);
			}
			if (PRINT_OPTIMIZATION_MESSAGES)
			{
				double len_dL = FunctionMinimization.euclideanNorm(dLdx);
				double len_x = FunctionMinimization.euclideanNorm(x);
				double rel_lendL = len_dL/len_x;

				System.out.println("#*EM.o "+iter+"\tLL "+LLnow+"\tdiff "+diff+"\trdiff "+delta
						+"\tx "+len_x+"\tdL "+len_dL+"\trdL "+rel_lendL
						+"\tmaxxd "+max_xd
						+"\tmaxdL "+max_dL); 
			}
			
			if (LLnow<LLprev) // dubious M step
			{
				System.out.println("#*EM.o done/decrease ("+iter+")"+"\t LL "+LLprev+"\tdiff "+(LLstart-LLprev));
				setSurvivalParameters(xprev);
				break;
			}
			LLprev = LLnow;
			xprev = x;
			if (max_xd<=FunctionMinimization.DFP_TOLX || max_dL<eps)
			{
				System.out.println("#*EM.o done/converged ("+iter+")"+"\t LL "+LLprev+"\tdiff "+(LLstart-LLprev));
				break;
			}
			
			iter++;
		}
		if (iter==maxiter && maxiter != 0)
			System.out.println("#*EM.o done/iterations ("+iter+")"+"\t LL "+LLprev+"\tdiff "+(LLstart-LLprev));
		return -LLprev;
	}
	
	private double[] invertSurvivalParameters(double[] survival_params)
	{
		
		final IndexedTree tree = gradient.factory.tree;
		final TreeWithRates rates = gradient.factory.rates;
		int num_nodes = tree.getNumNodes();
		double[] params = new double[survival_params.length];
		assert (params.length == 3*num_nodes);
		
		int node = num_nodes;
		while (node>0)
		{
			--node;
			double e;
			if (tree.isLeaf(node))
			{
				e = 0.0;
			} else
			{
				e = 1.0;
				int nc = tree.getNumChildren(node);
				for (int c=0; c<nc; c++)
				{
					int child = tree.getChild(node, c);
					double pc = survival_params[3*child+PARAMETER_LOSS];
					e *= pc;
				}
			}
			double p = survival_params[3*node+PARAMETER_LOSS];
			double q = survival_params[3*node+PARAMETER_DUPLICATION];
			double r = survival_params[3*node+PARAMETER_GAIN];
			double eq1 = e*(1.0-q);
			double denom = 1.0-eq1;
			p = (p-eq1)/denom;
			q = q/denom;
			
//			if (PRINT_OPTIMIZATION_MESSAGES && (p<0.0 || q>1.0))
//				System.out.println("#**EM.iSP node "+node+"\tp "+p+"\t/"+survival_params[3*node+PARAMETER_LOSS]+"\tq "+q+"\tr "+r
//						+"\te "+e
//						+"\teq1 "+eq1
//						+"\t// "+gradient.factory.rates.toString(node));	
			
			double min_p = Double.min(rates.getLossParameter(node), 1.0-MAX_PROB_NOT1);
//			p = Double.max(p, min_p);
			if (p<=0.0) p =min_p; // not ideal ... but if likelihood increases, then the iterations stop
			double min_q = Double.min( rates.getDuplicationParameter(node), (1.0-MAX_PROB_NOT1)*(1.0-MAX_PROB_NOT1));
//			q = Double.max(q,min_q);
			if (q<=0.0) q = min_q;
			
			assert (p<=1.0);
			assert (q<1.0);
			assert (0.0<=p);
			assert (0.0<=q);
			
			
			params[3*node+PARAMETER_LOSS] = p;
			params[3*node+PARAMETER_DUPLICATION] = q;
			params[3*node+PARAMETER_GAIN] = r;
			
		}
		return params;
		
	}
	
	private void MStep(EStatistics ES)
	{
		
		final IndexedTree tree = gradient.factory.tree;
		int num_nodes = tree.getNumNodes();
//		int nF = getTotalFamilyCount();
		double[] survival_params=new double[3*num_nodes];
		
		// M-estimate of loss:
		int num_leaves = tree.getNumLeaves();
		int u = num_nodes;
		while (num_leaves<u)
		{
			--u;
			int nc = tree.getNumChildren(u);
			assert (nc == 2);
			int v = tree.getChild(u, 0);
			int w = tree.getChild(u, 1);
			NodeStatistics vs = ES.node_stats[v];
			NodeStatistics ws = ES.node_stats[w];
//			double Nu_Sv = mean(vs.Nu_Sv);
//			double Nu_Sw = mean(ws.Nu_Sv);
			double Sv = mean(vs.Sv);
			double Sw = mean(ws.Sv);
			double Nu = mean(ES.node_stats[u].Nv);
			double[] tNu = tail(ES.node_stats[u].Nv);
			double[] tSv = tail(vs.Sv);
			double[] tSw = tail(ws.Sv);
			double[] tNu_Sv = tNu.clone();
			for (int i=0; i<tSv.length; i++)
			{
				tNu_Sv[i] = Double.max(0.0,tNu[i]-tSv[i]); // shouldn't be <0.0 but numerical errors 
			}
			double[] tNu_Sw = tNu.clone();
			for (int i=0; i<tSw.length; i++)
			{
				tNu_Sw[i] = Double.max(0.0,tNu[i]-tSw[i]);
			}
			double Nu_Sv = sum(tNu_Sv); // mean
			double Nu_Sw = sum(tNu_Sw); // mean
			
			if (Nu<Sv || Nu<Sw) 
			{
				// debug
				double mNu = sum(tNu);
				double mSv = sum(tSv);
				double mSw = sum(tSw);
				
				System.out.println("#**EM.M loss "+u
						+"\tNu "+Nu+"\t/"+mNu
						+"\tSv "+Sv+"\t/"+mSv
						+"\tSw "+Sw+"\t/"+mSw
						+"\tNu-Sv "+(Nu-Sv)+"\t/s "+Nu_Sv+"\tm "+(mNu-mSv)
						+"\tNu-Sw "+(Nu-Sw)+"\t/s "+Nu_Sw+"\tm "+(mNu-mSw)
						);
			}
			
			double pv;
			double pw;
			
			if (optimize_node[v])
			{
				
				if (optimize_node[w])
				{
					
					pv = Nu_Sv/Sw;
					pw = Nu_Sw/Sv;
				} else
				{
					pw = gradient.factory.getLossParameter(w); // survival loss parameter
					Nu_Sw = Nu-pw*Sv;
					pv = Nu_Sv/Nu_Sw;
				}
			} else if (optimize_node[w])
			{
				pv = gradient.factory.getLossParameter(v); // survival loss parameter
				Nu_Sv = Nu-pv*Sw;
				pw = Nu_Sw/Nu_Sv;
			} else
			{
				pv = gradient.factory.getLossParameter(v);
				pw = gradient.factory.getLossParameter(w);
			}
			survival_params[3*v+PARAMETER_LOSS] = pv;
			survival_params[3*w+PARAMETER_LOSS] = pw;
			
			if (PRINT_OPTIMIZATION_MESSAGES)
			{
//				
//				System.out.println("#*EM.M loss v,w "+v+","+w
//						+"\tpv "+pv
//						+"\tpw "+pw
//						+"\tNuSv "+Nu_Sv+" ("+Nu_Sv2+")"
//						+"\tNuSw "+Nu_Sw+" ("+Nu_Sw2+")"
//						+"\tNu "+Nu
//						+"\tSv "+Sv
//						+"\tSw "+Sw
//						+"\t//v "+factory.rates.toString(v)
//						+"\t//w "+factory.rates.toString(w));
			
			}
			
			if (pv<0.0)
			{
				System.out.println("#**EM.M node "+v+"\tpv "+pv+"\tNSw "+Nu_Sw+"\tNSv "+Nu_Sv
						+"\tNu "+Nu+"\tSv "+Sv+"\tSw "+Sw
						+"\tov "+optimize_node[v]+"\tow "+optimize_node[w]
						+"\t// v "+tree.toString(v)+"\tu "+tree.toString(u));
			}
			
			assert (0.0<=pv);
			assert (0.0<=pw);
			
		}
		
		survival_params[3*tree.getRoot()+PARAMETER_LOSS]=1.0; 
		
		// now set q and kappa
		
		double kappa_tol = 1.0/(1<<24);
		
		u = 0;
		double Lunobs = gradient.getUnobservedLL();
		double denom =  -Math.expm1(Lunobs); // ==(1-Pr(unobs))
		final double ncorr = gradient.getTotalFamilyCount()/denom;
		while (u<num_nodes)
		{
			if (optimize_node[u])
			{
				
				
				NodeStatistics stats = ES.node_stats[u];
				
				//   q = (N-S)/(N+k*F')
				//   1-q=(S+k*F')/(N+k*F') with F' = nF / (1-p(unobs))
//				final double Nu_Su = mean(stats.Nv_Sv);
				final double Nu = mean(stats.Nv);
				final double Su = mean(stats.Sv);
				double Nu_Su2 = Nu-Su;
	//			if (PRINT_OPTIMIZATION_MESSAGES)
	//				System.out.println("#*EM.M node "+u+"\tNuSu "+Nu_Su+"\tNuSu2 "+Nu_Su2);
										
				double[] tailN = tail(stats.Nv);
				double[] tailS = tail(stats.Sv);
				
				//final double[] tail = tail(stats.Nv_Sv);
				double[] tail = tailN.clone();
				for (int i=0; i<tailS.length; i++)
				{
					tail[i] -= tailS[i];
				}
				
				DoubleFunction<Double> qpar =
						new DoubleFunction<>()
						{
							@Override
							public Double apply(double κ)
							{
								double q = Nu_Su2/(Nu+κ *ncorr);
								return q;
							}
						};
				
				DoubleFunction<Double> dL =
						new DoubleFunction<>()
						{
							@Override
							public Double apply(double κ)
							{
								++calls_kappa;
								
								double sum_tails=0.0;
								for (int i=0; i<tail.length; i++)
								{
									sum_tails += tail[i]/(κ+i);
								}
								double q = qpar.apply(κ);
								double qterm = ncorr*Math.log1p(-q);
								double dLdκ = sum_tails + qterm;
								
	//							if (PRINT_OPTIMIZATION_MESSAGES)
	//								System.out.println("#*EM.M "+calls_kappa+"\t"+κ+"\tdL "+dLdκ+"\tst "+sum_tails+"\tqt "+qterm+"\tq "+q);							
								
								return dLdκ ; 
							}
						};
				
				//	bracketing
				double kappa1;
				kappa1 = 1.0/128.0;
				
				double d1 = dL.apply(kappa1);
				double kappa2, d2;
				
				int iter=0;
				
				if (0.0<d1)
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
	//					System.out.println("#**EM.M node "+u+"\tsetting max kappa "+kappa1);
					}
					
				} else if (d1==0.0)
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
					}
				}
				double κ;
				
				if (kappa1==kappa2)
				{
					κ = kappa1;
				} else
				{
					κ = FunctionMinimization.zbrent(dL, kappa1, kappa2, kappa_tol);
				}
				double q = qpar.apply(κ);
				survival_params[3*u+PARAMETER_GAIN]=κ;
				survival_params[3*u+PARAMETER_DUPLICATION]=q;
				
	//			if (PRINT_OPTIMIZATION_MESSAGES)
	//				System.out.println("#*EM.M node "+u
	//						+"\tp "+survival_params[3*u+PARAMETER_LOSS]
	//						+"\tq "+survival_params[3*u+PARAMETER_DUPLICATION]
	//						+"\tr "+survival_params[3*u+PARAMETER_GAIN]								
	//						+"\t// "+ factory.rates.toString(u)
	//						+"\t// "+factory.tree.toString(u)
	//					);
			} else
			{
				survival_params[3*u+PARAMETER_GAIN]=gradient.factory.getGainParameter(u);
				survival_params[3*u+PARAMETER_DUPLICATION]=gradient.factory.getDuplicationParameter(u);
			}
			u++;
		}
		
		double[] distribution_params = invertSurvivalParameters(survival_params);
		while (0<u)
		{
			--u;
			double r = distribution_params[3*u+PARAMETER_GAIN]	;
			double p = distribution_params[3*u+PARAMETER_LOSS]	;
			double q = distribution_params[3*u+PARAMETER_DUPLICATION]	;
			
//			q=factory.rates.getDuplicationParameter(u);
//			r=factory.rates.getGainParameter(u);
			
			gradient.factory.rates.setParameters(u, r, p, q);		
		}
		gradient.computeParameters();
	}
	
//	
//	private final NodeStatistics[] node_stats;
	
	private EStatistics getSampleStats()
	{
		int nF = gradient.factory.table.getFamilyCount();
		final int unit_task;
		ForkJoinPool thread_pool = Count.threadPool();
		unit_task = Count.unitTask(nF);
//		if (thread_pool != null)
//		{
//			unit_task = THREAD_UNIT_TASK;
//		}
//		else
//		{
//			unit_task = nF; // do not fork
//		}

		class PartialE extends RecursiveTask<EStatistics>
		{
			/**
			 * First family index.
			 */
			private final int minF;
			/**
			 * Last family index, exclusive.
			 */
			private final int maxF;

			PartialE(int min, int max)
			{
				this.minF = min;
				this.maxF = max;
			}
			@Override
			protected EStatistics compute() 
			{
				try
				{
					EStatistics ES;
					if (maxF-minF>unit_task)
					{
						int medF = (minF+maxF)/2;
						PartialE left = new PartialE(minF, medF);
						PartialE right = new PartialE(medF, maxF);
						left.fork();
						ES = right.compute();
						EStatistics ES2 = left.join();
						ES.addStats(ES2);
					} else
					{
						ES  = new EStatistics();
						for (int f = minF; f<maxF; f++)
						{
							Gradient.Profile G = gradient.getGradient(f);
							int mul = gradient.getMultiplicity(f);
							ES.addProfile(G.post, mul);
						}
					}
					return ES;
				} catch (Throwable t)
				{
					throw t;
				}
			}
		}

		PartialE bigjob = new PartialE(0, nF);
		EStatistics ES;
		try
		{
			if (nF>unit_task)
			{
				ES = thread_pool.invoke(bigjob);
			} else
			{
				ES = bigjob.compute();
			}
		} catch (Throwable t)
		{
			throw new RuntimeException(t);
		}
		return ES;
	}
	
	
	public static void main(String[] args) throws Exception
	{
		PRINT_OPTIMIZATION_MESSAGES = true;
		
		count.io.CommandLine cli = new count.io.CommandLine(args, EM.class);

		PrintStream out = System.out;
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(EM.class));
    	    out.println(CommandLine.getStandardRuntimeInfo(EM.class, args));
    	}
    	
    	GammaInvariant model = null; 
    	TreeWithRates starting_rates;
    	if (cli.getModel()==null)
    	{
    		Random RND = cli.getOptionRND(out);
//    		if (cli.getOptionValue(OPT_RND)!=null)
//    		{
//    			int rnd_seed = cli.getOptionInt(OPT_RND, 0);
//    			RND = (rnd_seed==0?new Random():new Random(rnd_seed));
//    			out.println(CommandLine.getStandardHeader("Random initialization: -"+OPT_RND+" "+rnd_seed));    			
//    		}
    		model = GammaInvariant.nullModel(cli.getTree(), RND);
    		starting_rates = model.getBaseModel();
//    		model = GammaInvariant.nullModel(cli.getTree());
//    		starting_rates = model.getBaseModel();
    		
    	} else
    	{
    		starting_rates = cli.getRates(); 
    		model = cli.getModel();
    	}

    	{
    		count.matek.DiscreteDistribution root_prior = cli.getOptionDistribution(OPT_MODEL_ROOT_PRIOR, null);
    		if (root_prior != null)
    		{
    			starting_rates.setRootDistribution(root_prior);
    			out.println(CommandLine.getStandardHeader("Root prior set: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
    		}
			out.println(CommandLine.getStandardHeader("(Root prior get: "+starting_rates.getRootDistribution()+")"));
    	}
    	
    	AnnotatedTable table = cli.getTable();
    	
    	EM M = new EM(starting_rates, table);//new UniqueProfileTable(table));
    	
        int min_copies = Integer.min(2, table.minCopies());
        min_copies = cli.getOptionInt(OPT_MINCOPY, min_copies);
		M.setMinimumObservedCopies(min_copies);
		out.println(CommandLine.getStandardHeader("Minimum observed copies: -"+OPT_MINCOPY+" "+min_copies));
		
		int maxiter = cli.getOptionInt(OPT_ROUNDS, Integer.MAX_VALUE);
		double eps = cli.getOptionDouble(OPT_EPS, 0.00125);
        
        out.println(CommandLine.getStandardHeader("Convergence: -"+OPT_ROUNDS+" "+maxiter+" -"+OPT_EPS+" "+eps));
    	
        double score =
        		M.optimize(eps, maxiter);
        double ascore = score/table.getFamilyCount();

		out.println(count.io.RateVariationParser.printRates(model));
        out.println("#SCORE "+score);
        out.println("#AVGSCORE "+ascore);
        
        
		if (out != System.out)
		{
			if (out.checkError())
				throw new java.io.IOException("Write failed.");
			out.close();
		}
        
		
	}
}
