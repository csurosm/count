package count.model;
/*
 * Copyright 2025 Mikl&oacute;s Cs&#369;r&ouml;s.
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

import static count.io.CommandLine.OPT_ANCESTRAL;
import static count.io.CommandLine.OPT_OUTPUT;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;

import count.Count;
import count.ds.AnnotatedTable;
import count.ds.ProfileTable;
import count.io.CommandLine;
import count.matek.Logarithms;
import count.matek.PseudoRandom;


/**
 * Viterbi reconstruction for most likely copy number history.
 * 
 * 
 */
public class ViterbiHistory extends LikelihoodParametrized {
	private static final boolean DEBUG_CALCULATIONS = false;
	
	public ViterbiHistory(TreeWithRates rates, ProfileTable table) {
		super(rates, table);
	}
	
	class HistoryProfile {
		HistoryProfile(int family_idx){
			this.family_idx = family_idx;
			int num_nodes = tree.getNumNodes();
			this.node_likelihoods = new double[num_nodes][];
			this.edge_likelihoods = new double[num_nodes][];
			this.node_trace = new int[num_nodes][][];
			this.edge_trace = new int[num_nodes][];
			for (int node=0; node<num_nodes; node++)
				this.setCalculationWidth(node);
			
			computeAll();
		}
		protected final int family_idx;
		private final double[][] node_likelihoods;
		private final int[][][] node_trace;
		private final double[][] edge_likelihoods;
		private final int[][] edge_trace;
		
		private int[] node_history;
		private int[] edge_history;
		
		private double prob_history;
		
		
		private int setCalculationWidth(int node){
			int w;
			if (tree.isLeaf(node)){
				int n = table.getFamilyProfile(family_idx)[node];
				w=1+Integer.max(-1, n);
			} else {
				int num_children = tree.getNumChildren(node);
				int n=0;
				int ambi=0;
				for (int ci=0; ci<num_children; ci++)
				{
					int child  = tree.getChild(node, ci);
					int cn = node_likelihoods[child].length-1;
					if (cn<0) ambi++;
					else n += cn;
				}
				if (ambi == num_children) w=0;
				else
					w=n+1;
			}
			node_likelihoods[node]=new double[w];
			return  w;
		}
		
		private void computeAll() {
			for (int node=0; node<tree.getNumNodes(); node++){
				computeNode(node);
				computeEdge(node);
			}
			computeHistory();
		}
		
		private void computeNode(int u) {
			double[] C = node_likelihoods[u];
			int m = C.length;
			if (m==0) return;
			
			if (tree.isLeaf(u)){
				Arrays.fill(C, Double.NEGATIVE_INFINITY);
				C[m-1]=0.0;
			} else {
				int nc = tree.getNumChildren(u);
				assert (nc==2);
				int v = tree.getChild(u, 0);
				int w = tree.getChild(u, 1);
				double log1_pv = getLogLossComplement(v);
				double log1_pw = getLogLossComplement(w);
				double logx10 = log1_pv + getLogLossParameter(w);
				double logx11 = log1_pv + log1_pw;
				double logx01 = getLogLossParameter(v)+log1_pw;
				double logy = Logarithms.add(logx11, Logarithms.add(logx01, logx10));
				double logp10 = logx10-logy; // survival on left only
				double logp11 = logx11-logy; // survival in both
				double logp01 = logx01-logy; // survival on right only
				
				double[] Kv = edge_likelihoods[v];
				double[] Kw = edge_likelihoods[w];
				
				
				
				node_trace[u]=new int[m][2];
				
				// brutal triple loop 
				for (int n=0; n<m ; n++) { 
					C[n]=Double.NEGATIVE_INFINITY;
					node_trace[u][n][0] = node_trace[u][n][1] = -1; // trace to v and w 
					
					int smax = Integer.min(n,Kv.length-1);
					int tmax = Integer.min(n,Kw.length-1);
					
					for (int s=0; s<=smax; s++) { // inherited on left
						for (int t=n-s; t<=tmax; t++) { // inherited on right
							int k = s+t-n; // both inherit
							//	inheritance n -> 
							// 		(n-t) only on left,
							//		s+t-n in both
							// 		(n-s) only on right
							double binom = factorials.factln(n)-factorials.factln(n-t)-factorials.factln(k)-factorials.factln(n-s);
							double logt10 = n==t?0.0:(n-t)*logp10;
							double logt11 = k==0?0.0:k*logp11;
							double logt01 = n==s?0.0:(n-s)*logp01;
							double z = Kv[s]+Kw[t]+logt10+logt11+logt01+binom;
							if (C[n]<z) {
								C[n]=z;
								node_trace[u][n][0]=s;
								node_trace[u][n][1]=t;
							}
						}
					}
				}
				if (DEBUG_CALCULATIONS) {
					System.out.print("#**VH.cN "+family_idx+"/"+u);
					for (int n=0; n<m; n++) {
						System.out.printf("\t%.2g(%d,%d)", node_likelihoods[u][n],node_trace[u][n][0],node_trace[u][n][1]);
					}
					System.out.println();
				}
			}
		}
		
		private void computeEdge(int u) {
			double[] K; // the edge likelihoods that will be set 
			double[] C = node_likelihoods[u];
			int m = C.length;
			assert (m!=0);
			
			if (tree.isRoot(u)){
				if (getLogLossComplement(u)==Double.NEGATIVE_INFINITY) //   logp==0.0) 
					K = new double[1];
				else 
					K = new double[2];
			} else {
				K = new double[m];
			}
			edge_likelihoods[u] = K; 
			
			edge_trace[u] = new int[K.length];
			
			double logq = getLogDuplicationParameter(u);
			double log1_q = getLogDuplicationComplement(u);
			
			if (logq==Double.NEGATIVE_INFINITY){
				// Poisson distribution with parameter r 
				double logr = getLogGainParameter(u);
				double r = Math.exp(logr);
				if (logr == Double.NEGATIVE_INFINITY){
					// no gain, no duplication ??
					for (int s=0; s<K.length; s++) {
						K[s] = C[s];
						edge_trace[u][s]=s;
					}
				} else {
					for (int s=0; s<K.length; s++) {
						K[s] = Double.NEGATIVE_INFINITY;
						edge_trace[u][s]=-1;
						
						int t=0, n=s;
						while (n<C.length) {
							double t_logr = (t==0?0.0:t*logr);
							double z = C[n]-r + t_logr -factorials.factln(t) ;
							
							if (K[s]<z) {
								K[s]=z;
								edge_trace[u][s]=n;
							}
							++n;
							++t;
						}
					}
				}
			} else {
				// Polya distribution with parameters kappa+s and q 
				//double κ = g;
				double log_κ = getLogGainParameter(u);
				assert (log_κ  != Double.NEGATIVE_INFINITY);
				double loglog1_q = Logarithms.logitToLogLogComplement(logq-log1_q);
				double log_kappa_log1_q = log_κ+loglog1_q;
				double kappa_log1_q = Math.exp(log_kappa_log1_q);
				
				for (int s=0; s<K.length; s++){
					K[s] = Double.NEGATIVE_INFINITY;
					edge_trace[u][s]=-1;
	
					double factln_s = gain_factorials[u].factln(s);
					int t=0, n=s+t;
					while (n<C.length) {
						double binom = gain_factorials[u].factln(n)-factln_s-factorials.factln(t);
						double t_logq = (t==0?0.0:t*logq);
						double ks_log1_q = s*log1_q-kappa_log1_q;
						double z = C[n]+binom + ks_log1_q + t_logq;
						
						if (K[s]<z) {
							K[s]=z;
							edge_trace[u][s]=n;
						}
						++n;
						++t;
					}
				} // for
			} // else 
			if (DEBUG_CALCULATIONS) {
				System.out.print("#**VH.cE "+family_idx+"/"+u);
				for (int s=0; s<K.length; s++) {
					System.out.printf("\t%.2g(%d)", K[s], edge_trace[u][s]);
				}
				System.out.println();
			}
		} // computeEdge
		
		
		/**
		 * Traceback for optimal reconstruction
		 */
		private void computeHistory() {
			node_history = new int[tree.getNumNodes()];
			edge_history = new int[node_history.length];
			Arrays.fill(node_history, -1);
			Arrays.fill(edge_history, -1);
			int u=tree.getRoot();
			{
				double[] K = edge_likelihoods[u];
				int s = (K.length==1||K[1]<=K[0])?0:1;
				edge_history[u]=s;
				node_history[u]=edge_trace[u][s];
				
				Likelihood.Profile lp = ViterbiHistory.this.getProfileLikelihood(family_idx);
				double LL = lp.getLogLikelihood();
				this.prob_history = Math.exp(K[s]-LL);
				if (Count.THREAD_PARALLELISM==1)
					System.out.println("#**VH.HP.cH\t"+family_idx+"\tprob "+this.prob_history+"\t("+K[s]+","+LL+")\t// @ "+Thread.currentThread());
			}
			assert (tree.getNumChildren(u)==2);
			for (int c=0; c<2; c++) {
				int v = tree.getChild(u, c);
				int n = node_history[u];
				edge_history[v] = node_trace[u][n][c];
			}
			// trace back in prefix traversal 
			while (0<u){
				--u;
				int s = edge_history[u];
				int n = node_history[u] = edge_trace[u][s];
				int nc = tree.getNumChildren(u);
				assert (nc==2 ||nc==0);
				for (int c=0; c<nc; c++) {
					int v = tree.getChild(u, c);
					edge_history[v] = node_trace[u][n][c];
				}
			}
			
			

		}
		
		public int getBirth(int u) {
			int n = node_history[u];
			int s = edge_history[u];
			return n-s;
		}
		public int getDeath(int v) {
			int s = edge_history[v];
			int n;
			if (tree.isRoot(v)) {
				double[] K=edge_likelihoods[v];
				n = K.length-1;
			} else {
				int u = tree.getParent(v);
				n = node_history[u];
			}
			return n-s;
		}
		public int getNodeCount(int u) {
			return node_history[u];
		}
		public int getEdgeCount(int u) {
			return edge_history[u];
		}
		
		/**
		 * Whether all copies are inherited exclusively left or right
		 * @param u
		 * @return
		 */
		public boolean isFamilySplit(int u) {
			int inherited_copies=0;
			int inherited_lineages=0;
			for (int c=0; c<tree.getNumChildren(u); c++) {
				int v = tree.getChild(u, c);
				int s = edge_history[v];
				inherited_copies += s;
				inherited_lineages += s==0?0:1;
			}
			return 1<inherited_lineages && inherited_copies == node_history[u];
		}
	} // Profile 
	
	class LineageStatistics {
		final int[] lineage_family_present;
		final int[] lineage_family_multi;
		final int[] lineage_family_gain ;
		final int[] lineage_family_loss;
		final int[] lineage_family_expand;
		final int[] lineage_family_contract;
		final int[] lineage_copies_count;
		final int[] lineage_copies_birth;
		final int[] lineage_copies_death; 
		final int[] lineage_family_splits;
		
		LineageStatistics(){
			int num_nodes = tree.getNumNodes();
			lineage_family_present =  new int[num_nodes];
			lineage_family_multi =  new int[num_nodes];
			lineage_family_gain =  new int[num_nodes];
			lineage_family_loss =  new int[num_nodes];
			lineage_family_expand =  new int[num_nodes];
			lineage_family_contract =  new int[num_nodes];
			lineage_copies_birth =  new int[num_nodes];
			lineage_copies_death =  new int[num_nodes];
			lineage_copies_count =  new int[num_nodes];
			lineage_family_splits = new int[num_nodes];
		}
		
		
		void add(int family_idx, int mul) {
//			System.out.println("#**VH.add\t"+family_idx+"\tstart\t// @ "+Thread.currentThread());
			HistoryProfile H = new HistoryProfile(family_idx);
			int v = tree.getRoot();
			do {
				int n = H.getNodeCount(v);
				int s = H.getEdgeCount(v);
				int b = H.getBirth(v);
				int d = H.getDeath(v);
				
				int fp = 0<n?1:0;
				int fm = 1<n?1:0;
				int fg = (0<n && s==0)?1:0;
				int fl, fe, fc;
				if (tree.isRoot(v)) {
					fl = 0;
					fe = (1<n && s==1)?1:0;
					fc = 0;
				} else {
					int u = tree.getParent(v);
					int nu = H.getNodeCount(u);
					fl = (0<nu && s==0)?1:0;
					fc = (1<nu && 1==n && 1==s)?1:0;
					fe = (1==nu && 1==s && 1<n)?1:0;
				}
				int sp = H.isFamilySplit(v)?1:0;
				lineage_family_present[v] 	+= fp*mul;
				lineage_family_multi[v] 	+= fm*mul;
				lineage_family_gain[v]   	+= fg*mul;
				lineage_family_loss[v]  	+= fl*mul;
				lineage_family_expand[v]	+= fe*mul;
				lineage_family_contract[v]	+= fc*mul;
				lineage_copies_count[v]  	+= n*mul;
				lineage_copies_birth[v] 	+= b*mul;
				lineage_copies_death[v]  	+= d*mul;
				lineage_family_splits[v]	+= sp*mul;
				
				--v;
			} while (0<=v);
//			System.out.println("#**VH.add\t"+family_idx+"\tdone\t// @ "+Thread.currentThread());
		}
		void add(LineageStatistics that) {
			int u = tree.getRoot();
			while (0<u) {
				this.lineage_family_present[u] += that.lineage_family_present[u];
				this.lineage_family_multi[u] += that.lineage_family_multi[u];
				this.lineage_family_gain[u] += that.lineage_family_gain[u];
				this.lineage_family_loss[u] += that.lineage_family_loss[u];
				this.lineage_family_expand[u] += that.lineage_family_expand[u];
				this.lineage_family_contract[u] += that.lineage_family_contract[u];
				this.lineage_copies_count[u] += that.lineage_copies_count[u];
				this.lineage_copies_birth[u] += that.lineage_copies_birth[u];
				this.lineage_copies_death[u] += that.lineage_copies_death[u];
				this.lineage_family_splits[u] += that.lineage_family_splits[u];
				--u;
			}
		}
	}
	
	
	protected LineageStatistics computeLineageStatistics() {
		int nF = table.getFamilyCount();
		final int unit_task = Count.unitTask(nF);
		ForkJoinPool thread_pool = Count.threadPool();
		
		// random order to mix easy and difficult profiles 
		final int[] familyOrder;
		if (nF <= unit_task) {
			familyOrder = new int[nF];
			for (int f=0; f<nF; f++)
				familyOrder[f]=f;
		} else 
			familyOrder = PseudoRandom.randomPermutation(nF);
	
		
		class PartialS extends RecursiveTask<LineageStatistics>{

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
			protected LineageStatistics compute() {
				try
				{
					LineageStatistics compute;
					if (maxF-minF > unit_task)
					{
						int medF = (minF+maxF)/2;
						PartialS left = new PartialS(minF, medF);
						PartialS right = new PartialS(medF, maxF);
						left.fork();
						compute = right.compute();
						LineageStatistics c2 = left.join();
						compute.add(c2);
					} else
					{
						compute = new LineageStatistics();
						for (int k=minF; k<maxF; k++)
						{
							int f = familyOrder[k];
							
							compute.add(f,1);
						}
					}
					return compute;
				} catch (Throwable t)
				{
					throw t;
				}
			}
		} // RecursiveTask
		
		PartialS bigjob = new PartialS(0,nF);
		LineageStatistics getLineageStatistics;
		try {
			if (unit_task<nF) {
				getLineageStatistics = thread_pool.invoke(bigjob);
			} else {
				getLineageStatistics = bigjob.compute();
			}
			return getLineageStatistics;
		} catch (Throwable oops) {
			throw new RuntimeException(oops);
		}
	}
	
	
	private void printLineageStatistics(PrintStream out) {
		LineageStatistics stats = computeLineageStatistics();
		String prefix_ancestral = OPT_ANCESTRAL.toUpperCase();
		out.println(prefix_ancestral+"\tnode\tpresence:p\tmulti:m\tcorrected:p.\tgain:g\tloss:l\texpand:++\tcontract:--\tscopies:n\ttcopies:n.\tsplits:sp\tbirth:+.\tdeath:-.");
		int root = tree.getRoot();
		for (int node=0; node<=root; node++) {
			int fp = stats.lineage_family_present[node];
			int fm = stats.lineage_family_multi[node];

			int fgain = stats.lineage_family_gain[node];
			int floss = stats.lineage_family_loss[node];
			int fexpand = stats.lineage_family_expand[node];
			int fcontract = stats.lineage_family_contract[node];
			
			int fsurviving = stats.lineage_copies_count[node];
			
			int fbirth = stats.lineage_copies_birth[node];
			int fdeath = stats.lineage_copies_death[node];
			
			int fsplits = stats.lineage_family_splits[node]; 
			
			out.printf("%s\t%s", prefix_ancestral, tree.getIdent(node));
			out.printf("\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" 
						, fp, fm, fp
						, fgain, floss, fexpand, fcontract
						, fsurviving, fsurviving
						, fsplits
						, fbirth, fdeath);
			
		}
		
	}
	
	public static void main(String[] args) throws Exception {
		Class<?> these = java.lang.invoke.MethodHandles.lookup().lookupClass();
		CommandLine cli = new CommandLine(args, these);
        AnnotatedTable table = cli.getTable();
        TreeWithRates rates = cli.getRates();
		ViterbiHistory VH = new ViterbiHistory(rates, table);
		PrintStream out = System.out;
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(these));
    	    out.println(CommandLine.getStandardRuntimeInfo(these, args));
    	}
		
		VH.printLineageStatistics(out);
		
		if (out_file!=null) out.close();
	}
	
}
