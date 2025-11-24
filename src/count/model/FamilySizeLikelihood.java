package count.model;

import static count.io.CommandLine.OPT_OUTPUT;

import java.io.PrintStream;
import java.util.Arrays;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.Phylogeny;
import count.ds.ProfileTable;
import count.io.CommandLine;
import count.matek.Logarithms;

/**
 * Likelihood by family size for unobserved copies
 * 
 */
public class FamilySizeLikelihood extends LikelihoodParametrized {
	
	public FamilySizeLikelihood(TreeWithLogisticParameters rates, int maxsize) {
		super(rates, maxCopyTable(rates.getTree(), maxsize));
		this.sizeProfiles = new SizeProfile[maxsize+1];
		initDataStructures();
	}
	
	/**
	 * Decoy table to preset max family size at all nodes 
	 * from instantiation. The constructed table has one profile 
	 * for each leaf with max copies. 
	 * 
	 * @param tree
	 * @param max
	 * @return
	 */
	private static ProfileTable maxCopyTable(IndexedTree tree, int max) {
		ProfileTable maxCopyTable = new ProfileTable() {

			@Override
			public int getFamilyCount() {
				return tree.getNumLeaves();
			}

			@Override
			public int getTaxonCount() {
				return tree.getNumLeaves();
			}

			@Override
			public int[] getFamilyProfile(int fam) {
				int[] profile = new int[getTaxonCount()];
				profile[fam]=max;
				return profile;
			}
			@Override
			public String[] getTaxonNames()
			{
				return tree.getLeafNames();
			}
		};
		return maxCopyTable;
	}
	
	private final SizeProfile[] sizeProfiles;
	
	public int maxCopies() { return sizeProfiles.length-1;}
	
	private void initDataStructures() {
		for (int m=0; m<sizeProfiles.length; m++) {
			sizeProfiles[m] = new SizeProfile(m);
			//System.out.println("#**FSL.iDS allocate SP["+m+"]");
		}
	}
	
	private void computeAll() {
		computeParameters();
		
		int u=0;
		while (u<tree.getNumNodes()) {
			for (int m=0; m<sizeProfiles.length; m++) {
				SizeProfile spm = sizeProfiles[m];
				spm.computeNode(u);
				spm.computeEdge(u);
			}
			//System.out.println("#**FSL.cA done/inner "+u+"\t"+tree.toString(u));
			
//			if (!tree.isLeaf(u)) { // free memory 
//				int nc = tree.getNumChildren(u);
//				for (int ci=0; ci<nc; ci++) {
//					int v = tree.getChild(u, ci);
//					for (int m=0; m<sizeProfiles.length; m++) {
//						SizeProfile spm = sizeProfiles[m];
//						//spm.destroyLikelihoods(v);
//					}
//				}
//			}
			
			++u;
		}
		// now outer likelihoods in prefix order 
		while (0<u) {
			--u;
			for (int m=0; m<sizeProfiles.length; m++) {
				SizeProfile spm = sizeProfiles[m];
				spm.computeOutside(u);
			}
			//System.out.println("#**FSL.cA done/outer "+u+"\t"+tree.toString(u));
		}
		
		// replace outer likelihoods by cumulative values on outer profile size 
		// in order to speed up posterior computations
		for (int m=1; m<sizeProfiles.length; m++) {
			SizeProfile prev = sizeProfiles[m-1];
			SizeProfile cumul = sizeProfiles[m];
			for (u=0; u<tree.getNumNodes(); u++) {
				double[] px = prev.edge_outside[u];
				double[] cx = cumul.edge_outside[u];
				assert (px.length<=cx.length);
				for (int s=0; s<px.length; s++) 
					cx[s] = Logarithms.add(cx[s],px[s]);
				px = prev.node_outside[u];
				cx = cumul.node_outside[u];
				assert (px.length<=cx.length);
				for (int n=0; n<px.length; n++) 
					cx[n] = Logarithms.add(cx[n],px[n]);
				double[][] pa = prev.edge_transitions[u];
				double[][] ca = cumul.edge_transitions[u];
				assert (pa.length<=ca.length);
				for (int n=0; n<pa.length; n++) {
					px = pa[n];
					cx = ca[n];
					assert (px.length<=cx.length);
					for (int s=0; s<px.length; s++) 
						cx[s] = Logarithms.add(cx[s],px[s]);
				}
				pa = prev.node_transitions[u];
				ca = cumul.node_transitions[u];
				assert (pa.length<=ca.length);
				for (int n=0; n<pa.length; n++) {
					px = pa[n];
					cx = ca[n];
					assert (px.length<=cx.length);
					for (int s=0; s<px.length; s++) 
						cx[s] = Logarithms.add(cx[s],px[s]);
				}
			}
		}
	}

	
	/**
	 * Computations for size-dependent likelihoods.
	 * {@link #computeEdge(int)} and {@link #computeNode(int)} are good as they are. 
	 */
	protected class SizeProfile extends Profile {
		
		protected SizeProfile(int size) {
			super(size);
			int num_nodes = tree.getNumNodes();
			this.pairing_likelihoods = new double[num_nodes][][];
			node_outside = new double[num_nodes][];
			node_transitions = new double[num_nodes][][];
			edge_outside = new double[num_nodes][];
			edge_transitions = new double[num_nodes][][];
		}
		
		/**
		 * Repurposing family index as family size
		 * 
		 * @return
		 */
		private int size() { return family_idx;}
		
		/**
		 * Calculation width is the same everywhere
		 */
		@Override
		protected int setCalculationWidth(int node) {
			setCalculationWidth(node, size()+1);
			return getCalculationWidth(node);
		}
		
		private double[][][] pairing_likelihoods;
		
		private final double[][] node_outside;
		private final double[][] edge_outside;
		
		private final double[][][] node_transitions;
		private final double[][][] edge_transitions;
		

//		protected void computePairingLikelihoods(int w) {
//			this.pairing_likelihoods[w] = getPairingLikelihoods(w);
//		}
		private double[][] getPairingLikelihoods(int w) {
			if (this.pairing_likelihoods[w]!=null) return this.pairing_likelihoods[w];
			
			double[][] Lw = new double[sizeProfiles.length][];
			double[] Kw = getEdgeLikelihoods(w);
			// log1_pw == log1_e
			final double log_pw = FamilySizeLikelihood.this.getLogLossParameter(w);
			final double log1_pw = FamilySizeLikelihood.this.getLogLossComplement(w);
			
			int m = size();
			
			for (int ell=0; ell<Lw.length; ell++) {
				int t =ell;
				Lw[ell] = new double[m+1];
				double x;
				if (t<=m) {
					x =  Lw[ell][t] = Kw[t];
				} else {
					x  = Double.NEGATIVE_INFINITY;
					t = m+1;
				}
				while (0<t) {
					--t;
					x += log1_pw;
					double y = Lw[ell-1][t]+log_pw;
					x = Lw[ell][t] = Logarithms.add(x,y);
				}
			}
			this.pairing_likelihoods[w] = Lw;
			return Lw;
		}
		
		
		
		@Override
		protected double[] computeNodeFromChildren(int u){
			if (tree.getNumChildren(u)!=2) throw new IllegalArgumentException("Only for binary nodes; "+u+" has "+tree.getNumChildren(u)+" children");
			
			int v = tree.getChild(u, 0);
			int w = tree.getChild(u, 1);
			
			final int m = size();
			double[] C = new double[m+1];
			
			final double logit_pv = FamilySizeLikelihood.this.getLogitSurvivalParameter(v, PARAMETER_LOSS);
			final double log1_e = FamilySizeLikelihood.this.getLogLossComplement(w);
			
			final double logit_p = logit_pv + log1_e;
			double log_p = Logarithms.logitToLogValue(logit_p);
			double log1_p = Logarithms.logitToLogComplement(logit_p);
			
			double[] terms = new double[m+1];
			
			for (int ell=0; ell<C.length; ell++) {
				double log_ellfact = factorials.factln(ell);
				int s = ell;
				int t = ell-s; // == 0
				
				while (0<=s) {
					double KvLw = Double.NEGATIVE_INFINITY;
					// need s<=m-mw and t<=mw (we have s+t=ell<=m)
					for (int mw=t; mw<=m-s; mw++) { // alas, cubic time bc of this inner loop 
						double[] Kv = sizeProfiles[m-mw].getEdgeLikelihoods(v);
						double[][] Lmw = sizeProfiles[mw].getPairingLikelihoods(w);
						double[] Lw = Lmw[ell];
						KvLw = Logarithms.add(Kv[s]+Lw[t], KvLw);
					}
					
					double binom = log_ellfact - factorials.factln(s)-factorials.factln(t); // ell chose s 
					double slog1_p = s==0?0.0:s*log1_p;
					double tlog_p = t==0?0.0:t*log_p;
					
					terms[t] = KvLw + binom + slog1_p + tlog_p;
					
					--s;
					++t;
				}
				C[ell] = Logarithms.sum(terms, t);
			}
			
			
			return C;
			
			
		}		
		
		
		
		@Override
		public void computeLikelihoods() {
//			FamilySizeDistribution.this.computeAll();
			throw new UnsupportedOperationException();
		}
		
		@Override
		protected void destroyLikelihoods(int node) {
			super.destroyLikelihoods(node);
			this.pairing_likelihoods[node]=null;
		}
		
		
		protected void computeOutside(int  node){
			int M=sizeProfiles.length-1;
			double[][] Jns = computeOuterEdgeTransitions(node);	
			double[] J = new double[M+1];
			for (int s=0; s<J.length; s++){
				double sum = Double.NEGATIVE_INFINITY;
				for (int n=s; n<Jns.length; n++){
					if (s<Jns[n].length)
						sum = Logarithms.add(sum, Jns[n][s]);
				}
				J[s] = sum;
				assert (J[s]<=0.0);// since it is a probability
			}
			edge_transitions[node] = Jns;
			edge_outside[node] = J;
			
			double[][] Bns = computeOuterNodeTransitions(node);
			double[] B = new double[Bns.length]; 
			for (int n=0; n<Bns.length; n++)
			{
				B[n] = Logarithms.sum(Bns[n], Bns[n].length);
				assert (B[n]<=0.0);// since it is a probability
			}
			node_transitions[node] = Bns;
			node_outside[node] = B; // computeNode(node);		}
		}
		
		/**
		 * Algorithm for computing the outside log-likelihood on an edge. 
		 * @param v node called in preorder 
		 * @return
		 */
		private double[][] computeOuterEdgeTransitions(int v)
		{
			int m = size();
			int M = sizeProfiles.length-1; // 0,1...M
			double[][] Jns;  // return  value
			//double log_pv = getLogLossParameter(v);
			if (tree.isRoot(v))
			{
				double log1_pv = getLogLossComplement(v);
				if (log1_pv==Double.NEGATIVE_INFINITY) // (proot == 1.0)
				{
					Jns = new double[1][M+1];
					Arrays.fill(Jns[0], Double.NEGATIVE_INFINITY);
					if (m==0)
						Jns[0][0]=0.0;
				} else { 
					throw new UnsupportedOperationException("Model with loss at root is  not supported.");
//					// ROOTLOSS
//					Jns = new double[2][];
//					Jns[0] = new double[1];
//					Jns[0][0] = Double.NEGATIVE_INFINITY;
//					Jns[1] = new double[m==0?1:2];
//					Jns[1][0] = log_pv;
//					if (0<m)
//						Jns[1][1] = log1_pv;
				}				
			} else {
				int u = tree.getParent(v);
				int w = tree.getSibling(v);
				final double logit_pv = FamilySizeLikelihood.this.getLogitSurvivalParameter(v, PARAMETER_LOSS);
				final double log1_e = FamilySizeLikelihood.this.getLogLossComplement(w);
				
				final double logit_p = logit_pv + log1_e;
				double log_p = Logarithms.logitToLogValue(logit_p);
				double log1_p = Logarithms.logitToLogComplement(logit_p);
				
				Jns = new double[M+1][];
				for (int n=0; n<Jns.length; n++) {
					double log_nfact = factorials.factln(n);

					Jns[n]=new double[n+1];
					Arrays.fill(Jns[n], Double.NEGATIVE_INFINITY);
					for (int s=0; s<=n; s++) {
						double BuLw = Double.NEGATIVE_INFINITY;
						int t = n-s;
						for (int mw=0; mw<=m; mw++) {
							double[] Bu = sizeProfiles[m-mw].node_outside[u];
							double[][] Lmw = sizeProfiles[mw].getPairingLikelihoods(w);
							double[] Lw = Lmw[n];
							if (n<Bu.length && t<Lw.length)
								BuLw = Logarithms.add(Bu[n]+Lw[t], BuLw);
						}
						double binom = log_nfact - factorials.factln(s)-factorials.factln(t); // n chosse s 
						double slog1_p = s==0?0.0:s*log1_p;
						double tlog_p = t==0?0.0:t*log_p;
						
						Jns[n][s] = BuLw + binom + slog1_p + tlog_p;
					}
				}
			}
			return Jns;
		}
		
		/**
		 * Algorithm for computing the outside log-likelihood at a node. The profile 
		 * must not be ambiguous at the root. 
		 * 
		 * @param u node called in preorder, after setting edge outer likelihood here, and node outer at the parent
		 * @return
		 */
		private double[][] computeOuterNodeTransitions(int u){
			int M = sizeProfiles.length-1; // 0,1...M
			int m = size();
			double[][] Bns; // return value outside transitions
			Bns = new double[M+1][];
			double logq = getLogDuplicationParameter(u); 
			double log1_q = getLogDuplicationComplement(u); 

			double log_gain = getLogGainParameter(u);
			
			double[] J = edge_outside[u]; // already set
			
			if  (logq==Double.NEGATIVE_INFINITY)  { 
				// Poisson
				//double r = factory.getGainParameter(node);
				double logr = log_gain;
				double r = Math.exp(logr);

				for (int n=0; n<=M; n++){
					double[] terms = new double[n+1];
					Arrays.fill(terms,Double.NEGATIVE_INFINITY);
					int t=n,s=0;
					while (s<J.length && s<=n){
						double t_logr = (t==0?0.0:t*logr);
						terms[s] = J[s] - r + t_logr - factorials.factln(t);
						++s;
						--t;
					}
					Bns[n] = terms;
				}				
			} else {
				// Pólya
				double log_kappa = log_gain;
				if (log_kappa==Double.NEGATIVE_INFINITY) {// untested
					double[] terms =new double[1];
					// negative binomial with s and q 
					terms[0] = J[0]; // no conservation
					assert (J[0]<=0.0);
					Bns[0]=terms;
					for (int n=1; n<Bns.length; n++){
						terms = new double[n+1];
						Arrays.fill(terms, Double.NEGATIVE_INFINITY);
						double factln_n = factorials.factln(n-1); // n-1 and not n
						
						int t=n-1, sm1=0, s=1; // no contribution from s=0
						//assert (J.length<=ell+1);
						while (s<J.length && s<=n) // sm1==s-1; t+s=ell
						{
							double binom = factln_n - factorials.factln(sm1) - factorials.factln(t); // s-1 and not s
							terms[sm1] = J[s]+binom + s*log1_q + t*logq;
							assert (terms[sm1]<=0.0);
							
							--t;
							++s;
							++sm1;
						}
						Bns[n] = terms;
					}
				} else {// negative binomial with 0<kappa
					double loglog1_q = Logarithms.logitToLogLogComplement(logq-log1_q);
					for (int n=0; n<Bns.length; n++)
					{
						double[] terms = new double[n+1];
						Arrays.fill(terms, Double.NEGATIVE_INFINITY);
						double log_nfact = gain_factorials[u].factln(n);
						int t=n, s=0; 
						while (s<J.length && s<=n)
						{
							double log_kappa_s = Logarithms.add(log_kappa, Math.log(s));
							double log_kappa_slog1_q = log_kappa_s + loglog1_q;
							double kappa_s_log1_q = -Math.exp(log_kappa_slog1_q); // we want -(-ln(1-q))*(kappa+s)
							double t_logq= t==0?0.0:t*logq;
							double binom = log_nfact - gain_factorials[u].factln(s) - factorials.factln(t);
							terms[s] = J[s] + binom + kappa_s_log1_q + t_logq;
							assert (terms[s]<=0.0); // since it is a probability
							++s;
							--t;
						}
						Bns[n] = terms;
					} // for n
				} // positive kappa				
			} // polya 
			return Bns;
		}		
		
		
	}
	
	/**
	 * A "profile" for families with given maximum total copies; 
	 * with getters adapted to size-dependent distributions 
	 * 
	 * @param rates
	 * @param maxUnobserved
	 * @return
	 */
	public static Posteriors.Profile getUnobservedProfile(TreeWithRates rates, int maxUnobserved){
		TreeWithLogisticParameters lrates = (rates instanceof TreeWithLogisticParameters)
				?(TreeWithLogisticParameters)rates
				:new TreeWithLogisticParameters(rates,false);
		FamilySizeLikelihood zis = new FamilySizeLikelihood(lrates, maxUnobserved);
		return zis.getUnobservedProfile();
	}
	
	
	public Posteriors.Profile getUnobservedProfile(){
		MockPosteriors MP = new MockPosteriors();
		return MP.new UnobservedProfile();
	}
	
	/**
	 * Used for enclosing a Posteriors.Profile object that serves up the posteriors
	 */
	private class MockPosteriors extends Posteriors {
		
		MockPosteriors(){
			super(FamilySizeLikelihood.this);
		}
		
		/**
		 * Posteriors.Profile interface to posterior 
		 * statistics over an unobserved profile.
		 * 
		 */
		class UnobservedProfile extends Posteriors.Profile{
			/**
			 * Instantiation triggers the computing of all 
			 * inner and outer likelihoods 
			 */
			UnobservedProfile(){
				super(0); // immaterial family index 
				this.computeLikelihoods();
			}
			
			/**
			 * Calls {@link FamilySizeLikelihood#computeAll()}
			 */
			@Override 
			public void computeLikelihoods() {
				FamilySizeLikelihood.this.computeAll();
			}
			
			/**
			 * Posterior probabilities for loss changes
			 * on logarithmic scale : post[<var>n</var>][<var>s</var>] is the 
			 * logarithm of the posterior probability
			 * for <var>n</var>&rarr;<var>s</var> transition  
			 * to have <var>s</var> surviving copies at node
			 * with <var>n</var>  at parent.
			 * 
			 * @param v node
			 * @return
			 */
			@Override
			public double[][] getLogEdgeTransitionPosteriors(int v){
				int M = sizeProfiles.length-1;
				double LL = Double.NEGATIVE_INFINITY;
				for (int s=0; s<=M; s++) {
					for (int m=s; m<=M; m++) {
						double [] J = sizeProfiles[M-m].edge_outside[v];
						double [] K = sizeProfiles[m].getEdgeLikelihoods(v);
						if (s<K.length) {
							double t = J[s]+K[s];
							LL = Logarithms.add(LL, t);						
						}
					}
				}
				//System.out.println("#**FSL.gLETP "+v+"\tLL "+LL);
				double[][] lpost = new double[M+1][]; // return value
				for (int n=0; n<=M; n++) {
					lpost[n] = new double[n+1];
					for (int s=0; s<=n; s++) {
						double t = Double.NEGATIVE_INFINITY;
						for (int m=s; m<=M; m++) {
							double[][] Jns = sizeProfiles[M-m].edge_transitions[v];
							double[] K = sizeProfiles[m].getEdgeLikelihoods(v);
							if (s<K.length && n<Jns.length)
								t = Logarithms.add(t,  Jns[n][s]+K[s]);
						}
						lpost[n][s] = t-LL;
					}
				}
				return lpost;
			}
			
			@Override
			public double getLogLikelihood() {
				int root = tree.getRoot();
				int M = sizeProfiles.length-1;
				double LL = Double.NEGATIVE_INFINITY;
				for (int m=0; m<=M; m++) {
					double [] J = sizeProfiles[M-m].edge_outside[root];
					double [] K = sizeProfiles[m].getEdgeLikelihoods(root);
					for (int s=0; s<K.length; s++) {
						double t = J[s]+K[s];
						LL = Logarithms.add(LL, t);						
					}
				}
				return LL;
			}
			
			/**
			 * Posterior probabilities for gain+duplication changes
			 * on logarithmic scale : post[<var>n</var>][<var>s</var>] is the 
			 * logarithm of the posterior probability
			 * for <var>s</var>&rarr;<var>n</var> transition  
			 * to have <var>n</var> surviving copies at node
			 * with <var>s</var> inherited among them. 
			 * 
			 * @param v
			 * @return
			 */
			@Override
			public double[][] getLogNodeTransitionPosteriors(int v){
				int M = sizeProfiles.length-1;
				double LL = Double.NEGATIVE_INFINITY;
				for (int n=0; n<=M; n++) {
					for (int m=n; m<=M;m++) {
						double[] B = sizeProfiles[M-m].node_outside[v];
						double[] C = sizeProfiles[m].getNodeLikelihoods(v);
						double t = B[n]+C[n];
						LL = Logarithms.add(LL, t);
					}
				}
				//System.out.println("#**FSL.gLNTP "+v+"\tLL "+LL);
				double[][] lpost = new double[M+1][];  // return value
				for (int n=0; n<=M; n++) {
					lpost[n] = new double[n+1];
					for (int s=0; s<=n; s++) {
						double t = Double.NEGATIVE_INFINITY;
						for (int m=n; m<=M; m++) {
							double[][] Bns = sizeProfiles[M-m].node_transitions[v];
							double[] C = sizeProfiles[m].getNodeLikelihoods(v);
							if (s<Bns[n].length)
								t = Logarithms.add(Bns[n][s]+C[n], t);
						}
						lpost[n][s] = t-LL;
					}
				}
				return lpost;
			}
			@Override
			public double[] getLogEdgePosteriors(int v){
				int M = sizeProfiles.length-1;
				double LL = Double.NEGATIVE_INFINITY;
				double[] lpost = new double[M+1]; // return value
				for (int s=0; s<=M; s++) {
					double z = Double.NEGATIVE_INFINITY;
					for (int m=s; m<=M; m++) {
						double [] J = sizeProfiles[M-m].edge_outside[v];
						double [] K = sizeProfiles[m].getEdgeLikelihoods(v);
						
						if (s<K.length) {
							double t = J[s]+K[s];
							z = Logarithms.add(z, t);
						}
					}
					lpost[s] = z;
					LL = Logarithms.add(LL, z);						
				}
				//System.out.println("#**FSL.gLEP "+v+"\tLL "+LL);
				
				for (int s=0; s<=M; s++) {
					lpost[s]-=LL;
				}
				return lpost;
			}
			@Override
			public double[] getLogNodePosteriors(int v){
				int M = sizeProfiles.length-1;
				double LL = Double.NEGATIVE_INFINITY;
				double[] lpost = new double[M+1];
				for (int n=0; n<=M; n++) {
					double z= Double.NEGATIVE_INFINITY;
					for (int m=n; m<=M;m++) {
						double[] B = sizeProfiles[M-m].node_outside[v];
						double[] C = sizeProfiles[m].getNodeLikelihoods(v);
						double t = B[n]+C[n];
						z = Logarithms.add(z, t);
					}
					lpost[n] = z;
					LL = Logarithms.add(LL, z);
				}
				//System.out.println("#**FSL.gLNP "+v+"\tLL "+LL);
				for (int n=0; n<=M; n++) {
					lpost[n] -= LL;
				}
				return lpost;
			}
			

			@Override 
			public double[] getEdgePosteriors(int v) {
				double[] lpost = this.getLogEdgePosteriors(v);
				return exp(lpost);
			}
			
			@Override
			public double[] getNodePosteriors(int v) {
				double[] lpost = this.getLogNodePosteriors(v);
				return exp(lpost);
			}
			/*
			 * Unsupported operations
			 */
			
			/**
			 * Should not be called because {@link #computeLikelihoods()} is redefined
			 * @throws UnsupportedOperationException always
			 */
			@Override
			protected void computeOutside(int v) 
			{throw new UnsupportedOperationException();}
			
			/**
			 * @throws UnsupportedOperationException always
			 */
			@Override
			protected void computeOutsideFromRoot(int v)
			{throw new UnsupportedOperationException();}
			
			
			/**
			 * @throws UnsupportedOperationException always
			 */
			@Override
			protected double[] getNodeOutside(int v)
			{throw new UnsupportedOperationException();}
			
			/**
			 * @throws UnsupportedOperationException always
			 */
			@Override
			protected double[] getEdgeOutside(int node)
			{throw new UnsupportedOperationException();}
			
			
			/**
			 * @throws UnsupportedOperationException always
			 */
			@Override
			public double[] getNodeAncestorPosteriors(int node) 
			{throw new UnsupportedOperationException();}
			/**
			 * @throws UnsupportedOperationException always
			 */
			@Override
			public double[] getFamilyEventPosteriors(int node)
			{throw new UnsupportedOperationException();}

			/**
			 * @throws UnsupportedOperationException always
			 */
			@Override
			public double[] getEdgeAncestorPosteriors(int node)
			{throw new UnsupportedOperationException();}
			
		} // UnobservedProfile
	} // MockPosteriors
		
	/**
	 * Point-mass-function on logarithmic scale 
	 * 
	 * @return
	 */
	public double[] getLogPMF() {
		double[] getLogPMF = new double[sizeProfiles.length];
		for (int m=0; m<sizeProfiles.length; m++)
			getLogPMF[m] = sizeProfiles[m].getLogLikelihood();
		return getLogPMF;
	}
	
	private void reportStats(PrintStream out, ProfileTable table) {
		LogGradient LG = new LogGradient((TreeWithLogisticParameters)this.rates, table);
		LG.setMinimumObservedCopies(1);
		LogGradient.PosteriorStatistics stat0 = LG.getUnobservedStatistics();
		LG.setMinimumObservedCopies(2);
		LogGradient.PosteriorStatistics stat1 = LG.getUnobservedStatistics();

		LogGradient.PosteriorStatistics statunobs = LG.getProfileStatistics(getUnobservedProfile());
		
		int[] ncopies = new int[sizeProfiles.length+1]; // last cell counts max< profiles
		int nF = table.getFamilyCount();
		for (int f=0; f<nF; f++) {
			int m = table.getMemberCount(f);
			ncopies[Integer.min(ncopies.length-1, m)]++;
		}
		
		//this.computeAll();
		double[] log_pm = getLogPMF();
		
		double cdf=0.0;
		double tail=1.0;
		int mincopy = table.minCopies();
		double p0 = 0.0;
		double p_not0 = 1.0;
		for (int m=0; m<mincopy; m++) {
			double lpm = m<log_pm.length?log_pm[m]:Double.NEGATIVE_INFINITY;
			double pm = Math.exp(lpm);
			p0 += pm;
			p_not0 -= pm;
		}
		
		boolean observed = false;
		
		out.println("#SIZE: distribution of number of copies (observed: min "+mincopy+")");
		out.println("#SIZE\tm\tlogp\tp\tcdf\ttail\tobs(n)\texp(n)\test(cdf)");
		int sum_ncopies = 0;
		int nFtail = nF;
		for (int m=0; m<ncopies.length; m++) {
			double lpm = m<log_pm.length?log_pm[m]:Double.NEGATIVE_INFINITY;
			double pm = Math.exp(lpm);
			cdf += pm;
			int n = ncopies[m];
			
			double est_cdf;
			if (m==0)
				est_cdf = Math.exp(stat0.getLogLikelihood());
			else if (m==1)
				est_cdf = Math.exp(stat1.getLogLikelihood());
			else 
				est_cdf = 0.0;
//			out.println(m+"\t"+log_pm[m]+"\t"+pm+"\t"+cdf+"\t"+tail
//					+"\t"+n
//					+"\t"+(nF*pm)
//					+"\t"+est_cdf);
			
			double exp_n;
			if (m<mincopy) {
				exp_n = nF*pm / p_not0;
				assert (n==0);
			} else {
				// nFtail are the families with at least m copies 
				// tail is the probability of at least m copies 

				exp_n = nFtail * pm/tail;
			}
			
//			double exp_n;
//			if (log_pm.length <= m)
//				exp_n = tail*nF;
//			else 
//				exp_n = nF*pm;
//			if (observed || n!=0) {
//				exp_n /= p_not0;
//				observed = true;
//			} else {
//				p0 += pm;
//				p_not0 -= pm;
//			}
			
			out.printf("#SIZE\t%d\t%g\t%g\t%g\t%g\t%d\t%g\t%g\t// nFtail = %d pm/tail=%g\n", 
					m, lpm, pm, cdf, tail, n, exp_n, est_cdf, nFtail, pm/tail);
			nFtail -= n;
			sum_ncopies += n;
			tail -= pm;
		}
		LogGradient.PosteriorStatistics statexp = (maxCopies()==0?stat0:stat1);
		out.printf("#POST: posteriors; LL %f exp %f\n", statunobs.getLogLikelihood(), statexp.getLogLikelihood());

		for (int u=0; u<tree.getNumNodes(); u++) {
			double[] sUnobs = statunobs.getLogEdgePosteriors(u);
			double[] sExp = statexp.getLogEdgePosteriors(u);
			
			out.printf("#POST\t%d\tS", u);
			for (int i=0; i<sUnobs.length || i<sExp.length; i++) {
				double su = i<sUnobs.length?sUnobs[i]:Double.NEGATIVE_INFINITY;
				double se = i<sExp.length?sExp[i]:Double.NEGATIVE_INFINITY;
				out.printf("\t%f/%f", su,se);
			}
			out.println();
		}			
		for (int u=0; u<tree.getNumNodes(); u++) {
			double[] sUnobs = statunobs.getLogBirthTails(u);
			double[] sExp = statexp.getLogBirthTails(u);
			out.printf("#POST\t%d\tbirth", u);
			for (int i=0; i<sUnobs.length || i<sExp.length; i++) {
				double su = i<sUnobs.length?sUnobs[i]:Double.NEGATIVE_INFINITY;
				double se = i<sExp.length?sExp[i]:Double.NEGATIVE_INFINITY;
				out.printf("\t%f/%f", su,se);
			}
			out.println();
		}
		for (int u=0; u<tree.getNumNodes(); u++) {
			double[] sUnobs = statunobs.getLogDeathTails(u);
			double[] sExp = statexp.getLogDeathTails(u);
			out.printf("#POST\t%d\tdeath", u);
			for (int i=0; i<sUnobs.length || i<sExp.length; i++) {
				double su = i<sUnobs.length?sUnobs[i]:Double.NEGATIVE_INFINITY;
				double se = i<sExp.length?sExp[i]:Double.NEGATIVE_INFINITY;
				out.printf("\t%f/%f", su,se);
			}
			out.println();
		}
		
	}
	
	public static void main(String[] args) throws Exception
	{
		Class<?> our_class = java.lang.invoke.MethodHandles.lookup().lookupClass();
		count.io.CommandLine cli = new count.io.CommandLine(args, our_class);

		PrintStream out = System.out;
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(our_class));
    	    out.println(CommandLine.getStandardRuntimeInfo(our_class, args));
    	}
    	
    	RateVariationModel model = cli.getVariationModel();
		TreeWithRates rates = model.getBaseModel();    
		TreeWithLogisticParameters lrates = rates instanceof TreeWithLogisticParameters?
					(TreeWithLogisticParameters)rates:new TreeWithLogisticParameters(rates);
		
    	Phylogeny phylo = cli.getTree();
		AnnotatedTable table = cli.getTable();
		int max_input = table==null?0:table.getMaxFamilySizes(phylo)[phylo.getRoot()];
		int max_copy = cli.getOptionInt(CommandLine.OPT_MAXCOPY, max_input);		
		
		out.println(CommandLine.getStandardHeader("Maximum copies: "+CommandLine.OPT_MAXCOPY+" "+max_copy));
		
		FamilySizeLikelihood D = new FamilySizeLikelihood(lrates, max_copy);
		D.reportStats(out, table);
		
		if (out_file!=null)
			out.close();
		
		
	}	
	
	
	
}
