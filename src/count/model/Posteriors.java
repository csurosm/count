/*
 * Copyright 2021 Mikl&oacute;s Cs&#369;r&ouml;s.
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

import count.ds.ProfileTable;
import count.ds.TreeTraversal;
import count.io.CommandLine;
import count.matek.Logarithms;

import java.io.PrintStream;
import java.util.Arrays;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;


/**
 * Computing prior probabilities top-down, and expectations.
 * 
 *
 *
 */
public class Posteriors
{
	public Posteriors(TreeWithRates rates, ProfileTable table)
	{
		//this (new Likelihood(rates,table));
		this(new LikelihoodParametrized(rates, table)); // LikelihoodParametrized works on an extreme range of parameter values 
		//System.out.println("#**P() with "+factory.getClass());
	}
	
	/**
	 * Dummy instantiation not for computing with.
	 */
	private Posteriors() {this(null);}
	
	public Posteriors(Likelihood factory) 
	{ 
		this.factory = factory;
	}
	
	protected Profile getPosteriors(int family_idx)
	{
		
		Profile P = new Profile(family_idx);
//		for (int node=factory.tree.getNumNodes(); node>0; )
//		{
//			--node;
//			System.out.println("#*P.getP family "+family_idx+"\tedge "+node+"\tJ[] = "
//					+Arrays.toString(P.getEdgeOutside(node)));
//			System.out.println("#*P.getP family "+family_idx+"\tnode "+node+"\tB[] = "
//					+Arrays.toString(P.getNodeOutside(node)));
//		}
    	
		return P;
	}
	
	protected final Likelihood factory;	
	protected Posteriors[] ancestor_posteriors;
	
	private double ancestor_deviation = 0.0; // set to 0 for don't care
	
	public void initAncestorPosteriors(int ancestor_width)
	{
		int num_nodes = factory.tree.getNumNodes();
		
		ancestor_posteriors = new Posteriors[num_nodes];
		for (int node=factory.tree.getNumLeaves(); node<num_nodes; node++) // only for inner nodes 
		{
			Ancestor anc = new Ancestor(factory, node);
			anc.setCalculationWidth(ancestor_width);
			Posteriors P = new Posteriors(anc);
			ancestor_posteriors[node] = P;
		}
	}
	
	public void setCalculationWidthThresholds(int absolute, double relative)
	{
		factory.setCalculationWidthThresholds(absolute, relative);
	}
	
	public int getCalculationWidthAbsolute()
	{
		return factory.getCalculationWidthAbsolute(); 
	}
	
	public double getCalculationWidthRelative()
	{
		return factory.getCalculationWidthRelative();
	}
	
	public void setAncestorWidthThreshold(double deviation)
	{
		this.ancestor_deviation = deviation;
	}
	
	
//	private Profile getAncestorPosteriors(int node, int family_idx)
//	{
//		if (ancestor_posteriors == null)
//			initAncestorPosteriors();
//		return ancestor_posteriors[node].getPosteriors(family_idx);
//	}
	
	
	public enum FamilyEvent { GAIN, LOSS, EXPAND, CONTRACT};

	static double[] computeCumulative(double[] Pr)
	{
		double prev = Pr[0];
		for (int i=1; i<Pr.length; i++)
		{
			prev = Pr[i] = Pr[i]+prev;
		}
		return Pr;
	}
	
	static double[] computeTail(double[] p)
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
	
	static double computeExpectation(double[] Pr)
	{
		double m = 0.0;
		if (Pr==null) return m;
		int j = Pr.length;
		while (1<j)
		{
			--j;
			m = Math.fma(j, Pr[j], m); // Math.fma(a,b,c) = a*b+c
		}
		return m;
//		
//		double expect = 0.0;
//		for (int ell=1; ell<Pr.length; ell++)
//			expect += ell*Pr[ell];
//		return expect;
	}
	
	protected static double[] computeLogPosteriors(double[] outside, double[] inside)
	{
		double[] posteriors = outside.clone();
		if (inside.length!=0)
		{
			for (int ell=0; ell<posteriors.length; ell++)
			{
				posteriors[ell]+=inside[ell];
			}
		}
		double LL = Logarithms.sum(posteriors, posteriors.length);
		for (int ell=0; ell<posteriors.length; ell++)
		{
			posteriors[ell] = posteriors[ell] - LL; 
		}
		return posteriors;		
	}
	
	protected static double computeLogMean(double[] log_p)
	{
		// compute tails
		double log_tail = Double.NEGATIVE_INFINITY;
		double log_mean = Double.NEGATIVE_INFINITY;
		int j=log_p.length;
		while (0<j)
		{
			--j;
			log_mean = Logarithms.add(log_mean,  log_tail);
			double x = log_p[j];
			log_tail = Logarithms.add(log_tail, x);
		}
		return log_mean;
	}
	
	
	
	/**
	 * 
	 * @param outside outside log-likelihoods
	 * @param inside inside log-likelihoods
	 * @return posterior probabilities normalized by sum of cell-by-cell products
	 */
	protected static double[] computePosteriors(double[] outside, double[] inside)
	{
		double[] posteriors = outside.clone();
		if (inside.length!=0)
		{
			// not ambiguous
			for (int ell=0; ell<posteriors.length; ell++)
			{
//				if (Double.isNaN(posteriors[ell]) || Double.isNaN(inside[ell]))
//				{
//					System.out.println("#**P.P.cP "+ell+"\tpost "+posteriors[ell]+"\tins "+inside[ell]);
//			
//				}
				posteriors[ell]+=inside[ell];
			}
		}
		double LL = Logarithms.sum(posteriors, posteriors.length);
		for (int ell=0; ell<posteriors.length; ell++)
		{
			posteriors[ell]=Math.exp(posteriors[ell]-LL);
//			if (Double.isNaN(posteriors[ell]))
//			{
//				System.out.println("#**P.P.cP "+ell+"\tpost "+posteriors[ell]+"\tLL "+LL);
//			}
//			assert (!Double.isNaN(posteriors[ell]));
		}
		return posteriors;
	}
	
	protected static double computeLogLikelihood(double[] outside, double[] inside)
	{
		//double[] posteriors = outside.clone();
		double LL = Double.NEGATIVE_INFINITY;
		if (inside.length!=0)
		{
			// not ambiguous
			for (int ell=0; ell<outside.length; ell++)
			{
				LL = Logarithms.add(LL, outside[ell]+inside[ell]);
			}
		}
		return LL;
	}
	

	/**
	 * If set to true, then it is verified if the 
	 * Ancestor calculations in {@link Profile#getNodeAncestorPosteriors(int)} 
	 * give believable 
	 * posteriors for ancestor copies: if greater than 1, then 
	 * exception is thrown and debug messages are written about the calculations. 
	 * Otherwise, greater than 1 is dismissed as a numerical error, and 
	 * replaced by 1.0. 
	 * (For very large trees, a small 1e-8 relative error may occur: cannot attribute it to bad code).   
	 */
	private static final boolean DEBUG_EVENT_POSTERIORS = true;

	public class Profile 
	{
		private Profile(int family_idx)
		{
			this.inside = factory.getProfileLikelihood(family_idx);
			int num_nodes = factory.tree.getNumNodes();
			node_outside = new double[num_nodes][];
			node_transitions = new double[num_nodes][][];
			edge_outside = new double[num_nodes][];
			edge_transitions = new double[num_nodes][][];
		}
		private final double[][] node_outside;
		private final double[][] edge_outside;
		
		private final double[][][] node_transitions;
		private final double[][][] edge_transitions;
		
		final Likelihood.Profile inside;
		
		private Profile[] ancestor=null;
		
		public void computeLikelihoods()
		{
			inside.computeLikelihoods();
			int root = factory.tree.getRoot();
			for (int node=root; node>=0; --node)
			{
				computeOutside(node);
//				edge_outside[node] = computeEdge(node);
//				node_outside[node] = computeNode(node);
			}
		}
		
		protected void computeOutside(int  node)
		{
			double[][] Jls = computeEdgeTransitions(node);
			double[] J = new double[inside.getEdgeLikelihoods(node).length];
			for (int s=0; s<J.length; s++)
			{
				double sum = Double.NEGATIVE_INFINITY;
				for (int ell=s; ell<Jls.length; ell++)
				{
					if (s<Jls[ell].length)
						sum = Logarithms.add(sum, Jls[ell][s]);
					
//					if (!(sum<=0.0))
//					{
//						System.out.println("#***P.cO node "+node+"\tJ[s="+s+"] "+sum+"\tell "+ell+"\tJls "+Jls[ell][s]);
//					}
					
				}
				J[s] = sum;
//				if (!(J[s]<=0.0))
//				{
//					System.out.println("#***P.cO node "+node+"\tJ[s="+s+"] "+J[s]+"\t// "+factory.rates.toString(node)+"\t// "+factory.toString(node));
//				}
				assert (J[s]<=0.0);// since it is a probability
			}
			edge_transitions[node] = Jls;
			edge_outside[node] = J;
//			edge_outside[node] = computeEdge(node);
			
			double[][] Bls = computeNodeTransitions(node);
			double[] B = new double[Bls.length]; 
			for (int ell=0; ell<Bls.length; ell++)
			{
				B[ell] = Logarithms.sum(Bls[ell], Bls[ell].length);
				
//				if (!(B[ell]<=0.0))
//				{
//					System.out.println("#***P.cO node "+node+"\tB[ell="+ell+"] "+B[ell]+"\tBls "+Arrays.toString(Bls[ell]));
//				}
				assert (B[ell]<=0.0);// since it is a probability
			}
			node_transitions[node] = Bls;
			node_outside[node] = B; // computeNode(node);
		}
				
		protected void computeAncestors()
		{
			synchronized(Posteriors.this) //  
			{
				if (ancestor_posteriors == null)
					initAncestorPosteriors(2);
			}

			ancestor = new Profile[factory.tree.getNumNodes()];
			for (int a=factory.tree.getNumLeaves(); a<ancestor.length; a++)
			{
				
				Profile A 
				= ancestor[a] 
				= ancestor_posteriors[a].getPosteriors(inside.family_idx);
				
				
				if (ancestor_deviation != 0.0)
				{
					int m = factory.table.maxCopies(inside.getFamilyIdx());
					int delta =(int) Math.ceil(ancestor_deviation*Math.sqrt(m+1.0));
					Ancestor anc = (Ancestor) A.inside.getOwner();
					anc.setCalculationWidth(m+delta+1); 
				}
				
				
				// will be consulted only outside of a's subtree and path to root,
				// but we calculate it everywhere (filtering out a's subtree 
				// is only a few nodes most of the time).
				
				A.inside.computeLikelihoods(); 
				
				// outside likelihoods are needed only between the root and this node a
				A.computeOutsideFromRoot(a);
//				int[] path = TreeTraversal.getPathToRoot(factory.tree, a);
//				for (int i=path.length-1; i>=0; i--)
//				{
//					int node=path[i];
//					A.edge_outside[node] = A.computeEdge(node);
//					A.node_outside[node] = A.computeNode(node);
//				}
			}
			
			
		}

		/**
		 * Computes outside likelihoods along the path from the root to a node
		 * (noth included). 
		 */
		protected void computeOutsideFromRoot(int node)
		{
			int[] path = TreeTraversal.getPathToRoot(factory.tree, node);
			for (int i=path.length-1; i>=0; i--)
			{
				int a=path[i];
//				edge_outside[a] = computeEdge(a);
//				node_outside[a] = computeNode(a);
				computeOutside(a);
			}
		}
		
		
		public double[] getNodeAncestorPosteriors(int node) 
		{
			double[] B = getNodeAncestor(node);
			double[] C = inside.getNodeLikelihoods(node);
			double epsi = factory.getExtinction(node);

			double[] getNodeAncestorPosteriors = computeAncestorPosteriors(B, C, epsi);
			
			if (getNodeAncestorPosteriors==null) // numerical problem: DEBUG code
			{
				IndexedTree tree = inside.getOwner().tree;
				int root = tree.getRoot();
				double[] K = inside.getEdgeLikelihoods(root);
				
				double[] Bsurv = getNodeOutside(node);
				double LL = Double.NEGATIVE_INFINITY;	
				for (int ell=0; ell<Bsurv.length; ell++)
				{
					LL = Logarithms.add(LL,  Bsurv[ell]+C[ell]);
				}
				System.out.println("#***P.P.gNAP fam "+inside.family_idx+"\tnode "+node
						+"\tB "+Arrays.toString(B)
						+"\tC "+Arrays.toString(C)
						+"\tBsurv "+Arrays.toString(Bsurv)
						+"\tpost "+Arrays.toString(computePosteriors(Bsurv, C))
						+"\tLL "+LL
						+"\tKroot "+Arrays.toString(K)
						+"\tLLempty "+inside.getOwner().getEmptyLL()
						+"\t// "+inside.toString());
				
				for (int v=0; v<=root; v++)
				{
					Bsurv = getNodeOutside(v);
					C = inside.getNodeLikelihoods(v);
					double LLn = computeLogLikelihood(Bsurv,C);
					double[] Jsurv = getEdgeOutside(v);
					K = inside.getEdgeLikelihoods(v); 
					double LLe = computeLogLikelihood(Jsurv,K);
					System.out.println("#***P.P.gNAP fam "+inside.family_idx+"\tnode "+v
							+"\tB "+Arrays.toString(Bsurv)
							+"\tC "+Arrays.toString(C)	
							+"\tJ "+Arrays.toString(Jsurv)
							+"\tK "+Arrays.toString(K)
							+"\tLLn "+LLn
							+"\tLLe "+LLe
							+"\tkappa "+factory.rates.getGainParameter(v)
							+"\tlog1_q "+factory.rates.getLogDuplicationComplement(v)
							+"\tmul "+(factory.rates.getGainParameter(v)*factory.rates.getLogDuplicationComplement(v))
							+"\t// "+inside.getOwner().rates.toString(v));
				}
				
				
				
			}
			assert getNodeAncestorPosteriors!=null;
			return getNodeAncestorPosteriors;
		}
		
		public double[] getEdgeAncestorPosteriors(int node)
		{
			double[] J = getEdgeAncestor(node);
			double[] K = inside.getEdgeLikelihoods(node);
			double epsi = factory.getExtinction(node) * factory.getDuplicationParameterComplement(node);
					// *(1.0-factory.getDuplicationParameter(node)); 
			return computeAncestorPosteriors(J, K, epsi);
		}
		
		/**
		 * Posterior distribution of ancestral copy counts. 
		 * 
		 * @param node
		 * @return
		 */
		public double[] getNodePosteriors(int node)
		{
			double[] B = getNodeOutside(node); // sets up inside likelihoods if necessary
			double[] C = inside.getNodeLikelihoods(node);
//			System.out.println("#**P.P.gNP "+node+"\tB "+Arrays.toString(B)+"\tC "+Arrays.toString(C));
			
			return computePosteriors(B, C);
		}
		
		/**
		 * Posterior probabilities for gain+duplication changes
		 * on logarithmic scale : post[<var>n</var>][<var>s</var>] is the 
		 * logarithm of the posterior probability
		 * for <var>s</var>&rarr;<var>n</var> transition  
		 * to have <var>n</var> surviving copies at node
		 * with <var>s</var> inherited among them. 
		 * 
		 * @param node
		 * @return
		 */
		public double[][] getLogNodeTransitionPosteriors(int node)
		{
			double[] B = getNodeOutside(node); // sets up inside likelihoods if necessary
			double[] C = inside.getNodeLikelihoods(node);
			
			assert (C.length != 0);
			
			double LL = Double.NEGATIVE_INFINITY;
			for (int ell=0; ell<B.length; ell++)
			{
				double t = B[ell]+C[ell];
				LL = Logarithms.add(LL, t);
			}
			double[][] post = new double[B.length][];
			double[][] Bls = node_transitions[node];
			for (int ell=0; ell<B.length; ell++)
			{
				post[ell] = new double[ell+1];
				for (int s=0; s<=ell; s++)
					post[ell][s] = Bls[ell][s] + C[ell] - LL;
			}			
			return post;
		}
		
		/**
		 * Posterior probabilities for loss changes
		 * on logarithmic scale : post[<var>n</var>][<var>s</var>] is the 
		 * logarithm of the posterior probability
		 * for <var>n</var>&rarr;<var>s</var> transition  
		 * to have <var>s</var> surviving copies at node
		 * with <var>n</var>  at parent.
		 * 
		 * @param node
		 * @return
		 */
		public double[][] getLogEdgeTransitionPosteriors(int node)
		{
			double[] J = getEdgeOutside(node); // sets up inside likelihoods if necessary
			double[] K = inside.getEdgeLikelihoods(node);
			double LL = Double.NEGATIVE_INFINITY;
			for (int s=0; s<J.length; s++)
			{
				double t = J[s]+K[s];
				LL = Logarithms.add(LL, t);
			}
			double[][] Jls = edge_transitions[node];
			double[][] post = new double[Jls.length][];
			for (int ell=0; ell<post.length; ell++)
			{
				post[ell] = new double[ell+1];
				Arrays.fill(post[ell], Double.NEGATIVE_INFINITY);
//				if (Jls[ell].length < Math.min(K.length, ell+1)) // debug
//				{
//					System.out.println("#***P.gLETP "+node+"\tell "+ell+"\tJls.len "+Jls[ell].length+"\tKlen "+K.length
//							+"\t// Jls "+Arrays.toString(Jls[ell])+"\tK "+Arrays.toString(K)
//							+"\t// "+factory.toString(node)
//							+"\t// "+factory.rates.toString(node));
//				}
				assert Math.min(K.length, ell+1)<=Jls[ell].length;
				
				for (int s=0; s<=ell && s<K.length; s++)
				{
					post[ell][s] = Jls[ell][s]+K[s]-LL;
				}
			}
			return post;
		}
		
		public double getLogEdgeMean(int node)
		{
			double[] J = getEdgeOutside(node); // sets up inside likelihoods if necessary
			double[] K = inside.getEdgeLikelihoods(node);
			double[] log_p  = computeLogPosteriors(J, K);
			// compute tails
			double log_tail = Double.NEGATIVE_INFINITY;
			int j=log_p.length;
			while (0<j)
			{
				--j;
				double x = log_p[j];
				log_p[j] = log_tail;
				log_tail = Logarithms.add(log_tail, x);
			}
			double log_mean = Logarithms.sum(log_p, log_p.length);
			return log_mean;
		}
		
		public double[] getLogEdgePosteriors(int node)
		{
			double[] J = getEdgeOutside(node); // sets up inside likelihoods if necessary
			double[] K = inside.getEdgeLikelihoods(node);
			double[] log_p  = computeLogPosteriors(J, K);
			return log_p;
		}
		
		public double[] getLogNodePosteriors(int node)
		{
			double[] B = getNodeOutside(node); // sets up inside likelihoods if necessary
			double[] C = inside.getNodeLikelihoods(node);
			double[] log_p = computeLogPosteriors(B, C);
			
			return log_p;
		}
		
		public double getLogNodeMean(int node)
		{
			double[] log_p = getLogNodePosteriors(node);
			double log_mean = computeLogMean(log_p);
			return log_mean;
		}
		
		/**
		 * Array of log-tail-difference for copy births. 
		 * Summing the probabilities in this array gives the 
		 * posterior expectation log EXP(N-S) for N=node copies, S=edge copies.
		 * 
		 * @param node
		 * @return [log (Prob(node copies&gt;i)-Prob(edge copies&gt;i))]
		 */
		public double[] getLogNodePosteriorIncrease(int  node)
		{
			double[][] Pls = getLogNodeTransitionPosteriors(node);
			double[] N_S = logTailDifference(Pls);
			return N_S;
		}
		
		public double[] getNodeBirthTails(int node)
		{
			double[] N_S = this.getLogNodePosteriorIncrease(node);
			for (int j=0;j<N_S.length; j++)
				N_S[j]=Math.exp(N_S[j]);
			return N_S;
		}
		
		public double getLogNodeIncrease(int node)
		{
			double[][] Pls = getLogNodeTransitionPosteriors(node);
			double[] N_S = logTailDifference(Pls);
			return Logarithms.sum(N_S, N_S.length);
		}
		
		
		
		public double getLogEdgeDecrease(int node)
		{
			double[][] Pls = getLogEdgeTransitionPosteriors(node);
			double[] N_S = logTailDifference(Pls);
			return Logarithms.sum(N_S, N_S.length);
		}
		
		/**
		 * Array of log-tail-difference for copy deaths. 
		 * Summing the probabilities in this array gives the 
		 * posterior expectation log EXP(N-S) for N=parent copies, S=edge copies.
		 * 
		 * @param node
		 * @return [log (Prob(parent copies&gt;i)-Prob(edge copies&gt;i))]
		 */
		public double[] getLogEdgePosteriorDecrease(int node)
		{
			double[][] Pls = getLogEdgeTransitionPosteriors(node);
			double[] N_S = logTailDifference(Pls);
			return N_S;
		}

		
		public double[] getEdgeDeathTails(int node)
		{
			double[] N_S = this.getLogEdgePosteriorDecrease(node);
			for (int j=0;j<N_S.length; j++)
				N_S[j]=Math.exp(N_S[j]);
			return N_S;
		}
		
		
		/**
		 * N_S[k] = logsum (s=0..k; ell=k+1..max) Tls 
		 * 
		 * @param Tls
		 * @return
		 */
		private double[] logTailDifference(double[][] Tls)
		{
			double[] N_S = new double[Tls.length];
			Arrays.fill(N_S, Double.NEGATIVE_INFINITY);
			for (int s=0; s<N_S.length; s++)
			{
				double tail = Double.NEGATIVE_INFINITY;
				for (int ell=Tls.length-1; ell>=s; ell--)
				{
					// tails sums ell=s+1..max
					N_S[ell] = Logarithms.add(N_S[ell], tail); 
					
//					if (Tls[ell].length<=s) // DEBUG
//					{
//						System.out.println("#***P.lTD ell "+ell+"\ts "+s+"\tTell "+Arrays.toString(Tls[ell]));
//					}
					assert (s<Tls[ell].length);
					tail = Logarithms.add(tail, Tls[ell][s]);
				}
			}
			return N_S;
		}
		
		
		public double[] getEdgePosteriors(int node)
		{
			double[] J = getEdgeOutside(node); // sets up inside likelihoods if necessary
			double[] K = inside.getEdgeLikelihoods(node);
			return computePosteriors(J, K);
		}
		
		public double getNodeMean(int node)
		{
			double[] Pr = getNodePosteriors(node);
			return computeExpectation(Pr);
		}

		public double[] getNodeMeans()
		{
			int num_nodes = factory.tree.getNumNodes();
			double[] N = new double[num_nodes]; // return value
			for (int node=0; node<N.length; node++)
				N[node] = getNodeMean(node);
			return N;
		}
		
		public double getEdgeMean(int node)
		{
			double[] Pr = getEdgePosteriors(node);
			return computeExpectation(Pr);
		}

		public double[] getEdgeMeans()
		{
			int num_nodes = factory.tree.getNumNodes();
			double[] S = new double[num_nodes]; // return value
			for (int node=0; node<S.length; node++)
			{				
				S[node] = getEdgeMean(node);
			}
			return S;
		}
		
		
		public double[] getNodeCDF(int node)
		{
			return computeCumulative(getNodePosteriors(node));
		}
		
		public double[] getEdgeCDF(int node)
		{
			return computeCumulative(getEdgePosteriors(node));
		}
		
		public double[] getNodeTail(int node)
		{
			return computeTail(getNodePosteriors(node));
		}
		
		public double[] getEdgeTail(int node)
		{
			return computeTail(getEdgePosteriors(node));
		}
		
		
//		public double getNodeTails(int node)
//		{
//			double[] Ncdf = getNodeCDF(node);
//			double[] Ntail = new double[Ncdf.length];
//			for (int ell=0; ell<Ncdf.length-1; ell++) // last entry is 1.0
//				Ntail[ell]=1.0-Ncdf[ell];
//			
//		}
		
		/**
		 * 
		 * Family event posteriors calculated for conserved ancestral copies.
		 * 
		 * @param node
		 * @return
		 */
		public double[] getFamilyLogSurvivalEventPosteriors(int node)
		{
			double[][] log_death_ls = this.getLogEdgeTransitionPosteriors(node);
			double[][] log_birth_ls = this.getLogNodeTransitionPosteriors(node);

			// family gain for sum_n>0 log_birth[n][0]
			// family loss for sum_l>0 log_death[l][0]
			// family contraction sum_l>1 log_death[l][1]
			// family expansion sum_n>1 log_birth[n][1]
			
			double log_family_gain = Double.NEGATIVE_INFINITY;
			double log_family_loss = Double.NEGATIVE_INFINITY;
			double log_family_expand = Double.NEGATIVE_INFINITY;
			double log_family_contract = Double.NEGATIVE_INFINITY;
			
			for (int ell=1; ell<log_death_ls.length; ell++)
			{
				
				log_family_loss = Logarithms.add(log_family_loss, log_death_ls[ell][0]);
				if (1<ell)
					log_family_contract = Logarithms.add(log_family_contract, log_death_ls[ell][1]);
			}
			for (int ell=1; ell<log_birth_ls.length; ell++)
			{
				log_family_gain = Logarithms.add(log_family_gain, log_birth_ls[ell][0]);
				if (1<ell)
				{
					log_family_expand = Logarithms.add(log_family_expand, log_birth_ls[ell][1]);
				}
			}
			double[] log_event_probs = new double[FamilyEvent.values().length];
			log_event_probs[FamilyEvent.GAIN.ordinal()] = log_family_gain;
			log_event_probs[FamilyEvent.LOSS.ordinal()] = log_family_loss;
			log_event_probs[FamilyEvent.EXPAND.ordinal()] = log_family_expand;
			log_event_probs[FamilyEvent.CONTRACT.ordinal()] = log_family_contract;
			
			return log_event_probs;
		}
		
		
		/**
		 * Posterior probabilities for family-level events on the edge leading to a node.
		 * The probabilities are calculated by subtractions using complementary 
		 * events, but rounded to 0.0 if negative bc precision. 
		 * 
		 * @param node
		 * @return 4-element array indexed in {@link FamilyEvent} order
		 */
		public double[] getFamilyEventPosteriors(int node)
		{
			// TODO use logit parameters ? 
			
			double[] post_node = getNodeAncestorPosteriors(node);
			double[] post_inherit = getEdgeAncestorPosteriors(node);

			double[] event_probs = new double[FamilyEvent.values().length];
			event_probs[FamilyEvent.GAIN.ordinal()] = post_inherit[0]-post_node[0];
			
			if (factory.tree.isRoot(node))
			{
				return event_probs;
			}
			
			int parent = factory.tree.getParent(node);
			double[] post_parent = getNodeAncestorPosteriors(parent);
			
			event_probs[FamilyEvent.LOSS.ordinal()] = post_inherit[0]-post_parent[0];

			double q = factory.rates.getDuplicationParameter(node);
			
			if (q!=0.0) // otherwise expansion and contraction are meaningless 
			{
				// collect the siblings 
				int num_children = factory.tree.getNumChildren(parent);
				int[] siblings = new int[num_children-1];
				int sib_idx=0;
				for (int ci=0; ci<num_children; ci++)
				{
					int sister = factory.tree.getChild(parent, ci);
					if (sister != node)
						siblings[sib_idx++]=sister;
				}
				// now compute the supporting likelihoods across the siblings
				double sib0 = 0.0;
				double sib1 = 0.0;
				--sib_idx;
				while (sib_idx>=0)
				{
					int sister = siblings[sib_idx];
					double[] Kw = inside.getEdgeLikelihoods(sister);
					sib0 += Logarithms.add(sib0, Kw[0]);
					double pw = factory.getLossParameter(sister);
					double w0 = Math.log(pw)+Kw[0];
					double w = w0;
					if (Kw.length>1)
					{
						double w1 = Math.log1p(-pw)+Kw[1];
						w = Logarithms.add(w, w1);
					}
					sib1 += w;
					
					--sib_idx;
				}
				double[] Bu = getNodeAncestor(parent);
				assert Bu.length>1;
				
				double Bu0 = Bu[0]+sib0;
				double Bu1 = Bu[1]+sib1;
				// compute the inside likelihoods conditioned on ancestor copies 0 or 1 
				double epsi = factory.getExtinction(node);
				double[] Cv = inside.getNodeLikelihoods(node);
				double Cv0 = Cv[0];
				double Cv_1 = Cv.length==1?Double.NEGATIVE_INFINITY:Cv[1];
				double Cv1 = Logarithms.add(Math.log(epsi)+Cv[0], Math.log1p(-epsi)+Cv_1);
				
				// transitions from 0/1 copies at parent to 0/1 copies at node
				double p = factory.rates.getLossParameter(node);
				double κ = factory.rates.getGainParameter(node);
				double logp = Math.log(p);
				double log1_p = Math.log1p(-p);
				double logq = Math.log(q);
				double log1_q = Math.log1p(-q);
				double log_κ = Math.log(κ);

				double t1_0 = logp + log1_q*κ ;
				double t101 = logp + log_κ + logq;
				double t111 = log1_p+log1_q;
				double t1_1 = log1_q*κ + Logarithms.add(t101, t111);
				
				double LL = inside.getLogLikelihood(); 
				double Lnoexp = Bu1+Logarithms.add(t1_0+Cv0, t1_1+Cv1); 
				double pnoexp = Math.exp(Lnoexp-LL);
				
				//assert !Double.isNaN(pnoexp);
				
				event_probs[FamilyEvent.EXPAND.ordinal()] = post_parent[1]-pnoexp;
				double t0_1 = log1_q*κ+log_κ + logq;
				double Lnocon = Logarithms.add(Bu0+t0_1, Bu1+t1_1)+Cv1;
				double pnocon = Math.exp(Lnocon-LL);
				
				//assert !Double.isNaN(pnocon);
				
				event_probs[FamilyEvent.CONTRACT.ordinal()] = post_node[1]-pnocon;
			}
			for (int i=0; i<event_probs.length; i++)
			{
				double e = event_probs[i];
				event_probs[i]=Double.max(0.0, e);
			}
			return event_probs;
		}
		
		
		
		
		
		/**
		 * 
		 * @param outsideLL outside ancestor log-likelihoods
		 * @param insideLL inside log-likelihoods
		 * @param extinct extinction probability
		 * @return minimum 2-member array
		 */
		private double[] computeAncestorPosteriors(double[] outsideLL, double[] insideLL, double extinct)
		{
			double[] posteriors;
			if (extinct==0.0)
			{
				posteriors = computePosteriors(outsideLL, insideLL);
			} else
			{
				posteriors = outsideLL.clone();
				double LL = inside.getLogLikelihood();
				double log_e = Math.log(extinct);
				double log_1e = Math.log1p(-extinct);
				
				if (insideLL.length!=0)
				{
	
					for (int n=0; n<posteriors.length; n++)
					{
						int ell = 0;
						double nle = n==0?0.0:n*log_e;
						double ins = insideLL[ell]+n*log_e;
						++ ell;
						while ( ell<=n && ell<insideLL.length)
						{
							double binom = factory.factorials.factln(n)
									-factory.factorials.factln(ell)
									-factory.factorials.factln(n-ell);
							double l1e = ell*log_1e;
							nle = (n==ell?0.0:(n-ell)*log_e);
							double log_p = binom +l1e+nle;
							ins = Logarithms.add(ins, insideLL[ell]+log_p);
							++ ell;
						}
						posteriors[n] += ins;
					}
				}
				
				for (int n=0; n<posteriors.length; n++)
				{
					if (DEBUG_EVENT_POSTERIORS)
					{
						// TODO 
						// DEBUG
						if (LL<posteriors[n])
						{
							double d = posteriors[n]-LL;
							System.out.println("#***P.P.cAP fam "+inside.family_idx+"\text "+extinct+"\tLL "+LL+"\tn "+n+"\tp[n] "+posteriors[n]+"\tdiff "+d+"\t"+Arrays.toString(posteriors));
							int ell = 0;
							double ins = insideLL[ell]+n*log_e;	
							double ipluso = ins + outsideLL[n];
							System.out.println("#***P.P.cAP fam "+inside.family_idx+"\tell "+ell+"\tins "+ins+"\tiLL "+insideLL[ell]+"\toLL "+outsideLL[n]+"\tsum "+ipluso);
							++ ell;
							while ( ell<=n && ell<insideLL.length)
							{
								double binom = factory.factorials.factln(n)
										-factory.factorials.factln(ell)
										-factory.factorials.factln(n-ell);
								double log_p = binom + ell*log_1e+(n-ell)*log_e;
								ins = Logarithms.add(ins, insideLL[ell]+log_p);
								ipluso = ins + outsideLL[n];
								System.out.println("#***P.P.cAP fam "+inside.family_idx+"\tell "+ell+"\tins "+ins+"\tiLL "+insideLL[ell]+"\toLL "+outsideLL[n]+"\tsum "+ipluso);
								++ ell;
							}
							ipluso = ins + outsideLL[n];
							System.out.println("#***P.P.cAP fam "+inside.family_idx+"\tins "+ins+"\toLL "+outsideLL[n]+"\tsum "+ipluso);
							
							return null;
						}
					}

					posteriors[n]= Math.exp(posteriors[n]-LL);
					
					if (DEBUG_EVENT_POSTERIORS) assert (posteriors[n]<=1.0);
					else posteriors[n] = Double.min(posteriors[n],1.0); // numerical precision
//					
				}
			}
			if (posteriors.length==1)
			{
				posteriors = Arrays.copyOf(posteriors, 2);
			}
			
			return posteriors;
		}		
				
		protected double[] getNodeOutside(int node)
		{
			double[] B = node_outside[node];
			if (B==null) 
			{
				computeLikelihoods();
				B = node_outside[node];
			}
			return B;
		}
		
		protected double[] getEdgeOutside(int node)
		{
			double[] J = edge_outside[node];
			if (J==null)
			{
				computeLikelihoods();
				J = edge_outside[node];
			}
			return J;
		}
		
		private double[] getNodeAncestor(int node)
		{
			if (factory.tree.isLeaf(node))
			{
				return getNodeOutside(node);
			}
			
			if (ancestor==null)
				computeAncestors();
			
			Profile A = ancestor[node];
			double[] B = A.node_outside[node];
			return B;
		}
		
		private double[] getEdgeAncestor(int node)
		{
			if (factory.tree.isLeaf(node))
			{
				return getEdgeOutside(node);
			}

			if (ancestor==null)
				computeAncestors();
			Profile A = ancestor[node];
			double[] J = A.edge_outside[node];
			return J;
		}
		
		
//		private double[] computeNode(int node)
//		{
//			double[][] Binc = computeNodeTransitions(node);
//			return Binc[0];
//		}
		
		
		/**
		 * Algorithm for computing the outside log-likelihood at a node. The profile 
		 * must not be ambiguous at the root. 
		 * 
		 * @param node called in preorder, after setting edge outer likelihood here, and node outer at the parent
		 * @return
		 */
		private double[][] computeNodeTransitions(int node)
		{
//			if (factory.gain_factorials[node]==null && factory.rates.getDuplicationParameter(node)!=0.0)
//			{
//				if (inside.family_idx==0)
//					System.out.println("#***P.P.cNT fam "+inside.family_idx+"\tnode "+node+"\t"+factory.rates.toString(node));
//			}
			
			double[] C = inside.getNodeLikelihoods(node);
			//double[] B; // return value: outside node likelihoods
			double[][] Bls; // return value outside transitions
			double[] J = edge_outside[node]; // already set

			assert (C.length != 0);
//			if (C.length==0)
//			{
//				// ambiguous
//				assert (!factory.tree.isRoot(node)); 
//				int parent = factory.tree.getParent(node);
//				int n = node_outside[parent].length;
//				B = new double[n+1]; // all 1 
//
//			} else
//			{
//			 	B = new double[C.length];
				Bls = new double[C.length][];
//			}
//			double q = factory.getDuplicationParameter(node);
			double logq = factory.getLogDuplicationParameter(node); //   Math.log(q);
			double log1_q = factory.getLogDuplicationComplement(node); // Math.log1p(-q);// Math.log(1.0-q);

			double log_gain = factory.getLogGainParameter(node);
			
			if  (logq==Double.NEGATIVE_INFINITY)  // (q==0.0)
			{ 
				
				// Poisson
				//double r = factory.getGainParameter(node);
				double logr =log_gain;
				double r = Math.exp(logr);
				
				for (int ell=0; ell<Bls.length; ell++)
				{
					double[] terms = new double[ell+1];
					int t=ell,s=0;
					//assert (J.length<=ell+1);
					while (s<J.length && s<=ell)
					{
						double t_logr = (t==0?0.0:t*logr);
						terms[s] = J[s] - r + t_logr - factory.factorials.factln(t);
						
//						if (!(terms[s]<=0.0))
//						{
//							System.out.println("#***P.P.cNT node "+node+"\tBls[ell="+ell+"][s="+s+"] "+terms[s]+"\tJs "+J[s]+"\tr "+r+"\ttlogr "+t_logr+"\tln(t!) "+factory.factorials.factln(t));
//						}
						assert terms[s]<=0.0;
						
						
						++s;
						--t;
					}
					if (s<terms.length)
						Arrays.fill(terms,s,terms.length,Double.NEGATIVE_INFINITY);
//					B[ell] = Logarithms.sum(terms, s);
					Bls[ell] = terms;
				}
			} else
			{
				// Pólya
				double log_kappa = log_gain;
				
				//L double κ = factory.getGainParameter(node);
//				double logq = Math.log(q);
//				double log1_q = Math.log1p(-q);// Math.log(1.0-q);
				
				if (log_kappa==Double.NEGATIVE_INFINITY) // untested
				{
					double[] terms =new double[1];
					// negative binomial with s and q 
					terms[0] = J[0]; // no conservation
					assert (J[0]<=0.0);
//					B[0]=J[0]; // no conservation
					Bls[0]=terms;
//					Bplus[0]=Double.NEGATIVE_INFINITY;
					for (int ell=1; ell<Bls.length; ell++)
					{
						terms = new double[ell+1];
						double factln_ell = factory.factorials.factln(ell-1); // ell-1 and not ell
						
						int t=ell-1, sm1=0, s=1; // no contribution from s=0
						//assert (J.length<=ell+1);
						while (s<J.length && s<=ell) // sm1==s-1; t+s=ell
						{
							double binom = factln_ell - factory.factorials.factln(sm1) - factory.factorials.factln(t); // s-1 and not s
							// double binom = factln_ell - factory.factorials.factln(sm1)- factory.factorials.factln(t);
							
							terms[sm1] = J[s]+binom + s*log1_q + t*logq;
							
//							if (!(terms[sm1]<=0.0))
//							{
//								System.out.println("#***P.P.cNT node "+node+"\tBls[ell="+ell+"][s-1="+sm1+"] "+terms[s]+"\tt "+t+"\tJs "+J[s]+"\tbinom "+binom
//										+"\t("+factln_ell+","+factory.factorials.factln(sm1)+","+factory.factorials.factln(t)+")"
//										+"\t(k+s)logq "+(s*log1_q)+"\ttlogq "+(t*logq)+"\tln(t!) "+factory.factorials.factln(t));
//							}
							assert (terms[sm1]<=0.0);
							
							--t;
							++s;
							++sm1;
						}
						if (s<terms.length)
							Arrays.fill(terms,s,terms.length,Double.NEGATIVE_INFINITY);
//						Bplus[ell] = logCumulativeNodeIncrease(terms, sm1, ell);
//						B[ell] = Logarithms.sum(terms, sm1);
						Bls[ell] = terms;
					}
				} else // negative binomial with 0<kappa
				{
					
					double loglog1_q = Logarithms.logitToLogLogComplement(logq-log1_q);
					for (int ell=0; ell<Bls.length; ell++)
					{
						double[] terms = new double[ell+1];
						double log_ellfact = factory.gain_factorials[node].factln(ell);
						int t=ell, s=0; 
						while (s<J.length && s<=ell)
						{
							double log_kappa_s = Logarithms.add(log_kappa, Math.log(s));
							double log_kappa_slog1_q = log_kappa_s + loglog1_q;
							double kappa_s_log1_q = -Math.exp(log_kappa_slog1_q); // we want -(-ln(1-q))*(kappa+s)
							double t_logq= t*logq;
							double binom = log_ellfact - factory.gain_factorials[node].factln(s) - factory.factorials.factln(t);
							//L terms[s] = J[s] + binom + (κ+s)*log1_q + t * logq;
							terms[s] = J[s] + binom + kappa_s_log1_q + t_logq;
							
//							if (!(terms[s]<=0.0))
//							{
//								System.out.println("#***P.P.cNT node "+node+"\tBls[ell="+ell+"][s="+s+"] "+terms[s]+"\tt "+t+"\tJs "+J[s]+"\tbinom "+binom
//										+"\t("+log_ellfact+","+factory.gain_factorials[node].factln(s)+","+factory.factorials.factln(t)+")"
//										+"\t(k+s)logq "+kappa_s_log1_q+"\ttlogq "+t_logq+"\tln(t!) "+factory.factorials.factln(t));
//							}
//							
//							
//							// DEBUG
//							if (0.0<terms[s])
//							{
//								double f1 = log_kappa;
//								double f2 = f1+Logarithms.add(log_kappa,1.0);
//
//								System.out.println("#***P.cNT node "+node+"\ts,t,ell "+s+","+t+","+ell
//										+"\tJs "+J[s]
//										+"\tksl1q "+kappa_s_log1_q
//										+"\ttlq "+t_logq+"\tbin "+binom
//										+"\tlell! "+log_ellfact+"\tls! "+factory.gain_factorials[node].factln(s)+"\tlt! "+factory.factorials.factln(t)
//										+"\tf1 "+f1+"\tf2 "+f2
//										+"\t// "+factory.gain_factorials[node].toString());
//								
//								
//							}
								
							
							assert (terms[s]<=0.0); // since it is a probability
							
							++s;
							--t;
						}
//						Bplus[ell] = logCumulativeNodeIncrease(terms, s, ell);
//						B[ell] = Logarithms.sum(terms, s);
						if (s<terms.length)
							Arrays.fill(terms,s,terms.length,Double.NEGATIVE_INFINITY);
						Bls[ell] = terms;
					}
				}
			}
			return Bls;
		}
		
//		private double logCumulativeNodeIncrease(double[] Bterms, int nterms, int ell)
//		{
//			double tot = Double.NEGATIVE_INFINITY; // log(0)
//			double tottot = Double.NEGATIVE_INFINITY;
//			int s =0; 
//			while (s<ell && s<nterms)
//			{
//				tot = Logarithms.add(tot, Bterms[s]);
//				tottot = Logarithms.add(tottot, tot);
//				s++;
//			}
//			return tottot;
//		}
//		
//		
//
				
		/**
		 * Algorithm for computing the outside log-likelihood on an edge. 
		 * @param node called in preorder 
		 * @return
		 */
		private double[][] computeEdgeTransitions(int node)
		{
			double[] J; // return value
			double[][] Jls; 
			double logp = factory.getLogLossParameter(node);
			if (factory.tree.isRoot(node))
			{
				
//				double proot = factory.getLossParameter(node);
				if (factory.getLogLossComplement(node)==Double.NEGATIVE_INFINITY) // (proot == 1.0)
				{
					J = new double[1]; // and the value is log(1)
					Jls = new double[1][1];
				} else { // untested
					double log1_p = factory.getLogLossComplement(node);
					J = new double[2];
					J[0] = logp; //Math.log(proot);
					J[1] = log1_p; //  Math.log1p(-proot);
					Jls = new double[2][];
					Jls[0] = new double[1];
					Jls[0][0] = Double.NEGATIVE_INFINITY;
					Jls[1] = new double[2];
					Jls[1][0] = logp;
					Jls[1][1] = log1_p;
				}
			} else
			{
//				double p = factory.getLossParameter(node);
				int parent = factory.tree.getParent(node);
				double[] B = node_outside[parent];
				
				// TODO: test of allowing p=1.0 works here
//				if (factory.getLogLossComplement(node)==Double.NEGATIVE_INFINITY) //(logp == 0.0) //(p==1.0)
//				{
//					J = new double[1];// and the value is log(1)
//					Jls = new double[B.length][];
//					for (int ell=0; ell<Jls.length; ell++)
//					{
//						Jls[ell] = new double[ell+1];
//						Arrays.fill(Jls[ell], Double.NEGATIVE_INFINITY);
//						Jls[ell][0] = B[ell];
//						
//						// TODO fix this calculation: K2[ell] should be there?
//					}
//				} else
				{
					int n = (inside.getNodeLikelihoods(node)).length;
					if (n==0)
					{
						// ambiguous: copy from parent
						n = B.length;
					}
					J = new double[n];
					Jls = new double[B.length][];
					for (int ell=0; ell<Jls.length; ell++)
					{
						Jls[ell] = new double[ell+1];
						Arrays.fill(Jls[ell], Double.NEGATIVE_INFINITY);
					}

					// collect the siblings
					int num_children = factory.tree.getNumChildren(parent);
					int[] siblings = new int[num_children-1];
					int sib_idx=0;
					for (int ci=0; ci<num_children; ci++)
					{
						int sister = factory.tree.getChild(parent, ci);
						if (sister != node)
							siblings[sib_idx++]=sister;
					}
					// now compute the inside likelihood across the siblings
					--sib_idx;
					int sister=siblings[sib_idx];
//					double eps = factory.getLossParameter(sister); // extinction across siblings
					
					double loge = factory.getLogLossParameter(sister);
					double log1_e = factory.getLogLossComplement(sister);
					double logit_e = loge - log1_e;
					
					double[] K2 = inside.getEdgeLikelihoods(sister).clone();
					while (sib_idx>0)
					{
						--sib_idx;
						int junior = siblings[sib_idx];
//						K2 = inside.computeSibling(junior, K2, eps);
						
						K2 = inside.computeSiblingLogit(junior, K2, logit_e);
						
//						eps *= factory.getLossParameter(junior);
						double logit_pj = factory.getLogLossParameter(junior)-factory.getLogLossComplement(junior);
						logit_e = Logarithms.mulLogit(logit_e, logit_pj);
					}
					if (K2.length==0)
					{
						// all ambiguous on the other side 
						K2 = new double[B.length]; 
						// filled with probability == 1
					} 
					
					
					// combine siblings and the outside at the parent
					// get the numerical parameters for the formulas
					
//					double loga = Math.log1p(-p*eps); // log(1-pe)
//					double logp1 = Math.log1p(-p)-loga; // log((1-p)/(1-pe))
//					double logp2 = Math.log(p)+Math.log1p(-eps)-loga; // log(p*(1-e)/(1-pe))
//					loge = Math.log(eps); // log(e)
//					log1_e = Math.log1p(-eps); // log(1-e)

					// Logit-scale:
					if (0<=logit_e)
					{
						loge = Logarithms.logitToLogValue(logit_e);
						log1_e = loge-logit_e;
					} else
					{
						log1_e = Logarithms.logitToLogComplement(logit_e);
						loge = log1_e+logit_e;
					}
					double logit_p = logp - factory.getLogLossComplement(node);
					double logit_p2 = logit_p + log1_e;
					double logp1,logp2;
					if (0.0<=logit_p2)
					{
						logp2 = Logarithms.logitToLogValue(logit_p2);
						logp1 = logp2-logit_p2;
					} else
					{
						logp1 = Logarithms.logitToLogComplement(logit_p2);
						logp2 = logp1 + logit_p2;
					}
					
					
					
					
	
					
					for (int s=0; s<J.length; s++)
					{
						if (s>0)
						{
							// update K2
							// recursion: K2(s+t, t) = (1-e)*K2(s-1+t+1, t+1)+e*K2(s-1+t,t);
							{
								int t = K2.length-1;
								double y = K2[t];
								K2[t]+=loge;
								while (t>0)
								{
									--t;
									double x = y+log1_e;
									y = K2[t];
									double z = y+loge;
									K2[t] = Logarithms.add(x, z);
								}
							}
						}
						
						double[] terms = new double[B.length]; 
						
						int t=0;
						int ell=s+t;
						double log_sfact = factory.factorials.factln(s);
						while (t<terms.length 
								&& t<K2.length) // if not truncated, t<K2.length is enough
						{						
							double binom = factory.factorials.factln(ell)-factory.factorials.factln(t)-log_sfact;
	//						assert (ell<B.length); // if not truncated
							double B_ell = (ell<B.length?B[ell]:Double.NEGATIVE_INFINITY);
	//						assert (t<K2.length); // if not truncated
	//						assert (t<terms.length);
							double s_logp1 = s==0?0.0:s*logp1;
							double t_logp2 = t==0?0.0:t*logp2;
							terms[t] = B_ell + K2[t] + binom + s_logp1 + t_logp2;
							
//							if (!(terms[t]<=0.0))
//								System.out.println("#**P.P.cET node "+node+"\tterms[t="+t+"] "+terms[t]+"\ts "+s+"\tell "+ell+"\tBll "+B_ell+"\tK2 "+K2[t]+"\tbinom "+binom+"\tlogp1 "+logp1+"\tlogp2 "+logp2);

							assert (terms[t]<=0.0); // since it is a probability
							
							
							if (ell<Jls.length)
								Jls[ell][s] = terms[t];
							
							
							++t;
							++ell;
						}
						J[s] = Logarithms.sum(terms, t);
					}
				}
			}
			return Jls;
//			return J;
		}
	} // Profile inner class
	
	private static double mean(double[] p)
	{
		double m = 0.0;
		for (int i=1; i<p.length; i++)
			m += i*p[i];
		return m;
	}
	
	
	private double node_miss=0.0;
	private double edge_miss=0.0;
	
	
	/**
	 * Pretty print for inferred posteriors. 
	 * 
	 * @param family_idx family index 
	 * @param out where to print
	 */
	public void printPosteriors(PrintStream out, int family_idx)
	{
        int nNodes = factory.tree.getNumNodes();
        int nLeaves = factory.tree.getNumLeaves();
        
    	Profile P = getPosteriors(family_idx);

    	out.println("# Node conserved posteriors for family "+family_idx);
    	
    	for (int node=0; node<nNodes; node++ )
    	{
    		double[] pN = P.getNodePosteriors(node);
    		double mN = mean(pN);
    		out.print(family_idx+"\tNP\t"+node+"\t"+P.getNodeMean(node)+"\t("+mN+"/"+pN.length+")");
    		for (int n=0; n<pN.length; n++)
    		{
    			out.print("\t"+pN[n]);
    		}
    		out.println();
    	}
    	out.println("\n# Node raw posteriors for family "+family_idx);
    	
    	for (int node=nLeaves; node<nNodes; node++ )
    	{
    		double[] pN = P.getNodeAncestorPosteriors(node);
    		double mN = mean(pN);
    		out.print(family_idx+"\tNA\t"+node+"\t"+mN+"\t/"+pN.length);
    		double miss = 1.0;
    		for (int n=0; n<pN.length; n++)
    		{
    			out.print("\t"+pN[n]);
    			miss-=pN[n];
    		}
    		out.println("\t// missing "+miss);
    		node_miss = Double.max(node_miss, miss);
    	}

    	out.println("\n# Edge conserved posteriors for family "+family_idx);
    	for (int node=0; node<nNodes; node++ )
    	{
    		double[] pS = P.getEdgePosteriors(node);
    		double mS = mean(pS);
    		out.print(family_idx+"\tEP\t"+node+"\t"+P.getEdgeMean(node)+"\t("+mS+"/"+pS.length+")");
    		for (int s=0; s<pS.length; s++)
    		{
    			out.print("\t"+pS[s]);
    		}
    		out.println();
    	}
    	out.println("\n# Edge raw posteriors for family "+family_idx);
    	for (int node=0; node<nNodes; node++ )
    	{
    		double[] pS = P.getEdgeAncestorPosteriors(node);
    		double mS = mean(pS);
    		out.print(family_idx+"\tEA\t"+node+"\t"+mS+"\t/"+pS.length);
    		double miss=1.0;
    		for (int s=0; s<pS.length; s++)
    		{
    			out.print("\t"+pS[s]);
    			miss-=pS[s];
    		}
    		out.println("\t missing "+miss);
    		edge_miss=Double.max(edge_miss,miss);
    	}
    	out.println(family_idx+"\tLL\t*\t"+P.inside.getLogLikelihood());
//        out.println("\n# Log-likelihood for family "+family_idx+"\t"+P.inside.getLogLikelihood());
    }
	
	
	private void printPosteriors(PrintStream out)
	{
		
	}
	
	/**
	 * Testing from command-line.
	 * 
	 * @param args
	 * @throws Exception
	 */
	private static void mainmain(String[] args) throws Exception
	{
		CommandLine cli = new CommandLine(args, Posteriors.class);
//        int arg_idx = 0;
//        if (3+arg_idx != args.length)
//            throw new IllegalArgumentException("Call as java "+getClass().getName()+" tree table rates");
//
//        String tree_file = args[arg_idx++];
//        String table_file = args[arg_idx++];
//        String rates_file = args[arg_idx];
//		
//        Phylogeny tree = count.io.NewickParser.readTree(new java.io.FileReader(tree_file));
//        AnnotatedTable table = TableParser.readTable(tree.getLeafNames(), 
//        		GeneralizedFileReader.guessReaderForInput(table_file), true);
//        TreeWithRates rates = new TreeWithRates(tree);
//        java.io.BufferedReader R = new java.io.BufferedReader(GeneralizedFileReader.guessReaderForInput(rates_file));
//        RateVariationParser.initFromFile(rates, R);
//        R.close();
//
        PrintStream out = System.out;
//        
//        out.println(Executable.getStandardHeader(this.getClass()));
//        out.println(Executable.getStandardRuntimeInfo());
//        out.println(Executable.getStandardHeader("Tree file:  "+tree_file));
//        out.println(Executable.getStandardHeader("Table file: "+table_file));
//        out.println(Executable.getStandardHeader("Rates file: "+rates_file));
		
		TreeWithRates rates = cli.getRates();
		AnnotatedTable table = cli.getTable();
		double ancestor_width = 3.0;
		
//		out.println(RateVariationParser.printRates(cli.getModel()));
		out.println("# Empty profile ");
        Posteriors empty = new Posteriors(rates, ProfileTable.emptyProfile(rates.getTree()));
        empty.setAncestorWidthThreshold(ancestor_width);
        empty.printPosteriors(out, 0);
			
        Posteriors post = new Posteriors(rates, table);
        
    	int absolute = cli.getOptionTruncateAbsolute();
    	double relative = cli.getOptionTruncateRelative();
        if (cli.getOptionValue(CommandLine.OPT_TRUNCATE)!=null)
        {
    		out.println(CommandLine.getStandardHeader("Truncated computation (absolute,relative)="
        				+absolute+","+relative));
        	post.setCalculationWidthThresholds(absolute, relative);
		}
        
        post.setAncestorWidthThreshold(ancestor_width);
        
        int nFam = table.getFamilyCount();
        
        for (int f=0; f<nFam; f++)
        	post.printPosteriors(out, f);
        
        out.println("# Raw MISS max\t"+post.node_miss+"\t"+post.edge_miss);
        
        
	}	
	
	public static void main(String[] args) throws Exception
	{
		//(new Posteriors()).
		mainmain(args);
	}
	
	
}
