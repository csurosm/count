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

public class Posteriors
{
	public Posteriors(TreeWithRates rates, ProfileTable table)
	{
		this(new Likelihood(rates, table));
	}
	
	/**
	 * Dummy instantiation not for computing with.
	 */
	private Posteriors() {this(null);}
	
	public Posteriors(Likelihood factory) 
	{ 
		this.factory = factory;
	}
	
	Profile getPosteriors(int family_idx)
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
	
	final Likelihood factory;	
	Posteriors[] ancestor_posteriors;
	
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
	
	static double computeExpectation(double[] Pr)
	{
		double expect = 0.0;
		for (int ell=1; ell<Pr.length; ell++)
			expect += ell*Pr[ell];
		return expect;
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
	

	
	class Profile 
	{
		private Profile(int family_idx)
		{
			this.inside = factory.getProfileLikelihood(family_idx);
			int num_nodes = factory.tree.getNumNodes();
			node_outside = new double[num_nodes][];
			edge_outside = new double[num_nodes][];
		}
		private final double[][] node_outside;
		private final double[][] edge_outside;
		
		final Likelihood.Profile inside;
		
		private Profile[] ancestor=null;
		
		public void computeLikelihoods()
		{
			inside.computeLikelihoods();
			int root = factory.tree.getRoot();
			for (int node=root; node>=0; --node)
			{
				edge_outside[node] = computeEdge(node);
				node_outside[node] = computeNode(node);
			}
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
				edge_outside[a] = computeEdge(a);
				node_outside[a] = computeNode(a);
			}
		}
		
		
		public double[] getNodeAncestorPosteriors(int node) 
		{
			double[] B = getNodeAncestor(node);
			double[] C = inside.getNodeLikelihoods(node);
			double epsi = factory.extinction[node];

			return computeAncestorPosteriors(B, C, epsi);
		}
		
		public double[] getEdgeAncestorPosteriors(int node)
		{
			double[] J = getEdgeAncestor(node);
			double[] K = inside.getEdgeLikelihoods(node);
			double epsi = factory.extinction[node] * factory.getDuplicationParameterComplement(node);
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
		
		
		
//		public double getNodeTails(int node)
//		{
//			double[] Ncdf = getNodeCDF(node);
//			double[] Ntail = new double[Ncdf.length];
//			for (int ell=0; ell<Ncdf.length-1; ell++) // last entry is 1.0
//				Ntail[ell]=1.0-Ncdf[ell];
//			
//		}
		
		
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
				double epsi = factory.extinction[node];
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
				
				if (insideLL.length!=0)
				{
	
					double log_e = Math.log(extinct);
					double log_1e = Math.log1p(-extinct);
					
					
					for (int n=0; n<posteriors.length; n++)
					{
						int ell = 0;
						double ins = insideLL[ell]+n*log_e;
						++ ell;
						while ( ell<=n && ell<insideLL.length)
						{
							double binom = factory.factorials.factln(n)
									-factory.factorials.factln(ell)
									-factory.factorials.factln(n-ell);
							double log_p = binom + ell*log_1e+(n-ell)*log_e;
							ins = Logarithms.add(ins, insideLL[ell]+log_p);
							++ ell;
						}
						posteriors[n] += ins;
					}
				}
				for (int n=0; n<posteriors.length; n++)
				{
					posteriors[n]=Math.exp(posteriors[n]-LL);
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
		
		
		/**
		 * Algorithm for computing the outside log-likelihood at a node. The profile 
		 * must not be ambiguous at the root. 
		 * 
		 * @param node called in preorder, after setting edge outer likelihood here, and node outer at the parent
		 * @return
		 */
		private double[] computeNode(int node)
		{
			double[] C = inside.getNodeLikelihoods(node);
			double[] B; // return value
			double[] J = edge_outside[node]; // already set

			if (C.length==0)
			{
				// ambiguous
				assert (!factory.tree.isRoot(node)); 
				int parent = factory.tree.getParent(node);
				int n = node_outside[parent].length;
				B = new double[n+1];
			} else
			{
				B = new double[C.length];
			}
			double q = factory.getDuplicationParameter(node);
			if (q==0.0)
			{ // untested
				
				// Poisson
				double r = factory.getGainParameter(node);
				double logr = Math.log(r);
				
				double[] terms = new double[B.length];
				for (int ell=0; ell<B.length; ell++)
				{
					int t=ell,s=0;
					assert (J.length<=ell+1);
					while (s<J.length)
					{
						terms[s] = J[s] - r + t*logr - factory.factorials.factln(t);
						++s;
						--t;
					}
					B[ell] = Logarithms.sum(terms, s);
				}
			} else
			{
				// Pólya
				double κ = factory.getGainParameter(node);
				double logq = Math.log(q);
				double log1_q = Math.log1p(-q);// Math.log(1.0-q);
				double[] terms =new double[B.length];
				
				if (κ==0.0)
				{
					// negative binomial with s and q 
					B[0]=J[0]; // no conservation
					for (int ell=1; ell<B.length; ell++)
					{
						double factln_ell = factory.factorials.factln(ell-1); // ell-1 and not ell
						
						
						int t=ell-1, sm1=0, s=1; // no contribution from s=0
						assert (J.length<=ell+1);
						while (s<J.length) // sm1==s-1; t+s=ell
						{
							double binom = factln_ell - factory.gain_factorials[node].factln(sm1) - factory.factorials.factln(t); // s-1 and not s
							terms[sm1] = J[s]+binom + s*log1_q + t*logq;
							--t;
							++s;
							++sm1;
						}
						B[ell] = Logarithms.sum(terms, sm1);
					}
				} else
				{
					for (int ell=0; ell<B.length; ell++)
					{
						double log_ellfact = factory.gain_factorials[node].factln(ell);
						int t=ell, s=0; 
						while (s<J.length && s<=ell)
						{
							double binom = log_ellfact - factory.gain_factorials[node].factln(s) - factory.factorials.factln(t);
							terms[s] = J[s] + binom + (κ+s)*log1_q + t * logq;
							++s;
							--t;
						}
						B[ell] = Logarithms.sum(terms, s);
					}
				}
			}
			return B;
		}
		
		/**
		 * Algorithm for computing the outside log-likelihood on an edge. 
		 * @param node called in preorder 
		 * @return
		 */
		private double[] computeEdge(int node)
		{
			double[] J; // return value
			if (factory.tree.isRoot(node))
			{
				double proot = factory.getLossParameter(node);
				if (proot == 1.0)
				{
					J = new double[1];
					// and the value is log(1)
				} else {
					J = new double[2];
					J[0] = Math.log(proot);
					J[1] = Math.log1p(-proot);
				}
			} else
			{
				double p = factory.getLossParameter(node);
				if (p==1.0)
				{
					J = new double[1];
					// and the value is log(1)
				} else
				{
					int parent = factory.tree.getParent(node);
					double[] B = node_outside[parent];
	
					int n = (inside.getNodeLikelihoods(node)).length;
					if (n==0)
					{
						// ambiguous: copy from parent
						n = B.length;
					}
					J = new double[n];
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
					double eps = factory.getLossParameter(sister); // extinction across siblings
					double[] K2 = inside.getEdgeLikelihoods(sister).clone();
					while (sib_idx>0)
					{
						--sib_idx;
						int junior = siblings[sib_idx];
						K2 = inside.computeSibling(junior, K2, eps);
						eps *= factory.getLossParameter(junior);
					}
					if (K2.length==0)
					{
						// all ambiguous on the other side 
						K2 = new double[B.length]; 
						// filled with probability == 1
					} 
					
					
					// combine siblings and the outside at the parent
					// get the numerical parameters for the formulas
					
					double loga = Math.log1p(-p*eps); // log(1-pe)
					double logp1 = Math.log1p(-p)-loga; // log((1-p)/(1-pe))
					double logp2 = Math.log(p)+Math.log1p(-eps)-loga; // log(p*(1-e)/(1-pe))
					double loge = Math.log(eps); // log(e)
					double log1_e = Math.log1p(-eps); // log(1-e)
	
					double[] terms = new double[B.length]; // auxiliary reused 
					
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
							terms[t] = B_ell + K2[t] + binom + s*logp1 + t*logp2;
							++t;
							++ell;
						}
						J[s] = Logarithms.sum(terms, t);
					}
				}
			}
			return J;
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
    		out.print(family_idx+"\tNP\t"+node+"\t"+P.getNodeMean(node)+" ("+mN+"/"+pN.length+")");
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
    		out.print(family_idx+"\tNA\t"+node+"\t"+mN+"/"+pN.length);
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
    		out.print(family_idx+"\tEP\t"+node+"\t"+P.getEdgeMean(node)+" ("+mS+"/"+pS.length+")");
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
    		out.print(family_idx+"\tEA\t"+node+"\t"+mS+"/"+pS.length);
    		double miss=1.0;
    		for (int s=0; s<pS.length; s++)
    		{
    			out.print("\t"+pS[s]);
    			miss-=pS[s];
    		}
    		out.println("\t missing "+miss);
    		edge_miss=Double.max(edge_miss,miss);
    	}
        out.println("\n# Log-likelihood for family "+family_idx+"\t"+P.inside.getLogLikelihood());
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
