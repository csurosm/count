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
import count.io.GeneralizedFileReader;
import count.io.RateVariationParser;
import count.io.TableParser;
import count.matek.Logarithms;
import count.util.Executable;

import java.io.PrintStream;
import java.util.Arrays;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.Phylogeny;

public class Posteriors
{
	public Posteriors(RateModel.GLD rates, ProfileTable table)
	{
		this.factory = new Likelihood(rates, table);
	}
	
	/**
	 * Dummy instantiation not for computing with.
	 */
	private Posteriors() {this.factory=null;}
	
	Profile getPosteriors(int family_idx)
	{
		
		Profile P = new Profile(family_idx);
		for (int node=factory.tree.getNumNodes(); node>0; )
		{
			--node;
			System.out.println("#*P.getP family "+family_idx+"\tedge "+node+"\tJ[] = "
					+Arrays.toString(P.getEdgeOutside(node)));
			System.out.println("#*P.getP family "+family_idx+"\tnode "+node+"\tB[] = "
					+Arrays.toString(P.getNodeOutside(node)));
		}
    	
		return P;
	}
	
	final Likelihood factory;	
	
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
		
		private final Likelihood.Profile inside; 
		
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
		
		public double[] getNodePosteriors(int node)
		{
			double[] B = getNodeOutside(node); // sets up inside likelihoods if necessary
			double[] C = inside.getNodeLikelihoods(node);
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
		
		double[] computeCumulative(double[] Pr)
		{
			double prev = Pr[0];
			for (int i=1; i<Pr.length; i++)
			{
				prev = Pr[i] = Pr[i]+prev;
			}
			return Pr;
		}
		
		double computeExpectation(double[] Pr)
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
		 * @return
		 */
		private double[] computePosteriors(double[] outside, double[] inside)
		{
			double[] posteriors = outside.clone();
			if (inside.length!=0)
			{
				// not ambiguous
				for (int ell=0; ell<posteriors.length; ell++)
				{
					posteriors[ell]+=inside[ell];
				}
			}
			double LL = Logarithms.sum(posteriors, posteriors.length);
			for (int ell=0; ell<posteriors.length; ell++)
			{
				posteriors[ell]=Math.exp(posteriors[ell]-LL);
			}
			return posteriors;
		}
		
		private double[] getNodeOutside(int node)
		{
			double[] B = node_outside[node];
			if (B==null) 
			{
				computeLikelihoods();
				B = node_outside[node];
			}
			return B;
		}
		
		private double[] getEdgeOutside(int node)
		{
			double[] J = edge_outside[node];
			if (J==null)
			{
				computeLikelihoods();
				J = edge_outside[node];
			}
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
			{
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
				double log1_q = Math.log(1.0-q);
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
				double p = factory.getLossParameter(node);
				
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
					while (t<K2.length)
					{						
						double binom = factory.factorials.factln(ell)-factory.factorials.factln(t)-log_sfact;
						terms[t] = B[ell] + K2[t] + binom + s*logp1 + t*logp2;
						++t;
						++ell;
					}
					J[s] = Logarithms.sum(terms, t);
				}
			}
			return J;
		}
	} // Profile inner class
	
	/**
	 * Pretty print for inferred posteriors. 
	 * 
	 * @param family_idx family index 
	 * @param out where to print
	 */
	public void printPosteriors(PrintStream out, int family_idx)
	{
        int nNodes = factory.tree.getNumNodes();

    	out.println("# Node posteriors for family "+family_idx);
    	Profile P = getPosteriors(family_idx);
    	
    	for (int node=0; node<nNodes; node++ )
    	{
    		double[] pN = P.getNodePosteriors(node);
    		out.print(family_idx+"\t"+node+"\t"+P.getNodeMean(node));
    		for (int n=0; n<pN.length; n++)
    		{
    			out.print("\t"+pN[n]);
    		}
    		out.println();
    	}
    	out.println("# Edge posteriors for family "+family_idx);
    	for (int node=0; node<nNodes; node++ )
    	{
    		double[] pS = P.getEdgePosteriors(node);
    		out.print(family_idx+"\t"+node+"\t"+P.getEdgeMean(node));
    		for (int s=0; s<pS.length; s++)
    		{
    			out.print("\t"+pS[s]);
    		}
    		out.println();
    	}
    }
	


	
	/**
	 * Instance-linked main for testing from command-line.
	 * 
	 * @param args
	 * @throws Exception
	 */
	private void mainmain(String[] args) throws Exception
	{
        int arg_idx = 0;
        if (3+arg_idx != args.length)
            throw new IllegalArgumentException("Call as java "+getClass().getName()+" tree table rates");

        String tree_file = args[arg_idx++];
        String table_file = args[arg_idx++];
        String rates_file = args[arg_idx];
		
        Phylogeny tree = count.io.NewickParser.readTree(new java.io.FileReader(tree_file));
        AnnotatedTable table = TableParser.readTable(tree.getLeafNames(), 
        		GeneralizedFileReader.guessReaderForInput(table_file), true);
        TreeWithRates rates = new TreeWithRates(tree);
        java.io.BufferedReader R = new java.io.BufferedReader(GeneralizedFileReader.guessReaderForInput(rates_file));
        RateVariationParser.initFromFile(rates, R);
        R.close();

        PrintStream out = System.out;
        
        out.println(Executable.getStandardHeader(this.getClass()));
        out.println(Executable.getStandardRuntimeInfo());
        out.println(Executable.getStandardHeader("Tree file:  "+tree_file));
        out.println(Executable.getStandardHeader("Table file: "+table_file));
        out.println(Executable.getStandardHeader("Rates file: "+rates_file));
        
        Posteriors post = new Posteriors(rates, table);
        int nFam = table.getFamilyCount();
        
        for (int f=0; f<nFam; f++)
        	post.printPosteriors(out, f);
	}	
	
	public static void main(String[] args) throws Exception
	{
		(new Posteriors()).mainmain(args);
	}
	
	
}
