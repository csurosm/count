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

import java.io.PrintStream;
import java.util.Arrays;



import count.matek.Functions.RisingFactorial;
import count.util.Executable;
import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.TreeTraversal;
import count.ds.Phylogeny;
import count.ds.ProfileTable;
import count.io.GeneralizedFileReader;
import count.io.RateVariationParser;
import count.io.TableParser;
import count.matek.Logarithms;

public class Likelihood implements GLDParameters
{
	private static final boolean USE_CLASSIC_RECURSION = false; // in Count 12
	
	/**
	 * Dummy instantiation not for computing with.
	 */
	private  Likelihood()
	{
		this.rates = null;
		this.table = null;
		this.tree = null;
		this.factorials=null;
		this.max_family_size=null;
	}
	
	public Likelihood(RateModel.GLD rates, ProfileTable table)
	{
		this.rates = rates;
		this.table = table;
		this.tree = rates.getTree();
		this.max_family_size = table.getMaxFamilySizes(tree);
		int root = tree.getRoot();
		factorials = new RisingFactorial(1.0, 1+max_family_size[root]);
		computeParameters();
	}
	
	final RateModel.GLD rates;
	final ProfileTable table;
	final IndexedTree tree;
	final RisingFactorial factorials;
	final int[] max_family_size;

	RateModel getRates(){ return rates;}
	
	
	// common data structures 
	/**
	 * Extinction probability in the subtree below the node.
	 * Not used here but in extending classes. 
	 * 
	 */
	double[] extinction;
	/**
	 * 3 parameters per node
	 */
	private double[] survival_parameters;
	
	@Override
	public double getGainParameter(int node) {return survival_parameters[3*node+PARAMETER_GAIN];}
	
	@Override
	public double getLossParameter(int node) {return survival_parameters[3*node+PARAMETER_LOSS];}
	
	@Override
	public double getDuplicationParameter(int node) {return survival_parameters[3*node+PARAMETER_DUPLICATION];}
	
	/**
	 * Gain parameter for survival model. 
	 * 
	 * @param node
	 * @param g kappa (if Polya), or r (if Poisson) corrected for extinction
	 */
	private void setGainParameter(int node, double g) {survival_parameters[3*node+PARAMETER_GAIN]=g;}
	/**
	 * Loss parameter in survival model. 
	 * 
	 * @param node
	 * @param p extinction probability for a parental copy towards this lineage
	 */
	private void setLossParameter(int node, double p){survival_parameters[3*node+PARAMETER_LOSS]=p;}
	/**
	 * Duplication parameter for survival model. 
	 * 
	 * @param node
	 * @param q duplication parameter (for Polya) corrected for extinction
	 */
	private void setDuplicationParameter(int node, double q) {survival_parameters[3*node+PARAMETER_DUPLICATION]=q;}
	
	/**
	 * Precomputed rising factorials for Polya pmf. 
	 */
	RisingFactorial[] gain_factorials;
	
	/**
	 * Allocates the arrays for the parameters data structures. 
	 */
	private void initDataStructures()
	{
		int num_nodes = max_family_size.length;
		this.extinction = new double[num_nodes];
		this.survival_parameters = new double[3*num_nodes];
		this.gain_factorials = new RisingFactorial[num_nodes];
	}
	
	/**
	 * Sum of profile log-likelihoods across the tables. 
	 * 
	 * @return
	 */
	public double getLL()
	{
		double LL = 0.0;
		for (int f=0; f<table.getFamilyCount(); f++)
		{
			Profile P = new Profile(f);
			double Pll = P.getLogLikelihood();
			LL += Pll;
			
			for (int node=0; node<tree.getNumNodes(); node++)
			{
				System.out.println("#*L.getLL family "+f+"\tnode "+node+"\tC[] = "+Arrays.toString(P.getNodeLikelihoods(node))+"\t// p "+getLossParameter(node)+"\tq "+getDuplicationParameter(node)+"\tkappa "+getGainParameter(node));
				System.out.println("#*L.getLL family "+f+"\tedge "+node+"\tK[] = "+Arrays.toString(P.getEdgeLikelihoods(node)));
			}
			
		}
		return LL;
	}
	
	/**
	 * Corrected log-likelihood: conditioned on non-empty profiles. 
	 * 
	 * @return
	 */
	public double getCorrectedLL()
	{
		double LL = getLL();
		int F = table.getFamilyCount();
		double L0 = getEmptyLL();
		// LL-F*log(1-exp(L0))
		double x = -Math.expm1(L0);
		LL -= F*Math.log(x); 
		return LL;
	}
	
	/**
	 * Log-likelihood of the empty profile. 
	 * 
	 * @return
	 */
	public double getEmptyLL()
	{
		double LL = 0.0;
		for (int node=0; node<tree.getNumNodes(); node++)
		{
			double q = getDuplicationParameter(node);
			if (q==0.0)
			{
				// Poisson
				double r = getGainParameter(node);
				LL -= r;
			} else
			{
				// Pólya
				double κ = getGainParameter(node);
				LL += κ*Math.log1p(-q); // log(1-q)
			}
		}
		return LL;
	}
	
	/**
	 * Call if parameters change in the underlying rate model
	 * (automatically called at instantiation). Calculates the  
	 * {@link #survival_parameters} and the {@link #extinction}
	 */
	public void computeParameters()
	{
		initDataStructures();

		// compute extinction probabilities and survival parameters
		for (int node=0; node<tree.getNumNodes(); node++) // including root, where the prior is set
		{
			double p = rates.getLossParameter(node);
			double q = rates.getDuplicationParameter(node);
			double epsi; // exinction probability at node 
			if (tree.isLeaf(node))
			{
				epsi = 0.0;
			} else
			{
				int num_children = tree.getNumChildren(node);
				int ci = 0;
				int child = tree.getChild(node, ci);
				epsi = getLossParameter(child);
				++ci;
				while (ci<num_children)
				{
					child = tree.getChild(node, ci);
					epsi *= getLossParameter(child);
					++ci;					
				}
			}
			
			if (q==0.0)
			{
				double r = rates.getGainParameter(node);
				r *= (1.-epsi); // Poisson
				this.gain_factorials[node] = null;
//				System.out.println("#*L.cP "+node+"\tq "+q+"\tpoisson "+r);
				setGainParameter(node, r);
			}
			else 
			{
				double κ = rates.getGainParameter(node);
				if (κ!=0.0)
					this.gain_factorials[node] = new RisingFactorial(κ, 1+max_family_size[node]);
				setGainParameter(node, κ);
			}
			double a = 1.-q*epsi;
			q *= (1.-epsi)/a; // correction for survival
			p += (1.-p)*epsi*(1.-q); // OK with q=0.0
						
			this.extinction[node] = epsi;
			setLossParameter(node, p);
			setDuplicationParameter(node, q);
		}
	}
    
    Profile getProfileLikelihood(int family_idx)
    {
    	return new Profile(family_idx);
    }
    
	/**
	 * Class for storing profile-specific conditional likelihoods.
	 * 
	 * @author csuros
	 *
	 */
	class Profile
	{
		private Profile(int family_idx)
		{
			this.family_idx = family_idx;
			int num_nodes = tree.getNumNodes();
			this.node_likelihoods = new double[num_nodes][];
			this.edge_likelihoods = new double[num_nodes][];
			initDataStructures();
		}
		private final int family_idx;
		private final double[][] node_likelihoods;
		private final double[][] edge_likelihoods;
		
		/**
		 * Allocates the row arrays in {@link #node_likelihoods} to match the sum
		 * of copy numbers at the leaves in the subtree. 
		 */
		private void initDataStructures()
		{
			int[] profile = table.getFamilyProfile(family_idx);
			for (int node=0; node<node_likelihoods.length; node++)
			{
				if (tree.isLeaf(node))
				{
					int n = profile[node];
					if (n<0) // ambiguous
					{
						node_likelihoods[node]=new double[0];
					} else
					{
						node_likelihoods[node] = new double[n+1];
					}
				} else
				{
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
					if (ambi == num_children)
						node_likelihoods[node] = new double[0];
					else 
						node_likelihoods[node] = new double[n+1];
				}
			}
		}
		
		/**
		 * Computes the conservation log-likelihoods. Automatically called 
		 * by {@link #getEdgeconservationLL(int)} and {@link #getNodeLikelihoods(int)} 
		 * if they are not yet initalized. 
		 * 
		 * Call if  parameters change. 
		 */
		public void computeLikelihoods()
		{
			if (USE_CLASSIC_RECURSION) // old Count's way
			{
				int[] postOrder = TreeTraversal.postOrder(tree);
				for (int node:postOrder) // makes no difference but has the same execution trace
				{
					computeNode(node);
					computeEdge(node);
				}
			} else // current favorite for calculations
			{
				for (int node=0; node<tree.getNumNodes(); node++)
				{
					computeNode(node);
					computeEdge(node);
				}
			}
		}
	
		public double getLogLikelihood()
		{
			int root = tree.getRoot();
			double[] K = getEdgeLikelihoods(root);
			double LL = K[0];
			
			return LL;
		}
		
		
		/**
		 * Conditional log-likelihoods per ancestral surviving copies. 
		 * 
		 * @param node
		 * @return
		 */
		double[] getNodeLikelihoods(int node)
		{
			double[] C = node_likelihoods[node];
			if (C==null)
			{
				computeLikelihoods();
				C = node_likelihoods[node]; // now it's there
			}
			return C;
		}
		
		/**
		 * Conditional log-likelihoods per conserved copies from parent. 
		 * 
		 * @param node
		 * @return
		 */
		double[] getEdgeLikelihoods(int node)
		{
			double[] K = edge_likelihoods[node];
			if (K==null)
			{
				computeLikelihoods();
				K = edge_likelihoods[node];
			}
			return K;
		}
		
		/**
		 * Computes a row (= 1 node) for {@link #edge_likelihoods}
		 * using {@link #node_likelihoods}.
		 * 
		 * @param node must be called in postorder 
		 * @return array of edge (survival) log-likelihoods
		 */
		private void computeEdge(int node)
		{
			double[] K; // the edge likelihoods that will be set 
			double[] C = node_likelihoods[node];
			
			if (tree.isRoot(node))
			{
				K = new double[1];
				if (C.length==0)
				{	// ambiguous at the root? tsk tsk tsk
					edge_likelihoods[node] = K; // bc likelihood is 1
					return;
				}
			}
			else
				K=new double[C.length];
			
			double q = getDuplicationParameter(node);
			double g = getGainParameter(node);
			
			if (q==0.0) // test before the loop
			{
				// Poisson distribution with parameter r 
				double r = g;
				if (r == 0.0)
				{
					// no gain, no duplication
					for (int s=0; s<K.length; s++)
					{
						K[s] = C[s];
					}
				} else
				{
					double logr = Math.log(r);
					double[] terms = new double[C.length];
					for (int s=0; s<K.length; s++)
					{
						
						int t=0, ell=s; 				
						while (ell<C.length)
						{
							assert (ell == s+t); // loop invariant
							terms[t] = C[ell] -r + t*logr -factorials.factln(t) ;
							++t; ++ell;
						}
						K[s] = Logarithms.sum(terms, t);
					}				
				}
			} else
			{
				// Polya distribution with parameters kappa+s and q 
				double κ = g;
				double logq = Math.log(q);
				double log1_q = Math.log1p(-q);
				double[] terms =new double[C.length];
				if (κ == 0.0)
				{
					// no gain: negative binomial with s and q 
					K[0] = C[0]; // no conservation here
					for (int s=1; s<K.length; s++)
					{
						double factln_s = factorials.factln(s-1); // s-1 and not ell-1
						int t = 0, ell=s+t;
						while (ell < C.length)
						{
							double binom = factorials.factln(ell-1) // ell-1 and not ell
											-factln_s-factorials.factln(t);
							terms[t] = C[ell] +  binom + s*log1_q + t*logq;
							++t;
							++ell;
						}
						K[s] = Logarithms.sum(terms, t);
					}
				} else
				{
					for (int s=0; s<K.length; s++)
					{
						double factln_s = gain_factorials[node].factln(s);
						int t=0, ell=s+t;
						while (ell<C.length)
						{
							assert (ell == s+t); // loop invariant;
							double binom = gain_factorials[node].factln(ell)-factln_s-factorials.factln(t);
							double w = binom + (κ+s)*log1_q + t*logq;
							terms[t] = C[ell] + w  ;
							++t;
							++ell;
						}
						K[s]=Logarithms.sum(terms,t);
					}
				}
			}
			edge_likelihoods[node] = K;
		}
		
		/**
		 * Algorithm for computing the node likelihood, and 
		 * filling the values in {@link #node_likelihoods} at one row. The actual 
		 * calculations fare done by 
		 * {@link #computeSibling(int, double[], double, int)}, adding one child at a time
		 * to an increasing set of siblings.  
		 * 
		 * @param node called in postorder, after {@link #computeEdge(int)}
		 */
		private void computeNode(int node)
		{
			int n = node_likelihoods[node].length;
			if (n==0) return; // all ambiguous
			
			if (tree.isLeaf(node))
			{
				double[] C = node_likelihoods[node];
				Arrays.fill(C, Double.NEGATIVE_INFINITY);
				C[n-1]=0.0;
			} else // ancestral node
			{
				int num_children = tree.getNumChildren(node);
				// loop over the siblings
				double[] C=null;// will be reused
				double sib_extinct= 1.0; // extinction probability product across siblings
				for (int c=0; c<num_children; c++)
				{
					int junior = tree.getChild(node, c); // current sibling
					if (USE_CLASSIC_RECURSION) // old Count
					{
						C = computeSiblingClassic(junior, C, sib_extinct);
					} else
						C = computeSibling(junior, C, sib_extinct);
					sib_extinct *= getLossParameter(junior);
				}
				assert (C.length == node_likelihoods[node].length);
				//System.arraycopy(C, 0, node_likelihoods[node], 0, n+1); // no need to copy
				node_likelihoods[node]=C;
			}
		}
		
		/**
		 * Algorithm for adding one child to a set of siblings. Called from {@link #computeNode(int)}.
		 * (Summing by surviving copies in this lineage.)
		 * 
		 * @param junior new child to be added 
		 * @param C siblings' (survival) log-likelihoods (will be destroyed)
		 * @param eps siblings' extinction 
		 * @return survival log-likelihoods for junior+siblings 
		 */
		double[] computeSibling(int junior, double[] C, double eps)
		{
			double[] K = edge_likelihoods[junior];
			if (C==null || C.length==0)
			{
				// ambiguous on the right
				return K.clone();
			} else if (K.length==0)
			{
				// ambiguous on the left 
				// nothing to do C is OK as is
				return C;
			}
			int combined_family_size = (K.length-1)+(C.length-1);
			
			// numerical parameters for the formulas
			double p = getLossParameter(junior);
			
			double loga = Math.log1p(-p*eps); // log(1-pe)
			double logp1 = Math.log1p(-p)-loga; // log((1-p)/(1-pe))
			double logp2 = Math.log(p)+Math.log1p(-eps)-loga;// log(p*(1-e)/(1-pe));
			
			double loge = Math.log(eps);
			double log1_e = Math.log1p(-eps); // log(1-e)
			
			double[] C2 = new double[combined_family_size+1]; // return value
			double[] terms = new double[K.length]; // auxiliary array to store summing terms 
			for (int ell=0; ell<=combined_family_size; ell++)
			{
				// calculate C(ell, t): t is at most as large as the max copy number in C[]  
				int t = Math.min(ell,  C.length);
				double x= // keeps the value of K2(ell,t+1) in the loop
					(t==C.length?Double.NEGATIVE_INFINITY:C[t]); // = C(ell,ell) while ell<C.length
				// when t starts with t==ell, C(t, t) = C[t] still, does not get touched
				// when t starts with t==C.length, all those values are log(0)
				int s= ell-t;
				while (t>0 && s<K.length-1)
				{
					--t;
					++s;
					x += log1_e;
					double y = C[t]+loge; // at this point C[t] = C(ell-1, t)
					x = C[t] = Logarithms.add(x, y); // recursion: C(ell,t) = (1-e)*C(ell,t+1)+e*C(ell-1,t)
				}
				
				// calculate C2[ell]
				double log_ellfact = factorials.factln(ell);

				assert (t==ell-s);
				int tmin = t; 
				while (s>=0 && t<C.length) // t<=ell by loop design 
				{
					double binom = log_ellfact - factorials.factln(s)-factorials.factln(t); // ell chose s 
					terms[t-tmin] // fill lower cells
							= K[s] + C[t] + binom +  s*logp1 + t*logp2;
					--s;
					++t;
				} 
				C2[ell] = Logarithms.sum(terms, t-tmin);  
			}
			return C2;
		}
				
		/**
		 * Algorithm in old Count (summing by surviving copies in previous siblings). 
		 * 
		 * @param sister
		 * @param C
		 * @param eps
		 * @return
		 */
		private double[] computeSiblingClassic(int sister, double[] C, double eps)
		{
			double p = getLossParameter(sister);
			double log1_p = Math.log1p(-p); // log(1-p)
			double logp = Math.log(p);
			double[] D = edge_likelihoods[sister].clone();
			for (int s=0; s<D.length; s++)
				D[s] += s*log1_p; // these are valid now with t=0, ell=s --- but flawed if unary parent with this single child 

			if (C == null || C.length==0 || D.length==0) // first child or ambiguous 
				return D;
			
			int combined_size = (D.length-1) + (C.length-1); // for sister added to the set already in C[]
			double[] C2 = new double[combined_size+1]; // return value
			double[] terms = new double[D.length]; // auxiliary array reused in the loop on ell for summing; min(Dlen,Clen) capacity suffices
			
			// precomputed multipliers
			double loge = Math.log(eps);
//			double log1_e = Math.log1p(-eps); // log(1-eps)
			double log_a = Math.log1p(-p*eps); // log(1-pe)
			
			// loop over ell for filling the return entries C2[ell]; inner loops with s+t=ell 
			for (int ell=0; ell<C2.length; ell++)
			{
//				Arrays.fill(terms, Double.NEGATIVE_INFINITY); // log(0.0) for safety
				
				// D[ell] is set if ell < D.length
				int s = Math.min(ell, D.length); // maximum+1
				// initialize t, and keep synchronized (s+t=ell) during the loops over s
				int t = ell-s; 

				double y ; // saves D(ell, s+1) in the loop on s 
				if (s==D.length)
					y = Double.NEGATIVE_INFINITY;
				else
					y = D[s]; // don't touch it: D(s,s) is OK for t=0 
				
				// compute D(ell,s) for all lower values of s 
				while (s>0 && t<C.length-1) // bound on t bc we won't need smaller s anymore as ell increases
				{
					--s;
					++t;
					// recursion: D(ell,s) = D(ell,s+1)+p*D(ell-1,s)
					y = D[s] = Logarithms.add(y, D[s]+logp); 
				}
				
				int smin = s; 
				// t is at its maximum possible value 
				double factln_ell = factorials.factln(ell); // same for all s,t
				while (t>=0 && s<D.length)
				{
					double binom = factln_ell-factorials.factln(s)-factorials.factln(t);
					double w = binom + s*loge - ell*log_a; //+t*log1_e;
					double z = D[s] + C[t] + w; 
					terms[s-smin] // use only the lower cells in terms[] 
							= z;
					--t;
					++s;
				}
				C2[ell] = Logarithms.sum(terms, s-smin) ;
			}
			return C2;
		}
	} // Profile class
	
	
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
        if (USE_CLASSIC_RECURSION)
        {
            out.println(Executable.getStandardHeader("(using classic recursion for node likelihoods)"));
        	
        }
        
        Likelihood factory = new Likelihood(rates, table);
        double LL = factory.getLL();
        
        
        for (int node=0; node<tree.getNumNodes(); node++)
        {
        }
        
        
        
        double L0 = factory.getEmptyLL();
        double p0 = Math.exp(L0);
        double corrL = LL-table.getFamilyCount()*Math.log1p(-p0);
        out.println("Log-likelihood:      \t"+LL+"\t("+Math.exp(LL)+")");
        out.println("Empty log-likelihood:\t"+L0+"\t("+Math.exp(L0)+")");
        out.println("Corrected log-lik.:  \t"+corrL+"\t("+Math.exp(corrL)+")");
 	}
	
	
	public static void main(String[] args) throws Exception
	{
		Likelihood testL = new Likelihood();
		testL.mainmain(args);
	}

}
