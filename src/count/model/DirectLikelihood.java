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
import java.util.Arrays;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.io.CommandLine;
import count.matek.Logarithms;
import count.matek.Functions.RisingFactorial;


/**
 * Likelihood computation by direct application of 
 * Felsenstein's peeling method, conditioning on ancestor copies. 
 * 
 * Use {@link #setCalculationWidthThresholds(int, double)} 
 * for settng the maximum ancestor copies assumed. 
 */
public class DirectLikelihood implements GLDParameters
{
	private static final boolean TRACK_CALCULATIONS = false;
	
	public DirectLikelihood(TreeWithRates rates, ProfileTable table)
	{
		this(rates, table, false);
	}

	
	protected DirectLikelihood(TreeWithRates rates, ProfileTable table, boolean want_transition_posteriors)
	{
		this.rates = rates;
		this.table = table;
		this.tree = rates.getTree();
		this.max_family_size = table.getMaxFamilySizes(tree);
		this.transition_probs = new TransitionProbabilities[tree.getNumNodes()];
		this.want_copy_number_transitions = want_transition_posteriors;
		//		int root = tree.getRoot();
//		int m = getCalculationWidth(max_family_size[root]);
//		factorials = new RisingFactorial(1.0, 1+m); // capacity for precomputing all necessary values 
		computeParameters();
	}
	
	final TreeWithRates rates;
	final ProfileTable table;
	final IndexedTree tree;
	RisingFactorial factorials;
	final int[] max_family_size;
	
	
	/**
	 * Whether we need to precompute {@link TransitionProbabilities#copychange(int, int)}  in a 
	 * cubic-time procedure; if false, only {@link TransitionProbabilities#copybirth(int, int)} and {@link TransitionProbabilities#copysurvival(int, int)}
	 * can be used which are precomputed in quadratic-time .
	 */
	private final boolean want_copy_number_transitions;
	final TransitionProbabilities[] transition_probs;
	
	private int threshold_width_absolute = 3; //Integer.MAX_VALUE;
	private double threshold_width_relative = 3.0; //Double.POSITIVE_INFINITY;
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
		
//		int m = getCalculationWidth(max_family_size[tree.getRoot()]);
//		factorials = new RisingFactorial(1.0, 1+m); 	
		
//		System.out.println("#**DL.sCWT "+absolute+","+relative);
		
		this.computeParameters();
	}	
	
	protected void setSameCalculationWidthThresholds(DirectLikelihood factory)
	{
		this.setCalculationWidthThresholds(factory.threshold_width_absolute, factory.threshold_width_relative);
	}
	
	protected int getCalculationWidthAbsolute()
	{
		return this.threshold_width_absolute; 
	}
	
	protected double getCalculationWidthRelative()
	{
		return this.threshold_width_relative;
	}
	
	private int getCalculationWidth(int max_value)
	{
		double rel = Math.ceil(threshold_width_relative*Math.sqrt(max_value+1.0));
		double width =  max_value+rel+threshold_width_absolute;
		int cwidth = (int) Double.min(width,Integer.MAX_VALUE);
		return cwidth;
	}
	
	/**
	 * 3 parameters per node
	 */
	private double[] node_parameters;
	private double[] node_complements;
	
	private double[] log_parameters;
	private double[] log_complements;
	
	/**
	 * Precomputed rising factorials for Polya pmf. Calculation 
	 * of {@link RisingFactorial#factln(int)} is 
	 * correct even if too few members are allocated + precomputed, 
	 * but it's slower.  
	 */
	private RisingFactorial[] gain_factorials;
	
	
	/**
	 * Allocates the arrays for the parameters data structures. 
	 */
	private void initDataStructures()
	{
		int num_nodes = tree.getNumNodes();
		this.node_parameters = new double[3*num_nodes];
		this.node_complements = node_parameters.clone();
		this.log_parameters = node_parameters.clone();
		this.log_complements = node_parameters.clone();
		this.gain_factorials = new RisingFactorial[num_nodes];
		int m = getCalculationWidth(max_family_size[tree.getRoot()]);
		factorials = new RisingFactorial(1.0, 1+m); 		
	}
	
	
	public void computeParameters()
	{
		
		this.initDataStructures();
		for (int node=0; node<tree.getNumNodes(); node++)
			computeParameters(node);
		
	}
	
	public void computeParameters(int node)
	{
		double p = rates.getLossParameter(node);
		double q = rates.getDuplicationParameter(node);
		double p_1 = rates.getLossParameterComplement(node);
		double q_1 = rates.getDuplicationParameterComplement(node);
		setLossParameter(node, p);
		setLossComplement(node, p_1);
		setDuplicationParameter(node, q);
		setDuplicationComplement(node, q_1);
		
		double r = rates.getGainParameter(node);
		setGainParameter(node, r); // after duplication parameter is set
		
		this.transition_probs[node] = new TransitionProbabilities(node);
	}
	

	TreeWithRates getRates(){ return rates;}
	@Override
	public double getGainParameter(int node) {return node_parameters[3*node+PARAMETER_GAIN];}
	
	@Override
	public double getLossParameter(int node) {return node_parameters[3*node+PARAMETER_LOSS];}
	
	@Override
	public double getDuplicationParameter(int node) {return node_parameters[3*node+PARAMETER_DUPLICATION];}
	
	@Override
	public double getDuplicationParameterComplement(int node) { return node_complements[3*node+PARAMETER_DUPLICATION];}

	@Override
	public double getLossParameterComplement(int node) { return node_complements[3*node+PARAMETER_LOSS];}
	
	private double getLogLossParameter(int node) { return log_parameters[3*node+PARAMETER_LOSS];}
	
	private double getLogLossComplement(int node) { return log_complements[3*node+PARAMETER_LOSS];}
	
	private double getLogDuplicationParameter(int node) {return log_parameters[3*node+PARAMETER_DUPLICATION];}
	
	private double getLogDuplicationComplement(int node) {return log_complements[3*node+PARAMETER_DUPLICATION];}
	
	private double getLogGainParameter(int node) {return log_parameters[3*node+PARAMETER_GAIN];}
	/**
	 * Gain parameter (locally stored); recomputes the gain_factorials.
	 * 
	 * @param node
	 * @param g kappa (if Polya), or r (if Poisson) 
	 */
	protected void setGainParameter(int node, double g) 
	{
		int j = 3*node+PARAMETER_GAIN;
		node_parameters[j]=g;
		log_parameters[j]=Math.log(g);

		// and precalculate rising factorials for this kappa, if Polya 
		double q = getDuplicationParameter(node);
		if (q==0.0) // Poisson
		{
			this.gain_factorials[node] = null;
//			System.out.println("#**DL.cP "+node+"\tq=0"+"\tr "+r);
		} else if (g!=0.0) // Polya distribution
		{
			int m = max_family_size[tree.isLeaf(node)?node:tree.getRoot()]; // logic of allocation for node_likelihoods[]
			int cw = getCalculationWidth(m); 
			this.gain_factorials[node] = new RisingFactorial(g, 1+cw);	
		}
	}
	/**
	 * Loss parameter . 
	 * 
	 * @param node
	 * @param p loss probability for a parental copy towards this lineage
	 */
	protected void setLossParameter(int node, double p)
	{
		node_parameters[3*node+PARAMETER_LOSS]=p;
		log_parameters[3*node+PARAMETER_LOSS] = Math.log(p);
	}
	/**
	 * Duplication parameter . 
	 * 
	 * @param node
	 * @param q duplication parameter (for Polya) 
	 */
	protected void setDuplicationParameter(int node, double q) 
	{
		int j = 3*node+PARAMETER_DUPLICATION;
		node_parameters[j]=q;
		log_parameters[j]=Math.log(q);
	}


	protected void setLossComplement(int node, double one_minus_p) 
	{
		int j = 3*node+PARAMETER_LOSS;
		node_complements[j]=one_minus_p;
		log_complements[j] = Math.log(one_minus_p);
	}

	
	protected void setDuplicationComplement(int node, double one_minus_q) 
	{
		int j = 3*node+PARAMETER_DUPLICATION;
		node_complements[j]=one_minus_q;
		log_complements[j]=Math.log(one_minus_q);
	}
	
	/**
	 * Inverts the stored distribution parameters 
	 * to set the rates and edge lengths in the underlying 
	 * {@link #rates} model. Gives better precision 
	 * than {@link TreeWithRates#setParameters(int, double, double, double)}
	 * because we know the complements (which matters when p or q are 
	 * close to 1.0).
	 * 
	 * @param node
	 */
	protected void setRatesForParameters(int node)
	{
		double p = getLossParameter(node);
		double q = getDuplicationParameter(node);
		
		if (p==1.0)
		{
			double μ = rates.getLossRate(node);
			assert (μ!=0.0);
			assert (Double.isFinite(μ));
			assert (Double.isFinite(q));
			rates.setDuplicationRate(node, q*μ);
		} else // ordinary node 
		{
			double λ, μ;
			
			if (p==0.0)
			{
				μ = 0.0;
				λ = -getLogDuplicationComplement(node);
			} else if (q==0.0) // Poisson
			{
				μ = -getLogLossComplement(node);	
				λ = 0.0;
			} else // Pólya
			{
				double p_q = p-q;  // unavoidable subtraction
				
				if (p_q==0.0)
				{
					double p_1 = getLossParameterComplement(node);
					μ = λ = p/p_1;
				} else if (0<p_q) // q<p
				{
					double p_1 = getLossParameterComplement(node);
					double gap_t = Math.log1p(p_q/p_1); // log1p works better for small p_q
					μ = p/p_q*gap_t;
					λ = q/p_q*gap_t;
				} else // p<q
				{
					double q_1 = getDuplicationParameterComplement(node);
					double q_p = -p_q;
					assert (0<q_p); // argument to log1p is positive
					double gap_t = Math.log1p(q_p/q_1); // log1p works better for small q_p
					μ = p/q_p*gap_t;
					λ = q/q_p*gap_t;
				}
			}
			rates.setTotalRates(node, μ, λ);
		}
		
		// now set gain
		double r = getGainParameter(node);
		if (q==0.0 && p!=0.0) // Poisson
		{
			rates.setGainRate(node, r/p);
		} else
		{
			rates.setGainRate(node, r);
		}
//		{ // DEBUG
//		System.out.println("#**DL.sRP done "+node
//			+"\t"+toString(node)
//			+"\t// "+rates.toString(node));
//		}
		// and synchronize for precision
		this.computeParameters(node);
	}
	
	
	/**
	 * 
	 * @param outside outside log-likelihoods
	 * @param inside inside log-likelihoods
	 * @return array of posterior probabilities
	 */
	protected static double[] computePosteriors(double[] outside, double[] inside)
	{
		double[] posteriors = outside.clone();
		for (int ell=0; ell<posteriors.length; ell++)
		{
			posteriors[ell]+=(ell<inside.length?inside[ell]:Double.NEGATIVE_INFINITY);
		}
		double LL = Logarithms.sum(posteriors, posteriors.length);
		for (int ell=0; ell<posteriors.length; ell++)
		{
			posteriors[ell]=Math.exp(posteriors[ell]-LL);
		}
		return posteriors;
	}
	
	
	
	
	
	
	public Profile getProfile(int family)
	{
		return new Profile(family);
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
			this.edge_outside = new double[num_nodes][];
			this.node_outside = new double[num_nodes][];
			this.sibling_likelihoods = new double[num_nodes][];
			this.initDataStructures();
		}
		final int family_idx;
		private final double[][] node_likelihoods;
		private final double[][] edge_likelihoods;
		private final double[][] sibling_likelihoods;
		private final double[][] node_outside;
		private final double[][] edge_outside;
		
		
		private void initDataStructures()
		{
			int max_n = table.maxCopies(family_idx);
			int m = DirectLikelihood.this.getCalculationWidth(max_n);
			int[] profile = table.getFamilyProfile(family_idx);
			
			for (int  node=0; node<node_likelihoods.length; node++)
			{
				if (tree.isLeaf(node))
				{
					node_likelihoods[node] = new double[profile[node]+1];
				} else
				{
					node_likelihoods[node] = new double[m+1];
				}
				edge_likelihoods[node] = null;
				sibling_likelihoods[node] = null;
				node_outside[node] = null;
				edge_outside[node] = null;
			}
		}
		
		public void computeInside()
		{

			for (int node=0; node<tree.getNumNodes(); node++)
			{
				computeNodeInside(node);
				edge_likelihoods[node] = computeEdgeInside(node);
			}
		}
		
		/**
		 * Call after computeInside
		 */
		public void computeOutside()
		{
			int root = tree.getRoot();
			for (int node=root; node>=0; --node)
			{
				edge_outside[node] = computeEdgeOutside(node);
				node_outside[node] = computeNodeOutside(node);
			}
		}
		
		/**
		 * Allocates and computes {@link #edge_likelihoods}
		 * 
		 * @param node
		 * @return proper value for <code>edge_likelihoods[node]</code>
		 */
		private double[] computeEdgeInside(int node)
		{
			double[] K; // the edge likelihoods that will be set 
			double[] C = node_likelihoods[node];
			
			if (tree.isRoot(node))
			{
				double p = getLossParameter(node);
				if (p==1.0) 
					K = new double[1];
				else 
					K = new double[2];
			} else
			{
				K=new double[C.length];
			}
			double q = getDuplicationParameter(node);
			double r = getGainParameter(node);
			TransitionProbabilities Tn = transition_probs[node];
			
			if (q==0.0) // test before the loop
			{
				// Poisson distribution with parameter r ; untested
				if (r == 0.0)
				{
					// no gain, no duplication
					for (int s=0; s<K.length; s++)
					{
						K[s] = C[s];
					}
				} else
				{		
//					double logr = getLogGainParameter(node); // Math.log(r);
					double[] terms = new double[C.length];
					for (int s=0; s<K.length; s++)
					{
						int t=0, ell=s+t; 				
						while (ell<C.length)
						{
							assert (ell == s+t); // loop invariant
//							double t_logr = t*logr;
//							terms[t] = C[ell] -r + t_logr -factorials.factln(t) ;
							terms[t] = Tn.copybirth(s, ell) + C[ell];
							
							++t; ++ell;
						}
						K[s] = Logarithms.sum(terms, t);
					}
				}
			} else
			{
				// Polya distribution with parameters kappa+s and q 
				double κ = r;
//				double logq = getLogDuplicationParameter(node ); // Math.log(q);
//				double log1_q = getLogDuplicationComplement(node); //  Math.log(one_minus_q);
				double[] terms =new double[C.length];
				if (κ == 0.0)
				{
					// no gain: negative binomial with s and q 
					K[0] = C[0]; // no conservation here
					for (int s=1; s<K.length; s++)
					{
//						double factln_s = factorials.factln(s-1); // s-1 and not ell-1
						int t = 0, ell=s+t;
						while (ell < C.length)
						{
//							double binom = factorials.factln(ell-1) // ell-1 and not ell
//											-factln_s-factorials.factln(t);
//							double t_logq = t*logq;
//							terms[t] = C[ell] +  binom + s*log1_q + t_logq;
							terms[t] = Tn.copybirth(s, ell) + C[ell];
							++t;
							++ell;
						}
						K[s] = Logarithms.sum(terms, t);
					}
				} else
				{
//					if (gain_factorials[node]==null)
//					{ // DEBUG
//						System.out.println("#**DL.P.cEI f "+family_idx+"\tnode "+node+"\tq "+q+"\tkappa "+κ);
//					}
					for (int s=0; s<K.length; s++)
					{
//						double factln_s = gain_factorials[node].factln(s);
						int t=0, ell=s+t;
						while (ell<C.length)
						{
							assert (ell == s+t); // loop invariant;
//							double binom = gain_factorials[node].factln(ell)-factln_s-factorials.factln(t);
//							double t_logq = t*logq;
//							double w = binom + (κ+s)*log1_q + t_logq;
//							
//							double check = transition_probs[node].copybirth(s, ell);
//							if (check!=w)
//							{
//								System.out.println("#**DL.P.cEI "+family_idx+"\t"+node+"\ts "+s+"\tell "+ell+"\there "+w+"\ttrans "+check);
//								throw new RuntimeException();
//							}
//							
//							terms[t] = C[ell] + w  ;
							terms[t] = Tn.copybirth(s, ell) + C[ell]; // K[s] = sum for ell=s+t, t=0.. <C.length-s}; s fix 
							
							
							++t;
							++ell;
						}
						K[s]=Logarithms.sum(terms,t);
					}
					if (TRACK_CALCULATIONS)
						System.out.println("#**DL.P.CEI "+family_idx+"\t"+node+"\t("+K.length+")\t"+Arrays.toString(K));
						
				}
			}
			return K;
		}
		
		/**
		 * Computes {@link #node_likelihoods}; 
		 * allocates and computes {@link #sibling_likelihoods}
		 * for a node. 
		 * 
		 * @param node
		 */
		private void computeNodeInside(int node)
		{
			double[] C = node_likelihoods[node]; 
			int n=C.length;

			if (tree.isLeaf(node))
			{
				Arrays.fill(C, Double.NEGATIVE_INFINITY);
				C[n-1] = Math.log(1.0);
				
//				System.out.println("#**L.P.cN "+node+"\tleaf n="+n+"\tC "+Arrays.toString(C));
			} else // ancestral node
			{
				Arrays.fill(C, 0.0);
				
				// compute stem likelihoods 
				int num_children = tree.getNumChildren(node);
				
				double[][] sib = new double[num_children][];
				for (int ci=0; ci<num_children; ci++)
				{
					int child = tree.getChild(node, ci);
					TransitionProbabilities Tc = DirectLikelihood.this.transition_probs[child]; 

					double[] C2 = new double[C.length];
					
//					double logp = getLogLossParameter(child);
//					double log1_p = getLogLossComplement(child);
					
					double[] K = edge_likelihoods[child];
					double[] terms = new double[K.length]; // temporary, reused
					
					for (int ell=0; ell<C2.length; ell++)
					{
						int k=0; 
//						int d = ell-k; // 
//						double log_ellfact = factorials.factln(ell);
						while(k<K.length && k<=ell)
						{
//							double klog1_p = k==0?0.0:k*log1_p;
//							double dlogp = d==0?0.0:d*logp;
//							double binom = log_ellfact - factorials.factln(k)-factorials.factln(d);
//							terms[k] =  binom + klog1_p + dlogp; 
//							
//							double check = transition_probs[child].copysurvival(ell, k);
//							if (check!=terms[k])
//							{
//								System.out.println("#**DL.P.cNI "+family_idx+"\t"+child+"\tell "+ell+"\tk "+k+"\there "+terms[k]+"\ttrans "+check);
//								throw new RuntimeException();
//							}
							terms[k] = Tc.copysurvival(ell, k) + K[k]; // cNI C2[ell] = sum k=0..ell <K.length; ell fix 
							
							
//							d--;
							k++;
						}
						double z = Logarithms.sum(terms, k);
						C2[ell] = z;
						C[ell] += z; // product across children
					}
					sib[ci] = C2; // save it for later
				}
				if (TRACK_CALCULATIONS)
					System.out.println("#**DL.P.CNI "+family_idx+"\t"+node+"\t("+C.length+")\t"+Arrays.toString(C)+"\t// "+tree.toString(node)+"\t// "+DirectLikelihood.this.toString(node));
				
				// compute product of stem likelihoods across siblings 
				// -- in linear and not quadratic time 
				//  with arity  
				
				double[] lprod = new double[num_children];
				double[] rprod = new double[num_children];
				for (int ell=0; ell<C.length; ell++)
				{
					lprod[0]=0.0;
					int ci=0; 
					while (ci<num_children-1)
					{
						lprod[ci+1]= lprod[ci]+sib[ci][ell];
						++ci;
					}
					rprod[ci] = 0.0;
					while (0<ci)
					{
						rprod[ci-1] = rprod[ci]+sib[ci][ell];
						--ci;
					}
					while (ci<num_children)
					{
						sib[ci][ell] = lprod[ci]+rprod[ci];
						++ci;
					}
				}
				for (int ci=0; ci<num_children; ci++)
				{
					int child = tree.getChild(node, ci);
					sibling_likelihoods[child] = sib[ci]; // product of all siblings except child 
					if (TRACK_CALCULATIONS)
					{
						
						System.out.println("#**DL.P.CNI "+family_idx+"\tsib "+child+"\t"+Arrays.toString(sib[ci])+"\t// "+tree.toString(node));
					}

				}
			}
		}
		
		/**
		 * Allocates and computes {@link #node_outside}.
		 * 
		 * @param node
		 * @return
		 */
		private double[] computeNodeOutside(int node)
		{
			double[] C = node_likelihoods[node];
			double[] J = edge_outside[node]; // already set
			double[] B = new double[C.length]; // return value
			
			
			double q = getDuplicationParameter(node);
			TransitionProbabilities Tn = transition_probs[node];
			if (q==0.0)
			{ // untested
				// Poisson
//				double r = getGainParameter(node);
//				double logr = getLogGainParameter(node);
				
				double[] terms = new double[J.length];
				for (int ell=0; ell<B.length; ell++)
				{
					int t=ell,s=0;
					
					while (s<J.length && 0<=t) // s<=ell
					{
//						double tlogr = t==0?0.0:t*logr;
//						terms[s] = J[s] - r + tlogr - factorials.factln(t);
						terms[s] = Tn.copybirth(s, ell) + J[s]; // 
						++s;
						--t;
					}
					B[ell] = Logarithms.sum(terms, s);
				}
			}  else
			{
				// Pólya
				double κ = getGainParameter(node);
//				double logq = getLogDuplicationParameter(node);
//				double log1_q =getLogDuplicationComplement(node);
				double[] terms =new double[J.length];
				
				if (κ==0.0)
				{
					
					// negative binomial with s and q 
					B[0]=J[0]; // no conservation
					
					for (int ell=1; ell<B.length; ell++)
					{
//						double factln_ell = factorials.factln(ell-1); // ell-1 and not ell
						
						
						int t=ell-1, sm1=0, s=1; // no contribution from s=0
						while (s<J.length) // sm1==s-1; t+s=ell
						{
//							double binom = factln_ell - gain_factorials[node].factln(sm1) - factorials.factln(t); // s-1 and not s
//							terms[sm1] = binom + s*log1_q + t*logq;
//							double check = transition_probs[node].copybirth(s, ell);
//							if (check!=terms[sm1])
//							{
//								System.out.println("#**DL.P.cNO "+family_idx+"\t"+node+"\ts "+s+"\tell "+ell+"\there "+terms[sm1]+"\ttrans "+check);
//								throw new RuntimeException();
//							}
//							terms[sm1]+= J[s];
							terms[sm1] =  J[s]+Tn.copybirth(s, ell);
							
							
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
//						double log_ellfact = gain_factorials[node].factln(ell);
						int s=0; // t=ell 
						while (s<J.length && s<=ell) // 0<=t
						{
//							double binom = log_ellfact - gain_factorials[node].factln(s) - factorials.factln(t);
//							double tlogq = t*logq;
//							
//							terms[s] = J[s] + binom + (κ+s)*log1_q + tlogq;
							terms[s] = J[s] + Tn.copybirth(s, ell); // B[ell] = sum s=0..ell, < J.length fix ell
							
							++s;
//							--t;
						}
						B[ell] = Logarithms.sum(terms, s);
					}

				}
			}
			if (TRACK_CALCULATIONS)
				System.out.println("#**DL.P.CNO "+family_idx+"\t"+node+"\t("+B.length+")\t"+Arrays.toString(B)+"\tC "+Arrays.toString(C)+"\t// "+tree.toString(node));
			return B;
		}

		/**
		 * Allocates and computes {@link #edge_outside}.
		 * 
		 * @param node
		 * @return
		 */
		private double[] computeEdgeOutside(int node)
		{
			double[] J; // return value
			double p = getLossParameter(node);
			TransitionProbabilities Tn = transition_probs[node]; 
			
			if (tree.isRoot(node))
			{
				if (p == 1.0) // most likely, in standard models
				{
					J = new double[1];
					// and the value is log(1)
				} else {
					J = new double[2];
					J[0] = getLogLossParameter(node);
					J[1] = getLogLossComplement(node);
				}
			} else
			{ // descending from a parent, with siblings to consider 
				int parent = tree.getParent(node);
				double[] B = node_outside[parent];
				
				int n = node_likelihoods[node].length;
				J = new double[n];
				
//				double log1_p=getLogLossComplement(node);
//				double logp = getLogLossParameter(node);
				
				double[] terms = new double[B.length];
				
				for (int s=0; s<J.length; s++)
				{
//					double slog1_p =s==0?0.0:s*log1_p;
//					double logsfact = factorials.factln(s);
					
					int t = 0, ell = s+t;
					
					while (ell<B.length)
					{
//						double tlogp = t==0?0.0:t*logp;
//						double binom = factorials.factln(ell)-logsfact-factorials.factln(t);
//
//						double w = binom + slog1_p + tlogp ;
//						double check = transition_probs[node].copysurvival(ell, s);
//						if (check!=w)
//						{
//							System.out.println("#**DL.P.cEO "+family_idx+"\t"+node+"\ts "+s+"\tell "+ell+"\there "+w+"\ttrans "+check);
//							throw new RuntimeException();
//						}
//						
//						terms[t] = B[ell] + w + sibling_likelihoods[node][ell];
						terms[t] = B[ell] + sibling_likelihoods[node][ell] + Tn.copysurvival(ell, s) ; // cEO J[s] = sum ell = s.. < B.length; s fix

						t++;
						ell++;
					}
					J[s] = Logarithms.sum(terms, t);
				}
			}
			if (TRACK_CALCULATIONS)
				System.out.println("#**DL.P.CEO "+family_idx+"\t"+node+"\t("+J.length+")\t"+Arrays.toString(J)+"\tK "+Arrays.toString(getEdgeInside(node))+"\t// "+tree.toString(node));
			return J;
		}

	
		protected double[] getNodeInside(int node)
		{
			if (edge_likelihoods[node]==null) // check if computed them yet 
			{
				computeInside();
			}
			double[] C = node_likelihoods[node];
			return C;
		}
		
		protected double[] getEdgeInside(int node)
		{
			double[] K = edge_likelihoods[node];
			if (K==null)
			{
				computeInside();
				K = edge_likelihoods[node];
			}
			return K;
		}
		
		protected double[] getNodeOutside(int node)
		{
			double[] B = node_outside[node];
			if (B==null) 
			{
				computeOutside();
				B = node_outside[node];
			}
			return B;
			
		}
		
		protected double[] getEdgeOutside(int node)
		{
			double[] J = edge_outside[node];
			if (J==null)
			{
				computeOutside();
				J = edge_outside[node];
			}
			return J;
		}
		
		/**
		 * Posterior distribution of ancestor's copy count. 
		 * 
		 * @param node
		 * @return
		 */
		public double[] getNodePosteriors(int node)
		{
			double[] C = getNodeInside(node);
			double[] B = getNodeOutside(node); // sets up inside likelihoods if necessary
//			System.out.println("#**P.P.gNP "+node+"\tB "+Arrays.toString(B)+"\tC "+Arrays.toString(C));
			
			return computePosteriors(B, C);
		}
		
		/**
		 * Posterior distribution of ancestor's inherited copy count. 
		 * 
		 * @param node
		 * @return
		 */
		public double[] getEdgePosteriors(int node)
		{
			double[] K = getEdgeInside(node);
			double[] J = getEdgeOutside(node);
			
			return computePosteriors(J, K);
		}		
				
		/** 
		 * Likelihood of this profile
		 * @return
		 */
		public double getLogLikelihood()
		{
			int root = tree.getRoot();
			double[] K = getEdgeInside(root);
			
			double LL = K[0];
			double proot = getLossParameter(root);
			if (proot != 1.0 && K.length==2)
			{
				double logp = getLogLossParameter(root);
				LL += logp;
				//double p1root = DirectLikelihood.this.getLossParameterComplement(root);
				double log1_p = getLogLossComplement(root);
				LL = Logarithms.add(LL, K[1]+log1_p);
			}
			
			return LL;
		}
		
		/**
		 * Incorrect calculation that ignores the sibling's contribution.  
		 * 
		 * 
		 * @param node
		 * @return
		 */
		public double[][] getTransitionPosteriors(int node)
		{
			double[][] post; // return value
			if (tree.isRoot(node))
			{
				post = new double[1][];
				post[0] = getNodePosteriors(node);
			} else
			{
				int parent = tree.getParent(node);
				
				double[] B = getNodeOutside(parent);
				double[] C = getNodeInside(node);
				TransitionProbabilities T = transition_probs[node];
				double[] sib = sibling_likelihoods[node];
				
				post = new double[B.length][];
				double LL = Double.NEGATIVE_INFINITY; // log(0.0)
				
				double[] pC = getNodeInside(node);
				for (int n=0; n<post.length; n++)
				{
					double[] post_n = new double[C.length];
					double tsum = Double.NEGATIVE_INFINITY;
					
					for (int m=0; m<post_n.length; m++)
					{
						double term = T.copychange(n, m) +C[m] + sib[n];
						post_n[m] = term;
						tsum = Logarithms.add(tsum,  term);
					}
//					System.out.println("#**DL.P.gTP\t"+family_idx+"\tnode "+node
//							+"\tn "+n+"\ttsum "+tsum+"\tCn"+pC[n]);
					for (int m=0; m<post_n.length; m++)
						post_n[m] += B[n];
					tsum += B[n];
					
					
					
					post[n] = post_n;
					// we could normalize here 
					// using tsum to store m|n conditionals instead of joint n,m
					LL = Logarithms.add(LL, tsum);
				}
//				if (TRACK_CALCULATIONS)
//				System.out.println("#**DL.P.gTP\t"+family_idx+"\tnode "+node+"\tLL "+LL+"\treal "+getLogLikelihood());
				for (int n=0; n<post.length; n++)
				{
					double[] post_n = post[n];
					for (int m=0; m<post_n.length; m++)
					{
						post_n[m] = Math.exp(post_n[m]-LL);
					}
				}
				
//				double[] pP = getNodePosteriors(parent);
//				double[] C = getNodeInside(node);
//				TransitionProbabilities T = transition_probs[node];
//				post = new double[pP.length][C.length];
//				
//				// want B_n*C_m*sib_n*w[m|n]
//				for (int n=0; n<post.length; n++)
//				{
//					double pn = pP[n];
//					double[] post_n = post[n];
//					
//					// normalizer to get conditional probabilities (child = m| parent=n)
//					double log_sum = Double.NEGATIVE_INFINITY; // log(0.0)
//					
//					for (int m=0; m<C.length; m++)
//					{
//						double term = T.logw(n,m)+C[m]; // log-weight of n->m
//						log_sum = Logarithms.add(log_sum, term);
//						post_n[m] = term;
//					}
//					for (int m=0; m<C.length; m++)
//					{
//						post_n[m] = Math.exp(post_n[m]-log_sum)*pn;
//					}
//				}
				
			}
			
			return post;
		}
		
	
		
	} // Profile class
	
	public String toString(int node)
	{
		String node_name = tree.getIdent(node);

		StringBuilder sb = new StringBuilder(node_name);
		double p = getLossParameter(node);
		double q = getDuplicationParameter(node);
		double r = getGainParameter(node);
		sb.append("[p ").append(p)
			.append("/1-").append(getLossParameterComplement(node)).append(")")
			.append("; q ").append(q)
			.append("/1-").append(getDuplicationParameterComplement(node)).append(")")
			.append("; r ").append(r);
		
		if (q<p)
		{
			double n = q*r/(p-q);
			sb.append("; n ").append(n);
		} else
		{
			sb.append("; n ").append(Double.POSITIVE_INFINITY);
		}
		return sb.append("]").toString();
	}
	
	
	
	
	private static double mean(double[] p)
	{
		double m = 0.0;
		for (int i=1; i<p.length; i++)
			m += i*p[i];
		return m;
	}
	
	private class TransitionProbabilities
	{
		TransitionProbabilities(int node)
		{
			this.node = node;
			int nmax = DirectLikelihood.this.max_family_size[tree.isLeaf(node)?node:tree.getRoot()];
			nmax = getCalculationWidth(nmax);
			int pmax;
			if (tree.isRoot(node))
			{
				pmax = 0;
			} else
			{
				pmax = max_family_size[tree.getRoot()];
			}
			pmax = getCalculationWidth(pmax);
			this.log_w = new double[1+pmax][1+nmax];
			this.fromparent = new double[1+pmax][];
			this.tochild = new double[1+pmax][1+nmax];
			
//			System.out.println("#**DL.TP "+node+"\tpmax "+pmax+"\tnmax "+nmax);
			this.calculateValues();
		}
		
		private final int node;
		private double[][] log_w;
		
		private double[][] fromparent;
		private double[][] tochild;
		
		/**
		 * Transition probabilities are calculated by a cubic-time 
		 * algorithm (loops n,m,s) instead of the 
		 * quadratic-time recurrences bc the latter 
		 * involve subtraction when p+q &gt; 1.
		 */
		private void calculateValues()
		{
			double q = getDuplicationParameter(node);
			double r = getGainParameter(node);
			double logq = getLogDuplicationParameter(node);
			double loq1_q = getLogDuplicationComplement(node);
			double logr  = getLogGainParameter(node);
			

			int n = 0; // 0->m transitions

			double[] log_wn = log_w[n];
			double[] fromparent_n = fromparent[n] =new double[1];
			double[] tochild_n = tochild[0];
			fromparent_n[0]=Math.log(1.0);
			
			if (q == 0.0)
			{
				// Poisson
				// prob(0->m) = exp(-r)*r^m/m!
				for (int m=0; m<log_wn.length; m++)
				{
					double m_logr = m==0?0.0:m*logr;
					double term = -r + m_logr - factorials.factln(m);
					log_wn[m] = term;
					tochild_n[m] = term; 
				}
			} else
			{
				// Polya
				
				double rlog1_q = r*loq1_q;
				
				for (int m=0; m<log_wn.length; m++)
				{
					// n-> m transitions
					double binom = gain_factorials[node].factln(m)-factorials.factln(m);
					double m_logq = m==0.0?0.0:m*logq;
					double term = binom + rlog1_q + m_logq;
					log_wn[m] = term;
					tochild_n[m] = term;
				}
			}
			++n;
			// 
			if (!tree.isRoot(node))
			{
				//double p = getLossParameter(node);
				double logp = getLogLossParameter(node);
				double log1_p = getLogLossComplement(node);
				while(n<log_w.length)
				{
					log_wn = log_w[n];
					fromparent_n = fromparent[n] = new double[n+1];
					for (int k=0,d=n-k; k<fromparent_n.length; k++, d--)
					{
						double loss_binom = factorials.factln(n)-factorials.factln(k)-factorials.factln(d);
						double loss_term = loss_binom+k*log1_p+d*logp;
						fromparent_n[k] = loss_term;
					}
					{
						int k=n;
						Arrays.fill(tochild[k], 0, Integer.min(k, tochild[k].length), Double.NEGATIVE_INFINITY); // tochild[k][m] = log(0.0) when m<k 
						for (int j=0, m=k+j; m<log_wn.length; m++, j++)
						{
							double dup_term;
							if (q==0.0)
							{
								// Poisson
								// prob = exp(-r)*r^k/k!
								dup_term = -r + j*logr-factorials.factln(j);
							} else
							{
								// Polya
								double dup_binom = gain_factorials[node].factln(m)-gain_factorials[node].factln(k)
										-factorials.factln(j);
								dup_term = dup_binom + (r+k)*loq1_q + j*logq;
							}
							tochild[k][m] = dup_term;
						}
					}
					
					if (want_copy_number_transitions)
					{
						double[] terms = new double[log_wn.length];
						for (int m=0; m<log_wn.length; m++)
						{
							int k=0;
							int t=n-k;
							int j=m-k;
	//						double factln_m = gain_factorials[node].factln(m);
							while (0<=t && 0<=j)
							{
								//double loss_binom = factorials.factln(n)-factorials.factln(k)-factorials.factln(t);
								double loss_term = fromparent_n[k]; // loss_binom+k*log1_p+t*logp;
								double dup_term = tochild[k][m];
	//							if (q==0.0)
	//							{
	//								// Poisson
	//								// prob = exp(-r)*r^k/k!
	//								dup_term = -r + j*logr-factorials.factln(j);
	//							} else
	//							{
	//								// Polya
	//								double dup_binom = factln_m-gain_factorials[node].factln(k)
	//										-factorials.factln(j);
	//								dup_term = dup_binom + (r+k)*loq1_q + j*logq;
	//							}
								terms[k] = loss_term+dup_term;
	//							tochild[k][m] = dup_term;
								
								k++;
								t--;
								j--;
							}
							log_wn[m] = Logarithms.sum(terms, k);
						}
						log_w[n] = log_wn; // already same 
					}
					++n;
				}
			}
//			checkValues();
		}
		
		private void checkValues() // debug code
		{
			double[][] w = new double[log_w.length][];
			
			double p = getLossParameter(node);
			double p_1 = getLossParameterComplement(node);
			double q = getDuplicationParameter(node);
			double pq_1 = p_1-q; 
			System.out.println("#**DL.TP.checkV node "+node+"\tp "+p+"/1-"+p_1+"\tq "+q+"\t1-p-q "+pq_1);

			int n = 0;
			double[] log_wn = log_w[n];
			double[] w0 = new double[log_wn.length];
			for (int m=0; m<w0.length; m++)
			{
				w0[m] = Math.exp(log_wn[m]);
				System.out.println("#**DL.TP.checkV node "+node+"\tn "+n+"\tm "+m+"\twn "+w0[m]+"/"+Math.exp(log_wn[m])+"\tlog "+Math.log(w0[m])+"/"+log_wn[m]);
			}
			w[0]=w0;
			++n;
			while (n<w.length)
			{
				log_wn = log_w[n];
				double[] wn = new double[log_wn.length];
				int m=0;
				wn[m] = p*w[n-1][m];
				System.out.println("#**DL.TP.checkV node "+node+"\tn "+n+"\tm "+m+"\twn "+wn[m]+"/"+Math.exp(log_wn[m])+"\tlog "+Math.log(wn[m])+"/"+log_wn[m]);
				++m;
				while (m<wn.length)
				{
					wn[m] = p*w[n-1][m]+ q*wn[m-1] +pq_1*w[n-1][m-1];
					
					System.out.println("#**DL.TP.checkV node "+node+"\tn "+n+"\tm "+m+"\twn "+wn[m]+"/"+Math.exp(log_wn[m])+"\tlog "+Math.log(wn[m])+"/"+log_wn[m]);
					++m;
				}
				w[n] = wn;
				++n;
			}
		}
		
		/**
		 * Log-probability of n-to-m transition
		 * 
		 * @param n parent copy number (condition)
		 * @param m child copy number 
		 * @return
		 */
		double copychange(int n, int m)
		{
			assert (want_copy_number_transitions);
			return log_w[n][m];
		}
		
		/**
		 * Log-probability of surviving copies.
		 * 
		 * @param n_parent
		 * @param n_surviving_copies
		 * @return
		 */
		double copysurvival(int n_parent, int n_surviving_copies)
		{
			if (fromparent.length <= n_parent 
					|| fromparent[n_parent] == null
					|| fromparent[n_parent].length <= n_surviving_copies)
			{
				System.out.println("#**DL.TP.cs "+node+"\tnp "+n_parent+"/"+fromparent.length
						+"\tns "+n_surviving_copies
							+"/"+(fromparent[n_parent]==null?0:fromparent[n_parent].length)
							+"\t// "+tree.toString(node)
						);
			}
			return this.fromparent[n_parent][n_surviving_copies]; 
		}
//		terms[t] = B[ell] + sibling_likelihoods[node][ell] + Tn.copysurvival(ell, s) ; // cEO J[s] = sum ell = s.. < B.length; s fix
//		terms[k] = Tc.copysurvival(ell, k) + K[k]; // cNI C2[ell] = sum k=0..ell <K.length; ell fix 
		
		/**
		 * Log-probability of m copies given that k survived from parent
		 * 
		 * @param k_surviving_copies
		 * @param m_child_copies
		 * @return
		 */
		double copybirth(int k_surviving_copies, int m_child_copies)
		{
			return this.tochild[k_surviving_copies][m_child_copies];
		}
//		terms[t] = Tn.copybirth(s, ell) + C[ell]; // CEI K[s] = sum for ell=s+t, t=0.. <C.length-s}; s fix 
//		terms[s] = J[s] + Tn.copybirth(s, ell); // CNO B[ell] = sum s=0..ell, < J.length fix ell
		
		
	}
	
	
	// --- tests and main
	
	private void testCompareLikelihoods(PrintStream out, Posteriors post)
	{
		
		int nFam = table.getFamilyCount();
        int nNodes = tree.getNumNodes();
        int nLeaves = tree.getNumLeaves();
		
        double totdLL = 0.0;
        double totpLL = 0.0;
        double totdiff = 0.0;
        
        for (int f=-1; f<nFam; f++)		
        {
        	double dLL, pLL;
        	int m=0, ncopies=0, npresent=0, cwidth; 

        	if (f==-1)
        	{
        		pLL = post.factory.getEmptyLL();
        		ProfileTable empty = ProfileTable.emptyProfile(tree);
        		Posteriors EP = new Posteriors(rates, empty);
        		DirectLikelihood ED = new DirectLikelihood(rates, empty);
        		ED.setCalculationWidthThresholds(threshold_width_absolute, threshold_width_relative);
        		Profile DP = ED.getProfile(0);
        		dLL = DP.getLogLikelihood();
        		cwidth = ED.getCalculationWidth(0);
        		
        	} else
        	{
        		Profile DP = this.getProfile(f);
        		Posteriors.Profile PP = post.getPosteriors(f);
        		dLL = DP.getLogLikelihood();
        		pLL = PP.inside.getLogLikelihood();
        		m = table.maxCopies(f);
        		ncopies = table.getMemberCount(f);
        		npresent = table.getLineageCount(f);
        		cwidth = getCalculationWidth(m);
        	}
        	
            double diff = pLL-dLL;
            assert (diff >= 0.0);
            
            out.println("#FAMILY\t"+f
            		+"\t"+dLL+"\t"+pLL+"\t"+diff+"\t"+Math.abs(diff/pLL)
            		+"\t//\t"+m+"\t("+npresent+","+ncopies+")\t"+cwidth
            		);
            
            if (0<=f)
            {
	            totdLL += dLL;
	            totpLL += pLL;
	            totdiff += diff;
            }
        }
        
        out.println("#FAMILY\ttotals\t"+totdLL+"\t"+totpLL
        		+"\t"+totdiff+"\t"+Math.abs(totdiff/totpLL));
	}
	
	private void testComparePosteriors(PrintStream out, DirectLikelihood.Profile DP, Posteriors.Profile PP)
	{
        int nNodes = tree.getNumNodes();
        int nLeaves = tree.getNumLeaves();
		
        int family_idx = PP.inside.family_idx;
                
    	for (int node=nLeaves; node<nNodes; node++ )
    	{
    		double[] post_pN = PP.getNodeAncestorPosteriors(node);
    		double[] direct_pN = DP.getNodePosteriors(node);
    		
    		double post_mN = mean(post_pN);
    		double direct_mN = mean(direct_pN);
    		
    		out.print(family_idx+"\tN\t"+node+"\t"+direct_mN+"/"+post_mN);
    		
    		for (int n=0; n<post_pN.length || n<direct_pN.length; n++)
    		{
    			double pp = n<post_pN.length?post_pN[n]:0.0;
    			double dp = n<direct_pN.length?direct_pN[n]:0.0;
    			
    			out.print("\t"+dp+"/"+pp);
    		}
    		out.println();

    		double[] post_pS = PP.getEdgeAncestorPosteriors(node);
    		double[] direct_pS = DP.getEdgePosteriors(node);
    		
    		double post_mS = mean(post_pS);
    		double direct_mS = mean(direct_pS);
    		
    		out.print(family_idx+"\tS\t"+node+"\t"+direct_mS+"/"+post_mS);
    		for (int s=0; s<post_pS.length || s<direct_pS.length; s++)
    		{
    			double pp = s<post_pS.length?post_pS[s]:0.0;
    			double dp = s<direct_pS.length?direct_pS[s]:0.0;
    			
    			out.print("\t"+dp+"/"+pp);
    		}
    		{
        		double[][] trans = DP.getTransitionPosteriors(node);
            	double sum = 0.0;
	    		for (int n=0; n<trans.length; n++)
	    			for (int m=0; m<trans.length; m++)
	    			{
	    				out.printf("\t%d:%d %g", n,m,trans[n][m]);
	    				sum += trans[n][m];
	    			}
    		}
    		out.println();
    	}
	}
	
	private void testComparePosteriors(PrintStream out, Posteriors post)
	{
		
		
		int nFam = table.getFamilyCount();
		
        for (int f=0; f<nFam; f++)		
        {
        	Profile DP = new Profile(f);
        	Posteriors.Profile PP = post.getPosteriors(f);
        	testComparePosteriors(out, DP, PP);
        	
        }
		
		
	}
	
	public static void main(String[] args) throws Exception
	{
		CommandLine cli = new CommandLine(args, DirectLikelihood.class);
        AnnotatedTable table = cli.getTable();
        TreeWithRates rates = cli.getRates();

        PrintStream out = System.out;        
        
        DirectLikelihood DL = new DirectLikelihood(rates, table);
        
    	int absolute = 3;
    	double relative = 1.0;
        if (cli.getOptionValue(CommandLine.OPT_TRUNCATE)!=null)
        {
        	absolute = cli.getOptionTruncateAbsolute();
        	relative = cli.getOptionTruncateRelative();
		}
		out.println(CommandLine.getStandardHeader("Truncated computation (absolute,relative)="
				+absolute+","+relative));
		
		DL.setCalculationWidthThresholds(absolute, relative);
		
		Posteriors post = new Posteriors(rates, table);
		post.setCalculationWidthThresholds(12, 3.0);
		
//		DL.testCompareLikelihoods(out, post);
		DL.testComparePosteriors(out, post);
		
	}
	
}
