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

import static count.io.CommandLine.OPT_MINCOPY;

import java.io.PrintStream;
import java.util.Arrays;



import count.matek.Functions.RisingFactorial;
import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.ds.TreeTraversal;
import count.ds.UniqueProfileTable;
import count.io.CommandLine;
import count.matek.Logarithms;

/**
 * Algorithms for computing conditional likelihoods bottom-up.  
 * 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */
public class Likelihood implements GLDParameters
{
	private static final boolean USE_CLASSIC_RECURSION = false; // in Count 12
	private static final boolean ASSERT_CALCULATIONS = false; // debug numerical errors 
	
	
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
	
	public Likelihood(TreeWithRates rates, ProfileTable table)
	{
		this.rates = rates;
		this.table = table;
		this.tree = rates.getTree();
		this.max_family_size = table.getMaxFamilySizes(tree);
		int root = tree.getRoot();
		factorials = new RisingFactorial(1.0, 1+max_family_size[root]); // capacity for precomputing all necessary values 
		computeParameters();
	}
	
	public Likelihood(TreeWithRates rates)
	{
		this(rates, ProfileTable.emptyTable(rates.getTree()));
	}
	
	
	final TreeWithRates rates;
	final ProfileTable table;
	final IndexedTree tree;
	final RisingFactorial factorials;
	final int[] max_family_size;
	

	TreeWithRates getRates(){ return rates;}
	
	
	// common data structures 
	/**
	 * Extinction probability in the subtree below the node.
	 * Not used here but in extending classes. 
	 * 
	 */
	double[] extinction;
	double[] extinction_complement;
	
	/**
	 * 3 parameters per node
	 */
	private double[] survival_parameters;
	private double[] survival_complements;
	
	@Override
	public double getGainParameter(int node) {return survival_parameters[3*node+PARAMETER_GAIN];}
	
	@Override
	public double getLossParameter(int node) {return survival_parameters[3*node+PARAMETER_LOSS];}
	
	@Override
	public double getLossParameterComplement(int node) { return survival_complements[3*node+PARAMETER_LOSS];}
	
	@Override
	public double getDuplicationParameter(int node) {return survival_parameters[3*node+PARAMETER_DUPLICATION];}
	
	@Override
	public double getDuplicationParameterComplement(int node) { return survival_complements[3*node+PARAMETER_DUPLICATION];}
	
	public double getExtinction(int node) { return extinction[node];}

	public double getExtinctionComplement(int node)
	{
		double epsi = extinction[node];
		return epsi == 1.0?extinction_complement[node]:1.0-epsi;
	}
	
	private int threshold_width_absolute = Integer.MAX_VALUE;
	private double threshold_width_relative = Double.POSITIVE_INFINITY;
	
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
	}
	protected void setSameCalculationWidthThresholds(Likelihood factory)
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
	
	
	private void setLossComplement(int node, double p_1)
	{
		survival_complements[3*node+PARAMETER_LOSS]=p_1;
	}
	private void setDuplicationComplement(int node, double q_1)
	{
		survival_complements[3*node+PARAMETER_DUPLICATION]=q_1;
	}
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
		this.survival_complements = new double[survival_parameters.length];
		this.gain_factorials = new RisingFactorial[num_nodes];
		this.extinction_complement = new double[num_nodes];
	}
	
	/**
	 * Sum of profile log-likelihoods across the tables. 
	 * (Single-thread execution). 
	 * 
	 * @return
	 */
	public double getLL()
	{
		if (saved_log_likelihood < 0.0)
			return saved_log_likelihood;
		
		double LL = 0.0;
		for (int f=0; f<table.getFamilyCount(); f++)
		{
			Profile P = getProfileLikelihood(f);
			double Pll = P.getLogLikelihood();
			LL += Pll;
			
//			for (int node=0; node<tree.getNumNodes(); node++)
//			{
//				System.out.println("#*L.getLL family "+f+"\tnode "+node+"\tC[] = "+Arrays.toString(P.getNodeLikelihoods(node))+"\t// p "+getLossParameter(node)+"\tq "+getDuplicationParameter(node)+"\tkappa "+getGainParameter(node));
//				System.out.println("#*L.getLL family "+f+"\tedge "+node+"\tK[] = "+Arrays.toString(P.getEdgeLikelihoods(node)));
//			}
			
		}
		this.saved_log_likelihood = LL;
		return LL;
	}
	
	
	private double saved_log_likelihood = 0.0;
	
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
		double p0 = Math.exp(L0);
		// LL-F*log(1-exp(L0))
		//double x = 1.0-Math.exp(L0);
		double p_not0 = -Math.expm1(L0);
		double sub = F*Math.log(p_not0); 
//		System.out.println("#*L.gcLL "+LL+"\tL0 "+L0+"\tsub "+sub+"\tnot0 "+p_not0+"\tp0 "+p0);
		LL -= sub;
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
			} else if (q==1.0)
			{
				return Double.NEGATIVE_INFINITY;
			} else
			{
				// Pólya
				double κ = getGainParameter(node);
				double q_1 = getDuplicationParameterComplement(node);
				LL += κ* Math.log(q_1);  // Math.log1p(-q); //  Math.log(1.0-q); // log(1-q)
			}
		}
		double proot = getLossParameter(tree.getRoot());
		LL +=  Math.log(proot); // OK if proot=1 

		
		{ // DEBUG
			if (Double.isNaN(LL) || LL>0.0)
			{
				double sum = 0.0;
				for (int node=0; node<tree.getNumNodes(); node++)
				{
					double q = getDuplicationParameter(node);
					double q_1 = getDuplicationParameterComplement(node);
					double term;
					if (q==0.0)
					{
						// Poisson
						double r = getGainParameter(node);
						term = -r;
					} else if (q==1.0)
					{
						term = Double.NEGATIVE_INFINITY;
					} else
					{
						double κ = getGainParameter(node);
						term = κ* Math.log(q_1);  //Math.log1p(-q);
					}
					sum += term;
					System.out.println("#**L.gEL numerr node "+node+"\tsum "+sum+"\tterm "+term+"\tq "+q+"\t==1 "+(q==1.0)+"\tlog1p(-q) "+Math.log1p(-q)+"\t1_q "+q_1+"\tlog(q_1) "+Math.log(q_1)+"\t// rates "+rates.toString(node)+"\t// lik "+this.toString(node));
				}			
				System.out.println("#**L.gEL numerr "+LL+"\t("+sum+")\tproot "+proot+"/logp "+Math.log(proot));
			}
		}
			
		return LL;
	}
	
	private ProfileTable singletons = null;
	
	/**
	 * Log-probability of a profile with a single gene at a single leaf.  
	 * 
	 * @return
	 */
	public double getSingletonLL()
	{
		if (singletons == null)
			singletons = ProfileTable.singletonTable(tree);
		Likelihood S = new Likelihood(rates,singletons);
		double log_p1 = Double.NEGATIVE_INFINITY;
		for (int f=0; f<singletons.getFamilyCount(); f++)
		{
			Profile P = S.getProfileLikelihood(f);
			double log_pf = P.getLogLikelihood();
//			System.out.println("#**L.gSL leaf "+f+"\tLf "+log_pf+"\tpf "+Math.exp(log_pf)+"\tExp_nf "+Math.exp(log_pf)*table.getFamilyCount()+"\t// "+tree.toString(f));
			 
			log_p1 = Logarithms.add(log_p1, log_pf);
		 }
		 return log_p1;
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
			computeNodeParameters(node);
		}
		this.saved_log_likelihood = 0.0; // needs to be recalculated
	}
	
	public void computeNodeParameters(int node)
	{
		double epsi; // extinction probability at node 
		
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
		
		if (epsi>1.0)
		{
			System.out.println("#**L.cNP numerr epsi "+epsi+"\tnode "+node+"\t"+rates.toString(node) );
		}
		setNodeParameters(node, epsi);
	}
	
	private double computeExtinctionComplement(int node)
	{
		int num_children = tree.getNumChildren(node);
		
		double dp_ext0=1.0; // extinction in all children
		double dp_ext1=0.0; // non-extinction in at least one child
		
		for (int ci=0; ci<num_children; ci++)
		{
			int child = tree.getChild(node, ci);
			double p = getLossParameter(child);
			double p_1 = getLossParameterComplement(child); 
			
			dp_ext1 = dp_ext1 + dp_ext0*p_1;
			dp_ext0 = dp_ext0*p;
		}
		
		return dp_ext1;
	}
	
	
	/**
	 * Sets the node survival parameters for a given extinction probability. 
	 * @param node
	 * @param epsi
	 */
	void setNodeParameters(int node, double epsi)
	{
		this.extinction[node] = epsi;
		double epsi_1; // 1-epsilon: recalculate if epsilon is too close to 1.0
		if (epsi==1.0)
		{
			epsi_1 = this.extinction_complement[node] = 
				computeExtinctionComplement(node);
//			System.out.println("#**L.sNP "+node+"\tepsi "+epsi+"\te1 "+epsi_1);
		}
		else
			epsi_1 = this.extinction_complement[node] = 
				1.0-epsi;
		
		double p = rates.getLossParameter(node);
		

		double q = rates.getDuplicationParameter(node);

		if (ASSERT_CALCULATIONS) assert (!Double.isNaN(p));
		
		
		if (q==0.0)
		{
			double r = rates.getGainParameter(node);
			r *= epsi_1; // Poisson
			this.gain_factorials[node] = null;
//			System.out.println("#*L.cP "+node+"\tq "+q+"\tpoisson "+r);
			setGainParameter(node, r);
		}
		else 
		{
			double κ = rates.getGainParameter(node);
			if (κ!=0.0)
				this.gain_factorials[node] = new RisingFactorial(κ, 1+max_family_size[node]);
			setGainParameter(node, κ);
		}
		
		double p_1 = rates.getLossParameterComplement(node);
		double q_1 = rates.getDuplicationParameterComplement(node);
		
		// the numerical calculations avoid subtractions 1-x 
		// instead of a=1-q*e, we go by 
		// (1-q+q)*(1-e+e) = 1
		// (1-q)e+q(1-e)+(1-q)*(1-e) = 1-qe 
		// (1-q)e + (1-e)
		
		 // a = q_1*(1-epsi_1)+epsi_1 = q_1 -q_1*epsi_1+epsi_1 = q_1 + (1-q_1)*epsi_1
		
		double a = q_1*epsi+epsi_1; // 1.-q*epsi;
		
		double survq = q*epsi_1/a;
		// q *= epsi_1/a; // *=(1.-epsi)/a; // correction for survival
		// // a/epsi_1 =  q_1*epsi/epsi_1 + 1
		double survq_1 = q_1/a;
		// q_1 /= a;
		// // a/q_1 = epsi + epsi_1/q_1
		
		// double dp = 
		
		double survp = (p*epsi_1+epsi*q_1)/a; //p + p_1*epsi*survq_1;
		// p += p_1*epsi*q_1; 
		
		// p~ = (p+(1-p)*e*(1-q~) = (p(1-e)+e(1-q))/a
		// 1-p~ =( 1-qe -p+pe-e+eq)/a = ((1-p)*(1-e))/a = p_1*(epsi_1/a)

		// 1-p~ = 1-p-(1-p)e(1-q~) = (1-p)*(1-e(1-q~))
		// here we calculate 1-e(1-q~) = 1-e + q~e
		double b = epsi_1/a; //epsi_1 + survq*epsi; // 1.-epsi*q_1;
		double survp_1 = p_1*b;
		// p_1 *= b;
		
//		assert (survq<=1.0);
//		
//		if (survp>1.0)
//		{
//			System.out.println("#**L.sNP numerr p "+survp
//					+"\tnode "+node
//					+"\te "+epsi
//					+"\te1 "+epsi_1
//					+"\t"+rates.toString(node));
//			// double a = q_1*epsi+epsi_1;
//			 // a = q_1*(1-epsi_1)+epsi_1 = q_1 -q_1*epsi_1+epsi_1 = q_1 + (1-q_1)*epsi_1
//			// = q_1+epsi_1-q_1*epsi_1 
//			
//			double a2 = q_1+epsi_1-q_1*epsi_1;
//			double a3 = 1.0-q*epsi;
//			
//			// p~ = (p+(1-p)*e*(1-q~) = (p(1-e)+e(1-q))/a
//			double p2 = p + p_1*epsi*survq_1;
//			
//			// double b = epsi_1 + survq*epsi; 
//			// 1-p~ =( 1-qe -p+pe-e+eq)/a = ((1-p)*(1-e))/a = p_1*(epsi_1/a)
//			double b2 = epsi_1 + survq*epsi;
//			
//			System.out.println("#**L.sNP numerr p "+survp
//					+"\tnode "+node
//					+"\ta "+a+"\ta2 "+a2+"\ta3 "+a3
//					+"\tb "+b+"\tb2 "+b2
//					+"\tp1 "+survp_1
//					+"\tp2 "+p2
//					+"\t"+rates.toString(node));
//			
//		}
//		
//		assert (survp<=1.0);
//		
		if (ASSERT_CALCULATIONS) assert (!Double.isNaN(survp));
		setLossParameter(node, survp);
		setDuplicationParameter(node, survq);		
		setDuplicationComplement(node, survq_1);
		
		setLossComplement(node, survp_1);
		
//		if (epsi==1.0)
//			System.out.println("#**L.sNP "+node+"\tb "+b+"\t("+(1.-epsi*q_1)+")\tp1 "+p_1);
		
//		this.extinction[node] = epsi;
//		double epsi_1;
//		if (epsi==1.0)
//			epsi_1 = this.extinction_complement[node] = computeExtinctionComplement(node);
//		else
//			epsi_1 = this.extinction_complement[node] = 1.0-epsi;
//		
//		
//		double p = rates.getLossParameter(node);
//		
//
//		double q = rates.getDuplicationParameter(node);
//
//		if (ASSERT_CALCULATIONS) assert (!Double.isNaN(p));
//		
//		
//		if (q==0.0)
//		{
//			double r = rates.getGainParameter(node);
//			r *= epsi_1; // Poisson
//			this.gain_factorials[node] = null;
////			System.out.println("#*L.cP "+node+"\tq "+q+"\tpoisson "+r);
//			setGainParameter(node, r);
//		}
//		else 
//		{
//			double κ = rates.getGainParameter(node);
//			if (κ!=0.0)
//				this.gain_factorials[node] = new RisingFactorial(κ, 1+max_family_size[node]);
//			setGainParameter(node, κ);
//		}
//		
//		double p_1 = rates.getLossParameterComplement(node);
//		double q_1 = rates.getDuplicationParameterComplement(node);
//		
//		double a = 1.-q*epsi;
//		if (a==0.0)
//		{
//			// (1-q+q)*(1-e+e) = 1
//			// (1-q)e+q(1-e)+(1-q)*(1-e) = 1-qe 
//			a = q_1*epsi+epsi_1;
//			
//			// a/epsi_1 =  q_1*epsi/epsi_1 + 1
//			// a/q_1 = epsi + epsi_1/q_1
//		}
//		q *= epsi_1/a; // *=(1.-epsi)/a; // correction for survival
//		
//		
//		double dp = 
//		p += p_1*epsi*q_1; 
//		
//		if (ASSERT_CALCULATIONS) assert (!Double.isNaN(p));
//		setLossParameter(node, p);
//		q_1 /= a;
//		setDuplicationParameter(node, q);		
//		setDuplicationComplement(node, q_1);
//		setLossComplement(node, p_1*(1.-epsi*q_1));
//		
//		// 1-e(1-q) = aq

		if (ASSERT_CALCULATIONS && q==1.0 && q_1==.0)
		{
			System.out.println("#**L.sNP "+node+"\trq "+rates.getDuplicationParameter(node)
			+"\trq1 "+rates.getDuplicationParameterComplement(node)
			+"\teps "+epsi+"\t(1-e) "+(1.0-epsi)
			+"\ta "+a
			+"\tp " // +dp+"/"
			+p
			+"\tq "+q+"\tq1 "+q_1
			+"\t// rates "+rates.toString(node)
			+"\t// lik "+this.toString(node));
		}
			
			
		if (ASSERT_CALCULATIONS) assert (q!=1.0 || q_1!=0.0);

	}
	
	public String toString(int node)
	{
		String node_name = tree.getIdent(node);

		StringBuilder sb = new StringBuilder(node_name);
		double p = getLossParameter(node);
		double q = getDuplicationParameter(node);
		double r = getGainParameter(node);
		sb.append("[p~ ").append(p)
			.append(" (1-p~ ").append(getLossParameterComplement(node)).append(")")
			.append("; q~ ").append(q)
			.append(" (1-q~ ").append(getDuplicationParameterComplement(node)).append(")")
			.append("; r~ ").append(r)
			.append("; e ").append(extinction[node]);
		
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
		Profile(int family_idx)
		{
			this.family_idx = family_idx;
			int num_nodes = tree.getNumNodes();
			this.node_likelihoods = new double[num_nodes][];
			this.edge_likelihoods = new double[num_nodes][];
			initDataStructures();
		}
		final int family_idx;
		private final double[][] node_likelihoods;
		private final double[][] edge_likelihoods;
		
		/**
		 * Allocates the row arrays in {@link #node_likelihoods} to match the sum
		 * of copy numbers at the leaves in the subtree. 
		 */
		private void initDataStructures()
		{
			int[] profile = table.getFamilyProfile(family_idx);
			Arrays.fill(node_likelihoods, null);
			for (int node=0; node<node_likelihoods.length; node++)
			{
				int max_n = table.maxCopies(family_idx);
				double max_relative = 
						max_n==0?0.0:threshold_width_relative*max_n; // avoids problems with infinity
				assert (max_relative>=0.0);
				int relative = (max_relative>Integer.MAX_VALUE)?
						Integer.MAX_VALUE:(int)max_relative;
				int truncate_at = Integer.max(threshold_width_absolute, relative);
				if (tree.isLeaf(node))
				{
					int n = profile[node];
					if (n<0) // ambiguous
					{
						setCalculationWidth(node, 0);
						//node_likelihoods[node]=new double[0];
					} else
					{
						setCalculationWidth(node, Integer.min(n, truncate_at)+1);
						//node_likelihoods[node] = new double[n+1];
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
						setCalculationWidth(node, 0); 
						//node_likelihoods[node] = new double[0];
					else
						setCalculationWidth(node, Integer.min(truncate_at,n)+1);
						//node_likelihoods[node] = new double[n+1];
				}
			}
			
			
		}
		
		protected int computeCalculationWidth(int node)
		{
			int max_n = table.maxCopies(family_idx);
			double max_relative = 
					max_n==0?0.0:threshold_width_relative*max_n; // avoids problems with infinity
			assert (max_relative>=0.0);
			int relative = (max_relative>Integer.MAX_VALUE)?
					Integer.MAX_VALUE:(int)max_relative;
			int truncate_at = Integer.max(threshold_width_absolute, relative);
			
			if (tree.isLeaf(node))
			{
				int n = table.getFamilyProfile(family_idx)[node];
				if (n<0) // ambiguous
				{
					setCalculationWidth(node, 0);
					//node_likelihoods[node]=new double[0];
				} else
				{
					setCalculationWidth(node, Integer.min(n, truncate_at)+1);
					//node_likelihoods[node] = new double[n+1];
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
					setCalculationWidth(node, 0); 
					//node_likelihoods[node] = new double[0];
				else
					setCalculationWidth(node, Integer.min(truncate_at,n)+1);
					//node_likelihoods[node] = new double[n+1];
			}
			return getCalculationWidth(node);
		}			
		
		
		int getCalculationWidth(int node)
		{
			return node_likelihoods[node].length;
		}
		
		
		void setCalculationWidth(int node, int len)
		{
			node_likelihoods[node]=new double[len];
		}
		
		/** (Used with 
		 * AncestorPosteriors.)
		 * 
		 * @param that
		 */
		protected void copyLikelihoods(Profile that)
		{
			assert (this.node_likelihoods.length == that.node_likelihoods.length);
			for (int node=0; node<this.node_likelihoods.length; node++)
			{
				double[] thisC = this.node_likelihoods[node];
				double[] thatC = that.node_likelihoods[node];
				
				//assert (this.node_likelihoods[node].length>=that.node_likelihoods[node].length);
				
				double[] thatK = that.edge_likelihoods[node];
				double[] thisK = tree.isRoot(node)
						?new double[thatK.length]
						:new double[thisC.length];
				this.edge_likelihoods[node] = thisK;
				
				Arrays.fill(thisC, Double.NEGATIVE_INFINITY);
				Arrays.fill(thisK, Double.NEGATIVE_INFINITY);
				System.arraycopy(thatC, 0, thisC, 0, Integer.min(thatC.length, thisC.length));
				System.arraycopy(thatK, 0, thisK, 0, Integer.min(thatK.length, thisK.length));
			}
		}
		
		/**
		 * Computes the conservation log-likelihoods. Automatically called 
		 * by {@link #getEdgeLikelihoods(int)} and {@link #getNodeLikelihoods(int)} 
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
//					System.out.println("#**L.P.cL "+family_idx+"\t"+node
//							+"\tnode "+Arrays.toString(node_likelihoods[node])
//							+"\tedge "+Arrays.toString(node_likelihoods[node])
//							+"\t//r "+rates.toString(node)
//							+"\t// "+Likelihood.this.toString(node)
//							+"\t"+tree.toString(node));
				}
//				if (family_idx==0)
//					throw new RuntimeException();
			}
		}
	
		public double getLogLikelihood()
		{
			int root = tree.getRoot();
			double[] K = getEdgeLikelihoods(root);
			
			double LL = K[0];
			double proot = getLossParameter(root);
			if (proot != 1.0 && K.length==2)
			{
				double p1root = Likelihood.this.getLossParameterComplement(root);
				LL += Math.log(proot);
				LL = Logarithms.add(LL, K[1]+Math.log(p1root));
			}
			
//			if (Double.isNaN(LL))
//			{
//				System.out.println("#**L.P.gLL "+family_idx+"\tK "+Arrays.toString(K)+"\tpr "+proot+"\tprof "+Arrays.toString(table.getFamilyProfile(family_idx)));
//			}
			if (ASSERT_CALCULATIONS) assert (!Double.isNaN(LL));
			
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
		 * Allocates and computes a row for {@link #edge_likelihoods}
		 * using {@link #node_likelihoods}.
		 * 
		 * @param node must be called in postorder 
		 */
		private void computeEdge(int node)
		{
			double[] K; // the edge likelihoods that will be set 
			double[] C = node_likelihoods[node];
			
			int calc_width = getCalculationWidth(node);
			assert (C.length == calc_width);
			
			if (tree.isRoot(node))
			{
				double p = getLossParameter(node);
				if (p==1.0) 
					K = new double[1];
				else 
					K = new double[2];
				if (calc_width==0)
				{	// ambiguous at the root? tsk tsk tsk
					edge_likelihoods[node] = K; // bc likelihood is 1
					return;
				}
			}
			else
				K=new double[calc_width];
			
			double q = getDuplicationParameter(node);
			double one_minus_q = getDuplicationParameterComplement(node);
			if (ASSERT_CALCULATIONS) assert (q!=1.0 || one_minus_q != 0.0);
			double g = getGainParameter(node);
			
			assert Double.isFinite(g);
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
						if (ASSERT_CALCULATIONS) assert (!Double.isNaN(K[s]));
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
							double t_logr = (t==0?0.0:t*logr);
							terms[t] = C[ell] -r + t_logr -factorials.factln(t) ;
							++t; ++ell;
						}
						K[s] = Logarithms.sum(terms, t);
						if (ASSERT_CALCULATIONS) assert (!Double.isNaN(K[s]));
					}				
				}
			} else
			{
				// Polya distribution with parameters kappa+s and q 
				double κ = g;
				double logq = Math.log(q);
				double log1_q = Math.log(one_minus_q);
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
							double t_logq = (t==0?0.0:t*logq);
							terms[t] = C[ell] +  binom + s*log1_q + t_logq;
							++t;
							++ell;
						}
						K[s] = Logarithms.sum(terms, t);
						if (ASSERT_CALCULATIONS) assert (!Double.isNaN(K[s]));
					}
				} else
				{
					if (gain_factorials[node] == null)
						System.out.println("#**L.P.cE "+node+"\t"+tree.toString(node)+"\t"+rates.getDuplicationRate(node)+"\tq "+rates.getDuplicationParameter(node)+"\tq~ "+q);
					for (int s=0; s<K.length; s++)
					{
						double factln_s = gain_factorials[node].factln(s);
						int t=0, ell=s+t;
						while (ell<C.length)
						{
							assert (ell == s+t); // loop invariant;
							double binom = gain_factorials[node].factln(ell)-factln_s-factorials.factln(t);
							double t_logq = (t==0?0.0:t*logq);
							double w = binom + (κ+s)*log1_q + t_logq;
							if (ASSERT_CALCULATIONS) assert !Double.isNaN(binom);
							if (ASSERT_CALCULATIONS) assert !Double.isNaN(w);
							if (ASSERT_CALCULATIONS) assert !Double.isNaN(C[ell]);
							terms[t] = C[ell] + w  ;
							++t;
							++ell;
						}
						K[s]=Logarithms.sum(terms,t);
//						if (Double.isNaN(K[s]))
//						{
//							System.out.println("#**L.P.cE "+node+"/"+s+"\t"+tree.toString(node)+"\t"+rates.getDuplicationRate(node)+"\tq "+rates.getDuplicationParameter(node)+"\tq~ "+q+"\t"+Arrays.toString(terms)+"\tt="+t);
//						}
						if (ASSERT_CALCULATIONS) assert (!Double.isNaN(K[s]));
					}
				}
			}
			edge_likelihoods[node] = K;
		}
		
		/**
		 * Algorithm for computing the node likelihood, and 
		 * filling the values in {@link #node_likelihoods} at one row. The actual 
		 * calculations fare done by 
		 * {@link #computeSibling(int, double[], double)}, adding one child at a time
		 * to an increasing set of siblings.  
		 * 
		 * @param node called in postorder, after {@link #computeEdge(int)}
		 */
		private void computeNode(int node)
		{
			int n = getCalculationWidth(node);
			assert (n==node_likelihoods[node].length);
			
			if (n==0) return; // all ambiguous
			
			if (tree.isLeaf(node))
			{
				double[] C = node_likelihoods[node];
				Arrays.fill(C, Double.NEGATIVE_INFINITY);
				C[n-1]=0.0;
//				System.out.println("#**L.P.cN "+node+"\tleaf n="+n+"\tC "+Arrays.toString(C));
			} else // ancestral node
			{
				int num_children = tree.getNumChildren(node);
				// loop over the siblings
				double[] C=null;// will be reused
				if (extinction[node]==1.0) 					
				{ // deal with this before numerical complications
					C = new double[1];
					C[0]=0.0; // log(1.0)
					for (int c=0; c<num_children; c++)
					{
						int junior = tree.getChild(node, c); // current sibling
						double[] K = edge_likelihoods[junior];
						if (K.length>0)
							C[0] += K[0];
					}
				} else
				{
					double sib_extinct= 1.0; // extinction probability product across siblings
					for (int c=0; c<num_children; c++)
					{
						int junior = tree.getChild(node, c); // current sibling
	//					System.out.println("#**L.P.cN "+node+"\tjunior"+c+" "+junior+"\tC "+Arrays.toString(C)+"\tsibe "+sib_extinct);
						
						if (USE_CLASSIC_RECURSION) // old Count
						{
							C = computeSiblingClassic(junior, C, sib_extinct);
						} else
							C = computeSibling(junior, C, sib_extinct);
						
						double pj = getLossParameter(junior);
						sib_extinct *= pj;
					}
				}
				assert (threshold_width_absolute!=Integer.MAX_VALUE || C.length <= node_likelihoods[node].length);
				if (C.length <= node_likelihoods[node].length)
				{
					System.arraycopy(C, 0, node_likelihoods[node], 0, C.length); 
					for (int j=C.length; j<n; j++)
						node_likelihoods[node][j] = Double.NEGATIVE_INFINITY;
				} else 
				{
					// truncated computation; C.length>node_likelihoods[node].length
					System.arraycopy(C, 0, node_likelihoods[node], 0, node_likelihoods[node].length); 
				}
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
			double one_minus_p = getLossParameterComplement(junior);
			
			double loga = Math.log1p(-p*eps); // log(1-pe)
			double logp1 = Math.log(one_minus_p)-loga; // log((1-p)/(1-pe))
			double logp2 = Math.log(p)+Math.log1p(-eps)-loga;// log(p*(1-e)/(1-pe));
			
			if (ASSERT_CALCULATIONS) assert (!Double.isNaN(loga));
			if (ASSERT_CALCULATIONS) assert (!Double.isNaN(logp1));
			if (ASSERT_CALCULATIONS) assert (!Double.isNaN(logp2));
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
					if (ASSERT_CALCULATIONS) assert !Double.isNaN(x);
				}
				
				// calculate C2[ell]
				double log_ellfact = factorials.factln(ell);

				assert (t==ell-s);
				int tmin = t; 
				while (s>=0 && t<C.length) // t<=ell by loop design 
				{
					double binom = log_ellfact - factorials.factln(s)-factorials.factln(t); // ell chose s 
					double slogp1 = s==0?0.0:s*logp1; // avoids 0*Double.NEGATIVE_INFINITY
					double tlogp2 = t==0?0.0:t*logp2; // avoids 0*Double.NEGATIVE_INFINITY
					terms[t-tmin] // fill lower cells
							= K[s] + C[t] + binom +  slogp1 +tlogp2;
//					if (Double.isNaN(terms[t-tmin]))
//					{
//						System.out.println("#**L.P.cS "+family_idx+"/"+junior+"\tell "+ell+"\tC2 "+Arrays.toString(C2)+"\tt "+t+"\ttmin "+tmin
//									+"\ts "+s+"\tK "+K[s]+"\tC "+C[t]
//									+"\tlp1 "+logp1+"\tlp2 "+logp2
//									+"\tp "+p+"\t1p "+one_minus_p+"\tla "+loga+"\teps "+eps+"/"+log1_e);
//					}
					if (ASSERT_CALCULATIONS) assert (!Double.isNaN(binom));
					if (ASSERT_CALCULATIONS) assert (!Double.isNaN(K[s]));
					if (ASSERT_CALCULATIONS) assert (!Double.isNaN(C[t]));
					if (ASSERT_CALCULATIONS) assert (!Double.isNaN(terms[t-tmin]));
					--s;
					++t;
				} 
				C2[ell] = Logarithms.sum(terms, t-tmin);
//				if (Double.isNaN(C2[ell]))
//				{
//					System.out.println("#**L.P.cS "+family_idx+"/"+junior+"\tell "+ell+"\tC2 "+Arrays.toString(C2)+"\tt "+t+"\ttmin "+tmin+"// "+tree.toString(junior)+"\t// "+rates.toString(junior)
//							+"\tterms "+Arrays.toString(terms));
//				}
				if (ASSERT_CALCULATIONS) assert !Double.isNaN(C2[ell]);
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
		 * @deprecated
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
		
		@Override
		public int hashCode()
		{
			int[] profile = table.getFamilyProfile(family_idx);
			return Arrays.hashCode(profile);
		}
		
		public int[] get()
		{
			return table.getFamilyProfile(family_idx);
		}
		
		public int getFamilyIdx()
		{
			return family_idx;
		}
		
		private ProfileTable getTable()
		{
			return table;
		}
		
		public Likelihood getOwner() { return Likelihood.this;}
		
		@Override
		public boolean equals(Object o)
		{
			if (o!=null && o instanceof Profile)
			{
				Profile that = (Profile) o;
				int[] this_profile = this.get();
				assert (table == that.getTable());
				int[] that_profile = that.get();
				return Arrays.equals(this_profile, that_profile);
			} else
				return super.equals(o);
		}

		/**
		 * Shares node and edge likelihood tables between two equal profiles.
		 *  
		 * @param that
		 */
		void copyFromProfile(Profile that)
		{
			assert (equals(that));
			System.arraycopy(that.node_likelihoods, 0, this.node_likelihoods, 0, this.node_likelihoods.length);
			System.arraycopy(that.edge_likelihoods, 0, this.edge_likelihoods, 0, this.edge_likelihoods.length);
		}
	} // Profile class

	
	
	
	
	// ------- testing and main
	
	public void printParameters(java.io.PrintStream out)
	{
		for (int node=0; node<tree.getNumNodes(); node++)
		{
			if (tree.isRoot(node))
				out.print("#ROOTRATES ");
//			double len = rates.getEdgeLength(node);
//			double grate = rates.getGainRate(node);
//			double lrate = rates.getLossRate(node);
//			double drate = rates.getDuplicationRate(node);
			
//			// keeping legacy scaling
//			if (drate==0.0)
//			{
//				if (lrate!=0.0) grate *= lrate;  
//			} else 
//			{
//				grate *= drate;
//			}
//			out.printf("%9g\t%9g\t%9g\t%9g",len, drate, lrate, grate);
			out.println("// p~ "+getLossParameter(node)
					+"\tp "+rates.getLossParameter(node)
					+"\tq~ "+getDuplicationParameter(node)
					+"\tq "+rates.getDuplicationParameter(node)
					+ "\te "+extinction[node]
					+"\t"+tree.toString(node));
			}
	}
	
	
	private void printSingletons(java.io.PrintStream out)
	{
        double L1 = getSingletonLL();
        double p1 = Math.exp(L1);
        int nF = table.getFamilyCount();
        double Exp_n1 = nF*p1;
        int n1=0;
        for (int f=0; f<nF; f++)
        {
        	int[] copies = table.getFamilyProfile(f);
        	int tot_copies=0;
        	for (int leaf=0; leaf<copies.length && tot_copies<2; leaf++)
        	{
        		tot_copies += copies[leaf];
        	}
        	if (tot_copies==1)
        		n1++;
        }
        double L0 = getEmptyLL();
	    double p0 = Math.exp(L0);
        out.println("Singletons L1 "+L1+" (p1 "+p1+") (p0+p1 "+(p0+p1)+")\tn1 "+n1+"\tExp "+Exp_n1);
	
	}
	/**
	 * Instance-linked main for testing from command-line.
	 * 
	 * @param args
	 * @throws Exception
	 */
	private void mainmain(String[] args) throws Exception
	{
		CommandLine cli = new CommandLine(args, Likelihood.class);
        AnnotatedTable table = cli.getTable();
        TreeWithRates rates = cli.getRates();

        PrintStream out = System.out;        
        if (USE_CLASSIC_RECURSION)
        {
            out.println(CommandLine.getStandardHeader("(using classic recursion for node likelihoods)"));
        }
        
        boolean want_distribution =  cli.getOptionBoolean("distribution", false);
        Likelihood factory = new Likelihood(rates, table);
        
        if (cli.getOptionValue(CommandLine.OPT_TRUNCATE)!=null)
        {
        	int absolute = cli.getOptionTruncateAbsolute();
        	double relative = cli.getOptionTruncateRelative();
    		out.println(CommandLine.getStandardHeader("Truncated computation (absolute,relative)="
        				+absolute+","+relative));
        	factory.setCalculationWidthThresholds(absolute, relative);
		}
        int min_copies = cli.getOptionInt(OPT_MINCOPY, table.minCopies());
    	out.println(CommandLine.getStandardHeader("Minimum observed: -"+OPT_MINCOPY+" "
				+min_copies));

    	double LL = factory.getLL();
	    double L0 = factory.getEmptyLL();
	    double p0 = Math.exp(L0);
        double L1 = factory.getSingletonLL();
        double p1 = Math.exp(L1);
        
        double L_unobs = min_copies==0 ? Double.NEGATIVE_INFINITY : L0;
        if (min_copies == 2)
        {
        	L_unobs = Logarithms.add(L0, L1);
        }
        double p_unobs = Math.exp(L_unobs); 
        double p_obs = -Math.expm1(L_unobs);
        double L_obs = Math.log(p_obs);
        int nF = table.getFamilyCount();
    	double corrL = LL-nF*L_obs;
    	
    	double score = -corrL;
    	double ascore = score/nF;
	    out.println("#SCORE "+score);
	    out.println("#AVGSCORE "+ascore);
	    
	    if (want_distribution)
	    {
	    	out.println("# Uncorrected log-likelihood:      \t"+LL+"\t("+Math.exp(LL)+")");
		    out.println("# Empty log-likelihood:\t"+L0+"\t(p0 "+p0+")");
		    out.println("# Singleton log-likelihood:\t"+L1+"\t(p1 "+p1+")"+"\t(p0+p1 "+(p0+p1));
		    out.println("# Corrected log-likelihood:\t"+corrL+"\t("+Math.exp(corrL)+")\t// "+factory.getCorrectedLL());
		    
		    UniqueProfileTable uniqs = new UniqueProfileTable(table);
        	factory = new Likelihood(rates, uniqs);
        	int nU =uniqs.getFamilyCount();

		    out.println("# Number of families:\t"+nF+"\t; unique profiles:\t"+nU);
		    out.print("#\n# Profile\tMultiplicity\tnLin\tnMem\tuLL\tcLL");
		    String[] taxa = table.getTaxonNames();
		    for (int node=0; node<taxa.length; node++) out.print("\t"+taxa[node]+":n");
		    out.println();
        	
        	
        	for (int u=0; u<nU; u++)
        	{
        		int nLin = uniqs.getLineageCount(u);
        		int nMem = uniqs.getMemberCount(u);
        		int[] nCopies = uniqs.getFamilyProfile(u);
        		int mul = uniqs.getMultiplicity(u);
        		Likelihood.Profile P = factory.getProfileLikelihood(u);
        		double uLL = P.getLogLikelihood();
        		double cLL = uLL - L_obs;
        		out.printf("%d\t%d\t%d\t%d\t%f\t%f", u, mul, nLin, nMem, uLL, cLL);
        		for (int t=0; t<nCopies.length; t++)
        			out.print("\t"+nCopies[t]);
        		out.println();
        	}
	    } else
	    {
		    out.println("Log-likelihood:      \t"+LL+"\t("+Math.exp(LL)+")");
		    out.println("Empty log-likelihood:\t"+L0+"\t(p0 "+p0+")");
		    factory.printSingletons(out);
		    out.println("Corrected log-lik.:  \t"+corrL+"\t("+Math.exp(corrL)+")\t// "+factory.getCorrectedLL());
	    }
	}	
	
	public static void main(String[] args) throws Exception
	{
		Likelihood testL = new Likelihood();
		testL.mainmain(args);
	}

}
