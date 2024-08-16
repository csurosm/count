package count.model;
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

import static count.io.CommandLine.OPT_OUTPUT;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Random;

import count.ds.IndexedTree;
import count.ds.Phylogeny;
import count.io.CommandLine;
import count.matek.DiscreteDistribution;
import count.matek.Logarithms;
import count.matek.NegativeBinomial;
import count.matek.PointDistribution;
import count.matek.Poisson;
import count.matek.PseudoRandom;
import count.matek.ShiftedGeometric;

/**
 * Basic class for storing rate parameters on tree edges.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */
public class TreeWithRates implements GLDParameters
{
	/** 
	 * Initialization with random rates.
	 * 
	 * @param tree
	 * @param random_init constant values if null
	 */
	public TreeWithRates(IndexedTree tree, Random random_init)
	{
		this.tree = tree;
		int num_nodes = tree.getNumNodes();
		this.gain_rates = new double[num_nodes];
		this.loss_rates = new double[num_nodes]; 
		this.duplication_rates = new double[num_nodes];
		this.edge_lengths = new double[num_nodes];
		this.rate_gap = new double[num_nodes];
		
		final boolean has_init_length = tree instanceof Phylogeny;
		if (has_init_length)
		{
			scaleByPhylogeny();
		}
		if (random_init == null)
		{
			Arrays.fill(gain_rates, DEFAULT_GAIN_RATE);
			Arrays.fill(loss_rates, DEFAULT_LOSS_RATE);
			Arrays.fill(duplication_rates, DEFAULT_DUPLICATION_RATE);
			Arrays.fill(rate_gap, DEFAULT_LOSS_RATE-DEFAULT_DUPLICATION_RATE);
			if (! has_init_length)
				Arrays.fill(edge_lengths,DEFAULT_EDGE_LENGTH);
			edge_lengths[tree.getRoot()] = Double.POSITIVE_INFINITY;
			RND = null;
		} else
		{
			RND = new PseudoRandom(random_init);
			for (int node=0; node<num_nodes; node++)
			{
				loss_rates[node]= DEFAULT_LOSS_RATE;
				setGainRate(node, RND.nextExponential(1.0/DEFAULT_GAIN_RATE));
				setDuplicationRate(node, RND.nextUniform()*DEFAULT_DUPLICATION_RATE); //  RND.nextUniform(); // less than 1.0
				
				if (tree.isRoot(node))
				{
					setEdgeLength(node, Double.POSITIVE_INFINITY);
//					// DEBUG
//					System.out.println("#**TWR.init rnd root "+toString(node));
					
				} else if (!has_init_length)
				{
					setEdgeLength(node, RND.nextExponential(1.0/DEFAULT_EDGE_LENGTH));
				}
				// DEBUG
//				System.out.println("#**TWR.init rnd node "+toString(node));
			}
		}
	}
	
	public TreeWithRates(IndexedTree tree)
	{
		this(tree, null);
//		this.tree = tree;
//		int num_nodes = tree.getNumNodes();
//		this.gain_rates = new double[num_nodes];
//		this.loss_rates = new double[num_nodes]; 
//		this.duplication_rates = new double[num_nodes];
//		this.edge_lengths = new double[num_nodes];
//		Arrays.fill(gain_rates, DEFAULT_GAIN_RATE);
//		Arrays.fill(loss_rates, DEFAULT_LOSS_RATE);
//		Arrays.fill(duplication_rates, DEFAULT_DUPLICATION_RATE);
//		if (tree instanceof Phylogeny)
//			scaleByPhylogeny();
//		else
//			Arrays.fill(edge_lengths,DEFAULT_EDGE_LENGTH);
//		
//		

	}
		
	/**
	 * Creates a copy linked to the same phylogeny; rate parameters are copied. 
	 * @param same_rates base model
	 */
	public TreeWithRates(TreeWithRates same_rates)
	{
		this(same_rates, true);
	}
	
	/**
	 * 
	 * @param same_phylogeny
	 * @param copy if true, a copy is created for the rate parameters
	 */
	protected TreeWithRates(TreeWithRates same_phylogeny, boolean copy)
	{
		this.tree = same_phylogeny.getTree(); 
		if (copy)
		{
			int num_nodes = tree.getNumNodes();
			this.gain_rates = new double[num_nodes];
			this.loss_rates = new double[num_nodes]; 
			this.duplication_rates = new double[num_nodes];
			this.edge_lengths = new double[num_nodes];
			this.rate_gap = new double[num_nodes];
			for (int node=0; node<num_nodes; node++)
			{
				gain_rates[node] = same_phylogeny.getGainRate(node);
				loss_rates[node] = same_phylogeny.getLossRate(node);
				duplication_rates[node] = same_phylogeny.getDuplicationRate(node);
				edge_lengths[node] = same_phylogeny.getEdgeLength(node);
				rate_gap[node] = same_phylogeny.getRateGap(node);
			}  
		} else
		{
			this.gain_rates = same_phylogeny.gain_rates;
			this.loss_rates = same_phylogeny.loss_rates;
			this.duplication_rates = same_phylogeny.duplication_rates;
			this.edge_lengths = same_phylogeny.edge_lengths;
			this.rate_gap = same_phylogeny.rate_gap;
		}
		this.RND = same_phylogeny.RND; // share random generator 
	}
	


	
	private final IndexedTree tree;
	private final double[] gain_rates;
	private final double[] loss_rates;
	private final double[] duplication_rates;
	private final double[] edge_lengths;
	private final double[] rate_gap;
	
    public static final double DEFAULT_GAIN_RATE = 0.2;
    public static final double DEFAULT_LOSS_RATE = 1.0;
    public static final double DEFAULT_DUPLICATION_RATE = 0.5;
    public static final double DEFAULT_EDGE_LENGTH = 1.0;
    
    
    private PseudoRandom RND;
    public void setRandom(Random RND) 
    {
    	this.RND = RND==null?null:new PseudoRandom(RND);
    }
    
    public PseudoRandom getRandom()
    {
    	return this.RND;
    }
    
    /**
     * Powers of ten.
     */
    private static final double[] POW10 = {1.,10.,100.,1000.,10000.,100000.,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19,1e20,1e21,1e22,1e23,1e24,1e25,1e26,1e27,1e28,1e29,1e30,1e31};

    private void scaleByPhylogeny()
	{
    	Phylogeny phylo = (Phylogeny) tree;
    	// rescale edge lengths 
    	int num_nodes = phylo.getNumNodes(); 
    	int node_idx = num_nodes;
    	assert node_idx>0;
    	double d[] = new double[node_idx];// depths in the tree
    	--node_idx;
    	d[node_idx] = 0.0 ; // root
    	double dmax = 0.0;
    	while (node_idx>0)
    	{
    		--node_idx;
    		int parent_idx = tree.getParent(node_idx);
    		double parent_depth = d[parent_idx];
    		double node_depth = parent_depth + phylo.getLength(node_idx);
    		if (node_depth > dmax) dmax = node_depth;
    		d[node_idx]  = node_depth;
    	}
    	// find a normalizing factor
        double factor = 1.0;
        
        if (Double.isFinite(dmax))
        {
	        int dlog = (int)Math.log10(dmax);
	        if (dlog>0) 
	        {
	        	while (dlog>=POW10.length){factor *= 0.1;--dlog;}
	        	factor *= 1.0/POW10[dlog];
	        } else if (dlog<0)
	        {
	        	while (dlog<=-POW10.length) {factor *=10.0; ++dlog;}
	        	factor *= POW10[dlog];
	        }
        }
        while (node_idx<num_nodes)
        {
        	edge_lengths[node_idx] = phylo.getLength(node_idx)*factor;
        	++node_idx;
        }
	}
	
	/**
	 * The underlying tree. 
	 * 
	 * @return
	 */
    public IndexedTree getTree()
	{
		return tree;
	}

	/**
	 * Relative gain-duplication rate (<var>κ</var>) on the edge leading to a node, (or for root prior),
	 * or relative loss-duplication rate (<var>γ</var>) if no-duplication.
	 * 
	 * @param node_idx any node 
	 * @return non-negative gain rate 
	 */
	public double getGainRate(int node_idx)
	{
		return gain_rates[node_idx];
	}
	
	/**
	 * Logarithm of relative gain-duplication rate (<var>κ</var>) on the edge leading to a node, (or for root prior),
	 * or  of relative loss-duplication rate (<var>γ</var>) if no-duplication.
	 * 
	 * @param node_idx any node 
	 * @return non-negative logarithm of gain rate 
	 */
	public double getLogGainRate(int node_idx)
	{
		return Math.log(getGainRate(node_idx));
	}

	public void setGainRate(int node_idx, double gain_rate)
	{
		gain_rates[node_idx]=gain_rate;
	}

	/**
	 * Loss rate (<var>μ</var>) on the edge leading to a node.
	 * 
	 * @param node_idx non-leaf node
	 * @return non-negative loss rate 
	 */
	public double getLossRate(int node_idx)
	{
		return loss_rates[node_idx];
	}

	public void setLossRate(int node_idx, double loss_rate)
	{
		loss_rates[node_idx]=loss_rate;
		rate_gap[node_idx] = loss_rate-getDuplicationRate(node_idx);
	}

	/**
	 * Duplication rate (<var>λ</var>) on the edge leading to a node, or for root prior. 
	 * 
	 * @param node_idx any node 
	 * @return non-negative duplication rate 
	 */
	public double getDuplicationRate(int node_idx)
	{
		return duplication_rates[node_idx];
	}
	
	public void setDuplicationRate(int node_idx, double dup_rate)
	{
		duplication_rates[node_idx]=dup_rate;
		
		rate_gap[node_idx] = getLossRate(node_idx)-dup_rate;
		
//		{
//			double q = getDuplicationParameter(node_idx);
//			double q1 = getDuplicationParameterComplement(node_idx);
//			
//			
//			if (q==1.0 && q1 ==0.0)
//			{
//				System.out.println("#**TWR.sDR "+node_idx+"\tdr "+dup_rate+"\t"+toString(node_idx));
//			}
//		}
	}
	
	public void setDuplicationRate(int node_idx, double dup_rate, double loss_minus_dup)
	{
		duplication_rates[node_idx] = dup_rate;
		rate_gap[node_idx] = loss_minus_dup;
	}
	
    /**
     * Sets the rate parameters. 
     * 
     * @param node
     * @param length
     * @param dup_rate
     * @param loss_rate
     * @param gain_rate gain/loss if no-duplication (gamma), or gain/duplication (kappa) 
     */
    public void setRates(int node, double length, double gain_rate, double loss_rate, double dup_rate) 
    {
    	this.setEdgeLength(node, length);
    	this.setLossRate(node, loss_rate);
    	this.setDuplicationRate(node, dup_rate);
    	this.setGainRate(node, gain_rate);
    }
	
	
//	public void setRateGap(int node_idx, double delta)
//	{
//		rate_gap[node_idx] = delta;
//		duplication_rates[node_idx] = getLossRate(node_idx)-delta;
//	}
	
	
	
	
	public double getRateGap(int node_idx)
	{
		return rate_gap[node_idx]; 
	}
	
	/**
	 * Edge length (<var>t</var>) on the edge leading to a node, or for root prior. 
	 * 
	 * @param node_idx any node 
	 * @return non-negative gain rate 
	 */
	public double getEdgeLength(int node_idx)
	{
		return edge_lengths[node_idx];
	}
	
	public void setEdgeLength(int node_idx, double length)
	{
		
		//assert (length<100.0);
		edge_lengths[node_idx]=length;
	}
	
	/**
	 * (Random) reinitialization of node rate parameters. 
	 * Duplication rate will be surely inferior to loss rate. 
	 * 
	 * @param node
	 */
	public void initNodeParameters(int node)
	{
		final boolean has_init_length = tree instanceof Phylogeny;
		double mlen = medianLength();
		if (RND == null)
		{
			if (getGainRate(node)!=0.0)
				setGainRate(node, DEFAULT_GAIN_RATE);
			setLossRate(node, DEFAULT_LOSS_RATE);
			if (getDuplicationRate(node)!=0.0)
				setDuplicationRate(node,DEFAULT_DUPLICATION_RATE);
			if (node == tree.getRoot())
				setEdgeLength(node, Double.POSITIVE_INFINITY);
			else if (has_init_length)
				setEdgeLength(node, mlen);
			else
				setEdgeLength(node, DEFAULT_EDGE_LENGTH);
		} else
		{
			setLossRate(node, DEFAULT_LOSS_RATE);
			if (getGainRate(node)!=0.0)
				setGainRate(node, RND.nextExponential(1.0/DEFAULT_GAIN_RATE));
			
			if (getDuplicationRate(node)!=0.0)
				setDuplicationRate(node, 
						//tree.isRoot(node)?0.0: // Poisson root prior by default
						RND.nextUniform()*DEFAULT_DUPLICATION_RATE); //  RND.nextUniform(); // less than 1.0
				
			if (tree.isRoot(node) || Double.isInfinite(getEdgeLength(node)))
			{
				setEdgeLength(node, Double.POSITIVE_INFINITY);
//					// DEBUG
//					System.out.println("#**TWR.init rnd root "+toString(node));
			} else if (has_init_length)
			{
				setEdgeLength(node, RND.nextExponential(1.0/mlen));
			} else 
			{
				setEdgeLength(node, RND.nextExponential(1.0/DEFAULT_EDGE_LENGTH));
			}
		}		
	}
	
	private double medianLength()
	{
		double[] l = edge_lengths.clone();
		Arrays.sort(l);
		return l[l.length/2];
	}
	
	
	public boolean hasGain()
	{
		for (int node=0; node<gain_rates.length; node++)
			if (!tree.isRoot(node) && gain_rates[node]>0.0)
				return true;
		return false;
	}
	
	public boolean hasDuplication()
	{
		for (int node=0; node<duplication_rates.length; node++)
			if (!tree.isRoot(node) && duplication_rates[node]>0.0 )
				return true;
		return false;
	}
	
	public DiscreteDistribution getRootDistribution()
	{
		int root  = tree.getRoot();
		DiscreteDistribution root_distribution
		= (getLossParameter(root)==1.0)
				?getGainDistribution(root)
				:getDuplicationDistribution(root);
		return root_distribution;
	}
	
	public double getRootMean()
	{
		DiscreteDistribution D = getRootDistribution();
		int root  = tree.getRoot();
		double mean;
		if (D instanceof NegativeBinomial)
		{
			double r = getGainParameter(root);
			double q = getDuplicationParameter(root);
			double q1 = getDuplicationParameterComplement(root);
			mean = r*q/q1;
		} else if (D instanceof Poisson)
		{
			mean = getGainParameter(root);
		} else if (D instanceof ShiftedGeometric)
		{
			double p1 = getLossParameterComplement(root);
			double q1 = getDuplicationParameterComplement(root);
			mean = p1/q1;
		} else
		{
			// point distribution
			mean = getLossParameterComplement(root);
		}
		return mean;
		
	}
	
	public void setRootDistribution(DiscreteDistribution root_prior)
	{
        int root_idx = tree.getRoot();
        double[] params = root_prior.getParameters();
        if (root_prior instanceof Poisson)
        {
        	setParameters(root_idx, params[0], 1.0, 0.0); // loss prob=1 at root...
        } else if (root_prior instanceof NegativeBinomial)
        {
        	setParameters(root_idx, params[NegativeBinomial.GAIN_IDX], 1.0, params[NegativeBinomial.DUPLICATION_IDX]);
        } else if (root_prior instanceof ShiftedGeometric)
        {
        	setParameters(root_idx, 0.0, params[ShiftedGeometric.LOSS_IDX], params[ShiftedGeometric.DUPLICATION_IDX]);
        } else
        {
        	assert (root_prior instanceof PointDistribution);
        	setParameters(root_idx, 0.0, params[0], 0.0);
        }
		
	}
	
	public double getUniversalGainParameter(int node)
	{
		double r; // return value
		double gain_rate = getGainRate(node);
		double log1_q = getLogDuplicationComplement(node);
		if (log1_q == 0.0)
		{
			// q==0.0; Poisson
			r = gain_rate*getLossParameter(node);
		} else
		{
			r = -gain_rate*log1_q;
		}
		return r;
	}

	@Override
	public double getGainParameter(int node_idx)
	{
		double gain_rate = getGainRate(node_idx);
		double gainParameter;
		if (getLossRate(node_idx)==0.0 || getDuplicationRate(node_idx)!=0.0 || gain_rate==0.0)
			gainParameter = gain_rate;
		else
			gainParameter = gain_rate * getLossParameter(node_idx); // gain-loss-noduplication = Poisson
		
//		if (!Double.isFinite(gainParameter)) // DEBUG
//		{
//			String node_name = (tree.isLeaf(node_idx)?IndexedTree.LEAF_IDENT_PREFIX:IndexedTree.NODE_IDENT_PREFIX)+node_idx;
//
//			StringBuilder sb = new StringBuilder(node_name);
//			double p = getLossParameter(node_idx);
//			double q = getDuplicationParameter(node_idx);
//			sb.append("[p ").append(p)
//			.append("/1-").append(getLossParameterComplement(node_idx))
//			.append(", q ").append(q)
//			.append("/1-").append(getDuplicationParameterComplement(node_idx))
//			.append(", r ").append(gainParameter);
//			if (q<p)
//			{
//				double n = q*gainParameter/(p-q);
//				sb.append("; n ").append(n);
//			} else
//			{
//				sb.append("; n ").append(Double.POSITIVE_INFINITY);
//			}
//			sb.append("; gr ").append(getGainRate(node_idx))
//			.append(", lr ").append(getLossRate(node_idx))
//			.append(", dr ").append(getDuplicationRate(node_idx))
//			.append(", len ").append(getEdgeLength(node_idx))
//			.append("]");		
//			System.out.println("#***TWR.gGP badgain "+sb.toString());
//		}
		return gainParameter;
	}
	/**
	 * Logarithm of the gain distribution parameter: 
	 * gain intensity <var>r</var> if no-duplication, or gain-duplication rate <var>kappa</var>.
	 * 
	 * 
	 * @param node
	 * @return
	 */
	public double getLogGainParameter(int node)
	{
		return Math.log(getGainParameter(node));
	}
	
	
	@Override
	public double getLossParameter(int node_idx)
	{
		double loss_rate = getLossRate(node_idx);
		double lossParameter;
		if (loss_rate==0.0)
			lossParameter = 0.0;
		else
		{

			double dup_rate = getDuplicationRate(node_idx);
			double t = getEdgeLength(node_idx);
			if (loss_rate == dup_rate)
			{
				double mu_t = loss_rate*t;
				lossParameter = Double.isInfinite(mu_t)
						?1.0:mu_t/(1.0+mu_t);
//				assert (!Double.isNaN(lossParameter));
			} else if (dup_rate<loss_rate)
			{
				if (dup_rate == 0.0)
				{
					// p = 1-e^{-mu t}
					lossParameter = -Math.expm1(-loss_rate*t); // 1.0-Math.exp(-loss_rate*t);
					assert (!Double.isNaN(lossParameter));
				} else
				{
					double gap = getRateGap(node_idx);
					double d = gap*t; // ok if t== +infinity
					double E = Math.exp(-d);
					double E1 = -Math.expm1(-d);
					
					double delta = gap/loss_rate; //positive
					assert (0.0<delta);
					double denom;
					if (delta<0.5)
					{
						denom = loss_rate * (-Math.expm1(-d+Math.log1p(-delta)));
					} else
					{
						denom = loss_rate-dup_rate*E;
					}
					
					
					// p = (mu-mu*E)/(mu-lambda*E)
					lossParameter = loss_rate * E1/denom; //(loss_rate - dup_rate*E);
				}
			} else // separate dup_rate>loss_rate to use formulas with better numerical stability
			{
				if (!(loss_rate<dup_rate))
				{
					System.out.println("#**TWR.gLP "+node_idx+"\tlr "+loss_rate+"\tdr "+dup_rate);
				}
				
				assert (loss_rate<dup_rate);
				
				
				double gap = -getRateGap(node_idx);
				double d = gap*t; 
				double E = Math.exp(-d);
				double E1 = -Math.expm1(-d); // 1-exp(-d)
				
				double delta = gap/dup_rate;
				double denom;
				if (delta<0.5)
				{
					// lr/dr = 1+(lr-dr)/dr
					denom = dup_rate * (-Math.expm1(-d+Math.log1p(-delta)));
				} else
				{
					denom = dup_rate-loss_rate*E;
				}
				
				// p = (mu-mu*E)/(lambda-mu*E)
				
				lossParameter = loss_rate * E1/denom; //(dup_rate-loss_rate*E);
			}
		}
		return lossParameter;
	}
	
	@Override
	public double getLossParameterComplement(int node)
	{
		double loss_rate = getLossRate(node);
		double lossParameterC;
		if (loss_rate==0.0)
			lossParameterC = 1.0;
		else
		{

			double dup_rate = getDuplicationRate(node);
			double t = getEdgeLength(node);
			if (loss_rate == dup_rate)
			{
				double mu_t = loss_rate*t;
				lossParameterC = 1.0/(1.0+mu_t); // ok with t == +infinity
//				assert (!Double.isNaN(lossParameter));
			} else if (dup_rate<loss_rate)
			{
				if (dup_rate == 0.0)
				{
					// p = 1-e^{-mu t}
					double mu_t = loss_rate*t;
					lossParameterC = Math.exp(-mu_t); 
				} else 
				{
					double gap = getRateGap(node); //(loss_rate-dup_rate); // positive 
					double d = gap*t;
					double E = Math.exp(-d); // ok if t== +infinity
					// p = (mu-mu*E)/(mu-lambda*E)
					
					double delta = gap/loss_rate;
					double denom;
					if (delta<0.5)
					{
						denom = loss_rate * (-Math.expm1(-d+Math.log1p(-delta)));
					} else
					{
						denom = loss_rate-dup_rate*E;
					}					
					lossParameterC = gap * E/denom; //(loss_rate - dup_rate*E);
				}
			} else
			{
				assert (loss_rate<dup_rate);
				double gap = -getRateGap(node); // positive
				double d = gap*t; 
				double E = Math.exp(-d);
				// p = (mu-mu*E)/(lambda-mu*E)
				double delta = gap/dup_rate;
				double denom;
				if (delta<0.5)
				{
					// lr/dr = 1+(lr-dr)/dr
					denom = dup_rate * (-Math.expm1(-d+Math.log1p(-delta)));
				} else
				{
					denom = dup_rate-loss_rate*E;
				}
				lossParameterC = gap /denom; //(dup_rate-loss_rate*E);
			}
		}
		return lossParameterC;
	}
	
//	public double getLossDuplicationParameterDifference(int node)
//	{
//		double gap = getRateGap(node);
//		double p = getLossParameter(node);
//		double p_minus_q;
//		if (p==0.0)
//		{
//			double q = getDuplicationParameter(node);
//			if (q==0.0)
//				p_minus_q = 0.0;
//			else 
//			{
//				double delta = gap/q;
//				p_minus_q = q*delta;
//			}
//		} else
//		{
//			double delta = gap/p;
//			p_minus_q = p*delta;
//		}
//		return p_minus_q;
//	}
	
	@Override
	public double getDuplicationParameter(int node_idx)
	{
		double dup_rate = getDuplicationRate(node_idx);
		double t = getEdgeLength(node_idx);
		double duplicationParameter;
		if (dup_rate==0.0)
			duplicationParameter = 0.0;
		else 
		{
			double loss_rate = getLossRate(node_idx);
			if (loss_rate == dup_rate)
			{
				double lambda_t = dup_rate*t;
				duplicationParameter = Double.isInfinite(lambda_t)
						?1.0:lambda_t/(1.0+lambda_t);
			} else if (dup_rate<loss_rate)
			{
				double gap = getRateGap(node_idx);
				double d = gap*t; //  (loss_rate-dup_rate)*t;
				double E = Math.exp(-d);
				double E1 = -Math.expm1(-d);
				// q = (lambda-lambda*E)/(mu-lambda*E)
				
				double delta = gap/loss_rate; //positive
				assert (0.0<delta);
				double denom;
				if (delta<0.5)
				{
					// dr/lr = 1-(lr-dr)/lr
					denom = loss_rate * (-Math.expm1(-d+Math.log1p(-delta)));
				} else
				{
					denom = loss_rate-dup_rate*E;
				}
				duplicationParameter = dup_rate * E1/denom; //(loss_rate - dup_rate*E);
				
				
//				if (Double.isNaN(duplicationParameter))
//				{
//					System.out.println("#**TWR.gDP numerr node "+node_idx+"\td "+d+"\tE "+E+"\tE1 "+E1+"\tdr "+dup_rate+"\tlr "+loss_rate);
//				}
//				assert (!Double.isNaN(duplicationParameter));
			} else
			{
				assert (loss_rate<dup_rate);
				double gap = getRateGap(node_idx); // negative
				double d =  -gap*t; //(dup_rate-loss_rate)*t; 
				double E = Math.exp(-d);
				double E1 = -Math.expm1(-d); // 1-exp(-d)
				// p = (mu-mu*E)/(lambda-mu*E)
				
				double delta = gap/dup_rate;
				
				
				
				
					
				
				double denom;
				if (delta>-0.5)
				{
					// lr/dr = 1+(lr-dr)/dr
					denom = dup_rate * (-Math.expm1(-d+Math.log1p(delta)));
				} else
				{
					denom = dup_rate-loss_rate*E;
				}
				
				duplicationParameter = dup_rate * E1/denom; //(dup_rate-loss_rate*E);
				
				if (!(delta<0.0)) // DEBUG
				{
					System.out.println("#**TWR.gDP "+node_idx+"\tdenom "+denom+"\tor "+(loss_rate-dup_rate*E)+"\tgap "+gap+"\tdr "+dup_rate);
					System.out.println("#**TWR.gDP numerr node "+node_idx+"\td "+d+"\tE "+E+"\tE1 "+E1+"\tdr "+dup_rate+"\tlr "+loss_rate+"\tt "+t+"\tdelta "+delta
								+"\tdparam "+duplicationParameter+"\tdenom "+denom);
				}
				assert (delta<0.0);

				//				if (duplicationParameter>1.0 || duplicationParameter<0.0)
//				{
//					System.out.println("#**TWR.gDP numerr node "+node_idx+"\td "+d+"\tE "+E+"\tE1 "+E1+"\tdr "+dup_rate+"\tlr "+loss_rate+"\tt "+t);
//				}
				// assert (!Double.isNaN(duplicationParameter));
			}
		}
		return duplicationParameter;
	}
	
	

	@Override
	public double getDuplicationParameterComplement(int node)
	{
		double dup_rate = getDuplicationRate(node);
		double t = getEdgeLength(node);
		double duplicationParameterC;
		if (dup_rate==0.0)
			duplicationParameterC = 1.0;
		else 
		{
			double loss_rate = getLossRate(node);
			if (loss_rate == dup_rate)
			{
				double lambda_t = dup_rate*t;
				duplicationParameterC = 1.0/(1.0+lambda_t);
			} else if (dup_rate<loss_rate)
			{
				double gap = getRateGap(node); // (loss_rate-dup_rate);
				double d = gap*t;
				double E = Math.exp(-d);
				// q = (lambda-lambda*E)/(mu-lambda*E)
				double delta = gap/loss_rate;
				double denom;
				if (delta<0.5)
				{
					// dr/lr = 1-(lr-dr)/lr
					denom = loss_rate * (-Math.expm1(-d+Math.log1p(-delta)));
				} else
				{
					denom = loss_rate-dup_rate*E;
				}
				
				duplicationParameterC = gap/denom; //(loss_rate - dup_rate*E);
			} else
			{
				assert (loss_rate<dup_rate);
				double gap = -getRateGap(node);//  (dup_rate-loss_rate); 
				double d = gap*t;
				double E = Math.exp(-d);
				// p = (mu-mu*E)/(lambda-mu*E)
				double delta = gap/dup_rate;
				double denom;
				if (delta<0.5)
				{
					// lr/dr = 1+(lr-dr)/dr
					denom = dup_rate * (-Math.expm1(-d+Math.log1p(-delta)));
				} else
				{
					denom = dup_rate-loss_rate*E;
				}
				duplicationParameterC = gap * E/denom; //(dup_rate-loss_rate*E);
			}
		}
		return duplicationParameterC;
	}
	
	public double getLogLossParameter(int node){ return Math.log(getLossParameter(node));}
	public double getLogLossComplement(int node) {return Math.log(getLossParameterComplement(node));}
	
	public double getLogDuplicationParameter(int node){ return Math.log(getDuplicationParameter(node));}
	public double getLogDuplicationComplement(int node) {return Math.log(getDuplicationParameterComplement(node));}
	
	/**
	 * Calculates log(<var>q</var>/<var>p</var>) if (<var>p</var>&ge;<var>q</var>) or 
	 * log(<var>p</var>/<var>q</var>). 
	 * 
	 * @param node
	 * @return always negative
	 */
	public double getLogRelativeRate(int node)
	{
		return -Math.abs(getLogDuplicationParameter(node)-getLogLossParameter(node));
	}

//	private void invertParameters(int node, double p, double q)
//	{
//		assert (p<1.0);
//		
//		double tolerance = 1e-8;
//		
//		double λ = getDuplicationRate(node);
//		double μ = getLossRate(node);
//		
//		double dl_ratio = q/p;
//		
//		setDuplicationRate(node, dl_ratio*μ);
//		// construct a function for edge length
//		DoubleFunction<Double> pForLength
//		 	= t->{setEdgeLength(node, t); return getLossParameter(node)-p;};
//		
//		double max_t = getEdgeLength(node);
//		while (pForLength.apply(max_t)<0.0 && max_t<1e99)
//			max_t *= 2.0;
//		
//		if (pForLength.apply(max_t)<0.0)
//		{
//			System.out.println("#**TWR.iP "+node+" failed for p="+p+" q="+q+"\t// "+toString(node));
//			throw new RuntimeException("Cannot find edge length for parameters");
//		}
//		
////		System.out.println("#**TWR.iP start "+node+"\t"+p+","+q+"\tmax_t "+max_t+"\t// "+tree.getIdent(node));
//		double t = FunctionMinimization.zbrent(pForLength, 0.0, max_t, tolerance);
//		setEdgeLength(node, t);
////		double diff = getLossParameter(node)-p;
////		System.out.println("#**TWR.iP got "+node+"\t"+p+","+q+"\tt "+t+"\t"+toString(node)+"\t// diff "+diff);
//		
//	}
	
	/**
	 * Sets the gain, loss and duplication rates for an edge.
	 * 
	 * @param node edge leads to this child node  
	 * @param gain_param 
	 * @param loss_param loss probability p
	 * @param dup_param duplication probability q 
	 */
	public void setParameters(int node, double gain_param, double loss_param, double dup_param)
	{
		assert (loss_param>=0.0 && loss_param<=1.0);
		assert (dup_param>=0.0 && dup_param<=1.0);
		
		double loss_rate, dup_rate;
		
		if (loss_param==1.0)
		{
			setEdgeLength(node, Double.POSITIVE_INFINITY);
			double μ = getLossRate(node);
			assert (μ!=0.0);
			assert (Double.isFinite(μ));
			assert (Double.isFinite(dup_param));
			setDuplicationRate(node, dup_param*μ);
		} else
		{
			if (dup_param == loss_param)
			{ // includes ==0.0 
				dup_rate = loss_rate = loss_param/(1.0-loss_param);
			} else
			{		
				if (loss_param == 0.0)
				{
					loss_rate = 0.0;
					dup_rate = -Math.log1p(-dup_param);
				} else if (dup_param==0.0)
				{
					dup_rate = 0.0;
					loss_rate = -Math.log1p(-loss_param);
					
				} else
				{
//					invertParameters(node, loss_param, dup_param);
//					double t = getEdgeLength(node);
//					assert (Double.isFinite(t));
//					dup_rate = getDuplicationRate(node)*t;
//					loss_rate = getLossRate(node)*t;
//
//							
					// (1-q)/(1-p) = exp((1-q/p)*μ)
					// ln(1-q)-ln(1-p) = (1-q/p)*μ 
					// and ln ((1-q)/(1-p)) = ln(1+(p-q)/(1-p)) = ln (1+z)
					// we have 
					// delta = 1-q/p = (p-q)/p 
					// So, 
					// ln (1+delta * p/(1-p)) = delta*μ 
					
					// or 
					// (1-p)/(1-q) = exp((1-p/q)*λ)
					// ln(1-p)-ln(1-q) = (1-p/q)*λ
					assert (dup_param != 1.0);
					
					if (dup_param<loss_param)
					{
						double dl_ratio = dup_param/loss_param;
						double delta = (loss_param-dup_param)/loss_param; //1.0-dl_ratio;
						double z = (loss_param-dup_param)/(1.0-loss_param);
						loss_rate = Math.log1p(z)/delta; 
						if (loss_rate == 0.0)
						{
							System.out.println("#**TWR.sP "+node+"\tlr "+loss_rate+"\tz "+z+"\tl1z "+Math.log1p(z));
						}
									// (Math.log1p(-dup_param)-Math.log1p(-loss_param))/delta;
									// Math.log((1.0-dup_param)/(1.0-loss_param))/delta;
						dup_rate = loss_rate * dl_ratio;

					} else
					{
						double ld_ratio = loss_param/dup_param;
						double delta = (dup_param-loss_param)/dup_param; //  1.0-ld_ratio;
						double z = (dup_param-loss_param)/(1.0-dup_param);
						
						dup_rate = Math.log1p(z)/delta;
								// (Math.log1p(-loss_param)-Math.log1p(-dup_param))/delta;
								// Math.log((1.0-loss_param)/(1.0-dup_param))/delta;
						if (dup_rate == 0.0)
						{
							System.out.println("#**TWR.sP "+node+"\tdr "+dup_rate+"\tz "+z+"\tl1z "+Math.log1p(z));
						}
						loss_rate = dup_rate * ld_ratio;
					}
				}
			}
			setTotalRates(node, loss_rate, dup_rate);
//			double t = getEdgeLength(node);
//
//			if (getLossRate(node)==1.0)
//			{
//				assert Double.isFinite(loss_rate);
//				assert Double.isFinite(dup_rate);
//				if (loss_rate == 0.0)
//				{
//					System.out.println("#**TWR.sP "+node+"\tlr "+loss_rate
//							+"\tr "+dup_rate
//							+"\tp "+loss_param+"\tq "+dup_param+"\tr "+gain_param
//							+"\tt "+t
//							+"\t"+toString(node)+"\t// "+tree.toString(node));
//				}
//				assert (loss_rate != 0.0);
//				// mu = 1
//				// loss_rate = mu * length
//				// want dup_rate = lambda * length
//				setEdgeLength(node, loss_rate);
//				double λ = loss_rate==0.0?dup_rate:(dup_rate/loss_rate);
//				setDuplicationRate(node, λ);
//			} else
//			{
//				if (t==Double.POSITIVE_INFINITY)
//				{
//					setLossRate(node, DEFAULT_LOSS_RATE);
//					assert (loss_rate != 0.0); // cannot have 0 edge length
//					t = loss_rate/DEFAULT_LOSS_RATE;
//					setEdgeLength(node, t);
//					setDuplicationRate(node, dup_rate/t);
//				} else
//				{
//					assert Double.isFinite(loss_rate);
//					assert Double.isFinite(dup_rate);
//					assert (t!=0.0);
//					// keep edge length
//					setLossRate(node, loss_rate/t);
//					setDuplicationRate(node, dup_rate/t);
//				}
//			}
		}

		double gain_rate;
		if (gain_param==0.0)
		{
			gain_rate = 0.0;
		} else
		{
			if (dup_param == 0.0 && loss_param != 0.0)
				gain_rate = gain_param/loss_param;
			else
				gain_rate = gain_param;
		}
		
		setGainRate(node, gain_rate);
//		{ // DEBUG
//			System.out.println("#**TWR.sP done "+node
//				+"\tp "+loss_param
//				+"\tq "+dup_param
//				+"\tr "+gain_param
//				+"\t// "+toString(node));
//		}
						
	}
	
	/**
	 * Scales edge length, loss and duplication rates together. 
	 * 
	 * @param node
	 * @param loss_rate
	 * @param dup_rate
	 */
	
	public void setTotalRates(int node, double loss_rate, double dup_rate)
	{
		double t = getEdgeLength(node);

		if (getLossRate(node)==1.0) // standard scaling
		{
			assert Double.isFinite(loss_rate);
			assert Double.isFinite(dup_rate);
//			if (loss_rate == 0.0)
//			{
//				System.out.println("#**TWR.sP "+node+"\tlr "+loss_rate
//						+"\tr "+dup_rate
//						+"\tp "+loss_param+"\tq "+dup_param+"\tr "+gain_param
//						+"\tt "+t
//						+"\t"+toString(node)+"\t// "+tree.toString(node));
//			}
			assert (loss_rate != 0.0);
			// mu = 1
			// loss_rate = mu * length
			// want dup_rate = lambda * length
			setEdgeLength(node, loss_rate);
			double λ = loss_rate==0.0?dup_rate:(dup_rate/loss_rate);
			setDuplicationRate(node, λ);
		} else
		{
			if (t==Double.POSITIVE_INFINITY)
			{
				setLossRate(node, DEFAULT_LOSS_RATE);
				assert (loss_rate != 0.0); // cannot have 0 edge length
				t = loss_rate/DEFAULT_LOSS_RATE;
				setEdgeLength(node, t);
				setDuplicationRate(node, dup_rate/t);
			} else
			{
				assert Double.isFinite(loss_rate);
				assert Double.isFinite(dup_rate);
				assert (t!=0.0);
				// keep edge length
				setLossRate(node, loss_rate/t);
				setDuplicationRate(node, dup_rate/t);
			}
		}
	}
	
	

	
//	/**
//	 * log(1-x) if we know both x and 1.0-x at double precision. Uses the smaller 
//	 * value for calculation. 
//	 * 
//	 * @param x variable at double precision
//	 * @param xcomp complement at double precision 
//	 * @return
//	 */
//	private static final double log1m(double x, double xcomp)
//	{
//		return x<xcomp?Math.log1p(-x):Math.log(xcomp);
//	}
//	/**
//	 * High-precision rate setting.
//	 * 
//	 * @param node
//	 * @param p
//	 * @param p_1 1-p
//	 * @param q
//	 * @param q_1 1-q
//	 */
//	public void setDuplicationLossRates(int node, double p, double p_1, double q, double q_1)
//	{
//		if (p_1==0.0)
//		{
//			setEdgeLength(node, Double.POSITIVE_INFINITY);
//			double μ = getLossRate(node);
//			assert (μ!=0.0);
//			assert (Double.isFinite(μ));
//			assert (Double.isFinite(q));
//			setDuplicationRate(node, q*μ);
//		} else
//		{
//			double λ,μ;
//			double rate_gap;
//			if (p==q && p_1==q_1) // equal rates
//			{
//				λ=μ=p/p_1;
//				rate_gap = 0.0;
//			} else 
//			{
//				if (p == 0.0)
//				{
//					μ = 0.0;
//					
////					dup_rate = -Math.log1p(-dup_param);
////				} else if (dup_param==0.0)
////				{
////					dup_rate = 0.0;
////					loss_rate = -Math.log1p(-loss_param);
//					
//					
//					λ = -GLDParameters.log1m(q, q_1); //  -log1m(q,q_1);
//					rate_gap = -λ;
//				} else if (q==0.0)
//				{
//					λ = 0.0;
//					μ = -GLDParameters.log1m(p,p_1);
//					rate_gap = μ;
//				} else
//				{
//					double diff;
//					if (q<p || p_1<q_1)
//					{
//						diff = Double.max(p-q,q_1-p_1); // the larger value should have better precision when p or q is close to 0 or 1 (?)
//					} else
//					{
//						diff = -Double.max(q-p, p_1-q_1);
//					}
//					rate_gap = GLDParameters.log1m(q,q_1)-GLDParameters.log1m(p,p_1); // mu-lambda = ln(1-q)-ln(1-p)
//					μ = p/diff * rate_gap;
//					λ = q/diff * rate_gap;
//				}
//			}
//			double t = getEdgeLength(node);
//
//			if (getLossRate(node)==1.0)
//			{
//				assert Double.isFinite(λ);
//				assert Double.isFinite(μ);
//				assert (μ != 0.0);
//				// mu = 1
//				// loss_rate = mu * length
//				// want dup_rate = lambda * length
//				setEdgeLength(node, μ);
//				λ /= μ; rate_gap /= μ;
//				setDuplicationRate(node, λ, rate_gap);
//			} else
//			{
//				if (t==Double.POSITIVE_INFINITY)
//				{
//					setLossRate(node, DEFAULT_LOSS_RATE);
//					assert (μ != 0.0); // cannot have 0 edge length
//					t = μ/DEFAULT_LOSS_RATE;
//					setEdgeLength(node, t);
//					λ/=t; rate_gap/=t;
//					setDuplicationRate(node, λ, rate_gap);
//				} else
//				{
//					assert Double.isFinite(μ);
//					assert Double.isFinite(λ);
//					assert (t!=0.0);
//					// keep edge length
//					setLossRate(node, μ/t);
//					setDuplicationRate(node, λ/t, rate_gap/t);
//				}
//			}					
//		}
//	}
	
	/**
	 * Adjusts edge lengths while keeping the rates the same. 
	 * 
	 * @param node
	 * @param p loss parameter 
	 * @param one_minus_p ==1-p (redundancy for better numerical precision) 
	 */
	public void setEdgeLengthForLossParameter(int node, double p, double one_minus_p)
	{
		double μ = getLossRate(node);
		//double λ = getDuplicationRate(node);
		double diff = getRateGap(node); //μ-λ;
		//double dl_ratio = ;
		double delta = diff/μ ; // == 1.0-λ/μ; 
		
		double z = delta*p/one_minus_p; 
		double w = Math.log1p(z);
		double t = w/diff;
		setEdgeLength(node, t);
	}
	
	public void setGainRateForLossParameter(int node, double gain_param, double p)
	{
		double gain_rate=gain_param;
		double λ = getDuplicationRate(node);
		if (λ==0.0 && p!=0.0)
		{
			gain_rate /= p;
		} 
		setGainRate(node, gain_rate);
	}
	
	public String toString(int node)
	{
		String node_name = (tree.isLeaf(node)?IndexedTree.LEAF_IDENT_PREFIX:IndexedTree.NODE_IDENT_PREFIX)+node;

		StringBuilder sb = new StringBuilder(node_name);
		double p = getLossParameter(node);
		double q = getDuplicationParameter(node);
		double r = getGainParameter(node);
		double ru = getUniversalGainParameter(node);
		sb.append("[p ").append(p)
		.append("/1-").append(getLossParameterComplement(node))
		.append(", q ").append(q)
		.append("/1-").append(getDuplicationParameterComplement(node))
		.append(getLogDuplicationParameter(node)==Double.NEGATIVE_INFINITY?", r ":", kapa ").append(r)
		.append("/ru ").append(ru);
		if (q<p)
		{
			double n = q*r/(p-q);
			sb.append("; n ").append(n);
		} else
		{
			sb.append("; n ").append(Double.POSITIVE_INFINITY);
		}
		sb.append("; gr ").append(getGainRate(node))
		.append(", lr ").append(getLossRate(node))
		.append(", dr ").append(getDuplicationRate(node));
		sb.append(", len ").append(getEdgeLength(node));
		return sb.append("]").toString();
	}
	
	
	private void reportParameters(java.io.PrintStream out)
	{
		int num_nodes = tree.getNumNodes();
		for (int node=0; node<num_nodes; node++)
		{
			out.println("NODE\t"+node+"\t"+toString(node));
		}
	}
	
	
	public static void main(String[] args) throws Exception
	{
		Class<?> us = java.lang.invoke.MethodHandles.lookup().lookupClass();
		CommandLine cli = new CommandLine(args,  us);
		
        PrintStream out = System.out; 
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(us));
    	    out.println(CommandLine.getStandardRuntimeInfo(us, args));
    	}
		
        MixedRateModel zeb  = cli.getMixedrateModel();
        TreeWithRates rates = zeb.getBaseModel();
        
        
//        String maxgain = cli.getOptionValue("max"+CommandLine.OPT_GAIN);
//        if (maxgain == null)
//        {
        	rates.reportParameters(out);
//        } else
//        {
//        	double max = Double.parseDouble(maxgain);
//        	rates.setMaxGainDuplicationRate(max);
//        	
//        	out.println(count.io.RateVariationParser.printRates(zeb));        	
//        }
        
		if (out != System.out)
		{
			if (out.checkError())
				throw new java.io.IOException("Write failed.");
			out.close();
		}
        
	}
}
