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
import count.matek.Functions.RisingFactorial;
//import count.model.junkyard.DirectLikelihood;
import count.matek.Logarithms;

/**
 * Likelihood computation by direct application of 
 * Felsenstein's peeling method, conditioning on ancestor copies. 
 * Loss and duplication (<var>p</var> and <var>q</var>) distribution parameters 
 * are represented in logistic space via logit(<var>p</var>)=ln(<var>p</var>/(1-<var>p</var>)) and 
 * logit(<var>q</var>)=ln(<var>q</var>/(1-<var>q</var>)). Gain (<var>r</var>) parameter is 
 * represented directly.
 * 
 * Use {@link #setCalculationWidthThresholds(int, double)} 
 * for settng the maximum ancestor copies assumed. 
 */
public class StraightLikelihood implements GLDParameters
{
	private static final boolean TRACK_CALCULATIONS = false; // for debugging
	private static boolean USE_UNIVERSAL_GAIN_PARAMETER = true;
	
	/**
	 * Usual instantiation with a rate model. 
	 * 
	 * @param rates
	 * @param table
	 */
	public StraightLikelihood(TreeWithLogisticParameters rates, ProfileTable table)
	{
		this.rates = rates;
		this.table = table;
		this.tree = rates.getTree();
		this.max_node_copies = table.getMaxCopies(tree);//   table.getMaxFamilySizes(tree);
		this.transition_probs = new TransitionProbabilities[tree.getNumNodes()];
		this.parameter_factory = null;
		this.initDataStructures();
		this.copyParametersFromRates();
	}
	
	public StraightLikelihood(TreeWithRates rates, ProfileTable table)
	{
		this(
			(rates instanceof TreeWithLogisticParameters
					?(TreeWithLogisticParameters)rates
					:new TreeWithLogisticParameters(rates, false))
			, table);
	}
	
	
	/**
	 * Instantiation with parameters copied from another instance.
	 * 
	 * @param that
	 * @param table
	 */
	public StraightLikelihood(StraightLikelihood that, ProfileTable table)
	{
		this.rates = that.rates;
		this.table = table;
		this.tree = that.tree;
		this.max_node_copies = table.getMaxCopies(tree); // table.getMaxFamilySizes(tree);
		this.transition_probs = new TransitionProbabilities[tree.getNumNodes()];
		this.initDataStructures();
		this.parameter_factory = that;
		this.copyParametersFromFactory();
	}
	
	
	protected final TreeWithLogisticParameters rates;
	protected final ProfileTable table;
	protected final IndexedTree tree;
	protected RisingFactorial factorials;
	protected final int[] max_node_copies;
	protected final TransitionProbabilities[] transition_probs;
	protected final StraightLikelihood parameter_factory;
	
	/**
	 * Precomputed rising factorials for Polya pmf. Calculation 
	 * of {@link RisingFactorial#factln(int)} is 
	 * correct even if too few members are allocated + precomputed, 
	 * but it's slower.  
	 */
	private RisingFactorial[] gain_factorials;
	
	/**
	 * 3 parameters per node
	 */
	private double[] node_parameters;

	
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
		this.initDataStructures();
		if (parameter_factory==null)
			this.copyParametersFromRates();
		else 
			this.copyParametersFromFactory();
	}	
	
	/**
	 * Sets the parameters from underlying model or parameter factory
	 */
	public void computeParameters()
	{
		if (parameter_factory==null)
			this.copyParametersFromRates();
		else 
			this.copyParametersFromFactory();
	}
	
	/**
	 * Copies calculation  width from another instance.
	 * 
	 * @param factory
	 */
	protected void setSameCalculationWidthThresholds(StraightLikelihood factory)
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
	 * Allocates the arrays for the parameter data structures. 
	 */
	private void initDataStructures()
	{
		int num_nodes = tree.getNumNodes();
		this.node_parameters = new double[3*num_nodes];
		this.gain_factorials = new RisingFactorial[num_nodes];
		int m = getCalculationWidth(max_node_copies[tree.getRoot()]);
		factorials = new RisingFactorial(1.0, 1+m); 		
	}
	
	/**
	 * Updates the optimization parameters
	 * 
	 * @param node
	 * @param logit_p logit for loss parameter
	 * @param logit_q logit for duplication parameter
	 * @param r universal gain parameter
	 */
	public void setUniversalNodeParameters(int node, double logit_p, double logit_q, double r)
	{
		// TODO  use TreeWithLogisticParameters

		if (!USE_UNIVERSAL_GAIN_PARAMETER && logit_q != Double.NEGATIVE_INFINITY)
		{
			double log1_q = Logarithms.logitToLogComplement(logit_q);
			r /= -log1_q; // for kappa
		}
		
		setLogitLossParameter(node,logit_p);
		setLogitDuplicationParameter(node, logit_q);
		setGainParameter(node, r);
		
		this.transition_probs[node] = new TransitionProbabilities(node);			
	}
	
	/**
	 * Updates the optimization parameters
	 * 
	 * @param node
	 * @param logit_p
	 * @param logit_q
	 * @param r
	 */
	private void setNodeParameters(int node, double logit_p, double logit_q, double r)
	{
		// TODO  use TreeWithLogisticParameters

		setLogitLossParameter(node,logit_p);
		setLogitDuplicationParameter(node, logit_q);
		setGainParameter(node, r);
		
		this.transition_probs[node] = new TransitionProbabilities(node);
		
//		System.out.println("#**SL.sNP "+node+"\tlp "+logit_p+"\tlq "+logit_q+"\tr "+r);
	}

	public void setParametersAndCopy(double[] params)
	{
		assert (params.length == node_parameters.length);
		System.arraycopy(params, 0, node_parameters, 0, params.length);
		this.copyParametersToModel();
	}
	
	public double[] getParameters()
	{
		return node_parameters.clone();
	}
	
	public double getLogitLossParameter(int node) { return this.node_parameters[3*node+PARAMETER_LOSS];}
	public double getLogLossParameter(int node)
	{
		double logit_p = this.getLogitLossParameter(node);
		double log_p = Logarithms.logitToLogValue(logit_p);
		return log_p;
	}
	public double getLogLossComplement(int node)
	{
		double logit_p = this.getLogitLossParameter(node);
		double log1_p = Logarithms.logitToLogComplement(logit_p);
		return log1_p;
	}
	@Override
	public double getLossParameter(int node)
	{
		return Math.exp(this.getLogLossParameter(node));
	}
	@Override
	public double getLossParameterComplement(int node)
	{
		return Math.exp(this.getLogLossComplement(node));
	}
	
	public double getLogitDuplicationParameter(int node) {return node_parameters[3*node+PARAMETER_DUPLICATION];}
	
	public double getLogDuplicationParameter(int node)
	{
		double logit_q = this.getLogitDuplicationParameter(node);
		double log_q = Logarithms.logitToLogValue(logit_q);
		return log_q;
	}
	
	public double getLogDuplicationComplement(int node)
	{
		double logit_q = this.getLogitDuplicationParameter(node);
		double log1_q = Logarithms.logitToLogComplement(logit_q);
		return log1_q;
	}
	@Override 
	public double getDuplicationParameter(int node)
	{
		return Math.exp(this.getLogDuplicationParameter(node));
	}
	@Override 
	public double getDuplicationParameterComplement(int node)
	{
		return Math.exp(this.getLogDuplicationComplement(node));
	}
	
	@Override
	public double getGainParameter(int node) {return node_parameters[3*node+PARAMETER_GAIN];}
	
	public double getUniversalGainParameter(int node)
	{
		double r = node_parameters[3*node+PARAMETER_GAIN];
		if (!USE_UNIVERSAL_GAIN_PARAMETER)
		{
			double log1_q = getLogDuplicationComplement(node);
			if (getLogitDuplicationParameter(node)!=Double.NEGATIVE_INFINITY)  //(log1_q!=0.0)
			{
				// this is kappa
				double kappa = r;
				r = -kappa*log1_q;
			}
		}
		return r;
	}
	
	
	/**
	 * Sets gain parameter (locally stored); recomputes the gain_factorials.
	 * 
	 * @param node
	 * @param g kappa (if Polya), or r (if Poisson) ; r if using universal gain param
	 */
	private void setGainParameter(int node, double g) 
	{
		int j = 3*node+PARAMETER_GAIN;
		node_parameters[j]=g;

		
		// and precalculate rising factorials for this kappa, if Polya 
		double y = getLogitDuplicationParameter(node);
		if (y==Double.NEGATIVE_INFINITY) // Poisson with log(q/(1-q)) = log(0)
		{
			this.gain_factorials[node] = null;
//			System.out.println("#**DL.cP "+node+"\tq=0"+"\tr "+r);
		} else if (g!=0.0) // Polya distribution
		{
			int m = max_node_copies[tree.isLeaf(node)?node:tree.getRoot()]; // logic of allocation for node_likelihoods[]
			int cw = getCalculationWidth(m); 
			double kappa;
			if (USE_UNIVERSAL_GAIN_PARAMETER)
			{
				kappa = -g/getLogDuplicationComplement(node);
			} else
			{
				kappa = g;
			}
			
			this.gain_factorials[node] = new RisingFactorial(kappa, 1+cw);	
		}
		if (g==0.0)
		{
			System.out.println("#**SL.sGP "+node+"\tgain "+g+"\ty "+y+"(l1-q "+getLogDuplicationComplement(node)+")\tx "+getLogitLossParameter(node));
			this.gain_factorials[node]  = null;
		}

	
	}

	/**
	 * Loss parameter setting; must call {@link #setGainParameter(int, double)} after. 
	 * 
	 * @param node
	 * @param logit_p loss probability for a parental copy towards this lineage
	 */
	private void setLogitLossParameter(int node, double logit_p)
	{
		node_parameters[3*node+PARAMETER_LOSS]=logit_p;
	}
	/**
	 * Duplication parameter setting; must call {@link #setGainParameter(int, double)} after.  
	 * 
	 * @param node
	 * @param logit_q duplication parameter (for Polya) 
	 */
	private  void setLogitDuplicationParameter(int node, double logit_q) 
	{
		int j = 3*node+PARAMETER_DUPLICATION;
		node_parameters[j]=logit_q;
	}	
	
	
	public void copyParametersFromFactory()
	{
		for (int node=0; node<tree.getNumNodes(); node++)
		{
			double logit_p = parameter_factory.getLogitLossParameter(node);
			double logit_q = parameter_factory.getLogitDuplicationParameter(node);
			double r = parameter_factory.getGainParameter(node);
			setNodeParameters(node, logit_p, logit_q, r);
		}
		
	}
	
	public void copyParametersFromRates()
	{
		for (int node=0; node<tree.getNumNodes(); node++)
		{
			double logit_p = rates.getLogLossParameter(node)-rates.getLogLossComplement(node);
			double logit_q = rates.getLogDuplicationParameter(node)-rates.getLogDuplicationComplement(node);
			double r;
			if (USE_UNIVERSAL_GAIN_PARAMETER)
				r = rates.getUniversalGainParameter(node);
			else
				r = rates.getGainParameter(node);
			setNodeParameters(node, logit_p, logit_q, r);
		}
	}
	
	public void copyParametersToModel()
	{
		for (int node=0; node<tree.getNumNodes(); node++)
		{
			double x = this.getLogitLossParameter(node);
			double y = this.getLogitDuplicationParameter(node);			
			double r = getUniversalGainParameter(node);
			
			if (y==Double.NEGATIVE_INFINITY)
			{
				// Poisson
				rates.setLogitLossDuplication(node, x, y, r);
			} else
			{
				// Polya
				double log1_q = Logarithms.logitToLogComplement(y);
				double kappa = -r/log1_q;
				rates.setLogitLossDuplication(node, x, y, kappa);
			}

//			double log_p, log1_p;
//			if (0<=x)
//			{
//				log_p = Logarithms.logitToLogValue(x);
//				log1_p = log_p-x;
//			} else
//			{
//				log1_p = Logarithms.logitToLogComplement(x);
//				log_p = log1_p+x;
//			}
//			
//			double y = this.getLogitDuplicationParameter(node);			
//			double log_q,log1_q;
//			if (0<=y)
//			{
//				log_q = Logarithms.logitToLogValue(y);
//				log1_q = log_q-y;
//			} else
//			{
//				log1_q = Logarithms.logitToLogComplement(y);
//				log_q = log1_q+y;
//			}
			
//			if (x==Double.POSITIVE_INFINITY) // p=1; ln p/(1-p) = ln (1/0)
//			{
//				// loss p=1
//				double μ = rates.getLossRate(node);
//				assert (μ!=0.0);
//				assert (Double.isFinite(μ));
//				double mul = Math.exp(log_q);
//				// q/p=q
//				rates.setEdgeLength(node, Double.POSITIVE_INFINITY);
//				rates.setDuplicationRate(node, mul*μ);
//				
//				// DEBUG
//				System.out.println("#**SL.cPTM "+node+"\tx "+x+"\ty "+y+"/rates:"+(rates.getLogDuplicationParameter(node)-rates.getLogDuplicationComplement(node))
//						+"\t// "+toString(node)+"\t"+rates.getClass());
//				
//			} else if (y==Double.NEGATIVE_INFINITY)
//			{
//				// Poisson
//				// p=1-exp(-mu*t)
//				double μt = -log1_p;
//				double μ = rates.getLossRate(node);
//				rates.setDuplicationRate(node, 0.0);
//				rates.setEdgeLength(node, μt/μ);
//			} else
//			{
//				
//				
//				// not infinite length
//				if (x==y)
//				{
//					double μ = rates.getLossRate(node);
//					double μt = Math.exp(x);
//					rates.setDuplicationRate(node, μ);
//					rates.setEdgeLength(node, μt/μ);
//				} else if (y<x)
//				{
//					double log_E = log1_p-log1_q; // Logarithms.logitParameterRatio(-y, -x);
//					
//					double logit_λ = Logarithms.logitParameterRatio(x, y);
//					double log_λ, log1_λ;
//					if (0<=logit_λ)
//					{
//						log_λ = Logarithms.logitToLogValue(logit_λ);
//						log1_λ = log_λ-logit_λ;
//					} else
//					{
//						log1_λ = Logarithms.logitToLogComplement(logit_λ);
//						log_λ = log1_λ+ logit_λ;
//					}
//					double μ = rates.getLossRate(node);
//					
//					double μt = -log_E/Math.exp(log1_λ);
//					rates.setDuplicationRate(node, μ*Math.exp(log_λ));
//					rates.setEdgeLength(node, μt/μ);
//				} else
//				{
//					double log_E = log1_q-log1_p; // Logarithms.logitParameterRatio(-x, -y);
//					double logit_μ = Logarithms.logitParameterRatio(y, x);
//					double log_μ, log1_μ;
//					if (0<=logit_μ)
//					{
//						log_μ = Logarithms.logitToLogValue(logit_μ);
//						log1_μ = log_μ-logit_μ;
//					} else
//					{
//						log1_μ= Logarithms.logitToLogComplement(logit_μ);
//						log_μ = log1_μ+ logit_μ;
//					}
//					double λt = -log_E/Math.exp(log1_μ);
//					double μ = rates.getLossRate(node); // fixed
//					double ld_ratio = Math.exp(log_μ);
//					rates.setDuplicationRate(node, μ/ld_ratio);
//					rates.setEdgeLength(node, λt*ld_ratio/μ);
//				}
//			}
			
//			double r = getUniversalGainParameter(node);
//			
//			if (r==0.0)
//				rates.setGainRate(node, r);
//			else if (y==Double.NEGATIVE_INFINITY && x!=Double.NEGATIVE_INFINITY)
//			{
//				// Poisson w/ loss
//				double gamma = r*Math.exp(-log_p);
//				rates.setGainRate(node, gamma);
//				
////				if (!tree.isRoot(node))
////				System.out.println("#**SL.cPTM "+node+"\tx "+x+"("+Math.exp(log_p)+"/1-"+Math.exp(log1_p)+")"
////						+"\ty "+y+"("+Math.exp(log_q)+"/1-"+Math.exp(log1_q)+")"
////						+"\tr "+r+"\t// "+rates.toString(node));				
//			} else
//			{
//				double κ = -r/log1_q;
//				rates.setGainRate(node, κ);
//			}
//			if (tree.isRoot(node)) // DEBUG
//			{
//				System.out.println("#**SL.cPTM "+node+"\tx "+x
//						+"\ty "+y
//						+"\tr "+r+"\t// "+toString(node));	
//			}

		} // for each node
	}
	
	
	private class TransitionProbabilities
	{
		private final int node;
		
		final double[][] fromparent;
		final double[][] tochild;

		TransitionProbabilities(int node)
		{
			this.node = node;
			
			int nmax = StraightLikelihood.this.max_node_copies[tree.isLeaf(node)?node:tree.getRoot()];
			nmax = getCalculationWidth(nmax);
			int pmax;
			if (tree.isRoot(node))
			{
				pmax = 0;
			} else
			{
				pmax = max_node_copies[tree.getRoot()];
			}
			//System.out.println("#**SL.TP() node "+node+"\tpmax "+pmax+"\tcw "+getCalculationWidth(pmax));
			pmax = getCalculationWidth(pmax);

//			try // DEBUG out-of-memory 
//			{
				//System.out.println("#**SL.TP() node "+node+"\tpmax "+pmax+"\tnmax "+nmax);
				this.fromparent = new double[1+pmax][];
				this.tochild = new double[1+pmax][1+nmax];
				this.calculateValues();
//			} catch (OutOfMemoryError not_enought_memory)
//			{
//				throw not_enought_memory;
//			}
		}		
		
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
		/**
		 * Log-probability of surviving copies.
		 * 
		 * @param n_parent
		 * @param n_surviving_copies
		 * @return
		 */
		double copydeath(int n_parent, int n_surviving_copies)
		{
			return this.fromparent[n_parent][n_surviving_copies]; 
		}

		private void calculateValues()
		{
			double ru = getUniversalGainParameter(node);
			
			double logit_q = getLogitDuplicationParameter(node);
			double log_q, log1_q;
			if (0<=logit_q)
			{
				log_q = Logarithms.logitToLogValue(logit_q);
				log1_q = log_q - logit_q;
			} else
			{
				log1_q = Logarithms.logitToLogComplement(logit_q);
				log_q = log1_q + logit_q;
			}
			int n = 0; // 0->m transitions

			double[] fromparent_n = fromparent[n] =new double[1];
			double[] tochild_n = tochild[0];
			fromparent_n[0]=Math.log(1.0);
			if (logit_q == Double.NEGATIVE_INFINITY) // q==0.0;
			{
				double logr = Math.log(ru);
				// Poisson
				// prob(0->m) = exp(-r)*r^m/m!
				// ok even with r=0
				for (int m=0; m<tochild_n.length; m++)
				{
					double m_logr = m==0?0.0:m*logr;
					double term = -ru + m_logr - factorials.factln(m);
					tochild_n[m] = term; 
//					System.out.println("#***SL.TP.cV "+node+"\ttc "+n+"\t"+m+"\t"+tochild[n][m]);
				}
			} else	
			{
				// Polya
				if (ru==0.0)
				{
					// no-gain, so 0->m only with m=0
					tochild_n[0]=0.0;
					for (int m=1; m<tochild_n.length; m++)
					{
						tochild_n[m]=Double.NEGATIVE_INFINITY;
					}
				} else
				{
					//double rlog1_q = r*log1_q;
					double κlog1_q = -ru;
					for (int m=0; m<tochild_n.length; m++)
					{
						// n-> m transitions
						double binom = gain_factorials[node].factln(m)-factorials.factln(m);
						double m_logq = m==0.0?0.0:m*log_q;
						double term = binom + κlog1_q + m_logq;
						tochild_n[m] = term;
	//					System.out.println("#***SL.TP.cV "+node+"\ttc "+n+"\t"+m+"\t"+tochild[n][m]);
					}
				}
			}
//			System.out.println("#***SL.TP.cV "+node+"\tn "+n+"\t"+Arrays.toString(tochild[n]));
			
			
			++n;
			// 
			if (!tree.isRoot(node))
			{
				double logr = Math.log(ru);
				double κlog1_q = -ru;

				double logit_p = getLogitLossParameter(node);
				double log_p, log1_p;
				if (0<=logit_p)
				{
					log_p = Logarithms.logitToLogValue(logit_p);
					log1_p = log_p - logit_p;
				} else
				{
					log1_p = Logarithms.logitToLogComplement(logit_p);
					log_p = log1_p + logit_p;
				}
				while(n<fromparent.length)
				{
					fromparent_n = fromparent[n] = new double[n+1];
					tochild_n = tochild[n];
					for (int k=0,d=n-k; k<fromparent_n.length; k++, d--)
					{
						double loss_binom = factorials.factln(n)-factorials.factln(k)-factorials.factln(d);
						double klog1_p = k==0?0.0:k*log1_p;
						double dlog_p = d==0?0.0:d*log_p;
						
						double loss_term = loss_binom+klog1_p+dlog_p;
						fromparent_n[k] = loss_term;
					}
//					for (int m=0; m<n && m<tochild_n.length; m++)
//						tochild_n[m]=Double.NEGATIVE_INFINITY;
					Arrays.fill(tochild_n, 0, Integer.min(n, tochild_n.length), Double.NEGATIVE_INFINITY); // tochild[k][m] = log(0.0) when m<k 
					for (int j=0, m=n+j; m<tochild_n.length; m++, j++)
					{
						double dup_term;
						if (logit_q == Double.NEGATIVE_INFINITY)
						{
							// Poisson
							// prob = exp(-r)*r^k/k!
							double jlog_r = j==0?0.0:j*logr;
							dup_term = -ru + jlog_r-factorials.factln(j);
						} else
						{
							// Polya
							double dup_binom;
							if (ru==0.0)
							{
								// no-gain
								dup_binom = factorials.factln(m-1)-factorials.factln(n-1)-factorials.factln(j);
										// n*log1_q + j*log_q: j=n'-s', n=s' so n'=n+j=m s'=n n'-s'=j
										// want binom (n'-1,s'-1)=binom(m-1,n-1)
							} else
							{
								// Polya
								dup_binom = gain_factorials[node].factln(m)-gain_factorials[node].factln(n)
										-factorials.factln(j);
							}
//							= gain_factorials[node].factln(m)-gain_factorials[node].factln(n)
//									-factorials.factln(j);
							dup_term = dup_binom + κlog1_q + n*log1_q + j*log_q;
						}
						tochild[n][m] = dup_term;
//						System.out.println("#***SL.TP.cV "+node+"\ttc "+n+"\t"+m+"\t"+tochild[n][m]);
					} // for j
//					System.out.println("#***SL.TP.cV "+node+"\tn "+n+"\t"+Arrays.toString(tochild[n]));
					++n;
				} // while n
			} // at non-root nodes 
		}
	}
	
	public Profile getProfile(int family)
	{
		return new Profile(family);
	}
	
	/**
	 * Class for storing profile-specific conditional likelihoods.
	 */
	public class Profile
	{
		protected final int family_idx;
		/**
		 * Inside node log-likelihoods C[<var>u</var>][<var>n</var>] for 
		 * <var>n</var> node copies.
		 */
		private final double[][] node_likelihoods;
		/**
		 * Inside edge log-likelihoods K[<var>u</var>][<var>s</var>]
		 * for <var>s</var> eddge copies (inherited from parent).
		 */
		private final double[][] edge_likelihoods;
		
		/**
		 * Complementary sibling log-likelihoods [<var>u</var>][<var>n</var>] for 
		 * <var>n</var> parent copies.
		 */
		private final double[][] sibling_likelihoods;
		
		/**
		 * Outside node log-probabilties B[<var>u</var>][<var>n</var>] for 
		 * <var>n</var> node copies.
		 */
		private final double[][] node_outside;
		/**
		 * Outside edge log-probabilties J[<var>u</var>][<var>s</var>]
		 * for <var>s</var> eddge copies (inherited from parent).
		 */
		private final double[][] edge_outside;
		/**
		 * Inside edge log-likelihoods for birth transitions 
		 * K[<var>u</var>][<var>s</var>][<var>n</var>] for 
		 * <var>s</var> edge copies and <var>n</var> node copies.
		 */
		private final double[][][] edge_transition_likelihoods;
		/**
		 * Outside edge log-probabilties for death transitions 
		 * J[<var>u</var>][<var>s</var>][<var>n</var>] for 
		 * <var>s</var> edge copies and <var>n</var> parent copies.
		 */
		private final double[][][] edge_transition_outside;
		
		private Profile(int family_idx)
		{
			this.family_idx = family_idx;
			int num_nodes = tree.getNumNodes();
			this.node_likelihoods = new double[num_nodes][];
			this.edge_likelihoods = new double[num_nodes][];
			this.edge_outside = new double[num_nodes][];
			this.node_outside = new double[num_nodes][];
			this.sibling_likelihoods = new double[num_nodes][];
			this.edge_transition_likelihoods = new double[num_nodes][][];
			this.edge_transition_outside = new double[num_nodes][][];
			this.initDataStructures();
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
			double logp = getLogLossParameter(root);
			assert (logp==0.0);
			
			return LL;
		}		
		
		/**
		 * Array of posterior log-probabilities for death transitions.
		 * 
		 * @param node
		 * @return post[n][s] for n&rarr;s from parent
		 */
		private double[][] getLogDeathPosteriors(int node)
		{
			double[] K = getEdgeInside(node);
			double[] J = getEdgeOutside(node);
			
			double LL = Double.NEGATIVE_INFINITY;
			for (int s=0; s<J.length; s++)
			{
				double t = J[s]+K[s];
				LL = Logarithms.add(LL, t);
			}
			
			double[][] Jsn = edge_transition_outside[node];
			int Jlen = Jsn[0].length;
			
			double[][] post = new double[Jlen][];
//			double debug_sum = 0.0;
			
			for (int n=0; n<post.length; n++)
			{
				post[n] = new double[n+1];
				Arrays.fill(post[n], Double.NEGATIVE_INFINITY);
				for (int s=0; s<=n && s<K.length; s++)
				{
					double jsn = s<Jsn.length && n<Jsn[s].length?Jsn[s][n]:Double.NEGATIVE_INFINITY;
					post[n][s] = jsn + K[s]-LL;
//					debug_sum += Math.exp(post[n][s]);
				}
			}
//			System.out.println("#**SL.gLDP "+family_idx+"\t"+node+"\tsum "+debug_sum+"\tLL "+LL);
			return post;		
		}
		
		/**
		 * Array of posterior log-probabilities for birth transitions.
		 * 
		 * @param node
		 * @return post[n][s] for s&rarr;n towards child
		 */
		private double[][] getLogBirthPosteriors(int node)
		{
			double[] K = getEdgeInside(node);
			double[] J = getEdgeOutside(node);
			
			double LL = Double.NEGATIVE_INFINITY;
			for (int s=0; s<J.length; s++)
			{
				double t = J[s]+K[s];
				LL = Logarithms.add(LL, t);
			}
			int Clen = node_likelihoods[node].length;
			
			double[][] Ksn = edge_transition_likelihoods[node];
			
//			double debug_sum=0.0;
			
			double[][] post = new double[Clen][];
			for (int ell=0; ell<Clen; ell++)
			{
				post[ell] = new double[ell+1];		
				Arrays.fill(post[ell], Double.NEGATIVE_INFINITY);
				for (int s=0; s<=ell && s<J.length; s++)
				{
					post[ell][s] =  J[s] + Ksn[s][ell]-LL;
//					debug_sum += Math.exp(post[ell][s]);
				}
			}
//			System.out.println("#**SL.P.gLBP "+family_idx+"\t"+node+"\tsum "+debug_sum);
			
			return post;
		}
		
		/**
		 * Array of posterior log-probabilites for copy numbers at nodes.
		 * 
		 * @param node
		 * @return [log Prob(copies=i)]
		 * 
		 */
		public double[] getLogNodePosteriors(int node)
		{
			double[] C = getNodeInside(node);
			double[] B = getNodeOutside(node);
			
			assert (B.length==C.length);
			
			double LL = Double.NEGATIVE_INFINITY;
			for (int n=0; n<C.length; n++)
			{
				double t = B[n]+C[n];
				LL = Logarithms.add(LL, t);
			}
			double[] post = new double[C.length];
			for (int n=0; n<C.length; n++)
			{
				post[n] = B[n]+C[n]-LL;
			}
			return post;
		}
		/**
		 * Posterior log-expected value for node copies.
		 * 
		 * @param node
		 * @return log Exp(node copies)
		 */
		public double getLogNodeMean(int node)
		{
			double[] post = getLogNodePosteriors(node);
			// compute tails
			double log_tail = Double.NEGATIVE_INFINITY;
			int j=post.length;
			while (0<j)
			{
				--j;
				double x = post[j];
				post[j] = log_tail;
				log_tail = Logarithms.add(log_tail, x);
			}
			double log_mean = Logarithms.sum(post, post.length);
			return log_mean;
		}
		
		/**
		 * Posterior log-probabilities for edge copies: copies inherited from the parent. 
		 * 
		 * @param node
		 * @return [log Prob(edge copies=i)]
		 */
		public double[] getLogEdgePosteriors(int node)
		{
			double[] K = getEdgeInside(node);
			double[] J = getEdgeOutside(node);
			
			double LL = Double.NEGATIVE_INFINITY;
			for (int s=0; s<J.length; s++)
			{
				double t = J[s]+K[s];
				LL = Logarithms.add(LL, t);
			}
			double[] post = J.clone();
			for (int s=0;s<post.length; s++)
				post[s] += K[s]-LL;
			return post;
		}
		
		/**
		 * Posterior log-expected value for edge copies (inherited from parent).
		 * 
		 * @param node
		 * @return log Exp(edge copies)
		 */
		public double getLogEdgeMean(int node)
		{
			double[] post = getLogEdgePosteriors(node);
			// compute tails
			double log_tail = Double.NEGATIVE_INFINITY;
			int j=post.length;
			while (0<j)
			{
				--j;
				double x = post[j];
				post[j] = log_tail;
				log_tail = Logarithms.add(log_tail, x);
			}
			double log_mean = Logarithms.sum(post, post.length);
			return log_mean;
		}
		
		/**
		 * Posterior log-probabilities for copy number transitions: 
		 * T[n][m] = log Prob(parent copies =n, node copies=m)
		 * 
		 * @param node
		 * @return Tnm
		 */
		public double[][] getLogNodeTransitionPosteriors(int node)
		{
			double[] K = getEdgeInside(node);
			double[] J = getEdgeOutside(node);
			
			double LL = Double.NEGATIVE_INFINITY;
			for (int s=0; s<J.length; s++)
			{
				double t = J[s]+K[s];
				LL = Logarithms.add(LL, t);
			}
			
			
			double[][] Jsn = edge_transition_outside[node];
			double[][] Ksm = edge_transition_likelihoods[node];
			int Jlen = Jsn[0].length;
			int Klen = Ksm[0].length;
			
			double[][] Tnm = new double[Jlen][Klen];
			
			double debug_sum=0.0;
			
			for (int n=0; n<Jlen; n++)
			{
				for (int m=0; m<Klen; m++)
				{
					double t = Double.NEGATIVE_INFINITY;
					for (int s=0; s<=n && s<=m && s<Jsn.length && s<Ksm.length; s++)
					{
						Logarithms.add(t, Jsn[s][n]+Ksm[s][m]);
					}
					Tnm[n][m]=t-LL;
					debug_sum += Math.exp(Tnm[n][m]);
				}
			}
			
			return Tnm;
		}
		
		/**
		 * Array of log-tail-difference for copy births. 
		 * Summing the probabilities in this array gives the 
		 * posterior expectation log EXP(N-S) for N=node copies, S=edge copies.
		 * 
		 * @param node
		 * @return [log (Prob(node copies&gt;i)-Prob(edge copies&gt;i))]
		 */
		public double[] getLogBirthDifferenceTails(int node)
		{
			double[][] post = getLogBirthPosteriors(node);
			double[] N_S = logTailDifference(post);
			return N_S;
		}
		
		/**
		 * Array of log-tail-difference for copy deaths. 
		 * Summing the probabilities in this array gives the 
		 * posterior expectation log EXP(N-S) for N=parent copies, S=edge copies.
		 * 
		 * @param node
		 * @return [log (Prob(parent copies&gt;i)-Prob(edge copies&gt;i))]
		 */
		public double[] getLogDeathDifferenceTails(int node)
		{
			double[][] post = getLogDeathPosteriors(node);
			double[] N_S = logTailDifference(post);
			return N_S;
		}
		
		/**
		 * Calculating tail difference from posteriors 
		 * for birth or death. N denotes node copies, S denotes edge copies.
		 * we want N_S[k] = log Prob(N&gt;k and S&le;k) = log (Prob(N&gt;k)-Prob(S&gt;k)). 
		 * 
		 * N_S[k] = logsum (s=0..k; ell=k+1..max) Tls 
		 * 
		 * @param Tls
		 * @return [log (Prob(N&gt;k)-Prob(S&gt;k))]
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
					if (s<Tls[ell].length)
					{
						N_S[ell] = Logarithms.add(N_S[ell], tail); 
						tail = Logarithms.add(tail, Tls[ell][s]);
					}
				}
			}
			return N_S;
		}
		
		/**
		 * Intitialize or resets the data structures.
		 */
		private void initDataStructures()
		{
			int max_n = table.maxCopies(family_idx);
			int m = StraightLikelihood.this.getCalculationWidth(max_n);
			int[] profile = table.getFamilyProfile(family_idx);
			
//			System.out.println("#**SL.P.iDS "+family_idx+"\t"+Arrays.toString(profile));
			
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
				edge_transition_likelihoods[node] = null;
				edge_transition_outside[node] = null;
			}
		}		
		
		/**
		 * Retrieves the node likelihoods, and calculates all 
		 * inside likelihoods if not done yet.
		 * 
		 * @param node
		 * @return
		 */
		private double[] getNodeInside(int node)
		{
			if (edge_likelihoods[node]==null) // check if computed them yet 
			{
				computeInside();
			}
			double[] C = node_likelihoods[node];
			return C;
		}
		
		/**
		 * Retrieves the edge likelihoods, and calculates all 
		 * inside likelihoods if not done yet.
		 * 
		 */ 
		private double[] getEdgeInside(int node)
		{
			double[] K = edge_likelihoods[node];
			if (K==null)
			{
				computeInside();
				K = edge_likelihoods[node];
			}
			return K;
		}
		/**
		 * Retrieves the outside node likelihoods, and calculates all 
		 * outside likelihoods if not done yet.
		 * 
		 * @param node
		 * @return
		 */
		private double[] getNodeOutside(int node)
		{
			double[] B = node_outside[node];
			if (B==null) 
			{
				computeOutside();
				B = node_outside[node];
			}
			return B;
			
		}
		
		/**
		 * Retrieves the outside edge likelihoods, and calculates all 
		 * outside likelihoods if not done yet.
		 * 
		 * @param node
		 * @return
		 */
		private double[] getEdgeOutside(int node)
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
		 * Calculates the inside likelihoods in a postorder traversal of the nodes.
		 */
		private void computeInside()
		{
			for (int node=0; node<tree.getNumNodes(); node++)
			{
				this.computeNodeInside(node);
				edge_transition_likelihoods[node] = computeEdgeInside( node);
				double[] K = new double[edge_transition_likelihoods[node].length];
				for (int s=0; s<K.length; s++)
				{
					int n=s;
					double sum = edge_transition_likelihoods[node][s][n];
					++n;
					while (n<edge_transition_likelihoods[node][s].length)
					{
						sum = Logarithms.add(sum, edge_transition_likelihoods[node][s][n]);
						++n;
					}
					K[s] = sum;
				}
				if (TRACK_CALCULATIONS)
					System.out.println("#**SL.P.CI "+node+"\t("+K.length+")\tK "+Arrays.toString(K)+"\tKs "+Arrays.toString(edge_transition_likelihoods[node][0]));
				
				edge_likelihoods[node] = K;
			}
		}
		
		/**
		 * Calculates the outside likelihoods in a preorder traversal of the nodes.
		 */
		private void computeOutside()
		{
			int root = tree.getRoot();
			for (int node=root; node>=0; --node)
			{
				double[][] Jsn = edge_transition_outside[node] = computeEdgeOutside(node);
				double[] J = new double[Jsn.length];
				for (int s=0; s<J.length; s++)
				{
					int n = s;
					double sum = Jsn[s][n];
					++n;
					while (n<Jsn[s].length)
					{
						sum = Logarithms.add(sum, Jsn[s][n]);
						++n;
					}
					J[s] = sum;
				}
				edge_outside[node] = J;
				node_outside[node] = computeNodeOutside(node);
			}
		}
		
		/**
		 * Inside edge log-likelihoods, by copy birth transitions: 
		 * <var>s</var> edge copies (inherited from parent), <var>n</var>
		 * node copies.  
		 * 
		 * @param node
		 * @return K[s][n]
		 */
		private double[][] computeEdgeInside(int node)
		{
			double[][] Ksn; // the edge likelihoods that will be set 
			double[] C = node_likelihoods[node];
		
			if (tree.isRoot(node))
			{
				assert (getLogitLossParameter(node) == Double.POSITIVE_INFINITY);
				Ksn = new double[1][C.length];
			} else
			{
				Ksn=new double[C.length][C.length];
			}
			TransitionProbabilities Tn = transition_probs[node];
			
//			double r = getGainParameter(node);
//			assert (r!=0.0);
//			

			for (int s=0; s<Ksn.length; s++)
			{
				double[] Ks = Ksn[s];
				int n=s; 				
				Arrays.fill(Ks,0,n,Double.NEGATIVE_INFINITY);
				while (n<C.length)
				{
					Ks[n] = Tn.copybirth(s, n) + C[n];
//					System.out.println("#***SL.TP.cEI "+node+"\ts "+s+"\tn "+n+"\t"+Ksn[s][n]+"\tC "+C[n]+"\ttn "+Tn.copybirth(s, n));
					
					
					++n;
				}
				
//				if (s==0)
//				{
//					System.out.println("#***SL.P.cEI "+node+"\t"+Arrays.toString(Ks)+"\t"+Tn.copybirth(0, 0));
//				}
			}
			return Ksn;
		} // compute inside		

	
	
		/**
		 * Computes and sets inside node log-likelihoods {@link #node_likelihoods},
		 * as well as the complementary sibling log-likelihoods ({@link #sibling_likelihoods}
		 * for the children: 
		 * by <var>n</var>
		 * node copies.  
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
				
				if (TRACK_CALCULATIONS)
					System.out.println("#**SL.P.cNI "+node+"\tleaf n="+n+"\tC "+Arrays.toString(C));
			} else // ancestral node
			{
				Arrays.fill(C, 0.0);
				
				// compute stem likelihoods 
				int num_children = tree.getNumChildren(node);
				
				double[][] sib = new double[num_children][];
				for (int ci=0; ci<num_children; ci++)
				{
					int child = tree.getChild(node, ci);
					TransitionProbabilities Tc = StraightLikelihood.this.transition_probs[child];
					double[] C2 = new double[C.length];
					
					double[] K = edge_likelihoods[child];
					double[] terms = new double[K.length]; // temporary, reused
					for (int ell=0; ell<C2.length; ell++)
					{
						int k=0; 
						while(k<K.length && k<=ell)
						{
							terms[k] = Tc.copydeath(ell, k) + K[k]; // cNI C2[ell] = sum k=0..ell <K.length; ell fix 
							k++;
						}
						double z = Logarithms.sum(terms, k);
						C2[ell] = z;
						C[ell] += z; // product across children
					} // for ell
					sib[ci] = C2; // save it for later
				}
				
				if (TRACK_CALCULATIONS)
					System.out.println("#**DL.P.CNI "+family_idx+"\t"+node+"\t("+C.length+")\t"+Arrays.toString(C)+"\t// "+tree.toString(node)+"\t// "+StraightLikelihood.this.toString(node));
				
				
				
				// compute product of stem likelihoods across siblings 
				// -- in linear (and not quadratic) time with arity  
				
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
				}
			} // ancestral
		}

		/**
		 * Inside edge log-likelihoods, by copy death transitions: 
		 * <var>s</var> edge copies (inherited from parent), <var>n</var>
		 * parent copies.  
		 * 
		 * @param node
		 * @return J[<var>s</var>][<var>n</var>]
		 */
		private double[][] computeEdgeOutside(int node)
		{
			double[][] Jsn; // return value
			TransitionProbabilities Tn = transition_probs[node]; 
			if (tree.isRoot(node))
			{
				assert (getLogitLossParameter(node) == Double.POSITIVE_INFINITY);
				Jsn = new double[1][1];
				Jsn[0][0]=0.0; 
			} else
			{
				int parent = tree.getParent(node);
				double[] B = node_outside[parent];
				
				int Clen = node_likelihoods[node].length;
				Jsn = new double[Clen][B.length];
				for (int s=0; s<Jsn.length; s++)
				{
					double[] Js = Jsn[s];
					Arrays.fill(Js,Double.NEGATIVE_INFINITY);
					int n = s;
					while (n<B.length)
					{
						Js[n] = B[n] + sibling_likelihoods[node][n] + Tn.copydeath(n, s);
						++n;
					}
				}
			}
			return Jsn;
		}	

		/**
		 * Computes the outside node log-likelihoods {@link #node_outside}.
		 * by <var>n</var>
		 * node copies.  
		 * @param node
		 * @return B[<var>n</var>]
		 */
		private double[] computeNodeOutside(int node)
		{
			double[] C = node_likelihoods[node];
			double[] J = edge_outside[node]; // already set
			double[] B = new double[C.length]; // return value

			TransitionProbabilities Tn = transition_probs[node];
			
//			assert (getGainParameter(node)!=0.0);

			double[] terms =new double[J.length];
				
			for (int n=0; n<B.length; n++)
			{
				int s=0; 
				while (s<J.length && s<=n) 
				{
					terms[s] = J[s] + Tn.copybirth(s, n); 
					++s;
				}
				B[n] = Logarithms.sum(terms, s);
			}
			return B;
		}
	}
	
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
		sb.append("]");
		sb.append("::");
		sb.append(rates.toString(node));
		sb.append("::");
		sb.append("[lp ").append(getLogitLossParameter(node)).append("/rates:").append(rates.getLogLossParameter(node)-rates.getLogLossComplement(node));
		sb.append(",lq ").append(getLogitDuplicationParameter(node)).append("/rates:").append(rates.getLogDuplicationParameter(node)-rates.getLogDuplicationComplement(node));
		sb.append("]");
		
		
		return sb.toString();
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
        		StraightLikelihood SL = new StraightLikelihood(rates, empty);
        		SL.setCalculationWidthThresholds(threshold_width_absolute, threshold_width_relative);
        		Profile SP = SL.getProfile(0);
        		dLL = SP.getLogLikelihood();
        		cwidth = SL.getCalculationWidth(0);
        		
        	} else
        	{
        		Profile SP = this.getProfile(f);
        		Posteriors.Profile PP = post.getPosteriors(f);
        		dLL = SP.getLogLikelihood();
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
	
	
//	private void testComparePosteriors(PrintStream out, DirectLikelihood DL)
//	{
//		int nFam = table.getFamilyCount();
//		
//        for (int f=0; f<nFam; f++)		
//        {
//        	Profile SP = getProfile(f);
//        	DirectLikelihood.Profile DP = DL.getProfile(f);
//        	this.testComparePosteriors(out, SP, DP);
//        }
//		
//	}
//	
//	private void testComparePosteriors(PrintStream out, StraightLikelihood.Profile SP, DirectLikelihood.Profile DP)
//	{
//        int nNodes = tree.getNumNodes();
//        int nLeaves = tree.getNumLeaves();
//		
//        int family_idx = SP.family_idx;
//                
//    	for (int node=nLeaves; node<nNodes; node++ )
//    	{
//    		double[] direct_pN = DP.getNodePosteriors(node);
//    		double[] straight_pN = SP.getLogNodePosteriors(node);
//    		double straight_tN = 0.0;
//    		for (int n=0; n<straight_pN.length; n++)
//    		{
//    			straight_pN[n] = Math.exp(straight_pN[n]);
//    			straight_tN += straight_pN[n];
//    		}
//    		double direct_mN = DirectLikelihood.mean(direct_pN);
//    		double straight_mN = Math.exp(SP.getLogNodeMean(node));
//    		double diff_mN =  straight_mN-direct_mN;
//    		
//    		out.print(family_idx+"\tN\t"+node+"\t"+straight_mN+"("+DirectLikelihood.mean(straight_pN)+")/"+diff_mN+"\t"+direct_mN+"\ttot "+straight_tN);
//    		if (Math.abs(diff_mN)>1e-5)
//    		{
//	    		for (int n=0; n<straight_pN.length || n<direct_pN.length; n++)
//	    		{
//	    			double sp = n<straight_pN.length?straight_pN[n]:0.0;
//	    			double dp = n<direct_pN.length?direct_pN[n]:0.0;
//	    			
//	    			out.print("\t"+sp+"/"+dp);
//	    		}
//    		}
//    		out.println();
//    		
//    		double[] direct_pS = DP.getEdgePosteriors(node);
//    		double direct_mS = DirectLikelihood.mean(direct_pS);
//    		double[] straight_pS = SP.getLogEdgePosteriors(node);
//    		double straight_tS=0.0;
//    		for (int s=0; s<straight_pS.length; s++)
//    		{
//    			straight_pS[s] = Math.exp(straight_pS[s]);
//    			straight_tS += straight_pS[s];
//    		}
//    		double straight_mS = Math.exp(SP.getLogEdgeMean(node));
//    		double diff_mS = straight_mS-direct_mS;
//    		
//    		out.print(family_idx+"\tS\t"+node+"\t"+straight_mS+"("+DirectLikelihood.mean(straight_pS)+")/"+diff_mS+"\t"+direct_mS+"\ttot "+straight_tS);
//    		if (diff_mS>1e-5)
//    		{
//	    		for (int s=0; s<straight_pS.length || s<direct_pS.length; s++)
//	    		{
//	    			double sp = s<straight_pS.length?straight_pS[s]:0.0;
//	    			double dp = s<direct_pS.length?direct_pS[s]:0.0;
//	    			
//	    			out.print("\t"+sp+"/"+dp);
//	    		}
//    		}
//    		out.println();
//    		
//    		if (!tree.isRoot(node))
//    		{
//	    		double[] straight_tNu_Sv=SP.getLogDeathDifferenceTails(node);
//	    		double straight_mNu_Sv = 0.0;
//	    		for (int i=0; i<straight_tNu_Sv.length; i++)
//	    		{
//	    			double t = Math.exp(straight_tNu_Sv[i]);
//	    			straight_mNu_Sv += t;
//	    		}
//	    		double direct_Nu = DirectLikelihood.mean(DP.getNodePosteriors(tree.getParent(node)));
//	    		double direct_mNu_Sv = direct_Nu-direct_mS;
//	    		double diff_mNu_Sv = straight_mNu_Sv-direct_mNu_Sv;
//	    		out.println(family_idx+"\tNu_Sv\t"+node+"\t"+straight_mNu_Sv+"/"+diff_mNu_Sv+"\t"+direct_mNu_Sv);
//    		}
//    		double[] straight_tNv_Sv = SP.getLogBirthDifferenceTails(node);
//    		double straight_mNv_Sv = 0.0;
//    		for (int i=0; i<straight_tNv_Sv.length; i++)
//    		{
//    			double t = Math.exp(straight_tNv_Sv[i]);
//    			straight_mNv_Sv+=t;
//    		}
//    		double direct_mNv_Sv = direct_mN-direct_mS;
//    		double diff_mNv_Sv = straight_mNv_Sv-direct_mNv_Sv;
//    		out.println(family_idx+"\tNv_Sv\t"+node+"\t"+straight_mNv_Sv+"/"+diff_mNv_Sv+"\t"+direct_mNv_Sv);
//    	}
//	}
	
//	public static void main(String[] args) throws Exception
//	{
//		Class<?> our_class = java.lang.invoke.MethodHandles.lookup().lookupClass();
//		CommandLine cli = new CommandLine(args, our_class);
//        AnnotatedTable table = cli.getTable();
//        TreeWithRates rates = cli.getRates();
//
//        PrintStream out = System.out;        
//        StraightLikelihood SL = new StraightLikelihood(rates,table);
//    	int absolute = 3;
//    	double relative = 1.0;
//        if (cli.getOptionValue(CommandLine.OPT_TRUNCATE)!=null)
//        {
//        	absolute = cli.getOptionTruncateAbsolute();
//        	relative = cli.getOptionTruncateRelative();
//		}
//		out.println(CommandLine.getStandardHeader("Truncated computation (absolute,relative)="
//				+absolute+","+relative));
//		
//		SL.setCalculationWidthThresholds(absolute, relative);
//		
////		// checking likelihood calculations
////		Posteriors post = new Posteriors(rates, table);
////		post.setCalculationWidthThresholds(12, 3.0);
////		SL.testCompareLikelihoods(out, post);
//
//		// checking posterior calculations
//		DirectLikelihood DL = new DirectLikelihood(rates, table);
//		DL.setCalculationWidthThresholds(absolute, relative);
//		SL.testComparePosteriors(out, DL);
//		
//	}	
}
