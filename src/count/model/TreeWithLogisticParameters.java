package count.model;
/*
 * Copyright 2024 Mikl&oacute;s Cs&#369;r&ouml;s.
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

import count.io.CommandLine;
import count.matek.Logarithms;

/**
 * A layer above an underlying {@link count.model.TreeWithRates}; here we store 
 * distribution parameters on logistic and logarithmic scales. This is useful in numerical 
 * optimization where extreme parameter settings may be probed.  
 * 
 * 
 */
public class TreeWithLogisticParameters extends TreeWithRates
{
	private static final boolean FIX_INFINITE_GAIN=true;
	private static final double MINIMUM_LOGIT_LOSS = Logarithms.toLogit(1e-99); //  Logarithms.toLogit(Math.ulp(1.0)); // Double.NEGATIVE_INFINITY
	
	
	public TreeWithLogisticParameters(TreeWithRates rates)
	{
		this(rates, false);
	}
	
	public TreeWithLogisticParameters(TreeWithRates rates, boolean copy)
	{
		super(rates,copy);
		int n = rates.getTree().getNumNodes();
		logit_parameters = new double[n][4];
		log_parameters = new double[n][4];
		log_complements = new double[n][4];
		for (int node=0; node<n; node++)
			inferParametersFromRates(node);
	}	

	private final double[][] logit_parameters;
	private final double[][] log_parameters;
	private final double[][] log_complements;
	
//	public void initLogitParametersFromRates()
//	{
//		for (int node=0; node<logit_parameters.length; node++)
//			initLogitParametersFromRates(node);
//	}
	
	private void inferParametersFromRates(int node)
	{
		double logp = this.log_parameters[node][PARAMETER_LOSS] = Math.log(super.getLossParameter(node));
		double log1_p = this.log_complements[node][PARAMETER_LOSS]  =  Math.log(super.getLossParameterComplement(node));
		this.logit_parameters[node][PARAMETER_LOSS] = logp-log1_p;
		double logq = this.log_parameters[node][PARAMETER_DUPLICATION] =  Math.log(super.getDuplicationParameter(node));
		double log1_q = this.log_complements[node][PARAMETER_DUPLICATION]  =  Math.log(super.getDuplicationParameterComplement(node));
		this.logit_parameters[node][PARAMETER_DUPLICATION] = logq-log1_q;
		
		setLogitParameter(node,PARAMETER_RELATIVE_DUPLICATION, inferLogitRelativeRate(node));				
	}
	
	@Override
	public void initNodeParameters(int node)	
	{
		super.initNodeParameters(node);
		this.inferParametersFromRates(node);
	}
	
	@Override
	public double getLogLossParameter(int node)
	{
		return log_parameters[node][PARAMETER_LOSS];
	}
	
	@Override
	public double getLogLossComplement(int node)
	{
		return log_complements[node][PARAMETER_LOSS];
	}
	
	@Override
	public double getLogDuplicationParameter(int node)
	{
		return log_parameters[node][PARAMETER_DUPLICATION];
	}
	
	@Override
	public double getLogDuplicationComplement(int node)
	{
		return log_complements[node][PARAMETER_DUPLICATION];
	}
	
	@Override
	public double getLogRelativeRate(int node)
	{
		return log_parameters[node][PARAMETER_RELATIVE_DUPLICATION];
	}
	
	public double getLogRelativeComplement(int node)
	{
		return log_complements[node][PARAMETER_RELATIVE_DUPLICATION];
	}

	public double getLogitLossParameter(int node) { return logit_parameters[node][PARAMETER_LOSS];}
	public double getLogitDuplicationParameter(int node) { return logit_parameters[node][PARAMETER_DUPLICATION];} 
	public double getLogitRelativeRate(int node) {return logit_parameters[node][PARAMETER_RELATIVE_DUPLICATION];}
	
	@Override 
	public double getLossParameter(int node)
	{
		if (logit_parameters==null)
			return super.getLossParameter(node);
		else
			return Math.exp(getLogLossParameter(node));
	}
	
	@Override 
	public double getLossParameterComplement(int node)
	{
		if (logit_parameters==null)
			return super.getLossParameterComplement(node);
		else
			return Math.exp(getLogLossComplement(node));
	}
	
	@Override 
	public double getDuplicationParameter(int node)
	{
		if (logit_parameters==null)
			return super.getDuplicationParameter(node);
		else
			return Math.exp(getLogDuplicationParameter(node));
	}
	@Override 
	public double getDuplicationParameterComplement(int node)
	{
		if (logit_parameters==null)
			return super.getDuplicationParameterComplement(node);
		else
			return Math.exp(getLogDuplicationComplement(node));
	}
	
	public double getLogitTransientParameter(int node)
	{
		assert (getLogitDuplicationParameter(node)<=getLogitLossParameter(node));
		double log1_λ  = getLogRelativeComplement(node);
		double logit_σ = -getLogitLossParameter(node)-log1_λ;
		return logit_σ;
	}
	
	public double getLogTransientParameter(int node)		
	{
		assert (getLogitDuplicationParameter(node)<=getLogitLossParameter(node));
		return getLogLossComplement(node)-getLogDuplicationComplement(node);
	}
	
//	public double getLogLinearUniversalRatio(int node)
//	{
//		double log1_q = getLogDuplicationComplement(node);
//		
//		double log_q = getLogDuplicationParameter(node);
//		
//		// if q==0.0, then log_q=-infty, log1_q=0; log(-log1_q)=-infty
//		// if q==1.0, then log_q = 0, log1_q = -infty, log(-log1_q)=infty
//		
//		if (log1_q==0.0)
//		{
//			return 0.0; // =-ln(1-q) ~ q as q->0
//		} else
//		{
//			return log_q - Math.log(-log1_q);
//		}
//	}
	
	public double getLogGainParameter(int node, int gain_bound)
	{
		return this.getLogGainParameter(node, gain_bound, true);
	}
	/**
	 * Logarithm of the gain parameter by different policies 
	 * 
	 * 
	 * @param node
	 * @param gain_bound one of PARAMETER_LOSS, PARAMETER_DUPLICATION, PARAMETER_GAIN
	 * @param linear_gain whether linear or universal gain 
	 * @return
	 */
	public double getLogGainParameter(int node, int gain_bound, boolean linear_gain)
	{
		final double log_gain_param; // return value
		
		double log1_q = getLogDuplicationComplement(node);
//		double gr = getGainRate(node);
//		double log_gr = Math.log(gr);
		
		double log_gr = this.getLogGainRate(node);
		double gr = Math.exp(log_gr);
		
		if (getLogitDuplicationParameter(node)==Double.NEGATIVE_INFINITY) //(log1_q==0.0)
		{
			// q==0.0
			// gr=gamma 
			if (PARAMETER_LOSS == gain_bound)
			{
				// want gamma 
				log_gain_param = log_gr;
			} else
			{
				// want r 
				assert (PARAMETER_GAIN==gain_bound) || (PARAMETER_DUPLICATION == gain_bound);
				log_gain_param = log_gr + getLogLossParameter(node);
			}
		} else
		{
			// gr=kappa 
			 if (PARAMETER_DUPLICATION == gain_bound)
			{
				// want kappa
				log_gain_param = log_gr;
			} else 
			{
				assert (PARAMETER_GAIN==gain_bound) || (PARAMETER_LOSS == gain_bound);
				double log_r;
				
				if (linear_gain)
				{
					log_r = log_gr + getLogDuplicationParameter(node);
				} else
				{
					double loglog1_q =  Logarithms.logitToLogLogComplement(getLogitDuplicationParameter(node)); //  Math.log(-log1_q);
					log_r = log_gr+loglog1_q;
				}
				
				if (PARAMETER_GAIN==gain_bound) 
				{
					// want r 
					log_gain_param = log_r;
				} else
				{
					// want kappa = r/p 
					log_gain_param = log_r - getLogLossParameter(node);
					if (!linear_gain)
						throw new UnsupportedOperationException("Universal loss-bound gain is buggy.");
				}
			}
		}
		return log_gain_param;
	}
	
	
	private void setLogitParameter(int node, int parameter_idx, double x)
	{
		double logp, log1_p;
		if (0<=x)
		{
			logp = Logarithms.logitToLogValue(x);
			log1_p = logp - x;
		} else
		{
			log1_p = Logarithms.logitToLogComplement(x);
			logp = log1_p + x;
		}
		logit_parameters[node][parameter_idx] = x;
		log_parameters[node][parameter_idx] = logp;
		log_complements[node][parameter_idx] = log1_p;			
	}
	
	
	/**
	 * 
	 * @param node
	 * @param logit_p
	 * @param logit_q
	 * @param gain_param r or kappa
	 */
	public void setLogitLossDuplication(int  node, double logit_p, double logit_q, double gain_param)
	{
		setLogitParameter(node, PARAMETER_LOSS, logit_p);
		setLogitParameter(node, PARAMETER_DUPLICATION, logit_q);
		double logit_rel = inferLogitRelativeRate(node);
		setLogitParameter(node, PARAMETER_RELATIVE_DUPLICATION, logit_rel);
		copyParametersToRates(node, gain_param);
	}
	
    /**
     * Adjusts the edge length for a given loss probability, while keeping gain and duplication rates the same.
     * 
     * @param node
     * @param logit_p
     */
    public void setEdgeLengthForLogitLoss(int node, double logit_p)
    {
    	setLogitParameter(node, PARAMETER_LOSS, logit_p);
    	
    	double gain_param = getGainParameter(node);
    	copyParametersToRates(node, gain_param);
    }
	
	
	
	public void setLogitTransientRelativeDuplication(int node, double logit_σ, double logit_λ, double gain_param)
	{
		double log1_λ = Logarithms.logitToLogComplement(logit_λ);
		double x = -logit_σ-log1_λ;
		setLogitParameter(node, PARAMETER_LOSS, x);
		setLogitParameter(node, PARAMETER_RELATIVE_DUPLICATION, logit_λ);
		double log1_σ = Logarithms.logitToLogComplement(logit_σ); 
		double y = logit_λ + log1_σ;
		setLogitParameter(node, PARAMETER_DUPLICATION, y);
		
		if (x==Double.POSITIVE_INFINITY) // p=1; ln p/(1-p) = ln (1/0)
		{
			// loss p=1
			super.setEdgeLength(node, Double.POSITIVE_INFINITY);
			double μ = getLossRate(node);
			assert (μ!=0.0);
			assert (Double.isFinite(μ));
			double mul = Math.exp(log_parameters[node][PARAMETER_DUPLICATION]);
			// q/p=q
			super.setDuplicationRate(node, mul*μ);
		} else
		{ // not infinite length
			double delta;
			if (x==y) 
			{
				delta = 0.0;
			} else 
			{
				delta = Math.exp(log1_λ);
			}
			double μt;
			
			if (x==y) // logit_parameters[node][PARAMETER_RELATIVE_DUPLICATION]==Double.POSITIVE_INFINITY) // (x==y) // || delta==0.0)
			{
				// includes ==0.0 
				μt = Math.exp(x) ; // = p/(1-p)
			} else
			{
				if (y==Double.NEGATIVE_INFINITY) // q=0; ln q/(1-q)=ln (0/1)
				{			
					// gain-loss model	
					μt  = -getLogLossComplement(node);
				} else
				{
					double δμt = -Logarithms.logitToLogValue(logit_σ);
					μt  =  δμt/Math.exp(log1_λ);
				}
			}
			super.setEdgeLength(node, μt);
			super.setDuplicationRate(node, Math.exp(getLogRelativeRate(node)));
			
		}			
		// set gain rates
		if (gain_param==0.0)
		{
			super.setGainRate(node, 0.0);
		} else if (y==Double.NEGATIVE_INFINITY && x!=Double.NEGATIVE_INFINITY)
		{
			double gamma = gain_param*Math.exp(-getLogLossParameter(node));
			super.setGainRate(node, gamma);
		} else
		{
			super.setGainRate(node, gain_param);
		}
	}
	
	public void setLogitLossRelativeDuplication(int node, double logit_p, double logit_rel, double gain_param)
	{
		setLogitParameter(node, PARAMETER_LOSS, logit_p);
		setLogitParameter(node, PARAMETER_RELATIVE_DUPLICATION, logit_rel);
		double logit_q = Logarithms.mulLogit(logit_p, logit_rel);
		setLogitParameter(node, PARAMETER_DUPLICATION, logit_q);

		
		copyParametersToRates(node, gain_param);
	}
	public void setLogitLossRelativeDuplication(int v, double logit_p, double logit_λ, double log_gain,  int gain_bound)
	{
		this.setLogitLossRelativeDuplication(v, logit_p, logit_λ, log_gain, gain_bound, true);
	}

	public void setLogitLossRelativeDuplication(int v, double logit_p, double logit_λ, double log_gain,  int gain_bound, boolean linear_gain)
	{

		if (logit_p<MINIMUM_LOGIT_LOSS)
		{
			// the correction makes this rate model instance's parameters differ from the caller's intention
			// hopefully, just an intermediate step in line minimization
			
			double new_logit_p = MINIMUM_LOGIT_LOSS;

			// adjust lambda?
			if (gain_bound == PARAMETER_LOSS)
			{
				double log_r = log_gain + Logarithms.logitToLogValue(logit_p);
				double new_log_gain = log_r - Logarithms.logitToLogValue(new_logit_p);
				
				//System.out.println("#**TWLP.sLLRD "+v+"\tsmall logit_p "+logit_p+"\treset to "+ new_logit_p+"\tlogitlm "+logit_λ+"\tloggn "+log_gain+"\treset to "+new_log_gain);
				log_gain = new_log_gain;
			} else
			{
				//System.out.println("#**TWLP.sLLRD "+v+"\tsmall logit_p "+logit_p+"\treset to "+ new_logit_p+"\tlogitlm "+logit_λ+"\tloggn "+log_gain+"\tstays");
			}
			logit_p = new_logit_p;
		}
		
		double logit_q;
		
		if (logit_λ == Double.NEGATIVE_INFINITY) // Poisson
		{
			logit_q = Double.NEGATIVE_INFINITY;
			
			double log_r;
			if (gain_bound == PARAMETER_LOSS)
				log_r = log_gain + Logarithms.logitToLogValue(logit_p);
			else
				log_r = log_gain;
			
			if (!Double.isFinite(Math.exp(log_r))) // DEBUG
			{
				System.out.println("#**TWLP.sLLRD "+v+"\tlogr "+log_r+"\tloggn "+log_gain
				+"\tlogitp "+logit_p+"\tlogitq "+logit_q+"\tlogitlm "+logit_λ
				+"\tlog1_q "+Logarithms.logitToLogComplement(logit_q)
				+"\tlogq "+Logarithms.logitToLogValue(logit_q)
				+"\t// "+GLDParameters.paramName(gain_bound)+"/"+(linear_gain?"linear":"universal"));
			}
			
			
			
			log_gain = log_r;
			
			
		} else
		{
			logit_q = Logarithms.mulLogit(logit_p, logit_λ);

			// Polya
			double log_kappa;
			double log_r;
			if (gain_bound == PARAMETER_DUPLICATION)
			{
				log_kappa = log_gain;
				log_r = Double.NaN; // unused
			} else
			{
				if (gain_bound == PARAMETER_GAIN)
				{
					log_r = log_gain;
				}
				else 
				{
					assert (PARAMETER_LOSS == gain_bound);
					log_r = log_gain + Logarithms.logitToLogValue(logit_p);
				}
				
				if (linear_gain)
				{
					double log_q = Logarithms.logitToLogValue(logit_q);
					log_kappa = log_r - log_q;
					log_gain = log_kappa;
					
					
					double kappa = Math.exp(log_kappa);
						
					
					if (!Double.isFinite(kappa)) // DEBUG
					{
						System.out.println("#**TWLP.sLLRD "+v+"\tinfinite kappa\tlogr "+log_r+"\tlogkapa "+log_kappa+"\tloggn "+log_gain
						+"\tlogitp "+logit_p+"\tlogitq "+logit_q+"\tlogitlm "+logit_λ
						+"\tlog1_q "+Logarithms.logitToLogComplement(logit_q)
						+"\tlogq "+Logarithms.logitToLogValue(logit_q)
						+"\t// "+gain_bound+"/"+linear_gain);
					}
					if (FIX_INFINITE_GAIN && kappa==Double.POSITIVE_INFINITY && log_kappa!=Double.POSITIVE_INFINITY)
					{
						// set to a max finite value?
						kappa = 1024.0 / Math.ulp(1.0);
						log_kappa = Math.log(kappa);
						log_gain = log_kappa;
						
						// switch to Poisson
						logit_λ = Double.NEGATIVE_INFINITY;
						logit_q = Double.NEGATIVE_INFINITY;
						log_gain = log_r;
						
						// at this point, the caller and this instance have different parameters for this node 
						
						// DEBUG
						System.out.println("#**TWLP.sLLRD "+v+"\tinfinite kappa fixed: kappa "+kappa+"\tlogkapa "+log_kappa+"\tlogitq "+logit_q+"\tlog_gain "+log_gain+"/gain "+Math.exp(log_gain));
					} else if (kappa==Double.POSITIVE_INFINITY)
					{
						System.out.println("#**TWLP.sLLRD "+v+"\tinfinite kappa not fixed");
					}
					
				} else
				{
					// "universal gain" : dubious setting 
					double loglog1_q = Math.log(-Logarithms.logitToLogComplement(logit_q));
					if (loglog1_q==Double.NEGATIVE_INFINITY)
						log_kappa = log_r-Logarithms.logitToLogValue(logit_q); // same for small q
					else
						log_kappa = log_r - loglog1_q;
					if (!Double.isFinite(Math.exp(log_kappa))) // DEBUG
					{
						System.out.println("#**TWLP.sLLRD "+v+"\tlogr "+log_r+"\tlogkapa "+log_kappa+"\tloggn "+log_gain
						+"\tll1q "+loglog1_q
						+"\tlogitp "+logit_p+"\tlogitq "+logit_q+"\tlogitlm "+logit_λ
						+"\tlog1_q "+Logarithms.logitToLogComplement(logit_q)
						+"\tlogq "+Logarithms.logitToLogValue(logit_q)
						+"\t// "+gain_bound+"/"+linear_gain);
					}
					throw new UnsupportedOperationException("Universal loss-bound gain is buggy.");

				}
			}
			
		} // Polya
		
		setLogitParameter(v, PARAMETER_LOSS, logit_p);
		setLogitParameter(v, PARAMETER_RELATIVE_DUPLICATION, logit_λ);
		setLogitParameter(v, PARAMETER_DUPLICATION, logit_q);
		
		
		double gain_param = Math.exp(log_gain);

		if (!Double.isFinite(gain_param)) // DEBUG
		{
			System.out.println("#**TWLP.sLLRD "+v
			+"\tgain_param "+gain_param
			+"\tloggn "+log_gain
			+"\tlogitp "+logit_p+"\tlogitq "+logit_q+"\tlogitlm "+logit_λ
			+"\tlog1_q "+Logarithms.logitToLogComplement(logit_q)
			+"\tlogq "+Logarithms.logitToLogValue(logit_q)
			+"\t// "+gain_bound+"/"+linear_gain);
		}
		
		copyParametersToRates(v, gain_param);
	}
	
	/**
	 * 
	 * @param node
	 * @param gain_param r or kappa
	 */
	private void copyParametersToRates(int node, double gain_param)
	{
		double x= logit_parameters[node][PARAMETER_LOSS];
		double y = logit_parameters[node][PARAMETER_DUPLICATION];
		if (x==Double.POSITIVE_INFINITY) // p=1; ln p/(1-p) = ln (1/0)
		{
			// loss p=1
			super.setEdgeLength(node, Double.POSITIVE_INFINITY);
			double μ = getLossRate(node);
			assert (μ!=0.0);
			assert (Double.isFinite(μ));
			double mul = Math.exp(log_parameters[node][PARAMETER_DUPLICATION]);
			// q/p=q
			super.setDuplicationRate(node, mul*μ);
		} else
		{ // not infinite length
			
			
			double delta;
			double log_delta;
			if (x==y) 
			{
				delta = 0.0;
				log_delta = Math.log(delta);
			}
			else if (y<=x)
			{
//					double g = -Math.expm1(y-x);
//					if (0<=y)
//					{
//						double h = 1+Math.exp(-y);
//						delta = Math.exp(-y)*g/h;
//					} else
//					{
//						double h = 1+Math.exp(y);
//						delta = g/h;
//					}
				log_delta = Logarithms.logLogitRatioComplement(x, y); 
				delta = Math.exp(log_delta);
			} else // x<y, p<q
			{
				log_delta = Logarithms.logLogitRatioComplement(y, x); 
				delta = Math.exp(log_delta);
//					double g = -Math.expm1(x-y);
//					if (0<=x)
//					{
//						double h = 1+Math.exp(-x);
//						delta = Math.exp(-x)*g/h;
//					} else
//					{
//						double h = 1+Math.exp(x);
//						delta = g/h;
//					}
			}
			double μt,λt;
			
			if (x==y) // logit_parameters[node][PARAMETER_RELATIVE_DUPLICATION]==Double.POSITIVE_INFINITY) // (x==y) // || delta==0.0)
			{
//					System.out.println("#**LP.TWL.cPTR "+node+"\tx "+x+"\ty "+y+"\tlogd "+log_delta+"\tdelta "+delta
//								+"\tpard "+log_complements[node][PARAMETER_RELATIVE_DUPLICATION]+"\tparr "+log_parameters[node][PARAMETER_RELATIVE_DUPLICATION]);
				// includes ==0.0 
				μt = λt = Math.exp(x) ; // = p/(1-p)
			} else
			{
				if (y==Double.NEGATIVE_INFINITY) // q=0; ln q/(1-q)=ln (0/1)
				{
					// gain-loss model	
					μt  = -log_complements[node][PARAMETER_LOSS];
					λt = 0.0;
				} else
				{
					if (y<=x)
					{
						//
						// want (-ln (1-p)/(1-q))/(1-lambda)
						
						double δμt = Logarithms.logLogitComplementRatio(x, y);
						double mulλ = Math.exp(logit_parameters[node][PARAMETER_RELATIVE_DUPLICATION]);
						μt = δμt/delta;
						λt = δμt*mulλ; 
						
//						if (!Double.isFinite(μt) || !Double.isFinite(λt) || μt==0.0)
//						{
//							System.out.println("#**LP.TWL.cPTR "+node+"\tdmut "+δμt+"\tx "+x+"\ty "+y+"\tlogd "+log_delta+"\tmul "+mulλ+"\tlogitr "+logit_parameters[node][LikelihoodParametrized.PARAMETER_RELATIVE_RATE]);
//						}
					} else
					{
						double δλt = Logarithms.logLogitComplementRatio(y, x);
						double mulμ = Math.exp(logit_parameters[node][PARAMETER_RELATIVE_DUPLICATION]);
						λt  = δλt/delta;
						μt = λt*mulμ; 
					}
				}
			} // unequal rates 

			assert (getLossRate(node)==1.0);
			// standard scaling
			
			super.setEdgeLength(node, μt);
			if (y<=x)
			{
				// q<p 
				super.setDuplicationRate(node, Math.exp(log_parameters[node][PARAMETER_RELATIVE_DUPLICATION]));
			} else
			{
				super.setDuplicationRate(node, Math.exp(-log_parameters[node][PARAMETER_RELATIVE_DUPLICATION]));
			}
		} // not infinite length

		// set gain rates
		if (gain_param==0.0)
		{
			super.setGainRate(node, 0.0);
		} else if (y==Double.NEGATIVE_INFINITY && x!=Double.NEGATIVE_INFINITY)
		{
			double gamma = gain_param*Math.exp(-log_parameters[node][PARAMETER_LOSS]);
			
			if (!Double.isFinite(gamma))
			{
				System.out.println("#***TWLP.cPTR node "+node+"\tx "+x+"\tgainp "+gain_param+"\tlogp "+log_parameters[node][PARAMETER_LOSS]+"\tgm "+gamma);
			}
			
			
			assert (Double.isFinite(gamma));
			
			super.setGainRate(node, gamma);
		} else
		{
			assert (Double.isFinite(gain_param));
			
			super.setGainRate(node, gain_param);
		}
		
		
//			{ // DEBUG
//				double dp = Math.exp(this.getLogLossParameter(node))-super.getLossParameter(node);
//				double dq = Math.exp(this.getLogDuplicationParameter(node))-super.getDuplicationParameter(node);
//				System.out.println("#**LP.TWL.cPTR "+node
//						+"\tdp " +dp+"\tdq "+dq  
//						+"\tp "+super.getLossParameter(node)+"\tq "+super.getDuplicationParameter(node)+"\t"+this.toString(node));
//			}
		
	}

	private double inferLogitRelativeRate(int node)
	{
		double x = logit_parameters[node][PARAMETER_LOSS];
		double y = logit_parameters[node][PARAMETER_DUPLICATION];
		double logitR;
		
		if (y<=x)
		{
			logitR = Logarithms.logitParameterRatio(x, y);
		} else
		{
			logitR = Logarithms.logitParameterRatio(y, x);			
		}
		return logitR;
	}
	
	@Override
	public void setEdgeLengthForLossParameter(int node, double p, double one_minus_p)
	{
		super.setEdgeLengthForLossParameter(node, p, one_minus_p);
		this.inferParametersFromRates(node);
	}
	
	/**
	 * Infers also the logistic parameters from the rates
	 */
	@Override
	public void setEdgeLength(int node, double length)
	{
		super.setEdgeLength(node, length);
		this.inferParametersFromRates(node);
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
	@Override
    public void setRates(int node, double length, double gain_rate, double loss_rate, double dup_rate) 
    {
    	super.setEdgeLength(node, length);
    	super.setLossRate(node, loss_rate);
    	super.setDuplicationRate(node, dup_rate);
    	super.setGainRate(node, gain_rate);
		this.inferParametersFromRates(node);
    }
	
	
	public void setMaxGainDuplicationRate(double max_kappa)
	{
		double log_max =Math.log(max_kappa);
		System.out.println("#**TWLP.sMGDR maxkapa "+max_kappa+"\tlog "+log_max);
		int num_nodes = getTree().getNumNodes();
		for (int u=0; u<num_nodes; u++)
		{
			
			if (! Double.isInfinite(getLogitDuplicationParameter(u)) ) // including q=0, q=1 (the latter would be a problem anyway)
			{
				double log_kappa = this.getLogGainParameter(u);
				if (log_max < log_kappa)
				{
					// go to Poisson instead
					double log_q = this.getLogDuplicationParameter(u);
					double log_r = log_kappa + log_q;
					
					String old_node_str = this.toString(u);
					
					double logit_p = getLogitLossParameter(u);
					this.setLogitLossDuplication(u, logit_p, Double.NEGATIVE_INFINITY, Math.exp(log_r));
					System.out.println("#**TWLP.sMGDR node "+u+"\tlogkapa "+log_kappa+"\tlogq "+log_q+"\tto Poisson w/logr "+log_r+"\twas "+old_node_str+"\tnow "+this.toString(u));
				}
			}
		}
		
		
	}	
	
	
	@Override
	public String toString(int node)
	{
		StringBuilder sb = new StringBuilder(super.toString(node));
		sb.append("/logitp ").append(getLogitLossParameter(node))
		.append(", logitq ").append(getLogitDuplicationParameter(node)).append("(").append(Math.log(super.getDuplicationParameter(node))-Math.log(super.getDuplicationParameterComplement(node))).append(")")
		.append(", loggn ").append(this.getLogGainParameter(node, PARAMETER_DUPLICATION, true));
		return sb.toString();
	}
	
	/**
	 * Test code: convert large gain/small duplication to Poisson (with <tt>-maxgain</tt> <var>maxkappa</var>)
	 * 
	 * @param args
	 * @throws Exception
	 */
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
        TreeWithLogisticParameters lrates;
        if (rates instanceof TreeWithLogisticParameters)
        	lrates = (TreeWithLogisticParameters) rates;
        else
        	lrates = new TreeWithLogisticParameters(rates,false);
        
        
        String maxgain = cli.getOptionValue("max"+CommandLine.OPT_GAIN);
        if (maxgain == null)
        {
        } else
        {
        	double max = Double.parseDouble(maxgain);
        	lrates.setMaxGainDuplicationRate(max);
        	
        }
    	out.println(count.io.RateVariationParser.printRates(zeb));        	
        
		if (out != System.out)
		{
			if (out.checkError())
				throw new java.io.IOException("Write failed.");
			out.close();
		}
        
	}
	
}