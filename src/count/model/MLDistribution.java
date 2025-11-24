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

import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LOSS;
import static count.model.GLDParameters.PARAMETER_DUPLICATION_COMPLEMENT;
import static count.model.GLDParameters.PARAMETER_LOSS_COMPLEMENT;
//import static count.model.GLDParameters.PARAMETER_LENGTH;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.function.DoubleFunction;

import count.Count;
import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.ds.UniqueProfileTable;
import count.io.CommandLine;
import count.io.NewickParser;
import count.matek.FunctionMinimization;
import count.matek.Functions;
import count.matek.Logarithms;

import static count.io.CommandLine.OPT_EPS;
import static count.io.CommandLine.OPT_GAIN;
import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_OUTPUT;
import static count.io.CommandLine.OPT_ROUNDS;
import static count.io.CommandLine.OPT_TRUNCATE;

import static count.io.CommandLine.OPT_MODEL_ROOT_PRIOR;
import static count.io.CommandLine.OPT_REINIT;

/**
 * Maximum likelihood for one-class model 
 * using distribution parameters (loss probability, duplication probability and gain rate). 
 * 
 * @author csuros
 *
 */
public class MLDistribution extends ML //implements Count.UsesThreadpool // via Gradient
{
//	public static final int DEFAULT_OPT_ROUNDS = 100;
//	public static final double DEFAULT_OPT_EPS = 1e-2;
 
	
	public static boolean PRINT_OPTIMIZATION_MESSAGES = false;
	
	public static final double MAX_PROB_LOSS  =  MAX_PROB_NOT1; //0.92;
	public static final double MAX_GAIN_RATE = 33.0; // larger value is better for numerical optimization; with small values one may get stuck on the boundary 
	
	
	private static int REINIT_EXTREME_EDGES = 0; // =0 no; =1 edge only; =3 with parent & sibling; =5 with parent & sibling & children
	private static boolean DEBUG_GRADIENT = false;
	
	private static int DEFAULT_TRUNCATE_ABSOLUTE = 6;
	private static double DEFAULT_TRUNCATE_RELATIVE = 1.0;
	private static boolean AUTO_TRUNCATE = true;
	
	
	private static final boolean PURE_LOGISTIC = true;
	
	
	private static final boolean DEBUG_LOGISTIC = true;
	
	
	public MLDistribution(TreeWithRates rates, ProfileTable table)
	{
		this(rates, table, optimizableParams(rates));
		for (int node=0; node<rates.getTree().getNumNodes(); node++)
		{
			if (gradient.factory.getExtinction(node)==1.0)
			{
				fixGain(node,true);
				fixLoss(node,true);
				fixDuplication(node,true);
				// DEBUG
				System.out.println("#**MLD() fix "+node+"\t// "+rates.toString(node));
			}
		}
	}
	
	public MLDistribution(TreeWithRates rates, ProfileTable table, boolean[] optimize_parameter)
	{
		UniqueProfileTable utable;
		if (table instanceof UniqueProfileTable)
			utable = (UniqueProfileTable) table;
		else utable = new UniqueProfileTable(table);
		
		if (PURE_LOGISTIC && !(rates instanceof TreeWithLogisticParameters)) {
			rates = new TreeWithLogisticParameters(rates, false); // hard link
		}
		
		if (PURE_LOGISTIC) this.gradient = new Gradient(new LikelihoodParametrized(rates, utable));
		else this.gradient = new Gradient(rates,utable);
		
		this.rates = gradient.factory.rates; //.parameterCache();
		if (PURE_LOGISTIC) assert (rates instanceof TreeWithLogisticParameters); // high-precision base model

		this.optimize_parameters = optimize_parameter;
		this.node_parameters = new double[rates.getTree().getNumNodes()][];
		this.distribution_params = new ArrayList<>();
		
		//		this.table = table;
//		int parameter_count = 0; 
//		for (boolean opt: optimize_parameter)
//			if (opt) parameter_count++;
//		parameter_index = new int[parameter_count];
//		int i = optimize_parameter.length;
//		int j = parameter_count;
//		while (i>0)
//		{
//			--i;
//			if (optimize_parameter[i])
//			{
//				--j;
//				parameter_index[j] = i;
//			}
//		}
		
		
		
//		System.out.println("#*MLD.init "+gradient.factory.getCorrectedLL());
	}
	private final TreeWithRates rates;
//	private final ProfileTable table; 
	
	private final Gradient gradient;
	
	private boolean optimization_by_quasiNewton = true; // if false, use conjugate gradient (slower)
	
	private boolean do_gradient_descent = false; // if true, start with a few rounds of gradient descent (slower) 
	
	private boolean do_em = false;
	
	private boolean is_duprate_bounded = true; // if bounded, <1.0 is enforced
	
	/**
	 * Buggy with other settings 
	 */
	private int gain_parameter_type = PARAMETER_DUPLICATION;
	
	
	private boolean track_complements = false; // if true, high-precision calculation is attempted (numerically buggy) 
	
	
	private TreeWithLogisticParameters getRates()
	{
		assert (PURE_LOGISTIC);
		return (TreeWithLogisticParameters) rates;
	}
	
	public void setGainParameterType(int gain_type)
	{
		if (gain_type != PARAMETER_LOSS && gain_type != PARAMETER_DUPLICATION && gain_type != PARAMETER_GAIN)
			throw new IllegalArgumentException("setGainParameter called with unrecognized parameter type "+gain_type);
		this.gain_parameter_type = gain_type;
	}
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
	public void setCalculationWidth(int absolute, double relative)
	{
		gradient.setCalculationWidthThresholds(absolute, relative);
	}
	
	public void setMinimumObservedCopies(int m)
	{
		gradient.setMinimumObservedCopies(m);
	}
	
	
//	private int[] parameter_index;
//	
	private final double[][] node_parameters;
	
	private final List<ModelParameter> distribution_params;
//	private List<Double> optimization_history=null;
	private boolean[] optimize_parameters;
	
	public void fixGain(int node, boolean not_optimized)
	{
		optimize_parameters[3*node+PARAMETER_GAIN]=!not_optimized;
	}
	public void fixLoss(int node, boolean not_optimized)
	{
		optimize_parameters[3*node+PARAMETER_LOSS]=!not_optimized;
		// DEBUG
		System.out.println("#**MLD.fL "+node+"\topt "+optimize_parameters[3*node+PARAMETER_LOSS]);
	}
	public void fixDuplication(int node, boolean not_optimized)
	{
		optimize_parameters[3*node+PARAMETER_DUPLICATION]=!not_optimized;
	}
	
	@Override
	public void fixNodeParameters(int node, boolean do_not_optimize)
	{
//		// DEBUG
//		System.out.println("#**MLD.fNP "+node+"\t"+do_not_optimize);
		fixGain(node, do_not_optimize || rates.getGainParameter(node)==0.0);
		fixLoss(node, do_not_optimize || rates.getLossParameter(node)==1.0);
		fixDuplication(node, do_not_optimize || rates.getDuplicationParameter(node)==0.0);
	}
	
	
	private void initModelParameters()
	{
		int nfam = gradient.getTotalFamilyCount();
		double prob_small = 0.5/nfam;
//		double prob_big_threshold = 1.0-prob_small_threshold;
		
		int num_nodes = rates.getTree().getNumNodes();
		distribution_params.clear();
		for (int node=0; node<num_nodes; node++)
		{
			boolean optimize_gain = optimize_parameters[3*node+PARAMETER_GAIN];
			boolean optimize_loss = optimize_parameters[3*node+PARAMETER_LOSS];
			boolean optimize_duplication = optimize_parameters[3*node+PARAMETER_DUPLICATION];
			if (optimize_gain || optimize_loss || optimize_duplication)
			{
				double[] params = new double[5]; // not 3 , to keep complements
				
				
				params[PARAMETER_GAIN] = rates.getGainParameter(node);
				
//				if (PURE_LOGISTIC)
				{
					TreeWithLogisticParameters lrates = (TreeWithLogisticParameters )rates;
					double logit_p = params[PARAMETER_LOSS] = lrates.getLogitLossParameter(node);
					double logit_q = lrates.getLogitDuplicationParameter(node);
					if (is_duprate_bounded)
					{
						if (logit_p<logit_q)
						{
							double smaller_logit_q  = params[PARAMETER_LOSS]*MAX_PROB_NOT1;
							String old_node = lrates.toString(node);

							lrates.setLogitLossDuplication(node, params[PARAMETER_LOSS], smaller_logit_q, params[PARAMETER_GAIN]);
							System.out.println("#***MLD.iMP dup rate should be <= 1.0; reset to "+smaller_logit_q+"\tnow "+lrates.toString(node)+"\twas "+old_node);
							
							logit_q = smaller_logit_q;
						}
						params[PARAMETER_DUPLICATION] = lrates.getLogitRelativeRate(node);
					} else
					{
						params[PARAMETER_DUPLICATION] = logit_q;
					}
					// unused complements bc logistic scale has complete info 
					
					
					if (gain_parameter_type != PARAMETER_DUPLICATION)
					{
						// want to use gamma or r 
						double gain_param = Math.exp(lrates.getLogGainParameter(node, gain_parameter_type));
						params[PARAMETER_GAIN] = gain_param;
					}					
//				} else
//				{
//					params[PARAMETER_LOSS] = rates.getLossParameter(node);
//					params[PARAMETER_LOSS_COMPLEMENT] = rates.getLossParameterComplement(node);
//					if (is_duprate_bounded)
//					{
//						params[PARAMETER_DUPLICATION] = rates.getDuplicationRate(node);
//						params[PARAMETER_DUPLICATION_COMPLEMENT] = rates.getRateGap(node);
//					} else
//					{
//						params[PARAMETER_DUPLICATION] = rates.getDuplicationParameter(node);
//						params[PARAMETER_DUPLICATION_COMPLEMENT] = rates.getDuplicationParameterComplement(node);
//					}
				}
				node_parameters[node] = params;
				
				
				// add them 
				if (optimize_gain)
				{
					distribution_params.add(newLogistic(new GainParameter(node), MAX_GAIN_RATE)); //10.0*MAX_RATE));
				}
				if (optimize_duplication) 
				{	
					if (PURE_LOGISTIC)
					{
						if (is_duprate_bounded)
						{
							distribution_params.add(new LogisticDuplicationRate(node));
						} else
						{
							distribution_params.add(new LogisticDuplicationParameter(node));
						}
					} else
					{
						if (is_duprate_bounded)
						{
	//						distribution_params.add(newLogistic(new DuplicationRate(node), 1.0) );
							if (track_complements)
							{
								distribution_params.add(new BoundedLogistic(new DuplicationRate(node)));
							} else
							{
								distribution_params.add(newLogistic(new DuplicationRate(node), MAX_PROB_NOT1*rates.getLossRate(node)));
							}
						}
						else
						{
							if (track_complements)
							{
								distribution_params.add(new BoundedLogistic(new DuplicationParameter(node)));
							} else
							{
								distribution_params.add(newLogistic(new DuplicationParameter(node), MAX_PROB_NOT1));
	//							distribution_params.add(bracketedLogistic(new DuplicationRate(node), 1.0/8192.0));
							}
						}
					}
					
//					ModelParameter dpar = new Logistic(new DuplicationParameter(node), 1.0-1e-9);
//					distribution_params.add(dpar);
				}
				if (optimize_loss) 
				{
					ModelParameter lpar;
					if (PURE_LOGISTIC)
					{
						lpar = new LogisticLossParameter(node);
					} else
					{
						if (track_complements && !is_duprate_bounded)
						{
							lpar = new BoundedLogistic(new LossParameter(node));
						} else
						{
							lpar = newLogistic(new LossParameter(node), MAX_PROB_LOSS);
						}
					}
					// lpar = bracketedLogistic(new LossParameter(node),prob_small );
					distribution_params.add(lpar);
//					// DEBUG
//					System.out.println("#**MLD.iMP "+node+"\tloss "+lpar);
				}
			} else
			{
				node_parameters[node]=null;
			}
		}
		
//		copyNodeParametersToModel(); // since we may have changed them to squeeze into 0..max
		gradient.computeParameters();
//		setParameterValues(getParameterValues());
	}
	
	/**
	 * Sets the actual rates at all nodes.
	 */
	private void copyNodeParametersToModel()
	{
		for (int node=0; node<node_parameters.length; node++)
		{
			if (node_parameters[node]!=null)
			{
//				if (PURE_LOGISTIC)
				{
					double logit_loss_par = node_parameters[node][PARAMETER_LOSS];
					// logit_dup_par: relative rate or duplication rate 
					double logit_dup_par = node_parameters[node][PARAMETER_DUPLICATION];
					// gain_par: duplication-gain kappa or gain intensity r 
					double gain_par = node_parameters[node][PARAMETER_GAIN];
					
					if (gain_parameter_type == PARAMETER_GAIN)
					{
						if (logit_dup_par == Double.NEGATIVE_INFINITY)
						{
							// nothing to do, we have r already
						} else
						{
							double r = gain_par;
							// want kappa = r/q 
							double logit_q;
							if (is_duprate_bounded)
							{
								logit_q = Logarithms.mulLogit(logit_loss_par, logit_dup_par);
							} else
							{
								logit_q = logit_dup_par;
							}
							gain_par = r/Math.exp(Logarithms.logitToLogValue(logit_q));
						}
					} else if (gain_parameter_type == PARAMETER_LOSS)
					{
						double gamma = gain_par;
						if (logit_dup_par == Double.NEGATIVE_INFINITY)
						{
							// Poisson: want r = p*gamma
							gain_par = gamma * Math.exp(Logarithms.logitToLogValue(logit_loss_par));
						} else
						{
							// Polya want kappa = gamma * p/q 
							if (is_duprate_bounded)
							{
								double logit_lambda = logit_dup_par;
								gain_par = gamma/Math.exp(Logarithms.logitToLogValue(logit_lambda));
							} else
							{
								double log_lambda = Logarithms.logitToLogValue(logit_dup_par)-Logarithms.logitToLogValue(logit_loss_par);
								gain_par = gamma/Math.exp(log_lambda);
							}
						}
					} else
					{
						assert (gain_parameter_type == PARAMETER_DUPLICATION);
						// nothing to do 
					}
					
					TreeWithLogisticParameters lrates = getRates();
					if (is_duprate_bounded)
					{
						lrates.setLogitLossRelativeDuplication(node, logit_loss_par, logit_dup_par, gain_par);
					} else
					{
						lrates.setLogitLossDuplication(node, logit_loss_par, logit_dup_par, gain_par);
					}
//				} else
//				{
//					double κ = node_parameters[node][PARAMETER_GAIN];
//					double p = node_parameters[node][PARAMETER_LOSS];
//					
//					if (is_duprate_bounded)
//					{
//						double λ = node_parameters[node][PARAMETER_DUPLICATION];
//						if (track_complements)
//						{
//							double λcomp = node_parameters[node][PARAMETER_DUPLICATION_COMPLEMENT];
//							rates.setDuplicationRate(node, λ, λcomp);
//							double p_1 = node_parameters[node][PARAMETER_LOSS_COMPLEMENT];
//							rates.setEdgeLengthForLossParameter(node, p, p_1);
//						} else
//						{
//							rates.setDuplicationRate(node, λ);
//							rates.setEdgeLengthForLossParameter(node, p, 1.0-p);
//						}
//						rates.setGainRateForLossParameter(node, κ, p);
//					} else
//					{
//						double q = node_parameters[node][PARAMETER_DUPLICATION];
//						if (track_complements)
//						{
//							double p_1 = node_parameters[node][PARAMETER_LOSS_COMPLEMENT];
//							double q_1 = node_parameters[node][PARAMETER_DUPLICATION_COMPLEMENT];
//							if (q==1.0)
//							{
//								System.out.println("#**MLD.cNP node "+node+"\tp" +p+"\t1-"+p_1+"\tq "+q+"\t1-"+q_1+"\t"+rates.toString(node));
//							}
//							
//							rates.setDuplicationLossRates(node, p, p_1, q, q_1);
//							rates.setGainRateForLossParameter(node, κ, p);
//						} else
//						{
//							rates.setParameters(node, κ, p, q);	
//						}
//	//				System.out.println("#**MLD.cNPTM "+node+"\t"+Arrays.toString(node_parameters[node])+"\t"+rates.toString(node));
//					}
				}
			}
		}
	}
	
	@Override
	public int getModelParameterCount()
	{
		return distribution_params.size();
	}
	
	
	public double[] getParameterValues()
	{
		double[] x = new double[distribution_params.size()];
		for (int j=0;j<x.length; j++)
			x[j] = distribution_params.get(j).get();

		return x;
	}
		
	public void setParameterValues(double[] x)
	{
		for (int j=0; j<distribution_params.size(); j++)
		{
			ModelParameter P = distribution_params.get(j);
			P.set(x[j]);
		}
		copyNodeParametersToModel();
		gradient.computeParameters();
	}
	
	private int calls_optFunc=0;
	
	public double optDistrFunc(double[] x)
	{
		setParameterValues(x);
		double LL = gradient.getCorrectedLL();
		calls_optFunc++;
//		System.out.println("#*MLD.oDF "+calls_optFunc+"\t"+LL); //+"\t"+rates.toString(24));

		
		return -LL;
	}
	
	
	private int calls_optDiff = 0;
	
	private double[] parameterGradient(double[] distribution_gradient)
	{
		double[] D = new double[distribution_params.size()];
		for (int j=0; j<D.length; j++)
		{
			D[j] = -distribution_params.get(j).dL(distribution_gradient);
//			System.out.println("#**MLD.pG "+j+"\t"+distribution_params.get(j)+"\tD "+D[j]);
			
		}
		return D;
	}
	
	private double[] parameterGradient()
	{
		double[] dLL = gradient.getCorrectedGradient();
		dLL = gradient.getDistributionGradient(dLL);
		//calls_optDiff++;
		double[] D = parameterGradient(dLL);
		return D;
	}
	
	public double[] optDistrDiff(double[] x)
	{
		setParameterValues(x);
		
		double[] dLL = gradient.getCorrectedGradient();
		dLL = gradient.getDistributionGradient(dLL);
		calls_optDiff++;
		
		double[] D = parameterGradient(dLL);
//				new double[distribution_params.size()];
//		for (int j=0; j<D.length; j++)
//			D[j] = -distribution_params.get(j).dL(dLL);

		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#*MLD.oDD "+calls_optDiff+"\t"+gradient.getCorrectedLL()+"\tL0 "+gradient.getUnobservedLL()+"\tfcalls "+calls_optFunc+"\tgradient "+FunctionMinimization.euclideanNorm(dLL)+"\tdL "+FunctionMinimization.euclideanNorm(D));
		
		return D;
		
	}

	
	
	private class GainParameter implements ModelParameter
	{
		GainParameter(int node)
		{
			this.node = node;
//			double κ = rates.getGainParameter(node);
//			node_parameters[node][PARAMETER_GAIN]=κ;
			assert (PURE_LOGISTIC);
			double log_gpar = getRates().getLogGainParameter(node, gain_parameter_type);
			
			node_parameters[node][PARAMETER_GAIN]=Math.exp(log_gpar);
		}
		private final int node;
		
		@Override
		public double get()
		{
//			double κ = rates.getGainParameter(node);
//			node_parameters[node][PARAMETER_GAIN]=κ;
			double gpar = node_parameters[node][PARAMETER_GAIN];
			return gpar;
		}
		
		@Override
		public void set(double g)
		{
			if (!Double.isFinite(g)) 
			{
				System.out.println("#**MLD.GP.set gain parameter="+g+"\tnpars "+Arrays.toString(node_parameters[node])+"\t"+rates.toString(node)+"\t// "+gradient.factory.toString(node));
			}
			assert Double.isFinite(g);
			assert (0.0<=g);
			node_parameters[node][PARAMETER_GAIN]=g;
			// rates.setParameters(node, κ, node_parameters[node][PARAMETER_LOSS], node_parameters[node][PARAMETER_DUPLICATION]);
		}
		
		@Override 
		public double dL(double[] gradient)
		{
			TreeWithLogisticParameters lrates = getRates();
			// we use the node_parameters here, but getRates() is equally good since we call after parameters are copied to the model
			double dL_dgpar;
			if (gain_parameter_type == PARAMETER_DUPLICATION)
			{
				dL_dgpar = gradient[3*node+PARAMETER_GAIN];
			} else if (gain_parameter_type == PARAMETER_GAIN)
			{
				if (lrates.getLogitDuplicationParameter(node) == Double.NEGATIVE_INFINITY)
				{ // Poisson
					dL_dgpar = gradient[3*node+PARAMETER_GAIN];
				} else
				{
					double dL_dkappa = gradient[3*node+PARAMETER_GAIN];
					// dL_dr = dL_kappa * kappa/r = dL_kappa / q
					double log_q = lrates.getLogDuplicationParameter(node);
					dL_dgpar = dL_dkappa / Math.exp(log_q);
				}
			} else
			{
				assert (gain_parameter_type == PARAMETER_LOSS);
				if (lrates.getLogitDuplicationParameter(node) == Double.NEGATIVE_INFINITY)
				{
					// Poisson
					double dL_dr = gradient[3*node+PARAMETER_GAIN];
					// dL_dgamma = dL_dr * r/gamma = dL_dr*p
					double log_p = lrates.getLogLossParameter(node);
					dL_dgpar = dL_dr * Math.exp(log_p);				
				} else
				{
					// Polya
					double dL_dkappa = gradient[3*node+PARAMETER_GAIN];
					// q*kappa = p*gamma; gamma = kappa * q/p; gamma/kappa = q/p
					// dL_dgamma = dL_dkappa * kappa/gamma = dL_dkappa * p/q
					if (is_duprate_bounded)
					{
						dL_dgpar = dL_dkappa / Math.exp(lrates.getLogitRelativeRate(node));
					} else
					{
						double log_lambda = lrates .getLogDuplicationParameter(node)-lrates.getLogLossParameter(node);
						dL_dgpar = dL_dkappa / Math.exp(log_lambda);
					}
				}
			}
			
			return dL_dgpar;
		}
		@Override 
		public String toString()
		{
			String name;
			if (gain_parameter_type == PARAMETER_LOSS)
				name = "gamma";
			else if (gain_parameter_type == PARAMETER_GAIN)
				name = "r'";
			else 
			{
				assert (gain_parameter_type == PARAMETER_DUPLICATION);
				if (node_parameters[node][PARAMETER_DUPLICATION]==Double.NEGATIVE_INFINITY)
				{
					name = "r";
				} else
				{
					name = "kappa";
				}
			}
			return name+node+"="+get();
		}
	}
	
	private class LossParameter implements BoundedParameter
	{
		LossParameter(int node)
		{
			this.node = node;
			double p = rates.getLossParameter(node);
			node_parameters[node][PARAMETER_LOSS]=p;
		}
		private final int node;
		
		@Override
		public double get()
		{
//			double p = rates.getLossParameter(node);
//			node_parameters[node][PARAMETER_LOSS]=p;
			double p = node_parameters[node][PARAMETER_LOSS];
			return p;
		}
		
		@Override
		public double getComplement()
		{
			double p_1 = node_parameters[node][PARAMETER_LOSS_COMPLEMENT];
			return p_1;
		}
		
		@Override
		public void set(double p)
		{
			if (track_complements)
				set(p, 1.0-p);
			else
				node_parameters[node][PARAMETER_LOSS]=p;
//			rates.setParameters(node, node_parameters[node][PARAMETER_GAIN], p, node_parameters[node][PARAMETER_DUPLICATION]);
		}
		
		@Override 
		public void set(double p, double p_1)
		{
			node_parameters[node][PARAMETER_LOSS]=p;
			node_parameters[node][PARAMETER_LOSS_COMPLEMENT]=p_1;
		}
		
		@Override 
		public double dL(double[] gradient)
		{
			double dLdp = gradient[3*node+PARAMETER_LOSS];
			if (is_duprate_bounded) // transformation from (p, q) to (p, λ)
			{
				double dLdq = gradient[3*node+PARAMETER_DUPLICATION];
				double dl_ratio = rates.getDuplicationRate(node)/rates.getLossRate(node);
				dLdp += dLdq * dl_ratio;
			}
			return dLdp;
		}
		@Override 
		public String toString()
		{
			return "p"+node+"="+get();
		}
		
		@Override
		public double getBound()
		{
			return 1.0;
		}
	}
	
	private class LossComplementParameter implements BoundedParameter
	{
		LossComplementParameter(int node)
		{
			this.node = node;
			double p1 = rates.getLossParameterComplement(node);
			
			node_parameters[node][PARAMETER_LOSS_COMPLEMENT]=p1;
			node_parameters[node][PARAMETER_LOSS]=rates.getLossParameter(node);
		}
		private final int node;		
		@Override
		public double get()
		{
			return node_parameters[node][PARAMETER_LOSS_COMPLEMENT];
		}
		@Override
		public void set(double p1)
		{
			if (track_complements)
				set(p1, 1.0-p1);
			else
			{
				assert (Double.isFinite(p1));
				assert (0.0<=p1);
				assert (p1<1.0);
				node_parameters[node][PARAMETER_LOSS_COMPLEMENT]=p1;
			}
		}
		@Override
		public void set(double p1, double p)
		{
			node_parameters[node][PARAMETER_LOSS_COMPLEMENT]=p1;
			node_parameters[node][PARAMETER_LOSS]=p;
		}
		@Override
		public double getComplement()
		{
			double p = node_parameters[node][PARAMETER_LOSS];
			return p;
		}
		@Override 
		public double dL(double[] gradient)
		{
			double dLdp = gradient[3*node+PARAMETER_LOSS];
			if (is_duprate_bounded) // transformation from (p, q) to (p, λ)
			{
				double dLdq = gradient[3*node+PARAMETER_DUPLICATION];
				double dl_ratio = rates.getDuplicationRate(node)/rates.getLossRate(node);
				dLdp += dLdq * dl_ratio;
			}
			return -dLdp;
		}
		@Override 
		public String toString()
		{
			return "1p"+node+"="+get();
		}
	}
	
//	private class DuplicationComplementParameter implements ModelParameter
//	{
//		DuplicationComplementParameter(int node)
//		{
//			this.node = node;
//		}
//		private final int node;
//		
//		@Override
//		public double get()
//		{
//			double q1 = rates.getDuplicationParameterComplement(node);
//			node_parameters[node][PARAMETER_DUPLICATION]=rates.getDuplicationParameter(node);
//			return q1;
//		}
//		
//		@Override
//		public void set(double q1)
//		{
//			node_parameters[node][PARAMETER_DUPLICATION]=1.0-q1;
//		}
//		@Override 
//		public double dL(double[] gradient)
//		{
//			return -gradient[3*node+PARAMETER_DUPLICATION];
//		}
//		@Override 
//		public String toString()
//		{
//			return "cq"+node+"="+get();
//		}	
//	}
	
	private class DuplicationParameter implements BoundedParameter
	{
		DuplicationParameter(int node)
		{
			this.node = node;
			double q = rates.getDuplicationParameter(node);
//			System.out.println("#**MLD.DP "+node+"\tq "+q+"\trate "+rates.getDuplicationRate(node)+"\tlen "+rates.getEdgeLength(node));
//			node_parameters[node][PARAMETER_DUPLICATION]=q;
			double q_1 = rates.getDuplicationParameterComplement(node);
			set(q,q_1);
		}
		
		private final int node;
		
		@Override
		public double get()
		{
//			double q = rates.getDuplicationParameter(node);
//			node_parameters[node][PARAMETER_DUPLICATION]=q;
			double q=node_parameters[node][PARAMETER_DUPLICATION];
			return q;
		}
		
		@Override
		public void set(double q)
		{
			if (q==1.0)
				throw new IllegalArgumentException();
			if (!Double.isFinite(q))
			{
				System.out.println("#**MLD.DP.set q "+q+"\t"+rates.toString(node)+"\tnpars "+Arrays.toString(node_parameters[node])+"//\t "+gradient.factory.toString(node));
			}
			assert (Double.isFinite(q));
			assert (0.0<=q);
			assert (q<1.0);
			node_parameters[node][PARAMETER_DUPLICATION]=q;
//			rates.setParameters(node, node_parameters[node][PARAMETER_GAIN], node_parameters[node][PARAMETER_LOSS], q);
		}
		
		@Override
		public double getComplement()
		{
			double q_1 = node_parameters[node][PARAMETER_DUPLICATION_COMPLEMENT];
			return q_1;
		}
		
		@Override
		public void set(double q, double q_1)
		{
			node_parameters[node][PARAMETER_DUPLICATION]=q;
			node_parameters[node][PARAMETER_DUPLICATION_COMPLEMENT]=q_1;
		}
		
		@Override 
		public double dL(double[] gradient)
		{
			return gradient[3*node+PARAMETER_DUPLICATION];
		}
		@Override 
		public String toString()
		{
			return "q"+node+"="+get();
		}
		
		@Override
		public double getBound() { return 1.0;}
	}
	
	private class DuplicationRate implements BoundedParameter
	{
		private final int node;
		DuplicationRate(int node)
		{ 
			this.node = node;
		}

		@Override
		public double getBound()
		{
			return rates.getLossRate(node);
		}
		
		@Override
		public void set(double λ)
		{		
			if (track_complements && !is_duprate_bounded)
				set(λ, 1.0-λ);
			else
				node_parameters[node][PARAMETER_DUPLICATION]=λ;
//			rates.setDuplicationRate(node, λ);
		}
		@Override 
		public void set(double x, double xcomp)
		{
			node_parameters[node][PARAMETER_DUPLICATION] = x;
			node_parameters[node][PARAMETER_DUPLICATION_COMPLEMENT] = xcomp;
		}
		
		@Override
		public double get()
		{
			double λ = node_parameters[node][PARAMETER_DUPLICATION];
			return λ;
		}
		
		@Override
		public double getComplement()
		{
			return node_parameters[node][PARAMETER_DUPLICATION_COMPLEMENT];
		}
		
		@Override
		public double dL(double[] distribution_gradient)
		{
			return gradient.inferDuplicationRateGradient(node, distribution_gradient);
		}
		
		@Override 
		public String toString()
		{
			return "d"+node+"="+get();
		}
	}
	
	
	private class LogisticLossParameter implements ModelParameter
	{
		LogisticLossParameter(int node)
		{
			assert (PURE_LOGISTIC);
			
			this.node = node;
			TreeWithLogisticParameters lrates = (TreeWithLogisticParameters )rates;
			node_parameters[node][PARAMETER_LOSS] = lrates.getLogitLossParameter(node);
			
			if (DEBUG_LOGISTIC)
			{
				this.debug_param = new Logistic(new SimpleParameter() {
					@Override
					public double dL(double[] gradient)
					{
						double dLdp = gradient[3*node+PARAMETER_LOSS];
						if (is_duprate_bounded) // transformation from (p, q) to (p, λ)
						{
							double dLdq = gradient[3*node+PARAMETER_DUPLICATION];
							double dl_ratio = rates.getDuplicationRate(node)/rates.getLossRate(node);
							dLdp += dLdq * dl_ratio;
						}
						return dLdp;
					}
				});
				debug_param.set(this.get());
			} else
				this.debug_param = null;
			
		}
		private final int node;
		
		private final Logistic debug_param;
		
		@Override
		public double get()
		{
//			double p = rates.getLossParameter(node);
//			node_parameters[node][PARAMETER_LOSS]=p;
			double logit_p = node_parameters[node][PARAMETER_LOSS];
			return logit_p;
		}
		
		
		@Override
		public void set(double logit_p)
		{
			node_parameters[node][PARAMETER_LOSS]=logit_p;
			if (DEBUG_LOGISTIC) debug_param.set(logit_p);
		}
		
		@Override 
		public double dL(double[] dLdparams)
		{
			double dLdp = dLdparams[3*node+PARAMETER_LOSS];
			double dLdp0 = dLdp;
			TreeWithLogisticParameters lrates = getRates();
			double dpdx = Math.exp(lrates.getLogLossParameter(node)+lrates.getLogLossComplement(node));
			double dLdx = dLdp * dpdx;
			if (is_duprate_bounded) // transformation from (p, q) to (p, λ)
			{
				double dLdq = dLdparams[3*node+PARAMETER_DUPLICATION];
				// dL by logit(q) = dLdq * q*(1-q)
				// (dLdq * q*(1-q)) * (1-p)/(1-q) = dLdq*q*(1-p)
				double dqdx = Math.exp(lrates.getLogDuplicationParameter(node)+lrates.getLogLossComplement(node));
				dLdx += dLdq * dqdx;
				
				if (lrates.getLogitDuplicationParameter(node)==Double.NEGATIVE_INFINITY)
				{
					if (gain_parameter_type == PARAMETER_LOSS)
					{
						double dL_dlogr = dLdparams[3*node+PARAMETER_GAIN]*lrates.getGainParameter(node);
						double dlr_dx = Math.exp(lrates.getLogLossComplement(node));
						dLdx += dL_dlogr*dlr_dx;
					}
				} else
				{
					if (gain_parameter_type == PARAMETER_GAIN)
					{
						double dL_dlogkappa = dLdparams[3*node+PARAMETER_GAIN]*lrates.getGainParameter(node);
						double dlkappa_dx = -Math.exp(lrates.getLogLossComplement(node));
						dLdx += dL_dlogkappa*dlkappa_dx;
					} 
				}
			} else // logit-p, logit-q
			{
				if (gain_parameter_type == PARAMETER_LOSS)
				{
					double dL_dlogg = dLdparams[3*node+PARAMETER_GAIN]*lrates.getGainParameter(node);
					double dlg_dx = Math.exp(lrates.getLogLossComplement(node));
					dLdx += dL_dlogg*dlg_dx;
				} 
			}
//			System.out.println("#**MLD.LLP.dL "+node+"\tdldx "+dLdx+"\tdldp "+dLdp0
//					+"\tdebug "+debug_param.dL(dLdparams)+"\t"+this.toString()+"\t"+debug_param.toString()+"\t"+lrates.toString(node)+"\t"+gradient.factory.toString(node));
			return dLdx;
		}
		@Override 
		public String toString()
		{
			return "logitp"+node+"["
					+Math.exp(Logarithms.logitToLogValue(get()))
					+"/1-"+Math.exp(Logarithms.logitToLogComplement(get()))
					+"; get "+get()
					//+"; debug "+debug_param.toString();
					+"]";			
		}

	}
	
	private class LogisticDuplicationParameter implements ModelParameter
	{
		LogisticDuplicationParameter(int node)
		{
			assert (PURE_LOGISTIC);
			this.node = node;
			TreeWithLogisticParameters lrates = (TreeWithLogisticParameters )rates;
			node_parameters[node][PARAMETER_DUPLICATION] = lrates.getLogitDuplicationParameter(node);
			
			if (DEBUG_LOGISTIC)
			{
				this.debug_param = new Logistic(new SimpleParameter()
				{
					@Override 
					public double dL(double[] gradient)
					{
						return gradient[3*node+PARAMETER_DUPLICATION];
					}
	
				});
				debug_param.set(this.get());
			} else
				this.debug_param = null;
		}
		
		private final int node;
		private final Logistic debug_param;
		
		@Override
		public double get()
		{
			double q=node_parameters[node][PARAMETER_DUPLICATION];
			return q;
		}
		
		@Override
		public void set(double y)
		{
			node_parameters[node][PARAMETER_DUPLICATION]=y;
			if (DEBUG_LOGISTIC) debug_param.set(y);
		}
		
		@Override 
		public double dL(double[] distribution_gradient)
		{
			double dLdq = distribution_gradient[3*node+PARAMETER_DUPLICATION];
			TreeWithLogisticParameters lrates = getRates();
			double dqdy = Math.exp(lrates.getLogDuplicationParameter(node)+lrates.getLogDuplicationComplement(node));
			double dLdy = dLdq*dqdy;
			
//			System.out.println("#**MLD.LDP.dL "+node+"\tdldy "+dLdy+"\tdldq "+dLdq
//					+"\tdebug "+debug_param.dL(distribution_gradient)+"\t"+this.toString()+"\t"+debug_param.toString()+"\t"+lrates.toString(node)
//						+"\t"+gradient.factory.toString(node));	
			
			if (lrates.getLogitDuplicationParameter(node)!=Double.NEGATIVE_INFINITY)
			{
				if (gain_parameter_type == PARAMETER_GAIN || gain_parameter_type == PARAMETER_LOSS)
				{
					double dL_dlogkappa = distribution_gradient[3*node+PARAMETER_GAIN]*lrates.getGainParameter(node);
					double dlkappa_dy = -Math.exp(lrates.getLogDuplicationComplement(node));
					dLdy += dL_dlogkappa*dlkappa_dy;
				}
			}
			
			return dLdy;
		}
		@Override 
		public String toString()
		{
			return "logitq"+node+"["
					+Math.exp(Logarithms.logitToLogValue(get()))
					+"/1-"+Math.exp(Logarithms.logitToLogComplement(get()))
					+"; get "+get()
					+"]"
				; //+"; debug "+debug_param.toString();
		}
		
	}
	
	private class LogisticDuplicationRate implements ModelParameter
	{
		private final int node;
		private final Logistic debug_param;
		LogisticDuplicationRate(int node)
		{ 
			assert (PURE_LOGISTIC);
			this.node = node;
			TreeWithLogisticParameters lrates = (TreeWithLogisticParameters )rates;
			node_parameters[node][PARAMETER_DUPLICATION]=lrates.getLogitRelativeRate(node);
			if (DEBUG_LOGISTIC)
			{
				this.debug_param = new Logistic(new SimpleParameter()
				{
					@Override
					public double dL(double[] distribution_gradient)
					{
						return gradient.inferDuplicationRateGradient(node, distribution_gradient);
					}
				});
				debug_param.set(this.get());
			} else
				debug_param = null;
		}

		@Override
		public void set(double λ)
		{		
			node_parameters[node][PARAMETER_DUPLICATION]=λ;
			if (DEBUG_LOGISTIC) debug_param.set(λ);
		}
		
		@Override
		public double get()
		{
			double λ = node_parameters[node][PARAMETER_DUPLICATION];
			return λ;
		}
		
		
		@Override
		public double dL(double[] distribution_gradient)
		{
			TreeWithLogisticParameters lrates = getRates();
//			double dLdr = gradient.inferDuplicationRateGradient(node, distribution_gradient);
			double dLdq = distribution_gradient[3*node+PARAMETER_DUPLICATION];
			double dqdz = Math.exp(lrates.getLogDuplicationParameter(node)+lrates.getLogRelativeComplement(node));
			double dLdz = dLdq * dqdz;
			
			
			if (lrates.getLogitDuplicationParameter(node)!=Double.NEGATIVE_INFINITY)
			{
				if (gain_parameter_type == PARAMETER_LOSS || gain_parameter_type == PARAMETER_GAIN)
				{
					double dL_dlogkappa = distribution_gradient[3*node+PARAMETER_GAIN]*lrates.getGainParameter(node);
					double dlkappa_dz = -Math.exp(lrates.getLogRelativeComplement(node));
					dLdz += dL_dlogkappa*dlkappa_dz;
				}
			}
//			System.out.println("#**MLD.LDR.dL "+node+"\t"+dLdz+"\tdebug "+debug_param.dL(distribution_gradient)+"\t"+this.toString()+"\t"+debug_param.toString()+"\t"+lrates.toString(node));			
			
			return dLdz;
		}
		
		@Override 
		public String toString()
		{
			return "logitdr"+node+"["
					+Math.exp(Logarithms.logitToLogValue(get()))
					+"/1-"+Math.exp(Logarithms.logitToLogComplement(get()))
					+"; get "+get()
				+"]";			
		}
	}
	
	
	
	@Override
	public double optimize(double delta)
	{
		return this.optimize(delta, 1024);
	}
	
	@Override
	public double optimize(double delta, final int itmax)
	{
		if (PRINT_OPTIMIZATION_MESSAGES)
		{
			double LL =  gradient.getCorrectedLL();
			System.out.println("#*MLD.o before LL "+LL+"\t; optimization with dup rates < 1.0 "+is_duprate_bounded);
		}
		initModelParameters(); // may change the parameters to respect rate or parameter bounds 
		
		List<Double> history = getOptimizationHistory();
		if (history == null)
		{
			history = new ArrayList<>();
			setOptimizationHistory(history);
		}
		
		if (AUTO_TRUNCATE)
		{
			double adjust_eps = Double.min(delta,1e-10);
			this.adjustCalculationWidth(adjust_eps);
		}
		
		double LL =  gradient.getCorrectedLL();
		
		double[] x0 = getParameterValues();

		if (PRINT_OPTIMIZATION_MESSAGES)
		{
			System.out.println("#*MLD.o starting LL "+LL+"\tparams "+distribution_params.size()+"\teval "+optDistrFunc(x0)+"\tdL "+FunctionMinimization.euclideanNorm(optDistrDiff(x0)));
			// printParameters(System.out);		
		}
		
		
		
		double[] prev_x = x0.clone();

		int h0 = history.size();
		int step_count = history.size()-h0;
		int nepoch = 0;
		
		double[] xbest = x0.clone();
		double LLbest = Double.NEGATIVE_INFINITY; // not = LL because it may have changed already
		
		while (step_count<itmax)
		{
			double min;
			
			if (do_em)
			{
//				SurvivalEM O = new SurvivalEM(gradient);
//				min = O.optimize(delta, x0.length);
//				initModelParameters(); // recalculate from rates set by EM
//				x0 = getParameterValues(); // starting point for further optimization
//				//gradient.computeParameters();
			}
			if (do_gradient_descent)
			{
				//min = FunctionMinimization.powell(x0, delta, itmax, x->optDistrFunc(x), history);
				min = FunctionMinimization.gradientDescent(x0, Math.sqrt(delta), x0.length, x->optDistrFunc(x), x->optDistrDiff(x), history);
				x0 = getParameterValues();

				LL =  gradient.getCorrectedLL(); // == -min 
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#*MLD.o GD "+LL+"\t("+min+")"+"\tfcalls "+calls_optFunc+"\tdcalls "+calls_optDiff);
			}
			
			
			// not too many iterations per epoch, 3 epochs preferred
			final int preferred_steps = Integer.max(200,1+x0.length);
			final int opti_steps;
			if (itmax-step_count<12)
				opti_steps = itmax-step_count;
			else if (nepoch==0)
					opti_steps =  Integer.min(preferred_steps,(itmax-step_count)/4);
			else if (nepoch==1)
				opti_steps = Integer.min(preferred_steps,(itmax-step_count)/3);
			else
				opti_steps = Integer.min(preferred_steps,(itmax-step_count));
						
			if (optimization_by_quasiNewton)
			{					
				 min  = FunctionMinimization.dfpmin(x0, delta, opti_steps, x->optDistrFunc(x), x->optDistrDiff(x), history);
			} else
			{
				min = FunctionMinimization.frprmn(x0, delta, opti_steps, x->optDistrFunc(x), x->optDistrDiff(x), history);
			}
			setParameterValues(x0); // x0 is updated by optimization
//			x0=getParameterValues();

			LL =  gradient.getCorrectedLL(); // == -min 
			
			if (LL>LLbest) // both negative
			{
				xbest = x0.clone();
				LLbest = LL;
			}

			double max_dx = 0.0;
			for (int j=0; j<x0.length;j++)
				if (x0[j]!=prev_x[j])
				{
					double dx;
					if (Double.isInfinite(x0[j]))
					{
						dx=Double.isInfinite(prev_x[j])?2.0:1.0;
					} else if (Double.isInfinite(prev_x[j]))
					{
						dx = 1.0;
					} else
					{
						dx = Math.abs(x0[j]-prev_x[j])/Double.max(1.0, Math.abs(prev_x[j]));
					}
					max_dx = Double.max(dx,max_dx);
				}
			
			
			double D[] = optDistrDiff(x0);
			double gradient_length = FunctionMinimization.euclideanNorm(D);
//			double L0 = gradient.factory.getEmptyLL();
//			double L1 = gradient.factory.getSingletonLL();
			double rel_gradient = gradient_length/(-LL);
			
			boolean done_dx = (max_dx < FunctionMinimization.DFP_TOLX);
			boolean done_gradient = rel_gradient<delta;
			
			boolean done_converged = ( done_gradient|| done_dx);

			step_count = history.size()-h0;
			
//			debugGradient(System.out);
			
			Count.out.println("#*MLD.o epoch "+nepoch+"\tLL "+LL+"\t("+min+")"
						+"\tfcalls "+calls_optFunc+"\tdcalls "+calls_optDiff
						//+"\tL0 "+L0+"\tL1 "+L1
						+"\tdL "+gradient_length
						+"\trgrad "+rel_gradient+"\tmax_dx "+max_dx+"\t"+(done_converged?"DONE":"loop")
						+"\tsteps "+step_count
						+"\tdelta "+delta);
			
			if (DEBUG_GRADIENT) debugGradient(System.out);
			
			prev_x = x0.clone();
			
			// check if stuck in some extreme of too short or too long edges			
			if (step_count<itmax && 0<REINIT_EXTREME_EDGES
					&& nepoch % 2 == 1) // possibly more iterations
			{
				IndexedTree tree = rates.getTree();
				boolean updated = false;
				boolean[] updated_nodes=new boolean[tree.getNumNodes()];

				double small = Math.sqrt(1.0/gradient.factory.table.getFamilyCount());
				// 1. check root 
				int root = tree.getRoot();
				if (optimize_parameters[3*root+PARAMETER_DUPLICATION] || optimize_parameters[3*root+PARAMETER_GAIN])
				{
					double rm = rates.getRootMean();
					if (rm < small*small)
					{
						// reset to a random value 
						String oldnodestr = rates.toString(root);
						rates.initNodeParameters(root);
						if (PRINT_OPTIMIZATION_MESSAGES)
							System.out.println("#*MLD.o root mean too small "+rm+"\told "+oldnodestr
										+"\treset "+rates.toString(root));
						updated_nodes[root] = updated = true;
					} else
					{
						if (PRINT_OPTIMIZATION_MESSAGES)
							System.out.println("#*MLD.o root mean ok "+rm);
					}
				}
				// 2. check edges 
				for (int node = 0; node<root; node++)
				{
					if (optimize_parameters[3*node+PARAMETER_LOSS])
					{
						double p1  = rates.getLossParameterComplement(node);
						double p = rates.getLossParameter(node);
						if (p1<small)
						{
							String oldnodestr = rates.toString(node);
							rates.initNodeParameters(node);
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLD.o edge too long\told "+oldnodestr
											+"\treset "+rates.toString(node));
	
							updated_nodes[node] =  updated = true;
						} else if (p<small*small)
						{
							//String oldnodestr = rates.toString(node);
							
							rates.initNodeParameters(node);
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLD.o edge too short // oldp "+p //\told "+oldnodestr
											+"\treset "+rates.toString(node));
							updated_nodes[node] =  updated = true;
						}
					}
					if (optimize_parameters[3*node+PARAMETER_DUPLICATION]&& !updated_nodes[node])
					{
						double drate = rates.getDuplicationRate(node);
						if (drate>1.0) drate = 1.0/drate;
						boolean too_small = (drate<small*small) && false;
						boolean too_big = is_duprate_bounded && (1.0-small<drate)
								&& false;
						
						if (too_small || too_big)
						{
							String oldnodestr = rates.toString(node);
							if (optimize_parameters[3*node+PARAMETER_LOSS])
							{
								rates.initNodeParameters(node);
								
							} else
							{
								double p = rates.getLossParameter(node);
								double p1 = rates.getLossParameterComplement(node);
								rates.initNodeParameters(node);
								rates.setEdgeLengthForLossParameter(node, p, p1);
							}
							if (PRINT_OPTIMIZATION_MESSAGES)
							{
								if (too_small)
								{
									System.out.println("#*MLD.o edge duprate too small\told "+oldnodestr
												+"\treset "+rates.toString(node));
								} else 
								{
									assert (too_big);
									System.out.println("#*MLD.o edge duprate too big\told "+oldnodestr
											+"\treset "+rates.toString(node));
								}
							}
							updated_nodes[node] =  updated = true;
						}
					}
				}
				
				
				if (updated)
				{
					if (1<REINIT_EXTREME_EDGES)
					{
						boolean[] neighbor_nodes = new boolean[updated_nodes.length];
						for (int node = 0; node<updated_nodes.length; node++)
						{
							if (updated_nodes[node] && node != root)
							{
								int parent = tree.getParent(node);
								neighbor_nodes[parent] = true;
								if (tree.getNumChildren(parent)==2)
								{
									neighbor_nodes[tree.getChild(parent, 0)]=
									neighbor_nodes[tree.getChild(parent, 1)]=true;
									
								}
								if (3<REINIT_EXTREME_EDGES && tree.getNumChildren(node)==2)
								{
									neighbor_nodes[tree.getChild(node, 0)]=
									neighbor_nodes[tree.getChild( node, 1)]=true;
								}
							}
						}
						for (int node = 0; node<updated_nodes.length; node++)
							if (optimize_parameters[3*node+PARAMETER_DUPLICATION]
									|| optimize_parameters[3*node+PARAMETER_LOSS])
						{
							
							if (neighbor_nodes[node] && !updated_nodes[node])
							{
								
								String oldnodestr = rates.toString(node);
								if (optimize_parameters[3*node+PARAMETER_LOSS])
								{
									rates.initNodeParameters(node);
									
								} else
								{
									double p = rates.getLossParameter(node);
									double p1 = rates.getLossParameterComplement(node);
									rates.initNodeParameters(node);
									rates.setEdgeLengthForLossParameter(node, p, p1);
								}
								if (PRINT_OPTIMIZATION_MESSAGES)
									System.out.println("#*MLD.o neighbor\told "+oldnodestr
												+"\treset "+rates.toString(node));
							}
						}
					}
					
					initModelParameters(); // recalculates parameter values
					done_converged = false;
					x0 = getParameterValues();
				} 
				
				
				
			}
			if (done_converged) break;
			// back to more iterations 
			if (AUTO_TRUNCATE && step_count<itmax)
			{
				double adjust_eps = Double.min(delta,1e-10);
				this.adjustCalculationWidth(adjust_eps);
			}
						
			nepoch++;
		}
		
		if (LL<LLbest)
		{
			if (PRINT_OPTIMIZATION_MESSAGES)
				System.out.println("#*MLD.o returning to previous optimum\t"+LLbest+"\tinstead of "+LL);
			LL = LLbest;
			setParameterValues(xbest);
			x0 = getParameterValues();
			
		}
		
				
		return -LL;
	}
	
	/**
	 * Tries to increase/decrease absolute an relative calculation width parameters 
	 * for better estimation of true likelihood and gradient (with max 
	 * calculation width).
	 * 
	 * @param tol
	 * @return log-likelihood
	 */
	private double adjustCalculationWidth(double tol)
	{
		int absolute = gradient.getCalculationWidthAbsolute();
		double relative = gradient.getCalculationWidthRelative();
		
		double current_LL = gradient.getCorrectedLL();
		double current_dL = FunctionMinimization.euclideanNorm(parameterGradient());
		
		// need to increase?
		double step_size = Math.log(2.0)/3.0;
		
		
		boolean adjustCalculationWidth=false;
		
		int num_adjustments = 0;
		
		int dabs =0, drel = 0;
		
		final int maxiter = 8;
		boolean adjusted_in_iteration=true;
		while (adjusted_in_iteration && num_adjustments < maxiter)
		{
			adjusted_in_iteration = false;

			{ // try changing absolute
				
				
				int next_absolute = Integer.max(absolute+1,(int)Math.ceil(Math.exp(Math.log(absolute)+step_size)));
				this.setCalculationWidth(next_absolute, relative);
				
				double next_dL = FunctionMinimization.euclideanNorm(parameterGradient());
				double next_LL = gradient.getCorrectedLL();
				
				double next_delta = next_LL-current_LL;
				double next_rdiff = Math.abs(next_delta/current_LL); 
				
				double next_ddiff = Math.abs(next_dL-current_dL)/Math.max(1.0,current_dL);
	
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#**MLD.aCW ("+dabs+","+drel+") absolute "+absolute+"\tincrease "+next_absolute
							+"\trdiff "+next_rdiff+"\tddiff "+next_ddiff
							+"\twas "+current_LL+"\tnext "+next_LL
							+"\tdL "+current_dL+"\tnextdL "+next_dL+"\t(tol "+tol+")");
	
				if (tol < next_rdiff || tol < next_ddiff)
				{
					absolute = next_absolute;
					current_LL = next_LL;
					current_dL = next_dL;
					adjusted_in_iteration = adjustCalculationWidth = true;
					++num_adjustments;
					++dabs;
	
				} else if (DEFAULT_TRUNCATE_ABSOLUTE < absolute)
				{
					int prev_absolute = Integer.max(Integer.min((int)Math.ceil(Math.exp(Math.log(absolute)-step_size)), absolute-1), DEFAULT_TRUNCATE_ABSOLUTE);
					this.setCalculationWidth(prev_absolute, relative);
					
					double prev_LL = gradient.getCorrectedLL();
					double prev_dL = FunctionMinimization.euclideanNorm(parameterGradient());
					double prev_delta = prev_LL-current_LL;
					double prev_rdiff = Math.abs(prev_delta/current_LL);
					double prev_ddiff = Math.abs(prev_dL-current_dL)/Math.max(1.0,current_dL);
					
					if (PRINT_OPTIMIZATION_MESSAGES)
						System.out.println("#**MLD.aCW ("+dabs+","+drel+") absolute "+absolute+"\tdecrease "+prev_absolute
								+"\trdiff "+prev_rdiff+"\tddiff "+prev_ddiff
								+"\twas "+current_LL+"\tprev "+prev_LL
								+"\tdL "+current_dL+"\tprevdL "+prev_dL+"\t(tol "+tol+")");
					if (prev_rdiff < tol && prev_ddiff < tol)
					{
						absolute = prev_absolute;
						current_LL = prev_LL;
						current_dL = prev_dL;
						adjusted_in_iteration = adjustCalculationWidth = true;
						++num_adjustments;
						--dabs;
					}
				}
			}			
			{
				// try changing relative
				double next_rel = Math.exp(Math.log(relative)+step_size);
				this.setCalculationWidth(absolute, next_rel);
				
				double next_dL =  FunctionMinimization.euclideanNorm(parameterGradient());
				double next_LL = gradient.getCorrectedLL();
				
				double next_delta = next_LL-current_LL;
				double next_rdiff = Math.abs(next_delta/current_LL); 
				
				double next_ddiff = Math.abs(next_dL-current_dL)/Math.max(1.0,current_dL);
				
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#**MLD.aCW ("+dabs+","+drel+") relative "+relative+"\tincrease "+next_rel
							+"\trdiff "+next_rdiff+"\tddiff "+next_ddiff
							+"\twas "+current_LL+"\tnext "+next_LL
							+"\tdL "+current_dL+"\tnextdL "+next_dL+"\t(tol "+tol+")");
	
				if (tol < next_rdiff || tol < next_ddiff)
				{
					relative = next_rel;
					current_LL = next_LL;
					current_dL = next_dL;
					adjusted_in_iteration = adjustCalculationWidth = true;
					++num_adjustments;
					++drel;
				} else if (DEFAULT_TRUNCATE_RELATIVE < relative)
				{
					double prev_rel = Double.max(DEFAULT_TRUNCATE_RELATIVE, Math.exp(Math.log(relative)-step_size));
					this.setCalculationWidth(absolute, prev_rel);
					
					double prev_LL = gradient.getCorrectedLL();
					double prev_dL =  FunctionMinimization.euclideanNorm(parameterGradient());
					double prev_delta = prev_LL-current_LL;
					double prev_rdiff = Math.abs(prev_delta/current_LL);
					double prev_ddiff = Math.abs(prev_dL-current_dL)/Math.max(1.0,current_dL);
					
					if (PRINT_OPTIMIZATION_MESSAGES)
						System.out.println("#**MLD.aCW ("+dabs+","+drel+") relative "+relative+"\tdecrease "+prev_rel
								+"\trdiff "+prev_rdiff+"\tddiff "+prev_ddiff
								+"\twas "+current_LL+"\tprev "+prev_LL
								+"\tdL "+current_dL+"\tprevdL "+prev_dL+"\t(tol "+tol+")");								
					if (prev_rdiff < tol && prev_ddiff < tol)
					{
						relative = prev_rel;
						current_LL = prev_LL;
						current_dL = prev_dL;
						adjusted_in_iteration = adjustCalculationWidth = true;
						++num_adjustments;
						--drel;
					}
				}
			} 			
		} // while changes
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#*MLD.aCW setting "+absolute+","+relative);		
		this.setCalculationWidth(absolute, relative);
		
		return current_LL;
		
	}
	
	
	
	/**
	 * Which of the 3/node parameters should be optimized: gain rate if non-zero,
	 * loss rate if not 1.0, duplication if not zero.  
	 * 
	 * @param rates
	 * @return
	 */
	public static boolean[] optimizableParams(TreeWithRates rates)
	{
		IndexedTree tree = rates.getTree();
		int num_nodes = tree.getNumNodes();
		
		boolean[] opt_par = new boolean[3*num_nodes];
		for (int node=0; node<num_nodes; node++)
		{
//			if (tree.isRoot(node))
//			{
				opt_par[3*node+PARAMETER_GAIN] = rates.getLogGainParameter(node)!=Double.NEGATIVE_INFINITY;
				opt_par[3*node+PARAMETER_LOSS] = (rates.getLogLossComplement(node)!=Double.NEGATIVE_INFINITY);
				opt_par[3*node+PARAMETER_DUPLICATION] = rates.getLogDuplicationParameter(node)!=Double.NEGATIVE_INFINITY;
//				// DEBUG
//				System.out.println("#**MLD.oP root "+node
//						+"\tog "+opt_par[3*node+PARAMETER_GAIN]
//						+"\tod "+opt_par[3*node+PARAMETER_DUPLICATION]
//						+"\tol "+opt_par[3*node+PARAMETER_LOSS]
//						+"\t// "+rates.toString(node)
//						);
//			} else
//			{
//				opt_par[3*node + PARAMETER_GAIN] = rates.getGainParameter(node)>0.0;
//				opt_par[3*node+PARAMETER_LOSS] = (rates.getLossParameter(node)<1.0);
//				opt_par[3*node+PARAMETER_DUPLICATION] = rates.getDuplicationParameter(node)>0.0;
//
//			} 
		}
		
		return opt_par;
	}

//	private static final double RATES_DIFF_EPS = 1e-7; // cubic root of machine precision
//	
//	protected double[] estimateDistrDiff(double[] x)
//	{
//		double fx = optDistrFunc(x);
//		
//		double[] D = new double[x.length];
//		for (int p=0; p<x.length; p++)
//		{
//			double θ = x[p];
//			double h = Math.abs(θ*RATES_DIFF_EPS);
//			x[p] = θ+h;
//			double fd = optDistrFunc(x);
//			double delta = fd-fx;
//			double dfdθ = delta/h;
//			D[p] = dfdθ;
//			
//			x[p] = θ;
//		}
//		setParameterValues(x);
//		return D;
//		
//	}
	
	/*
	 * debug methods
	 */
	
	private void printParameters(PrintStream out)
	{
		for (int j=0; j<distribution_params.size(); j++)
		{
			out.println("#param\t"+distribution_params.get(j));
		}
	}

	private void debugGradient(PrintStream out)
	{
		//optimize(0,0); 
		double[] x = getParameterValues();
		
		double[] df = optDistrDiff(x);
		
		double[] df_est = FunctionMinimization.numericalGradient(θ->optDistrFunc(θ), x);// estimateDistrDiff(x); // FunctionMinimization.numericalGradient(θ->optDistrFunc(θ), x);
		
		double[] xafter = getParameterValues();

		boolean had_details = false;
		
		double sum_delta = 0.0;
		double sumsq_delta = 0.0;
		
		
		for (int i=0; i<x.length; i++)
		{
			ModelParameter P = distribution_params.get(i);
			
			double df_delta = df_est[i]-df[i];
			double rel_delta = Math.abs(df_delta/df[i]);
			
			sum_delta += df_delta;
			sumsq_delta += df_delta*df_delta;
			
			out.println("#**MLD.dG param "+i+"\tdfdx "+df[i]+"\test "+df_est[i]
					+"\tdiff "+df_delta+"\t("+rel_delta+")"
					+"\t"+P
					//+"\tx0 "+x[i]+"\tx1 "+xafter[i] // checking if parameters stay unchanged
				);
			
			boolean want_details = P instanceof LogisticLossParameter
					//&& df_est[i]*df[i]<0.0
					&& 1.0 < rel_delta 
					&& 1.0<Math.abs(df_est[i])+Math.abs(df[i]);
			
			if (want_details)
			{
				this.debugLossGradient(out, (LogisticLossParameter)P, 12);
				had_details = true;
			}
			
		}		
		double df_len = FunctionMinimization.euclideanNorm(df);
		double avg_delta = sum_delta/df_len;
		double sd_delta =  Math.sqrt(sumsq_delta/df_len);
				
		out.println("#**MLD.dG dflen "+df_len+"\tsumdelta "+sum_delta
				+"\tavgdelta "+avg_delta+"\tsd_delta "+sd_delta
				);
		
		
		
		if (had_details)
		{
			GammaInvariant model = new GammaInvariant(rates,1,1,1,1);
			out.println(count.io.RateVariationParser.printRates(model));
			System.exit(2024);
		}
		
	}
	
	private void debugLossGradient(PrintStream out,  LogisticLossParameter P ,  int steps)
	{
		// TODO
		double eps = FunctionMinimization.GRADIENT_DIFF_EPS;
		
		double x0 = P.get();
		double f0 = -gradient.getCorrectedLL();
		
		ModelParameter Pdup;
		if (is_duprate_bounded)
		{
			Pdup = new LogisticDuplicationRate(P.node);
		} else
		{
			Pdup = new LogisticDuplicationParameter(P.node);
		}
		ModelParameter Pgain = newLogistic(new GainParameter(P.node), MAX_GAIN_RATE);
		
		System.out.println("#*MLD.dLossG\tnode\tstep\tx\tfx"
				+"\tdx\tdf\tdf/dx"
				+"\tdfdlp\tdfdllm\tdfdlgn"
				+ "\tdsfdp\tdsfdq\tdfdp\tdfdq");
		
		double fminus=Double.NaN, fplus=Double.NaN;
		
		for (int d=-steps; d<=steps; d++)
		{
			double xdelta = Math.abs(x0)*(eps*d)/steps;
			double x = x0+xdelta;
			P.set(x);
			copyNodeParametersToModel();
			gradient.computeParameters();	
			
			double LL = gradient.getCorrectedLL();
			
			double[] dLL = gradient.getCorrectedGradient();
			
			double survival_dloss = -dLL[3*P.node+PARAMETER_LOSS]
					*Math.exp(gradient.factory.getLogLossParameter(P.node)+gradient.factory.getLogLossComplement(P.node));
			double survival_ddup = -dLL[3*P.node+PARAMETER_DUPLICATION]
					*Math.exp(gradient.factory.getLogLossParameter(P.node)+gradient.factory.getLogLossComplement(P.node));
			dLL = gradient.getDistributionGradient(dLL);
			
			double dfdp = -dLL[3*P.node+PARAMETER_LOSS]*Math.exp(rates.getLogLossParameter(P.node)+rates.getLogLossComplement(P.node));
			double dfdq	= -dLL[3*P.node+PARAMETER_DUPLICATION]*Math.exp(rates.getLogDuplicationParameter(P.node)+rates.getLogDuplicationComplement(P.node));
			
			double fx = -LL;
			double dloss = -P.dL(dLL);
			double ddup  = -Pdup.dL(dLL);
			double dgain  = -Pgain.dL(dLL);
			
			double fdelta = fx-f0;
			double dest = xdelta==0.0?0.0:fdelta/xdelta;
			
			if (d!=0)
			System.out.println("#*MLD.dLossG\t"+P.node
					+"\t"+d
					+"\t"+x
					+"\t"+fx
					+"\t"+xdelta
					+"\t"+fdelta
					+"\t"+dest
					+"\t"+dloss
					+"\t"+ddup
					+"\t"+dgain
					+"\t"+survival_dloss+"\t"+survival_ddup
					+"\t"+dfdp+"\t"+dfdq
					+"\t"+Math.exp(rates.getLogLossComplement(P.node)-rates.getLogDuplicationComplement(P.node))
					+"\t// "+P
					);
			
			if (d==-steps) fminus = fx;
			if (d==steps) fplus = fx;
		}
		double dest = (fplus-fminus)/(2.0*Math.abs(x0)*eps);
		System.out.println("#*MLD.dLossG\t"+P.node+"\t0\tx0\t"+f0+"\t"+dest);
		P.set(x0);
		copyNodeParametersToModel();
		gradient.computeParameters();	
	}
	
	
	
	
	private void gainLikelihoodInterval(PrintStream out, int node, double significance_alpha, double delta, int itmax)
	{
		double rtol = 1.0/(1<<20); 
		// want to find max gain and min gain s.t 
		// L(g') >= L(g*)* p 
		// log L(g') >= log L(g*)+log p
		// -log L(g') <= -log L(g*)-log p
		
		// -log L(g') - (-log L(g*)) <= -log p
		// -log L(g')- (-log L(g*)) + log p <= 0
		// starts with <0 as moving away from g*
		// 
		final double negLL0 = optimize(delta, 0);
		double[] model0 = getParameterValues();
		
		double r0 = rates.getGainParameter(node);
//		final double log_pregion = pregion<0.0?Math.log(-pregion):Math.log(pregion);
		
		final Map<Double,Double> diff_cache = new HashMap<>();
		diff_cache.put(r0, 0.0);
		final Map<Double, double[]> param_cache = new HashMap<>();
		param_cache.put(r0, model0);
		
		DoubleFunction<Double> diffLL  = new DoubleFunction<>()
		{
			@Override
			public Double apply(double r)
			{
				double diffLL;
				if (diff_cache.containsKey(r)) 
				{
					diffLL= diff_cache.get(r);
				} else
				{
					rates.initNodeParameters(node); // so that we start with random duplication and gain parameters 			
					double p = rates.getLossParameter(node);
					rates.setGainRateForLossParameter(node, r, p);
					fixGain(node, true);
					// gradient.computeParameters(); // not needed bc optimize() will do that
					double negLL = optimize(delta, itmax);
					diffLL = negLL-negLL0;
					
					double[] model = getParameterValues();
	
					fixGain(node, false);
					initModelParameters();
					setParameterValues(model0); // and reset
					
//					double dval = diff+log_pregion;
//					diff_cache.put(r, dval);
					diff_cache.put(r, diffLL);
					param_cache.put(r, model);
				
//					out.println("#**MLD.gLI\t"+r+"\tdiff "+diffLL+"\tr0 "+r0+"\tLL "+negLL+"\tLL0 "+negLL0);
				}
				// LRT test
				double chi_square_p = Functions.Chi_square_tail(1, 2.0*Math.max(0.0, diffLL)); 
				out.println("#**MLD.gLI\t"+r+"\tdiff "+diffLL+"\tr0 "+r0+"\tLL0 "+negLL0+"\tp "+chi_square_p);
				return Math.abs(significance_alpha) - chi_square_p;
			}
		};
		
		if (0.0<significance_alpha)
		{
			// bracket for rmax
			
			double r2 = r0;
			double d2 = significance_alpha-1.0; // negative 
			double r1,d1;
			
			assert (d2<0.0);
			do
			{
				r1 = r2;
				d1 = d2;
				r2 = 1.25*r1;
				d2 = diffLL.apply(r2);
			} while (d2<0.0 && r1*1.25<MAX_GAIN_RATE);
			
			double rmax;
			if (d2==0.0) // unlikely
			{
				// got rmax 
				rmax = r2;
			} else if (d2<0) // at max 
			{
				rmax = r2;
			} else 
			{
				// d2>0
				rmax = FunctionMinimization.zbrent(diffLL, r1, r2, rtol);
			}
			double dmax = diffLL.apply(rmax); // caches model for rmax
			out.println("#**MLD.gLI\trmax "+rmax);
			
			double[] model = param_cache.get(rmax);
			fixGain(node,true);
			double p = rates.getLossParameter(node);
			rates.setGainRateForLossParameter(node, rmax, p);
			initModelParameters();
			setParameterValues(model);
			fixGain(node,false);
			initModelParameters();
		} else
		{
			// bracket for rmin
			double r1 = r0;
			double d1 = -significance_alpha-1.0;
	
			double r2,d2;
			do
			{
				r2 = r1;
				d2 = d1;
				r1 = 0.8*r2;
				d1 = diffLL.apply(r1);		
			} while (d1<0.0 && r1>1e-6); // very small
			
			double rmin;
			if (d1<=0.0) rmin = r1;
			else rmin = FunctionMinimization.zbrent(diffLL, r1, r2, rtol);

			out.println("#**MLD.gLI\trmin "+rmin);
			double dmin = diffLL.apply(rmin);
			double[] model = param_cache.get(rmin);
			fixGain(node,true);
			double p = rates.getLossParameter(node);
			rates.setGainRateForLossParameter(node, rmin, p);
			initModelParameters();
			setParameterValues(model);
			fixGain(node,false);
			initModelParameters();
		}
	}
		
	
	private void testEdgeGain(PrintStream out, int node, double delta, int itmax)
	{
		double r = rates.getGainParameter(node);
		List<Double> rtest = new ArrayList<>();
		
		double sqr = Math.sqrt(r);
		double rmin = 1e-4;
		while (sqr<rmin)
		{
			rtest.add(sqr);
			sqr = Math.sqrt(sqr);
		}
		rtest.add(rmin);
		rtest.add(0.001);
		rtest.add(0.0025);
		rtest.add(0.005);
		rtest.add(0.01);
		rtest.add(0.02);
		rtest.add(0.04);
		rtest.add(0.05);
		rtest.add(0.06);
		rtest.add(0.08);
		rtest.add(0.1);
		rtest.add(0.15);
		rtest.add(0.2);
		rtest.add(0.25);
		rtest.add(0.32);
		rtest.add(0.5);
//		rtest.add(0.8);
		rtest.add(1.0);
		rtest.add(2.0);
		rtest.add(3.0);
		
		java.util.Collections.sort(rtest);
		final int ntest = rtest.size();

		out.println("#**MLD.tEG node\t"+node+"\tntest "+ntest+"\t"+rtest.toString()+"\t// "+ rates.toString(node));
		out.println("#NODEGAIN\tnode\tr\tLL\tgetr\t|optdL|\tdLdr\t|dL|\tnodeinfo");
		
		// starting values
		double negLL = optimize(delta, 0);
		double[] model0 = getParameterValues();
		double[] grad = gradient.getDistributionGradient(gradient.getCorrectedGradient());
		double glen = FunctionMinimization.euclideanNorm(grad); 
		double dLdr = grad[3*node + PARAMETER_GAIN];
		double[] pgrad = parameterGradient(grad);
		double dlen = FunctionMinimization.euclideanNorm(pgrad);
		
		out.println("#NODEGAIN0\t"+node+"\t"+r+"\t"+(-negLL)+"\t"+r+"\t"+dlen+"\t"+dLdr+"\t"+glen+"\t"+rates.toString(node)+"\t"+gradient.factory.toString(node));
		
		for (int round=0; round<ntest; round++)
		{
			double thisr = rtest.get(round);
			rates.initNodeParameters(node); // so that we start with random duplication and gain parameters 			

			double p = rates.getLossParameter(node);
			rates.setGainRateForLossParameter(node, thisr, p);
			
			out.println("#**MLD.tEL node\t"+node+"\tround "+round+"\tr "+thisr+"\t"+rates.toString(node));
			
			fixGain(node, true);
			// gradient.computeParameters(); // not needed bc optimize() will do that
			negLL = optimize(delta, itmax);

			r = rates.getGainParameter(node);
			
			grad = gradient.getDistributionGradient(gradient.getCorrectedGradient());
			dLdr = grad[3*node + PARAMETER_GAIN];
			pgrad = parameterGradient(grad);

			dlen = FunctionMinimization.euclideanNorm(pgrad);
			glen = FunctionMinimization.euclideanNorm(grad);
			
			r=rates.getGainParameter(node);
			out.println("#NODEGAIN\t"+node+"\t"+thisr+"\t"+(-negLL)+"\t"+r+"\t"+dlen+"\t"+dLdr+"\t"+glen+"\t"+rates.toString(node)+"\t"+gradient.factory.toString(node));

			fixGain(node, false);
			initModelParameters();
			setParameterValues(model0); // and reset
		}
	}
	
	
	private void testEdgeLoss(PrintStream out, int node, double delta, int itmax)
	{
		double p = rates.getLossParameter(node);
		double q = rates.getDuplicationParameter(node);
		double q1 = rates.getDuplicationParameterComplement(node);

		List<Double> ptest = new ArrayList<>();
		
//		ptest.add(p);
		double sqp = Math.sqrt(p);
		while (sqp<1e-6)
		{
			ptest.add(sqp);
			sqp = Math.sqrt(sqp);
		}
		ptest.add(1e-6);
		ptest.add(1e-5);
		ptest.add(2e-5);
		ptest.add(5e-5);
		ptest.add(0.0001);
		ptest.add(0.0005);
		ptest.add(0.001);
		ptest.add(0.0025);
		ptest.add(0.004);
		ptest.add(0.005);
		ptest.add(0.008);
		ptest.add(0.01);
		ptest.add(0.025);
		ptest.add(0.05);
		ptest.add(0.08);
		ptest.add(0.1);
		ptest.add(0.12);
		ptest.add(0.15);
		ptest.add(0.18);
		ptest.add(0.2);
		ptest.add(0.25);
		ptest.add(0.32);
		ptest.add(0.4);
		ptest.add(0.5);
		ptest.add(0.8);
		ptest.add(0.9);
		ptest.add(0.96);
		ptest.add(0.98);
//		ptest.add(0.99);
//		ptest.add(0.995);
//		ptest.add(0.999);
//		ptest.add(0.9999);
//		ptest.add(0.99999);
		
		java.util.Collections.sort(ptest);
		final int ntest = ptest.size();

		out.println("#**MLD.tEL node\t"+node+"\tntest "+ntest+"\t"+ptest.toString()+"\t// "+ rates.toString(node));
		out.println("#NODELOSS\tnode\tp\tLL\tlen\tp(len)\t|optdL|\tdLdp\t|dL|\tnodeinfo");
		
		// starting values
		double negLL = optimize(delta, 0);
		double[] model0 = getParameterValues();
		double[] grad = gradient.getDistributionGradient(gradient.getCorrectedGradient());
		double glen = FunctionMinimization.euclideanNorm(grad); 
		double dLdp = grad[3*node + PARAMETER_LOSS];
		double[] pgrad = parameterGradient(grad);
		double dlen = FunctionMinimization.euclideanNorm(pgrad);
//		double dLdθ = pgrad[3*node+PARAMETER_LOSS];
		
		double t = rates.getEdgeLength(node);
		out.println("#NODELOSS0\t"+node+"\t"+p+"\t"+(-negLL)+"\t"+t+"\t"+p+"\t"+dlen+"\t"+dLdp+"\t"+glen+"\t"+rates.toString(node)+"\t"+gradient.factory.toString(node));

		int parent = rates.getTree().getParent(node);
		int sib = (rates.getTree().isRoot(node) 
				|| rates.getTree().getNumChildren(parent)!= 2)?-1:rates.getTree().getSibling(node);
		
		if (!rates.getTree().isRoot(node) )
		{
			out.println("#NODEPARENT\t"+parent+"\t"+p+"\t\t"+rates.getEdgeLength(parent)+"\t"+rates.getLossParameter(parent)+"\t\t\t\t"+rates.toString(parent)+"\t"+gradient.factory.toString(parent));			
		}
		if (sib != -1)
		{
			out.println("#NODESIB\t"+sib+"\t"+p+"\t\t"+rates.getEdgeLength(sib)+"\t"+rates.getLossParameter(sib)+"\t\t"+grad[3*sib + PARAMETER_LOSS]+"\t\t"+rates.toString(sib)+"\t"+gradient.factory.toString(sib));			
		}
		
		for (int round=0; round<ntest; round++)
		{
			double thisp = ptest.get(round);
			rates.initNodeParameters(node); // so that we start with random duplication and gain parameters 			
			rates.setEdgeLengthForLossParameter(node, thisp, 1.0-thisp);
			t = rates.getEdgeLength(node);
			
			rates.setEdgeLength(node, t);
			p = rates.getLossParameter(node);
			out.println("#**MLD.tEL node\t"+node+"\tround "+round+"\tp "+thisp+"\t"+rates.toString(node));
			
			fixLoss(node, true);
			// gradient.computeParameters(); // not needed bc optimize() will do that
			negLL = optimize(delta, itmax);

			grad = gradient.getDistributionGradient(gradient.getCorrectedGradient());
			dLdp = grad[3*node + PARAMETER_LOSS];
			pgrad = parameterGradient(grad);

			dlen = FunctionMinimization.euclideanNorm(pgrad);
			glen = FunctionMinimization.euclideanNorm(grad);
			
			p=rates.getLossParameter(node);
			out.println("#NODELOSS\t"+node+"\t"+thisp+"\t"+(-negLL)+"\t"+t+"\t"+p+"\t"+dlen+"\t"+dLdp+"\t"+glen+"\t"+rates.toString(node)+"\t"+gradient.factory.toString(node));
			if (!rates.getTree().isRoot(node) )
			{
				out.println("#NODEPARENT\t"+parent+"\t"+thisp+"\t\t"+rates.getEdgeLength(parent)+"\t"+rates.getLossParameter(parent)+"\t\t\t\t"+rates.toString(parent)+"\t"+gradient.factory.toString(parent));			
			}
			if (sib != -1)
			{
				out.println("#NODESIB\t"+sib+"\t"+thisp+"\t\t"+rates.getEdgeLength(sib)+"\t"+rates.getLossParameter(sib)+"\t\t"+grad[3*sib + PARAMETER_LOSS]+"\t\t"+rates.toString(sib)+"\t"+gradient.factory.toString(sib));			
			}

			fixLoss(node, false);
			initModelParameters();
			setParameterValues(model0); // and reset
		}
	}
	
	
//	private void debugGradient()
//	{
//		double[] x=getParameterValues();
//		
//		double[] df = optDistrDiff(x);
//		
//		double[] df_est = FunctionMinimization.numericalGradient(θ->optDistrFunc(θ), x);
//		
//		final int n = x.length;
//		
//		double[] xafter = getParameterValues();
//		
//		for (int i=0; i<n; i++)
//		{
//			ModelParameter P = distribution_params.get(i);
//			double df_delta = df_est[i]-df[i];
//			double rel_delta = Math.abs(df_delta/df[i]);
//			
//			System.out.println("MLD.dG param "+i+"\t"+P+"\tdfdx "+df[i]+"\test "+df_est[i]
//					+"\tdiff "+df_delta+"\t("+rel_delta+")"+"\tx0 "+x[i]+"\tx1 "+xafter[i]);
//		}
//	}
	
	
	
	
	public static void main(String[] args) throws Exception
	{
		PRINT_OPTIMIZATION_MESSAGES = true;
		
		Class<?> us = java.lang.invoke.MethodHandles.lookup().lookupClass();
		count.io.CommandLine cli = new count.io.CommandLine(args, us);

		PrintStream out = System.out;
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(us));
    	    out.println(CommandLine.getStandardRuntimeInfo(us, args));
    	}
    	
    	MixedRateModel model = null; 
    	TreeWithRates starting_rates;
		count.matek.DiscreteDistribution root_prior = cli.getOptionDistribution(OPT_MODEL_ROOT_PRIOR, null);
    	if (cli.getMixedrateModel()==null)
    	{
    		
    		
    		Random RND = cli.getOptionRND(out);
//    		if (cli.getOptionValue(OPT_RND)!=null)
//    		{
//    			int rnd_seed = cli.getOptionInt(OPT_RND, 0);
//    			RND = (rnd_seed==0?new Random():new Random(rnd_seed));
//    			out.println(CommandLine.getStandardHeader("Random initialization: -"+OPT_RND+" "+rnd_seed));    			
//    		}
    		model =  GammaInvariant.nullModel(cli.getTree(), RND);
    		starting_rates = model.getBaseModel();
//    		model = GammaInvariant.nullModel(cli.getTree());
//    		starting_rates = model.getBaseModel();
    		if (root_prior != null)
    		{
    			starting_rates.setRootDistribution(root_prior);
    			out.println(CommandLine.getStandardHeader("Root prior preset: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
    		}
    		starting_rates.initNodeParameters(starting_rates.getTree().getRoot());
			out.println(CommandLine.getStandardHeader("(Root prior random: "+starting_rates.getRootDistribution()+")"));
    		
    	} else
    	{
    		starting_rates = cli.getRates(); 
    		assert (starting_rates != null);
    		RateVariationModel constant_rates = new RateVariationModel(starting_rates);
    		constant_rates.initConstantRates();
    		model = constant_rates;
//    		// set to a single class 
//    		if (1<model.getNumActiveClasses())
//    		{
//	    		model.setDuplicationForbidden(0.0);
//	    		model.setGainForbidden(0.0);
//	    		model.setClasses(1, 1, 1, 1);
//    		}
//    		assert starting_rates == model.getBaseModel();
    		Random RND = cli.getOptionRND(out);
    		if (RND!=null)
    			starting_rates.setRandom(RND);
    		if (root_prior != null)
    		{
    			starting_rates.setRootDistribution(root_prior);
    			out.println(CommandLine.getStandardHeader("Root prior preset: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
    		}
    	}
    	
    	AnnotatedTable table = cli.getTable();
		MLDistribution O = new MLDistribution(starting_rates, table);
		
		if (cli.getOptionValue(OPT_REINIT)!=null)
		{
			int reinit = cli.getOptionInt(OPT_REINIT, REINIT_EXTREME_EDGES);
			if (reinit == 0 || reinit == 1 || reinit == 3 || reinit == 5)
			{
				REINIT_EXTREME_EDGES = reinit;
			} else
			{
				throw new IllegalArgumentException("-"+OPT_REINIT+" allowed values: 0,1,3, and 5 but not "+reinit);
			}
			out.println(CommandLine.getStandardHeader("Extreme edge parameter reset neighborhood: -"+OPT_REINIT+" "+reinit)); 
		}
		
		
		String opt_debug_gradient = "debug.gradient";
		if (cli.getOptionValue(opt_debug_gradient)!=null)
		{
			boolean want_debug = cli.getOptionBoolean(opt_debug_gradient, false);
			DEBUG_GRADIENT = want_debug;
		}
		
		
		int absolute = AUTO_TRUNCATE?DEFAULT_TRUNCATE_ABSOLUTE:Integer.MAX_VALUE;
		double relative = AUTO_TRUNCATE?DEFAULT_TRUNCATE_RELATIVE:Double.POSITIVE_INFINITY;
        if (cli.getOptionValue(OPT_TRUNCATE)!=null)
        {
        	String truncate_val = cli.getOptionValue(OPT_TRUNCATE);
        	if (truncate_val.endsWith("auto"))
        	{
        		if ("auto".equals(truncate_val))
        		{
        			AUTO_TRUNCATE = true;
        			absolute = DEFAULT_TRUNCATE_ABSOLUTE;
        			relative = DEFAULT_TRUNCATE_RELATIVE;
        		} else if ("noauto".equals(truncate_val))
        		{
        			AUTO_TRUNCATE = false;
        			absolute = Integer.MAX_VALUE;
        			relative= Double.POSITIVE_INFINITY;
        		} else
        		{
        			throw new IllegalArgumentException("Use -"+OPT_TRUNCATE+" with auto or noauto [got "+truncate_val+"]");
        		}
        	} else
        	{
        		absolute = cli.getOptionTruncateAbsolute();
        		relative = cli.getOptionTruncateRelative();
	        	AUTO_TRUNCATE = false;
	        	DEFAULT_TRUNCATE_ABSOLUTE = Integer.min(absolute, DEFAULT_TRUNCATE_ABSOLUTE);
	        	DEFAULT_TRUNCATE_RELATIVE = Double.min(relative, DEFAULT_TRUNCATE_RELATIVE);
        	}
        } 
    	O.setCalculationWidth(absolute, relative);
        out.println(CommandLine.getStandardHeader("Truncated computation: -"
        		+OPT_TRUNCATE+" "+absolute+","+relative)
        		+" ; auto-truncation "+(AUTO_TRUNCATE?"on":"off"));

        int min_copies = Integer.min(2, cli.getTable().minCopies());
        min_copies = cli.getOptionInt(OPT_MINCOPY, min_copies);
		O.setMinimumObservedCopies(min_copies);
		out.println(CommandLine.getStandardHeader("Minimum observed copies: -"+OPT_MINCOPY+" "+min_copies));
		
		int maxiter = cli.getOptionInt(OPT_ROUNDS, Integer.MAX_VALUE);
		double eps = cli.getOptionDouble(OPT_EPS, 0.00125);
        
        out.println(CommandLine.getStandardHeader("Convergence: -"+OPT_ROUNDS+" "+maxiter+" -"+OPT_EPS+" "+eps));
        
        { // testing numerical optimization strategies 
	        O.optimization_by_quasiNewton = !(cli.getOptionBoolean("opt.cg", !O.optimization_by_quasiNewton));
	        O.do_gradient_descent = cli.getOptionBoolean("opt.gd", O.do_gradient_descent);
	        O.do_em = cli.getOptionBoolean("opt.EM", O.do_em);
	        out.println(CommandLine.getStandardHeader("Numerical optimization method: -opt.EM "+O.do_em+" -opt.gd "+O.do_gradient_descent+" -opt.cg "+!O.optimization_by_quasiNewton));
	        
	        boolean doarmijo = cli.getOptionBoolean("opt.armijo", FunctionMinimization.GD_ARMIJO);
	        if (doarmijo != FunctionMinimization.GD_ARMIJO && O.do_gradient_descent)
	        {
	        	FunctionMinimization.GD_ARMIJO = doarmijo;
	        	out.println(CommandLine.getStandardHeader("Gradient descent: -opt.armijo "+doarmijo));
	        }
	        boolean bfgsupdate = cli.getOptionBoolean("opt.bfgs", FunctionMinimization.DFP_BFGS_UPDATE);
	        if (bfgsupdate != FunctionMinimization.DFP_BFGS_UPDATE && O.optimization_by_quasiNewton)
	        {
	        	FunctionMinimization.DFP_BFGS_UPDATE = bfgsupdate;
	        	out.println(CommandLine.getStandardHeader("Quasi-Newton: -opt.bfgs "+bfgsupdate));
	        }
	        O.is_duprate_bounded =  cli.getOptionBoolean("opt.dupbound", O.is_duprate_bounded);
        }
        
        
        
//    	int gain_type = cli.getOptionInt(OPT_GAIN, PARAMETER_DUPLICATION);
//    	out.println(CommandLine.getStandardHeader("Gain parameter: -"+OPT_GAIN+" "+gain_type+" ("+ GLDParameters.paramName(gain_type)+")"));
//        O.setGainParameterType(gain_type);
        
        int testpnode = cli.getOptionInt("testp", -1);
        int testrnode = cli.getOptionInt("testr", -1);
        
        double root_gain_region = cli.getOptionDouble("lregion", 0.0);
        if (root_gain_region != 0.0)
        {
        	testrnode = starting_rates.getTree().getRoot();
        }
        
        double score = O.optimize(eps, (testpnode == -1&&testrnode==-1?maxiter:0));
        double bic_pty = 0.5*O.getModelParameterCount()*Math.log(table.getFamilyCount());
		
		out.println("#SCORE "+score+"\tBICpty "+bic_pty+"\tregularized "+(score+bic_pty));
        
		out.println("#TREE "+NewickParser.printTree(cli.getTree()));
		
		if (DEBUG_GRADIENT && maxiter==0) O.debugGradient(out);
		
//		if (PURE_LOGISTIC)
//		{
//			for (int node=0; node<O.rates.getTree().getNumNodes(); node++)
//				System.out.println("#**MLD.main rates "+O.rates.toString(node));
//		}
		

		if (testpnode != -1 || testrnode != -1)
		{
			if (testpnode != -1)
				O.testEdgeLoss(out, testpnode, eps, maxiter);
			
			if (testrnode != -1)
			{
				if (root_gain_region == 0.0)
					O.testEdgeGain(out, testrnode, eps, maxiter);
				else
				{
					//double r0 = starting_rates.getGainParameter(testrnode);
					O.gainLikelihoodInterval(out, testrnode, root_gain_region, eps, maxiter);
				}
			}
		}
		if ((testpnode == -1 && testrnode == -1) || root_gain_region != 0.0)
		{
			out.println(count.io.RateVariationParser.printRates(model));
		}

		
		if (out != System.out)
		{
			if (out.checkError())
				throw new java.io.IOException("Write failed.");
			out.close();
		}

	}
	
}
