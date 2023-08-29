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
import static count.model.GLDParameters.PARAMETER_LENGTH;
import static count.model.GLDParameters.PARAMETER_LOSS;
import static count.model.GLDParameters.PARAMETER_DUPLICATION_COMPLEMENT;
import static count.model.GLDParameters.PARAMETER_LOSS_COMPLEMENT;
//import static count.model.GLDParameters.PARAMETER_LENGTH;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import count.Count;
import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.io.CommandLine;
import count.io.NewickParser;
import count.matek.FunctionMinimization;

import static count.io.CommandLine.OPT_EPS;
import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_OUTPUT;
import static count.io.CommandLine.OPT_ROUNDS;
import static count.io.CommandLine.OPT_TRUNCATE;

import static count.io.CommandLine.OPT_MODEL_ROOT_PRIOR;

/**
 * Maximum likelihood for one-class model 
 * using distribution parameters (loss probability, duplication probability and gain rate). 
 * 
 * @author csuros
 *
 */
public class MLDistribution extends ML implements Count.UsesThreadpool // via Gradient
{
//	public static final int DEFAULT_OPT_ROUNDS = 100;
//	public static final double DEFAULT_OPT_EPS = 1e-2;
	
	
	public static boolean PRINT_OPTIMIZATION_MESSAGES = false;
	
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
		this.rates = rates; //.parameterCache();
		this.gradient = new Gradient(this.rates, table);
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
	
	private boolean is_duprate_bounded = false; // if bounded, <1.0 is enforced
	
	private boolean track_complements = false; // if true, high-precision calculation is attempted (numerically buggy) 
	
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
//		// DEBUG
//		System.out.println("#**MLD.fL "+node+"\t"+optimize_parameters[3*node+PARAMETER_LOSS]);
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
				params[PARAMETER_LOSS] = rates.getLossParameter(node);
				params[PARAMETER_LOSS_COMPLEMENT] = rates.getLossParameterComplement(node);
				if (is_duprate_bounded)
				{
					params[PARAMETER_DUPLICATION] = rates.getDuplicationRate(node);
					params[PARAMETER_DUPLICATION_COMPLEMENT] = rates.getRateGap(node);
				} else
				{
					params[PARAMETER_DUPLICATION] = rates.getDuplicationParameter(node);
					params[PARAMETER_DUPLICATION_COMPLEMENT] = rates.getDuplicationParameterComplement(node);
				}
				node_parameters[node] = params;
				// add them 
				if (optimize_gain)
				{
//					distribution_params.add(new Logarithmic(new GainParameter(node)));
					//addLogarithmic(new GainParameter(node));
					//addLogistic(new GainParameter(node), MAX_RATE);
					distribution_params.add(newLogistic(new GainParameter(node), MAX_RATE));
//					ModelParameter gpar = new Logistic(new GainParameter(node));
//					distribution_params.add(gpar);
				}
				if (optimize_duplication) 
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
						}
					}
					
//					ModelParameter dpar = new Logistic(new DuplicationParameter(node), 1.0-1e-9);
//					distribution_params.add(dpar);
				}
				if (optimize_loss) 
				{
					ModelParameter lpar;
					if (track_complements && !is_duprate_bounded)
					{
						lpar = new BoundedLogistic(new LossParameter(node));
					} else
					{
						lpar = newLogistic(new LossParameter(node), MAX_PROB_NOT1);
					}
					// lpar = bracketedLogistic(new LossParameter(node),prob_small );
					distribution_params.add(lpar);
					// DEBUG
					System.out.println("#**MLD.iMP "+node+"\tloss "+lpar);
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
				double κ = node_parameters[node][PARAMETER_GAIN];
				double p = node_parameters[node][PARAMETER_LOSS];
				
				if (is_duprate_bounded)
				{
					double λ = node_parameters[node][PARAMETER_DUPLICATION];
					if (track_complements)
					{
						double λcomp = node_parameters[node][PARAMETER_DUPLICATION_COMPLEMENT];
						rates.setDuplicationRate(node, λ, λcomp);
						double p_1 = node_parameters[node][PARAMETER_LOSS_COMPLEMENT];
						rates.setEdgeLengthForLossParameter(node, p, p_1);
					} else
					{
						rates.setDuplicationRate(node, λ);
						rates.setEdgeLengthForLossParameter(node, p, 1.0-p);
					}
					rates.setGainRateForLossParameter(node, κ, p);
				} else
				{
					double q = node_parameters[node][PARAMETER_DUPLICATION];
					if (track_complements)
					{
						double p_1 = node_parameters[node][PARAMETER_LOSS_COMPLEMENT];
						double q_1 = node_parameters[node][PARAMETER_DUPLICATION_COMPLEMENT];
						if (q==1.0)
						{
							System.out.println("#**MLD.cNP node "+node+"\tp" +p+"\t1-"+p_1+"\tq "+q+"\t1-"+q_1+"\t"+rates.toString(node));
						}
						
						rates.setDuplicationLossRates(node, p, p_1, q, q_1);
						rates.setGainRateForLossParameter(node, κ, p);
					} else
					{
						rates.setParameters(node, κ, p, q);	
					}
//				System.out.println("#**MLD.cNPTM "+node+"\t"+Arrays.toString(node_parameters[node])+"\t"+rates.toString(node));
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
//		System.out.println("#*MLD.oDF "+calls_optFunc+"\t"+LL);

		
		return -LL;
	}
	
	
	private int calls_optDiff = 0;
	
	private double[] parameterGradient(double[] distribution_gradient)
	{
		double[] D = new double[distribution_params.size()];
		for (int j=0; j<D.length; j++)
			D[j] = -distribution_params.get(j).dL(distribution_gradient);
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
			double κ = rates.getGainParameter(node);
			node_parameters[node][PARAMETER_GAIN]=κ;
		}
		private final int node;
		
		@Override
		public double get()
		{
//			double κ = rates.getGainParameter(node);
//			node_parameters[node][PARAMETER_GAIN]=κ;
			double κ = node_parameters[node][PARAMETER_GAIN];
			return κ;
		}
		
		@Override
		public void set(double κ)
		{
			if (!Double.isFinite(κ)) 
			{
				System.out.println("#**MLD.GP.set r="+κ+"\tnpars "+Arrays.toString(node_parameters[node])+"\t"+rates.toString(node)+"\t// "+gradient.factory.toString(node));
			}
			assert Double.isFinite(κ);
			assert (0.0<=κ);
			node_parameters[node][PARAMETER_GAIN]=κ;
			// rates.setParameters(node, κ, node_parameters[node][PARAMETER_LOSS], node_parameters[node][PARAMETER_DUPLICATION]);
		}
		
		@Override 
		public double dL(double[] gradient)
		{
			return gradient[3*node+PARAMETER_GAIN];
		}
		@Override 
		public String toString()
		{
			return "r"+node+"="+get();
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
	
	
//	public void Mstep()
//	{
//		int num_nodes = rates.getTree().getNumNodes();
//		double[] Nmeans = new double[]
//	}

	
//	public void set(double[] x)
//	{
//		int num_nodes = rates.getTree().getNumNodes();
//		boolean[] setnode = new boolean[num_nodes];
//		double[] node_params = new double[num_nodes*3];
//		
//		for (int j=0; j<parameter_index.length; j++)
//		{
//			int i = parameter_index[j];
//			int node = i/3;
//			int type = i % 3;
//			if (!setnode[node])
//			{
//				setnode[node]= true;
//				double κ = rates.getGainParameter(node);
//				double p = rates.getLossParameter(node);
//				double q = rates.getDuplicationParameter(node);
//				node_params[3*node+PARAMETER_GAIN] = κ;
//				node_params[3*node+PARAMETER_LOSS] = p;
//				node_params[3*node+PARAMETER_DUPLICATION] = q;
//			}
//			
//			if (type == PARAMETER_GAIN)
//			{
//				double κ = Math.exp(x[j]);
//				node_params[i] = κ;
//			} else if (type == PARAMETER_LOSS)
//			{
//				double p = 1.0/(1.0+Math.exp(-x[j]));
//				node_params[i] = p;
//			} else if (type == PARAMETER_DUPLICATION)
//			{
//				double q = 1.0/(1.0+Math.exp(-x[j]));
//				node_params[i] = q;
//			}
//		}
//		for (int node=0; node<num_nodes; node++)
//			if (setnode[node])
//			{
//				double κ = node_params[3*node+PARAMETER_GAIN];
//				double p = node_params[3*node + PARAMETER_LOSS];
//				double q = node_params[3*node + PARAMETER_DUPLICATION];
//				
//				rates.setParameters(node, κ, p, q);
//			}
//		
//		gradient.factory.computeParameters();
////		System.out.println(count.io.RateVariationParser.printRates(rates));
//	}
//	
//	public double[] get()
//	{
//		double[] x = new double[parameter_index.length];
//		for (int j=0; j<parameter_index.length; j++)
//		{
//			int i = parameter_index[j];
//			int node = i/3;
//			int type = i % 3;
//			if (type == PARAMETER_GAIN)
//			{
//				double κ = rates.getGainParameter(node); // logarithmic transformation
//				x[j] = Math.log(κ);
//			} else if (type == PARAMETER_LOSS)
//			{
//				double p = rates.getLossParameter(node);
//				x[j] = Math.log(p)-Math.log1p(-p); // logistic transformation
//			} else if (type == PARAMETER_DUPLICATION)
//			{
//				double q = rates.getDuplicationParameter(node);
//				x[j] = Math.log(q)-Math.log1p(-q);
//			}
//		}
//		return x;
//	}
//
	@Override
	public double optimize(double delta)
	{
		return this.optimize(delta, Integer.MAX_VALUE);
	}
	
	@Override
	public double optimize(double delta, int itmax)
	{
		if (PRINT_OPTIMIZATION_MESSAGES)
		{
			double LL =  gradient.getCorrectedLL();
			System.out.println("#*MLD.o before LL "+LL+"\t; optimization with dup rates < 1.0 "+is_duprate_bounded);
		}
		initModelParameters();
		List<Double> history = getOptimizationHistory();
		if (history == null)
		{
			history = new ArrayList<>();
			setOptimizationHistory(history);
		}
		
		int h0 = history.size();
		
		double LL =  gradient.getCorrectedLL();
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#*MLD.o starting LL "+LL+"\tparams "+distribution_params.size());
		
		while (itmax!=0 && history.size()<itmax+h0)
		{
			double[] x0 = getParameterValues();
			double[] old_x0 = x0.clone();
			
			double min;
			
			if (do_em)
			{
				EM O = new EM(gradient);
				min = O.optimize(delta, x0.length);
				initModelParameters(); // recalculate from rates set by EM
				x0 = getParameterValues(); // starting point for further optimization
				//gradient.computeParameters();
			}
			if (do_gradient_descent)
			{
				//min = FunctionMinimization.powell(x0, delta, itmax, x->optDistrFunc(x), history);
				min = FunctionMinimization.gradientDescent(x0, Math.sqrt(delta), x0.length, x->optDistrFunc(x), x->optDistrDiff(x), history);
				setParameterValues(x0); // x0 is updated by dfpmin

				LL =  gradient.getCorrectedLL(); // == -min 
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#*MLD.o GD "+LL+"\t("+min+")"+"\tfcalls "+calls_optFunc+"\tdcalls "+calls_optDiff);
				
			}
			
			if (optimization_by_quasiNewton)
			{					
				 min  = FunctionMinimization.dfpmin(x0, delta, itmax, x->optDistrFunc(x), x->optDistrDiff(x), history);
			}
			else
			{
				min = FunctionMinimization.frprmn(x0, delta, itmax, x->optDistrFunc(x), x->optDistrDiff(x), history);
			}
			
			setParameterValues(x0); // x0 is updated by dfpmin

			LL =  gradient.getCorrectedLL(); // == -min 

			double max_dx = 0.0;
			for (int j=0; j<x0.length;j++)
			{
				double dx = Math.abs(x0[j]-old_x0[j])/Double.max(1.0, Math.abs(x0[j]));
				max_dx = Double.max(dx,max_dx);
			}
			
			
			double D[] = optDistrDiff(x0);
			double gradient_length = FunctionMinimization.euclideanNorm(D);
//			double L0 = gradient.factory.getEmptyLL();
//			double L1 = gradient.factory.getSingletonLL();
			double rel_gradient = gradient_length/(-LL);
			
			boolean done_converged = (rel_gradient<delta || max_dx < FunctionMinimization.DFP_TOLX);

			int iter_done = history.size()-h0;
			
			Count.out.println("#*MLD.o iter LL "+LL+"\t("+min+")"+"\tfcalls "+calls_optFunc+"\tdcalls "+calls_optDiff
						//+"\tL0 "+L0+"\tL1 "+L1
						+"\tdL "+gradient_length
						+"\trgrad "+rel_gradient+"\tmax_dx "+max_dx+"\t"+(done_converged?"DONE":"loop")
						+"\titers "+iter_done);
			
			if (done_converged || true) break;
			//optimization_by_quasiNewton = !optimization_by_quasiNewton; // alternate between methods 
			
//			rates.setRates();
		}
		
				
		return -LL;
	}
	
	/**
	 * Which of the 3/node parameters should be optimized. 
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
			if (tree.isRoot(node))
			{
				opt_par[3*node+PARAMETER_GAIN] = rates.getGainParameter(node)>0.0;
				opt_par[3*node+PARAMETER_LOSS] = (rates.getLossParameter(node)<1.0);
				opt_par[3*node+PARAMETER_DUPLICATION] = rates.getDuplicationParameter(node)>0.0;
//				// DEBUG
//				System.out.println("#**MLD.oP root "+node
//						+"\tog "+opt_par[3*node+PARAMETER_GAIN]
//						+"\tod "+opt_par[3*node+PARAMETER_DUPLICATION]
//						+"\tol "+opt_par[3*node+PARAMETER_LOSS]
//						+"\t// "+rates.toString(node)
//						);
			} else
			{
				opt_par[3*node + PARAMETER_GAIN] = rates.getGainParameter(node)>0.0;
				opt_par[3*node+PARAMETER_LOSS] = (rates.getLossParameter(node)<1.0);
				opt_par[3*node+PARAMETER_DUPLICATION] = rates.getDuplicationParameter(node)>0.0;

			} 
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

	private void debugGradient(PrintStream out)
	{
		optimize(0,0); 
		double[] x = getParameterValues();
		
		double[] df = optDistrDiff(x);
		
		double[] df_est = FunctionMinimization.numericalGradient(θ->optDistrFunc(θ), x);// estimateDistrDiff(x); // FunctionMinimization.numericalGradient(θ->optDistrFunc(θ), x);
		
		double[] xafter = getParameterValues();

		for (int i=0; i<x.length; i++)
		{
			ModelParameter P = distribution_params.get(i);
			
			double df_delta = df_est[i]-df[i];
			double rel_delta = Math.abs(df_delta/df[i]);
			
			out.println("#**MLD.dG param "+i+"\t"+P+"\tdfdx "+df[i]+"\test "+df_est[i]
					+"\tdiff "+df_delta+"\t("+rel_delta+")"+"\tx0 "+x[i]+"\tx1 "+xafter[i]);
		}
		
		
	}
	
	private void testEdgeLoss(int node, double delta, int itmax)
	{
		double p = rates.getLossParameter(node);
		double q = rates.getDuplicationParameter(node);
		double q1 = rates.getDuplicationParameterComplement(node);

		List<Double> ptest = new ArrayList<>();
		
		ptest.add(p);
		double sqp = Math.sqrt(p);
		while (sqp<1e-5)
		{
			ptest.add(sqp);
			sqp = Math.sqrt(sqp);
		}
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
		ptest.add(0.4);
		ptest.add(0.5);
		ptest.add(0.8);
		ptest.add(0.9);
		ptest.add(0.96);
		ptest.add(0.98);
		ptest.add(0.99);
		ptest.add(0.995);
		ptest.add(0.999);
		ptest.add(0.9999);
		ptest.add(0.99999);
		
		java.util.Collections.sort(ptest);

		
		final int ntest = ptest.size();
		System.out.println("#**MLD.tEL node\t"+node+"\tntest "+ntest+"\t"+ptest.toString()+"\t// "+ rates.toString(node));
		
		double negLL = optimize(delta, itmax);
		double[] model0 = getParameterValues();
		double[] grad = gradient.getDistributionGradient(gradient.getCorrectedGradient());
		double glen = FunctionMinimization.euclideanNorm(grad); 
		double dLdp = grad[3*node + PARAMETER_LOSS];
		double[] pgrad = parameterGradient(grad);
		double dlen = FunctionMinimization.euclideanNorm(pgrad);
//		double dLdθ = pgrad[3*node+PARAMETER_LOSS];
		
		
		
		double t = rates.getEdgeLength(node);
		System.out.println("#NODELOSS0\t"+node+"\t"+p+"\t"+negLL+"\t"+t+"\t"+p+"\t"+dlen+"\t"+dLdp+"\t"+glen);
		
		for (int round=0; round<ntest; round++)
		{
			double thisp = ptest.get(round);
			rates.setEdgeLengthForLossParameter(node, thisp, 1.0-thisp);
			t = rates.getEdgeLength(node);
			double setp = rates.getLossParameter(node);
			
			fixLoss(node, true);
			// gradient.computeParameters(); // not needed cc optimize() will do that
			negLL = optimize(delta, itmax);

			grad = gradient.getDistributionGradient(gradient.getCorrectedGradient());
			dLdp = grad[3*node + PARAMETER_LOSS];
			pgrad = parameterGradient(grad);

			dlen = FunctionMinimization.euclideanNorm(pgrad);
			glen = FunctionMinimization.euclideanNorm(grad);
			
			System.out.println("#NODELOSS\t"+node+"\t"+thisp+"\t"+negLL+"\t"+t+"\t"+setp+"\t"+dlen+"\t"+dLdp+"\t"+glen);

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
		
		count.io.CommandLine cli = new count.io.CommandLine(args, MLDistribution.class);

		PrintStream out = System.out;
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(MLDistribution.class));
    	    out.println(CommandLine.getStandardRuntimeInfo(MLDistribution.class, args));
    	}
    	
    	GammaInvariant model = null; 
    	TreeWithRates starting_rates;
    	if (cli.getModel()==null)
    	{
    		Random RND = cli.getOptionRND(out);
//    		if (cli.getOptionValue(OPT_RND)!=null)
//    		{
//    			int rnd_seed = cli.getOptionInt(OPT_RND, 0);
//    			RND = (rnd_seed==0?new Random():new Random(rnd_seed));
//    			out.println(CommandLine.getStandardHeader("Random initialization: -"+OPT_RND+" "+rnd_seed));    			
//    		}
    		model = GammaInvariant.nullModel(cli.getTree(), RND);
    		starting_rates = model.getBaseModel();
//    		model = GammaInvariant.nullModel(cli.getTree());
//    		starting_rates = model.getBaseModel();
    		
    	} else
    	{
    		starting_rates = cli.getRates(); 
    		model = cli.getModel();
    	}

    	{
    		count.matek.DiscreteDistribution root_prior = cli.getOptionDistribution(OPT_MODEL_ROOT_PRIOR, null);
    		if (root_prior != null)
    		{
    			starting_rates.setRootDistribution(root_prior);
    			out.println(CommandLine.getStandardHeader("Root prior set: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
    		}
			out.println(CommandLine.getStandardHeader("(Root prior get: "+starting_rates.getRootDistribution()+")"));
    	}
    	
    	AnnotatedTable table = cli.getTable();
		MLDistribution O = new MLDistribution(starting_rates, table);
		
		int absolute = 12;
		double relative = 3.0;
        if (cli.getOptionValue(OPT_TRUNCATE)!=null)
        {
        	absolute = cli.getOptionTruncateAbsolute();
        	relative = cli.getOptionTruncateRelative();
        } 
    	O.setCalculationWidth(absolute, relative);
        out.println(CommandLine.getStandardHeader("Truncated computation: -"
        		+OPT_TRUNCATE+" "+absolute+","+relative));

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
        }
        
        double score = O.optimize(eps, maxiter);
        double bic_pty = 0.5*O.getModelParameterCount()*Math.log(table.getFamilyCount());
		
		out.println("#SCORE "+score+"\tBICpty "+bic_pty+"\tregularized "+(score+bic_pty));
        
		out.println("#TREE "+NewickParser.printTree(cli.getTree()));
		
		if (maxiter==0) O.debugGradient(out);
		
		out.println(count.io.RateVariationParser.printRates(model));
		if (out != System.out)
		{
			if (out.checkError())
				throw new java.io.IOException("Write failed.");
			out.close();
		}

		{
			int testnode = cli.getOptionInt("testp", -1);
			if (testnode != -1)
				O.testEdgeLoss(testnode, eps, maxiter);
		}
		
		
	}
	
}
