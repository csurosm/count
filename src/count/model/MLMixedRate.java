package count.model;
/*
 * Copyright 2022 Mikl&oacute;s Cs&#369;r&ouml;s.
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
import static count.model.GLDParameters.PARAMETER_LENGTH;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;
import java.util.function.DoubleFunction;
import java.util.List;


import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.matek.FunctionMinimization;

/**
 * Maximum likelihood for multi-class model with rate parameters and multipliers. 
 * 
 * @author csuros
 *
 */
public class MLMixedRate extends ML
{
	protected static boolean PRINT_OPTIMIZATION_MESSAGES = false;
	/**
	 * Bad experimental code.  
	 */
	private static final boolean OPTIMIZE_BY_DISTRIBUTION_PARAMETERS = false;
	
	public MLMixedRate(MixedRateModel.RateMultipliers mixed_model, ProfileTable table)
	{
		this(mixed_model, table, optimizableRates(mixed_model.getBaseModel()));
	}

	public MLMixedRate(MixedRateModel.RateMultipliers mixed_model, ProfileTable table, boolean[] optimize_rates)
	{
		this.init_mixed_model = mixed_model;
		this.init_table = table;

		this.base = mixed_model.getBaseModel();
		assert (optimize_rates.length == 4*base.getTree().getNumNodes());
		this.optimize_rates = optimize_rates;
		this.rate_params = new ArrayList<>();
		this.node_distribution_parameters = new double[base.getTree().getNumNodes()][];
		

		this.gradient = new MixedRateGradient(mixed_model, table);
	}
	
	private final MixedRateModel.RateMultipliers init_mixed_model;
	private final ProfileTable init_table;
	
	protected final MixedRateGradient gradient;
	private final TreeWithRates base;
	private final List<ModelParameter> rate_params;
	/**
	 * p,q,r parameters at each node for base model
	 */
	private final double[][] node_distribution_parameters;
	
	
	private boolean[] optimize_rates;
	private boolean uniform_gain=false;
	private boolean uniform_duplication=false;
	private boolean uniform_loss=false;
	private boolean stationary_root=false;
	
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
	@Override
	public void setCalculationWidth(int absolute, double relative)
	{
		gradient.setCalculationWidthThresholds(absolute, relative);
	}
	
	
	public void fixGain(int node)
	{
		optimize_rates[4*node+PARAMETER_GAIN]=false;
	}
	
	public void fixLoss(int node)
	{
		optimize_rates[4*node+PARAMETER_LOSS]=false;
	}
	public void fixDuplication(int node)
	{
		optimize_rates[4*node+PARAMETER_DUPLICATION]=false;
	}
	public void fixLength(int node)
	{
		optimize_rates[4*node+PARAMETER_LENGTH]=false;
	}
	
	@Override
	public void fixNodeParameters(int node, boolean do_not_optimize)
	{
		if (do_not_optimize)
		{
			fixGain(node);
			fixLoss(node);
			fixDuplication(node);
			fixLength(node);
		}
	}
	
	public void setUniformGain(boolean uniform) 
	{
		uniform_gain=uniform;
	}
	public void setUniformLoss(boolean uniform) 
	{
		uniform_loss=uniform;
	}
	public void setUniformDuplication(boolean uniform) 
	{
		uniform_duplication=uniform;
	}
	public void setStationaryRoot(boolean stationary)
	{
		stationary_root=stationary;
	}
	@Override
	public void setMinimumObservedCopies(int num_copies)
	{
		gradient.setMinimumObservedCopies(num_copies);
	}

//	private ModelParameter addLogistic(ModelParameter theta)
//	{
//		return addLogistic(theta, MAX_RATE);
//	}
//	private ModelParameter addLogistic(ModelParameter theta, double max_rate)
//	{
////		double θ = theta.get();
////		if (θ>=max_rate)
////		{
////			double newθ = (1.-1e-9)*max_rate;
////			System.out.println("#*MLMR.addLogist "+theta+"; resetting to "+newθ);
////			theta.set(newθ);
////		}
//		ModelParameter transformed = newLogistic(theta, max_rate);
//		rate_params.add(transformed); 
//		// new Logarithmic(theta); // new LogLogistic(theta);
//		return transformed;
//	}
	
//	private ModelParameter addLogarithmic(ModelParameter theta)
//	{
//		Logarithmic L = new Logarithmic(theta);
//		rate_params.add(L);
//		return L;
//	}
	
	public int getModelParameterCount()
	{
		return rate_params.size();
	}
	
	protected ModelParameter getRateParameter(int par_idx)
	{
		return rate_params.get(par_idx);
	}
	
	private void initModelParameters()
	{
		rate_params.clear();
		Arrays.fill(node_distribution_parameters, null);
		
		IndexedTree phylo = base.getTree();
		int root = phylo.getRoot();
		int num_nodes = phylo.getNumNodes();
		
		if (!OPTIMIZE_BY_DISTRIBUTION_PARAMETERS || uniform_gain || uniform_loss || uniform_duplication)
		{
			if (uniform_gain)
			{
				rate_params.add(new Logarithmic(new UniformGain())); //addLogistic(new UniformGain());
				if (!stationary_root && optimize_rates[4*root+PARAMETER_GAIN])
				{
					rate_params.add(new Logarithmic(new GainRate(root)));
					// addLogistic(new GainRate(root));
				}
			} else
			{
				for (int node = 0; node<num_nodes; node++)
					if (optimize_rates[4*node + PARAMETER_GAIN])
					{
//						addLogistic(new GainRate(node));
						rate_params.add(new Logarithmic(new GainRate(node)));
						// addLogarithmic(new GainRate(node));
					}
			}
			if (uniform_loss)
			{
				rate_params.add(newLogistic(new UniformLoss()));
				// addLogistic(new UniformLoss());
				if (optimize_rates[4*root+PARAMETER_LOSS])
				{
					rate_params.add(newLogistic(new LossRate(root))); //addLogistic(new LossRate(root));
				}
			} else
			{
				for (int node = 0; node<num_nodes; node++)
					if (optimize_rates[4*node + PARAMETER_LOSS])
					{
						rate_params.add(newLogistic(new LossRate(node)));
					}
			}
			if (uniform_duplication)
			{
				if (stationary_root)
				{
					rate_params.add(newLogistic(new UniformDuplication(), MAX_PROB_NOT1*base.getLossRate(0)));
					//addLogistic(new UniformDuplication(), (1.-1e-9)*base.getLossRate(0));
				} else
				{
					rate_params.add(newLogistic(new UniformDuplication()));
					//addLogistic(new UniformDuplication());
					if (optimize_rates[4*root+PARAMETER_DUPLICATION])
					{
						rate_params.add(newLogistic(new DuplicationRate(root), MAX_PROB_NOT1*base.getLossRate(root)));
						//addLogistic(new DuplicationRate(root), (1.-1e-9)*base.getLossRate(root));
					}
				}
			} else
			{
				for (int node = 0; node<num_nodes; node++)
					if (optimize_rates[4*node + PARAMETER_DUPLICATION])
					{
						if (phylo.isRoot(node) || base.getEdgeLength(node)==Double.POSITIVE_INFINITY)
						{
							rate_params.add(newLogistic(new DuplicationRate(node), MAX_PROB_NOT1*base.getLossRate(node)));
							//addLogistic(new DuplicationRate(node),(1.-1e-9)*base.getLossRate(node));
						} else
						{
							rate_params.add(newLogistic(new DuplicationRate(node)));
							//addLogistic(new DuplicationRate(node));
						}
					}
			}
			for (int node = 0; node<num_nodes; node++)
				if (optimize_rates[4*node + PARAMETER_LENGTH])
				{
					rate_params.add(newLogistic(new EdgeLength(node)));
					//addLogistic(new EdgeLength(node));
				}
		
		} else
		{ // experimental code: not by rates 
			for (int node = 0; node<num_nodes; node++)
			{
				boolean optimize_gain = optimize_rates[4*node+PARAMETER_GAIN];
				boolean optimize_loss = optimize_rates[4*node+PARAMETER_LOSS] || optimize_rates[4*node+PARAMETER_LENGTH];
				boolean optimize_duplication = optimize_rates[4*node+PARAMETER_DUPLICATION];

				if (optimize_gain || optimize_loss || optimize_duplication)
				{
					double[] params = new double[3];
					params[PARAMETER_GAIN] = base.getGainParameter(node);
					params[PARAMETER_LOSS] = base.getLossParameter(node);
					params[PARAMETER_DUPLICATION] = base.getDuplicationParameter(node);
					node_distribution_parameters[node] = params;
					// add them 
					if (optimize_gain)
					{
						rate_params.add(new Logarithmic(new GainParameter(node)));
					}
					if (optimize_loss)
					{
						rate_params.add(newLogistic(new LossParameter(node), MAX_PROB_NOT1));
					}
					if (optimize_duplication)
					{
						rate_params.add(newLogistic(new DuplicationParameter(node), MAX_PROB_NOT1));
					}
				}
			}
			// buggy: 
			throw new UnsupportedOperationException("Use rate parameters, not distribution parameters");
			
			// copyNodeParametersToModel();
		}
		
//		double[] x = getParameterValues();
////		setParameterValues(x);
//		for (int j=0; j<x.length; j++)
//		{
//			double diff = (rate_params.get(j).get()-x[j]);
//			System.out.println("#**MLMR.iMP "+j+"\t"+rate_params.get(j)
//				+"\t"+rate_params.get(j).get()+"\twas "+x[j]
//				+"\tdiff "+diff);
//		}

	}
	
	
	
	private class GainRate implements ModelParameter
	{
		private final int node;
		GainRate(int node)
		{
			this.node = node;
		}
		@Override
		public double get()
		{
			return base.getGainRate(node);
		}
		@Override
		public void set(double x)
		{
			base.setGainRate(node, x);
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dLdκ = rate_gradient[4*node+PARAMETER_GAIN];
			return dLdκ ;
		}
		@Override 
		public String toString()
		{
			return "g"+node+"="+get();
		}		
	}	
	
	private class DuplicationRate implements ModelParameter
	{
		private final int node;
		DuplicationRate(int node)
		{
			this.node = node;
		}
		@Override
		public double get()
		{
			double λ = base.getDuplicationRate(node);
			return λ;
		}
		@Override
		public void set(double x)
		{
			base.setDuplicationRate(node,x);
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dLdλ = rate_gradient[4*node+PARAMETER_DUPLICATION];
			return dLdλ ;
		}
		@Override 
		public String toString()
		{
			return "d"+node+"="+get();
		}
	}			
	
	private class EdgeLength implements ModelParameter
	{
		private final int node;
		EdgeLength(int node)
		{
			this.node = node;
		}
		@Override
		public double get()
		{
			return base.getEdgeLength(node);
		}
		@Override
		public void set(double x)
		{
			base.setEdgeLength(node, x);
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dLdt = rate_gradient[4*node+PARAMETER_LENGTH];
			return dLdt;
		}
		@Override 
		public String toString()
		{
			return "t"+node+"="+get();
		}
	}		

	private class LossRate implements ModelParameter
	{
		private final int node;
		LossRate(int node)
		{
			this.node = node;
		}
		@Override
		public double get()
		{
			return base.getLossRate(node);
		}
		@Override
		public void set(double x)
		{
			base.setLossRate(node, x);
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dLdμ = rate_gradient[4*node+PARAMETER_LOSS];
			return dLdμ ;
		}
		@Override 
		public String toString()
		{
			return "l"+node+"="+get();
		}
	}		
	
	private class UniformGain implements ModelParameter
	{
		@Override
		public double get()
		{
			return base.getGainRate(0);
		}
		@Override
		public void set(double x)
		{
			IndexedTree phylo = base.getTree();
			int num_nodes = phylo.getNumNodes();
			for (int node=0; node<num_nodes; node++)
			{
				if (stationary_root || !phylo.isRoot(node))
					base.setGainRate(node, x);
			}
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dL = 0.0;

			IndexedTree phylo = base.getTree();
			int num_nodes = phylo.getNumNodes();
			for (int node=0; node<num_nodes; node++)
				if (stationary_root || !phylo.isRoot(node))
				{
					double dLdκ = rate_gradient[4*node+PARAMETER_GAIN];
					dL += dLdκ;
				}
			return dL;
		}
		@Override 
		public String toString()
		{
			return "Ug="+get();
		}
	}
	
	private class UniformLoss implements ModelParameter
	{
		@Override
		public double get()
		{
			return base.getLossRate(0);
		}
		@Override
		public void set(double x)
		{
			IndexedTree phylo = base.getTree();
			int num_nodes = phylo.getNumNodes();
			for (int node=0; node<num_nodes; node++)
			{
				if (!phylo.isRoot(node))
					base.setLossRate(node, x);
			}
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dL = 0.0;

			IndexedTree phylo = base.getTree();
			int num_nodes = phylo.getNumNodes();
			for (int node=0; node<num_nodes; node++)
				if (!phylo.isRoot(node))
				{
					double dLdμ = rate_gradient[4*node+PARAMETER_LOSS];
					dL += dLdμ;
				}
			return dL;
		}
		@Override 
		public String toString()
		{
			return "Ul="+get();
		}
	}
	
	private class UniformDuplication implements ModelParameter
	{
		@Override
		public double get()
		{
			return base.getDuplicationRate(0);
		}
		@Override
		public void set(double x)
		{
			IndexedTree phylo = base.getTree();
			int num_nodes = phylo.getNumNodes();
			for (int node=0; node<num_nodes; node++)
			{
				if (stationary_root || !phylo.isRoot(node))
					base.setDuplicationRate(node, x);
			}
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dL = 0.0;

			IndexedTree phylo = base.getTree();
			int num_nodes = phylo.getNumNodes();
			for (int node=0; node<num_nodes; node++)
				if (stationary_root || !phylo.isRoot(node))
				{
					double dLdλ = rate_gradient[4*node+PARAMETER_DUPLICATION];
					dL += dLdλ;
				}
			return dL;
		}
		@Override 
		public String toString()
		{
			return "Ud="+get();
		}
	}
	
	private class GainParameter implements ModelParameter
	{
		GainParameter(int node)
		{
			this.node = node;
			double κ = base.getGainParameter(node);
			node_distribution_parameters[node][PARAMETER_GAIN]=κ;
		}
		private final int node;
		
		@Override
		public double get()
		{
//			double κ = rates.getGainParameter(node);
//			node_parameters[node][PARAMETER_GAIN]=κ;
			double κ = node_distribution_parameters[node][PARAMETER_GAIN];
			return κ;
		}
		
		@Override
		public void set(double κ)
		{
			if (!Double.isFinite(κ)) 
			{
				System.out.println("#**MLMR.GP.set r="+κ+"\tnpars "+Arrays.toString(node_distribution_parameters[node])+"\t"+base.toString(node));
			}
			assert Double.isFinite(κ);
			assert (0.0<=κ);
			node_distribution_parameters[node][PARAMETER_GAIN]=κ;
			base.setParameters(node, κ, node_distribution_parameters[node][PARAMETER_LOSS], node_distribution_parameters[node][PARAMETER_DUPLICATION]);
		}
		
		@Override 
		public double dL(double[] rateparams_gradient)
		{
			double dLdr = MLMixedRate.this.gradient.inferDistributionGradient(node, PARAMETER_LOSS, rateparams_gradient);
			return dLdr;
		}
		@Override 
		public String toString()
		{
			return "r"+node+"="+get();
		}
	}
	
	/**
	 * Parametrizing by the base model's loss distribution parameter.
	 * Call {@link MLMixedRate#copyNodeParametersToModel()} after {@link #set(double)}.  
	 * 
	 * @author csuros
	 *
	 */
	private class LossParameter implements ModelParameter
	{
		LossParameter(int node)
		{
			this.node = node;
			double p = base.getLossParameter(node);
			node_distribution_parameters[node][PARAMETER_LOSS]=p;
		}
		private final int node;
		
		@Override
		public double get()
		{
//			double p = rates.getLossParameter(node);
//			node_parameters[node][PARAMETER_LOSS]=p;
			double p = node_distribution_parameters[node][PARAMETER_LOSS];
			return p;
		}
		
		@Override
		public void set(double p)
		{
			node_distribution_parameters[node][PARAMETER_LOSS]=p;
//			rates.setParameters(node, node_parameters[node][PARAMETER_GAIN], p, node_parameters[node][PARAMETER_DUPLICATION]);
		}
		
		@Override 
		public double dL(double[] rateparams_gradient)
		{
			double dLdp = MLMixedRate.this.gradient.inferDistributionGradient(node, PARAMETER_LOSS, rateparams_gradient);
			return dLdp;
		}
		@Override 
		public String toString()
		{
			return "p"+node+"="+get();
		}
	}
	
	private class DuplicationParameter implements ModelParameter
	{
		DuplicationParameter(int node)
		{
			this.node = node;
			double q = base.getDuplicationParameter(node);
//			System.out.println("#**MLD.DP "+node+"\tq "+q+"\trate "+rates.getDuplicationRate(node)+"\tlen "+rates.getEdgeLength(node));
			node_distribution_parameters[node][PARAMETER_DUPLICATION]=q;
		}
		
		private final int node;
		
		@Override
		public double get()
		{
//			double q = rates.getDuplicationParameter(node);
//			node_parameters[node][PARAMETER_DUPLICATION]=q;
			double q=node_distribution_parameters[node][PARAMETER_DUPLICATION];
			return q;
		}
		
		@Override
		public void set(double q)
		{
			if (q==1.0)
				throw new IllegalArgumentException();
			if (!Double.isFinite(q))
			{
				System.out.println("#**MLMR.DP.set q "+q+"\t"+base.toString(node)+"\tnpars "+Arrays.toString(node_distribution_parameters[node]));
			}
			assert (Double.isFinite(q));
			assert (0.0<=q);
			assert (q<1.0);
			node_distribution_parameters[node][PARAMETER_DUPLICATION]=q;
//			rates.setParameters(node, node_parameters[node][PARAMETER_GAIN], node_parameters[node][PARAMETER_LOSS], q);
		}
		
		@Override 
		public double dL(double[] rateparams_gradient)
		{
			double dLdq = MLMixedRate.this.gradient.inferDistributionGradient(node, PARAMETER_DUPLICATION, rateparams_gradient);
			
			return dLdq; 	
		}
		@Override 
		public String toString()
		{
			return "q"+node+"="+get();
		}
	}
	
	/**
	 * Sets the node-specific distribution parameters.
	 */
	private void copyNodeParametersToModel()
	{
		for (int node=0; node<node_distribution_parameters.length; node++)
		{
			if (node_distribution_parameters[node]!=null)
			{
				double κ = node_distribution_parameters[node][PARAMETER_GAIN];
				double p = node_distribution_parameters[node][PARAMETER_LOSS];
				double q = node_distribution_parameters[node][PARAMETER_DUPLICATION];
				base.setParameters(node, κ, p, q);	
				//System.out.println("#**MLMR.cNPTM "+node+"\t"+Arrays.toString(node_distribution_parameters[node])+"\t"+base.toString(node));
				
			}
		}
	}
	
	
	public double[] getParameterValues()
	{
		double[] x = new double[rate_params.size()];
		for (int j=0; j<rate_params.size(); j++)
		{
			x[j] = rate_params.get(j).get();
		}
		return x;
	}
	
	public void setParameterValues(double[] x)
	{
		
		for (int j=0; j<rate_params.size(); j++)
		{
			ModelParameter P=rate_params.get(j);
			double old_x = P.get();
			P.set(x[j]);
//			System.out.println("#**MLMR.sPV "+j+"\t"+P+"\tx "+x[j]+"\txold "+old_x);
		}
		copyNodeParametersToModel();
		gradient.computeClasses();
//		System.out.println(count.io.RateVariationParser.printRates(base));
	}	
	
	private int calls_optFunc=0;
	
	/**
	 * Target function for minimization: negative log-likelihood.
	 * @param x
	 * @return
	 */
	public double optRatesFunc(double[] x)
	{
		setParameterValues(x);
		double LL = gradient.getCorrectedLL();
		
		calls_optFunc++;
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#*MLM.oRF "+calls_optFunc+"\t"+LL+"\tL0 "+gradient.getUnobservedLL());
		
		return -LL;
	}
	
	
	
	private int calls_optDiff = 0;
	
	/**
	/**
	 * Gradient for target function {@link #optRatesFunc(double[])}.
	 * @param x
	 * @return
	 */
	public double[] optRatesDiff(double [] x)
	{
		setParameterValues(x);
		calls_optDiff++; 

		
		double[] dL = gradient.getCorrectedGradient();
		double glen = FunctionMinimization.euclideanNorm(dL);
		
//		if (glen > 1e6)
//		{
//			gradient.reportClassParameters();
//			gradient.reportRatesGradient(dL);
//			
//
//		}
		double[] D = new double[rate_params.size()];
		for (int j=0; j<D.length; j++)
			D[j] = -rate_params.get(j).dL(dL);
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#*MLMR.oRD "+calls_optDiff+"\t"+gradient.getCorrectedLL()+"\tL0 "+gradient.getUnobservedLL()+"\tgrad "+glen+"\tdL "+FunctionMinimization.euclideanNorm(D));
		
		
		return D;
	}
	
	private static final double RATES_DIFF_EPS = 1e-7; // cubic root of machine precision
	
	protected double[] estimateRatesDiff(double [] x)
	{
		double fx = optRatesFunc(x);
		
		double[] D = new double[x.length];
		for (int p=0; p<x.length; p++)
		{
			double θ = x[p];
			double h = Math.abs(θ*RATES_DIFF_EPS);
			x[p] = θ+h;
			double fd = optRatesFunc(x);
			double delta = fd-fx;
			double dfdθ = delta/h;
			D[p] = dfdθ;
			
			x[p] = θ;
		}
		setParameterValues(x);
		return D;
	}
	
	
	protected void reportDistributionGradient()
	{
		double[] dL = gradient.getCorrectedGradient();
		int num_nodes = base.getTree().getNumNodes();
		assert(dL.length == 4*num_nodes);
		
		for (int node=0; node<num_nodes; node++)
		{
			double dLdp  = gradient.inferDistributionGradient(node, PARAMETER_LOSS, dL);
			double dLdq = gradient.inferDistributionGradient(node, PARAMETER_DUPLICATION, dL);
			double dLdr = gradient.inferDistributionGradient(node, PARAMETER_GAIN, dL);
			if (PRINT_OPTIMIZATION_MESSAGES)
				System.out.println("#**MLMR.rDG node "+node
						+"\tdLdr "+dLdr+"\t("+dL[4*node+PARAMETER_GAIN]+")"
						+"\tdLdp "+dLdp+"\t("+dL[4*node+PARAMETER_LOSS]+","+dL[4*node+PARAMETER_LENGTH]+")"
						+"\tdLdq "+dLdq+"\t("+dL[4*node+PARAMETER_DUPLICATION]+")"
					+"\t// "+base.toString(node));
		}
		
		
	}
	

	private class EdgeParameters implements Function<double[], Double>
	{
		private final int node;
		private final List<ModelParameter> edge_parameters;
		private final ModelParameter gain_parameter;
		private final ModelParameter loss_parameter;
		private final ModelParameter duplication_parameter;
		

		EdgeParameters(int node)
		{
			this.node = node;
			double[] params;
			if (node_distribution_parameters[node]==null)
			{
				params = node_distribution_parameters[node] = new double[3];
			} else
			{
				params = node_distribution_parameters[node];
			}
			params[PARAMETER_GAIN] = base.getGainParameter(node);
			params[PARAMETER_LOSS] = base.getLossParameter(node);
			params[PARAMETER_DUPLICATION] = base.getDuplicationParameter(node);
			
			
			this.edge_parameters = new ArrayList<>(3);
			boolean optimize_gain = optimize_rates[4*node+PARAMETER_GAIN];
			boolean optimize_loss = optimize_rates[4*node+PARAMETER_LOSS] || optimize_rates[4*node+PARAMETER_LENGTH];
			boolean optimize_duplication = optimize_rates[4*node+PARAMETER_DUPLICATION];
			
			if (optimize_gain)
			{
				gain_parameter = new GainParameter(node);
				edge_parameters.add(new Logarithmic(gain_parameter)); //newLogistic(new GainParameter(node), MAX_RATE)); 
			} else gain_parameter = null;
			if (optimize_loss)
			{
				loss_parameter = new LossParameter(node);
				edge_parameters.add(newLogistic(loss_parameter, MAX_PROB_NOT1));
			} else loss_parameter = null;
			if (optimize_duplication)
			{
				duplication_parameter = new DuplicationParameter(node);
				edge_parameters.add(newLogistic(duplication_parameter, MAX_PROB_NOT1));
			} else duplication_parameter = null;
		}
		
		/**
		 * 
		 * Negative log-likelihood. 
		 */
		@Override
		public Double apply(double[] x)
		{
			set(x);
			double LL = gradient.getCorrectedLL();
			calls_optFunc++;
			System.out.println("#*MLM.EP.a "+calls_optFunc+"\t"+LL+"\tL0 "+gradient.getUnobservedLL());
			
			return -LL;
		}
		
		public double[] get()
		{
			double[] x = new double[edge_parameters.size()];
			for (int i=0; i<x.length; i++)
			{
				x[i] = edge_parameters.get(i).get();
			}
			return x;
		}
		
		public void set(double[] x)
		{
			assert (x.length == edge_parameters.size());
			for (int i=0; i<x.length; i++)
			{
				edge_parameters.get(i).set(x[i]);
			}
			this.copyParametersToModel();
//			double κ = node_distribution_parameters[node][PARAMETER_GAIN];
//			double p = node_distribution_parameters[node][PARAMETER_LOSS];
//			double q = node_distribution_parameters[node][PARAMETER_DUPLICATION];
//			base.setParameters(node, κ, p, q);	
//			gradient.computeClasses();
		}
		
		private void copyParametersToModel()
		{
			double κ = node_distribution_parameters[node][PARAMETER_GAIN];
			double p = node_distribution_parameters[node][PARAMETER_LOSS];
			double q = node_distribution_parameters[node][PARAMETER_DUPLICATION];
			base.setParameters(node, κ, p, q);	
			gradient.computeClasses();
		}
		
		private ModelParameter getParameter(int param_type)
		{
			if (param_type == PARAMETER_GAIN)
			{
				return gain_parameter;
			} else if (param_type == PARAMETER_LOSS)
			{
				return loss_parameter;
			} else if (param_type == PARAMETER_DUPLICATION)
			{
				return duplication_parameter;
			} else
				throw new IllegalArgumentException(); // it shouldn't come to this
		}

		
		class PartialDerivative implements DoubleFunction<Double>
		{
			private final ModelParameter P;
			PartialDerivative(int param_type)
			{
				this.P = EdgeParameters.this.getParameter(param_type);
				// P should not be null
				// assert (P!=null);
			}
			
			
			/**
			 * Partial derivative of the negative log-likelihood wrt to this edge parameter.
			 * Updates the underlying model. 
			 * 
			 * @param θ parameter value; underlying model is set 
			 * @return derivative at point θ 
			 */
			@Override
			public Double apply(double θ) 
			{
				set(θ);
				
				double[] dL = gradient.getCorrectedGradient();
				double dLdθ = P.dL(dL);
				
				return -dLdθ;
			}
			
			public double get()
			{
				return P.get();
			}
			
			public void set(double θ)
			{
				P.set(θ);
				EdgeParameters.this.copyParametersToModel();
			}
			
			/**
			 * Negative log-likelihood, to be minimized as a function of this single parameter 
			 * 
			 * @return
			 */
			public DoubleFunction<Double> optParameterFunc()
			{
				return new DoubleFunction<Double>() {
					@Override 
					public Double apply(double θ)
					{
						set(θ);
						
						double LL = gradient.getCorrectedLL();
						
						calls_optFunc++;
						if (PRINT_OPTIMIZATION_MESSAGES)
							System.out.println("#*MLMR.EP.PD "+calls_optFunc+"\t"+LL+"\tL0 "+gradient.getUnobservedLL());
						
						return -LL;
					}
				};
			}
		}
	}
	
	
	
	/**
	 * Slow method optimizing one parameter at a time
	 * @param delta
	 * @param itmax
	 * @return
	 */
	protected double edgewiseOptimize(double delta, int itmax)
	{
		double LL = gradient.getCorrectedLL();
		double min = -LL;
		
		List<Double> optimization_history = getOptimizationHistory();
		if (!uniform_duplication && !uniform_gain && !uniform_loss)
		{
			int num_nodes = base.getTree().getNumNodes();
			
			int iter = 0;
			double interval_threshold = FunctionMinimization.POWELL_TOL;

			while (iter<itmax) 
			{
				double delta_iter = 0.0; // will be negative
				for (int node=0; node<num_nodes; node++)
				{
					double emin = min;
					EdgeParameters EP = new EdgeParameters(node);
	
					
					double delta_edge = 0.0; // will be negative
					EdgeParameters.PartialDerivative D = EP.new PartialDerivative(PARAMETER_LOSS);
				
					if (D.P!=null )
					{
						assert (Double.isFinite(base.getEdgeLength(node)));
						
						
						double x0 = D.get(); 
						DoubleFunction<Double> func = D.optParameterFunc();
						
						double x1 =1.0-MAX_PROB_NOT1 ; // Double.max(base.getDuplicationParameter(node),  );
						
						double x2 = MAX_PROB_NOT1;
						
//						double f2 = func.apply(x2);
//						double f1 = func.apply(x1);
//						double f0 = func.apply(x0);
						
						
						
						if (x0<=x1 || x2<=x0) // does not work with mnbrak
						{
							double[] braket = FunctionMinimization.mnbrak(x1, x2, func);
							// a,b, c, func(a), func(b), func(c)
							
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLMR.eO bracket loss "+Arrays.toString(braket));
							x1 = braket[0];
							x0 = braket[1];
							x2 = braket[2];
						}
						
						double[] xy = 
									FunctionMinimization.dbrent(x1,x0,x2,func, D, interval_threshold);
									// FunctionMinimization.brent(x1, x0, x2, func, interval_threshold);
									// FunctionMinimization.golden(x1, x0, x2, func, delta);
						// xy = {xm,ym} for minimum
						double xm = xy[0];
						emin = xy[1];
						
						if (optimization_history!=null) optimization_history.add(emin);
					
						
						
						double dLdxm = D.apply(xm); // so that the parameters are copied into model
						double LLxm = gradient.getCorrectedLL();
						double deltaLL = LLxm-LL;
						if (PRINT_OPTIMIZATION_MESSAGES)
							System.out.println("#**MLMR.eO node "+node+"\t"+LLxm+"\t("+deltaLL+")\tp "+xm+"\t("+x0+")\tdL "+dLdxm+"\t// mn "+Arrays.toString(xy));
								// +"\tf1 "+f1+"\tf2 "+f2+"\tf0 "+f0);
					
						LL = LLxm;
						delta_edge += deltaLL;
					}
					D = EP.new PartialDerivative(PARAMETER_DUPLICATION);
					if (D.P!=null)
					{
						double x0 = D.get(); 
						double x1 = 1.0-MAX_PROB_NOT1;
						double x2 = MAX_PROB_NOT1;
						if (Double.isInfinite(base.getEdgeLength(node)))
								x2 = base.getLossParameter(node);
						DoubleFunction<Double> func = D.optParameterFunc();
//						double f2 = func.apply(x2);
//						double f1 = func.apply(x1);
//						double f0 = func.apply(x0);
						
						if (x0<=x1 || x2<=x0) // does not work with mnbrak
						{
							double[] braket = FunctionMinimization.mnbrak(x1, x2, func);
							// a,b, c, func(a), func(b), func(c)
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLMR.eO bracket dup "+Arrays.toString(braket));
							x1 = braket[0];
							x0 = braket[1];
							x2 = braket[2];
						}
						
						double[] xy = FunctionMinimization.dbrent(x1,x0,x2,func, D, interval_threshold);
							// FunctionMinimization.brent(x1, x0, x2, func, interval_threshold);
							// FunctionMinimization.golden(x1, x0, x2, func, delta);
							// FunctionMinimization.mnbrak(x1, x2, D.optParameterFunc());
						double xm = xy[0];
						emin = xy[1];
						if (optimization_history!=null) optimization_history.add(emin);
					
						double dLdxm = D.apply(xm); // so that the parameters are copied into model
						double LLxm = gradient.getCorrectedLL();
						double deltaLL = LLxm-LL;
						if (PRINT_OPTIMIZATION_MESSAGES)
							System.out.println("#**MLMR.eO node "+node+"\t"+LLxm+"\t("+deltaLL+")\tq "+xm+"\t("+x0+")\tdL "+dLdxm+"\t// x12 "+x1+", "+x2+"\t; mn "+Arrays.toString(xy));
									//+"\tf1 "+f1+"\tf2 "+f2+"\tf0 "+f0);

						LL = LLxm;
						delta_edge += deltaLL;
					}
					D = EP.new PartialDerivative(PARAMETER_GAIN);
					if (D.P!=null)
					{
						double x0 = D.get(); 
						double x1 = 1.0-MAX_PROB_NOT1;
						double x2 = MAX_RATE;
						
						DoubleFunction<Double> func = D.optParameterFunc();
//						double f2 = func.apply(x2);
//						double f1 = func.apply(x1);
//						double f0 = func.apply(x0);
						
						if (x0<=x1 || x2<=x0) // does not work with mnbrak
						{
							double[] braket = FunctionMinimization.mnbrak(x1, x2, func);
							// a,b, c, func(a), func(b), func(c)
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLMR.eO bracket gain "+Arrays.toString(braket));
							x1 = braket[0];
							x0 = braket[1];
							x2 = braket[2];
						}
						
	
						double[] xy = FunctionMinimization.dbrent(x1,x0,x2,func, D, interval_threshold);
								// FunctionMinimization.brent(x1, x0, x2, func, interval_threshold);
								// FunctionMinimization.golden(x1, x0, x2, func, delta);
						double xm = xy[0];
						emin = xy[1];
						if (optimization_history!=null) optimization_history.add(emin);
	
						double dLdxm = D.apply(xm); // so that the parameters are copied into model
						double LLxm = gradient.getCorrectedLL();
						double deltaLL = LLxm-LL;
						if (PRINT_OPTIMIZATION_MESSAGES)
							System.out.println("#**MLMR.eO node "+node+"\t"+LLxm+"\t("+deltaLL+")\tr "+xm+"\t("+x0+")\tdL "+dLdxm+"\t// mn "+Arrays.toString(xy));
							//+"\tf1 "+f1+"\tf2 "+f2+"\tf0 "+f0);
						
						LL = LLxm;
						delta_edge += deltaLL;
						
					}
					delta_iter += delta_edge;
					System.out.println("#**MLMR.eO node "+node+"\tdelta "+delta_edge);

					System.out.println("#**MLM.eO node "+node+"\temin "+emin+"//\tmin "+min);
					assert (emin<=min);
					min=emin;
				} // for node
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#**MLM.eO iter "+iter+"\tdelta "+delta_iter);
				iter ++;
				if (delta_iter>-delta) break;
			} // while ()
		}
		return min;
	}
	
	@Override
	public double optimize(double delta)
	{
		return optimize(delta, 1000);
	}
	
	/** 
	 * Optimizes the base rates. 
	 */
	@Override
	public double optimize(double delta, int itmax)
	{
		if (PRINT_OPTIMIZATION_MESSAGES)
		{
			double LL =  gradient.getCorrectedLL();
			System.out.println("#*MLMR.o before LL "+LL);
		}
		initModelParameters();
		
		double LL =  gradient.getCorrectedLL();
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#**MLMR.o starting LL "+LL+"\tparams "+getModelParameterCount());
		
		gradient.computeClasses();
//		{
//			double cLL =  gradient.getCorrectedLL();
//			System.out.println("#*MLMR.o recompute LL "+cLL);
//		}
			//+"\t// "+rate_params.toString());
		
		double[] x0 = getParameterValues();
		setParameterValues(x0);
		if (PRINT_OPTIMIZATION_MESSAGES)
		{
			double cLL =  gradient.getCorrectedLL();
			System.out.println("#*MLMR.o reset LL "+cLL);
		}
		
		
		List<Double> optimization_history = getOptimizationHistory();
		
		
		if (itmax>0)
		{

			
//			{ // debug
//				double[] dLL = optRatesDiff(x0);
//				for (int pidx=0; pidx<rate_params.size(); pidx++)
//				{
//					System.out.println("#*MLMR.o param "+pidx+"\t"+rate_params.get(pidx)+"\tdLdx "+dLL[pidx]);
//				}
//			}
			double min;
			
//			min = FunctionMinimization.gradientDescent(x0, delta, itmax, x->optRatesFunc(x), x->optRatesDiff(x), optimization_history);
			
			
			
			if (x0.length==2 )
			{
				min = FunctionMinimization.powell(x0, delta, itmax, x->optRatesFunc(x), optimization_history);
			} else
			{
//				{
//					// init step
//					double[] dL = optRatesDiff(x0);
////					for (int pidx=0; pidx<dL.length; pidx++) dL[pidx]=-dL[pidx];
////					double initmin = FunctionMinimization.linmin(x0, dL, x->optRatesFunc(x));
//					System.out.println("#**MLMR.o init dL "+FunctionMinimization.euclideanNorm(dL));
//				}
				
				
				
				min  = FunctionMinimization.dfpmin(x0, delta, itmax, x->optRatesFunc(x), x->optRatesDiff(x), optimization_history);
						  // 
			}
	        if (Thread.currentThread().isInterrupted()) // keeps interrupt status
	        {
	            return Double.NaN;
	        }

	        setParameterValues(x0);
		}
		
		
		gradient.getCorrectedClassGradients();

        LL =  gradient.getCorrectedLL();
        
        
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#**MLM.oR final LL "+LL);
		
		
//		setParameterValues(x0);
//		LL =  gradient.getCorrectedLL();
//		System.out.println("#**MLM.oR by dfpmin "+LL);
		
		return -LL;
	}

	/**
	 * Default setting for rate optimization: not loss, but all other non-0  
	 * 
	 * @param rates
	 * @return
	 */
	private static boolean[] optimizableRates(TreeWithRates rates)
	{
		IndexedTree tree = rates.getTree();
		int num_nodes = tree.getNumNodes();
		
		boolean[] opt_par = new boolean[4*num_nodes];
		for (int node=0; node<num_nodes; node++)
		{
			if (tree.isRoot(node))
			{
				opt_par[4*node+PARAMETER_GAIN] = rates.getGainParameter(node)>0.0;
				opt_par[4*node+PARAMETER_LOSS] = false;
				opt_par[4*node+PARAMETER_DUPLICATION] = rates.getDuplicationParameter(node)>0.0;
				double len = rates.getEdgeLength(node);
				opt_par[4*node+PARAMETER_LENGTH] = Double.isFinite(len) && len>0.0;				
			} else
			{
				opt_par[4*node+PARAMETER_GAIN] = rates.getGainParameter(node)>0.0;
				opt_par[4*node+PARAMETER_LOSS]=false;
				opt_par[4*node+PARAMETER_DUPLICATION] = rates.getDuplicationParameter(node)>0.0;
				double len = rates.getEdgeLength(node);
				opt_par[4*node+PARAMETER_LENGTH] = Double.isFinite(len); // && len>0.0;				
			} 
		}
		
		return opt_par;
	}	
	
	private void debugGradient()
	{
		optimize(0,0); 
		double[] x = getParameterValues();
		
		double[] df = optRatesDiff(x);
		
		double[] df_est = estimateRatesDiff(x);
		
		for (int i=0; i<x.length; i++)
		{
			ModelParameter P = rate_params.get(i);
			
			double df_delta = df_est[i]-df[i];
			double rel_delta = Math.abs(df_delta/df[i]);
			
			System.out.println("#**MLMR.dG param "+i+"\t"+P+"\tdfdx "+df[i]+"\test "+df_est[i]
					+"\tdiff "+df_delta+"\t("+rel_delta+")");
		}
		
		
	}
	
	private void debugOptimization()
	{
		optimize(1e-5);
	}
		
	public static void main(String[] args) throws Exception
	{
		PRINT_OPTIMIZATION_MESSAGES = true;

		count.io.CommandLine cli = new count.io.CommandLine(args, MLMixedRate.class);
		MLMixedRate O = new MLMixedRate(cli.getModel(), cli.getTable());
		
		O.debugGradient();

		//System.out.println(count.io.RateVariationParser.printRates(cli.getModel()));

	}
}
