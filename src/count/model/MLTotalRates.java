package count.model;

import java.util.ArrayList;
import java.util.List;


import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LOSS;

import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.matek.FunctionMinimization;
import count.model.ML.ModelParameter;

/** 
 * Optimization of rate parameters.
 * 
 * @author csuros
 * @deprecated
 */
public class MLTotalRates extends ML
{
	private static final boolean OPTIMIZE_COMBINED_GAIN_RATE = true;
	public MLTotalRates(TreeWithRates rates, ProfileTable table)
	{
		this(rates, table, optimizableParams(rates));
	}
	
	private static final double MAX_RATE = 8.0;
	
	public MLTotalRates(TreeWithRates rates, ProfileTable table, boolean[] optimize_parameter)
	{
		this.rates = rates;
		this.table = table;
		assert (optimize_parameter.length == 3*rates.getTree().getNumNodes());
		
		this.optimize_parameter = optimize_parameter;
		this.params = new ArrayList<>();
		this.gradient = new Gradient(rates, table);
		System.out.println("#*MLT.init "+gradient.factory.getCorrectedLL());
	}
	
	private final TreeWithRates rates;
	private final ProfileTable table;
	private boolean[] optimize_parameter;
	private final List<ModelParameter> params;
	private final Gradient gradient;
	
	
	private List<Double> optimization_history=null;
	/**
	 * Sets the tracking for iterations during optimization. 
	 * 
	 * @param history will contain the successive function values calculated during optimization
	 */
	public void setOptimizationHistory(List<Double> history)
	{
		this.optimization_history = history;
	}
	
	@Override
	public void setCalculationWidth(int absolute, double relative)	
	{
		gradient.setCalculationWidthThresholds(absolute, relative);
	}
	
	@Override 
	public void setMinimumObservedCopies(int m)
	{
		gradient.setMinimumObservedCopies(m);
	}
	
	
	private void initModelParameters()
	{
		params.clear();
		
		final IndexedTree phylo = rates.getTree();
		final int num_nodes = phylo.getNumNodes();
		
		for (int node=0; node<num_nodes; node++)
		{
			if (optimize_parameter[3*node + PARAMETER_GAIN])
			{
				ModelParameter θ = new GainRate(node); 
				addLogistic(θ);
				//System.out.println("#**MLTR.iMP "+θ);
			}
			if (optimize_parameter[3*node+PARAMETER_LOSS])
			{
				ModelParameter θ = new LossRate(node); 
				addLogistic(θ);
				//System.out.println("#**MLTR.iMP "+θ);
			}
			if (optimize_parameter[3*node + PARAMETER_DUPLICATION])
			{
				ModelParameter θ = new DuplicationRate(node); 
				
				if (phylo.isRoot(node) || rates.getEdgeLength(node)==Double.POSITIVE_INFINITY)
				{
					addLogistic(θ,(1.-1e-9)*rates.getLossRate(node));
				} else
				{
					addLogistic(θ);
				}
				//System.out.println("#**MLTR.iMP "+θ);
			}
		}
		
//		
//		int i = 0;
//		int j = 0;
//		while (i<optimize_parameter.length)
//		{
//			if (optimize_parameter[i]){ j++;}
//			i++;
//		}
//		this.params = new ModelParameter[j];
//		while (j>0)
//		{
//			i--;
//			if (optimize_parameter[i])
//			{
//				j--;
//				int node = i/3;
//				int type = i%3;
//				ModelParameter theta=null;
//				if (type == PARAMETER_GAIN)
//					theta = new GainRate(node);
//				else if (type == PARAMETER_LOSS)
//					theta = new LossRate(node);
//				else if (type == PARAMETER_DUPLICATION)
//					theta = new DuplicationRate(node);
//				
//				double θ = theta.get();
//				if (θ>=MAX_RATE)
//				{
//					System.out.println("#*MLT.init "+theta+"; resetting. //"
//							+ " len "+rates.getEdgeLength(node)
//							+ "\tdup "+rates.getDuplicationRate(node));
//					theta.set((1.-1e-9)*MAX_RATE);
//				}
//				params[j] = new Logistic(theta, MAX_RATE); // new Logarithmic(theta); // new LogLogistic(theta);
////				params[j] = new Logarithmic(theta);
//			}
//		}
//		
	}

	private void addLogistic(ModelParameter theta)
	{
		addLogistic(theta, MAX_RATE);
	}
	private void addLogistic(ModelParameter theta, double max_rate)
	{
		double θ = theta.get();
		if (θ>=max_rate)
		{
			double newθ = (1.-1e-9)*max_rate;
			System.out.println("#*MLTR.init "+theta+"; resetting to "+newθ);
			theta.set(newθ);
		}
		params.add(new Logistic(theta, max_rate)); 
		// new Logarithmic(theta); // new LogLogistic(theta);
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
			return rates.getGainRate(node);
		}
		@Override
		public void set(double x)
		{
			rates.setGainRate(node, x);
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dLdκ = rate_gradient[3*node+PARAMETER_GAIN];
			return dLdκ ;
		}
		@Override 
		public String toString()
		{
			return "g"+node+"="+get();
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
			double μ = rates.getLossRate(node); 
			double t = rates.getEdgeLength(node);
			if (Double.isInfinite(t))
				return μ;
			else
				return μ*t;
		}
		@Override
		public void set(double x)
		{
//			assert (!Double.isNaN(x));
			double μ = rates.getLossRate(node); 
			double λ = rates.getDuplicationRate(node);
			double t = rates.getEdgeLength(node);
			if (Double.isInfinite(t))
				rates.setLossRate(node, x);
			else
			{ // keep the same loss rate and total duplication; change edge length
				double tt = x/μ;
				rates.setEdgeLength(node, tt);
				// x/mu * lm = lambda*t
				// lm = lambda*t*mu/x
				rates.setDuplicationRate(node, λ*t/tt);
			}
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dLdμ = rate_gradient[3*node+PARAMETER_LOSS];
			return dLdμ ;
		}
		@Override 
		public String toString()
		{
			return "l"+node+"="+get();
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
			double λ = rates.getDuplicationRate(node);
			double t = rates.getEdgeLength(node);
			if (Double.isInfinite(t))
				return λ;
			else
				return λ*t;
		}
		@Override
		public void set(double x)
		{
			double t = rates.getEdgeLength(node);
			if (Double.isInfinite(t))
			{
				rates.setDuplicationRate(node,x);
			} else
			{
				rates.setDuplicationRate(node, x/t);
			}
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dLdλ = rate_gradient[3*node+PARAMETER_DUPLICATION];
			return dLdλ ;
		}
		@Override 
		public String toString()
		{
			return "d"+node+"="+get();
		}
	}			

	public double[] get()
	{
		double[] x = new double[params.size()];
		for (int j=0; j<params.size(); j++)
		{
			x[j] = params.get(j).get();
		}
		return x;
	}
	
	public void set(double[] x)
	{
		for (int j=0; j<params.size(); j++)
			params.get(j).set(x[j]);
		gradient.computeParameters();
//		gradient.factory.printParameters(System.out);
//		System.out.println(count.io.RateVariationParser.printRates(rates));
	}
	
	private int calls_optDiff=0;
	private int calls_optFunc=0;
	
	public double optFunc(double[] x)
	{
//		System.out.println("#*MLR.oF "+java.util.Arrays.toString(x));
		this.set(x);
		double LL = gradient.getCorrectedLL();
		calls_optFunc++;
		//System.out.println("#*MLT.oF "+calls_optFunc+"\t"+LL+"\tL0 "+gradient.getUnobservedLL());
		return -LL;
	}
	
	public double[] optDiff(double[] x)
	{
		set(x);
		double LL = gradient.getCorrectedLL();
		calls_optDiff++;
		System.out.println("#*MLT.oD "+calls_optDiff+"\t"+LL);

		double[] survival_gradient = gradient.getCorrectedGradient();
		double[] distribution_gradient = gradient.getDistributionGradient(survival_gradient);
		double[] rate_gradient = gradient.getTotalRateGradient(distribution_gradient);
		
		double[] D = new double[params.size()];
		for (int j=0; j<D.length; j++)
		{
			D[j] = -params.get(j).dL(rate_gradient);
		}
		return D;
	}
	
	private static final double RATES_DIFF_EPS = 1e-7; // cubic root of machine precision
	protected double[] estDiff(double[] x)
	{
		double fx = optFunc(x);
		
		double[] D = new double[x.length];
		for (int p=0; p<x.length; p++)
		{
			double θ = x[p];
			double h = Math.abs(θ*RATES_DIFF_EPS);
			x[p] = θ+h;
			double fd = optFunc(x);
			double delta = fd-fx;
			double dfdθ = delta/h;
			D[p] = dfdθ;
			
			x[p] = θ;
		}
		set(x);
		return D;
	}
	
	@Override
	public double optimize(double delta)
	{
		return optimize(delta, 10000);
	}
	
	@Override
	public double optimize(double delta, int max_iter)
	{
		double LL =  gradient.getCorrectedLL();

		this.initModelParameters();
		System.out.println("#**MLTR.o starting LL "+LL+"\tparams "+params.size());
		
		double[] x0 = get();
		
		if (max_iter > 0)
		{
			double min  = FunctionMinimization.dfpmin(x0, delta, max_iter, x->optFunc(x), x->optDiff(x), optimization_history);
		
			LL =  gradient.getCorrectedLL();
			System.out.println("#**MLTR.o final LL "+LL+"\t(min "+min+")"+"\tfcalls "+calls_optFunc+"\tdcalls "+calls_optDiff);
//		set(x0);
		}
//		LL =  gradient.getCorrectedLL();
//		System.out.println("#**MLT.o by dfpmin "+LL);
		
		return -LL;
	}
	
	private static boolean[] optimizableParams(TreeWithRates rates)
	{
		IndexedTree tree = rates.getTree();
		int num_nodes = tree.getNumNodes();
		
		boolean[] opt_par = new boolean[3*num_nodes];
		for (int node=0; node<num_nodes; node++)
		{
			if (tree.isRoot(node))
			{
				opt_par[3*node + PARAMETER_GAIN] = rates.getGainParameter(node)>0.0;
				opt_par[3*node+PARAMETER_LOSS] = !Double.isInfinite(rates.getEdgeLength(node));
				opt_par[3*node+PARAMETER_DUPLICATION] = rates.getDuplicationParameter(node)>0.0;
			} else
			{
				opt_par[3*node+PARAMETER_GAIN] = rates.getGainParameter(node)>0.0;
				opt_par[3*node+PARAMETER_LOSS]=rates.getEdgeLength(node)>0;
				opt_par[3*node+PARAMETER_DUPLICATION] = rates.getDuplicationParameter(node)>0.0;
			} 
		}
		
		return opt_par;
	}	
		
	private void debugGradient()
	{
		optimize(0,0); 
		double[] x = get();
		
		double[] df = optDiff(x);
		
		double[] df_est = estDiff(x);
		
		for (int i=0; i<x.length; i++)
		{
			ModelParameter P = params.get(i);
			
			double df_delta = df_est[i]-df[i];
			double rel_delta = Math.abs(df_delta/df[i]);
			
			System.out.println("#**MLTR.dG param "+i+"\t"+P+"\tdfdx "+df[i]+"\test "+df_est[i]
					+"\tdiff "+df_delta+"\t("+rel_delta+")");
		}
	}
	
	private void debugOptimization()
	{
		optimize(1e-9);
	}
	
	public static void main(String[] args) throws Exception
	{
		count.io.CommandLine cli = new count.io.CommandLine(args, MLRates.class);
    	GammaInvariant model = null; 
    	TreeWithRates base_rates;
    	if (cli.getModel()==null)
    	{
    		model = GammaInvariant.nullModel(cli.getTree());
    		base_rates = model.getBaseModel();
    		
    	} else
    	{
    		base_rates = cli.getRates(); 
    		model = cli.getModel();
    	}
		
		MLTotalRates O = new MLTotalRates(base_rates, cli.getTable());
		
		
		O.debugGradient();
		
		//System.out.println(count.io.RateVariationParser.printRates(model));
	}
	
}
