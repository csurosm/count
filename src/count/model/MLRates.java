package count.model;

import static count.io.CommandLine.OPT_EPS;
import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_OUTPUT;
import static count.io.CommandLine.OPT_ROUNDS;
import static count.io.CommandLine.OPT_TRUNCATE;
import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LENGTH;
import static count.model.GLDParameters.PARAMETER_LOSS;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.io.CommandLine;
import count.io.NewickParser;
import count.matek.FunctionMinimization;

/**
 * 
 * 
 * @author csuros
 *
 */
public class MLRates extends ML
{
	public static boolean PRINT_OPTIMIZATION_MESSAGES = false;

	private static final boolean OPTIMIZE_COMBINED_GAIN_RATE = false; // true is untested
	private boolean do_gradient_descent = true; // if true, start with a few rounds of gradient descent (slower) 
	private boolean optimization_by_quasiNewton = true; // if false, use conjugate gradient (slower)
	private static final boolean CALCULATE_FROM_TOTAL_RATES = false;
	private static final boolean TRACK_RATE_GAP = false; // true is too much precision, leads to numerical errors
	
	
	public MLRates(TreeWithRates rates, ProfileTable table)
	{
		this(rates, table, optimizableParams(rates));
	}

	public MLRates(TreeWithRates rates, ProfileTable table, boolean[] optimize_parameter)
	{
		this.rates = rates;
		this.table = table;
		assert (optimize_parameter.length == 4*rates.getTree().getNumNodes());
		
		this.params = new ArrayList<>();
		this.optimize_parameter = optimize_parameter;
//		
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
//				int node = i/4;
//				int type = i%4;
//				ModelParameter theta=null;
//				if (type == PARAMETER_GAIN)
//					theta = new GainRate(node);
//				else if (type == PARAMETER_LOSS)
//					theta = new LossRate(node);
//				else if (type == PARAMETER_DUPLICATION)
//					theta = new DuplicationRate(node);
//				else if (type == PARAMETER_LENGTH)
//					theta = new EdgeLength(node);
//				params[j] = new Logistic(theta); // new Logarithmic(theta); // new LogLogistic(theta);
//				//if (type == PARAMETER_DUPLICATION) params[j] = new Logarithmic(theta,0.001);
//			}
//		}
		this.gradient = new Gradient(rates, table);
		//System.out.println("#*MLR.init "+gradient.factory.getCorrectedLL());
		this.likelihood_scaling = 1.0; //rates.getTree().getNumLeaves()*table.getFamilyCount()*1e-3; 
	}
	
	private final TreeWithRates rates;
	private final ProfileTable table;
	private final List<ModelParameter> params;
	private final boolean[] optimize_parameter;
//	private List<Double> optimization_history=null;
	private final Gradient gradient;
	private double likelihood_scaling;

//	public void setOptimizationHistory(List<Double> history)
//	{
//		this.optimization_history = history;
//	}
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
	
//	private void addLogistic(ModelParameter theta)
//	{
//		addLogistic(theta, MAX_RATE);
//	}
//	private void addLogistic(ModelParameter theta, double max_rate)
//	{
//		double θ = theta.get();
//		if (θ>=max_rate)
//		{
//			double newθ = (1.-1e-9)*max_rate;
//			System.out.println("#*MLR.init "+theta+"; resetting to "+newθ);
//			theta.set(newθ);
//		}
//		params.add(new Logistic(theta, max_rate)); 
//		// new Logarithmic(theta); // new LogLogistic(theta);
//	}
//	
//	private void addLogarithmic(ModelParameter theta)
//	{
//		params.add(new Logarithmic(theta));
//	}
	public boolean fixGain(int node, boolean not_optimized)
	{
		boolean b = optimize_parameter[4*node+PARAMETER_GAIN];
		optimize_parameter[4*node+PARAMETER_GAIN]=!not_optimized;
		return b;
	}
	public boolean fixLoss(int node, boolean not_optimized)
	{
		boolean b = optimize_parameter[4*node+PARAMETER_LOSS];
		optimize_parameter[4*node+PARAMETER_LOSS]=!not_optimized;
		return !b;
	}
	public boolean fixDuplication(int node, boolean not_optimized)
	{
		boolean b = optimize_parameter[4*node+PARAMETER_DUPLICATION];
		optimize_parameter[4*node+PARAMETER_DUPLICATION]=!not_optimized;
		return b;
	}
	public boolean fixLength(int node, boolean not_optimized)
	{
		boolean b = optimize_parameter[4*node+PARAMETER_LENGTH];
		optimize_parameter[4*node+PARAMETER_LENGTH]=!not_optimized;
		return !b;
	}
	
	@Override
	public void fixNodeParameters(int node, boolean do_not_optimize)
	{
		fixGain(node, do_not_optimize);
		fixLoss(node, do_not_optimize);
		fixDuplication(node, do_not_optimize);
		fixLength(node, do_not_optimize);
	}
	
	private void initModelParameters()
	{
		params.clear();
		
		final IndexedTree phylo = rates.getTree();
		final int num_nodes = phylo.getNumNodes();
		
		for (int node=0; node<num_nodes; node++)
		{
			if (optimize_parameter[4*node + PARAMETER_GAIN])
			{
				ModelParameter θ = new GainRate(node); 
				//addLogistic(θ);
				//addLogarithmic(θ);
//				params.add(new Logarithmic(θ));
				params.add(newLogistic(θ, MAX_RATE));
//				System.out.println("#**MLR.iMP "+θ);
			}
			if (optimize_parameter[4*node+PARAMETER_LOSS])
			{
				ModelParameter θ = new LossRate(node); 
				//addLogistic(θ);
				params.add(newLogistic(θ, MAX_RATE));
//				System.out.println("#**MLR.iMP "+θ);
			}
			if (optimize_parameter[4*node + PARAMETER_DUPLICATION])
			{
				BoundedParameter θ = new DuplicationRate(node); 
				
//				if (phylo.isRoot(node) || rates.getEdgeLength(node)==Double.POSITIVE_INFINITY)
//				{
//					params.add(newLogistic(θ,MAX_PROB_NOT1*rates.getLossRate(node)));
				if (TRACK_RATE_GAP)
					params.add(new BoundedLogistic(θ));
				else
					params.add(newLogistic(θ,1.0*rates.getLossRate(node)));
//				} else
//				{
//					params.add(newLogistic(θ));
//				}
//				System.out.println("#**MLR.iMP "+θ);
			}
			if (optimize_parameter[4*node + PARAMETER_LENGTH])
			{
				ModelParameter θ = new EdgeLength(node); 
				//addLogistic(θ);
				params.add(newLogistic(θ, MAX_RATE));
//				System.out.println("#**MLR.iMP "+θ);
			}
		}
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
			double κ = rates.getGainRate(node);
			if (OPTIMIZE_COMBINED_GAIN_RATE)
			{
				double λ = rates.getDuplicationRate(node);
				return ( λ==0.0 ? κ : κ*λ );
			} else
				return κ;
		}
		@Override
		public void set(double x)
		{
			if (OPTIMIZE_COMBINED_GAIN_RATE)
			{
				double λ = rates.getDuplicationRate(node);
				if (λ == 0.0)
					rates.setGainRate(node, x);
				else
					rates.setGainRate(node, x/λ);
			} else
				rates.setGainRate(node, x);
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dLdκ = CALCULATE_FROM_TOTAL_RATES ? rate_gradient[3*node+PARAMETER_GAIN] : rate_gradient[4*node+PARAMETER_GAIN];
			if (OPTIMIZE_COMBINED_GAIN_RATE)
			{
				double λ = rates.getDuplicationRate(node);
				if (λ == 0.0)
				{
					return dLdκ;
				} else
				{
					double κ = rates.getGainRate(node);
					double dLdλ = rate_gradient[3*node+PARAMETER_DUPLICATION];
					return dLdκ*λ + dLdλ * κ;
				}
			}
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
			return rates.getLossRate(node);
		}
		@Override
		public void set(double x)
		{
			rates.setLossRate(node, x);
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dLdμ;
			if (CALCULATE_FROM_TOTAL_RATES)
			{
				dLdμ = rate_gradient[3*node+PARAMETER_LOSS];
				double t = rates.getEdgeLength(node);
				if (!rates.getTree().isRoot(node) || !Double.isInfinite(t))
					dLdμ *= t;
			} else
			{
				dLdμ = rate_gradient[4*node+PARAMETER_LOSS];
			}
			return dLdμ ;
				
		}
		@Override 
		public String toString()
		{
			return "l"+node+"="+get();
		}
	}		
	private class DuplicationRate implements BoundedParameter
	{
		private final int node;
		DuplicationRate(int node)
		{
			this.node = node;
		}
		@Override
		public double get()
		{
			return rates.getDuplicationRate(node);
		}
		@Override
		public double getComplement()
		{
			return rates.getRateGap(node);
		}
		@Override
		public void set(double x)
		{
			rates.setDuplicationRate(node,x);
			//System.out.println("#**MLR.DR.set "+x+"\t"+this);
		}
		@Override 
		public void set(double x, double xcomp)
		{
			rates.setDuplicationRate(node, x, xcomp);
			if (x==getBound())
			{
				System.out.println("#**MLR.DR.set2 "+toString());
			}
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dLdλ;
			if (CALCULATE_FROM_TOTAL_RATES)
			{
				dLdλ = rate_gradient[3*node+PARAMETER_DUPLICATION];
				double t = rates.getEdgeLength(node);
				if (!rates.getTree().isRoot(node) || !Double.isInfinite(t))
					dLdλ *= t;
			} else
			{
				dLdλ = rate_gradient[4*node+PARAMETER_DUPLICATION];
			}
			return dLdλ ;
		}
		@Override 
		public String toString()
		{
			return "d"+node+"="+get()+"/"+getBound()+"-"+getComplement();
		}
		
		@Override
		public double getBound()
		{
			return rates.getLossRate(node);
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
			return rates.getEdgeLength(node);
		}
		@Override
		public void set(double x)
		{
			rates.setEdgeLength(node, x);
		}
		@Override 
		public double dL(double[] rate_gradient)
		{
			double dLdt;
			if (CALCULATE_FROM_TOTAL_RATES)
			{
				double dLdμ = rate_gradient[3*node+PARAMETER_LOSS];
				double dLdλ = rate_gradient[3*node+PARAMETER_DUPLICATION];
				double t = rates.getEdgeLength(node);
				if (rates.getTree().isRoot(node) && Double.isInfinite(t))
					dLdt = 0.0;
				else
				{
					double μ = rates.getLossRate(node);
					double λ = rates.getDuplicationRate(node);
					
					dLdt =  μ * dLdμ + λ * dLdλ;
				}
			} else
			{
				dLdt = rate_gradient[4*node+PARAMETER_LENGTH];
			}
			return dLdt;
		}
		@Override 
		public String toString()
		{
			return "t"+node+"="+get();
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
	}
	
	@Override
	public int getModelParameterCount()
	{
		return params.size();
	}
	
	
	private int calls_optDiff=0;
	private int calls_optFunc=0;
	
	public double optFunc(double[] x)
	{
//		System.out.println("#*MLR.oF "+java.util.Arrays.toString(x));
		this.set(x);
		double LL = gradient.getCorrectedLL();
		calls_optFunc++;
		//System.out.println("#*MLR.oF "+calls_optFunc+"\t"+LL);
		//assert (!Double.isInfinite(LL));
		return -LL/likelihood_scaling;
	}
	
	
	private double[] parameterGradient(double[] rate_gradient)
	{
		double[] D = new double[params.size()];
		for (int j=0; j<params.size(); j++)
		{
			D[j] = -params.get(j).dL(rate_gradient)/likelihood_scaling;
		}
		return D;
	}	
	
	public double[] optDiff(double[] x)
	{
		set(x);
		double LL = gradient.getCorrectedLL();
		calls_optDiff++;

		double[] survival_gradient = gradient.getCorrectedGradient();
		double[] distribution_gradient = gradient.getDistributionGradient(survival_gradient);
		double[] rate_gradient = CALCULATE_FROM_TOTAL_RATES ? gradient.getTotalRateGradient(distribution_gradient)
									: gradient.getRateComponentGradient(distribution_gradient);
		
		double[] D = parameterGradient(rate_gradient) ; 
//		new double[params.size()];
//		for (int j=0; j<params.size(); j++)
//		{
//			D[j] = -params.get(j).dL(rate_gradient)/likelihood_scaling;
//		}
		if (PRINT_OPTIMIZATION_MESSAGES)
		{
			System.out.println("#*MLR.oD "+calls_optDiff+"\tfcalls "+calls_optFunc
					+"\t"+LL+"\trgrad "+FunctionMinimization.euclideanNorm(rate_gradient)
					+"\tdL "+FunctionMinimization.euclideanNorm(D));
		}

		return D;
	}
	
//	private static final double RATES_DIFF_EPS = 1e-7; // cubic root of machine precision
//	protected double[] estDiff(double[] x)
//	{
//		double fx = optFunc(x);
//		
//		double[] D = new double[x.length];
//		for (int p=0; p<x.length; p++)
//		{
//			double θ = x[p];
//			double h = Math.abs(θ*RATES_DIFF_EPS);
//			x[p] = θ+h;
//			double fd = optFunc(x);
//			double delta = fd-fx;
//			double dfdθ = delta/h;
//			D[p] = dfdθ;
//			
//			x[p] = θ;
//		}
//		set(x);
//		return D;
//		
//	}
	
	public double optimize(double delta)
	{
		return this.optimize(delta, Integer.MAX_VALUE);
	}
	public double optimize(double delta, int max_iter)
	{
		double LL =  gradient.getCorrectedLL();
		
		
		
		this.initModelParameters();
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#**MLR.o starting LL "+LL+"\tparams "+params.size());
		
		List<Double> history = getOptimizationHistory();
		if (history == null)
		{
			history = new ArrayList<>();
			setOptimizationHistory(history);
		}		
		
		if (max_iter>0)
		{
			double[] x0 = get();
			double min;   
			
			if (do_gradient_descent)
			{
				min = FunctionMinimization.gradientDescent(x0, Math.sqrt(delta), x0.length, x->optFunc(x), x->optDiff(x), history);
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#**MLR.o gd "+min);
			}
			if (optimization_by_quasiNewton)
			{
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#**MLR.o bfgs start");
				min = FunctionMinimization.dfpmin(x0, delta, max_iter, x->optFunc(x), x->optDiff(x), history);
			}
			else
			{
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#**MLR.o cg start");
				min = FunctionMinimization.frprmn(x0, delta, max_iter, x->optFunc(x), x->optDiff(x), history);
			}

			// min = FunctionMinimization.dfpmin(x0, delta, max_iter, x->optFunc(x), x->optDiff(x), history);
			
			
		}
		
		LL =  gradient.getCorrectedLL();
		double L0 = gradient.factory.getEmptyLL();
		double L1 = gradient.factory.getSingletonLL();
		
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#**MLR.o final LL "+LL+"\tfcalls "+calls_optFunc+"\tdcalls "+calls_optDiff+"\tL0 "+L0+"\tL1 "+L1);
//		set(x0);
//		LL =  gradient.factory.getCorrectedLL();
//		System.out.println("#** By dfpmin "+LL);
//		
		return -LL;
	}
	
	private static boolean[] optimizableParams(TreeWithRates rates)
	{
		IndexedTree tree = rates.getTree();
		int num_nodes = tree.getNumNodes();
		
		boolean[] opt_par = new boolean[4*num_nodes];
		for (int node=0; node<num_nodes; node++)
		{
			if (tree.isRoot(node))
			{
				opt_par[4*node + PARAMETER_GAIN] = rates.getGainParameter(node)>0.0;
				opt_par[4*node+PARAMETER_LOSS] = (rates.getLossParameter(node)<1.0);
				opt_par[4*node+PARAMETER_DUPLICATION] = rates.getDuplicationParameter(node)>0.0;
				opt_par[4*node+PARAMETER_LENGTH]=!Double.isInfinite(rates.getEdgeLength(node));
			} else
			{
				opt_par[4*node+PARAMETER_GAIN] = rates.getGainParameter(node)>0.0;
				opt_par[4*node+PARAMETER_LOSS]=false;
				opt_par[4*node+PARAMETER_DUPLICATION] = rates.getDuplicationParameter(node)>0.0;
				double t = rates.getEdgeLength(node);
				opt_par[4*node+PARAMETER_LENGTH] = t>0.0 && !Double.isInfinite(t);
			} 
		}
		
		return opt_par;
	}	
	
	private void debugGradient(PrintStream out)
	{
		optimize(0.0,0); 
		double[] x = get();
		
		double[] survival_gradient = gradient.getCorrectedGradient();
		double[] distribution_gradient = gradient.getDistributionGradient(survival_gradient.clone());
		
		for (int node=0; node<rates.getTree().getNumNodes(); node++)
		{
			double dLdp = distribution_gradient[3*node+PARAMETER_LOSS];
			double dLdq = distribution_gradient[3*node+PARAMETER_DUPLICATION];
			out.println("#GRADIENT\t"+node //+"\tdp~ "+survival_gradient[3*node+PARAMETER_LOSS]+"\tdq~ "+survival_gradient[3*node+PARAMETER_DUPLICATION]
					+"\tdp "+dLdp+"\tdq "+dLdq
						+"\t"+rates.toString(node));
			double μ = rates.getLossRate(node);
			double λ = rates.getDuplicationRate(node);
			double t = rates.getEdgeLength(node);	
			double μt = μ*t;
			double λt = λ*t;
			double d = (μ-λ)*t;
			double E = Math.exp(-d);
			double E1 = -Math.expm1(-d); // 1.0-E;
			
			double denom = μt-λt*E;     
			double divby = denom*denom;
			
			
			double dpdμ = E*(μt*d-λt*E1)/divby;
			double dpdλ = E*μt*(E1-d)/divby;
			double dqdμ = λt*(E*d-E1)/divby;
			double dqdλ = (E1*μt-λt*d*E)/divby;
			double zdt = E*d*d/divby;

			out.println("#GRADIENTr\t"+node+"\tE "+E+"\t(1-E) "+E1+"\tdenom "+denom+"\tdpdl "+dpdμ+"\tdqdl "+dqdμ+"\tdpdd "+dpdλ+"\tdqdd "+dqdλ
					+"\tzdt "+zdt
					+"\tdiv "+divby+"\tdLdl="+dLdp+"*"+dpdμ+"+"+dLdq+"*"+dqdμ);
		}
		
		
		double[] df_est = FunctionMinimization.numericalGradient(θ->optFunc(θ), x); //   estDiff(x);
		double[] df = optDiff(x);
		
		for (int i=0; i<x.length; i++)
		{
			ModelParameter P = params.get(i);
			
			double df_delta = df_est[i]-df[i];
			double rel_delta = Math.abs(df_delta/df[i]);
			
			out.println("#GRADIENT param "+i+"\t"+P+"\tdfdx "+df[i]+"\test "+df_est[i]
					+"\tdiff "+df_delta+"\t("+rel_delta+")");
		}
				
	}
	
//	private void shortenEdges(double scale)
//	{
//		for (int node=0; node<rates.getTree().getNumNodes(); node++)
//		{
//			double t = rates.getEdgeLength(node);
//			t *= (1.-scale);
//			rates.setEdgeLength(node, t);
//			double duprate = rates.getDuplicationRate(node);
//			rates.setDuplicationRate(node, Math.sqrt(duprate));
//		}
//	}
	
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
		ptest.add(0.033);
		ptest.add(0.05);
		ptest.add(0.08);
		ptest.add(0.1);
		ptest.add(0.12);
		ptest.add(0.15);
		ptest.add(0.18);
		ptest.add(0.2);
		ptest.add(0.25);
		ptest.add(0.3);
		ptest.add(0.35);
		ptest.add(0.4);
		ptest.add(0.45);
		ptest.add(0.5);
		ptest.add(0.6);
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
		System.out.println("#**MLR.tEL node\t"+node+"\tntest "+ntest+"\t"+ptest.toString()+"\t// "+ rates.toString(node));
		
		double negLL = optimize(delta, itmax);
		double[] model0 = get();
		double[] grad = gradient.getDistributionGradient(gradient.getCorrectedGradient());
		double[] rgrad = CALCULATE_FROM_TOTAL_RATES ? gradient.getTotalRateGradient(grad)
				: gradient.getRateComponentGradient(grad);
		double glen = FunctionMinimization.euclideanNorm(rgrad); 
		double dLdp = grad[3*node + PARAMETER_LOSS];
		
		double[] pgrad = parameterGradient(rgrad);
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
			
			boolean floss = fixLoss(node, true);
			boolean flength = fixLength(node, true);
			// gradient.computeParameters(); // not needed bc optimize() will do that
			negLL = optimize(delta, itmax);

			grad = gradient.getDistributionGradient(gradient.getCorrectedGradient());
			rgrad = CALCULATE_FROM_TOTAL_RATES ? gradient.getTotalRateGradient(grad)
					: gradient.getRateComponentGradient(grad);
			dLdp = grad[3*node + PARAMETER_LOSS];
			pgrad = parameterGradient(rgrad);

			dlen = FunctionMinimization.euclideanNorm(pgrad);
			glen = FunctionMinimization.euclideanNorm(rgrad);
			
			System.out.println("#NODELOSS\t"+node+"\t"+thisp+"\t"+negLL+"\t"+t+"\t"+setp+"\t"+dlen+"\t"+dLdp+"\t"+glen);

			fixLoss(node, floss);
			fixLength(node, flength);
			initModelParameters();
			set(model0); // and reset
		}
	}
	
	
	
	public static void main(String[] args) throws Exception
	{
		PRINT_OPTIMIZATION_MESSAGES = true;
		count.io.CommandLine cli = new count.io.CommandLine(args, MLRates.class);
		
		PrintStream out = System.out;
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(MLRates.class));
    	    out.println(CommandLine.getStandardRuntimeInfo(MLRates.class, args));
    	}
    	
    	
    	out.println(CommandLine.getStandardHeader("Gain rate optimization: "+(OPTIMIZE_COMBINED_GAIN_RATE?"combined":"separate")));
    	
    	GammaInvariant model = null; 
    	TreeWithRates base_rates;
    	if (cli.getModel()==null)
    	{
    		Random RND = cli.getOptionRND(out);
    		model = GammaInvariant.nullModel(cli.getTree(), RND);
    		base_rates = model.getBaseModel();
    		
    	} else
    	{
    		base_rates = cli.getRates(); 
    		model = cli.getModel();
    	}
		
		MLRates O = new MLRates(base_rates, cli.getTable());
		
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
		double eps = cli.getOptionDouble(OPT_EPS, 1e-4);
        
        out.println(CommandLine.getStandardHeader("Convergence: -"+OPT_ROUNDS+" "+maxiter+" -"+OPT_EPS+" "+eps));
        
		
//		if (cli.getOptionDouble("shorten", 0.0)!=0.0)
//		{
//			double scale = cli.getOptionDouble("shorten", 0.0);
//			System.out.println("#**MLR.main shortening by "+scale);
//			O.shortenEdges(scale);
//		}
        
        
		double score = O.optimize(eps,maxiter);
		
		
		out.println("#SCORE "+score);
		
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
