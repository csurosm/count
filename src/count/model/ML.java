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

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.beans.PropertyChangeSupport;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeEvent;

import count.ds.IndexedTree;

//import static count.matek.FunctionMinimization.EPS;
/**
 * Root class for model parameter setting by likelihood maximization, 
 * defining the common interface for extending classes.   
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */
public abstract class ML 
{
	/**
	 * Smallest power of 0.5 such that 
	 * 1.0+EPS &gt; 1.0 but 1.0+EPS/2.0==1.0
	 * at double precision
	 */
	public static final double EPS = Math.ulp(1.0);
	
//	static
//	{
//		double eps = 0.25;
//		double sqeps = 0.5;
//		while (1.0+eps/2.0 > 1.0)
//		{
//			eps = eps/2.0;
//			if (1.0+eps/2.0>1.0)
//			{
//				eps = eps/2.0;
//				sqeps = sqeps/2.0;
//			}
//		}
//		assert eps + (1.0-eps) == 1.0;
//		EPS = eps;
//	}
	
	public static final double MAX_RATE = 4.0;
	
	/**
	 * Constant to bound probability parameters away from 1.0
	 */
	public static final double MAX_PROB_NOT1 = 1.0-1.0/(1<<30); 
	
	
	
	public abstract double optimize(double delta);
	
	/**
	 * Default optimization ignores max_iterations.
	 * 
	 * @param delta
	 * @param max_iterations ignored here
	 * @return same as {@link #optimize(double)}
	 */
	public double optimize(double delta, int max_iterations)
	{
		return optimize(delta);
	}
	
	public abstract void setCalculationWidth(int absolute, double relative);
	
	public abstract void setMinimumObservedCopies(int m);
	
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
	
	public List<Double> getOptimizationHistory()
	{
		return  this.optimization_history;
	}
	
	public abstract void fixNodeParameters(int node, boolean do_not_optimize);

	public abstract int getModelParameterCount();
	
	
	/*
	 * Propetry change support
	 */
	public static final String PROPERTY_OPTIMIZATION_PHASE = "phase";
	
	private PropertyChangeSupport property_support=null;
	private Map<String,Object> property_value;
	
	private PropertyChangeSupport getSupport()
	{
		if (property_support==null)
		{
			property_support = new PropertyChangeSupport(this);
			property_value = new HashMap<>();
		}
		return property_support;
	}
	
	public synchronized void addPropertyChangeListener(String propertyName, PropertyChangeListener listener)
	{
		//getSupport().addPropertyChangeListener(propertyName, listener);
	}
	public synchronized void addPropertyChangeListener(PropertyChangeListener listener)
	{
		//getSupport().addPropertyChangeListener(listener);
	}
	protected synchronized void firePropertyChange​(String propertyName, Object oldValue, Object newValue)
	{
		if (property_support!=null) 
		{
			property_support.firePropertyChange(propertyName, oldValue, newValue);
			property_value.put(propertyName, newValue);
		}
	}
	protected synchronized void firePropertyChange​(String propertyName, Object newValue)
	{
		if (property_support!=null) 
		{
			Object oldValue = property_value.get(propertyName);
			if (!Objects.equals(newValue, oldValue))
			{
				property_support.firePropertyChange(propertyName, oldValue, newValue);
				property_value.put(propertyName, newValue);
			}
		}
	}
	protected synchronized void firePropertyChange​(PropertyChangeEvent e)
	{
		if (property_support!=null) 
		{
			property_support.firePropertyChange(e);
			property_value.put(e.getPropertyName(), e.getNewValue());
		}
	}
	public void removePropertyChangeListener(PropertyChangeListener listener)
	{
		if (property_support!=null) 
		{
			property_support.removePropertyChangeListener(listener);
			// set to null if empty?
		}
	}
	public void removePropertyChangeListener(String propertyName, PropertyChangeListener listener)
	{
		if (property_support!=null) 
		{
			property_support.removePropertyChangeListener(propertyName, listener);
			// set to null if empty?
		}
	}
	
	/* convenience methods  */
	
	
	/**
	 * Which of the 3/node parameters can be optimized. By default, 
	 * loss == 1.0 (at any node) is not optimized, and 
	 * duplication == 0.0 is not optimized at the root.
	 * 
	 * @param rates
	 * @param relative_dup whether relative duplication rate (q/p) or the duplication parameter (q) is to be optimized
	 * @return
	 */
	public static boolean[] optimizableParams(TreeWithRates rates, boolean relative_dup)
	{
		IndexedTree tree = rates.getTree();
		int num_nodes = tree.getNumNodes();
		
		boolean[] opt_par = new boolean[3*num_nodes];
		for (int node=0; node<num_nodes; node++)
		{
			double log_p = rates.getLogLossParameter(node);
			double log_q = rates.getLogDuplicationParameter(node);
			
			
			opt_par[3*node+PARAMETER_GAIN] = true; //rates.getGainParameter(node)>0.0;
			opt_par[3*node+PARAMETER_LOSS] = log_p<0.0;
			opt_par[3*node+PARAMETER_DUPLICATION] 
					= !tree.isRoot(node) || log_q != Double.NEGATIVE_INFINITY;
					//(log_q != Double.NEGATIVE_INFINITY && (!relative_dup || log_q<log_p));
		}
		
		return opt_par;
	}
	
	public static double[] getParameterValues(List<? extends ModelParameter> params)
	{
		int npar = params.size();
		double[] x = new double[npar];
		for (int j=0; j<npar; j++)
		{
			ModelParameter P = params.get(j);
			x[j] = P.get();
		}
		return x;
	}
	
	public static void setParameterValues(double[] x, List<? extends ModelParameter> params)
	{
		int npar = params.size();
		assert (npar<=x.length);
		for (int j=0; j<npar; j++)
		{
			ModelParameter P = params.get(j);
			P.set(x[j]);
			
		}
	}
	
	
	
	/**
	 * Model parameters. 
	 * 
	 * @author csuros
	 *
	 */
 	public static interface ModelParameter
	{
		/**
		 * Parameter value (<var>x</var>).
		 * 
		 * @return
		 */
		abstract double get();
		/**
		 * Parameter setting.
		 * @param x new value for parameter
		 */
		abstract void set(double x);
		/**
		 * Differentiation: implementation 
		 * calculates using total rate gradient,
		 * like from {@link Gradient#getTotalRateGradient(double[])}
		 * 
		 * @param rate_gradient
		 * @return
		 */
		default double dL(double[] rate_gradient)
		{
			throw new UnsupportedOperationException();
		}
	}
 	
 	/**
 	 * Dummy implementation hoklding a parameter value
 	 * 
 	 */
 	public static class SimpleParameter implements ModelParameter 
 	{
 		private double x;
 		
 		@Override
 		public double get() { return x;}
 		
 		@Override
 		public void set(double x) {this.x = x;}
 		
 		@Override
 		public String toString()
 		{
 			return "x="+x;
 		}
 	}
	
 	/**
 	 * Model parameter with simultaneous tracking of its complement. 
 	 * @author csuros
 	 *
 	 */
	public static interface  BoundedParameter extends ModelParameter
	{
		default double getBound()
		{
			return 1.0;
		}
		default double getComplement()
		{
			return getBound()-get();
		}
		default void set(double x, double xcomplement)
		{
			set(x);
		}
	}
	
	
	
	/**
	 * Sets up a parameter in the range epsilon ... (1-epsilon)
	 * 
	 * @param theta
	 * @param epsilon
	 * @return parameter for unconstrained optimization
	 */
	protected static Logistic bracketedLogistic(ModelParameter theta, double epsilon)
	{
		Bracketed bracket = new Bracketed(theta, epsilon);
		double y = bracket.get();
		if (y>1.0)
		{
			double oldθ = theta.get();
			bracket.set(MAX_PROB_NOT1);
			double newθ = theta.get();
			System.out.println("#*ML.bL "+theta+" too large "+oldθ+"; resetting");
		} else if (y<0.0)
		{
			double oldθ = theta.get();
			bracket.set(1.0-MAX_PROB_NOT1);
			double newθ = theta.get();
			System.out.println("#*ML.bL "+theta+" too small "+oldθ+"; resetting");
		}
		return new Logistic(bracket);
	}
	
	/**
	 * Transforms parameter  ε  ≤ θ ≤ (1-ε) to [0,1].
	 * 
	 * @author csuros
	 *
	 */
	protected static class Bracketed implements ModelParameter
	{
		protected Bracketed (ModelParameter theta, double eps) 
		{
			this.offset = eps;
			this.range = (1.0-2.0*eps);
			this.theta = theta;
		}
		private final ModelParameter theta;
		private final double offset;
		private final double range;

		@Override
		public double get()
		{
			double θ = theta.get();
			double y = (θ-offset)/range;
			if (y>1.0 || y<0.0) // DEBUG
			{
				System.out.println("#**ML.B.get "+theta+"\ty " +y);
			}
			return y;
		}
		
		@Override
		public void set(double y)
		{
			double θ = y*range + offset;
			theta.set(θ);
			if (y>1.0 || y<0.0) // DEBUG
			{
				System.out.println("#**ML.B.set "+theta+"\ty" +y);
				
			}
		}
		@Override 
		public double dL(double[] gradient)
		{
			double d = theta.dL(gradient);
			return d*range;
		}

		@Override 
		public String toString()
		{
			return this.getClass().getSimpleName()+"["+theta.toString()+"]-"+offset+"/"+range;
		}
	}
	
	
//	/**
//	 * Logistic transformations for range [ε,1-ε]
//	 *
//	 */
//	public static class BracketedLogistic implements ModelParameter
//	{
//		BracketedLogistic(ModelParameter theta, double eps, double max_value)
//		{
//			this.eps = eps;
//			this.scale = (1.0-2.0*eps)*max_value;
//			this.ltheta = new Logistic()
//		}
//		
//		final double eps;
//		final double scale;
//		final Logistic ltheta;
//	}
	
	/**
	 * Logistic with maximum MAX_RATE
	 * @param theta
	 * @return
	 */
	protected Logistic newLogistic(ModelParameter theta)
	{
		return newLogistic(theta, MAX_RATE);
	}
	
	protected Logistic newLogistic(ModelParameter theta, double max_value)
	{
		double θ = theta.get();
		if (θ>=max_value)
		{
			double newθ =MAX_PROB_NOT1*max_value;
			System.out.println("#*ML.newLogist "+theta+"; resetting to "+newθ);
			theta.set(newθ);
		}
		Logistic newLogistic = new Logistic(theta, max_value);
		return newLogistic;
	}
	
	/**
	 * Logistic transformation for parameter 0 ≤ θ ≤ A.
	 * 
	 */
	protected static class Logistic implements ModelParameter
	{
		protected Logistic(ModelParameter theta, double max_value)
		{
			this.theta = theta;
			this.max_value = max_value;
		}
		
		protected Logistic(ModelParameter theta)
		{
			this(theta, 1.0);
		}
		
		private final ModelParameter theta;
		private final double max_value;
		
		public double getOriginal() { return theta.get();}
		
		@Override
		public double get()
		{
			double t = theta.get()/max_value;
			
//			if (t>1.0)
//			{
//				System.out.println("#*ML.L.get "+this+"\tt "+t);
//			}
			
			assert (t <= 1.0);
			double x = Math.log(t)-Math.log1p(-t);
			return x;
		}
		
		@Override
		public void set(double x)
		{
			double t;
			if (x>=0.0)
			{
				t = 1.0/(1.0 + Math.exp(-x));
			} else
			{
				double e =Math.exp(x);
				t = e/(1.0+e);
			}
			if (Double.isFinite(x))
			{
//				if (t<EPS || t>1.0-EPS)
//				{
//					System.out.println("#*ML.L.set "+this+"\tx "+x+"\tt "+t+"\teps "+EPS+"\t1-e "+(1.0-EPS));
//				}
				t = Double.min(Double.max(EPS, t), 1.0-EPS);
				
			}
			if (!(Double.isInfinite(x) || (0.0<t && t<1.0))) //(t<=0.0 || t>=1.0) && !Double.isInfinite(x))
			{
				System.out.println("#*ML.L.set "+this+"\tx "+x+"\tt "+t+"\t"+(0.0<t)+","+(t<1.0));
				// 
			}
//			assert (Double.isInfinite(x) || (0.0<t && t<1.0)); // the formulas should never give exactly 0 or 1 for finite x 
			
			
			theta.set(max_value*t);
		}
		
		@Override 
		public double dL(double[] rate_gradient)
		{
			double d = theta.dL(rate_gradient);
			double θ = theta.get();
			d *= θ*(1.0-θ/max_value);
			
			return d;
		}

		@Override 
		public String toString()
		{
			return this.getClass().getSimpleName()+"["+theta.toString()+"; get "+get()
				+(max_value==1.0?"":"; max "+max_value)+"]";
		}
	
	}
	

	/**
	 * Logistic transformation for range[0,max]. 
	 * @author csuros
	 *
	 */
	protected static class BoundedLogistic implements ModelParameter
	{
		protected BoundedLogistic(BoundedParameter theta)
		{
			this.theta = theta;
			this.scale = theta.getBound();
			double θ = theta.get();
			if (θ>=scale)
			{
				double newθ =MAX_PROB_NOT1*scale;
				System.out.println("#*ML.BoundedLogistic "+theta+"; resetting to "+newθ);
				theta.set(newθ);
			}
		}
		private final BoundedParameter theta;
		private final double scale;
		
		@Override
		public double get()
		{
			double t = theta.get();
			double t_1 = theta.getComplement();
			
			double x = Math.log(t)-Math.log(t_1);
			
			return x;
		}
		
		@Override
		public void set(double x)
		{
			double t, t_1;
			if (x>=0)
			{
				double e = Math.exp(-x);
				t = 1.0/(1.0 + e);
				t_1 = e*t;
			} else
			{
				double e =Math.exp(x);
				t_1 = 1.0/(1.0+e);
				t = e*t_1;
			}
			t *= scale;
			t_1 *= scale;
			theta.set(t, t_1);
		}
		
		@Override 
		public double dL(double[] rate_gradient)
		{
			double d = theta.dL(rate_gradient);
			double t = theta.get();
			double t_1 = theta.getComplement();
			
			d *= t*t_1/scale;
			return d;
		}
	}
	

	
	
	
	
	
//	protected static class EpsilonLogistic implements ModelParameter
//	{
//		protected EpsilonLogistic(ModelParameter theta)
//		{
//			this.theta = theta;
//		}
//		
//		private final ModelParameter theta;
//		
//		@Override
//		public double get()
//		{
//			double t = theta.get();
//			
//			assert (t <= 1.0);
//			
//			assert (t>=EPS/2.0);
//			t-=EPS/2.0;
//			double x = Math.log(t)-Math.log1p(-t);
//			return x;
//		}
//		
//		@Override
//		public void set(double x)
//		{
//			double t;
//			if (x>=0.0)
//			{
//				t = 1.0/(1.0 + Math.exp(-x));
//			} else
//			{
//				double e =Math.exp(x);
//				t = e/(1.0+e);
//			}
////			if (t==0.0 || t==1.0 || true)
////				System.out.println("#*ML.L.set "+this+"\tx "+x+"\tt "+t);
//			theta.set(t+EPS/2.0);
//		}
//		
//		@Override 
//		public double dL(double[] rate_gradient)
//		{
//			double d = theta.dL(rate_gradient);
//			double θ = theta.get();
//			d *= θ*(1.0-θ);
//			return d;
//		}
//
//		@Override 
//		public String toString()
//		{
//			return this.getClass().getSimpleName()+"["+theta.toString()+"]";
//		}
//	}

	
	
//	/**
//	 * Transformation for parameter 0 ≤ θ; 
//	 * transformed <var>x</var> can be arbitrary.
//	 *
//	 *
//	 */
//	protected static class LogLogistic implements ModelParameter
//	{
//		static final double TOO_BIG;
//		static 
//		{
//			double big = 10;
//			double exp = Math.exp(big);
//			while (exp>exp-1.0)
//			{
//				big += 1.0;
//				exp = Math.exp(big);
//			}
//			TOO_BIG = big;
//		}
//
//		protected LogLogistic(ModelParameter theta)
//		{
//			this.theta = theta;
//		}
//		
//		private final ModelParameter theta;
//		
//		@Override
//		public double get()
//		{
//			
//			
//			double θ = theta.get();
//			double x = θ>TOO_BIG
//					?θ:
//					Math.log(Math.expm1(θ)); // ln(e^t-1)
////			System.out.println("#*MLR.LL.get "+this+"\tu "+θ+"\tx "+x);
//			
//			return x;
//		}
//		
//		@Override
//		public void set(double x)
//		{
//			double θ = x>TOO_BIG
//					?x
//					:Math.log1p(Math.exp(x)); // ln(1+e^x)
////			System.out.println("#*MLR.LL.set "+this+"\tu "+θ+"\tx "+x);
//			theta.set(θ);
//		}
//		
//		@Override 
//		public double dL(double[] rate_gradient)
//		{
//			double d = theta.dL(rate_gradient);
//			double θ = theta.get();
//			d *= -Math.expm1(-θ); // d*(1-exp(-t))
//			return d;
//		}		
//		@Override 
//		public String toString()
//		{
//			return this.getClass().getSimpleName()+"["+theta.toString()+"]";
//		}
//	}
	
	/**
	 * Logarithmic transformation for parameter 0 ≤ θ
	 * 
	 */
	protected static class Logarithmic implements ModelParameter
	{
		protected Logarithmic(ModelParameter theta, double scaling_factor)
		{
			this.theta = theta;
			this.scaling_factor = scaling_factor;
		}
		protected Logarithmic(ModelParameter theta)
		{
			this(theta, 1.0);
		}
		
		protected ModelParameter ogParameter() { return theta;}
		
		
		private final ModelParameter theta;
		private final double scaling_factor;
		
		@Override
		public double get()
		{
			double t = theta.get();
			double y = Math.log(t);
			return y/scaling_factor;
		}
		
		@Override
		public void set(double x)
		{
			double y = scaling_factor*x;
			double t = Math.exp(y);
			if (!Double.isFinite(t))
				System.out.println("#***ML.L.set big argument "+this+"\tx "+x+"\te^y=t "+t+"\ty "+y);

			assert (!Double.isNaN(x));
			theta.set(t);
		}
		
		@Override 
		public double dL(double[] rate_gradient)
		{
			double d = theta.dL(rate_gradient);
			double t = theta.get();
			d *= t;
			return d*scaling_factor;
		}		
		@Override 
		public String toString()
		{
			return this.getClass().getSimpleName()+"["+theta.toString()+"]";
		}
		
	}
}
