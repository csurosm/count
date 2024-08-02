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


import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LOSS;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;

import count.ds.IndexedTree;
import count.matek.Logarithms;

/**
 * Rate variation model with experimental category modifiers. 
 * Model in main BFGS implementation {@link MLRateVariation}.
 *
 */
public class RateVariationModel implements MixedRateModel, Iterable<RateVariationModel.Category>
{
	private static boolean DEBUG_NUMERICAL = true;

	public RateVariationModel(TreeWithLogisticParameters base_rates)
	{
		this.main_rates = base_rates;
		this.rate_categories = new ArrayList<>();
	}
	
	public RateVariationModel(TreeWithRates base_rates)
	{
		this(new TreeWithLogisticParameters(base_rates, false)); // hard link 
	}
	
	private final TreeWithLogisticParameters main_rates; 
	/**
	 * Rate-variation categories
	 */
	private final List<Category> rate_categories;

	private Class<? extends Category> default_category_type = LogisticShift.class;

	
	/*
	 * Common gain
	 * 
	 */
	private int common_gain_by = PARAMETER_GAIN; //PARAMETER_DUPLICATION;
	private boolean use_universal_gain = false;
	
	/*
	 * Setters and getters
	 */
	@Override
	public TreeWithLogisticParameters getBaseModel()
	{
		return this.main_rates;
	}

	/**
	 * Sets the policy for common gain rate: it is bound  
	 * either to PARAMETER_LOSS (<var>γ</var>), or to PARAMETER_GAIN (<var>r</var>), or
	 * to PARAMETER_DUPLICATION (<var>κ</var>/<var>r</var> if no dup); 
	 * <var>γ</var>=<var>r</var>/<var>p</var>, and <var>r</var> 
	 * computed as universal or linear gain: 
	 * <ul>
	 * <li>linear gain: <var>r</var>=<var>κ</var><var>q</var>
	 * <li>universal gain: <var>r</var>=<var>κ</var>(-ln(1-<var>q</var>))
	 * </ul> 
	 * 
	 * @param common_gain_by gain bound to which parameter
	 * @param universal_gain whether universal or linear gain
	 */
	public void setCommonGain(int common_gain_by, boolean universal_gain)
	{
		this.common_gain_by = common_gain_by;
		this.use_universal_gain = universal_gain;
	}

	public void setCommonGain(int common_gain_by)
	{
		this.setCommonGain(common_gain_by, false);
	}
	
	
	public int getCommonGainType()
	{
		return this.common_gain_by;
	}
	
	public boolean isUniversalGain()
	{
		return this.use_universal_gain;
	}
	
	/*
	 * Rate categories
	 */
	@Override
	public int getNumClasses()
	{
		return rate_categories.size();
	}
	
	public Category getCategory(int cat)
	{
		return rate_categories.get(cat);
	}
	
	@Override
	public double getClassProbability(int cat)
	{
		return Math.exp(getCategory(cat).getLogCatProbability());
	}
	
	@Override
	public TreeWithRates getClassModel(int cat)
	{
		return getCategory(cat).getRates();
	}
	
	public <C extends Category> void setDefaultType(Class<C> default_type)
	{
		this.default_category_type = default_type;
	}

	public void initCategories(double[] p_c, double[] mod_length, double[] mod_duplication)
	{
		int ncat = p_c.length;
		assert (mod_length.length == ncat);
		assert (mod_duplication.length == ncat);
		
		rate_categories.clear(); 
		
		for (int k=0; k<ncat; k++)
		{
			if (default_category_type == LogisticShift.class)
			{
				rate_categories.add(new LogisticShift(p_c[k], mod_length[k], mod_duplication[k]));
			} else
			{
				assert (default_category_type == Multiplier.class);
				rate_categories.add(new Multiplier(p_c[k], mod_length[k], mod_duplication[k]));
			}
		}
	}
	
	/**
	 * A single category.
	 */
	public void initConstantRates()
	{
		double[] p_c = {1.0};
		double[] mod_length = {0.0};
		double[] mod_dup = {0.0};
		this.initCategories(p_c, mod_length, mod_dup);
	}
	
	public void initRandomCategories(int num_categories, Random rnd)
	{
		if (rnd == null) rnd = new Random();
		rate_categories.clear();
		for (int k=0; k<num_categories; k++)
		{
			double mod_len, mod_dup;
			if (k==0)
			{
				mod_len = 0.0;
				mod_dup = 0.0;
			} else
			{
				double u = rnd.nextDouble(); // Uniform(0,1)
				double x = -Math.log(u); // Exponential(1)
				mod_len = Math.log(x);
				mod_dup = Math.log(-Math.log(rnd.nextDouble()));
			}
			Category C = this.addCategory(default_category_type, mod_len, mod_dup);
			// DEBUG
			System.out.println("#**RVM.iRC "+k+"\t"+C);
		}
	}
	
	public <Type extends Category> void addRandomCategory(Class<Type> type, Random rnd)
	{
		if (rnd == null) rnd = new Random();
		double mod_len = Math.log(-Math.log(rnd.nextDouble()));
		double mod_dup = Math.log(-Math.log(rnd.nextDouble()));

		int k = this.getNumClasses();
		Category C = this.addCategory(type, mod_len, mod_dup);
		// DEBUG
		System.out.println("#**RVM.aRC "+k+"\t"+C);
	}	
		
	
	/**
	 * Adds a new category, and reweighs the probabilities of existing ones
	 * so that they sum to 1.
	 * 
	 * @param C
	 */
	public Category addCategory(Category C)
	{
		for (Category known: rate_categories)
			if (known.equals(C))
				return known;

		double log_p = C.getLogCatProbability();
		double log1_p = Logarithms.logToLogComplement(log_p);
		for (Category previous: rate_categories)
		{
			previous.setLogCatProbability(previous.getLogCatProbability()+log1_p);
		}
		rate_categories.add(C);
		return C;
	}
	
	
	public <Type extends Category> Type addCategory(Class<Type> type, double mod_length, double mod_duplication)
	{
		int ncat = rate_categories.size();
		double pcat = 1.0/(1.0+ncat);
		Category C;
		if (type == LogisticShift.class)
		{
			C = new LogisticShift(pcat, mod_length, mod_duplication);
		} else if (type == Multiplier.class)
		{
			C = new Multiplier(pcat, mod_length, mod_duplication);
		} else
			throw new IllegalArgumentException("Type "+type+" is unknown here.");

		C = this.addCategory(C);
		return type.cast(C);
	}
	
	public static RateVariationModel convert(MixedRateModel.RateMultipliers model)
	{
		RateVariationModel convert=new RateVariationModel(model.getBaseModel());
		convert.initCategories(model);
		return convert;
	}
	
	/**
	 * An approximation of length- and duplication categories in a multiplier model.
	 * 
	 * @param gamma_invariant
	 */
	public void initCategories(MixedRateModel.RateMultipliers gamma_invariant)
	{
		rate_categories.clear();
		Map<Category,Double> p_cat = new HashMap<>();
		for (int k=0; k<gamma_invariant.getNumClasses(); k++)
		{
			double p =gamma_invariant.getClassProbability(k);  
			if (p!=0.0)
			{
				double mul_len = gamma_invariant.getEdgeLengthMultiplier(k);
				double mul_dup = gamma_invariant.getDuplicationRateMultiplier(k);
				Category C = addCategory(LogisticShift.class, Math.log(mul_len), Math.log(mul_dup));
				if (p_cat.containsKey(C))
				{
					p_cat.put(C, p_cat.get(C)+p);
				} else
				{
					p_cat.put(C, p);
				}
			}
		}
		for (Category C: rate_categories)
		{
			C.setLogCatProbability(Math.log(p_cat.get(C)));
		}
	}
	
	/**
	 * Resets the base model to the highest-probability category, and shifts 
	 * the category modifiers accordingly.
	 * 
	 * @return true if base model was reset
	 */
	public boolean recenterCategories()
	{
		if (rate_categories.size()<2) return false;
		
		Category center = null;
		double max_log_p=Double.NEGATIVE_INFINITY;
		int max_idx = -1;

		for (int k=0; k<rate_categories.size(); k++)
		{
			Category C = rate_categories.get(k);
			double log_p = C.log_pclass;
			double mod_len = C.getModLength();
			double mod_dup = C.getModDuplication();
			if (max_log_p<log_p
					|| (max_log_p == log_p && mod_len == 0.0 && mod_dup == 0.0))
			{
				max_log_p = log_p;
				center = C;
				max_idx = k;
			}
		}
		double shift_mod_len = center.getModLength();
		double shift_mod_dup = center.getModDuplication();
		
		boolean recenter = shift_mod_len != 0.0 || shift_mod_dup != 0.0;
		
		if (recenter)
		{
			// reset main rates to center's model
			int num_nodes = main_rates.getTree().getNumNodes();
			for (int v=0; v<num_nodes; v++)
			{
				double logit_p = center.lrates.getLogitLossParameter(v);
				double logit_lm = center.lrates.getLogitRelativeRate(v);
				
				double log_common_gain = center.lrates.getLogGainParameter(v, getCommonGainType(), !isUniversalGain());
				
				main_rates.setLogitLossRelativeDuplication(v, logit_p, logit_lm, log_common_gain, getCommonGainType(), !isUniversalGain());
				
			}
		
			for (Category C: rate_categories)
			{
				assert (center.getClass().equals(C.getClass())); // same type of mods 
				
				double mod_len = C.getModLength();
				double mod_dup = C.getModDuplication();
				C.setModifiers(mod_len-shift_mod_len, mod_dup-shift_mod_dup);
				C.computeParameters(); // although category's rates (C.lrates) don't change if calculations are exact
			}
			
//			if (max_idx != 0)
//			{
//				// swap so that category 0
//				// has maximum probability
//				Category headC = rate_categories.get(0);
//				rate_categories.set(max_idx, headC);
//				rate_categories.set(0, center);
//			}
		}
		return recenter;
	}
	
	/**
	 * Deep copy of the model
	 * 
	 * @return
	 */
	public RateVariationModel copy()
	{
		RateVariationModel copy = new RateVariationModel(new TreeWithLogisticParameters(main_rates,true));
		copy.setCommonGain(common_gain_by, use_universal_gain);
		int ncat = 0;
		for (Category C: rate_categories)
		{
			copy.addCategory(C.getClass(), C.getModLength(), C.getModDuplication());
			++ncat;
		}
		for (int k=0; k<ncat; k++)
		{
			Category C = this.getCategory(k);
			copy.getCategory(k).setLogCatProbability(C.getLogCatProbability());
		}
		return copy;
	}
	
	private static int NEXT_CATID=1;
	public abstract class Category 
	{
		private final int cat_id;

		private double mod_length=0.0;
		private double mod_duplication=0.0;
		private double log_pclass;

		protected TreeWithLogisticParameters lrates; 
//		private LogGradient gradient;

		private Category(double p, double mod_length, double mod_duplication)
		{
			this.cat_id = RateVariationModel.this.NEXT_CATID++;
			this.log_pclass = Math.log(p);
			this.mod_length = mod_length;
			this.mod_duplication = mod_duplication;
			this.lrates = new TreeWithLogisticParameters(main_rates, true); // make a copy
		}

		public double getLogCatProbability()
		{
			return this.log_pclass;
		}
		public void setLogCatProbability(double log_p)
		{
			assert (log_p <= 0.0);
			this.log_pclass = log_p;
		}
		
		public double getModLength()
		{
			return this.mod_length;
		}

		public double getModDuplication()
		{
			return this.mod_duplication;
		}

		public void setModifiers(double mod_len, double mod_dup)
		{
			if (this.mod_length != mod_len|| this.mod_duplication != mod_dup)
			{
				this.mod_length = mod_len;
				this.mod_duplication = mod_dup;
				
//				System.out.println("#**RVM.C.sM "+this.toString()); // DEBUG
			}
		}
		
		/**
		 * Updates the category's rate model
		 */
		public void computeParameters()
		{
			int num_nodes = main_rates.getTree().getNumNodes();
			for (int v=0; v<num_nodes; v++)
			{
				this.updateNodeParameters(v);
			}
		}
		
		/**
		 * Updates the category's rate model, and returns access to it.
		 * 
		 * @return
		 */
		public TreeWithLogisticParameters getRates()
		{
			this.computeParameters();
			return this.lrates;
		}
		
		/**
		 * Computes the parameters in the underlying model {@link #lrates}.
		 * Extending classes call {@link #setLogitLossRelativeDuplication(int, double, double)}
		 * at the node.
		 * 
		 * @param node
		 */
		protected abstract void updateNodeParameters(int node);
		
		/**
		 * Sets node parameters in underlying model {@link #lrates}, respecting the common gain policy.
		 * Called by extending classes from {@link #updateNodeParameters(int)}.
		 * 
		 * @param v node 
		 * @param cat_logit_p
		 * @param cat_logit_λ
		 */
		protected void setLogitLossRelativeDuplication(int v, double cat_logit_p, double cat_logit_λ)
		{
			double gr = main_rates.getGainRate(v);
			double log_gr = Math.log(gr);
			
			if (cat_logit_λ == Double.NEGATIVE_INFINITY) // Poisson
			{
				double cat_log_r;
				if (RateVariationModel.this.common_gain_by == PARAMETER_LOSS)
				{
					// same gamma
					cat_log_r = log_gr + Logarithms.logitToLogValue(cat_logit_p);
				} else
				{
					// same r 
					cat_log_r =  log_gr + main_rates.getLogLossParameter(v);
				}
				
				double r = Math.exp(cat_log_r);
				
				if (DEBUG_NUMERICAL && !Double.isFinite(r)) // DEBUG
				{
					System.out.println("#**RVM.C.sLLRD "+v+"\tr "+r+"\tlogr "+cat_log_r+"\tlogitp "+cat_logit_p+"\tlogitlm "+cat_logit_λ+"\t"+main_rates.toString(v));
				}
				if (DEBUG_NUMERICAL) assert(Double.isFinite(r));
				
				lrates.setLogitLossRelativeDuplication(v, cat_logit_p, cat_logit_λ, Math.exp(cat_log_r));
			} else
			{
				// Polya
				if (cat_logit_λ != Double.POSITIVE_INFINITY)
				{				
					cat_logit_p += main_rates.getLogRelativeComplement(v)-Logarithms.logitToLogComplement(cat_logit_λ);				
				}
				
				double cat_log_kappa;

				if (RateVariationModel.this.common_gain_by == PARAMETER_DUPLICATION)
				{ // same kappa
					cat_log_kappa = log_gr;
				} else 
				{
					// same r 
					assert (PARAMETER_GAIN==common_gain_by) || (PARAMETER_LOSS == common_gain_by);
					double log_r;
					
					if (RateVariationModel.this.use_universal_gain)
					{
						double loglog1_q = Logarithms.logitToLogLogComplement(main_rates.getLogitDuplicationParameter(v));
								//Math.log(-main_rates.getLogDuplicationComplement(v));
						log_r = log_gr+loglog1_q;
					} else
					{
						log_r = log_gr + main_rates.getLogDuplicationParameter(v);
					} 
					
					double cat_logit_q = Logarithms.mulLogit(cat_logit_p, cat_logit_λ);
					if (RateVariationModel.this.common_gain_by == PARAMETER_GAIN)
					{
						// same r 
						
						if (RateVariationModel.this.use_universal_gain)
						{
							double cat_loglog1_q = Logarithms.logitToLogLogComplement(cat_logit_q);
							if (DEBUG_NUMERICAL)
							{
								// log(-log(1-q)) == -infty only if log(1-q)=0 or q=0; log(1-q)==0.0 in double precision if log(q) is smaller than log(min double)
								assert (cat_loglog1_q != Double.NEGATIVE_INFINITY || Logarithms.logitToLogValue(cat_logit_q)==Double.NEGATIVE_INFINITY);
							}
							cat_log_kappa = log_r -  cat_loglog1_q;
						} else
						{
							cat_log_kappa = log_r - Logarithms.logitToLogValue(cat_logit_q);
						}
					} else
					{
						assert (PARAMETER_LOSS == common_gain_by);
						// same gamma
						double log_gamma = log_r - main_rates.getLogLossParameter(v);
						if (RateVariationModel.this.use_universal_gain)
						{
							double cat_log_r = log_gamma+Logarithms.logitToLogValue(cat_logit_p);
							double cat_loglog1_q = Logarithms.logitToLogLogComplement(cat_logit_q);
							if (DEBUG_NUMERICAL)
							{
//								System.out.println("#**RVM.C.sLLRD smallq "
//										+"\tnode "+v
//										+ "\tclogq "+ Logarithms.logitToLogValue(cat_logit_q)
//										+ "\tclogit "+cat_logit_q
//										+"\tloglog1q "+cat_loglog1_q
//										+"\tcat "+this.toString()
//										+"\t// "+main_rates.toString(v)
//										+"\t// "+lrates.toString(v)
//										);
//								
								// log(-log(1-q)) == -infty only if log(1-q)=0 or q=0; log(1-q)==0.0 in double precision if log(q) is smaller than log(min double)
								// so log(-log(1-q)) == -infty must imply log(p)=-infty 
								assert (cat_loglog1_q != Double.NEGATIVE_INFINITY || Logarithms.logitToLogValue(cat_logit_q)==Double.NEGATIVE_INFINITY);
								
								
								
							}
							cat_log_kappa = cat_log_r -  cat_loglog1_q;
						} else
						{
							cat_log_kappa = log_gamma - Logarithms.logitToLogValue(cat_logit_λ);
						}
					}
				}
				double kappa = Math.exp(cat_log_kappa);
				if (DEBUG_NUMERICAL && !Double.isFinite(kappa)) // DEBUG
				{
					System.out.println("#**RVM.C.sLLRD "+v+"\tkappa "+kappa+"\tlogitp "+cat_logit_p+"\tlogitlm "+cat_logit_λ+"\t"+main_rates.toString(v));
				}
				if (DEBUG_NUMERICAL) assert Double.isFinite(kappa);
				
				lrates.setLogitLossRelativeDuplication(v, cat_logit_p, cat_logit_λ,kappa);
			} // Polya
		}
		

		@Override
		public String toString()
		{
			StringBuilder sb = new StringBuilder(getClass().getSimpleName());
			sb.append("#")
			.append(this.cat_id)
			.append("[")
			.append(this.mod_length);
			if (mod_length != 0.0)
				sb.append("(").append(mod_length<0.0?"/"+Math.exp(-mod_length):"*"+Math.exp(mod_length)).append(")");
			sb.append(",")
			.append(this.mod_duplication);
			if (mod_duplication!=0.0)
				sb.append("(").append(mod_duplication<0.0?"/"+Math.exp(-mod_duplication):"*"+Math.exp(mod_duplication)).append(")");
			sb.append("; p=").append(Math.exp(log_pclass)).append("/logp=").append(log_pclass);
			sb.append("]");
			return sb.toString();
		}
		
		@Override 
		public int hashCode()
		{
			int hclass = getClass().hashCode();
			int hmod = Double.hashCode(mod_length)*17+Double.hashCode(mod_duplication);
			return hclass*31 + hmod;
		}
		
		@Override
		public boolean equals(Object o)
		{
			if (o == null) return false;
			if (o instanceof Category)
			{
				Category that= (Category) o;
				return this.getClass().equals(that.getClass())
						&& this.mod_length == that.mod_length
						&& this.mod_duplication == that.mod_duplication;
			} else 
			{
				return super.equals(o);
			}
		}
		
//		/**
//		 * Calculates gradient by base model's logit-p, logit-lambda, log-gain; 
//		 * and rate modifiers. The last two entries of the returned array
//		 * are the derivatives by the {@link #mod_length} and {@link #mod_duplication}
//		 * with indexes 3<var>n</var>+{@link RateVariationLogGradient#PARAMETER_MOD_LENGTH}
//		 * and 3<var>n</var>+ {@link RateVariationLogGradient#PARAMETER_MOD_DUPLICATION}.
//		 * 
//		 * Default implementation just adds then last two cells with 0 gradient; useful for 
//		 * initialization.
//		 * 
//		 * @param logCatGradient category gradient from {@link #convertToLogCatGradient(double[][])}; unaffected
//		 * 	
//		 * @return array of length 3<var>n</var>+2 for <var>n</var> nodes in the tree
//		 */
//		protected double[][] convertToLogMainGradient(double[][] logCatGradient)
//		{
//			int num_nodes = lrates.getTree().getNumNodes();
//			double[][] logMD = new double[3*num_nodes+2][];
//			int j=0; 
//			while (j<3*num_nodes)
//			{
//				logMD[j]=logCatGradient[j].clone();
//				j++;
//			}
//			logMD[j++] = Logarithms.ldiff();
//			logMD[j++] = Logarithms.ldiff();			
//			return logMD;
//		}
//		
//		/**
//		 * Log-gradient by base model's logit(p), logit(lambda), log(common gain)
//		 * and rate modifiers. The last two entries of the returned array
//		 * are the derivatives by the {@link #mod_length} and {@link #mod_duplication}
//		 * with indexes 3<var>n</var>+{@link RateVariationLogGradient#PARAMETER_MOD_LENGTH}
//		 * and 3<var>n</var>+ {@link RateVariationLogGradient#PARAMETER_MOD_DUPLICATION}.
//		 * 
//		 * 
//		 * @param logSD log-gradient by survival parameters 
//		 * @return array of length 3<var>n</var>+2 for <var>n</var> nodes in the tree
//		 */
//		protected double[][] getLogMainGradientCommonGain(double[][] logSD)
//		{
//			double[][] log_Dp = gradient.convertToLogDistributionGradient(logSD);
//			
//			double[][] log_Dpλ = LogGradient.convertToLogRelativeRateGradient(lrates, log_Dp);
//			log_Dpλ = LogGradient.convertToLogGainGradient(lrates, log_Dpλ, RateVariationLogGradient.this.common_gain_by, RateVariationLogGradient.this.use_universal_gain);
//			
//			double[][] logD = this.convertToLogMainGradient(log_Dpλ);
//			
//			return logD;
//		}
		
	
	} // Category class	

	/**
	 * This rate variation model shifts the logit parameters, indirectly 
	 * changing edge length (short edges as if a multiplier, long edges 
	 * less so).  
	 * The loss-gradients is calculated 
	 * in one step with a well-scaled Jacobian for converting from 
	 * category-specific parameters to base-model parameters.  
	 */
	public class LogisticShift extends Category
	{
		private LogisticShift(double p, double mod_length, double mod_duplication)
		{
			super(p, mod_length, mod_duplication);
		}
		
		@Override
		protected void updateNodeParameters(int v)
		{
			double logit_λ = main_rates.getLogitRelativeRate(v);
			double logit_p = main_rates.getLogitLossParameter(v);
			
			double cat_logit_p = logit_p + getModLength(); //mod_length;
			
			if (logit_λ == Double.NEGATIVE_INFINITY) // Poisson
			{
				this.setLogitLossRelativeDuplication(v, cat_logit_p, logit_λ);
			} else
			{
				// Polya
				double cat_logit_λ = logit_λ + getModDuplication();// mod_duplication;
				if (logit_λ != Double.POSITIVE_INFINITY)
				{				
					// add log(1-lambda)-log(1-cat_lambda)
					cat_logit_p += main_rates.getLogRelativeComplement(v)-Logarithms.logitToLogComplement(cat_logit_λ);				
				}
				this.setLogitLossRelativeDuplication(v, cat_logit_p, cat_logit_λ);
			} // Polya
		} // updateNodeParameters
	} // LogisticShift
		
	/**
	 * This rate variation model multiplies the effective edge length -log(transient) 
	 * directly. The loss-gradient is calculated in two steps of converting to 
	 * and from transient-gradients.  
	 */
	public class Multiplier extends Category
	{
		private Multiplier(double p, double mod_length, double mod_duplication)
		{
			super(p, mod_length, mod_duplication);
		}
		
		@Override
		protected void updateNodeParameters(int v)
		{
			IndexedTree tree = lrates.getTree();
			double logit_λ = main_rates.getLogitRelativeRate(v);
			double logit_p = main_rates.getLogitLossParameter(v);
			
			double mod_length = getModLength();
			double mod_duplication = getModDuplication();
			
			if (tree.isRoot(v) || logit_p == Double.POSITIVE_INFINITY)
			{
				// same p 
				double cat_logit_p = logit_p;
				if (logit_λ == Double.NEGATIVE_INFINITY) // Poisson
				{
					this.setLogitLossRelativeDuplication(v, cat_logit_p, logit_λ);
				} else
				{
					// Polya
					double cat_logit_λ = logit_λ + mod_duplication;
					this.setLogitLossRelativeDuplication(v, cat_logit_p, cat_logit_λ);
				}
			} else
			{ // p<1.0
				if (logit_λ == Double.POSITIVE_INFINITY) // p==q
				{
					double cat_logit_p = logit_p + mod_length;
					this.setLogitLossRelativeDuplication(v, cat_logit_p, logit_λ);
				} else
				{
					double cat_logit_λ = logit_λ + mod_duplication;
				
					double log_σ = main_rates.getLogTransientParameter(v);
				
					double log_length_mul;
					if (logit_λ == Double.NEGATIVE_INFINITY)
						log_length_mul = mod_length; // Poisson
					else 
					{
						log_length_mul = mod_length 
								+ Logarithms.logitToLogComplement(cat_logit_λ)
								- main_rates.getLogRelativeComplement(v);
					}
					double cat_log_σ = log_σ * Math.exp(log_length_mul);
					double cat_logit_σ = Logarithms.logToLogit(cat_log_σ);
					
					if (logit_λ == Double.NEGATIVE_INFINITY) // Poisson
					{
						double cat_logit_p = -cat_logit_σ;
						this.setLogitLossRelativeDuplication(v, cat_logit_p, logit_λ);
					} else // Polya 
					{
						double cat_logit_p = -cat_logit_σ-Logarithms.logitToLogComplement(cat_logit_λ);
						this.setLogitLossRelativeDuplication(v, cat_logit_p, cat_logit_λ);
					}
				}
			}

		}
	}

	/**
	 * Iterator over the rate categories.
	 */
	@Override
	public Iterator<Category> iterator() 
	{
		return Collections.unmodifiableList(rate_categories).iterator();
	}
	
}
