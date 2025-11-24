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

import static count.io.CommandLine.OPT_EPS;
import static count.io.CommandLine.OPT_GAIN;
import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_MODEL_CATEGORIES;
import static count.io.CommandLine.OPT_MODEL_DUPLICATION_CATEGORIES;
import static count.io.CommandLine.OPT_MODEL_LENGTH_CATEGORIES;
import static count.io.CommandLine.OPT_MODEL_ROOT_PRIOR;
import static count.io.CommandLine.OPT_OUTPUT;
import static count.io.CommandLine.OPT_PVALUE;
import static count.io.CommandLine.OPT_ROUNDS;
import static count.io.CommandLine.OPT_TRUNCATE;
import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LOSS;
import static count.model.VariationGradientFactory.PARAMETER_MOD_DUPLICATION;
import static count.model.VariationGradientFactory.PARAMETER_MOD_LENGTH;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.function.DoubleFunction;
import java.util.function.Function;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.ds.UniqueProfileTable;
import count.io.CommandLine;
import count.io.NewickParser;
import count.io.RateVariationParser;
import count.matek.FunctionMinimization;
import count.matek.FunctionMinimization.CGVariant;
import count.matek.Functions;
import count.matek.Logarithms;


/**
 * 
 * Numerical optimization for {@link count.model.RateVariationModel}. 
 * Current favorite for BFGS. 
 */
public class MLRateVariation extends ML
{
	/**
	 * Verbosity for optimization steps; set to true if launched from command line 
	 */
	public static boolean PRINT_OPTIMIZATION_MESSAGES = false;
	

	/**
	 * The model to be optimized
	 */
	protected final RateVariationModel variation_model;

	/**
	 * Sample profiles (uniques with multiplicity)
	 */
	protected final UniqueProfileTable utable;

	/**
	 * Computations of log-likelihood, posteriors and gradient across the sample
	 */
	protected final VariationGradientFactory gradient_factory;
	
	
	/**
	 * Initialized for node parameters, then refitted by {@link #initDataStructures()}
	 */
	protected boolean[] do_optimize_parameters;
	
	/**
	 * Array for storing node parameter values
	 * (logit-p, logit-lambda, log-gain) ; must be synchronized with model 
	 */
	protected final double[][] node_parameters;
	
	/**
	 * Array for storing category parameter values ; must be synchronized with model 
	 */
	protected double[][] category_parameters;

	/**
	 * List of all parameter instances 
	 */
	protected final List<MinGradient> full_parameters;
		
	/**
	 * List of adjustable/optimizable parameter instances 
	 */
	protected final List<MinGradient> adjustable_parameters;
	
	
	public MLRateVariation(RateVariationModel model, ProfileTable table)
	{
		if (table instanceof UniqueProfileTable)
			utable = (UniqueProfileTable) table;
		else
			utable = new UniqueProfileTable(table);
		
		this.variation_model = model;
		this.gradient_factory = new VariationGradientFactory(model, utable);
		//this.lrates = model.getBaseModel();
		TreeWithLogisticParameters lrates = model.getBaseModel();
		this.do_optimize_parameters = ML.optimizableParams(lrates, false); // want to optimize even p==q model
		int num_nodes = lrates.getTree().getNumNodes();
		this.node_parameters = new double[num_nodes][];
		this.adjustable_parameters = new ArrayList<>();
		this.full_parameters = new ArrayList<>();
		this.initParameters();
	}
	
	
	@Override
	public void setMinimumObservedCopies(int m) 
	{
		gradient_factory.setMinimumObservedCopies(m);
	}
	
	@Override
	public void setCalculationWidth(int absolute, double relative) 
	{
		gradient_factory.setCalculationWidthThresholds(absolute, relative);
		gradient_factory.computeClasses();
	}
	
	private static CGVariant USE_CONJUGATE_GRADIENT = null; // CGVariant.DY_HS; // CGVariant.DaiYuan;//CGVariant.HestenesStiefel;//  CGVariant.FR_PR; //   CGVariant.PolakRibiere;
	
	private static final double MAX_KAPPA = Double.POSITIVE_INFINITY; // 128.0; //99.0*Math.log(2.0);//   33.0; //Double.POSITIVE_INFINITY;
	private double max_modifier = 33.0;
	private static int DEFAULT_TRUNCATE_ABSOLUTE = 6;
	private static double DEFAULT_TRUNCATE_RELATIVE = 1.0;
	private boolean auto_truncate = true;
	private static boolean DEBUG_GRADIENT = false; 
	
	/**
	 * Unbounded loss rate parameter (using log) is not numerically stable; so 
	 * we bound it and optimize it as logit(scaled gain).
	 */
	private final static boolean ALWAYS_UNBOUNDED_GAIN_LOSS = false; // keep false - do not change	
	private static double MAX_DUPLICATION_MARGIN = 0.0; //do not change
	private static boolean USE_LOGISTIC_GAIN = true;
	private static final boolean USE_BRACKETS = false;
	private static final boolean REGULARIZE_GAIN = true;
	
	
//	public void setGainParameterBound(double max_gain)
//	{
//		this.max_gain = max_gain;
//		this.initParameters();
//	}
//	
//	public void setModifierBound(double max_mod)
//	{
//		this.max_modifier = max_mod;
//		this.initParameters();
//	}
	
	public void setWantAutoTruncation(boolean auto_truncate)
	{
		this.auto_truncate = auto_truncate;
		if (auto_truncate)
		{
			this.setCalculationWidth(DEFAULT_TRUNCATE_ABSOLUTE, DEFAULT_TRUNCATE_RELATIVE);
		}
	}
	
	public int getCommonGainType() {
		return variation_model.getCommonGainType();
	}
	
	public boolean isDuplicationBounded() {return variation_model.isDuplicationBounded();}
	
	public void setDuplicationBounded(boolean isBounded) {
		if (variation_model.isDuplicationBounded() != isBounded) {
			variation_model.setBoundedDuplication(isBounded);
		}
		this.initParameters();
	}
	
	public boolean isGainBounded() {return USE_LOGISTIC_GAIN;}
	
	public void setGainBounded(boolean isBounded) {
		USE_LOGISTIC_GAIN = isBounded;
		this.initParameters();
	}
	
	
	
	
	/*
	 * 
	 * Synchronization between locally stored parameters and model parameters 
	 * 
	 */
	
	/**
	 * Called from instantiation, and must be called if the 
	 * underlying {@link #variation_model} structure changes.
	 * Checks the ranges of node parameters, populates 
	 * the {@link #full_parameters} list and {@link #do_optimize_parameters} array, finally resets the 
	 * gradient factory. Calls {@link #initOptimizableCategoryParameters()} and {@link #initCategoryParameters()} 
	 * to add category parameters.
	 * 
	 */
	protected void initParameters()
	{
		int num_nodes = node_parameters.length;
//		// DEBUG
//		System.out.println("#**MLRV.iP n "+num_nodes);
		
		TreeWithLogisticParameters lrates = variation_model.getBaseModel();
		int common_gain = variation_model.getCommonGainType();

		boolean have_adjusted_parameters = false;
		boolean bounded_duplication = variation_model.isDuplicationBounded();
		
		for (int v=0; v<num_nodes; v++)
		{
			double logit_p = lrates.getLogitLossParameter(v);
			double log_gain = lrates.getLogGainParameter(v, common_gain, !variation_model.isUniversalGain());
			
			
			if (bounded_duplication) {
				// check duplication rate
				double logit_q = lrates.getLogitDuplicationParameter(v);
				if (logit_p<logit_q) { // thus q is not 0
					double new_logit_q = Logarithms.mulLogit(logit_p, 0.0); // q=p/2
					double log_q = lrates.getLogDuplicationParameter(v);
					double new_log_q = Logarithms.logitToLogValue(new_logit_q);
					if (common_gain == PARAMETER_DUPLICATION || common_gain == PARAMETER_GAIN) {
						double new_log_gain = log_gain + log_q - new_log_q;
						if (PRINT_OPTIMIZATION_MESSAGES) {
							if (common_gain==PARAMETER_DUPLICATION)
								System.out.println("#*MLRV.iP "+v+" reset q log="+log_q+"\tto "+new_log_q+"\tkappa log="+log_gain+"\tto "+new_log_gain+"\t// "+lrates.toString(v));
							else
								System.out.println("#*MLRV.iP "+v+" reset q log="+log_q+"\tto "+new_log_q+"\tr log="+log_gain+"\tto "+new_log_gain+"\t// "+lrates.toString(v));
						}
						lrates.setLogitLossDuplication(v, logit_p, new_logit_q, new_log_gain, common_gain, !variation_model.isUniversalGain());
						log_gain = new_log_gain;
					} else {
						if (PRINT_OPTIMIZATION_MESSAGES) 
							System.out.println("#*MLRV.iP "+v+" reset q log="+log_q+"\tto "+new_log_q+"\tgamma log="+log_gain+"\tstays\t// "+lrates.toString(v));
						lrates.setLogitLossDuplication(v, logit_p, new_logit_q, log_gain, common_gain, !variation_model.isUniversalGain());
					}
					logit_q = new_logit_q;
					have_adjusted_parameters = true;
				}
				double logit_lambda = lrates.getLogitRelativeRate(v);
				// check large gain 
				double log_max_gain = Math.log(
						common_gain == PARAMETER_DUPLICATION?
						MAX_KAPPA:1.0);
				if (log_max_gain<log_gain)
				{
					double new_log_gain = log_max_gain + Math.log(1023.0/1024.0);
					if (common_gain == PARAMETER_LOSS) 
					{
						if (!USE_LOGISTIC_GAIN && ALWAYS_UNBOUNDED_GAIN_LOSS)
						{
							// nothing to do: max_gain is ignored
						} else
						{
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLRV.iP "+v+" reset gamma log="+log_gain+"\tto "+new_log_gain+"\t// "+lrates.toString(v));
							lrates.setLogitLossRelativeDuplication(v, logit_p, logit_lambda, new_log_gain, common_gain, !variation_model.isUniversalGain());
							log_gain = new_log_gain;
							have_adjusted_parameters = true;
						}
					} else if (common_gain == PARAMETER_DUPLICATION && logit_lambda != Double.NEGATIVE_INFINITY)
					{
						// try to adjust lambda so that we have same r=q*kappa and same gamma=lambda*kappa
						double log_lambda = lrates.getLogRelativeRate(v);
						double new_log_lambda = log_lambda + log_gain - new_log_gain;
						if (new_log_lambda < 0.0)
						{
							logit_lambda = Logarithms.logToLogit(new_log_lambda);
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLRV.iP "+v+" reset kappa log="+log_gain+"\tto "+new_log_gain+"\tlambda log="+log_lambda+"\tto "+new_log_lambda+"\t// "+lrates.toString(v));
							lrates.setLogitLossRelativeDuplication(v, logit_p, logit_lambda, new_log_gain, common_gain, !variation_model.isUniversalGain());
						} else
						{
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLRV.iP "+v+" reset kappa log="+log_gain+"\tto "+new_log_gain+"\tlambda stays log="+log_lambda+"\t// "+lrates.toString(v));
							lrates.setLogitLossRelativeDuplication(v, logit_p, logit_lambda, new_log_gain, common_gain, !variation_model.isUniversalGain());						
						}
						log_gain = new_log_gain;
						have_adjusted_parameters = true;
					} else
					{
						double log_lambda = lrates.getLogRelativeRate(v);
						if (common_gain == PARAMETER_GAIN || (common_gain==PARAMETER_DUPLICATION && logit_lambda==Double.NEGATIVE_INFINITY))
						{
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLRV.iP "+v+" reset r log="+log_gain+"\tto "+new_log_gain+"\tlambda stays log="+log_lambda+"\t// "+lrates.toString(v));
						} else 
						{
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLRV.iP "+v+" reset kappa log="+log_gain+"\tto "+new_log_gain+"\tlambda stays log="+log_lambda+"\t// "+lrates.toString(v));
						}
						lrates.setLogitLossRelativeDuplication(v, logit_p, logit_lambda, new_log_gain, common_gain, !variation_model.isUniversalGain());						
						log_gain = new_log_gain;
						have_adjusted_parameters = true;
					}
				}
				
				// check lambda == 1.0
				double reset_lambda;
				double small = 1.0/(1L<<20); // this will maybe decrease the likelihood
				if (logit_lambda == Double.POSITIVE_INFINITY)
				{
					reset_lambda = -Logarithms.logToLogit(Math.log(small));
				} else
					reset_lambda = logit_lambda;
				if (reset_lambda != logit_lambda)
				{
					String old_node_str = lrates.toString(v);
					lrates.setLogitLossRelativeDuplication(v, logit_p, reset_lambda, log_gain, common_gain, !variation_model.isUniversalGain());
					have_adjusted_parameters=true;
					
					// DEBUG
					if (PRINT_OPTIMIZATION_MESSAGES)
						System.out.println("#*MLRV.iP "+v+" reset lambda "+logit_lambda+"\tto "+reset_lambda+"\t// "+lrates.toString(v)+"\t// was "+old_node_str);
					
					logit_lambda = reset_lambda;
				}			
				
				boolean first_init = node_parameters[v]==null;
				
				double[] params = new double[3];
				params[PARAMETER_LOSS] 		  = logit_p;
				params[PARAMETER_DUPLICATION] = logit_lambda;
				
				if (USE_LOGISTIC_GAIN)
					params[PARAMETER_GAIN] = Logarithms.logToLogit(log_gain-log_max_gain);
				else
					params[PARAMETER_GAIN] = log_gain;
				
				this.node_parameters[v] = params;
//				// DEBUG
//				System.out.println("#**MLRV.iP "+v+"\t"+Arrays.toString(params)+"\t// "+lrates.toString(v));
				
				
				if (first_init)
				{
					fixLoss(v, logit_p == Double.POSITIVE_INFINITY); // || lrates.getTree().isRoot(v)
					fixDuplication(v, logit_lambda == Double.NEGATIVE_INFINITY);
					fixGain(v, log_gain == Double.NEGATIVE_INFINITY);
				}
			} else {
				double logit_q = lrates.getLogitDuplicationParameter(v);
				// check large gain 
				double log_max_gain = Math.log(MAX_KAPPA);
				if (log_max_gain<log_gain)
				{
					double new_log_gain = log_max_gain + Math.log(1023.0/1024.0);
					if (common_gain == PARAMETER_LOSS) 
					{
						if (ALWAYS_UNBOUNDED_GAIN_LOSS)
						{
							// nothing to do: max_gain is ignored
						} else
						{
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLRV.iP "+v+" reset gamma log="+log_gain+"\tto "+new_log_gain+"\t// "+lrates.toString(v));
							lrates.setLogitLossDuplication(v, logit_p, logit_q, new_log_gain, common_gain, !variation_model.isUniversalGain());
							log_gain = new_log_gain;
							have_adjusted_parameters = true;
						}
					} else if (common_gain == PARAMETER_DUPLICATION && logit_q != Double.NEGATIVE_INFINITY)
					{
						double log_q = lrates.getLogDuplicationParameter(v);
						double new_log_q = log_q + log_gain-new_log_gain;
						
						// q needs to grow to have same r=q*kappa and same gamma=lambda*kappa		
						if (new_log_q<0.0) {
							logit_q = Logarithms.logToLogit(new_log_q);
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLRV.iP "+v+" reset kappa log="+log_gain+"\tto "+new_log_gain+"\tq log="+log_q+"\tto "+new_log_q+"\t// "+lrates.toString(v));
							lrates.setLogitLossDuplication(v, logit_p, logit_q, new_log_gain, common_gain, !variation_model.isUniversalGain());
						} else
						{
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLRV.iP "+v+" reset kappa log="+log_gain+"\tto "+new_log_gain+"\tq stays log="+log_q+"\t// "+lrates.toString(v));
							lrates.setLogitLossDuplication(v, logit_p, logit_q, new_log_gain, common_gain, !variation_model.isUniversalGain());						
						}
						log_gain = new_log_gain;
						have_adjusted_parameters = true;
					} else
					{
						double log_q = lrates.getLogDuplicationParameter(v);
						if (common_gain == PARAMETER_GAIN || (common_gain==PARAMETER_DUPLICATION && logit_q==Double.NEGATIVE_INFINITY))
						{
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLRV.iP "+v+" reset r log="+log_gain+"\tto "+new_log_gain+"\tq stays log="+log_q+"\t// "+lrates.toString(v));
						} else 
						{
							if (PRINT_OPTIMIZATION_MESSAGES)
								System.out.println("#*MLRV.iP "+v+" reset kappa log="+log_gain+"\tto "+new_log_gain+"\tq stays log="+log_q+"\t// "+lrates.toString(v));
						}
						lrates.setLogitLossDuplication(v, logit_p, logit_q, new_log_gain, common_gain, !variation_model.isUniversalGain());						
						log_gain = new_log_gain;
						have_adjusted_parameters = true;
					}
				}
				boolean first_init = node_parameters[v]==null;
				
				double[] params = new double[3];
				params[PARAMETER_LOSS] 		  = logit_p;
				params[PARAMETER_DUPLICATION] = logit_q;
				if (USE_LOGISTIC_GAIN)
					params[PARAMETER_GAIN] = Logarithms.logToLogit(log_gain-log_max_gain);
				else
					params[PARAMETER_GAIN] = log_gain;
				
				this.node_parameters[v] = params;
					
				if (first_init)
				{
					fixLoss(v, logit_p == Double.POSITIVE_INFINITY); // lrates.getTree().isRoot(v) || 
					fixDuplication(v, logit_q == Double.NEGATIVE_INFINITY);
					fixGain(v, log_gain == Double.NEGATIVE_INFINITY);
				}
			}
		} 

		this.full_parameters.clear();
		for (int v=0; v<num_nodes; v++)
		{
			full_parameters.add(defaultNodeParameter(v,PARAMETER_GAIN));
			full_parameters.add(defaultNodeParameter(v,PARAMETER_DUPLICATION));
			full_parameters.add(defaultNodeParameter(v,PARAMETER_LOSS));
		}
		this.initCategoryParameters();
		this.initOptimizableCategoryParameters();
		gradient_factory.computeClasses();
	}
	
	/**
	 * Adds category parameters to the {@link #full_parameters} (2 parameters) and {@link #category_parameters}
	 * list. 
	 * 
	 */
	protected void initCategoryParameters() {
		int ncat = variation_model.getNumClasses();
		this.category_parameters = new double[ncat][];
		for (int k=0; k<ncat; k++)
		{
			RateVariationModel.Category C = variation_model.getCategory(k);
			double[] mods = new double[2];
			mods[PARAMETER_MOD_LENGTH] = C.getModLength();
			mods[PARAMETER_MOD_DUPLICATION] = C.getModDuplication();
			this.category_parameters[k] = mods;
		}
		
		for (int k=0; k<ncat; k++)
		{
			full_parameters.add(defaultCategoryParameter(k, PARAMETER_MOD_LENGTH));
			full_parameters.add(defaultCategoryParameter(k, PARAMETER_MOD_DUPLICATION));
		}
		//return 2*ncat;
	}	

//	private boolean initNodeParameters(int v)
//	{
//		TreeWithLogisticParameters lrates = variation_model.getBaseModel();
//		int common_gain = variation_model.getCommonGainType();
//
//		double logit_p = lrates.getLogitLossParameter(v);
//		double logit_lambda = lrates.getLogitRelativeRate(v);
//		double log_gain = lrates.getLogGainParameter(v, common_gain, !variation_model.isUniversalGain());
//		
//		boolean have_adjusted_parameters = false;
//		
//		
//		
//		
//		
//		return have_adjusted_parameters;
//	}
	
	
	/**
	 * Sets initial values for in {@link #do_optimize_parameters} 
	 * for category mods. The first 3n cells are 
	 * used for the parameters of n nodes.  
	 */
	protected void initOptimizableCategoryParameters()
	{
		int num_nodes = node_parameters.length;
		int ncat = variation_model.getNumClasses();
		this.do_optimize_parameters = Arrays.copyOf(do_optimize_parameters, 3*num_nodes+2*ncat);

		for (int k=0; k<ncat; k++)
		{
			RateVariationModel.Category C = variation_model.getCategory(k);
			
			int j_mod_len = 3*num_nodes + 2*k + PARAMETER_MOD_LENGTH;
			int j_mod_dup = 3*num_nodes + 2*k + PARAMETER_MOD_DUPLICATION;
			if ( C.getLogCatProbability()==Double.NEGATIVE_INFINITY)
			{
				do_optimize_parameters[j_mod_len] =
				do_optimize_parameters[j_mod_dup] = false;		
			} else
			{
				do_optimize_parameters[j_mod_len] = C.getModLength()!=0.0
						&& Double.isFinite(C.getModLength());
				do_optimize_parameters[j_mod_dup] = C.getModDuplication()!=0.0
						&& C.getModDuplication() != Double.NEGATIVE_INFINITY;
			}
		}
	}
	
	public void fixGain(int node, boolean not_optimized)
	{
		do_optimize_parameters[3*node+PARAMETER_GAIN]=!not_optimized;
	}
	public void fixLoss(int node, boolean not_optimized)
	{
		do_optimize_parameters[3*node+PARAMETER_LOSS]=!not_optimized;
	}
	public void fixDuplication(int node, boolean not_optimized)
	{
		do_optimize_parameters[3*node+PARAMETER_DUPLICATION]=!not_optimized;
		
	}	

	@Override
	public void fixNodeParameters(int node, boolean do_not_optimize)
	{
		fixGain(node, do_not_optimize || node_parameters[node][PARAMETER_GAIN]==Double.NEGATIVE_INFINITY);
		
		fixLoss(node, do_not_optimize || node_parameters[node][PARAMETER_LOSS]==Double.POSITIVE_INFINITY);
		fixDuplication(node, do_not_optimize || node_parameters[node][PARAMETER_DUPLICATION]==Double.NEGATIVE_INFINITY);
	}
	
	/**
	 * Fills in 
	 * the adjustable parameters (in {@link #adjustable_parameters}.
	 * Calls {@link #copyCategoryParametersFromModel()} for 
	 * adding the adjustable category parameters. 
	 * 	
	 */
	private void copyParametersFromModel()
	{
		gradient_factory.computeClasses();
		
		int num_nodes = node_parameters.length;
		
		// now collect adjustable parameters
		adjustable_parameters.clear();
		for (int v=0; v<num_nodes; v++)
		{
			double[] params = node_parameters[v];
			if (params != null)
			{
				boolean want_gain = this.do_optimize_parameters[3*v+PARAMETER_GAIN];
				boolean want_loss = this.do_optimize_parameters[3*v+PARAMETER_LOSS];
				boolean want_dup = this.do_optimize_parameters[3*v+PARAMETER_DUPLICATION];
				
				if (want_gain)
				{
					MinGradient P = defaultNodeParameter(v, PARAMETER_GAIN);
					if (!Double.isInfinite(P.get()))
						adjustable_parameters.add(P);
				}
				MinGradient dup_par = defaultNodeParameter(v, PARAMETER_DUPLICATION);  // new LogisticDuplicationRate(v);
				if (want_dup)
				{
					double logit_lm = dup_par.get();
					
					if (!Double.isInfinite(logit_lm)) {
						adjustable_parameters.add(dup_par);
						
//						if (v==num_nodes-1) { // DEBUG
//							System.out.println("#**MLRVL.cPFM node "+v+"\tadd dup "+dup_par);
//						}
					} else {
//						if (v==num_nodes-1) { // DEBUG
//							System.out.println("#**MLRVL.cPFM node "+v+"\tnoadd dup "+dup_par);
//						}
					}
				} else {
//					if (v==num_nodes-1) { // DEBUG
//						System.out.println("#**MLRVL.cPFM node "+v+"\tnoadd dup "+dup_par);
//					}
				}
				
				
				if (want_loss)
				{
					MinGradient P = defaultNodeParameter(v, PARAMETER_LOSS);
					if (!Double.isInfinite(P.get()))
						adjustable_parameters.add(P);
				}
			}
		}
		this.copyCategoryParametersFromModel();
	}
		
	/**
	 * 	Adds adjustable category parameters.
	 */
	protected void copyCategoryParametersFromModel() {	
		int num_nodes = node_parameters.length;
		int ncat = variation_model.getNumClasses(); 
		for (int k=0; k<ncat; k++)
		{
			double[] mods = category_parameters[k];
			if (mods!=null)
			{
				boolean want_len = this.do_optimize_parameters[3*num_nodes+2*k+PARAMETER_MOD_LENGTH];
				boolean want_dup = this.do_optimize_parameters[3*num_nodes+2*k+PARAMETER_MOD_DUPLICATION];
				
				if (want_len)
				{
					MinGradient M = defaultCategoryParameter(k, PARAMETER_MOD_LENGTH);
					adjustable_parameters.add(M);
				}
				if (want_dup)
				{
					MinGradient M = defaultCategoryParameter(k, PARAMETER_MOD_DUPLICATION);
					adjustable_parameters.add(M);
				}
			}
		}
	}
	
	/**
	 * Copies the locally stored parameter values to 
	 * the underlying {@link #variation_model}.
	 * Calls {@link #copyCategoryParametersToModel()} 
	 * for setting category parameters, and 
	 * resets the underlying gradient factory. 
	 */
	protected final void copyParametersToModel()
	{
		TreeWithLogisticParameters lrates = variation_model.getBaseModel();
		int num_nodes = node_parameters.length;
		boolean bounded_duplication = variation_model.isDuplicationBounded();
		for (int v=0; v<num_nodes; v++)
		{
			double[] params = node_parameters[v];
			if (params!=null)
			{
				double logit_p =  params[PARAMETER_LOSS];
				double log_gain;
				if (USE_LOGISTIC_GAIN) {
					double log_max_gain = Math.log(
							getCommonGainType() == PARAMETER_DUPLICATION?
									MAX_KAPPA:1.0);
					if (log_max_gain == Double.POSITIVE_INFINITY)
						log_gain = params[PARAMETER_GAIN];
					else {
						double logit_gn = params[PARAMETER_GAIN];
						log_gain = Logarithms.logitToLogValue(logit_gn)
								+log_max_gain;
//						// DEBUG
//						if ((log_gain==Double.NEGATIVE_INFINITY || log_gain == log_max_gain)) { //  && !Double.isInfinite(logit_gn)
//							System.out.println("#**MLRV.cPTM node "+v+"/gain\tlogit "+logit_gn+"\tlog "+log_gain+"\tlmax "+log_max_gain+"+("+Logarithms.logitToLogValue(logit_gn)+")"
//									+"\tlogitp "+logit_p
//									+(bounded_duplication?"\tlogitq ":"\tlogitlm ")
//										+params[PARAMETER_DUPLICATION]
//									);
//						}
					}
				} else 
					log_gain = params[PARAMETER_GAIN];
				
				if (bounded_duplication) {
					double logit_lambda = params[PARAMETER_DUPLICATION];
					lrates.setLogitLossRelativeDuplication(v, logit_p, logit_lambda, log_gain, variation_model.getCommonGainType());
				}
				else {
					double logit_q = params[PARAMETER_DUPLICATION];
					lrates.setLogitLossDuplication(v, logit_p, logit_q, log_gain, variation_model.getCommonGainType());
				}
			}
		}
		this.copyCategoryParametersToModel();
		gradient_factory.computeClasses();
	}
	
	protected void copyCategoryParametersToModel() {
		int ncat = category_parameters.length;
		for (int k=0; k<ncat; k++)
		{
			double[] params = category_parameters[k];
			RateVariationModel.Category C = variation_model.getCategory(k);
//			if (params != null)
//			{
				C.setModifiers(params[PARAMETER_MOD_LENGTH], params[PARAMETER_MOD_DUPLICATION]);
				// parameters automatically recomputed
//			} else
//				C.computeParameters();
		}
	}	
	
	
	/*
	 * 
	 * Model adjustments 
	 * 
	 */
	/**
	 * Replaces 1.0 duplication rates to facilitate optimization
	 * in the base rates . 
	 * 
	 * The same functionality is assumed by {@link #initParameters()}.
	 * 
	 * @return true if at least one lambda was reset
	 */
	private boolean reduceInfiniteDuplicationRates(double small)
	{
		boolean have_adjusted = false;
		TreeWithLogisticParameters lrates = variation_model.getBaseModel();
		IndexedTree tree = lrates.getTree();
		int num_nodes = tree.getNumNodes();
		for (int v=0; v<num_nodes; v++)
		{
			if (do_optimize_parameters[3*v+PARAMETER_DUPLICATION])
			{
				double logit_lambda = lrates.getLogitRelativeRate(v);
				
				double reset_lambda;
				
				
//					if (logit_lambda == Double.NEGATIVE_INFINITY)
//					{
//						reset_lambda = Logarithms.logToLogit(Math.log(small));
//					} else 
				if (logit_lambda == Double.POSITIVE_INFINITY)
				{
					reset_lambda = -Logarithms.logToLogit(Math.log(small));
				} else
					reset_lambda = logit_lambda;
				if (reset_lambda != logit_lambda)
				{
					double logit_p = lrates.getLogitLossParameter(v);
					double log_gain = lrates.getLogGainParameter(v, variation_model.getCommonGainType(), !variation_model.isUniversalGain());
					
					String old_node_str = lrates.toString(v);
					lrates.setLogitLossRelativeDuplication(v, logit_p, reset_lambda, log_gain, variation_model.getCommonGainType(), !variation_model.isUniversalGain());
//					this.fixDuplication(v, false);
					
					
					
					have_adjusted=true;
					
					// DEBUG
					if (PRINT_OPTIMIZATION_MESSAGES)
						System.out.println("#*MLRV.rIDR "+v+" lambda resets "+reset_lambda+"\t// "+lrates.toString(v)+"\t// was "+old_node_str);
				}
			}
			
		}
		
		if (have_adjusted) // bug fix 2024/07/26: need to update if rates change 
		{
			gradient_factory.computeClasses();
			
		}
		
		return have_adjusted;
	}

	/**
	 * Shifts the main category to the one with maximum prior probability
	 * 
	 * @return true if model structure changed (still same log-likelihood, but the gradients are different as parameters change)
	 */
	protected boolean recenterBaseModel()
	{
		boolean recenter = variation_model.recenterCategories();
		if (recenter)
		{
			int ncat = variation_model.getNumClasses();
			int nnodes = node_parameters.length;
			for (int k=0; k<ncat; k++)
			{
				RateVariationModel.Category C = variation_model.getCategory(k);
				System.out.println("#*MLRV.rBM before cat "+k+" "+C.toString()
					+"(stored "
						+full_parameters.get(3*nnodes+2*k)+"/"+do_optimize_parameters[3*nnodes+2*k]
					+","+full_parameters.get(3*nnodes+2*k+1)+"/"+do_optimize_parameters[3*nnodes+2*k+1]
					+")"					
					);
			}
			
			this.initParameters();
//			this.initOptimizableCategoryParameters();
			this.copyParametersFromModel();
			
			
			for (int k=0; k<ncat; k++)
			{
				RateVariationModel.Category C = variation_model.getCategory(k);
				System.out.println("#*MLRV.rBM after cat "+k+" "+C.toString()
					+"(stored "
						+full_parameters.get(3*nnodes+2*k)+"/"+do_optimize_parameters[3*nnodes+2*k]
					+","+full_parameters.get(3*nnodes+2*k+1)+"/"+do_optimize_parameters[3*nnodes+2*k+1]
					+")"					
					);
			}
		}
		return recenter;
	}
	
	
	
	/*
	 * 
	 * Saving and restoring optimization states 
	 * 
	 */
	/**
	 * Stores a complete model setting with associated likelihood, posteriors and gradient 
	 * 
	 */
	protected class OptimizationState 
	{
		final double[] allx;
		final VariationGradientFactory.SampleGradient G;
		final double[] log_prior;
		final long computing_time_nano;
		
		OptimizationState()
		{
			this.allx = ML.getParameterValues(MLRateVariation.this.full_parameters);
			long t0 = System.nanoTime();
			boolean want_dup = !variation_model.isDuplicationBounded();
			this.G = gradient_factory.computeSampleGradient(want_dup); // we want by lambda and not by dup parameter q 
			G.correctForUnobserved(want_dup);
			this.computing_time_nano = System.nanoTime()-t0;

			int ncat = variation_model.getNumClasses();
			this.log_prior = new double[ncat];
			for (int k=0; k<ncat; k++)
			{
				log_prior[k] = variation_model.getCategory(k).getLogCatProbability();
			}
		}
		
		void resetModel()
		{
			ML.setParameterValues(allx, MLRateVariation.this.full_parameters);
			int ncat = variation_model.getNumClasses();
			for (int k=0; k<ncat; k++)
			{
				variation_model.getCategory(k).setLogCatProbability(log_prior[k]);
			}
			MLRateVariation.this.initOptimizableCategoryParameters(); // maybe mods==0.0 changed 
			MLRateVariation.this.copyParametersToModel();
		}
		
		/**
		 * Max normalized change of parameter values, 
		 * including category prior probabilities
		 * 
		 * @param that
		 * @return
		 */
		double maxRelDelta(OptimizationState that)
		{
			double maxRelDelta = 0.0;
			assert (this.allx.length == that.allx.length);
			
			for (int j=0; j<this.allx.length;j++)
			{
				if (this.allx[j]!=that.allx[j])
				{
					double dx;
					if (Double.isInfinite(this.allx[j]))
					{
						dx = Double.isInfinite(that.allx[j])?2.0:1.0;
						// if both infinite, then with the opposite signs
					} else if (Double.isInfinite(that.allx[j]))
					{
						dx = 1.0;
					} else
					{
						// both finite
						dx = Math.abs(this.allx[j]-that.allx[j])/Double.max(1.0, Math.abs(this.allx[j]));
					}
					maxRelDelta = Double.max(dx, maxRelDelta);
				}
			}
			
			assert (this.log_prior.length == that.log_prior.length);
			
			for (int k=0; k<log_prior.length; k++)
			{
				double thisp = Math.exp(this.log_prior[k]);
				double thatp = Math.exp(that.log_prior[k]);
				double dx =  Math.abs(thisp-thatp);
				maxRelDelta = Double.max(dx, maxRelDelta);
			}
			return maxRelDelta;
		}
		
		/**
		 * Maximum absolute displacement (L-infinity norm for parameter change vector)
		 * 
		 * @param that
		 * @return
		 */
		double maxDelta(OptimizationState that) {
			double maxDelta = 0.0;
			assert (this.allx.length == that.allx.length);
			
			for (int j=0; j<this.allx.length;j++)
			{
				if (this.allx[j]!=that.allx[j])
				{
					double dx;
					if (Double.isFinite(this.allx[j]) && Double.isFinite(that.allx[j])){
						maxDelta = Double.max(maxDelta, Math.abs(this.allx[j]-that.allx[j]));
					}
				}
			}
			assert (this.log_prior.length == that.log_prior.length);
			
			for (int k=0; k<log_prior.length; k++)
			{
				double thisp = Math.exp(this.log_prior[k]);
				double thatp = Math.exp(that.log_prior[k]);
				maxDelta = Double.max(maxDelta, Math.abs(thisp-thatp));
			}
			return maxDelta;
			
		}
		
		
		double getGradientL2(boolean want_node_params, boolean want_category_params)
		{
			double[] D = getAdjustableParameterGradient(G.get(), want_node_params, want_category_params);
			return  FunctionMinimization.euclideanNorm(D);
		}
		
		
		double getGradientDiffL2(OptimizationState that, boolean want_node_params, boolean want_category_params)
		{
			double[] thisD = getAdjustableParameterGradient(this.G.get(), want_node_params, want_category_params);
			double[] thatD = getAdjustableParameterGradient(that.G.get(), want_node_params, want_category_params);
			double getGradientDiffL2 = FunctionMinimization.euclideanDistance(thisD, thatD);
			return getGradientDiffL2;
		}
		
		
		
		double getGradientLmax(boolean want_node_params, boolean want_category_params)
		{
			double[] D = getAdjustableParameterGradient(G.get(), want_node_params, want_category_params);
			return  FunctionMinimization.maxNorm(D);
		}
		
		
		double getLL()
		{
			return G.getLogLikelihood();
		}
		
		double[] getGradient() {
			return G.get();
		}
		
		long nanoTime()
		{
			return this.computing_time_nano;
		}
		
		public boolean equals (OptimizationState that) {
			return Arrays.equals(this.allx, that.allx)
					&& Arrays.equals(this.log_prior, that.log_prior);
		}
		
		@Override 
		public boolean equals(Object o) {
			if (o instanceof OptimizationState) {
				OptimizationState that = (OptimizationState) o;
				return this.equals(that);
			} else return false;
		}
		
		@Override 
		public int hashCode() {
			return (Arrays.hashCode(this.allx)*37)^Arrays.hashCode(this.log_prior);
		}
	
		
		@Override
		public String toString()
		{
			StringBuilder sb = new StringBuilder("{");
			double LL = G.getLogLikelihood();
			sb.append("LL ").append(LL);
			sb.append("; pcat {");
			int ncat = variation_model.getNumClasses();
			for (int k=0; k<ncat; k++)
			{
				if (0<k)
					sb.append(", ");
				sb.append(Math.exp(log_prior[k]));
			}
			sb.append("}");
			double Dlen = this.getGradientL2(true,true);
			double rlen = Dlen/(-LL);
			
			sb.append("; |dL| ").append(Dlen)
				.append("/").append(rlen);
			
			sb.append("}");
			return sb.toString();
		}
		
		
	}

	private OptimizationState currentState() { return new OptimizationState();}	
	
	
	private OptimizationState adjustCalculationWidth(OptimizationState current, double tol)
	{
		
		int absolute = gradient_factory.getCalculationWidthAbsolute();
		double relative = gradient_factory.getCalculationWidthRelative();
		if (absolute == Integer.MAX_VALUE || relative == Double.POSITIVE_INFINITY)
			return current;
		
		firePropertyChange​(PROPERTY_OPTIMIZATION_PHASE, "Adjusting truncation parameters");

		double current_LL = current.getLL();
		double current_dL = current.getGradientL2(true, false);
		
		// need to increase?
		double step_size = Math.log(2.0)/3.0;
		
		
		
		
		int dabs =0, drel = 0;
		
		double ftol = tol; // tolerance on function value
		double dtol = tol; //Math.sqrt(tol); // tolerance on derivative
		
		// computing the true values 
		this.setCalculationWidth(Integer.MAX_VALUE, Double.POSITIVE_INFINITY);
		OptimizationState trueS = currentState();
		double true_LL = trueS.getLL();
		double true_dL = trueS.getGradientL2(true, false);
		final double dLscale = Double.max(1.0,true_dL);
		double current_delta = current_LL-true_LL;
		double current_rdiff = Math.abs(current_delta/true_LL);
		double current_dist = Math.abs(trueS.getGradientDiffL2(current, true,  false));
		
		double current_ddiff = current_dist/dLscale;
		
					
					
		
		boolean have_approximation = current_rdiff<=ftol && current_ddiff <= dtol ;
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#**MLRV.aCW abs "+absolute+"\trel "+relative
					+"\trdiff "+current_rdiff+"\tddiff "+current_ddiff
					+"\tcur "+current_LL+"\tdL "+current_dL
					+"\ttru "+true_LL+"/"+true_dL
					+"\t(tol "+ftol+","+dtol+")"
					+"\tapprox? "+have_approximation);
		
		boolean adjustCalculationWidth = true;
		int num_adjustments = 0; // avoid too many loops; likely never attains 33
		while (adjustCalculationWidth && num_adjustments < 33)
		{
			adjustCalculationWidth = false;
			if (!have_approximation)
			{ // try changing absolute
				int next_absolute = Integer.max(absolute+1,(int)Math.ceil(Math.exp(Math.log(absolute)+step_size)));
				this.setCalculationWidth(next_absolute, relative);
				
				OptimizationState nextS = currentState();
				double next_dL = nextS.getGradientL2(true, false);
				double next_LL = nextS.getLL();
				
				double next_delta = next_LL-true_LL;
				double next_rdiff = Math.abs(next_delta/true_LL); 
				double next_dist = Math.abs(trueS.getGradientDiffL2(nextS, true,  false));				
				double next_ddiff = next_dist/dLscale;
	
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#**MLRV.aCW ("+dabs+","+drel+") absolute "+absolute+"\tincrease "+next_absolute
							+"\trdiff "+next_rdiff+"\tddiff "+next_ddiff
							+"\twas "+current_LL+"\tnext "+next_LL
							+"\tdL "+current_dL+"\tnextdL "+next_dL
							+"\ttru "+true_LL+"/"+true_dL
							+"\t(tol "+ftol+","+dtol+")");
	
					absolute = next_absolute;
					current_LL = next_LL;
					current_dL = next_dL;
					current = nextS;
					adjustCalculationWidth = true;
					++num_adjustments;
					++dabs;
					have_approximation = next_rdiff<=ftol && next_ddiff<=dtol;
			} else if (DEFAULT_TRUNCATE_ABSOLUTE < absolute)
			{
				int prev_absolute = Integer.max(Integer.min((int)Math.ceil(Math.exp(Math.log(absolute)-step_size)), absolute-1), DEFAULT_TRUNCATE_ABSOLUTE);
				this.setCalculationWidth(prev_absolute, relative);
				
				OptimizationState prevS = currentState();
				double prev_LL = prevS.getLL();
				double prev_dL = prevS.getGradientL2(true, false);
				double prev_delta = prev_LL-true_LL;
				double prev_rdiff = Math.abs(prev_delta/true_LL);
				double prev_dist = Math.abs(trueS.getGradientDiffL2(prevS, true,  false));
				double prev_ddiff = prev_dist/dLscale;
				
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#**MLRV.aCW ("+dabs+","+drel+") absolute "+absolute+"\tdecrease "+prev_absolute
							+"\trdiff "+prev_rdiff+"\tddiff "+prev_ddiff
							+"\twas "+current_LL+"\tprev "+prev_LL
							+"\tdL "+current_dL+"\tprevdL "+prev_dL
							+"\ttru "+true_LL+"/"+true_dL
							+"\t(tol "+ftol+","+dtol+")");
				if (prev_rdiff <= ftol && prev_ddiff <= dtol)
				{
					absolute = prev_absolute;
					current_LL = prev_LL;
					current_dL = prev_dL;
					current = prevS;
					adjustCalculationWidth = true;
					++num_adjustments;
					--dabs;
					have_approximation = true;
				}
			}
			
			if (!have_approximation)
			{
				// try changing relative
				double next_rel = Math.exp(Math.log(relative)+step_size);
				this.setCalculationWidth(absolute, next_rel);
				
				OptimizationState nextS = currentState();
				double next_dL = nextS.getGradientL2(true, false);
				double next_LL = nextS.getLL();
				
				double next_delta = next_LL-true_LL;
				double next_rdiff = Math.abs(next_delta/true_LL); 
				double next_dist = Math.abs(trueS.getGradientDiffL2(nextS, true,  false));
				double next_ddiff = next_dist/dLscale;
				
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#**MRLV.aCW ("+dabs+","+drel+") relative "+relative+"\tincrease "+next_rel
							+"\trdiff "+next_rdiff+"\tddiff "+next_ddiff
							+"\twas "+current_LL+"\tnext "+next_LL
							+"\tdL "+current_dL+"\tnextdL "+next_dL
							+"\ttru "+true_LL+"/"+true_dL
							+"\t(tol "+ftol+","+dtol+")");
	
				relative = next_rel;
				current_LL = next_LL;
				current_dL = next_dL;
				current = nextS;
				adjustCalculationWidth = true;
				++num_adjustments;
				++drel;
				have_approximation = next_rdiff<=ftol && next_ddiff<=dtol;
			} else if (DEFAULT_TRUNCATE_RELATIVE < relative)
			{
				double prev_rel = Double.max(DEFAULT_TRUNCATE_RELATIVE, Math.exp(Math.log(relative)-step_size));
				this.setCalculationWidth(absolute, prev_rel);
				
				OptimizationState prevS = currentState();
				double prev_LL = prevS.getLL();
				double prev_dL = prevS.getGradientL2(true, false);
				double prev_delta = prev_LL-true_LL;
				double prev_rdiff = Math.abs(prev_delta/true_LL);
				double prev_dist = Math.abs(trueS.getGradientDiffL2(prevS, true,  false));
				double prev_ddiff = prev_dist/dLscale;
				
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#**MLRV.aCW ("+dabs+","+drel+") relative "+relative+"\tdecrease "+prev_rel
							+"\trdiff "+prev_rdiff+"\tddiff "+prev_ddiff
							+"\twas "+current_LL+"\tprev "+prev_LL
							+"\tdL "+current_dL+"\tprevdL "+prev_dL
							+"\ttru "+true_LL+"/"+true_dL
							+"\t(tol "+ftol+","+dtol+")");								
				if (prev_rdiff <= ftol && prev_ddiff <= dtol)
				{
					relative = prev_rel;
					current_LL = prev_LL;
					current_dL = prev_dL;
					current = prevS;
					adjustCalculationWidth = true;
					++num_adjustments;
					--drel;
					have_approximation = true;
				}
			}

			firePropertyChange​(PROPERTY_OPTIMIZATION_PHASE, "Adjusting truncation: "+num_adjustments+"("+absolute+","+((int)(relative*10.0+0.5))/10.0+")");
		} 
		
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#**MLRV.aCW setting "+absolute+","+relative
					+"\tapproxOK? "+have_approximation);	
		if (absolute != gradient_factory.getCalculationWidthAbsolute()
				|| relative != gradient_factory.getCalculationWidthRelative())
		{
			this.setCalculationWidth(absolute, relative);
			current = currentState();
		}
		return current;
	}
	
	
	/*
	 * 
	 * Adjustable parameters
	 * 
	 */
	
	
	/**
	 * Sets adjustable parameter values 
	 * 
	 * @param x
	 * @param want_model_parameters
	 * @param want_category_parameters
	 */
	private void setAdjustableParameterValues(double[] x, boolean want_model_parameters, boolean want_category_parameters)
	{
		int j=0;
		for (MinGradient P: adjustable_parameters)
		{
			if (!P.isCategoryParameter())
			{
				if (want_model_parameters) {
					P.set(x[j]); 
//					if (!Double.isFinite(P.get())) { // DEBUG
//						System.out.println("#**MLRV.sAP "+P+"\tx "+x[j]+"\tg "+P.get());
//					}
//					assert (Double.isFinite(P.get()));
					j++;
				}
			} else
			{
				if (want_category_parameters) {
					P.set(x[j]);
					assert (x[j]==P.get());
					j++;
				}
			}
		}
		this.copyParametersToModel();
	}

	/**
	 * Number of adjustable parameter values 
	 * 
	 * @param want_model_parameters
	 * @param want_category_parameters
	 * @return
	 */
	private int numAdjustableParameters(boolean want_model_parameters, boolean want_category_parameters)
	{
		int num_base_pars = 0;
		int num_cat_pars = 0;
		for (MinGradient P : adjustable_parameters)
		{
			if (want_model_parameters && !P.isCategoryParameter())
				num_base_pars++;
			if (want_category_parameters && P.isCategoryParameter())
				num_cat_pars++;
		}
		return num_base_pars + num_cat_pars;
	}
	
	private double[] getAdjustableParameterValues(boolean want_model_parameters, boolean want_category_parameters)
	{
		int npar = numAdjustableParameters(want_model_parameters, want_category_parameters);
		double[] x = new double[npar];
		int j=0;
		for (MinGradient P : adjustable_parameters)
		{
			if (!P.isCategoryParameter())
			{
				if (want_model_parameters)
					x[j++] = P.get();
			} else
			{
				if (want_category_parameters)
					x[j++] = P.get();
			}
		}
		return x;
	}

	private double[] getAdjustableParameterGradient(double[] dLL, boolean want_model_parameters, boolean want_category_parameters)
	{
		int npar = numAdjustableParameters(want_model_parameters, want_category_parameters);
		double[] D = new double[npar];
		int j=0;
		for (MinGradient G: adjustable_parameters)
		{
			if ((G.isCategoryParameter() && want_category_parameters)
					|| (!G.isCategoryParameter() && want_model_parameters))
			{
				D[j++]=G.dL(dLL);
			}
		}
		return D;
	}	
	
	/*
	 * Function evaluation and gradients 
	 */
	private int calls_optFunc=0;
	
//	/**
//	 * Negative log-likelihood 
//	 * 
//	 * @param x
//	 * @param want_model_parameters
//	 * @param want_category_parameters
//	 * @return
//	 */
//	private double optFunc(double [] x, boolean want_model_parameters, boolean want_category_parameters)
//	{
//		setAdjustableParameterValues(x, want_model_parameters, want_category_parameters);
//		double LL = gradient_factory.getCorrectedLL();
//		
////		if (PRINT_OPTIMIZATION_MESSAGES)
////			System.out.println("#**MLRV.oF/"+want_model_parameters+","+want_category_parameters+"\t"+calls_optFunc+"\t"+LL);
//		this.calls_optFunc++;
//
//		return -LL;
//	}
	
	private int calls_optDiff = 0;
//	private double[] optDiff(double[] x, boolean want_model_parameters, boolean want_category_parameters)
//	{
//		setAdjustableParameterValues(x, want_model_parameters, want_category_parameters);
//		double[] dLL = gradient_factory.getCorrectedGradient(!variation_model.isDuplicationBounded());
//		double[] D = getAdjustableParameterGradient(dLL, want_model_parameters, want_category_parameters);
//		
//		this.calls_optDiff++;
//		if (PRINT_OPTIMIZATION_MESSAGES)
//			System.out.println("#**MLRV.oD "+calls_optDiff+"\t"+gradient_factory.getCorrectedLL()
//					+"\tfcalls "+calls_optFunc
//					+"\t|dL| "+FunctionMinimization.euclideanNorm(D));
//		return D;
//		
//	}
	
	/**
	 * Evaluation for computing likelihood and gradient in
	 * one passage with caching. 
	 */
	private class GradientEvaluation {
		/**
		 * Cached parameter values at last evaluation
		 */
		private double[] allx;
		/**
		 * Cached class probabilities at last evaluation
		 */
		private double[] logp;
		/**
		 * Cached sample gradient
		 */
		private VariationGradientFactory.SampleGradient sgrad; 

		GradientEvaluation() {
			this.allx=null;
			this.logp=null; 
			this.sgrad = null;
		}
		
		private boolean isCached(boolean enforceCache) {
			double[] x = ML.getParameterValues(MLRateVariation.this.full_parameters);
			int ncat = variation_model.getNumClasses();
			boolean samex = Arrays.equals(x, allx) && logp.length==ncat; // if x==allx, then logp is not null either
			if (samex) {
				int ci=0;
				while (ci<ncat-1 && samex) {
					samex = this.logp[ci] == variation_model.getCategory(ci).getLogCatProbability();
					++ci;
				}
			}
			if (!samex && enforceCache) {
				this.allx = x;
				this.logp=new double[ncat];
				for (int ci=0; ci<ncat; ci++) {
					this.logp[ci] = variation_model.getCategory(ci).getLogCatProbability();
				}
				this.sgrad = gradient_factory.computeSampleGradient(!variation_model.isDuplicationBounded());
				sgrad.correctForUnobserved(!variation_model.isDuplicationBounded());
	
				calls_optDiff++;
				samex=true;
			}
			return samex;
		}
		
		/**
		 * Function value at a point x; if called after {@link #gradient(double[], boolean, boolean)}
		 * at the same x, then it is already computed, otherwise it computes 
		 * the likelihood on the spot but not the gradient.
		 * 
		 * @param x
		 * @param want_model_parameters
		 * @param want_category_parameters
		 * @return
		 */
		double apply(double[] x, boolean want_model_parameters, boolean want_category_parameters) {
			setAdjustableParameterValues(x, want_model_parameters, want_category_parameters);
			double LL;
			if (!isCached(false)) {
				LL = gradient_factory.getCorrectedLL(); // not calculating gradient 
				MLRateVariation.this.calls_optFunc++;
			} else {
				LL = sgrad.getLogLikelihood();
			}
			return -LL;
		}
		
		/**
		 * Gradient of the log-likelihood at a given point x; call before {@link #apply(double[], boolean, boolean)}.
		 * Subsequent calls with same point retrieve the cached value.
		 * 
		 * @param x
		 * @param want_model_parameters
		 * @param want_category_parameters
		 * @return gradient vector
		 */
		double[] gradient(double[] x, boolean want_model_parameters, boolean want_category_parameters) {
			setAdjustableParameterValues(x, want_model_parameters, want_category_parameters);
			isCached(true); // fills cache if different x
			double[] dLL = sgrad.get();
			double[] D = getAdjustableParameterGradient(dLL, want_model_parameters, want_category_parameters);
			if (PRINT_OPTIMIZATION_MESSAGES)
				System.out.println("#**MLRV.oD "+calls_optDiff+"\t"+gradient_factory.getCorrectedLL()
						+"\tfcalls "+calls_optFunc
						+"\t|dL| "+FunctionMinimization.euclideanNorm(D));
			return D;
		}
		
	}
	
	
	
	/*
	 * Helper classes for manipulating model parameters 
	 */
	protected MinGradient defaultNodeParameter(int node, int param_type)
	{
		final MinGradient P;
		if (param_type == PARAMETER_GAIN)
		{
			if (USE_LOGISTIC_GAIN && (MAX_KAPPA!=Double.POSITIVE_INFINITY || getCommonGainType() != PARAMETER_DUPLICATION)) {
				
				//P = new LogisticGain(node);
				MinGradient LG = new LogisticGain(node);
				if (USE_BRACKETS)
					P = new BracketedFromLogistic(LG);	
				else
					P=LG;
			} else {
				
				MinGradient LG = new LogGain(node);
				P = LG;
//				if (MAX_KAPPA==Double.POSITIVE_INFINITY || (ALWAYS_UNBOUNDED_GAIN_LOSS && variation_model.getCommonGainType()==PARAMETER_LOSS))
//				{
//					P = LG;
//				} else
//				{
//					// maximum set at 1.0 for loss- and gain-linked common gain 
//					double max_gain = getCommonGainType() == PARAMETER_DUPLICATION
//							?MAX_KAPPA:1.0;
//					
//					P = new LogisticFromLogarithmic(LG, 
//							max_gain) ;
//				}
			}
		} else if (param_type == PARAMETER_LOSS) 
		{
			P = new LogisticLoss(node);
		} else 
		{
			assert (param_type == PARAMETER_DUPLICATION);
			MinGradient LDR = new LogisticDuplicationRate(node);
			
			if (Double.isFinite(LDR.get()))
			{
				if (0.0<MAX_DUPLICATION_MARGIN)
				{
					P = new MarginsFromLogistic (LDR, MAX_DUPLICATION_MARGIN);
				} else
				{
					if (USE_BRACKETS)
						P = new BracketedFromLogistic(LDR);
					else
						P = LDR;
				}
			} else
				P = LDR;
		} 
		
//		System.out.println("#**MLRV.dNP "+node+"/"+GLDParameters.paramName(param_type)+"\t"+P);
		return P;
	}
	
	private MinGradient defaultCategoryParameter(int cat, int param_type)
	{
		MinGradient M;
		if (param_type == PARAMETER_MOD_LENGTH)
		{
			M = new ModLength(cat);
		} else
		{
			assert (param_type == PARAMETER_MOD_DUPLICATION);
			M = new ModDuplication(cat);
		}
		MinGradient catP;
		if (max_modifier == Double.POSITIVE_INFINITY)
		{
			catP = M;
		} else
		{
			catP = new LogisticFromLogarithmic(M,max_modifier);
		}
		
		return catP;
	}
	
	
	abstract class MinGradient implements ModelParameter
	{
		MinGradient(int par_idx, int emt_idx, boolean is_cat_param)
		{
			this.pidx = par_idx;
			this.is_category_param = is_cat_param;
			this.emt_idx = emt_idx;
		}
		/**
		 * Parameter index in full log-gradient from {@link VariationGradientFactory#getCorrectedGradient(boolean)}
		 */
		protected final int pidx;
		/**
		 * Whether category or node parameter (can be inferred from  {@link #pidx} &ge; 3*number of nodes 
		 */
		protected final boolean is_category_param;
		/**
		 * Node index or category index
		 */
		protected final int emt_idx; 
		
		/**
		 * Derivative of the negative log-likelihood (to be minimized)
		 * 
		 * @param dLarray gradient of the log-likelihood 
		 */
		@Override
		public double dL(double[] dLarray)
		{
			return -dLarray[this.pidx];
		}
		
		boolean isCategoryParameter()
		{
			return this.is_category_param;
		}
		
	}
	

	private class LogisticFromLogarithmic extends MinGradient 
	{
		private final MinGradient LG;
		private final double log_max_value;
		LogisticFromLogarithmic(MinGradient log_par, double max_value)
		{
			super(log_par.pidx, log_par.emt_idx, log_par.is_category_param);
			this.LG = log_par;
			this.log_max_value = Math.log(max_value);
			double log_gn = LG.get();
			if (log_max_value < log_gn)
			{
				double newval = log_max_value + Math.log(1023.0/1024.0);
				System.out.println("#**MLRV.LogisticFromLogarithmic: too large "+log_gn+"\t; resetting "+newval+"\t"+variation_model.getBaseModel().toString(LG.emt_idx));
				LG.set(newval);
			}
		}
		@Override
		public void set(double logit_t)
		{
			double log_t = Logarithms.logitToLogValue(logit_t);
			
			assert (Double.isFinite(logit_t));
			if (log_t==0.0 || !Double.isFinite(log_t)) {
				System.out.println("#**MLRV.LFL.set log "+log_t+"\tlogit "+logit_t);
			}
			
			double log_gn = log_t + log_max_value;
			LG.set(log_gn);
		}
		
		@Override 
		public double get()
		{
			double log_gn = LG.get();
			double log_t = log_gn - log_max_value;
			double logit_t = Logarithms.logToLogit(log_t);
			
			if (Double.isFinite(log_gn) && !Double.isFinite(logit_t)) { // DEBUG
				System.out.println("#**MLRV.LFL.get "+LG+"\tloggn "+log_gn+"\tlogt "+log_t+"\tlogit "+logit_t+"\tcomp "+Logarithms.logitToLogLogComplement(log_t));
			}
			
			return logit_t;
		}
		
		@Override 
		public double dL(double[] gradient)
		{
			double log_gn = LG.get();
			double log_t = log_gn - log_max_value;
			//double one_minus_t = -Math.expm1(log_t);  // 1-exp(log_t)
			double one_minus_t = Math.exp(Logarithms.logToLogComplement(log_t));
			
			double dL = LG.dL(gradient)*one_minus_t;
			return dL;
		}
		
		@Override 
		public String toString()
		{
			return "logit"+"["+LG
					+"; get "+get()
					+"; max "+Math.exp(log_max_value)
			+"]";			
		}	
	}
	
	/**
	 * 
	 * Transformation from a probability bounded away from 0 and 1 
	 * as <var>ε</var> &le; <var>p</var> &ge; (1-<var>ε</var>),
	 * represented on logistic scale by <var>y</var>=ln(<var>p</var>)-ln(1-<var>p</var>) 
	 * to unbounded parameter <var>x</var>=ln(<var>p</var>-<var>ε</var>)-ln(1-<var>p</var>-<var>ε</var>) 
	 * 
	 */
	class BracketedFromLogistic extends MinGradient
	{
		private final MinGradient MG;
		private final double log_eps;
		private final double log1_2eps;
		
		BracketedFromLogistic(MinGradient logistic_par)
		{
			this(logistic_par, 1.0/(1L<<54));
		}		
		
		BracketedFromLogistic(MinGradient logistic_par, double epsilon)
		{
			super(logistic_par.pidx, logistic_par.emt_idx, logistic_par.is_category_param);
			this.MG = logistic_par;
			this.log_eps = Math.log(epsilon);
			this.log1_2eps = Math.log1p(-2.0*epsilon);
			double log1_eps = Math.log1p(-epsilon);
			
			
			double logit_val = MG.get();
			double new_offset = 1.0/1024.0;
			if (Logarithms.logitToLogValue(logit_val)<log_eps)
			{
				double new_logval = log_eps + Math.log1p(new_offset);
				double new_logit = Logarithms.logToLogit(new_logval);
				 // 1-e*(1+o) // (1+o)*(1-e) = 1+o - (1+o)e
				
				
				System.out.println("#**MLRV.BracketedFromLogistic(): too small "+logit_val+"\t; resetting "+new_logit+"\t// "+MG);
				MG.set(new_logit);
			}
			if (log1_eps < Logarithms.logitToLogValue(logit_val))
			{
				double new_logcmp = log_eps + Math.log1p(new_offset); // log-complement 
				double new_logit = -Logarithms.logToLogit(new_logcmp);
				System.out.println("#**MLRV.BracketedFromLogistic(): too large "+logit_val+"\t; resetting "+new_logit+"\t// "+MG);
				MG.set(new_logit);
			}
		}
		
		
		@Override
		public void set(double logit_t)
		{
			double log_t = Logarithms.logitToLogValue(logit_t);
			double log1_t = Logarithms.logitToLogComplement(logit_t);
			
			double log_p = Logarithms.add(log_eps, log1_2eps+log_t);
			double log1_p = Logarithms.add(log_eps, log1_2eps+log1_t);
			
			double logit_p = log_p - log1_p; 
			
			//assert (Double.isFinite(logit_t)==);
			
			assert Double.isFinite(logit_p);
			MG.set(logit_p);
		}
		
		@Override
		public double get()
		{
			double logit_p = MG.get();
			double log_p = Logarithms.logitToLogValue(logit_p);
			double log1_p = Logarithms.logitToLogComplement(logit_p);
			double[] ld_numer = Logarithms.ldiff(log_p, log_eps);
			double[] ld_denom = Logarithms.ldiff(log1_p, log_eps);
			
			double logit_t = Logarithms.ldiffLogValue(ld_numer)-Logarithms.ldiffLogValue(ld_denom);
			
//			// DEBUG
//			if (Double.isFinite(logit_t)!=Double.isFinite(logit_p)) {
//				System.out.println("#**MLRV.BFL.get/"+MG+"\tlt "+logit_t+"\tlp "+logit_p+"\tlogp "+log_p+"\tlog1p "+log1_p+"\tle "+log_eps
//						+"\tld "+Logarithms.ldiffLogValue(ld_numer)+"-"+Logarithms.ldiffLogValue(ld_denom));
//			}
//			
//			assert (Double.isFinite(logit_t)==Double.isFinite(logit_p));

			//			System.out.println("#**MLRV.BFL.get "+logit_t+"\t// "+this.toString());
			
			return logit_t;
		}
		
		
		
		@Override 
		public double dL(double[] gradient)
		{
			double logit_p = MG.get();
			
			assert Double.isFinite(logit_p);
			double log_p = Logarithms.logitToLogValue(logit_p);
			double log1_p = Logarithms.logitToLogComplement(logit_p);
			double[] ld_numer = Logarithms.ldiff(log_p, log_eps);
			double[] ld_denom = Logarithms.ldiff(log1_p, log_eps);
			
			double logit_t = Logarithms.ldiffLogValue(ld_numer)-Logarithms.ldiffLogValue(ld_denom);
			
			double log_t = Logarithms.logitToLogValue(logit_t);
			double log1_t = Logarithms.logitToLogComplement(logit_t);
			
			double log_mul = (log_t+log1_t)-(log_p+log1_p)+log1_2eps;
			double dL_dlogitp = MG.dL(gradient);
			double dL = dL_dlogitp * Math.exp(log_mul);
			
//			System.out.println("#**MLRV.BFL.dL "+dL
//					+"\tlogmul "+log_mul+"/*"+Math.exp(log_mul)
//					+"\tdLlp "+dL_dlogitp
//					+"\tlt "+logit_t+"("+log_t+"/"+log1_t+")"
//					+"\tlp "+logit_p+"("+log_p+"/"+log1_p
//						+", "+Logarithms.ldiffLogValue(ld_denom)+")"
//						+"\tl12e "+log1_2eps
//					+"\t// "+this.toString());
			return dL;
		}		
		@Override 
		public String toString()
		{
			return "brk["+MG.toString()
			+"]";			
		}	
	}
	
	
	private class MarginsFromLogistic extends MinGradient
	{
		private final MinGradient MG;
		
		private final double log_eps_lo;
		private final double log_eps_hi;
		private final double log_range;
		
		MarginsFromLogistic(MinGradient logistic_par, double eps_hi)
		{
			this(logistic_par, 0.25*Math.ulp(1.0), eps_hi);
		}
		
		MarginsFromLogistic(MinGradient logistic_par, double eps_lo, double eps_hi)
		{
			super(logistic_par.pidx, logistic_par.emt_idx, logistic_par.is_category_param);
			this.MG = logistic_par;
			this.log_eps_lo = Math.log(eps_lo);
			this.log_eps_hi = Math.log(eps_hi);
			this.log_range = Math.log1p(-(eps_lo+eps_hi)); 
			double log1_eps = Math.log1p(-eps_hi);
			
			double logit_val = MG.get();
			double new_offset = 1.0/1024.0;
			if (Logarithms.logitToLogValue(logit_val)<log_eps_lo)
			{
				//double logit_eps = Logarithms.toLogit(eps_lo);
				double new_logval = log_eps_lo + Math.log1p(new_offset);
				double new_logit = Logarithms.logToLogit(new_logval);
				 // 1-e*(1+o) // (1+o)*(1-e) = 1+o - (1+o)e
				
				
				System.out.println("#**MLRV.MarginsFromLogistic(): too small "+logit_val+"\t; resetting "+new_logit+"\t// "+MG);
				MG.set(new_logit);
			}
			if (log1_eps < Logarithms.logitToLogValue(logit_val))
			{
				//double logit_eps = Logarithms.toLogit(eps_hi);
				double new_logcmp = log_eps_hi + Math.log1p(new_offset); // log-complement 
				double new_logit = -Logarithms.logToLogit(new_logcmp);
				System.out.println("#**MLRV.MarginsFromLogistic(): too large "+logit_val+"\t; resetting "+new_logit
						+"\tlog1_eps "+log1_eps
						+"\tlogval "+Logarithms.logitToLogValue(logit_val)
						+"\tnewlogit "+new_logit+"/"+Math.exp(Logarithms.logitToLogValue(new_logit))
						+"\t// "+MG);
				MG.set(new_logit);
			}			
		}		
		@Override
		public void set(double logit_t)
		{
			double log_t = Logarithms.logitToLogValue(logit_t);
			double log1_t = Logarithms.logitToLogComplement(logit_t);
			
			double log_p = Logarithms.add(log_eps_lo, log_range+log_t);
			double log1_p = Logarithms.add(log_eps_hi, log_range+log1_t);
			
			double logit_p = log_p - log1_p; 
			MG.set(logit_p);
		}
		@Override
		public double get()
		{
			double logit_p = MG.get();
			double log_p = Logarithms.logitToLogValue(logit_p);
			double log1_p = Logarithms.logitToLogComplement(logit_p);
			double[] ld_numer = Logarithms.ldiff(log_p, log_eps_lo);
			double[] ld_denom = Logarithms.ldiff(log1_p, log_eps_hi);
			
			double logit_t = Logarithms.ldiffLogValue(ld_numer)-Logarithms.ldiffLogValue(ld_denom);
			
//			System.out.println("#**MLRV.BFL.get "+logit_t+"\t// "+this.toString());
			
			return logit_t;
		}
		public double dL(double[] gradient)
		{
			double logit_p = MG.get();
			double log_p = Logarithms.logitToLogValue(logit_p);
			double log1_p = Logarithms.logitToLogComplement(logit_p);
			double[] ld_numer = Logarithms.ldiff(log_p, log_eps_lo);
			double[] ld_denom = Logarithms.ldiff(log1_p, log_eps_hi);
			
			double logit_t = Logarithms.ldiffLogValue(ld_numer)-Logarithms.ldiffLogValue(ld_denom);
			
			double log_t = Logarithms.logitToLogValue(logit_t);
			double log1_t = Logarithms.logitToLogComplement(logit_t);
			
			double log_mul = (log_t+log1_t)-(log_p+log1_p)+log_range;
			double dL_dlogitp = MG.dL(gradient);
			double dL = dL_dlogitp * Math.exp(log_mul);
			
//			System.out.println("#**MLRV.BFL.dL "+dL
//					+"\tlogmul "+log_mul+"/*"+Math.exp(log_mul)
//					+"\tdLlp "+dL_dlogitp
//					+"\tlt "+logit_t+"("+log_t+"/"+log1_t+")"
//					+"\tlp "+logit_p+"("+log_p+"/"+log1_p
//						+", "+Logarithms.ldiffLogValue(ld_denom)+")"
//						+"\tl12e "+log1_2eps
//					+"\t// "+this.toString());
			return dL;
		}			
		@Override 
		public String toString()
		{
			return "mrg["+MG.toString()
			+"]";			
		}	
		
		
		
		
	}
	
	private class LogisticLoss extends MinGradient
	{
		LogisticLoss(int node)
		{
			super(3*node+PARAMETER_LOSS, node,  false);
			this.node = node;
		}
		private final int node;
		
		@Override
		public double get()
		{
			double logit_p = node_parameters[node][PARAMETER_LOSS];
			return logit_p;
		}
		
		
		@Override
		public void set(double logit_p)
		{
			node_parameters[node][PARAMETER_LOSS]=logit_p;
			
			// DEBUG
			if (variation_model.getBaseModel().getTree().isRoot(node)) {
				System.out.println("#**MLRV.LL.set root "+node+"\t"+this.toString()+"\t// "+variation_model.getBaseModel().toString(node));
			}
		}
		@Override
		public double dL(double[] dLarray)
		{
			double dL = super.dL(dLarray);
			// DEBUG
			if (variation_model.getBaseModel().getTree().isRoot(node)) {
				System.out.println("#**MLRV.LL.dL root "+node+"\tdLLdlogitp "+(-dL)+"\t"+this.toString()); //+"\t// "+variation_model.getBaseModel().toString(node));
			}
			return dL;
		
		}
		
		@Override 
		public String toString()
		{
			return "logitp"+node+"["
					+Math.exp(Logarithms.logitToLogValue(get()))
					+"/1-"+Math.exp(Logarithms.logitToLogComplement(get()))
					+"; get "+get()
			+"]";			
		}
	}
	
	/**
	 * Logistic transformation for relative duplication
	 * rate lambda = p/q, bounded at 1.0. 
	 */
	private class LogisticDuplicationRate extends MinGradient
	{
		LogisticDuplicationRate(int node)
		{
			super(3*node+PARAMETER_DUPLICATION, node, false);
			this.node = node;
		}
		private final int node;
		
		@Override
		public double get()
		{
			double logit_lm = node_parameters[node][PARAMETER_DUPLICATION];
			
			return logit_lm;
		}
		
		
		@Override
		public void set(double logit_lambda)
		{
			node_parameters[node][PARAMETER_DUPLICATION]=logit_lambda;
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
	
	
	
	
	private class LogGain extends MinGradient
	{
		LogGain(int node)
		{
			super(3*node+PARAMETER_GAIN, node, false);
			this.node = node;
		}
		private final int node;
		
		@Override
		public double get()
		{
			double log_gn = node_parameters[node][PARAMETER_GAIN];
			return log_gn;
		}
		
		
		@Override
		public void set(double log_gn)
		{
			node_parameters[node][PARAMETER_GAIN]=log_gn;
		}
		@Override 
		public String toString()
		{
			TreeWithLogisticParameters lrates = variation_model.getBaseModel();
			return "loggn"+node+"["+get()
			+", gn "+Math.exp(get())
			+"/grate "+lrates.getGainRate(node)
			+"]";			
		}	
		
	}
	
	
	private class LogisticGain extends MinGradient {
		LogisticGain(int node)
		{
			super(3*node+PARAMETER_GAIN, node, false);
			assert (USE_LOGISTIC_GAIN);
			this.node = node;
		}
		private final int node;
		
		@Override
		public double get()
		{
			double get = node_parameters[node][PARAMETER_GAIN];
			return get;
		}
		
		@Override
		public void set(double logit_gn)
		{
			node_parameters[node][PARAMETER_GAIN]=logit_gn;
		}
		
		@Override 
		public double dL(double[] gradient)
		{
			double logit_gn = get();
			double log_1g = Logarithms.logitToLogComplement(logit_gn);
			
			double dL = super.dL(gradient);
			dL *= Math.exp(log_1g);
			return dL;
		}
		
		
		@Override 
		public String toString()
		{
			return "logitgn"+node+"["
					+Math.exp(Logarithms.logitToLogValue(get()))
					+"/1-"+Math.exp(Logarithms.logitToLogComplement(get()))
					+"; get "+get()
					+"]";			
		}	
		
		
	}

	private class ModLength extends MinGradient
	{
		ModLength(int cat)
		{
			super(3*node_parameters.length+2*cat+PARAMETER_MOD_LENGTH, cat, true);
			this.cat = cat;
		}
		private final int cat;
		
		@Override
		public double get() 
		{
			double mod_len = category_parameters[cat][PARAMETER_MOD_LENGTH];
			return mod_len;
		}
		@Override
		public void set(double mod_len) 
		{
			category_parameters[cat][PARAMETER_MOD_LENGTH]=mod_len;
		}
		@Override
		public String toString()
		{
			return "modlen"+cat+"["+get()+", mul "+Math.exp(get())+"]";
		}
	}
	
	private class ModDuplication extends MinGradient
	{
		ModDuplication(int cat)
		{
			super(3*node_parameters.length+2*cat+PARAMETER_MOD_DUPLICATION, cat, true);
			this.cat = cat;
		}
		private final int cat;
		
		@Override
		public double get() 
		{
			double mod_dup = category_parameters[cat][PARAMETER_MOD_DUPLICATION];
			return mod_dup;
		}
		@Override
		public void set(double mod_dup) 
		{
			category_parameters[cat][PARAMETER_MOD_DUPLICATION]=mod_dup;
		}
		@Override
		public String toString()
		{
			return "moddup"+cat+"["+get()+", mul "+Math.exp(get())+"]";
		}
	}
	
	
	/*
	 * 
	 * Optimization
	 * 
	 */
	
	/**
	 * Model parameters: up to 3 node parameters per node, 
	 * 	category parameters, and category probabilities
	 */
	@Override
	public int getModelParameterCount()
	{
//		System.out.println("#**MLRV.gMPC adjustable "+adjustable_parameters.size());
		
		int getModelParameterCount = 0;
		for (int j=0; j<do_optimize_parameters.length; j++)
		{
			if (do_optimize_parameters[j])
				++getModelParameterCount;
		}
		
		int num_prob_params = variation_model.getNumClasses()-1;
		
		getModelParameterCount += num_prob_params;
		
		return getModelParameterCount;
	}

	
	
	/**
	 * Expectation-maximization for category prior probabilities
	 * @param current
	 * @param eps small change in probability that shortcuts the EM iterations
	 * @return new optimization state 
	 */
	protected OptimizationState optimizeCategoryProbabilities(OptimizationState current, double eps, List<Double> history)
	{

		final int EMiter = 96;
		final double EMdelta = 0.1; 
		
		int ncat = category_parameters.length;
		
		if (1<ncat)
		{
			OptimizationState prev;
			for (int iter=0; iter<EMiter; iter++)
			{
				prev = current;
				double max_dp = 0.0;
				for (int k=0; k<ncat; k++)
				{
					RateVariationModel.Category C = variation_model.getCategory(k);
					double prior_p = Math.exp(C.getLogCatProbability());
					double log_post_p = current.G.getLogPosteriorFrequency(k);
					double post_p = Math.exp(log_post_p);
					
					double dp = Math.abs(post_p-prior_p);
					max_dp = Double.max(dp, max_dp);
					
					C.setLogCatProbability(log_post_p);
				}
				current = currentState();
				double LLdelta = current.getLL()-prev.getLL();
				
				if (PRINT_OPTIMIZATION_MESSAGES)
				{
					System.out.println("#*MLRV.oCProb "+iter+"\tcur "+current+"\tmax_dp "+max_dp+"\tLLinc "+LLdelta);
				}
				if (history != null)
				{
					history.add(-current.getLL());
				}
				if (max_dp<eps || LLdelta<EMdelta) break;
			}
			for (int k=0; k<ncat; k++)
			{
				RateVariationModel.Category C = variation_model.getCategory(k);
				System.out.println("#*MLRV.oCProb cat "+k+" "+C.toString());
			}
			
		}
		return current;
	}
	
	/**
	 * Optimization for category parameters without using the gradient 
	 * 
	 * @param delta
	 * @param itmax
	 * @param history
	 */
	protected void optimizeCategoryParams(double delta, int itmax, List<Double> history)
	{
		final List<MinGradient> cat_param_list = new ArrayList<>();
		
		for (MinGradient par: adjustable_parameters)
		{
			if (par.isCategoryParameter())
				cat_param_list.add(par);
		}
		
		if (0<cat_param_list.size())
		{
			double[] x0 = ML.getParameterValues(cat_param_list);
			
			Function<double[], Double> cat_param_func
			= new Function<>()
				{
					@Override
					public Double apply(double[] x) 
					{
						ML.setParameterValues(x, cat_param_list);
						MLRateVariation.this.copyParametersToModel();
						double LL = gradient_factory.getCorrectedLL();
						MLRateVariation.this.calls_optFunc++;
//						if (PRINT_OPTIMIZATION_MESSAGES)
//						{
//							System.out.println("#**MLV.oCPar "+LL);
//						}
						return -LL;
					}
				};
			double opt = FunctionMinimization.powell(x0, delta, itmax, cat_param_func, history);
			// set optimum parameter values
			ML.setParameterValues(x0, cat_param_list);
			MLRateVariation.this.copyParametersToModel();
		}
	}
	
	
	
	
	
	

	@Override
	public double optimize(double delta)
	{
		return this.optimize(delta, Integer.MAX_VALUE);
	}
	
	@Override
	public double optimize(double delta, final int itmax)
	{
		List<Double> history = getOptimizationHistory();
		if (history == null)
		{
			history = new ArrayList<>();
			setOptimizationHistory(history);
		}
		this.copyParametersFromModel();

		OptimizationState current_state = currentState();
		
		double LL =  current_state.getLL(); //gradient.getCorrectedLL();

		if (PRINT_OPTIMIZATION_MESSAGES)
		{
			System.out.println("#*MRLV.o model "+variation_model.getClass().getSimpleName()+"/"+variation_model.getNumClasses()+" classes"
					+"; cgain "+GLDParameters.paramName(variation_model.getCommonGainType())+", dupbound "+variation_model.isDuplicationBounded()
					+"; gainbound "+USE_LOGISTIC_GAIN+", bracketed "+USE_BRACKETS
					+"; optimization "+(USE_CONJUGATE_GRADIENT==null?"BFGS":"CG."+USE_CONJUGATE_GRADIENT.toString()));
			//FunctionMinimization.DEBUG_ZLNSRCH = true;
		}
		
		
		
		if (REGULARIZE_GAIN) {
			TreeWithLogisticParameters lrates = variation_model.getBaseModel();
			boolean have_adjusted = lrates.regularizeGainRates();
			this.initParameters();
			this.copyParametersFromModel();
			if (have_adjusted) {
				OptimizationState adjusted_state = currentState();
				double adjustedLL = adjusted_state.getLL();
				if (PRINT_OPTIMIZATION_MESSAGES) {
					System.out.println("#*MRLV.o regGain/init\tLL "+adjustedLL+"\twas "+LL);
				}
				LL = adjustedLL;
				current_state = adjusted_state;
			}
		}
		
		
//		if (PRINT_OPTIMIZATION_MESSAGES)
//		{
//			double nano = 1e-9;
//			// get timing info
//			double time_gradient = current_state.nanoTime()*nano;
//			long t0 = System.nanoTime();
//			double forgettable = gradient_factory.getCorrectedLL();
//			double time_eval = (System.nanoTime()-t0)*nano;
//			
//			System.out.printf("#*MRLV.o starting %s\tevaltime %.3fs\tgradienttime %.3fs\tparameter count %d\n",
//					current_state.toString(), time_eval, time_gradient, getModelParameterCount());
////			printAdjustableParameters(System.out);		
//		}
		int h0 = history.size();
		int step_count = history.size()-h0;
		int nepoch = 0;
		
		final double truncate_precision = Double.min(delta/8.0,1.0/(1L<<33));
		if (auto_truncate)
		{
			current_state = this.adjustCalculationWidth(current_state, truncate_precision);
		}
		
		if (PRINT_OPTIMIZATION_MESSAGES)
		{
			double nano = 1e-9;
			// get timing info
			double time_gradient = current_state.nanoTime()*nano;
			long t0 = System.nanoTime();
			double forgettable = gradient_factory.getCorrectedLL();
			double time_eval = (System.nanoTime()-t0)*nano;
			
			System.out.printf("#*MRLV.o starting %s\tevaltime %.3fs\tgradienttime %.3fs\tnum.parameters %d\tadjustable %d\tfactoryLL %f\tmodel %s\n",
					current_state.toString(), time_eval, time_gradient, getModelParameterCount(), numAdjustableParameters(true,true), forgettable, variation_model.toString());
		}
		
		OptimizationState best = current_state;
		LL = current_state.getLL();
		double LLbest = LL;
		
		if (step_count < itmax && variation_model.isDuplicationBounded())
		{
			boolean have_adjusted = this.reduceInfiniteDuplicationRates(1.0/(1L<<40));
			current_state = currentState();
			if (have_adjusted && PRINT_OPTIMIZATION_MESSAGES)
			{
				System.out.println("#*MLRV.o  adjusted ==1.0 duplication rates "+current_state);
			}
		}
		OptimizationState prev_state = current_state;	
		
		final GradientEvaluation eval = new GradientEvaluation();
		
		while (step_count<itmax) // optimization loop until convergence or max iterations
		{
			
			long epochT0 = System.nanoTime();
			// timing
			/*
			 *  I. adjust category parameters, if not constant rates 
			 */
			
			
			if (1<variation_model.getNumClasses() ) 
			{
				firePropertyChange​(PROPERTY_OPTIMIZATION_PHASE, "Setting rate categories");

				// adjust category parameters
				int opti_cat = Integer.min((itmax-step_count)/2, 32);
				double cat_delta = Double.max(delta, 1.0/(1L<<24));
				
				this.optimizeCategoryParams(cat_delta, opti_cat, history);
				current_state = currentState();
				
				if (PRINT_OPTIMIZATION_MESSAGES)
				{
					System.out.println("#*MLRV.o optcat\t"+current_state);
					for (int k=0; k<variation_model.getNumClasses(); k++)
					{
						RateVariationModel.Category C = variation_model.getCategory(k);
						System.out.println("#*MLRV.o\tcat"+k+" "+C.toString());
					}
				}			
				// adjust category probabilities
				double eps_p = Double.max(delta,1.0/(1L<<16)); // likelihood is not as sensitive to class probs
				current_state = optimizeCategoryProbabilities(current_state, eps_p, history);
				
				// category priors changed: recenter?
				boolean recenter = recenterBaseModel();
				if (recenter)
				{
					OptimizationState centered = currentState();
					if (PRINT_OPTIMIZATION_MESSAGES)
					{ 
						System.out.println("#*MLRV.main recenter: switched base category to max. probability: "+centered
								+"\twas "+current_state);
					}
					current_state = centered;
				}
				
			}
			
			/*
			 * II. numerical optimization for base model parameters with fixed categories 
			 */
			firePropertyChange​(PROPERTY_OPTIMIZATION_PHASE, "Epoch "+(1+nepoch)+" ("+this.getClass().getSimpleName()+")");

			step_count = history.size()-h0;
			double[] x0 = getAdjustableParameterValues(true,false);
			
			for (int i=0; i<x0.length; i++) // DEBUG
				assert Double.isFinite(x0[i]);
			
			// we want at least a few number times n iterations of bfgs to build up 
			// and exploit a solid inverse Hessian approximation
			// but we need to update class parameters, so not too many 
			int muln = nepoch<2?2:Integer.min(nepoch,4);
			
			int preferred_steps = Integer.max(128,1+muln*x0.length);
			//if (nepoch < 12) preferred_steps = (preferred_steps+2)/3;
			
			int opti_steps;
			if (itmax-step_count<100)
				opti_steps = itmax-step_count;
			else if (nepoch==0 && (1<variation_model.getNumClasses() || auto_truncate))
				opti_steps =  Integer.min((preferred_steps+2)/3,(itmax-step_count)/4);
			else if (nepoch<=1)
				opti_steps = Integer.min((preferred_steps+2)/3,(itmax-step_count)/2);
			else
				opti_steps = Integer.min(preferred_steps,(itmax-step_count));

			double min;
			if (USE_CONJUGATE_GRADIENT==null) {
				min = FunctionMinimization.quasiNewtonMin(x0, delta, opti_steps, x->eval.apply(x, true, false), x->eval.gradient(x, true, false), history);
				// min = FunctionMinimization.dfpmin(x0, delta, opti_steps, x->eval.apply(x, true, false), x->eval.gradient(x, true, false), history);
			} else {
				min = FunctionMinimization.conjugateGradientMin(x0, delta, opti_steps, x->eval.apply(x,  true,  false), x->eval.gradient(x, true, false), history, USE_CONJUGATE_GRADIENT);
			}
			
			if (Thread.interrupted()) break; // clear status
			
			for (int i=0; i<x0.length; i++) // DEBUG
				assert Double.isFinite(x0[i]);
			
			setAdjustableParameterValues(x0,true,false);
			x0 = getAdjustableParameterValues(true,false);
			
			current_state = currentState();
			LL = current_state.getLL();

			// hopefully this is a better state 
			
			if (LL>LLbest) // both negative
			{
				best = current_state;
				LLbest = LL;
			}
			
			// check convergence 
			double max_dx = current_state.maxRelDelta(prev_state);
			double gradient_length = current_state.getGradientLmax(true,false); //current_state.getGradientL2(true,false);
			double rel_gradient = gradient_length/(-LL);

			boolean done_dx = (max_dx < FunctionMinimization.DFP_TOLX);
			boolean done_gradient =  gradient_length<delta; //rel_gradient<delta; //
			boolean done_converged = ( done_gradient|| done_dx);
			
			step_count = history.size()-h0;
			
			// timing
			double epoch_seconds = (System.nanoTime()-epochT0)*1e-9;
			
			if (PRINT_OPTIMIZATION_MESSAGES)
			{
				System.out.println("#*MLRV.o epoch "+nepoch
						+"\tLL "+LL+"\t(min "+min+")"
						+"\tfcalls "+calls_optFunc+"\tdcalls "+calls_optDiff
						+"\t|dL| "+gradient_length
						+"\trgrad "+rel_gradient+"/"+done_gradient
						+"\tmax_dx "+max_dx+"/"+done_dx+"\t"+(done_converged?"DONE":"loop")
						+"\tsteps "+step_count
						+"\tdelta "+delta
						+"\ttiming "+epoch_seconds+" s");
			}

			if (DEBUG_GRADIENT) debugGradient(System.out);
			
			if (done_converged)
			{
				if (0<nepoch && done_dx) break; // avoid restarting again 
				
				if (done_gradient) break;
			}
			
			boolean has_infty = false;
			for (int i=0; i<x0.length && !has_infty; i++) // DEBUG
				has_infty = Double.isInfinite(x0[i]);
			
			
			if (REGULARIZE_GAIN) {
				TreeWithLogisticParameters lrates = variation_model.getBaseModel();
				boolean have_adjusted = lrates.regularizeGainRates();
				if (have_adjusted) {
					this.initParameters();
					this.copyParametersFromModel();
					OptimizationState adjusted_state = currentState();
					double adjustedLL = adjusted_state.getLL();
					if (PRINT_OPTIMIZATION_MESSAGES) {
						System.out.println("#*MRLV.o regGain/loop\tLL "+adjustedLL+"\twas "+LL+"\t(epoch "+nepoch+")");
					}
					if (LL==LLbest) {// likely
						LLbest = adjustedLL; // do not change it back later
						best = adjusted_state;
					}
					LL = adjustedLL;
					current_state = adjusted_state;
				}
			}				
			
			if (auto_truncate && step_count<itmax)
			{
				current_state = this.adjustCalculationWidth(current_state, truncate_precision);
			}
			
			if (step_count < itmax && variation_model.isDuplicationBounded())
			{
				boolean have_adjusted = this.reduceInfiniteDuplicationRates(1.0/(1L<<20));
				current_state = currentState();
				if (have_adjusted && PRINT_OPTIMIZATION_MESSAGES)
				{
					System.out.println("#*MLRV.o  reboost ==1.0 duplication rates "+current_state);
				}
			}
			
			nepoch++;
			
			prev_state = current_state;
		} // end of optimization loop 
		
		
		if (LL<LLbest-1.0 || !Double.isFinite(LL))
		{
			if (PRINT_OPTIMIZATION_MESSAGES)
				System.out.println("#*MLRV.o returning to previous optimum\t"+best+"\tinstead of "+current_state);
			LL = LLbest;
			current_state = best;
			current_state.resetModel();
		} else {
			if (REGULARIZE_GAIN) {
				TreeWithLogisticParameters lrates = variation_model.getBaseModel();
				boolean have_adjusted = lrates.regularizeGainRates();
				this.initParameters();
				this.copyParametersFromModel();
				if (have_adjusted) {
					OptimizationState adjusted_state = currentState();
					double adjustedLL = adjusted_state.getLL();
					if (PRINT_OPTIMIZATION_MESSAGES) {
						System.out.println("#*MRLV.o regGain/final\tLL "+adjustedLL+"\twas "+LL+"\t(final)");
					}
					LL = adjustedLL;
					current_state = adjusted_state;
				}
			}					
		}
		
		
		return -LL;
		
		
	}

	
	/*
	 * 
	 * Tests and main
	 * 
	 */
	private void debugGradient(PrintStream out)
	{
		boolean wantCategories = false; // category gradients are not used 
		
		double[] x = getAdjustableParameterValues(true,wantCategories);
		
		final GradientEvaluation eval = new GradientEvaluation();
		
		double[] df = eval.gradient(x,true,wantCategories);
		
		double[] df_est = FunctionMinimization.numericalGradientTwoPoint(θ->eval.apply(θ, true, wantCategories), x); //  optFunc(θ,true,wantCategories), x);
		
//		double[] xafter = getAdjustableParameterValues(true,true);

		double sum_delta = 0.0;
		double sumsq_delta = 0.0;
		
		for (int i=0; i<x.length; i++)
		{
			MinGradient P = adjustable_parameters.get(i);
			
			double df_delta = df_est[i]-df[i];
			double rel_delta = Math.abs(df_delta/df[i]);
			sum_delta += df_delta;
			sumsq_delta += df_delta*df_delta;
			
			out.println("#**MLRV.dG param "+i
							+"\tdfdx "+df[i]+"\test "+df_est[i]
							+"\tdiff "+df_delta+"\t("+rel_delta+")"
							+"\t"+P
							//+"\tx0 "+x[i]+"\tx1 "+xafter[i] // checking if parameters stay unchanged
						);	
		}		
		double df_len = FunctionMinimization.euclideanNorm(df);
		double avg_delta = sum_delta/df_len;
		double sd_delta =  Math.sqrt(sumsq_delta/df_len);
				
		out.println("#**MLRV.dG dflen "+df_len+"\tsumdelta "+sum_delta
				+"\tavgdelta "+avg_delta+"\tsd_delta "+sd_delta
				);
	}
	
	
	
	private void debugModel(PrintStream out)
	{
		out.println("#**MLRV.dM model same: "+(variation_model==gradient_factory.getVariationModel()));
		
	}
	
	
	
	private static final boolean LIKELIHOOD_INTERVAL_OPTIMIZE_START = false;

	
	/**
	 * LRT test
	 * 
	 * @param node
	 * @param param_type
	 * @param significance_alpha chi-square-test significance level: negative for lower tail; positive for upper tail 
	 * @param delta
	 * @param itmax
	 */
	private double findLikelihoodInterval(int node, int param_type, double significance_alpha, double delta, int itmax)
	{
		PrintStream out = System.out;
				
		final double negLL0 = optimize(delta, LIKELIHOOD_INTERVAL_OPTIMIZE_START?itmax:0);
		
		TreeWithLogisticParameters lrates = variation_model.getBaseModel();
		double logit_p0 = lrates.getLogitLossParameter(node);
		double logit_lambda0 =lrates.getLogitRelativeRate(node);
		double log_gain0 = lrates.getLogGainParameter(node, variation_model.getCommonGainType());
		
		double x0;
		if (PARAMETER_LOSS == param_type)
		{
			x0 = logit_p0;
		} else if (PARAMETER_DUPLICATION == param_type)
		{
			x0 = logit_lambda0;
		} else 
		{
			assert (PARAMETER_GAIN == param_type);
			x0 = log_gain0;
		}
		out.println("#**MLRV.fLI start\t"+node+"/"+param_type
				+"\tx0 "+x0
				+"\tlp "+logit_p0
				+"\tllm "+logit_lambda0  
				+"\tlgn " +log_gain0
				+"\tnegLL0 "+negLL0
				+"\t// "+lrates.toString(node));
		// diff_cache stores difference in log-likelihood for x settings  
		final Map<Double,Double> diff_cache = new HashMap<>();
		diff_cache.put(x0, 0.0);
		// param_cache stores best model parameters for given x settings
		final Map<Double, OptimizationState> param_cache = new HashMap<>();
		final OptimizationState state0 = currentState();
		param_cache.put(x0, state0);
		
		DoubleFunction<Double> diffLL  = new DoubleFunction<>()
		{
			@Override
			public Double apply(double x)
			{
				double diffLL;
				if (diff_cache.containsKey(x)) 
				{
					diffLL= diff_cache.get(x);
				} else
				{
					
					double logit_p, logit_lambda, log_gain;
					if (PARAMETER_LOSS == param_type)
					{
						logit_p = x;
					} else
					{
						logit_p = logit_p0;
					}
					if (PARAMETER_DUPLICATION == param_type)
					{
						logit_lambda = x;
					} else
					{
						logit_lambda = logit_lambda0;
					}
					if (PARAMETER_GAIN == param_type)
					{
						log_gain = x;
					} else
					{
						log_gain = log_gain0;
					}
					
					lrates.setLogitLossRelativeDuplication(node, logit_p, logit_lambda, log_gain, variation_model.getCommonGainType());
					
					initParameters(); // bc we need to put the model parameters into our nodeparameters

					if (PARAMETER_LOSS == param_type)
					{
						fixLoss(node,true);
					} else if (PARAMETER_DUPLICATION == param_type)
					{
						fixDuplication(node, true);
					} else
					{
						fixGain(node, true);
					}
					
					copyParametersFromModel(); // after fix
					
					
					double negLL = optimize(delta, itmax);
					diffLL = negLL-negLL0;
					OptimizationState current = currentState();
					diff_cache.put(x, diffLL);
					param_cache.put(x, current);
					
					out.println("#**MLRV.fLI.dL put\t"+x
							+"\tdiff "+diffLL+"\tx0 "+x0+"\tLL "+negLL+"\tLL0 "+negLL0
							+"\tlp "+lrates.getLogitLossParameter(node)
							+"\tllm "+lrates.getLogitRelativeRate(node)  
							+"\tlgn " +lrates.getLogGainParameter(node, variation_model.getCommonGainType())
							+"\t// "+lrates.toString(node));
					if (diffLL<0.0) // DEBUG
					{
						out.println("#**MLRV.fLI.dL decrease "+diffLL+"\tnegLL "+negLL+"\tnegLL0 "+negLL0);
						
						RateVariationModel vmodel = new RateVariationModel(lrates);
						vmodel.initConstantRates();
						out.println(RateVariationParser.printRates(vmodel));
						System.exit(99);
					}
					
					// reset
					if (PARAMETER_LOSS == param_type)
					{
						fixLoss(node,false);
					} else if (PARAMETER_DUPLICATION == param_type)
					{
						fixDuplication(node, false);
					} else
					{
						fixGain(node, false);
					}
					state0.resetModel();
				}
				// LRT test
				double chi_square_p = Functions.Chi_square_tail(1, 2.0*Math.max(0.0, diffLL)); 
				out.println("#**MLRV.fLI.dL\t"+x+"\tdiff "+diffLL+"\tp "+chi_square_p);
				return Math.abs(significance_alpha) - chi_square_p;
			}
		};	
		double step = (PARAMETER_GAIN == param_type?1.25:1.0);			
		int bracket_iter = 30;
		double xtol = 1.0/(1<<20); 

		OptimizationState state_interval;
		
		if (0.0<significance_alpha)
		{
			// bracket for xmax
			
			double x2,d2;
			if (x0==Double.NEGATIVE_INFINITY)
			{
				x2 = -Math.log(1L<<52);
				d2 = diffLL.apply(x2);
			} else
			{
				x2 = x0;
				d2 = significance_alpha-1.0;
			}
			assert (d2<0.0);

			double x1,d1; // left bracket 
			int iter=0;
			do
			{
				x1 = x2;
				d1 = d2;
				if (PARAMETER_GAIN == param_type)
					x2 = x1*step;
				else
					x2 = x1+step;
				d2 = diffLL.apply(x2);
				++iter;
			} while (d2<0.0 && iter<bracket_iter);
			
			double xmax;
			if (d2==0.0) // unlikely
			{
				// got rmax 
				xmax = x2;
			} else if (d2<0) // at maxiter 
			{
				xmax = x2;
			} else 
			{
				// d2>0
				out.println("#**MLRV.fLI brackethi\t"+x1+"\t"+x2);
				xmax = FunctionMinimization.zbrent(diffLL, x1, x2, xtol);
			}
			
			if (PRINT_OPTIMIZATION_MESSAGES)
			{
				Double[] cached_x =  diff_cache.keySet().toArray(new Double[0]);
				Arrays.sort(cached_x);
				for (double x: cached_x)
				{
					out.println("#*MLRV.fLI\t"+node+"/"+param_type+"\t"+x
							+"\t"+diff_cache.get(x)
							+"\t"+(Math.abs(significance_alpha)-diffLL.apply(x))
							);
				}
				
			}
			
			
			
			double dmax = diffLL.apply(xmax); // caches model for rmax
			state_interval = param_cache.get(xmax);
			
			state_interval.resetModel();
			out.println("#**MLRV.fLI\txmax "+xmax+"\tmax "+dmax+"\tx0 "+x0+"\t// "+lrates.toString(node));
		} else
		{
			// bracket for xmin
			
			double x1,d1;
			if (x0==Double.POSITIVE_INFINITY)
			{
				x1 = Math.log(1L<<52);
				d1 = diffLL.apply(x1); 
			} else
			{
				x1 = x0;
				d1 = -significance_alpha-1.0;
			}
			double x2,d2;
			int iter = 0;
			do
			{
				x2 = x1;
				d2 = d1;
				if (PARAMETER_GAIN == param_type)
					x1 = x2/step;
				else
					x1 = x2-step;
				d1 = diffLL.apply(x1);
				
				if (diff_cache.get(x1)<0.0) break; // what?
				
				++iter;
			} while (d1<0.0 && iter<bracket_iter); 
			
			double xmin;
			if (d1<=0.0) xmin = x1;
			else
			{
				out.println("#**MLRV.fLI bracketlow\t"+x1+"\t"+x2);
				xmin = FunctionMinimization.zbrent(diffLL, x1, x2, xtol);
			}

			
			if (PRINT_OPTIMIZATION_MESSAGES)
			{
				Double[] cached_x =  diff_cache.keySet().toArray(new Double[0]);
				Arrays.sort(cached_x);
				for (double x: cached_x)
				{
					out.println("#*MLRV.fLI\t"+node+"/"+param_type+"\t"+x
							+"\t"+diff_cache.get(x)
							+"\t"+(Math.abs(significance_alpha)-diffLL.apply(x))
							);
				}
				
			}
			
			
			double dmin = diffLL.apply(xmin);
			state_interval = param_cache.get(xmin);
			state_interval.resetModel();
			out.println("#**MLRV.fLI\txmin "+xmin+"\tmin "+dmin+"\tx0 "+x0+"\t// "+lrates.toString(node));
		}
		
		
		// sanity check nobody should gave a negative difference : that would be a better model than the starting one 
		double xbest=x0;
		double dbest=0.0;
		
		for (double x: diff_cache.keySet())
		{
			if (diff_cache.get(x)<dbest)
			{
				dbest = diff_cache.get(x);
				xbest = x;
			}
		}
		
		if (dbest<0.0) // hopefully not 
		{
			state_interval = param_cache.get(xbest);
			state_interval.resetModel();
			out.println("#*MLRV.fLI unexpected better model\txbest "+xbest+"\tbest "+dbest+"\tx0 "+x0+"\t// "+lrates.toString(node));
		}
				
		
		return -state_interval.getLL();
	}	
	
	/**
	 * Parses and sets truncation and min. observed copies 
	 * 
	 * @param out
	 * @param cli
	 */
	protected void parseComputeParameters(PrintStream out, CommandLine cli) {
		
    	boolean dup_bounded = cli.getOptionBoolean(CommandLine.OPT_DUP_BOUNDED, isDuplicationBounded()) ;
    	if (dup_bounded != isDuplicationBounded()) {
    		this.setDuplicationBounded(dup_bounded);
			out.println(CommandLine.getStandardHeader(
					"Bounded duplication: -"
					+CommandLine.OPT_DUP_BOUNDED+" "+dup_bounded
					));
    		
    	}
    	boolean gain_bounded = cli.getOptionBoolean(CommandLine.OPT_GAIN_BOUNDED, isGainBounded());
    	if (gain_bounded != isGainBounded()) {
    		this.setGainBounded(gain_bounded);
			out.println(CommandLine.getStandardHeader(
					"Bounded gain: -"
					+CommandLine.OPT_GAIN_BOUNDED+" "+gain_bounded
					));
    		
    	}
    	

		
		int absolute = Integer.MAX_VALUE;
		double relative = Double.POSITIVE_INFINITY;
        if (cli.getOptionValue(OPT_TRUNCATE)!=null)
        {
        	String truncate_val = cli.getOptionValue(OPT_TRUNCATE);
        	if (truncate_val.endsWith("auto"))
        	{
        		if ("auto".equals(truncate_val))
        		{
        			this.setWantAutoTruncation(true);
        			absolute = DEFAULT_TRUNCATE_ABSOLUTE;
        			relative = DEFAULT_TRUNCATE_RELATIVE;
        		} else if ("noauto".equals(truncate_val))
        		{
        			this.setWantAutoTruncation(false);
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
	        	this.setWantAutoTruncation(false);
        	}
        } 
        
        if (!this.auto_truncate)
        {
        	this.setCalculationWidth(absolute, relative);
        }
        
        out.println(CommandLine.getStandardHeader("Truncated computation: -"
        		+OPT_TRUNCATE+" "+absolute+","+relative)
        		+" ; auto-truncation "+(this.auto_truncate?"on":"off"));

        int min_copies = Integer.min(2, cli.getTable().minCopies());
        min_copies = cli.getOptionInt(OPT_MINCOPY, min_copies);
		this.setMinimumObservedCopies(min_copies);
		out.println(CommandLine.getStandardHeader("Minimum observed copies: -"+OPT_MINCOPY+" "+min_copies));
		boolean gotalgo = setAlgorithm(cli.getOptionValue(CommandLine.OPT_ALGORITHM));
		if (gotalgo)
			out.println(CommandLine.getStandardHeader("Numerical optimization algorithm: -"+CommandLine.OPT_ALGORITHM+" "
					+(USE_CONJUGATE_GRADIENT==null?"BFGS":USE_CONJUGATE_GRADIENT.toString())));
	}
	
	
    private static boolean setAlgorithm(String algo) {
//    	System.out.println("#**MLRV.sA "+algo);
    	boolean setAlgorithm=true;
		if (algo != null) {
			if ("BFGS".equalsIgnoreCase(algo)) {
				USE_CONJUGATE_GRADIENT=null;
			} else {
				CGVariant requested = CGVariant.valueOf(algo);
				if (requested == null) {
					if ("FR".equalsIgnoreCase(algo))
						requested = CGVariant.FletcherReeves;
					else if ("PR".equalsIgnoreCase(algo))
						requested = CGVariant.PolakRibiere;
					else if ("HS".equalsIgnoreCase(algo))
						requested = CGVariant.HestenesStiefel;
					else if ("DY".equalsIgnoreCase(algo))
						requested = CGVariant.DaiYuan;
				}
				if (setAlgorithm = (requested != null))
					USE_CONJUGATE_GRADIENT = requested;
			}
		}    	
		return setAlgorithm;
    }
	
	
	/**
	 * Main entry 
	 */
	public static void main(String[] args) throws Exception
	{
		PRINT_OPTIMIZATION_MESSAGES = true;
		Class<?> our_class = java.lang.invoke.MethodHandles.lookup().lookupClass();
		count.io.CommandLine cli = new count.io.CommandLine(args, our_class);

		PrintStream out = System.out;
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(our_class));
    	    out.println(CommandLine.getStandardRuntimeInfo(our_class, args));
    	}
    	
    	RateVariationModel model = null; 
    	TreeWithRates starting_rates;
		count.matek.DiscreteDistribution root_prior = cli.getOptionDistribution(OPT_MODEL_ROOT_PRIOR, null);
		
		Random RND = cli.getOptionRND(out);
		/*
		 * rate variation
		 */
		int dup_k = cli.getOptionInt(OPT_MODEL_DUPLICATION_CATEGORIES, 1);
		int length_k = cli.getOptionInt(OPT_MODEL_LENGTH_CATEGORIES, 1);
		int cat_k = cli.getOptionInt(OPT_MODEL_CATEGORIES, 1);
		
		boolean cold_start; 
    	if (cold_start = (cli.getVariationModel()==null))
    	{
    		starting_rates = new TreeWithRates(cli.getTree(), RND);
    		if (root_prior != null)
    		{
    			starting_rates.setRootDistribution(root_prior);
    			out.println(CommandLine.getStandardHeader("Root prior preset: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
    		}
    		starting_rates.initNodeParameters(starting_rates.getTree().getRoot());
			out.println(CommandLine.getStandardHeader("(Root prior random: "+starting_rates.getRootDistribution()+")"));
    		model = new RateVariationModel(starting_rates);
    		model.initConstantRates();
    		
    		starting_rates = model.getBaseModel();
    	} else
    	{
    		model = cli.getVariationModel();
    		starting_rates = model.getBaseModel();
    		if (root_prior != null)
    		{
    			starting_rates.setRootDistribution(root_prior);
    			out.println(CommandLine.getStandardHeader("Root prior preset: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
    		}
    		if (RND!=null)
    		{
    			starting_rates.setRandom(RND);
    		} else {
    			boolean reinit = cli.getOptionBoolean(CommandLine.OPT_REINIT, false);
    			if (reinit) {
    				TreeWithLogisticParameters lrates = (TreeWithLogisticParameters) starting_rates;
    				boolean adjusted = lrates.increaseZeroDuplicationRates();
    				if (adjusted) {
    					out.println(CommandLine.getStandardHeader("Zero-duplication edges reinitialized: -"+CommandLine.OPT_REINIT+" "+reinit));
    				}
    			}
    		}
    	}
    	
    	
		double proot = cli.getOptionDouble("opt.proot", 1.0);
		if (proot<1.0) {
			TreeWithLogisticParameters lrates = (TreeWithLogisticParameters) starting_rates;
			int root = cli.getTree().getRoot();
			double logitlm = lrates.getLogitRelativeRate(root);
			double log_gamma = lrates.getLogGainParameter(root, PARAMETER_LOSS);
			double logitp = Logarithms.toLogit(proot);
			lrates.setLogitLossRelativeDuplication(root, logitp, logitlm, log_gamma, PARAMETER_LOSS);
			
			out.println(CommandLine.getStandardHeader("Loss at root: -opt.proot "+proot));
		}
    	
		String opt_common_gain = CommandLine.OPT_COMMON_GAIN; 
		String common_str = cli.getOptionValue(opt_common_gain);
		int common_gain_by = CommandLine.parseOptionParameterType(common_str);
		if (common_gain_by!=-1)
		{
			model.setCommonGain(common_gain_by);
			out.println(CommandLine.getStandardHeader(
					"Common gain: -"
					+opt_common_gain+" "+common_str
					));
		}		
		if (model.isUniversalGain())
		{
			out.println(CommandLine.getStandardHeader("Model reset to linear gain (universal is not supported because it's buggy)"));
			model.setCommonGain(model.getCommonGainType(), false);
		}
		
		// bounded duplication rate (away from 1.0)
//		String opt_max_dup = "opt.dupbound";
//		double max_dup = cli.getOptionDouble(opt_max_dup, 1.0-MAX_DUPLICATION_MARGIN);
//		if (max_dup != 1.0-MAX_DUPLICATION_MARGIN)
//		{
//			MAX_DUPLICATION_MARGIN = 1.0-max_dup;
//		}
//		out.println(CommandLine.getStandardHeader("Bound on relative duplication rate: -"+opt_max_dup+" "+max_dup+"\t; (margin "+MAX_DUPLICATION_MARGIN+")"));
			
		
		// debugging
		String opt_debug_gradient = "debug.gradient";
		if (cli.getOptionValue(opt_debug_gradient)!=null)
		{
			boolean want_debug = cli.getOptionBoolean(opt_debug_gradient, false);
			DEBUG_GRADIENT = want_debug;
		}		
		
    	AnnotatedTable table = cli.getTable();
    	
//    	double avg_copy_number = table.getMeanCopies(true);
//    	System.out.println("#**MLRV.main table mean copies "+avg_copy_number);
    	
    	
    	MLRateVariation O = new MLRateVariation(model, table);
    	O.parseComputeParameters(out, cli);
    	
    	out.println(CommandLine.getStandardHeader("Unique profiles: "+O.utable.tableStatistics(starting_rates.getTree())));
    			
		int maxiter = cli.getOptionInt(OPT_ROUNDS, Integer.MAX_VALUE);
		double eps = cli.getOptionDouble(OPT_EPS, 1.0/(1L<<26));
        
        out.println(CommandLine.getStandardHeader("Convergence: -"+OPT_ROUNDS+" "+maxiter+" -"+OPT_EPS+" "+eps));

		String opt_dup0 = "opt.dup0";
		boolean has_poisson=cli.getOptionBoolean(opt_dup0, false);
		out.println(CommandLine.getStandardHeader("No-duplication class: -"+opt_dup0+" "+has_poisson));
        
//    	boolean dup_bounded = cli.getOptionBoolean(CommandLine.OPT_DUP_BOUNDED, O.isDuplicationBounded()) ;
//    	if (dup_bounded != O.isDuplicationBounded()) {
//    		O.setDuplicationBounded(dup_bounded);
//			out.println(CommandLine.getStandardHeader(
//					"Bounded duplication: -"
//					+CommandLine.OPT_DUP_BOUNDED+" "+dup_bounded
//					));
//    		
//    	}
        
		/*
		 * 
		 * Start optimization
		 */
    	int pretrain = cli.getOptionInt(CommandLine.OPT_PRETRAIN, 0);
    	cold_start = cold_start && (0<pretrain);

		if (cold_start 
				&& (1<cat_k || 1<dup_k*length_k))
		{
			out.println(CommandLine.getStandardHeader("Pretraining: -"+CommandLine.OPT_PRETRAIN+" "+pretrain));
			O.optimize(1.0/(1L<<20), Integer.min(pretrain, maxiter)); // just a little bit of adjustment for the base model before introducing multiple categories
		}
		
		String opt_ctype = "opt.ctype";
		String optClass = cli.getOptionValue(opt_ctype);
		if (optClass != null) {
			if (optClass.equals(RateVariationModel.Multiplier.class.getSimpleName())) {
				dup_k = 1;
				model.setDefaultType(RateVariationModel.Multiplier.class);
			} else if (optClass.equals(RateVariationModel.LogisticShift.class.getSimpleName())) {
				model.setDefaultType(RateVariationModel.LogisticShift.class);
			} else optClass = optClass+" (unrecognized: stay with default type)"; 
			
			out.println(CommandLine.getStandardHeader("Category type: -"+opt_ctype+" "+optClass));
		}
		
		
		
		
		// 
		// rate variation?
		// 
		if (1<dup_k*length_k)
		{
			
			
			
			out.println(CommandLine.getStandardHeader("Lattice init categories: -"+OPT_MODEL_LENGTH_CATEGORIES+" "+length_k+" -"+OPT_MODEL_DUPLICATION_CATEGORIES+" "+dup_k));			
			
			// 0.5 .. 2 
			double max_mod = 2.0;
			double mod_range = 2.0*Math.log(max_mod);
			
			// mid = floor(nK/2)
			// d = range / (nK-1)
			// mod(k) = (k-nK/2)*d;
			
			
			double delta_len = length_k==1?0.0:mod_range/(length_k-1.0);
			double delta_dup = dup_k==1?0.0:mod_range/(dup_k-1.0);
			
			// 
			int ncat=Integer.max(length_k, dup_k);
			
			
			double[] cat_mod_len = new double[ncat];
			double[] cat_mod_dup = new double[ncat];
			
			//		0-, 10, 2-, 30, 4-, 50, 6-, 70 
			//		0-, 10, 2+, 3-, 40, 5+, 6-, 70, + 
			//	2,3	     *
			//	4,5	        **
			//	6,7	            **
			//	8,9	                **
			
			
			int di=0;
			int li=0;
			for (int k=0; k<ncat; )
			{
				double mod_len = (li-length_k/2)*delta_len;
				double mod_dup = (di-dup_k/2)*delta_dup;
				cat_mod_len[k] = mod_len;
				cat_mod_dup[k] = mod_dup;
				
				// DEBUG
				k++;
				di = (di+1) % dup_k;
				li = (li+1) % length_k;
			}
			double[] cat_p = new  double[ncat];
			Arrays.fill(cat_p, 1.0/ncat);
			
			if (has_poisson) {
				int pcat = ncat-1;
				while (cat_mod_len[pcat]==0.0 && cat_mod_dup[pcat]==0.0) --pcat;
				cat_mod_dup[pcat]=Double.NEGATIVE_INFINITY;
				if (cat_mod_len[pcat]==0.0) cat_mod_len[pcat] = -delta_len;
			}
			
			model.initCategories(cat_p, cat_mod_len, cat_mod_dup);
			for (int k=0; k<ncat; k++)
			{
				System.out.println("#**MLRV.main initcat "+k+"\t"+model.getCategory(k)); 
			}
			O.initParameters();
		} else if (1<cat_k)
		{
			int ncat = model.getNumClasses();
			assert (1<=ncat);
			
			if (cat_k <= ncat)
			{
				model.initRandomCategories(cat_k, RND);
			} else 
			{
				for (int k=ncat; k<cat_k; k++)
					model.addRandomCategory(RateVariationModel.LogisticShift.class, RND);
			}
			O.initParameters();
		}
		
		
        int testpnode = cli.getOptionInt("testp", -1);
        int testrnode = cli.getOptionInt("testr", -1);
        int testqnode = cli.getOptionInt("testq", -1);
        boolean  want_parameter_test = (testpnode!=-1 || testrnode != -1 ||  testqnode != -1);
        double score;
        
        if (want_parameter_test)
        {
        	double pvalue = cli.getOptionDouble(OPT_PVALUE, 0.05);
        	int node, param_type;
        	if (testpnode != -1)
        	{
        		node = testpnode;
        		param_type = PARAMETER_LOSS;
        	} else if (testqnode != -1)
        	{
        		node = testqnode;
        		param_type = PARAMETER_DUPLICATION;
        	} else
        	{
        		node = testrnode;
        		param_type = PARAMETER_GAIN;
        	}
        	
        	
        	
        	double score0 = O.optimize(eps,0);
        	
        	System.out.println("#**MLRV.main likelihood interval "+pvalue+"\tstartscore "+score0);
        	
        	score = O.findLikelihoodInterval(node, param_type, pvalue, eps, maxiter);
        } else		
        {
        	score = O.optimize(eps, maxiter);
        }
        
        // DEBUG info
		if (maxiter==0 && DEBUG_GRADIENT)
		{
			O.debugGradient(out);
		}
        
        
        
		// model fit 

		double ascore = score/table.getFamilyCount();
		
		int npars = O.getModelParameterCount();
        double bic_pty = 0.5*npars*Math.log(table.getFamilyCount());

		out.println("#TREE "+NewickParser.printTree(cli.getTree()));
        out.println("#SCORE "+score+"\tBICpty "+bic_pty+"\tregularized "+(score+bic_pty)
        			+"\tnum.parameters "+npars+"\tsamplesize "+table.getFamilyCount());
        out.println("#AVGSCORE "+ascore);
        		
		// save model
		if (!want_parameter_test || out != System.out)
			out.println(count.io.RateVariationParser.printRates(model));

		if (out.checkError()) throw new java.io.IOException("Write failed.");
		out.close();
	}
	
}
