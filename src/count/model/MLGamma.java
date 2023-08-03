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


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import count.ds.AnnotatedTable;
import count.ds.ProfileTable;
import count.io.CommandLine;
import count.io.NewickParser;
import count.matek.FunctionMinimization;

import static count.io.CommandLine.OPT_EPS;
import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_MODEL_DUPLICATION_CATEGORIES;
import static count.io.CommandLine.OPT_MODEL_FORBIDDEN_DUPLICATION;
import static count.io.CommandLine.OPT_MODEL_FORBIDDEN_GAIN;
import static count.io.CommandLine.OPT_MODEL_GAIN_CATEGORIES;
import static count.io.CommandLine.OPT_MODEL_LENGTH_CATEGORIES;
import static count.io.CommandLine.OPT_MODEL_ROOT_PRIOR;
import static count.io.CommandLine.OPT_MODEL_UNIFORM_DUPLICATION;
import static count.io.CommandLine.OPT_MODEL_UNIFORM_GAIN;
import static count.io.CommandLine.OPT_OUTPUT;
import static count.io.CommandLine.OPT_RND;
import static count.io.CommandLine.OPT_ROUNDS;
import static count.io.CommandLine.OPT_TRUNCATE;

/**
 * Maximum likelihood for GammaInvariant model 
 * using rate and shape parameters. Not really robust: 
 * may croak with very large rate classes implied by skewed 
 * multipliers. 
 */  

public class MLGamma extends MLMixedRate
{
	private static final boolean POWELL_WITH_ALL_PARAMS = false;
	private static final boolean EDGEWISE_OPTIMIZE = false;
	private static final double MAX_ALPHA = 50.0;
	public MLGamma(GammaInvariant model, ProfileTable table)
	{
		super(model, table);
		this.model = model;
		this.class_params = new ArrayList<>();
	}
//	public MLGamma(GammaInvariant model, ProfileTable table, boolean[] optimize_rates)
//	{
//		super(model, table, optimize_rates);
//		this.model = model;
//		this.class_params = new ArrayList<>();
//	}
	
	private final GammaInvariant model;
	
	private final List<ModelParameter> class_params;
	
	private boolean fixed_gain=false;
	private boolean fixed_loss=false;
	private boolean fixed_duplication=false;
	private boolean fixed_length=false;
	private boolean fixed_forbidden_gain=false;
	private boolean fixed_forbidden_duplication=false;
	
	public GammaInvariant getModel()
	{
		return model;
	}
	public void fixGainAlpha()
	{
		this.fixed_gain = true;
	}
	public void fixLossAlpha()
	{
		this.fixed_loss = true;
	}
	public void fixDuplicationAlpha()
	{
		this.fixed_duplication = true;
	}
	public void fixLengthAlpha()
	{
		this.fixed_length=true;
	}
	public void fixForbiddenGain()
	{
		this.fixed_forbidden_gain=true;
	}
	public void fixForbiddenDuplication()
	{
		this.fixed_forbidden_duplication=true;
	}

	private void initModelParameters()
	{
		class_params.clear();
		
		if (!fixed_gain && model.getNumGainGammaCategories()>1)
			class_params.add(new Logistic(new GainAlpha(), MAX_ALPHA));
		if (!fixed_loss && model.getNumLossGammaCategories()>1)
			class_params.add(new Logistic(new LossAlpha(), MAX_ALPHA));
		if (!fixed_duplication && model.getNumDuplicationGammaCategories()>1)
			class_params.add(new Logistic(new DuplicationAlpha(), MAX_ALPHA));
		if (!fixed_length && model.getNumLengthGammaCategories()>1)
			class_params.add(new Logistic(new LengthAlpha(), MAX_ALPHA));
		if (!fixed_forbidden_gain && model.getGainForbidden()>0.0)
			class_params.add(new Logistic(new GainForbidden()));
		if (!fixed_forbidden_duplication && model.getDuplicationForbidden()>0.0)
			class_params.add(new Logistic(new DuplicationForbidden()));
		
		
	}
	
	@Override 
	public int getModelParameterCount()
	{
		return super.getModelParameterCount() + class_params.size();
	}
	
	private class GainAlpha implements ModelParameter
	{
		@Override
		public double get()
		{
			return model.getGainAlpha();
		}
		
		@Override 
		public void set(double a)
		{
			model.setGainAlpha(a);
		}
		@Override
		public String toString()
		{
			return "gA="+get();
		}
	}
	
	private class LossAlpha implements ModelParameter
	{
		@Override
		public double get()
		{
			return model.getLossAlpha();
		}
		
		@Override 
		public void set(double a)
		{
			model.setLossAlpha(a);
		}
		@Override
		public String toString()
		{
			return "lA="+get();
		}
	}
	private class DuplicationAlpha implements ModelParameter
	{
		@Override
		public double get()
		{
			return model.getDuplicationAlpha();
		}
		
		@Override 
		public void set(double a)
		{
			model.setDuplicationAlpha(a);
		}
		@Override
		public String toString()
		{
			return "dA="+get();
		}
	}
	private class LengthAlpha implements ModelParameter
	{
		@Override
		public double get()
		{
			return model.getLengthAlpha();
		}
		
		@Override 
		public void set(double a)
		{
			model.setLengthAlpha(a);
		}
		@Override
		public String toString()
		{
			return "tA="+get();
		}
	}
	private class GainForbidden implements ModelParameter
	{
		@Override
		public double get()
		{
			return model.getGainForbidden();
		}
		
		@Override 
		public void set(double f)
		{
			model.setGainForbidden(f);
		}
		@Override
		public String toString()
		{
			return "g0="+get();
		}
	}
//	private class LossForbidden implements ModelParameter
//	{
//		@Override
//		public double get()
//		{
//			return model.getLossForbidden();
//		}
//		
//		@Override 
//		public void set(double f)
//		{
//			model.setLossForbidden(f);
//		}
//	}
	private class DuplicationForbidden implements ModelParameter
	{
		@Override
		public double get()
		{
			return model.getDuplicationForbidden();
		}
		
		@Override 
		public void set(double f)
		{
			model.setDuplicationForbidden(f);
		}
		@Override
		public String toString()
		{
			return "d0="+get();
		}
	}
	
	public double[] get()
	{
		double[] x = new double[class_params.size()];
		for (int j=0; j<x.length; j++)
		{
			x[j] = class_params.get(j).get();
		}
		return x;
	}
	
	public void set(double[] x)
	{
		for (int j=0; j<x.length; j++)
		{
			class_params.get(j).set(x[j]);
		}
		gradient.computeClasses();
	}
	
	private int calls_optFunc=0;
	public double optClassFunc(double[] x)
	{
		set(x);
		double LL = gradient.getCorrectedLL();
		calls_optFunc++;
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#*MLG.oCF "+calls_optFunc+"\t"+LL+"\tL0 "+gradient.getUnobservedLL());
		
		return -LL;
	}
	
	/**
	 * Vector of rate + class parameters. 
	 * 
	 * @return array of length {@link #getModelParameterCount()}
	 */
	public double[] getAllParameterValues()
	{
		double[] x = new double[getModelParameterCount()];
		double[] rx = super.getParameterValues();
		double[] cx = this.get();
		assert (x.length == rx.length + cx.length);
		System.arraycopy(rx, 0, x, 0, rx.length);
		System.arraycopy(cx, 0, x, rx.length, cx.length);
		return x;
	}
	
	public void setAllParameterValues(double[] x)
	{
		int rn = super.getModelParameterCount();
		double[] rx = new double[rn];
		System.arraycopy(x, 0, rx, 0, rn);
		double[] cx = new double[x.length-rn];
		assert (cx.length == class_params.size());
		System.arraycopy(x, rn, cx, 0, cx.length);
		this.set(cx);
		super.setParameterValues(rx);
	}
	
	/**
	 * Negative log-likelihood (to be minimized) as a 
	 * function of all parameter values.
	 * 
	 * @param x vector of length {@link #getModelParameterCount()}
	 * @return
	 */
	public double optAllParameterFunc(double[] x)
	{
		int rn = super.getModelParameterCount();
		double[] rx = new double[rn];
		System.arraycopy(x, 0, rx, 0, rn);
		double[] cx = new double[x.length-rn];
		assert (cx.length == class_params.size());
		System.arraycopy(x, rn, cx, 0, cx.length);
		this.set(cx);
		return super.optRatesFunc(rx);
	}
	
	private static final int MAX_ITER=1000;
	
	@Override
	public double optimize(double delta)
	{
		return this.optimize(delta, MAX_ITER);
	}
	
	@Override
	public double optimize(double delta, int it_max)
	{
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#**MLG.o start "+delta+"\t"+it_max);
		
		initModelParameters();
		
		if (class_params.size()==0)
			return super.optimize(delta, it_max);

		super.optimize(0.0, 0); // setup parameters

		List<Double> history = getOptimizationHistory();
		if (history == null)
		{
			history = new ArrayList<>();
			setOptimizationHistory(history);
		}
		
		double LL =  gradient.getCorrectedLL();
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#**MLG.o starting LL "+LL+"\tclass params:"+class_params.size());

		double Ldiff = 1.0;
		double[] x0 ;
		double Ls = -LL;
		int iter = 0;
		
		
		if (it_max>0)
		{
			boolean first = true;
			do
			{
				double Lr;
				{
					// optimize only class parameters
					x0 = get();
					Lr = FunctionMinimization.powell(x0, delta, 3*class_params.size(), // quick adjustment of shape params 
							x->optClassFunc(x), history);
					set(x0);
					if (PRINT_OPTIMIZATION_MESSAGES)
						System.out.println("#*MLG.o "+iter+"\tpowell/classes "+Lr);
				}
				if (POWELL_WITH_ALL_PARAMS && first)
				{
					x0 = getAllParameterValues();
					Lr = FunctionMinimization.powell(x0, Math.sqrt(delta), 6, // adjustment of all params 
							x->optAllParameterFunc(x), history);
					setAllParameterValues(x0);
					if (PRINT_OPTIMIZATION_MESSAGES)
						System.out.println("#*MLG.o "+iter+"\tpowell/all "+Lr);
				}		
				if (EDGEWISE_OPTIMIZE && first)
				{
					double Le = super.edgewiseOptimize(Math.sqrt(delta), 2); // go through the edge parameters at most twice
					if (PRINT_OPTIMIZATION_MESSAGES)
						System.out.println("#*MLG.o "+iter+"\tedgewise "+Le);
					Lr = Le;
				}
				first = false;

				double Lc = super.optimize(delta, it_max-history.size()); // may reset params to respect bracketing
	
				if (Thread.currentThread().isInterrupted()) // keeps interrupt status
	            {
	                return Double.NaN;
	            }
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#*MLG.o "+iter+"\tdfpmin "+Lc);
				
//				x0 = get();
//				if (history.size()+2<it_max)
//				{
//					Lc = FunctionMinimization.powell(x0, delta, it_max-history.size(), 
//						x->optClassFunc(x), history);
//		            if (Thread.currentThread().isInterrupted()) // keeps interrupt status
//		            {
//		                return Double.NaN;
//		            }
//					System.out.println("#*MLG.o "+iter+"\tpowell "+Lc);
//				}
	
				++iter;
	
				Ldiff = (Ls-Lc); // /Ls;
				if (PRINT_OPTIMIZATION_MESSAGES)
					System.out.println("#*MLG.o "+iter+"\tLd "+Ldiff+"\tLs "+Ls+"\tLr "+Lr); //+"\tLc "+Lc);
				//assert (Ldiff >=  -delta); //0.0);
				Ls = Lc;
			} while (Ldiff>delta && history.size()+2<it_max && iter<it_max);
	
			LL =  gradient.getCorrectedLL();
			// debug
			{
//				gradient.reportClassParameters();
//				reportDistributionGradient();
//				
//				double x[] = getParameterValues();
//				double D[] = optRatesDiff(x);
//				double gradient_length = FunctionMinimization.euclideanNorm(D);
//				double param_length = FunctionMinimization.euclideanNorm(x);
//				double rel_length = gradient_length/param_length;
//				
//				// decreasing order by absolute value of gradient 
//				List<Integer> sorted_params = new ArrayList<>();
//				for (int pidx=0; pidx<super.getModelParameterCount(); pidx++) // super only has rate params
//					sorted_params.add(pidx);
//		        sorted_params.sort(new java.util.Comparator<>(){
//		        	@Override
//		        	public int compare(Integer p1, Integer p2)
//		        	{
//		        		double d1 = Math.abs(D[p1]);
//		        		double d2 = Math.abs(D[p2]);
//		        		return Double.compare(d2, d1);
//		        	}
//		        });
//		        StringBuilder grad_sb = new StringBuilder();
//		        for (int pidx=0; pidx<3; pidx++)
//		        {
//		        	int j = sorted_params.get(pidx);
//		        	double rel = Math.abs(D[j]/gradient_length);
//		        	grad_sb.append(";\t")
//		        	.append(getRateParameter(j))
//		        	.append(": dLL/dx ")
//		        	.append(D[j])
//		        	.append(" (")
//		        	.append(rel)
//		        	.append(")");
//		        }
//		        System.out.println("#*MLD.o	gradient\t"+gradient_length+"\t("+rel_length+")"+grad_sb);        
//		        
//		        
//			
				if (PRINT_OPTIMIZATION_MESSAGES)
				{
					double L0 = gradient.getEmptyLL();
					double L1 = gradient.getSingletonLL();
					double Lo = gradient.getObservedLL();
					
					System.out.println("#**MLG.o final LL "+LL+"\tL0 "+L0+"\tL1 "+L1+"\tLobs "+Lo
							+"\t// "+Math.exp(L0)+"\t"+Math.exp(L1)+"\t"+Math.exp(Lo));
				}
			}
		} 
		return -LL;
	}
	
	public static void main(String[] args) throws Exception
	{
		PRINT_OPTIMIZATION_MESSAGES=true;

		count.io.CommandLine cli = new count.io.CommandLine(args, MLGamma.class);
		
		PrintStream out = System.out;
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(MLGamma.class));
    	    out.println(CommandLine.getStandardRuntimeInfo(MLGamma.class, args));
    	}
    	AnnotatedTable table = cli.getTable();

		int absolute = 12;
		double relative = 3.0;
        if (cli.getOptionValue(OPT_TRUNCATE)!=null)
        {
        	absolute = cli.getOptionTruncateAbsolute();
        	relative = cli.getOptionTruncateRelative();
        } 
        out.println(CommandLine.getStandardHeader("Truncated computation: -"
        		+OPT_TRUNCATE+" "+absolute+","+relative));
        int min_copies = Integer.min(2, table.minCopies());
        min_copies = cli.getOptionInt(OPT_MINCOPY, min_copies);
		out.println(CommandLine.getStandardHeader("Minimum observed copies: -"+OPT_MINCOPY+" "+min_copies));
		int maxiter = cli.getOptionInt(OPT_ROUNDS, Integer.MAX_VALUE);
		double eps = cli.getOptionDouble(OPT_EPS, 1e-4);
        
        out.println(CommandLine.getStandardHeader("Convergence: -"+OPT_ROUNDS+" "+maxiter+" -"+OPT_EPS+" "+eps));

		
    	
    	GammaInvariant model = null; 
    	TreeWithRates base_rates;
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
    		base_rates = model.getBaseModel();
    	} else
    	{
    		base_rates = cli.getRates(); 
    		model = cli.getModel();
    	}
    	
    	int gain_k = cli.getOptionInt(OPT_MODEL_GAIN_CATEGORIES, model.getNumGainGammaCategories());
    	int length_k = cli.getOptionInt(OPT_MODEL_LENGTH_CATEGORIES, model.getNumLengthGammaCategories());
    	int dup_k = cli.getOptionInt(OPT_MODEL_DUPLICATION_CATEGORIES, model.getNumDuplicationGammaCategories());
    	boolean forbidden_gain = cli.getOptionBoolean(OPT_MODEL_FORBIDDEN_GAIN, model.getGainForbidden()>0.0);
    	boolean forbidden_duplication = cli.getOptionBoolean(OPT_MODEL_FORBIDDEN_DUPLICATION, model.getDuplicationForbidden()>0.0);
    	boolean uniform_gain = cli.getOptionBoolean(OPT_MODEL_UNIFORM_GAIN, false);
    	boolean uniform_duplication = cli.getOptionBoolean(OPT_MODEL_UNIFORM_DUPLICATION, false);
    	out.println(CommandLine.getStandardHeader("Categories: -"
    			+OPT_MODEL_GAIN_CATEGORIES+" "+gain_k+" -"
    			+OPT_MODEL_DUPLICATION_CATEGORIES+" "+dup_k+" -"
    			+OPT_MODEL_LENGTH_CATEGORIES+" "+length_k+" -"
    			+OPT_MODEL_FORBIDDEN_GAIN+" "+forbidden_gain+" -"
    			+OPT_MODEL_FORBIDDEN_DUPLICATION + " " + forbidden_duplication
    			+" -"+OPT_MODEL_UNIFORM_GAIN+" "+uniform_gain
    			+" -"+OPT_MODEL_UNIFORM_DUPLICATION+" "+uniform_duplication
    			));
		count.matek.DiscreteDistribution root_prior = cli.getOptionDistribution(OPT_MODEL_ROOT_PRIOR, null);
		if (root_prior != null)
		{
			base_rates.setRootDistribution(root_prior);
			out.println(CommandLine.getStandardHeader("Root prior set: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
		}
		out.println(CommandLine.getStandardHeader("(Root prior get: "+base_rates.getRootDistribution()+")"));
    	
    	
    	if (cli.getModel()==null)
    	{
    		// optimize the base model first
    		ML Obase;
    		if (!uniform_duplication && !uniform_gain)
    		{
	    		Obase = new MLDistribution(base_rates, table);
    		} else
    		{
    			MLGamma Og = new MLGamma(model, table);
    			Og.setUniformDuplication(uniform_duplication);
    			Og.setUniformGain(uniform_gain);
    			Obase = Og;
    		}
    		Obase.setMinimumObservedCopies(min_copies);
    		Obase.setCalculationWidth(absolute, relative);
    		out.println("#** MLG.main: optimizing base model with "+Obase.getClass().getName());
    		double bscore = Obase.optimize(eps, maxiter);
    	}
    	
    	model.setClasses(gain_k, model.getNumLossGammaCategories(), dup_k, length_k);
    	{ // 0-rate classes 
    		double g0 = model.getGainForbidden();
    		if (forbidden_gain && g0==0.0)
    			model.setGainForbidden(1e-5); // little bit bigger than 0
    		else if (!forbidden_gain)
    			model.setGainForbidden(0.0);
    		double d0 = model.getDuplicationForbidden();
    		if (forbidden_duplication && d0 == 0.0)
    			model.setDuplicationForbidden(1e-5);
    		else if (!forbidden_duplication)
    			model.setDuplicationForbidden(0.0);
    	}
    	
		MLGamma O = new MLGamma(model, table);
		O.setUniformDuplication(uniform_duplication);
		O.setUniformGain(uniform_gain);
    	O.setCalculationWidth(absolute, relative);
		O.setMinimumObservedCopies(min_copies);
		
		
		
//		O.setMinimumObservedCopies(0);
//		O.setMinimumObservedCopies(1);
		double score = O.optimize(eps, maxiter);
        double bic_pty = 0.5*O.getModelParameterCount()*Math.log(table.getFamilyCount());

		out.println("#SCORE "+score+"\tBICpty "+bic_pty+"\tregularized "+(score+bic_pty));
        
		out.println("#TREE "+NewickParser.printTree(cli.getTree()));
		out.println(count.io.RateVariationParser.printRates(model));
		if (out != System.out)
		{
			if (out.checkError())
				throw new java.io.IOException("Write failed.");
			out.close();
		}
	}
	
	

}
