package count.model;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.DoubleConsumer;
import java.util.function.DoubleFunction;
import java.util.function.DoubleSupplier;

import count.ds.AnnotatedTable;
import count.ds.ProfileTable;
import count.io.CommandLine;
import count.io.NewickParser;
import count.matek.FunctionMinimization;

import static count.io.CommandLine.OPT_EPS;
import static count.io.CommandLine.OPT_MODEL_DUPLICATION_CATEGORIES;
import static count.io.CommandLine.OPT_MODEL_LENGTH_CATEGORIES;
import static count.io.CommandLine.OPT_MODEL_ROOT_PRIOR;
import static count.io.CommandLine.OPT_OUTPUT;
import static count.io.CommandLine.OPT_ROUNDS;
import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LOSS;
import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_LENGTH;
import static count.model.VariationGradientFactory.PARAMETER_MOD_DUPLICATION;
import static count.model.VariationGradientFactory.PARAMETER_MOD_LENGTH;

/**
 * Numerical optimization for Gamma model. 
 */
public class MLGamma extends MLRateVariation {
	public MLGamma(GammaInvariant model, ProfileTable table){
		super(RateVariationModel.convert(RateVariationModel.Multiplier.class, model), table);
		this.gamma_model = model;
		this.gamma_parameters = new GammaParameter[4];
		this.initParameters();
	}
	
	private final GammaInvariant gamma_model;
	
	private final GammaParameter[] gamma_parameters;
	
	private static final double MAX_ALPHA = 50.0;
	
	/**
	 * Whether Brent one-by-one or Powell with all params at once. 
	 * If only length variation, then brent is better.
	 */
	private static final boolean SEPARATE_GAMMA_OPTIMIZATION = true;
	
	
	private abstract class GammaParameter extends MinGradient {
		protected GammaParameter(int par, DoubleSupplier modelGetter, DoubleConsumer modelSetter, boolean doOptimize){
			super(3*node_parameters.length+par, par, true);
			this.modelGetter = modelGetter;
			this.modelSetter = modelSetter;
			this.isOptimized = doOptimize;
			copyFromModel();
		}
		
		protected final DoubleSupplier modelGetter;
		protected final DoubleConsumer modelSetter;
		protected boolean isOptimized; 
		
		protected double a; 
		
		@Override
		public double get() {
			return a; //fromLogistic(logit_a, MAX_ALPHA);
		}

		@Override
		public void set(double a) {
			this.a = a; //toLogistic(alpha, MAX_ALPHA);
		}
		
		protected abstract void copyFromModel();
		protected abstract void copyToModel();
		protected final boolean isOptimized() {return isOptimized;}
		
		/**
		 * No derivatives for gamma parameters
		 */
		@Override
		public double dL(double[] ignored) {
			return 0.0;
		}
		
	}
	
	/**
	 * Alpha shape parameter for Gamma rate multipliers: represented as logistic for bounding
	 */
	private class LogisticGammaParameter extends GammaParameter {
		LogisticGammaParameter(int par, DoubleSupplier modelGetter, DoubleConsumer modelSetter, boolean doOptimize){
			super(par,modelGetter,modelSetter,doOptimize);
		}
		
		@Override
		protected void copyFromModel() {
			double alpha = modelGetter.getAsDouble();
			set(toLogistic(alpha, MAX_ALPHA));
		}
		
		@Override
		protected void copyToModel() {
			double alpha = fromLogistic(get(), MAX_ALPHA);
			modelSetter.accept(alpha);
		}
		
		@Override 
		public String toString() {
			return "LGP["+GLDParameters.paramName(this.emt_idx)+"; logit "+a+"; alpha "+fromLogistic(a, MAX_ALPHA)+"; adjust "+isOptimized()+"]";
		}
		
	}
	
	private class DirectGammaParameter extends GammaParameter {
		DirectGammaParameter(int par, DoubleSupplier modelGetter, DoubleConsumer modelSetter, boolean doOptimize){
			super(par,modelGetter,modelSetter,doOptimize);
		}

		@Override
		protected void copyFromModel() {
			double alpha = modelGetter.getAsDouble();
			set(alpha);
		}
		@Override
		protected void copyToModel() {
			double alpha = get();
			modelSetter.accept(alpha);
		}
		
		@Override 
		public String toString() {
			return "DGP["+GLDParameters.paramName(this.emt_idx)+"; alpha "+get()+"; adjust "+isOptimized()+"]";
		}		
	}
	
	@Override
	protected void initCategoryParameters() {
		super.initCategoryParameters();
		
		int ncat = variation_model.getNumClasses();
		int delete = 2*ncat; // so many added by super's initCategoryParameters
		while (0<delete) {
			full_parameters.remove(full_parameters.size()-1); 
			--delete;
		}
		GammaParameter gainAlpha;
		GammaParameter lossAlpha;
		GammaParameter duplicationAlpha;
		GammaParameter lengthAlpha;
		
		if (SEPARATE_GAMMA_OPTIMIZATION) {
			gainAlpha = new DirectGammaParameter(PARAMETER_GAIN
					, ( )-> gamma_model.getGainAlpha()
					, (a)-> gamma_model.setGainAlpha(a)
					, false);
				lossAlpha = new DirectGammaParameter(PARAMETER_LOSS
					, ( )-> gamma_model.getLossAlpha()
					, (a)-> gamma_model.setLossAlpha(a)
					, false);
				duplicationAlpha = new DirectGammaParameter(PARAMETER_DUPLICATION
					, ( )-> gamma_model.getDuplicationAlpha()
					, (a)-> gamma_model.setDuplicationAlpha(a)
					, 1<gamma_model.getNumDuplicationGammaCategories());
				lengthAlpha = new DirectGammaParameter(PARAMETER_LENGTH
					, ( )-> gamma_model.getLengthAlpha()
					, (a)-> gamma_model.setLengthAlpha(a)
					, 1<gamma_model.getNumLengthGammaCategories());
			
		} else {
			gainAlpha = new LogisticGammaParameter(PARAMETER_GAIN
				, ( )-> gamma_model.getGainAlpha()
				, (a)-> gamma_model.setGainAlpha(a)
				, false);
			lossAlpha = new LogisticGammaParameter(PARAMETER_LOSS
				, ( )-> gamma_model.getLossAlpha()
				, (a)-> gamma_model.setLossAlpha(a)
				, false);
			duplicationAlpha = new LogisticGammaParameter(PARAMETER_DUPLICATION
				, ( )-> gamma_model.getDuplicationAlpha()
				, (a)-> gamma_model.setDuplicationAlpha(a)
				, 1<gamma_model.getNumDuplicationGammaCategories());
			lengthAlpha = new LogisticGammaParameter(PARAMETER_LENGTH
				, ( )-> gamma_model.getLengthAlpha()
				, (a)-> gamma_model.setLengthAlpha(a)
				, 1<gamma_model.getNumLengthGammaCategories());
		}
		
		gamma_parameters[PARAMETER_GAIN] = gainAlpha;
		gamma_parameters[PARAMETER_LOSS] = lossAlpha;
		gamma_parameters[PARAMETER_DUPLICATION] = duplicationAlpha;
		gamma_parameters[PARAMETER_LENGTH] = lengthAlpha;
		
		full_parameters.add(gainAlpha);
		full_parameters.add(lossAlpha);
		full_parameters.add(duplicationAlpha);
		full_parameters.add(lengthAlpha);
		
		// DEBUG
		System.out.println("#**MLG.iCP len "+lengthAlpha+"\tdup "+duplicationAlpha);
	}
	
	@Override 
	protected void initParameters() {
		if (gamma_model==null) return; // called from super's instantiation
		else super.initParameters();
	}
	
	@Override
	protected void initOptimizableCategoryParameters()
	{
		int num_nodes = node_parameters.length;
		this.do_optimize_parameters = Arrays.copyOf(do_optimize_parameters, 3*num_nodes+4);
		for (GammaParameter a : gamma_parameters) {
			this.do_optimize_parameters[a.pidx] = a.isOptimized();
		}
	}
	
	@Override
	protected void copyCategoryParametersFromModel() {
		GammaParameter da = gamma_parameters[PARAMETER_DUPLICATION];
		GammaParameter la = gamma_parameters[PARAMETER_LENGTH];
		
		boolean want_dup = this.do_optimize_parameters[da.pidx];
		boolean want_len = this.do_optimize_parameters[la.pidx];
		if (want_dup) {
			adjustable_parameters.add(da);
		}
		if (want_len) {
			System.out.println("#**MLG.cCPFM add length/"+la); // DEBUG
			adjustable_parameters.add(la);
		}
	}
	
	@Override
	protected void copyCategoryParametersToModel() {
		for (GammaParameter a: gamma_parameters) {
			a.copyToModel();
		}
		variation_model.updateCategories(gamma_model);
		// we also update the stored category parameter values
		int ncat = category_parameters.length;
		for (int k=0; k<ncat; k++)
		{
			RateVariationModel.Category C = variation_model.getCategory(k);
			double[] mods = category_parameters[k];	
			mods[PARAMETER_MOD_LENGTH] = C.getModLength();
			mods[PARAMETER_MOD_DUPLICATION] = C.getModDuplication();
		}
	}
	
	/**
	 * No recentering in Gamma model
	 */
	@Override
	protected boolean recenterBaseModel() {return false;}
	
	/**
	 * No update to class probabilities in Gamma model
	 */
	@Override
	protected OptimizationState optimizeCategoryProbabilities(OptimizationState current, double eps, List<Double> history) {
		if (PRINT_OPTIMIZATION_MESSAGES) {
			double max_dp = 0.0;
			int ncat = category_parameters.length;
			for (int k=0; k<ncat; k++)
			{
				RateVariationModel.Category C = variation_model.getCategory(k);
				double prior_p = Math.exp(C.getLogCatProbability());
				double log_post_p = current.G.getLogPosteriorFrequency(k);
				double post_p = Math.exp(log_post_p);
				
				double dp = Math.abs(post_p-prior_p);
				max_dp = Double.max(dp, max_dp);
				
				System.out.printf("#**MLG.pCProb cat %d\tprior %.4f\tpost %.4f\n", 
						k, prior_p, post_p);
			}
			System.out.println("#**MLG.pCProb max_dp "+max_dp); 
		}
		
		
		
		return current;
	}
	
	
	@Override
	protected void optimizeCategoryParams(double delta, int itmax, List<Double> history)
	{
		if (SEPARATE_GAMMA_OPTIMIZATION) {
			for (GammaParameter gpar: gamma_parameters) {
				if (gpar.isOptimized) {
					double a0 = gpar.get();
					DoubleFunction<Double> func = (x)->{
						gpar.set(x); 
						copyParametersToModel();
						double LL = gradient_factory.getCorrectedLL();
						
						if (PRINT_OPTIMIZATION_MESSAGES)
							System.out.println("#**MLG.oCP.f() "+gpar+"\t"+LL);
						return -LL;
					};
					double[] xmin_ymin = FunctionMinimization.brent(0.0, a0, MAX_ALPHA, 
							func, Double.max(Math.sqrt(Math.ulp(1.0)), delta));
					double amin = xmin_ymin[0];
					gpar.set(amin);
					copyParametersToModel();
					if (PRINT_OPTIMIZATION_MESSAGES)
						System.out.println("#**MLG.oCP "+gpar+"\tbest "+a0+" -> "+amin+"\tLL " +(-xmin_ymin[1])+"\tmodel "+gpar.modelGetter.getAsDouble());
					if (history!=null) history.add(xmin_ymin[1]);
				}
			}
		} else {
			super.optimizeCategoryParams(delta, itmax, history); // via Powell
		}
	}
			
	
	@Override
	public int getModelParameterCount() {
		int getModelParameterCount = super.getModelParameterCount();
		getModelParameterCount -= (variation_model.getNumClasses()-1); // fixed class prior probs
		return getModelParameterCount;
	}
	
	public void setCommonGainType(int common_gain) {
		variation_model.setCommonGain(common_gain,false);
		this.initParameters();
	}
	
	
	
	public void fixAlpha(int parameter_type, boolean do_not_optimize) {
		int num_nodes = node_parameters.length;
		this.do_optimize_parameters[3*num_nodes+parameter_type]=!do_not_optimize;
		gamma_parameters[parameter_type].isOptimized=!do_not_optimize;
		System.out.println("#**MLG.fixA "+GLDParameters.paramName(parameter_type)+"\t"+do_not_optimize+"\t"+gamma_parameters[parameter_type]);
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
    	RateVariationModel rvm_model = cli.getVariationModel(); 
    	GammaInvariant gamma_model = cli.getGammaModel();
    	TreeWithRates starting_rates;
		count.matek.DiscreteDistribution root_prior = cli.getOptionDistribution(OPT_MODEL_ROOT_PRIOR, null);
		
		Random RND = cli.getOptionRND(out);
		/*
		 * rate variation
		 */
		int dup_k = cli.getOptionInt(OPT_MODEL_DUPLICATION_CATEGORIES, 0);
		int length_k = cli.getOptionInt(OPT_MODEL_LENGTH_CATEGORIES, 0);
    	
		boolean cold_start; 
    	if (cold_start = (rvm_model==null))
    	{
    		starting_rates = new TreeWithRates(cli.getTree());
    		if (root_prior != null)
    		{
    			starting_rates.setRootDistribution(root_prior);
    			out.println(CommandLine.getStandardHeader("Root prior preset: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
    		}
			starting_rates.setRandom(RND);
    		if (RND != null) {
    			starting_rates.initNodeParameters();
    			out.println(CommandLine.getStandardHeader("(Root prior reset random: "+starting_rates.getRootDistribution()+")"));
    		}
			gamma_model = new GammaInvariant(starting_rates,1,1,Integer.max(dup_k,1),Integer.max(length_k,1));
    	} else {
    		starting_rates = rvm_model.getBaseModel();
    		if (gamma_model==null) {
        		starting_rates = rvm_model.getBaseModel();
        		if (RND!=null){
        			starting_rates.setRandom(RND);
        		}
        		if (root_prior != null){
        			starting_rates.setRootDistribution(root_prior);
        			out.println(CommandLine.getStandardHeader("Root prior preset: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
        		}
        		gamma_model = new GammaInvariant(starting_rates,1,1,Integer.max(dup_k,1),Integer.max(length_k,1));
    		}
    	}
    	if ((0<length_k && gamma_model.getNumLengthGammaCategories()!=length_k)
    		|| (0<dup_k && gamma_model.getNumDuplicationGammaCategories()!=dup_k)) {
    		
    		if (length_k==0) length_k = gamma_model.getNumLengthGammaCategories();
    		if (dup_k==0) dup_k = gamma_model.getNumDuplicationGammaCategories();
    		gamma_model.setClasses(1, 1, dup_k, length_k);
    	}
    	
    	
    	AnnotatedTable table = cli.getTable();
    	
    	MLGamma opt = new MLGamma(gamma_model, table);
    	opt.parseComputeParameters(out, cli);
    	
		String opt_common_gain = CommandLine.OPT_COMMON_GAIN; 
		String common_str = cli.getOptionValue(opt_common_gain);
		int common_gain_by = CommandLine.parseOptionParameterType(common_str);
		if (common_gain_by!=-1)
		{
			opt.setCommonGainType(common_gain_by);
			out.println(CommandLine.getStandardHeader(
					"Common gain: -"
					+opt_common_gain+" "+common_str
					));
		} else {
			out.println(CommandLine.getStandardHeader(
					"Common gain: -"
					+opt_common_gain+" "+GLDParameters.paramName(opt.variation_model.getCommonGainType())
					+"\t[unset]"
					));
			
		}
    	boolean dup_bounded = cli.getOptionBoolean(CommandLine.OPT_DUP_BOUNDED, opt.isDuplicationBounded()) ;
    	if (dup_bounded != opt.isDuplicationBounded()) {
    		opt.setDuplicationBounded(dup_bounded);
			out.println(CommandLine.getStandardHeader(
					"Bounded duplication: -"
					+CommandLine.OPT_DUP_BOUNDED+" "+dup_bounded
					));
    		
    	}
		
		
		int maxiter = cli.getOptionInt(OPT_ROUNDS, Integer.MAX_VALUE);
		double eps = cli.getOptionDouble(OPT_EPS, 1.0/(1L<<26));
        
        out.println(CommandLine.getStandardHeader("Convergence: -"+OPT_ROUNDS+" "+maxiter+" -"+OPT_EPS+" "+eps));

        
		double score = opt.optimize(eps, maxiter);
        
		// model fit 
		double ascore = score/table.getFamilyCount();
		
		int npars = opt.getModelParameterCount();
        double bic_pty = 0.5*npars*Math.log(table.getFamilyCount());

		out.println("#TREE "+NewickParser.printTree(cli.getTree()));
        out.println("#SCORE "+score+"\tBICpty "+bic_pty+"\tregularized "+(score+bic_pty)
        			+"\tnum.parameters "+npars+"\tsamplesize "+table.getFamilyCount());
        out.println("#AVGSCORE "+ascore);
        		
		// save model
        if (0<maxiter)
        	out.println(count.io.RateVariationParser.printRates(gamma_model));
		
		if (out.checkError())
			throw new java.io.IOException("Write failed.");
		out.close();    	    	
    	
    	
	}
}
