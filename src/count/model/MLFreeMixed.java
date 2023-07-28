package count.model;

import count.ds.AnnotatedTable;
import count.ds.Phylogeny;
import count.ds.ProfileTable;
import count.ds.UniqueProfileTable;
import count.io.CommandLine;
import count.io.RateVariationParser;
import count.matek.FunctionMinimization;

import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LOSS;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static count.io.CommandLine.OPT_EPS;
import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_MODEL_CATEGORIES;
import static count.io.CommandLine.OPT_OUTPUT;
import static count.io.CommandLine.OPT_RND;
import static count.io.CommandLine.OPT_ROUNDS;
import static count.io.CommandLine.OPT_TRUNCATE;

/**
 * 
 * @author csuros
 * @deprecated
 */
public class MLFreeMixed extends ML
{
	public MLFreeMixed(MixedRateModel.RateMultipliers
			mixed_model, ProfileTable table)	
	{
		this(mixed_model, table, MLDistribution.optimizableParams(mixed_model.getBaseModel()));
	}
	public MLFreeMixed(MixedRateModel
			mixed_model, ProfileTable table)	
	{
		this(mixed_model, table, MLDistribution.optimizableParams(mixed_model.getClassModel(0)));
	}
	public MLFreeMixed(MixedRateModel mixed_model, ProfileTable table, boolean[] optimize_params)	
	{
		if (mixed_model instanceof FreeMixedModel) 
		{
			this.model = (FreeMixedModel) mixed_model;
		}
		else
		{
			this.model = new FreeMixedModel(mixed_model);
		}
			
		if (table instanceof UniqueProfileTable)
		{
			this.utable = (UniqueProfileTable)table;
		} else
		{
			this.utable = new UniqueProfileTable(table);
		}
		this.gradient = new MixedDistributionGradient(model, utable);

		this.optimize_parameters = optimize_params;
		this.distribution_params = new ArrayList<>();
		this.class_params = new ArrayList<>();
		this.min_copies = Integer.min(2, table.minCopies());
	}
	
	
	
	private final FreeMixedModel model;
	private final UniqueProfileTable utable;
	
	public FreeMixedModel getModel()
	{
		return model;
	}

	private MixedDistributionGradient gradient;
	
	private int min_copies;
	
	public void setMinimumObservedCopies(int min_copies)
	{
		if (min_copies<0 || min_copies>2)
			throw new UnsupportedOperationException("Minimum copies must be 0,1 or 2.");
		this.min_copies = min_copies;
	}
	
	private int threshold_width_absolute = Integer.MAX_VALUE;
	private double threshold_width_relative = Double.POSITIVE_INFINITY;
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
		this.threshold_width_absolute = absolute;
		this.threshold_width_relative = relative;
	}
	
	
	private double[][] node_parameters;
	
	private final List<ModelParameter> distribution_params;
	private final List<ClassModelParameter> class_params;
	private boolean[] optimize_parameters;
	
	
	public void fixGain(int node, boolean not_optimized)
	{
		optimize_parameters[3*node+PARAMETER_GAIN]=!not_optimized;
	}
	public void fixLoss(int node, boolean not_optimized)
	{
		optimize_parameters[3*node+PARAMETER_LOSS]=!not_optimized;
	}
	public void fixDuplication(int node, boolean not_optimized)
	{
		optimize_parameters[3*node+PARAMETER_DUPLICATION]=!not_optimized;
	}	
	
	private boolean wantEdge(int node)
	{
		int j = 3*node;
		return optimize_parameters[j+PARAMETER_GAIN]
				|| optimize_parameters[j+PARAMETER_LOSS]
				|| optimize_parameters[j+PARAMETER_DUPLICATION];
	}
	
	private class ClassModelParameter implements ModelParameter
	{
		ClassModelParameter(int cat, int node, int par_idx)
		{
			this.cat = cat;
			this.node = node;
			this.offset = 3*node+par_idx;
		}
		final int cat;
		final int node;
		private final int offset;
		
		@Override
		public double get() 
		{
			return node_parameters[cat][offset];
		}
		@Override
		public void set(double x) 
		{
			node_parameters[cat][offset] = x;
		}
		@Override
		public double dL(double[] gradient)
		{
			return gradient[offset];
		}		
	}
	
	
	private class GainParameter extends ClassModelParameter
	{
		GainParameter(int cat, int node)
		{
			super(cat, node, PARAMETER_GAIN);
			set(model.getClassModel(cat).getGainParameter(node));
		}
		@Override
		public String toString()
		{
			return "r"+node+"/"+cat+"="+get();
		}
	}
	private class LossParameter extends ClassModelParameter
	{
		LossParameter(int cat, int node)
		{
			super(cat, node, PARAMETER_LOSS);
			set(model.getClassModel(cat).getLossParameter(node));
		}
		@Override
		public String toString()
		{
			return "p"+node+"/"+cat+"="+get();
		}
	}
	private class DuplicationParameter extends ClassModelParameter
	{
		DuplicationParameter(int cat, int node)
		{
			super(cat, node, PARAMETER_DUPLICATION);
			TreeWithRates rates = model.getClassModel(cat);
			set(rates.getDuplicationParameter(node));
		}
		@Override
		public String toString()
		{
			return "q"+node+"/"+cat+"="+get();
		}
	}
	
	private Logistic addLogistic( ClassModelParameter theta, double max)
	{
//		double θ = theta.get();
//		if (θ>max)
//		{
//			double newθ = max;
//			System.out.println("#*MLMD.addL "+theta+"; resetting to "+newθ);
//			theta.set(newθ);
//		}
		Logistic L = newLogistic(theta, max);
		distribution_params.add(L);
		class_params.add(theta);
		return L;
	}	
	
	private Logarithmic addLogarithmic(ClassModelParameter theta)
	{
		Logarithmic L = new Logarithmic(theta);
		distribution_params.add(L);
		class_params.add(theta);
		return L;
	}
	
	private void initModelParameters()
	{
		distribution_params.clear();
		class_params.clear();
		int num_nodes = -1;
		node_parameters = new double[model.getNumClasses()][];
		
		for (int c=0; c<model.getNumClasses(); c++)
		{
			TreeWithRates rates = model.getClassModel(c);
			if (num_nodes == -1)
				num_nodes = rates.getTree().getNumNodes();
			else
			{
				assert(num_nodes == rates.getTree().getNumNodes()); 
			}
			node_parameters[c] = new double[3*num_nodes];
			for (int node=0; node<num_nodes; node++)
			{
				int j = 3*node;
				if (wantEdge(node))
				{
					boolean want_gain = optimize_parameters[j+PARAMETER_GAIN];
					boolean want_loss = optimize_parameters[j+PARAMETER_LOSS];
					boolean want_duplication = optimize_parameters[j+PARAMETER_DUPLICATION];
					// instantiate all 3 parameters to set node_parameters[c][j+..]
					GainParameter gain = new GainParameter(c, node);
					LossParameter loss = new LossParameter(c, node);
					DuplicationParameter dup = new DuplicationParameter(c, node);
					//
					if (want_gain)
						addLogarithmic(gain);// 
					if (want_loss)
						addLogistic(loss, MAX_PROB_NOT1);
					if (want_duplication)
						addLogistic(dup, MAX_PROB_NOT1);
				}
			}
		}
//		for (int j=0; j<distribution_params.size(); j++)
//		{
//			System.out.println("#**MLFM.iMP "+j+"\t"+distribution_params.get(j)+"\tcp "+class_params.get(j));
//		}
		gradient.setCalculationWidthThresholds(threshold_width_absolute, threshold_width_relative);;
		gradient.setMinimumObservedCopies(min_copies);
		gradient.computeClasses();
	}
	
	/**
	 * Sets the actual rates at all nodes.
	 */
	private void copyNodeParametersToModel()
	{
		for (int cat=0; cat<model.getNumClasses(); cat++)
		{
			TreeWithRates rates = model.getClassModel(cat); 
			int num_nodes = rates.getTree().getNumNodes();
			for (int node=0; node<num_nodes; node++)
			{
				if (wantEdge(node))
				{
					double κ = node_parameters[cat][3*node+PARAMETER_GAIN];
					double p = node_parameters[cat][3*node+PARAMETER_LOSS];
					double q = node_parameters[cat][3*node+PARAMETER_DUPLICATION];
					rates.setParameters(node, κ, p, q);	
//					System.out.println("#**MLD.cNPTM "+node+"/"+cat+"\t"+
//							Arrays.toString(Arrays.copyOfRange(node_parameters[cat],3*node, 3*node+3))
//							+"\t"+rates.toString(node)
//							+"\t// "+κ+", "+p+", "+q);
				}
			}
		}
	}	
	
	public double[] getDistributionParameterValues()
	{
		double[] x = new double[distribution_params.size()];
		for (int j=0; j<x.length; j++)
		{
			x[j] = distribution_params.get(j).get();
		}
		return x;
	}
	
	public void setDistributionParameterValues(double[] x)
	{
		for (int j=0; j<distribution_params.size(); j++)
		{
			ModelParameter P = distribution_params.get(j);
			P.set(x[j]);
		}
		copyNodeParametersToModel();
		gradient.computeClasses();
	}
	
	
	private int calls_optFunc=0;
	private double optDistributionFunc(double[] x)
	{
		setDistributionParameterValues(x);
		double LL = gradient.getCorrectedLL();
		++calls_optFunc;
		return -LL;
	}
	
	
	private int calls_optDiff = 0;
	private double[] optDistributionGradient(double[] x)
	{
		setDistributionParameterValues(x);
		double[][] classD = gradient.getCorrectedGradient();
		double[] dLL = new double[x.length];
		for (int j=0; j<x.length; j++)
		{
			int cat = class_params.get(j).cat;
			dLL[j] = -distribution_params.get(j).dL(classD[cat]);
		}
		++calls_optDiff;
		System.out.println("#*MLMR.oRD "+calls_optDiff+"\tcallF "+calls_optFunc+"\tLL "+gradient.getCorrectedLL()
				+"\tdLL "+FunctionMinimization.euclideanNorm(dLL));
//		if (calls_optDiff % 10==1)
//		{
//			for (int j=0; j<x.length; j++)
//			{
//				System.out.println("#**MLMR.oRD "+calls_optDiff
//						+"\t"+distribution_params.get(j)
//						+"\t"+dLL[j]);
//			}
//		}
		
		return dLL;
	}
	
	private double optEMClassProbabilities()
	{
		double[] post = gradient.getEmpiricalClassProbabilities();
		model.setClassProbabilities(post);
		gradient.computeClasses();
		double LL = gradient.getCorrectedLL();
		return -LL;
	}
	
	public int getParameterCount()
	{
		return distribution_params.size()+model.getNumClasses()-1;
	}
	
	@Override
	public double optimize(double delta, int itmax) 
	{
		
		initModelParameters();
		
		List<Double> opt_history = new ArrayList<>();
		
		double LL =  gradient.getCorrectedLL();
		System.out.println("#**MLFM.o starting LL "+LL+"\tparams "+getParameterCount());
		opt_history.add(-LL);
		
		
		while (itmax>0)
		{
			double[] x0 = getDistributionParameterValues();
			int histr_len = opt_history.size();
			double min = FunctionMinimization.dfpmin(x0, delta, itmax, x->optDistributionFunc(x), x->optDistributionGradient(x), opt_history);
			int dfp_steps = opt_history.size()-histr_len;
			itmax -= dfp_steps;
			
			LL = gradient.getCorrectedLL();
			System.out.println("#**MLFM.o opt/distr LL "+LL+"\tmin "+min);
			
			double min2 = optEMClassProbabilities();
			LL = gradient.getCorrectedLL();
			System.out.println("#**MLFM.o opt/distr LL "+LL+"\tmin2 "+min2);
			System.out.println("#**MLFM.o classes: "+Arrays.toString(model.getClassProbabilities()));
			
			if (min2>min-Math.sqrt(delta)) break;
		}
		
		return -LL;
	}
	
	@Override
	public double optimize(double delta)
	{
		return this.optimize(delta, Integer.MAX_VALUE);
	}
	
	public static void main(String[] args) throws Exception
	{

		CommandLine cli = new CommandLine(args, MLFreeMixed.class);
		AnnotatedTable table = cli.getTable();
		FreeMixedModel input_model = cli.getFreeModel();
		
		PrintStream out = System.out;
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(MLFreeMixed.class));
    	    out.println(CommandLine.getStandardRuntimeInfo(MLFreeMixed.class, args));
    	}
		
		if (cli.getOptionValue(OPT_MODEL_CATEGORIES)!=null)
		{
			int k = cli.getOptionInt(OPT_MODEL_CATEGORIES, 2);
			TreeWithRates rates = cli.getRates();
			
			if (cli.getOptionValue(OPT_RND)==null)
			{
				
				int dup_k = 1;
				int len_k = k;
				while (len_k>4*dup_k && (len_k % 2)==0)
				{
					dup_k *= 2;
					len_k /= 2;
				}
				
				//len_k = k; dup_k=1;
				System.out.println("#**MLFM.main init k "+k+"\tdupk "+dup_k+"\tlenk "+len_k);
				
				GammaInvariant model_k = new GammaInvariant(rates, 1, 1, dup_k, len_k);
				input_model = new FreeMixedModel(model_k);
			} else
			{
				Random RND = cli.getOptionRND(out);
//				long seed = cli.getOptionLong(OPT_RND, 0L);
//				if (seed == 0L) RND = new Random();
//				else RND = new Random(seed);
				
				TreeWithRates[] class_rates = new TreeWithRates[k];
				for (int c=0; c<k; c++)
				{
					if (c==0) class_rates[c] = rates;
					else
					{
						TreeWithRates mod_rates = new TreeWithRates(rates);
						// random permutation
						int[] perm = new int[mod_rates.getTree().getNumNodes()-1];
						for (int i=0; i<perm.length; i++) perm[i] = i;
						for (int i=0; i<perm.length-1;i++)
						{
							int j = i+RND.nextInt(perm.length-i);
							if (i != j)
							{
								int pj = perm[j];
								perm[j] = perm[i];
								perm[i]=pj;
							}
						}
						for (int node=0; node<=perm.length; node++)
						{
							double p =rates.getLossParameter(node);
							double q = rates.getDuplicationParameter(node);
							double r = rates.getGainParameter(node);
							
							double modp = p==1.0?p:p*RND.nextDouble();
							double modq = q*RND.nextDouble();
							double modr = r*RND.nextDouble();
							
							int mnode = node==perm.length?node:perm[node];
							mod_rates.setParameters(mnode, modr, modp, modq);
						}
						class_rates[c] = mod_rates;
					}
				}
				input_model = new FreeMixedModel(class_rates, null);		
			}
		}
		
		MLFreeMixed M = new MLFreeMixed(input_model, table);
        if (cli.getOptionValue(OPT_TRUNCATE)!=null)
        {
        	int absolute = cli.getOptionTruncateAbsolute();
        	double relative = cli.getOptionTruncateRelative();
        	M.setCalculationWidth(absolute, relative);
        	out.println(CommandLine.getStandardHeader("Truncated computation: -"+OPT_TRUNCATE+" "
        				+absolute+","+relative));
        }
        if (cli.getOptionValue(OPT_MINCOPY)!=null)
        {
        	int min_copies = cli.getOptionInt(OPT_MINCOPY, table.minCopies());
        	M.setMinimumObservedCopies(min_copies);
        	out.println(CommandLine.getStandardHeader("Minimum observed: -"+OPT_MINCOPY+" "
    				+min_copies));
        }
		int maxiter = cli.getOptionInt(OPT_ROUNDS, Integer.MAX_VALUE);
		double eps = cli.getOptionDouble(OPT_EPS, 0.00125);
        
        out.println(CommandLine.getStandardHeader("Convergence: -"+OPT_ROUNDS+" "+maxiter+" -"+OPT_EPS+" "+eps));

        double score = M.optimize(eps, maxiter);
        double bic_pty = 0.5*M.getParameterCount()*Math.log(table.getFamilyCount());
        double bic1 = 0.5*(3.0*cli.getTree().getNumNodes()-1.0)*Math.log(table.getFamilyCount());
		System.out.println("#SCORE "+score+"\tBIC/1class "+(score+(bic_pty-bic1))+"\tbic "+bic_pty+"\tbic1 "+bic1+"\tregularized "+(score+bic_pty));
		System.out.println("#**MLFM.main elapsed time: "+cli.getMillisSinceStart()/1000.0+" sec");
        
        out.print(RateVariationParser.printRates(input_model));
        
		if (out != System.out)
		{
			if (out.checkError())
				throw new java.io.IOException("Write failed.");
			out.close();
		}
		
	}
	
	
}
