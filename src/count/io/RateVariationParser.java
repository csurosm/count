/*
 * Copyright 2021 Mikl&oacute;s Cs&#369;r&ouml;s.
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
package count.io;

import java.io.Reader;
import java.util.ArrayList;
import java.util.List;
import java.io.BufferedReader;
import java.io.IOException;

import count.ds.TreeTraversal;
import count.ds.IndexedTree;
import count.matek.DiscreteDistribution;
import count.matek.NegativeBinomial;
import count.matek.PointDistribution;
import count.matek.Poisson;
import count.matek.ShiftedGeometric;
import count.model.FreeMixedModel;
import count.model.GammaInvariant;
import count.model.TreeWithRates;

/**
 * Static methods for reading and writing rate variation models. 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class RateVariationParser 
{
    /**
     * The strings used in the rate file to mark the rate variation parameters
     */
    static String RATE_VARIATION_PREFIX = "|variation";
    static String ROOT_PRIOR_PREFIX = "|root";
    static String MODEL_END = "|endrates"; // for free mixed model
    
	private RateVariationParser() {}
	
//	public static void initFromFile(TreeWithRates rates, BufferedReader input) throws IOException
//	{
//		initFromFile(rates, input, true);
//	}
//	private static void initFromFile(TreeWithRates rates, BufferedReader input, boolean compatibility) throws IOException
//	{
//		IndexedTree tree = rates.getTree();
//		int num_nodes = tree.getNumNodes();
//		final int[] postorder = TreeTraversal.postOrder(tree); // legacy logic: input file lists the rates in postorder
//		
//		int stop_offset = compatibility?1:0; // not for root in legacy format
//		
//		for (int i=0; i<num_nodes-stop_offset; i++) 
//		{
//			String line = null;
//	        do
//	        {
//	            line = input.readLine();
//	        } while (line != null && line.startsWith("#")); // skip comments
//	
//	        int node_idx = postorder[i];
//	        String[] fields=line.split("\\s+");
//	        double len = Double.parseDouble(fields[0]);
//	        double drate =  Double.parseDouble(fields[1]);
//	        double lrate = Double.parseDouble(fields[2]);
//	        double grate = Double.parseDouble(fields[3]);
//	        
//	        rates.setEdgeLength(node_idx, len);
//	        rates.setDuplicationRate(node_idx, drate);
//	        rates.setLossRate(node_idx, lrate);
//	        if (drate==0.0)
//	        {
//	        	if (lrate==0.0)
//	    	        rates.setGainRate(node_idx, grate);
//	        	else
//	        		rates.setGainRate(node_idx, grate/lrate); // legacy logic: store κ*μ
//	        } else
//	        {
//	        	rates.setGainRate(node_idx, grate/drate); // legacy logic: store κ*λ
//	        }
//		}
//	}
	
	private static String last_rate_line= null; 
	
	public static GammaInvariant readRates(BufferedReader input, IndexedTree tree) throws FileFormatException, IOException	
	{
		final TreeWithRates base_rates = new TreeWithRates(tree);
		final GammaInvariant readRates = new GammaInvariant(base_rates);
		
		
		// parse the root prior and the gamma/forbidden rate categories
		int gain_categories = 1;
		int loss_categories = 1;
		int duplication_categories = 1;
		int length_categories = 1;
		
		double gain_alpha = 1.0;
		double loss_alpha = 1.0;
		double duplication_alpha = 1.0;
		double length_alpha = 1.0;
		
		double gain_forbidden = 0.0;
		double loss_forbidden = 0.0;
		double duplication_forbidden = 0.0;
		
		
		int root = tree.getRoot();
		DiscreteDistribution root_prior=null;
		
		final int[] postorder = TreeTraversal.postOrder(tree); // legacy logic: input file lists the rates in postorder
		int node_idx = 0;

		String line = null;
        do
        {
            line=input.readLine();
            if (line != null)
            {
                if (line.startsWith(MODEL_END))
                    break;
                line = line.trim();
                if (line.length()==0 || line.startsWith("#"))
                    continue;
                
                if (line.startsWith(RATE_VARIATION_PREFIX))
                {
                    String variation_data = line.substring(RATE_VARIATION_PREFIX.length()+1);
                    String[] fields = variation_data.split("\\s+");
                    if (fields.length<3)
                        throw new FileFormatException("Rate variation line has bad syntax: "+line);
                            
                    int num_categories = Integer.parseInt(fields[1]);
                    double alpha = Double.parseDouble(fields[2]);
                    double zero = 0.0;
                    if (fields.length>3)
                    {
                        zero = Double.parseDouble(fields[3]);
                    }
                    if ("transfer".equals(fields[0]) || "gain".equals(fields[0]))
                    {
                    	gain_categories = num_categories;
                    	gain_alpha = alpha;
                    	gain_forbidden = zero;
                    } else if ("loss".equals(fields[0]))
                    {
                    	loss_categories = num_categories;
                    	loss_alpha = alpha;
                    	loss_forbidden = zero;
                    } else if ("duplication".equals(fields[0]))
                    {
                    	duplication_categories = num_categories;
                    	duplication_alpha = alpha;
                    	duplication_forbidden = zero;
                    } else if ("length".equals(fields[0]))
                    {
                    	length_categories = num_categories;
                    	length_alpha = alpha;
                    } else 
                    {
                        throw new FileFormatException("Variation type '"+fields[0]+"' is not recognized in the line '"+line+"'");
                    }
                } else if (line.startsWith(ROOT_PRIOR_PREFIX))
                {
                    String prior_data = line.substring(ROOT_PRIOR_PREFIX.length()+1);
                    String[] fields = prior_data.split("\\s+");
                    double[] params = new double[fields.length-1];
                    for (int i=0; i<params.length; i++)
                        params[i] = Double.parseDouble(fields[i+1]);
                    if (fields[0].equals(Poisson.class.getSimpleName()))
                        root_prior = new Poisson(params[0]);
                    else if (fields[0].equals(NegativeBinomial.class.getSimpleName()))
                        root_prior = new NegativeBinomial(params[0], params[1]);
                    else if (fields[0].equals(PointDistribution.class.getSimpleName()))
                        root_prior = new PointDistribution(params[0]);
                    else if (fields[0].equals(ShiftedGeometric.class.getSimpleName()))
                        root_prior = new ShiftedGeometric(params[0], params[1]);
                    else
                        throw new FileFormatException("Root prior distribution '"+fields[0]+"' is unknown in line '"+line+"'");
                } else
                {
                	// regular node with edge length and rates
                	int node = postorder[node_idx];
                	
        	        String[] fields=line.split("\\s+");
        	        // legacy order: length, dup, loss, gain
        	        double len = Double.parseDouble(fields[0]);
        	        double drate =  Double.parseDouble(fields[1]);
        	        double lrate = Double.parseDouble(fields[2]);
        	        double grate = Double.parseDouble(fields[3]);
        	        
        	        base_rates.setEdgeLength(node, len);
        	        base_rates.setDuplicationRate(node, drate);
        	        base_rates.setLossRate(node, lrate);
        	        if (drate==0.0)
        	        {
        	        	if (lrate==0.0)
        	    	        base_rates.setGainRate(node, grate);
        	        	else
        	        		base_rates.setGainRate(node, grate/lrate); // legacy logic: store κ*μ
        	        } else
        	        {
        	        	base_rates.setGainRate(node, grate/drate); // legacy logic: store κ*λ
        	        }
        	        
        	        node_idx++;
                }
            }
        } while (line != null);
        
        last_rate_line = line;
        
        if (node_idx<=tree.getNumNodes()-1 && root_prior==null)
        	return null;
        
        // set the root prior
        if (root_prior == null)
        { // root rates were listed in postorder 
			root_prior = base_rates.getLossParameter(root)==1.0
				?base_rates.getGainDistribution(root)
				:base_rates.getDuplicationDistribution(root);
        } else
        {
	        int root_idx = tree.getRoot();
	        double[] params = root_prior.getParameters();
	        if (root_prior instanceof Poisson)
	        {
	        	base_rates.setParameters(root_idx, params[0], 1.0, 0.0); // loss prob=1 at root...
	        } else if (root_prior instanceof NegativeBinomial)
	        {
	        	base_rates.setParameters(root_idx, params[NegativeBinomial.GAIN_IDX], 1.0, params[NegativeBinomial.DUPLICATION_IDX]);
	        } else if (root_prior instanceof ShiftedGeometric)
	        {
	        	base_rates.setParameters(root_idx, 0.0, params[ShiftedGeometric.LOSS_IDX], params[ShiftedGeometric.DUPLICATION_IDX]);
	        } else
	        {
	        	assert (root_prior instanceof PointDistribution);
	        	base_rates.setParameters(root_idx, 0.0, params[0], 0.0);
	        }
        }
        
        // set the rate categories 
        readRates.setClasses(gain_categories, loss_categories, duplication_categories, length_categories);
        readRates.setGainForbidden(gain_forbidden);
        readRates.setLossForbidden(loss_forbidden);
        readRates.setDuplicationForbidden(duplication_forbidden);
        readRates.setGainAlpha(gain_alpha);
        readRates.setLossAlpha(loss_alpha);
        readRates.setDuplicationAlpha(duplication_alpha);
        readRates.setLengthAlpha(length_alpha);
		
		return readRates;
	}
	
	public static FreeMixedModel readFreeRates(BufferedReader input, IndexedTree tree, GammaInvariant class_model) throws FileFormatException, IOException	
	{
		List<TreeWithRates> class_rates_list = new ArrayList<>();
		List<Double> class_probs_list = new ArrayList<>();
		String line = null;
		if (class_model == null)
			class_model = readRates(input, tree);		
		
		while(class_model != null)
		{
			line = last_rate_line;
			
//			System.out.println("#**RVP.rFR line "+line+"\t// "+class_rates_list.size());
			
			if (line != null && line.startsWith(MODEL_END))
			{
				String class_data = line.substring(MODEL_END.length()+1);
                String[] fields = class_data.split("\\s+");			
                int class_idx = Integer.parseInt(fields[0]);
                double pc = Double.parseDouble(fields[1]);
                while (class_idx>=class_rates_list.size())
                {
                	class_rates_list.add(null);
                	class_probs_list.add(null);
                }
                class_rates_list.set(class_idx, class_model.getBaseModel());
                class_probs_list.set(class_idx, pc);
                
                class_model = readRates(input, tree);
			} else if (line == null)
			{
				for (int c=0; c<class_model.getNumClasses(); c++)
				{
					double pc = class_model.getClassProbability(c);
					if (pc != 0.0)
					{
						class_rates_list.add(class_model.getClassModel(c));
						class_probs_list.add(pc);
					}
				}
				class_model = null; // no more 
			} else
			{
				 throw new FileFormatException("Unrecognized variation line: "+line);				
			}
		} 
		
		TreeWithRates[] class_rates = class_rates_list.toArray(new TreeWithRates[0]);
		double[] class_weights = new double[class_rates.length];
		for (int c=0; c<class_weights.length; c++)
		{
			class_weights[c] = class_probs_list.get(c);
		}
		FreeMixedModel free_model = new FreeMixedModel(class_rates, class_weights);
		return free_model;
	}
	
	
	
	public static String printRates(TreeWithRates base_model)
	{
		IndexedTree tree = base_model.getTree();
		int[] all_nodes = TreeTraversal.postOrder(tree); // must be listed in this order bc that's how it is expected on input
		
		StringBuilder sb = new StringBuilder();
		
		for (int node:all_nodes)
		{
			if (tree.isRoot(node))
				sb.append("#ROOTRATES ");
			if (!tree.isRoot(node) || true) // root is last entry, but commented out
			{
				double len = base_model.getEdgeLength(node);
				double grate = base_model.getGainRate(node);
				double lrate = base_model.getLossRate(node);
				double drate = base_model.getDuplicationRate(node);
				
				// keeping legacy scaling
				if (drate==0.0)
				{
					if (lrate!=0.0) grate *= lrate;  
				} else 
				{
					grate *= drate;
				}
				double p = base_model.getLossParameter(node);
				double q = base_model.getDuplicationParameter(node);
				double r = base_model.getGainParameter(node);
				sb.append(len)
				.append("\t").append(drate)
				.append("\t").append(lrate)
				.append("\t").append(grate);
				sb.append("\t// params")
				.append("\t").append(p)
				.append("\t").append(q)
				.append("\t").append(r)
				.append("\t// ").append(tree.toString(node));
				sb.append("\n");
			}
		}
		return sb.toString();
	}
	
	public static String printRootPrior(TreeWithRates base_model)
	{
		IndexedTree tree = base_model.getTree();
		StringBuilder sb = new StringBuilder();
        sb.append(ROOT_PRIOR_PREFIX);
        sb.append('\t');
        
        int root = tree.getRoot();
        DiscreteDistribution root_prior;
        if (base_model.getGainRate(root)==0.0)
        	root_prior = base_model.getDuplicationDistribution(root);
        else 
        	root_prior = base_model.getGainDistribution(root);
        sb.append(root_prior.getClass().getSimpleName());
        double[] params = root_prior.getParameters();
        for (int i=0; i<params.length; i++)
        {
            sb.append('\t');
            sb.append(params[i]);
        }
        sb.append("\n");
//      sb.append(MODEL_END);
//      sb.append("\n");
        return sb.toString();		
	}
	
	public static String printRates(GammaInvariant rates)
	{
		TreeWithRates base_model = rates.getBaseModel();
		IndexedTree tree = base_model.getTree();
//		int[] all_nodes = TreeTraversal.postOrder(tree); // must be listed in this order bc that's how it is expected on input
		
		StringBuilder sb = new StringBuilder(printRates(base_model));
		
//		for (int node:all_nodes)
//		{
//			if (tree.isRoot(node))
//				sb.append("#ROOTRATES ");
//			if (!tree.isRoot(node) || true) // root is last entry
//			{
//				double len = base_model.getEdgeLength(node);
//				double grate = base_model.getGainRate(node);
//				double lrate = base_model.getLossRate(node);
//				double drate = base_model.getDuplicationRate(node);
//				
//				// keeping legacy scaling
//				if (drate==0.0)
//				{
//					if (lrate!=0.0) grate *= lrate;  
//				} else 
//				{
//					grate *= drate;
//				}
//				sb.append(len+"\t"+drate+"\t"+lrate+"\t"+grate+"\t// "+tree.toString(node)+"\n");
//			}
//		}
        sb.append(RATE_VARIATION_PREFIX);
        sb.append("\tduplication\t");
        sb.append(rates.getNumDuplicationGammaCategories());
        sb.append('\t');
        sb.append(rates.getDuplicationAlpha());
        sb.append('\t');
        sb.append(rates.getDuplicationForbidden());
        sb.append("\n");
        sb.append(RATE_VARIATION_PREFIX);
        sb.append("\tloss\t");
        sb.append(rates.getNumLossGammaCategories());
        sb.append('\t');
        sb.append(rates.getLossAlpha());
        sb.append('\t');
        sb.append(rates.getLossForbidden());
        sb.append("\n");
        sb.append(RATE_VARIATION_PREFIX);
        sb.append("\ttransfer\t"); // = gain --- backward compatibility for output
        sb.append(rates.getNumGainGammaCategories());
        sb.append('\t');
        sb.append(rates.getGainAlpha());
        sb.append('\t');
        sb.append(rates.getGainForbidden());
        sb.append("\n");
        sb.append(RATE_VARIATION_PREFIX);
        sb.append("\tlength\t");
        sb.append(rates.getNumLengthGammaCategories());
        sb.append('\t');
        sb.append(rates.getLengthAlpha());
        sb.append("\n");
        sb.append(ROOT_PRIOR_PREFIX);
        sb.append('\t');
        
        int root = tree.getRoot();
        DiscreteDistribution root_prior;
        if (base_model.getGainRate(root)==0.0)
        	root_prior = base_model.getDuplicationDistribution(root);
        else 
        	root_prior = base_model.getGainDistribution(root);
        sb.append(root_prior.getClass().getSimpleName());
        double[] params = root_prior.getParameters();
        for (int i=0; i<params.length; i++)
        {
            sb.append('\t');
            sb.append(params[i]);
        }
        sb.append("\n");
//      sb.append(MODEL_END);
//      sb.append("\n");
        return sb.toString();
	}

	public static String printRates(FreeMixedModel model)
	{
		StringBuilder sb = new StringBuilder();
		for (int c=0; c<model.getNumClasses(); c++)
		{
			TreeWithRates class_rates = model.getClassModel(c);
			sb.append(printRates(class_rates));
			sb.append(printRootPrior(class_rates));

			double pc = model.getClassProbability(c);
			sb.append(MODEL_END);
			sb.append("\t").append(c);
			sb.append("\t").append(pc);
			sb.append("\n");
		}
		return sb.toString();
	}
	
    /**
     * Our own exception type.
     */
    public static class FileFormatException extends IOException
    {
        private FileFormatException(String msg)
        {
            super(msg);
        }
    }
    
    private void mainmain(String[] args) throws Exception
    {
        
        int arg_idx = 0;
        String tree_file = args[arg_idx++];
        String rate_file = args[arg_idx++];
        
        IndexedTree tree = NewickParser.readTree(new count.io.GeneralizedFileReader(tree_file));
        GammaInvariant zeb = readRates(new count.io.GeneralizedFileReader(rate_file), tree);
        
        System.out.println("Read model:");
        System.out.println(printRates(zeb));
        
    }
    
    public static void main(String[] args) throws Exception
    {
        RateVariationParser O = new RateVariationParser();
        O.mainmain(args);
        
    }
    
	

}
