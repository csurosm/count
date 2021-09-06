package ca.umontreal.iro.evolution.genecontent;

/**
 * Class for computing the posterior probability for 
 * gains on each edge
 * 
 * @author csuros
 */

import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.Verbose;

public class PosteriorGain extends BasicExecutable
{
    private static int MAX_PARALOGS = 1000000;
    private static int MIN_PRESENT_LINEAGES = 1;

    private PosteriorGain(){}
            
    
    private void go(String[] args) throws Exception
    {
        if (args.length != 3)
        {
            System.err.println("Call as java "+getClass().getCanonicalName()+" tree table rates");
            System.exit(2008);
        }
        reportLaunch(args);
        reportOtherArguments("Max. paralogs "+MAX_PARALOGS+", min. lineages "+MIN_PRESENT_LINEAGES);
        
        String tree_file = args[0];
        String table_file = args[1];
        String rates_file = args[2];

        TreeWithRates input_tree = new TreeWithRates(ca.umontreal.iro.evolution.Parser.readNewick(new java.io.FileReader(tree_file)));
        OccurrenceTable input_table = new OccurrenceTable(input_tree.getLeaves());
        input_table.readTable(new java.io.FileReader(table_file));
        OccurrenceTable filtered_table = input_table.filterByMaximumSize(MAX_PARALOGS,MIN_PRESENT_LINEAGES);
        int num_families = filtered_table.getNumFamilies();

        RateVariation input_model = RateVariation.read(new java.io.FileReader(rates_file),input_tree);
        
        NodeWithRates[] nodes = input_tree.getDFT();
        
        double[][] likelihood = new double[nodes.length][num_families];
        double[] absent_probability = new double[nodes.length];

        // first compute normal likelihoods, stored under the root's index
        int root_idx = nodes.length-1;
        {
            StableComputation SC = new StableComputation(filtered_table, input_model);
            absent_probability[root_idx] = SC.getAbsentProfileProbability(1);
            for (int pidx=0; pidx<num_families; pidx++)
            {
                double p = SC.getLikelihood(pidx);
                likelihood[root_idx][pidx] = p;
            }
        }
        for (int node_idx=0; node_idx<nodes.length;node_idx++)
        {
            NodeWithRates N = nodes[node_idx];
            if (!N.isRoot())
            {
                double g = N.getTransferRate();
                N.setTransferRate(0.0);
                // 
                RateVariation R = input_model.sameModelForDifferentTree(input_tree);
                StableComputation SC = new StableComputation(filtered_table, R);
                absent_probability[node_idx] = SC.getAbsentProfileProbability(1);
                for (int pidx=0; pidx<num_families; pidx++)
                {
                    double p = SC.getLikelihood(pidx);
                    likelihood[node_idx][pidx] = p;
                }
                
                N.setTransferRate(g);// reset gain rate
            }
        }
        
        System.out.print("Family\tPattern\tLikelihood\tCorrected");
        for (int node_idx=0; node_idx<nodes.length;node_idx++)
        {
            NodeWithRates N = nodes[node_idx];
            if (!N.isRoot())
            {
                String name = N.newickName();
                System.out.print("\t"+name);
            }
        }
        System.out.println();
        for (int pidx=0; pidx<num_families; pidx++)
        {
            PhyleticProfile profile = filtered_table.getProfile(pidx);
            System.out.print(filtered_table.getFamilyName(pidx)+"\t"+profile.getPatternString());
            double lik = likelihood[root_idx][pidx];
            double p0 = absent_probability[root_idx];
            double p = lik/(1.-p0);
            System.out.print("\t"+lik+"\t"+p);
            for (int node_idx=0; node_idx<nodes.length;node_idx++)
            {
                NodeWithRates N = nodes[node_idx];
                if (!N.isRoot())
                {
                    double likN = likelihood[node_idx][pidx]/(1.0-absent_probability[node_idx]);
                    double p_gain = 1.0-(likN / p);
                    System.out.print("\t"+p_gain);
                }
            }
            System.out.println();
        }
        
    }
    
    public static void main(String[] args)
    {
        Verbose.setVerbose(false);
        
        PosteriorGain O = new PosteriorGain(); //InclusionExclusionComputation());
        
        int num_switches = 0;
        try {
            while (args.length>2*num_switches && args[2*num_switches].startsWith("-"))
            {
                String arg_switch = args[2*num_switches].substring(1);
                if (arg_switch.equals("h"))
                    O.go(new String[0]); // will throw an Exception
                if (args.length==2*num_switches+1)
                    throw new IllegalArgumentException("Missing argument for switch "+args[2*num_switches]);
                String arg_value = args[2*num_switches+1];
                if (arg_switch.equals("v"))
                {
                    Verbose.setVerbose(arg_value.equals("true"));
                } else if (arg_switch.equals("max_paralogs"))
                {
                    MAX_PARALOGS = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("min_lineages"))
                {
                    MIN_PRESENT_LINEAGES = Integer.parseInt(arg_value);
                } 
                else 
                    throw new IllegalArgumentException("Switch not recognized: '"+args[2*num_switches]+"'");
                    
                num_switches++;
            }
            
            String[] rest=new String[args.length-2*num_switches];
            for (int j=0; j<rest.length; j++)
                rest[j]=args[2*num_switches+j];
            O.go(rest);
        } catch (Exception E)
        {
            die(E);
        }        
        
    }

}
