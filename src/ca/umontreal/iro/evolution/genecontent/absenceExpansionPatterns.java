/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ca.umontreal.iro.evolution.genecontent;

/**
 *
 * @author csuros
 */

import java.util.Arrays;

import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.Verbose;

import ca.umontreal.iro.matek.SpearmanRankOrderCorrelation;

public class absenceExpansionPatterns extends BasicExecutable
{

    private void go(String[] args) throws Exception
    {
        if (args.length != 2)
        {
            System.err.println("Call as java "+getClass().getCanonicalName()+" tree table");
            System.exit(2008);
        }
        reportLaunch(args);
        reportOtherArguments("Max. paralogs "+MAX_PARALOGS+", min. lineages "+MIN_PRESENT_LINEAGES);
        
        String tree_file = args[0];
        String table_file = args[1];

        TreeWithRates input_tree = new TreeWithRates(ca.umontreal.iro.evolution.Parser.readNewick(new java.io.FileReader(tree_file)));
        OccurrenceTable input_table = new OccurrenceTable(input_tree.getLeaves());
        input_table.readTable(new java.io.FileReader(table_file));
        OccurrenceTable filtered_table = input_table.filterByMaximumSize(MAX_PARALOGS,MIN_PRESENT_LINEAGES);
        
        int num_families = filtered_table.getNumFamilies();
        System.out.println("Family\tPattern\tZeros\tDifferent positive values");
        double[] nz = new double[num_families];
        double[] nd = new double[num_families];
        for (int family_idx=0; family_idx<num_families; family_idx++)
        {
            PhyleticProfile PP = filtered_table.getProfile(family_idx);
            int[] pattern = PP.getProfile().clone();
            Arrays.sort(pattern);
            int num_zeros = 0;
            int num_different = 0;
            for (int i=0; i<pattern.length; i++)
            {
                if (pattern[i]==0)
                    num_zeros++;
                else
                    if (i==0 || pattern[i] != pattern[i-1])
                        num_different ++;
            }
            System.out.println(filtered_table.getFamilyName(family_idx)+"\t"+PP.getPatternString()+"\t"+num_zeros+"\t"+num_different);
            nz[family_idx]=num_zeros;
            nd[family_idx]=num_different;
        } // for families
        
        SpearmanRankOrderCorrelation Corr = new SpearmanRankOrderCorrelation(nz,nd);
        System.out.println("#SPEARMAN "+Corr);
        
    }
    
    private static int MAX_PARALOGS = 150;
    private static int MIN_PRESENT_LINEAGES = 1;
    
    /**
     * Outputs the posterior probabilities in a current_table format to standard output
     * 
     * @param args Command-line arguments
     */
    public static void main(String[] args)
    {
        Verbose.setVerbose(false);
                
        absenceExpansionPatterns O = new absenceExpansionPatterns(); //InclusionExclusionComputation());
        
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
