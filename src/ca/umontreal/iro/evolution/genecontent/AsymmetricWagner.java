/*
 * AsymmetricWagner.java
 *
 * Created on June 16, 2008, 10:01 AM
 */

package ca.umontreal.iro.evolution.genecontent;

import java.io.BufferedReader;

import java.util.Hashtable;

import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.Verbose;


/**
 *
 * @author  csuros
 */
public class AsymmetricWagner extends BasicExecutable
{
    public static final String VERSION_INFO = "11.0502";
    
    public static double GAIN_PENALTY = 5.0;
    
    /**
     * Minimum number of lineages in which families are present in the input table.
     * Must be 0,1, or 2.
     */
    public static int MIN_PRESENT_LINEAGES = 0;
    
    /**
     * Maximum total number of paralogs in the input tables
     */
    public static int MAX_PARALOGS = Integer.MAX_VALUE;
    
    public static String FUNCTIONAL_CATEGORIES = null;
    
    /** Creates a new instance of AsymmetricWagner */
    private AsymmetricWagner() 
    {
    }
    
    @Override
    protected String executableInfo()
    {
        return getClass().getName()+" "+VERSION_INFO;
    }

    private void go(String[] args) throws Exception
    {
        if (args.length < 2 || args.length>3)
        {
            System.err.println("Call as java "+this.getClass().getCanonicalName()+" [-gain x] [-max_paralogs n] [-min_lineages] phylogeny table");
            System.exit(2008);
        }
        reportLaunch(args);
        
        String tree_file = args[0];
        String table_file = args[1];
        
        reportOtherArguments("Gain penalty "+GAIN_PENALTY);
        

        TreeWithRates main_tree = new TreeWithRates(ca.umontreal.iro.evolution.Parser.readNewick(new java.io.FileReader(tree_file)));
        OccurrenceTable table = new OccurrenceTable(main_tree.getLeaves());
        table.readTable(new java.io.FileReader(table_file));
        OccurrenceTable filtered_table = table.filterByMaximumSize(MAX_PARALOGS,MIN_PRESENT_LINEAGES);
        
        int num_profiles = filtered_table.getNumFamilies();
        //num_families = 1;
        
        NodeWithRates[] nodes = main_tree.getDFT();
        int[] num_family_gains = new int[nodes.length-1];
        int[] num_family_losses = new int[nodes.length-1];
        int[] num_gene_gains = new int[nodes.length-1];
        int[] num_gene_losses = new int[nodes.length-1];
        int[] num_gene_duplications = new int[nodes.length-1];
        int[] num_family_duplications = new int[nodes.length-1];
        int[] num_family_contractions = new int[nodes.length-1];
        int[] num_genes = new int[nodes.length];
        int[] num_families = new int[nodes.length];
        int[] num_single_families = new int[nodes.length];
        
        Hashtable<String,Integer> category_index = new Hashtable<String,Integer>();
        Hashtable<String,String> family_category = new Hashtable<String,String>();
        int[][] category_gain=null;
        int[][] category_expansion=null;
        int[][] category_loss=null;
        int[][] category_reduction=null;
        
        if (FUNCTIONAL_CATEGORIES != null)
        {
            BufferedReader BR = new BufferedReader(new java.io.FileReader(FUNCTIONAL_CATEGORIES));
            
            String line=null;
            do 
            {
                line = BR.readLine();
                if (line == null)
                    break;
                
                String[] fields = line.split("\\,");
                String name = fields[0];
                String category = fields[1];
                String description = fields[2];
                
                if (!category_index.containsKey(category))
                {
                    int idx = category_index.size();
                    category_index.put(category,new Integer(idx));
                }
                family_category.put(name,category); 
            } while (line!=null);
            BR.close();
            
            int num_categories = category_index.size();
            
            category_gain = new int[nodes.length-1][num_categories];
            category_expansion = new int[nodes.length-1][num_categories];
            category_loss = new int[nodes.length-1][num_categories];
            category_reduction = new int[nodes.length-1][num_categories];
        }
        
        
        
        System.out.print("# FAMILY\tname");
        for (int node_idx=0;node_idx<nodes.length; node_idx++)
            System.out.print("\t"+nodes[node_idx].getTaxonName());
        System.out.println("\tGains\tLosses\tExpansions\tReductions");
        for (int profile_idx=0; profile_idx<num_profiles; profile_idx++)
        {
            PhyleticProfile PP = filtered_table.getProfile(profile_idx);
            int[] counts = PP.computeWagnerParsimony(main_tree,GAIN_PENALTY);
            int category_idx = -1;
            if (FUNCTIONAL_CATEGORIES!=null)
            {
                String cat = family_category.get(filtered_table.getFamilyName(profile_idx));
                category_idx = category_index.get(cat).intValue();
            }
            int family_gains = 0;
            int family_losses = 0;
            int family_expansions = 0;
            int family_reductions = 0;
            for (int node_idx=0; node_idx<counts.length; node_idx++)
            {
                num_genes[node_idx]+=counts[node_idx];
                if (counts[node_idx]>0)
                    num_families[node_idx]++;
                if (counts[node_idx]==1)
                    num_single_families[node_idx]++;
                if (node_idx != nodes.length-1) // root
                {
                    int parent_idx = main_tree.getParentIndex(node_idx);
                    int parent_count = counts[parent_idx];
                    if (parent_count>counts[node_idx])
                    {
                        if (counts[node_idx]==0)
                        {
                            family_losses++;
                            num_family_losses[node_idx]++;
                            if (category_idx!=-1)
                                category_loss[node_idx][category_idx]++;
                        } else 
                        {
                            family_reductions++;
                            num_family_contractions[node_idx]++;
                            if (category_idx != -1)
                                category_reduction[node_idx][category_idx]++;
                        }
                        num_gene_losses[node_idx]+=(parent_count-counts[node_idx]);
                        
                    } else if (parent_count<counts[node_idx])
                    {
                        if (parent_count==0)
                        {
                            family_gains++;
                            num_family_gains[node_idx]++;
                            if (category_idx != -1)
                                category_gain[node_idx][category_idx]++;
                        }
                        else
                        {
                            family_expansions++;
                            num_family_duplications[node_idx]++;
                            if (category_idx != -1)
                                category_expansion[node_idx][category_idx]++;
                            num_gene_duplications[node_idx]+=(counts[node_idx]-parent_count);
                        }
                        num_gene_gains[node_idx]+=(counts[node_idx]-parent_count);
                    }
                }
            }
            
            System.out.print("# FAMILY\t"+filtered_table.getFamilyName(profile_idx));
            
            for (int node_idx=0; node_idx<counts.length; node_idx++)
            {
                System.out.print("\t"+counts[node_idx]);
            }
            System.out.println("\t"+family_gains+"\t"+family_losses+"\t"+family_expansions+"\t"+family_reductions);
        }
        
        System.out.println("# PRESENT\tnode\tgenes\tfamilies\tsingle-member families");
        for (int node_idx=0; node_idx<nodes.length; node_idx++)
        {
            System.out.println("# PRESENT\t"+main_tree.getNode(node_idx).getTaxonName()+"\t"+num_genes[node_idx]+"\t"+num_families[node_idx]+"\t"+num_single_families[node_idx]);
        }

        int tot_family_gains=0;
        int tot_gene_gains=0;
        int tot_family_losses=0;
        int tot_gene_losses=0;

        System.out.println("# CHANGE\tnode\tgene_gain\tfamily_gain\tgene_duplications\tfamily_expansions\tgene_loss\tfamily_loss\tfamily_contractions");
        for (int node_idx=0; node_idx<nodes.length-1; node_idx++)
        {
            System.out.println("# CHANGE\t"+main_tree.getNode(node_idx).getTaxonName()
                +"\t"+num_gene_gains[node_idx]+"\t"+num_family_gains[node_idx]
                +"\t"+num_gene_duplications[node_idx]+"\t"+num_family_duplications[node_idx]
                +"\t"+num_gene_losses[node_idx]+"\t"+num_family_losses[node_idx]
                +"\t"+num_family_contractions[node_idx]);
            tot_family_gains += num_family_gains[node_idx];
            tot_gene_gains += num_gene_gains[node_idx];
            tot_family_losses += num_family_losses[node_idx];
            tot_gene_losses += num_gene_losses[node_idx];
        }
        
        double pg = tot_gene_gains/(num_profiles*(nodes.length-1.0));
        double pl = tot_gene_losses/(num_profiles*(nodes.length-1.0));
        double p0 = 1.0-pg-pl;
        System.out.println("# CHANGE\ttotal"
                +"\t"+tot_gene_gains+"\t"+tot_family_gains
                +"\t\t"
                +"\t"+tot_gene_losses+"\t"+tot_family_losses
                +"\t"
                +"\t// N="+num_profiles+"\tlog(gain)/log(loss)="+Math.log(pg/p0)/Math.log(pl/p0));
        
        if (FUNCTIONAL_CATEGORIES != null)
        {
            String[] category_list = category_index.keySet().toArray(new String[0]);
            java.util.Arrays.sort(category_list);
            System.out.print("# CATEGORY\tedge");
            for (int i=0; i<category_list.length; i++)
                System.out.print("\tgain-"+category_list[i]);
            for (int i=0; i<category_list.length; i++)
                System.out.print("\texpansion-"+category_list[i]);
            for (int i=0; i<category_list.length; i++)
                System.out.print("\tloss-"+category_list[i]);
            for (int i=0; i<category_list.length; i++)
                System.out.print("\treduction-"+category_list[i]);
            System.out.println();
            for (int node_idx=0; node_idx<nodes.length-1; node_idx++)
            {
                System.out.print("# CATEGORY\t"+nodes[node_idx].getTaxonName());
                for (int i=0; i<category_list.length; i++)
                {
                    int cat = category_index.get(category_list[i]).intValue();
                    System.out.print("\t"+category_gain[node_idx][cat]);
                }
                for (int i=0; i<category_list.length; i++)
                {
                    int cat = category_index.get(category_list[i]).intValue();
                    System.out.print("\t"+category_expansion[node_idx][cat]);
                }
                for (int i=0; i<category_list.length; i++)
                {
                    int cat = category_index.get(category_list[i]).intValue();
                    System.out.print("\t"+category_loss[node_idx][cat]);
                }
                for (int i=0; i<category_list.length; i++)
                {
                    int cat = category_index.get(category_list[i]).intValue();
                    System.out.print("\t"+category_reduction[node_idx][cat]);
                }
                System.out.println();
            } // for each node
        } // if categories
        
        
    }
    
    public static void main(String[] args) 
    {
        Verbose.setVerbose(false);
                
        AsymmetricWagner O = new AsymmetricWagner();
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
                } else if (arg_switch.equals("gain"))
                {
                    GAIN_PENALTY = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("max_paralogs"))
                {
                    MAX_PARALOGS = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("min_lineages"))
                {
                    MIN_PRESENT_LINEAGES = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("func"))
                {
                    FUNCTIONAL_CATEGORIES = new String(arg_value);
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
