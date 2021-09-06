
package ca.umontreal.iro.evolution.genecontent;

/**
 * Class for producing bootstrap samples.
 * 
 * @author csuros
 */

import java.util.Random;
import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.Verbose;


public class BootstrapTable extends BasicExecutable
{
    public BootstrapTable(OccurrenceTable table)
    {
        setTable(table);
    }
    
    private BootstrapTable(){}
    
    private void setTable(OccurrenceTable table)
    {
        this.table = table;
        this.RND = new Random();
    }
    
    private OccurrenceTable table;
    private Random RND;
    

    /**
     * Picks a random family from the tabe and returns its profile
     * 
     * @param RND random number generator
     * @return a profile picked randomly from the table
     */
    public RandomProfile randomProfile(Random RND)
    {
        int rnd_idx = RND.nextInt(table.getNumFamilies());
        return new RandomProfile(rnd_idx);
    }
    
    /**
     * Picks a random family from the tabe and returns its profile
     * 
     * @return a profile picked randomly from the table
     */
    public RandomProfile randomProfile()
    {
        return randomProfile(RND);
    }
    
    /**
     * Generates a bootstrap sample
     * 
     * @return a new table constructed from the one used at instantiation
     */
    public OccurrenceTable bootstrapSample()
    {
        int num_families = table.getNumFamilies();
        int[][] bootstrap_table = new int[num_families][];
        String[] family_names = new String[num_families];
        for (int family_idx=0; family_idx<num_families; family_idx++)
        {
            
            RandomProfile P = randomProfile();
            bootstrap_table[family_idx] = P.getProfile();
            family_names[family_idx] = "F"+(1+family_idx)+"/"+table.getFamilyName(P.getIndex());
        }
        String[] leaves = table.getTerminalTaxonNames();
        OccurrenceTable O = new OccurrenceTable(leaves);
        
        O.setTable(bootstrap_table,family_names);
        return O;
    }
    
    private void go(String[] args) throws Exception
    {
        if (args.length != 2)
        {
            System.err.println("Call as "+getClass().getCanonicalName()+" tree table");
            System.exit(1107);
        }
        
        this.reportLaunch(args);
        
        String tree_file = args[0];
        String table_file = args[1];
        
        

        TreeWithRates main_tree = new TreeWithRates(ca.umontreal.iro.evolution.Parser.readNewick(new java.io.FileReader(tree_file)));
        OccurrenceTable input_table = new OccurrenceTable(main_tree.getLeaves());
        input_table.readTable(new java.io.FileReader(table_file));
        OccurrenceTable filtered_table = input_table.filterByMaximumSize(MAX_PARALOGS,MIN_PRESENT_LINEAGES);
        Verbose.message("BT.go "+filtered_table.getNumFamilies()+" families");
        setTable(filtered_table);
        OccurrenceTable bootstrap_table = bootstrapSample();
        System.out.println(bootstrap_table.getFormattedTable());
   }
    
   public class RandomProfile extends PhyleticProfile
   {
       private RandomProfile(int profile_index)
       {
           super(table.getProfile(profile_index).getProfile());
           this.profile_index = profile_index;
       }
       private int profile_index;
       
       public int[] getProfile()
       {
           return table.getProfile(profile_index).getProfile();
       }
       
       public int getIndex()
       {
           return profile_index;
       }
   }
    
   private static int MAX_PARALOGS=1000000;
   private static int MIN_PRESENT_LINEAGES=1;
    
    public static void main(String[] args)
    {
        Verbose.setVerbose(false);
                
        BootstrapTable O = new BootstrapTable();
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
