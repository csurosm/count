package ca.umontreal.iro.evolution.genecontent;


import java.io.Reader;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Vector;

import ca.umontreal.iro.evolution.TreeNode;

/**
 * Stores a table of family sizes across a number of organisms. 
 *
 * A table is always linked to a list of 
 * terminal taxon names. This list may be set 
 * at instantiation, or else set by readTable().
 * 
 * @author csuros
 */

public class SimpleOccurrenceTable 
{
    /**
     * Instantiates a table with the names of these taxa.
     * @param terminal_taxa array of terminal taxon names
     */
    public SimpleOccurrenceTable(String[] terminal_taxa) 
    {
        init(terminal_taxa);
    }
    
    /**
     * Instantiates a table with a given array of taxon names.
     * @param terminal_taxa array of terminal nodes; their names are used in this class  
     */
    public SimpleOccurrenceTable(TreeNode[] terminal_taxa)
    {
        String[] names = new String[terminal_taxa.length];
        for (int j=0; j<names.length; j++)
            names[j] = terminal_taxa[j].getName();
        init(names);
    }
    
    /**
     * Initializes the table with a null set of terminal taxa.
     * Terminal taxa will be determined when the table is set.
     */
    public SimpleOccurrenceTable(){ this((String[])null);}
    
    
    private String[] terminal_taxa;
    
    private int[][] table;
    
    private String[] family_names;
    

    /**
     * Returns the number of families
     * @return number of families 
     */
    public int getNumFamilies()
    {
        return family_names.length;
    }
    
    /**
     * Returns the array of terminal taxon names used at instantiation
     * @return array of terminal taxon names 
     */
    public String[] getTerminalTaxonNames()
    {
        return terminal_taxa;
    }

    /**
     * Size distribution of the given family across the terminal taxa
     * @param family_idx index of the family
     * @return array of family sizes for the terminal taxa
     */
    public int[] getSizes(int family_idx)
    {
        return table[family_idx];
    }
    
    /**
     * Name of the given family
     * @param family_idx index of the family
     * @return name of the family with the given index
     */
    public String getFamilyName(int family_idx)
    {
        return family_names[family_idx];
    }
    
    
    /**
     * Called by instantiation methods 
     *
     * @param taxa may be null
     */
    private void init(String[] taxa)
    {
        //System.out.println("#**IT.i "+taxa.length);
        this.terminal_taxa = taxa;
        table = null;
        family_names = new String[0];
    }
    
    /**
     * Sets the table from externally constructed data. 
     * 
     * @param table family sizes: table[i][j] is the size of the i-th family in the j-th species (this latter indexed as in terminal_taxa)
     * @param family_names names for the families in the order of the indexes into the int[] array; if null, then they are set to are set to "F1,F2,..."
     */
    public void setTable(int[][] table, String[] family_names)
    {
        this.table = table;
        int nfam = table.length;
        this.family_names = family_names;
        if (family_names == null)
        {
            this.family_names = new String[nfam];
            for (int i=0; i<nfam; i++)
                this.family_names[i] = "F"+Integer.toString(i+1);
        }
        
    }
    
    /**
     * Sets the table from externally constructed data. Family names are set to are set to "F1,F2,..."
     * 
     * @param table family sizes: table[i][j] is the size of the i-th family in the j-th species (this latter indexed as in terminal_taxa)
     */ 
    public void setTable(int[][] table)
    {
        setTable(table, null);
    }    
    
    /**
     * Adds a new profile to the table. The new profile is inserted at index 0.
     * 
     * @param new_profile family size across terminal taxa 
     * @param new_profile_name name of the new profile
     */
    public void addProfile(int[] new_profile, String new_profile_name)
    {
        int[][] new_table = new int[table.length+1][];
        System.arraycopy(table, 0, new_table, 1, table.length);
        new_table[0] = new_profile;
        String[] new_family_names = new String[family_names.length+1];
        System.arraycopy(family_names, 0, new_family_names, 1, family_names.length);
        new_family_names[0] = new_profile_name;
        setTable(new_table, new_family_names);
    }
    
    /*
     * Reads in an OccurrenceTable using a reader
     * 
     * @param reader the input reader from which the lines are read
     * 
     * @throws java.io.IOException if there is something wrong with I/O 
     */
    public void readTable(Reader reader) throws IOException
    {
        // my list of terminal taxa may be in a different order than in the file ...
        Hashtable<String,Integer> taxon_indexes = new Hashtable<String,Integer>();
        for (int i=0; i<terminal_taxa.length; i++)
            taxon_indexes.put(terminal_taxa[i], new Integer(i));
        
        // here we collect the names of the families
        Vector<String> read_names = new Vector<String>();
        Vector<int[]> read_sizes = new Vector<int[]>();

        // here we figure out the proper permutation of the columns 
        int[] taxon_order = null;

        BufferedReader R=new BufferedReader(reader);
        String line=null;
        int num_lines=0;
        do 
        {
            line = R.readLine();
            if (line == null || line.trim().length()==0 || line.startsWith("#"))
                continue;
            num_lines++;
            String[] fields=line.split("\t");
            if (num_lines == 1)
            {
                // header line
                taxon_order=new int[fields.length-1];
                for (int j=1; j<fields.length; j++)
                {
                    if (taxon_indexes.containsKey(fields[j]))
                    {
                        int taxon_idx = ((Integer)taxon_indexes.get(fields[j])).intValue();
                        taxon_order[j-1]=taxon_idx;
                    } else 
                    {
                        taxon_order[j-1]=-1;
                    }
                }
            } else
            {
                // data lines
                String family = fields[0]; 
                int[] copies = new int[terminal_taxa.length];
                for (int j=1; j<fields.length; j++)
                {
                    int idx = taxon_order[j-1];
                    if (idx != -1)
                    {
                        if ("?".equals(fields[j]))
                        {
                            copies[idx]=-1;
                        } else
                        {
                            int m=Integer.parseInt(fields[j]);
                            copies[idx]=m;
                        }
                    } 
                }
                read_names.add(family);
                read_sizes.add(copies);
            }
        } while (line !=null);
        
        R.close();
        
        setTable(read_sizes.toArray(new int[0][]), read_names.toArray(new String[0]));
    }
    
    /**
     * Selects a subset of families
     * 
     * @param family_ok array of booleans: family_ok[i] specifies if tyyhe i-th family should be included
     * 
     * @return a new table that contains only the selected families
     */
    public SimpleOccurrenceTable tableForFamilies(boolean[] family_ok)
    {
        int num_profiles = getNumFamilies();
        int num_filtered_families=0;
        for (int pidx=0; pidx<num_profiles; pidx++)
            if (family_ok[pidx])
                num_filtered_families++;
        int[][] new_table = new int[num_filtered_families][];
        String[] new_names = new String[num_filtered_families];
        {
            int family_idx = 0;
            for (int pidx=0; pidx<num_profiles; pidx++)
                if (family_ok[pidx])
                {
                    new_table[family_idx]=table[pidx];
                    new_names[family_idx]=family_names[pidx];
                    family_idx++;
                }
        }
        SimpleOccurrenceTable T = new SimpleOccurrenceTable();
        T.terminal_taxa = this.terminal_taxa;
        T.setTable(new_table, new_names);
        // share all propeties
        
        return T;
    }
            
    /**
     * Convenience method for selecting a range of families
     * 
     * @param first_family_idx first family index to be included
     * @param last_family_idx last family index to be included 
     * @return a table for families with indexes first_family_idx..last_family_idx (both ends inclusive)
     */
    public SimpleOccurrenceTable tableForFamilies(int first_family_idx, int last_family_idx)
    {
        boolean[] family_ok = new boolean[getNumFamilies()];
        for (int i=first_family_idx; i<=last_family_idx; i++)
            family_ok[i]=true;
        return tableForFamilies(family_ok);
    }
     
}
