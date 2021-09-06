/*
 * OccurrenceTable.java
 *
 * Created on April 17, 2008, 3:23 PM
 */

package ca.umontreal.iro.evolution.genecontent;


import java.util.Hashtable;
import java.util.Properties;
import java.util.Vector;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

import ca.umontreal.iro.banality.StringSplit;
import ca.umontreal.iro.banality.Verbose;

import ca.umontreal.iro.evolution.TreeNode;
import ca.umontreal.iro.evolution.Parser;

/**
 * Stores a table of family sizes across a number of organisms. 
 *
 * A table is always linked to a list of 
 * terminal taxon names. This list may be initialized 
 * at instantiation, or else set by {@link #readTable(java.io.Reader) }
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */


public class OccurrenceTable 
{
    
    /**
     * Instantiates a table with the names of these taxa.
     * @param terminal_taxa array of terminal taxon names
     */
    public OccurrenceTable(String[] terminal_taxa) 
    {
        init(terminal_taxa);
    }
    
    /**
     * Instantiates a table with a given array of taxon names.
     * @param terminal_taxa array of terminal nodes; their names are used in this class  
     */
    public OccurrenceTable(TreeNode[] terminal_taxa)
    {
        String[] names = new String[terminal_taxa.length];
        for (int j=0; j<names.length; j++)
            names[j] = terminal_taxa[j].getName();
        init(names);
    }
    
    /**
     * Instantiates a table with a given array of taxon names.
     * @param terminal_taxa array of terminal nodes; their names are used in this class  
     */
    public OccurrenceTable(NodeWithRates[] terminal_taxa)
    {
        this((TreeNode[])terminal_taxa);
    }
    

    /**
     * Initializes the table with a null set of terminal taxa.
     * Terminal taxa will be determined when the table is set.
     */
    public OccurrenceTable(){ this((String[])null);}
    
    /**
     * Returns the array of terminal taxon names used at instantiation
     * @return array of terminal taxon names 
     */
    public String[] getTerminalTaxonNames()
    {
        return terminal_taxa;
    }

    /**
     * Terminal taxon names are stored here
     */
    private String[] terminal_taxa;

    /**
     * Family profiles are stored here (row=family, column=taxon)
     */
    private int[][] table;

    /**
     * Family names are stored here
     */
    private String[] family_names;

    /**
     * Family properties are stored here (one entry per family).
     */
    private Properties[] family_properties;

    /**
     * Property names are stored here.      
     * Property 0 is Family.
     */
    private String[] property_names;

    private boolean has_missing_entries; // set by setTable()

    /**
     * Whether this table has missing entries.
     *
     * @return
     */
    public boolean hasMissingEntries()
    {
        return has_missing_entries;
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
        family_properties = new Properties[0];
        property_names = new String[0];
    }

    private void copyTerminalTaxa(OccurrenceTable template)
    {
        this.terminal_taxa = template.terminal_taxa;
    }
    
    /**
     * Returns the number of families
     * @return number of families 
     */
    public int getNumFamilies()
    {
        return family_names.length;
    }
    

    /**
     * Returns the number of families present at a terminal taxon
     * 
     * @param taxon_idx index of the terminal taxon (same order as getTerminalTaxonNames)
     * @return number of families with non-zero size
     */
    public int getNumFamiliesPresent(int taxon_idx)
    {
        int s = 0;
        for (int i=0; i<table.length; i++)
            if (table[i][taxon_idx]>0)
                s++;
        return s;
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
     * Computes an array of profiles for this occurrence table
     * 
     * @return array of phyletic profiles
     */
    public PhyleticProfile[] getProfiles()
    {
        int num_profiles = getNumFamilies();
        
        PhyleticProfile[] profiles = new PhyleticProfile[num_profiles];
        for (int profile_idx=0; profile_idx<num_profiles; profile_idx++)
        {
            profiles[profile_idx] =getProfile(profile_idx);
        }
        return profiles;
    }
    
    /**
     * Retrieves one of the phyletic profiles in the data table.
     * 
     * @param profile_idx index of family for which the profile is needed
     * @return a phyletic profile describing family sizes across the terminal taxa
     */
    public PhyleticProfile getProfile(int profile_idx)
    {
        return new PhyleticProfile(table[profile_idx]);
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
     * Each family is associated with a set of properties.
     * Initially, the single known property is the family name.
     * Further properties can be added later. 
     * 
     * @return array of known property names
     */
    public String[] getKnownProperties()
    {
        return property_names.clone();
    }
    
    /**
     * Each family is associated with a set of properties.
     * Initially, the single known property is the family name.
     * Further properties can be added later. 
     * 
     * @param property_idx 0-based index of a known property
     * 
     * @return name of the given property
     */
    public String getPropertyName(int property_idx)
    {
        return property_names[property_idx];
    }
    
    /**
     * Each family is associated with a set of properties.
     * Initially, the single known property is the family name.
     * Further properties can be added later. 
     * 
     * @return number of known properties (at least 1)
     */
    public int getKnownPropertiesCount()
    {
        return property_names.length;
    }

    /**
     * Each family is associated with a set of properties.
     * Initially, the single known property is the family name.
     * Further properties can be added later. 
     * 
     * @param family_idx 0-based index of the family
     * @param property_name name of a known property
     * @return the value of the given property, or null if the property is unknown
     */
    public String getFamilyProperty(int family_idx, String property_name)
    {
        return family_properties[family_idx].getProperty(property_name);
    }

    /**
     * Each family is associated with a set of properties.
     * Initially, the single known property is the family name.
     * Further properties can be added later. 
     * 
     * @param family_idx 0-based index of the family
     * @param property_idx 0-based index of known properties
     * @return the value of the given property
     */
    public String getFamilyProperty(int family_idx, int property_idx)
    {
        return family_properties[family_idx].getProperty(property_names[property_idx]);
    }
    
    public void setFamilyProperty(int family_idx, String property_name, String property_value)
    {
        family_properties[family_idx].setProperty(property_name, property_value);
    }
    
    public void setFamilyProperty(int family_idx, int property_idx, String property_value)
    {
        setFamilyProperty(family_idx,property_names[property_idx],property_value);
    }
    
    /**
     * Registers a new property that applies to every family. 
     * 
     * @param property_name name of the new property
     * @return index of the new property
     */
    public int registerProperty(String property_name)
    {
        if (property_name==null)
            property_name = "P"+Integer.toString(property_names.length);
        String[] new_props = new String[property_names.length+1];
        for (int j=0; j<property_names.length; j++)
        {
            if (property_name.equals(property_names[j]))
                throw new IllegalArgumentException("Cannot add more than one property with the same name `"+property_name+"'");
            new_props[j]=property_names[j];
        }
        //System.out.println("#*OT.rP '"+property_name+"'\tidx "+property_names.length);
        new_props[property_names.length]=property_name;
        property_names = new_props;
        return property_names.length-1;
    }

    /**
     * Adds a new profile to the table. The new profile is inserted
     * at index 0, and all otehr family indices are shifted up by one.
     *
     * @param new_profile the new profile
     * @param new_profile_name name of the new family
     */
    public String[] addProfile(int[] new_profile, String new_profile_name)
    {
        int[][] new_table = new int[table.length+1][];
        System.arraycopy(table, 0, new_table, 1, table.length);
        new_table[0] = new_profile;
        String[] new_family_names = new String[family_names.length+1];
        System.arraycopy(family_names, 0, new_family_names, 1, family_names.length);
        new_family_names[0] = new_profile_name;
        setTable(new_table, new_family_names);
        return new_family_names;
    }
    
    
    /**
     * Sets the table from externally constructed data. 
     * 
     * @param table family sizes: table[i][j] is the size of the i-th family in the j-th species (this latter indexed as in terminal_taxa)
     * @param family_names names for the families in the order of the indexes into the int[] array; if null, then they are set to are set to "F1,F2,..."
     */
    public void setTable(int[][] table, String[] family_names)
    {
        setTableOnly(table);

        int nfam = table.length;
        this.family_names = family_names;
        if (family_names == null)
        {
            this.family_names = new String[nfam];
            for (int i=0; i<nfam; i++)
                this.family_names[i] = "F"+Integer.toString(i+1);
        }
        
        family_properties = new Properties[nfam];
        property_names = new String[1];
        property_names[0] = "Family";
        for (int family_idx = 0; family_idx<nfam; family_idx++)
        {
            Properties prop = new Properties();
            //if (family_names != null)
                prop.setProperty(property_names[0], this.family_names[family_idx]);
            family_properties[family_idx] = prop;
        }
        checkMissingEntries();
    }
    
    protected void setTableOnly(int[][] table)
    {
        this.table = table;
    }

    protected void setNamesOnly(String[] family_names)
    {
        this.family_names = family_names;
    }

    protected void checkMissingEntries()
    {
        this.has_missing_entries = false;
        for (int i=0; i<table.length && !has_missing_entries; i++)
            for (int j=0; j<table[i].length && !has_missing_entries; j++)
                has_missing_entries = (table[i][j]<0) ;
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

    /*
     * @param reader the input reader from which the lines are read
     * 
     * @throws java.io.IOException if there is something wrong with I/O 
     */
    public void readTable(Reader reader) throws IOException
    {
        readTable(reader, false);
    }
    
    /*
     * Reads in an OccurrenceTable using a reader
     * 
     * @param reader the input reader from which the lines are read
     * @param includes_properties whether columns that do not correspond to terminal taxa should be considered as property columns
     * 
     * @throws java.io.IOException if there is something wrong with I/O 
     */
    public void readTable(Reader reader, boolean includes_properties) throws IOException
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
        Vector<String> property_names_in_input = new Vector<String>();
        Vector<Vector<String>> family_properties = new Vector<Vector<String>>();

        BufferedReader R=new BufferedReader(reader);
        String line=null;
        int num_lines=0;
        do 
        {
            line = R.readLine();
            if (line == null || line.trim().length()==0 || line.startsWith("#"))
                continue;
            num_lines++;
            String[] fields=StringSplit.splitAtTabs(line);
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
                        if (includes_properties)
                            property_names_in_input.add(fields[j]);
                    }
                }
            } else
            {
                Vector<String> properties_in_line = new Vector<String>();
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
                    } else if (includes_properties)
                    {
                        properties_in_line.add(fields[j]);
                    }
                }
                read_names.add(family);
                read_sizes.add(copies);
                if (includes_properties)
                {
                    family_properties.add(properties_in_line);
                }
            }
        } while (line !=null);
        
        R.close();
        
        setTable(read_sizes.toArray(new int[0][]), read_names.toArray(new String[0]));
        if (includes_properties)
        {
            for (int j=0; j<property_names_in_input.size(); j++)
                registerProperty(property_names_in_input.get(j));
            int num_families = getNumFamilies();
            int num_properties = property_names_in_input.size();
            for (int family_idx = 0; family_idx<num_families; family_idx++)
            {
                Vector<String> V = family_properties.get(family_idx);
                for (int j=0; j<num_properties; j++)
                    setFamilyProperty(family_idx, j+1, V.get(j)); // property 0 is family name that is already there
            }
        }
    }
    
    /**
     * Returns a String that can be written directly into a file, in the same 
     * format as expected by readTable()
     * 
     * @return a multiline String description of the table
     */
    public String getFormattedTable()
    {
        return getFormattedTable(false);
    }
    
    public String getFormattedTable(boolean include_properties)
    {
        StringBuffer sb = new StringBuffer();
        sb.append("Family");
        if (include_properties)
            for (int prop_idx=1; prop_idx<getKnownPropertiesCount(); prop_idx++)
            {
                sb.append('\t');
                sb.append(getPropertyName(prop_idx));
            }
        for (int j=0; j<terminal_taxa.length; j++)
        {
            sb.append("\t");
            sb.append(terminal_taxa[j]);
        }
        sb.append("\n");
        for (int i=0; i<table.length; i++)
        {
            sb.append(getFamilyName(i));
            if (include_properties)
            {
                int num_properties = getKnownPropertiesCount();
                for (int prop_idx=1; prop_idx<num_properties; prop_idx++)
                {
                    sb.append('\t');
                    sb.append(getFamilyProperty(i,prop_idx));
                }
            }
            for (int j=0; j<terminal_taxa.length; j++)
            {
                sb.append("\t");
                sb.append(Integer.toString(table[i][j]));
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    /**
     * Creates an instance of this table, in which
     * profiles are binary. Names and properties are shard
     * with the new instance.
     *
     * @return a binary version of the same table
     */
    public TransformedTable binaryTable()
    {
        TransformedTable ot = new TransformedTable(this);
        int num_families  = getNumFamilies();
        int[][] tbl = new int[num_families][];
        for (int family_idx=0; family_idx<num_families; family_idx++)
        {
            PhyleticProfile pp = getProfile(family_idx);
            tbl[family_idx] = pp.binaryProfile().getProfile();
        }
        ot.setTable(tbl);
//        ot.table = tbl;
//        ot.family_names = family_names;
//        ot.family_properties = family_properties;
//        ot.property_names = property_names;
//        ot.terminal_taxa = terminal_taxa;

        return ot;
    }

    /**
     * Produces a table that includes only those entries where the sum of family sizes is not too large,
     * and that have at least 1 occurrence in somewhere 
     * 
     * @param maximum_family_size size threshold (inclusive)
     * @return a new table composed of families not surpassing the size threshold
     */
    public DerivedTable filterByMaximumSize(int maximum_family_size)
    {
        return filterByMaximumSize(maximum_family_size,1);
    }
    
    /**
     * Produces a table that includes only those entries where the sum of family sizes is not too large,
     * and that have at least 1 occurrence in somewhere 
     * 
     * @param maximum_family_size size threshold
     * @param min_lineages lineage count threshold
     * @return a new table composed of families present in at least as many lineages as the lineage count threshold, and at most as many total members as the size threshold
     */
    public DerivedTable filterByMaximumSize(int maximum_family_size, int min_lineages)
    {
        int num_profiles = getNumFamilies();
//        System.out.println("#*OT.DT.fBMS "+num_profiles+" families\tmax family size "+maximum_family_size+"\tmin lineages "+min_lineages);
        boolean[] family_ok = new boolean[num_profiles];
        int num_filtered_families=0;
        for (int pidx=0; pidx<num_profiles; pidx++)
        {
            int[] t = table[pidx];
            int s=0;
            int lg = 0;
//            System.out.println("#*OT.DT.fBMS iter "+pidx+"\ttlen "+t.length);
            for (int j=0; j<t.length; j++)
            {
//                System.out.println("#*OT.DT.fBMS "+pidx+":"+j+"/"+t.length+"\t= "+t[j]+"\ts "+s+"\tlg "+lg);
                s+=t[j];
                if (t[j]!=0)
                    lg++;
            }
            family_ok[pidx] = (lg>=min_lineages && s<=maximum_family_size);
//            System.out.println("#*OT.DT.fBMS "+family_ok[pidx]);

//            Verbose.message("OT.DT.fBMS (<="+maximum_family_size+") \t"+pidx+"/"+family_names[pidx]+"\t"+s+":"+lg+"\t"+family_ok[pidx]+"\t"+num_filtered_families+"\t"+getProfile(pidx).getPatternString());
            if (family_ok[pidx])
                num_filtered_families++;
        }

//        System.out.println("#*OT.DT.fBMS done: "+num_filtered_families+" families kept");

        return tableForFamilies(family_ok);
    }

    /**
     * Selects a subset of families
     * 
     * @param family_ok array of booleans: family_ok[i] specifies if tyyhe i-th family should be included
     * 
     * @return a new table that contains only the selected families
     */
    public DerivedTable tableForFamilies(boolean[] family_ok)
    {
        int num_profiles = getNumFamilies();
        int num_filtered_families=0;
        for (int pidx=0; pidx<num_profiles; pidx++)
            if (family_ok[pidx])
                num_filtered_families++;
        int[][] new_table = new int[num_filtered_families][];
//        String[] new_names = new String[num_filtered_families];
        {
            int family_idx = 0;
            for (int pidx=0; pidx<num_profiles; pidx++)
                if (family_ok[pidx])
                {
                    new_table[family_idx]=table[pidx];
//                    new_names[family_idx]=family_names[pidx];
                    family_idx++;
                }
        }
        DerivedTable T = new DerivedTable(this, family_ok);
        T.setTable(new_table);
        
        return T;
    }
            
    /**
     * Convenience method for selecting a range of families
     * 
     * @param first_family_idx first family index to be included
     * @param last_family_idx last family index to be included 
     * @return a table for families with indexes first_family_idx..last_family_idx (both ends inclusive)
     */
    public OccurrenceTable tableForFamilies(int first_family_idx, int last_family_idx)
    {
        boolean[] family_ok = new boolean[getNumFamilies()];
        for (int i=first_family_idx; i<=last_family_idx; i++)
            family_ok[i]=true;
        return tableForFamilies(family_ok);
    }
    /**
     * Produces a table containing one profile for each family.
     *
     * @return array of tables (one for each family)
     */
    public OccurrenceTable[] allTablesForFamilies()
    {
        int num_families  = getNumFamilies();
        OccurrenceTable[] family_tables = new OccurrenceTable[num_families];
        for (int family_idx=0; family_idx<num_families; family_idx++)
        {
            int[][] new_table = new int[1][];
//            String[] new_names = new String[1];
            new_table[0] = table[family_idx];
//            new_names[0] = family_names[family_idx];
            SingleFamilyTable T = new SingleFamilyTable(this, family_idx);
            T.setTable(new_table);
            // share all propeties
//            T.property_names = property_names;
//            T.family_properties[0] = family_properties[family_idx];
            family_tables[family_idx] = T;
        }
        return family_tables;
    }

    /**
     * A class for different profile entries over a subset of (or all)
     * families and taxa. A transformed table has the same
     * metadata (family names, properties) as the parent table, but has its
     * own data table.
     */
    public static class TransformedTable extends OccurrenceTable
    {
        public TransformedTable(OccurrenceTable parent)
        {
            super();
            this.parent = parent;
            ((OccurrenceTable)this).copyTerminalTaxa(parent);
        }
        private OccurrenceTable parent;

        protected final OccurrenceTable getParent(){return parent;}

        /**
         * Hook for subclasses
         *
         * @param family_idx family index in this filtered table
         * @return family index at the parent table
         */
        protected int getFamilyIndexAtParent(int family_idx)
        {
            return family_idx;
        }

        /**
         * Terminal taxon names are copied from parent
         *
         * @return
         */
        @Override
        public String[] getTerminalTaxonNames()
        {
            return getParent().getTerminalTaxonNames();
        }

        /**
         * Family properties are copied from parent table.
         *
         * @param family_idx
         * @param property_name
         * @return
         */
        @Override
        public String getFamilyProperty(int family_idx, String property_name)
        {
            int parent_idx = getFamilyIndexAtParent(family_idx);
            if (parent_idx<0)
                return super.getFamilyProperty(family_idx, property_name);
            else
                return getParent().getFamilyProperty(parent_idx, property_name);
        }

        /**
         * Family properties are copied from parent table.
         *
         * @param family_idx
         * @param property_idx
         * @return
         */
        @Override
        public String getFamilyProperty(int family_idx, int property_idx)
        {
            int parent_idx = getFamilyIndexAtParent(family_idx);
            if (parent_idx<0)
                return super.getFamilyProperty(family_idx, property_idx);
            else
                return getParent().getFamilyProperty(parent_idx, property_idx);
        }

        /**
         * Family properties are copied from parent table.
         *
         * @return
         */
        @Override
        public String[] getKnownProperties()
        {
            return getParent().getKnownProperties();
        }

        /**
         * Family properties are copied from parent table.
         *
         * @return
         */
        @Override
        public int getKnownPropertiesCount()
        {
            return getParent().getKnownPropertiesCount();
        }

        /**
         * Family properties are copied from parent table.
         *
         * @param property_idx
         * @return
         */
        @Override
        public String getPropertyName(int property_idx)
        {
            return getParent().getPropertyName(property_idx);
        }

        /**
         * Family properties are registered at the parent table.
         *
         * @param property_name
         * @return
         */
        @Override
        public int registerProperty(String property_name)
        {
            return getParent().registerProperty(property_name);
        }

        /**
         * Family properties are set at the parent table.
         *
         * @param family_idx
         * @param property_name
         * @param property_value
         */
        @Override
        public void setFamilyProperty(int family_idx, String property_name, String property_value)
        {
            int parent_idx = getFamilyIndexAtParent(family_idx);
            if (parent_idx<0)
                super.setFamilyProperty(family_idx, property_name,property_value);
            else
                getParent().setFamilyProperty(parent_idx, property_name, property_value);
        }


        /**
         * Family properties are set at the parent table.
         * 
         * @param family_idx
         * @param property_idx
         * @param property_value
         */
        @Override
        public void setFamilyProperty(int family_idx, int property_idx, String property_value)
        {
            int parent_idx = getFamilyIndexAtParent(family_idx);
            if (parent_idx<0)
                super.setFamilyProperty(family_idx, property_idx, property_value);
            else
                getParent().setFamilyProperty(parent_idx, property_idx, property_value);
        }

        /**
         * Family names are copied from the parent table.
         *
         * @param family_idx
         * @return
         */
        @Override
        public String getFamilyName(int family_idx)
        {
            int parent_idx = getFamilyIndexAtParent(family_idx);
            if (parent_idx<0)
                return super.getFamilyName(family_idx);
            else
                return getParent().getFamilyName(parent_idx);
        }

        /**
         * Family names and properties are not initialized here,
         * @param table data table: one entry for each row
         * @param family_names ignored
         *
         */
        @Override
        public void setTable(int[][] table, String[] family_names)
        {
            setTableOnly(table);
            checkMissingEntries();
        }

        @Override
        public int getNumFamilies()
        {
            return parent.getNumFamilies();
        }
    }

    /**
     * Class for a descendant table with a single family.
     */
    public static class SingleFamilyTable extends TransformedTable
    {
        private int family_index_at_parent;
        public SingleFamilyTable(OccurrenceTable parent, int parent_family)
        {
            super(parent);
            this.family_index_at_parent = parent_family;
        }

        @Override
        protected int getFamilyIndexAtParent(int family_idx)
        {
            assert (family_idx==0);
            return family_index_at_parent;
        }

        @Override
        public int getNumFamilies()
        {
            return 1;
        }
    }

    /**
     * A class for a data set over a restricted set of families.
     */
    public static class DerivedTable extends TransformedTable
    {
        private int[] family_index_at_parent;
        private int num_families;
        private int num_added_later;

        public DerivedTable(OccurrenceTable parent, boolean[] family_included)
        {
            super(parent);
            setMapping(family_included);
        }

        @Override
        protected int getFamilyIndexAtParent(int family_idx)
        {
            if (family_idx<num_added_later)
                return -1;
            else
                return family_index_at_parent[family_idx-num_added_later];
        }

        private void setMapping(boolean[] family_ok)
        {
            int num_profiles = getParent().getNumFamilies();
            int num_filtered_families=0;
            for (int pidx=0; pidx<num_profiles; pidx++)
                if (family_ok[pidx])
                    num_filtered_families++;
            family_index_at_parent = new int[num_filtered_families];

            int family_idx = 0;
            for (int pidx=0; pidx<num_profiles; pidx++)
                if (family_ok[pidx])
                    family_index_at_parent[family_idx++]=pidx;
            this.num_families = family_idx;
            this.num_added_later = 0;
//            System.out.println("#*OT.DT.sM "+num_families+" included in derived table");
        }

        @Override
        public int getNumFamilies()
        {
            return this.num_families;
        }

        @Override
        public String[] addProfile(int[] new_profile, String new_profile_name)
        {
            String[] new_profile_names = super.addProfile(new_profile, new_profile_name);
            setNamesOnly(new_profile_names);
            ++num_added_later;
            assert (new_profile_names.length == num_added_later);
            ++num_families;
//            System.out.println("#*OT.DT.aP later added "+num_added_later+" of "+num_families+"\t"+PhyleticProfile.getPatternString(new_profile));

            return new_profile_names;
        }
    }
    
    /**
     * Computes for each t how many family profiles appear exactly t times in the data.
     * 
     * @return array of multiplicities; [t] gives the number of phyletic patterns that occur t times  
     */
    public int[] getProfileMultiplicityDistribution()
    {
        Hashtable<String,Integer> profile_count=new Hashtable<String,Integer>();
        PhyleticProfile[] profiles = getProfiles();
        int max_cnt = 0;
        for (int i=0;i<profiles.length; i++)
        {
            String pattern = profiles[i].getPatternString();
            Verbose.message("OT.gPMD "+pattern+"\t"+profile_count.get(pattern)+"\t"+max_cnt);
            if (profile_count.containsKey(pattern))
            {
                int cnt = profile_count.get(pattern).intValue();
                cnt++;
                if (cnt>max_cnt)
                    max_cnt = cnt;
                profile_count.put(pattern, new Integer(cnt));
            } else
            {
                profile_count.put(pattern, new Integer(1));
                if (max_cnt==0)
                    max_cnt=1;
            }
        }
        int[] multiplicity_distribution = new int[max_cnt+1];
        java.util.Enumeration<String> different_profiles = profile_count.keys();
        while (different_profiles.hasMoreElements())
        {
            String pattern = different_profiles.nextElement();
            int cnt = profile_count.get(pattern).intValue();
            multiplicity_distribution[cnt]++;
        }
        return multiplicity_distribution;
    }

    public double getUnconstrainedLogLikelihood()
    {
        int[] multiplicity_distribution = getProfileMultiplicityDistribution();
        int num_families = 0;
        double support = 0.0;
        for (int t=1; t<multiplicity_distribution.length; t++)
            if (multiplicity_distribution[t]!=0)
            {
                num_families += t*multiplicity_distribution[t];
                double z = t*Math.log(t);
                support += z*multiplicity_distribution[t];
            }
        double E = num_families*Math.log(num_families);
        Verbose.message("OT.gUL nf "+num_families+" sum n_i log n_i = "+support+"\tE "+E);
        return support - E;
    }
    
    /**
     * Minimum number of lineages in which families are present in the input table.
     * Must be 0,1, or 2.
     */
    private static int MIN_PRESENT_LINEAGES = 2;
    
    /**
     * Maximum total number of paralogs in the input tables
     */
    private static int MAX_PARALOGS = 100; //Integer.MAX_VALUE;
    
    /**
     * Test code --- reads a phylogeny and a table, and then writes them to stdout.
     * @param args command line arguments
     * @throws IOException if no file
     * @throws Parser.parseException if bad syntax
     */
    public static void main(String[] args) throws java.io.IOException, Parser.ParseException
    {
        Verbose.setVerbose(false);
        

        int num_switches = 0;
        while (args.length>2*num_switches && args[2*num_switches].startsWith("-"))
        {
            String arg_switch = args[2*num_switches].substring(1);
            if (args.length==2*num_switches+1)
                throw new IllegalArgumentException("Missing argument for switch "+args[2*num_switches]);
            String arg_value = args[2*num_switches+1];
            if (arg_switch.equals("max_paralogs"))
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
        
        if (rest.length != 2)
        {
            System.err.println("Call as java "+OccurrenceTable.class.getCanonicalName()+" phylogeny table");
            System.exit(2008);
        }
        
        
        String tree_file = rest[0];
        String table_file = rest[1];
        
        Reader TR = new java.io.FileReader(tree_file);
        TreeNode root = Parser.readNewick(TR);
        TR.close();

        
        OccurrenceTable T0 = new OccurrenceTable(root.getTraversal().getLeaves());
        Reader BR = new java.io.FileReader(table_file);
        T0.readTable(BR);
        BR.close();
        OccurrenceTable OT = T0.filterByMaximumSize(MAX_PARALOGS, MIN_PRESENT_LINEAGES);
        
        
        System.out.println(OT.getFormattedTable());
        
        //TreeWithRates tree = new TreeWithRates(root);
        //SurvivalProbabilityFactory factory = new SurvivalProbabilityFactory(tree, OT);
        //SurvivalProbabilityFactory.Rates R = factory.getSurvivalProbability(tree);

        //NodeWithRates[] nodes = tree.getDFT();
        //for (int node_idx=0; node_idx<nodes.length; node_idx++)
        //{
        //    NodeWithRates N = nodes[node_idx];
        //    System.out.println("------ node "+N.getId()+"/"+N.getTaxonName());
        //    System.out.println("\textinction "+R.getExtinction(node_idx));
        //    if (!N.isRoot())
        //    {
        //        int[] sizes = factory.getExtremeSizes(node_idx);
        //        for (int n=0; n<sizes.length; n++)
        //        {
        //            System.out.print("\t"+n+"\t0.."+sizes[n]);
        //            double[] sp = R.getSurvival(node_idx,n);
        //            for (int m=0; m<Math.min(sp.length,5); m++)
        //                System.out.print("\t"+sp[m]);
        //            System.out.println();
        //        }
        //    }
        //}
        
        TreeWithRates tree = new TreeWithRates(root);
        NodeWithRates[] leaves = tree.getLeaves();
        PhyleticProfile[] profiles = OT.getProfiles();
        int[] num_singletons = new int[leaves.length];
        int[] num_families = new int[leaves.length];
        int[] num_genes = new int[leaves.length];
        for (int pidx=0; pidx<profiles.length; pidx++)
        {
            int[] pattern = profiles[pidx].getProfile();
            int lineage = -1;
            int num_present=0;
            for (int i=0;i<pattern.length;i++)
                if (pattern[i]>0)
                {
                    lineage = i;
                    num_present++;
                    num_families[i]++;
                    num_genes[i]+=pattern[i];
                }
            if (num_present == 1)
                num_singletons[lineage]++;
        }
        System.out.println("# LS node\ttotal families\tlineage-specific families\tgenes");
        for (int i=0; i<leaves.length; i++)
        {
            NodeWithRates N = leaves[i];
            System.out.println("# LS "+i+"/"+N.newickName()+"\t"+num_families[i]+"\t"+num_singletons[i]+"\t"+num_genes[i]);
        }
        
        //    System.out.println("====== Profile "+pidx);
        //    PhyleticProfile PP = profiles[pidx];
        //    PhyleticProfile.LowerConditionals C = PP.getLowerConditional(R);
        //    for (int node_idx=0; node_idx<nodes.length; node_idx++)
        //    {
        //        NodeWithRates N = nodes[node_idx];
        //        System.out.println("\tnode "+N.getId()+"/"+N.getTaxonName());
        //        double[] L = C.getLikelihoods(node_idx);
        //        System.out.print("\t");
        //        for (int n=0; n<Math.min(L.length,4); n++)
        //            System.out.print("\t"+L[n]);
        //        System.out.println();
        //    }
        //}
        
        int[] freq = OT.getProfileMultiplicityDistribution();
        for (int t=0;t<freq.length; t++ )
            if (freq[t]!=0)
            {
                System.out.println(".: "+freq[t]
                        +"\tprofile"+(freq[t]==1?" ":"s")
                        +" occur"+(freq[t]==1?"s ":"  ")+t+" time"+(t==1?"":"s"));
            }
        
        System.out.println(".:.... unconstrained log-likelihood "+OT.getUnconstrainedLogLikelihood());
        
                
    }
    
    
}
