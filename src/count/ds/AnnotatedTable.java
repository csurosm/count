package count.ds;
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

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Random;
import java.util.function.Predicate;

//import org.json.JSONObject;

import count.io.CommandLine;
import count.io.TableParser;

import static count.io.CommandLine.OPT_RND;

/**
 * Stores a table of family sizes across a number of organisms. 
 *
 * A table is always linked to a list of 
 * terminal taxon names. This list may be initialized 
 * at instantiation. Families may have associated properties beyond then copy numbers and the family name. 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 * @since April 17, 2008, 3:23 PM
 */
public class AnnotatedTable implements ProfileTable
{
    private AnnotatedTable()
    {
        this.terminal_taxa = new HashMap<>();
    }
    
    /**
     * Instantiates a table with the names of these taxa.
     * @param terminal_taxon_names array of terminal taxon names
     */
    public AnnotatedTable(String[] terminal_taxon_names) 
    {
        this();
        initNames(terminal_taxon_names);
    }

    /**
     * Terminal taxon names are stored here
     */
    private final Map<String,Integer> terminal_taxa;
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
     * Taxon names in the order the profiles.
     * (Inverted index for taxon names.)
     * 
     * @return 
     */
    public String[] getTaxonNames()
    {
    	String[] taxon_names = new String[getTaxonCount()];
    	for (Map.Entry<String, Integer> pair: terminal_taxa.entrySet())
			taxon_names[pair.getValue()] = pair.getKey();
    	return taxon_names;
    }
    
    /**
     * Called by instantiation method.
     *
     * @param taxon_names may not be null
     */
    private void initNames(String[] taxon_names)
    {
        terminal_taxa.clear();
        for (int i=0; i<taxon_names.length; i++)
        {
            if (terminal_taxa.containsKey(taxon_names[i]))
                throw new IllegalArgumentException("Repeated taxon name "+taxon_names[i]);
            terminal_taxa.put(taxon_names[i], i);
        }
        
        //System.out.println("#**IT.i "+taxa.length);
        table = null;
        family_names = new String[0];
        family_properties = new Properties[0];
        property_names = new String[0];
    }
    
    
    
    /**
     * Whether this table has missing entries.
     *
     * @return
     */
    public boolean hasMissingEntries()
    {
        return has_missing_entries;
    }    
    
//    public String[] getTaxonNames()
//    {
//        String[] names = new String[terminal_taxa.size()];
//        for (String node_name: terminal_taxa.keySet())
//        {
//            int column_idx = terminal_taxa.get(node_name);
//            names[column_idx] = node_name;
//        }
//        return names;
//    }

    @Override
    public int getTaxonCount()
    {
    	return terminal_taxa.size();
    }

    @Override
    public int getFamilyCount(){ return this.family_names.length;}    

    @Override
    public int[] getFamilyProfile(int family_idx)
    {
        return table[family_idx];
    }
    
     /**
     * Name of the given family
     * @param family_idx index of the family
     * @return name of the family with the given index
     */
    @Override
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
     * Returns the number of families present at a terminal taxon
     * 
     * @param taxon_idx index of the terminal taxon (same order as getTerminalTaxonNames)
     * @return number of families with non-zero size
     */
    public int getNumFamiliesPresent(int taxon_idx)
    {
        int s = 0;
        for (int i=0; i<getFamilyCount(); i++)
        {
        	int[] profile = getFamilyProfile(i);
            if (profile[taxon_idx]>0)
                s++;
        }
        return s;
    }
    
    
    
//    /**
//     * Reads in a table.
//     * @param reader the input reader from which the lines are read
//     * @throws java.io.IOException if there is something wrong with I/O 
//     */
//    public void readTable(Reader reader) throws IOException
//    {
//        readTable(reader, false);
//    }
//    
    
    /**
     * Column header and property for family name. 
     */
    public static final String FAMILY_NAME_PROPERTY = "Family";
    /**
     * Sets the table from externally constructed data. 
     * 
     * @param table family sizes: table[i][j] is the size of the i-th family in the j-th species (this latter indexed as in terminal_taxa)
     * @param family_names names for the families in the order of the indexes into the int[] array; if null, then they are set to are set to "F1,F2,..."
     */
    public void setTable(int[][] table, String[] family_names)
    {
        this.table=table;

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
        property_names[0] = FAMILY_NAME_PROPERTY;
        for (int family_idx = 0; family_idx<nfam; family_idx++)
        {
            Properties prop = new Properties();
            //if (family_names != null)
                prop.setProperty(property_names[0], this.family_names[family_idx]);
            family_properties[family_idx] = prop;
        }
        
        checkMissingEntries();
    }    
    
    public void mapTo(IndexedTree tree)
    {
    	int[] column_order = new int[getTaxonCount()];
    	String[] tree_leaves = tree.getLeafNames();
    	
    	assert (column_order.length == tree_leaves.length);
    	for (int col=0; col<column_order.length; col++)
    	{
    		int old_col = terminal_taxa.get(tree_leaves[col]);
    		column_order[old_col] = col;
    	}
    	for (int f = 0; f<table.length; f++)
    	{
    		int[] old_row = table[f];
    		int[] new_row = new int[old_row.length];
    		for (int col=0; col<column_order.length; col++)
    		{
    			new_row[column_order[col]]=old_row[col];
    			table[f] = new_row;
    		}
    	}
    	for (int leaf=0; leaf<tree_leaves.length; leaf++)
    	{
    		terminal_taxa.put(tree.getName(leaf), leaf);
    	}
    }
    
    private void checkMissingEntries()
    {
        this.has_missing_entries = false;
        for (int i=0; i<getFamilyCount() && !has_missing_entries; i++)
        {
        	int[] row = getFamilyProfile(i);
            for (int j=0; j<row.length && !has_missing_entries; j++)
            	has_missing_entries = (row[j]<0) ;
        }
    }
    
    /**
     * String representation of the pattern: numbers with at least two 
     * digits are enclosed in parentheses, one-digit numbers are simply listed.
     * Example: 017(12)4
     * 
     * @param profile array of sizes 
     * @return a String representation of the extended phyletic pattern
     */
    public static String getPatternString(int[] profile)
    {
        StringBuilder sb = new StringBuilder();
        for (int i=0; i<profile.length; i++)
        {
            if (profile[i]<0)
                sb.append('?');
            else if (profile[i]<10)
                sb.append(Integer.toString(profile[i]));
            else
            {
                sb.append('(');
                sb.append(Integer.toString(profile[i]));
                sb.append(')');
            }
        }
        return sb.toString();
    }
    
    public PhyleticProfile getPhyleticProfile(int family)
    {
    	return new PhyleticProfile(family);
    }
    
    /**
    *
    * Computations over profiles (size distribution across tree leaves)
    * 
    * @since April 19, 2008, 12:00 AM
    * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
    */
   public class PhyleticProfile 
   {
       public PhyleticProfile(int family) 
       {
    	   this(family, getFamilyProfile(family));
       }
       
       private PhyleticProfile(int family, int[] profile)
       {
    	   this.family_idx = family;
    	   this.profile = profile;
    	   
       }
       
       private final int family_idx;
       private final int[] profile;

       /**
        * Profile with which this instance was created
        * 
        * @return the array for the family with which this instance was created
        */
       public int[] getProfile()
       {
           return profile;
       }
       
       /**
        * Profile entry (family size) for a particular node.
        * 
        * @param leaf_idx index of a leaf (0..length of profile-1)
        * @return profile value at leaf_idx
        */
       public int get(int leaf_idx)
       {
           return profile[leaf_idx];
       }
       
       public IntegerOrMissing getValue(int leaf_idx)
       {
           return new IntegerOrMissing(get(leaf_idx));
       }
       
       /**
        * String representation of the pattern: numbers with at least two 
        * digits are enclosed in parentheses, one-digit numbers are simply listed.
        * Example: 017(12)4
        * 
        * @return a String representation of the extended phyletic pattern
        */
       public String getPatternString()
       {
           return AnnotatedTable.getPatternString(getProfile());
       }
       
       @Override 
       public int hashCode()
       {
    	   return Arrays.hashCode(profile);
       }
       
       @Override
       public boolean equals(Object o)
       {
    	   if (o!= null && o instanceof PhyleticProfile)
    	   {
    		   PhyleticProfile other = (PhyleticProfile) o;
    		   return Arrays.equals(this.profile, other.profile);
    	   } else
    		   return super.equals(o);
       }
       
   }
   
   
   private class DerivedTable extends AnnotatedTable
   {
	   DerivedTable()
	   {
		   super(); 
	   }
	   
	   
	   	AnnotatedTable getParentTable()
	   	{
	   		return AnnotatedTable.this;
	   	}
	   
  		@Override
   	    public String[] getTaxonNames()
   	    {
   			return AnnotatedTable.this.getTaxonNames();
   	    }
   	    @Override
   	    public int getTaxonCount()
   	    {
   	    	return AnnotatedTable.this.getTaxonCount();
   	    }
   	    @Override
   	    public int getFamilyCount()
   	    { 
   	    	return AnnotatedTable.this.getFamilyCount();
   	    }    
   	    @Override
   	    public int[] getFamilyProfile(int f)
   	    {
   	    	return AnnotatedTable.this.getFamilyProfile(f);
   	    }
   	    @Override
   	    public String getFamilyName(int f)
   	    {
   	    	return AnnotatedTable.this.getFamilyName(f);
   	    }
   	    @Override
   	    public String[] getKnownProperties()
   	    {
   	    	return AnnotatedTable.this.getKnownProperties();
   	    }
   	    @Override
   	    public String getPropertyName(int property_idx)
   	    {
   	    	return AnnotatedTable.this.getPropertyName(property_idx);
   	    }
   	    @Override
   	    public int getKnownPropertiesCount()
   	    {
   	    	return AnnotatedTable.this.getKnownPropertiesCount();
   	    }
   	    @Override
   	    public String getFamilyProperty(int f, String property_name)
   	    {
   	    	return AnnotatedTable.this.getFamilyProperty(f, property_name);
   	    }
   	    @Override
   	    public String getFamilyProperty(int f, int property_idx)
   	    {
   	    	return AnnotatedTable.this.getFamilyProperty(f, property_idx);
   	    }
   	    @Override
   	    public void setFamilyProperty(int f, String property_name, String property_value)
   	    {
   	    	AnnotatedTable.this.setFamilyProperty(f, property_name, property_value);
   	    }
   	    @Override
   	    public void setFamilyProperty(int f, int property_idx, String property_value)
   	    {
   	    	AnnotatedTable.this.setFamilyProperty(f, property_idx, property_value);
   	    }
   	    @Override
   	    public int registerProperty(String property_name)
   	    {
   	    	return AnnotatedTable.this.registerProperty(property_name);
   	    }
   	    @Override
   	    public int getNumFamiliesPresent(int taxon_idx)
   	    {
   	    	return AnnotatedTable.this.getNumFamiliesPresent(taxon_idx);

   	    }
   	    @Override
   	    public void setTable(int[][] table, String[] family_names)
   	    {
   	    	throw new UnsupportedOperationException();
   	    }
   	    
   	    
   }
   
   public AnnotatedTable filteredTable(Predicate<Integer> family_selected)
   {
	   int nF = getFamilyCount();
	   boolean[] selected = new boolean[nF];
	   int f=0;
	   int s=0;
	   while (f<nF)
	   {
		   if (selected[f] = family_selected.test(f)) // test evaluated only once 
			   s++;
		   f++;
	   }
	   int[] selected_rows = new int[s];
	   assert (f==nF);
	   while (s>0) // selected_rows[s..] is filled
	   {
		   do --f; while (!selected[f]);
		   --s;
		   selected_rows[s] = f;
	   }
	   return filteredTable(selected_rows);
   }
   
   public AnnotatedTable filteredTable(int[] selected_rows)
   {
   	
   	class FilteredTable extends DerivedTable
   	{
   		int[] parent_rows;
   		FilteredTable()
   		{
   			this.parent_rows = selected_rows.clone();
   			checkMissingEntries();
   		}
   	    @Override
   	    public int getFamilyCount()
   	    { 
   	    	return parent_rows.length;
   	    }    
   	    @Override
   	    public int[] getFamilyProfile(int family_idx)
   	    {
   	    	int f = parent_rows[family_idx];
   	    	return AnnotatedTable.this.getFamilyProfile(f);
   	    }
   	    @Override
   	    public String getFamilyName(int family_idx)
   	    {
   	    	int f = parent_rows[family_idx];
   	    	return AnnotatedTable.this.getFamilyName(f);
   	    }
   	    @Override
   	    public String getFamilyProperty(int family_idx, String property_name)
   	    {
   	    	int f = parent_rows[family_idx];
   	    	return AnnotatedTable.this.getFamilyProperty(f, property_name);
   	    }
   	    @Override
   	    public String getFamilyProperty(int family_idx, int property_idx)
   	    {
   	    	int f = parent_rows[family_idx];
   	    	return AnnotatedTable.this.getFamilyProperty(f, property_idx);
   	    }
   	    @Override
   	    public void setFamilyProperty(int family_idx, String property_name, String property_value)
   	    {
   	    	int f = parent_rows[family_idx];
   	    	AnnotatedTable.this.setFamilyProperty(f, property_name, property_value);
   	    }
   	    @Override
   	    public void setFamilyProperty(int family_idx, int property_idx, String property_value)
   	    {
   	    	int f = parent_rows[family_idx];
   	    	AnnotatedTable.this.setFamilyProperty(f, property_idx, property_value);
   	    }
   	}
   	
   	
   	
   	return new FilteredTable();
   }
   
   public AnnotatedTable mappedToTree(IndexedTree tree)
   {
	   String[] table_taxon_names = getTaxonNames();
	   Map<String, Integer> table_column_index = new HashMap<>();
	   for (int t=0; t<table_taxon_names.length; t++)
		   table_column_index.put(table_taxon_names[t], t);
	   String[] tree_taxon_names = tree.getLeafNames();
//	   System.out.println("#**AT.mTT table "+Arrays.toString(table_taxon_names)
//	   	+"\ttree "+Arrays.toString(tree_taxon_names));

	   
	   //	   if (tree_taxon_names.length != table_taxon_names.length)
//	   {
//		   throw new IllegalArgumentException("Cannot match all columnns: requested "
//				   	+Arrays.toString(tree_taxon_names)
//				   	+";\tknown "+Arrays.toString(table_taxon_names));
//	   }
	   int[] column_order = new int[tree_taxon_names.length];
	   Arrays.fill(column_order, -1);
	   boolean has_unmapped_leaves = false;
	   for (int c=0; c<tree_taxon_names.length; c++)
	   {
		   String name = tree_taxon_names[c];
		   if (table_column_index.containsKey(name))
		   {	
			   if (column_order[c]!=-1)
				   has_unmapped_leaves=true;
			   column_order[c]=table_column_index.get(name);
		   }
		   else
		   {
			   has_unmapped_leaves = true;
			   column_order[c]=-1;
		   }
	   }
	   if (has_unmapped_leaves)
	   {
		   System.err.println("Cannot match all leaf names to columnns: requested "
				   	+Arrays.toString(tree_taxon_names)
				   	+";\tknown "+Arrays.toString(table_taxon_names));
		   throw new IllegalArgumentException();
	   }
	   
	   
	   class MappedTable extends DerivedTable
	   {
		   private final int[] column_order;
		   /**
		    * 
		    * 
		    * @param column_order every entry must be a column in the parent table
		    */
		   MappedTable(int[] column_order)
		   {
			   boolean same_order = true;
			   for (int c=0; c<column_order.length && same_order; c++)
				   same_order = column_order[c]==c;
			   this.column_order = same_order?null:column_order;
			   checkMissingEntries();
		   }
		   @Override
		   public String[] getTaxonNames()
		   {
			   String[] orig_taxons = AnnotatedTable.this.getTaxonNames();
			   String[] our_taxons;
			   if (column_order==null)
				   our_taxons = orig_taxons;
			   else
			   {
				   our_taxons = new String[column_order.length];
				   for (int c=0; c<column_order.length; c++)
					   our_taxons[c] = orig_taxons[column_order[c]];
			   }
			   return our_taxons;
		   }
		   
		   @Override
		   public int getTaxonCount()
		   {
			   return column_order==null?AnnotatedTable.this.getTaxonCount():column_order.length;
			   
		   }
		   @Override
		   public int[] getFamilyProfile(int f)
		   {
			   int[] orig_profile = AnnotatedTable.this.getFamilyProfile(f);
			   int[] our_profile = null;
			   if (column_order==null)
				   our_profile = orig_profile;
			   else
			   {
				   our_profile = new int[column_order.length];
				   for (int c=0; c<column_order.length; c++)
					   our_profile[c]=orig_profile[column_order[c]];
			   }
			   return our_profile;
		   }
		   
		   @Override 
		   public boolean isBinaryTable()
		   {
			   return AnnotatedTable.this.isBinaryTable();
		   }
	   } // class
	   MappedTable mapped = new MappedTable(column_order);
	   
//	   System.out.println("#**AT.mTT order "+Arrays.toString(column_order)
//	   	+"\tsize "+mapped.getFamilyCount()+"*"+mapped.getTaxonCount());
	   
	   return mapped;
   }

   public AnnotatedTable binaryTable()
   {
	   class BinaryTable extends DerivedTable
	   {
		   /**
		    * Precomputed table so that fetching 
		    * a profile is faster. 
		    */
		   private final int[][] binary_table;
		   
		   BinaryTable()
		   {
			   binary_table = new int[getFamilyCount()][];
			   fillTable();
		   }
		   
		   /**
		    * Fills {@link #binary_table}
		    * by detecting identical binary profiles. 
		    */
		   private void fillTable()
		   {
			   Map<PhyleticProfile, Integer> binary_profiles = new HashMap<>();
			   int nF = getFamilyCount();
			   for (int f=0; f<nF; f++)
			   {
				   int[] original_profile = AnnotatedTable.this.getFamilyProfile(f);
				   int[] profile01 = new int[original_profile.length];
				   for (int c=0; c<original_profile.length; c++)
					   profile01[c]=Integer.min(1,original_profile[c]);
				   PhyleticProfile this_row = new PhyleticProfile(f,profile01);
				   if (binary_profiles.containsKey(this_row))
				   {
					   int same_f = binary_profiles.get(this_row);
					   binary_table[f]=binary_table[same_f];
				   } else
				   {
					   binary_table[f]=profile01;
					   binary_profiles.put(this_row, f);
				   }
			   }
		   }
		   
	   	    @Override
	   	    public int[] getFamilyProfile(int f)
	   	    {
	   	    	return binary_table[f];
	   	    }
		   
	   	    @Override
	   	    public boolean isBinaryTable()
	   	    {
	   	    	return true;
	   	    }
	   	    
	   }
	   
	   return new BinaryTable();
   }
   

   public AnnotatedTable bootstrap(Random RND)
   {
	   int num_families = getFamilyCount();
	   int[] bootstrap_rows = new int[num_families];
	   for (int f=0; f<num_families; f++)
	   {
		   int row = RND.nextInt(num_families);
		   bootstrap_rows[f]=row;
	   }
	   return filteredTable(bootstrap_rows);
   }
   
//   public JSONObject toJSON()
//   {
//	   
//   }

   public static void main(String[] args) throws Exception
   {
		PrintStream out = System.out;
	    out.println(CommandLine.getStandardHeader(AnnotatedTable.class));
	    out.println(CommandLine.getStandardRuntimeInfo(AnnotatedTable.class, args));
	    CommandLine cli = new CommandLine(args,AnnotatedTable.class, 2) ;
	    AnnotatedTable table = cli.getTable();
	    Phylogeny phylo = cli.getTree();
	   
	    if (cli.getOptionValue(OPT_RND)!=null)
	    {
//	    	int seed = cli.getOptionInt(OPT_RND, 0);
	    	Random RND = cli.getOptionRND(out);
//	    			seed==0?new Random():new Random(seed);
		    table = table.bootstrap(RND);
	    }
	    out.println(TableParser.getFormattedTable(table, true));
	    
	    int[] profile_distribution = new int[phylo.getNumLeaves()+1];
	    for (int f=0; f<table.getFamilyCount(); f++)
	    {
	    	int num_lineages = table.getLineageCount(f);
	    	profile_distribution[num_lineages]++;
	    }
		out.println("#DISTRIBUTION\tsize\tn");
	    for (int size=0; size<profile_distribution.length; size++)
	    {
	    	out.println("#DISTRIBUTION\t"+size+"\t"+profile_distribution[size]);
	    }
	    
	    Arrays.fill(profile_distribution, 0);
	    for (int f=0; f<table.getFamilyCount(); f++)
	    {
	    	int[] copies = table.getFamilyProfile(f);
	    	for (int leaf=0; leaf<phylo.getNumLeaves(); leaf++)
	    	{
	    		int c = copies[leaf];
	    		if (c>0) profile_distribution[leaf]++;
	    	}
	    }
		out.println("#FAMILIES\tleaf\tn");
    	for (int leaf=0; leaf<phylo.getNumLeaves(); leaf++)
    	{
    		out.println("#FAMILIES\t"+leaf+"\t"+profile_distribution[leaf]+"\t"+phylo.getIdent(leaf));
    	}
   }
}