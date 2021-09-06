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
package count.ds;

import java.util.HashMap;
import java.util.Map;
import java.util.Properties;

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
     * @param taxa may be null
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
        for (int i=0; i<table.length; i++)
            if (table[i][taxon_idx]>0)
                s++;
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
    
    private void checkMissingEntries()
    {
        this.has_missing_entries = false;
        for (int i=0; i<table.length && !has_missing_entries; i++)
            for (int j=0; j<table[i].length && !has_missing_entries; j++)
                has_missing_entries = (table[i][j]<0) ;
    }
    
}