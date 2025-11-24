package count.io;

/*
 * Copyright 2023 Mikl&oacute;s Cs&#369;r&ouml;s.
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


/**
 * Defines the XML tags and attributes used in Count 
 *
 */
public interface CountXML 
{
	public static final String PREFIX_TREEID="T";
	public static final String PREFIX_RATESID="R";
	public static final String PREFIX_DATAID="D";
	
	public static final String EMT_ROOT = "CountML";
	
    public static final String EMT_SESSION = "session";
	/**
	 * Section tag for phylogenies
	 */
    public static final String EMT_TREE = "tree";
    /**
     * Section tag for data tables and analysis tables 
     */
    public static final String EMT_DATA = "data";
    /**
     * Section tag for rate models
     */
    public static final String EMT_RATES = "rates";

    /**
     * Tag for one rate model
     */
    public static final String EMT_RATEMODEL = "model";

    /**
     * Tag for one table
     */
    public static final String EMT_TABLE = "table";
    
    /**
     * Tag for non-table items (ancestral reconstructions) in data section
     */
    public static final String EMT_DATAITEM = "item";
    /**
     * Legacy tag for data item (Count 2023 recomputes all values)
     */
    public static final String EMT_POSTERIORS = "posteriors";
    
    
    public static final String ATT_TYPE = "type";
    public static final String ATT_ID = "id";
    public static final String ATT_PARENT = "parent";
    public static final String ATT_NAME = "name";
    public static final String ATT_RATES = EMT_RATES; 
    public static final String ATT_TREE = EMT_TREE;
    public static final String ATT_TABLE =  EMT_TABLE;
    // history 
    public static final String ATT_SELECTION = "selection";
    public static final String ATT_CONTROLS = "controls";
    
    // in numerical parsimony
    public static final String ATT_GAIN = CommandLine.OPT_GAIN;
    public static final String ATT_LOSS = CommandLine.OPT_LOSS;
    public static final String ATT_DUPLICATION = CommandLine.OPT_DUPLICATION;
    // in Dollo
    public static final String ATT_FIX_ROOT = "fixroot";
    // in simulation
    public static final String ATT_RND = CommandLine.OPT_RND;
    public static final String ATT_ROWCOUNT = CommandLine.OPT_N;
    
    public static final String ATT_MINCOPY = CommandLine.OPT_MINCOPY;
    
    public static final String ATT_BINARY = "isbinary";
    
    // general setup
    public static final String ATT_THREADS = CommandLine.OPT_THREADS;
    public static final String ATT_TRUNCATE = CommandLine.OPT_TRUNCATE;
    
    public static final String ATT_DATE = "date";
    public static final String ATT_VERSION = "version";
    
    public static final String XML_DECLARATION = "<?xml version=\"1.0\" encoding=\"us-ascii\" ?>";
    
}
