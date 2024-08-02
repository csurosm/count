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

import static count.io.CommandLine.OPT_OUTPUT;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;
import java.io.StringReader;
import java.io.ByteArrayOutputStream;

import java.util.Map;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import count.ds.AnnotatedTable;
import count.ds.Phylogeny;
import count.ds.TreeComparator;

/**
 * 
 * Parsing and writing methods for copy-number table data
 * (stored as tab-delimited text). 
 * 
 * @author csuros
 *
 */
public class TableParser 
{
	private TableParser()
	{
		
	}
	
    public static String MISSING_ENTRY = "?";

    /*
     * Reads in an AnnotatedTable using a reader.
     * 
     * @param taxon_names column headers recognized as containing family size 
     * @param reader the input reader from which the lines are read
     * @param includes_properties whether columns that do not correspond to terminal taxa should be considered as property columns
     * 
     * @throws java.io.IOException if there is something wrong with I/O 
     */
	public static AnnotatedTable readTable(String[] taxon_names,
			Reader reader, boolean includes_properties) throws IOException
	{
        AnnotatedTable table=null;
        final Map<String,Integer> terminal_taxa = new HashMap<>();
        if (taxon_names != null)
        {
        	table = new AnnotatedTable(taxon_names);
	        for (int i=0; i<taxon_names.length; i++)
	        {
	            if (terminal_taxa.containsKey(taxon_names[i]))
	                throw new IllegalArgumentException("Repeated taxon name "+taxon_names[i]);
	            terminal_taxa.put(taxon_names[i], i);
	        }
        }
		
        List<String> property_names_in_input = new ArrayList<>();
        List<int[]> parsed_sizes = new ArrayList<>();
        List<String> parsed_family_names = new ArrayList<>();
        List<List<String>> parsed_family_properties = new ArrayList<>();

        BufferedReader R	= (reader instanceof BufferedReader)
        					?(BufferedReader)reader
        					:new BufferedReader(reader);
        String line;
        int[] taxon_order = null; // how to permute the columns to match terminal_taxa: -1 denotes unknown
        do 
        {
            line = R.readLine();
            if (line == null || line.trim().length()==0 || line.startsWith("#"))
                continue;
            
            String[] fields=line.split("\\t");
            if (taxon_order == null)
            {
                // header line
                taxon_order=new int[fields.length-1]; // column 0 is for family name 
                
                if (taxon_names ==null)
                {
                	for (int j=1; j<fields.length; j++)
	                {
                		String taxon = fields[j];
                		if (terminal_taxa.containsKey(taxon))
        	                throw new IllegalArgumentException("Repeated taxon name "+taxon);
                		terminal_taxa.put(taxon,  j-1);
	                }
                	taxon_names = new String[terminal_taxa.size()];
                	for (String taxon: terminal_taxa.keySet())
                	{
                		int j = terminal_taxa.get(taxon);
                		taxon_names[j] = taxon;
                		taxon_order[j]=j;
                	}
                	table = new AnnotatedTable(taxon_names);
                } else // order/filtering supplied on input
                {
	                for (int j=1; j<fields.length; j++)
	                {
	                    if (terminal_taxa.containsKey(fields[j])) // we know this guy
	                    {
	                        int taxon_idx = terminal_taxa.get(fields[j]);
	                        taxon_order[j-1]=taxon_idx;
	                    } else 
	                    {
	                        taxon_order[j-1]=-1;
	                        if (includes_properties)
	                        {	
	                        	String prop_name = "".equals(fields[j])?"column"+Integer.toString(j):fields[j];
	                            property_names_in_input.add(prop_name);
	                        }
	                    }
	                }
                }
            } else
            {
                // data lines
                List<String> properties_in_line = (includes_properties?new ArrayList<>():null);
                
                String family = fields[0]; 
                int[] copies = new int[terminal_taxa.size()];
                for (int j=1; j<fields.length; j++)
                {
                    int idx = taxon_order[j-1];
                    if (idx != -1)
                    {
                    	if (MISSING_ENTRY.equals(fields[j])
                                || fields[j].trim().equals("")) // or all whitespace
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
                parsed_family_names.add(family);
                parsed_sizes.add(copies);
                if (includes_properties)
                    parsed_family_properties.add(properties_in_line);
            }
        } while (line !=null);
        
        R.close();
        table.setTable(parsed_sizes.toArray(new int[0][]),parsed_family_names.toArray(new String[0]));
        
//        // DEBUG
//        System.out.println("#**TP.rT nF "+parsed_sizes.size()+"\tnames "+parsed_family_names.size()
//        		+"\tprops "+property_names_in_input.size());
        
        if (includes_properties)
        {
            for (int j=0; j<property_names_in_input.size(); j++)
                table.registerProperty(property_names_in_input.get(j));
            int num_properties = property_names_in_input.size();
            int num_families = parsed_family_properties.size();
            for (int family_idx = 0; family_idx<num_families; family_idx++)
            {
                List<String> V = parsed_family_properties.get(family_idx);
                for (int j=0; j<num_properties; j++)
                    table.setFamilyProperty(family_idx, j+1, V.get(j)); // property 0 is family name that is already there
            }
        }        
	    return table;
	}
	
	public static Map<String,String> readTwoColumnMap(BufferedReader reader) throws IOException
	{
		Map<String,String> map = new HashMap<>();
		
		String separator = null;
		String line;
		do
		{
			line = reader.readLine();
            if (line == null || line.trim().length()==0 || line.startsWith("#"))
                continue;
            String[] fields;
            if (separator == null)
            {
            	String[] tfields = line.split("\t");
            	String[] cfields = line.split(",");
            	if (tfields.length<cfields.length)
            	{
            		separator = ",";
            		fields = cfields;
            	} else if (1<tfields.length)
            	{
            		separator = "\t";
            		fields = tfields;
            	} else
            	{
            		throw new IOException("Cannot determine table format: neither tab-, nor comma-separated.");
            	}
            } else
            {
            	fields = line.split(separator);
            }
            String key = fields[0];
            String value = fields[1];
            if (!"".equals(key) && !"".equals(value))
            	map.put(key, value);
		} while (line != null);
		return map;
	}

	/**
	 * Reads membership data from COG.csv files.
	 * 
	 * @param table_reader
	 * @param genome_id
	 * @param unique
	 * @return
	 * @throws IOException
	 */
	public static AnnotatedTable readMembershipData(BufferedReader table_reader
			, Map<String,String> genome_id, boolean unique) throws IOException
	{
		return readMembershipData(table_reader, 0, 1, 6, 7, 3, genome_id, null, unique);
	}
	
	
	
	/**
	 * Generic parsing for COG-style data table: 1 gene-family membership per line.  
	 * 
	 * @param table_reader
	 * @param gene_column
	 * @param genome_column
	 * @param family_column
	 * @param membership_column column for membership type (0-3), set to -1 if unused
	 * @param max_membership maximum membership type for inclusion
	 * @param genome_id mapping for genome ids to taxon ids (shared between table and phylogeny); null if identity mapping
	 * @param unique whether only first gene annotation is kept
	 * @return
	 * @throws IOException
	 */
	public static AnnotatedTable readMembershipData(BufferedReader table_reader
			, int gene_column, int genome_column, int family_column
			, int membership_column
			, int max_membership
			, Map<String,String> genome_id
			, String[] selected_taxa
			, boolean unique) throws IOException
	{
		// TODO add membership column 
		
		
		boolean ident_genome_id = genome_id == null;
		final Map<String,Integer> terminal_taxa = new HashMap<>();
		int num_genomes = 0;
		if (ident_genome_id)
		{
			genome_id = new HashMap<>(); // we'll fill it up
		} else if (selected_taxa==null)
		{
			for (String txid: genome_id.values())
			{
				terminal_taxa.put(txid, num_genomes++);
			}
		}
		if (selected_taxa != null)
		{
			num_genomes = selected_taxa.length;
			for (int t=0; t<num_genomes; t++)
			{
				terminal_taxa.put(selected_taxa[t], t);
			}
		}
		
		// family copy numbers
		Map<String, int[]> family_profiles = new HashMap<>();
		Map<String,String> gene_family_memberships = new HashMap<>(); // tracking multiple annotations for same gene
		
		int num_orphan_genes = 0;
    	ByteArrayOutputStream bstream = new ByteArrayOutputStream(256);
    	PrintStream orphan_name_formatter = new PrintStream(bstream);
		
		String separator = null;
		String line;
		do
		{
			line = table_reader.readLine();
            if (line == null || line.trim().length()==0 || line.startsWith("#"))
                continue;
            String[] fields;
            if (separator == null)
            {
            	String[] tfields = line.split("\t");
            	String[] cfields = line.split(",");
            	if (tfields.length<cfields.length)
            	{
            		separator = ",";
            		fields = cfields;
            	} else if (1<tfields.length)
            	{
            		separator = "\t";
            		fields = tfields;
            	} else
            	{
            		throw new IOException("Cannot determine table format: neither tab-, nor comma-separated.");
            	}
            } else
            {
            	fields = line.split(separator);
            }
            String gene_id = fields[gene_column];
            String genome = fields[genome_column];
            String family = null;
            if (family_column<fields.length)
            	family=fields[family_column];
            int member = 0;
            if (0<= membership_column && membership_column<fields.length)
            {
            	member = Integer.parseInt(fields[membership_column]);
            }
            
            if (unique && gene_family_memberships.containsKey(gene_id))
            	continue;
            if (max_membership < member)
            	continue;

            if (family == null || "".equals(family))
            {
            	num_orphan_genes++;
            	bstream.reset();
            	orphan_name_formatter.printf("orphan%05d", num_orphan_genes);
            	family = bstream.toString();
            }
            int[] copies;
        	if (family_profiles.containsKey(family))
        	{
        		copies = family_profiles.get(family);
        	} else
        	{
        		copies = new int[num_genomes];
        	}
        	gene_family_memberships.put(gene_id, family);
            
            String txid;
            if (ident_genome_id)
            {
            	txid = genome;
            	if (!genome_id.containsKey(genome))
            	{
            		genome_id.put(genome, txid);
            		if (selected_taxa==null)
            		{
            			if (terminal_taxa.containsKey(txid))
            				throw new IllegalArgumentException("Repeated taxon identifier "+txid);
            			terminal_taxa.put(txid, num_genomes++);
            		}
            	}
            } else
            {
            	txid = genome_id.get(genome);
            }
            
//            System.out.println("#**TP.rMD family "+family+"\tgene "+gene_id+"\tgenome "+genome+"\ttxid "+txid+"\torder "+terminal_taxa.get(txid));
            
            
            if (txid!=null && terminal_taxa.containsKey(txid))
            {
            	int taxon_order = terminal_taxa.get(txid);
            	if (family_profiles.containsKey(family))
            	{
            		copies = family_profiles.get(family);
            	} else
            	{
            		copies = new int[num_genomes];
            	}
            	if (copies.length<num_genomes)
            		copies = Arrays.copyOf(copies, num_genomes);
            	copies[taxon_order]++;
            }
        	family_profiles.put(family, copies);
        	
        	
		} while (line != null);
		
		// prepare the table data
		String[] taxon_names = new String[num_genomes];
		for (String txid: terminal_taxa.keySet())
		{
			int taxon_order = terminal_taxa.get(txid);
			taxon_names[taxon_order] = txid;
		}
		AnnotatedTable readMembershipData = new AnnotatedTable(taxon_names);
		
		int num_families = family_profiles.size();
		int [][] table_entries = new int[num_families][];
		String[] family_names =  family_profiles.keySet().toArray(new String[0]) ; //Arrays.sort();
		Arrays.sort(family_names);
		
		// make sure the profiles have the same length
		for (int f=0; f<family_names.length; f++)
		{
			table_entries[f] = Arrays.copyOf(family_profiles.get(family_names[f]), num_genomes);
		}
		readMembershipData.setTable(table_entries, family_names);
		
		return readMembershipData;
	}
	

	
	/**
	 * 
	 * 
	 * Match class 0: footprint covers most of the protein and most of domain
	 * Match class 1: footprint covers most of the domain and part of the protein
	 * Match class 2: footprint covers most of the protein and part of the domain
	 * Match class 3: partial match on both protein and domain
	 * 
	 * Match class 0 and 1 are immediate evidence of domain occurrence, each increment the copy number 
	 * Match class 2 is an evidence of copy number 1 
	 * Match class 3 is partial: multiple hits are combined into a single occurrence but only if  
	 * 
	 * @param table_reader
	 * @param gene_column
	 * @param genome_column
	 * @param family_column
	 * @param match_class_column
	 * @param genome_id
	 * @param selected_taxa
	 * @param unique
	 * @return
	 * @throws IOException
	 */
	
//	public static AnnotatedTable readCOGData(BufferedReader table_reader
//			, int gene_column, int genome_column, int family_column, int match_class_column
//			, Map<String,String> genome_id
//			, String[] selected_taxa
//			, boolean unique) throws IOException
//	{
//		boolean ident_genome_id = genome_id == null;
//		final Map<String,Integer> terminal_taxa = new HashMap<>();
//		int num_genomes = 0;
//		if (ident_genome_id)
//		{
//			genome_id = new HashMap<>(); // we'll fill it up
//		} else if (selected_taxa==null)
//		{
//			for (String txid: genome_id.values())
//			{
//				terminal_taxa.put(txid, num_genomes++);
//			}
//		}
//		if (selected_taxa != null)
//		{
//			num_genomes = selected_taxa.length;
//			for (int t=0; t<num_genomes; t++)
//			{
//				terminal_taxa.put(selected_taxa[t], t);
//			}
//		}
//		
//		// family copy numbers
//		Map<String, int[]> family_profiles = new HashMap<>();
//		Map<String,String> gene_family_memberships = new HashMap<>(); // tracking multiple annotations for same gene
//		
//		// naming orphan genes
//		int num_orphan_genes = 0;
//    	ByteArrayOutputStream bstream = new ByteArrayOutputStream(256);
//    	PrintStream orphan_name_formatter = new PrintStream(bstream);
//
//		String separator = null;
//		String line;
//		do
//		{
//			line = table_reader.readLine();
//            if (line == null || line.trim().length()==0 || line.startsWith("#"))
//                continue;
//            String[] fields;
//            if (separator == null)
//            {
//            	String[] tfields = line.split("\t");
//            	String[] cfields = line.split(",");
//            	if (tfields.length<cfields.length)
//            	{
//            		separator = ",";
//            		fields = cfields;
//            	} else if (1<tfields.length)
//            	{
//            		separator = "\t";
//            		fields = tfields;
//            	} else
//            	{
//            		throw new IOException("Cannot determine table format: neither tab-, nor comma-separated.");
//            	}
//            } else
//            {
//            	fields = line.split(separator);
//            }
//            String gene_id = fields[gene_column];
//            String genome = fields[genome_column];
//        	int match_class = 0;
//        	if (0<=match_class_column && match_class_column<fields.length)
//        	{
//        		match_class = Integer.parseInt(fields[match_class_column]);
//        	}
//            String family = null;
//            if (family_column<fields.length)
//            	family=fields[family_column];
//            if (family == null || "".equals(family))
//            {
//            	num_orphan_genes++;
//            	bstream.reset();
//            	orphan_name_formatter.printf("orphan%05d", num_orphan_genes);
//            	family = bstream.toString();
//            }
//            int[] copies;
//        	if (family_profiles.containsKey(family))
//        	{
//        		copies = family_profiles.get(family);
//        	} else
//        	{
//        		copies = new int[num_genomes];
//        	}
//        	gene_family_memberships.put(gene_id, family);
//        	family_profiles.put(family, copies);
//        	
//            
//            
//		} while (line != null); 	
//    	
//	
//	}
	
	/**
	 * Parses MCL-style data of clustering results: enumerated members in each line 
	 * 
	 * @param reader
	 * @param gene_taxon_mapping
	 * @param selected_taxa
	 * @throws IOException
	 */
	public static AnnotatedTable readClusteringData(BufferedReader reader
				, Map<String,String> gene_taxon_mapping
				, String[] selected_taxa)
				throws IOException
	{
		final Map<String,Integer> terminal_taxa = new HashMap<>();
		int num_genomes;
		if (selected_taxa != null)
		{
			num_genomes = selected_taxa.length;
			for (int t=0; t<num_genomes; t++)
			{
				terminal_taxa.put(selected_taxa[t], t);
			}
		} else
		{
			num_genomes=0;
			for (String txid: gene_taxon_mapping.values())
			{
				if (!terminal_taxa.containsKey(txid))
					terminal_taxa.put(txid, num_genomes++);
			}
			selected_taxa = new String[num_genomes];
			for (String txid: terminal_taxa.keySet())
				selected_taxa[terminal_taxa.get(txid)]=txid;
		}
		assert (num_genomes == terminal_taxa.size());			

		List<int[]> family_profiles = new ArrayList<>();
		String separator = null;
		String line;
		do
		{
    		line = reader.readLine();
            if (line == null || line.trim().length()==0 || line.startsWith("#"))
                continue;
            String[] fields;
            if (separator == null)
            {
            	String[] tfields = line.split("\t");
            	String[] cfields = line.split(",");
            	if (tfields.length<cfields.length)
            	{
            		separator = ",";
            		fields = cfields;
            	} else if (1<tfields.length)
            	{
            		separator = "\t";
            		fields = tfields;
            	} else
            	{
            		throw new IOException("Cannot determine table format: neither tab-, nor comma-separated.");
            	}
            } else
            {
            	fields = line.split(separator);
            }
            int[] copies = new int[num_genomes];
            for (int i=0; i<fields.length; i++)
            {
            	String gene = fields[i];
            	String txid = gene_taxon_mapping.get(gene);
            	if (txid!=null && terminal_taxa.get(txid)!=null)
            	{
            		int taxon_order = terminal_taxa.get(txid);
            		copies[taxon_order]++;
            	}
            }
            family_profiles.add(copies);
		} while (line!=null);
		
		// construct family names
		List<String> family_names = new ArrayList<>();
		int ndigits = 2;
		int max_nfam = 100;
		while (max_nfam<family_profiles.size())
		{
			ndigits++;
			max_nfam *= 10;
		}
		String counter_format = "%0"+ndigits+"d";
		int og_cnt = 0;
		int lse_cnt = 0;
		int single_cnt = 0;
		ByteArrayOutputStream bstream = new ByteArrayOutputStream(256);
    	PrintStream name_formatter = new PrintStream(bstream);
    	Map<String, Integer> original_index=new HashMap<>();
    	
		for (int f=0; f<family_profiles.size(); f++)
		{
			int[] copies = family_profiles.get(f);
			int num_lineages = 0;
			int num_copies = 0;
			for (int t=0; t<copies.length && num_lineages<2; t++)
			{
				if (0<copies[t]) num_lineages++;
				num_copies += copies[t];
			}
			int cnt;
			String type;
			
			if (1<num_lineages)
			{
				type = "og";
				cnt = ++og_cnt;
			} else if (1<num_copies)
			{
				type = "lse";
				cnt = ++lse_cnt;
			} else
			{
				type = "single";
				cnt = ++single_cnt;
			}
			String name_format = type + counter_format;
			bstream.reset();
			name_formatter.printf(name_format, cnt);
			String family = bstream.toString();
			family_names.add(family);
			original_index.put(family, f);
		}
		Collections.sort(family_names);
		int[][] table_entries = new int[family_profiles.size()][];
		for (int f=0; f<table_entries.length; f++)
		{
			String family = family_names.get(f);
			table_entries[f] = family_profiles.get(original_index.get(family));
		}
		
		AnnotatedTable readClusteringData = new AnnotatedTable(selected_taxa);
		readClusteringData.setTable(table_entries,family_names.toArray(new String[0]));
		return readClusteringData;
	}
	
	
	
	
	/**
	 * Reads EggNOG data (version 5). 
	 * 
	 * @param members_reader
	 * @param annotations_reader
	 * @param trees_reader
	 * @param gene_trees trees for families with exactly 1 copy at each leaf
	 * @return
	 * @throws IOException
	 */
    public static AnnotatedTable eggNOG(
    		Reader members_reader, 
    		Reader annotations_reader,
    		Reader trees_reader,
    		List<DataFile<Phylogeny>> gene_trees
    		) throws IOException
    {
    	BufferedReader R = new BufferedReader(members_reader);
        final Map<String,Integer> terminal_taxa = new HashMap<>();
        int num_terminals=0;
    	List<int[]> table_rows = new ArrayList<>();
    	List<String> row_names = new ArrayList<>();
    	Map<String, Integer> family_indexes = new HashMap<>();
    	
    	char geneid_separator = '.';
    	String the_taxon_id=null;
    	
    	String line;
    	do
    	{
    		line = R.readLine();
            if (line == null || line.trim().length()==0 || line.startsWith("#"))
                continue;
            
            String[] fields=line.split("\\t");
    		
            int field_idx = 0;
            String taxon_id = fields[field_idx++];
            if (the_taxon_id==null)
            	the_taxon_id = taxon_id;
            String family_name = fields[field_idx++];
            int num_members = Integer.parseInt(fields[field_idx++]);
            int num_lineages = Integer.parseInt(fields[field_idx++]);
            String[] members = fields[field_idx++].split(",");
            assert (members.length==num_members);

            int[] lineages = new int[members.length];
            for (int m=0; m<members.length; m++)
            {
            	String mem = members[m];
            	int split_at = mem.indexOf(geneid_separator);
            	String otu = mem.substring(0,split_at);
            	int lineage_id;
            	if (terminal_taxa.containsKey(otu))
            	{
            		lineage_id = terminal_taxa.get(otu);
            	} else
            	{
            		lineage_id = num_terminals++;
            		terminal_taxa.put(otu, lineage_id);
            	}
            	lineages[m]=lineage_id;
            }
            int[] row = new int[num_terminals];
            int n=0;
            for (int m=0; m<members.length; m++)
            {
            	int leaf = lineages[m];
            	if (row[leaf]==0) n++;
            	row[leaf]++;
            }
            assert (num_lineages==n);
            table_rows.add(row);
            family_indexes.put(family_name, row_names.size());
            row_names.add(family_name);
    	} while (line != null);
    	
    	R.close();
    	
    	// gather the different lineages encountered 
    	String[] taxon_names = new String[num_terminals];
    	for (String name: terminal_taxa.keySet())
    	{
    		taxon_names[terminal_taxa.get(name)]=name;
    	}
    	// create the profile table set; pad the short rows 
    	int[][] table_profiles = new int[table_rows.size()][];
    	for (int f=0; f<table_profiles.length; f++)
    	{
    		table_profiles[f] = Arrays.copyOf(table_rows.get(f), num_terminals);
    	}
    	AnnotatedTable table = new AnnotatedTable(taxon_names);

    	table.setTable(table_profiles, row_names.toArray(new String[0]));
    	
    	if (annotations_reader!=null)
    	{
    		int cat_col = table.registerProperty("Category");
    		int desc_col = table.registerProperty("Description");
    		
    		R = new BufferedReader(annotations_reader);
    		int row = 0;
    		do
    		{
    			line = R.readLine();
                if (line == null || line.trim().length()==0 || line.startsWith("#"))
                    continue;
                
                String[] fields=line.split("\\t");
        		
                int field_idx = 0;
                String taxon_id = fields[field_idx++];
                String family_name = fields[field_idx++];
                String category = fields[field_idx++];
                String description = field_idx<fields.length?fields[field_idx++]:"";

                assert (table.getFamilyName(row).equals(family_name));
                table.setFamilyProperty(row, cat_col, category);
                table.setFamilyProperty(row, desc_col, description);;
                
                row++;
    		} while (line != null);
    		
    		R.close();
    	}

    	if (trees_reader != null)
    	{
        	String root_id=null;
    		int method_col = table.registerProperty("Method");
    		int trees_col = table.registerProperty("Trees");
    		R = new BufferedReader(trees_reader);
    		int row = 0;
    		Map<Integer, DataFile<Phylogeny>> family_trees = new HashMap<>();
    		do
    		{
    			line = R.readLine();
                if (line == null || line.trim().length()==0 || line.startsWith("#"))
                    continue;
                
                String[] fields=line.split("\\t");
        		
                int field_idx = 0;
                String taxon_id = fields[field_idx++];
                if (root_id==null) root_id=taxon_id;
                String family_name = fields[field_idx++];
    			
                String method = fields[field_idx++];
                String tree = fields[field_idx++];
                
                if (!table.getFamilyName(row).equals(family_name))
                {
                	System.out.println("#**TP.eN different family for trees in row "+row+": have "+table.getFamilyName(row)
                		+", read "+family_name+"\t// resynchronizing");
                	row = family_indexes.get(family_name);
                }
                
                assert (table.getFamilyName(row).equals(family_name));
                if (table.getMemberCount(row)==table.getLineageCount(row)
                		&& table.getLineageCount(row)==terminal_taxa.size())
                {
                	Phylogeny phylo = NewickParser.readTree(new StringReader(tree));
                	int num_nodes = phylo.getNumNodes();
                	int num_fixed_edges = phylo.fixZeroEdges();
                	for (int leaf=0; leaf<phylo.getNumLeaves(); leaf++)
                	{
                		String gene = phylo.getName(leaf);
                		int split_at = gene.indexOf(geneid_separator);
                		String lineage_name = gene.substring(0,split_at);
                		assert (terminal_taxa.containsKey(lineage_name));
                		phylo.getNode(leaf).setName(lineage_name);
                	}
                	File gene_file = new File(root_id, table.getFamilyName(row));
					DataFile<Phylogeny> phylo_data = new DataFile<>(phylo, gene_file);
					if (num_fixed_edges>0)
					{
						phylo_data.setNote(Integer.toString(num_fixed_edges)
								+"("+
								Integer.toString(num_nodes)+"->"+phylo.getNumNodes()
								+")");
					}
					phylo_data.setDirty(true);	
                	family_trees.put(row, phylo_data);
	                table.setFamilyProperty(row, method_col, method);
	                table.setFamilyProperty(row, trees_col, NewickParser.printTree(phylo));
                } else
                {
                	table.setFamilyProperty(row, method_col, "");
                	table.setFamilyProperty(row,  trees_col, "");
                }
    			row++;
    		} while (line != null);
    		R.close();

			System.out.println("#**TP.eN family trees "+family_trees.size());

    		if (family_trees.isEmpty())
    		{
    			Phylogeny star = Phylogeny.starTree(taxon_names);
    			DataFile<Phylogeny> star_data = 
    			new DataFile<>(star, new File((File)null, the_taxon_id));
    			star_data.setNote("star tree because there are no usable gene trees");
    			family_trees.put(-1,star_data);
    			System.out.println("#**TP.eN added star tree");
    		}
    		
    		if (!family_trees.isEmpty())
    		{
    			List<DataFile<Phylogeny>> uniq_trees = new ArrayList<>();
	    		List<TreeComparator> comparators = new ArrayList<>();
	    		for (int f: family_trees.keySet())
	    		{
	    			DataFile<Phylogeny> phylo_data = family_trees.get(f);
	    			Phylogeny phylo = phylo_data.getContent();
	    			boolean uniq=true;
	    			for (TreeComparator cmp: comparators)
	    			{
	    				TreeComparator.NodeMap map = cmp.map(phylo);
	    				int[] mapped = map.fromReference();
	    				boolean diff=false;
	    				for (int node=0; node<mapped.length && !diff; node++)
	    					diff = mapped[node]==-1; // loop stops when an nmapped node is found
	    				uniq = diff;
	    				if (!uniq) break;
	    			}
    				if (uniq)
    				{
    					comparators.add(new TreeComparator(phylo));
    					uniq_trees.add(phylo_data);
    					if (gene_trees != null)
    						gene_trees.add(phylo_data);
    				} else
    				{
    					System.out.println("#**TP.eNOG tree skipped for "+f);
    				}
	    		}
	    		Phylogeny main_phylo = uniq_trees.get(0).getContent();
	    		int[] column_order = new int[terminal_taxa.size()];
	    		for (int leaf=0; leaf<main_phylo.getNumLeaves(); leaf++)
	    		{
	    			String name = main_phylo.getName(leaf);
	    			int col = terminal_taxa.get(name);
	    			column_order[col]=leaf;
	    		}
	        	for (int f=0; f<table_profiles.length; f++)
	        	{
	        		int[] old_row = table_profiles[f];
	        		int[] new_row = new int[old_row.length];
	        		for (int col=0; col<column_order.length; col++)
	        		{
	        			new_row[column_order[col]] = old_row[col];
	        		}
	        		table_profiles[f] = new_row;
	        	}
	        	table.mapTo(main_phylo);
    		}
    		
    	}
    	
    	return table;
    }
    
	public static void printFormattedTable(PrintStream out, AnnotatedTable table, boolean include_properties)
	{
		StringBuilder sb = appendFormattedTableRow(null, -1, table, include_properties); // header
		out.println(sb.toString());
    	for (int i=0; i<table.getFamilyCount(); i++)
    	{
    		sb.setLength(0); // clear
    		sb = appendFormattedTableRow(sb, i, table, include_properties);
    		out.println(sb.toString());
    	}
	}

    private static StringBuilder appendFormattedTableRow(StringBuilder sb, int row, AnnotatedTable table, boolean include_properties)
    {
    	if (sb==null) sb = new StringBuilder();
    	if (row==-1)
    	{
            sb.append(AnnotatedTable.FAMILY_NAME_PROPERTY);
            if (include_properties)
            {
                for (int prop_idx=1; prop_idx<table.getKnownPropertiesCount(); prop_idx++)
                {
                    sb.append('\t');
                    sb.append(table.getPropertyName(prop_idx));
                }
            }
            String[] taxon_names = table.getTaxonNames();
            
            for (int j=0; j<taxon_names.length; j++)
            {
                sb.append("\t");
                sb.append(taxon_names[j]);
            }
    	} else
    	{
    		assert (0 <= row);
    		int i = row;
            int num_properties = table.getKnownPropertiesCount();
            sb.append(table.getFamilyName(i));
            if (include_properties)
            {
                for (int prop_idx=1; prop_idx<num_properties; prop_idx++)
                {
                    sb.append('\t');
                    sb.append(table.getFamilyProperty(i,prop_idx));
                }
            }
            final int[] profile = table.getFamilyProfile(i);
            
            for (int j=0; j<profile.length; j++)
            {
                sb.append("\t");
                if (profile[j]<0)
                    sb.append(MISSING_ENTRY);
                else
                    sb.append(Integer.toString(profile[j]));
            }
    	}
    	return sb;
    }
	
    /**
     * Calculates a string representation of the table data.
     * 
     * @param include_properties
     * @return 
     */
    public static String getFormattedTable(AnnotatedTable table, boolean include_properties)
    {

//        StringBuilder sb = new StringBuilder();
//        sb.append(AnnotatedTable.FAMILY_NAME_PROPERTY);
//        if (include_properties)
//        {
//            for (int prop_idx=1; prop_idx<table.getKnownPropertiesCount(); prop_idx++)
//            {
//                sb.append('\t');
//                sb.append(table.getPropertyName(prop_idx));
//            }
//        }
//        String[] taxon_names = table.getTaxonNames();
//        for (int j=0; j<taxon_names.length; j++)
//        {
//            sb.append("\t");
//            sb.append(taxon_names[j]);
//        }
//        sb.append("\n");
//        int num_properties = table.getKnownPropertiesCount();
    	
    	StringBuilder sb = appendFormattedTableRow(null, -1, table, include_properties); // header
    	for (int i=0; i<table.getFamilyCount(); i++)
    	{
    		
    		sb = appendFormattedTableRow(sb.append("\n"), i, table, include_properties);
    	}
    	sb.append("\n");
//    	
//        for (int i=0; i<table.getFamilyCount(); i++)
//        {
//            sb.append(table.getFamilyName(i));
//            if (include_properties)
//            {
//                for (int prop_idx=1; prop_idx<num_properties; prop_idx++)
//                {
//                    sb.append('\t');
//                    sb.append(table.getFamilyProperty(i,prop_idx));
//                }
//            }
//            final int[] profile = table.getFamilyProfile(i);
//            
//            for (int j=0; j<taxon_names.length; j++)
//            {
//                sb.append("\t");
//                if (profile[j]<0)
//                    sb.append(MISSING_ENTRY);
//                else
//                    sb.append(Integer.toString(profile[j]));
//            }
//            sb.append("\n");
//        }
        return sb.toString();
    }
    
//    public static String getTransposedTable(AnnotatedTable annotated_table)
//    {
//        StringBuilder sb = new StringBuilder();
//        String[] terminal_taxa = annotated_table.getTaxonNames();
//        int[][] table =  new int[annotated_table.getFamilyCount()][];
//        for (int f=0; f<table.length; f++)
//        	table[f] = annotated_table.getFamilyProfile(f);
//        
//        for (int j=0; j<terminal_taxa.length; j++)
//        {
//            sb.append(terminal_taxa[j]);
//            sb.append('\t');
//            for (int i=0; i<table.length; i++)
//            {
//                int n = table[i][j];
//                if (n>9)
//                {
//                    sb.append('(');
//                    sb.append(n);
//                    sb.append(')');
//                } else if (n>=0)
//                {
//                    sb.append(n);
//                }
//            }
//            sb.append('\n');
//        }
//        for (int i=0; i<table.length; i++)
//        {
//            sb.append("#MAP\t");
//            sb.append(i);
//            sb.append('\t');
//            sb.append(annotated_table.getFamilyName(i));
//            sb.append("\t0\n"); // phase in Malin input files 
//        }
//        return sb.toString();
//    	
//    }
    
    private void testEggNOGMembers(String[] args) throws Exception
    {
        int arg_idx = 0;
        if (1+arg_idx > args.length)
        {
            throw new IllegalArgumentException("Call as java "+getClass().getName()+" membes.tsv[.gz]");
        }
    	String members_file = args[arg_idx++];
    	Reader members_reader = GeneralizedFileReader.guessReaderForInput(members_file);
    	Reader annotations_reader = null;
    	if (arg_idx<args.length)
    	{
    		String annotations_file = args[arg_idx++];
    		if (!".".equals(annotations_file))
    			annotations_reader = GeneralizedFileReader.guessReaderForInput(annotations_file);
    	}
    	Reader trees_reader=null;
    	if (arg_idx<args.length)
    	{
    		String trees_file = args[arg_idx++];
    		if (!".".equals(trees_file))
    			trees_reader = GeneralizedFileReader.guessReaderForInput(trees_file);
    	}
    	List<DataFile<Phylogeny>> uniq_trees = new ArrayList<>();
    	AnnotatedTable table = TableParser.eggNOG(members_reader,annotations_reader, trees_reader, uniq_trees);
    	
        PrintStream out = System.out;
        
        out.println(CommandLine.getStandardHeader(this.getClass()));
        out.println(CommandLine.getStandardRuntimeInfo());
        out.println(CommandLine.getStandardHeader("Members file:"+members_file));
        out.println(getFormattedTable(table, true));
        for (int t=0; t<uniq_trees.size(); t++)
        {
        	Phylogeny phylo = uniq_trees.get(t).getContent();
        	out.println("# Tree "+(1+t)+"\t"+uniq_trees.get(t).getFile()+"\t"+NewickParser.printTree(phylo));
        }
        
    }
	
    /**
     * Test code --- reads a phylogeny and a table, and then writes them to stdout.
     */
    private void testTableAndTree(String[] args) throws Exception
    {
        int arg_idx = 0;
        if (2+arg_idx != args.length)
        {
            
            throw new IllegalArgumentException("Call as java "+getClass().getName()+" tree table");
        }
        
        String tree_file = args[arg_idx++];
        String table_file = args[arg_idx++];
        
        Phylogeny tree = NewickParser.readTree(new java.io.FileReader(tree_file));
        AnnotatedTable table = TableParser.readTable(tree.getLeafNames(), 
        		GeneralizedFileReader.guessReaderForInput(table_file), true);
        
        PrintStream out = System.out;
        
        out.println(CommandLine.getStandardHeader(this.getClass()));
        out.println(CommandLine.getStandardRuntimeInfo());
        out.println(CommandLine.getStandardHeader("Tree file: "+tree_file));
        out.println(CommandLine.getStandardHeader("Table file:"+table_file));
        
        out.println(getFormattedTable(table, true));
     }

//    /**
//     * Test code --- reads a phylogeny and a table, and then writes them to stdout.
//     * @param args command line arguments
//     */
//    public static void main(String[] args) throws Exception
//    {
//    	TableParser T = new TableParser();
////        T.testTableAndTree(args);
//        T.testEggNOGMembers(args);
//    }

	public static void main(String[] args) throws Exception
	{
		Class<?> our_class = java.lang.invoke.MethodHandles.lookup().lookupClass();
		count.io.CommandLine cli = new count.io.CommandLine(args, our_class, 1);
		
		PrintStream out = System.out;
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(our_class));
    	    out.println(CommandLine.getStandardRuntimeInfo(our_class, args));
    	}

        Phylogeny tree = cli.getTree();
    	
		if (out != System.out)
		{
			if (out.checkError())
				throw new java.io.IOException("Write failed.");
			out.close();
		}
		
	}
    
    
}
