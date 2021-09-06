package count.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;
import java.util.Map;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import count.ds.AnnotatedTable;
import count.ds.Phylogeny;
import count.util.Executable;
public class TableParser 
{
	private TableParser()
	{
		
	}
	
    public static String MISSING_ENTRY = "?";

    /*
     * Reads in an FamilyTable using a reader.
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
        AnnotatedTable table = new AnnotatedTable(taxon_names); // return value

        final Map<String,Integer> terminal_taxa = new HashMap<>();
        for (int i=0; i<taxon_names.length; i++)
        {
            if (terminal_taxa.containsKey(taxon_names[i]))
                throw new IllegalArgumentException("Repeated taxon name "+taxon_names[i]);
            terminal_taxa.put(taxon_names[i], i);
        }
		
        List<String> property_names_in_input = new ArrayList<>();
        List<int[]> parsed_sizes = new ArrayList<>();
        List<String> parsed_family_names = new ArrayList<>();
        List<List<String>> parsed_family_properties = new ArrayList<>();

        BufferedReader R=new BufferedReader(reader);
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
                            property_names_in_input.add(fields[j]);
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

    /**
     * Calculates a string representation of the table data.
     * 
     * @param include_properties
     * @return 
     */
    public String getFormattedTable( AnnotatedTable table, boolean include_properties)
    {
        StringBuilder sb = new StringBuilder();
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
        sb.append("\n");
        int num_properties = table.getKnownPropertiesCount();
        for (int i=0; i<table.getFamilyCount(); i++)
        {
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
            
            for (int j=0; j<taxon_names.length; j++)
            {
                sb.append("\t");
                if (profile[j]<0)
                    sb.append(MISSING_ENTRY);
                else
                    sb.append(Integer.toString(profile[j]));
            }
            sb.append("\n");
        }
        return sb.toString();
    }
    	
	
    /**
     * Test code --- reads a phylogeny and a table, and then writes them to stdout.
     */
    private void mainmain(String[] args) throws Exception
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
        
        out.println(Executable.getStandardHeader(this.getClass()));
        out.println(Executable.getStandardRuntimeInfo());
        out.println(Executable.getStandardHeader("Tree file: "+tree_file));
        out.println(Executable.getStandardHeader("Table file:"+table_file));
        
        out.println(getFormattedTable(table, true));
     }

    /**
     * Test code --- reads a phylogeny and a table, and then writes them to stdout.
     * @param args command line arguments
     */
    public static void main(String[] args) throws Exception
    {
    	TableParser T = new TableParser();
        T.mainmain(args);
    }

}
