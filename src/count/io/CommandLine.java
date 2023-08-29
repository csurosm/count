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
package count.io;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import java.util.Random;

import count.Count;
import count.ds.AnnotatedTable;
import count.ds.Phylogeny;
import count.matek.DiscreteDistribution;
import count.matek.Poisson;
import count.matek.NegativeBinomial;
import count.matek.ShiftedGeometric;
import count.model.FreeMixedModel;
import count.model.GammaInvariant;
import count.model.TreeWithRates;

/**
 * Common initialization of problem instance from the command line. 
 * 
 * @author Miklós Csűrös
 *
 */
public class CommandLine 
{
    public static final String HDR_PREFIX = "#| ";

    public static final String NO_FILE = ".";
    
    private long time_start;
	
	public CommandLine(String[] args, Class<?> for_class) throws IOException
	{
		this(args, for_class, 3);
	}
	
	
	
	/**
	 * 
	 * @param args 3-element array of tree, table and rate files
	 * @param for_class parsing for this class
	 * @throws java.io.IOException if files cannot be read
	 */
	public CommandLine(String[] args, Class<?> for_class, int num_mandatory_arguments) throws IOException
	{
		this.cli_arguments=args;
	    this.cli_options = new Properties();
	    this.cli_for_class = for_class;
	    int arg_idx = 0;
	    
	    while (arg_idx<args.length && args[arg_idx].startsWith("-"))
	    {
	    	String opt;
	    	String val="";
	    	if (args[arg_idx].startsWith("--"))
	    	{
	    		int eq = args[arg_idx].indexOf('=');
	    		if (eq>=0)
	    		{
	    			opt = args[arg_idx].substring(2, eq);
	    			val = args[arg_idx].substring(eq+1);
	    			arg_idx++;
	    		} else
	    		{
	    			opt = args[arg_idx].substring(2);
	    			arg_idx++;
	    			if (arg_idx<args.length && !args[arg_idx].startsWith("-"))
	    				val = args[arg_idx++];
	    		}
	    	} else
	    	{
	    		opt = args[arg_idx].substring(1);
    			arg_idx++;
    			if (arg_idx<args.length && !args[arg_idx].startsWith("-"))
    				val = args[arg_idx++];
	    	}
	    	cli_options.setProperty(opt, val);
//	    	System.out.println("#**CL() got '"+opt+"':'"+val+"'");
	    }
	    
	    
	    if (num_mandatory_arguments+arg_idx > args.length)
	        throw new IllegalArgumentException("Call as java "+for_class.getCanonicalName()+" tree"
	        			+(num_mandatory_arguments==1?"":" table"
	        				+(num_mandatory_arguments<3?"":" rates")));
	

	    setCommonOptions();
	    
	    
	    String tree_file=null;
	    if (cli_options.containsKey(OPT_TREE))
	    	tree_file = cli_options.getProperty(OPT_TREE);
	    else if (arg_idx<args.length)
		    tree_file = args[arg_idx++];
	    
	    
	    String table_file=null;
	    if (cli_options.containsKey(OPT_TABLE))
	    	table_file = cli_options.getProperty(OPT_TABLE);
	    else if (arg_idx<args.length)
	    	table_file = args[arg_idx++];
	    
	    String rates_file=null;
	    if (cli_options.containsKey(OPT_RATES))
	    	rates_file = cli_options.getProperty(OPT_RATES);
	    else if (arg_idx<args.length)
	    	rates_file = args[arg_idx++];
	    
	    extra_arguments = new ArrayList<>();
	    while (arg_idx<args.length)
	    	extra_arguments.add(args[arg_idx++]);
	    	
		
	    Phylogeny tree;
	    if (tree_file != null && !NO_FILE.equals(tree_file))
	    {
	    	tree = count.io.NewickParser.readTree(new java.io.FileReader(tree_file));
	    	int num_nodes = tree.getNumNodes();
	    	int num_fixed_edges = tree.fixZeroEdges();
        	if (num_fixed_edges>0)
        	{
        		System.out.println("#**CL.gTFA() adjusted "
            			+num_fixed_edges
            			+"\tzero-length edges to small positive length"
        			+"\tnodes/leaves "+num_nodes+"->"+tree.getNumNodes()+"/"+tree.getNumLeaves());
        	}
	    	tree_data = new DataFile<>(tree, new File(tree_file));
	    }
	    else
	    {
	    	tree = null;
	    	tree_data = null;
	    }
	    
	    if (table_file != null && !NO_FILE.equals(table_file))
	    {
	    	
	    	String[] taxon_names = (tree==null)?null:tree.getLeafNames();
		    AnnotatedTable table = TableParser.readTable(taxon_names, 
		    		GeneralizedFileReader.guessReaderForInput(table_file), false);
		    table_data = new DataFile<>(table, new File(table_file));
	    } else
	    {
	    	table_data = null;
	    }
	    if (rates_file != null && !NO_FILE.equals(rates_file) && tree!=null)
	    {
		    BufferedReader B = new BufferedReader(GeneralizedFileReader.guessReaderForInput(rates_file));
		    GammaInvariant input_model = RateVariationParser.readRates(B, tree);
		    FreeMixedModel free_model = RateVariationParser.readFreeRates(B, tree, input_model);
		    B.close();
		    model_data = new DataFile<>(input_model, new File(rates_file));
		    free_model_data = new DataFile<>(free_model, new File(rates_file));
	    } else
	    {
	    	model_data = null;
	    	free_model_data = null;
	    }
	
	    if (num_mandatory_arguments>0)
	    {
		    PrintStream out = System.out;
		    
		    out.println(getStandardHeader(for_class));
		    //out.println(getStandardRuntimeInfo());
		    out.println(getStandardRuntimeInfo(for_class, args));
		    if (tree_file != null)
		    {
		    	out.println(getStandardHeader("Tree file:  "+tree_file));
		    }
	    	if (table_file != null)
	    	{
	    		out.println(getStandardHeader("Table file: "+table_file+
	    				(table_data==null
	    					?""
	    					:"\t(hash "+table_data.getContent().tableHashCode()+")")));
	    	}
    		if (rates_file != null)
    			out.println(getStandardHeader("Rates file: "+rates_file));
	    }
	    
	    this.time_start = System.currentTimeMillis();
	    
	}
	
	public DataFile<Phylogeny> getTreeFromArgument(String tree_file) throws IOException
	{
		DataFile<Phylogeny> getTree = null;
		if (tree_file != null && !NO_FILE.equals(tree_file))
	    {
	    	Phylogeny tree = count.io.NewickParser.readTree(new java.io.FileReader(tree_file));
	    	int num_nodes = tree.getNumNodes();
	    	int num_fixed_edges = tree.fixZeroEdges();
        	if (num_fixed_edges>0)
        	{
        		System.out.println("#**CL.gTFA() adjusted "
        			+num_fixed_edges
        			+"\tzero-length edges to small positive length"
        			+"\tnodes/leaves "+num_nodes+"->"+tree.getNumNodes()+"/"+tree.getNumLeaves());
        	}
        	getTree = new DataFile<>(tree, new File(tree_file));
	    }	
		return getTree;
	}
	
	public DataFile<GammaInvariant> getModelFromArgument(Phylogeny tree, String rates_file) throws IOException
	{
		DataFile<GammaInvariant> getModel;
	    if (rates_file != null && !NO_FILE.equals(rates_file) && tree!=null)
	    {
		    BufferedReader B = new BufferedReader(GeneralizedFileReader.guessReaderForInput(rates_file));
		    GammaInvariant input_model = RateVariationParser.readRates(B, tree);
		    FreeMixedModel free_model = RateVariationParser.readFreeRates(B, tree, input_model);
		    B.close();
		    getModel = new DataFile<>(input_model, new File(rates_file));
	    } else
	    {
	    	getModel = null;
	    }
	    return getModel;
		
	}
	
	
	private final Properties cli_options;
	private final List<String> extra_arguments;
	private final String[] cli_arguments;
	private final Class<?> cli_for_class;

	
	public TreeWithRates getRates(){ return getModel().getBaseModel();}
	public AnnotatedTable getTable() { return DataFile.getContent(table_data);}
	public DataFile<AnnotatedTable> getTableData(){ return table_data;}
	public Phylogeny getTree() { return DataFile.getContent(tree_data);}
	public DataFile<Phylogeny> getTreeData(){ return tree_data;}
	public GammaInvariant getModel() { return DataFile.getContent(model_data);}
	public FreeMixedModel getFreeModel() { return DataFile.getContent(free_model_data);}
	public DataFile<GammaInvariant> getModelData(){ return model_data;}
	
	private final DataFile<Phylogeny> tree_data;
	private final DataFile<AnnotatedTable> table_data;
	private final DataFile<GammaInvariant> model_data;
	
	private final DataFile<FreeMixedModel> free_model_data;
	
	
	public static final String OPT_TABLE = "table";
	public static final String OPT_TREE = "tree";
	public static final String OPT_RATES = "rates";
	
	public static final String OPT_OUTPUT = "o";
	public static final String OPT_OUTPUT_TREE = "otree";
//	public static final String OPT_OUTPUT_SAVEALL = "saveall";
	public static final String OPT_LOAD = "load";
	public static final String OPT_SAVE = "save";
	public static final String OPT_SNAPSHOT = "snapshots";
	
	public static final String OPT_TRUNCATE = "truncate";
	public static final String OPT_THREADS = "threads";
	public static final String OPT_TASKS_UNIT = "unittask"; 
	
	public static final String OPT_ROUNDS = "maxiter";
	public static final String OPT_EPS = "convergence";
	
	public static final String OPT_TOP = "top";
	
	public static final String OPT_MINCOPY = "mincopy";
	
	
	public static final String OPT_RELABEL = "relabel";
	public static final String OPT_FILTER = "filter";
	
	public static final String OPT_SPR = "walk";
	public static final String OPT_NNI = "nni";
	public static final String OPT_BUILD = "build";
	public static final String OPT_REROOT = "reroot";
	public static final String OPT_CONTRACT = "contract";
	public static final String OPT_POLLARD = "pollard";
	public static final String OPT_TOUR = "tour";
	public static final String OPT_SEARCH = "search";
	public static final String OPT_PLACE = "place";
	public static final String OPT_NEAR = "near";
	
	public static final String OPT_BOOTSTRAP = "bootstrap";

	public static final String OPT_PARSIMONY = "parsimony";
	public static final String OPT_PARSIMONY_FIT = OPT_PARSIMONY+".fit";
	
	
    public static final String OPT_GAIN = "gain";
    public static final String OPT_DUPLICATION = "duplication";
    public static final String OPT_LOSS = "loss";
    
	public final static String OPT_ROWCOUNT = "n";
	
	public final static String OPT_RND = "rnd";

	
	public static final String OPT_MODEL_CATEGORIES = "k";
	public static final String OPT_MODEL_GAIN_CATEGORIES = "gain_k";
	public static final String OPT_MODEL_LENGTH_CATEGORIES = "length_k";
	public static final String OPT_MODEL_DUPLICATION_CATEGORIES = "duplication_k";
	public static final String OPT_MODEL_FORBIDDEN_GAIN = "forbidden_gain";
	public static final String OPT_MODEL_FORBIDDEN_DUPLICATION = "forbidden_duplication";
	public static final String OPT_MODEL_ROOT_PRIOR = "root"; 
	
	public static final String OPT_MODEL_UNIFORM_DUPLICATION = "uniform_duplication";
	public static final String OPT_MODEL_UNIFORM_GAIN = "uniform_gain";
	

	/**
	 * 
	 * @param opt optional argument name (without starting dashes) 
	 * @return null if was not specified 
	 */
	public String getOptionValue(String opt)
	{
		return cli_options.getProperty(opt);
	}
	public double getOptionDouble(String opt, double default_value)
	{
		return cli_options.containsKey(opt)?Double.parseDouble(getOptionValue(opt)):default_value;
	}
	public int getOptionInt(String opt, int default_value)
	{
		return cli_options.containsKey(opt)?Integer.parseInt(getOptionValue(opt)):default_value;
	}
	
	public long getOptionLong(String opt, long default_value)
	{
		return cli_options.containsKey(opt)?Long.parseLong(getOptionValue(opt)):default_value;
	}
	
	public boolean getOptionBoolean(String opt, boolean default_value)
	{
		
		String val  = cli_options.getProperty(opt, Boolean.toString(default_value));
		boolean b = Boolean.parseBoolean(val);
		if (!b && !"false".equalsIgnoreCase(val))
			throw new IllegalArgumentException("Argument value must be false or true (case ignored) for -"+opt+" [default value "+default_value+"]");
		return b;
	}
	
	
	public DiscreteDistribution getOptionDistribution(String opt, DiscreteDistribution default_value)
	{
		DiscreteDistribution getOptionDistribution = default_value; // return value
		if (cli_options.containsKey(opt))
		{
			String val = cli_options.getProperty(opt);
	        String[] fields = val.split(",");
	        double[] params = new double[fields.length-1];
	        for (int i=0; i<params.length; i++)
	            params[i] = Double.parseDouble(fields[i+1]);
	        if (Poisson.class.getSimpleName().equals(fields[0]))
	        {
	        	getOptionDistribution = new Poisson(params[0]);
	        } else if (NegativeBinomial.class.getSimpleName().equals(fields[0]))
	        {
        		getOptionDistribution = new NegativeBinomial(params[0], params[1]);
	        } else if (ShiftedGeometric.class.getSimpleName().equals(fields[0]))
	        {
	        	getOptionDistribution = new ShiftedGeometric(params[0], params[1]);
	        } else
	        	throw new java.lang.IllegalArgumentException("Prior distribution -"+opt+" '"+fields[0]+"' is unknown.");
			
		} 
		return getOptionDistribution;
	}
	
	public static String encodeDistributionOption(DiscreteDistribution D)
	{
		StringBuilder sb = new StringBuilder(D.getClass().getSimpleName());
		double[] params = D.getParameters();
		for (int pidx=0; pidx<params.length; pidx++) 
		{
			sb.append(",").append(params[pidx]);
		}
		return sb.toString();
	}
	
	public PrintStream getOutput(String opt, PrintStream default_value) throws java.io.FileNotFoundException
	{
		PrintStream out = default_value;
		String out_file = getOptionValue(opt);
    	if (out_file!=null)
    	{
    		if ("-".equals(out_file))
    		{
    			out = System.out;
    		} else 
    			out = new PrintStream(out_file);
    		if (out != default_value)
    		{
    			out.println(CommandLine.getStandardHeader(cli_for_class));
    			out.println(CommandLine.getStandardRuntimeInfo(cli_for_class, cli_arguments));
    		}
    	}
    	return out;
	}
	
	
	/**
	 * Sets up initialization with the {@link #OPT_RND} parameter}.
	 * 
	 * @return null if no random init
	 */
	public Random getOptionRND(PrintStream out)
	{
		Random RND = null;
		if (getOptionValue(OPT_RND)!=null)
		{
			int rnd_seed = getOptionInt(OPT_RND, 0);
			RND = (rnd_seed==0?new Random():new Random(rnd_seed));
			out.println(getStandardHeader("Random initialization: -"+OPT_RND+" "+rnd_seed));    			
		}
		return RND;
		
	}
	
	public void setOptionValue(String opt, String val)
	{
		cli_options.setProperty(opt, val);
	}
	
	public int getExtraArgumentCount()
	{
		return extra_arguments.size();
	}
	
	public String getExtraArgument(int arg_idx)
	{
		return extra_arguments.get(arg_idx);
	}
	
	public long getMillisSinceStart()
	{
		return System.currentTimeMillis()-time_start;
	}
	
	public int getOptionTruncateAbsolute()
	{
        String opt_truncate = getOptionValue(OPT_TRUNCATE);
        return parseTruncateAbsolute(opt_truncate);
//        int absolute = Integer.MAX_VALUE;
//        if (opt_truncate != null)
//        {
//        	int comma_at = opt_truncate.indexOf(',');
//        	if (comma_at<0)
//        	{
//        		throw new IllegalArgumentException(OPT_TRUNCATE+": give comma-separated values absolute,relative");
//        	} 
//        	else
//        	{
//        		absolute = Integer.parseInt(opt_truncate.substring(0,comma_at));
//        	}
//        }
//        return absolute;
	}
	
	public double getOptionTruncateRelative()
	{
        String opt_truncate = getOptionValue(OPT_TRUNCATE);
        return parseTruncateRelative(opt_truncate);
//        double relative = Double.MAX_VALUE;
//        if (opt_truncate != null)
//        {
//        	int comma_at = opt_truncate.indexOf(',');
//        	if (comma_at<0)
//        		throw new IllegalArgumentException(OPT_TRUNCATE+": give comma-separated values absolute,relative");
//    		relative = Double.parseDouble(opt_truncate.substring(comma_at+1));
//        }
//        return relative;	
	}
	
	public static int parseTruncateAbsolute(String opt_truncate)
	{
        int absolute = Integer.MAX_VALUE;
        if (opt_truncate != null)
        {
        	int comma_at = opt_truncate.indexOf(',');
        	if (comma_at<0)
        	{
        		throw new IllegalArgumentException(OPT_TRUNCATE+": give comma-separated values absolute,relative");
        	} 
        	else
        	{
        		absolute = Integer.parseInt(opt_truncate.substring(0,comma_at));
        	}
        }
        return absolute;
	}
	
	public static double parseTruncateRelative(String opt_truncate)
	{
        double relative = Double.MAX_VALUE;
        if (opt_truncate != null)
        {
        	int comma_at = opt_truncate.indexOf(',');
        	if (comma_at<0)
        		throw new IllegalArgumentException(OPT_TRUNCATE+": give comma-separated values absolute,relative");
    		relative = Double.parseDouble(opt_truncate.substring(comma_at+1));
        }
        return relative;	
	}
	
	
	private void setCommonOptions()
	{
		Count.THREAD_PARALLELISM = Integer.parseInt(cli_options.getProperty(OPT_THREADS,
					Integer.toString(Count.THREAD_PARALLELISM)));
		Count.THREAD_UNIT_TASK = Integer.parseInt(cli_options.getProperty(OPT_TASKS_UNIT,
					Integer.toString(Count.THREAD_UNIT_TASK)));
		
	}
    
    public static String getStandardHeader(Class<?> C)
    {
        return getStandardHeader(Count.getAppFullName()+"::"+C.getName());
    }

    public static String getStandardHeader(String info)
    {
        return HDR_PREFIX+info;
    }
    
    public static String getStandardRuntimeInfo()
    {
        java.util.Date now = java.util.Calendar.getInstance().getTime();
        
        java.util.Properties Props=System.getProperties();

        String system = Props.getProperty("os.name", "[unknown OS]")+" "+Props.getProperty("os.version","[Unknown version]")+" "+Props.getProperty("os.arch","[Unknown architecture]");
        String java = Props.getProperty("java.vm.name","[Unknown VM]")+" "+Props.getProperty("java.vm.version","[Unknown version]")+" ("+Props.getProperty("java.runtime.version","[Unknwon runtime]")+") "+Props.getProperty("java.vm.info","")+", "+Props.getProperty("java.vm.vendor","[Uknown vendor]");
        String cwd = Props.getProperty("user.dir","[unknown]");

        StringBuilder message = new StringBuilder();
        message.append(getStandardHeader("System: ")).append(system);
        message.append("\n").append(getStandardHeader("Java engine:")).append(java);
        message.append("\n").append(getStandardHeader("Date: ")).append(now);
        message.append("\n").append(getStandardHeader("Current directory: ")).append(cwd);
        return message.toString();
    }
    
    public static String getStandardRuntimeInfo(Class<?> C, String[] args)
    {
    	StringBuilder sb = new StringBuilder(getStandardRuntimeInfo());
    	sb.append("\n").append(getStandardHeader("Command-line arguments:"))
    	.append(" ").append(C.getCanonicalName());
    	for (String a: args)
    	{
    		sb.append(" ").append(a);
    	}
    	return sb.toString();
    }
	
	
	
}
