package count.model;

import static count.io.CommandLine.OPT_MODEL_ROOT_PRIOR;
import static count.io.CommandLine.OPT_OUTPUT;

import java.io.BufferedReader;
import java.io.PrintStream;
import java.util.Random;

import count.ds.AnnotatedTable;
import count.ds.Phylogeny;
import count.ds.TreeTraversal;
import count.io.CommandLine;
import count.io.GeneralizedFileReader;
import count.io.RateVariationParser;
import count.io.TableParser;

import static count.model.GLDParameters.PARAMETER_LOSS;



/**
 * Ad-hoc class for initializing a model by joining models trained on the root's child trees. 
 * Edge parameters are copied from input 
 * 
 */
public class JoinRatesModels {
	private JoinRatesModels(Phylogeny main_phylo, TreeWithLogisticParameters rates) {
		this.main_phylo = main_phylo;
		this.rates = rates;
	}
	
	private final Phylogeny main_phylo;
	private final TreeWithLogisticParameters rates;
	
	private void setSubtreeRates(int subtree_root, TreeWithLogisticParameters subtree_rates) {
		
		int[] subtree_nodes = TreeTraversal.postOrder(subtree_rates.getTree());
		
		//double child_loss_p = 0.5; 
		
		int[] our_nodes = TreeTraversal.postOrder(main_phylo, subtree_root);
		for (int i=0; i<subtree_nodes.length; i++) {
			int snode = subtree_nodes[i];
			int node = our_nodes[i];
			
			System.out.println("#**JRM.sSR/"+subtree_root+":"+i+"\tfrom "+subtree_rates.getTree().toString(snode)
					+"\tto "+main_phylo.toString(node)
					+"\t// "+subtree_rates.toString(snode));
			
			double logitp = subtree_rates.getLogitLossParameter(snode);
			double logitq = subtree_rates.getLogitDuplicationParameter(snode);
			double log_gamma = subtree_rates.getLogGainParameter(snode, PARAMETER_LOSS);
			if (logitq <= logitp) {
				double logitλ = subtree_rates.getLogitRelativeRate(snode);
				if ( node==subtree_root || logitp == Double.POSITIVE_INFINITY) {
					logitp = rates.getLogitLossParameter(node);
					logitp = 0.0; 
				}
				rates.setLogitLossRelativeDuplication(node, logitp, logitλ, log_gamma, PARAMETER_LOSS);
			} else {
				if ( node==subtree_root || logitp == Double.POSITIVE_INFINITY) {
					logitp = rates.getLogitLossParameter(node);
					logitp = 0.0;
					double logitλ = logitq; 
					rates.setLogitLossRelativeDuplication(node, logitp, logitλ, log_gamma, PARAMETER_LOSS);
				} else 
					rates.setLogitLossDuplication(node, logitp, logitq, log_gamma, PARAMETER_LOSS);
			}
		}
	}
	
	
	
	/**
	 * Main entry : 
	 */
	public static void main(String[] args) throws Exception
	{
		Class<?> our_class = java.lang.invoke.MethodHandles.lookup().lookupClass();
		count.io.CommandLine cli = new count.io.CommandLine(args, our_class, 1);
    	Phylogeny phylo = cli.getTree();
    	
    	int xtra = cli.getExtraArgumentCount();
    	if (xtra != 1 + phylo.getRootNode().getNumChildren()) {
    		throw new IllegalArgumentException("Call as java ... "+our_class+" [-rnd n -root distribution,paramlist] phylo table rates1 rates2 ... ratesk for k="+phylo.getRootNode().getNumChildren()+" root children" );
    	}

		
		
		PrintStream out = System.out;
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(our_class));
    	    out.println(CommandLine.getStandardRuntimeInfo(our_class, args));
    	}
    	
		Random RND = cli.getOptionRND(out);
    	
		count.matek.DiscreteDistribution root_prior = cli.getOptionDistribution(OPT_MODEL_ROOT_PRIOR, null);
		
		TreeWithRates starting_rates = new TreeWithRates(cli.getTree(), RND);
		if (root_prior != null)
		{
			starting_rates.setRootDistribution(root_prior);
			out.println(CommandLine.getStandardHeader("Root prior preset: -"+OPT_MODEL_ROOT_PRIOR+" "+cli.getOptionValue(OPT_MODEL_ROOT_PRIOR)));
		}
		starting_rates.initNodeParameters(starting_rates.getTree().getRoot());
		out.println(CommandLine.getStandardHeader("(Root prior random: "+starting_rates.getRootDistribution()+")"));
		RateVariationModel model = new RateVariationModel(starting_rates);
		model.initConstantRates();
    	JoinRatesModels joint = new JoinRatesModels(phylo, model.getBaseModel());
    	
    	int xargi = 0;
    	AnnotatedTable table;
    	if (cli.getTable()==null) {
    		String table_file =  cli.getExtraArgument(xargi++);
    		if (CommandLine.NO_FILE.equals(table_file)) {
    			table = null;
    		} else {
    			table = TableParser.readTable(phylo.getLeafNames(), 
	    		GeneralizedFileReader.guessReaderForInput(table_file), true);
    		}
    	} else 
    		table = cli.getTable();
    	
    	int child_idx = 0;
    	while (xargi<xtra) {
    		String rates_file = cli.getExtraArgument(xargi++);
    		int child = phylo.getChild(phylo.getRoot(), child_idx++);
    		Phylogeny child_tree = new Phylogeny(phylo, child);
    		BufferedReader B = new BufferedReader(GeneralizedFileReader.guessReaderForInput(rates_file));
    		MixedRateModel M = RateVariationParser.readModel(B, child_tree);
    		assert (M instanceof RateVariationModel);
    		RateVariationModel rvm = (RateVariationModel) M;
    		joint.setSubtreeRates(child, rvm.getBaseModel());
    	}
    	
		out.println(count.io.RateVariationParser.printRates(model));
		if (out_file!=null) out.close();
	}	

}
