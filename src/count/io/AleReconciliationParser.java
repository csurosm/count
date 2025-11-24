package count.io;
/*
 * Copyright 2022 Mikl&oacute;s Cs&#369;r&ouml;s.
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


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import count.ds.IndexedTree;
import count.ds.Phylogeny;

import static count.io.CommandLine.OPT_ANCESTRAL;
import static count.io.CommandLine.OPT_HISTORY;
import static count.io.CommandLine.OPT_OUTPUT;

/**
 * Parses a reconciliation output from ALE 
 */
public class AleReconciliationParser {
	public AleReconciliationParser(IndexedTree tree) {
		this.tree = tree;
		this.max_lineages = tree.getNumNodes();
		this.nodeIndex = new HashMap<>();
		initDataStructures();
	}
	
	private void initDataStructures() {
		for (int u=0; u<tree.getNumNodes(); u++) {
			String uname = tree.getName(u);
			nodeIndex.put(uname, u);
			System.out.println("#**aRP.iDS "+u+"\t"+uname+"\t// "+tree.toString(u));
		}
	}
	
	private final IndexedTree tree;
	private final Map<String, Integer> nodeIndex;
	private final String prefix_families = "#"+OPT_HISTORY.toUpperCase();
	private final String prefix_ancestral =  "#"+OPT_ANCESTRAL.toUpperCase();
	private final String prefix_famstat= "#FAMSTAT";
	private int min_lineages = 0;
	private int max_lineages;
	private int min_copies = 0;
	private int max_copies = Integer.MAX_VALUE;
	private boolean want_families=true;
	
	
	
	private static final int STATS_NDUP = 0;
	private static final int STATS_NTRA = 1;
	private static final int STATS_NLOSS = 2;
	private static final int STATS_NCOPY = 3;
	private static final int STATS_PPRES = 4;
	
	private void setLineageRange(int min, int max) {
		this.min_lineages = min;
		this.max_lineages = max;
	}
	
	
	private void printHeader(PrintStream out) {
		out.print("Family\t");
		for (int leaf = 0; leaf<tree.getNumLeaves(); leaf++) {
			out.printf("\t%s", tree.getName(leaf));
		}		
		out.println();
		
		if (want_families)
			out.println(prefix_families+"\tfamilyidx\tnodeidx\tpresence:p\tignored:m\tignored:p.\tgain:g\tloss:l\tdup:d\torig:+*\tscopies:n\tignored:n.");
		
		out.println(prefix_famstat+"\tnlin\tnmem\tnloss\tngain\tnacquis\tndup\n");
	}
	
	
	private void printTrailer(PrintStream out, double[][] stats) {
		out.println(prefix_ancestral+"\tnode\tpresence:p\tignored:m\tignored:p.\ttransfer:g\tloss:l\tdup:d\toriginations:ori\tscopies:n\tignored:n.");
		int num_nodes = tree.getNumNodes();
		for (int u=0; u<num_nodes; u++) {
			out.printf("%s\t%s", prefix_ancestral, tree.getIdent(u));
			double ppres = stats[u][STATS_PPRES];
			double ndup = stats[u][STATS_NDUP];
			double ntra = stats[u][STATS_NTRA];
			double nloss = stats[u][STATS_NLOSS];
			double ncopy = stats[u][STATS_NCOPY];
			
			double pcopy = tree.isRoot(u)?0.0:stats[tree.getParent(u)][STATS_NCOPY];
			double norig = tree.isRoot(u)?ncopy
					:ncopy-(pcopy+ndup+ntra-nloss);
			
			out.printf("\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" 
						, ppres , 0.0, ppres
						, ntra, nloss, ndup, norig
						, ncopy, ncopy);
			
		}
	}
	
	/**
	 * 
	 * @param out
	 * @param f
	 * @param ale_reconciliation_file
	 * @return success
	 * @throws IOException
	 */
	private boolean printFamily(PrintStream out, int f, File ale_reconciliation_file, double[][] stats) throws IOException {
		BufferedReader R = new BufferedReader(new FileReader(ale_reconciliation_file));

		String family_name = ale_reconciliation_file.getName();
		
		double fam_drate = -1.0;
		double fam_trate = -1.0;
		double fam_lrate = -1.0;
		
		int num_nodes = tree.getNumNodes();
		
		double[] ndup = new double[num_nodes];
		double[] ntra = new double[num_nodes];
		double[] nloss = new double[num_nodes];
		double[] ncopy = new double[num_nodes];
		double[] norig = new double[num_nodes];

		int col_dup=-1;
		int col_transfer = -1;
		int col_loss= -1;
		int col_copy = -1;
		int col_orig = -1;
		
		Arrays.fill(ndup,0.0);
		Arrays.fill(ntra, 0.0);
		Arrays.fill(nloss, 0.0);
		Arrays.fill(ncopy, 0.0);
		Arrays.fill(norig, Double.NaN);
		
		int num_data_lines = 0;
		String line;
		
		int tot_copies= 0;
		int tot_taxa = 0;
		
		boolean fam_rates_follow = false;
		do {
			line = R.readLine();
			if (line != null) {
				String[] fields = line.split("\t");
				if ("rate of".equals(fields[0])) {
					fam_rates_follow = true;
				} else if (fam_rates_follow){ // rate of	 Duplications	Transfers	Losses
					fam_drate = Double.parseDouble(fields[1]);
					fam_trate = Double.parseDouble(fields[2]);
					fam_lrate = Double.parseDouble(fields[3]);
					fam_rates_follow = false;
				} else if ("# of".equals(fields[0])) { // # of	 Duplications	Transfers	Losses	Originations	copies	singletons	extinction_prob
					// header to posteriors table 
					for (int fidx=1; fidx<fields.length; fidx++) {
						fields[fidx]=fields[fidx].trim();
						if ("Duplications".equals(fields[fidx])) col_dup = 1+fidx;
						if ("Transfers".equals(fields[fidx])) col_transfer = 1+fidx;
						if ("Losses".equals(fields[fidx])) col_loss = 1+fidx;
						if ("Originations".equals(fields[fidx])) col_orig = 1+fidx;
						if ("copies".equalsIgnoreCase(fields[fidx])) col_copy = 1+fidx;
						// System.out.println("#**ARP.pF got columns dup "+col_dup+" tra "+col_transfer+" loss "+col_loss+" copy "+col_copy);
					}
				} else if ("S_terminal_branch".equals(fields[0]) || "S_internal_branch".equals(fields[0])) {
					boolean leaf_data = "S_terminal_branch".equals(fields[0]);
					int fidx = 1 ;
					String node_name = fields[fidx++];
					if (leaf_data) {
						int paren_at = node_name.indexOf('(');
						if (0<=paren_at) {
							node_name = node_name.substring(0, paren_at);
						}
					}
					if (!nodeIndex.containsKey(node_name)) {
						System.out.println("#**ARP.pF unknown node "+node_name+" in line "+line+" of file "+ale_reconciliation_file);
						continue;
					}
					
					++num_data_lines;
					int u = nodeIndex.get(node_name);
					ndup[u] = Double.parseDouble(fields[col_dup]);
					ntra[u] = Double.parseDouble(fields[col_transfer]);
					nloss[u] = Double.parseDouble(fields[col_loss]);
					ncopy[u] = Double.parseDouble(fields[col_copy]);
					
					if (col_orig!=-1) {
						norig[u] = Double.parseDouble(fields[col_orig]);
					}
				}
			}
		} while (line != null);
		
		R.close();
		
		boolean got_them = (0 < num_data_lines);

		String profile = "...";

		boolean want_them = min_copies<=tot_copies && tot_copies<=max_copies
					&& min_lineages<=tot_taxa && tot_taxa <= max_lineages;
		
		if (got_them && want_them) {
			if (f==0) printHeader(out);
			
			int ext_at = family_name.indexOf('.');
			if (0<=ext_at) {
				family_name = family_name.substring(0, ext_at);
			}
			out.printf("%s", family_name);
					//"\t%.5f\t%.5f\t%.5f", fam_drate, fam_trate, fam_lrate);
			
			profile = "";
			
			for (int leaf = 0; leaf<tree.getNumLeaves(); leaf++) {
				int n =(int) ncopy[leaf];
				out.printf("\t%d", n);
				
				if (9<n) {
					profile = profile +"("+n+")"; 
				} else {
					profile = profile + Integer.toString(n);
				}
			}
			out.println();
			
			
			double tot_loss=0.0;
			double tot_acquis = 0.0;
			double tot_gain = 0.0;
			double tot_dup= 0.0;
			
			
			for (int node=0; node<num_nodes; node++) {	
				double fsurviving = ncopy[node];
				double fp = Double.min(1.0, fsurviving);
				double ftruecopies = fsurviving;
				double fpc = Double.min(1.0, ftruecopies);
				double fm = 1.<=fsurviving?Double.min(fsurviving-1.0,1.0):0.0; // repurposed for non-singletons (?)
				double fgain=ntra[node];// repurposed for gains
				double floss = nloss[node]; // repurposed for deaths
				double fexpand = ndup[node]; // repurposed for dups 
				
				double forig = norig[node];
				if (Double.isNaN(forig)) {
					double pcopy=0.0;
					if (!tree.isRoot(node)) {
						int parent = tree.getParent(node);
						pcopy = ncopy[parent];
					}
					// ncopy = pcopy-loss + gain + dup + orig
					forig = ncopy[node]-(pcopy-nloss[node])-ntra[node]-ndup[node];
				}
				if (want_families) {
					out.printf("%s\t%d\t%d", prefix_families, f, node);
					
					out.printf("\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n"
							, fp, fm, fpc
							, fgain, floss, fexpand, forig
							, fsurviving, ftruecopies);
				}
				if (stats != null) {
					stats[node][STATS_NDUP]+=ndup[node];
					stats[node][STATS_NTRA]+=ntra[node];
					stats[node][STATS_NLOSS]+=nloss[node];
					stats[node][STATS_NCOPY]+=ncopy[node];
					stats[node][STATS_PPRES]+=fp;
				}
				tot_loss += nloss[node];
				tot_acquis += forig+fgain;
				tot_gain += fgain;
				tot_dup += ndup[node];
			}
			
			
			out.printf("%s\t%d\t%d\t%g\t%g\t%g\t%g\n", prefix_famstat, tot_taxa, tot_copies, tot_loss, tot_gain,tot_acquis,tot_dup);
		}
//		System.out.println("#FAMILYNAME\t"+(f+1)+"\t"+family_name
//				+"\t"+(got_them?"DATA":"IGNORE")
//				+"\t"+profile
//				);
		return got_them && want_them;
	}
	
	/**
	 * Object-linked main with command-line argument parsing
	 * @param cli
	 * @param args
	 * @throws Exception
	 */
	private void main(CommandLine cli, String[] args) throws Exception {
		Class<?> this_class = getClass();
		String directory = cli.getOptionValue("dir");
		String file = cli.getOptionValue("file");
		int n = cli.getOptionInt(CommandLine.OPT_N, 0);
		if (directory==null && file == null) {
			throw new IllegalArgumentException(
					"Call as java -cp ... "
							+this_class
							+" -dir directory tree.newick");
		}
		
		this.min_lineages = cli.getOptionInt(CommandLine.OPT_MINTAXA, this.min_lineages);
		this.max_lineages = cli.getOptionInt(CommandLine.OPT_MAXTAXA, this.max_lineages);
		this.min_copies = cli.getOptionInt(CommandLine.OPT_MINCOPY, this.min_copies);
		this.max_copies = cli.getOptionInt(CommandLine.OPT_MAXCOPY, this.max_copies);
		this.want_families = cli.getOptionBoolean(OPT_HISTORY, this.want_families);

		PrintStream out = System.out;
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	    out.println(CommandLine.getStandardHeader(this_class));
    	    out.println(CommandLine.getStandardRuntimeInfo(this_class, args));
    	}
    	
    	//
    	// Now process the output
    	// 
    	if (directory != null) {
    		if (file != null) {
    			out.println("#** Switch -file "+file+" ignored wiht -dir option");
    		}
    		File d = new File(directory);
    		int f=n;
        	double[][] lineage_stats = new double[tree.getNumNodes()][5];
    		for (File fam_file: d.listFiles()) {
    			if (!fam_file.isDirectory()) {
    				if (printFamily(out, f, fam_file, lineage_stats)) f++;
    			}
    		}
    		printTrailer(out, lineage_stats);
    	} else {
    		printFamily(out, n, new File(file), null);
    	}
    	
    	if (out_file != null)
    		out.close();
	}
	
	public static void main(String[] args) throws Exception {
		Class<?> these = java.lang.invoke.MethodHandles.lookup().lookupClass();
		CommandLine cli = new CommandLine(args, these, 1);
		Phylogeny tree = cli.getTree();
		
		AleReconciliationParser P = new AleReconciliationParser(tree);
		P.main(cli, args);
		
		
		
		
	}
}
