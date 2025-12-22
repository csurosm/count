package count.io;
/*
 * Copyright 2025 Mikl&oacute;s Cs&#369;r&ouml;s.
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

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.Phylogeny;
import count.ds.Reconciliation;

import static count.io.CommandLine.OPT_ANCESTRAL;
import static count.io.CommandLine.OPT_HISTORY;
import static count.io.CommandLine.OPT_OUTPUT;

/**
 * Parses reconciliation output 
 */
public class ReconciliationParser {
	public ReconciliationParser(IndexedTree tree) {
		this.tree = tree;
		this.max_lineages = tree.getNumNodes();
		this.nodeIndex = new HashMap<>();
		initDataStructures();
	}
	
	private void initDataStructures() {
		for (int u=0; u<tree.getNumNodes(); u++) {
			String uname = tree.getName(u);
			if (uname != null) {
				uname.replace('_', ' '); 
				nodeIndex.put(uname, u);
			}
			if (want_histories)
				System.out.println("#**RP.iDS NODE\t"+u+"\t"+uname+"\t// "+tree.toString(u));
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
	private static boolean want_histories=true;
	
	
	
	private static final int STATS_NDUP = 0;
	private static final int STATS_NTRA = 1;
	private static final int STATS_NLOSS = 2;
	private static final int STATS_NCOPY = 3;
	private static final int STATS_PPRES = 4;
	
	private void setLineageRange(int min, int max) {
		this.min_lineages = min;
		this.max_lineages = max;
	}
    /**
     * Our own exception type.
     */
    public static class FileFormatException extends IOException
    {
        private FileFormatException(String msg)
        {
            super("Error in input file:\n"+msg);
        }
    }
	
    /**
     * Parses a master table 
     * 
     * @param r
     * @param tree
     * @return
     * @throws IOException
     */
	public static Reconciliation<Double> parseTable(BufferedReader r, IndexedTree tree) throws IOException {
		Reconciliation<Double> table = new Reconciliation<>(tree, Double.class);		
		
		int nlines = -1;
		int c_node= -1;
		int c_fam = -1;
		int c_ncopy = -1;
		int c_ndup = -1;
		int c_ntra = -1;
		int c_nloss = -1;
		int c_single = -1; // inherited copies
		int c_pres = -1;
		int c_orig = -1;
		int c_cog = -1;
		String coglabel = null;
		String[] cogannot = CogAnnotator.getAvailableProperties();
		
		String line;
		Reconciliation<Double>.HistoryProfile current= null;
		do {
			line = r.readLine();
			if (line != null) {
				String[] fields = line.split("\t");
				if (nlines == -1) {
					for (int f=0; f<fields.length; f++) {
						String ff = fields[f].toLowerCase();
						if (ff.startsWith("family")) {
							if (c_fam!=-1)
								System.out.println("#**RP.pMT Family header already seen "+fields[c_fam]+" and now "+fields[f]);
							c_fam = f;
						} else if ("BranchName".equals(fields[f]) 
								|| fields[f].startsWith("species")) {
							if (c_node!=-1)
								System.out.println("#**RP.pMT Node header already seen "+fields[c_node]+" and now "+fields[f]);
							c_node = f;
						} else if ("duplications".equals(ff)) {
							if (c_ndup != -1)
								System.out.println("#**RP.pMT Duplications header already seen "+fields[c_ndup]+" and now "+fields[f]);
							c_ndup = f;
						} else if ("transfers".equals(ff)) {
							if (c_ntra != -1)
								System.out.println("#**RP.pMT Transfers header already seen "+fields[c_ntra]+" and now "+fields[f]);
							c_ntra= f;
						} else if ("losses".equals(ff)) {
							if (c_nloss != -1)
								System.out.println("#**RP.pMT Losses header already seen "+fields[c_nloss]+" and now "+fields[f]);
							c_nloss= f;
						} else if (ff.startsWith("origination")) {
							if (c_orig != -1)
								System.out.println("#**RP.pMT Originations header already seen "+fields[c_orig]+" and now "+fields[f]);
							c_orig= f;
						} else if ("presence".equals(ff)) {
							if (c_pres != -1)
								System.out.println("#**RP.pMT Presence header already seen "+fields[c_pres]+" and now "+fields[f]);
							c_pres= f;
						} else if ("singletons".equals(ff)) {
							if (c_single != -1)
								System.out.println("#**RP.pMT Singletons header already seen "+fields[c_single]+" and now "+fields[f]);
							c_single= f;
						} else if ("copies".equals(ff)) {
							if (c_ncopy != -1)
								System.out.println("#**RP.pMT Copies header already seen "+fields[c_ncopy]+" and now "+fields[f]);
							c_ncopy= f;
						} else if (ff.endsWith("cog")) {
							if (c_cog != -1)
								System.out.println("#**RP.pMT COG header already seen "+fields[c_ncopy]+" and now "+fields[f]);
							c_cog= f;
							coglabel = fields[f];
						}
					}
					
					StringBuilder error = null;
					
					if (c_fam==-1) {
						error = new StringBuilder();
						error.append("Family column not found");
					}
					if (c_node==-1) {
						if (error ==null) error = new StringBuilder();
						else error.append(", ");
						error.append("Node column not found");
					}
					if (c_ncopy==-1) {
						if (error ==null) error = new StringBuilder();
						else error.append(", ");
						error.append("Copies column not found");
					}
					if (want_histories && c_ndup==-1) {
						if (error ==null) error = new StringBuilder();
						else error.append(", ");
						error.append("Duplications column not found");
					}
					if (want_histories && c_ntra==-1) {
						if (error ==null) error = new StringBuilder();
						else error.append(", ");
						error.append("Transfers column not found");
					}
					if (want_histories && c_nloss==-1) {
						if (error ==null) error = new StringBuilder();
						else error.append(", ");
						error.append("Losses column not found");
					}
					if (c_single==-1)
						System.out.println("#**RP.pMT Singleton column not found");
					if (c_pres==-1)
						System.out.println("#**RP.pMT Presence column not found");
					if (c_orig==-1)
						System.out.println("#**RP.pMT Origination column not found");
					if (c_cog == -1) {
						System.out.println("#**RP.pMT COG column not found");
					} else {
						table.registerProperty(coglabel);
						if (cogannot != null) 
							for(String prop: cogannot)
								table.registerProperty(prop);
					}
					
					if (error != null) {
						throw new FileFormatException(error.toString());
					}
					
				} else {
					String name = fields[c_fam];
					String branch = fields[c_node].replace('_', ' ');
					double ncopy = Double.parseDouble(fields[c_ncopy]);
					double ndup = c_ndup==-1?0.0:Double.parseDouble(fields[c_ndup]);
					double ntra = c_ntra==-1?0.0:Double.parseDouble(fields[c_ntra]);
					double nloss = c_nloss==-1?0.0:Double.parseDouble(fields[c_nloss]);
					double norig;
					if (c_orig==-1) norig=0.0;
					else norig = Double.parseDouble(fields[c_orig]);
					double nsingle;
					if (c_single == -1)
						nsingle = ncopy-(ndup+ntra+norig);
					else 
						nsingle = Double.parseDouble(fields[c_single]);
					double present;
					if (c_pres==-1) 
						present = Double.min(1.0, ncopy);
					else 
						present = Double.parseDouble(fields[c_pres]);
					
					if (current == null || !current.getName().equals(name)) {
						if (want_histories && nlines % 10 == 9 && nlines < 10000) {
							String error_message =  current.checkDstats();
							if (error_message!=null) {
								throw new FileFormatException(error_message);
							}
						}
						
						current = table.addFamily(name);
						if (coglabel != null) {
							String og = fields[c_cog];
							current.setProperty(coglabel, og);
							if (cogannot!=null) {
								String[] cogdesc = CogAnnotator.getAnnotation(og);
								if (cogdesc != null) {
									for (int i=0; i<cogdesc.length; i++) {
										current.setProperty(cogannot[i], cogdesc[i]);
									}
								}
							}
						}
					}
					int node = table.getNode(branch);
					if (node<0) {
						// unknown node
						if (want_histories)
							throw new FileFormatException("Node name "+branch+" not recognized");
					} else {
						current.setStats(node, ncopy, ndup, ntra+norig, nloss, present, nsingle);
					}
				}
				nlines++;
			}
		} while (line != null);
		
		return table;
	}
	
	private void printHeader(PrintStream out) {
		out.print("Family\t");
		for (int leaf = 0; leaf<tree.getNumLeaves(); leaf++) {
			out.printf("\t%s", tree.getName(leaf));
		}		
		out.println();
		
		if (want_histories)
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
						// System.out.println("#**RP.pF got columns dup "+col_dup+" tra "+col_transfer+" loss "+col_loss+" copy "+col_copy);
					}
				} else if ("S_terminal_branch".equals(fields[0]) || "S_internal_branch".equals(fields[0])) {
					boolean leaf_data = "S_terminal_branch".equals(fields[0]);
					int fidx = 1 ;
					String node_name = fields[fidx++];
					if (leaf_data) {
						int paren_at = node_name.indexOf('('); // node index in parens 
						if (0<=paren_at) {
							node_name = node_name.substring(0, paren_at);
						}
					}
					if (!nodeIndex.containsKey(node_name)) {
						System.out.println("#**RP.pF unknown node "+node_name+" in line "+line+" of file "+ale_reconciliation_file);
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
					
					
					if (leaf_data) {
						assert (tree.isLeaf(u));
						if (0<ncopy[u])
							tot_taxa++;
						tot_copies += ncopy[u];
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
			
//			int chop_at = family_name.indexOf("FAM");
//			if (0<=chop_at) {
//				// c, c+1, c+2 ...
//				int chop_until = chop_at+3;
//				while (Character.isDigit(family_name.charAt(chop_until)))
//					++chop_until;
//				family_name = family_name.substring(chop_at, chop_until);
//			} else {
				int ext_at = family_name.indexOf('.');
				if (0<=ext_at) {
					family_name = family_name.substring(0, ext_at);
				}
//			}
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
				if (want_histories) {
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
			
			
			out.printf("%s\t%d\t%d\t%g\t%g\t%g\t%g\t%s\n", prefix_famstat, tot_taxa, tot_copies, tot_loss, tot_gain,tot_acquis,tot_dup,family_name);
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
							+this_class.getName()
							+"[-"+OPT_HISTORY+" true|false]"+" [-dir directory|-file file] tree.newick"
						    +"\none of -file and -dir must be specified");
		}
		
		this.min_lineages = cli.getOptionInt(CommandLine.OPT_MINTAXA, this.min_lineages);
		this.max_lineages = cli.getOptionInt(CommandLine.OPT_MAXTAXA, this.max_lineages);
		this.min_copies = cli.getOptionInt(CommandLine.OPT_MINCOPY, this.min_copies);
		this.max_copies = cli.getOptionInt(CommandLine.OPT_MAXCOPY, this.max_copies);
		

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
    	if (want_histories) {
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
    	} else {
    		if (file == null)
    			throw new IllegalArgumentException("Specify -file with -"+OPT_HISTORY+" true");
    		
    		BufferedReader R = GeneralizedFileReader.guessBufferedReaderForInput(file);
    		Reconciliation<Double> rec = parseTable(R, this.tree);
    		
    		AnnotatedTable table = rec.toTable();
    		out.println(TableParser.getFormattedTable(table, true));
    		
    		R.close();
    	}
    	
    	if (out_file != null)
    		out.close();
	}
	
	public static void main(String[] args) throws Exception {
		Class<?> these = java.lang.invoke.MethodHandles.lookup().lookupClass();
		CommandLine cli = new CommandLine(args, these, 1);
		Phylogeny tree = cli.getTree();
		
		want_histories = cli.getOptionBoolean(OPT_HISTORY, want_histories);
		
		ReconciliationParser recP = new ReconciliationParser(tree);
		recP.main(cli, args);
		
	}
}
