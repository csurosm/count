package count.ds;
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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A table of arbitrary reconstructions with integer 
 * or floating-point precision.
 * 
 * @param <T> type of ancestral copy statistics 
 */
public class Reconciliation<T extends Number> {
	public Reconciliation(IndexedTree tree, Class<T> type) {
		this.tree = tree;
		if (!Integer.class.equals(type) && !Double.class.equals(type)) {
			throw new IllegalArgumentException("typem must be Integer or Double.");
		}
		this.type = type;
		this.nodeIndex = new HashMap<>();
		this.families = new ArrayList<>();
		this.propIndex = new HashMap<>(); 
		initDataStructures();
	}
	
	private void initDataStructures() {
		for (int u=0; u<tree.getNumNodes(); u++) {
			String uname = tree.getName(u);
			if (uname != null)	{
				uname = uname.replace('_', ' ');
				nodeIndex.put(uname, u);
				System.out.println("#**Rec.iDS\tNODE\t"+u+"\t"+uname+"\t// "+tree.toString(u));
			}
		}
	}
	
	private final IndexedTree tree;
	private Class<T> type;
	private final Map<String, Integer> nodeIndex;
	private final List<HistoryProfile> families;
	
	private final Map<String, Integer> propIndex;
	
	
	/*
	 * Getters
	 */
	public IndexedTree getTree() {return tree;}
	
	public boolean hasIntegers() { return Integer.class.equals(type);}
	
	public void registerProperty(String name) {
		if (!propIndex.containsKey(name)) {
			propIndex.put(name,  propIndex.size());
		}
	}
	
	public int getFamilyCount() {
		return families.size();
	}
	public HistoryProfile getFamily(int fam) {
		return families.get(fam);
	}
	
	/**
	 * Node with this name (all internal nodes must be named).
	 * 
	 * @param name
	 * @return -1 if no such name
	 */
	public int getNode(String name) {
		String uname = name.replace('_', ' ');
		return nodeIndex.containsKey(uname)?nodeIndex.get(uname):-1;
	}
	
	/**
	 * Creates a corresponding profile table with properties and leaf copies.
	 * @return
	 */
	public AnnotatedTable toTable() {
		AnnotatedTable toTable = new AnnotatedTable(tree.getLeafNames());
		
		
		int nfam = families.size();
		int[][] table = new int[nfam][];
		String[] family_names = new String[nfam];
		for (int f=0; f<nfam; f++) {
			HistoryProfile hp = families.get(f);
			table[f] = hp.profile;
			family_names[f] = hp.name;
		}
		toTable.setTable(table, family_names);
		
		String[] props = new String[propIndex.size()];
		
		for (String p: propIndex.keySet())
			props[propIndex.get(p)] = p;
		for (String p: props) {
			toTable.registerProperty(p);
		}
		for (int i=0; i<props.length; i++) {
			String p = props[i];
			for (HistoryProfile hp : families) {
				String value = hp.properties[i];
				if (value != null)
					toTable.setFamilyProperty(hp.fam, p, value);
			}
		}
		
		return toTable;
	}
	
	
	
	/* *****************
	 * HistoryProfiles
	 */
	
	
	/*
	 * Constants for column indices in HistoryProfile tables
	 */
	/**
	 * copy number
	 */
	private static final int HIST_NCOPY = 0;
	/*
	 * DTL event counts 
	 */
	/**
	 * Copy duplications
	 */
	private static final int HIST_DUPLICATION = 1;
	/**
	 * Transfers, including donorless originations
	 */
	private static final int HIST_TRANSFER = 2;
	/**
	 * Copy losses
	 */
	private static final int HIST_LOSS = 3;
	
	/**
	 * Family presence
	 */
	private static final int HIST_P = 4;
	/**
	 * Edge copies / founding inherited copies / vertical inheritance
	 */
	private static final int HIST_SCOPY = 5;
	
	private static final int HIST_STATS = HIST_SCOPY+1; 
	
	/**
	 * Initializing a new family and adding it to the table
	 */
	public HistoryProfile addFamily(String name) {
		HistoryProfile addFamily = this.new HistoryProfile(families.size(), name);
		families.add(addFamily);
		return addFamily;
	}
	
	/**
	 * Storage for family history reconstructions
	 */
	public class HistoryProfile {
		
		private HistoryProfile(int fam, String name) {
			this.name = name;
			int num_leaves = tree.getNumLeaves();
			int num_nodes = tree.getNumNodes();
			this.profile = new int[num_leaves];
			if (Integer.class.equals(type)) {
				dhistory = null;
				ihistory = new int[num_nodes][HIST_STATS];
			} else {
				dhistory = new double[num_nodes][HIST_STATS];
				ihistory = null;
			}
			this.properties = new String[propIndex.size()];
			this.fam = fam;
		}
		/**
		 * Index in {@link Reconciliation#families}
		 */
		private final int fam;
		private final String name;
		private int[] profile;
		/*
		 * one or the other
		 */
		private final double[][] dhistory;
		private final int[][] ihistory;
		
		/**
		 * Family proprties
		 */
		private final String[] properties;
		
		public void setStats(String nodename, double ncopy, double ndup, double ntransfer, double nloss, double presence, double scopy) {
			int node = nodeIndex.get(nodename.replace('_', ' '));
			this.setStats(node, ncopy, ndup, ntransfer, nloss, presence, scopy);
		}
		
		public void setStats(int node, double ncopy, double ndup, double ntransfer, double nloss, double presence, double scopy) {
			assert (Double.class.equals(type));
			if (tree.isLeaf(node)) {
				profile[node] = (int) ncopy;
			}
			dhistory[node][HIST_NCOPY] = ncopy;
			dhistory[node][HIST_DUPLICATION] = ndup;
			dhistory[node][HIST_TRANSFER] = ntransfer;
			dhistory[node][HIST_LOSS] = nloss;
			dhistory[node][HIST_P] = presence;
			dhistory[node][HIST_SCOPY] = scopy;
		}

		public void setStats(String nodename, int ncopy, int ndup, int ntransfer, int nloss, int presence, int scopy) {
			int node = nodeIndex.get(nodename.replace('_', ' '));
			this.setStats(node, ncopy, ndup, ntransfer, nloss, presence, scopy);
		}
		
		public void setStats(int node, int ncopy, int ndup, int ntransfer, int nloss, int presence, int scopy) {
			assert (Integer.class.equals(type));
			if (tree.isLeaf(node)) {
				profile[node] = ncopy;
			}
			ihistory[node][HIST_NCOPY] = ncopy;
			ihistory[node][HIST_DUPLICATION] = ndup;
			ihistory[node][HIST_TRANSFER] = ntransfer;
			ihistory[node][HIST_LOSS] = nloss;
			ihistory[node][HIST_P] = presence;
			ihistory[node][HIST_SCOPY] = scopy;
		}
		
		public void setProperty(String propertyName, String propertyValue) {
			if (!propIndex.containsKey(propertyName))
				throw new IllegalArgumentException("No registered property with name "+propertyName);
			int i = propIndex.get(propertyName);
			properties[i] = propertyValue;
		}
		
		public double dNcopy(int node) {return dhistory[node][HIST_NCOPY];}
		public double dNdup(int node) {return dhistory[node][HIST_DUPLICATION];}
		public double dNtransfer(int node) {return dhistory[node][HIST_TRANSFER];}
		public double dNloss(int node) {return dhistory[node][HIST_LOSS];}
		public double dPresent(int node) {return dhistory[node][HIST_P];}
		public double dScopy(int node) {return dhistory[node][HIST_SCOPY];}
		
		public String getName() { return name;}
		public int getOrder() { return fam;}
		
		
		public String checkDstats() {
			StringBuilder msg = null;
			double diff=0.0;
			for (int v=0; v<tree.getNumNodes(); v++) {
				if (!tree.isRoot(v)) {
					int u = tree.getParent(v);
					double est_ncopy = dNcopy(u)+dNdup(v)+dNtransfer(v)-dNloss(v);
					
					double dc = Math.abs(est_ncopy-dNcopy(v)); 
					if (0.05<dc) {
						if (msg==null) msg=new StringBuilder();
						else msg.append("\n");
						msg.append("Family ").append(fam).append("/").append(name);
						msg.append(", edge ").append(u).append("->").append(v)
						  .append(": copies ").append(dNcopy(u)).append("->").append(dNcopy(v));
						msg.append(", change ").append(dNdup(v)).append("+").append(dNtransfer(v)).append("-").append(dNloss(v));
						msg.append("; gives ").append(est_ncopy);
						msg.append("\t// child ").append(tree.toString(v))	
							.append(", parent ").append(tree.toString(u));
//						+"\t// "++"\t"+tree.toString(u));
					}
					diff += dc;
				}
			}
			if (1.0<diff) {
				if (msg==null) msg=new StringBuilder();
				else msg.append("\n");
				msg.append("Family ").append(fam).append("/").append(name)
					.append(" total difference from expected ").append(diff);
				
			}
			return msg==null?null:msg.toString();
		}
		
	}
}
