package count.gui;
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

import java.awt.Color;
import java.awt.Component;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JEditorPane;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.ToolTipManager;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import count.io.CommandLine;
import count.io.CountXML;
import count.io.DataFile;
import count.io.ModelBundle;
import count.io.Removable;
import count.io.SavableData;
import count.io.ModelBundle.Entry;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.Phylogeny;
import count.ds.TreeComparator;
import count.ds.TreeTraversal;
import count.gui.kit.BoxIcon;
import count.gui.kit.DiamondIcon;
import count.gui.kit.DiscIcon;
import count.gui.kit.PointIcon;
import count.gui.kit.TableIcon;
import count.gui.kit.TreeIcon;

import count.model.MixedRateModel;

/**
 * GUI component (JTree) for displaying a {@link count.io.ModelBundle}. 
 * Tree {@link Node}s have an associated JComponent to show. 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */
public class BundleTree extends JTree
{
	private static boolean REPORT_COLORS_STDOUT = false;
	
	/**
	 * Instantiation with an empty or non-empty bundle. For a non-empty bundle, 
	 * follow up with {@link #copyBundleStructure()}.
	 * 
	 * @param bundle
	 * @param top_item which children should be displayed at the root ({@link #TREES} or {@link #TABLES}
	 */
	protected BundleTree(ModelBundle bundle, Predicate<ModelBundle.Entry> top_item)
	{
		super();
		this.bundle = bundle;
		this.tree_root =  new Node(bundle.getRoot());
		
		this.top_selector = top_item;
//		this.entry_display = new HashMap<>();
		
		tree_model = new DefaultTreeModel(tree_root);
		setModel(tree_model);

        setRootVisible(false);
        setShowsRootHandles(true);
        setEditable(false);
        getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        
        setCellRenderer(new EntryRenderer());
        ToolTipManager.sharedInstance().registerComponent(this); // no tooltips by default
	}
	

	/**
	 * Completes the instantiation if a non-empty bundle is used. 
	 */
	protected void copyBundleStructure()
	{
        copyStructure(tree_root);
        setComponents(tree_root);
        
//        Node lastLeaf = tree_root;
//        while (!lastLeaf.isLeaf())
//        {
//			lastLeaf = (Node)lastLeaf.getChildAt(lastLeaf.getChildCount()-1);
//        }
//        if (lastLeaf == tree_root)
//        	setSelectedNode(null);
//    	else
//    		setSelectedNode(lastLeaf);
	}
	
//	private Node lastNode(Node subtree_root)
//	{
//		if (subtree_root.isLeaf())
//			return subtree_root;
//		else 
//		{
//			Node last_child = (Node)subtree_root.getChildAt(subtree_root.getChildCount()-1);
//			return lastNode(last_child);
//		}
//	}
	
	
	private final Node tree_root;
	private final ModelBundle bundle;
	private final Predicate<Entry> top_selector;
	private final DefaultTreeModel tree_model;
	
	/**
	 * Selection by {@link Entry#isTreeEntry()}, meaning that
	 * {@link Entry#getTreeData()} is not null.
	 */
	public static final Predicate<Entry> TREES = E->E.isTreeEntry();
	
	/**
	 * Selection by {@link Entry#isTableEntry()}, meaning that
	 * {@link Entry#getTableData()} is not null.
	 */
	public static final Predicate<Entry> DATA = E->E.isTableEntry();
	/**
	 * Selection by {@link Entry#isRatesEntry()}, meaning that
	 * {@link Entry#getRatesData()} is not null.
	 */
	public static final Predicate<Entry> RATES = E->E.isRatesEntry();
	/**
	 * Selection for {@link Entry#isTableEntry()} with a non-null
	 * {@link AnnotatedTable} content ({@link DataFile#getContent()}) 
	 * in {@link Entry#getTableData()}.  
	 * 
	 */
	public static final Predicate<Entry> TABLES = E->(E.isTableEntry() && E.getTableData().getContent()!=null);
	
	public static final Predicate<Entry> ALL = E->true; 
	
//	public Predicate<Entry> topSelector()
//	{
//		return top_selector;
//	}
	
	protected Node getRoot()
	{
		return tree_root;
	}
	
	protected Entry getMainTree() {return bundle.getMainTree();}
	
	protected Phylogeny getMainPhylogeny() 
	{
		Phylogeny getMainPhylogeny; // return value

//		System.out.println("#**BT.gMP ");
		
		if (main_phylo == null)
		{
			Entry first_tree = getMainTree();
			getMainPhylogeny = first_tree==null?null:first_tree.getTreeData().getContent();
			main_phylo = getMainPhylogeny;
		} else
		{
			getMainPhylogeny = main_phylo;
		}
		return getMainPhylogeny;
	}
	
    protected Node findItem(JComponent item)
    {
    	Node findItem = null;
    	Enumeration<TreeNode> node_enum = getRoot().postorderEnumeration();
    	while (node_enum.hasMoreElements() && findItem==null)
    	{
    		Node node = (Node) node_enum.nextElement();
    		if (node.getComponent()==item)
    			findItem = node;
    	}
    	return findItem;
    }
    
    /**
     * Finds a node where the entry has this name for file name. 
     * 
     * @param name
     * @return
     */
    protected Node findTree(String name)
    {
    	Node findTree = null;
    	Enumeration<TreeNode> node_enum = getRoot().postorderEnumeration();
    	while (node_enum.hasMoreElements() && findTree==null)
    	{
    		Node node = (Node) node_enum.nextElement();
    		Entry enode = node.getEntry();
    		if (enode.isTreeEntry())
    		{
    			File f = enode.getTreeData().getFile();
    			if (name.equals(f.getName()))
    				findTree = node;
    		}
    	}
    	return findTree;
    }
    
    protected Node findEntry(String entry_id)
    {
    	Node findEntry = null;
    	Enumeration<TreeNode> node_enum = getRoot().postorderEnumeration();
    	while (node_enum.hasMoreElements() && findEntry==null)
    	{
    		Node node = (Node) node_enum.nextElement();
    		Entry enode = node.getEntry();
    		if (enode.getId().equals(entry_id))
    			findEntry = node;
    	}
    	return findEntry;
    }
	
    /**
     * Hues for the terminal nodes, in the order of the first tree.
     */
    private float[] terminal_hues;
    
    private static final float BRIGHTNESS = 0.8f;
    
    private static final float MIN_HUE = 0f;
    private static final float MAX_HUE = 5f/6f; //keep gap between red and purple
	
    private TreeComparator tree_comparator=null;
    
    
    
    
//    private Map<String,Integer> leaf_name_mapping=null;
    
    private Phylogeny main_phylo = null;
	
    /**
     * Copies the structure of the bundle to the tree.
     * 
     * @param node
     */
	private void copyStructure(Node node)
	{
		Entry enode = node.getEntry();
		for (int ci=0; ci<enode.getNumChildren(); ci++)
		{
			Entry cnode = enode.getChild(ci);
			if (!enode.isRoot() || top_selector.test(cnode))
			{
				Node child = new Node(cnode);
				tree_model.insertNodeInto(child, node, node.getChildCount());
				// DEBUG
//				System.out.println("#**BT.cS "+cnode+"\tchild "+child+"\tenode "+enode);
				
//				node.insert(child, node.getChildCount());
				copyStructure(child);
			}
		}
	}
	
	protected Node getSelectedNode()
	{
		Node getSelectedNode = null; // return value
		TreePath selected = this.getSelectionPath();
		if (selected != null)
			getSelectedNode = (Node) selected.getLastPathComponent();
		return getSelectedNode;
	}
	
    protected boolean removeSelectedNode()
    {
    	Node selected_node = getSelectedNode();
        boolean removeSelectedNode = selected_node != null;
        if (removeSelectedNode && selected_node.getComponent() instanceof Removable)
    	{
    		Removable removable = (Removable) selected_node.getComponent();
    		removeSelectedNode = removable.remove();
    	}
        if (removeSelectedNode)
        {
        	Entry selected_entry = selected_node.getEntry();
        	selected_entry.removeFromParent();
        	
        	clearSelection();
        	tree_model.removeNodeFromParent(selected_node);
        	Node root = getRoot();
        	for (int ci=0; ci<root.getChildCount(); ci++)
        	{
        		Node child = (Node)root.getChildAt(ci);
        		if (top_selector.test(child.getEntry()))
        			setSelectedNode(child);
        	}
        }
        return removeSelectedNode;
    }
	
	
	protected Node getSelectedNode(Predicate<Entry> wanted)
	{
		Node getSelectedNode = null; // return value
		TreePath selected = this.getSelectionPath();
		if (selected != null)
		{
			int step = selected.getPathCount();
			while (0<step && getSelectedNode==null)
			{
				--step;
				Node node = (Node)selected.getPathComponent(step);
				Entry enode = node.getEntry();
				if (wanted.test(enode))
					getSelectedNode = node;
			}
		}
		return getSelectedNode;
	}
	
	protected void setSelectedNode(Node node)
	{
		Node prev = getSelectedNode();
		if (node==null)
			this.clearSelection();
		else
		{
			TreePath view_path = new TreePath(node.getPath());
			this.scrollPathToVisible(view_path);
			if (prev != node)
				this.setSelectionPath(view_path);
		}
		
//		System.out.println("#**BT.setSN "+node+"\tprev "+prev);
	}
	
	
	/**
	 * Recursive function for initializing the GUI components for each 
	 * Entry in the bundle. Used when initialized with a non-empty bundle.
	 * 
	 * @param node subtree root
	 */
	private void setComponents(Node node)
	{
		Entry enode = node.getEntry();
		// DEBUG
//		System.out.println("#**BT.sComp "+enode+"\tnode "+node);
		
		if (enode.isTreeEntry())
		{
			node.setComponent(node.createTreeComponent());
		} else if (enode.isRatesEntry())
		{
			node.setComponent(node.createRatesComponent());
		} else if (enode.isTableEntry())
		{
			DataFile<AnnotatedTable> table_data = enode.getTableData();
			if (table_data.getContent() == null)
			{
				String type = enode.getType();
				Entry dnode = enode.getClosestAncestor(TABLES);
				DataFile<AnnotatedTable> orig_data = dnode.getTableData();
				
				if (ParsimonyView.class.getName().equals(type))
				{
					// numerical parsimony ancestral reconstruction

					String sgain = enode.getAttributeValue(CountXML.ATT_GAIN);
					double pty_gain = Double.parseDouble(sgain);
					String sloss = enode.getAttributeValue(CountXML.ATT_LOSS);
					double pty_loss = Double.parseDouble(sloss);
					String sdup = enode.getAttributeValue(CountXML.ATT_DUPLICATION);
					double pty_dup = Double.parseDouble(sdup);
					String tree_id = enode.getAttributeValue(CountXML.ATT_TREE);
					Entry tnode = bundle.getEntry(tree_id);
					assert (tnode != null);
					assert (tnode.isTreeEntry());
					
					DataFile<Phylogeny> tree_data = tnode.getTreeData();
					Phylogeny selected_phylo = tree_data.getContent();
					
					AnnotatedTable mapped_table = orig_data.getContent().mappedToTree(selected_phylo);
					DataFile<AnnotatedTable> mapped_data = new DataFile<>(mapped_table, orig_data.getFile());
					
					ParsimonyView PV = new ParsimonyView(tree_data, mapped_data);
					PV.setPenalties(pty_gain, pty_loss, pty_dup);
					node.setHistoryComponent(PV);
//					decorateByMainTree(PV.getTreePanel());
//					PV.computeAll();
				} else if (PosteriorsView.class.getName().equals(type))
				{
					// posteriors ancestral reconstruction
					String rates_id = enode.getAttributeValue(CountXML.ATT_RATES);
					Entry rnode = bundle.getEntry(rates_id);
					assert (rnode != null);
					assert (rnode.isRatesEntry());
					DataFile<MixedRateModel> rates_data = rnode.getRatesData();
					Entry tnode = rnode.getClosestAncestor(TREES);
					DataFile<Phylogeny> tree_data = tnode.getTreeData();
					Phylogeny selected_phylo = tree_data.getContent();
					
					AnnotatedTable mapped_table = orig_data.getContent().mappedToTree(selected_phylo);
					DataFile<AnnotatedTable> mapped_data = new DataFile<>(mapped_table, orig_data.getFile());
					
					PosteriorsView PV = new PosteriorsView(rates_data, mapped_data);
					if (enode.getAttributeValue(CountXML.ATT_TRUNCATE)!= null)
					{
						String truncate = enode.getAttributeValue(CountXML.ATT_TRUNCATE);
						int absolute = CommandLine.parseTruncateAbsolute(truncate);
						double relative = CommandLine.parseTruncateRelative(truncate);
						PV.setCalculcationWidthThresholds(absolute, relative);
					}
					node.setHistoryComponent(PV);
//					System.out.println("#**BT.sC node "+node+"\tpv "+PV
//							+"\ttruncate "+PV.getCalculationWidthAbsolute()+","+PV.getCalculationWidthRelative()
//							+"\ttree "+tree_data
//							+"\ttable "+mapped_data
//							+"\trates "+rates_data
//					); 
					
//					PV.computeAll();
				} else if (DolloView.class.getName().equals(type))
				{
					// Dollo reconstruction
					String tree_id = enode.getAttributeValue(CountXML.ATT_TREE);
					Entry tnode = bundle.getEntry(tree_id);
					assert (tnode != null);
					assert (tnode.isTreeEntry());
					
					DataFile<Phylogeny> tree_data = tnode.getTreeData();
					Phylogeny selected_phylo = tree_data.getContent();
					
					assert orig_data.getContent().isBinaryTable();
					
					AnnotatedTable mapped_table = orig_data.getContent().mappedToTree(selected_phylo);
					DataFile<AnnotatedTable> mapped_data = new DataFile<>(mapped_table, orig_data.getFile());
					
					DolloView D = new DolloView(tree_data, mapped_data);
					node.setHistoryComponent(D);
//					decorateByMainTree(D.getTreePanel());
//					D.computeAll();
				} else if (SimulationView.class.getName().equals(type))
				{
					String rates_id = enode.getAttributeValue(CountXML.ATT_RATES);
					Entry rnode = bundle.getEntry(rates_id);
					assert (rnode != null);
					assert (rnode.isRatesEntry());
					DataFile<MixedRateModel> rates_data = rnode.getRatesData();
					
					String srnd = enode.getAttributeValue(CountXML.ATT_RND);
					long rnd_seed = Long.parseLong(srnd);
					String srowcount = enode.getAttributeValue(CountXML.ATT_ROWCOUNT);
					int rowcount = Integer.parseInt(srowcount);
					String smincopy = enode.getAttributeValue(CountXML.ATT_MINCOPY);
					int mincopy = Integer.parseInt(smincopy);
					
					SimulationView S = new SimulationView(rates_data, rnd_seed, rowcount, mincopy); 
					node.setHistoryComponent(S);
//					// DEBUG
//					System.out.println("#**BT.sC node "+node+"\tsimview "+S);
			
					
//					S.computeAll();
					// replace parent as well
					AnnotatedTablePanel Stable = S.asTable();
					DataFile<AnnotatedTable> Stable_data = Stable.getDataFile();
					Node parent = (Node) node.getParent();
					parent.setComponent(Stable);
					assert (parent.getEntry() == dnode);
					dnode.setTableData(Stable_data);
				} else
				{
					// what is this?
				}
				HistoryView H = (HistoryView) node.getComponent();
				String list_control_string = enode.getAttributeValue(CountXML.ATT_CONTROLS);
				if (list_control_string != null) {
					String[] wanted_controls = list_control_string.split(",");
					H.getTreeControl().getTreePanel().setControlsByName(wanted_controls);
				}
				String list_selected = enode.getAttributeValue(CountXML.ATT_SELECTION);
				if (list_selected!=null) {
					int[] selected_rows;
					if ("all".equals(list_selected)) {
						int nFam = H.getTableScroll().getDataTable().getRowCount();
						selected_rows = new int[nFam];
						for (int j=0; j<selected_rows.length; j++)
							selected_rows[j] = j;
					} else {
						String[] selected_rows_str = list_selected.split(",");
						selected_rows = new int[selected_rows_str.length];
						for (int j=0; j<selected_rows.length; j++)
							selected_rows[j]= Integer.parseInt(selected_rows_str[j]);
					}
					H.setSelectedFamilies(selected_rows);
				}
				
				
				//((HistoryView) node.getComponent()).setPhyleticProfileColoring(getLeafColors());
			} else // data not null
			{
				if ("true".equals(enode.getAttributeValue(CountXML.ATT_BINARY)))
				{
					// descends from a table
					Entry tnode = enode.getClosestAncestor(TABLES);
					if (tnode != null && !tnode.getTableData().getContent().isBinaryTable())
					{
						// replace by binary conversion
						AnnotatedTable binary = tnode.getTableData().getContent().binaryTable();
						enode.getTableData().setData(binary);
					}
				}
				node.setComponent(node.createTableComponent());
			}
		} else
		{
			// not any of the three types: this is the root
		}
		setSelectedNode(node);
		
		// DEBUG
//		System.out.println("#**BT.sComp "+enode+"\tnode "+node+"\tcomp "+node.getComponent()
//			+"\tctype "+(node.getComponent()==null?null:node.getComponent().getClass().getName()));
		
		
		
		for (int ci=0; ci<node.getChildCount(); ci++)
		{
			Node child = (Node)node.getChildAt(ci);
			setComponents(child);
		}
	}
	
	/**
	 * Sets the attributes of History items to their current value. 
	 */
	public void setAllItemAttributes()
	{
		setAttributes(getRoot());
	}

	private void setAttributes(Node node)
	{
		Entry enode = node.getEntry();
		if (!enode.isRoot() && !enode.isTableEntry()) return; // none of its descendants are, either
		
		JComponent item = node.getComponent();
		if (item instanceof HistoryView)
		{
			HistoryView Hitem = (HistoryView) item;	
			int[] selected_rows = Hitem.getSelectedFamilies();
			StringBuffer selected_list = new StringBuffer();
			if (selected_rows.length == Hitem.getTableScroll().getDataTable().getRowCount()) {
				selected_list.append("all");
			} else {
				for (int s=0; s<selected_rows.length; s++) {
					if (0<s) selected_list.append(",");
					selected_list.append(selected_rows[s]);
				}
			}
			enode.setAttribute(CountXML.ATT_SELECTION, selected_list.toString());
			
			StringBuffer selected_charts = new StringBuffer();
			List<JCheckBox> Hcontrols = Hitem.getTreePanel().getAllControls();
			
			for (JCheckBox cb: Hcontrols) {
				if (cb.isSelected()) {
					if (selected_charts.length()!=0) selected_charts.append(",");
					selected_charts.append(cb.getName());
				}
			}
			
			
			enode.setAttribute(CountXML.ATT_CONTROLS, selected_charts.toString());
			
			
			if (item instanceof ParsimonyView)
			{
				ParsimonyView Pitem = (ParsimonyView) item;
				enode.setAttribute(CountXML.ATT_GAIN, Double.toString(Pitem.getGainPenalty()));
				enode.setAttribute(CountXML.ATT_LOSS, Double.toString(Pitem.getLossPenalty()));
				enode.setAttribute(CountXML.ATT_DUPLICATION, Double.toString(Pitem.getDuplicationPenalty()));
			} else if (item instanceof DolloView)
			{
				DolloView Ditem = (DolloView) item;
				enode.setAttribute(CountXML.ATT_FIX_ROOT, Boolean.toString(Ditem.presentAtRoot()));
			} else if (item instanceof PosteriorsView)
			{	
				PosteriorsView Pitem = (PosteriorsView) item;
				enode.setAttribute(CountXML.ATT_MINCOPY, Integer.toString(Pitem.getMinimumObserved()));
				double relative = Pitem.getCalculationWidthRelative();
				int absolute = Pitem.getCalculationWidthAbsolute();
				if (!Double.isInfinite(relative))
					enode.setAttribute(CountXML.ATT_TRUNCATE, Integer.toString(absolute)+","+Double.toString(relative));
				
				
			} else if (item instanceof SimulationView)
			{
				SimulationView Sitem = (SimulationView) item;
				enode.setAttribute(CountXML.ATT_RND, Long.toString(Sitem.getRandomSeed()));
				enode.setAttribute(CountXML.ATT_ROWCOUNT, Integer.toString(enode.getParent().getTableData().getContent().getFamilyCount()));
				enode.setAttribute(CountXML.ATT_MINCOPY, Integer.toString(Sitem.getMinimumObserved()));
			}
			
		}
		
		
		
		for (int ci = 0; ci<node.getChildCount(); ci++)
		{
			Node child = (Node)node.getChildAt(ci);
			setAttributes(child);
		}
	}
	
	/**
     * Initializes hues for the terminal nodes, splitting the color char by 
     * given tree topology. 
     * 
     * @param node subtree root
     * @param min_hue min hue for subtree
     * @param max_hue max hue for subtree
     * @param sizes calculation needs this node height array; reset at first call with tree root 
     * @param colors terminal hues; reset at first call with tree root
     */
    private static float[] initHues(Phylogeny tree, int node, float min_hue, float max_hue, int[] sizes, float[] colors)
    {
    	if (tree.isRoot(node)) // first call
    	{
    		colors = new float[tree.getNumLeaves()];
			sizes = TreeTraversal.getSubtreeSizes(tree);// TreeTraversal.getHeights(tree);
    	}
    	if (tree.isLeaf(node))
    	{
    		colors[node] = (min_hue+max_hue)/2.0f;
    	} else 
    	{
    		int nc = tree.getNumChildren(node);
    		float[] child_wt = new float[nc];
    		float tot_wt=0f;
    		int ci =0;
    		while (ci < nc)
    		{
    			int child = tree.getChild(node, ci);
    			float wt= sizes[child]; // weighing by subtree size
    			//wt = 1f; // better separation for clade colors 
    			wt = (float)Math.sqrt(wt+1f); 
    			tot_wt += 
    			child_wt[ci] = wt;
    			ci++;
    		}
    		while (ci>0)
    		{
    			--ci;
    			child_wt[ci] /= tot_wt;
    		}
    		tot_wt = 0f;
    		float minc = min_hue;
    		while (ci<nc)
    		{
    			tot_wt += child_wt[ci];
    			float maxc = min_hue * (1f-tot_wt)+max_hue*tot_wt;
    			int child = tree.getChild(node, ci);
    			colors = initHues(tree, child, minc, maxc, sizes,colors);
    			minc = maxc;
    			ci++;
    		}
    	}
    	return colors;
    }	
    
    private class Zuum extends Zoom<TreePanel> implements Removable, SavableData<Phylogeny>
    {
		Zuum(TreePanel panel, DataFile<Phylogeny> tree_data)
		{
			super(panel);
	        Box control_bar =  this.getControlBar();
	        control_bar.removeAll();
	        control_bar.add(panel.createLayoutChooser());
	        control_bar.add(Box.createHorizontalGlue());
	        control_bar.add(panel.createSaveImageButton()); // save as PNG
	        control_bar.add(getSpinner());
	        this.tree_data = tree_data;
		}
		private final DataFile<Phylogeny> tree_data;
		@Override
		public boolean remove()
		{
			if (getTreePanel().getTreeData().getContent() == main_phylo)
				return false;
			// check if there are any dependent rates / history reconstructions
			
//			Session sesh = Session.getSession(this);
//			return sesh.getRatesBrowser().getSelectedPrimaryItem()==null;
			return true;
		}
		@Override
		public DataFile<Phylogeny> getDataFile() 
		{
			return tree_data;
		}
		@Override
		public void saveData(File f) throws IOException 
		{
			getTreePanel().saveData(f);
		}
    }
    
    
    /**
     * Checks if a new tree has the same set of terminal nodes as the main one.
     * 
     * 
     * @param tree new tree to be added
     * @return false if not the same terminal nodes
     */
    private boolean hasSameLeaves(IndexedTree tree)
    {
        Phylogeny main_tree = getMainPhylogeny();
        Map<String,Integer> original_leaves = new HashMap<>();
        for (int leaf=0; leaf<main_tree.getNumLeaves(); leaf++)
        {
        	String leaf_name = main_tree.getName(leaf);
        	original_leaves.put(leaf_name, leaf);
        }        
        
        Map<String, Integer> mapped_leaves = new HashMap<>(); // index mapping in this tree
        StringBuilder sb_problem = new StringBuilder();
        
        for (int leaf=0; leaf<tree.getNumLeaves(); leaf++)
        {
            String leaf_name = tree.getName(leaf);
            if (mapped_leaves.containsKey(leaf_name))
            {
                sb_problem.append("<li>Leaf ").append(leaf_name).append(" appears more than once in the tree.</li>\n");
            } else
            {
            	mapped_leaves.put(leaf_name, leaf);
            }
        	if (!original_leaves.containsKey(leaf_name))
        	{
                sb_problem.append("<li>Leaf ").append(leaf_name).append(" does not appear in main tree.</li>\n");
            }
        }
        for (int leaf=0; leaf<main_tree.getNumLeaves(); leaf++)
        {
        	String leaf_name = main_tree.getName(leaf);
        	if (!mapped_leaves.containsKey(leaf_name))
        	{
                sb_problem.append("<li>Leaf ").append(leaf_name).append(" is missing in this tree.</li>\n");
        	}
        }
        boolean tree_is_correct = (sb_problem.length()==0);
        
        if (!tree_is_correct)
        {
            StringBuilder page_text = new StringBuilder("<h1>Cannot add this tree</h1>");
            page_text.append("<p><em>(But you can start a new session with it...)</em></p>");
            page_text.append("<ul>").append(sb_problem).append("</ul>");
            JEditorPane problems_pane = new JEditorPane("text/html", page_text.toString());
            problems_pane.setEditable(false);
            problems_pane.setBackground(AppFrame.WARNING_COLOR);
            JScrollPane problems_scroll = new JScrollPane(problems_pane);
            problems_scroll.setMaximumSize(new java.awt.Dimension(500,600));
            problems_scroll.setPreferredSize(problems_scroll.getMaximumSize());
            JOptionPane.showMessageDialog(Session.getApp(BundleTree.this),
                    problems_scroll,
                    "Is your tree file correct?",
                    JOptionPane.WARNING_MESSAGE
                    );
        }
        return tree_is_correct;
    }
    
    
    
    
	private JComponent newTreeDisplay(final DataFile<Phylogeny> tree_data)
	{
		if (tree_comparator == null)
		{
			initDataStructures(tree_data.getContent());
		}
	    final TreePanel panel = coloredTreePanel(tree_data);
	    
//	    class Zuum extends Zoom<TreePanel> implements Removable, SavableData<Phylogeny>
//		{
//			Zuum()
//			{
//				super(panel);
//		        Box control_bar =  this.getControlBar();
//		        control_bar.removeAll();
//		        control_bar.add(panel.createLayoutChooser());
//		        control_bar.add(Box.createHorizontalGlue());
//		        control_bar.add(getSpinner());
//			}
//			
//			@Override
//			public boolean remove()
//			{
//				if (getTreePanel().getTreeData().getContent() == main_phylo)
//					return false;
//				// check if there are any dependent rates / history reconstructions
//				
////				Session sesh = Session.getSession(this);
////				return sesh.getRatesBrowser().getSelectedPrimaryItem()==null;
//				return true;
//			}
//			@Override
//			public DataFile<Phylogeny> getDataFile() 
//			{
//				return tree_data;
//			}
//			@Override
//			public void saveData(File f) throws IOException 
//			{
//				getTreePanel().saveData(f);
//			}
//		}	
	    
	    Zuum addTree = new Zuum(panel, tree_data);
	    return addTree;
	}
	
	private TreePanel coloredTreePanel(DataFile<Phylogeny> tree_data)
	{
		TreePanel panel = new TreePanel(tree_data);
		panel.setNormalFontSize((panel.getNormalFontSize()*7)/6);
		decorateByMainTree(panel);
        return panel;		
	}
	
//	/**
//	 * Returns the main tree {@link #main_phylo} 
//	 * for this bundle that defines the default ordering of leaf names. 
//	 * 
//	 * @return
//	 */
//	protected Phylogeny getMainTree()
//	{
//		Entry root_entry = tree_root.getEntry();
//    	Entry first_tree = root_entry.getChild(0, BundleTree.TREES);
//    	Phylogeny getMainTree = first_tree==null?null:first_tree.getTreeData().getContent();
//    	return getMainTree;
//	}
	
	public static void printHsb(float hue, float saturation, float value) {

	    int h = (int)(hue * 6);
	    float f = hue * 6 - h;
	    float p = value * (1 - saturation);
	    float q = value * (1 - f * saturation);
	    float t = value * (1 - (1 - f) * saturation);

	    float r,g,b;
	    
	    switch (h) {
	      case 0: r=value; g=t; b=p; break;
	    	  //return rgbToString(value, t, p);
	      case 1: r=p; g=value; b=p; break;
	    	  //return rgbToString(q, value, p);
	      case 2: r=p; g=value; b=t; break;
	    	  //return rgbToString(p, value, t);
	      case 3: r=p; g=q; b=value; break;
	    	  //return rgbToString(p, q, value);
	      case 4: r=t; g=p; b=value; break;
	    	  //return rgbToString(t, p, value);
	      case 5: r=value; g=p; b=q; break;
	    	  // return rgbToString(value, p, q);
	      default: r=0f; g=0f; b=0f; break;
	    }
	}	
	
	public Color[] decorateByMainTree(TreePanel panel)
	{
		if (tree_comparator == null)
		{
			initDataStructures(bundle.getMainTree().getTreeData().getContent());
		}
		
        int pt_size = panel.getTreePointSize();
        IndexedTree tree = panel.getTreeData().getContent();
        
        int num_nodes = tree.getNumNodes();
        float[] node_hues = new float[num_nodes];

        final int[] height = TreeTraversal.getHeights(tree);
        float root_height = height[height.length-1]; 
        int num_leaves = tree.getNumLeaves();
        
        TreeComparator.NodeMap node_map = tree_comparator.map(tree);
        int[] original_index = node_map.toReference();
        int node=0;
        
        Color[] decorateByMainTree = new Color[num_leaves];
        while (node < num_leaves)
        {
        	int j = original_index[node];
        	assert (j>=0); // there is a leaf with the same name 
        	
            float hue = node_hues[node] = terminal_hues[j];
            
            Color col = Color.getHSBColor(hue, 1.0f, BRIGHTNESS);
            decorateByMainTree[node] = col;	
            
            
            BoxIcon leaf_icon = new BoxIcon(pt_size, true); // filled
            leaf_icon.setDrawColor(TreePanel.TREE_UNSELECTED_LEAF_COLOR);
            leaf_icon.setFillColor(col);
            
            BoxIcon selected_leaf_icon = new BoxIcon(pt_size, true);
            selected_leaf_icon.setDrawColor(TreePanel.TREE_SELECTED_LEAF_COLOR);
            selected_leaf_icon.setFillColor(col);
            selected_leaf_icon.setCrossing(TreePanel.TREE_SELECTED_LEAF_COLOR);

            TreePanel.DisplayedNode D = panel.getNode(node);
            D.setIcon(true, selected_leaf_icon);
            D.setIcon(false, leaf_icon);
        	
            if (REPORT_COLORS_STDOUT)
	        	System.out.printf("#BT.dBMT.COLOR\t[ %.3f %.3f %.3f]\t%% %s\n"
	        			, col.getRed()/255.0
	        			, col.getGreen()/255.0
	        			, col.getBlue()/255.0
	        			, tree.toString(node)
	        			);
            
        	node++;
        }
        while (node < num_nodes)
        {
        	assert (!tree.isLeaf(node));
        	int num_children=tree.getNumChildren(node);
            float hue = 0.0f;
            
            for (int ci=0; ci<num_children; ci++)
            {
                int child = tree.getChild(node, ci);
                float child_hue = node_hues[child];
                hue += child_hue;
            }
            hue /= num_children;
            node_hues[node] = hue;
            
            float sat = 1.0f-height[node]/root_height;
            Color col = Color.getHSBColor(hue, sat, BRIGHTNESS);
            
            int j = original_index[node];
            PointIcon node_unsel;
            PointIcon node_selected;
            if (j<0)
            {
            	// no mapping
            	node_selected = new DiscIcon(pt_size, true);
            	node_unsel = new DiscIcon(pt_size, true);
            } else
            {
            	node_selected = new DiamondIcon(pt_size, true);
            	node_unsel = new DiamondIcon(pt_size, true);
            	// copy name if not set
            	String tname = tree.getName(node);
            	if (tree instanceof Phylogeny && 
            			(tname==null || "".equals(tname)))
            	{
            		Phylogeny phylo = (Phylogeny) tree;
            		phylo.getNode(node).setName(main_phylo.getName(j));
            	} else
            	{
            		// keep that name
            	}
            }
            node_unsel.setDrawColor(TreePanel.TREE_UNSELECTED_NODE_COLOR);
            node_unsel.setFillColor(col);
            node_selected.setDrawColor(TreePanel.TREE_SELECTED_NODE_COLOR);
            node_selected.setFillColor(col);
            
            TreePanel.DisplayedNode D = panel.getNode(node);
            D.setIcon(true, node_selected);
            D.setIcon(false, node_unsel);
            
            
            if (REPORT_COLORS_STDOUT)
	        	System.out.printf("#BT.dBMT.COLOR\t[ %.3f %.3f %.3f]\t%% %s\n"
	        			, col.getRed()/255.0
	        			, col.getGreen()/255.0
	        			, col.getBlue()/255.0
	        			, tree.toString(node)
	        			);
            
            
            ++node;
        }
        return decorateByMainTree;
	}
		
	
	public Color[] getLeafColors()
	{
		JComponent selected_tree_embedding_panel =  getSelectedNode(BundleTree.TREES).getComponent();
		if (selected_tree_embedding_panel ==null || !Zuum.class.isInstance(selected_tree_embedding_panel))
		{
			
			throw new UnsupportedOperationException("Looking for "+Zuum.class+"; got "	
					+(selected_tree_embedding_panel==null?"null":selected_tree_embedding_panel.getClass().getName()));
		}
		Zuum zuum = (Zuum) selected_tree_embedding_panel;
		TreePanel tree_panel = zuum.getTreePanel();
		int num_leaves = tree_panel.getTreeData().getContent().getNumLeaves();
				
		Color[] leaf_colors = new Color[num_leaves];
		for (int leaf=0; leaf<num_leaves; leaf++)
		{
			leaf_colors[leaf] = tree_panel.getNode(leaf).getIcon(true).getFillColor();
		}
		return leaf_colors;
	}
		
	
	/**
	 * Sets the main phylogeny, and initializes {@link #terminal_hues} and {@link #tree_comparator}.
	 * 
	 * @param main_phylo
	 */
	private void initDataStructures(Phylogeny main_phylo)
	{
		this.main_phylo = main_phylo;
        terminal_hues = initHues(main_phylo, main_phylo.getRoot(), MIN_HUE, MAX_HUE, null, null);
//		leaf_name_mapping = new HashMap<>();
//        for (int leaf_idx=0; leaf_idx<main_phylo.getNumLeaves(); leaf_idx++)
//        	leaf_name_mapping.put(main_phylo.getName(leaf_idx),leaf_idx);
		tree_comparator = new TreeComparator(main_phylo);
	}
	
	protected class Node extends DefaultMutableTreeNode
	{
		private Node(ModelBundle.Entry bundle_node)
		{
			super(bundle_node);
		}
		
		private JComponent component;
		
		protected ModelBundle.Entry getEntry()
		{
			return (ModelBundle.Entry)this.getUserObject();
		}
		
		protected void setComponent(JComponent component)
		{
			this.component = component;
		}
		
		/**
		 * Sets the associated component to be a HistoryView : 
		 * also decorates the tree and launches the computations.
		 * @param history_view
		 */
		protected void setHistoryComponent(HistoryView history_view)
		{
			Color[] leafColors = decorateByMainTree(history_view.getTreePanel());
			history_view.setPhyleticProfileColoring(leafColors);
			setComponent(history_view);
			history_view.repaint();
//			System.out.println("#***BT.N.sHC "+toString()+"\thv "+history_view);
			try
			{
				history_view.computeAll();
			} catch (Throwable t)
			{
				System.out.println("#***BT.N.sHC "+toString()+"\tcaught "+t);
			}
		}
		
		protected JComponent getComponent()
		{
			return this.component;
		}
		
		private RateVariationPanel createRatesComponent()
		{
			ModelBundle.Entry enode = getEntry();
			assert enode.isRatesEntry();
			DataFile<MixedRateModel> rates_data = enode.getRatesData();
			RateVariationPanel rates_panel = new RateVariationPanel(rates_data);
			return rates_panel;
		}
		
		private AnnotatedTablePanel createTableComponent()
		{
			ModelBundle.Entry enode = getEntry();
			DataFile<AnnotatedTable> table_data = enode.getTableData();
			assert (enode.isTableEntry() && table_data != null);
			AnnotatedTablePanel table_panel = new AnnotatedTablePanel(table_data);
			return table_panel;
		}
		
		private JComponent createTreeComponent()
		{
			ModelBundle.Entry enode = getEntry();
			assert enode.isTreeEntry();			
			DataFile<Phylogeny> tree_data = enode.getTreeData();
			JComponent tree_panel = newTreeDisplay(tree_data);
			return tree_panel;
		}
		
		/**
		 * Attempts to add a new tree but cheks if it has the same leaf set.
		 * 
		 * @param tree_data
		 * @return null if invalid tree
		 */
		protected Node addTree(DataFile<Phylogeny> tree_data)
		{
			Node child;
			if ( main_phylo != null //tree_comparator != null
					&& !hasSameLeaves(tree_data.getContent()))
			{
				child = null;
			} else
			{
				Entry enode = getEntry();
				Entry etree = enode.addTree(tree_data);
				child = new Node(etree);
	//			this.insert(child, this.getChildCount());
				tree_model.insertNodeInto(child, this, this.getChildCount());
				child.setComponent(child.createTreeComponent());
			}
			return child;
		}
		
		protected Node addTable(DataFile<AnnotatedTable> table_data, AnnotatedTablePanel table_panel)
		{
			Entry enode = getEntry();
			Entry etable = enode.addTable(table_data);
			Node child = new Node(etable);
//			this.insert(child, this.getChildCount());
			tree_model.insertNodeInto(child, this, this.getChildCount());
			if (table_panel == null)
				table_panel = child.createTableComponent();
			child.setComponent(table_panel);
//			enode.setAttribute(CountXML.ATT_TYPE, table_panel.getClass().getName());
			
			
			return child;
		}
		
		protected Node addRates(DataFile<MixedRateModel> rates_data, RateVariationPanel rates_panel)
		{
			Entry enode = getEntry();
			Entry erates = enode.addRates(rates_data);
			Node child = new Node(erates);
//			this.insert(child,  this.getChildCount());
			tree_model.insertNodeInto(child, this, this.getChildCount());
			if (rates_panel==null)
				rates_panel = child.createRatesComponent();
			child.setComponent(rates_panel);
			return child;
		}
		
		protected Node addHistoryItem(HistoryView comp, Entry related_tree, Entry related_rates)
		{
			DataFile<AnnotatedTable> comp_file = new DataFile<>(null, new File(comp.toString()));
			Entry enode = getEntry();
			Entry eitem = enode.addDataItem(comp_file, comp.getClass().getCanonicalName());
			if (related_tree != null)
			{
				eitem.setAttribute(CountXML.ATT_TREE, related_tree.getId());
			}
			if (related_rates != null)
			{
				eitem.setAttribute(CountXML.ATT_RATES, related_rates.getId());
			}
			Node child = new Node(eitem);
//			this.insert(child, this.getChildCount());
			tree_model.insertNodeInto(child, this, this.getChildCount());
			child.setHistoryComponent(comp);
			return child;
		}

		@Override
		public String toString()
		{
			Entry enode = getEntry();
			if (enode.isRoot()) { return "*";}
			else if (enode.isTreeEntry())
			{
				return DataFile.chopDirectory(enode.getTreeData().getFile().getName());
			} else if (enode.isRatesEntry())
			{
				return DataFile.chopDirectory(enode.getRatesData().getFile().getName());
			} else if (TABLES.test(enode))
			{
				return DataFile.chopDirectory(enode.getTableData().getFile().getName());
			} else if (getComponent()!=null)
			{
				return getComponent().toString();
			} else
				return super.toString();
		}
		
		
		
	}
	
	/**
	 * Our special rendering by ModelBundle.Entry types
	 *
	 */
	private class EntryRenderer extends DefaultTreeCellRenderer
	{
		EntryRenderer()
		{
			super();
			table_icon = new TableIcon(24, false);
			filtered_table_icon = new TableIcon(24, true);
			history_icon = TableIcon.historyIcon();
			tree_icon = new TreeIcon();
			rates_icon = TreeIcon.ratesIcon();
			binary_table_icon = new TableIcon(24, false);
			binary_table_icon.setBinary(true);
			binary_filtered_table_icon = new TableIcon(24, true);
			binary_filtered_table_icon.setBinary(true);
		}
		private final TableIcon table_icon;
		private final TableIcon filtered_table_icon;
		private final TableIcon history_icon;
		private final TreeIcon tree_icon;
		private final TreeIcon rates_icon;
		private final TableIcon binary_table_icon;
		private final TableIcon binary_filtered_table_icon;
		
        @Override
        public Component getTreeCellRendererComponent(
                        JTree tree,
                        Object value,
                        boolean sel,
                        boolean expanded,
                        boolean leaf,
                        int row,
                        boolean hasFocus) 
        {
            super.getTreeCellRendererComponent(
                            tree, value, sel,
                            expanded, leaf, row,
                            hasFocus);
            
            Node node = (Node) value;
            Entry enode = node.getEntry();
            if (enode.isTreeEntry())
            {
            	setIcon(tree_icon);
            	setToolTipText("Phylogeny");
            } else if (enode.isRatesEntry())
            {
            	setIcon(rates_icon);
            	setToolTipText("Probabilistic model");
            } else if (enode.isTableEntry())
            {
            	DataFile<AnnotatedTable> table_data = enode.getTableData();
            	if (table_data.getContent()==null)
            	{
            		setIcon(history_icon);
            		setToolTipText("Ancestral reconstruction");
            	} else 
            	{
            		Entry eparent = enode.getParent();
            		boolean isbinary = table_data.getContent().isBinaryTable();
            		boolean isfiltered = eparent.isTableEntry();
            		if (isbinary && isfiltered)
            		{
            			isfiltered = eparent.getTableData().getContent()==null || eparent.getTableData().getContent().isBinaryTable();
            		}
            		if (isbinary)
            		{
            			setIcon(isfiltered?binary_filtered_table_icon:binary_table_icon);
                		setToolTipText("Presence-absence table");
            		} else
            		{
            			setIcon(isfiltered?filtered_table_icon:table_icon);
                		setToolTipText(isfiltered?"Filtered phylogenetic profiles":"Table of integer-valued phylogenetic profiles");
            		}
            	} 
            }
            Color tree_bg = BundleTree.this.getBackground();
            this.setBackground(tree_bg);
            this.setBackgroundNonSelectionColor(tree_bg);
            this.setBackgroundSelectionColor(tree_bg==Color.WHITE?Color.BLUE:Color.DARK_GRAY);
            return this;
        }
		
	}
}
