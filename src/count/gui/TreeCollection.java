package count.gui;
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

import java.awt.Color;
import java.awt.GridLayout;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.Box;
import javax.swing.JComponent;
import javax.swing.JEditorPane;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.event.ChangeListener;

import count.ds.IndexedTree;
import count.ds.Phylogeny;
import count.ds.TreeComparator;
import count.ds.TreeTraversal;
import count.gui.kit.BoxIcon;
import count.gui.kit.DiamondIcon;
import count.gui.kit.DiscIcon;
import count.gui.kit.PointIcon;
import count.io.DataFile;
import count.io.Removable;
import count.io.SavableData;

/**
*
* Collection of trees: multiple phylogenies with the same terminal node set.
* 
* @author Mikl&oacute;s Cs&#369;r&ouml;s
* @deprecated
*/
public class TreeCollection extends JPanel
{
	private final Browser<Zuum> tree_browser;
//    private final JTabbedPane tree_panels;
    /**
     * Common set of terminal names, and their indexes in the trees.
     */
    private final List<Map<String, Integer>> map_names;
    /**
     * Hues for the terminal nodes, in the order of the first tree.
     */
    private float[] terminal_hues;
    
    private static final float BRIGHTNESS = 0.8f;
    
    private static final float MIN_HUE = 0f;
    private static final float MAX_HUE = 5f/6f; //keep gap between red and purple
    
    private final Phylogeny main_tree;
    private final TreeComparator main_tree_mapper;

    public TreeCollection(DataFile<Phylogeny> main_tree_data)
    {
    	this(main_tree_data,null);
    }
    /**
     * 
     * 
     * @param main_tree_data
     * @param embed_pane this is a standalone Panel if null
     */
    public TreeCollection(DataFile<Phylogeny> main_tree_data, JSplitPane embed_pane)
    {
        super(new GridLayout(1,1));
        if (embed_pane==null)
        {
        	this.tree_browser = new Browser(Zuum.class);
        }
        else
        {
        	this.tree_browser = new Browser<>(Zuum.class, embed_pane);
        }
        this.main_tree = main_tree_data.getContent();
        this.main_tree_mapper = new TreeComparator(main_tree);
        this.map_names = new ArrayList<>();
        initHues(main_tree.getRoot(), MIN_HUE, MAX_HUE, TreeTraversal.getHeights(main_tree));
//        this.tree_panels= new JTabbedPane(JTabbedPane.LEFT);
//        super.add(tree_panels);
        addTopTree(main_tree_data);
        if (embed_pane == null)
        	add(tree_browser);
    }
    
    /** 
     * TreePanel with Savable phylogenies (Newick format). 
     * @author csuros
     *
     */
    private static class PhylogenyPanel extends TreePanel
		implements SavableData<Phylogeny>
    {
		private final DataFile<Phylogeny> phylo_data;
		PhylogenyPanel(DataFile<Phylogeny> phylo_data) 
		{
			super(phylo_data);
			this.phylo_data = phylo_data;
	        setNormalFontSize((super.getNormalFontSize()*7)/6); // 12->14, 11->12.8, 10->11.7
		}
	
		@Override
		public DataFile<Phylogeny> getDataFile() 
		{
			return phylo_data;
		}
		
    }
    
    /**
     * Our class for properly typing the {@link #tree_browser},
     * and for adding the layout chooser. 
     * 
     * @author csuros
     *
     */
    private class Zuum extends Zoom<PhylogenyPanel>
    		implements Removable, SavableData<Phylogeny>
    {
    	Zuum(PhylogenyPanel panel)
    	{
    		super(panel);
            Box control_bar =  this.getControlBar();
            control_bar.removeAll();
            control_bar.add(panel.createLayoutChooser());
            control_bar.add(Box.createHorizontalGlue());
            control_bar.add(getSpinner());
    	}
    	
    	@Override
    	public boolean remove()
    	{
    		if (getTreePanel().getTreeData().getContent() == main_tree)
    			return false;
    		Session sesh = Session.getSession(this);
			return sesh.getRatesBrowser().getSelectedPrimaryItem()==null;
    	}

		@Override
		public void saveData(File f) throws IOException 
		{
			getTreePanel().saveData(f);
		}

		@Override
		public DataFile<Phylogeny> getDataFile() 
		{
			return getTreePanel().getDataFile();
		}
    	
    }

    public TreePanel addTopTree(DataFile<Phylogeny> tree_data)
    {
    	return addTreePanel(new PhylogenyPanel(tree_data), true);
    }
    
    public TreePanel addTree(DataFile<Phylogeny> tree_data)
    {
    	return addTreePanel(new PhylogenyPanel(tree_data), false);
    }
    
    private TreePanel addTreePanel(PhylogenyPanel panel, boolean at_top)
    {
        Map<String, Integer> leaf_name_mapping = mapTerminals(panel);
        if (leaf_name_mapping==null) // bad tree
            return null;

        map_names.add(leaf_name_mapping); 
        if (map_names.size()==1)
        {
        	// first, the main tree is added
        	calculateNodeColors(panel, leaf_name_mapping);
        } else
        {
        	decorateByMainTree(panel);
        }

        Zuum zum = new Zuum(panel);
    	if (at_top)
    	{
    		tree_browser.addTopItem(zum);
    	}
    	else
    	{
    		tree_browser.addItem(zum);
    	}
    	return panel;
    }
    
//    public TreePanel addTree(TreePanel panel)
//    {
//        Map<String, Integer> leaf_name_mapping = mapTerminals(panel);
//        if (leaf_name_mapping==null) // bad tree
//            return null;
//
//        map_names.add(leaf_name_mapping); 
//        calculateNodeColors(panel, leaf_name_mapping);
//                
//
//        Zoom<TreePanel> zoom = new Zoom<>(panel);
//        Box control_bar =  zoom.getControlBar();
//        control_bar.removeAll();
//        control_bar.add(panel.createLayoutChooser());
//        control_bar.add(Box.createHorizontalGlue());
//        control_bar.add(zoom.getSpinner());
//
//        // add zoom and set short name
//        tree_panels.add(DataFile.chopFileExtension(panel.getTreeData().getFile().getName()), zoom);
//        tree_panels.setSelectedIndex(tree_panels.getTabCount()-1);
//        
////        System.out.println("#**TC.aT add tree "+panel.getTreeData().getContent());
//        
//        return panel;
//    }
        
    
    public TreePanel getSelectedTree()
    {
    	Zuum zum = tree_browser.getSelectedPrimaryItem();
    	if (zum==null) return null;
    	else return zum.getTreePanel();
    }
//    	this.getSelectedItem()
//        int active_idx = tree_panels.getSelectedIndex();
//        if (active_idx==-1)
//            return null;
//        else
//        {
//            Zoom<TreePanel> zum = (Zoom<TreePanel>)tree_panels.getComponentAt(active_idx); // warning but <TreePanel> is sure
//            return zum.getTreePanel();
//        }
//    }
    
//    private TreePanel getTree(int idx)
//    {
//        return ((Zoom<TreePanel>) tree_panels.getComponentAt(idx)).getTreePanel(); // warning but we only add Zoomed TreePanel components here 
//    }
    
    public TreePanel getTree(String tree_name)
    {
    	for (Zuum zum: tree_browser.primaryItemList())
    	{
    		PhylogenyPanel panel = zum.getTreePanel();
    		String panel_name = panel.getTreeName();
    		if (tree_name.equals(panel_name))
    			return panel;
    	}
    	return null;
    }
    
    public TreePanel selectTree(IndexedTree tree)
    {
    	TreePanel TP = getSelectedTree();
    	if (TP.getTreeData().getContent() == tree)
    	{
    		return TP;
    	} else
    	{
    		for (Zuum zum: tree_browser.primaryItemList())
        	{
    			TP = zum.getTreePanel();
	    		IndexedTree TPtree = TP.getTreeData().getContent();
	    		if (TPtree == tree) // same object
	    		{
	    			tree_browser.selectItem(zum);
	    			return TP;
	    		}
	    	}
	    	return null;
    	}
    }

    
    
    
    /**
     * Adds a change listener to tab selection.
     * 
     * @param listener
     */
    public void addChangeListener(ChangeListener listener)
    {
    	tree_browser.addSelectionChangeListener(listener);
    }
    
    public JComponent getBrowserComponent()
    {
    	return tree_browser.getBrowserComponent();
    }
    
    /**
     * Names at the leaves, in the order of the first (main) tree.
     * 
     * @return
     */
    public String[] getLeafNames()
    {
    	return main_tree.getLeafNames();
    }
    
    
    public IndexedTree getMainTree()
    {
    	return main_tree;
    }
    
    /**
     * Initializes hues for the terminal nodes, splitting tye color chart by 
     * main tree topology. 
     * 
     * @param node subtree root
     * @param min_hue min hue for subtree
     * @param max_hue max hue for subtree
     * @param heights calculation needs this node height array 
     */
    private void initHues(int node, float min_hue, float max_hue, int[] heights)
    {
    	if (main_tree.isRoot(node)) // first call
    		terminal_hues = new float[main_tree.getNumLeaves()];
    		
    	if (main_tree.isLeaf(node))
    	{
    		terminal_hues[node] = (min_hue+max_hue)/2.0f;
    	} else // weighing by subtree heights
    	{
    		int nc = main_tree.getNumChildren(node);
    		float[] child_wt = new float[nc];
    		float tot_wt=0f;
    		int ci =0;
    		while (ci < nc)
    		{
    			int child = main_tree.getChild(node, ci);
    			tot_wt += 
    			child_wt[ci] = heights[child]+1f;
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
    			int child = main_tree.getChild(node, ci);
    			initHues(child, minc, maxc, heights);
    			minc = maxc;
    			ci++;
    		}
    	}
    }
    
    /**
     * Checks if a new tree has the same set of terminal nodes as the first one,
     * and computes the mapping.
     * 
     * @param panel new tree panel to be added
     * @return null if not the same terminal nodes
     */
    private Map<String, Integer> mapTerminals(TreePanel panel)
    {
        
        IndexedTree tree = panel.getTreeData().getContent();
        
        if (map_names.isEmpty())
        { // first tree to be added
            Map<String, Integer> terminal_names =new HashMap<>();
            for (int leaf_idx=0; leaf_idx<tree.getNumLeaves(); leaf_idx++)
                terminal_names.put(tree.getName(leaf_idx),leaf_idx);
            
            return terminal_names;
        } else
        {
            Map<String, Integer> terminal_names = map_names.get(0); // index mapping in first tree
            Map<String, Integer> mapped_leaves = new HashMap<>(); // index mapping in this tree
            StringBuilder sb_problem = new StringBuilder();
            
            for (int leaf_idx=0; leaf_idx<tree.getNumLeaves(); leaf_idx++)
            {
                String leaf_name = tree.getName(leaf_idx);
                
                if (terminal_names.containsKey(leaf_name))
                {
                    if (mapped_leaves.containsKey(leaf_name))
                    {
                        sb_problem.append("<li>Leaf ").append(leaf_name).append(" appears more than once in the tree.</li>\n");
                    } else
                    {
                        int orig_idx = terminal_names.get(leaf_name);
                        String name_str  = main_tree.getName(orig_idx); // same String reused
                        mapped_leaves.put(name_str, leaf_idx);
                    }
                } else
                {
                    sb_problem.append("<li>Leaf ").append(leaf_name).append(" does not appear in first tree.</li>\n");
                }
            }
            for (String leaf_name: terminal_names.keySet())
            {
                if (!mapped_leaves.containsKey(leaf_name))
                {
                    sb_problem.append("<li>Leaf ").append(leaf_name).append(" from the first tree is missing from this tree.</li>\n");
                }
            }
            
            boolean tree_is_correct = (sb_problem.length()==0);
            
            if (tree_is_correct)
            {
                return mapped_leaves; 
            } else
            {
                StringBuilder page_text = new StringBuilder("<h1>Cannot add this tree</h1>");
                page_text.append("<p><em>(But you can start a new session with it...)</em></p>");
                page_text.append("<ul>").append(sb_problem).append("</ul>");
                JEditorPane problems_pane = new JEditorPane("text/html", page_text.toString());
                problems_pane.setEditable(false);
                problems_pane.setBackground(AppFrame.WARNING_COLOR);
                JScrollPane problems_scroll = new JScrollPane(problems_pane);
                problems_scroll.setMaximumSize(new java.awt.Dimension(500,400));
                problems_scroll.setPreferredSize(problems_scroll.getMaximumSize());
                JOptionPane.showMessageDialog(tree_browser,
                        problems_scroll,
                        "Is your tree file correct?",
                        JOptionPane.WARNING_MESSAGE
                        );
                return null;
            }
        }
    }

    public void decorateByMainTree(TreePanel panel)
    {
        IndexedTree tree = panel.getTreeData().getContent();
        TreeComparator.NodeMap node_map = main_tree_mapper.map(tree);
//        
//    	
//    	
//        Map<String, Integer> leaf_index = mapTerminals(panel);
//        if (leaf_index != null)
//        {
//            this.calculateNodeColors(panel, leaf_index);
//        }
        this.decorateByMainTree(panel, node_map);
    }    
    
    
    private void decorateByMainTree(TreePanel panel, TreeComparator.NodeMap node_map)
    {
        int pt_size = panel.getTreePointSize();

        IndexedTree tree = panel.getTreeData().getContent();
        int num_nodes = tree.getNumNodes();
        float[] node_hues = new float[num_nodes];

        final int[] height = TreeTraversal.getHeights(tree);
        float root_height = height[height.length-1]; 
        int num_leaves = tree.getNumLeaves();
        
        int[] original_index = node_map.toReference();
        int node=0;
        while (node < num_leaves)
        {
        	int j = original_index[node];
        	assert (j>=0); // there is a leaf with the same name 
        	
            float hue = node_hues[node] = terminal_hues[j];
            
            Color col = Color.getHSBColor(hue, 1.0f, BRIGHTNESS);
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
            		phylo.getNode(node).setName(main_tree.getName(j));
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
            
            ++node;
        }
    }
    
    /**
     * Colors the tree nodes by using the first tree as reference.
     * 
     * @param panel tree panel where the colors are set
     * @param leaf_index mapping from leaf names 
     */
    private void calculateNodeColors(TreePanel panel, Map<String, Integer> leaf_index)
    {
        IndexedTree tree = panel.getTreeData().getContent();
        int num_leaves = leaf_index.size();
        int num_nodes = tree.getNumNodes();
        
        float[] node_hues = new float[num_nodes];
        
        final int[] height = TreeTraversal.getHeights(tree);
        float root_height = height[height.length-1]; 

        int pt_size = panel.getTreePointSize();
        Map<String, Integer> original_index = map_names.get(0);
        int node_idx = 0;
        while (node_idx<num_leaves)
        {
            int j = original_index.get(tree.getName(node_idx));
            float hue = node_hues[node_idx] = terminal_hues[j];
            
            Color col = Color.getHSBColor(hue, 1.0f, BRIGHTNESS);
            BoxIcon leaf_icon = new BoxIcon(pt_size, true); // filled
            leaf_icon.setDrawColor(TreePanel.TREE_UNSELECTED_LEAF_COLOR);
            leaf_icon.setFillColor(col);
            
            BoxIcon selected_leaf_icon = new BoxIcon(pt_size, true);
            selected_leaf_icon.setDrawColor(TreePanel.TREE_SELECTED_LEAF_COLOR);
            selected_leaf_icon.setFillColor(col);
            selected_leaf_icon.setCrossing(Color.RED);

            TreePanel.DisplayedNode D = panel.getNode(node_idx);
            D.setIcon(true, selected_leaf_icon);
            D.setIcon(false, leaf_icon);
            node_idx++;
        }
        while (node_idx<tree.getNumNodes())
        {
            float hue = 0.0f;
            for (int ci=0; ci<tree.getNumChildren(node_idx); ci++)
            {
                int child_idx = tree.getChild(node_idx, ci);
                float child_hue = node_hues[child_idx];
                hue += child_hue;
            }
            hue /= tree.getNumChildren(node_idx);
            node_hues[node_idx] = hue;
            
//            System.out.println("#*TC.cNC "+node_idx+"\t"+hue+"\t"+tree.getNode(node_idx));
            
            float sat = 1.0f-height[node_idx]/root_height;
            Color col = Color.getHSBColor(hue, sat, BRIGHTNESS);
//            DiamondIcon node_empty = new DiamondIcon(pt_size, false);
            DiamondIcon node_full = new DiamondIcon(pt_size, true);
            
            node_full.setDrawColor(TreePanel.TREE_UNSELECTED_NODE_COLOR);
            node_full.setFillColor(col);
            
            DiamondIcon selected_node_icon = new DiamondIcon(pt_size, true);
            selected_node_icon.setDrawColor(TreePanel.TREE_SELECTED_NODE_COLOR);
            selected_node_icon.setFillColor(col);
            
            TreePanel.DisplayedNode D = panel.getNode(node_idx);
            D.setIcon(true, selected_node_icon);
            D.setIcon(false, node_full);

            node_idx++;
        }
    }    
}
