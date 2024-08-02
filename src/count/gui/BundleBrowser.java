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

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Font;

import java.util.Enumeration;
import java.util.HashSet;
import java.util.Set;
import java.util.function.Predicate;

import javax.swing.JComponent;
import javax.swing.JEditorPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.TreeModelEvent;
import javax.swing.event.TreeModelListener;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;

import count.ds.AnnotatedTable;
import count.ds.Phylogeny;
import count.gui.BundleTree.Node;
import count.io.CountXML;
import count.io.DataFile;
import count.io.ModelBundle;
import count.io.ModelBundle.Entry;
import count.model.MixedRateModel;

/**
 * GUI component (JPanel with JSplitPane) for displaying a 
 * {@link BundleTree} and following the node selection. 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */
public class BundleBrowser extends JPanel
{
    public static JSplitPane defaultSplitPane()
    {
    	JSplitPane content = new JSplitPane();
        content = new JSplitPane();
        content.setDividerLocation(200+content.getInsets().top);
        content.setOrientation(JSplitPane.HORIZONTAL_SPLIT);
        content.setDebugGraphicsOptions(javax.swing.DebugGraphics.NONE_OPTION);
        content.setBorder(null);
        content.setResizeWeight(0.0);
        content.setOneTouchExpandable(true);
        return content;
    }
    
    /**
     * Instantiation for a standalone component. 
     * 
     * @param bundle
     * @param top_selector
     */
    public BundleBrowser(ModelBundle bundle, Predicate<Entry> top_selector)
    {
    	this(bundle, top_selector, defaultSplitPane());
    	setLayout(new BorderLayout());
        add(content, BorderLayout.CENTER);
    }

    /**
     * Instantiation for placing the browser into a JSplitPane (left component 
     * is the BundleTree, right component is the selected node's associated JComponent.
     * 
     * @param bundle
     * @param top_selector
     * @param container
     */
    public BundleBrowser(ModelBundle bundle, Predicate<Entry> top_selector, JSplitPane container)
    {
    	super();
    	this.item_browser = new BundleTree(bundle, top_selector);
    	initComponents(container);
    }
    
//    protected static JSplitPane embedBrowser(ModelBundle bundle, Predicate<Entry> top_selector, JSplitPane embedPane)
//    {
//    	BundleBrowser browser = new BundleBrowser(bundle, top_selector, embedPane);
//    	return browser.content;
//    }
    
    /**
     * Main element in the panel. 
     * Left-hand side contains the browser, the right-hand side 
     * shows the information associated with the selected 
     * item. 
     */
    private JSplitPane content;
    
    /**
     * The browsing is implemented via JTree.
     */
    private final BundleTree item_browser;
    /**
     * The text field where messages are written when there is no selected item.
     */
    private JEditorPane no_selection_text;
    
    private void initComponents(JSplitPane container)
    {
    	int div = container.getDividerLocation();
    	this.content = container;

        item_browser.addTreeSelectionListener(e->showSelectedItem());
    	
        item_browser.setFont(new java.awt.Font("Serif", Font.BOLD,14));
        item_browser.setFocusable(true); // is true by default on JTree
        item_browser.setVisible(true);
        
        JScrollPane item_browser_scroll = new JScrollPane();
        item_browser_scroll.setViewportView(item_browser);
        item_browser_scroll.setMinimumSize(new Dimension(200,200));
        content.setLeftComponent(item_browser_scroll);        

        no_selection_text = new JEditorPane();
        no_selection_text.setContentType("text/html");
        no_selection_text.setBackground(getBackground().brighter());
        no_selection_text.setEditable(false);
        content.setRightComponent(no_selection_text);

        content.setDividerLocation(div);    	
    }
    
    public BundleTree getBrowserComponent()
    {
    	return item_browser;
    }
    
    
    public Phylogeny getMainPhylogeny()
    {
    	return item_browser.getMainPhylogeny();
//    	Entry root_entry = item_browser.getRoot().getEntry();
//    	Entry first_tree = root_entry.getChild(0, BundleTree.TREES);
//    	Phylogeny getMainTree = first_tree==null?null:first_tree.getTreeData().getContent();
//    	return getMainTree;
    }
    
    public boolean isMainTreeSelected()
    {
//    	Entry root_entry = item_browser.getRoot().getEntry();
//    	Entry first_tree = root_entry.getChild(0, BundleTree.TREES);
    	Entry first_tree = item_browser.getMainTree();
    	return getSelectedEntry() == first_tree;
    	
    }
    
    /**
     * Adds a change listener that responds to changes of 
     * node selection or structure change. 
     * 
     * @param listener
     */
    public void addSelectionChangeListener(ChangeListener listener)
    {
        class ChangeAdapter implements TreeSelectionListener, TreeModelListener
        {
            private final ChangeListener listener;
            
            ChangeAdapter(ChangeListener L)
            {
                listener = L;
            }

            private void structureChanged()
            {
                JComponent src = getSelectedItem();
                if (src == null)
                    src = BundleBrowser.this;
//                System.out.println("#**BB.aSCL.CA.sC "+src+"\tsnode "+item_browser.getSelectedNode());
                
                listener.stateChanged(new ChangeEvent(src));
            }
            
            @Override
            public void treeNodesChanged(TreeModelEvent e)
            {
                // nothing
            }
            
            @Override
            public void treeNodesInserted(TreeModelEvent e)
            {
                structureChanged();
            }
            @Override
            public void treeNodesRemoved(TreeModelEvent e)
            {
                structureChanged();
            }
            @Override
            public void treeStructureChanged(TreeModelEvent e)
            {
                // nothing
            }
            
            @Override
            public void valueChanged(TreeSelectionEvent treeSelectionEvent) 
            {
                JComponent daddy = getSelectedItem();
                // DEBUG
//                System.out.println("#**BB.aSCL.CA.vC "+daddy+"\tsnode "+item_browser.getSelectedNode());
                if (daddy != null)
                    listener.stateChanged(new ChangeEvent(daddy));
            }    	
        }
        ChangeAdapter pisl = new ChangeAdapter(listener);
        item_browser.addTreeSelectionListener(pisl);
        item_browser.getModel().addTreeModelListener(pisl);
    }
    
//	public JComponent getSelectedItem()
//	{
//		Node selected_node = getSelectedNode();
//		JComponent getSelectedItem = (selected_node==null?null:selected_node.getComponent());
//		return getSelectedItem;
//	}
//	
//	public Entry getSelectedEntry()
//	{
//		Node selected_node = getSelectedNode();
//		Entry getSelectedEntry = selected_node==null?null:selected_node.getEntry();
//		return getSelectedEntry;
//		
//	}
//	
//	public JComponent getSelectedPrimaryItem(Predicate<Entry> primary_selector)
//	{
//		Node selected_node = getSelectedNode(primary_selector);
//		JComponent getSelectedItem = (selected_node==null?null:selected_node.getComponent());
//		return getSelectedItem;
//	}
    
    public JComponent getSelectedItem()
    {
//		System.out.println("#**BB.getSI "+item_browser.getSelectedNode());
    	return getSelectedPrimaryItem(e->true);
    }
    
    public JComponent getSelectedPrimaryItem(Predicate<Entry> primary_type)
    {
		Node selected_node = item_browser.getSelectedNode(primary_type);
		
//		System.out.println("#**BB.getSP "+selected_node+"\ttreesel "+item_browser.getSelectedNode());
		
		JComponent getSelectedItem = (selected_node==null?null:selected_node.getComponent());
		return getSelectedItem;
    }
    
    public ModelBundle.Entry getSelectedEntry()
    {
    	return getSelectedPrimaryEntry(e->true);
    }
    
    public ModelBundle.Entry getSelectedPrimaryEntry(Predicate<Entry> primary_type)
    {
		Node selected_node = item_browser.getSelectedNode(primary_type);
		Entry getSelectedEntry = selected_node==null?null:selected_node.getEntry();
		return getSelectedEntry;
    }
    
    public ModelBundle.Entry getSelectedRatesEntry()
    {
    	return getSelectedPrimaryEntry(BundleTree.RATES);    	
    }
    public ModelBundle.Entry getSelectedTreeEntry()
    {
    	return getSelectedPrimaryEntry(BundleTree.TREES);    	
    }
    public ModelBundle.Entry getSelectedTableEntry()
    {
    	return getSelectedPrimaryEntry(BundleTree.TABLES);
    }
    
    public Entry getTree(String name)
    {
    	Node node = item_browser.findTree(name);
    	Entry getTree = node==null?null:node.getEntry();
    	return getTree;
    }
    
    public boolean selectItem(JComponent item)
    {
    	Node selected_node = item_browser.findItem(item);
    	boolean selectItem = selected_node != null;
    	if (selectItem)
    	{
    		item_browser.setSelectedNode(selected_node);
    	}
    	
    	return selectItem;
    }
    
    public boolean selectEntry(String id)
    {
    	Node selected_node = item_browser.findEntry(id);
    	boolean selectEntry = selected_node != null;
    	if (selectEntry)
    		item_browser.setSelectedNode(selected_node);
    	return selectEntry;
    }
    
    /**
     * Calculates the set of ids that are in the subtree of the
     * selected node. 
     * 
     * @return non-null set
     */
    public Set<String> getAllSelectedEntries()
    {
    	Set<String> getAllSelected = new HashSet<>();
    	Node selected_node = item_browser.getSelectedNode();
    	Enumeration<?> selected_subtree = selected_node.depthFirstEnumeration();
    	while (selected_subtree.hasMoreElements())
    	{
    		Node node = (Node) selected_subtree.nextElement();
    		Entry entry = node.getEntry();
    		getAllSelected.add(entry.getId());
    	}
    	return getAllSelected;
    }
    
    /**
     * Checks if a non-removed node contains a tree or model reference 
     * to the removed entries within the entire {@link #item_browser} tree. 
     * 
     * @param removed_entries set of ids marked for deletion
     * @return true if there are no reference to the specified set
     */
    public boolean mayRemove(Set<String> removed_entries)
    {
    	Enumeration<?> all_nodes = item_browser.getRoot().depthFirstEnumeration();
    	boolean mayRemove = true;
    	while (all_nodes.hasMoreElements() && mayRemove)
    	{
    		Node node = (Node) all_nodes.nextElement();
    		Entry entry = node.getEntry();
    		boolean to_be_removed = removed_entries.contains(entry.getId());
    		if (!to_be_removed && entry.isTableEntry() && entry.getTableData().getContent()==null
    				&& !to_be_removed)
    		{
	    		if (entry.getAttributeValue(CountXML.ATT_TREE) != null)
	    		{
	    			String tree_id = entry.getAttributeValue(CountXML.ATT_TREE);
	    			boolean tree_removed = removed_entries.contains(tree_id);
	    			mayRemove = mayRemove && !tree_removed;
	    		}
	    		if (entry.getAttributeValue(CountXML.ATT_RATES) != null)
	    		{
	    			String rates_id = entry.getAttributeValue(CountXML.ATT_RATES);
	    			boolean rates_removed = removed_entries.contains(rates_id);
	    			mayRemove = mayRemove && !rates_removed;
	    		}
    		}
    	}
    	return mayRemove;
    }    
    
    public boolean removeSelectedItem()
    {
    	boolean removeSelectedItem = item_browser.removeSelectedNode();
    	
    	if (removeSelectedItem)
    	{
    		showSelectedItem();
    	}
    	return removeSelectedItem;
    }
    
//	public void addItem(JComponent component)
//	{
//		// 1. (Session.addRates) RateVariationPanel from rates_data
//		// 2. (Session.doDollo) DolloView / HistoryView [launched]
//		// 3. (Session.doOptimize) RateOptimizationPanel / RateVariationPanel [launched]
//		// 4. (Session.doParsimony) ParsimonyView / HistoryView [launched]
//		// 5. (Session.doPosteriors) PosteriorsView / HistoryView [launched]
//		// 6. (Session.doSimulation) SimulationView / HistoryView  [launched]
//		// 7. (Session.showBinaryProfiles) AnnotatedTablePanel
//		// 8. (Session.showFilteredFamilies) AnnotatedTablePanel 
//		
//		
//		// top1: (Session.addDataSet) AnnotatedTablePanel
//		// top2. (Session.addRates) RateVariationPanel from rates_data
//		// top2: (Session.doSimulation) AnnotatedTablePanel [launched]
//		// top4. (Session.doOptimize) RateOptimizationPanel / RateVariationPanel [launched]
//		
//	}
    
    
	public JComponent addTable(DataFile<AnnotatedTable> table_data, boolean at_root)
	{
		return addTable(table_data, null, at_root, null, null);
	}
	
	public JComponent addTable(AnnotatedTablePanel table_panel, boolean at_root)
	{
		return addTable(null, table_panel, at_root, null, null);
	}
	
	public JComponent addTopTable(AnnotatedTablePanel table_panel, Entry related_rates)
	{
		return addTable(null, table_panel, true, null, related_rates);
	}
	
	private JComponent addTable(DataFile<AnnotatedTable> table_data, AnnotatedTablePanel table_panel, boolean at_root, Entry related_table, Entry related_rates)
	{
		if (table_data==null)
			table_data = table_panel.getDataFile();
		Node parent = at_root?item_browser.getRoot():item_browser.getSelectedNode(BundleTree.TABLES);
		if (parent == null)
			return null;
		Node addTable = parent.addTable(table_data, table_panel);
		if (related_table != null)
		{
			addTable.getEntry().setAttribute(CountXML.ATT_TABLE, related_table.getId());
		}
		if (related_rates != null)
		{
			addTable.getEntry().setAttribute(CountXML.ATT_RATES, related_rates.getId());
		}
		item_browser.setSelectedNode(addTable);
		item_browser.requestFocus();
		return addTable.getComponent();		
	}	
	public JComponent addTree(DataFile<Phylogeny> tree_data, boolean at_root)
	{
		Node parent = at_root?item_browser.getRoot():item_browser.getSelectedNode(BundleTree.TREES);
		if (parent == null)
			return null;
		Node addTree = parent.addTree(tree_data);
		item_browser.setSelectedNode(addTree);
		item_browser.requestFocus();
		
		
		return addTree.getComponent();
	}
	
	public JComponent addRates(DataFile<MixedRateModel> rates_data, boolean at_root)
	{
		return addRates(rates_data, null, at_root);
//		Node parent = at_root?getSelectedNode(TREES):getSelectedNode(RATES);
//		if (parent == null)
//			return null;
//		Node addRates = parent.addRates(rates_data, null);
//		setSelectedNode(addRates);
//		return addRates.getComponent();
	}
	
	public JComponent addRates(RateVariationPanel rates_panel, boolean at_root)
	{
		return addRates(null, rates_panel, at_root);
//		DataFile<GammaInvariant> rates_data = rates_panel.getDataFile();
//		Node parent = at_root?getSelectedNode(TREES):getSelectedNode(RATES);
//		if (parent == null)
//			return null;
//		Node addRates = parent.addRates(rates_data, rates_panel);
//		setSelectedNode(addRates);
//		return addRates.getComponent();
	}
	
	private JComponent addRates(DataFile<MixedRateModel> rates_data, RateVariationPanel rates_panel, boolean at_root)
	{
		if (rates_data==null)
			rates_data = rates_panel.getDataFile();
		Node parent = at_root?item_browser.getSelectedNode(BundleTree.TREES):item_browser.getSelectedNode(BundleTree.RATES);
		if (parent == null)
			return null;
		Node addRates = parent.addRates(rates_data, rates_panel);
		item_browser.setSelectedNode(addRates);
		item_browser.requestFocus();
		return addRates.getComponent();
	}
	
	public JComponent addHistoryItem(HistoryView view_panel, Entry related_phylo_entry, Entry related_rates_entry)
	{
		Node parent = item_browser.getSelectedNode(BundleTree.TABLES);
		if (parent == null)
			return null;
		Node addHistoryItem = parent.addHistoryItem(view_panel, related_phylo_entry, related_rates_entry);

//		decorateByMainTree(view_panel.getTreePanel());
		
		item_browser.setSelectedNode(addHistoryItem);
		item_browser.requestFocus();
		
		
		return addHistoryItem.getComponent();
	}
	
//	public void decorateByMainTree(TreePanel panel)
//	{
//		item_browser.decorateByMainTree(panel);
//	}
//    
    
    
    /**
     * The message displayed when there are no available items.
     */
    private static final String NO_AVAILABLE_ITEMS = "<em>No available items.</em>";
    
    /**
     * The message displayed when there is no selection.
     */
    private static final String SELECT_AN_ITEM = "<em>Select an item.</em>";
    

    /**
     * Called whenever item selection changes
     */
    private void showSelectedItem()
    {
        JComponent selected = getSelectedItem();
        int div = content.getDividerLocation(); // we will reset so that the divider doesn't change 
        if (selected == null)
        {
            if (item_browser.getRoot().getChildCount()==0)
                no_selection_text.setText(NO_AVAILABLE_ITEMS);
            else
                no_selection_text.setText(SELECT_AN_ITEM);
            content.setRightComponent(no_selection_text);
            content.setDividerLocation(div);
        } else 
        {
            content.setRightComponent(selected);
            content.setDividerLocation(div);
        }        
        
//        System.out.println("#**BB.showSI "+selected+"\tnode "+item_browser.getSelectedNode());
    }        
    
    
    
}
