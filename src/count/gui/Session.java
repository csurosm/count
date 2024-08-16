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
import java.awt.Component;
import java.awt.Container;
import java.awt.FileDialog;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.util.Set;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.SwingWorker;
import javax.swing.event.ChangeListener;
//import javax.swing.plaf.metal.MetalProgressBarUI;

import count.io.CountXML;
import count.io.DataFile;
import count.io.ModelBundle;
import count.io.SavableData;
import count.io.ExportableData;
import count.model.MixedRateModel;
import count.ds.AnnotatedTable;
import count.ds.Phylogeny;


/**
*
* Main component (JSplitPane) to show the 
* browsers (models and data) associated with a session. 
* 
* @author Mikl&oacute;s Cs&#369;r&ouml;s 
*/
public class Session extends JSplitPane
{
    public static final Color WORKTAB_TREES_COLOR = new Color(24,117,52);
    public static final Color WORKTAB_DATA_COLOR = new Color(108,24,76);
    public static final Color WORKTAB_RATES_COLOR = new Color(250,78,25);

    
    public static final Color WORKTAB_TREES_BACKGROUND = new Color(204,204,204); //  new Color(240,240,240); //   Color(240,255,240);
    public static final Color WORKTAB_DATA_BACKGROUND = WORKTAB_TREES_BACKGROUND;  //new Color(255,240,255);
    public static final Color WORKTAB_RATES_BACKGROUND = WORKTAB_TREES_BACKGROUND;//   new Color(255,240,240);
    
    public Session(AppFrame app, DataFile<Phylogeny> main_tree)
	{
    	this(app, null, main_tree);
	}
    
    public Session(AppFrame app, ModelBundle bundle)
    {
    	this(app, bundle, null);
    }
    
    private Session(AppFrame app, ModelBundle bundle, DataFile<Phylogeny> main_tree)
    {
    	super(HORIZONTAL_SPLIT);
        setDividerLocation(200+getInsets().top);
        setDebugGraphicsOptions(javax.swing.DebugGraphics.NONE_OPTION);
        setBorder(null);
        setResizeWeight(0.0);
        setOneTouchExpandable(true);

        work_area = new JTabbedPane(JTabbedPane.LEFT);
		browser_area = new Box(BoxLayout.PAGE_AXIS);

		this.app = app;
		
		if (bundle==null)
		{
	        this.session_key = getSessionKey(main_tree.getFile().getName());
	
	// BUNDLE		
	//        trees_browser = new TreeCollection(main_tree, new TabbedSplitPane(TREE_PANEL));
	//		
	//        data_browser = new Browser<>(AnnotatedTablePanel.class, new TabbedSplitPane(DATA_PANEL));
	//
	//        rates_browser = new Browser<>(RateVariationPanel.class, new TabbedSplitPane(RATES_PANEL));
	        
			session_bundle = new ModelBundle(session_key);
			session_bundle.getRoot().setAttribute(CountXML.ATT_TYPE, this.getClass().getName());
		} else
		{
			this.session_bundle = bundle;
			this.session_key = bundle.getRoot().getId();
		}
		
		model_browser = new BundleBrowser(session_bundle, BundleTree.TREES, new TabbedSplitPane(MODEL_PANEL));
		data_browser = new BundleBrowser(session_bundle, BundleTree.DATA, new TabbedSplitPane(DATA_PANEL));
		
        setLeftComponent(browser_area);
        setRightComponent(work_area);


        initComponents();
        
        if (bundle == null)
        	model_browser.addTree(main_tree, true);
        else
        {
        	model_browser.getBrowserComponent().copyBundleStructure();
        	data_browser.getBrowserComponent().copyBundleStructure();
        }
        colorBrowserBackgrounds();
	}
	
    private final ModelBundle session_bundle;
	private final AppFrame app;
	private final String session_key;
	private final JTabbedPane work_area;
	private final JComponent browser_area;

	// BUNDLE    
//    /**
//     * Title for the Rates tab
//     */
//    public static final String RATES_PANEL = "Rates";
    
    /**
     * Title for the Data tab
     */
    public static final String DATA_PANEL = "Data";

// BUNDLE    
//    /**
//     * Title for the Tree tab
//     */
//    public static final String TREE_PANEL = "Tree"; 
    
    /**
     * Title for the Model tab
     */
    public static final String MODEL_PANEL = "Models";
    
// BUNDLE        
//    /**
//     * Browser in the Data tab
//     */
//    private final Browser<AnnotatedTablePanel> data_browser;
//    
//    /**
//     * Browser in the Rate tab
//     */
//    private final Browser<RateVariationPanel> rates_browser;    
//    
//    /**
//     * Trees
//     */
//    private final TreeCollection trees_browser;
    
	/**
	 * Browser in the Data tab
	 */
    private final BundleBrowser data_browser;
	/**
	 * Browser in the Models tab
	 */
    private final BundleBrowser model_browser;
    
    /**
     * Our own display for multiple JSplitPanes: their left components ({@link BundleTree} with node selection)
     * are arranged on the left in the {@link #browser_area}, their right components 
     * (displaying the selected node) go into the tabbed pane of {@link #work_area}.
     *
     */
    private class TabbedSplitPane extends JSplitPane
    {
    	TabbedSplitPane(String tab_name)
    	{
    		super(JSplitPane.HORIZONTAL_SPLIT);
    		this.tab_name = tab_name;
    		this.tab_idx = work_area.getTabCount();
    		work_area.addTab(tab_name, new JPanel()); // for now 
    	}
    	private final String tab_name;
    	private int browser_idx = -1;
    	private final int tab_idx;
    	
    	@Override
    	public String toString() { return tab_name;}
    	
    	@Override 
    	public void setLeftComponent(Component C)
    	{
    		if (C instanceof JScrollPane) // with JTree in it 
    		{
    			JScrollPane scroll = (JScrollPane) C;
//    			scroll.getViewport().addPropertyChangeListener(chg->System.out.println("#**S.TSP.sLC vp "+chg+"\tfor "+scroll));
//    			scroll.getViewport().getView().addPropertyChangeListener(chg->System.out.println("#**S.TSP.sLC vw "+chg+"\tfor "+scroll));
//    			C.addPropertyChangeListener(chg->System.out.println("#**S.TSP.sLC "+chg+"\tfor "+scroll));
    			
	    		int num_browsers = browser_area.getComponentCount();
	    		if (browser_idx==-1)
	    		{
	    			browser_idx = browser_area.getComponentCount();
	    			browser_area.add(C);
//	    			System.out.println("#**S.TSP.sL add "+C+"\tidx "+browser_idx
//	    					+"\tbcount "+browser_area.getComponentCount()
//	    					+"\t//"+Thread.currentThread());
	    		} else 
	    		{
	    			Component[] browsers = browser_area.getComponents();
	    			browser_area.removeAll();
	    			int bi=0;
	    			for (; bi<browser_idx; bi++)
	    			{
	    				browser_area.add(browsers[bi]);
	    			}
	    			browser_area.add(C);
	    			for (; bi<browsers.length; bi++)
	    			{
	    				browser_area.add(browsers[bi]);
	    			}
//	    			System.out.println("#**S.TSP.sL set "+C+"\tidx "+browser_idx
//	    					+"\tbcount "+browser_area.getComponentCount()
//	    					+"\t//"+Thread.currentThread());
	    		}
        		browser_area.repaint();
    		} else 
    			super.setLeftComponent(C);
    	}
    	
    	@Override
    	public void setRightComponent(Component C)
    	{
    		if (tab_name == null) // premature call from super's instantiation
    			super.setRightComponent(C);
    		else
    		{
    			work_area.setComponentAt(tab_idx, C);
    			work_area.repaint();
    		}
    	}
    	
    	@Override
    	public Component getRightComponent()
    	{
    		if (tab_name == null)
    			return super.getRightComponent();
    		else 
    			return work_area.getComponentAt(tab_idx);
    	}
    }
	
    public static final String getSessionKey(String file_name)
    {
    	return DataFile.chopFileExtension(file_name);
    }
    
    public String getSessionKey()
    {
    	return session_key;
    }
	
    public String getSessionTitle()
    { 
    	StringBuilder session_name = new StringBuilder();
// BUNDLE
//    	session_name.append(trees_browser.getSelectedTree().toString());
//    	RateVariationPanel rates_panel = rates_browser.getSelectedPrimaryItem();
//    	if (rates_panel != null)
//    	{
//    		session_name.append(":").append(rates_panel.toString());
//    	}
//    	AnnotatedTablePanel data_panel = data_browser.getSelectedPrimaryItem();
    	JComponent tree_panel =  model_browser.getSelectedPrimaryItem(BundleTree.TREES);
    	if (tree_panel != null)
    	{
	    	String tree_panel_name = model_browser.getSelectedPrimaryItem(BundleTree.TREES).toString();
	    	session_name.append(tree_panel_name);
    	}
    	JComponent rates_panel = model_browser.getSelectedPrimaryItem(BundleTree.RATES);
    	if (rates_panel != null)
    	{
    		session_name.append(":").append(rates_panel.toString());
    	}
    	JComponent data_panel = data_browser.getSelectedPrimaryItem(BundleTree.TABLES);
// BUNDLE
    	if (data_panel != null)
    	{
    		session_name.append("@").append(data_panel.toString());
    	}
    	if (session_name.length()==0)
    		session_name.append(getSessionKey());
    	return session_name.toString();
    }
    
    public ModelBundle getModelBundle()
    {
    	data_browser.getBrowserComponent().setAllItemAttributes();
    	// model_browser.getBrowserComponent().setAllItemAttributes();
    	return session_bundle;
    }
	
	
    /**
     * Adds a new (possibly first) tree to the tree panel
     * @param tree_data 
     */

    public void addTree(DataFile<Phylogeny> tree_data)
    {
// BUNDLE
//        trees_browser.addTree(tree_data);
//        work_area.setSelectedIndex(work_area.indexOfTab(TREE_PANEL));
    	model_browser.addTree(tree_data, false);
    	work_area.setSelectedIndex(work_area.indexOfTab(MODEL_PANEL));
    }
    
    public void addTopTree(DataFile<Phylogeny> tree_data)
    {
// BUNDLE
//    	trees_browser.addTopTree(tree_data);
//        work_area.setSelectedIndex(work_area.indexOfTab(TREE_PANEL));
    	model_browser.addTree(tree_data, true);
    	work_area.setSelectedIndex(work_area.indexOfTab(MODEL_PANEL));
    }

    /**
     * Adds a new data set to the Data tab's browser and switches to it.
     * @param table_data
     */ 
    public void addDataSet(DataFile<AnnotatedTable> table_data)
    {
// BUNDLE
//    	AnnotatedTablePanel TP = new AnnotatedTablePanel(table_data);
//      data_browser.addTopItem(TP);
    	
//		// color for the profiles
//		JComponent selected_tree_embedding_panel =  getModelBrowser().getBrowserComponent().getSelectedNode(BundleTree.TREES).getComponent();
//    	if (selected_tree_embedding_panel instanceof Zoom<?>)
//    	{
//    		Zoom<? extends TreePanel> zum = (Zoom<? extends TreePanel>) selected_tree_embedding_panel;
//    		TreePanel tree_panel = zum.getTreePanel();
//    		int num_leaves = model_browser.getMainPhylogeny().getNumLeaves();
//    		for (int leaf=0; leaf<num_leaves; leaf++)
//    		{
//    			
//    		}
//    	}
		
		
		
    	JComponent table_component = data_browser.addTable(table_data, true);
    	
    	
    	
    	if (table_component instanceof AnnotatedTablePanel)
    	{
    		AnnotatedTablePanel table_panel = (AnnotatedTablePanel) table_component;
    		table_panel.setPhyleticProfileColoring(model_browser.getBrowserComponent().getLeafColors());
    	} else
    	{
    		throw new UnsupportedOperationException("Expecting "+AnnotatedTablePanel.class+"; got "+table_component.getClass());
    	}
        work_area.setSelectedIndex(work_area.indexOfTab(DATA_PANEL));
    }
    
    /**
     * Adds a new GammaInvariant model to the rates browser.
     */
    public void addRates(DataFile<MixedRateModel> rates_data, boolean add_at_top)
    {
// BUNDLE
//        RateVariationPanel RP = new RateVariationPanel(rates_data);
//        if (add_at_top)
//            rates_browser.addTopItem(RP);
//        else
//            rates_browser.addItem(RP);
        
    	model_browser.addRates(rates_data, add_at_top);
        work_area.setSelectedIndex(work_area.indexOfTab(MODEL_PANEL));
    }
    
    /**
     * Retrieves the closest ancestor that 
     * belongs to a given class. 
     * 
     * @param <A> ancestor class
     * @param component
     * @return null if component is null, or none of its ancestors belong to A
     */
    public static <A extends Container> A getAncestor(Component component, Class<A> ancestor_class)
    {
    	Container anc = (component==null?null:component.getParent());
    	
        while (anc != null)
        {
            if (ancestor_class.isInstance(anc))
            {
            	return ancestor_class.cast(anc);
            } 
            anc = anc.getParent();
        }
        // there is no such ancestor
    	return null;
    }
    
    /**
     * Retrieves the Session ancestor for a Component. 
     */
    public static Session getSession(Component C)
    {
    	return getAncestor(C, Session.class);
    }
    
    
    /**
     * Finds the enclosing AppFrame.
     * 
     * @param C A component within Count.
     * @return
     */
    public static AppFrame getApp(Component C)
    {
    	Session sesh = getSession(C);
    	if (sesh==null)
    		return null;
    	else
    		return sesh.app;
    }
    
    /**
     * Adds a change listener to tab selection.
     * 
     * @param listener
     */
    public void addChangeListener(ChangeListener listener)
    {
    	work_area.addChangeListener(listener);
    }
    
    public void addSelectionChangeListener(ChangeListener listener)
    {
    	this.addChangeListener(listener);
    	getModelBrowser().addSelectionChangeListener(listener);
    	getDataBrowser().addSelectionChangeListener(listener);
    }
    
    public BundleBrowser getModelBrowser() { return model_browser;}
    public BundleBrowser getDataBrowser() { return data_browser;}
    
// BUNDLE
//    public TreeCollection getTreesBrowser()
//    {
//    	return trees_browser;
//    }
//    
//    public Browser<AnnotatedTablePanel> getDataBrowser()
//    {
//    	return data_browser;
//    }
//    
//    public Browser<RateVariationPanel> getRatesBrowser()
//    {
//    	return rates_browser;
//    }
    
    /**
     * Constructs a version of the 
     * selected data table with column indexing mapped 
     * to the selected phylogeny. 
     * 
     * @return original table if main tree is selected, otherwise mapped to the selected tree
     */
    public DataFile<AnnotatedTable> getSelectedData()
    {
    	DataFile<AnnotatedTable> selected_data;
// BUNDLE
//    	AnnotatedTablePanel table_panel = data_browser.getSelectedPrimaryItem();
//    	if (table_panel == null)
//    	{
//    		selected_data = null;
//    	}
//    	else
//    	{
//    		DataFile<AnnotatedTable> orig_data = table_panel.getDataFile();
//    		TreePanel tree_panel = trees_browser.getSelectedTree();;
//    		IndexedTree selected_tree = tree_panel.getTreeData().getContent();
//    		if (selected_tree == trees_browser.getMainTree())
//    		{
////    			System.out.println("#**S.gSD for "+tree_panel.toString()+"\tsame");
//    			selected_data = orig_data;
//    		} else
//    		{ 
//	    		AnnotatedTable mapped_table = orig_data.getContent().mappedToTree(selected_tree);
//	    		selected_data = new DataFile<>(mapped_table,orig_data.getFile());
//    		}
//    	}
// BUNDLE
    	ModelBundle.Entry table_entry = data_browser.getSelectedTableEntry();
//    	System.out.println("#**S.gSD "+table_entry.getTableData().getFile());
    	if (table_entry==null)
    	{
    		selected_data = null;
    	} else
    	{
    		DataFile<AnnotatedTable> orig_data = table_entry.getTableData();
    		ModelBundle.Entry tree_entry = model_browser.getSelectedPrimaryEntry(BundleTree.TREES);
    		Phylogeny selected_tree =  tree_entry.getTreeData().getContent();
    		if (selected_tree == model_browser.getMainPhylogeny())
    		{
//    			System.out.println("#**S.gSD for "+tree_panel.toString()+"\tsame");
    			selected_data = orig_data;
    		} else
    		{ 
//    			System.out.println("#**S.gSD for "+tree_entry.getTreeData().getFile().toString()+"\tmapping "+orig_data.getFile()+"\tcontent "+orig_data.getContent());
    			
	    		AnnotatedTable mapped_table = orig_data.getContent().mappedToTree(selected_tree);
	    		selected_data = new DataFile<>(mapped_table,orig_data.getFile());
    		}
    	}
    	
    	
    	return selected_data;
    }
    
    
    public JComponent getSelectedComponent()
    {
    	return (JComponent) work_area.getSelectedComponent();
    }
    
    private void initComponents()
    {
// BUNDLE
//        work_area.setToolTipTextAt(work_area.indexOfTab(TREE_PANEL), "Phylogeny ");
////        work_area.setForegroundAt(work_area.indexOfTab(TREE_PANEL), WORKTAB_TREES_COLOR);
////        work_area.setBackgroundAt(work_area.indexOfTab(TREE_PANEL), WORKTAB_TREES_BACKGROUND);
////        work_area.add(DATA_PANEL, data_browser);
//        work_area.setToolTipTextAt(work_area.indexOfTab(DATA_PANEL), "Data set(s)");
////        work_area.setForegroundAt(work_area.indexOfTab(DATA_PANEL), WORKTAB_DATA_COLOR);
////        work_area.add(RATES_PANEL, rates_browser);
//        work_area.setToolTipTextAt(work_area.indexOfTab(RATES_PANEL), "Rate model(s)");
////        work_area.setForegroundAt(work_area.indexOfTab(RATES_PANEL), WORKTAB_RATES_COLOR);
//        
// BUNDLE
    	work_area.setToolTipTextAt(work_area.indexOfTab(DATA_PANEL), "Data set(s)");
    	work_area.setToolTipTextAt(work_area.indexOfTab(MODEL_PANEL), "Model(s)");
// BUNDLE

    	
//    	data_browser.getBrowserComponent().setBackground(WORKTAB_DATA_BACKGROUND);
        
// BUNDLE
//        trees_browser.getBrowserComponent().setBackground(WORKTAB_TREES_BACKGROUND);
//        trees_browser.addChangeListener(e->   // switch to a rate model that fits this phylogeny
//        			{
//        				TreePanel TP = trees_browser.getSelectedTree();
//        				IndexedTree TPtree = TP.getTreeData().getContent();
//        				
//        				RateVariationPanel RP = rates_browser.getSelectedPrimaryItem();
//        				boolean select_rates;
//        				if (RP==null)
//        					select_rates = true;
//        				else
//        				{
//        					GammaInvariant model = RP.getDataFile().getContent();
//        					IndexedTree RPtree = model.getBaseModel().getTree();
//        					select_rates = (RPtree != TPtree);
//        				}
//        				if (select_rates)
//        				{
//        					int nr = rates_browser.getTopItemCount();
//            				for (int i=0; i<nr; i++)
//        					{
//        						RP = rates_browser.getTopItem(i);
//        						GammaInvariant model = RP.getDataFile().getContent();
//            					IndexedTree RPtree = model.getBaseModel().getTree();
//            					if (RPtree==TPtree)
//            					{
//            						rates_browser.selectItem(RP);
//                    				return;
//            					} 
//            				}
//        					rates_browser.clearSelection(); // no such rates 
//        				}
//        			});
//        
//        rates_browser.getBrowserComponent().setBackground(WORKTAB_RATES_BACKGROUND);
//        rates_browser.addSelectionChangeListener(e-> // switch to this model's phylogeny
//        			{
//        				RateVariationPanel RP = rates_browser.getSelectedPrimaryItem();
//        				if (RP != null)
//        				{
//							GammaInvariant model = RP.getDataFile().getContent();
//	    					IndexedTree RPtree = model.getBaseModel().getTree();
//	
//	    					trees_browser.selectTree(RPtree);
//        				}
//        			});
// BUNDLE
//      model_browser.getBrowserComponent().setBackground(WORKTAB_TREES_BACKGROUND);
// BUNDLE
    	data_browser.addSelectionChangeListener(e-> // coupled selection changes in history views
    	{
    		ModelBundle.Entry selected_entry = data_browser.getSelectedEntry();
    		if (selected_entry != null)
    		{
				String rates_id = selected_entry.getAttributeValue(CountXML.ATT_RATES);
    			if (rates_id!=null)
    			{
    				model_browser.selectEntry(rates_id);
    			} else
    			{
    				String tree_id = selected_entry.getAttributeValue(CountXML.ATT_TREE);
    				if (tree_id!=null)
    				{
    					model_browser.selectEntry(tree_id);
    				}
    			}
    		}
    	});
    	
    	data_browser.addSelectionChangeListener(chg->work_area.setSelectedIndex(work_area.indexOfTab(DATA_PANEL)));
    	model_browser.addSelectionChangeListener(chg->work_area.setSelectedIndex(work_area.indexOfTab(MODEL_PANEL)));
    	
    	// add change selection listener to  model_browser: move up in data_browser, then 
    	// depth-first traversal to find matching rates or tree id
    	
        work_area.addChangeListener(chg->colorBrowserBackgrounds());
//					        {
//					        	int tab_idx= work_area.getSelectedIndex();
//					        	for (int ci=0; ci<browser_area.getComponentCount(); ci++)
//					        	{
//					        		JComponent browser = (JComponent)((JScrollPane ) browser_area.getComponent(ci)).getViewport().getView();		        		
//					        		if (ci==tab_idx)
//					        			browser.setBackground(Color.WHITE);
//					        		else
//										browser.setBackground(WORKTAB_TREES_BACKGROUND);
//					        	}
//					        	
//					        	//					        	scroll.requestFocus();
//					        	
////					        	Component W = work_area.getComponentAt(tab_idx);
////					        	Color bg = browse.getBackground();
//////					        	if (tab_idx==work_area.indexOfTab(DATA_PANEL))
//////					        		bg = WORKTAB_DATA_BACKGROUND;
//////					        	else if (tab_idx==work_area.indexOfTab(TREE_PANEL))
//////					        		bg = WORKTAB_TREES_BACKGROUND;
//////					        	else
//////					        	{
//////					        		assert (tab_idx == work_area.indexOfTab(RATES_PANEL));
//////					        		bg = WORKTAB_RATES_BACKGROUND;
//////					        	}
////					        	W.setBackground(bg);
//					        });
        
        
        
        
        for (int ci=0; ci<browser_area.getComponentCount(); ci++)
        {
        	Component C = browser_area.getComponent(ci);
        	JScrollPane scroll = (JScrollPane ) C;
        	JComponent browse = (JComponent) scroll.getViewport().getView();
        	browse.setOpaque(true);
        	class BrowserListener implements // FocusListener, 
        		MouseListener
        	{
        		private final int ci;
        		private final JComponent browser;
        		BrowserListener(int ci)
        		{
        			this.ci = ci;
        			this.browser = (JComponent)((JScrollPane ) browser_area.getComponent(ci)).getViewport().getView();
        		}
//				@Override
//				public void focusGained(FocusEvent e)
//				{
////					browser.setBackground(Color.WHITE);
////					int selected=  work_area.getSelectedIndex();
////					if (selected>=0 && selected != ci) // avoids infinite loop via selection change listener that sets focs here 
////					{
////						work_area.setSelectedIndex(ci);
////					}
////					for (int cj=0; cj<browser_area.getComponentCount(); cj++)
////					{
////						if (ci!=cj)
////						{
////							JComponent other_browser = (JComponent)((JScrollPane ) browser_area.getComponent(cj)).getViewport().getView(); 
////							other_browser.setOpaque(true);
////							other_browser.setBackground(WORKTAB_TREES_BACKGROUND);
////						}
////					}
//				}
//
//				@Override
//				public void focusLost(FocusEvent e) 
//				{
////					browser.setOpaque(true);
////					browser.setBackground(WORKTAB_TREES_BACKGROUND);
//				}
				
				@Override
				public void mouseClicked(MouseEvent e) 
				{
					int selected=  work_area.getSelectedIndex();
					if (selected>=0 && selected != ci) // avoids infinite loop via selection change listener that sets focs here 
					{
						work_area.setSelectedIndex(ci);
					}
				}

				@Override
				public void mousePressed(MouseEvent e) 
				{
					browser.setBackground(Color.WHITE);
				}

				@Override
				public void mouseReleased(MouseEvent e) {
					// TODO Auto-generated method stub
					
				}

				@Override
				public void mouseEntered(MouseEvent e) {
					// TODO Auto-generated method stub
					
				}

				@Override
				public void mouseExited(MouseEvent e) {
					// TODO Auto-generated method stub
					
				}				
        		
        	}
        	BrowserListener listener = new BrowserListener(ci);
//        	browse.addFocusListener(listener);
        	browse.addMouseListener(listener);
        	
        	
        	
        	
        	
//        	browse.requestFocus();
//        	Color tabcol = work_area.getForegroundAt(ci);
//        	
//        	int light = 2;
//        	Color bgcol = new Color((tabcol.getRed()+light*255)/(light+1),
//        							(tabcol.getGreen()+light*255)/(light+1),
//        							(tabcol.getBlue()+light*255)/(light+1));
////        	browse.setBackground(bgcol);
//        	browse.setForeground(tabcol);
        	
        }
    }		
    
    /**
     * Called when selection changes between model and data. 
     */
    private void colorBrowserBackgrounds()
    {
		int tab_idx= work_area.getSelectedIndex();
		for (int ci=0; ci<browser_area.getComponentCount(); ci++)
		{
			BundleTree browser = (BundleTree)((JScrollPane ) browser_area.getComponent(ci)).getViewport().getView();		        		
			if (ci==tab_idx)
				browser.setBackground(Color.WHITE);
			else
				browser.setBackground(WORKTAB_TREES_BACKGROUND);
			browser.repaint();
		}
    }
    
    
    /**
     * Required interface for allowing the "extract filtered families" action
     * 
     *
     */
    public static interface FamilySelection
    {
    	public abstract int[] getSelectedFamilies();
    }

    /**
     * Adds a data tabkle witjh the filtered families, Current selection in 
     * {@link #data_browser} must have a {@link FamilySelection} associated component. 
     */
    public void showFilteredFamilies()
    {
    	JComponent comp = data_browser.getSelectedItem();
    	if (comp instanceof FamilySelection)
    	{
        	FamilySelection selection_panel = (FamilySelection) comp;
        	int[] selected = selection_panel.getSelectedFamilies();
// BUNDLE
//        	AnnotatedTablePanel table_panel = data_browser.getSelectedPrimaryItem();
//        	DataFile<AnnotatedTable> table_data = table_panel.getDataFile();
        	ModelBundle.Entry table_entry = data_browser.getSelectedPrimaryEntry(BundleTree.TABLES);
        	DataFile<AnnotatedTable> table_data = table_entry.getTableData();
// BUNDLE
        	AnnotatedTable filtered_table = table_data.getContent().filteredTable(selected);
        	File filtered_file = new File((File)null, "filt:"+table_data.getFile().getName());
        	DataFile<AnnotatedTable> filtered_data = new DataFile<>(filtered_table, filtered_file);
// BUNDLE
//        	AnnotatedTablePanel filtered_panel =  new AnnotatedTablePanel(filtered_data);
//        	data_browser.addItem(filtered_panel);
        	data_browser.addTable(filtered_data, false);
// BUNDLE
	    	work_area.setSelectedIndex(work_area.indexOfTab(DATA_PANEL));
    	} else
    	{
    		throw new IllegalArgumentException();
    	}
    }
    
    public void showSplitTable()
    {
    	throw new UnsupportedOperationException();
    }
    
    public void showBinaryProfiles()
    {
// BUNDLE
//    	AnnotatedTablePanel table_panel = data_browser.getSelectedPrimaryItem();
//    	if (table_panel==null)
//    	{
//    		throw new IllegalArgumentException();
//    	}
//    	DataFile<AnnotatedTable> table_data = table_panel.getDataFile();
    	ModelBundle.Entry table_entry = data_browser.getSelectedTableEntry(); //data_browser.getSelectedPrimaryEntry(BundleTree.TABLES);
    	
    	if (table_entry==null)
    	{
    		throw new IllegalArgumentException();
    	}
    	DataFile<AnnotatedTable> table_data = table_entry.getTableData();
// BUNDLE
    	AnnotatedTable binary = table_data.getContent().binaryTable();
    	File binary_file = new File((File)null, "01:"+table_data.getFile().getName());
    	DataFile<AnnotatedTable> binary_data = new DataFile<>(binary, binary_file);
// BUNDLE
//    	AnnotatedTablePanel binary_panel = new AnnotatedTablePanel(binary_data);
//    	data_browser.addItem(binary_panel);
    	data_browser.addTable(binary_data, false);
// BUNDLE
    	work_area.setSelectedIndex(work_area.indexOfTab(DATA_PANEL));
    }
    
    /**
     * Brings up a model selection/optimization parameters dialog.
     */
    public void doOptimize()
    {
    	ModelSelectionDialog optD = new ModelSelectionDialog(app, "Rate model optimization");
    	RateOptimizationPanel optP = optD.showOptimization(); 
    	if (optP != null)
    	{
// BUNDLE
//	    	if (optP.isDescendantModel())
//	    	{
//	    		rates_browser.addItem(optP);
//	    	} else
//	    	{
//	    		rates_browser.addTopItem(optP);
//	    	}
//	    	optP.addPropertyChangeListener(RateOptimizationPanel.PROPERTY_OPTIMIZATION_DONE,
//	    			change->{
//			    			rates_browser.selectItem(optP);
//			    	    	work_area.setSelectedIndex(work_area.indexOfTab(RATES_PANEL));
//			    	    	// updates menus 
//			    	});
//	    	work_area.setSelectedIndex(work_area.indexOfTab(RATES_PANEL));
    		model_browser.addRates(optP, !optP.isDescendantModel());
	    	optP.addPropertyChangeListener(RateOptimizationPanel.PROPERTY_OPTIMIZATION_DONE,
	    			change->{
			    			model_browser.selectItem(optP);
			    	    	work_area.setSelectedIndex(work_area.indexOfTab(MODEL_PANEL));
			    	    	// updates menus 
			    	});
	    	work_area.setSelectedIndex(work_area.indexOfTab(MODEL_PANEL));
    		
// BUNDLE
    	}
    }
    
    public void doDollo()
    {
// BUNDLE
//    	DataFile<? extends IndexedTree> phylo_data = trees_browser.getSelectedTree().getTreeData();
//        DataFile<AnnotatedTable> table_data = getSelectedData();
//    	DolloView D = new DolloView(phylo_data, table_data);
//    	trees_browser.decorateByMainTree(D.getTreePanel());
//    	D.computeAll();
//    	data_browser.addItem(D);
    	ModelBundle.Entry phylo_entry =  model_browser.getSelectedTreeEntry();
    	DataFile<Phylogeny> phylo_data = phylo_entry.getTreeData();
    	DataFile<AnnotatedTable> table_data = this.getSelectedData();
    	DolloView D = new DolloView(phylo_data, table_data);
//    	model_browser.decorateByMainTree(D.getTreePanel());
    	data_browser.addHistoryItem(D, phylo_entry, null);
//    	D.computeAll();
    	// BUNDLE
    	D.setPhyleticProfileColoring(model_browser.getBrowserComponent().getLeafColors());    	

    	work_area.setSelectedIndex(work_area.indexOfTab(DATA_PANEL));
    }
    
    public void doParsimony()
    {
// BUNDLE
//    	DataFile<? extends IndexedTree> phylo_data = trees_browser.getSelectedTree().getTreeData();
//        DataFile<AnnotatedTable> table_data = getSelectedData();
//    	ParsimonyView P = new ParsimonyView(phylo_data, table_data);
//    	
//    	trees_browser.decorateByMainTree(P.getTreePanel());
//    	P.computeAll();
//    	data_browser.addItem(P);
    	ModelBundle.Entry phylo_entry =  model_browser.getSelectedPrimaryEntry(BundleTree.TREES);
    	DataFile<Phylogeny> phylo_data = phylo_entry.getTreeData();
    	DataFile<AnnotatedTable> table_data = this.getSelectedData();
    	ParsimonyView P = new ParsimonyView(phylo_data, table_data);
//    	model_browser.decorateByMainTree(P.getTreePanel());
    	data_browser.addHistoryItem(P, phylo_entry, null);
//    	P.computeAll();
// BUNDLE
    	P.setPhyleticProfileColoring(model_browser.getBrowserComponent().getLeafColors());    	
    	work_area.setSelectedIndex(work_area.indexOfTab(DATA_PANEL));
    }
    
    public void doPosteriors()
    {
    	ModelSelectionDialog optD = new ModelSelectionDialog(app, "Rate model and algorithm selection");
    	PosteriorsView P = optD.showPosteriors();
    	if (P!=null)
    	{
// BUNDLE
//	    	trees_browser.decorateByMainTree(P.getTreePanel());
//	    	data_browser.addItem(P);
//    		model_browser.decorateByMainTree(P.getTreePanel());
    		ModelBundle.Entry rate_entry = model_browser.getSelectedPrimaryEntry(BundleTree.RATES);
    		data_browser.addHistoryItem(P, null, rate_entry);
//    		P.computeAll();
    		
// BUNDLE
    		P.setPhyleticProfileColoring(model_browser.getBrowserComponent().getLeafColors());
	    	work_area.setSelectedIndex(work_area.indexOfTab(DATA_PANEL));
    	}
    	
//    	RateVariationPanel rate_panel = rates_browser.getSelectedPrimaryItem();
//    	AnnotatedTablePanel table_panel = data_browser.getSelectedPrimaryItem();
//    	DataFile<AnnotatedTable> data_table = getSelectedData();
//    	PosteriorsView P = new PosteriorsView(rate_panel.getDataFile(), data_table);
//    	trees_browser.calculateNodeColors(P.getTreePanel());
//    	data_browser.addItem(P);
//    	P.computeAll();
//    	work_area.setSelectedIndex(work_area.indexOfTab(DATA_PANEL));
    }
    
    public void doSimulation()
    {
    	ModelSelectionDialog optD = new ModelSelectionDialog(app, "Rate model and sample size");
    	SimulationView S = optD.showSimulation();
    	if (S!=null)
    	{
// BUNDLE
//    		data_browser.addTopItem(S.asTable());
//    		trees_browser.decorateByMainTree(S.getTreePanel());
//    		data_browser.addItem(S);
    		ModelBundle.Entry rate_entry = model_browser.getSelectedPrimaryEntry(BundleTree.RATES);
    		data_browser.addTopTable(S.asTable(), rate_entry);
//    		model_browser.decorateByMainTree(S.getTreePanel());
    		data_browser.addHistoryItem(S, null, rate_entry);
//    		S.computeAll();
// BUNDLE
        	S.setPhyleticProfileColoring(model_browser.getBrowserComponent().getLeafColors());    	
    		work_area.setSelectedIndex(work_area.indexOfTab(DATA_PANEL));
    	}
    }
    

    public void doSaveData()
    {
    	Object selected = getSelectedComponent();
    	final String action_name;
    	String file_name = null;
    	String directory = null;
    	final SavableData<?> selected_savable;
    	final ExportableData selected_exportable;
    	
    	if (selected instanceof SavableData)
    	{
    		selected_exportable = null;
    		selected_savable = (SavableData)selected;
    		File data_file = selected_savable.getDataFile().getFile();
    		action_name = "Save";
    		
    		
    		file_name = data_file.getName();
    		directory = data_file.getParent();
    	} else if (selected instanceof ExportableData)
    	{
    		selected_exportable = (ExportableData) selected;
    		selected_savable = null;
    		action_name = "Export";
    	} else
    	{
    		assert(false); // should never get here, but final variables need this branch 
    		action_name = null;
    		selected_exportable = null;
    		selected_savable = null;
    	}
//		System.out.println("#**S.dSV start "+action_name+"\tsav "+selected_savable+"\texp "+selected_exportable);
    	
    	
    	if (action_name != null)
    	{
    		
    		// TODO : use progress monitor on save
    		FileDialog dialog = new FileDialog(app, action_name+" \""+selected.toString()+"\"" , FileDialog.SAVE);
			if (file_name!=null)
				dialog.setName(file_name);
    		if (directory !=null)
    			dialog.setDirectory(directory);
            dialog.setVisible(true);
            file_name = dialog.getFile();
            directory = dialog.getDirectory();                
    		if (file_name != null)
    		{
                final File save_file = new File(directory, file_name); 
    			
    			
                final JDialog wait_until_saved = new JDialog(app,dialog.getTitle(),true);
                BoxLayout saving_layout = new BoxLayout(wait_until_saved.getContentPane(),BoxLayout.PAGE_AXIS);
                wait_until_saved.getContentPane().setLayout(saving_layout);
                wait_until_saved.add(Box.createVerticalGlue());
                wait_until_saved.add(new JLabel("Writing data to disk: please wait ..."));
                JProgressBar progress = new JProgressBar(0,100);
                
                
                progress.setUI(new javax.swing.plaf.metal.MetalProgressBarUI()); // native proressbar does not animate on MacOS
                progress.setIndeterminate(true); 
            	progress.setVisible(true);
            	progress.setStringPainted(true);

                Box progress_box = new Box(BoxLayout.LINE_AXIS);
                progress_box.add(Box.createHorizontalStrut(12));
                progress_box.add(progress);
                progress_box.add(Box.createHorizontalStrut(12));
                wait_until_saved.add(progress_box);
                wait_until_saved.add(Box.createVerticalGlue());

                SwingWorker<Void,Void> save_worker = new SwingWorker<Void,Void>()
                {
                    @Override
                    public Void doInBackground()
                    {
                        try 
                        {
                        	progress.setValue(0);                        	
    		    			if (selected instanceof SavableData)
    		    			{
    		                    selected_savable.saveData(save_file);
    		    			} else if (selected instanceof ExportableData)
    		    			{
    		    				selected_exportable.saveData(save_file);
    	                    } 
                        }
                        catch (java.io.IOException E)
                        {
                        	app.getExceptionHandler().handle(E, "I/O error", "File error while "+action_name+"ing data.");
                        } catch (Exception E)
                        {
                        	app.getExceptionHandler().handle(E, "A bug maybe?", "Error while "+action_name+"ing data.");
                        }
                        return null;
                    }

                    @Override
                    public void done()
                    {
                        wait_until_saved.setVisible(false);
                        app.setTitle(getSessionTitle());   
                    }    
                }; 
                wait_until_saved.pack();
                java.awt.Dimension frameD = app.getSize();
                wait_until_saved.setBounds(frameD.width/3, frameD.height/3, frameD.width/3, frameD.height/6);
                save_worker.execute();
                wait_until_saved.setVisible(true);             			

//                try 
//                {
//	    			if (selected instanceof SavableData)
//	    			{
//	                    File save_file = new File(directory, file_name); 
//	                    selected_savable.saveData(save_file);
//	    			} else if (selected instanceof ExportableData)
//	    			{
//	    				File save_file = new File(directory, file_name);
//	    				selected_exportable.saveData(save_file);
//                    } 
//	    			
////	    			System.out.println("#**S.dSV action "+action_name+"\tfile "+file_name);
//                } catch (IOException trouble)
//                {
//                    app.getExceptionHandler().handle(trouble, "I/O error", "Error while writing data.");
//                }
    		}
//			System.out.println("#**S.dSV done "+action_name+"\tfile "+file_name);
    	}
    }
    
    public boolean isSelectedRemovable()
    {
    	BundleBrowser selected_browser = getSelectedBrowser();
    	boolean isSelectedRemovable = false;
    	if (selected_browser != null)
    	{
    		isSelectedRemovable = selected_browser.getSelectedItem()!=null
    				 && !selected_browser.isMainTreeSelected();
    	}
    	return isSelectedRemovable;
    }
    
    public BundleBrowser getSelectedBrowser()
    {
    	BundleBrowser selected_browser = null;
    	if (work_area.getSelectedIndex() == work_area.indexOfTab(DATA_PANEL))
    	{
    		selected_browser = data_browser;
    	} else if (work_area.getSelectedIndex() == work_area.indexOfTab(MODEL_PANEL))
    	{
    		selected_browser = model_browser;
    	}
    	return selected_browser;

    	
    }
    
    public void doRemoveItem()
    {
    	BundleBrowser selected_browser = getSelectedBrowser();
    	
    	if (selected_browser != null)
    	{
        	Object selected = getSelectedComponent();
        	if (selected != selected_browser.getSelectedItem())
        	{ // DEBUG
        		System.out.println("#**S.dRI selected "+selected+"\titem "+selected_browser.getSelectedItem());
        	}
        	assert (selected == selected_browser.getSelectedItem());
        	if (selected_browser.isMainTreeSelected())
        	{
    			JOptionPane.showMessageDialog(app, "Cannot remove the main tree for the session. Close the session instead." , "Removal impossible", JOptionPane.ERROR_MESSAGE);
        		
        	} else
        	{
	        	Set<String> selected_ids = selected_browser.getAllSelectedEntries();
	    		if (data_browser.mayRemove(selected_ids) && model_browser.mayRemove(selected_ids))
	    		{
		        	String msg = "Please confirm: do you really want to remove \""+selected.toString()+"?\""
		        			+" This operation will also remove all descendant views, if any.";
		        	int confirm = JOptionPane.showConfirmDialog(app, msg, "Removing "+selected.toString(), JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);
		        	if (confirm == JOptionPane.YES_OPTION)
		        	{
		        		selected_browser.removeSelectedItem();
		        	}
	    		} else
	    		{
	    			JOptionPane.showMessageDialog(app, "Cannot remove the selected views because there are dependencies (ancestral reconstruction / trees, models) outside the selection. Delete the dependent views first." , "Removal impossible", JOptionPane.ERROR_MESSAGE);
	    		}
        	}
    	}
    }
}
