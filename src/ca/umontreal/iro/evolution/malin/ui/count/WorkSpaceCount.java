/*
 * WorkSpaceCount.java
 *
 * Created on March 2, 2008, 4:16 PM
 *
 */

package ca.umontreal.iro.evolution.malin.ui.count;

import java.util.Map;
import java.util.Hashtable;
import java.util.List;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.StringReader;

import java.util.zip.GZIPInputStream;

import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.helpers.XMLReaderFactory;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;

import javax.swing.JComponent;
import javax.swing.JFrame;

import javax.swing.tree.DefaultMutableTreeNode;

import ca.umontreal.iro.evolution.Parser;
import ca.umontreal.iro.evolution.TreeNode;

import ca.umontreal.iro.evolution.malin.DataFile;
import ca.umontreal.iro.evolution.malin.IndexedPoint;

import ca.umontreal.iro.evolution.malin.ui.Browser;
import ca.umontreal.iro.evolution.malin.ui.Dealer;
import ca.umontreal.iro.evolution.malin.ui.DrawString;
import ca.umontreal.iro.evolution.malin.ui.EmbellishedTreePanel;
import ca.umontreal.iro.evolution.malin.ui.WorkSpace;
import ca.umontreal.iro.evolution.malin.ui.ZoomableTreePanel;

import ca.umontreal.iro.evolution.genecontent.NodeWithRates;
import ca.umontreal.iro.evolution.genecontent.OccurrenceTable;
import ca.umontreal.iro.evolution.genecontent.RateVariation;
import ca.umontreal.iro.evolution.genecontent.TreeWithRates;

/**
 * Main workspace class with three tabs: tree, data, rates. 
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class WorkSpaceCount extends WorkSpace
{
    
    public WorkSpaceCount(Dealer D, TreeNode root, File path_to_tree_file) 
    {
        super(D);
        setPhylogeny(root, path_to_tree_file);
    }

    /**
     * Title for the Rates tab
     */
    public static final String RATES_PANEL = "Rates";
    
    /**
     * Title for the Data tab
     */
    public static final String DATA_PANEL = "Data";

    /**
     * Tree with which this work space was initialized
     */
    private TreeNode original_tree_root;
    
    /**
     * The tree displayed in the Tree tab
     */
    private TreeWithRates main_tree;
    

    private Map<NodeWithRates, TreeNode> node_mapping;
    /**
     * Browser in the Data tab
     */
    protected Browser rates_panel;
    
    /**
     * Browser in the Rate tab
     */
    protected Browser data_panel;    


    /**
     * Display for the tree
     */
    private ZoomableTreePanel zoomed_tree;
    
    /**
     * Finds the DealerCount parent of this guy.
     */
    public DealerCount getDealerCount()
    {
        return (DealerCount)super.getDealer();
    }
    
    public static DealerCount getDealerCount(Component C)
    {
        return (DealerCount)WorkSpace.getDealer(C);
    }
    
    public static WorkSpaceCount getWorkSpaceCount(Component C)
    {
        return (WorkSpaceCount)WorkSpace.getWorkSpace(C);
    }


    /**
     * Returns the view/primary item/tree panel that the user is looking at now
     */
    public Component getDisplayedComponent()
    {
        if (work_area.getTabCount()==0)
            return null;
        else
        {
            int tab_idx = work_area.getSelectedIndex();
            String tab_title =  work_area.getTitleAt(tab_idx);
            if (DATA_PANEL.equals(tab_title))
            {
                Component C = (Component) getDatasetsBrowser().getSelectedItem();
                if (C==null)
                    return getDatasetsBrowser();
                else
                    return C;
            }
            else if (RATES_PANEL.equals(tab_title))
            {
                Component C = (Component) getRatesBrowser().getSelectedItem();
                if (C==null)
                    return getRatesBrowser();
                else
                    return C;
            }
            else if (TREE_PANEL.equals(tab_title))
                return tree_panel;
            return work_area.getComponentAt(tab_idx);
        }
    }

    /**
     * Returns the browser associated with the Rates tab.
     */
    public Browser getRatesBrowser()
    {
        return rates_panel;
    }
    
    /**
     * Returns the browser associated with the Data tab.
     */
    public Browser getDatasetsBrowser()
    {
        return data_panel;
    }    

    
    public FamilySizeTableDisplay getSelectedTableDisplay()
    {
        return (FamilySizeTableDisplay)
                data_panel.getSelectedPrimaryItem();
    }

    public FamilySizeTableDisplay getSelectedRootTableDisplay()
    {
        return (FamilySizeTableDisplay)
                data_panel.getSelectedPrimaryRootItem();
    }
    
    public Object getSelectedPrimaryDataDisplay()
    {
        return data_panel.getSelectedPrimaryItem();
    }
            
    public DataFile<OccurrenceTable> getSelectedDataset()
    {
        FamilySizeTableDisplay TD = getSelectedTableDisplay();
        if (TD==null)
            return null;
        else 
            return TD.getData();
    }
    
    public RateModelDisplay getSelectedRateModelDisplay()
    {
        return (RateModelDisplay) rates_panel.getSelectedPrimaryItem();
    }
    
    public DataFile<RateVariation> getSelectedRateModel()
    {
        RateModelDisplay RD = getSelectedRateModelDisplay();
        if (RD==null)
            return null;
        else
            return RD.getRateModel();
    }
    
    /**
     * Initializes the main three tabs; displays the tree on the Tree tab.
     * The origial edge lengths are ignored in the main tree.
     */
    protected void setPhylogeny(TreeNode root, File path_to_tree_file)
    {
        original_tree_root = root;
        main_tree = new TreeWithRates(root);
        computeNodeMapping();
        
        work_area.removeAll();
        
        tree_panel = new TRPanel();
        tree_panel.setBackground(java.awt.Color.WHITE);
        tree_panel.setAssociatedFile(path_to_tree_file);
        tree_panel.setToolTipText("Main phylogeny ("+path_to_tree_file.getName()+")");
        //main_panel_scroll = new JScrollPane(tree_panel);
        //main_panel_scroll.getViewport().setBackground(tree_panel.getBackground());
        
        zoomed_tree = new ZoomableTreePanel((EmbellishedTreePanel)tree_panel);
        addToWorkArea(TREE_PANEL, zoomed_tree);
        
        data_panel = new Browser(Browser.PrimaryItem.class);
        addToWorkArea(DATA_PANEL,data_panel);
        rates_panel = new Browser(RateModelDisplay.class);
        //rates_panel = new Browser(RateModelDisplay.class);
        addToWorkArea(RATES_PANEL, rates_panel);
        
        work_area.setToolTipTextAt(work_area.indexOfTab(TREE_PANEL), "Outline of the main phylogeny ("+path_to_tree_file.getName()+")");
        work_area.setToolTipTextAt(work_area.indexOfTab(DATA_PANEL), "Data set(s)");
        work_area.setToolTipTextAt(work_area.indexOfTab(RATES_PANEL), "Rate model(s)");
        if (this.getDealer().isMac())
        {
            work_area.setForegroundAt(work_area.indexOfTab(TREE_PANEL), new Color(24,117,52));
            work_area.setForegroundAt(work_area.indexOfTab(DATA_PANEL), new Color(108,24,76));
            work_area.setForegroundAt(work_area.indexOfTab(RATES_PANEL), new Color(250,78,25));
        }

        checkPhylogeny();
    }
    
    /**
     * Initializes the node_mapping hash table
     */
    private void computeNodeMapping()
    {
        node_mapping = main_tree.getTreeNodeMapping(original_tree_root);
//
//
//                new Hashtable<NodeWithRates,TreeNode>();
//        // all nodes
//        TreeNode[] original_nodes = original_tree_root.getTraversal().getDFT();
//        NodeWithRates[] new_nodes = main_tree.getDFT();
//        // map the leaves
//        Hashtable<String,TreeNode> original_leaves = new Hashtable<String,TreeNode>();
//        for (int node_idx=0; node_idx<original_nodes.length; node_idx++)
//        {
//            TreeNode N = original_nodes[node_idx];
//            if (N.isLeaf())
//                original_leaves.put(N.getName(), N);
//        }
//        for (int node_idx=0; node_idx<new_nodes.length; node_idx++)
//        {
//            NodeWithRates N = new_nodes[node_idx];
//            if (N.isLeaf())
//            {
//                String name = N.getName();
//                TreeNode oN = original_leaves.get(name);
//                node_mapping.put(N, oN);
//            } else
//            {
//                // an inner node
//                // the trees have the same topology: it's enough to look at one child
//                NodeWithRates C = (NodeWithRates)N.getChild(0);
//                TreeNode oC = node_mapping.get(C);
//                TreeNode oN = oC.getParent();
//                node_mapping.put(N, oN);
//            }
//        }
    }

    /**
     * The original tree node that corresponds to a node in the main tree.
     * 
     * @param node_index node index in the main tree
     * @return the corresponding original node
     */
    public TreeNode getOriginalNode(int node_index)
    {
        return node_mapping.get(main_tree.getNode(node_index));
    }
    
    public double getOriginalEdgeLength(int node_index)
    {
        return getOriginalNode(node_index).getLength();
    }
    
    public double getOriginalEdgeLength(NodeWithRates N)
    {
        return node_mapping.get(N).getLength();
    }
    
    /**
     * Adds a new data set to the Data tab's browser and switches to it.
     * 
     * @param ot table of family sizes defining the data set 
     * @param path file associated with this data set (where it was loaded from)
     */
    public void addDataSet(OccurrenceTable ot, File path)
    {
        FamilySizeTableDisplay TD = new FamilySizeTableDisplay(new DataFile<OccurrenceTable>(ot, path));
        
        data_panel.addTopItem(TD);
        selectWorkAreaElement(DATA_PANEL);
    }
    
    public boolean selectDataItem(JComponent item)
    {
        boolean found_item = data_panel.selectItem(item);
        if (found_item)
            selectWorkAreaElement(DATA_PANEL);
        
        return found_item;
    }
    
    public boolean selectRatesItem(JComponent item)
    {
        boolean found_item = rates_panel.selectItem(item);
        if (found_item)
            selectWorkAreaElement(RATES_PANEL);
            
        
        return found_item;
    }
    /**
     * Adds a new rate variation model to the Rates tab's browser as aq top element, and switches the focus to it.
     * 
     * @param rates the rate variation model
     * @param path file associated with the rate variation model (e.g., where it was loaded from)
     */
    public void addRates(RateVariation rates, File path)
    {
        addRates(rates, path, true);
    }
    /**
     * Adds a new rate variation model to the Rates tab's browser and switches the focus to it.
     * 
     * @param rates the rate variation model
     * @param path file associated with the rate variation model (e.g., where it was loaded from)
     * @param add_at_top if this is to be a top item, or a derived rates model
     */
    public void addRates(RateVariation rates, File path, boolean add_at_top)
    {
        RateModelDisplay RD = new RateModelDisplay(new DataFile<RateVariation>(rates, path));
        if (add_at_top)
            rates_panel.addTopItem(RD);
        else
            rates_panel.addItem(RD);
        selectWorkAreaElement(RATES_PANEL);
    }
    
    public void showFilteredFamilies()
    {
        Object d_selected = data_panel.getSelectedItem();
        if (d_selected !=null && (d_selected instanceof FilterableFamilies))
        {
            FilterableFamilies f_selected = (FilterableFamilies) d_selected;
            JComponent fd = f_selected.newDisplayWithSelectedFamilies();
            if (fd != null)
            {
                data_panel.addItem(fd);
                selectWorkAreaElement(DATA_PANEL);
            }
        } else
        {
            java.awt.Toolkit.getDefaultToolkit().beep();
        }
            
    }
    
    public void showBinaryProfiles()
    {
        FamilySizeTableDisplay td = getSelectedTableDisplay();
        if (td != null)
        {
            FamilySizeTableDisplay binary_td = td.createBinaryTableDisplay();
            data_panel.addItem(binary_td);
            selectWorkAreaElement(DATA_PANEL);
        } else
            java.awt.Toolkit.getDefaultToolkit().beep();
        
    }
    
    public void showFamilySharing()
    {
        
    }
    
    public void showDollo()
    {
        DataFile<OccurrenceTable> T = getSelectedDataset();
        DolloDisplay dd = new DolloDisplay(main_tree, T, this);
        getSelectedTableDisplay().addTableModelListener(dd.getFamilyTableScroll());
        data_panel.addItem(dd);
        selectWorkAreaElement(DATA_PANEL);
    }
    
    public void showWagner()
    {
        showWagner(Double.NaN);
    }
    
    private void showWagner(double pty)
    {
        WagnerDisplay wd = null;
        if (Double.isNaN(pty))
            wd = new WagnerDisplay(main_tree, getSelectedTableDisplay().getData(), this);//,1.1);
        else 
            wd = new WagnerDisplay(main_tree, getSelectedTableDisplay().getData(), pty, this);
        getSelectedTableDisplay().addTableModelListener(wd.getFamilyTableScroll());
        data_panel.addItem(wd);
        selectWorkAreaElement(DATA_PANEL);
    }
    
    public void showPosteriors()
    {
        DataFile<RateVariation> R = getSelectedRateModel();
        DataFile<OccurrenceTable> T = getSelectedDataset();
        PosteriorDisplay pd = new PosteriorDisplay(T, R, this);
        getSelectedTableDisplay().addTableModelListener(pd.getFamilyTableScroll());
        data_panel.addItem(pd);
        selectWorkAreaElement(DATA_PANEL);
    }
    
    public void showPGL()
    {
        DataFile<OccurrenceTable> T = getSelectedDataset();
        PGLDisplay pd = new PGLDisplay(main_tree, T, this);
        getSelectedTableDisplay().addTableModelListener(pd.getFamilyTableScroll());
        data_panel.addItem(pd);
        selectWorkAreaElement(DATA_PANEL);
    }
        
    public void optimizeRates(JFrame frame)
    {
        DataFile<RateVariation> R = getSelectedRateModel();
        DataFile<OccurrenceTable> T = getSelectedDataset();
        RateOptimizationUI optimization_ui=new RateOptimizationUI(main_tree, T,R);
        final RateOptimizationUI.ModelDisplay optimized_model_display = optimization_ui.optimizeModel(getDealerCount());
        if (optimized_model_display != null)
        {
            optimized_model_display.addPropertyChangeListener("done",new PropertyChangeListener()
            {
                @Override
                public void propertyChange(PropertyChangeEvent E)
                {
                    if (E.getNewValue() == Boolean.TRUE)
                    {
                        rates_panel.selectItem(optimized_model_display);
                        selectWorkAreaElement(RATES_PANEL);
                        
                        getDealerCount().initMenu();
                    }
                }
            });

            if (optimized_model_display.isDerivedModel())
                rates_panel.addItem(optimized_model_display);
            else
                rates_panel.addTopItem(optimized_model_display);
            selectWorkAreaElement(RATES_PANEL);
        }
    }
    
    private String session_id;
    
    public void saveSession(PrintStream out, String session_id)
    {
        this.session_id = session_id;
        out.println("<"+EMT_SESSION+" id=\""+session_id+"\" >");
        saveTree(out);
        Hashtable<DataFile<RateVariation>, String> rate_id = saveRates(out);
        saveData(rate_id, out);
        
        out.println("</"+EMT_SESSION+">");
    }
    
    private void saveTree(PrintStream out)
    {
        out.println("<"+EMT_TREE+" id=\""+session_id+".T \" name=\""+getAssociatedFile().getName()+"\" >");
        String newick = original_tree_root.newickTree(true, true, false, false);
        out.println("<![CDATA[\n"+newick+"\n]]>");
        out.println("</"+EMT_TREE+">");
    }
    
    private Hashtable<DataFile<RateVariation>, String> saveRates(PrintStream out)
    {
        out.println("<"+EMT_RATES+" id=\""+session_id+".R\">");
        
        int model_idx = 0;
        
        Hashtable<DataFile<RateVariation>, String> rate_id = new Hashtable<DataFile<RateVariation>, String>();
        Browser.ItemEnumeration traversal = rates_panel.preorderEnumeration();
        while (traversal.hasMoreElements())
        {
            Browser.EnumeratedItem thingy = traversal.nextElement();
            if (thingy.getItem() instanceof RateModelDisplay )
            {
                RateModelDisplay display = (RateModelDisplay) thingy.getItem();
                
                DataFile<RateVariation> model = display.getRateModel();
                if (model == null) // e.g., optimization is still running
                    continue;

                String current_id = session_id+".R"+Integer.toString(model_idx);
                model_idx++;
                rate_id.put(model, current_id);
                
                File current_file = model.getFile();
                String parent_info = session_id+".R"; // default parent is the root
                { // let's find if there is a valid ancestor
                    DefaultMutableTreeNode parent_node = thingy.getNode();
                    do 
                    {
                        parent_node = (DefaultMutableTreeNode)parent_node.getParent();
                        if (parent_node != null)
                        {
                            Object parent_item = parent_node.getUserObject();
                            if (parent_item != null && parent_item instanceof RateModelDisplay)
                            {
                                RateModelDisplay parent_display = (RateModelDisplay) parent_item;
                                DataFile<RateVariation> parent_model = parent_display.getRateModel();
                                if (parent_model != null)
                                    parent_info = rate_id.get(parent_model);

                            }
                        }
                    } while (parent_node != null && "".equals(parent_info));
                }
                
                out.println("<"+EMT_RATEMODEL+" id=\""+current_id+"\" " +
                        "parent=\""+parent_info+"\" "+
                        "name=\""+current_file.getName()+"\" " +
                        ">");
                String model_description = model.getData().tableRates();
                out.println("<![CDATA[");
                out.println(model_description);
                out.println("]]>");
                out.println("</"+EMT_RATEMODEL+">");
            }
            
        }
        out.println("</"+EMT_RATES+">");
        
        return rate_id;
    }
    private void saveData(Hashtable<DataFile<RateVariation>,String> rate_id, PrintStream out) 
    {
        out.println("<"+EMT_DATA+" id=\""+session_id+".D\">");
        Browser.ItemEnumeration traversal = data_panel.preorderEnumeration();
        Hashtable<Object, String> item_identifiers = new Hashtable<Object, String>();
        int current_item_idx = 0;
        while (traversal.hasMoreElements())
        {
            Browser.EnumeratedItem thingy = traversal.nextElement();
            Object display_item = thingy.getItem();
            String current_id = session_id+".D"+Integer.toString(current_item_idx);
            current_item_idx++;
            item_identifiers.put(display_item, current_id);
            
            Object parent_item = thingy.getParent();
            String parent_info = session_id+".D"; // default parent is the root
            if (parent_item != null)
                parent_info = item_identifiers.get(parent_item);
            
            if (display_item instanceof FamilySizeTableDisplay)
            {
                FamilySizeTableDisplay TD = (FamilySizeTableDisplay) display_item;
                DataFile<OccurrenceTable> table = TD.getData();
 
                OccurrenceTable T = table.getData();
                out.println("<"+EMT_TABLE+" " +
                        "id=\""+current_id+"\" " +
                        "parent=\""+parent_info+"\" "+
                        "name=\""+table.getFile().getName()+"\" " +
                        //"propertiesCount=\""+T.getKnownPropertiesCount()+"\" "+
                        ">");
                String table_content = table.getData().getFormattedTable(true);
                out.println("<![CDATA[");                
                out.println(table_content);
                out.println("]]>");
                out.println("</"+EMT_TABLE+">");
            } else if (display_item instanceof PosteriorDisplay)
            {
                PosteriorDisplay PD = (PosteriorDisplay) display_item;
                DataFile<RateVariation> model = PD.getRateModel();
                
                String ri = rate_id.get(model);
                
                out.println("<"+EMT_POSTERIORS+" " +
                        "id=\""+current_id+"\" " +
                        "parent=\""+parent_info+"\" "+
                        "rates=\""+ri+"\" " +
                        ">");
                out.println("<![CDATA[");                
                PD.savePosteriors(out);
                out.println("]]>");
                out.println("</"+EMT_POSTERIORS+">");
            } else 
            { // everything else will be recomputed
                String attributes = "";
                if (display_item instanceof WagnerDisplay)
                {
                    WagnerDisplay WD = (WagnerDisplay) display_item;
                    attributes = "gain=\""+WD.getGainPenalty()+"\"";
                } 
                out.println("<"+EMT_DATAITEM+" " +
                        "id=\""+current_id+"\" " +
                        "type=\""+display_item.getClass().getCanonicalName()+"\" " +
                        "parent=\""+parent_info+"\" "+
                        attributes+" "+
                        "/>");
            }
        }
        
        
        out.println("</"+EMT_DATA+">");
    }
    //public static final String DTD_FPID = "\"-//Csuros//DTD count 1.01//EN\"";
    //public static final String DTD_URI = "\"http://www.iro.umontreal.ca/~csuros/count/count1.dtd\"";
    public static final String EMT_TOP = "CountML";
    public static final String EMT_SESSION = "session";
    public static final String EMT_TREE = "tree";
    public static final String EMT_DATA = "data";
    public static final String EMT_RATES = "rates";
    public static final String EMT_RATEMODEL = "model";
    public static final String EMT_TABLE = "table";
    public static final String EMT_POSTERIORS = "posteriors";
    public static final String EMT_DATAITEM = "item";
    
    public static void loadSessions(File f, DealerCount dealer) throws IOException, SAXException
    {
        // check if compressed file ...
        String file_name = f.getName();
        int dot_position = file_name.lastIndexOf('.');
        BufferedReader BR = null;
        if (dot_position != -1)
        {
            String extension = file_name.substring(dot_position+1);
            if ("gz".equals(extension))
                BR = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f))));
        }
        
        if (BR==null)
            BR = new BufferedReader(new FileReader(f));
        
        
        XMLReader R = XMLReaderFactory.createXMLReader();
        R.setContentHandler(new CountMLHandler(dealer));
        R.parse(new InputSource(BR));
        BR.close();
    }
    
    
    /**
     * SAX ContentHandler for CountML documents
     */
    public static class CountMLHandler extends DefaultHandler
    {
        public CountMLHandler(DealerCount dealer)
        {
            this.dealer = dealer;
            sessions = new Hashtable<String,WorkSpaceCount>();
            tables = new Hashtable<String, FamilySizeTableDisplay>();
            rates = new Hashtable<String, RateModelDisplay>();
        }
        
        private DealerCount dealer;

        private boolean reading_tree;
        private String current_tree_name;
        private StringBuffer tree_sb;

        private boolean reading_table;
        private StringBuffer table_sb;
        private Hashtable<String,FamilySizeTableDisplay> tables;
        private String current_table_id;
        private String current_table_parent;
        private String current_table_name;

        private boolean reading_rate_model;
        private StringBuffer rate_sb;
        private Hashtable<String, RateModelDisplay> rates;
        private String current_rate_id;
        private String current_rate_parent;
        private String current_rate_name;
        
        private Attributes current_item_attributes;

        private boolean reading_posteriors;
        private StringBuffer posteriors_sb;
        private String current_posteriors_table;
        private String current_posteriors_rates;
        
        private String active_session;
        private Hashtable<String,WorkSpaceCount> sessions;
        private String current_session_id;
        
        
        
        @Override
        public void startElement(String uri,
                  String localName,
                  String qName,
                  Attributes attributes)
                  throws SAXException      
        {
            String elementName=qName;
            if (!"".equals(uri)) elementName=localName;
            
            if (EMT_TOP.equals(elementName))
                startTop(attributes);
            else if (EMT_SESSION.equals(elementName))
                startSession(attributes);
            else if (EMT_TREE.equals(elementName))
                startTree(attributes);
            else if (EMT_DATA.equals(elementName))
                startData(attributes);
            else if (EMT_RATES.equals(elementName))
                startRates(attributes);
            else if (EMT_TABLE.equals(elementName))
                startTable(attributes);
            else if (EMT_RATEMODEL.equals(elementName))
                startRateModel(attributes);
            else if (EMT_POSTERIORS.equals(elementName))
                startPosteriors(attributes);
            else if (EMT_DATAITEM.equals(elementName))
                startDataItem(attributes);
            else
                throw new SAXException("Undefined start tag for "+elementName);
        }
        
        @Override
        public void endElement(String uri,
                String localName,
                String qName)
                throws SAXException
        {
            String elementName=qName;
            if (!"".equals(uri)) elementName=localName;

            try 
            {
                if (EMT_TOP.equals(elementName))
                    endTop();
                else if (EMT_SESSION.equals(elementName))
                    endSession();
                else if (EMT_TREE.equals(elementName))
                    endTree();
                else if (EMT_DATA.equals(elementName))
                    endData();
                else if (EMT_RATES.equals(elementName))
                    endRates();
                else if (EMT_TABLE.equals(elementName))
                    endTable();
                else if (EMT_RATEMODEL.equals(elementName))
                    endRateModel();
                else if (EMT_POSTERIORS.equals(elementName))
                    endPosteriors();
                else if (EMT_DATAITEM.equals(elementName))
                    endDataItem();
                else
                    throw new SAXException("Undefined end tag for "+elementName);
            } catch (IOException E)
            {
                throw new SAXException(E);
            } catch (Parser.ParseException E)
            {
                throw new SAXException(E);
            }
        }
        
        @Override
        public void characters(char[] ch,
                int start,
                int length)
                throws SAXException
        {
            if (reading_tree || reading_table || reading_rate_model || reading_posteriors)
            {
                
                String reading_what = (reading_tree?"tree ":"");
                reading_what = reading_what + (reading_table?"table ":"");
                reading_what = reading_what + (reading_rate_model?"rates ":"");
                reading_what = reading_what + (reading_posteriors?"posteriors":"");

                //System.out.println("#*WSC.CMLH.ch reading "+reading_what);
                
                if (reading_tree)
                    tree_sb.append(ch, start, length);
                else if (reading_table)
                    table_sb.append(ch, start, length);
                else if (reading_rate_model)
                    rate_sb.append(ch,start,length);
                else if (reading_posteriors)
                    posteriors_sb.append(ch,start,length);
            } else
            {
                //System.out.println("#*WSC.CMLH.ch "+length+" characters skipped");
            }
        }
        
        /**
         * Start tag for top element
         * 
         * @param attributes associated attributes
         */
        private void startTop(Attributes attributes)
        {
            active_session = attributes.getValue("activesession");
            //System.out.println("#*WSC.CMLH.sT active "+active_session);
        }
        
        /**
         * Start tag for session
         * 
         * @param attributes associated attributes
         */
        private void startSession(Attributes attributes)
        {
            current_session_id = attributes.getValue("id");
        }
        
        /**
         * Start tag for tree
         * 
         * @param attributes associated attributes
         */
        private void startTree(Attributes attributes)
        {
            reading_tree = true;
            tree_sb = new StringBuffer();
            current_tree_name = attributes.getValue("name");
        }
        
        /**
         * Start tag for data
         * 
         * @param attributes associated attributes
         */
        private void startData(Attributes attributes)
        {
            WorkSpaceCount ws = sessions.get(current_session_id);
            ws.selectWorkAreaElement(DATA_PANEL);
        }
        
        
        /**
         * Start tag for rates
         * 
         * @param attributes associated attributes
         */
        private void startRates(Attributes attributes)
        {
            WorkSpaceCount ws = sessions.get(current_session_id);
            ws.selectWorkAreaElement(RATES_PANEL);
        }

        /**
         * Start tag for table
         * 
         * @param attributes associated attributes
         */
        private void startTable(Attributes attributes)
        {
            reading_table = true;
            table_sb = new StringBuffer();
            current_table_id = attributes.getValue("id");
            current_table_parent = attributes.getValue("parent");
            current_table_name = attributes.getValue("name");
        }
        
        /**
         * Start tag for rate model
         * 
         * @param attributes associated attributes
         */
        private void startRateModel(Attributes attributes)
        {
            reading_rate_model = true;
            rate_sb = new StringBuffer();
            current_rate_id = attributes.getValue("id");
            current_rate_parent = attributes.getValue("parent");
            current_rate_name = attributes.getValue("name");
        }
        
        /**
         * Start tag for posteriors
         * 
         * @param attributes associated attributes
         */
        private void startPosteriors(Attributes attributes)
        {
            reading_posteriors = true;
            posteriors_sb = new StringBuffer();
            current_posteriors_table = attributes.getValue("parent");
            current_posteriors_rates = attributes.getValue("rates");
        }
        
        /**
         * Start tag for data item
         * 
         * @param attributes associated attributes
         */
        private void startDataItem(Attributes attributes)
        {
            current_item_attributes = attributes;
        }
        
        /**
         * End tag for top element
         */
        private void endTop()
        {
        }
        
        
        /**
         * End tag for session
         */
        private void endSession()
        {
            current_session_id=null;
        }
        
        /**
         * End tag for tree
         */
        private void endTree() throws IOException, Parser.ParseException
        {
            StringReader R = new StringReader(tree_sb.toString());
            TreeNode main_phylo = Parser.readNewick(R,false,true,false);
            File file = new File((String)null, current_tree_name);
            WorkSpaceCount ws = (WorkSpaceCount)dealer.newWorkSpace(main_phylo, file);
            dealer.addWorkSpace(ws);
            sessions.put(current_session_id, ws);
            reading_tree=false;
            tree_sb = null;
        }
        
        /**
         * End tag for data
         */
        private void endData()
        {
            
        }
        
        /**
         * End tag for rates
         */
        private void endRates()
        {
            
        }

        /**
         * End tag for table
         */
        private void endTable() throws IOException
        {
            StringReader R = new StringReader(table_sb.toString());
            WorkSpaceCount ws = sessions.get(current_session_id);
            OccurrenceTable T = new OccurrenceTable(ws.main_tree.getLeaves());
            T.readTable(R, true);
            File f = new File((String)null, current_table_name);
            FamilySizeTableDisplay F = new FamilySizeTableDisplay(new DataFile<OccurrenceTable>(T, f));
            if (tables.containsKey(current_table_parent))
            {
                // could verify if parent really is a key in sessions, but we assume that it is
                FamilySizeTableDisplay mommy = tables.get(current_table_parent);
                ws.data_panel.selectItem(mommy);
                ws.data_panel.addItem(F);
            } else
            {
                ws.data_panel.addTopItem(F);
            }
            tables.put(current_table_id, F);
            
            reading_table = false;
            table_sb = null;
        }
    
        /**
         * End tag for rate model
         */
        private void endRateModel() throws IOException
        {
            StringReader R = new StringReader(rate_sb.toString());
            WorkSpaceCount ws = sessions.get(current_session_id);

            RateVariation model = RateVariation.read(R, ws.main_tree);
            File f = new File((String)null, current_rate_name);
            RateModelDisplay RD = new RateModelDisplay(new DataFile<RateVariation>(model, f));

            if (rates.containsKey(current_rate_parent))
            {
                RateModelDisplay mommy = rates.get(current_rate_parent);
                ws.rates_panel.selectItem(mommy);
                ws.rates_panel.addItem(mommy);
            } else
            {
                ws.rates_panel.addTopItem(RD);
            }
            rates.put(current_rate_id, RD);
            
            reading_rate_model = false;
            rate_sb = null;
            
        }
        /**
         * End tag for posteriors
         */
        private void endPosteriors() throws IOException
        {
            WorkSpaceCount ws = sessions.get(current_session_id);

            FamilySizeTableDisplay FD = tables.get(current_posteriors_table);
            ws.data_panel.selectItem(FD);            
            PosteriorDisplay pd = new PosteriorDisplay(ws.main_tree, FD.getData(), ws);

            RateModelDisplay       RD = rates.get(current_posteriors_rates);
            StringReader reader = new StringReader(posteriors_sb.toString());
            pd.loadPosteriors(reader, RD.getRateModel());
            FD.addTableModelListener(pd.getFamilyTableScroll());
            ws.data_panel.addItem(pd);

            reading_posteriors = false;
            posteriors_sb = null;
        }
        
        /**
         * End tag for data item
         */
        private void endDataItem() throws SAXException
        {
            WorkSpaceCount ws = sessions.get(current_session_id);

            String id = current_item_attributes.getValue("id");
            String parent = current_item_attributes.getValue("parent");
            String type = current_item_attributes.getValue("type");
            
            FamilySizeTableDisplay F = tables.get(parent);
            ws.data_panel.selectItem(F);
            
            if (DolloDisplay.class.getCanonicalName().equals(type))
            {
                ws.showDollo();
            } else if (WagnerDisplay.class.getCanonicalName().equals(type))
            {
                double gain_penalty = Double.parseDouble(current_item_attributes.getValue("gain"));
                ws.showWagner(gain_penalty);
            } else if (PGLDisplay.class.getCanonicalName().equals(type))
            {
                ws.showPGL();
            } else throw new SAXException("Data item with unrecognized type '"+type+"'");
        }
        
    }
    
    
/*
<!-- XML DTD for Count sessions
	Miklos Csuros csuros@iro.umontreal.ca
	
	Pubid:
	"-//Csuros//DTD Count-session 1.01//EN"
-->

<!ENTITY % ID          "id ID #REQUIRED" >

<!-- top-level element -->
<!ELEMENT count-session (tree data rates)>
*/    
    
    private class TRPanel extends EmbellishedTreePanel
    {
        
        private TRPanel()
        {
            super(main_tree.getRoot(),false);
            setLayout(LayoutStyle.PHENOGRAM);
            super.setPadding(new java.awt.Insets(10,10,10,10*label_font_size+10));
        }

        @Override
        public void computeNodeLabelBoundingBoxes(Graphics g)
        {
            Graphics2D g2 = (Graphics2D)g.create();
            
            FontMetrics leaf_fm = g.getFontMetrics(new Font("Serif", Font.PLAIN, label_font_size));
            FontMetrics inner_fm = g.getFontMetrics(new Font("Serif", Font.ITALIC, label_font_size));
            for (int node_idx=0; node_idx<node.length; node_idx++)
            {
                TreeNode N = node[node_idx];
                
                if (N.isLeaf())
                {
                    Rectangle2D R = DrawString.getBoundingBoxForRotatedString(g2, leaf_fm, N.getName(), 0, -POINT_SIZE/2-3, 0.0, 0.5f);
                    R.add(new Rectangle2D.Double(-POINT_SIZE,-POINT_SIZE,POINT_SIZE*2, POINT_SIZE*2)); // space for the little squares/diamonds
                    setNodeLabelBoundingBox(N, R);
                } else
                {
                    String label = (getMagnification()<1.0?" ":"mmmmm");
                    Rectangle2D R = DrawString.getBoundingBoxForRotatedString(g2, inner_fm, label, POINT_SIZE/2+5, label_font_size/2-2, 0.0, 0.f);
                    R.setRect(R.getX()-1.0, R.getY()-2.0, R.getWidth()+2, R.getHeight()+2);
                    R.add(new Rectangle2D.Double(-POINT_SIZE,-POINT_SIZE,2*POINT_SIZE, 2*POINT_SIZE)); // space for the little squares/diamonds
                    setNodeLabelBoundingBox(N,R);
                }
            }
        }
        
        @Override
        public void computeEdgeLabelBoundingBoxes(Graphics g)
        {
            for (int node_idx=0; node_idx<node.length; node_idx++)
            {
                TreeNode N = node[node_idx];
                
                if (!N.isRoot())
                {
                    setEdgeLabelBoundingBox(N, -1, 3, 1.7*label_font_size);
                }
            }
        }
        
        /**
         * Displays all the node names
         */
        protected void plotNodeNames(Graphics g)
        {
            Font leaf_font = new Font("Serif", Font.PLAIN, label_font_size);
            Font inner_font = new Font("Serif", Font.ITALIC, label_font_size);
            Color label_background = new Color(180,180,180,180);
            Color old_color = g.getColor();

            int inner_node_idx = 1;
            for (int nidx=0; nidx<node.length; nidx++)
            {
                TreeNode N = node[nidx];
                
                int display_idx = getDisplayNodeIndex(N);

                Font label_font=(N.isLeaf()?leaf_font:inner_font);
                FontMetrics label_fm = g.getFontMetrics(label_font);
                
                
                if (N.isLeaf())
                {
                    String label_text = N.getName();
                    int w = label_fm.stringWidth(label_text);
                    //int h = label_font_size+2;
                    
                    int x = (int)displayed_node_location[display_idx].x-w/2;
                    int y = (int)displayed_node_location[display_idx].y-POINT_SIZE/2-3;
                    g.setFont(label_font);
                    g.setColor(old_color);
                    g.drawString(label_text, x, y);              
                    //g.translate((int)displayed_node_location[display_idx].x,(int)displayed_node_location[display_idx].y);
                    //Rectangle2D R = node_label_bounding_box[nidx];
                    //System.out.println("#**WSC.pNN "+nidx+"\t"+R+"\t//"+N);
                    //((Graphics2D)g).draw(R);
                    //g.translate(-(int)displayed_node_location[display_idx].x,-(int)displayed_node_location[display_idx].y);
                } else
                {
                    int h = label_font_size+2;
                    int x = (int)displayed_node_location[display_idx].x+POINT_SIZE/2+5;
                    int y = (int)displayed_node_location[display_idx].y+label_font_size/2-2;
                    String label_text = Integer.toString(inner_node_idx);

                    inner_node_idx++;
                    
                    if (N.getName()!=null && !N.getName().equals(""))
                    {
                        String full_name = N.getName();
                        for (int j=full_name.length(); j>=0; j--)
                        {
                            String dn = full_name.substring(0,j);
                            String cand_label = label_text+" ["+dn;
                            if (j<full_name.length())
                                cand_label = cand_label+"...";
                            cand_label = cand_label+"]";
                            int w = label_fm.stringWidth(cand_label)+8;
                            Rectangle covered_by_label = new Rectangle(x-1, y-h+3,w,h);
                            List<IndexedPoint> nodes_covered = point_set.withinRectangle(covered_by_label);

                            if (j==0 || nodes_covered.size()==0)
                            {
                                label_text = cand_label;
                                break;
                            }
                        }
                    }
                    int w = label_fm.stringWidth(label_text)+2;
                    Rectangle covered_by_label = new Rectangle(x-1, y-h+3,w,h);
                    g.setColor(label_background);//new Color(174,249,63,50)); //new Color(24,24,200,80)); //
                    g.fillRoundRect(covered_by_label.x,covered_by_label.y,covered_by_label.width, covered_by_label.height, 10, 10);
                    g.setFont(label_font);
                    g.setColor(unselected_node_color);
                    g.drawString(label_text, x, y);
                    //g.translate((int)displayed_node_location[display_idx].x,(int)displayed_node_location[display_idx].y);
                    //Rectangle2D R = node_label_bounding_box[nidx];
                    //System.out.println("#**WSC.pNN "+nidx+"\t"+R+"\t//"+N);
                    //((Graphics2D)g).draw(R);
                    //g.translate(-(int)displayed_node_location[display_idx].x,-(int)displayed_node_location[display_idx].y);
                } // not a leaf 
                g.setColor(old_color);
            }
        }
    }    
}
