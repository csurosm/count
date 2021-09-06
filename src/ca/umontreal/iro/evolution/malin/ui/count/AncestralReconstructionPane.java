
package ca.umontreal.iro.evolution.malin.ui.count;


import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Stroke;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;

import java.awt.geom.Rectangle2D;

import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeEvent;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

//import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.SwingWorker;

import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;

import javax.swing.table.AbstractTableModel;
//import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableRowSorter;

import ca.umontreal.iro.evolution.TreeNode;

import ca.umontreal.iro.evolution.genecontent.NodeWithRates;
import ca.umontreal.iro.evolution.genecontent.OccurrenceTable;
import ca.umontreal.iro.evolution.genecontent.TreeWithRates;

import ca.umontreal.iro.evolution.malin.DataFile;
import ca.umontreal.iro.evolution.malin.DoubleRoundedForDisplay;       
import ca.umontreal.iro.evolution.malin.Exportable;
import ca.umontreal.iro.evolution.malin.IndexedPoint;

import ca.umontreal.iro.evolution.malin.ui.ColoredValueRenderer;
import ca.umontreal.iro.evolution.malin.ui.DrawString;
import ca.umontreal.iro.evolution.malin.ui.EmbellishedTreePanel;
import ca.umontreal.iro.evolution.malin.ui.FrozenColumnsTable;
import ca.umontreal.iro.evolution.malin.ui.ZoomableTreePanel;

/**
 * Common methods for displaying the results of ancestral reconstruction
 *
 * @author csuros
 */
public class AncestralReconstructionPane extends JSplitPane 
        implements Exportable,ListSelectionListener, FilterableFamilies, PropertyChangeListener
{
    /**
     * Instantiates a new ancestral reconstruction browser. 
     * This method only sets the basic class variables (the two 
     * input arguments). init() must be called separately by the 
     * subclass in order to set up the whole pane.  
     * 
     * @param main_tree the underlying phylogeny
     * @param data_file the occurrence table structure on which the ancestral reconstruction will be carried out
     * @param W enclosing work space  
     */
    public AncestralReconstructionPane(TreeWithRates main_tree, DataFile<OccurrenceTable> data_file, WorkSpaceCount W)
    {
        this.main_tree = main_tree;
        this.data_file = data_file;
        this.work_space = W;
    }
    
    /*
     * Initialization of the pane. 
     * 
     * The following methods are called (in the 
     * same order).
     * <ol>
     * <li><code>initGeneralVariables()</code>: initia</li>
     * <li><code>initReconstructionVariables()</code></li>
     * <li><code>initSelectionVariables()</code></li>
     * <li><code>initComponents()</code></li>
     * <li><code>computeAll()</code></li>
     * </ol>
     */
    protected void init()
    {
        initGeneralVariables();
        initReconstructionVariables();
        initSelectionVariables();
        initComponents();
        computeAll();
    }
    
    
    // -----------------------------------------------------------------------
    // ------------
    // ------------ General data structure variables
    // ------------
    // -----------------------------------------------------------------------

    
    protected WorkSpaceCount work_space;
    
    /**
     * Phylogeny for the terminal taxa: only the fixed node traversal is used. 
     */
    protected TreeWithRates main_tree;

    /**
     * All the tree nodes in the defined postorder traversal order
     */
    protected NodeWithRates[] tree_nodes;
    
    /**
     * Data file containing family sizes across the terminal taxa.
     */
    protected DataFile<OccurrenceTable> data_file;
    
    /**
     *  Names of families
     */
    protected String[] families;
    
    /**
     *  Number of lineages in which a famiy is present
     */
    protected int[] present_leaves;
    
    /**
     * Number of genes in the family
     */
    protected int[] member_count;
    
    /**
     * Data about family size at the leaves.
     */
    protected int[][] family_size;

    
    // -----------------------------------------------------------------------
    // ------------
    // ------------ Variables for reconstruction
    // ------------
    // -----------------------------------------------------------------------
    
    /**
     * Probability for at least one representative in the family at the nodes
     */
    protected double[][] reconstructed_family_count;
    
    /**
     * Probability for more than one representative in the family at the nodes
     */
    private double[][] reconstructed_family_multi;
    

    // -----------------------------------------------------------------------
    // ------------
    // ------------ Variables for aggregate information in selected families
    // ------------
    // -----------------------------------------------------------------------

    protected double[] selected_family_count;
    protected double[] selected_family_multi;
    protected double[] selected_family_gain;
    protected double[] selected_family_loss;
    protected double[] selected_family_expansion;
    protected double[] selected_family_reduction;
    
    
    /**
     * Stores last row selection in order to skip 
     * updating the selected_[] arrays when it is unnecessary.
     */
    protected int[] last_row_selection;
    
    /**
     * Background thread for computing ancestral reconstruction
     */
    protected SwingWorker<Void,Void> compute_reconstruction_task = null;    
    
    
    // -----------------------------------------------------------------------
    // ------------
    // ------------ Swing components
    // ------------
    // -----------------------------------------------------------------------

    /**
     * Main table for displaying the ancestra reconstruction:
     * rows correspond to families, and columns to 
     * events (gains, losses, presence etc) on nodes or edges. 
     */
    protected JTable family_table;
    /**
     * The enclosing scroller for the family_table.
     */
    protected FrozenColumnsTable<OccurrenceTableModel> family_table_scroll;

    /**
     * Secondary table for displaying information about selected familes. 
     * Every row corresponds to a node/edge; columns correspond to 
     * events  (gains, losses, presence etc).
     */
    protected JTable lineage_table;
    
    /**
     * The enclosing scroller for the lineage_table
     */
    protected JScrollPane lineage_table_scroll;
    
    /**
     * The tree display in the bottom part. 
     */
    protected HistoryTreePanel tree_panel;
    
    /**
     * The enclosing panel for the tree, with the zoom 
     * and other buttons in the bottom bar.
     */
    protected ZoomableTreePanel tree_zoom;

    /**
     * A label about what families are selected: placed in the bottom bar.
     */
    protected JLabel selected_rows_information; 

    /**
     * Progress bar for computing ancestral reconstruction.
     */
    protected JProgressBar computation_progress;    
    
    /**
     * A semaphore variabel for managing the displayed 
     * information about the current family selection. 
     * Families can be selected by usual table selection, 
     * or by a double-click in a cell. In the latter case, 
     * the information can be set by the selectSimilarFamilies()
     * method. The default information is displayed by a ListSelectionListener
     * for the family table. This variable is set to false in 
     * order to disable the update by the selection listener 
     * when the family selection is done by a computation.
     * 
     */
    protected boolean update_selected_rows_information = true;
    
    
    public DataFile<OccurrenceTable> getData()
    {
        return data_file;
    }
    
    public FrozenColumnsTable<OccurrenceTableModel> getFamilyTableScroll()
    {
        return family_table_scroll;
    }
    
    // -----------------------------------------------------------------------
    // ------------
    // ------------ Initialization of data structures
    // ------------
    // -----------------------------------------------------------------------
    
    /**
     * Initializes general information about families: their names and sizes. 
     */
    protected void initGeneralVariables()
    {
        tree_nodes = main_tree.getDFT();

        OccurrenceTable data_table = data_file.getData();
        int num_families = data_table.getNumFamilies();
        families = new String[num_families];
        for (int family_idx=0; family_idx<num_families; family_idx++)
            families[family_idx] = data_table.getFamilyName(family_idx);
        
        int num_leaves = main_tree.getNumLeaves();

        present_leaves=new int[num_families];
        member_count = new int[num_families];
        
        family_size = new int[num_families][];
        
        for(int family_idx=0; family_idx<num_families; family_idx++)
        {
            int[] pattern = data_table.getSizes(family_idx);
            family_size[family_idx] = pattern;
            
            present_leaves[family_idx]=0;
            member_count[family_idx]=0;
            for (int leaf_idx=0; leaf_idx<num_leaves; leaf_idx++)
            {
                if (pattern[leaf_idx]>=0) // skip missing data
                {
                    member_count[family_idx]+=pattern[leaf_idx];
                    if (pattern[leaf_idx]!=0)
                    {
                        present_leaves[family_idx]++;
                    }
                }
            }
        }
    }
    /**
     * Allocates the reconstructed_[] arrays. Subclasses 
     * may use this method to initialize any other
     * data structure that is needed.
     */
    protected void initReconstructionVariables()
    {
        int num_nodes = tree_nodes.length;
        int num_edges = tree_nodes.length-1;
        int num_families = families.length;
                
        reconstructed_family_count = new double[num_families][num_nodes];
        reconstructed_family_multi = new double[num_families][num_nodes];
    }
    
    /**
     * Allocates the selected_[] arrays
     */
    protected void initSelectionVariables()
    {
        int num_nodes = tree_nodes.length;
        int num_edges = tree_nodes.length-1;
                
        selected_family_multi = new double[num_nodes];
        selected_family_count = new double[num_nodes];
        selected_family_gain = new double[num_edges];
        selected_family_loss = new double[num_edges];
        selected_family_expansion = new double[num_edges];
        selected_family_reduction = new double[num_edges];
    }
    
    
    // -----------------------------------------------------------------------
    // ------------
    // ------------ Executing the ancestral reconstruction  
    // ------------
    // -----------------------------------------------------------------------

    protected double getReconstructedFamilyCount(int family_idx, int node_idx)
    {
        if (reconstructed_family_count[family_idx]!=null)
            return reconstructed_family_count[family_idx][node_idx];
        else
            return Double.NaN;
    }
    
    protected double getReconstructedFamilyMulti(int family_idx, int node_idx)
    {
        return reconstructed_family_multi[family_idx][node_idx];
    }

    /**
     * Default implementation just assumes gain at every leaf.
     * 
     * @param family_idx index of gene family
     * @param node_idx index of a non-root node (by TreeWithRates)
     * @return probability that the family has gained a first member on the edge leading to the node
     */
    protected double getReconstructedFamilyGain(int family_idx, int node_idx)
    {
        return getReconstructedFamilyCount(family_idx, node_idx);
    }

    
    /**
     * Default implementation just assumes no loss ever.
     * 
     * @param family_idx index of gene family
     * @param node_idx index of a non-root node (by TreeWithRates)
     * @return probability that the family has lost all its members on the edge leading to the node
     */
    protected double getReconstructedFamilyLoss(int family_idx, int node_idx)
    {
        return 0.0;
    }
    

    /**
     * Default implementation just assumes no expansion ever.
     * 
     * @param family_idx index of gene family
     * @param node_idx index of a non-root node (by TreeWithRates)
     * @return probability that the family has expanded (from size 1) on the edge leading to the node
     */
    protected double getReconstructedFamilyExpansion(int family_idx, int node_idx)
    {
        return 0.0;
    }
    
    /**
     * Default implementation just assumes no reduction ever.
     * 
     * @param family_idx index of gene family
     * @param node_idx index of a non-root node (by TreeWithRates)
     * @return probability that the family has contracted (to size 1) on the edge leading to the node
     */
    protected double getReconstructedFamilyReduction(int family_idx, int node_idx)
    {
        return 0.0;
    }

    
    /**
     * Computes ancestral reconstruction
     * using a background process. Should be 
     * invoked in the event thread only.
     */
    protected void computeAll()
    {
        if (compute_reconstruction_task != null) // a background task is still running 
            compute_reconstruction_task.cancel(false); // no interrupt because we will check isCanceled() anyway...
        initComputeAll();
        
        for (int j=0; j<reconstructed_family_count.length ; j++)
            reconstructed_family_count[j] = null;

        compute_reconstruction_task = newComputingTask();
        // set up the tracking of the progress 
        compute_reconstruction_task.addPropertyChangeListener(this);
        // launch in background and return
        compute_reconstruction_task.execute();        
    }    
    
    protected ComputingTask newComputingTask()
    {
        return new ComputingTask();
    }
    
    /**
     * Displays progress in the progress bar
     * 
     * @param evt better be for the <q>progress</q> property
     */
    @Override
    public void propertyChange(PropertyChangeEvent evt)
    {
        //System.out.println("#*ARP.pC "+evt.getPropertyName()+"\t// evt "+evt);
        
        if ("progress".equals(evt.getPropertyName()))
        {
            int pct = (Integer) evt.getNewValue();
            computation_progress.setIndeterminate(false);
            computation_progress.setValue(pct);
            computation_progress.setString("Reconstruction: "+Integer.toString(pct)+"%");
        }
    }

    protected class ComputingTask extends SwingWorker<Void,Void>
    {
        @Override
        public Void doInBackground() 
        {
            try 
            {
                prepareComputation();
//                System.out.println("#*ARP.CT.dIB prepareComputation() done.");
                for (int j=0; j<reconstructed_family_count.length && !isCancelled(); j++)
                {
                    computeAncestralReconstruction(j);
                    reportFamilyDone(j);
                    //try // slow down for debugging 
                    //{
                    //    Thread.sleep(10L);
                    //} catch (InterruptedException ignore){}
                }
            } catch (Exception E)
            {
                WorkSpaceCount.getDealer(AncestralReconstructionPane.this).exceptionCaught(E, "Computation of ancestral reconstruction failed.", "This is probably due to a programming bug.");
            }


            return null;
        }
        
        // executed on event dispatch thread
        @Override
        public void done()
        {
            finishComputeAll();
            compute_reconstruction_task = null;
        }
        

        /**
         * Prepares the ancestral reconstruction in all the families. 
         * This method is activated in a background process: it should not be 
         * invoked in the event dispatch thread. The default implementation 
         * does nothing: subclasses may override this to reform costly 
         * initializations while the progress bar is set to indeterminate.
         * After this methd is done, progress is set to 0, which is captured 
         * by a property listener to the computation task, and 
         * that is when the progress bar switches to determinate mode. 
         */
        protected void prepareComputation()
        {

        }

        /**
         * Called by the worker after computeAncestralReconstruction().
         * Here, it updates the progress (<q>setProgress</q>).
         * 
         * @param family_idx the family for which everything is computed
         */
        protected void reportFamilyDone(int family_idx)
        {
            int pct = Math.min(100*family_idx/reconstructed_family_count.length,99);
            setProgress(pct);
        }
        
    
        /**
         * Computes the ancestral reconstruction for a family. 
         * The default implementation is for debug purposes only: subclasses 
         * are expected to override this method. In particular, 
         * the given row of the reconstructed_[] arrays needs to be 
         * set. This method is called by a background process only; it should
         * no be invoked n the event dispatch thread. 
         * 
         * @param family_idx 0-based index of the family within the occurrence table
         */
        protected void computeAncestralReconstruction(int family_idx)
        {
            int num_nodes = tree_nodes.length;
            //int num_edges = tree_nodes.length-1;

            reconstructed_family_count[family_idx] = new double[num_nodes];
            reconstructed_family_multi[family_idx] = new double[num_nodes];

            for (int node_idx=0; node_idx<num_nodes; node_idx++)
            {
                NodeWithRates N = tree_nodes[node_idx];
                if (N.isLeaf())
                {
                    int n = family_size[family_idx][node_idx];
                    if (n>0)
                    {
                        reconstructed_family_count[family_idx][node_idx]=1.0;
                        if (n>1)
                            reconstructed_family_multi[family_idx][node_idx]=1.0;
                    }
                }
            }

        }    
    }
    
    

    /**
     * Updates the user interface prior to starting a 
     * computation task in the background. 
     * Here, it updates the progress bar. 
     * Executed in the same thread as computeAll().
     */
    protected void initComputeAll()
    {
        computation_progress.setValue(0);
        computation_progress.setString("Recomputing ..."); // this won't show with Indeterminate ...
        computation_progress.setIndeterminate(true);
        computation_progress.setVisible(true);
    }

    /**
     * Updates the user interface after finishing a 
     * computation task in the background. 
     * Here, it updates the progress bar, fires a 
     * table data update event on the family table, and 
     * resets the original selection.
     * This procedure is executed in the event thread.
     */
    protected void finishComputeAll()
    {
        computation_progress.setString("");
        computation_progress.setIndeterminate(true);
        int[] selected_rows = family_table_scroll.getSelectedModelRows();
        int selected_lineage = lineage_table.getSelectedRow();
        family_table_scroll.getModel().fireTableDataChanged();
        // when it comes back, we redo row selection
        last_row_selection = null; // force recomputing
        family_table_scroll.setSelectedModelRows(selected_rows);
        if (selected_lineage != -1)
            lineage_table.setRowSelectionInterval(selected_lineage,selected_lineage);

        computation_progress.setVisible(false);
    }
    
    
    
    // -----------------------------------------------------------------------
    // ------------
    // ------------ Summing posteriors in the selected rows of the family table
    // ------------
    // -----------------------------------------------------------------------
    
    
    /**
     * Sets the values for the selected...[] arrays 
     * as cumulative for selected columns, but only if this row selection differs from the last one.
     * 
     * @return whether the current row selection differs from the last one used in the computation
     */
     protected boolean validateSelectedRows()
     {
        int[] selected_rows = family_table_scroll.getSelectedModelRows();
        
        if (last_row_selection != null && selected_rows.length == last_row_selection.length)
        {
            boolean same = true;
            for (int i=0; i<selected_rows.length && same; i++)
                if (selected_rows[i]!=last_row_selection[i])
                    same = false;
            if (same)
                return false;
        } 

        sumInRows(selected_rows);
        last_row_selection = selected_rows;

        return true;
    }    
    
    /**
     * Updates the values in selected_... arrays.
     * 
     * @param selected_rows rows in which the summing is carried out
     */
    protected void sumInRows(int[] selected_rows)
    {
        for (int node_idx =0; node_idx<tree_nodes.length; node_idx++)
        {
            selected_family_count[node_idx] = 0.;
            selected_family_multi[node_idx] = 0.;
            for (int j=0; j<selected_rows.length; j++)
            {
                int idx = selected_rows[j];
                if (reconstructed_family_count[idx]!=null)
                {
                    selected_family_count[node_idx] += getReconstructedFamilyCount(idx,node_idx);
                    selected_family_multi[node_idx] += getReconstructedFamilyMulti(idx,node_idx);
                }
            }
        }
        for (int edge_idx=0; edge_idx<tree_nodes.length-1; edge_idx++)
        {
            selected_family_gain[edge_idx] = 0.;
            selected_family_loss[edge_idx] = 0.;
            selected_family_expansion[edge_idx] = 0.;
            selected_family_reduction[edge_idx] = 0.;
            for (int j=0; j<selected_rows.length; j++)
            {
                int idx = selected_rows[j];
                if (reconstructed_family_count[idx]!=null)
                {
                    selected_family_gain[edge_idx] += getReconstructedFamilyGain(idx,edge_idx);
                    selected_family_loss[edge_idx] += getReconstructedFamilyLoss(idx,edge_idx);
                    selected_family_expansion[edge_idx] += getReconstructedFamilyExpansion(idx,edge_idx);
                    selected_family_reduction[edge_idx] += getReconstructedFamilyReduction(idx,edge_idx);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // ------------
    // ------------ Setting up Swing components
    // ------------
    // -----------------------------------------------------------------------
    
    /**
     * Initializes the Swing components in the pane.
     */
    protected void initComponents()
    {
        setBackground(Color.WHITE);
        setDividerLocation(300+getInsets().top);
        setBorder(null);
        
        setResizeWeight(0.5);
        setOrientation(JSplitPane.VERTICAL_SPLIT);
        this.setOneTouchExpandable(true);
        
        initFamilyTable();
        initLineageTable();
        
        JSplitPane table_pane = new JSplitPane();
        table_pane.setBorder(null);
        table_pane.setOrientation(JSplitPane.HORIZONTAL_SPLIT);
        table_pane.setLeftComponent(family_table_scroll);
        table_pane.setRightComponent(lineage_table_scroll);
        table_pane.setResizeWeight(0.5);
        setTopComponent(table_pane);
        
        Font tp_font_rm = new Font("Serif", Font.PLAIN, LookAndFeel.TREE_PANEL_FONT_SIZE);
        Font tp_font_it = tp_font_rm.deriveFont(Font.ITALIC);

        tree_panel = createTreePanel();
        tree_panel.setBackground(getBackground());
            
        computation_progress = new JProgressBar(0,100);
        
        computation_progress.setIndeterminate(false);
        computation_progress.setStringPainted(true);
        computation_progress.setString("");
        computation_progress.setMaximumSize(new Dimension(280,20));
        computation_progress.setFont(tp_font_rm);
        //selected_stats_progress.setEnabled(false);
        computation_progress.setVisible(false);
        computation_progress.setBorderPainted(false);
        
        selected_rows_information = new JLabel(":");
        selected_rows_information.setFont(tp_font_it.deriveFont(LookAndFeel.TABLE_FONT_SIZE*0.8f));
        selected_rows_information.setOpaque(false);
        selected_rows_information.setLabelFor(family_table);
        selected_rows_information.setMaximumSize(new Dimension(520,30));
        selected_rows_information.setMinimumSize(selected_rows_information.getMaximumSize());
        selected_rows_information.setPreferredSize(selected_rows_information.getMaximumSize());
        //selected_rows_information.setBorder(BorderFactory.createEtchedBorder());

        tree_zoom = createZoomableTreePanel();
        setZoomableTreePanel();

        family_table.getSelectionModel().addListSelectionListener(this);
        family_table.setRowSelectionInterval(0,0);      
        
        lineage_table.getModel().addTableModelListener(tree_panel);
        lineage_table.getSelectionModel().addListSelectionListener(tree_panel);
    }
    
    /**
     * This is called from within initComponents() to add the tree in the bottom panel. 
     * Here, it adds tree_zoom as a bottom component in the split pane. 
     */
    protected void setZoomableTreePanel()
    {
        super.setBottomComponent(tree_zoom);
    }
    
    protected HistoryTreePanel createTreePanel()
    {
        return new HistoryTreePanel();        
    }
    
    /**
     * Called from initComponents to obtain the proper 
     * ZoomableTreePanel. Subclasses can use this method to
     * override the default bottom bar setup. 
     * 
     * @return a ZoomableTreePanel for tree_panel 
     */
    protected ZoomableTreePanel createZoomableTreePanel()
    {
        return new ARZoomablePane();
    }
    
    /**
     * Called with the default implementation of the zoomable tree panel.
     * The default implementation puts the label for selected families
     * (selected_row_information) on the left, a progress bar in the middle,
     * and the zoom spinner on the right. 
     * 
     * This method adds the selected_row_information label to the bottom bar. 
     * 
     * @param bb the component for the bottom bar
     */
    protected void createBottomBarLeft(Box bb)
    {
        bb.add(selected_rows_information);
        bb.add(Box.createHorizontalGlue());
    }
    
    /**
     * Called with the default implementation of the zoomable tree panel.
     * The default implementation puts the label for selected families
     * (selected_row_information) on the left, a progress bar in the middle,
     * and the zoom spinner on the right. 
     * 
     * This method adds just a glue to the bottom bar. 
     * 
     * @param bb the component for the bottom bar
     */
    protected void createBottomBarRight(Box bb)
    {
        bb.add(Box.createHorizontalGlue());
    }
    

    /**
     * Called when the family table is initialized.
     * 
     * @return the data model for the famiyl table
     */
    protected OccurrenceTableModel createFamilyTableModel()
    {
        return new OccurrenceTableModel(data_file.getData(), false);
    }
    

    /**
     * Sets up the family_table and family_table_scroll
     * No listeners are attached here.
     */
    private void initFamilyTable()
    {
        OccurrenceTableModel model = createFamilyTableModel();
        
        family_table_scroll = createFamilyTableScroll(model);
        
        family_table = family_table_scroll.getDataTable();
        family_table.setFont(new Font("Serif",Font.PLAIN,LookAndFeel.TABLE_FONT_SIZE));
        
        //int largest_membership = 0;
        //{
        //    OccurrenceTable data_table = data_file.getData();
        //    int num_leaves = main_tree.getNumLeaves();
        //    for(int family_idx=0; family_idx<families.length; family_idx++)
        //    {
        //        int[] pattern = data_table.getSizes(family_idx);
        //        for (int leaf_idx=0; leaf_idx<num_leaves; leaf_idx++)
        //            if (pattern[leaf_idx]>largest_membership)
        //                largest_membership = pattern[leaf_idx];
        //    }
        //}
        
        //DefaultTableCellRenderer family_size_renderer = new DefaultTableCellRenderer()
        //{
        //    @Override
        //    public Component getTableCellRendererComponent(JTable table,
        //                                           Object value,
        //                                           boolean isSelected,
        //                                           boolean hasFocus,
        //                                           int row,
        //                                           int column)  
        //    {
        //        Component retval = super.getTableCellRendererComponent(table, value, false, hasFocus, row, column);
        //        if (row != -1)// && !isSelected)
        //        {
        //            if (value instanceof Number)
        //            {
        //                int x = ((Number)value).intValue();
        //                Color C = LookAndFeel.ABSENCE_COLOR;
        //                if (x==1)
        //                    C=LookAndFeel.SINGLE_PRESENCE_COLOR;
        //                if (x>1)
        //                     C = LookAndFeel.MULTI_PRESENCE_COLOR;
        //                retval.setBackground(C);
        //
        //                if (isSelected)
        //                {
        //                    JLabel label = (JLabel)retval;
        //                    label.setBorder(BorderFactory.createLineBorder(LookAndFeel.MULTI_PRESENCE_COLOR, 1));
        //                }
        //            } else
        //                retval.setBackground(Color.RED); // error!
        //        }
        //        return retval;
        //    }            
        //};
        //for (int j=0; j<main_tree.getNumLeaves(); j++)
        //    family_table.getColumnModel().getColumn(2+j).setCellRenderer(family_size_renderer);

        
        family_table.setDefaultRenderer(DoubleRoundedForDisplay.class, new DoubleRoundedForDisplay.Renderer());
        model.setDefaultRenderers(family_table);
        
        // row header data_table
        JTable row_header = family_table_scroll.getHeaderTable();
        row_header.setFont(family_table.getFont());
        
        family_table_scroll.getViewport().setBackground(getBackground());

        TableRowSorter<OccurrenceTableModel> sorter = new TableRowSorter<OccurrenceTableModel>(family_table_scroll.getModel());
        family_table_scroll.setRowSorter(sorter);
    }
    
    protected FrozenColumnsTable<OccurrenceTableModel> createFamilyTableScroll(OccurrenceTableModel model)
    {
        return new FamilyScroll(model);
    }
        
    /**
     * Produces the lineage table model when the lineage table is initialized. 
     * 
     * @return
     */
    protected LineageTableModel createLineageTableModel()
    {
        return new LineageTableModel();
    }
    
    /**
     * Initializes lineage_table and lineage_table_scroll. Should be called after initFamilyTable(), because
     * parameters of family_table are used here to impose similar appearance.
     * The lineage table's model is attached to the family table as a ListSelectionListener. 
     */
    private void initLineageTable()
    {
        final LineageTableModel table_model = createLineageTableModel();
        
        final TableColumnModel column_model = new DefaultTableColumnModel();
        lineage_table = new JTable(table_model, column_model)
        {
            @Override
            protected JTableHeader createDefaultTableHeader() 
            {
                return new JTableHeader(column_model) 
                {
                    @Override
                    public String getToolTipText(MouseEvent e) 
                    {
                        Point p = e.getPoint();
                        int index = column_model.getColumnIndexAtX(p.x);
                        int realIndex = column_model.getColumn(index).getModelIndex();                        
                        return table_model.getHeaderToolTip(realIndex);
                    }
                };
            }        
        };
        
        lineage_table.setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
        lineage_table.setColumnSelectionAllowed(false);
        lineage_table.setRowSelectionAllowed(true);
        lineage_table.setColumnSelectionAllowed(false);
        lineage_table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        lineage_table.setFont(family_table.getFont());
        lineage_table.createDefaultColumnsFromModel();
        column_model.getColumn(0).setPreferredWidth(180); // first column should be wide (node name)
        for (int j=1; j<table_model.getColumnCount(); j++)
            column_model.getColumn(j).setPreferredWidth(45);

        lineage_table.setDefaultRenderer(DoubleRoundedForDisplay.class, new DoubleRoundedForDisplay.Renderer());       
        
        lineage_table_scroll = new JScrollPane(lineage_table);
        // scroll bars are always shown (MacOS X philosophy)
        lineage_table_scroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        lineage_table_scroll.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS); 
        lineage_table_scroll.getViewport().setBackground(getBackground());
        
        family_table.getSelectionModel().addListSelectionListener(table_model);
    }    
    
    
    /**
     * Constructs a new FamilySizeTableDisplay with the selected families.
     * 
     * @return a FamilySizeTableDisplay for the selected families only
     */
    @Override
    public JComponent newDisplayWithSelectedFamilies()
    {
        WorkSpaceCount ws = (WorkSpaceCount)WorkSpaceCount.getWorkSpace(this);
        FamilySizeTableDisplay td = ws.getSelectedTableDisplay();
        if (td!=null)
        {
            OccurrenceTable data_table = data_file.getData();
            int num_families = data_table.getNumFamilies();
            boolean[] family_is_selected = new boolean[num_families];
            int[] selected_rows = family_table.getSelectedRows();
            for (int j=0; j<selected_rows.length; j++)
            {
                int row_idx = selected_rows[j];
                int model_row_idx = family_table.convertRowIndexToModel(row_idx);
                family_is_selected[model_row_idx] = true;
            }
            OccurrenceTable filtered_table = data_table.tableForFamilies(family_is_selected);

            File filtered_file = new File((File)null, "filt:"+data_file.getFile().getName());
            DataFile<OccurrenceTable> filtered_data = new DataFile<OccurrenceTable>(filtered_table,filtered_file);
            FamilySizeTableDisplay filtered_display = td.newDisplayWithSelectedFamilies(family_is_selected);
            
            //PosteriorDisplay filtered_display = new PosteriorDisplay(filtered_data,rate_file);
            return filtered_display;
        } else
            return null;
    }    
    
    
    // -----------------------------------------------------------------------
    // ------------
    // ------------ User interaction
    // ------------
    // -----------------------------------------------------------------------
    
    
    /**
     * We are listening to selection changes in the family table. 
     * When selection changes, the selection inormation in the bottom 
     * bar is updated. 
     * 
     * @param e a list selection event: the info is updated only if the event is not <q>adjusting</q>
     */
    @Override
    public void valueChanged(ListSelectionEvent e)
    {
        if (!e.getValueIsAdjusting())
        {
            if (update_selected_rows_information)
            {

                int num_selected = family_table.getSelectedRowCount();


                String selection_info_text = "";
                if (num_selected == 0)
                    selection_info_text = "No families selected.";
                else 
                {
                    if (num_selected==1)
                        selection_info_text = "One family selected";
                    else
                        selection_info_text = Integer.toString(num_selected)+" families selected";
                    if (num_selected<TREE_PANEL_MAX_INDIVIDUAL_PLOT)
                    {
                        selection_info_text += " (";
                        int[] selected_rows = family_table_scroll.getSelectedModelRows();
                        for (int j=0; j<selected_rows.length; j++)
                        {
                            if (j!=0) selection_info_text += ", ";
                            selection_info_text += families[selected_rows[j]];
                        }
                        selection_info_text += ")";
                    }
                }

                displaySelectionInfo(selection_info_text);
            }
            else
                update_selected_rows_information = true; // next time
        }    
    }
    
    /**
     * Displays a message about the current family selection. 
     * 
     * @param text what is to be shown 
     */
    private void displaySelectionInfo(String text)
    {
        if (selected_rows_information != null)
        {
            selected_rows_information.setText(text);
        } 
    }

    
    protected javax.swing.JComponent createModelBasedChoicesPane(final int row_idx, final int column_idx)
    {
        String column_name = family_table.getModel().getColumnName(column_idx);
        Object cell_value = family_table.getModel().getValueAt(row_idx, column_idx);
        

        JPanel choices = new JPanel();
        BoxLayout layout = new BoxLayout(choices, BoxLayout.PAGE_AXIS);

        StringBuffer description = new StringBuffer("Multiple selection by column \"");
        description.append(column_name);
        description.append("\" using reference family ");
        description.append(family_table.getModel().getColumnName(0));
        description.append(" (");
        description.append(family_table.getModel().getColumnName(1)); //family name
        description.append(")");
        JLabel whats_happening = new JLabel(description.toString());
        whats_happening.setFont(new Font("Serif", Font.PLAIN, LookAndFeel.TABLE_FONT_SIZE));
        choices.add(whats_happening);
        
        JRadioButton[] choice_buttons = null;
        final ButtonGroup choice_button_group= new ButtonGroup();
        
        if (column_idx > data_file.getData().getKnownPropertiesCount()) // column 0 is id, column 1 ... num_props are properties, others are numeric
        {
            choice_buttons = new JRadioButton[5];
            choice_buttons[0] = new JRadioButton("All families with \""+column_name+"\" < "+cell_value);
            choice_buttons[1] = new JRadioButton("All families with \""+column_name+"\" <= "+cell_value);
            choice_buttons[2] = new JRadioButton("All families with \""+column_name+"\" = "+cell_value);
            choice_buttons[3] = new JRadioButton("All families with \""+column_name+"\" >= "+cell_value);
            choice_buttons[4] = new JRadioButton("All families with \""+column_name+"\" > "+cell_value);
            
            choice_buttons[0].setActionCommand("lt");
            choice_buttons[1].setActionCommand("le");
            choice_buttons[2].setActionCommand("eq");
            choice_buttons[3].setActionCommand("ge");
            choice_buttons[4].setActionCommand("gt");

            for (int i=0; i<choice_buttons.length; i++)
                choice_button_group.add(choice_buttons[i]);
            choice_buttons[2].setSelected(true);
        } else
        {
            choice_buttons = new JRadioButton[2];
            choice_buttons[0] = new JRadioButton("All families with \""+column_name+"\" equal to "+cell_value);
            choice_buttons[1] = new JRadioButton("All families with \""+column_name+"\" containing "+cell_value);
            choice_buttons[0].setActionCommand("eq");
            choice_buttons[1].setActionCommand("contains");
            for (int i=0; i<choice_buttons.length; i++)
                choice_button_group.add(choice_buttons[i]);
            choice_buttons[0].setSelected(true);
        }
        
        
        for (int i=0; i<choice_buttons.length; i++)
            choices.add(choice_buttons[i]);
        
        JPanel OKCancel = new JPanel();
        BoxLayout bottom_layout = new BoxLayout(OKCancel, BoxLayout.LINE_AXIS);
        JButton ok_button = new JButton("Perform selection");
        ok_button.addActionListener(new ActionListener()
        {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                String command = choice_button_group.getSelection().getActionCommand();
                family_table_scroll.selectSimilarFamilies(row_idx, column_idx, command);
            }
        });
        JButton cancel_button = new JButton("Cancel");

        
        OKCancel.add(cancel_button);
        OKCancel.add(ok_button);
        
        return choices;
        
    }
    
    // -----------------------------------------------------------------------
    // ------------
    // ------------ Exportable
    // ------------
    // -----------------------------------------------------------------------
    
    /**
     * Exportable interface: writes the FamilyTable using its model
     * in the default implementation.
     * 
     * @param f File into which the output is written
     * 
     * @throws java.io.IOException when something goes wrong 
     */
    @Override
    public void saveData(File f) throws IOException
    {
        PrintStream PS = new PrintStream(f);
        PS.println(WorkSpaceCount.getDealerCount(this).getStandardHeader(getClass()));
        
        OccurrenceTableModel table_model = family_table_scroll.getModel();
        int num_columns = table_model.getColumnCount();
        
        // header
        PS.print("# ");
        for (int col=0; col<num_columns; col++)
        {
            String col_name = table_model.getColumnName(col);
            if (col!=0)
                PS.print("\t");
            PS.print(col_name);
        }
        PS.println();
        
        int num_rows = table_model.getRowCount();
        for (int row=0; row<num_rows; row++)
        {
            for (int col=0; col<num_columns; col++)
            {
                if (col!=0)
                    PS.print("\t");
                Object V = table_model.getValueAt(row, col);
                if (V instanceof DoubleRoundedForDisplay)
                {
                    DoubleRoundedForDisplay D = (DoubleRoundedForDisplay) V;
                    PS.print(D.doubleValue());
                } else
                    PS.print(V);
            }
            PS.println();
        }
        
        PS.close();
    }
    
    // -----------------------------------------------------------------------
    // ------------
    // ------------ Convenience methods
    // ------------
    // -----------------------------------------------------------------------

    /**
     * Short description of this instance
     * 
     * @return a String that can be used in the browser
     */
    @Override
    public String toString()
    {
        return "Ancestral reconstruction";
    }

    /**
     * Constructs a short name for the node.
     * 
     * @param idx index of tree node (into the traversal array)
     * @return leaf name if leaf, or a positive integer for non-leaf nodes
     */
    protected String getShortNodeName(int idx)
    {
        return LookAndFeel.getShortNodeName(main_tree, idx);
    }
    
    /**
     * Constructs a long name for the node.
     * 
     * @param idx index of tree node (into the traversal array)
     * @return leaf name if leaf, or String of style "int [name]" for non-leaf nodes
     */
    protected String getLongNodeName(int idx)
    {
        return LookAndFeel.getLongNodeName(main_tree, idx);
    }
    

    /**
     * Usually invoked when selection changes in the family table.
     * validateSelectedRows() is called, and a table data change event 
     * is fired, when this is really a new selection.
     */
    public void recomputeLineageTableForCurrentSelection()
    {
        //System.out.println("#*ARP.rLTFCS "+family_table.getSelectedRowCount());
        if (family_table.getSelectedRowCount()!=0)
        {
            if (validateSelectedRows())
            {
                int selected_row = lineage_table.getSelectedRow();
                ((AbstractTableModel)lineage_table.getModel()).fireTableDataChanged();
                if (selected_row != -1)
                    lineage_table.setRowSelectionInterval(selected_row, selected_row);
            }
        }
    }
    
    /**
     * Constants defining the order of columns
     */
    protected static final int LINEAGE_TABLE_COLUMN_FAMILY_COUNT = 1;
    protected static final int LINEAGE_TABLE_COLUMN_FAMILY_MULTI = 2;
    protected static final int LINEAGE_TABLE_COLUMN_FAMILY_GAINS = 3;
    protected static final int LINEAGE_TABLE_COLUMN_FAMILY_LOSSES = 4;
    protected static final int LINEAGE_TABLE_COLUMN_FAMILY_EXPANSIONS = 5;
    protected static final int LINEAGE_TABLE_COLUMN_FAMILY_REDUCTIONS = 6;
    
    /**
     * Model for the lineage history table.
     * 
     * Columns:
     * 0 node
     * 1 family count
     * 2 gene count
     * 3 family gains
     * 4 family losses
     * 5 family expansions
     * 6 family reductions
     */
    protected class LineageTableModel 
            extends AbstractTableModel 
            implements ListSelectionListener
    {
        public LineageTableModel()
        {
            super();
        }

        @Override
        public int getColumnCount()
        {
            return 7;
        }
        
        @Override
        public int getRowCount()
        {
            return tree_nodes.length;
        }
        
        @Override
        public Object getValueAt(int row_idx, int column_idx)
        {
            if (column_idx==0)
                return getLongNodeName(row_idx);
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_COUNT)
                return new DoubleRoundedForDisplay(selected_family_count[row_idx]);
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_MULTI)
                return new DoubleRoundedForDisplay(selected_family_multi[row_idx]);

            if (tree_nodes[row_idx].isRoot())
                return null;

            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_GAINS)
                return new DoubleRoundedForDisplay(selected_family_gain[row_idx]);
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_LOSSES)
                return new DoubleRoundedForDisplay(selected_family_loss[row_idx]);
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_EXPANSIONS)
                return new DoubleRoundedForDisplay(selected_family_expansion[row_idx]);
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_REDUCTIONS)
                return new DoubleRoundedForDisplay(selected_family_reduction[row_idx]);
            
            // should never get here
            return null;
        }

        /**
         * Specify column classes to use the correct comparator in sorting
         * 
         * @param column_idx index of the column
         * @return what class the colum belongs to
         */
        @Override
        public Class getColumnClass(int column_idx)
        {
            if (column_idx == 0) 
                return String.class; 
            else
                return DoubleRoundedForDisplay.class;
        }
        
        @Override
        public String getColumnName(int column_idx)
        {
            if (column_idx==0)
                return "Node";
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_COUNT)
                return "Families";
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_MULTI)
                return ":m (Multi-member)";
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_GAINS)
                return ":g (Gains)";
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_LOSSES)
                return ":l (Losses)";
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_EXPANSIONS)
                return "++ (Expansions)";
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_REDUCTIONS)
                return "-- (Contractions)";
            // should never get here
            return null;
        }
        
        public String getHeaderToolTip(int column_idx)
        {
            if (column_idx==0)
                return "Node name";
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_COUNT)
                return "Number of selected families present at the node";
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_MULTI)
                return "Number of selected families with more than one member at the node";
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_GAINS)
                return "Number of selected families acquired in the lineage leading to the node";
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_LOSSES)
                return "Number of selected families lost in the lineage leading to the node";
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_EXPANSIONS)
                return "Number of selected families expanded (from size 1) in the lineage leading to the node";
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_REDUCTIONS)
                return "Number of selected families contracted (to size 1) in the lineage leading to the node";
            // should never get here
            return null;
            
        }
        
        /**
         * This is called when the selection changes in the family table. 
         * validateSelectedRows() is called, and a table data change event 
         * is fired, when this is really a new selection.
         * 
         * @param e
         */
        @Override
        public void valueChanged(ListSelectionEvent e)
        {
            //System.out.println("#**ARP.LTM.vC selection changed in family_table "+e);
            if (!e.getValueIsAdjusting())
            {
                recomputeLineageTableForCurrentSelection();
            }
        }
        
    }    

    protected static int TREE_PANEL_EDGE_RESOLUTION = 16;
    
    protected static int TREE_PANEL_MAX_INDIVIDUAL_PLOT = 7;
    
    /**
     * Default graphic display of the ancestral reconstruction.
     */
    protected class HistoryTreePanel 
            extends EmbellishedTreePanel
            implements TableModelListener, ListSelectionListener
    {
        protected double getTreePanelPresenceRadius(int label_font_size)
        {
         
            return 1.4*Math.max((double)label_font_size, ((double)label_font_size)*label_font_size/this.normal_label_font_size);
        }
        
        protected double getTreePanelPresenceWidth(int label_font_size)
        {
            double d = getMagnification()*normal_label_font_size;
            return 3.5*d;
        }
        
        protected double getTreePanelSelectedInfoWidth(int label_font_size)
        {
            return 14.*label_font_size;
        }

        //protected static double TREE_PANEL_PRESENCE_RADIUS = 1.4*LookAndFeel.TREE_PANEL_FONT_SIZE;

        protected double getTreePanelGainLossWidth(int label_font_size)
        {
            return 2.5*getTreePanelPresenceRadius(label_font_size);
        }
        protected double getTreePanelGainLossHeight(int label_font_size)
        {
            return 2.0*getTreePanelPresenceRadius(label_font_size);
        }

        protected int getTreePanelPadding(int label_font_size)
        {
            return label_font_size/2;
        }
        
        @Override
        public void setMagnification(double mag)
        {
            int sep = (int)(mag*6.0+0.5);
            if (sep<1) sep=1;
            if (sep>20) sep=20;
            bounding_box_separation = sep; 
            super.setMagnification(mag);
        }
        
        private boolean still_adjusting; //for scroll scynchronization
        
        public HistoryTreePanel()
        {
            this(main_tree.getRoot(),EmbellishedTreePanel.LayoutStyle.PHENOGRAM);
        }
        
        protected HistoryTreePanel(NodeWithRates root, EmbellishedTreePanel.LayoutStyle layout)
        {
            super(root,layout,true,false);
            super.setPadding(new java.awt.Insets(10,10,10,10));
            super.setBoundingBoxSeparation(getTreePanelPadding(LookAndFeel.TREE_PANEL_FONT_SIZE));
            still_adjusting = false;
            this.setCloseRadius(getTreePanelPresenceRadius(LookAndFeel.TREE_PANEL_FONT_SIZE));
            this.setNormalFontSize(LookAndFeel.TREE_PANEL_FONT_SIZE);
        }
        
        /**
         * Called when table model changes (better be the lineage table), which in turn happens when selection changes in history table
         * 
         * @param ignored ignored argument
         */
        @Override
        public void tableChanged(TableModelEvent ignored)
        {
            repaint();
        }
        
        protected void removeTreeSelection()
        {
            super.removeSelection();
        }
        
        @Override
        protected void removeSelection()
        {
            removeTreeSelection();
            lineage_table.clearSelection();
        }

        /**
         * selects a given node and synchronizes table selection if necessary
         * 
         * @param num_mouse_clicks number of clicks (single, double, etc)
         * @param isCtrlDown whether the Ctrl key is pressed at the same time
         */
        @Override
        protected synchronized void selectPoint(IndexedPoint P, int num_mouse_clicks, boolean isCtrlDown)
        {
            if (!still_adjusting)
            {
                still_adjusting = true;
                removeTreeSelection();
                int display_idx = P.getIndex();
                selected[display_idx]=true;
                //super.selectPoint(P, num_mouse_clicks, isCtrlDown);
                TreeNode N= node[display_idx];
                int row_idx = N.getId();
                lineage_table.getSelectionModel().setSelectionInterval(row_idx,row_idx);
                lineage_table.scrollRectToVisible(lineage_table.getCellRect(row_idx,0,true));
                still_adjusting = false;
            }
        }
        
        /**
         * Called when selection changes in the lineage_table
         * 
         * @param E the event that led us here
         */
        @Override
        public void valueChanged(ListSelectionEvent E)
        {
            if (!still_adjusting)
            {
                still_adjusting = true;
                removeTreeSelection();
                int row_idx = lineage_table.getSelectedRow();
                if (row_idx != -1)
                {
                    int point_idx = getDisplayNodeIndex(tree_nodes[row_idx]);
                    
                    selected[point_idx]=true;
                    double nx = displayed_node_location[point_idx].getX();
                    double ny = displayed_node_location[point_idx].getY();
                    
                    double w = Math.max(getTreePanelPresenceWidth(label_font_size), getTreePanelSelectedInfoWidth(label_font_size));
                    
                    Rectangle nRect = new Rectangle((int)(nx-w/2.0),(int)(ny-label_font_size),
                            (int)(w),7*label_font_size);
                    tree_panel.scrollRectToVisible(nRect);
                    //System.out.println("#**PD.PTP.vC row "+row_idx+",\ttn "+tree_nodes[row_idx]+"\tpidx "+point_idx+"\tn "+node[point_idx]+"\trect "+nRect);
                    //selectSubtree(tree_nodes[row_idx]);
                }
                still_adjusting = false;
                repaint();
            }
        }        
        
        
        private float[] dash = {2.0f,4.0f}; 
        private Stroke dashed_stroke  = new BasicStroke(1.1f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 10.0f, dash, 0);
        protected Stroke thick_stroke = new BasicStroke(4.0f);
        protected Stroke medium_stroke = new BasicStroke(2.0f);
        protected Stroke thin_stroke = new BasicStroke(1.0f);
        private float[] dot = {1.0f, 1.0f};
        private Stroke dotted_stroke  = new BasicStroke(1.1f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 10.0f, dot, 0);
        
        private Color rectangle_fill = new Color(220,220,220,220);
        //private Color disc_fill = new Color(255,200,255);//,220);    
        

//        @Override
//        public void paintComponent(Graphics g) 
//        {
//            //computeNodeLabelBoundingBoxes(g);
//            //computeEdgeLabelBoundingBoxes(g);
//            
//            super.paintComponent(g);

//            calculateDisplayedNodeLocations();
//            rebuildPointSet();

//            plotEdges(g);
//            plotNodes(g);
//            plotNodeNames(g);
//        }     
        
        protected Font getNodeCountFont()
        {
            return new Font("Serif",Font.BOLD,label_font_size);
        }
        
        protected Font getLeafLabelFont()
        {
            return new Font("Serif",Font.PLAIN,label_font_size);
        }
        
        protected Font getEdgeLabelFont()
        {
            return new Font("Serif",Font.PLAIN,label_font_size);            
        }
        
        @Override
        public void paintComponent(Graphics g)
        {
            super.paintComponent(g);
            plotSelectedNode(g);
            //g.setColor(Color.MAGENTA);
            //g.drawString(Double.toString(getDisplayTransform().getScaleY()), 100, 100);
        }
        
        /**
         * String that corresponds to a double value
         * 
         * @param d value for displaying
         * @return what should be shown
         */
        private String infoValue(double d)
        {
            if (d==0.0)
                return "-";
            if (d==(int)d)
                return Integer.toString((int)d);
            double rounding = 100000.0;
            d = ((int)(d*rounding+0.5))/rounding;
            return Double.toString(d);
        }
        
        protected void plotSelectedNode(Graphics g)
        {
            Graphics2D g2 = (Graphics2D) g.create();
            int small_font_size = label_font_size*4/5;
            Font small_font = new Font("Serif",Font.PLAIN,small_font_size);
            Font small_title_font = small_font.deriveFont(Font.BOLD);//.deriveFont(Font.ITALIC);

            int selected_idx = lineage_table.getSelectedRow();
            if (selected_idx != -1)
            {
                int[] selected_families = family_table_scroll.getSelectedModelRows();
                int num_selected_families = selected_families.length;
                boolean plot_individually = (num_selected_families<TREE_PANEL_MAX_INDIVIDUAL_PLOT);
                NodeWithRates N = tree_nodes[selected_idx];
                int Nidx = getDisplayNodeIndex(N);
                double nx = displayed_node_location[Nidx].getX();
                double ny = displayed_node_location[Nidx].getY();

                double max_width = getTreePanelSelectedInfoWidth(label_font_size);
                int dw = (int)(max_width+0.5);
                int title_width = g2.getFontMetrics(small_title_font).stringWidth(getLongNodeName(selected_idx)+"  ");
                if (title_width>dw)
                    dw = title_width;

                Color plotbg =  LookAndFeel.TREE_PANEL_SELECTED_NODE_INFO_BACKGROUND;
                Color textc = Color.BLACK;

                int dsep = 5;
                int num_lines = 3; // title, presence, multi
                if (!N.isRoot() || plot_individually)
                    num_lines += 4; // gain, loss, ...
                int dh = num_lines*(small_font_size+1) ;
                int dx = (int)nx;
                int dy = (int)ny;
                if (dx-dw/2<0)
                    dx=dw/2;
                if (dx+dw/2>=getWidth())
                    dx = getWidth()-1-dw/2;
                {   // draw background
                    g2.setColor(plotbg);

                    int[] triangle_x = new int[3];
                    int[] triangle_y = new int[3];
                    triangle_x[0] = (int)nx;
                    triangle_y[0] = (int)ny; // top vertex
                    triangle_x[1] = dx-dw/2+dw; // right vertex
                    triangle_y[1] = dy+dsep; 
                    triangle_x[2] = dx-dw/2; // left vertex
                    triangle_y[2] = triangle_y[1];
                    g2.fillPolygon(triangle_x, triangle_y, 3);
                    g2.fillRect(triangle_x[2], triangle_y[2], dw, dh);
                    g2.setColor(Color.RED);
                    int r = 4;
                    g2.fillOval((int)nx-r, (int)ny-r, 2*r, 2*r);

                }
                {   // draw node name
                    g2.setColor(textc);
                    g2.setFont(small_title_font);
                    DrawString.drawCenteredString(g2, getLongNodeName(selected_idx), dx, dy+dsep+small_font_size);
                }

                int bw = dw/5;
                int bar_left = dx-dw/2+dsep;
                int bar_y =  dy+2*dsep+2*small_font_size;
                
                // draw presence info
                if (plot_individually)
                {
                        bar_y -= small_font_size/2;
                        g2.setStroke(thin_stroke);

                        int individual_padding = 1;
                        int individual_width = (int)(bw/(TREE_PANEL_MAX_INDIVIDUAL_PLOT-1)+0.5)-individual_padding;

                        int total_width = selected_families.length*(individual_padding+individual_width)-individual_padding;
                        int plot_x = dx-dw/2+dsep;
                        int label_x = plot_x + total_width + dsep;

                        g2.setFont(small_font);
                        for (int j=0; j<num_selected_families; j++)
                        {
                            int x = plot_x +j*(individual_padding+individual_width);
                            int family_idx = selected_families[j];
                            g2.setStroke(thin_stroke);
                            plotIndividualPresence(g2, x, bar_y, individual_width, small_font_size, family_idx, selected_idx);
                            g2.setColor(LookAndFeel.MULTI_PRESENCE_COLOR);
                            g2.setStroke(dotted_stroke);
                            int label_y = bar_y+(selected_families.length-1-j)*small_font_size;
                            int from_x = x+individual_width/2;
                            if (j==selected_families.length-1)
                                from_x = x+individual_width;
                            if (j!=selected_families.length-1)
                                g2.drawLine(from_x, bar_y+small_font_size/2, from_x, label_y);
                            g2.drawLine(from_x, label_y, label_x, label_y);
                            g2.setColor(textc);
                            String family_info = families[family_idx];
                            if (reconstructed_family_count[family_idx]!=null)
                            {
                                double cnt = ((int)(100.0*getReconstructedFamilyCount(family_idx, selected_idx)+0.5))/100.0;
                                family_info += ": "+infoValue(cnt);
                                double threshold = 0.05;
                                if (getReconstructedFamilyMulti(family_idx, selected_idx)>threshold)
                                {
                                    family_info += " m";
                                    double m = ((int)(100.0*getReconstructedFamilyMulti(family_idx, selected_idx)+0.5))/100.0;
                                    if (m<1.0)
                                        family_info += "("+infoValue(m)+")";
                                }
                                if (!N.isRoot())
                                {
                                    if (getReconstructedFamilyLoss(family_idx, selected_idx)>threshold)
                                    {
                                        family_info += " loss";
                                        double ls = ((int)(100.0*getReconstructedFamilyLoss(family_idx, selected_idx)+0.5))/100.0;
                                        if (ls<1.)
                                            family_info += "("+infoValue(ls)+")";
                                    }
                                    if (getReconstructedFamilyGain(family_idx, selected_idx)>threshold)
                                    {
                                        family_info += " gain";
                                        double gn = ((int)(100.0*getReconstructedFamilyGain(family_idx, selected_idx)+0.5))/100.0;
                                        if (gn<1.)
                                            family_info += "("+infoValue(gn)+")";
                                        
                                    }
                                    if (getReconstructedFamilyExpansion(family_idx, selected_idx)>threshold)
                                    {
                                        family_info += " ++";
                                        double exp = ((int)(100.0*getReconstructedFamilyExpansion(family_idx, selected_idx)+0.5))/100.0;
                                        if (exp<1.)
                                            family_info += "("+infoValue(exp)+")";
                                    }
                                    if (getReconstructedFamilyReduction(family_idx, selected_idx)>threshold)
                                    {
                                        family_info += " --";
                                        double red = ((int)(100.0*getReconstructedFamilyReduction(family_idx, selected_idx)+0.5))/100.0;
                                        if (red<1.)
                                            family_info += "("+infoValue(red)+")";
                                    }
                                }
                            }
                            g2.drawString(family_info, label_x, label_y+small_font_size/2);
                        }
                } else
                { // multiple selection
                    int mbw = plotPresenceBars(g2,bar_left+bw/2, bar_y, bw, small_font_size/2, null, selected_idx);
                    
                    double threshold = 0.05;

                    int multi_x = bar_left+mbw;
                    int count_x = bar_left+bw;
                    int count_y = bar_y+small_font_size/4;
                    int multi_y = bar_y-small_font_size/4;

                    g2.setColor(textc);
                    g2.setFont(small_font);

                    int label_count_x = count_x+2*dsep;
                    int label_multi_x = label_count_x;
                    int label_count_y = (count_y+multi_y)/2+dsep/2+small_font_size/2;
                    int label_multi_y = label_count_y - small_font_size;

                    g2.setStroke(dotted_stroke);
                    g2.drawLine(count_x, count_y, label_count_x, label_count_y-small_font_size/2);
                    g2.setStroke(medium_stroke);
                    g2.drawString("Present: "+infoValue(selected_family_count[selected_idx]), label_count_x, label_count_y);
                    if (selected_family_multi[selected_idx]>threshold)
                    {

                        g2.setStroke(dotted_stroke);
                        g2.drawLine(multi_x, multi_y, label_multi_x, label_multi_y-small_font_size/2);
                        g2.setStroke(medium_stroke);
                        g2.drawString("Multi: "+infoValue(selected_family_multi[selected_idx]), label_multi_x, label_multi_y);
                    }

                    if (!N.isRoot())
                    {
                        int line_height = small_font_size;
                        int gain_y = label_count_y + line_height;
                        int loss_y = gain_y + line_height;
                        int expansion_y = loss_y + line_height;
                        int reduction_y = expansion_y + line_height;
                        
                        int change_bw = bw/2;
                        int lineage_x = bar_left+change_bw;
                        int change_x = label_count_x;
                        g2.setStroke(thick_stroke);
                        g2.setColor(EDGE_COLOR);
                        g2.drawLine(lineage_x,gain_y-small_font_size/2, lineage_x, reduction_y);
                        g2.setStroke(medium_stroke);
                        if (selected_family_gain[selected_idx]>threshold)
                        {
                            g2.setColor(LookAndFeel.GAIN_COLOR);
                            g2.fillRect(lineage_x, gain_y-3*small_font_size/4, change_bw, small_font_size/2);
                            g2.setColor(textc);
                            g2.drawString("Gains: "+infoValue(selected_family_gain[selected_idx]),change_x, gain_y);
                        }
                        if (selected_family_loss[selected_idx]>threshold)
                        {
                            g2.setColor(LookAndFeel.LOSS_COLOR);
                            g2.fillRect(lineage_x-change_bw, loss_y-3*small_font_size/4, change_bw, small_font_size/2);
                            g2.setColor(textc);
                            g2.drawString("Losses: "+infoValue(selected_family_loss[selected_idx]),change_x, loss_y);
                        }
                        if (selected_family_expansion[selected_idx]>threshold)
                        {
                            g2.setColor(LookAndFeel.EXPANSION_COLOR);
                            g2.fillRect(lineage_x, expansion_y-3*small_font_size/4, change_bw, small_font_size/2);
                            g2.setColor(textc);
                            g2.drawString("Expansions: "+infoValue(selected_family_expansion[selected_idx]),change_x, expansion_y);
                        }
                        if (selected_family_reduction[selected_idx]>threshold)
                        {
                            g2.setColor(LookAndFeel.REDUCTION_COLOR);
                            
                            g2.fillRect(lineage_x-change_bw, reduction_y-3*small_font_size/4, change_bw, small_font_size/2);
                            g2.setColor(textc);
                            g2.drawString("Contractions: "+infoValue(selected_family_reduction[selected_idx]),change_x, reduction_y);
                        }
                    }
                }
            } // if at least one selected 
        }
        
        protected void plotIndividualPresence(Graphics2D g2, int x, int y, int width, int height, int family_idx, int node_idx)
        {
            if (reconstructed_family_count[family_idx]!=null)
            {
                {
                    Color C = g2.getColor();
                    int half_plot_size = height/2;
                    double p = getReconstructedFamilyCount(family_idx,node_idx);
                    double pm = getReconstructedFamilyMulti(family_idx,node_idx);
                    double ps = p-pm;
                    int plotted_width = (int)(width*p+0.5);
                    int plotted_width_multi = (int)(width * pm+0.5);
                    if (plotted_width>0)
                    {
                        g2.setColor(LookAndFeel.SINGLE_PRESENCE_COLOR);
                        //g2.fillRect(x, y+label_font_size/2-hs, width, hs);
                        g2.fillRect(x, y, plotted_width, half_plot_size);
                    }
                    if (plotted_width_multi>0)
                    {
                        g2.setColor(LookAndFeel.SINGLE_PRESENCE_COLOR);
                        g2.fillRect(x, y-half_plot_size, plotted_width_multi, half_plot_size);
                    }
                    g2.setColor(LookAndFeel.MULTI_PRESENCE_COLOR);
                    g2.drawRect(x-1, y-half_plot_size-1, width+1, 2*half_plot_size+1);
                    g2.setColor(C);
                }
                if (false)
                {
                    Color C = g2.getColor();
                    int half_plot_size = label_font_size/2;
                    double p = getReconstructedFamilyCount(family_idx,node_idx);
                    double pm = getReconstructedFamilyMulti(family_idx,node_idx);
                    double ps = p-pm;
                    int plotted_height = (int)(half_plot_size*p+0.5);
                    int plotted_height_single = (int)(half_plot_size * ps +0.5);
                    int plotted_height_multi = plotted_height-plotted_height_single;
                    if (plotted_height_single>0)
                    {
                        g2.setColor(LookAndFeel.SINGLE_PRESENCE_COLOR);
                        //g2.fillRect(x, y+label_font_size/2-hs, width, hs);
                        g2.fillRect(x, y, width, plotted_height_single);
                    }
                    if (plotted_height_multi>0)
                    {
                        g2.setColor(LookAndFeel.MULTI_PRESENCE_COLOR);
                        g2.fillRect(x, y-plotted_height_multi, width, plotted_height_multi);
                    }
                    g2.setColor(LookAndFeel.MULTI_PRESENCE_COLOR);
                    g2.drawRect(x-1, y-half_plot_size-1, width+1, 2*half_plot_size+1);
                    g2.setColor(C);
                }
            }
        }
        
        private double max_present;
        private double max_multi;
        
        private int plotPresenceBars(Graphics2D g2, int bar_x, int bar_y, int bar_width, int bar_height, Color bgcolor, int node_idx)
        {
            int bar_width_multi = (int)(bar_width*max_multi /max_present);
            if (bgcolor != null)
            {
                g2.setColor(bgcolor);
                g2.fillRect(bar_x-bar_width/2, bar_y-bar_height, bar_width_multi, bar_height);
                g2.fillRect(bar_x-bar_width/2, bar_y, bar_width, bar_height);
            }

            int mw = (int)(bar_width*selected_family_multi[node_idx]/max_present+0.5);
            int pw = (int)(bar_width*selected_family_count[node_idx]/max_present+0.5);

            //System.out.println("#*ARP.HTP.pN node "+i+"/"+getShortNodeName(i)+"\tcnt "+selected_family_count[i]+"\tmulti "+selected_family_multi[i]);

            g2.setColor(LookAndFeel.SINGLE_PRESENCE_COLOR);
            if (mw != 0)
                g2.fillRect(bar_x-bar_width/2, bar_y-bar_height, mw, bar_height);
            g2.fillRect(bar_x-bar_width/2, bar_y, pw, bar_height);

            g2.setStroke(thin_stroke);
            g2.setColor(LookAndFeel.MULTI_PRESENCE_COLOR);
            if (bar_width_multi != 0)
            {
                g2.drawLine(bar_x-bar_width/2-1, bar_y-bar_height-1, bar_x-bar_width/2+bar_width_multi, bar_y-bar_height-1);
                g2.drawLine(bar_x-bar_width/2+bar_width_multi, bar_y-bar_height-1,bar_x-bar_width/2+bar_width_multi, bar_y-1);
            }
            g2.drawLine(bar_x-bar_width/2+bar_width_multi, bar_y-1, bar_x-bar_width/2+bar_width, bar_y-1);
            g2.drawLine(bar_x-bar_width/2+bar_width, bar_y-1,bar_x-bar_width/2+bar_width, bar_y+bar_height);
            g2.drawLine(bar_x-bar_width/2+bar_width, bar_y+bar_height,bar_x-bar_width/2-1, bar_y+bar_height);
            g2.drawLine(bar_x-bar_width/2-1, bar_y+bar_height,bar_x-bar_width/2-1, bar_y-1);
            if (bar_width_multi != 0)
                g2.drawLine(bar_x-bar_width/2-1, bar_y-1,bar_x-bar_width/2-1, bar_y-bar_height-1);
            
            return bar_width_multi;
        }
        
        @Override
        protected void plotNodes(Graphics g)
        {
            if (family_table.getSelectedRowCount()==0 )
                return; // nothing to do

            Graphics2D g2 = (Graphics2D) g.create();

            max_present = 1.0;
            max_multi = 0.0;
            for (int i=0; i<tree_nodes.length; i++)
            {
                if (selected_family_count[i]>max_present)
                    max_present = selected_family_count[i];
                if (selected_family_multi[i]>max_multi)
                    max_multi = selected_family_multi[i];
            }
            double max_width = getTreePanelPresenceWidth(label_font_size);
            int bar_width = (int)(max_width+0.5);
            
            g2.setFont(getNodeCountFont());
            
            for (int i=0; i<tree_nodes.length; i++)
            {
                NodeWithRates N = tree_nodes[i];
                int Nidx = getDisplayNodeIndex(N);
                double nx = displayed_node_location[Nidx].getX();
                double ny = displayed_node_location[Nidx].getY();
                
                int bar_x = (int)(nx+0.5); //
                int bar_y = (int)(ny-label_font_size/2+0.5);

                if (compute_reconstruction_task==null || compute_reconstruction_task.isDone())
                {
                    if (family_table.getSelectedRowCount()<TREE_PANEL_MAX_INDIVIDUAL_PLOT)
                    { // presence shown individually
                        g2.setStroke(thin_stroke);
                        int[] selected_families = family_table_scroll.getSelectedModelRows();

                        int individual_padding = 1;
                        int individual_width = (int)(max_width/(TREE_PANEL_MAX_INDIVIDUAL_PLOT-1)+0.5)-individual_padding;

                        int total_width = selected_families.length*(individual_padding+individual_width)-individual_padding;
                        int x_offset = total_width/2;

                        g2.setColor(getBackground());
                        g2.fillRect(bar_x-x_offset, bar_y-label_font_size/2, total_width, label_font_size);

                        for (int j=0; j<selected_families.length; j++)
                        {
                            int x = bar_x - x_offset+j*(individual_padding+individual_width);
                            int family_idx = selected_families[j];
                            plotIndividualPresence(g2, x, bar_y, individual_width, label_font_size, family_idx, i);
                        }
                    } else
                    { // presence shown with horizontal bars
                        plotPresenceBars(g2, bar_x, bar_y, bar_width, label_font_size/2, getBackground(), i);

                        g2.setColor(Color.BLACK);

                        //g2.drawString(Integer.toString(selected_family_count[i]), bar_x+2, bar_y+(int)(label_font_size-0.5)+1);
                        int cnt = (int)(selected_family_count[i]+0.5);

                        DrawString.drawRotatedStringRight(g2, Integer.toString(cnt)//+"/"+Integer.toString(selected_gene_count[i])
                                ,bar_x+bar_width/2,bar_y-3, 0.0);

                    }
                }
                
                if (false){
                    Rectangle2D bb = this.node_label_bounding_box[Nidx];
                    g2.setColor(Color.MAGENTA);
                    g2.translate(nx,ny);
                    g2.draw(bb);
                    g2.translate(-nx,-ny);
                    g2.setColor(Color.BLACK);
                }
            }            
        }                
        
        @Override
        public void computeNodeLabelBoundingBoxes(Graphics g)
        {
            Graphics2D g2 = (Graphics2D)g.create();
            FontMetrics count_fm = g.getFontMetrics(getNodeCountFont());
            double max_width = getTreePanelPresenceWidth(label_font_size);
            
            for (int i=0; i<tree_nodes.length; i++)
            {
                NodeWithRates N = tree_nodes[i]; // traversal by TreeWithRates
                int Nidx = getDisplayNodeIndex(N);
                
                double width  = max_width;
                double height = label_font_size/2;
                Rectangle2D R = new Rectangle2D.Double(-width/2.0-1.0,-2.0*height-1.0,width+1.0,label_font_size+1.0);
                {
                // counts displayed
                    int count = (int)(selected_family_count[i]+0.5);
                    Rectangle2D Rcnt = DrawString.getBoundingBoxForRotatedString(g2, count_fm,Integer.toString(count),(int)(max_width/2.0), -label_font_size/2-3,0.0,1.0f);
                    R.add(Rcnt);
                }
                if (N.isLeaf())
                { // leaf name displayed
                    String label_text = getShortNodeName(i);
                    Rectangle2D Rlab = DrawString.getBoundingBoxForRotatedString(g2,label_text, 0, (int)(-2.0*label_font_size),-Math.PI*0.4,0.f);
                    R.add(Rlab);
                }
                if (N.isRoot())
                {
                    int dw = (int)(getTreePanelSelectedInfoWidth(label_font_size));
                    Rectangle Rselected_node_info = new Rectangle(-dw/2,0,dw,7*(label_font_size*2/3+1));
                    R.add(Rselected_node_info);
                }
                setNodeLabelBoundingBox(N, R);
                
            }
        }
        
        @Override
        public void computeEdgeLabelBoundingBoxes(Graphics g)
        {
            double max_width = getTreePanelPresenceWidth(label_font_size);
            int bar_width = 1+(int)(max_width-0.5);
            int tree_panel_padding = getTreePanelPadding(label_font_size);
            
            for (int node_idx=0; node_idx<node.length; node_idx++)
            {
                NodeWithRates N = tree_nodes[node_idx];
                if (!N.isRoot())
                {
                    Rectangle2D R = new Rectangle2D.Double(-bar_width/2, 0, bar_width, 3*label_font_size/2+2*tree_panel_padding);
                    setEdgeLabelBoundingBox(N, R.getX(), R.getWidth(), R.getHeight());
                }
            }
        }
        
        /**
         * Plots bars for the lineage-specific changes (gains and losses)
         * 
         * @param g2 the graphics context
         * @param shade_edge_presence whether there is a smooth shading for presence changes on the edges 
         */
        protected void plotEdgeChanges(Graphics2D g2, boolean shade_edge_presence)
        {
            g2.setStroke(thin_stroke);
            g2.setFont(getEdgeLabelFont());
            int tree_panel_padding = getTreePanelPadding(label_font_size);
            
            double max_width = getTreePanelPresenceWidth(label_font_size);
            int yoffset_from_child = label_font_size + tree_panel_padding;
            
            
            if (family_table.getSelectedRowCount()!=0)
            {
                double max_presence = 0.0;
                for (int i=0; i<tree_nodes.length; i++)
                {
                    if (selected_family_count[i]>max_presence)
                        max_presence = selected_family_count[i];
                }
                if (max_presence != 0.0 && 
                             (compute_reconstruction_task==null || compute_reconstruction_task.isDone()))
                {
                    if (family_table.getSelectedRowCount()<TREE_PANEL_MAX_INDIVIDUAL_PLOT)
                    {
                        max_presence=TREE_PANEL_MAX_INDIVIDUAL_PLOT-1;
                    }
                    for (int i=0; i<tree_nodes.length-1; i++)
                    {
                        NodeWithRates N = tree_nodes[i];

                        if (!N.isRoot())
                        {
                            int Nidx = getDisplayNodeIndex(N);
                            int Pidx = getDisplayNodeIndex(N.getParent());

                            double nx = displayed_node_location[Nidx].getX();
                            double ny = displayed_node_location[Nidx].getY();
                            double px = displayed_node_location[Pidx].getX();
                            double py = displayed_node_location[Pidx].getY();

                            if (shade_edge_presence)
                            {

                                g2.setStroke(medium_stroke);

                                int Pi = main_tree.getParentIndex(i);
                                //System.out.println("#*ARP.HTP.pE i "+i+"\tcnt "+selected_family_count[i]+"\tprn "+selected_family_count[Pi]+"\tmax "+max_presence);
                                Color Pcolor = ColoredValueRenderer.intermediateColor(
                                        Color.WHITE,
                                        Color.BLACK,
                                        selected_family_count[Pi]/max_presence); // rounding is correct, since max_presence is double
                                Color Ncolor =ColoredValueRenderer.intermediateColor(
                                        Color.WHITE,
                                        Color.BLACK,
                                        selected_family_count[i]/max_presence); 

                                g2.setColor(Pcolor);
                                g2.drawLine((int)px,(int)py,(int)nx,(int)py);

                                for (int j=0; j<TREE_PANEL_EDGE_RESOLUTION; j++)
                                {
                                    int y_prev = (int)(py + (ny-py)*j/((double)TREE_PANEL_EDGE_RESOLUTION));
                                    int y_current = (int)(py + (ny-py)*(j+1.0)/((double)TREE_PANEL_EDGE_RESOLUTION));
                                    Color col = ColoredValueRenderer.intermediateColor(Pcolor,Ncolor,(j+1.0)/TREE_PANEL_EDGE_RESOLUTION);
                                    g2.setColor(col);
                                    g2.drawLine((int)nx,y_prev,(int)nx,y_current);
                                }
                            }
                            
                            g2.setStroke(thin_stroke);
                            
                            //double labelx = nx - (tree_panel_gainloss_width+10.0)/2.;
                            //double labely = ny + tree_panel_presence_radius+tree_panel_padding;
                            
                            double barx = nx;
                            double bary = ny+yoffset_from_child; 
                            int bar_height = label_font_size/2;
                            
                            if (selected_family_loss[i]>0 || selected_family_gain[i]>0)
                            {
                                g2.setStroke(thin_stroke);
                                int zl = (int)(max_width*selected_family_loss[i]/max_presence+0.5);
                                int zg = (int)(max_width*selected_family_gain[i]/max_presence+0.5);
                                if (zg>zl )
                                {
                                    g2.setColor(LookAndFeel.GAIN_COLOR);
                                    int w = zg-zl;
                                    g2.fillRect((int)barx, (int)bary, w, bar_height);
                                } else if (zl>zg)
                                {
                                    g2.setColor(LookAndFeel.LOSS_COLOR);
                                    int w =zl-zg;
                                    g2.fillRect((int)barx-w, (int)bary, w, bar_height);
                                }
                                if (zl>0)
                                {
                                    g2.setColor(LookAndFeel.LOSS_COLOR);
                                    g2.drawRect((int)barx-zl, (int) bary, zl, bar_height);
                                }
                                if (zg>0)
                                {
                                    g2.setColor(LookAndFeel.GAIN_COLOR);
                                    g2.drawRect((int)barx, (int)bary, zg, bar_height);
                                }
                            }
                            
                            if (false){
                                Rectangle2D bb = this.edge_label_bounding_box[Nidx];
                                g2.setColor(Color.RED);
                                //g2.translate(nx,labely);
                                g2.draw(bb);
                                //g2.translate(-nx,-labely);
                                g2.setColor(Color.BLACK);
                            }
                            
                            
                        }
                    }
                } // max_presence > 0
            }            
        }
        
        @Override
        protected void plotEdges(Graphics g)
        {
            Graphics2D g2 = (Graphics2D) g.create();
            super.EDGE_COLOR = new Color(180,180,180); // 240
            g2.setStroke(thick_stroke);
            super.plotEdges(g2); // draws thick edges with gray
            
            plotEdgeChanges(g2, true);
        }         
        
        
        @Override
        protected void plotNodeNames(Graphics g)
        {

            Graphics2D g2 = (Graphics2D) g.create();
            g2.setFont(getLeafLabelFont());
            g2.setColor(rectangle_fill.darker());
            
            double max_width = getTreePanelPresenceWidth(label_font_size);
            int bar_width = 1+(int)(max_width-0.5);

            //int selected_idx = lineage_table.getSelectedRow();

            for (int i=0; i<tree_nodes.length; i++)
            {
                NodeWithRates N = tree_nodes[i];
                int Nidx = getDisplayNodeIndex(N);                    
                double nx = displayed_node_location[Nidx].getX();
                double ny = displayed_node_location[Nidx].getY();
                int bar_x = (int)(nx+0.5); //
                int bar_y = (int)(ny-label_font_size/2+0.5);
                
                if (N.isLeaf())
                {
                    String label_text = getShortNodeName(i);
                    DrawString.drawRotatedStringLeft(g2,label_text, bar_x, bar_y-2*label_font_size,-Math.PI*0.4);
                } 
                //else if (i==selected_idx)
                //{
                //    Color C = g2.getColor();
                //    g2.setColor(Color.RED);
                //    if (!N.isLeaf())
                //        DrawString.drawRotatedStringLeft(g2, getLongNodeName(i), bar_x-bar_width/2-12, bar_y-10, -Math.PI/6.0);
                //    g2.setColor(C);
                //}
            }
        }        
        
    }
        
    protected class ARZoomablePane extends ZoomableTreePanel
    {
        protected ARZoomablePane()
        {
            super(tree_panel);
        }
        
        /**
         * Called with the default implementation of the zoomable tree panel.
         * The default implementation puts the label for selected families
         * (selected_row_information) on the left, a progress bar in the middle,
         * and the zoom spinner on the right. 
         * Subclasses can use the hooks createBottomBarLeft() and createBottomBarRight()
         * to insert on the left or on the right of the progress bar.
         *  
         * @param bb the component for the bottom bar
         */
        @Override
        protected void initBottomBarElements(Box bb)
        {
            createBottomBarLeft(bb);
            
            JPanel progress_panel = new JPanel();
            progress_panel.add(computation_progress);
            progress_panel.setPreferredSize(computation_progress.getPreferredSize());
            progress_panel.setMaximumSize(computation_progress.getMaximumSize());
            //{
            //    Dimension smaller = computation_progress.getMaximumSize();
            //    //java.awt.Insets margin = progress_panel.getInsets();
            //    smaller.setSize(smaller.width-80, smaller.height-6);
            //    computation_progress.setMaximumSize(smaller);
            //    computation_progress.setPreferredSize(smaller); 
            //}
            //progress_panel.setBorder(BorderFactory.createEtchedBorder());
            bb.add(progress_panel);
            
            createBottomBarRight(bb);

            super.initBottomBarElements(bb);
        }
    }
    
    protected class FamilyScroll extends FrozenColumnsTable<OccurrenceTableModel>
    {
        protected FamilyScroll(OccurrenceTableModel model)
        {
            super(model, 2);
        }

        @Override 
        protected int getPreferredColumnWidth(int idx)
        {
            if (idx == 1 ) // family name
                return 120;
            else if (idx == data_file.getData().getKnownPropertiesCount()+3) // profile outline
                return getModel().getPreferredRendererWidth();
            else
                
                return super.getPreferredColumnWidth(idx);
        }

        @Override
        public int selectSimilarFamilies(int family_idx, int col, String command)
        {
            if (command == null)
                return 0;
            else
            {
                int num_selected = super.selectSimilarFamilies(family_idx, col, command);
                String info = "";
                if (num_selected == 1)
                {
                    info = "One row";
                } else
                {
                    info = Integer.toString(num_selected)+" rows";
                }
                info += " selected with "+model.getColumnName(col);
                if (command.equals("eq"))
                    info += "=";
                else if (command.equals("le"))
                    info += "\u2264";
                else if (command.equals("ge"))
                    info += "\u2265";
                info += getModel().getValueAt(family_idx, col);
                displaySelectionInfo(info);
                update_selected_rows_information = false;

                return num_selected;
            }
        }

        @Override
        protected String getHeaderToolTip(int column_idx)
        {
//            System.out.println("#*ARP.FS.gHTT "+column_idx+"\t"+getModel().getColumnName(column_idx));
            return //family_table_scroll.
                    getModel().getColumnHeaderToolTip(column_idx);
        }
        
        @Override
        protected String getCellToolTip(int row, int col)
        {
            return //family_table_scroll.
                    getModel().getCellToolTip(row, col);
        }
        
        
        @Override
        public String getRowName(int family_idx)
        {
            return families[family_idx];
        }
    }
    

}
