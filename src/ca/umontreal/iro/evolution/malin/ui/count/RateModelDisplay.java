
package ca.umontreal.iro.evolution.malin.ui.count;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import java.util.Hashtable;


import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Stroke;
import java.awt.BasicStroke;


import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;

import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeEvent;

import javax.swing.Box;
import javax.swing.Icon;
import javax.swing.JLayeredPane;
import javax.swing.JPanel;
import javax.swing.JCheckBox;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;

import javax.swing.event.ListSelectionListener;

import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.TableColumnModel;

//import ca.umontreal.iro.banality.Verbose;

import ca.umontreal.iro.banality.Functions;
import ca.umontreal.iro.matek.DiscreteDistribution;

import ca.umontreal.iro.evolution.DiscreteGamma;
import ca.umontreal.iro.evolution.TreeNode;

import ca.umontreal.iro.evolution.genecontent.NodeWithRates;
import ca.umontreal.iro.evolution.genecontent.RateVariation;
import ca.umontreal.iro.evolution.genecontent.TreeWithRates;

import ca.umontreal.iro.evolution.malin.DataFile;
import ca.umontreal.iro.evolution.malin.DoubleRoundedForDisplay;
import ca.umontreal.iro.evolution.malin.IndexedPoint;
import ca.umontreal.iro.evolution.malin.Saveable;

import ca.umontreal.iro.evolution.malin.ui.DiscreteDistributionPlot;
import ca.umontreal.iro.evolution.malin.ui.DrawString;
import ca.umontreal.iro.evolution.malin.ui.EmbellishedTreePanel;
import ca.umontreal.iro.evolution.malin.ui.LayeredBorderLayout;
import ca.umontreal.iro.evolution.malin.ui.ZoomableTreePanel;
        
/**
 * Swing component for displaying a RateVariation model.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class RateModelDisplay extends JSplitPane implements Saveable
{
    
    public RateModelDisplay(DataFile<RateVariation>rate_file)
    {
        this.rate_file = rate_file;
        initDataStructures();
        initComponents();
    }

    private JTable table;
    private JScrollPane table_scroll;

    private RateVariationPanel rate_variation_panel;
    private JScrollPane rate_variation_scroll;
    
    private JSplitPane top_pane;

    protected EmbellishedRateTreePanel tree_panel;
    /**
     * The enclosing panel for the tree, with the zoom 
     * and other buttons in the bottom bar.
     */
    protected ZoomableTreePanel tree_zoom;

    
    private JLayeredPane layers;
    //private JPanel tree_layer;
    //private JPanel legend_layer;
    private RateLegend rate_legend;
        
    private DataFile<RateVariation> rate_file;

    private TreeWithRates main_tree;
    private int num_leaves;
    private Hashtable<TreeNode,Integer> tree_node_index;
    private NodeWithRates[] tree_nodes;
    
    public DataFile<RateVariation> getRateModel()
    {
        return rate_file;
    }
    
    public boolean hasAssociatedFile()
    {
        return rate_file.getFile().getParent()!=null;
    }
    
    public void saveData() throws IOException
    {
        saveData(rate_file.getFile());
    }
    
    public void saveData(File f) throws IOException
    {
        PrintStream PS = new PrintStream(f);
        PS.println(WorkSpaceCount.getDealerCount(this).getStandardHeader(getClass()));
        RateVariation rate_model = rate_file.getData();
        String model_description = rate_model.tableRates();
        
        PS.println("# tree "+main_tree.getRoot().newickTree());
        PS.println(model_description);

        if (PS.checkError()) // also flushes
        {
            throw new IOException("Cannot write the table.");

        }
        PS.close();
        
        rate_file.setFile(f);
        setDirty(false);
    }
    
    public boolean isDirty()
    {
        return rate_file.isDirty();
    }
    
    public void setDirty(boolean dirty)
    {
        rate_file.setDirty(dirty);
    }
        
    @Override
    public String toString()
    {
        StringBuffer sb = new StringBuffer();
        sb.append(rate_file.getFile().getName());
        return sb.toString();
    }       
    
    private void initDataStructures()
    {
        RateVariation rate_model = rate_file.getData();
        int num_classes = rate_model.getNumClasses();
        main_tree = rate_model.getMainTree();
        
        NodeWithRates leaf[] = main_tree.getLeaves();
        num_leaves = leaf.length;
        
        tree_nodes = main_tree.getDFT();
        tree_node_index = new Hashtable<TreeNode, Integer>();
        for (int j=0; j<tree_nodes.length; j++)
        {
            tree_node_index.put(tree_nodes[j],new Integer(j));
            //System.out.println("#*RMD.iDS node "+j+"\t"+tree_nodes[j]);
        }
    }

    
    private void initComponents()
    {
        setBackground(Color.WHITE);
        setBorder(null);
        setDividerLocation(300+getInsets().top);
        setResizeWeight(0.5);
        setOrientation(JSplitPane.VERTICAL_SPLIT);
        this.setOneTouchExpandable(true);
        
        table = createTable();
        table_scroll = new JScrollPane(table);
        table_scroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        table_scroll.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        table_scroll.getViewport().setBackground(getBackground());
        
        rate_variation_panel = createRateVariationPanel();
        rate_variation_panel.setBackground(getBackground());
        rate_variation_scroll = new JScrollPane(rate_variation_panel);
        rate_variation_scroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        rate_variation_scroll.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        

        top_pane = new JSplitPane();
        top_pane.setBorder(null);
        top_pane.setOrientation(JSplitPane.HORIZONTAL_SPLIT);
        top_pane.setLeftComponent(table_scroll);
        top_pane.setRightComponent(rate_variation_scroll);
        top_pane.setResizeWeight(0.1);
        top_pane.setOneTouchExpandable(true);
        setTopComponent(top_pane);
        
        tree_panel = createRateTreePanel();
        tree_panel.setBackground(this.getBackground());
        table.getSelectionModel().addListSelectionListener(tree_panel);
        
        tree_zoom = createZoomableTreePanel();
        
        layers = new JLayeredPane();
        LayeredBorderLayout layers_layout = new LayeredBorderLayout();
        layers.setLayout(layers_layout);
        layers_layout.addLayoutComponent(LayeredBorderLayout.CENTER, tree_zoom);
        layers.add(tree_zoom, JLayeredPane.DEFAULT_LAYER);

        rate_legend = createRateLegend();
        layers_layout.addLayoutComponent(LayeredBorderLayout.WEST, rate_legend);
        layers.add(rate_legend, JLayeredPane.PALETTE_LAYER);

        setBottomComponent(layers);       
        
        tree_panel.recomputeTreeLayout();
    }
    

    protected ZoomableTreePanel createZoomableTreePanel()
    {
        return new RateZoomableTreePanel();
    }
    
    private JTable createTable()
    {
        RateTableModel M = new RateTableModel();
        TableColumnModel TCM = new DefaultTableColumnModel();
        
        JTable T = new JTable(M,TCM);
        T.setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
        T.setColumnSelectionAllowed(false);
        T.setRowSelectionAllowed(true);
        T.setColumnSelectionAllowed(false);
        T.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        T.setFont(new Font("Serif",Font.PLAIN,LookAndFeel.TABLE_FONT_SIZE));
        T.createDefaultColumnsFromModel();
        TCM.getColumn(0).setPreferredWidth(180); // first column should be wide (node name)
        T.setDefaultRenderer(DoubleRoundedForDisplay.class, new DoubleRoundedForDisplay.Renderer());
        return T;    
    }
    
    private RateVariationPanel createRateVariationPanel()
    {
        RateVariationPanel P = new RateVariationPanel();
        return P;
    }
    
    protected EmbellishedRateTreePanel createRateTreePanel()
            
    {
        NodeWithRates root = rate_file.getData().getMainTree().getRoot();
        EmbellishedRateTreePanel P = new EmbellishedRateTreePanel(root);
        
        return P;
    }
    
    private RateLegend createRateLegend()
    {
        RateLegend P = new RateLegend();
        //P.setBorder(javax.swing.BorderFactory.createLineBorder(Color.RED));
        return P;
    }
    
    protected void createBottomBarLeft(Box bb)
    {
        Font bb_font = new Font("Serif", Font.PLAIN, LookAndFeel.TREE_PANEL_FONT_SIZE);
        
        JCheckBox legendB = new JCheckBox("Legend");
        legendB.setFont(bb_font);
        legendB.setSelected(true);
        legendB.addItemListener(new ItemListener()
            {
                public void itemStateChanged(ItemEvent E)
                {
                    int ch = E.getStateChange();
                    
                    if (ch == ItemEvent.DESELECTED)
                    {
                        tree_panel.setPadding(new java.awt.Insets(10,10+(int)(3.*tree_panel.getNodeRadius(LookAndFeel.TREE_PANEL_FONT_SIZE)),10,10));
                        rate_legend.setVisible(false);
                    }
                    else if (ch==ItemEvent.SELECTED)
                    {
                        tree_panel.setPadding(new java.awt.Insets(10,10+getPreferredRateLegendDimension().width+(int)(3.*tree_panel.getNodeRadius(LookAndFeel.TREE_PANEL_FONT_SIZE)),10,10));
                        rate_legend.setVisible(true);
                    }
                    tree_panel.repaint();
                }
            });
        bb.add(legendB);
        
        boolean has_duplication = main_tree.hasDuplication();
        boolean has_gain = main_tree.hasGain();

        JCheckBox lossB = new JCheckBox("Loss");
        lossB.setBackground(LookAndFeel.LOSS_COLOR);
        lossB.setFont(bb_font);
        lossB.setSelected(true);
        lossB.addItemListener(new ItemListener()
            {
                public void itemStateChanged(ItemEvent E)
                {
                    int ch = E.getStateChange();
                    if (ch == ItemEvent.DESELECTED)
                        tree_panel.setDrawLoss(false);
                    else if (ch==ItemEvent.SELECTED)
                        tree_panel.setDrawLoss(true);
                    tree_panel.repaint();
                }
            });
        bb.add(lossB);
        if (has_duplication)
        {
            JCheckBox dupB = new JCheckBox("Duplication");
            dupB.setBackground(LookAndFeel.DUPLICATION_COLOR);
            dupB.setFont(bb_font);
            dupB.setSelected(true);
            dupB.addItemListener(new ItemListener()
                {
                    public void itemStateChanged(ItemEvent E)
                    {
                        int ch = E.getStateChange();
                        if (ch == ItemEvent.DESELECTED)
                            tree_panel.setDrawDuplication(false);
                        else if (ch==ItemEvent.SELECTED)
                            tree_panel.setDrawDuplication(true);
                        tree_panel.repaint();
                    }
                });
            bb.add(dupB);
        }

        if (has_gain)
        {
            JCheckBox gainB = new JCheckBox("Gain");
            gainB.setBackground(LookAndFeel.GAIN_COLOR);
            gainB.setFont(bb_font);
            gainB.setSelected(true);
            gainB.addItemListener(new ItemListener()
                {
                    public void itemStateChanged(ItemEvent E)
                    {
                        int ch = E.getStateChange();
                        if (ch == ItemEvent.DESELECTED)
                            tree_panel.setDrawGain(false);
                        else if (ch==ItemEvent.SELECTED)
                            tree_panel.setDrawGain(true);
                        tree_panel.repaint();
                    }
                });
            bb.add(gainB);
        }
    }
    
    public static double roundToMostSignificantDigit(double x, int[] exact_values)
    {
        return roundToMostSignificantDigit(x, exact_values, false);
    }
    /**
     * Rounds a double to one significant digit (called msd). 
     * 
     * @param x the value to be rounded
     * @param exact_values if not null, then (msd, exponent) pair is put there 
     * @param round_to_zero whether we round towards 0 (convrsion to int), or in the usual manner
     * @return the rounded value
     */
    public static double roundToMostSignificantDigit(double x, int[] exact_values, boolean round_to_zero)
    {
        
        if (x==0.0)
        {
            if (exact_values!=null)
            {
                exact_values[0]=0;
                exact_values[1]=0;
            }
            return x;
        } 
        int sign = 1;
        if (x<0.0)
        {
            sign = -1;
            x=-x;
        }
        
        double z = Math.log10(x/(round_to_zero?1.0:0.95));
        int exponent = (int)z;
        if (z<0.) exponent--;
        
        
        double pow10 = Math.pow(10.0, exponent);
        int msd = (int)(x/pow10+(round_to_zero?0.0:0.5));
        msd = sign * msd;
        if (exact_values != null)
        {
            exact_values[0] = msd;
            exact_values[1] = exponent;
        }
        double v = pow10 * msd;
        //System.out.println("#*RMD.rTMSD x "+x+"\t v"+v+"\t"+msd+"E"+exponent);
        
        return v;
    }
    
    
    private class RateTableModel extends AbstractTableModel
    {
        public int getRowCount()
        {
            return main_tree.getNumNodes();
        }
        
        public int getColumnCount()
        {
            return 4;
        }
        
        public Object getValueAt(int row, int column)
        {
            NodeWithRates N = tree_nodes[row];
            if (column == 0)
            {
                return LookAndFeel.getLongNodeName(main_tree, row);
            } else if (column == 1)
            {
                if (N.isRoot())
                    return "";
                else
                {
                    double x = N.getLossRate()*N.getLength();
                    return new DoubleRoundedForDisplay(x);
                }
            } else if (column == 2)
            {
                if (N.isRoot())
                    return "";
                else
                {
                    double x = N.getDuplicationRate()*N.getLength();
                    return new DoubleRoundedForDisplay(x);
                }
            } else if (column == 3)
            {
                if (N.isRoot())
                    return "";
                else
                {
                    double x = N.getTransferRate()*N.getLength();
                    return new DoubleRoundedForDisplay(x);
                }
            } else 
                throw new IllegalArgumentException("RateTableModel.getValueAt --- column index must be 0..3 [got "+column+"]");
            
        }
        
        @Override
        public String getColumnName(int column)
        {
            if (column==0)
                return "Node";
            else if (column==1)
                return "Loss rate";
            else if (column == 2)
                return "Duplication rate";
            else if (column == 3)
                return "Gain rate";
            else 
                throw new IllegalArgumentException("RateTableModel.getColumnName --- column index must be 0..3 [got "+column+"]");
            
        }
        
        @Override
        public Class getColumnClass(int column)
        {
            if (column==0)
                return String.class;
            else return DoubleRoundedForDisplay.class;
        }
                
        
    }
    
    /**
     * Produces a HTML text describing the rate variation model.
     * 
     * @return String of MIME-type <q>text/html</q>
     */
    public String getRateDisplayInfo()
    {
        RateVariation rates = rate_file.getData();
        StringBuffer content = new StringBuffer("<html>\n");
        content.append("<h1>");
        content.append("Birth-death model: rate variation");
        content.append("</h1>\n");
        {
            content.append("<p><strong>Edge length:</strong> ");
            int n = rates.getNumEdgeLengthGammaCategories();
            if (n==1)
            {
                content.append(" identical length");
            } else 
            {
                content.append(Integer.toString(n));
                content.append(" discrete Gamma categories with alpha =");
                content.append(Double.toString(rates.getEdgeLengthAlpha()));
            }
            
            content.append("</p>\n");
        }

        {
            int n = rates.getNumDuplicationRateGammaCategories();
            content.append("<p><strong>Duplication rate:</strong> ");
            if (n==1)
            {
                content.append(" identical rates");
            } else 
            {
                content.append(Integer.toString(n));
                content.append(" discrete Gamma categories with alpha =");
                content.append(Double.toString(rates.getDuplicationRateAlpha()));
            }
            double pz = rates.getDuplicationForbiddenProportion();
            if (pz==0.0)
                content.append(", no forbidden duplications");
            else
            {
                content.append(", ");
                content.append(Double.toString(pz));
                content.append(" proportion of families with forbidden duplications");
            }
            content.append("</p>\n");
        }

        {
            content.append("<p><strong>Loss rate:</strong> ");
            int n = rates.getNumLossRateGammaCategories();

            if (n==1)
            {
                content.append("identical rates");
            } else 
            {
                content.append(Integer.toString(n));
                content.append(" discrete Gamma categories with alpha =");
                content.append(Double.toString(rates.getLossRateAlpha()));
            }
            double pz = rates.getLossForbiddenProportion();
            if (pz==0.0)
                content.append(", no forbidden losses");
            else
            {
                content.append(", ");
                content.append(Double.toString(pz));
                content.append(" proportion of families with forbidden losses");
            }
            content.append("</p>\n");
        }

        {
            content.append("<p><strong>Transfer rate:</strong> ");
            int n = rates.getNumTransferRateGammaCategories();
            if (n==1)
            {
                content.append("identical rates");
            } else 
            {
                content.append(Integer.toString(n));
                content.append(" discrete Gamma categories with alpha =");
                content.append(Double.toString(rates.getTransferRateAlpha()));
            }
            double pz = rates.getTransferForbiddenProportion();
            if (pz==0.0)
                content.append(", no forbidden transfers");
            else
            {
                content.append(", ");
                content.append(Double.toString(pz));
                content.append(" proportion of families with forbidden transfers");
            }
            content.append("</p>\n");
        }
        
        content.append("</html>");
        return content.toString();
    }
    
    
    private class RateVariationPanel extends JPanel
    {
        private RateVariationPanel()
        {
            super();
            initComponents();
        }
        
        private DiscreteGamma D_edge_length;
        private DiscreteGamma D_loss;
        private DiscreteGamma D_duplication;
        private DiscreteGamma D_gain;
        
        private DiscreteGammaPlot P_edge_length;
        private DiscreteGammaPlot P_loss;
        private DiscreteGammaPlot P_duplication;
        private DiscreteGammaPlot P_gain;
                

        private void initComponents()
        {
            setBackground(Color.WHITE);
            //JEditorPane rate_legend_text = new javax.swing.JEditorPane();
            //rate_legend_text.setContentType("text/html");
            //rate_legend_text.setFont(new Font("Serif",Font.PLAIN,LookAndFeel.TABLE_FONT_SIZE));
            //rate_legend_text.setEditable(false);
            //rate_legend_text.setText(getRateDisplayInfo());
            //add(rate_legend_text, BorderLayout.CENTER);
            
        }   
        
        private DiscreteGammaPlot createPlot(DiscreteGamma distribution, int num_categories)
        {
            DiscreteGammaPlot plot = new DiscreteGammaPlot(distribution, num_categories, 
                            RATE_VARIATION_PLOT_WIDTH, RATE_VARIATION_PLOT_HEIGHT);
            plot.setPlottingColor(LookAndFeel.RATE_VARIATION_DISTRIBUTION_COLOR);
            plot.setCategoryColor(LookAndFeel.RATE_VARIATION_CATEGORY_COLOR);
            return plot;
        }
        
        @Override 
        public Dimension getPreferredSize()
        {
            int w = RATE_VARIATION_PLOT_WIDTH + 2*RATE_VARIATION_PADDING + 160;
            int h = LookAndFeel.TABLE_FONT_SIZE*6/5+2*RATE_VARIATION_PADDING+4*(RATE_VARIATION_PADDING+RATE_VARIATION_PLOT_HEIGHT);
            return new Dimension(w,h);
        }
        
        @Override 
        public Dimension getMinimumSize()
        {
            return getPreferredSize();
        }

        @Override
        public void paintComponent(Graphics g)
        {
            RateVariation rates = rate_file.getData();

            Graphics2D g2 = (Graphics2D) g.create();
            int w = getWidth();
            int h = getHeight();
            
            Font normal_font = new Font("Serif",Font.PLAIN,LookAndFeel.TABLE_FONT_SIZE);
            Font title_font  = new Font("Serif", Font.BOLD, LookAndFeel.TABLE_FONT_SIZE*6/5);
            Font subtitle_font = normal_font.deriveFont(Font.BOLD);
            Font it_font = normal_font.deriveFont(Font.ITALIC);
            Font plot_font = new Font("Serif", Font.PLAIN, LookAndFeel.TABLE_FONT_SIZE*4/5);
            
            int current_y = title_font.getSize()+RATE_VARIATION_PADDING;
            g2.setFont(title_font);
            g2.drawString("Rate variation across families", RATE_VARIATION_PADDING, current_y);
            
            int subtitle_length = g2.getFontMetrics(subtitle_font).stringWidth("Duplication rate:");
            
            {
                current_y+=RATE_VARIATION_PADDING;
                g2.setFont(subtitle_font);
                g2.drawString("Edge length:", RATE_VARIATION_PADDING, current_y+subtitle_font.getSize());
                if (P_edge_length == null)
                {
                    int n = rates.getNumEdgeLengthGammaCategories();
                    if (n==1)
                    {
                        g2.setFont(it_font);
                        g2.drawString("no variation", RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING, current_y+subtitle_font.getSize());
                    } else
                    {
                        D_edge_length = new DiscreteGamma(rates.getEdgeLengthAlpha());
                        P_edge_length = new DiscreteGammaPlot(D_edge_length, n);
                        Color c = Color.DARK_GRAY;
                        P_edge_length.setCategoryColor(c);
                        P_edge_length.setPlottingColor(c.brighter().brighter());
                    }
                }
                if (P_edge_length != null)
                {
                    int n = rates.getNumEdgeLengthGammaCategories();
                    if (n==1)
                    {
                        D_edge_length = null;
                        P_edge_length = null;
                    } else
                    {
                        double alpha = rates.getEdgeLengthAlpha();
                        D_edge_length.setAlpha(alpha);
                        P_edge_length.setNumCategories(n);
                        g2.setFont(plot_font);
                        P_edge_length.paintIcon(this, g2, RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING, current_y);
                    }
                }
                current_y += RATE_VARIATION_PLOT_HEIGHT;
            }
            {
                current_y+=RATE_VARIATION_PADDING;
                g2.setFont(subtitle_font);
                g2.setColor(LookAndFeel.LOSS_COLOR);
                g2.drawString("Loss rates:", RATE_VARIATION_PADDING, current_y+subtitle_font.getSize());
                if (P_loss == null)
                {
                    int n = rates.getNumLossRateGammaCategories();
                    if (n==1)
                    {
                        g2.setFont(it_font);
                        g2.drawString("no variation", RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING, current_y+subtitle_font.getSize());
                    } else
                    {
                        D_loss = new DiscreteGamma(rates.getLossRateAlpha());
                        P_loss = new DiscreteGammaPlot(D_loss, n);
                        Color c = LookAndFeel.LOSS_COLOR;
                        P_loss.setCategoryColor(c);
                        P_loss.setPlottingColor(c.brighter().brighter());
                    }
                }
                if (P_loss != null)
                {
                    int n = rates.getNumLossRateGammaCategories();
                    if (n==1)
                    {
                        D_loss = null;
                        P_loss = null;
                    } else
                    {
                        double alpha = rates.getLossRateAlpha();
                        D_loss.setAlpha(alpha);
                        P_loss.setNumCategories(n);
                        g2.setFont(plot_font);
                        P_loss.paintIcon(this, g2, RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING, current_y);
                    }
                }
                current_y += RATE_VARIATION_PLOT_HEIGHT;
            }
            {
                current_y+=RATE_VARIATION_PADDING;
                g2.setFont(subtitle_font);
                g2.setColor(LookAndFeel.DUPLICATION_COLOR);
                g2.drawString("Duplication rates:", RATE_VARIATION_PADDING, current_y+subtitle_font.getSize());
                if (P_duplication == null)
                {
                    int n = rates.getNumDuplicationRateGammaCategories();
                    if (n==1)
                    {
                        double p0 = rates.getDuplicationForbiddenProportion();
                        int txt_x = RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING;
                        int txt_y = current_y+subtitle_font.getSize();
                        if (p0==0.0)
                        {
                            g2.setFont(it_font);
                            g2.drawString("no variation", txt_x, txt_y);
                        } else
                        {
                            g2.setFont(plot_font);
                            g2.drawString("P(no-duplication)="+Double.toString(p0), txt_x, txt_y);
                        }
                    } else
                    {
                        D_duplication = new DiscreteGamma(rates.getDuplicationRateAlpha());
                        P_duplication = new DiscreteGammaPlot(D_duplication, n);
                        Color c = LookAndFeel.DUPLICATION_COLOR;
                        P_duplication.setCategoryColor(c);
                        P_duplication.setPlottingColor(c.brighter().brighter());
                    }
                }
                if (P_duplication != null)
                {
                    int n = rates.getNumDuplicationRateGammaCategories();
                    if (n==1)
                    {
                        D_duplication = null;
                        P_duplication = null;
                        double p0 = rates.getDuplicationForbiddenProportion();
                        if (p0!=0.0)
                        {
                            int txt_x = RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING;
                            int txt_y = current_y+subtitle_font.getSize();
                            g2.setFont(plot_font);
                            g2.drawString("P(no-duplication)="+Double.toString(p0), txt_x, txt_y);
                        }
                    } else
                    {
                        double alpha = rates.getDuplicationRateAlpha();
                        D_duplication.setAlpha(alpha);
                        P_duplication.setNumCategories(n);
                        g2.setFont(plot_font);
                        P_duplication.paintIcon(this, g2, RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING, current_y);
                        double p0 = rates.getDuplicationForbiddenProportion();
                        if (p0!=0.0)
                        {
                            int txt_x = RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING+RATE_VARIATION_PLOT_WIDTH+RATE_VARIATION_PADDING;
                            int txt_y = current_y+subtitle_font.getSize();
                            g2.setFont(plot_font);
                            g2.drawString("P(no-duplication)="+Double.toString(p0), txt_x, txt_y);
                        }
                    }
                }
                current_y += RATE_VARIATION_PLOT_HEIGHT;
            }
            {
                current_y+=RATE_VARIATION_PADDING;
                g2.setFont(subtitle_font);
                g2.setColor(LookAndFeel.GAIN_COLOR);
                g2.drawString("Gain rates:", RATE_VARIATION_PADDING, current_y+subtitle_font.getSize());
                if (P_gain == null)
                {
                    int n = rates.getNumTransferRateGammaCategories();
                    if (n==1)
                    {
                        double p0 = rates.getTransferForbiddenProportion();
                        int txt_x = RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING;
                        int txt_y = current_y+subtitle_font.getSize();
                        if (p0==0.0)
                        {
                            g2.setFont(it_font);
                            g2.drawString("no variation", txt_x, txt_y);
                        } else
                        {
                            g2.setFont(plot_font);
                            g2.drawString("P(no-gain)="+Double.toString(p0), txt_x, txt_y);
                        }
                    } else
                    {
                        D_gain = new DiscreteGamma(rates.getTransferRateAlpha());
                        P_gain = new DiscreteGammaPlot(D_gain, n);
                        Color c = LookAndFeel.GAIN_COLOR;
                        P_gain.setCategoryColor(c);
                        P_gain.setPlottingColor(c.brighter().brighter());
                    }
                }
                if (P_gain != null)
                {
                    int n = rates.getNumTransferRateGammaCategories();
                    if (n==1)
                    {
                        D_gain = null;
                        P_gain = null;
                        double p0 = rates.getTransferForbiddenProportion();
                        if (p0!=0.0)
                        {
                            int txt_x = RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING;
                            int txt_y = current_y+subtitle_font.getSize();
                            g2.setFont(plot_font);
                            g2.drawString("P(no-gain)="+Double.toString(p0), txt_x, txt_y);
                        }
                    } else
                    {
                        double alpha = rates.getTransferRateAlpha();
                        D_gain.setAlpha(alpha);
                        P_gain.setNumCategories(n);
                        g2.setFont(plot_font);
                        P_gain.paintIcon(this, g2, RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING, current_y);
                        double p0 = rates.getTransferForbiddenProportion();
                        if (p0!=0.0)
                        {
                            int txt_x = RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING+RATE_VARIATION_PLOT_WIDTH+RATE_VARIATION_PADDING;
                            int txt_y = current_y+subtitle_font.getSize();
                            g2.setFont(plot_font);
                            g2.drawString("P(no-gain)="+Double.toString(p0), txt_x, txt_y);
                        }
                    }
                }
                current_y += RATE_VARIATION_PLOT_HEIGHT;
            }
        }
    }


    private static double TREE_PANEL_SHORT_EDGE_THRESHOLD = 0.05;
    private static int TREE_PANEL_ROOT_DISTRIBUTION_WIDTH = 120;
    private static int TREE_PANEL_ROOT_DISTRIBUTION_HEIGHT = 40;
    
    private static int RATE_VARIATION_PLOT_WIDTH = 320;
    private static int RATE_VARIATION_PLOT_HEIGHT = 60;
    private static int RATE_VARIATION_PADDING = 5;
    
    /**
     * Given an array of sorted edge lengths x[0]<=x[1]<=...<=x[n-1], 
     * finds an L and the corresponding 
     * quantiles x[0..t-1] and x[n-t..n-1] such that 
     * when the maximum solid edge length is L, 
     * only short edges (i<t) and long edges (i>=n-t)
     * are displayed with a relative solid edge length below epsilon.
     * 
     * @param sorted_values edge lenth values in increasing order
     * @param epsilon (0&lt;epsilon&lt;1 is assumed) the small fractional edge length below which there is no clear visual distinction in the display
     * @return appropriate setting for L, the maximum edge length that is displayed normally (may be 0 in some cases!)
     */
    public static double getDisplayLengthThreshold(double[] sorted_values, double epsilon)
    {
        int n = sorted_values.length;
        double f = epsilon; // epsilon*epsilon / (1.-epsilon);
        double retval = 9e99; // such a weird return value would indicate a problem
        for (int t=1; t-1<=n-t; t++)
        {
            double low = sorted_values[t-1]; // constraint L<= low/epsilon
            double high  = sorted_values[n-t]; // constraint L>=high*epsilon/(1-epsilon)
            if (low >= high*f)
            { // found a good range 
                retval = high; // Math.sqrt(low*high)/Math.sqrt(1.-epsilon);
                break;
            }
        }

        if (retval == 0.0)
        {
            for (int t=0; t<n; t++)
                if (sorted_values[t]!=0.0)
                {
                    retval = sorted_values[(t+n-1)/2];
                    break;
                }
        }

        //retval = 2.0*sorted_values[n/2];

        return retval;
    }
    
    /*
     * Draws the vertical part for an edge: long edges have a proportional 
     * dashed region. The plotted length of the edge is determined by the difference 
     * between the end vertices' coordinates.
     * 
     * @param g2 graphics context: stroke may be altered on return 
     * @param solid_stroke stroke used when edge length is below the maximum display edge length
     * @param dashed_stroke stroke used for part of long edges
     * @param value the edge length to be displayed
     * @param dvalue maximum edge length displayed with solid stroke
     * @param x_i position of the edge
     * @param ny child's vertical coordinate
     * @param ny parent's vertical coordinate
     */ 
    public static void drawEdge(Graphics2D g2, Stroke solid_stroke, Stroke dashed_stroke, double value, double dlength, int x_i, double ny, double py)
    {
        //int x_i = (int)(x+0.5);
        int py_i = (int)(py+0.5);

        if (value>dlength)
        {
            // long edge
            double solid_fraction = dlength/value;
            double my = (ny+(py-ny)*solid_fraction);
            int my_i = (int)(my+0.5);
            int ny_i = (int)(ny+0.5);
            g2.setStroke(solid_stroke);
            g2.drawLine(x_i, ny_i, x_i, my_i);
            g2.setStroke(dashed_stroke);
            g2.drawLine(x_i, my_i, x_i, py_i);
        } else
        {
            double ey = py + (ny-py) * value / dlength;
            int ey_i = (int)(ey+0.5);
            g2.setStroke(solid_stroke);
            g2.drawLine(x_i, py_i, x_i, ey_i);
        }
    }        
    
    
    /**
     * Our own local class for displaying gain/loss/duplication rates on each edge
     */
    protected class EmbellishedRateTreePanel extends EmbellishedTreePanel implements ListSelectionListener
    {
        
        protected EmbellishedRateTreePanel(NodeWithRates root)
        {
            super(root,EmbellishedTreePanel.LayoutStyle.SCALED_PHENOGRAM,true,false);    

            root_distribution = rate_file.getData().getRootPrior();
            
            super.setPadding(new java.awt.Insets(10,10+getPreferredRateLegendDimension().width+(int)(3.*getNodeRadius(LookAndFeel.TREE_PANEL_FONT_SIZE)),10,10));
            super.setRootStemLength(0); // we will have a window instead

            super.EDGE_COLOR = new Color(180,180,180); 
            
            super.setBoundingBoxSeparation(getPadding(LookAndFeel.TREE_PANEL_FONT_SIZE));
            still_adjusting = false;
            dash  = new float[2]; dash[0]=2.0f; dash[1]=4.0f;
            dashed_stroke = new BasicStroke(2.f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 10.0f, dash, 0);
            solid_stroke = new BasicStroke(2.f);
            this.setCloseRadius(getNodeRadius(LookAndFeel.TREE_PANEL_FONT_SIZE));
            this.setNormalFontSize(LookAndFeel.TREE_PANEL_FONT_SIZE);
        }
        
        @Override
        protected void initDataStructures()
        {
            draw_loss = true;
            draw_gain = main_tree.hasGain();
            draw_duplication = main_tree.hasDuplication();
            super.initDataStructures();
        }
                
        
        private DiscreteDistribution root_distribution; 
        
        private boolean still_adjusting;
        private float[] dash;
        private Stroke dashed_stroke;
        private Stroke solid_stroke ;
        
        protected boolean draw_loss;
        protected boolean draw_gain;
        protected boolean draw_duplication;
        
        public void recomputeTreeLayout()
        {
            setupNodeLocations();
            tree_zoom.setValidBoundingBoxes(false);
        }
        
        private void setDrawLoss(boolean b)
        {
            this.draw_loss = b;
            recomputeTreeLayout();
        }
        private void setDrawGain(boolean b)
        {
            this.draw_gain = b;
            recomputeTreeLayout();
        }
        private void setDrawDuplication(boolean b)
        {
            this.draw_duplication = b;
            recomputeTreeLayout();
        }
        
        protected int getPadding(int label_font_size)
        {
            return label_font_size/2;
        }
        
        protected double getNodeRadius(int label_font_size)
        {
         
            return 1.4*Math.max((double)label_font_size, ((double)label_font_size)*label_font_size/this.normal_label_font_size);
        }
        
        protected double getEdgeDisplayHeight(int label_font_size)
        {
            return label_font_size;// * 10.0;
        }

        protected Font getLeafLabelFont()
        {
            return new Font("Serif",Font.PLAIN,label_font_size);
        }
        
        /**
         * We intercept the call to also set bounding box separation appropriately 
         * 
         * @param mag new magnification factor
         */
        @Override
        public void setMagnification(double mag)
        {
            int sep = (int)(mag*6.0+0.5);
            if (sep<1) sep=1;
            if (sep>20) sep=20;
            setBoundingBoxSeparation(sep); 
            super.setMagnification(mag);
            if (rate_legend != null)
                rate_legend.revalidate();
        }
        
        /**
         * Intercepted to calculate proper normalization. 
         */
        @Override
        protected void setupNodeLocations()
        {
            setMaximumDisplayLengths();
            super.setupNodeLocations();
        }
        
        
        protected void setMaximumDisplayLengths()
        {
            int num_edges = main_tree.getNumEdges();
            double[] gain_rate        = new double[num_edges];
            double[] duplication_rate = new double[num_edges];
            double[] loss_rate       = new double[num_edges];
            for (int edge_idx = 0; edge_idx<num_edges; edge_idx++)
            {
                NodeWithRates bottom_node = tree_nodes[edge_idx];
                double len = bottom_node.getLength();
                double gn  = len*bottom_node.getTransferRate();
                double dp  = len*bottom_node.getDuplicationRate();
                double ls  = len*bottom_node.getLossRate();
                gain_rate[edge_idx]        = gn;
                duplication_rate[edge_idx] = dp;
                loss_rate[edge_idx]        = ls;
                //System.out.println("#*RMD.ERTP.sMDL "+edge_idx+"/"+tree_nodes[edge_idx].newickName()+"\tl "+ls+"\tg "+gn+"\td "+dp);
            }
            java.util.Arrays.sort(gain_rate);
            java.util.Arrays.sort(loss_rate);
            java.util.Arrays.sort(duplication_rate);
            
            double epsilon = TREE_PANEL_SHORT_EDGE_THRESHOLD;

            if (draw_gain)
            {
                max_gain_display = getDisplayLengthThreshold(gain_rate, epsilon );
                max_gain_display_exact = new int[2];
                max_gain_display = 2.0*roundToMostSignificantDigit(max_gain_display/2.0, max_gain_display_exact);
                max_gain_display_exact[0]*=2;
                if (max_gain_display_exact[0]>=10)
                {
                    max_gain_display_exact[0]*=0.1;
                    max_gain_display_exact[1]++;
                }
                if (max_gain_display == 0.0)
                    System.out.println("#*RMD.ERTP.sMDL gd "+getDisplayLengthThreshold(gain_rate, epsilon )+"\tmax "+max_gain_display);
            }
            if (draw_duplication)
            {
                max_duplication_display = getDisplayLengthThreshold(duplication_rate, epsilon );
                max_duplication_display_exact = new int[2];
                max_duplication_display = 2.0*roundToMostSignificantDigit(max_duplication_display/2.0, max_duplication_display_exact);
                max_duplication_display_exact[0] *= 2;
                if (max_duplication_display_exact[0]>=10)
                {
                    max_duplication_display_exact[0]*=0.1;
                    max_duplication_display_exact[1]++;
                }
            }
            
            max_loss_display = getDisplayLengthThreshold(loss_rate, epsilon );
            max_loss_display_exact = new int[2];
            max_loss_display = 2.*roundToMostSignificantDigit(max_loss_display/2, max_loss_display_exact);
            max_loss_display_exact[0] *= 2;
            if (max_loss_display_exact[0]>=10)
            {
                max_loss_display_exact[0]*=0.1;
                max_loss_display_exact[1]++;
            }

            //System.out.println("#*RMD.ERTP.sMDL rounded\tmax_l "+max_loss_display+"\tmax_g "+max_gain_display+"\tmax_d "+max_duplication_display);

        }
        
        private double max_gain_display;
        private int[] max_gain_display_exact;
        private double max_duplication_display;
        private int[] max_duplication_display_exact;
        private double max_loss_display;
        private int[] max_loss_display_exact;

        /**
         * Intercept the call to clear the selection in the table
         */
        @Override
        protected void removeSelection()
        {
            super.removeSelection();
            if (!still_adjusting)
                table.getSelectionModel().clearSelection();
        }

        /**
         * Intercept the call to select the corresponding table row and to scroll the table to it
         * @param P the newly selected point (corresponding to a tree node)
         * @param num_mouse_clicks number of clicks on the point (ignored)
         * @param ignored (whether Ctrl was pressed at the same time)
         */
        @Override
        protected void selectPoint(IndexedPoint P, int num_mouse_clicks, boolean ignored)
        {
            if (!still_adjusting)
            {
                still_adjusting = true;
                super.selectPoint(P, num_mouse_clicks, false);
                int display_idx = P.getIndex();
                TreeNode N= node[display_idx];
                int row_idx = ((Integer)tree_node_index.get(N)).intValue();
                table.getSelectionModel().setSelectionInterval(row_idx,row_idx);
                table.scrollRectToVisible(table.getCellRect(row_idx,0,true));
                still_adjusting = false;
            }
        }

        /**
         * When selection changes in the table, we scroll to the corresponding 
         * node in the tree.
         * 
         * @param E the same action is performed for all values of E
         */
        public void valueChanged(javax.swing.event.ListSelectionEvent E)
        {
            if (!still_adjusting)
            {
                still_adjusting = true;
                removeSelection();
                int row_idx = table.getSelectedRow();
                if (row_idx != -1)
                {
                    int point_idx = getDisplayNodeIndex(tree_nodes[row_idx]);
                    selected[point_idx]=true;
                    double nx = displayed_node_location[point_idx].getX();
                    double ny = displayed_node_location[point_idx].getY();
                    int r = (int)(getNodeRadius(LookAndFeel.TREE_PANEL_FONT_SIZE)+0.5);
                    Rectangle nRect = new Rectangle((int)(nx-r),(int)(ny-r),2*r,7*r);
                    tree_panel.scrollRectToVisible(nRect);
                }
                repaint();
                still_adjusting = false;
            }
        }
        
        @Override
        public void paintComponent(Graphics g)
        {
            super.paintComponent(g);
            plotSelectedNode(g);
            //g.setColor(Color.MAGENTA);
            //g.drawString(Double.toString(getDisplayTransform().getScaleY()), 100, 100);
        }
        
        protected void plotSelectedNode(Graphics g)
        {
            Graphics2D g2 = (Graphics2D) g.create();
            int small_font_size = label_font_size*4/5;
            Font small_font = new Font("Serif",Font.PLAIN,small_font_size);
            Font small_title_font = small_font.deriveFont(Font.BOLD);//.deriveFont(Font.ITALIC);

            int selected_idx = table.getSelectedRow();
            if (selected_idx != -1)
            {
                
                NodeWithRates N = tree_nodes[selected_idx];
                int Nidx = getDisplayNodeIndex(N);
                double nx = displayed_node_location[Nidx].getX();
                double ny = displayed_node_location[Nidx].getY();
                
                //super.plotNode(g2, N);
                g2.setColor(Color.RED);
                g2.setStroke(solid_stroke);
                double r = (int)(getNodeRadius(label_font_size)+0.5);
                int d = (int)(2.*r+0.5);
                
                if (N.isRoot())
                {
                    g2.drawOval((int)(nx-r),(int)(ny-r),d,d);
                } else
                {
                    int dh = (int)(1.5*d);
                    int dw = dh*5/2;
                    int dsep = 5;

                    int title_width = g2.getFontMetrics(small_title_font).stringWidth(LookAndFeel.getLongNodeName(main_tree,selected_idx)+"  ");
                    if (title_width>dw)
                        dw = title_width;
                    
                    int dx = (int)nx;
                    int dy = (int)ny;
                    if (N.isLeaf() || true) // always put it below the node
                        dy += r+3;
                    else
                        dy -= r+2*dh+3*dsep;
                    
                    if (dy<0) dy=0;
                    if (dx<dw/2+dsep) dx=dw/2+dsep;
                    if (dy+2*dh+3*dsep>getHeight())
                        dy=getHeight()-2*dh-3*dsep;
                    if (dx+dw/2+dsep>getWidth())
                        dx = getWidth()-dw/2-dsep;
                    

                    DiscreteDistribution pdistr = N.getDuplicationDistribution();
                    DiscreteDistribution xdistr = N.getTransferDistribution();
                    DiscreteDistributionPlot pplot = new DiscreteDistributionPlot(pdistr, dw, dh);
                    Color plotbg =  LookAndFeel.TREE_PANEL_SELECTED_NODE_INFO_BACKGROUND;
                    Color plotc = Color.RED;
                    Color legendc = Color.BLACK;
                    g2.setColor(plotbg);
                    
                    int[] triangle_x = new int[3];
                    int[] triangle_y = new int[3];
                    triangle_x[0] = (int)nx;
                    triangle_y[0] = (int)ny;
                    triangle_x[1] = dx-dw/2+dw+dsep;
                    triangle_y[1] = dy-dsep;
                    triangle_x[2] = dx-dw/2;
                    triangle_y[2] = triangle_y[1];
                    g2.fillPolygon(triangle_x, triangle_y, 3);

                    g2.fillRect(dx-dw/2-dsep, dy-dsep, dw+2*dsep, 2*dh+3*dsep);
                    g2.setColor(Color.RED);
                    int pici_r = 4;
                    g2.fillOval((int)nx-pici_r, (int)ny-pici_r, 2*pici_r, 2*pici_r);
                    
                    {   // draw node name
                        g2.setColor(legendc);
                        g2.setFont(small_title_font);
                        DrawString.drawCenteredString(g2, LookAndFeel.getLongNodeName(main_tree, selected_idx), dx, dy);
                    }

                    pplot.setBackground(null);
                    pplot.setLegendColor(legendc);
                    pplot.setPlottingColor(plotc);
                    g2.setFont(small_font);
                    pplot.paintIcon(this, g2, dx-dw/2, dy);
                    g2.setColor(legendc);
                    g2.setFont(small_font.deriveFont(Font.BOLD));
                    g2.drawString("Inparalogs", dx-dw/2+6, dy+small_font.getSize());
                    DiscreteDistributionPlot xplot = new DiscreteDistributionPlot(xdistr, dw, dh);
                    xplot.setBackground(null);
                    xplot.setLegendColor(legendc);
                    xplot.setPlottingColor(plotc);
                    g2.setFont(small_font);
                    xplot.paintIcon(this, g2, dx-dw/2, dy+dh+dsep);
                    g2.setColor(legendc);
                    g2.setFont(small_font.deriveFont(Font.BOLD));
                    g2.drawString("Xenologs", dx-dw/2+6, dy+dh+dsep+small_font.getSize());
                }
            }
        }
        
        @Override
        protected void plotNodes(Graphics g)
        {
            
            
            Graphics2D g2 = (Graphics2D) g.create();
            Font small_font = new Font("Serif",Font.PLAIN,LookAndFeel.TREE_PANEL_FONT_SIZE*4/5);
            if (root_distribution != null)
            {
                NodeWithRates root = main_tree.getRoot();
                int root_idx = getDisplayNodeIndex(root);
                int nx = (int)(displayed_node_location[root_idx].getX()+0.5);
                int ny = (int)(displayed_node_location[root_idx].getY()+0.5);
                
                DiscreteDistributionPlot root_plot = new DiscreteDistributionPlot(root_distribution, TREE_PANEL_ROOT_DISTRIBUTION_WIDTH, TREE_PANEL_ROOT_DISTRIBUTION_HEIGHT);
                root_plot.setBackground(LookAndFeel.DISTRIBUTION_PLOT_BACKGROUND);
                root_plot.setLegendColor(Color.BLACK);
                root_plot.setPlottingColor(LookAndFeel.MULTI_PRESENCE_COLOR);
                //root_plot.setRange(1, 10);
                g2.setFont(small_font);
                
                root_plot.paintIcon(this, g2, nx-TREE_PANEL_ROOT_DISTRIBUTION_WIDTH/2, ny+LookAndFeel.TREE_PANEL_FONT_SIZE);
            }
            
            
            if (false){
                for (int i=0; i<tree_nodes.length; i++)
                {
                    NodeWithRates N = tree_nodes[i];
                    int Nidx = getDisplayNodeIndex(N);
                    double nx = displayed_node_location[Nidx].getX();
                    double ny = displayed_node_location[Nidx].getY();
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
        protected void plotNodeNames(Graphics g)
        {
            Graphics2D g2 = (Graphics2D) g.create();
            Font label_font = getLeafLabelFont();
            
            for (int node_idx=0; node_idx<node.length; node_idx++)
            {
                NodeWithRates N = tree_nodes[node_idx];
                
                if (N.isLeaf())
                {
                    // our own display for leaves

                    int Nidx = getDisplayNodeIndex(N);                    
                    double nx = displayed_node_location[Nidx].getX();
                    double ny = displayed_node_location[Nidx].getY();
                    
                    int label_x = (int)(nx+0.5);
                    int label_y = (int) (ny+0.5);
                    
                    g2.setFont(label_font);
                    g2.setColor(Color.DARK_GRAY);
                    String label_text = LookAndFeel.getShortNodeName(main_tree,node_idx);
                    DrawString.drawRotatedStringLeft(g2,label_text, label_x, label_y,-Math.PI*0.4);                    
                } else
                    super.plotNodeLabel(g,N,LookAndFeel.getShortNodeName(main_tree,node_idx));        
            }
        }
        
        @Override
        public void computeNodeLabelBoundingBoxes(Graphics g)
        {
            Graphics2D g2 = (Graphics2D)g.create();
            double r = getNodeRadius(label_font_size);
            
            
            for (int i=0; i<tree_nodes.length; i++)
            {
                NodeWithRates N = tree_nodes[i]; // traversal by TreeWithRates
                int Nidx = getDisplayNodeIndex(N);

                Rectangle2D R = new Rectangle2D.Double(-r/2,-r/2,r+1, r+1);
                
                
                if (N.isLeaf())
                { // leaf name displayed
                    String label_text = LookAndFeel.getShortNodeName(main_tree, i);
                    Rectangle2D Rlab = DrawString.getBoundingBoxForRotatedString(g2,label_text, 0, 0,-Math.PI*0.4,0.f);
                    R.add(Rlab);
                    
                }
                if (N.isRoot())
                {
                    Rectangle2D Rplot = new Rectangle2D.Double(-TREE_PANEL_ROOT_DISTRIBUTION_WIDTH/2, LookAndFeel.TREE_PANEL_FONT_SIZE, TREE_PANEL_ROOT_DISTRIBUTION_WIDTH, TREE_PANEL_ROOT_DISTRIBUTION_HEIGHT);
                    R.add(Rplot);
                }
                setNodeLabelBoundingBox(N, R);
            }
        }        
        
        @Override
        public void computeEdgeLabelBoundingBoxes(Graphics g)
        {
            super.computeEdgeLabelBoundingBoxes(g);
            //int hgt = (int)(getEdgeDisplayHeight(label_font_size)+0.5);
            //
            //for (int node_idx=0; node_idx<node.length; node_idx++)
            //{
            //    NodeWithRates N = tree_nodes[node_idx];
            //    if (!N.isRoot())
            //    {
            //        setEdgeLabelBoundingBox(N, -label_font_size, 2*label_font_size, hgt);
            //    }
            //}
        }
        
        private void drawEdge(Graphics2D g2, double value, double dlength, int x_i, double ny, double py)
        {
            RateModelDisplay.drawEdge(g2, solid_stroke, dashed_stroke, value, dlength, x_i, ny, py);
        }

        @Override
        public void plotEdges(Graphics g)
        {
            Graphics2D g2 = (Graphics2D) g.create();
            
            boolean has_duplication = main_tree.hasDuplication();
            boolean has_gain = main_tree.hasGain();
            
            for (int node_idx=0; node_idx<tree_nodes.length; node_idx++)
            {
                NodeWithRates N = tree_nodes[node_idx];
                
                if (N.isRoot()) // just a little stem
                {
                    int Nidx = getDisplayNodeIndex(N);
                    int nx = (int)(displayed_node_location[Nidx].getX()+0.5);
                    int ny = (int)(displayed_node_location[Nidx].getY()+0.5);
                    g2.setStroke(solid_stroke);
                    g2.setColor(EDGE_COLOR); 
                    g2.drawLine(nx,ny,nx,ny+LookAndFeel.TREE_PANEL_FONT_SIZE);
                }
                else
                {
                    int Nidx = getDisplayNodeIndex(N);
                    int Pidx = getDisplayNodeIndex(N.getParent());
                    
                    int nx = (int)(displayed_node_location[Nidx].getX()+0.5);
                    double ny = displayed_node_location[Nidx].getY();
                    int px = (int)(displayed_node_location[Pidx].getX()+0.5);
                    double py = displayed_node_location[Pidx].getY();
                    int py_i = (int)(py+0.5);
                    
                    int num_bars=0;
                    if (draw_loss) num_bars++;
                    if (draw_gain) num_bars++;
                    if (draw_duplication) num_bars++;
                   
                    int bar_x1, bar_x2, bar_x3;
                    if (num_bars==3)
                    {
                        bar_x1=nx-4;
                        bar_x2=nx;
                        bar_x3=nx+4;
                    } else if (num_bars==2) 
                    {
                        bar_x1 = nx-2;
                        bar_x2 = bar_x3 = nx+2;
                    } else 
                    {
                        bar_x1=bar_x2=bar_x3=nx;
                    }
                    
                    //int loss_x = nx - 4;
                    //int gain_x = nx + 4;
                    //int dup_x = nx;
                            
                    g2.setStroke(solid_stroke);
                    g2.setColor(EDGE_COLOR);
                    
                    if (px<nx)
                    {
                        g2.drawLine(px, py_i, bar_x3, py_i);
                    } else
                    {
                        g2.drawLine(px, py_i, bar_x1, py_i);
                    }
                    double length = N.getLength();

                    double loss = N.getLossRate()*length;
                    double gain = N.getTransferRate()*length;
                    double dup = N.getDuplicationRate()*length;
                    
                    if (max_loss_display != 0.0)
                        loss /= max_loss_display;
                    if (draw_gain && max_gain_display != 0.0)
                        gain /= max_gain_display;
                    if (draw_duplication && max_duplication_display != 0.0)
                        dup /= max_duplication_display;
                    
                    double dlength = getDisplayEdgeLength(N);
                    
                    if (num_bars == 0)
                    {
                        g2.setColor(EDGE_COLOR);
                        drawEdge(g2, dlength, dlength, nx, ny, py);
                    } else
                    {
                        int bar_idx = 0;
                        if (draw_loss)
                        {
                            bar_idx++;
                            g2.setColor(LookAndFeel.LOSS_COLOR);
                            drawEdge(g2, loss, dlength, bar_x1, ny, py);
                            //System.out.println("#*RMD.ERTP.pE "+node_idx+"\tloss "+loss+"\tdl "+dlength+"\tny "+ny+"\tpy "+py+"\t// node "+node_location[Nidx]+"\tparent "+node_location[Pidx]);
                        }
                        if (draw_gain)
                        {
                            g2.setColor(LookAndFeel.GAIN_COLOR);
                            drawEdge(g2, gain, dlength, bar_x3, ny, py);
                            //System.out.println("#*RMD.ERTP.pE "+node_idx+"\tgain "+gain+"\tdl "+dlength+"\tny "+ny+"\tpy "+py+"\t// node "+node_location[Nidx]+"\tparent "+node_location[Pidx]);
                        }
                        if (draw_duplication)
                        {
                            int bar_x = (bar_idx==0?bar_x1:bar_x2);
                            bar_idx++;
                            g2.setColor(LookAndFeel.DUPLICATION_COLOR);
                            drawEdge(g2, dup, dlength, bar_x, ny, py);
                            //System.out.println("#*RMD.ERTP.pE "+node_idx+"\tdup  "+dup+"\tdl "+dlength+"\tny "+ny+"\tpy "+py+"\t// node "+node_location[Nidx]+"\tparent "+node_location[Pidx]);
                        }
                    }
                    
                }
            }
        }                

        
        /**
         * Length of the plot for each edge is at most for edge length 1.0  
         * 
         * @param N child node
         * @return maximum edge length or this edge length, whichever is smaller
         */
        @Override
        protected double getDisplayEdgeLength(TreeNode N)
        {
            double lN = getNormalizedEdgeLength((NodeWithRates)N);
            double retval = Math.max(TREE_PANEL_SHORT_EDGE_THRESHOLD, Math.min(lN, 1.0));
            
            //System.out.println("#*RMD.gDEL length "+retval+"\tnode "+N);
            return retval;
        }
        
        /**
         * Computes the normalized edge length: maximum of rate/mr, where rate=loss, gain, duplication and mr 
         * is the corresponding max_rate_display
         * 
         * @param N one of the tree nodes
         * @return normalized edge length
         */
        protected double getNormalizedEdgeLength(NodeWithRates N)
        {
            double length = N.getLength();
            double loss = length*N.getLossRate();
            double duplication = length*N.getDuplicationRate();
            double gain = length*N.getTransferRate();
            
            
            if (max_loss_display!=0.)
                loss /= max_loss_display;
            if (duplication != 0. && max_duplication_display!=0.)
                duplication /= max_duplication_display;
            if (gain != 0. && max_gain_display!=0.)
                gain /= max_gain_display;
            
            double m = Double.NEGATIVE_INFINITY;
            if (draw_loss || draw_gain || draw_duplication)
            {
                if (draw_loss)
                    m = loss;
                if (draw_gain)
                    m = Math.max(m, gain);
                if (draw_duplication)
                    m = Math.max(m, duplication);
            } else
                m = 1.0;
            
            //System.out.println("#**RMD.gNEL "+m+"\tnode "+N+"\t"+draw_loss);
            
            return m;
        }        
        
        /**
         * We extend the acces level to public.
         * 
         * @return the appropriate display transformation from node positions to display positions
         */
        @Override
        public AffineTransform getDisplayTransform()
        {
            return super.getDisplayTransform();
        }
        
    }
    
    private Dimension getPreferredRateLegendDimension()
    {
        int f = LookAndFeel.TREE_PANEL_FONT_SIZE; // line height for text
        int d = LookAndFeel.TREE_LEGEND_INSET; // separation between elements
        int lw = LookAndFeel.TREE_PANEL_FONT_SIZE * 5/3; // width of individual bars
        int bar_height = LookAndFeel.TREE_LEGEND_BAR_HEIGHT; // max. height of the bars
        int w = lw*4+d*5; // total width
        int num_rates = 1; // loss
        if (main_tree.hasDuplication())
            num_rates ++;
        if (main_tree.hasGain())
            num_rates++;
        
        int h = num_rates*(2*d+2*f+bar_height)+d; // total height
        Dimension D = new Dimension (w,h);
        return D;
    }
    
    private class RateLegend extends JPanel implements PropertyChangeListener
    {
        private RateLegend()
        {
            super();
            setBackground(LookAndFeel.SMOKY_BACKGROUND);
            setBorder(javax.swing.BorderFactory.createRaisedBevelBorder());
            user_preferred = null;
        }
        
        private Dimension user_preferred;
        
        public void propertyChange(PropertyChangeEvent E)
        {
            if (EmbellishedTreePanel.MAGNIFICATION_PROPERTY.equals(E.getPropertyName())
                    && E.getSource() instanceof EmbellishedTreePanel)
            {
                // our kind of event!
                double old_mag = ((Double)E.getOldValue()).doubleValue();
                double new_mag = ((Double)E.getNewValue()).doubleValue();
                if (old_mag != new_mag)
                {
                    //revalidate();
                }
            }
        }
        
        @Override
        public void paintComponent(Graphics g)
        {
            //System.out.println("#*RMD.RL.pC start");
            super.paintComponent(g);
            
            int num_rates = 0; // loss
            if (tree_panel.draw_loss)
                num_rates++;
            if (tree_panel.draw_gain)
                num_rates++;
            if (tree_panel.draw_duplication)
                num_rates++;
            
            if (num_rates>0)
            {
                Graphics2D g2 = (Graphics2D) g.create();
                Color old_color = g2.getColor();

                int w = getWidth();
                int h = getHeight();


                int font_size = LookAndFeel.TREE_PANEL_FONT_SIZE*4/5;
                Font legend_font = new Font("Serif", Font.PLAIN, font_size);
                Font legend_title_font = legend_font.deriveFont(Font.BOLD);
                int skip = LookAndFeel.TREE_LEGEND_INSET; // this is how much space is left between elements
                int plot_height = (h - skip)/num_rates;            
                int bar_height = plot_height - 3*skip - 2*font_size;  
                int bar_width = (w-2*skip)/4;
                double scale_y = Math.abs(tree_panel.getDisplayTransform().getScaleY());

                int legend_y = skip;
                for (int rate_idx=0; rate_idx<3; rate_idx++)
                {
                    String legend_subtitle = "";
                    double max_display;
                    int[] max_exact;
                    Color rate_color;

                    boolean rate_shown = false;


                    if (rate_idx==0)
                    { // loss
                        legend_subtitle = "Loss rates";
                        max_display = tree_panel.max_loss_display;
                        max_exact = tree_panel.max_loss_display_exact;
                        rate_color = LookAndFeel.LOSS_COLOR;
                        rate_shown= tree_panel.draw_loss;
                    } else if (rate_idx==1)
                    {
                        legend_subtitle = "Duplication rates";
                        max_display = tree_panel.max_duplication_display;
                        max_exact = tree_panel.max_duplication_display_exact;
                        rate_color = LookAndFeel.DUPLICATION_COLOR;
                        rate_shown= tree_panel.draw_duplication;
                    } else 
                    {
                        legend_subtitle = "Gain rates";
                        max_display = tree_panel.max_gain_display;
                        max_exact = tree_panel.max_gain_display_exact;
                        rate_color = LookAndFeel.GAIN_COLOR;
                        rate_shown= tree_panel.draw_gain;
                    }

                    if (rate_shown)
                    {
                        g2.setFont(legend_title_font);
                        g2.drawString(legend_subtitle, skip, legend_y+font_size);
                        g2.setFont(legend_font);

                        double rate_ref = bar_height / scale_y;
                        if (rate_ref>=1.0)
                        {
                            double level0 = legend_y+font_size+skip+bar_height;
                            double level1 = level0 - scale_y;
                            g2.setColor(rate_color);
                            drawLegend(g2, font_size+skip, skip+bar_width/2, level0, level1, 0.5, max_exact);
                            g2.setColor(rate_color);
                            drawLegend(g2, font_size+skip, skip+3*bar_width/2, level0, level1, 1.0, max_exact);
                            g2.setColor(rate_color);
                            drawLegend(g2, font_size+skip, skip+5*bar_width/2, level0, level1, 2.0, max_exact);
                            g2.setColor(rate_color);
                            drawLegend(g2, font_size+skip, skip+7*bar_width/2, level0, level1, 4.0, max_exact);
                        } else
                        {
                            // we can accomodate a bar corresponding to a displayed length of rate_ref at most
                            double max_legend_length = rate_ref * max_display; // true length corresponding to rate_ref
                            int [] rounded_half_exact = new int[2];
                            double rounded_half = roundToMostSignificantDigit(max_legend_length/2.0, rounded_half_exact, true);
                            double rounded_half_rel = rounded_half/max_display;
                            double level0 = legend_y+font_size+skip+bar_height;
                            double level1 = level0 - 2.0*rounded_half_rel * scale_y;
                            rounded_half_exact[0] *= 2.0;

                            //System.out.println("#*RMD.RL.pC rate "+rate_idx+"\tbh "+bar_height+"\tref "+rate_ref
                            //        +"\tmax "+max_display+"\tscale "+scale_y+"\tmll "+max_legend_length
                            //        +"\trh "+rounded_half+"\tl0 "+level0+"\tl1 "+level1+"// "+legend_subtitle);

                            g2.setColor(rate_color);
                            drawLegend(g2, font_size+skip, skip+bar_width/2, level0, level1, 0.25, rounded_half_exact);
                            g2.setColor(rate_color);
                            drawLegend(g2, font_size+skip, skip+3*bar_width/2, level0, level1, 0.5, rounded_half_exact);
                            g2.setColor(rate_color);
                            drawLegend(g2, font_size+skip, skip+5*bar_width/2, level0, level1, 0.75, rounded_half_exact);
                            g2.setColor(rate_color);
                            drawLegend(g2, font_size+skip, skip+7*bar_width/2, level0, level1, 1.0, rounded_half_exact);
                        }
                        legend_y += plot_height;
                    } else
                    {
                        //g2.setFont(legend_title_font);
                        //g2.drawString("["+legend_subtitle+"]", skip, legend_y+font_size);
                    }
                }
            }
            
            //g2.drawString(Double.toString(tree_panel.getDisplayTransform().getScaleY()), 2, 20);
        }

        private void drawLegend(Graphics2D g2, double y_unit, int x, double level0, double level1, double factor, int[] scale)
        {
            tree_panel.drawEdge(g2, factor, 1.0, x, level1, level0);
            g2.setColor(Color.BLACK);
            DrawString.drawCenteredString(g2, getScaledValueDisplay(factor, scale), x, (int)(level0+y_unit+0.5));
        }
        
        private String getScaledValueDisplay(double scaling_factor, final int[] scale)
        {
            double mantissa = scaling_factor*scale[0];
            int exponent = scale[1];
            if (mantissa<0.1)
            {
                mantissa *= 10.0;
                exponent--;
            } else if (mantissa>=10.0)
            {
                mantissa /= 10.0;
                exponent++;
            }
            // check if we can get away without scientific notation
            if (exponent<=0)
            {
                if (exponent>=-2)
                    while (exponent<0){exponent++; mantissa*=0.1;}
            } else if (exponent<=2)
                while (exponent>0){exponent--; mantissa*=10.0;}
            
            // get rid of floating-point errors
            //mantissa = ((int)(mantissa*1000000.0+0.5))/1000000.0;
            
            String retval = "";
            if (mantissa == (int)mantissa)
            {
                retval = Integer.toString((int)mantissa);
            } else
            {
                retval = Float.toString((float)mantissa);
            }
            if (exponent!=0)
            {
                retval = retval + "e" + Integer.toString(exponent);
            }
            //System.out.println("#*RMD.RL.gSVD factor "+scaling_factor+"\tretval `"+retval+"'\tm "+mantissa+"\te "+exponent+"\tsc "+scale[0]+"\t"+scale[1]);
            
            return retval;
        }

       /**
         * Sets a user-preferred size, which overrules the locally computed values.
         * Reset to default by calling with null argument.
         */
        @Override
        public void setPreferredSize(Dimension size)
        {
            this.user_preferred = size;
            //System.out.println("#**ISD.sPS "+size);
        }
        

        @Override
        public Dimension getPreferredSize()
        {
            Dimension D = user_preferred;
            if (D ==null)
            {
                D = getPreferredRateLegendDimension();
            }
            //System.out.println("#*RMD.RL.gPS "+D+"\t"+tree_panel.getDisplayTransform().getScaleY());
            return D;
        }
        
        @Override
        public Dimension getMinimumSize()
        {
            return getPreferredSize();
        }        
    }
    
    
    protected class RateZoomableTreePanel extends ZoomableTreePanel
    {
        public RateZoomableTreePanel()
        {
            super(tree_panel);
        }
        
        @Override
        protected void initBottomBarElements(Box bb)
        {
            createBottomBarLeft(bb);

            super.initBottomBarElements(bb);
        }
    }
    
    public static class DiscreteGammaPlot implements Icon
    {
        private static final int DEFAULT_PARTITIONS = 8;

        /**
         * Instantiation with default size and default partition. 
         * 
         * @param distribution a discretized Gamma distribution that we want to plot
         */
        public DiscreteGammaPlot(DiscreteGamma distribution)
        {
            this(distribution, DEFAULT_PARTITIONS);
        }

        public DiscreteGammaPlot(DiscreteGamma distribution, int num_categories)
        {
            this(distribution, num_categories, RATE_VARIATION_PLOT_WIDTH, RATE_VARIATION_PLOT_HEIGHT);
        }

        public DiscreteGammaPlot(DiscreteGamma distribution, int num_categories, int width, int height)
        {
            setDistribution(distribution);
            setNumCategories(num_categories);
            setWidth(width);
            setHeight(height);
            background_color=legend_color=plotting_color= null;
        }

        private DiscreteGamma distribution;
        private int num_categories;

        private int width;
        private int height;

        /**
         * Background color
         */
        private Color background_color;
        /**
         * Color used for drawing the axes and parameter info
         */
        private Color legend_color;
        /**
         * Color used for plotting the distribution
         */
        private Color plotting_color;
        
        /**
         * Color used for plotting the categories
         */
        private Color category_color;

        public void setWidth(int width){ this.width = width;}
        public void setHeight(int height){this.height = height;}
        
        public int getIconWidth(){return width;}
        public int getIconHeight(){return height;}


        /**
         * Sets the bakground color for the plot. By default, this is null, meaning that 
         * no background is painted (transparent icon).
         * 
         * @param c a color (null means default behavior)
         */
        public void setBackground(Color c)
        {
            this.background_color = c;
        }

        /**
         * Sets the color for the axes, and parameter info. By 
         * default, this is the normal color for the Graphics 
         * context when paintComponent() is called.  
         * 
         * @param c a color (null means default behavior)
         */
        public void setLegendColor(Color c)
        {
            this.legend_color = c;
        }

        /**
         * Sets the plotting color for the distribution. By 
         * default, this is the same as the legend color   
         * 
         * @param c
         */
        public void setPlottingColor(Color c)
        {
            this.plotting_color = c;
        }
        
        /**
         * Sets the plotting color for the categories. By 
         * default, this is the same as the plotting color
         * 
         * @param c plotting color
         */
        public void setCategoryColor(Color c)
        {
            this.category_color = c;
        }
        
        /**
         * Sets the distribution that is to be plotted here.
         * 
         * @param D a discrete Gamma distribution 
         */
        public void setDistribution(DiscreteGamma D){distribution = D;}    

        /**
         * Sets the number of partitions used in the discretization
         * 
         * @param num_partitions number of partitions 
         */
        public void setNumCategories(int num_categories){this.num_categories=num_categories;}

        /**
         * Plots the distribution at a given location. 
         * 
         * @param ignored this parameter is ignored, may be even null (required parameter for Icon interface)
         * @param g graphics context (should be convertible to Graphics2D)
         * @param x coordinate for the left-hand side
         * @param y coordinate for the top
         */
        public void paintIcon(java.awt.Component ignored,
                   Graphics g,
                   int x,
                   int y)
        {
            Graphics2D g2 = (Graphics2D)g.create(); // we will work with a copy
            Color original_color = g2.getColor();
            Stroke original_stroke=g2.getStroke();
            // draw background
            if (background_color != null)
            {
                g2.setColor(background_color);
                g2.fillRect(x, y, width, height);
                g2.setColor(original_color);
            }

            // now comes the legend: X, Y axes and tics
            if (legend_color != null)
                g2.setColor(legend_color);
            
            // parameters for the plot
            int font_size = g2.getFont().getSize();
            int tic_length = 3;
            
            // how much pace is needed outside the axes
            int axis_offset_x = tic_length+g2.getFontMetrics().stringWidth("9e-9")+1;
            int axis_offset_y = tic_length+font_size*6/5;

            // position of the origin
            int origin_x = x+axis_offset_x+1;
            int origin_y = y+height-axis_offset_y+1;
            
            // size of the the actual plot area 
            int plot_width  = width - axis_offset_x - tic_length-1;
            int plot_height = height- axis_offset_y - tic_length-1;
            
            // draw the axes
            g2.drawLine(origin_x-1, origin_y+1, origin_x-1, origin_y-plot_height);
            g2.drawLine(origin_x-1, origin_y+1, origin_x+plot_width, origin_y+1);

            // get the discretization for the distribution 
            double[] partition_boundary = new double[num_categories-1];
            double[] partition_center = distribution.getPartitionMeans(num_categories, partition_boundary);

            // horizontal mapping for plot area 
            double max_value_x = partition_center[num_categories-1]*1.05; // leave a little bit after the last quantile's mean
            double scale_x = plot_width/max_value_x; // x is plotted at coordinate origin_x+scale_x*x
            
            // draw the tics
            double unit_x = roundToMostSignificantDigit(max_value_x/6.0,null, true);  // just a few tics 

            for (int i=1; i*unit_x<max_value_x; i++)
            {
                double r = i*unit_x; // this is the position of this tic
                int tic_x = (int)(r*scale_x); // transform to display coordinates
                g2.drawLine(origin_x+tic_x, origin_y, origin_x+tic_x, origin_y+tic_length);
                if (i==1)
                    DrawString.drawCenteredString(g2, Float.toString((float)r), origin_x+tic_x, origin_y+tic_length+font_size);
            }
            
            
            // mode is at x=(alpha-1)/alpha, but here the actual plotted maximum is computed instead
            double max_value_y = 0.0; // will be positive
            
            // computed the points to be plotted
            double[] pdf = new double[plot_width];
            {
                double alpha = distribution.getAlpha();
                double d = alpha * Math.log(alpha)-Functions.gammln(alpha);
                for (int i=0; i<plot_width; i++)
                {
                    // conversion from display coordinate to x
                    double r = (i+1.0)/scale_x;
                    double ln_f = d-alpha*r+(alpha-1.0)*Math.log(r);
                    double f = Math.exp(ln_f); // this is the value of the density function at r
                    pdf[i] = f;
                    //System.out.println("#*RMD.DGP pI "+i+"\tr "+r+"\tf "+f);
                    if (f>max_value_y) max_value_y=f;
                }
            }
            double scale_y = plot_height/max_value_y;

            
            // plot the distributions
            if (plotting_color != null)
                g2.setColor(plotting_color);
            g2.setStroke(new BasicStroke(1.f));
            for (int i=0; i<plot_width; i++)
            {
                int bar_height = (int)(scale_y*pdf[i]);
                //System.out.println("#*RMD.DGP pI "+i+"\tf "+pdf[i]+"\tbh "+bar_height+"\t// scale "+scale_y+"\tph "+plot_height+"\tmaxy "+max_value_y);
                g2.drawLine(origin_x+i, origin_y, origin_x+i, origin_y-bar_height);
            }
            
            // plot categories
            if (category_color != null)
                g2.setColor(category_color);
            
            g2.setStroke(new BasicStroke(2.f));
            for (int i=0; i<partition_center.length; i++)
            {
                int loc_x = (int)(partition_center[i]*scale_x);
                g2.drawLine(origin_x+loc_x, origin_y, origin_x+loc_x, origin_y-2*plot_height/4);
            }
            
            g2.setStroke(original_stroke);
            g2.setColor(original_color);
            
            if (legend_color != null)
                g2.setColor(legend_color);
            
            g2.drawString(num_categories+" categories, \u03b1="+distribution.getAlpha() //DoubleRoundedForDisplay.toString(distribution.getAlpha())+"\u22ef"
                    , origin_x+50, origin_y-plot_height+font_size);
        }
    }
    
    

    
  
}
