

package ca.umontreal.iro.evolution.malin.ui.count;


import java.awt.geom.AffineTransform;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;

import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JLayeredPane;
import javax.swing.JPanel;

import ca.umontreal.iro.evolution.TreeNode;

import ca.umontreal.iro.evolution.genecontent.NodeWithRates;
import ca.umontreal.iro.evolution.genecontent.OccurrenceTable;
import ca.umontreal.iro.evolution.genecontent.TreeWithRates;

import ca.umontreal.iro.evolution.malin.DataFile;
import ca.umontreal.iro.evolution.malin.DoubleRoundedForDisplay;

import ca.umontreal.iro.evolution.malin.ui.DrawString;
import ca.umontreal.iro.evolution.malin.ui.EmbellishedTreePanel;
import ca.umontreal.iro.evolution.malin.ui.LayeredBorderLayout;

/**
 *
 * Propensity for Gene Loss (PGL), introduced by
 * Krylov et al [Krylov, Wolf, Rogozin, Koonin. <q>Gene loss,
 * protein sequence divergence, gene dispensability,
 * expression level, and interactivity are
 * correlated in eukaryotic evolution</q>,
 * <em>Genome Research</em>, <b>13</b>:2229-2235, 2003].
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class PGLDisplay extends DolloDisplay
{
    /**
     * Instantiates a new PGL browser. 
     * 
     * @param main_tree the underlying phylogeny
     * @param data_file the occurrence table structure on which the ancestral reconstruction will be carried out
     * @param daddy the work space where this guy will go (need it to have original edge lengths)
     */
    public PGLDisplay(TreeWithRates main_tree, DataFile<OccurrenceTable> data_file, WorkSpaceCount daddy)
    {
        super(main_tree, data_file, daddy);
    }
    
    private double[] original_edge_length;
    private double[] pgl;
    private double[] original_tree_length;
    
    private LengthLegend length_legend;
    private JLayeredPane layers;
    
    
    @Override
    public String toString()
    {
        return "PGL";
    }
    
    @Override
    protected void init()
    {
        int num_edges = main_tree.getNumEdges();
        int num_nodes = main_tree.getNumNodes();

        original_edge_length = new double[num_edges];
        original_tree_length = new double[num_nodes];
        
        NodeWithRates[] nodes = main_tree.getDFT();
        for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
        {
            double len = work_space.getOriginalEdgeLength(edge_idx);
            original_edge_length[edge_idx]=len;
        }
        for (int node_idx=0; node_idx<num_nodes; node_idx++)
        {
            NodeWithRates N = nodes[node_idx];
            if (N.isLeaf())
                original_tree_length[node_idx]=0.0;
            else
            {
                double s = 0.0;
                int num_children = N.getNumChildren();
                for (int ci=0; ci<num_children; ci++)
                {
                    int child_idx = main_tree.getChildIndex(node_idx, ci);
                    s+=(original_edge_length[child_idx]+original_tree_length[child_idx]);
                }
                original_tree_length[node_idx]=s;
            }
        }
        
        super.init();
    }
    
    private Dimension getPreferredLengthLegendDimension()
    {
        int f = LookAndFeel.TREE_PANEL_FONT_SIZE; // line height for text
        int d = LookAndFeel.TREE_LEGEND_INSET; // separation between elements
        int lw = LookAndFeel.TREE_PANEL_FONT_SIZE * 5/3; // width of individual bars
        int bar_height = LookAndFeel.TREE_LEGEND_BAR_HEIGHT; // max. height of the bars
        int w = lw*4+d*5; // total width
        int num_rates = 1; // length
        
        int h = num_rates*(2*d+2*f+bar_height)+d; // total height
        Dimension D = new Dimension (w,h);
        return D;
    }
    

    @Override
    protected void createBottomBarLeft(Box bb)
    {
        super.createBottomBarLeft(bb);
        Font bb_font = new Font("Serif", Font.PLAIN, LookAndFeel.TREE_PANEL_FONT_SIZE);
        
        JCheckBox legendB = new JCheckBox("Legend");
        legendB.setFont(bb_font);
        legendB.setSelected(true);
        legendB.addItemListener(new ItemListener()
            {
                @Override
                public void itemStateChanged(ItemEvent E)
                {
                    int ch = E.getStateChange();
                    
                    if (ch == ItemEvent.DESELECTED)
                    {
                        tree_panel.setPadding(new java.awt.Insets(10,10,10,10));
                        length_legend.setVisible(false);
                    }
                    else if (ch==ItemEvent.SELECTED)
                    {
                        tree_panel.setPadding(new java.awt.Insets(10,10+getPreferredLengthLegendDimension().width,10,10));
                        length_legend.setVisible(true);
                    }
                    tree_panel.repaint();
                }
            });
        bb.add(legendB);
    }

    
    @Override
    protected void setZoomableTreePanel()
    {
        layers = new JLayeredPane();
        LayeredBorderLayout layers_layout = new LayeredBorderLayout();
        layers.setLayout(layers_layout);
        layers_layout.addLayoutComponent(LayeredBorderLayout.CENTER, tree_zoom);
        layers.add(tree_zoom, JLayeredPane.DEFAULT_LAYER);

        length_legend = createLengthLegend();
        layers_layout.addLayoutComponent(LayeredBorderLayout.WEST, length_legend);
        layers.add(length_legend, JLayeredPane.PALETTE_LAYER);

        setBottomComponent(layers);       
    }
    
    private LengthLegend createLengthLegend()
    {
        LengthLegend P = new LengthLegend();
        //P.setBorder(javax.swing.BorderFactory.createLineBorder(Color.RED));
        return P;
    }
    
    
    private class LengthLegend extends JPanel implements PropertyChangeListener
    {
        private LengthLegend()
        {
            super();
            setBackground(LookAndFeel.SMOKY_BACKGROUND);
            setBorder(javax.swing.BorderFactory.createRaisedBevelBorder());
            user_preferred = null;
            this.tree_panel = (PGLTreePane) PGLDisplay.this.tree_panel;
        }
        
        private PGLTreePane tree_panel;
        
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
            
            Graphics2D g2 = (Graphics2D) g.create();
            Color old_color = g2.getColor();

            int w = getWidth();
            int h = getHeight();


            int font_size = LookAndFeel.TREE_PANEL_FONT_SIZE*4/5;
            Font legend_font = new Font("Serif", Font.PLAIN, font_size);
            Font legend_title_font = legend_font.deriveFont(Font.BOLD);
            int skip = LookAndFeel.TREE_LEGEND_INSET; // this is how much space is left between elements
            int plot_height = (h - skip);            
            int bar_height = plot_height - 3*skip - 2*font_size;  
            int bar_width = (w-2*skip)/4;
            double scale_y = Math.abs(tree_panel.getDisplayTransform().getScaleY());

            int legend_y = skip;
            String legend_subtitle = "";
            double max_display;
            int[] max_exact;
            Color rate_color;

            legend_subtitle = "Edge lengths";
            max_display = tree_panel.display_edge_length_threshold;

            max_exact = tree_panel.display_edge_length_threshold_exact;
            rate_color = Color.BLACK;

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
                double rounded_half = RateModelDisplay.roundToMostSignificantDigit(max_legend_length/2.0, rounded_half_exact, true);
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
                D = getPreferredLengthLegendDimension();
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
        
    
    @Override
    protected HistoryTreePanel createTreePanel()
    {
        return new PGLTreePane();
    }
    
    private static final double TREE_PANEL_SHORT_EDGE_THRESHOLD = 1./12;
    
    protected class PGLTreePane extends DolloTreePane
    {
        protected PGLTreePane()
        {
            super(main_tree.getRoot(), EmbellishedTreePanel.LayoutStyle.SCALED_PHENOGRAM);
            super.setPadding(new java.awt.Insets(10,10+getPreferredLengthLegendDimension().width,10,10));
            super.setRootStemLength(0); // we will have a window instead

            float[] dash  = new float[2]; dash[0]=2.0f; dash[1]=4.0f;
            dashed_stroke = new BasicStroke(2.f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 10.0f, dash, 0);
            solid_stroke = new BasicStroke(2.f);
        }
        
        private BasicStroke solid_stroke;
        private BasicStroke dashed_stroke;

        double display_edge_length_threshold;
        int[] display_edge_length_threshold_exact;
        
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
        
        
        /**
         * Intercepted to calculate proper display thresholds. 
         */
        @Override
        protected void setupNodeLocations()
        {
            setEdgeLengthThreshold();
            super.setupNodeLocations();
        }
        
        private void setEdgeLengthThreshold()
        {
            double epsilon = TREE_PANEL_SHORT_EDGE_THRESHOLD;
            
            int num_edges = main_tree.getNumEdges();
            double[] len = new double[num_edges];
            for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
            {
                NodeWithRates N = main_tree.getNode(edge_idx);
                len[edge_idx] = work_space.getOriginalEdgeLength(N);
            }
            java.util.Arrays.sort(len);
            display_edge_length_threshold = RateModelDisplay.getDisplayLengthThreshold(len, epsilon);

            display_edge_length_threshold_exact = new int[2];
            display_edge_length_threshold = 2.0*RateModelDisplay.roundToMostSignificantDigit(display_edge_length_threshold/2.0, display_edge_length_threshold_exact);
            display_edge_length_threshold_exact[0]*=2;
            if (display_edge_length_threshold_exact[0]>=10)
            {
                display_edge_length_threshold_exact[0]*=0.1;
                display_edge_length_threshold_exact[1]++;
            }
            if (display_edge_length_threshold == 0.0 )
                System.out.println("#*PGLD.PGLTP.sELT gd "+RateModelDisplay.getDisplayLengthThreshold(len, epsilon )+"\tmax "+display_edge_length_threshold);
            
                    
        }
        
        protected double getNormalizedEdgeLength(NodeWithRates N)
        {
            double len = work_space.getOriginalEdgeLength((NodeWithRates)N);
            double retval = len;

            if (display_edge_length_threshold!=0.)
                retval /= display_edge_length_threshold;
            return retval;
        }

        
        @Override
        protected double getDisplayEdgeLength(TreeNode N)
        {
            double lN = getNormalizedEdgeLength((NodeWithRates)N);
            double retval = Math.max(TREE_PANEL_SHORT_EDGE_THRESHOLD, Math.min(lN, 1.0));
            
            //System.out.println("#*RMD.gDEL length "+retval+"\tnode "+N);
            return retval;
        }
    
        private void drawEdge(Graphics2D g2, double value, double dlength, int x_i, double ny, double py)
        {
            //System.out.println("#*PGLD.PGLTP.dE "+value+"/"+dlength+"\tx "+x_i+"\tny "+ny+"\tpy "+py);
            RateModelDisplay.drawEdge(g2, solid_stroke, dashed_stroke, value, dlength, x_i, ny, py);
        }            

        @Override
        public void plotEdges(Graphics g)
        {
            Graphics2D g2 = (Graphics2D) g.create();
            for (int node_idx=0; node_idx<tree_nodes.length; node_idx++)
            {
                NodeWithRates N = tree_nodes[node_idx];
                
                if (!N.isRoot()) // just a little stem
                {
                    int Nidx = getDisplayNodeIndex(N);
                    int Pidx = getDisplayNodeIndex(N.getParent());
                    
                    int nx = (int)(displayed_node_location[Nidx].getX()+0.5);
                    double ny = displayed_node_location[Nidx].getY();
                    int px = (int)(displayed_node_location[Pidx].getX()+0.5);
                    double py = displayed_node_location[Pidx].getY();
                    int py_i = (int)(py+0.5);
                                                
                    double length = getNormalizedEdgeLength(N);
                    double dlength = getDisplayEdgeLength(N);
                    
                    g2.setStroke(solid_stroke);
                    g2.setColor(EDGE_COLOR);
                    g2.drawLine(px, py_i, nx, py_i);
                    drawEdge(g2, length, dlength, nx, ny, py);
                }
            }
            plotEdgeChanges(g2,false);
        }                        

        @Override
        public void computeEdgeLabelBoundingBoxes(Graphics g)
        {
            // nuthin'
        }
    
    
    }

    @Override
    protected void initReconstructionVariables()
    {
        super.initReconstructionVariables();
        int num_families = families.length;
        pgl = new double[num_families];
    }
    
    @Override
    protected AncestralReconstructionPane.ComputingTask newComputingTask()
    {
        return new PGLComputingTask();
    }
    
    protected class PGLComputingTask extends DolloComputingTask
    {
        @Override
        protected void computeAncestralReconstruction(int family_idx)
        {
            super.computeAncestralReconstruction(family_idx);
            
            int num_edges = main_tree.getNumEdges();
            double tot_loss = 0.0;
            for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
            {
                double ls = getReconstructedFamilyLoss(family_idx, edge_idx);
                if (ls == 1.0)
                    tot_loss += original_edge_length[edge_idx];
            }
            int first = family_first_at_node[family_idx];
            
            double subtree_size = (first==-1?0.0:original_tree_length[first]);
            if (subtree_size == 0.0)
                pgl[family_idx]=0.0;
            else
                pgl[family_idx] = tot_loss/subtree_size;
            //System.out.println("#*PGLD.CT.cAR "+family_idx+"/"+families[family_idx]
            //        +"\tfirst "+family_first_at_node[family_idx]+"/"+getLongNodeName(family_first_at_node[family_idx])
            //        +"\tls "+tot_loss+"\ttree "+subtree_size);
        }
    }
    
    @Override
    protected OccurrenceTableModel createFamilyTableModel()
    {
        return new FamilyTableModel();
    }
    
    protected class FamilyTableModel extends OccurrenceTableModel
    {     
        protected FamilyTableModel()
        {
            super(data_file.getData(), false);
        }
        
        @Override
        public int getColumnCount()
        {
            return super.getColumnCount()
                    +1; // PGL
        }

        @Override
        public String getColumnName(int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            if (column_idx < num_usual_columns)
                return super.getColumnName(column_idx);
            else
                column_idx-=num_usual_columns;
            if (column_idx==0)
                return "PGL";
            else
                return "[No such column]";
        }

        @Override
        public Class getColumnClass(int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            if (column_idx<num_usual_columns)
                return super.getColumnClass(column_idx);
            else 
                column_idx-= num_usual_columns;
            if (column_idx==0)
                return DoubleRoundedForDisplay.class;
            else
                return Object.class;
        }
        
        @Override
        public String getColumnHeaderToolTip(int column_idx)        
        {
            int num_usual_columns = super.getColumnCount();
            //System.out.println("#**PGLD.FTM.gCHTT "+column_idx+"/"+num_usual_columns+"\tprops "+data_table.getKnownPropertiesCount());
            if (column_idx < num_usual_columns)
                return super.getColumnHeaderToolTip(column_idx);
            else
                column_idx-=num_usual_columns;
            
            if (column_idx==0)
                return "PGL: Propensity for Gene Loss (Krylov et al. 2003)";
            else
                return "[No such column]";
        }
        
        
        @Override
        public Object getValueAt(int row_idx, int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            if (column_idx<num_usual_columns)
                return super.getValueAt(row_idx, column_idx);
            else 
                column_idx-= num_usual_columns;
            if (column_idx==0)
                return new DoubleRoundedForDisplay(pgl[row_idx]);
            else
                return "[No such column]";
        }
        
        @Override 
        public String getCellToolTip(int row_idx, int column_idx)
        {
            Object value = getValueAt(row_idx,column_idx);
            int num_usual_columns = super.getColumnCount();
            if (column_idx < num_usual_columns)
                return super.getCellToolTip(row_idx, column_idx);
            else
                column_idx -= num_usual_columns;
            if (column_idx==0)
                return Double.toString(pgl[row_idx]);
            else
                return "[No such column]";
        }
        
    }

    @Override
    protected AncestralReconstructionPane.LineageTableModel createLineageTableModel()
    {
        return new LineageTableModel();
    }
    
    protected class LineageTableModel extends DolloDisplay.LineageTableModel
    {
        @Override
        public int getColumnCount()
        {
            int num_usual_columns = super.getColumnCount();
            return num_usual_columns
                    +1;
        }

        @Override
        public String getColumnName(int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            if (column_idx<num_usual_columns)
                return super.getColumnName(column_idx);
            else
                column_idx -= num_usual_columns;
            if (column_idx==0)
                return "Length";
            else 
                return "[No such column]";
        }
        
        @Override
        public String getHeaderToolTip(int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            if (column_idx<num_usual_columns)
                return super.getHeaderToolTip(column_idx);
            else
                column_idx -= num_usual_columns;
            if (column_idx==0)
                return "Edge length (from input tree in "+data_file.getFile().getName()+")";
            else 
                return "[No such column]";
        }

        @Override
        public Object getValueAt(int row_idx, int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            if (column_idx<num_usual_columns)
                return super.getValueAt(row_idx, column_idx);
            else
                column_idx -= num_usual_columns;
            if (column_idx==0)
            {
                if (main_tree.getNode(row_idx).isRoot())
                    return null;
                else
                    return new DoubleRoundedForDisplay(original_edge_length[row_idx]);
            }
            else 
                return "[No such column]";
        }   
        
        @Override
        public Class getColumnClass(int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            if (column_idx<num_usual_columns)
                return super.getColumnClass(column_idx);
            else
                return DoubleRoundedForDisplay.class;
        }
            
    }
}
