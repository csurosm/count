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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Rectangle2D;
import java.util.Arrays;

import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.JRadioButton;
import javax.swing.ListSelectionModel;

import count.ds.IndexedTree;
import count.gui.kit.DiscreteDistributionPlot;
import count.gui.kit.DrawString;
import count.io.DataFile;
import count.matek.DiscreteDistribution;
import count.matek.Functions;
import count.model.TreeWithRates;


import count.model.GLDParameters;

/**
*
* A {@link TreePanel} that shows the model's rate parameters on the phylogeny. 
* 
* @author Mikl&oacute;s Cs&#369;r&ouml;s 
*/
public class RatesTreePanel extends TreePanel 
{
	
	public RatesTreePanel(DataFile<TreeWithRates> rates_data)
	{
		super(rates_data.getContent().getTree(), TreePanel.LayoutStyle.PHENOGRAM, ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		this.rates_data = rates_data;
		initGraphicalElements();
	}

	public static final double RATE_TREE_SHORT_EDGE = 0.05;
    public static final int RATE_TREE_ROOT_DISTRIBUTION_WIDTH = 120;
    public static final int  RATE_TREE_ROOT_DISTRIBUTION_HEIGHT = 40;
    public static final Color RATE_TREE_DISTRIBUTION_PLOT_COLOR = Color.RED;

    public static final Color RATES_LOSS_COLOR = new Color(255,128,0); // Tangerine
    public static final Color RATES_GAIN_COLOR = new Color(64,128,0); // Fern
    public static final Color RATES_DUPLICATION_COLOR = new Color(104,118,231); // Evening Blue
    public static final Color RATES_EDGELENGTH_COLOR =  Color.DARK_GRAY;
    public static final Color RATE_TREE_EDGE_COLOR =  new Color(180,180,180);
    
	private final DataFile<TreeWithRates> rates_data;

    private final float[] dash = {2.0f, 4.0f};
    private final Stroke dashed_stroke = new BasicStroke(2.f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 10.0f, dash, 0);
    private final Stroke solid_stroke = new BasicStroke(2.f);
	
    private GLDParameters.Type displayed_parameter = null;
    
//    private double max_length_display;
//    private int[] max_length_exact;
    
    private void initGraphicalElements()
    {
        IndexedTree main_tree = getTree();
        for (int leaf_idx=0; leaf_idx<main_tree.getNumLeaves(); leaf_idx++)
        {
        	DisplayedNode L = getNode(leaf_idx);
//        	L.setPosition(0, 0, -Math.PI*0.4, 0f);
            L.setColor(Color.DARK_GRAY);
        }
        setRootStem(0); // we will have a window instead
        setBoundingBoxSeparation(getNormalFontSize()/2);
        setDisplayedParameter(GLDParameters.Type.LOSS);
    }
    
    
    public DataFile<TreeWithRates> getRatesData()
    {
 	   return this.rates_data;
    }    
    
    public void setDisplayedParameter(GLDParameters.Type param)
    {
    	this.displayed_parameter = param;
    	setTreeLayoutStyle(LayoutStyle.PHYLOGRAM); // recalculate layout and repaint 
    }
    
    /**
     * Creates a Box with three radio buttons for choosing which parameter to show. 
     * 
     * @return
     */
    public Box createParameterChooser()
    {
        TreeWithRates rates_model = getModel();

        String apple_button_size = "small";
        
        JRadioButton gainRB = new JRadioButton("gain");
        gainRB.putClientProperty("JComponent.sizeVariant", apple_button_size);        
        gainRB.setSelected(GLDParameters.Type.GAIN.equals(displayed_parameter));
        gainRB.addActionListener(button_pushed->setDisplayedParameter(GLDParameters.Type.GAIN));
        gainRB.setForeground(RATES_GAIN_COLOR);
        gainRB.setEnabled(rates_model.hasGain());

        JRadioButton lossRB = new JRadioButton("loss");
        lossRB.putClientProperty("JComponent.sizeVariant", apple_button_size);        
        lossRB.setSelected(GLDParameters.Type.LOSS.equals(displayed_parameter));
        lossRB.addActionListener(button_pushed->setDisplayedParameter(GLDParameters.Type.LOSS));
        lossRB.setForeground(RATES_LOSS_COLOR);
        
        JRadioButton dupRB = new JRadioButton("duplication");
        dupRB.putClientProperty("JComponent.sizeVariant", apple_button_size);        
        dupRB.setSelected(GLDParameters.Type.DUPLICATION.equals(displayed_parameter));
        dupRB.addActionListener(button_pushed->setDisplayedParameter(GLDParameters.Type.DUPLICATION));
        dupRB.setForeground(RATES_DUPLICATION_COLOR);
        dupRB.setEnabled(rates_model.hasDuplication());
        
        ButtonGroup parameter_choices = new ButtonGroup();
        parameter_choices.add(gainRB);
        parameter_choices.add(lossRB);
        parameter_choices.add(dupRB);
        
        
        
        Box B = Box.createHorizontalBox();
        B.add(lossRB);
        B.add(dupRB);
        B.add(gainRB);
        
        return B;
    }
    
    
    
    /**
     * Sets magnification level (and adjusts bounding box separation accordingly)
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
    }    
            
    
//    @Override
//    protected double getDisplayEdgeLength(int node)
//    {
//    	if (displayed_parameter==null) // called during super's instantiation
//    	{
//    		return super.getDisplayEdgeLength(node); // init with phylogeny's edge lengths
//    	}
//    	double r = getRate(node);
//    	if (max_length_display != 0.0)
//    		r /= max_length_display;
//        double retval = Math.max(RATE_TREE_SHORT_EDGE, Math.min(r, 1.0));
//    	return retval;
//    }
    
    @Override
    protected void calculateBoundingBoxes(Graphics g)
    {
        if (g==null)
            super.calculateBoundingBoxes(g);
        else
        {
            Graphics2D g2 = (Graphics2D)g.create();
            double r = getNodeRadius();
            IndexedTree tree = getTree();
            int font_size = getNormalFontSize();
            for (int node_idx=0; node_idx<tree.getNumNodes(); node_idx++)
            {
            	DisplayedNode N = getNode(node_idx);
                Rectangle2D R = new Rectangle2D.Double(-r/2,-r/2,r+1, r+1);
                R.add(N.getLabelBoundingBox(g2));
                if (tree.isRoot(node_idx))
                {
                    {
                        Rectangle2D Rplot = new Rectangle2D.Double(-RATE_TREE_ROOT_DISTRIBUTION_WIDTH/2, font_size, RATE_TREE_ROOT_DISTRIBUTION_WIDTH, DiscreteDistributionPlot.DISTRIBUTION_PLOT_HEIGHT);
                        R.add(Rplot);
                    }
                }
                N.setNodeBoundingBox(R);
            }
        }
    }    
    
    @Override
    protected void paintComponent(Graphics g)
    {
        super.paintComponent(g);
        this.paintSelectedNodes(g);
    }
    
    @Override
    protected void paintNodes(Graphics g)
    {
//        Graphics2D g2 = (Graphics2D) g.create();
//        Font small_font = new Font("Serif",Font.PLAIN,getTreeNormalFontSize(getTreePointSize())*4/5);
//        
//        TreeWithRates model = getModel();
//        int root = model.getTree().getRoot();
//        DiscreteDistribution root_distribution
//        	= model.getLossParameter(root)==1.0
//        		?model.getGainDistribution(root)
//        		:model.getDuplicationDistribution(root);
//
//        DisplayedNode Rpt = getNode(root);
//        int nx = (int)(Rpt.getX()+0.5);
//        int ny = (int)(Rpt.getY()+0.5);
//
//        DiscreteDistributionPlot root_plot = new DiscreteDistributionPlot(root_distribution,RATE_TREE_ROOT_DISTRIBUTION_WIDTH, RATE_TREE_ROOT_DISTRIBUTION_HEIGHT);
//        root_plot.setBackground(DiscreteDistributionPlot.DISTRIBUTION_PLOT_BACKGROUND);
//        root_plot.setLegendColor(Color.BLACK);
//        root_plot.setPlottingColor(RATE_TREE_DISTRIBUTION_PLOT_COLOR);
//        g2.setFont(small_font);
//        root_plot.paintIcon(this, g2, nx-RATE_TREE_ROOT_DISTRIBUTION_WIDTH/2, ny+getTreeNormalFontSize(getTreePointSize()));

    }

    /**
     * Plots the inparalog and xenolog distributions at selected nodes.
     * 
     * @param g
     */
    private void paintSelectedNodes(Graphics g)
    {
        ListSelectionModel selection_model = getSelectionModel();
        if (!selection_model.isSelectionEmpty())
        {
            Graphics2D g2 = (Graphics2D) g.create();
            int small_font_size = Math.max(getLabelFontSize(), getNormalFontSize()); //*4/5;
            Font small_font = new Font("Serif",Font.PLAIN,small_font_size);
            Font small_title_font = small_font.deriveFont(Font.BOLD);//.deriveFont(Font.ITALIC);

            IndexedTree main_tree = getTree();
            
            int i0 = selection_model.getMinSelectionIndex();
            if (i0 != -1)
            {
                int i1 = selection_model.getMaxSelectionIndex();
                for (int idx=i0; idx<=i1; idx++)
                    if (selection_model.isSelectedIndex(idx))
                    {
                    	DisplayedNode DN = getNode(idx);
                        double nx = DN.getX();
                        double ny = DN.getY();

                        //super.plotNode(g2, N);
                        g2.setColor(Color.RED);
                        g2.setStroke(solid_stroke);
                        double r = (int)(getNodeRadius()+0.5);
                        int d = (int)(2.*r+0.5);
                        
                        

                        if (main_tree.isRoot(idx))
                        {
//                            g2.drawOval((int)(nx-r),(int)(ny-r),d,d); // bof
                            TreeWithRates model = getModel();
                            DiscreteDistribution root_distribution = model.getRootDistribution();

                        	DiscreteDistributionPlot root_plot = new DiscreteDistributionPlot(root_distribution,RATE_TREE_ROOT_DISTRIBUTION_WIDTH, RATE_TREE_ROOT_DISTRIBUTION_HEIGHT);
                            root_plot.setBackground(TREE_NODE_INFO_BACKGROUND);
                            root_plot.setLegendColor(Color.BLACK);
                            root_plot.setPlottingColor(RATE_TREE_DISTRIBUTION_PLOT_COLOR);
                            g2.setFont(small_font);
                            root_plot.paintIcon(this, g2, (int)nx-RATE_TREE_ROOT_DISTRIBUTION_WIDTH/2, (int)ny+getTreeNormalFontSize(getTreePointSize()));
                        } else
                        {
                            d=Math.max(d, 4*small_font_size);
                            
                            int dh = (int)(1.5*d);
                            int dw = dh*5/2;
                            int dsep = 5;

                            int title_width = g2.getFontMetrics(small_title_font).stringWidth(main_tree.getName(idx)+"  ");
                            if (title_width>dw)
                                dw = title_width;

                            int dx = (int)nx;
                            int dy = (int)ny;
                            if (main_tree.isLeaf(idx) || true) // always put it below the node
                                dy += r+3;
                            else
                                dy -= r+2*dh+3*dsep;

                            if (dy<0) dy=0;
                            if (dx<dw/2+dsep) dx=dw/2+dsep;
                            if (dy+2*dh+3*dsep>getHeight())
                                dy=getHeight()-2*dh-3*dsep;
                            if (dx+dw/2+dsep>getWidth())
                                dx = getWidth()-dw/2-dsep;

                            TreeWithRates rates = getModel();
                            DiscreteDistribution xdistr = rates.getGainDistribution(idx);
                            DiscreteDistribution pdistr = rates.getDuplicationDistribution(idx);
                            DiscreteDistributionPlot pplot = new DiscreteDistributionPlot(pdistr, dw, dh);
                            Color plotbg =  TREE_NODE_INFO_BACKGROUND;
                            Color plotc = RATE_TREE_DISTRIBUTION_PLOT_COLOR;
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
                            	String node_name = main_tree.getName(idx);
                            	if (node_name == null) node_name = main_tree.getIdent(idx);
                                g2.setColor(legendc);
                                g2.setFont(small_title_font);
                                DrawString.drawCentered(g2, node_name, dx, dy);
                            }

                            pplot.setBackground(null);
                            pplot.setLegendColor(legendc);
                            pplot.setPlottingColor(plotc);
                            g2.setFont(small_font);
                            pplot.paintIcon(this, g2, dx-dw/2, dy);
                            g2.setColor(legendc);
                            g2.setFont(small_font.deriveFont(Font.BOLD));
                            g2.drawString("Inparalogs", dx-dw/2+6, dy+small_font_size);
                            DiscreteDistributionPlot xplot = new DiscreteDistributionPlot(xdistr, dw, dh);
                            xplot.setBackground(null);
                            xplot.setLegendColor(legendc);
                            xplot.setPlottingColor(plotc);
                            g2.setFont(small_font);
                            xplot.paintIcon(this, g2, dx-dw/2, dy+dh+dsep);
                            g2.setColor(legendc);
                            g2.setFont(small_font.deriveFont(Font.BOLD));
                            g2.drawString("Xenologs", dx-dw/2+6, dy+dh+dsep+small_font_size);
                        }
                    } // if selected
                // for
            } // at least 1 selected
        } // can be selected
    }

    @Override
    protected void paintNodeLabels(Graphics g)
    {
        
        Font leaf_font = new Font("Serif", Font.PLAIN, getLabelFontSize());
        Font inner_font = new Font("Serif", Font.ITALIC, getLabelFontSize());
        
        Graphics2D g2 = (Graphics2D) g.create();
        IndexedTree tree = getTree();
        
        for (int nidx=0; nidx<tree.getNumNodes(); nidx++)
        {
            if (!this.isSelected(nidx) || tree.isRoot(nidx))
            {
                Font label_font=(tree.isLeaf(nidx)?leaf_font:inner_font);
                g2.setFont(label_font);
                getNode(nidx).paintLabel(g2);
            }
        }
    }
    
    
    @Override
    protected void paintEdges(Graphics g)
    {
    	Color old_color = g.getColor();
    	switch (displayed_parameter)
        {
        	case GAIN: g.setColor(RATES_GAIN_COLOR); break;
        	case LOSS: g.setColor(RATES_LOSS_COLOR); break;
        	case DUPLICATION: g.setColor(RATES_DUPLICATION_COLOR); break;
        	default: g.setColor(RATES_EDGELENGTH_COLOR);
        }
    	super.paintEdges(g);
    	g.setColor(old_color);
    }    
    
//    /*
//     * Draws the vertical part for an edge: long edges have a proportional 
//     * dashed region. The plotted length of the edge is determined by the difference 
//     * between the end vertices' coordinates.
//     * 
//     * @param g2 graphics context: stroke may be altered on return 
//     * @param value the edge length to be displayed
//     * @param dlength maximum edge length displayed with solid stroke
//     * @param x_i position of the edge
//     * @param ny child's vertical coordinate
//     * @param py parent's vertical coordinate
//     */ 
//    private void drawEdge(Graphics2D g2, double value, double dlength, int x_i, double ny, double py)
//    {
//        int py_i = (int)(py+0.5);
//
//        if (value>dlength)
//        {
//            // long edge
//            double solid_fraction = dlength/value;
//            double my = (ny+(py-ny)*solid_fraction);
//            int my_i = (int)(my+0.5);
//            int ny_i = (int)(ny+0.5);
//            g2.setStroke(solid_stroke);
//            g2.drawLine(x_i, ny_i, x_i, my_i);
//            g2.setStroke(dashed_stroke);
//            g2.drawLine(x_i, my_i, x_i, py_i);
//        } else
//        {
//            double ey = py + (ny-py) * value / dlength;
//            int ey_i = (int)(ey+0.5);
//            g2.setStroke(solid_stroke);
//            g2.drawLine(x_i, py_i, x_i, ey_i);
//        }
//    }
     
    
//    /**
//     * Given an array of sorted edge lengths x[0]<=x[1]<=...<=x[n-1], 
//     * finds an L and the corresponding 
//     * quantiles x[0..t-1] and x[n-t..n-1] such that 
//     * when the maximum solid edge length is L, 
//     * only short edges (i<t) and long edges (i>=n-t)
//     * are displayed with a relative solid edge length below epsilon.
//     * 
//     * @param sorted_values edge lenth values in increasing order
//     * @param epsilon (0&lt;epsilon&lt;1 is assumed) the small fractional edge length below which there is no clear visual distinction in the display
//     * @return appropriate setting for L, the maximum edge length that is displayed normally (may be 0 in some cases!)
//     */
//    private static double getDisplayLengthThreshold(double[] sorted_values, double epsilon)
//    {
//        int n = sorted_values.length;
//        double f = epsilon; // epsilon*epsilon / (1.-epsilon);
//        double retval = 9e99; // such a weird return value would indicate a problem
//        for (int t=1; t-1<=n-t; t++)
//        {
//            double low = sorted_values[t-1]; // constraint L<= low/epsilon
//            double high  = sorted_values[n-t]; // constraint L>=high*epsilon/(1-epsilon)
//            if (low >= high*f)
//            { // found a good range 
//                retval = high; // Math.sqrt(low*high)/Math.sqrt(1.-epsilon);
//                break;
//            }
//        }
//
//        if (retval == 0.0)
//        {
//            for (int t=0; t<n; t++)
//                if (sorted_values[t]!=0.0)
//                {
//                    retval = sorted_values[(t+n-1)/2];
//                    break;
//                }
//        }
//
//        //retval = 2.0*sorted_values[n/2];
//
//        return retval;
//    }
    
    private TreeWithRates getModel() { return rates_data.getContent();}
    private IndexedTree getTree() { return getTreeData().getContent();}
    
    @Override 
    protected double getTrueEdgeLength(int node)
    {
    	TreeWithRates model = rates_data.getContent();
    	double r=1.0;
    	switch (displayed_parameter)
    	{
    		case GAIN:
    			r = model.getGainRate(node); 
    		case DUPLICATION:
    			r *= model.getDuplicationRate(node)*model.getEdgeLength(node); break;
    		case LOSS:
    			r = model.getLossRate(node)*model.getEdgeLength(node); break;
    		default:
    			r = 1.0;
    	}
    	return r;
    }
    
    @Override
    protected double getDisplayEdgeRangeFactor() { return 1.0/64.0;}
    
//    private void calculateMaxLengthDisplay()
//    {
//    	TreeWithRates model = rates_data.getContent();
//    	IndexedTree tree= model.getTree();
//    	int num_nodes = tree.getNumNodes();
//    	
//    	double[] rates = new double[num_nodes];
//    	for (int node=0; node<num_nodes; node++)
//    		rates[node] = getRate(node);
//    	Arrays.sort(rates);
//    	this.max_length_display = getDisplayLengthThreshold(rates, RATE_TREE_SHORT_EDGE);
//        this.max_length_exact = new int[2];
//        max_length_display = 2.*Functions.roundToMostSignificantDigit(max_length_display/2.0, max_length_exact);
//        max_length_exact[0] *= 2;
//        if (max_length_exact[0]>=10)
//        {
//        	max_length_exact[0]*=0.1;
//        	max_length_exact[1]++;
//        }
//    	
//    }

    @Override
    public String toString()
    {
    	return DataFile.chopFileExtension(rates_data.getFile().getName());
    }
}
