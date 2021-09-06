
package ca.umontreal.iro.evolution.malin.ui;


import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Insets;

import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;

import ca.umontreal.iro.evolution.TreeNode;

/**
 * This is a TreePanel for including <q>labels</q>
 * at nodes and edges.
 *
 * @author csuros
 */
public class EmbellishedTreePanel extends TreePanel
{
    public EmbellishedTreePanel(TreeNode root, LayoutStyle layout_style, boolean node_selection_allowed, boolean area_selection_allowed)
    {
        super(root, layout_style, node_selection_allowed, area_selection_allowed);
        setPadding(new Insets(5,5,5,5));
        initBoundingBoxes();
    }
    
    
    public EmbellishedTreePanel(TreeNode root, boolean node_selection_allowed, boolean area_selection_allowed)
    {
        this(root, LayoutStyle.CLADOGRAM, node_selection_allowed, area_selection_allowed);
    }

    public EmbellishedTreePanel(TreeNode root, boolean node_selection_allowed)
    {
        this(root, node_selection_allowed, node_selection_allowed);
    }
    
    public EmbellishedTreePanel(TreeNode root)
    {
        this(root,true,true);
    }
    
    /**
     * Padding around the tree in the panel. 
     * 
     * @param I the new padding values
     */
    public void setPadding(Insets I)
    {
        this.padding = I;
    }
    
    protected Insets padding;
    
    @Override
    protected AffineTransform getDisplayTransform()
    {
        int w=getWidth();
        int h=getHeight();
        
        int last_leaf_idx = node.length-1;
        while (!node[last_leaf_idx].isLeaf()) last_leaf_idx--;
        
        double x_stretch = (w-
                     padding.left 
                    + node_label_bounding_box[0].getX()
                    - node_label_bounding_box[last_leaf_idx].getX()
                    - node_label_bounding_box[last_leaf_idx].getWidth()
                    - padding.right)/(number_of_leaves-1.0);
        
        double stem = Math.max(root_stem_length, node_label_bounding_box[node.length-1].getY()+node_label_bounding_box[node.length-1].getHeight());
            
        double y_stretch = (h-padding.top-padding.bottom
                    - stem
                    + node_label_bounding_box[0].getY())/tree_size[DEPTH];
        AffineTransform plot_transform=new AffineTransform();        
        plot_transform.setToIdentity();

        plot_transform.translate(padding.left-node_label_bounding_box[0].getX(), h-padding.bottom-stem);
        plot_transform.scale(x_stretch,-y_stretch);
        return plot_transform;
    }
    
    
    /**
     * Each node has a bounding box (BB) for its label. The BB is assumed 
     * to cover at least the symbol used to plot the node.
     * 
     * @param N the node
     * @param dx relative X position of the BB with respect to the node cordinates
     * @param dy relative Y position of the BB with respect to the node cordinates
     * @param w width of the BB 
     * @param h height of the BB
     */
    public void setNodeLabelBoundingBox(TreeNode N, double dx, double dy, double w, double h)
    {
        Rectangle2D R = new Rectangle2D.Double(dx,dy,w,h);
        setNodeLabelBoundingBox(N, R);
    }

    /**
     * Each node has a bounding box (BB) for its label. The BB is assumed 
     * to cover at least the symbol used to plot the node.
     * 
     * @param N the node
     * @param R the bounding box for the node 
     */
    public void setNodeLabelBoundingBox(TreeNode N, Rectangle2D R)
    {
        int node_idx = getDisplayNodeIndex(N);
        node_label_bounding_box[node_idx] = R;
    }    
    
    /**
     * Each edge has a bounding box (BB) for its label. The BB is assumed 
     * to cover at least the edge thickness.
     * 
     * @param N the node to which the edge leads
     * @param dx relative X position of the BB with respect to the child node cordinates
     * @param w width of the BB 
     * @param h height of the BB
     */
    public void setEdgeLabelBoundingBox(TreeNode N, double dx, double w, double h)
    {
        Rectangle2D R = new Rectangle2D.Double(dx, 0, w, h);
        int node_idx = getDisplayNodeIndex(N);
        edge_label_bounding_box[node_idx] = R;
    }
    
    protected Rectangle2D[] node_label_bounding_box;
    
    protected Rectangle2D[] edge_label_bounding_box;
    
    /**
     * Hook for recomputing the size of the bounding boxes in a 
     * manner that depends on the graphics context. Does not do anything 
     * in the default implementation.
     * 
     * @param g graphics context
     */
    public void computeNodeLabelBoundingBoxes(Graphics g)
    {
        
    }
    
    /**
     * Hook for recomputing the size of the bounding boxes in a 
     * manner that depends on the graphics context.
     * It does not do anything in the default implementation. 
     * 
     * @param g graphics context
     */
    public void computeEdgeLabelBoundingBoxes(Graphics g)
    {
        
    }

    /**
     * Allocates up the bounding box arrays, and initializes the bounding boxes. 
     */
    private void initBoundingBoxes()
    {
        node_label_bounding_box = new Rectangle2D[node.length];
        edge_label_bounding_box = new Rectangle2D[node.length];
        for (int node_idx=0; node_idx<node.length; node_idx++)
        {
            TreeNode N = node[node_idx];
            setNodeLabelBoundingBox(N, -POINT_SIZE/2, -POINT_SIZE/2, POINT_SIZE, POINT_SIZE);
            setEdgeLabelBoundingBox(N, -1, 3, 0);
        }
    }
    
    protected int root_stem_length; 
    
    /**
     * Sets the length of the little stem leading to the root
     * @param len stem length in pixels
     */
    public void setRootStemLength(int len)
    {
        root_stem_length = len;
    }

    protected int bounding_box_separation = 12;
    
    public void setBoundingBoxSeparation(int d)
    {
        this.bounding_box_separation = d;
    }
    
    protected double magnification = 1.0;
    public static final String MAGNIFICATION_PROPERTY = "magnification";
    
    public void setMagnification(double mag)
    {
        double old_mag = this.magnification;
        this.magnification = mag;
        label_font_size = (int)(normal_label_font_size * Math.sqrt(magnification)+0.5);
        if (label_font_size<6)
            label_font_size = 6;
        if (label_font_size>32)
            label_font_size = 32;
        this.firePropertyChange(MAGNIFICATION_PROPERTY, old_mag, mag);
    }
    
    public void setNormalFontSize(int d)
    {
        this.normal_label_font_size = d;
        setMagnification(magnification);
    }
    
    protected int normal_label_font_size = 14;
    
    public double getMagnification()
    {
        return this.magnification;
    }
    
    /**
     * Computes the preferred size of this panel, unless 
     * it was set previously via setPreferredSize()
     * 
     * @return the preferred size
     */
    @Override
    public Dimension getPreferredSize()
    {
        if (user_preferred==null)
        {
            // figure out minimum stretch for X coordinates
            double min_x_stretch = 0.0;
            int last_leaf_idx = -1;
            {
                int previous_leaf_idx = -1;
                for (int node_idx=0; node_idx<node.length; node_idx++)
                {
                    TreeNode N = node[node_idx];
                    if (N.isLeaf())
                    {
                        if (previous_leaf_idx!=-1)
                        {
                            Rectangle2D leaf_BB = node_label_bounding_box[previous_leaf_idx];
                            double previous_leaf_ends = leaf_BB.getX()+leaf_BB.getWidth();
                            Rectangle2D edge_BB = edge_label_bounding_box[previous_leaf_idx];
                            double previous_edge_ends = edge_BB.getX()+edge_BB.getWidth();
                            double previous_ends = Math.max(previous_leaf_ends, previous_edge_ends);
                            
                            double current_leaf_starts = node_label_bounding_box[node_idx].getX();
                            double current_edge_starts = edge_label_bounding_box[node_idx].getX();
                            double current_starts = Math.min(current_leaf_starts, current_edge_starts);
                            
                            double current_min_stretch = previous_ends - current_starts;
                            if (current_min_stretch>min_x_stretch)
                                min_x_stretch = current_min_stretch;
                        }
                        previous_leaf_idx = node_idx;
                    }
                } // for all nodes
                last_leaf_idx = previous_leaf_idx;
            }
            min_x_stretch += bounding_box_separation;
            
            //min_x_stretch *= magnification;
            
            double minimum_width 
                    = padding.left 
                    - node_label_bounding_box[0].getX()
                    + min_x_stretch * (number_of_leaves-1.0)
                    + node_label_bounding_box[last_leaf_idx].getX()
                    + node_label_bounding_box[last_leaf_idx].getWidth()
                    + padding.right;
            
            double min_y_stretch = 0.0;
            {
                for (int node_idx=0; node_idx<node.length; node_idx++)
                {
                    TreeNode N = node[node_idx];
                    if (!N.isRoot())
                    {
                        TreeNode P = N.getParent();
                        int parent_idx = getDisplayNodeIndex(P);
                        Rectangle2D cR = node_label_bounding_box[node_idx];
                        Rectangle2D eR = edge_label_bounding_box[node_idx];
                        Rectangle2D pR = node_label_bounding_box[parent_idx];
                        
                        double eh = eR.getHeight();
                        double current_stretch = (cR.getY()+cR.getHeight()+eh-pR.getY());

                        int x = (eh==0.?0:1);
                        if (eh!=0. && pR.getHeight()!=0.)
                        {
                            x++;
                        }
                        current_stretch += x * bounding_box_separation;

                        double diff_y = node_location[node_idx].getY()-node_location[parent_idx].getY();
                        current_stretch /= diff_y;
                        
                        //System.out.println("#*ETP.gPS "+node_idx+"\tstretch "+current_stretch+"\tNbb "+cR.getY()+"\tPbb "+pR.getY()+"\tdiff "+diff_y+"\tx "+x+"\tmin "+min_y_stretch);
                        
                        if (current_stretch > min_y_stretch)
                            min_y_stretch = current_stretch;
                    }
                }
            }
            
            min_y_stretch *= magnification;
            
            double farthest_label = 0.0;
            {
                for (int node_idx=0; node_idx<node.length; node_idx++)
                {
                    Rectangle2D cR = node_label_bounding_box[node_idx];
                        
                    double label_y = node_location[node_idx].getY()*min_y_stretch-cR.getY();//+cR.getWidth();
                    if (label_y > farthest_label)
                        farthest_label = label_y;
                }
                
            }
            
            double minimum_height = padding.top + root_stem_length + padding.bottom 
                    + farthest_label;
            
            
            double reasonable_height = 300*magnification;
            //int w = (int)(Math.max(minimum_height,reasonable_height)/minimum_height*minimum_width);
            int h = (int) Math.max(minimum_height, reasonable_height);
            
            int w = (int)Math.max(minimum_width,12*node.length*magnification);
            //int h = (int)minimum_height;
            
            //System.out.println("#*ETP.gPS "+w+"*"+h+"\txs "+min_x_stretch+"\tys "+min_y_stretch+"\tsep "+bounding_box_separation+"\tfarthest "+farthest_label);
            
            return new Dimension(Math.max(w,80),Math.max(h,80)); // werfuse to draw something smaller than 80x80. no matter the magnification
        } else
            return user_preferred;
    }
    
}
