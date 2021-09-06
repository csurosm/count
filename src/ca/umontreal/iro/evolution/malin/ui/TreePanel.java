/*
 * TreePanel.java
 *
 * Created on November 14, 2007, 12:17 AM
 */

package ca.umontreal.iro.evolution.malin.ui;


import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Rectangle;

import java.awt.geom.AffineTransform;

import javax.swing.Scrollable;

import java.util.List;
import java.util.Hashtable;

import java.io.File;

import ca.umontreal.iro.evolution.TreeNode;

import ca.umontreal.iro.evolution.malin.IndexedPoint;
import ca.umontreal.iro.evolution.malin.SimplePointSet;

/**
 * JPanel capable of displaying a phylogenetic tree.
 * Each tree node has an assigned IndexedPoint, which is used
 * to map points selected on the screen to the nodes. The indexing of the tree
 * nodes here, used for the display/selection purposes, is independent
 * from the TreeNode objects own getId().
 *
 * @author  csuros
 */
public class TreePanel extends PointSetPanel implements Scrollable
{
    /**
     * root node for underlying main phylogeny
     */
    private TreeNode root;
  
    /**
     * Number of tree nodes, computed once for convenience.
     */
    protected int number_of_leaves;

    ///**
    // * Number of terminal tree nodes, computed once for convenience.
    // */
    //protected int number_of_nodes;
    
    /**
     * Array of tree nodes. This array gives a post-order traversal of the tree. 
     */
    protected TreeNode[] node;

    /**
     * Inverse index from tree nodes to their indices in the node[] array
     */
    private Hashtable<TreeNode,Integer> display_node_index; // may or may not be the same as getId()
    
    /**
     * A boolean value for each tree node to track selections; 
     * indexing as in <var>node[]</var>.
     */
    protected boolean[] selected;

    /*
     * Array of locations in <q>tree layout space</q>, 
     * an arbitrary fixed coordinate system 
     * (defined by the tree layout algorithm).
     * Indexing as in <var>node[]</var>.
     */
    protected IndexedPoint[] node_location;

    /**
     * Array of displayed locations, obtained by mapping from <var>node_location</var> to the display.
     * Indexing as in <var>node[]</var>.
     */
    protected IndexedPoint[] displayed_node_location;

    /**
     * Quick lookup of styles for selected/unselected nodes. 
     * displayed_node_style[idx] is a two-element array of selected/unselected styles.
     *
     * Indexing as in <var>node[]</var>.
     */
    protected PointDisplay[][] displayed_node_style;
    
    /**
     * Constant for accessing <var>displayed_node_style</var>.
     */
    protected static final int UNSELECTED_NODE=0;
    /**
     * Constant for accessing <var>displayed_node_style</var>.
     */
    protected static final int SELECTED_NODE=1;
    
    /**
     * Style for selected leaf.
     */
    private PointDisplay selected_leaf_display;
    /**
     * Style for unselected leaf.
     */
    private PointDisplay leaf_display;
    /**
     * Style for selected inner node.
     */
    private PointDisplay selected_node_display;
    
    /**
     * Style for unselected inner node.
     */
    private PointDisplay node_display;
    
    /**
     * Size of the node displays.
     */
    public static int POINT_SIZE = 8;

    public static Color unselected_leaf_color = Color.getHSBColor((float)(92.0/360.0), 0.82f, 0.64f); // Fern
    public static Color unselected_node_color = Color.getHSBColor((float)(138.0/360.0), 0.79f, 0.48f); // Pine
    public static Color selected_leaf_color = Color.getHSBColor((float)(92.0/360.0), 0.82f, 0.64f);  // Fern
    public static Color selected_node_color = Color.getHSBColor((float)(138.0/360.0), 0.79f, 0.48f); // Pine
  
    public  Color EDGE_COLOR = Color.getHSBColor((float)(28.0/360.0), 0.62f, 0.36f); // Raw Sienna
  
//    public double SUNNY_MIN_ANGLE = 0; // degrees
//    public double SUNNY_MAX_ANGLE = 180.0; // degrees

    /**
     * Dimensions of the tree layout space.
     */
    protected double[] tree_size;
    
    /**
     * Indexing constant (longest path from the root to a leaf) used with <var>tree_size</var>
     */
    protected static final int DEPTH=0;

    /**
     * Indexing constants used with <var>tree_size</var>
     */
    protected static final int BREADTH=1;
    
    public enum LayoutStyle {CLADOGRAM, PHENOGRAM, SCALED_PHENOGRAM};

    /**
     * layout style used with this tree
     */
    private LayoutStyle layout_style=null;
    
    /**
     * Most general 
     * constructor with complete specification of selection models and layout.
     */
    public TreePanel(TreeNode root, LayoutStyle layout_style, boolean node_selection_allowed, boolean area_selection_allowed)
    {
        super(node_selection_allowed, area_selection_allowed);
        this.layout_style=layout_style;
        setRoot(root);
        setLayout(layout_style);
    }

    /**
     * Constructor with complete specification of selection models.
     */
    public TreePanel(TreeNode root, boolean node_selection_allowed, boolean area_selection_allowed)
    {
        this(root, LayoutStyle.CLADOGRAM, node_selection_allowed, area_selection_allowed);
    }

    /**
     * Constructor with the same node and area selection setting.
     */
    public TreePanel(TreeNode root, boolean node_selection_allowed)
    {
        this(root, node_selection_allowed, node_selection_allowed);
    }
    
    /**
     * Instantiates with the underlying phylogeny (specified by its root),
     * with selection allowed.
     *
     * @param root the root node of the phylogeny
     */
    public TreePanel(TreeNode root)
    {
        this(root,true,true);
    }
  
    /**
     * Sets the phylogeny. Setting the root of course triggers 
     * the complete reorganization of the underlying data structures, and 
     * the display.
     *
     * @param root the root node of the phylogeny
     */
    public void setRoot(TreeNode root)
    {
        this.root=root;
        initDataStructures();
        initComponents();
        setupNodeLocations();
    }
    
    public void setLayout(LayoutStyle layout_style)
    {
        this.layout_style=layout_style;
        setupNodeLocations();
    }

    /**
     * Returns the current tree's root
     */
    public TreeNode getRoot()
    {
        return root;
    }
  
    /**
     * Each TreePanel has an associated  
     * file, which is the place from which the phylogeny was loaded.
     */ 
    private File phylo_path;
    
    /**
     * Sets the associated file 
     */
    public void setAssociatedFile(File F)
    {
        this.phylo_path = F;
    }
    
    /**
     * Returns the associated file
     */
    public File getAssociatedFile()
    {
        return phylo_path;
    }
    
    /**
     * Allocates space for the used data structures.
     */
    protected void initDataStructures()
    {
        // init node[] 
        TreeNode[] dft = root.getTraversal().getDFT();
        number_of_leaves = root.getTraversal().getLeaves().length;
        node = dft;

        // init display_node_index
        display_node_index = new Hashtable<TreeNode, Integer>();
        for (int i=0; i<node.length; i++)
        {
            display_node_index.put(node[i], new Integer(i));
        }
        
        // init node_location[]
        node_location=new IndexedPoint[node.length];
        displayed_node_location=new IndexedPoint[node.length];
        for (int i=0; i<node.length; i++)
        {
            node_location[i]=new IndexedPoint(i);
            displayed_node_location[i]=new IndexedPoint(i);
        }
        
        // init selected[]
        selected=new boolean[node.length];
        
        // init subtree weights 
        subtree_weight = new int[node.length]; // subtree_weight is used to set X coordinates
        for (int i=0; i<node.length; i++)
        {
            TreeNode N = node[i];
            subtree_weight[i]=1;
            if (!N.isLeaf())
            {
                for (int child_idx=0; child_idx<N.getNumChildren(); child_idx++)
                {
                    TreeNode C = N.getChild(child_idx);
                    int Cidx = getDisplayNodeIndex(C);
                    subtree_weight[i] += subtree_weight[Cidx];
                }
            }
        } // for each node
        
        tree_size=new double[2];
    }

    /**
     * Initializes the few graphical objects used by TreePanel.
     */
    private void initComponents()
    {
        setSelectionAreaColor(Color.pink);
        
        setCloseRadius(POINT_SIZE);
        
        leaf_display=new BoxIcon(POINT_SIZE,true);
        leaf_display.setDrawColor(unselected_leaf_color);
        leaf_display.setFillColor(getBackground());
        selected_leaf_display=new BoxIcon(POINT_SIZE,true);
        selected_leaf_display.setDrawColor(selected_leaf_color);
        selected_leaf_display.setFillColor(selected_leaf_color);
        node_display=new DiamondIcon(POINT_SIZE,true);
        node_display.setDrawColor(unselected_node_color);
        node_display.setFillColor(getBackground());
        selected_node_display=new DiamondIcon(POINT_SIZE,true);
        selected_node_display.setDrawColor(selected_node_color);
        selected_node_display.setFillColor(selected_node_color);

        displayed_node_style=new PointDisplay[node.length][2];
        for (int i=0; i<node.length; i++)
        {
            int didx = getDisplayNodeIndex(node[i]);
            if (node[i].isLeaf())
            {
                displayed_node_style[didx][SELECTED_NODE]=selected_leaf_display;
                displayed_node_style[didx][UNSELECTED_NODE]=leaf_display;
            } else 
            {
                displayed_node_style[didx][SELECTED_NODE]=selected_node_display;
                displayed_node_style[didx][UNSELECTED_NODE]=node_display;
            }
        } // for i
    }

    /**
     * Returns the index of a tree node into <var>node[]</node> array, which 
     * is the same index for all auxiliary arrays.
     */
    protected int getDisplayNodeIndex(TreeNode N)
    {
        return ((Integer)display_node_index.get(N)).intValue();
    }

  
    /**
     * Computes the node coordinates in <q>tree space</q>: for display,
     * calculateDisplayedNodeLocations() is called, which figures out the coordinates 
     * in graphics space.
     *
     * Does nothing if current layout style is null.
     */
    protected void setupNodeLocations()
    {
        tree_size[DEPTH]=tree_size[BREADTH]=0.;
        if (layout_style == LayoutStyle.CLADOGRAM)
            cladogram();
        else if (layout_style == LayoutStyle.PHENOGRAM)
            phenogram();
        else if (layout_style == LayoutStyle.SCALED_PHENOGRAM)
            scaled_phenogram();
    }
  
    /**
    * Needed in tree layout.
    */
    private int[] subtree_weight;
  
    private static final double CLADOGRAM_SLANT = 0.75;

    /**
     * Tree layout with equispaced leaf placement in cladogram style.
     * X cordinates are placed proportionally by subtree size.
     * Y coordinates are placed so that root is at 0, children always have a larger Y
     * cordinate than their parents, and all leaves are 
     * at the same Y coordinate. Y coordinate of bifurcating inner nodes is set
     * so that all branches have the same slant.
     */
    private void cladogram()
    {
    
        // initial node placement: leaves are at Y=0, parents at Y>0
        int leaf_idx=0;
        for (int i=0; i<node.length; i++)
        {
            TreeNode N = node[i];
            if (N.isLeaf())
            {
                node_location[i].setLocation(leaf_idx,0.0);
                tree_size[BREADTH]=leaf_idx;
                leaf_idx++;
            } else
            {
                double nx = 0.0;
                double ny = 0.0;
                if (N.getNumChildren()==2)
                {
                    TreeNode C_left = N.getChild(0);
                    TreeNode C_right = N.getChild(1);
                    int Cidx_left = getDisplayNodeIndex(C_left);
                    int Cidx_right = getDisplayNodeIndex(C_right);
                    double cx_left = node_location[Cidx_left].getX();
                    double cy_left = node_location[Cidx_left].getY();
                    double cx_right = node_location[Cidx_right].getX();
                    double cy_right = node_location[Cidx_right].getY();
                    
                    double xdiff = Math.abs(cx_left-cx_right);
                    double ydiff = Math.abs(cy_left-cy_right);
                    
                    double h = 0.5*(xdiff-CLADOGRAM_SLANT*ydiff)/CLADOGRAM_SLANT;
                    
                    ny = Math.max(cy_left,cy_right)+h;
                    nx = cx_left + (ny-cy_left)*CLADOGRAM_SLANT;
                } else
                {
                    for (int child_idx=0; child_idx<N.getNumChildren(); child_idx++)
                    {
                        TreeNode C = N.getChild(child_idx);
                        int Cidx = getDisplayNodeIndex(C);
                        double cx = node_location[Cidx].getX();
                        double cy = node_location[Cidx].getY();
                        nx += subtree_weight[Cidx]*cx;
                        ny = Math.max(cy+C.getLength(), ny);
                    }
                    nx /= subtree_weight[i]-1.0;
                }
                node_location[i].setLocation(nx,ny);
                tree_size[DEPTH] = Math.max(tree_size[DEPTH],ny);
            }
        }
    

        // Reposition the Y cordinates so that root is at Y=0 and children are at Y>0 
        for (int i=0; i<node.length; i++)
        {
            TreeNode N = node[i];
            double ny = node_location[i].getY();
            double nx = node_location[i].getX();
            node_location[i].setLocation(nx, node_location[node.length-1].getY()-ny);
        }
    }
    
    /**
     * Tree layout with equispaced leaf placement in phenogram style.
     * X cordinates are placed proportionally by subtree size.
     * Y coordinates are placed so that root is at 0, children always have a larger Y
     * cordinate than their parents, and all leaves are 
     * at the same Y coordinate.
     */
    private void phenogram()
    {
        // initial node placement: leaves are at Y=0, parents at Y>0
        int leaf_idx=0;
        for (int i=0; i<node.length; i++)
        {
            TreeNode N = node[i];
            if (N.isLeaf())
            {
                node_location[i].setLocation(leaf_idx,0.0);
                tree_size[BREADTH]=leaf_idx;
                leaf_idx++;
            } else
            {
                // X coordinate at the median
                int cm1 = (N.getNumChildren()-1)/2;
                int cm2 = N.getNumChildren()/2;
                double nx = 0.5*(node_location[getDisplayNodeIndex(N.getChild(cm1))].getX()
                            + node_location[getDisplayNodeIndex(N.getChild(cm2))].getX());
                // Y coordinate is one larger than the children's max
                double ny = 0.0;
                for (int child_idx=0; child_idx<N.getNumChildren(); child_idx++)
                {
                    TreeNode C = N.getChild(child_idx);
                    int Cidx = getDisplayNodeIndex(C);
                    //double cx = node_location[Cidx].getX();
                    double cy = node_location[Cidx].getY();
                    //nx += cx;
                    ny = Math.max(cy+1.0, ny);
                }
                //nx /= N.getNumChildren();
                
                
                
                //if (N.getNumChildren()%2==1)
                //{
                //    // odd number of children: put x coordinate at the middle child
                //    nx = node_location[getDisplayNodeIndex(N.getChild(N.getNumChildren()/2))].getX();
                //}
                node_location[i].setLocation(nx,ny);
                tree_size[DEPTH] = Math.max(tree_size[DEPTH],ny);
            } // not leaf
        } // for all nodes

        // Reposition the Y cordinates so that root is at Y=0 and children are at Y>0 
        for (int i=0; i<node.length; i++)
        {
            TreeNode N = node[i];
            double ny = node_location[i].getY();
            double nx = node_location[i].getX();
            node_location[i].setLocation(nx, node_location[node.length-1].getY()-ny);
        }
    }
    

    /**
     * Default implementation of edge length calculations: 
     * Math.max(N.getLength(),0.).
     */
    protected double getDisplayEdgeLength(TreeNode N)
    {
        //System.out.println("#**TP.gDEL "+N.newickName()+"\t"+N.getLength());
        return Math.max(0.0,N.getLength());
    }
    
    /**
     * Tree layout with equispaced leaf placement in phenogram style.
     * X cordinates are placed proportionally by subtree size.
     * Y coordinates are placed so that root is at 0, children always have a larger Y
     * cordinate than their parents, and the difference between
     * a node N's and its parent's Y coordinate equals the edge length getDisplayEdgeLength(N). 
     */
    private void scaled_phenogram()
    {
        // initial node placement: leaves are at Y=0, parents at Y>0
        int leaf_idx=0;
        for (int i=0; i<node.length; i++)
        {
            TreeNode N = node[i];
            if (N.isLeaf())
            {
                node_location[i].setLocation(leaf_idx,0.0);
                tree_size[BREADTH]=leaf_idx;
                leaf_idx++;
            } else
            {
                //double nx = 0.0;
                double ny = 0.0;
                for (int child_idx=0; child_idx<N.getNumChildren(); child_idx++)
                {
                    TreeNode C = N.getChild(child_idx);
                    int Cidx = getDisplayNodeIndex(C);
                    //double cx = node_location[Cidx].getX();
                    double cy = node_location[Cidx].getY();
                    //nx += cx;
                    ny = Math.max(cy+getDisplayEdgeLength(C), ny);
                    //System.out.println("#**TP.sp/set "+i+"/"+N.newickName()+"\tchild "+child_idx+"/"+C.newickName()+"\tcy "+cy+"\tlen "+getDisplayEdgeLength(C)+"\tny "+ny+"\t"+C);
                    
                }
                //nx /= N.getNumChildren();
                double nx = 0.0;
                if (N.getNumChildren()%2==1)
                {
                    // odd number of children: put x coordinate at the middle child, but shift it a little bit 
                    // or else it looks too confusing with edge labels and fancy colorings
                    nx = node_location[getDisplayNodeIndex(N.getChild(N.getNumChildren()/2))].getX();
                    if (N.getNumChildren()!=1)
                    { // i.e., at least 3
                        int neighbor_idx = N.getNumChildren()/2+(((int)ny) % 2==0
                                   ? 1 : -1);
                        double nx_shifted = (0.6*nx+0.4*node_location[getDisplayNodeIndex(N.getChild(neighbor_idx))].getX());
                        //System.out.println("#*TP.p "+i+"/"+node[i].newickName()+"\t"+nx_shifted+" ["+nx+"]\t"+ny);
                        nx = nx_shifted;
                    }
                } else
                {
                    int cm1 = (N.getNumChildren()-1)/2;
                    int cm2 = N.getNumChildren()/2;
                    nx = 0.5*(node_location[getDisplayNodeIndex(N.getChild(cm1))].getX()
                                + node_location[getDisplayNodeIndex(N.getChild(cm2))].getX());
                }
                //System.out.println("#**TP.sp/initial "+i+"/"+N.newickName()+"\t("+nx+","+ny+")\t"+node_location[i]);
                node_location[i].setLocation(nx,ny);
                tree_size[DEPTH] = Math.max(tree_size[DEPTH],ny);
            } // not leaf
        } // for all nodes

        // Reposition the Y cordinates so that root is at Y=0 and children are at Y>0 
        for (int i=node.length-1; i>=0; i--)
        {
            TreeNode N = node[i];
            if (N.isRoot())
            {
                double ny = 0.0;
                double nx = node_location[i].getX();
                node_location[i].setLocation(nx,ny);
            } else
            {
                TreeNode P = N.getParent();
                int Pidx = getDisplayNodeIndex(P);
                double ny = node_location[Pidx].getY()+getDisplayEdgeLength(N);
                double nx = node_location[i].getX();
                //System.out.println("#*TP.sp/final "+i+"/"+N.newickName()+"\t("+nx+","+ny+")\t"+node_location[i]+"\tlen "+getDisplayEdgeLength(N));
                node_location[i].setLocation(nx,ny);
            }
            
        }
        
        
    }    
    
    /**
     * paintComponent for a tree calls calculateDisplayedNodeLocations(), 
     * plotEdges(), plotNodes() and plotNodeNames(), in this order.
     */
    public void paintComponent(Graphics g) 
    {
        super.paintComponent(g);
        //if (window_resized || points_moved)
        //{
        calculateDisplayedNodeLocations();
        rebuildPointSet();
        //}
        plotEdges(g);
        plotNodes(g);
        plotNodeNames(g);
    }
    
    /**
     * Called by the superclass when multiple nodes are selected by the mouse.
     * The normal behavior is that all nodes in the subtrees are set to <q>selected</q>.
     *
     * @param L list of points, may be null. 
     */
    protected void selectPoints(List<IndexedPoint> L, boolean is_Ctrl_down) 
    {
        if (!is_Ctrl_down)
            removeSelection();
        if (L!=null)
        {
            for (int i=0; i<L.size(); i++)
                selectSubtree(node[((IndexedPoint)L.get(i)).getIndex()]);
        }
        repaint();
    }
    
    /**
     * Called by the superclass when a single tree node is selected by the mouse.
     * The normal behavior is that all nodes in the subtree are set to <q>selected</q>.
     */
    protected void selectPoint(IndexedPoint P, int num_mouse_clicks, boolean isCtrlDown) 
    {
        if (!isCtrlDown)
            removeSelection();
        selectSubtree(node[P.getIndex()]);
        repaint();
    }

    /**
     * Sets all tree nodes to <q>unselected</q> and repaints the tree
     */
    protected void removeSelection() 
    {
        //System.out.println("#**TP.rS");
        for(int i=0; i<node.length; i++)
            selected[i]=false;
        repaint();
    }
    
    /**
     * Sets the nodes in the subtree to <q>selected</q>.
     */
    protected void selectSubtree(TreeNode subtree_root)
    {
        TreeNode[] dft = subtree_root.getTraversal().getDFT();
        for (int j=0; j<dft.length; j++)
            selected[getDisplayNodeIndex(dft[j])]=true;
    }
    
    /**
     * Called by the superclass to display the appropriate 
     * tool tip.
     *
     * @param x X-coordinate of the mouse
     * @param y Y-coordinate of the mouse
     * @param p closest point in the set (may be null)
     */
    protected String getToolTipTextForPoint(int x, int y, IndexedPoint p) 
    {
        if (p==null)
            return getToolTipText();
        else 
        {
            int i=p.getIndex();
            StringBuffer retval=new StringBuffer();
            retval.append(node[i].newickName());
            //retval.append(" // loc ");
            //retval.append(node_location[i]);
            //retval.append(" // disp ");
            //retval.append(displayed_node_location[i]);
            return retval.toString();
        }
    }

    /**
     * Called by superclass to initialize the point set.
     */
    protected void initPointSet() 
    {
        point_set=new SimplePointSet();
    }
    
    /**
     * Recomputes the display locations for the tree nodes after something moved or otherwise changed.
     */ 
    protected void rebuildPointSet() 
    {
        point_set.clear();
        for (int i=0; i<node.length; i++)
            point_set.add(displayed_node_location[i]);
    }

    /**
     * This is the transformation used to convert internal coordinates to display coordinates.
     */
    //private AffineTransform plot_transform=new AffineTransform();
  
    /**
     * White space around the tree in the JPanel
     */
    private double tree_padding = 5.0;
    
    /**
     * Font size for displaying node names.
     */
    protected int label_font_size = 14;
    
    protected AffineTransform getDisplayTransform()
    {
        int w=getWidth();
        int h=getHeight();
        double x_stretch = (w-2.*tree_padding-label_font_size*20)/tree_size[BREADTH];
        double y_stretch = (h-2*tree_padding-4*label_font_size)/tree_size[DEPTH];
        AffineTransform plot_transform=new AffineTransform();        
        plot_transform.setToIdentity();

        plot_transform.translate(tree_padding+10*label_font_size, h-tree_padding-2*label_font_size);
        plot_transform.scale(x_stretch,-y_stretch);
        return plot_transform;
    }
  
    /**
    * Calculates displayed_node_location[] and sets plot_transform.
    */
    protected void calculateDisplayedNodeLocations()
    {
        AffineTransform plot_transform=getDisplayTransform();
        plot_transform.transform(node_location,0,displayed_node_location,0,node.length);
    }

    /**
    * Draws the edges of the tree.
    */
    protected void plotEdges(Graphics g)
    {
        Color old_color = g.getColor();
        g.setColor(EDGE_COLOR);  
        for (int i=0; i<node.length; i++)
            if (!node[i].isRoot())
            {
                IndexedPoint parent_pt = displayed_node_location[getDisplayNodeIndex(node[i].getParent())];
                if (layout_style == LayoutStyle.CLADOGRAM)
                    drawThickLine(g,displayed_node_location[i],parent_pt);
                else if (layout_style == LayoutStyle.PHENOGRAM || layout_style == LayoutStyle.SCALED_PHENOGRAM)
                    drawBentLine(g,displayed_node_location[i],parent_pt);
            }
        g.setColor(old_color);
    }

//    if (false) //treeStyle == LEFT_TO_RIGHT_DENDOGRAM || treeStyle == TOP_TO_BOTTOM_DENDOGRAM)
//    {
//      IndexedPoint mid_point=new IndexedPoint(0);
//      IndexedPoint displayed_mid_point=new IndexedPoint(0);
//      for (int i=0; i<number_of_nodes; i++)
//        if (!node[i].isRoot())
//        {
//          int parent=node[i].getParent().getId();
//          mid_point.x=node_location[parent].x;
//          mid_point.y=node_location[i].y;
//          plot_transform.transform(mid_point,displayed_mid_point);
//          drawLine(g,displayed_node_location[parent],displayed_mid_point);
//          drawLine(g,displayed_mid_point,displayed_node_location[i]);
//        }
//    } else 
//    {
//    }

    /**
     * Line drawing for edges (straight thick lines)
     */
    protected void drawThickLine(Graphics g, IndexedPoint p1, IndexedPoint p2)
    {
        g.drawLine((int)p1.x,(int)p1.y,(int)p2.x,(int)p2.y);
        g.drawLine((int)p1.x-1, (int)p1.y,(int)p2.x-1,(int)p2.y);
        g.drawLine((int)p1.x+1, (int)p1.y,(int)p2.x+1,(int)p2.y);
    }

    /**
     * Line drawing for edges (bent lines)
     */
    protected void drawBentLine(Graphics g, IndexedPoint p1, IndexedPoint p2)
    {
        g.drawLine((int)p1.x,(int)p1.y,(int)p1.x,(int)p2.y);
        g.drawLine((int)p1.x,(int)p2.y,(int)p2.x,(int)p2.y);
    }

    
    /**
     * Draws a single node of the tree
     */
    protected void plotNode(Graphics g, TreeNode N)
    {
        int i = getDisplayNodeIndex(N);
        if (selected[i])
            displayed_node_style[i][SELECTED_NODE].paint(g,
              (int)displayed_node_location[i].x,
              (int)displayed_node_location[i].y);
        else
            displayed_node_style[i][UNSELECTED_NODE].paint(g,
              (int)displayed_node_location[i].x,
              (int)displayed_node_location[i].y);
    }

    /**
    * Draws the tree nodes: calls plotNode() for each node
    */
    protected void plotNodes(Graphics g)
    {
        for (int i=0; i<node.length; i++) plotNode(g, node[i]);
    }
  
    /**
     * Draws the the label of a single tree: calls plotNodeLabel with newickName()
     */ 
    protected void plotNodeLabel(Graphics g, TreeNode N)
    {
        plotNodeLabel(g,N,N.newickName());
    }
  
    /**
     * Draws the label of a tree node.
     * For leaves, the label is right above the display location, centered; 
     * for inner nodes, the label is placed to the right of the displayed node, 
     * if there is enough space (i.e., the label does not cover any displayed node.
     */
    protected void plotNodeLabel(Graphics g, TreeNode N, String label_text)
    {
        Color old_color = g.getColor();
        int node_idx = getDisplayNodeIndex(N);
        //TreeNode N = node[node_idx];
        Font label_font=null;
        if (N.isLeaf())
            label_font=new Font("Serif", Font.PLAIN, label_font_size);
        else
            label_font=new Font("Serif", Font.ITALIC, label_font_size);

        FontMetrics label_fm = g.getFontMetrics(label_font);
        int x=0, y=0;
        int w = label_fm.stringWidth(label_text);
        int h = label_font_size+2;
        if (N.isLeaf())
        {
            x = (int)displayed_node_location[node_idx].x-w/2;
            y = (int)displayed_node_location[node_idx].y-POINT_SIZE/2-3;
            g.setFont(label_font);
            g.setColor(old_color);
            g.drawString(label_text, x, y);              
        } else
        {
            w+=2;
            x = (int)displayed_node_location[node_idx].x-w/2;
            y = (int)displayed_node_location[node_idx].y+label_font_size/2-2;
            Rectangle covered_by_label = new Rectangle(x-1, y-h+3,w,h);
            List<IndexedPoint> nodes_covered = point_set.withinRectangle(covered_by_label);

            if (nodes_covered.size()==1) // this very node is there ...
            {
                g.setColor(new Color(180,180,180,160));//new Color(174,249,63,50)); //new Color(24,24,200,80)); //
                g.fillRoundRect(covered_by_label.x,covered_by_label.y,covered_by_label.width, covered_by_label.height, 10, 10);
                g.setFont(label_font);
                g.setColor(unselected_node_color);
                g.drawString(label_text, x, y);
            } 
        }
        g.setColor(old_color);
    }
  
    /**
     * Displays all the node names: calls plotNodeLabel() with each tree node.
     */
    protected void plotNodeNames(Graphics g)
    {
        for (int node_idx=0; node_idx<node.length; node_idx++)
            plotNodeLabel(g,node[node_idx]);
    }
  
    protected Dimension user_preferred=null;

    /**
     * Sets a user-preferred size, which overrules the locally computed values.
     * Reset to default by calling with null argument.
     */
    public void setPreferredSize(Dimension size)
    {
        this.user_preferred = size;
        //System.out.println("#**ISD.sPS "+size);
    }    
    
    /**
     * Computes the preferred size of this panel, unless 
     * it was set previously via setPreferredSize()
     */
    public Dimension getPreferredSize()
    {
        
        if (user_preferred==null)
        {
            int max_name_length=0;
            for (int i=0; i<node.length; i++)
            {
                TreeNode N = node[i];
                if (N.isLeaf())
                {
                    int name_length = N.getName().length();
                    if (name_length>max_name_length)
                        max_name_length = name_length;
                }
            }
            int w = Math.max(800,4*(max_name_length)*number_of_leaves*label_font_size/5);
            int h = Math.max(600,(int)(tree_size[DEPTH]*w/tree_size[BREADTH]));
            return new Dimension(w,h);
        } else
            return user_preferred;
    }

  
//  /**
//   * Node placement using an <q>organic</q> procedure.
//   */
//  private void organicTreeLayout()
//  {
//      // compute path length to farthest leaf
//      double[] node_height = new double[number_of_nodes];
//      subtree_weight = new int[number_of_nodes];
//      
//      for (int node_idx=0; node_idx<number_of_nodes; node_idx++)
//      {
//          subtree_weight[node_idx]=1;
//          TreeNode N = node[node_idx];
//          if (N.isLeaf())
//          {
//              node_height[node_idx]=0.0;
//              subtree_weight[node_idx]=1;
//          }
//          else
//          {
//              for (int child_idx=0; child_idx<N.getNumChildren(); child_idx++)
//              {
//                  TreeNode C = N.getChild(child_idx);
//                  int Cidx = C.getId();
//                  double ch = node_height[Cidx]+C.getLength();
//                  node_height[node_idx]=Math.max(node_height[node_idx],ch);
//                  subtree_weight[node_idx]+=subtree_weight[Cidx];
//              }
//          }
//          //System.out.println("#**TP.oTL weight "+node_idx+"\t"+N.newickName()+"\t"+subtree_weight[node_idx]+"\t"+node_height[node_idx]);
//      }
//      
//      // compute path length from root
//      //double[] 
//      node_depth = new double[number_of_nodes];
//      for (int node_idx=number_of_nodes-2; node_idx>=0; node_idx--) // root is at 0.0
//      {
//          TreeNode N = node[node_idx];
//          int parent_idx = N.getParent().getId();
//          node_depth[node_idx]=node_depth[parent_idx]+N.getLength();
//      }
//      
//      // layout for subtrees
//      
//      subtree_sector = new CircularSector[number_of_nodes];
//      for (int node_idx=number_of_nodes-1; node_idx>=0; node_idx--)
//      {
//          subtree_sector[node_idx]=new CircularSector();
//          subtree_sector[node_idx].setRadius(node_height[node_idx]);
//      }
//
//      // set up data structure for left and right neighbors
//      int [] left_neighbor= new int[number_of_nodes];
//      int [] right_neighbor = new int[number_of_nodes];
//      
//      // constraints for root
//      CircularSector tree_right_border = new CircularSector(0.0,0.0,node_height[number_of_nodes-1],SUNNY_MIN_ANGLE,SUNNY_MIN_ANGLE); 
//      CircularSector tree_left_border = new CircularSector(0.0,0.0,node_height[number_of_nodes-1],SUNNY_MAX_ANGLE,SUNNY_MAX_ANGLE); 
//      left_neighbor[number_of_nodes-1] = -1;
//      right_neighbor[number_of_nodes-1]= -1;
//      
//      Heap Queue = new Heap(
//        new Comparator()
//        {
//            public int compare(Object o1, Object o2)
//            {
//                TreeNode N1 = (TreeNode)o1;
//                TreeNode N2 = (TreeNode)o2;
//                double y1 = subtree_sector[N1.getId()].getCenterY();
//                double y2 = subtree_sector[N2.getId()].getCenterY();
//                return Double.compare(y1,y2);
//                //double d1 = node_depth[N1.getId()];
//                //double d2 = node_depth[N2.getId()];
//                //return Double.compare(d1,d2);
//            }
//      });
//      Queue.add(node[number_of_nodes-1]);
//      
//      while (!Queue.isEmpty())
//      {
//          TreeNode N = (TreeNode) Queue.deleteMin();
//          
//          int node_idx = N.getId();
//          System.out.println("#**TP.oTL node "+node_idx+"/"+N.newickName()+"\t"+node_depth[node_idx]+"\t"+subtree_sector[node_idx]);
//
//          // subtree sector's radius is set at this point
//          // subtree sector's center is set at this point
//          // min and max angles need to be determined based on left and right constraints
//          double min_a = SUNNY_MIN_ANGLE;
//          do 
//          {
//              CircularSector right_constraint = tree_right_border;
//              if (right_neighbor[node_idx]!=-1)
//              {
//                  right_constraint = subtree_sector[right_neighbor[node_idx]];
//              }
//              double [] min_relative_position 
//                = subtree_sector[node_idx].centerDirectionFrom(right_constraint.getCenterX(),right_constraint.getCenterY());
//              double min_l = min_relative_position[0];
//              double min_angle = min_relative_position[1];
//              min_a = subtree_sector[node_idx].intersectionSector(min_l, min_angle, right_constraint);
//              if (Double.isNaN(min_a) || min_a<SUNNY_MIN_ANGLE || min_a>180.0)
//                min_a = SUNNY_MIN_ANGLE;
//              System.out.println("#**TP.oTL min "+node_idx+"/"+N.newickName()+"\tmin_a "+min_a+"\t"+right_neighbor[node_idx]+"\t"+right_constraint);
//
//              if (right_neighbor[node_idx]!=-1 && min_a == SUNNY_MIN_ANGLE)
//              {
//                  if (subtree_sector[node_idx].getRadius()>min_l)
//                  {
//                      int second_neighbor = right_neighbor[right_neighbor[node_idx]];
//                      right_neighbor[node_idx]=second_neighbor;
//                  } else
//                      break;
//              }
//          } while (right_neighbor[node_idx] != -1 && min_a == SUNNY_MIN_ANGLE);
//          subtree_sector[node_idx].setAngleStart(min_a);
//
//          double max_a = SUNNY_MAX_ANGLE;
//          do 
//          {
//              CircularSector left_constraint = tree_left_border;
//              if (left_neighbor[node_idx]!=-1)
//                  left_constraint = subtree_sector[left_neighbor[node_idx]];
//              
//              double [] max_relative_position 
//                = subtree_sector[node_idx].centerDirectionFrom(left_constraint.getCenterX(),left_constraint.getCenterY());
//              double max_l = max_relative_position[0];
//              double max_angle = max_relative_position[1];
//              max_a = subtree_sector[node_idx].intersectionSector(max_l, max_angle, left_constraint);
//              if (Double.isNaN(max_a) || max_a>SUNNY_MAX_ANGLE || max_a<0.0)
//                max_a = SUNNY_MAX_ANGLE;
//              System.out.println("#**TP.oTL max "+node_idx+"/"+N.newickName()+"\tmax_a "+max_a+"\t"+left_neighbor[node_idx]+"\t"+left_constraint);
//
//              if (left_neighbor[node_idx] != -1 && max_a == SUNNY_MAX_ANGLE)
//              {
//                  if (subtree_sector[node_idx].getRadius()>max_l)
//                  {
//                      int second_neighbor = left_neighbor[left_neighbor[node_idx]];
//                      left_neighbor[node_idx]=second_neighbor;
//                  } else
//                      break;
//              }
//          } while (left_neighbor[node_idx]!=-1 && max_a == SUNNY_MAX_ANGLE);
//          
//          subtree_sector[node_idx].setAngleEnd(max_a);
//
//          if (!N.isLeaf())
//          {
//              // set up children's initial layout
//              this.organicLayoutChildren(node_idx, subtree_sector);
//              
//              // set left and right neighbors for children; add them to queue
//              for (int child_idx=0; child_idx<N.getNumChildren(); child_idx++)
//              {
//                  TreeNode C = N.getChild(child_idx);
//                  int Cidx = C.getId();
//                  
//                  if (child_idx==0)
//                      right_neighbor[Cidx] = right_neighbor[node_idx];
//                  else
//                      right_neighbor[Cidx] = N.getChild(child_idx-1).getId();
//                  if (child_idx==N.getNumChildren()-1)
//                      left_neighbor[Cidx] = left_neighbor[node_idx];
//                  else
//                      left_neighbor[Cidx] = N.getChild(child_idx+1).getId();
//                  
//                  Queue.add(C);
//                  System.out.println("#**TP.oTL queue add "+node_idx+"+"+child_idx+"\t"+Cidx+"/"+C.newickName()+"\t"+Queue.size());
//              } // for all children of N
//          } // not a leaf
//      } // Queue
//      
//      /*
//      // compute the playing field
//      for (int node_idx=number_of_nodes-1; node_idx>=0; node_idx--)
//      {
//          TreeNode N = node[node_idx];
//          if (node_idx==number_of_nodes-1)
//          {
//              // root
//              subtree_sector[node_idx].setCircularSector(0.,0.,node_height[node_idx],SUNNY_MIN_ANGLE,SUNNY_MAX_ANGLE);
//              node_location[node_idx].setLocation(0.0, 0.0);
//          } 
//          if (!N.isLeaf())
//          {
//              this.organicLayoutChildren(node_idx, subtree_sector);
//          } // not a leaf
//      } // node_idx
//       */
//
//      // shift nodes into first quadrant
//      double farthest_left = 0.0;
//      for (int node_idx=0; node_idx<number_of_nodes; node_idx++)
//      {
//          double nx = subtree_sector[node_idx].getCenterX();
//          double ny = subtree_sector[node_idx].getCenterY();
//
//          tree_size[DEPTH]=Math.max(tree_size[DEPTH],ny);
//          tree_size[BREADTH]=Math.max(tree_size[BREADTH],nx);
//          farthest_left = Math.min(farthest_left, nx);
//          node_location[node_idx].setLocation(nx,ny);
//      }
//      
//      tree_size[BREADTH]-=farthest_left;
//      
//      for (int node_idx=0; node_idx<number_of_nodes;node_idx++)
//      {
//          TreeNode N = node[node_idx];
//          double ny = node_location[node_idx].getY();
//          double nx = node_location[node_idx].getX();
//          node_location[node_idx].setLocation(nx-farthest_left,ny);
//      }
//  }
  
//  /**
//   * Computes the appropriate angles for the children of a node.
//   * Needs to be called in a descending order from root towards the leaves.
//   */
//  private void organicLayoutChildren(int node_idx, CircularSector[] subtree_sector)
//  {
//      TreeNode N = node[node_idx];
//      double nx = subtree_sector[node_idx].getCenterX();
//      double ny = subtree_sector[node_idx].getCenterY();
//              
//      // set locations for children
//      double diff = subtree_sector[node_idx].getAngleExtent(); // subtree_sector[node_idx] is set already
//      subtree_sector[N.getChild(0).getId()].setAngleStart(subtree_sector[node_idx].getAngleStart());
//      int cumulative_weight = 0;
//      for (int child_idx=0; child_idx<N.getNumChildren(); child_idx++)
//      {
//          TreeNode C = N.getChild(child_idx);
//          int Cidx = C.getId();
//          cumulative_weight += subtree_weight[Cidx]; 
//          double m = subtree_sector[node_idx].getAngleStart()+diff*cumulative_weight/(subtree_weight[node_idx]-1.0);    //(1.0+child_idx)/N.getNumChildren(); 
//          subtree_sector[Cidx].setAngleEnd(m);
//          if (child_idx!=N.getNumChildren()-1)
//          {
//              TreeNode S = N.getChild(child_idx+1);
//              int Sidx = S. getId();
//              subtree_sector[Sidx].setAngleStart(m);
//          }
//          // set location in the middle
//          double mid_angle = 0.5*(subtree_sector[Cidx].getAngleStart()+subtree_sector[Cidx].getAngleEnd());
//          double cx = nx+C.getLength()*CircularSector.cos(mid_angle);
//          double cy = ny+C.getLength()*CircularSector.sin(mid_angle);
//                  
//          subtree_sector[Cidx].setCenterX(cx);
//          subtree_sector[Cidx].setCenterY(cy);
//                  
//          //System.out.println("#**TP.TL loc "+Cidx+"/"+C.newickName()+"\t("+cx+", "+cy+") // "+mid_angle+" "+subtree_sector[Cidx]+"\t"+diff+"\t"+subtree_sector[node_idx]);
//
//          // OK, we have a playing field for the children
//          // given as slices based at the root
//          // refine the min-max angle values if necessary
//          if (!C.isLeaf())
//          {
//              double l = C.getLength();
//
//              double[] min_fan = subtree_sector[Cidx].intersectionHalfLine(l, mid_angle, subtree_sector[Cidx].getAngleStart());
//              if (min_fan.length==0)
//                  subtree_sector[Cidx].setAngleStart(SUNNY_MIN_ANGLE);
//              else 
//                  {
//                      double min_a = min_fan[min_fan.length-1];
//                      
//                      if (min_a<SUNNY_MIN_ANGLE || min_a>180.0)
//                          min_a = SUNNY_MIN_ANGLE;
//                      subtree_sector[Cidx].setAngleStart(min_a);
//                  }
//              double[] max_fan = subtree_sector[Cidx].intersectionHalfLine(l, mid_angle, subtree_sector[Cidx].getAngleEnd());
//              if (max_fan.length==0)
//                  subtree_sector[Cidx].setAngleEnd(SUNNY_MAX_ANGLE);
//              else
//              {
//                  double max_a = max_fan[max_fan.length-1];
//                  if (max_a>SUNNY_MAX_ANGLE || max_a<0.0)
//                      max_a = SUNNY_MAX_ANGLE;
//                  subtree_sector[Cidx].setAngleEnd(max_a);
//              }
//          } // child is not a leaf
//          //System.out.println("#**TP.oTL "+Cidx+"/"+C.newickName()+"\t"+subtree_sector[Cidx]);
//      } // for all children      
//  }

    
        public boolean getScrollableTracksViewportHeight()
        {
            return false;
        }
        
        public boolean getScrollableTracksViewportWidth()
        {
            return false;
        }
        
        public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction)
        {
            return 10;
        }
        
        public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation, int direction)
        {
            return 60;
        }
        
        /**
         * This will return whatever getPreferredSize() says.
         */
        public Dimension getPreferredScrollableViewportSize()
        {
            Dimension D= getPreferredSize();
            //System.out.println("#**ISD.gPrefSVS "+D);
            return D;
        }

    
    
//    /**
//     * Class used to store tree layout styles
//     */
//    public static class LayoutStyle
//    {
//        private LayoutStyle(){}
//    }
            
}