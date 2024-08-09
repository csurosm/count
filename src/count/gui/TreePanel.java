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
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.Rectangle;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.JRadioButton;
import javax.swing.ListSelectionModel;
import javax.swing.Scrollable;
import javax.swing.SwingConstants;
import javax.swing.event.ListSelectionEvent;

import count.ds.IndexedTree;
import count.ds.TreeTraversal;
import count.gui.kit.BoxIcon;
import count.gui.kit.DiamondIcon;
import count.gui.kit.DiscIcon;
import count.gui.kit.DrawString;
import count.gui.kit.IndexedPoint;
import count.gui.kit.MagnificationSpinner;
import count.gui.kit.PointIcon;
import count.gui.kit.PointSetPanel;
import count.io.CommandLine;
import count.io.DataFile;
import count.io.ExportableData;
import count.io.NewickParser;
import count.matek.Functions;



/** 
 * JPanel capable of displaying a phylogenetic tree.
 * Each tree node has an assigned IndexedPoint, which is used
 * to map points selected on the screen to the nodes. The indexing of the tree
 * nodes here, used for the display/selection purposes, is independent
 * from the tree nodes' own index within the tree.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 * @since November 14, 2007, 12:17 AM
 */
public class TreePanel 
					extends PointSetPanel 
					implements Scrollable, ExportableData
{
    /**
     * Color for {@link count.gui.TreePanel}.
     */
    public static Color TREE_UNSELECTED_LEAF_COLOR = Color.getHSBColor((float)(92.0/360.0), 0.82f, 0.64f); // Fern
    /**
     * Color for {@link count.gui.TreePanel}.
     */
    public static Color TREE_UNSELECTED_NODE_COLOR = Color.BLACK; //Color.getHSBColor((float)(138.0/360.0), 0.79f, 0.48f); // Pine
    /**
     * Color for {@link count.gui.TreePanel}.
     */
    public static Color TREE_SELECTED_LEAF_COLOR = Color.getHSBColor((float)(92.0/360.0), 0.82f, 0.64f);  // Fern
    /**
     * Color for {@link count.gui.TreePanel}.
     */
    public static Color TREE_SELECTED_NODE_COLOR = Color.getHSBColor((float)(138.0/360.0), 0.79f, 0.48f); // Pine
    /**
     * Color for {@link count.gui.TreePanel}.
     */
    public static Color TREE_EDGE_COLOR = Color.getHSBColor((float)(28.0/360.0), 0.62f, 0.36f); // Raw Sienna
    /**
     * Size of the nodes in {@link count.gui.TreePanel}.
     */
    public static int TREE_POINT_SIZE = 8;
    
    public static Color AREA_SELECTION_COLOR = new Color(210, 240, 255, 128);
    
    public static double TREE_MAGNIFICATION_MIN = 0.1;
    public static double TREE_MAGNIFICATION_MAX = 10.;


    public static final Color TREE_NODE_LABEL_BACKGROUND = new Color(180,180,180,180);
    public static final Color TREE_NODE_INFO_BACKGROUND = new Color(255,220,220,240);

    public static final Color TREE_EDGE_BACKGROUND = TREE_NODE_LABEL_BACKGROUND;
    
    public static int FONT_SIZE_MIN = 7; 
    public static int FONT_SIZE_MAX = 32;
    
    /**
     * Thick edge width.
     */
    public static int TREE_THICK_EDGE = 2;
	
   /**
    * Slant for the cladogram drawing.
    */
   private static final double CLADOGRAM_SLANT = 0.75;
   /**
    * Layout styles. 
    * May be (slanted) cladogram, phenogram (=rectangular cladogram) or (rectangular) phylogram. 
    * In a cladogram and a phenogram, all leaves are at the same level. In a phylogam, displayed edge lengths 
    * are informative. 
    */
   public enum LayoutStyle 
   {
	   /**
	    * Slanted edges leading to children; placement by topology.
	    */
	   CLADOGRAM, 
	   /**
	    * Bent edges leading to children; placement buy topology
	    */
	   PHENOGRAM, 
	   /**
	    * Bent edges leading to children; placement by {@link TreePanel#getDisplayEdgeLength(int)}
	    */
	   PHYLOGRAM,
	   /**
	    * Curved edges leading to children; placement by topology
	    */
	   PITCHFORK,
	   /**
	    * Node list with vertical drawing
	    */
	   NODE_TABLE
   };
   /**
    * Layout style used with this tree.
    */
   private LayoutStyle layout_style;
   
   /**
    * The tree that is displayed here. 
    */
   private final IndexedTree tree; 
   private final DataFile<? extends IndexedTree> tree_data;
   /**
    * Node placement in tree space. 
    */
   private final IndexedPoint[] node_locations;
   
   private double tree_space_width;
   private double tree_space_height;
   
   /**
    * Node placement in screen coordinates.
    */
   private final DisplayedNode[] displayed_node;

//   private final boolean[] selected;
   /**
    * Number of nodes within each subtree - calculated once on instantiation.
    */
   private final int[] subtree_sizes;
   /**
    * Preferred user dimensions, if set.
    */
   private Dimension user_preferred=null;
   
   /**
    * Flag to indicate that node placement needs to be recalculated.
    */
   private boolean valid_bounding_boxes;
   
   
   /**
    * Magnification on the display. 
    */
   private double magnification = 1.0;
   
   /**
    * Property change fired when magnification changes. 
    */
   public static final String MAGNIFICATION_PROPERTY = "magnification";
   

   /** 
    * Font size for 100% magnification.  
    */
   private int normal_label_font_size;
   
   /** 
    * Font size for displaying node names.
    */
   private int label_font_size;

   
   /**
    * Padding around the tree. 
    */
   private Insets padding;
   
   private int root_stem_length;
   
   private int bounding_box_separation;
   
   
   private final static float[] dash = {2.0f, 4.0f};
   private final static Stroke dashed_stroke = new BasicStroke(.5f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 10.0f, dash, 0);
   private final static Stroke thick_dashed_stroke = new BasicStroke(TREE_THICK_EDGE, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 10.0f, dash, 0);
   private final static Stroke thick_stroke =  new BasicStroke(TREE_THICK_EDGE);
   
   
   /**
    * Instantiation.
    * 
    * 
    * @param tree_data Tree to be displayed. Edge lengths will be accessed through the {@link IndexedTree#getLength(int) } interface. 
    * @param layout_style
    * @param selection_mode one of the mode constants from {@link javax.swing.ListSelectionModel}
    */
   public TreePanel(DataFile<? extends IndexedTree> tree_data, LayoutStyle layout_style, int selection_mode)
   {
       super(selection_mode, true, true); // with area selection
       this.tree_data = tree_data;
       this.tree = tree_data.getContent();
       
//       this.node_selection = new javax.swing.DefaultListSelectionModel();
//       node_selection.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
       
       this.layout_style = layout_style;
       int num_nodes = tree.getNumNodes();
       node_locations = new IndexedPoint[num_nodes];
       subtree_sizes = TreeTraversal.getSubtreeSizes(tree);
       displayed_node = new DisplayedNode[num_nodes];

       int tree_point_size = getTreePointSize();
       setNormalFontSize(getTreeNormalFontSize(tree_point_size));       
       int pad = 2*getTreeDefaultPadding(tree_point_size);
       setPadding(new Insets(pad, pad, pad, LayoutStyle.NODE_TABLE==layout_style?0:10*label_font_size));
       setRootStem(getTreeDefaultStem(tree_point_size));
       setBoundingBoxSeparation(LayoutStyle.NODE_TABLE==layout_style?1:getTreeBBoxSperation(tree_point_size));
       
       calculateNodeLocations();
       initializeGraphicalElements();
       setValidBoundingBoxes(false);
   }

   public TreePanel(DataFile<? extends IndexedTree> tree_data)
   {
       this(tree_data, 
    		   (tree_data.getContent().getNumLeaves()<60?
    				   LayoutStyle.PITCHFORK:LayoutStyle.NODE_TABLE)
//    				   (tree_data.getContent().hasLength()
//    				   ?LayoutStyle.PHYLOGRAM
//    						   :LayoutStyle.PITCHFORK):LayoutStyle.NODE_TABLE), 
    		   , ListSelectionModel.SINGLE_SELECTION);
   }
   
   public TreePanel(IndexedTree phylo, LayoutStyle layout_style, int selection_mode)
   {
	   this(new DataFile<>(phylo),layout_style,selection_mode);
   }
   public DataFile<? extends IndexedTree> getTreeData()
   {
	   return this.tree_data;
   }
   
   /**
    * To be shown while browsing TreePanels.
    * 
    * @return
    */
   public String getTreeName()
   {
	   return DataFile.chopFileExtension(tree_data.getFile().getName());
   }
   
   public LayoutStyle getTreeLayoutStyle()
   {
	   return layout_style;
   }
   
   /**
    * Resets to default node positions and recalculates 
    * label bounding boxes. 
    * 
    * @param layout
    */
   public void setTreeLayoutStyle(LayoutStyle layout)
   {
       this.layout_style = layout;
       if (LayoutStyle.NODE_TABLE==layout_style)
       {
    	   setBoundingBoxSeparation(1);
    	   padding.right = padding.left;
       } else
       {
           setBoundingBoxSeparation(getTreeBBoxSperation(getTreePointSize()));
           padding.right = padding.left + 10*label_font_size;
       }
       for (int node=0; node<tree.getNumNodes(); node++)
    	   displayed_node[node].setDefaultPosition();
       this.calculateNodeLocations();
       this.revalidate();
       this.repaint();
   }
   
   @Override
   public void setPreferredSize(Dimension D)
   {
	   this.user_preferred = D;
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
    	   if (LayoutStyle.NODE_TABLE == layout_style)
    	   {
    		   double phylo_width = tree_space_height*getTreePointSize();
    		   
    		   // rightmost ending node label
    		   double label_width = 0.0;
    		   for (int node=0; node<tree.getNumNodes(); node++)
    		   {
    			   Rectangle2D nodeBB = displayed_node[node].node_bounding_box;
    			   double node_ends = nodeBB.getMaxX();
    			   if (node_ends>label_width) label_width=node_ends;
    		   }
    		   // Y stretch
    		   double min_y_stretch = label_font_size;
    		   int[] postorder = TreeTraversal.postOrder(tree);
    		   assert (postorder[0]==0);
    		   Rectangle2D prevBB = displayed_node[0].node_bounding_box;
    		   for (int i=1; i<tree.getNumNodes(); i++)
    		   {
    			   int node = postorder[i];
    			   Rectangle2D currentBB = displayed_node[node].node_bounding_box;
    			   double prev_ends = prevBB.getMaxY();
    			   double current_starts = currentBB.getMinY();
    			   double current_stretch = prev_ends-current_starts;
    			   if (current_stretch>min_y_stretch)
    				   min_y_stretch = current_stretch;
    			   prevBB = currentBB;
    		   }
    		   min_y_stretch += bounding_box_separation;
    		   
//    		   System.out.println("#**TP.gPS "+layout_style+"\tystretch "+min_y_stretch+"\tlfont "+label_font_size+"\t");
    		   
    		   double minimum_height 
    		   			 = padding.left 
//    		   			 - displayed_node[0].node_bounding_box.getMinY()
    		   			 + min_y_stretch * (tree.getNumNodes()-1.0)
//    		   			 + prevBB.getMaxY()
    		   			 + padding.right;
    		   double minimum_width 
    		   			= padding.top 
    		   				+ root_stem_length 
    		   				+ phylo_width
    		   				+ bounding_box_separation
    		   				+ label_width
    		   				+ padding.bottom;
    		   
    		   int w = (int) minimum_width;
    		   int h = (int) minimum_height;

//    		   System.out.println("#**TP.gPS "+layout_style+"\tystretch "+min_y_stretch+"\tlfont "+label_font_size+"\tw "+w+"\th "+h);
    		   
	           int tiny = 2*(padding.top + root_stem_length + padding.bottom);
	           return new Dimension(Math.max(w,tiny),Math.max(h,tiny)); // refuse to draw something too small, no matter the magnification
    	   } else
    	   {
	           // figure out minimum stretch for X coordinates
	           double min_x_stretch = 0.0;
	           int num_leaves = tree.getNumLeaves();
	           for (int leaf=1; leaf<num_leaves; leaf++)
	           {
	               final Rectangle2D leaf_BB = displayed_node[leaf-1].node_bounding_box;
	               double previous_leaf_ends = leaf_BB.getX()+label_font_size; //leaf_BB.getWidth();
	               final Rectangle2D edge_BB = displayed_node[leaf-1].edge_bounding_box;
	               double previous_edge_ends = edge_BB.getX()+edge_BB.getWidth();
	               double previous_ends = Math.max(previous_leaf_ends, previous_edge_ends);
	
	               double current_leaf_starts = displayed_node[leaf].node_bounding_box.getX();
	               double current_edge_starts = displayed_node[leaf].edge_bounding_box.getX();
	               double current_starts = Math.min(current_leaf_starts, current_edge_starts);
	
	               double current_min_stretch = previous_ends - current_starts;
	               if (current_min_stretch>min_x_stretch)
	                   min_x_stretch = current_min_stretch;
	           }
	           min_x_stretch += bounding_box_separation;
	           
	           //min_x_stretch *= magnification;
	           
	           final int last_leaf_idx = num_leaves -1;
	           double minimum_width 
	                   = padding.left 
	                   - displayed_node[0].node_bounding_box.getX()
	                   + min_x_stretch * (num_leaves-1.0)
	                   + displayed_node[last_leaf_idx].node_bounding_box.getX()
	                   + displayed_node[last_leaf_idx].node_bounding_box.getWidth()
	                   + padding.right;
	           double reasonable_width
	           		= minimum_width + 2*label_font_size * (num_leaves-1.0);
	           
	           int num_nodes = tree.getNumNodes();
	           double min_y_stretch = 0.0;
	           for (int node=0; node<num_nodes-1; // skip root at index num_nodes-1
	                   node++) 
	           {
	               int parent = tree.getParent(node);
	               Rectangle2D cR = displayed_node[node].node_bounding_box;
	               Rectangle2D eR = displayed_node[node].edge_bounding_box;
	               Rectangle2D pR = displayed_node[parent].node_bounding_box;
	                       
	               double eh = eR.getHeight();
	               double current_stretch = (cR.getY()+cR.getHeight()+eh-pR.getY());
	
	               int x = (eh==0.?0:1); // counting how many non-empty bounding boxes need to be separated here
	               if (eh!=0. && pR.getHeight()!=0.)
	               {
	                   x++;
	               }
	               current_stretch += x * bounding_box_separation;
	
	               double diff_y = node_locations[node].getY()-node_locations[parent].getY();
	               
//	               if (diff_y == 0.0)
//	               {
////	            	   System.out.println("#**TP.gPS zero edge "+parent+"->"+node+"\t"
////	            			   +node_locations[node]+"\t"+node_locations[parent]
////	            			   +"\tlen "+tree.getLength(node)+"/"+getDisplayEdgeLength(node));
//	               }
	               
	               current_stretch /= diff_y;
	                       
	               if (current_stretch > min_y_stretch)
	                   min_y_stretch = current_stretch;
	           }
	           
	           min_y_stretch *= getMagnification();
	           
	           double farthest_label = 0.0;
	           {
	               for (int node=0; node<num_nodes; node++)
	               {
	                   Rectangle2D cR = displayed_node[node].node_bounding_box;
	                       
	                   double label_y = node_locations[node].getY()*min_y_stretch-cR.getY();//+cR.getWidth();
	                   if (label_y > farthest_label)
	                       farthest_label = label_y;
	               }
	           }
	           double minimum_height = padding.top + root_stem_length + padding.bottom 
	                   + farthest_label;
	
	           double reasonable_height = 300*magnification;
	           //int w = (int)(Math.max(minimum_height,reasonable_height)/minimum_height*minimum_width);
	           int h = (int) Math.max(minimum_height, reasonable_height);
	           
	           int w = (int)Math.max(minimum_width,3*label_font_size*num_leaves);
	           //int h = (int)minimum_height;
	           
//	           System.out.println("#*TP.gPS "+w+"*"+h+"\txs "+min_x_stretch+"\tys "+min_y_stretch
//	        		   +"\tsep "+bounding_box_separation+"\tfarthest "+farthest_label
//	        		   +"\ttsp "+tree_space_width+"*"+tree_space_height
//	        		   //+"\t"+Arrays.toString(node_locations)
//	        		   );
	           
	           int tiny = 2*(padding.top + root_stem_length + padding.bottom);
	           return new Dimension(Math.max(w,tiny),Math.max(h,tiny)); // refuse to draw something too small, no matter the magnification
    	   }
       } else
       return user_preferred;
   }
      

   @Override
   protected void paintComponent(Graphics g) 
   {
       super.paintComponent(g); // as PointSetPanel, draws selection area 
       
//       System.out.println("#*TP.pC "+g);
       
       if (!valid_bounding_boxes)
       {
    	   this.calculateBoundingBoxes(g);
           setValidBoundingBoxes(true);
           revalidate(); // and repaint after
       } 
       AffineTransform plot_transform = getDisplayTransform();
       plot_transform.transform(node_locations,0,displayed_node,0,node_locations.length);
       
       Color old_color = g.getColor();
       g.setColor(TREE_EDGE_COLOR);
       paintEdges(g);
       g.setColor(old_color);
       paintNodes(g);
       paintNodeLabels(g);
   }    
   
   
   protected void paintNodes(Graphics g)
   {
       Graphics2D g2 = (Graphics2D)g.create();
       for (int node=0; node<tree.getNumNodes(); node++)
		   displayed_node[node].paintNode(g2);
   }
   
   protected void paintEdges(Graphics g)
   {
       Graphics2D g2 = (Graphics2D)g.create();
       for (int node=tree.getNumLeaves(); node<tree.getNumNodes(); node++)
       {
		   displayed_node[node].paintEdges(g2);
       }
   }
   
   protected void paintNodeLabels(Graphics g)
   {
       Graphics2D g2 = (Graphics2D)g.create();
       for (int node=0; node<tree.getNumNodes(); node++)
    	   if (!isSelected(node))
    		   displayed_node[node].paintLabel(g2);
       
       for (int node=0; node<tree.getNumNodes(); node++)
    	   if (isSelected(node))
    		   displayed_node[node].paintLabel(g2);
	   
   }
   
   protected int getTreePointSize() {return TREE_POINT_SIZE;}
   
   /**
    * Default white space amount around the tree. 
    */
   public static int getTreeDefaultPadding(int tree_point_size){return (5*tree_point_size)/8;}
   
   public static int getTreeBBoxSperation(int tree_point_size){ return (1*tree_point_size)/2;}
   
   public static int getTreeNormalFontSize(int tree_point_size){ return (3*tree_point_size)/2;}
   
   public static int getTreeNormalFontSize() {return getTreeNormalFontSize(TREE_POINT_SIZE);}
   
   public static int getTreeDefaultStem(int tree_point_size){ return 0;}
   
   /**
    * Padding around the tree in the panel. 
    * 
    * @param I the new padding values
    */
   protected final void setPadding(Insets I)
   {
       this.padding = I;
   }
   
   protected final Insets getPadding()
   {
	   return this.padding;
   }
   
   /**
    * Sets the length of the little stem leading to the root.
    * 
    * @param len stem length in pixels
    */
   protected final void setRootStem(int len)
   {
       root_stem_length = len;
   }
   
   /**
    * Separation between bounding boxes for nodes. 
    * @param d new value
    */
   protected final void setBoundingBoxSeparation(int d)
   {
       this.bounding_box_separation = d;
   }
   
   protected final int getTreeBoundingBoxSeparation()
   {
	   return this.bounding_box_separation;
   }
   
   protected final double getTreeHeight()
   {
	   return tree_space_height;
   }
   
   protected final double getTreeWidth()
   {
	   return tree_space_width;
   }
   
   @Override
   protected String getTooltipText(int x, int y, IndexedPoint p) 
   {
       if (p==null)
           return getToolTipText();
       else 
       {
           int i=p.getIndex();
           return tree.getIdent(i);
       }
   }
   
   /**
    * A Spinner attached to the magnification value: 
    * changing the spinner calls {@link #setMagnification(double) }.
    * 
    * @return 
    */
   public MagnificationSpinner createMagnificationSpinner()
   {
	   MagnificationSpinner magnification_spin = new MagnificationSpinner();
       magnification_spin.addChangeListener(event->{
               double value = ((Number)magnification_spin.getValue()).doubleValue();
               if (value!=getMagnification())
               {
                   setMagnification(value);
               }
           });
   
       addPropertyChangeListener(event->{
           if (MAGNIFICATION_PROPERTY.equals(event.getPropertyName()))
           {
               double current_mag = getMagnification();
               double spin_value =((Number)magnification_spin.getValue()).doubleValue();
               if (spin_value != current_mag)
               {
                   magnification_spin.setValue(current_mag);
               }
           }
       });
   
       return magnification_spin;
   }
   
   @Override
   public void valueChanged(ListSelectionEvent e)
   {
       if (!e.getValueIsAdjusting())
       {
           int idx = getSelectionModel().getLeadSelectionIndex();
           if (idx != -1)
           {
               double nx = displayed_node[idx].getX();
               double ny = displayed_node[idx].getY();
               int r = (int)(getNodeRadius()+0.5);
               Rectangle nRect = new Rectangle((int)(nx-r),(int)(ny-r),2*r,7*r);
               scrollRectToVisible(nRect);
           }
       }
       super.valueChanged(e);
   }
   

   public int getNodeRadius()
   {
       int font_size = getLabelFontSize();
       int r = Math.max(7*font_size, 7*font_size*font_size/getNormalFontSize())/5;
       return r;
   }
   
   
   /**
    * A little box for choosing the layout, with the attached listeners to update the tree layout.
    * 
    * @return 
    */
   public Box createLayoutChooser()
   {
       String apple_button_size = "small";
       
       JRadioButton pitchfork = new JRadioButton("curved");
       pitchfork.putClientProperty("JComponent.sizeVariant", apple_button_size);        
       pitchfork.setSelected(TreePanel.this.layout_style == LayoutStyle.PITCHFORK);
       
       JRadioButton clado = new JRadioButton("cladogram");
       clado.setSelected(TreePanel.this.layout_style == LayoutStyle.CLADOGRAM);
       clado.putClientProperty("JComponent.sizeVariant", apple_button_size);        

       JRadioButton phylo = new JRadioButton("phylogram");
       phylo.setSelected(TreePanel.this.layout_style == LayoutStyle.PHYLOGRAM);
       phylo.setEnabled(TreePanel.this.tree.hasLength());
       phylo.putClientProperty("JComponent.sizeVariant", apple_button_size);        

       JRadioButton pheno = new JRadioButton("phenogram");
       pheno.setSelected(TreePanel.this.layout_style == LayoutStyle.PHENOGRAM);
       pheno.putClientProperty("JComponent.sizeVariant", apple_button_size);        
       
       JRadioButton table = new JRadioButton("table");
       table.setSelected(TreePanel.this.layout_style == LayoutStyle.PHENOGRAM);
       table.putClientProperty("JComponent.sizeVariant", apple_button_size);        
       
       ButtonGroup layout_choices = new ButtonGroup();
       layout_choices.add(pitchfork);
       layout_choices.add(clado);
       layout_choices.add(pheno);
       layout_choices.add(phylo);
       layout_choices.add(table);
       
       pitchfork.setSelected(layout_style == LayoutStyle.PITCHFORK);
       clado.setSelected(layout_style == LayoutStyle.CLADOGRAM);
       phylo.setSelected(layout_style == LayoutStyle.PHYLOGRAM);
       pheno.setSelected(layout_style == LayoutStyle.PHENOGRAM);
       table.setSelected(layout_style == LayoutStyle.NODE_TABLE);
       
       pitchfork.addActionListener(a->setTreeLayoutStyle(LayoutStyle.PITCHFORK));
       clado.addActionListener(a->setTreeLayoutStyle(LayoutStyle.CLADOGRAM));
       phylo.addActionListener(a->setTreeLayoutStyle(LayoutStyle.PHYLOGRAM));
       pheno.addActionListener(a->setTreeLayoutStyle(LayoutStyle.PHENOGRAM));
       table.addActionListener(a->setTreeLayoutStyle(LayoutStyle.NODE_TABLE));
       
       Box B = Box.createHorizontalBox();
       B.add(pheno);
       B.add(clado);
       B.add(phylo);
       B.add(pitchfork);
       B.add(table);
       
       return B;
   }   
   
   
   /**
    * Places the nodes in tree space ({@link #node_locations}), 
    * by {@link #layout_style} logic.
    */
   protected void calculateNodeLocations()
   {
       int num_nodes = tree.getNumNodes();
       int num_leaves = tree.getNumLeaves();

       

       if (LayoutStyle.NODE_TABLE == layout_style)
	   {
		    placeNodeTable();
	   } else
	   {
	       for (int leaf_idx=0; leaf_idx<num_leaves; leaf_idx++)
	       {
	           assert (tree.isLeaf(leaf_idx)); // by IndexedTree specification
	           node_locations[leaf_idx]=new IndexedPoint(leaf_idx, leaf_idx, 0.0);
	       }
	       
	       if (LayoutStyle.PHYLOGRAM == layout_style)
	       {
	    	   calculateDisplayEdgeRange();
	       }
	       
	       for (int node_idx=num_leaves; node_idx<num_nodes; node_idx++)
	       {
	           assert (!tree.isLeaf(node_idx)); // by IndexedTree specification
	           switch (layout_style)
	           {
	               case CLADOGRAM: node_locations[node_idx]=placeNodeCladogram(node_idx); break;
	               case PHENOGRAM: 
	               case PITCHFORK:     
	                   node_locations[node_idx]=placeNodePhenogram(node_idx); break;
	               case PHYLOGRAM: node_locations[node_idx]=placeNodePhylogram(node_idx); break;
	               default: // no such
	           }
	       }
	   }
       tree_space_width = tree_space_height = 0.0;
       for (int node=0; node<num_nodes; node++)
       {
    	   tree_space_width = Double.max(tree_space_width, node_locations[node].getX());
    	   tree_space_height = Double.max(tree_space_height, node_locations[node].getY());
       }
       // actually, there is no need to check all nodes, we know which ones define the bounding box
       assert (num_leaves == 0 || tree_space_width==(layout_style==LayoutStyle.NODE_TABLE?num_nodes-1:num_leaves-1));
       assert (num_nodes==0 || tree_space_height==node_locations[tree.getRoot()].getY());
	       
       // flip the Y coordinates so that root is at Y=0 and children are at Y>0 
       switch(layout_style)
       {
           case CLADOGRAM:
           case PHENOGRAM:
           case PITCHFORK: 
           case NODE_TABLE:
               for (int node_idx=0; node_idx<num_nodes; node_idx++)
               {
                   double ny = node_locations[node_idx].getY();
                   double nx = node_locations[node_idx].getX();
                   node_locations[node_idx].setLocation(nx, tree_space_height-ny);
               }
               break;
           case PHYLOGRAM:
               int node_idx = num_nodes-1;
               node_locations[node_idx].setLocation(node_locations[node_idx].getX(), 0.0);
               while(node_idx>0)
               {
                   node_idx--;
                   int Pidx = tree.getParent(node_idx);
                   double nx = node_locations[node_idx].getX();
                   double ny = node_locations[Pidx].getY()+getDisplayEdgeLength(node_idx);
                   node_locations[node_idx].setLocation(nx, ny);
               }
       }
       setValidBoundingBoxes(false);
   }

   /**
    * Tree layout in cladogram style.
    * 
    * X cordinates are placed proportionally by subtree size.
    * Y coordinate of bifurcating inner nodes is set
    * so that all branches have the same slant.
    */
   private IndexedPoint placeNodeCladogram(int node)
   {
       double nx = 0.0;
       double ny = 0.0;
       int num_children = tree.getNumChildren(node);
       if (num_children==2)
       {
           int C_left = tree.getChild(node, 0);
           double cx_left = node_locations[C_left].getX();
           double cy_left = node_locations[C_left].getY();
           int C_right = tree.getChild(node, 1);
           double cx_right = node_locations[C_right].getX();
           double cy_right = node_locations[C_right].getY();
                   
           double xdiff = Math.abs(cx_left-cx_right);
           double ydiff = Math.abs(cy_left-cy_right);
                   
           double h = 0.5*(xdiff-CLADOGRAM_SLANT*ydiff)/CLADOGRAM_SLANT;
                   
           ny = Math.max(cy_left,cy_right)+h;
           nx = cx_left + (ny-cy_left)*CLADOGRAM_SLANT;
       } else
       {
           for (int child_idx=0; child_idx<num_children; child_idx++)
           {
               int C = tree.getChild(node, child_idx);
               double cx = node_locations[C].getX();
               double cy = node_locations[C].getY();
               nx += subtree_sizes[C]*cx;
               ny = Math.max(cy+1.0, ny);
//               ny = Math.max(cy+tree.getLength(Cidx), ny);
           }
           nx /= subtree_sizes[node]-1.0;
       }
       return new IndexedPoint(node,nx,ny);
   }
   
   
   /**
    * Edge length for tree display, used with {@link LayoutStyle#PHYLOGRAM}.
    * (Default implementation uses {@link IndexedTree#getLength(int) }.
    * 
    * @param node_idx
    * @return 
    */
   protected double getDisplayEdgeLength(int node_idx)
   {
       return Double.max(display_edge_shortest, 
    		   Double.min(display_edge_longest, getTrueEdgeLength(node_idx)));
   }
   
   protected double getTrueEdgeLength(int node)
   {
	   return tree.getLength(node);
   }

   protected double display_edge_shortest=0.0;
   protected double display_edge_longest=0.0;
   protected int[] display_edge_exact;
   
   private static final double DISPLAY_EDGE_RANGE = 1.0/256.0;

   protected double getDisplayEdgeRangeFactor() { return DISPLAY_EDGE_RANGE;}
   private void calculateDisplayEdgeRange()
   {
	   int num_nodes = tree.getNumNodes();
	   double[] edges = new double[num_nodes];
	   int num_edges = 0;
	   for (int node=0; node<num_nodes; node++)
	   {
		   if (!tree.isRoot(node))
		   {
			   double len = getTrueEdgeLength(node);
			   if (len != 0.0 && len != Double.POSITIVE_INFINITY)
			   {
				   edges[num_edges++]=len;
			   }
		   }
	   }
	   edges = Arrays.copyOf(edges, num_edges);
	   Arrays.sort(edges);
	   double[] length_range = calculateDisplayLengthRange(edges, getDisplayEdgeRangeFactor());
	   if (length_range[0]==length_range[1]) // 0 or 1 finite edges 
	   {
		   assert Double.isFinite(length_range[0]);
		   length_range[1] += 1.0;
	   }
	   display_edge_exact = new int[2]; // used for legend 
	   display_edge_longest = length_range[1];
	   Functions.roundToMostSignificantDigit(display_edge_longest, display_edge_exact); // but don't keep the value
	   display_edge_shortest = length_range[0];
	   
//	   System.out.println("#**TP.cDER "
//	   	+"\tsh "+display_edge_shortest+"\tln "+display_edge_longest
//	   	+"\ttsh "+tree.shortestEdgeLength()
//	   	+"\ttln "+tree.longestEdgeLength()+"\t// "+Arrays.toString(edges));
   }
   
   
   /**
    * Given an array of sorted edge lengths x[0]&le;x[1]&le;...&le;x[n-1], 
    * finds an L and the corresponding range x[i]..x[j] 
    * s.t. x[i] &ge; x[j] * epsilon.
    * 
    * @param edges edge lenth values in increasing order
    * @param epsilon (0&lt;epsilon&lt;1 is assumed) the small fractional edge length below which there is no clear visual distinction in the display
    * @return min and max settings.
    */
   protected static double[] calculateDisplayLengthRange(double[] edges, double epsilon)
   {
       int n = edges.length;
       double[] lrange = new double[2];
       if (n!=0)
       {
	       lrange[0] = edges[0];
	       lrange[1] = edges[edges.length-1];
       }
       if (lrange[0]< lrange[1]/epsilon)
       {
	       double f = epsilon; // epsilon*epsilon / (1.-epsilon);
	       double f2 = Math.sqrt(f);
	       
	       
	       
	       int m = n/2;
	       // try centering around the median
	       double mid = edges[m]; // median
	       
	       int lo = m;
	       while (lo>0 && edges[lo-1]>= mid*f2)
	    	   lo--;       
	       int hi = m;
	       while (hi+1<edges.length && edges[hi+1]<= mid/f2)
	    	   hi++;
	       
	       int imin = lo;
	       int imax = hi;
	       while (imax+1<edges.length && edges[imax+1]<=edges[imin]/f)
	    	   imax++;
	       
	       int[] irange = new int[2];
	       irange[0] = imin;
	       irange[1] = imax;
	       int max_range = imax-imin; // don't bother with the +1
	       
	       
//	       System.out.println("#**TP.cDLR lo "+lo+"\tmid "+mid+"\thi "+hi+"\timin "+imin+"\timax "+imax+"\telen "+edges.length
//	    		   	+"\t"+Arrays.toString(edges));
	       if (max_range<edges.length-1)
	       {
		       imax = edges.length-1;
		       imin = imax;
		       
		       do // step one-by-one with imax, extend imin  
		       {
//			       System.out.println("#**TP.cDLR iter imin "+imin+"\timax "+imax);
			       while (imin>0 && edges[imin-1]>=edges[imax]*f)
			    	   --imin;
//			       System.out.println("#**TP.cDLR fini imin "+imin+"\timax "+imax);
			       assert (imin<=imax);
		    	   if ((imax-imin)>max_range)
		    	   {
		    		   irange[0]=imin;
		    		   irange[1]=imax;
		    		   max_range = imax-imin;
		    	   }
		    	   if (imin==0) break; // avoids imax==-1
		    	   --imax;
		       } while (max_range<edges.length-1);
	       }
	       lrange[0] = edges[irange[0]];
	       lrange[1] = edges[irange[1]];
       }

       
       return lrange;
   }
   
   
   /**
    * Node placement in phenogram style.
    * 
    * @param node
    * @return
    */
   private IndexedPoint placeNodePhenogram(int node)
   {
       int num_children=tree.getNumChildren(node);
       // X coordinate at the median
       int cm1 = (num_children-1)/2;
       int cm2 = num_children/2;
       double nx = 0.5*(node_locations[tree.getChild(node, cm1)].getX()
                       + node_locations[tree.getChild(node, cm2)].getX());
       // Y coordinate is one larger than the children's max
       double ny = 0.0;
       for (int child_idx=0; child_idx<num_children; child_idx++)
       {
           int C = tree.getChild(node,child_idx);
           double cy = node_locations[C].getY();
           ny = Math.max(cy+1.0, ny);
       }
       return new IndexedPoint(node,nx,ny);
   }

   /** 
    * Placement of parent node phylogram style --- 
    * using edge lengths ({@link #getDisplayEdgeLength(int)}) to children. 
    * 
    * @param node
    * @return 
    */
   private IndexedPoint placeNodePhylogram(int node)
   {
       int num_children = tree.getNumChildren(node);
       double ny = 0.0;
       for (int child_idx=0; child_idx<num_children; child_idx++)
       {
           int C = tree.getChild(node,child_idx);
           double cy = node_locations[C].getY();
           
           ny = Math.max(cy+getDisplayEdgeLength(C), ny);
       }
       
//       {
//    	   double max_cy = 0.0;
//    	   for (int ci=0; ci<num_children; ci++)
//    	   {
//    		   int child = tree.getChild(node, ci);
//    		   double cy = node_locations[child].getY();
//    		   if (ny==cy)
//    		   {
//    			   System.out.println("#**TP.pNP "+node+"->"+child
//    					   +"\tlen "+tree.getLength(child)+"/"+getDisplayEdgeLength(child)
//    					   +"\tny "+ny+"\tcy "+cy);
//    		   }
//    		   max_cy = Double.max(max_cy, cy);
//    		   assert (ny>cy);
//    	   }
//    	   assert (ny>max_cy);
//       }
       

       double nx;
       if (num_children%2==1)
       {
           // odd number of children: put x coordinate at the middle child, but shift it a little bit 
           // or else it looks too confusing with edge labels and fancy colorings
           nx = node_locations[tree.getChild(node,num_children/2)].getX();
           if (num_children!=1)
           { // i.e., at least 3
               int neighbor = num_children/2+(((int)ny) % 2==0 ? 1 : -1); // somewhat randomly left or right neighbor 
               double nx_shifted = (0.6*nx+0.4*node_locations[tree.getChild(node,neighbor)].getX());
               //System.out.println("#*TP.p "+i+"/"+node[i].newickName()+"\t"+nx_shifted+" ["+nx+"]\t"+ny);
               nx = nx_shifted;
           }
       } else
       {
           int cm1 = (num_children-1)/2;
           int cm2 = num_children/2;
           nx = 0.5*(node_locations[tree.getChild(node,cm1)].getX()
                       + node_locations[tree.getChild(node,cm2)].getX());
       }
       return new IndexedPoint(node,nx,ny);
   }
   
   /**
    * Computes {@link #node_locations} by table logic; root is at maximum y, x is by post order.
    */
   private void placeNodeTable()
   {
	   int[] node_heights = TreeTraversal.getHeights(tree);
	   int[] node_order = TreeTraversal.postOrder(tree);
	   
	   for (int i=0; i<node_heights.length; i++)
	   {
		   int node = node_order[i];
		   double y = node_heights[node];
		   double x = i;
		   
		   node_locations[node] = new IndexedPoint(node, x, y);
//		   System.out.println("#**TP.pNT "+i+"\tnode "+node+"\t"+node_locations[node]+"\t"+tree.toString(node));
	   }
   }
   
   
   private void initializeGraphicalElements()
   {
       setSelectionAreaColor(AREA_SELECTION_COLOR);
       int tree_point_size = getTreePointSize();
       setCloseRadius(tree_point_size);
       PointIcon leaf_icon=new BoxIcon(tree_point_size,true);
       leaf_icon.setDrawColor(TREE_UNSELECTED_LEAF_COLOR);
       leaf_icon.setFillColor(getBackground());
       PointIcon selected_leaf_icon=new BoxIcon(tree_point_size,true);
       selected_leaf_icon.setDrawColor(TREE_SELECTED_LEAF_COLOR);
       selected_leaf_icon.setFillColor(TREE_SELECTED_LEAF_COLOR);
       PointIcon node_icon=new DiamondIcon(tree_point_size,true);
       node_icon.setDrawColor(TREE_UNSELECTED_NODE_COLOR);
       node_icon.setFillColor(getBackground());
       PointIcon selected_node_icon=new DiamondIcon(tree_point_size,true);
       selected_node_icon.setDrawColor(TREE_SELECTED_NODE_COLOR);
       selected_node_icon.setFillColor(TREE_SELECTED_NODE_COLOR);
       
       int num_nodes = tree.getNumNodes();
       for (int node=0; node<num_nodes; node++)
       {
    	   DisplayedNode node_display;
    	   if (tree.isLeaf(node))
    	   {
    		   node_display 
    		   = new DisplayedNode(node); //, tree.getName(node), null);
    		   node_display.setIcon(true, selected_leaf_icon);
    		   node_display.setIcon(false, leaf_icon);
           } else 
           {
               String label_text = IndexedTree.NODE_IDENT_PREFIX+Integer.toString(node);
        	   node_display 
        	   = new DisplayedNode(node); //, label_text, tree.getName(node));
        	   node_display.setIcon(false, node_icon);
        	   node_display.setIcon(true, selected_node_icon);
               node_display.setLabelBackgroundFill(TREE_NODE_LABEL_BACKGROUND);
           }
		   displayed_node[node] = node_display;
    	   this.addPoint(node_display);
       }
       this.setBackground(Color.WHITE);
   }
   
   
   protected void calculateBoundingBoxes(Graphics g)
   {
	   int num_nodes = tree.getNumNodes();
	   for (int node =0; node<num_nodes; node++)
	   {
		   displayed_node[node].setNodeBoundingBox(g);
		   if (!tree.isRoot(node))
			   displayed_node[node].setEdgeBoundingBox(g);
	   }
   }
   
   protected final void setNormalFontSize(int d)
   {
       this.normal_label_font_size = d;
       setMagnification(magnification);
   }
   
   /**
    * Font size used at this magnification level. 
    * @return 
    */
   protected final int getLabelFontSize()
   {
       return this.label_font_size;
   }
   
   /**
    * Font size used with 100% magnification.
    * 
    * @return 
    */
   protected final int getNormalFontSize()
   {
       return this.normal_label_font_size;
   }

   protected final double getMagnification()
   {
       return this.magnification;
   }
   
   /**
    * Sets the new magnification value. 
    * Setting this value updates the font size {@link #label_font_size} and calls
    * {@link #repaint() }. Also fires a property change (named {@link #MAGNIFICATION_PROPERTY}). 
    * 
    * @param mag desired scale (positive, normal size for 1.0)
    */
   public void setMagnification(double mag)
   {
       double old_mag = this.magnification;
       this.magnification = mag;
       label_font_size = (int)(normal_label_font_size * Math.sqrt(magnification)+0.5);
       //label_font_size = Math.max(label_font_size, Mien.FONT_SIZE_MIN);
       label_font_size = Math.min(label_font_size, FONT_SIZE_MAX);
       this.firePropertyChange(MAGNIFICATION_PROPERTY, old_mag, mag);
       setValidBoundingBoxes(false);
       repaint();
   }
   
   protected final void setValidBoundingBoxes(boolean is_valid)
   {
       this.valid_bounding_boxes = is_valid;
   }
   
   
   /**
    *  @return affine transformation to find display coordinates
    */
  protected AffineTransform getDisplayTransform()
  {
      int w=getWidth();
      int h=getHeight();
      
      AffineTransform plot_transform=new AffineTransform();        
      plot_transform.setToIdentity();
         
      if (LayoutStyle.NODE_TABLE==layout_style)
      {
    	  int root = tree.getRoot();
    	  double node_stretch 
    	  		= (h
    			  -padding.left
    			  +displayed_node[0].node_bounding_box.getY()
    			  -displayed_node[root].node_bounding_box.getMaxY()
    			  -padding.right)/tree_space_width;
    	  double label_width=0.0;
    	  for (int node=0; node<tree.getNumNodes(); node++)
    	  {
    		  Rectangle2D nodeBB = displayed_node[node].node_bounding_box;
    		  double node_ends = nodeBB.getMaxX();
    		  if (node_ends>label_width)
    			  label_width = node_ends;
    	  }
    	  // not correct here ?
    	  label_width+=bounding_box_separation;
    	  double level_stretch
    	  		= (w
    	  			-padding.top-padding.bottom
    	  			-bounding_box_separation
    	  			-label_width)/tree_space_height;
    	  
    	  
//    	  System.out.println("#**TP.gDT ns "+node_stretch+" ls "+level_stretch+"\tw "+w+"\th "+h);
    	  
    	  
	      double tx = padding.top+bounding_box_separation;
	      double ty = padding.left-displayed_node[0].node_bounding_box.getMinY();
	    		  // h; //-padding.left-displayed_node[0].node_bounding_box.getMinY();
	      plot_transform.translate(tx,ty); // wrt original XY coordinates 
	      plot_transform.scale(level_stretch,-node_stretch);
	      plot_transform.quadrantRotate(-1); // -90==270 degrees
	      
      } else
      {
	      int num_leaves = tree.getNumLeaves();
	      int last_leaf_idx = num_leaves-1;
	      double x_scale = (w-
	                   padding.left 
	                  + displayed_node[0].node_bounding_box.getMinX()
	                  - displayed_node[last_leaf_idx].node_bounding_box.getMaxX()
	                  - padding.right)/tree_space_width;
	
	      int root_idx = tree.getRoot();
	      
	      double stem = Math.max(root_stem_length, 
	    		  displayed_node[root_idx].node_bounding_box.getMaxY());
	          
	      double y_scale = (h-padding.top
	                  - stem
	                  + displayed_node[0].node_bounding_box.getMinY()
	                  -padding.bottom
	                  )/tree_space_height;
	
	      plot_transform.translate(padding.left-displayed_node[0].node_bounding_box.getMinX(), h-padding.bottom-stem);
	      plot_transform.scale(x_scale,-y_scale);
      }
      return plot_transform;
  }
   
   protected DisplayedNode getNode(int node)
   {
	   return displayed_node[node];
   }
   
   @Override
   public IndexedPoint getClosestPoint(double x, double y)
   {
	   IndexedPoint closest = super.getClosestPoint(x, y);
	   if (closest == null && layout_style == LayoutStyle.NODE_TABLE)
	   {
		   	double mindist = 2.0*getTreePointSize(); //   Double.POSITIVE_INFINITY;
		   	for (DisplayedNode N: displayed_node)
		   	{
		   		if (N.x<=x) // to the right of the displayed node 
		   		{
			   		double dist = Math.abs(N.y-y);
			   		if ( dist<mindist)
			   		{
			   			mindist = dist;
			   			closest = N;
			   		}
		   		}
		   	}
	   }
	   return closest;
   }
   
   private final static double DEFAULT_LEAF_ROTATION=-0.4*Math.PI;
   private double getDefaultTheta(int node)
   {
	   double theta;
	   if (LayoutStyle.NODE_TABLE==layout_style)
	   {
		   theta= 0.0;
	   } else 
	   {
		   if (tree.isLeaf(node))
		   {
			   theta = DEFAULT_LEAF_ROTATION;
		   } else 
		   {
    		   theta = 0.0;
		   }
	   }	   
	   return theta;
   }
   /** 
    * Class for displaying node labels.
    * 
    * @author csuros
    *
    */
   protected class DisplayedNode extends IndexedPoint
   {
       private int dx;
       private int dy;
//       private String label_text;
//       private String full_text;
       private double theta;
       private float anchor_point;
       private PointIcon selected_icon=null;
       private PointIcon unselected_icon=null;
       private Rectangle2D node_bounding_box;
       private Rectangle2D edge_bounding_box;
       
//       protected DisplayedNode(int node, String label_text, String full_text, int dx, int dy, double theta, float anchor_point)
//       {
//    	   super(node);
//           this.label_text = label_text;
//           this.full_text  = full_text;
//           
//           setPosition(dx,dy,theta,anchor_point);
//           
//           this.setNodeBoundingBox((Graphics)null);
//           if (!tree.isRoot(node))
//        	   this.setEdgeBoundingBox(null);
//       }
       
       protected DisplayedNode(int node)
       {
    	   super(node);
//           this.label_text = label_text;
//           this.full_text  = full_text;
           setDefaultPosition();
           this.setNodeBoundingBox((Graphics)null);
           if (!tree.isRoot(node))
        	   this.setEdgeBoundingBox(null);
       }
       
       protected final int getOffsetX(){ return dx;}
       protected final int getOffsetY(){ return dy;}
       protected final String getLabelText()
       { 
    	   int node = getIndex();
    	   String label_text;
    	   if (tree.isLeaf(node))
		   {
    		   label_text = tree.getName(node);
		   } else
    	   {
    		   label_text = IndexedTree.NODE_IDENT_PREFIX+Integer.toString(node); 		   
    	   }
    	   return label_text;
       }
       protected final String getFullText()
       {
    	   int node = getIndex();
    	   if (tree.isLeaf(node)) return null;
    	   else return tree.getName(node);
       }
//       protected final void setLabelText(String txt){this.label_text=txt;} 
       protected final double getRotation(){ return theta;}
       protected final float getAnchorOffset(){ return anchor_point;}
       
       protected final void setPosition( int dx, int dy, double theta, float anchor_point)
       {
    	   this.dx=dx; this.dy=dy;
    	   this.theta=theta;
    	   this.anchor_point = anchor_point;
       }
       
       protected final void setDefaultPosition()
       {
    	   int dx, dy;
    	   float anchor_point;
    	   double theta = getDefaultTheta(getIndex());
    	   if (LayoutStyle.NODE_TABLE==layout_style)
    	   {
    		   dx = getTreePointSize()/2+5;
    		   dy = getLabelFontSize()/2;
    		   anchor_point = 0f;
    	   } else 
    	   {
			   int point_size = getTreePointSize();
    		   if (tree.isLeaf(getIndex()))
    		   {
    			   dx = point_size/2;
    			   dy = -point_size/2;
    			   anchor_point = 0f;
    		   } else 
    		   {
        		   dx = point_size/2+5;
        		   dy = point_size/2;
        		   anchor_point = 0f;
    		   }
    	   }
    	   setPosition(dx, dy, theta, anchor_point);
       }

//       @Override
//       public double getX(){ return super.getX()+dx;}
//       @Override 
//       public double getY(){ return super.getY()+dy;}
       
       private Color label_background_fill;
       private Color label_stroke_color;
//       private Color label_background_stroke;
       
       protected void setLabelBackgroundFill(Color c)
       {
           this.label_background_fill=c;
       }
       
       protected void setColor(Color c)
       {
           this.label_stroke_color = c;
       }
       
       protected Color getColor()
       {
    	   return label_stroke_color;
       }
       
//       public void setLabelBackgroundStroke(Color c)
//       {
//           this.label_background_stroke=c;
//       }
       protected final int getFontHeight(Graphics2D g2)
       {
           return g2.getFontMetrics().getHeight();
       }
       
       /**
        * Draws the node label string.
        * 
        * @param g2
        */
       protected void paintLabel(Graphics2D g2)
       {
           Color old_color = g2.getColor();
    	   int node = getIndex();
    	   Font label_font;
    	   if (tree.isLeaf(node))
    		   label_font = new Font("Serif", Font.PLAIN, label_font_size);
    	   else
    		   label_font = new Font("Serif", Font.ITALIC, label_font_size);
    	   
    	   if (TreePanel.this.isSelected(node))
    	   {
               float selected_size = Math.max(label_font_size, normal_label_font_size);
               label_font = label_font.deriveFont(selected_size).deriveFont(label_font.getStyle() + Font.BOLD);
    	   }
    	   g2.setFont(label_font);

    	   
    	   int hgt = label_font.getSize(); // getFontHeight(g2);
           if (hgt<FONT_SIZE_MIN) // too small to show 
               return;
           
           int x = (int) getX();
           if (LayoutStyle.NODE_TABLE  == layout_style)
           {
        	   x = (int) displayed_node[0].getX(); // everybody the same
           }
           
           int y = (int) getY();
           
           int labelx = x+dx;
           int labely = y+dy;
           
           //int h = label_fm.getHeight()+3;
           String label_text = getLabelText();
           String full_text = getFullText();
           String written_label = label_text;
           if (full_text!=null && !full_text.equals(""))
           {               
               for (int j=full_text.length(); j>=0; j--)
               {
                   String dn = full_text.substring(0,j);
                   String cand_label = label_text+" ["+dn;
                   if (j<full_text.length())
                       cand_label = cand_label+"...";
                   cand_label = cand_label+"]";
                   
                   Rectangle2D covered_by_label = DrawString.getBoundingBoxForRotatedString(g2, cand_label, labelx, labely, theta, anchor_point);
                   int num_nodes_covered = coveredPoints(covered_by_label).size();
                   if (covered_by_label.contains(this))
                   {
                       num_nodes_covered --;
                   }
                   
                   if (j==0 || num_nodes_covered==0 || TreePanel.this.isSelected(node))
                   {
                       written_label = cand_label;
                       break;
                   }
               }
           }
           
           if (label_background_fill != null && layout_style != LayoutStyle.NODE_TABLE)
           {
               Color label_color = g2.getColor();
               g2.setColor(label_background_fill);
               Rectangle2D labelR = DrawString.getBoundingBoxForRotatedString(g2, written_label, labelx, labely, theta, anchor_point);
               
               g2.fillRoundRect((int)labelR.getX(), (int)labelR.getY(), (int)labelR.getWidth(), (int)labelR.getHeight(), hgt/2, hgt/2);
               g2.setColor(label_color);
           }
           if (label_stroke_color != null)
               g2.setColor(label_stroke_color);
           DrawString.drawRotatedString(g2,written_label, labelx, labely, theta, anchor_point); 
           if (LayoutStyle.NODE_TABLE==layout_style && !tree.isLeaf(node))
           {
        	   Stroke old_stroke = g2.getStroke();
        	   if (isSelected(node))
        	   {
        		   g2.setStroke(thick_dashed_stroke);
        		   g2.setColor(Color.BLACK);
        	   } else
        	   {
	        	   g2.setStroke(dashed_stroke);
	        	   g2.setColor(Color.DARK_GRAY);
        	   }
        	   g2.drawLine((int)getX()+dx, labely-dy, (int)(labelx /*-dx+node_bounding_box.getMinX()*/), labely-dy);
        	   g2.setStroke(old_stroke);
           }
           g2.setColor(old_color);
       }
       
       /**
        * Bounding box relative to the node position.
        * 
        * @param g2
        * @return 
        */
       protected Rectangle2D getLabelBoundingBox(Graphics2D g2)
       {
           int hgt = getFontHeight(g2);
           if (hgt<FONT_SIZE_MIN)
           {
               return new Rectangle2D.Double(-hgt/2.0,0,hgt,hgt);
           } else
           {
        	   String label_text = getLabelText();
        	   String full_text = getFullText();
               int x = getOffsetX();
               int y = getOffsetY();
               String shown_text = label_text;
               if (LayoutStyle.NODE_TABLE==layout_style && full_text!=null)
               {
            	   shown_text=label_text+" ["+full_text+"]"; // match displayed String in paintLabel
               }
               return DrawString.getBoundingBoxForRotatedString(g2, shown_text, x, y, theta, anchor_point);
           }
       }
       
       /**
        * Bounding box relative to node position, covering label and node icon. 
        * 
        * @param g if null, then node labels are not included 
        * @return the bounding box
        */
       protected Rectangle2D setNodeBoundingBox(Graphics g)
       {
		   Rectangle2D R; // return value
		   int tree_point_size = getTreePointSize();
    	   if (g==null)
    	   {
    		   R = new Rectangle2D.Double(-tree_point_size/2, -tree_point_size/2, tree_point_size, tree_point_size);
    	   } else
    	   {
               Graphics2D g2 = (Graphics2D)g.create();
               int node = getIndex();
               if (tree.isLeaf(node))
               {
                   Font leaf_font = new Font("Serif", Font.PLAIN, getLabelFontSize());
                   g2.setFont(leaf_font);
                   R = getLabelBoundingBox(g2);
               } else
               {
                   Font inner_font = new Font("Serif", Font.ITALIC, getLabelFontSize());
                   g2.setFont(inner_font);
                   R = getLabelBoundingBox(g2);
                   // enlarge a bit
                   R.setRect(R.getX()-1, R.getY()-2, R.getWidth()+2, R.getHeight()+2);
               }
//               if (layout_style!=LayoutStyle.NODE_TABLE)
//               {
            	   R.add(new Rectangle2D.Double(-tree_point_size,-tree_point_size,tree_point_size*2, tree_point_size*2)); // space for the little square/diamonds icons
//               }
//               System.out.println("#*TP.DN.sNBB "+R+"\t"+this);
    	   }
    	   return node_bounding_box = R;
       }
       
       /**
        * Sets the node bounding box directly.
        * 
        * @param R
        */
       protected final void setNodeBoundingBox(Rectangle2D R)
       {
    	   node_bounding_box = R;
       }
       
       protected final Rectangle2D getNodeBoundingBox()
       {
    	   return node_bounding_box;
       }
       
       protected Rectangle2D setEdgeBoundingBox(Graphics g)
       {
    	   return edge_bounding_box = new Rectangle2D.Double(-(TREE_THICK_EDGE/2), 0, TREE_THICK_EDGE, 0);
       }
       
       protected void setIcon(boolean selected, PointIcon icon)
       {
    	   if (selected) this.selected_icon = icon;
    	   else this.unselected_icon = icon;
       }
       
       protected PointIcon getIcon(boolean for_selected)
       {
    	   return (for_selected?this.selected_icon:this.unselected_icon);
       }
       
       protected void paintNode(Graphics g)
       {
    	   int nx = (int)this.getX();
    	   int ny = (int)this.getY();
           if (TreePanel.this.isSelected(getIndex()) && selected_icon != null)
               selected_icon.paint(g, nx, ny);
           else if (label_font_size>=FONT_SIZE_MIN && unselected_icon != null)
               unselected_icon.paint(g, nx, ny);
       }
       
       /**
        * Draws the edges from this node to its children (or nothing if 
        * this is a leaf). 
        * 
        * @param g2
        */
       protected void paintEdges(Graphics2D g2)
       {
    	   int node = getIndex();
    	   if (tree.isLeaf(node))
    		   return; 
    	   
    	   Color old_color = g2.getColor();
           if (layout_style == LayoutStyle.PITCHFORK)
           {
               double parent_y = getY();
               int cidx=0; 
               int child = tree.getChild(node, cidx);
               double min_diff = Math.abs(parent_y-displayed_node[child].getY());
               int min_child = child;
               cidx++;
               while (cidx<tree.getNumChildren(node))
               {
                   child = tree.getChild(node, cidx);
                   double diff = Math.abs(parent_y-displayed_node[child].getY());
                   if (diff<min_diff)
                   {
                       min_diff = diff;
                       min_child = child;
                   }
                   cidx++;
               }
               double ref_Y = displayed_node[min_child].getY();

               do
               {
                   --cidx;
                   child = tree.getChild(node, cidx);
                   DisplayedNode child_node = getNode(child);
                   g2.setColor(child_node.selected_icon.getFillColor());
                   drawCurvedLine(g2, child_node, this, ref_Y);
                   //System.err.println("#*TP.pE "+node_idx+"->"+child+"\t"+parent_pt+"\tchl "+displayed_node_locations[child]+"\t// "+ref_Y+"\t// "+Mien.getLongNodeName(tree, node_idx)+"; "+Mien.getLongNodeName(tree, child));
               } while (cidx>0);
           } else if (layout_style == LayoutStyle.CLADOGRAM)
           {
               for (int cidx=0; cidx<tree.getNumChildren(node); cidx++)
               {
                   int child = tree.getChild(node, cidx);
                   DisplayedNode child_node = getNode(child);
                   g2.setColor(child_node.selected_icon.getFillColor());
                   drawThickLine(g2,child_node,this);
               }
           } else if (layout_style == LayoutStyle.PHENOGRAM )
           {
               for (int cidx=0; cidx<tree.getNumChildren(node); cidx++)
               {
                   int child = tree.getChild(node, cidx);
                   DisplayedNode child_node = getNode(child);
                   g2.setColor(child_node.selected_icon.getFillColor());
                   drawBentLineVertical(g2,child_node,this);
               }
           } else if (layout_style == LayoutStyle.PHYLOGRAM)
           {
               for (int cidx=0; cidx<tree.getNumChildren(node); cidx++)
               {
                   int child = tree.getChild(node, cidx);
                   double tlen = getTrueEdgeLength(child);
                   double dlen = getDisplayEdgeLength(child);
                   drawBentRangedVertical(g2,tlen,dlen,displayed_node[child],this);
               }
        	   
           }           
           else if (layout_style == LayoutStyle.NODE_TABLE)
           {
               for (int cidx=0; cidx<tree.getNumChildren(node); cidx++)
               {
                   int child = tree.getChild(node, cidx);
                   DisplayedNode child_node = getNode(child);
                   g2.setColor(child_node.selected_icon.getFillColor());
                   drawBentLineHorizontal(g2,child_node,this);
               }
           }
           g2.setColor(old_color);
       }
   }
   /**
    * Line drawing for edges (straight thick lines).
    * 
    * @param g graphics context
    * @param p1 one endpoint
    * @param p2 other endpoint
    */
   private static void drawThickLine(Graphics g, Point2D p1, Point2D p2)
   {
       final int x1 = (int)p1.getX();
       final int y1 = (int)p1.getY();
       final int x2 = (int)p2.getX();
       final int y2 = (int)p2.getY();
       final int x_offset_min = -TREE_THICK_EDGE/2;
       final int x_offset_max = TREE_THICK_EDGE + x_offset_min;
       for (int dx=x_offset_min; dx<x_offset_max; dx++)
           g.drawLine(x1+dx, y1, x2+dx, y2);
//       g.drawLine((int)p1.getX(),(int)p1.getY(),(int)p2.getX(),(int)p2.getY());
//       g.drawLine((int)p1.getX()-1, (int)p1.getY(),(int)p2.getX()-1,(int)p2.getY());
//       g.drawLine((int)p1.getX()+1, (int)p1.getY(),(int)p2.getX()+1,(int)p2.getY());
   }

   private static final double CURVE_CONTROL1 = 1.0;
   private static final double CURVE_CONTROL2 = 0.0;
   
   public static void drawCurvedLine(Graphics2D g2, Point2D p_to, Point2D p_from, double mid_Y)
   {
       double y1 = p_from.getY();
       double x1 = p_from.getX();
       double y2 = p_to.getY();
       double x2 = p_to.getX();
       
       double bx1 = x1;
       double by1 = y1*(1.-CURVE_CONTROL1)+mid_Y*CURVE_CONTROL1;
       double bx2 = x2;
       double by2 = y1*(1.-CURVE_CONTROL2)+y2*CURVE_CONTROL2;
       
       Path2D link = new Path2D.Double();
       link.moveTo(x1, y1);
       link.curveTo(bx1, by1, bx2, by2, x2, y2);
       g2.draw(link);
   }
   
   /**
    * Line drawing for edges (bent lines).
    * 
    * @param g graphics context
    * @param child one endpoint
    * @param parent other endpoint
    */
   public static void drawBentLineVertical(Graphics g, Point2D child, Point2D parent)
   {
	   Color old_color = g.getColor();
	   int px = (int)parent.getX();
	   int py = (int)parent.getY();
	   int cx =  (int) child.getX();
	   int cy = (int) child.getY();
	   
	   g.drawLine(cx, py, cx, cy); 
	   g.setColor(TREE_EDGE_BACKGROUND);
	   g.drawLine(px, py, cx, py); 
	   g.setColor(old_color);
   }
   
   /**
    * 
    * @param g2 Graphics context
    * @param tlength true length
    * @param dlength display length
    * @param child child node
    * @param parent parent node
    */
   private static void drawBentRangedVertical(Graphics2D g2, double tlength, double dlength, DisplayedNode child, DisplayedNode parent)
   {
	   Stroke old_stroke = g2.getStroke();
	   Color old_color = g2.getColor();
	   int px = (int)parent.getX();
	   int py = (int)parent.getY();
	   int cx =  (int) child.getX();
	   int cy = (int) child.getY();
	   
	   g2.setStroke(thick_stroke);
	   g2.setColor(TREE_EDGE_BACKGROUND);
	   g2.drawLine(px, py, cx, py); 
	   //g2.drawLine(cx, py, cx, cy);// background for dashed line
	   
	   if (tlength>dlength)
       {
		   // long edge
    	   g2.setColor(old_color);
		   double solid_fraction = dlength/tlength;
		   int my = (int)(cy + (py-cy)*solid_fraction);
		   g2.drawLine(cx, my, cx, cy);
		   g2.setStroke(dashed_stroke);
		   g2.drawLine(cx, py, cx, my);
       } else
       {
    	   // short edge
//    	   g2.drawLine(cx, py, cx, cy);// background for dashed line
    	   g2.setColor(old_color);
    	   int ey = (int)(py + (cy-py) * tlength / dlength);
    	   g2.drawLine(cx,py, cx,ey);
       }
	   g2.setStroke(old_stroke);
	   g2.setColor(old_color);
   }
      
   
   /**
    * Line drawing for edges (bent lines).
    * 
    * @param g graphics context
    * @param child one endpoint
    * @param parent other endpoint
    */
   public void drawBentLineHorizontal(Graphics g, Point2D child, Point2D parent)
   {
	   Color old_color = g.getColor();
	   int px = (int)parent.getX();
	   int py = (int)parent.getY();
	   int cx =  (int) child.getX();
	   int cy = (int) child.getY();
	   
	   g.drawLine(px, cy, cx, cy);
	   g.setColor(TREE_EDGE_BACKGROUND);
	   g.drawLine(px, py, px, cy);
	   g.setColor(old_color);
   }
   


   /*
    * Scrollable interface
    */
   @Override
   public Dimension getPreferredScrollableViewportSize() 
   {
	   return getPreferredSize();
   }
	
	@Override
	public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction) 
	{
		return 4*label_font_size;
	}
	
	@Override
	public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation, int direction) 
	{
		int increment;
		int bw = 9*visibleRect.width/10;
		int bh = 9*visibleRect.height/10;
		
		if (orientation == SwingConstants.VERTICAL)
		{
			increment = bh;
		} else
		{
			increment = bw;
		}
		
//		if (LayoutStyle.NODE_TABLE == layout_style)
//		{
//			if (orientation == SwingConstants.VERTICAL)
//			{
//				// skip between nodes 
//				DisplayedNode last_before = null;
//				DisplayedNode first_after = null;
//				for (DisplayedNode N: displayed_node)
//				{
//					if (N.y>=visibleRect.getMinY())
//					{
//						if (N.y<=visibleRect.getMaxY())
//						{
//							// within 
//						} else // after
//							if (first_after==null || first_after.y>N.y)
//								first_after = N;
//					} else // before
//						if (last_before==null || last_before.y<N.y)
//							last_before = N;
//				}
//				if (direction<0)
//				{
//					if (last_before != null)
//					{
//						double d = visibleRect.getCenterY()-last_before.y;
//						assert d>0;
//						bh = Integer.min(bh, (int)d);
//					}
//				} else
//				{
//					if (first_after != null)
//					{
//						double d = first_after.y-visibleRect.getCenterY();
//						assert d>0;
//						bh = Integer.min(bh, (int)d);
//					}
//				}
//				increment = bh;
//			} else
//			{
//				increment = bw;
//			}
//		} else // vertical layout
//		{
//			if (orientation == SwingConstants.HORIZONTAL)
//			{
//				// skip between nodes 
//				DisplayedNode last_before = null;
//				DisplayedNode first_after = null;
//				for (DisplayedNode N: displayed_node)
//				{
//					if (N.x>=visibleRect.getMinX())
//					{
//						if (N.y<=visibleRect.getMaxX())
//						{
//							// within 
//						} else // after
//							if (first_after==null || first_after.x>N.x)
//								first_after = N;
//					} else // before
//						if (last_before==null || last_before.x<N.x)
//							last_before = N;
//				}
//				if (direction<0)
//				{
//					if (last_before != null)
//					{
//						double d = visibleRect.getCenterX()-last_before.x;
//						assert d>0;
//						bw = Integer.min(bw, (int)d);
//					}
//				} else
//				{
//					if (first_after != null)
//					{
//						double d = first_after.x-visibleRect.getCenterX();
//						assert d>0;
//						bw = Integer.min(bw, (int) d);
//					}
//				}
//				increment = bw;
//			} else
//			{
//				increment = bh;
//			}
//		}
		return increment;
	}
	
	@Override
	public boolean getScrollableTracksViewportWidth()
	{
//		return false;
//		return true;
		return layout_style == LayoutStyle.NODE_TABLE;
	}	
	@Override
	public boolean getScrollableTracksViewportHeight()
	{
		return false;
//		return layout_style != LayoutStyle.NODE_TABLE;
	}

	@Override
	public void saveData(File f) throws IOException 
	{
        PrintStream PS = new PrintStream(f);
        PS.println(CommandLine.getStandardHeader(getClass()));
        PS.println(CommandLine.getStandardRuntimeInfo());
        PS.println(NewickParser.printTree(tree));
        if (PS.checkError()) // also flushes
        {
            PS.close();
            throw new IOException("Cannot write the table.");
        }
        PS.close();		
        tree_data.setFile(f);
        tree_data.setDirty(false);
	}
	
	@Override
	public String toString()
	{
		return getTreeName();
	}
   
}