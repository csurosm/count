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

import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Point2D;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.ListSelectionModel;

import count.ds.IndexedTree;
import count.ds.Phylogeny;
import count.ds.TreeTraversal;
import count.gui.kit.ConvexHull;
import count.gui.kit.IndexedPoint;
import count.io.DataFile;

/**
 * A JDialog for editing a tree. 
 */ 
public class TreeEditDialog extends JDialog
{
    private static final String EDIT_PRUNE = "Prune";
    private static final String EDIT_GRAFT_ABOVE = "Graft as sibling";
    private static final String EDIT_GRAFT_BELOW = "Graft as a new child";
    private static final String EDIT_REROOT = "Root here";
    private static final String EDIT_FUSE = "Fuse into parent";
    private static final String EDIT_CANCEL = "Cancel selection";
    
    private final AppFrame main_application;
    private final DataFile<? extends IndexedTree> initial_data;
    private final List<Zoom<TreeManipulator>> edited_trees;
    private int currently_edited_tree;
    private CardLayout layout;
    
    /**
     * Initialization with the currently selected tree in the active session.
     * 
     * @param daddy_o
     * @param title
     */
    public TreeEditDialog(AppFrame daddy_o, String title)
    {
        super(daddy_o, title, true); // modal
        this.main_application = daddy_o;
        Session sesh = main_application.getActiveSession();
// BUNDLE
        initial_data = sesh.getModelBrowser().getSelectedTreeEntry().getTreeData();
        this.edited_trees = new ArrayList<>();
        initComponents();
    }
    
    /**
     * Final edited tree after user pressed cancel or ok.
     * 
     * @return null if cancelled, or tree was not edited 
     */
    public DataFile<Phylogeny> getEditedTree()
    {
        if (currently_edited_tree==0)            
            return null;
        else
        {
        	TreeManipulator current_tree = edited_trees.get(currently_edited_tree).getTreePanel();
        	DataFile<Phylogeny> edited_data = current_tree.getTreeData();
        	// remove inner node names 
        	Phylogeny phylo = edited_data.getContent();
        	for (int node=0; node<phylo.getNumNodes(); node++)
        	{
        		Phylogeny.Node N = phylo.getNode(node);
        		if (!N.isLeaf() && !N.isRoot())
        			N.setName(null);
        	}
        	return edited_data;
        }
    }
    
    
    private void initComponents()
    {
        layout = new CardLayout();
        this.getContentPane().setLayout(layout);
        addCopy(createCopy(initial_data));
    }

    private DataFile<Phylogeny> createCopy(DataFile<? extends IndexedTree> template)
    {
        Phylogeny phylo = new Phylogeny(template.getContent());
        String orig_name = (template.getFile()==null || template.getFile().getName()==null)?"":template.getFile().getName();
        String orig_parent = (template.getFile()==null)?null:template.getFile().getParent();
        String copy_name = edited_trees.isEmpty()?orig_name:DataFile.createIdentifier(edited_trees.size()-1); //orig_name+"_"+DataFile.anyIdentifier();
        
        File copy_file = new File(orig_parent, copy_name);
        DataFile<Phylogeny> copy = new DataFile<>(phylo,copy_file);
        
        return copy;
    }
    
    private TreeManipulator addCopy(DataFile<Phylogeny> copy)
    {
        while (edited_trees.size()>currently_edited_tree+1)
        {
            int idx = edited_trees.size()-1;
            Zoom<TreeManipulator> deleted_tree_panel = edited_trees.get(idx);
            layout.removeLayoutComponent(deleted_tree_panel);
            this.getContentPane().remove(deleted_tree_panel);
            edited_trees.remove(idx);
            
            assert (edited_trees.size() == idx);
        }
        TreePanel.LayoutStyle tree_layout =edited_trees.isEmpty()
                        ?TreePanel.LayoutStyle.PHENOGRAM
                        :edited_trees.get(currently_edited_tree).getTreePanel().getTreeLayoutStyle();
        
        this.currently_edited_tree = edited_trees.size();

        TreeManipulator edit_this = new TreeManipulator(copy, tree_layout);
        
// BUNDLE        
        main_application.getActiveSession().getModelBrowser().decorateByMainTree(edit_this);
        
        Zoom<TreeManipulator> zumm = new Zoom<>(edit_this);
        edit_this.setupControlBar(zumm.getControlBar());

        String zumm_name = DataFile.createIdentifier(currently_edited_tree);
        this.getContentPane().add(zumm, zumm_name);
        layout.show(this.getContentPane(), zumm_name);
        edited_trees.add(zumm);

        return edit_this;
    }
    
    
        
    private void undoEditedTrees()
    {
        assert (this.currently_edited_tree>0);
        
        this.currently_edited_tree--;
        Zoom<TreeManipulator> zumm = edited_trees.get(currently_edited_tree);
        String zumm_name = DataFile.createIdentifier(currently_edited_tree);
        this.getContentPane().add(zumm, zumm_name);
        layout.show(this.getContentPane(), zumm_name);
    }
    
    
    private class TreeManipulator 
    	extends TreePanel 
    	implements MouseListener, MouseMotionListener // we shall listen to ourselves
    {
        private final Phylogeny starting_phylo;

        TreeManipulator(DataFile<Phylogeny> starting_data, TreePanel.LayoutStyle layout)
        {
            super(starting_data, layout, ListSelectionModel.SINGLE_SELECTION);
            this.starting_phylo = (Phylogeny) starting_data.getContent();
            this.setAreaSelectionEnabled(false);
            this.setPointSelectionEnabled(false);
            
            initListeners();
            
        }
        private IndexedPoint prune_position;
        private int[] subtree_node_indices;
        private Polygon pruned_tree_hull;

        private Point regraft_position;

        @Override
        public DataFile<Phylogeny> getTreeData()
        {
        	return (DataFile<Phylogeny>) super.getTreeData(); // warning, but this must be instantiated with a Phylogeny
        }
        
        /**
         * Makes a small drawing of the selected subtree. 
         * 
         * @param g 
         */
        private void paintGraftingSubtree(Graphics g)
        {
            Graphics2D g2 = (Graphics2D)g.create();

            
            double FIXED_SIZE = 25.0*TREE_POINT_SIZE;

            double w = Double.max(pruned_tree_hull.getBounds2D().getWidth(),1.5*FIXED_SIZE);
            double h = Double.max(pruned_tree_hull.getBounds2D().getHeight(),1.5*FIXED_SIZE);
            
            double scala = Math.min(FIXED_SIZE/(w+1e-9), FIXED_SIZE/(h+1e-9));
            
            g2.scale(scala, scala);
            g2.translate(regraft_position.getX()/scala-prune_position.getX(), regraft_position.getY()/scala-prune_position.getY());
            g2.setColor(AppFrame.SMOKY_BACKGROUND);
            if (pruned_tree_hull.npoints>2)
            	g2.fillPolygon(pruned_tree_hull);
            g2.setColor(TREE_EDGE_COLOR);
            for (int pidx: subtree_node_indices)
            {
                if (!starting_phylo.isLeaf(pidx))
                	getNode(pidx).paintEdges(g2);
            }
            for (int pidx: subtree_node_indices)
            {
            	getNode(pidx).paintNode(g2);
            }
        }
        
        @Override
        protected void paintComponent(Graphics g)
        {
            super.paintComponent(g);
            if (prune_position != null )
            {
                Color col = g.getColor();
                if (pruned_tree_hull.npoints>2)
                {
                    g.setColor(AREA_SELECTION_COLOR);
	                g.fillPolygon(pruned_tree_hull);
	                g.setColor(Color.RED);
	                g.drawPolygon(pruned_tree_hull);
                } else
                	g.setColor(Color.RED);
                
                g.drawLine((int) prune_position.getX(), (int)prune_position.getY(), (int)regraft_position.getX(), (int)regraft_position.getY());
                paintGraftingSubtree(g);
                g.setColor(col);
            }
        }
        
        
        private void initListeners()
        {
            this.addMouseListener(this);
            this.addMouseMotionListener(this);
        }
      
        private void setupControlBar(Box bar)
        {
    		bar.removeAll();
	        bar.add(createLayoutChooser());
	        bar.add(Box.createHorizontalGlue());
	        
	        final JButton undo_button = new JButton("Undo");
	        final JButton ok_button = new JButton("Done editing");
	        final JButton cancel_button = new JButton("Cancel editing");
	
	        bar.add(cancel_button);            
	        cancel_button.addActionListener(e-> 
	              {
	                  TreeEditDialog.this.currently_edited_tree = 0;
	                  TreeEditDialog.this.dispose();
	              });
	
	        if (TreeEditDialog.this.currently_edited_tree>0)
	        {
	        	bar.add(undo_button);
	            undo_button.addActionListener(e->
	                  {
	                      undoEditedTrees();
	                  });
	        }
	        bar.add(ok_button);
	        ok_button.addActionListener(e->
	              {
	                  TreeEditDialog.this.dispose();
	              });
	
	        bar.add(Box.createHorizontalGlue());
	        bar.add(createMagnificationSpinner());
	    }

        private void clearPruningNode()
        {
            this.prune_position = null;
            this.pruned_tree_hull = null;
            this.subtree_node_indices = null;
//            getSelectionModel().clearSelection();
        }
        
        private void setPruningNode(int node_index)
        {
            this.prune_position = getNode(node_index);
            this.regraft_position = new Point((int)prune_position.getX(), (int)prune_position.getY());
            
            final List<Point2D> point_list = new ArrayList<>();
            final List<DisplayedNode> subtree_nodes = new ArrayList<>();
            
            TreeTraversal TT = new TreeTraversal(starting_phylo);
            TT.traverse(node_index, node->{
                DisplayedNode node_loc = getNode(node);
                point_list.add(node_loc);
                subtree_nodes.add(node_loc);
            	if (node != node_index)
                { // if not the subtree root, then add a point for the convex hull
                    int parent = starting_phylo.getParent(node);
                    DisplayedNode parent_loc = getNode(parent);
                    Point2D edge_bend = new Point2D.Double(node_loc.getX(), parent_loc.getY());
                    point_list.add(edge_bend);
                }
            }, null);
            
            
            subtree_node_indices = new int[subtree_nodes.size()];
            for (int i=0;i<subtree_nodes.size(); i++)
            {
            	DisplayedNode P = subtree_nodes.get(i);
                subtree_node_indices[i] = P.getIndex();
            }
            if (subtree_nodes.size()==1)
            {
            	assert (starting_phylo.isLeaf(node_index));
            	int pt = getTreePointSize();
            	point_list.add(new Point2D.Double(prune_position.x-pt, prune_position.y-pt));
            	point_list.add(new Point2D.Double(prune_position.x+pt, prune_position.y-pt));
            	point_list.add(new Point2D.Double(prune_position.x+pt, prune_position.y-pt));
            	point_list.add(new Point2D.Double(prune_position.x+pt, prune_position.y+pt));
            }
            
            
            ConvexHull GS = new ConvexHull(point_list.toArray(new Point2D[]{}));
            java.awt.geom.Point2D[] subtree_convex_hull = GS.getHull();
            this.pruned_tree_hull = new Polygon();
            for (java.awt.geom.Point2D P: subtree_convex_hull)
            {
                this.pruned_tree_hull.addPoint((int)P.getX(), (int)P.getY());
            }
        }

        @Override
        public void mouseMoved(MouseEvent evt)
        {
            if (this.prune_position != null)
            {
                this.regraft_position = evt.getPoint();
                repaint();
            }            
        }

        @Override
        public void mouseClicked(MouseEvent e) 
        {
            int x=e.getX();
            int y=e.getY();

            IndexedPoint p=getClosestPoint(x,y);
            if (p==null)
            {
                if (prune_position != null)
                {
                    clearPruningNode();
                    //repaint();
                } 
                this.getSelectionModel().clearSelection(); // triggers repaint
            } else
            {
                int pidx = p.getIndex();
                this.getSelectionModel().setSelectionInterval(pidx, pidx);
                
                NodeActions node_actions = new NodeActions(pidx);
                node_actions.popup();
            }
        }

        @Override
        public void mousePressed(MouseEvent e) 
        {
            // ignore
        }

        @Override
        public void mouseReleased(MouseEvent e) 
        {
            // ignore
        }

        @Override
        public void mouseEntered(MouseEvent e) 
        {
            // ignore
        }

        @Override
        public void mouseExited(MouseEvent e) 
        {
            // ignore
        }

        @Override
        public void mouseDragged(MouseEvent e) 
        {
            // ignore
        }

        private class NodeActions implements ActionListener
        {
            private final int node;
            
            NodeActions(int node)
            {
                this.node = node;
            }
            
            @Override 
            public void actionPerformed(ActionEvent e)
            {
                JMenuItem src = (JMenuItem) e.getSource();
                String menu_text = src.getText();
//                System.out.println("#*TED.TM.NA.aP "+node+"\t"+menu_text+"\te "+e);
                
                if (EDIT_PRUNE.equals(menu_text))
                {
                    setPruningNode(node);
                } else if (EDIT_GRAFT_ABOVE.equals(menu_text))
                {

                    DataFile<Phylogeny> modified_tree = createCopy(TreeManipulator.this.getTreeData());
                    Phylogeny copied_phylo = modified_tree.getContent();
                    copied_phylo.pruneAndRegraft(copied_phylo.getNode(prune_position.getIndex()), copied_phylo.getNode(node), false);
                    addCopy(modified_tree);
                } else if (EDIT_GRAFT_BELOW.equals(menu_text))
                {
                    DataFile<Phylogeny> modified_tree = createCopy(TreeManipulator.this.getTreeData());
                    Phylogeny copied_phylo = modified_tree.getContent();
                    copied_phylo.pruneAndRegraft(copied_phylo.getNode(prune_position.getIndex()), copied_phylo.getNode(node), true);
                    addCopy(modified_tree);
                } else if (EDIT_REROOT.equals(menu_text))
                {
                    DataFile<Phylogeny> modified_tree = createCopy(TreeManipulator.this.getTreeData());
                    Phylogeny copied_phylo = modified_tree.getContent();
                    copied_phylo.reroot(copied_phylo.getNode(node));
                    addCopy(modified_tree);
                } else if (EDIT_FUSE.equals(menu_text))
                {
                    DataFile<Phylogeny> modified_tree = createCopy(TreeManipulator.this.getTreeData());
                    Phylogeny copied_phylo = modified_tree.getContent();
                    copied_phylo.fuseIntoParent(copied_phylo.getNode(node));
//                    
//                    System.out.println("#*TED.TM.NA.aP "+TreeManipulator.this.getData().getContent().hasLength()+" -> "+copied_phylo.hasLength());
                    addCopy(modified_tree);
                } else if (EDIT_CANCEL.equals(menu_text))
                {
                    // do nothing
                    clearPruningNode();
                    getSelectionModel().clearSelection();
                }
            }
            
            void popup()
            {
                JPopupMenu popup_options = new JPopupMenu();
                popup_options.setBorder(BorderFactory.createTitledBorder(starting_phylo.getName(node)));
                
                if (prune_position == null)
                {
                    if (!starting_phylo.isRoot(node))
                    {
                        JMenuItem prune = new JMenuItem(EDIT_PRUNE);
                        prune.addActionListener(this);
                        popup_options.add(prune);
                        
                        if (!starting_phylo.isLeaf(node))
                        {
                            JMenuItem reroot = new JMenuItem(EDIT_REROOT);
                            reroot.addActionListener(this);
                            popup_options.add(reroot);

                            JMenuItem fuse = new JMenuItem(EDIT_FUSE);
                            fuse.addActionListener(this);
                            popup_options.add(fuse);
                        }
                    }                
                } else // regrafting options
                {
                    JMenuItem graft_above = new JMenuItem(EDIT_GRAFT_ABOVE);
                    graft_above.addActionListener(this);
                    popup_options.add(graft_above);
                    if (!starting_phylo.isLeaf(node))
                    {
                        JMenuItem graft_below = new JMenuItem(EDIT_GRAFT_BELOW);
                        graft_below.addActionListener(this);
                        popup_options.add(graft_below);
                    }

                }
                JMenuItem cancel = new JMenuItem(EDIT_CANCEL);
                cancel.addActionListener(this);
                popup_options.add(cancel);
                
                TreePanel.DisplayedNode point = getNode(node);
                
                popup_options.show(TreeManipulator.this, (int)point.getX(), (int)point.getY());
            }
        }        
    
    } // TreeManipulator
}
