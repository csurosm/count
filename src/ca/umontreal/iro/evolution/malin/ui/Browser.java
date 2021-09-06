/*
 * Browser.java
 *
 * Created on December 30, 2007, 1:21 AM
 *
 */

package ca.umontreal.iro.evolution.malin.ui;


import java.io.IOException;
import java.io.File;

import java.util.Enumeration;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.Rectangle;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.BorderFactory;
import javax.swing.JComponent;
import javax.swing.JEditorPane;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTree;

import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.event.TreeModelListener;
import javax.swing.event.TreeModelEvent;

import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreeSelectionModel;
import javax.swing.tree.TreePath;
import javax.swing.tree.DefaultTreeCellRenderer;

import ca.umontreal.iro.evolution.malin.Exportable;
import ca.umontreal.iro.evolution.malin.Removable;
import ca.umontreal.iro.evolution.malin.Saveable;

/**
 *
 * Class for browsing among primary items (data sets, rate models) and
 * derived views (posterior predictions etc).
 *
 * @author csuros
 */
public class Browser extends JPanel 
{
    public Browser() 
    {
        this(PrimaryItem.class);
    }
    
    
    /**
     * Creates a new instance of Browser
     */
    public Browser(Class primary_item_class) 
    {
        this.primary_item_class = primary_item_class;
        initDataStructures();
        initComponents();
    }

    private static final String REMOVE_NODE_TEXT = "Remove";
    private static final String EXPORT_NODE_TEXT = "Export ...";
    private static final String SAVE_NODE_TEXT   = "Save";
    private static final String SAVE_AS_NODE_TEXT = "Save as ...";
    
    protected Class primary_item_class;
    
    /**
     * A hook for anything that needs to be done at initialization,
     * but before the components are set up by initComponents().
     * Here, it does nothing.
     */
    protected void initDataStructures()
    {
        
    }

    /**
     * Main element in the panel. 
     * Left-hand side containts the browser, the right-hand side 
     * shows the information associated with the selected 
     * item. 
     */
    protected JSplitPane content;

    /**
     * The browser is implemented as a JTree.
     */
    protected JTree item_browser;
    
    /**
     * Root node for the browser.
     */
    protected DefaultMutableTreeNode item_browser_root;
    
    /**
     * Data model for the browser.
     */
    protected DefaultTreeModel item_browser_model;    
    
    /**
     * The scroll for the browser.
     */
    protected JScrollPane item_browser_scroll;
        
    /**
     * The scroll for the item info
     */
    protected JScrollPane item_info_scroll;
    
    /**
     * The text field where messages are written when there is no selected item.
     */
    protected JEditorPane no_selection_text;
    
    protected void initComponents()
    {
        content = new JSplitPane();
        content.setDividerLocation(200+content.getInsets().top);
        content.setOrientation(JSplitPane.HORIZONTAL_SPLIT);
        content.setDebugGraphicsOptions(javax.swing.DebugGraphics.NONE_OPTION);
        content.setBorder(null);
        content.setResizeWeight(0.0);
        content.setOneTouchExpandable(true);
        
        item_browser_root = new DefaultMutableTreeNode("browser root");
        item_browser_model = new DefaultTreeModel(item_browser_root);
        item_browser = new JTree(item_browser_model);
        
        ItemInfoRenderer renderer = new ItemInfoRenderer();
        item_browser.setCellRenderer(renderer);
        
        item_browser.setRootVisible(false);
        item_browser.setShowsRootHandles(true);
        item_browser.setEditable(false);
        item_browser.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        item_browser.addTreeSelectionListener(new TreeSelectionListener() 
        {
            public void valueChanged(TreeSelectionEvent treeSelectionEvent) 
            {
                setItemInfo();
            }
        });
        item_browser.addMouseListener(new MouseAdapter()
            {
                public void mousePressed(MouseEvent e)
                {
                        popupViewMenu(e);
                }
                
                public void mouseReleased(MouseEvent e)
                {
                        popupViewMenu(e);
                }
                
            });
        
        item_browser.setFont(new java.awt.Font("Serif", Font.BOLD,14));
        //rate_chooser.setSelectionBackground(new Color(255,148,86)); // Grapefruit
        item_browser.setVisible(true);
        
        item_browser_scroll = new JScrollPane();
        item_browser_scroll.setViewportView(item_browser);
        item_browser_scroll.setMinimumSize(new Dimension(200,400));
        content.setLeftComponent(item_browser_scroll);        
        
        no_selection_text = new JEditorPane();
        no_selection_text.setContentType("text/html");
        no_selection_text.setBackground(getBackground().brighter());
        no_selection_text.setEditable(false);
        content.setRightComponent(no_selection_text);

        setItemInfo();

        setLayout(new BorderLayout());
        add(content, BorderLayout.CENTER);
        
    }
    
    /**
     * Finds the closest ancestor node that contains a primary item.
     * @return null, if none of the ancestors are of primary class.
     */
    public DefaultMutableTreeNode getSelectedPrimaryNode()
    {
        TreePath selected = item_browser.getSelectionPath();
        if (selected == null)
            return null;
        
        for (int j=selected.getPathCount()-1; j>=0; j--)
        {
            DefaultMutableTreeNode selected_node = (DefaultMutableTreeNode)selected.getPathComponent(j);
            Object selected_node_object = selected_node.getUserObject();
            if (primary_item_class.isInstance(selected_node_object))
                    return selected_node;
        }
        
        return null;
    }

    
    /**
     * Returns the selected item
     * @return null if no selection
     */
    public Object getSelectedItem()
    {
        TreePath selected = item_browser.getSelectionPath();
        if (selected == null)
            return null;
        else
            return ((DefaultMutableTreeNode)selected.getLastPathComponent()).getUserObject();
    }
    
    public Object getSelectedPrimaryItem()
    {
        DefaultMutableTreeNode node = getSelectedPrimaryNode();
        if (node == null)
            return null;
        else 
            return node.getUserObject();
    }

    public Object getSelectedPrimaryRootItem()
    {
        TreePath selected = item_browser.getSelectionPath();
        if (selected == null)
            return null;

        for (int j=0; j<selected.getPathCount(); ++j)
        {
            DefaultMutableTreeNode selected_node = (DefaultMutableTreeNode)selected.getPathComponent(j);
            Object selected_node_object = selected_node.getUserObject();
            if (primary_item_class.isInstance(selected_node_object))
                    return selected_node.getUserObject();
        }

        return null;
    }
        
    /**
     * Adds a listener for changes in the primary item selection and browser structure.
     */
    public void addBrowserSelectionListener(BrowserSelectionListener listener)
    {
        PrimaryItemSelectionListener pisl = new PrimaryItemSelectionListener(listener);
        item_browser.addTreeSelectionListener(pisl);
        item_browser_model.addTreeModelListener(pisl);
    }

    /**
     * The message displayed when there are no available items.
     */
    protected String NO_AVAILABLE_ITEMS = "<em>No available items.</em>";
    
    /**
     * The message displayed when there is no selection.
     */
    protected String SELECT_AN_ITEM = "<em>Select an item.</em>";
    
    /**
     * Called whenever item selection changes
     */
    protected void setItemInfo()
    {
        TreePath selected = item_browser.getSelectionPath();
        int div = content.getDividerLocation(); // so that the divider doesn't change 
        if (selected == null)
        {
            if (item_browser_model.getChildCount(item_browser_root)==0)
                no_selection_text.setText(NO_AVAILABLE_ITEMS);
            else
                no_selection_text.setText(SELECT_AN_ITEM);
            content.setRightComponent(no_selection_text);
            content.setDividerLocation(div);
        } else 
        {
            DefaultMutableTreeNode selected_node = (DefaultMutableTreeNode)selected.getLastPathComponent();
            content.setRightComponent((javax.swing.JComponent) selected_node.getUserObject());
            content.setDividerLocation(div);
        }        
    }

    /**
     * Adds a new item at the top level.
     *
     * @throws IllegalArgumentException if the item does not belong to the primary item class.
     */
    public void addTopItem(JComponent item)
    {
        if (!primary_item_class.isInstance(item))
            throw new IllegalArgumentException("Top level items must be of class "+primary_item_class);
        addItem(item, item_browser_root);
    }

    /** 
     * Adds a new item as the child of the selected primary item.
     */
    public void addItem(JComponent item)
    {
        DefaultMutableTreeNode parent = getSelectedPrimaryNode();
        if (parent != null)
        {
            addItem(item, parent);
        }
    }
    
    /**
     * Adds a new item as the child of a given node.
     */
    protected void addItem(JComponent item, DefaultMutableTreeNode parent)
    {
        DefaultMutableTreeNode view_node = new DefaultMutableTreeNode(item);
        item_browser_model.insertNodeInto(view_node,parent,parent.getChildCount());
        TreePath view_path = new TreePath(view_node.getPath());

        item_browser.scrollPathToVisible(view_path);
        item_browser.setSelectionPath(view_path);
    }

    
    /**
     * Called when a mouse popup trigger is fired on a tree node
     */
    protected void popupViewMenu(MouseEvent e)
    {
        if (e.isPopupTrigger())
        {
            TreePath clicked_node_path = item_browser.getClosestPathForLocation(e.getX(),e.getY());
            if (clicked_node_path == null)
                return;
            Rectangle clicked_node_R = item_browser.getPathBounds(clicked_node_path);
            if (clicked_node_R.contains(e.getX(),e.getY()))
            {
                item_browser.setSelectionPath(clicked_node_path);
                DefaultMutableTreeNode selected_node = (DefaultMutableTreeNode)clicked_node_path.getLastPathComponent();
                JComponent view = (JComponent) selected_node.getUserObject();
                //dataset_chooser_model.removeNodeFromParent(selected_node);

                JPopupMenu popup_options = new JPopupMenu();
                popup_options.setBorder(BorderFactory.createRaisedBevelBorder());
                PopupActionsOnNode popup_actions = new PopupActionsOnNode(selected_node);

                JMenuItem close = new JMenuItem(REMOVE_NODE_TEXT);
                close.addActionListener(popup_actions);

                JLabel tlabel = new JLabel(view.toString());
                tlabel.setForeground(tlabel.getForeground().brighter());
                //tlabel.setBorder(BorderFactory.createLineBorder(Color.BLACK));
                //tlabel.setBackground(popup_options.getForeground());
                popup_options.add(tlabel);
                popup_options.add(new javax.swing.JSeparator());
                popup_options.add(close);

                if (view instanceof Saveable)
                {
                    Saveable Sview = (Saveable)view;
                    if (Sview.isDirty())
                    {
                        tlabel.setText("¥ "+view.toString());
                    } else 
                    {
                        tlabel.setText("  "+view.toString());
                    }
                    JMenuItem save = new JMenuItem(SAVE_NODE_TEXT);
                    save.setEnabled(Sview.hasAssociatedFile() && Sview.isDirty());
                    save.addActionListener(popup_actions);
                    JMenuItem save_as = new JMenuItem(SAVE_AS_NODE_TEXT);
                    save_as.addActionListener(popup_actions);

                    popup_options.addSeparator();
                    popup_options.add(save);
                    popup_options.add(save_as);
                } else if (view instanceof Exportable)
                {
                    Exportable Eview = (Exportable) view;
                    JMenuItem export = new JMenuItem(EXPORT_NODE_TEXT);
                    export.addActionListener(popup_actions);
                    popup_options.addSeparator();
                    popup_options.add(export);
                }

                popup_options.show(this,e.getX(),e.getY());
            }
        }
    }
    
    /**
     * Removes an item from the browser by finding its node through tree traversal.
     */
    public void removeItem(JComponent item)
    {
        java.util.Enumeration traversal = item_browser_root.depthFirstEnumeration();
        while(traversal.hasMoreElements())
        {
            DefaultMutableTreeNode node = (DefaultMutableTreeNode)traversal.nextElement();
            if (item == node.getUserObject())
            {
                removeNode(node);
                break;
            }
        }
    }
    
    private void removeNode(DefaultMutableTreeNode node)
    {
        boolean doRemove = true;
        java.util.Enumeration enum_subtree = node.depthFirstEnumeration();
        int subtree_size = 0;
        while (enum_subtree.hasMoreElements())
        {
            subtree_size ++;
            DefaultMutableTreeNode descendant = (DefaultMutableTreeNode) enum_subtree.nextElement();
        }
        if (subtree_size != 1)
        {
            // verify if the user really wants to close this 
            // even if there are descendant views.
            String msg = null; 
            if (subtree_size == 2)
            {
                msg = "There is a descendant view, which will also be removed along with this one.\nIs that OK?";
            } else
            {
                msg = "There are "+(subtree_size-1)+" descendant views, which will also be removed along with this one.\nIs that OK?";
            }
            int ans = JOptionPane.showConfirmDialog(item_browser,msg,"Remove from browser",JOptionPane.YES_NO_OPTION,JOptionPane.WARNING_MESSAGE);
            doRemove = (ans == JOptionPane.YES_OPTION);
        }
        
        if (doRemove && node.getUserObject() instanceof Removable)
        {
            Removable R = (Removable) node.getUserObject();
            doRemove = R.remove();
        }

        if (doRemove)
        {
            item_browser.clearSelection();
            item_browser_model.removeNodeFromParent(node);
            setItemInfo();
        }
    }
    
    
    /**
     * Selects a given item in the browser. 
     * 
     * @param item an item that was added to the browser before
     * @return true if the item was found, or false if not
     */
    public boolean selectItem(JComponent item)
    {
        Enumeration item_enum = item_browser_root.postorderEnumeration();
        
        while (item_enum.hasMoreElements())
        {
            Object visited_item = item_enum.nextElement();
            if (visited_item instanceof DefaultMutableTreeNode) // it is ...
            {
                DefaultMutableTreeNode visited_node = (DefaultMutableTreeNode)visited_item;
                Object visited_node_content = visited_node.getUserObject();
                if (visited_node_content==item)
                {
                    // found the right guy
                    TreePath view_path = new TreePath(visited_node.getPath());

                    item_browser.scrollPathToVisible(view_path);
                    item_browser.setSelectionPath(view_path);
                    //System.out.println("#*B.sI found "+item+"\tnode "+visited_node);
                    
                    return true;
                }
            }
        }
        //System.out.println("#*B.sI not found "+item);
        return false;
    }

    private class PopupActionsOnNode implements ActionListener
    {
        private PopupActionsOnNode(DefaultMutableTreeNode node)
        {
            this.node=node;
        }
        
        public void actionPerformed(ActionEvent E)
        {
            JMenuItem src = (JMenuItem) E.getSource();
            JComponent view = (JComponent)node.getUserObject();
            
            String menu_text = src.getText();
            if (menu_text.equals(REMOVE_NODE_TEXT))
            {
                removeNode(node);
            } // delete
            else if (EXPORT_NODE_TEXT.equals(menu_text))
            {
                FileDialog dialog =  new FileDialog((java.awt.Frame)getTopLevelAncestor(),"Export \""+view.toString()+"\"",FileDialog.SAVE);
                dialog.setVisible(true);
                String file_name = dialog.getFile();
                String directory = dialog.getDirectory();                
                if (file_name!=null)
                {
                    try 
                    {
                        File export_file = new File(directory, file_name); // also checks if overwriting is OK
                        ((Exportable)view).saveData(export_file);
                    } catch (IOException exc)
                    {
                        Dealer dealer = WorkSpace.getDealer(Browser.this);
                        dealer.exceptionCaught(exc, "I/O error", "Error while exporting data view.");
                    }
                }
            } else if (SAVE_NODE_TEXT.equals(menu_text))
            {
                try {
                    Saveable Sview = (Saveable)view;
                    Sview.saveData();
                } catch (IOException exc)
                {
                    Dealer dealer = WorkSpace.getDealer(Browser.this);
                    dealer.exceptionCaught(exc, "I/O error", "Error while saving data.");
                }
            } else if (SAVE_AS_NODE_TEXT.equals(menu_text))
            {
                FileDialog dialog =  new FileDialog((java.awt.Frame)getTopLevelAncestor(),"Save \""+view.toString()+"\"",FileDialog.SAVE);
                dialog.setVisible(true);
                String file_name = dialog.getFile();
                String directory = dialog.getDirectory();                
                if (file_name!=null)
                {
                    try {
                        Saveable Sview = (Saveable)view;
                        File save_file = new File(directory, file_name); 
                        Sview.saveData(save_file);
                        Browser.this.repaint(); // since name may have changed
                    } catch (IOException exc)
                    {
                        Dealer dealer = WorkSpace.getDealer(Browser.this);
                        dealer.exceptionCaught(exc, "I/O error", "Error while saving data.");
                    }
                }
                
            }
            
        } 
        
        private DefaultMutableTreeNode node;

    }    

    /**
     * Our own rendering for tree nodes: primary items as folders, others as documents.
     */
    protected class ItemInfoRenderer extends DefaultTreeCellRenderer
    {
        
        
        public java.awt.Component getTreeCellRendererComponent(
                        JTree tree,
                        Object value,
                        boolean sel,
                        boolean expanded,
                        boolean leaf,
                        int row,
                        boolean hasFocus) {
            super.getTreeCellRendererComponent(
                            tree, value, sel,
                            expanded, leaf, row,
                            hasFocus);
            
            DefaultMutableTreeNode node = (DefaultMutableTreeNode) value;
            if (primary_item_class.isInstance(node.getUserObject()))
            {
                if (expanded)
                    setIcon(openIcon);
                else
                    setIcon(closedIcon);
            } 
            return this;
        }
    }    
    
    public static class PrimaryItemSelectionEvent extends java.util.EventObject
    {
        public PrimaryItemSelectionEvent(Object src)
        {
            super(src);
        }
    }
    
    protected class PrimaryItemSelectionListener implements TreeSelectionListener, TreeModelListener
    {
        
        BrowserSelectionListener listener;
        
        protected PrimaryItemSelectionListener(BrowserSelectionListener L)
        {
            listener = L;
        }

        private void structureChanged()
        {
            Object src = getSelectedPrimaryNode();
            if (src == null)
                src = Browser.this;
            else
                src = ((DefaultMutableTreeNode)src).getUserObject();
            listener.valueChanged(new PrimaryItemSelectionEvent(src));
        }
        
        public void treeNodesChanged(TreeModelEvent e)
        {
            // nothing
        }
        
        public void treeNodesInserted(TreeModelEvent e)
        {
            structureChanged();
        }
        public void treeNodesRemoved(TreeModelEvent e)
        {
            structureChanged();
        }
        public void treeStructureChanged(TreeModelEvent e)
        {
            // nothing
        }
        
        public void valueChanged(TreeSelectionEvent treeSelectionEvent) 
        {
            DefaultMutableTreeNode daddy = getSelectedPrimaryNode();
            if (daddy != null)
                listener.valueChanged(new PrimaryItemSelectionEvent(daddy.getUserObject()));
        }
    }
    
    /**
     * A common superclass for primary items.
     */
    public static class PrimaryItem extends JPanel
    {
        
    }
    
    
    public Enumeration preorderNodeEnumeration()
    {
        Enumeration E = item_browser_root.preorderEnumeration();
        E.nextElement() ; // skip root
        return E;
    }
    
    public ItemEnumeration preorderEnumeration(){return new ItemEnumeration();}

    /**
     * Class for enumerating items in the browser
     */
    public class ItemEnumeration 
    {
        private Enumeration traversal;
        private ItemEnumeration()
        {
            traversal = item_browser_root.preorderEnumeration();
            traversal.nextElement(); // skip root
        }
        
        /**
         * Whether there are more elements there.
         * 
         * @return true if the enumeration is not done yet
         */
        public boolean hasMoreElements()
        {
            return traversal.hasMoreElements();
        }

        /**
         * Next item of the enumeration 
         * @return next item or null
         */
        public EnumeratedItem nextElement()
        {
            DefaultMutableTreeNode node = (DefaultMutableTreeNode)traversal.nextElement();
            Object o = node.getUserObject();
            
            DefaultMutableTreeNode parent = (DefaultMutableTreeNode)node.getParent();
            Object p = null;
            if (parent != item_browser_root)
            {
                p = parent.getUserObject();
            }
            EnumeratedItem item = new EnumeratedItem(node,p);
            return item;
        }
    }
    
    /**
     * Class used in the enumeration of items
     */
    public static class EnumeratedItem 
    {
        private Object parent;
        private DefaultMutableTreeNode node;
        
        private EnumeratedItem(DefaultMutableTreeNode node, Object parent)
        {
            this.node = node;
            this.parent = parent;
        }
        /**
         * The item in the browser
         * 
         * @return item
         */
        public Object getItem(){return node.getUserObject();}
        
        public DefaultMutableTreeNode getNode(){return node;}
        /**
         * Parent item 
         * @return parent item or null
         */
        public Object getParent(){return parent;}
    }
    
    
    
}
