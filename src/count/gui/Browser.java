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

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;

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
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.TreeModelEvent;
import javax.swing.event.TreeModelListener;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import count.io.DataFile;
import count.io.ExportableData;
import count.io.Removable;
import count.io.SavableData;

/**
 * First implementation of browsing GUI.
 * 
 *
 * @param <TOP>
 * @deprecated
 */
public class Browser<TOP extends JComponent> extends JPanel 
{ 
    /**
     * Creates a new instance of Browser as a standalone graphical element.  
     * @param primary_item_class
     */
    public Browser(Class<TOP> primary_item_class) 
    {
    	this(primary_item_class, defaultSplitPane());

    	setLayout(new BorderLayout());
        add(content, BorderLayout.CENTER);
    }
    
    /**
     * Puts the two components into an externally constructed JSplitPane
     * 
     * @param primary_item_class
     * @param container 
     */
    public Browser(Class<TOP> primary_item_class, JSplitPane container)
    {
    	
        this.primary_item_class = primary_item_class;
        initComponents(container);
    }

    /** 
     * Embeds the browser in a given JSplitPane, without creating 
     * a standalone graphical element: the embedding pane's right component is  
     * set by tree selection in left component. 
     * 
     * @param <P> the type of primary elements in this browser
     * @param primary_item_class 
     * @param embedPane 
     * @return left component is a JScrollPane around a JTree; right component changes as items are added 
     */
    protected static <P extends Browser.PrimaryItem> JSplitPane embedBrowser(Class<P> primary_item_class, JSplitPane embedPane)
    {
    	Browser<P> browser = new Browser<>(primary_item_class, embedPane);
    	return browser.content;
    }
    
    private static final String REMOVE_NODE_TEXT = "Remove";
    private static final String EXPORT_NODE_TEXT = "Export ...";
    private static final String SAVE_NODE_TEXT   = "Save";
    private static final String SAVE_AS_NODE_TEXT = "Save as ...";
    
    private final Class<TOP> primary_item_class;

    /**
     * Main element in the panel. 
     * Left-hand side contains the browser, the right-hand side 
     * shows the information associated with the selected 
     * item. 
     */
    private JSplitPane content;

    /**
     * The browser is implemented as a JTree.
     */
    private JTree item_browser;
    
    /**
     * Root node for the browser.
     */
    private DefaultMutableTreeNode item_browser_root;
    
    /**
     * Data model for the browser.
     */
    private DefaultTreeModel item_browser_model;    
    
    /**
     * The text field where messages are written when there is no selected item.
     */
    private JEditorPane no_selection_text;
    
    public static JSplitPane defaultSplitPane()
    {
    	JSplitPane content = new JSplitPane();
        content = new JSplitPane();
        content.setDividerLocation(200+content.getInsets().top);
        content.setOrientation(JSplitPane.HORIZONTAL_SPLIT);
        content.setDebugGraphicsOptions(javax.swing.DebugGraphics.NONE_OPTION);
        content.setBorder(null);
        content.setResizeWeight(0.0);
        content.setOneTouchExpandable(true);
        return content;
    }
    
    private void initComponents(JSplitPane container)
    {
    	int div = container.getDividerLocation();
    	this.content = container;

        item_browser_root = new DefaultMutableTreeNode("browser root");
        item_browser_model = new DefaultTreeModel(item_browser_root);
        item_browser = new JTree(item_browser_model);
                
        DefaultTreeCellRenderer renderer = new DefaultTreeCellRenderer() 
        {
            @Override
            public java.awt.Component getTreeCellRendererComponent(
                            JTree tree,
                            Object value,
                            boolean sel,
                            boolean expanded,
                            boolean leaf,
                            int row,
                            boolean hasFocus) 
            {
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
        	
        };
        item_browser.setCellRenderer(renderer);
        
        item_browser.setRootVisible(false);
        item_browser.setShowsRootHandles(true);
        item_browser.setEditable(false);
        item_browser.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);

        item_browser.addTreeSelectionListener(e->showSelectedItem());
//        		
//        		new TreeSelectionListener() // when tree node selection changes, we update the content panel
//        {
//            @Override
//            public void valueChanged(TreeSelectionEvent treeSelectionEvent) 
//            {
//                showSelectedItem();
//            }
//        });
        item_browser.addMouseListener(new MouseAdapter() // make the nodes respond to popup mouse clicks
            {
                @Override
                public void mousePressed(MouseEvent e)
                {
                        popupViewMenu(e);
                }
                
                @Override
                public void mouseReleased(MouseEvent e)
                {
                        popupViewMenu(e);
                }
                
            });
        
        item_browser.setFont(new java.awt.Font("Serif", Font.BOLD,14));
        item_browser.setFocusable(true); // is true by default on JTree
        item_browser.setVisible(true);
        
        JScrollPane item_browser_scroll = new JScrollPane();
        item_browser_scroll.setViewportView(item_browser);
        item_browser_scroll.setMinimumSize(new Dimension(200,200));
        content.setLeftComponent(item_browser_scroll);        

        no_selection_text = new JEditorPane();
        no_selection_text.setContentType("text/html");
        no_selection_text.setBackground(getBackground().brighter());
        no_selection_text.setEditable(false);
        content.setRightComponent(no_selection_text);

        content.setDividerLocation(div);
        showSelectedItem();

        
    }

    /** 
     * Finds the selected tree node.
     * 
     * @return
     */
    private DefaultMutableTreeNode getSelectedNode()
    {
        TreePath selected = item_browser.getSelectionPath();
        if (selected == null)
            return null;
        else 
        {
            DefaultMutableTreeNode node = (DefaultMutableTreeNode)selected.getLastPathComponent();
            return node;
        }
    }

    /**
     * Finds the closest ancestor node that contains a primary item.
     * @return null, if none of the ancestors are of primary class.
     */
    private DefaultMutableTreeNode getSelectedPrimaryNode()
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
    
    public JComponent getBrowserComponent()
    {
    	return item_browser;
    }
    
    /**
     * Returns the selected item (associated with {@link #getSelectedNode()}).
     * 
     * @return null if no selection
     */
    public JComponent getSelectedItem()
    {
        DefaultMutableTreeNode selected_node = getSelectedNode();
        if (selected_node==null)
            return null;
        else
        {
            return (JComponent)selected_node.getUserObject();
        }
    }
    
    public TOP getSelectedPrimaryItem()
    {
    	DefaultMutableTreeNode primary_node = getSelectedPrimaryNode();
    	if (primary_node==null)
    		return null;
    	else
    		return primary_item_class.cast(primary_node.getUserObject());
    }
    
    /**
     * All members of {@link #primary_item_class} in the tree, in preorder enumeration.
     * 
     * @return
     */
    public List<TOP> primaryItemList()
    {
    	ArrayList<TOP> item_list = new ArrayList<>();
    	Enumeration<TreeNode> node_enum = item_browser_root.preorderEnumeration();
    	while (node_enum.hasMoreElements())
    	{
    		DefaultMutableTreeNode node = (DefaultMutableTreeNode) node_enum.nextElement();
        	Object payload = node.getUserObject();
        	if (payload!=null && primary_item_class.isInstance(payload))
    	    		item_list.add(primary_item_class.cast(payload));
    	}
    	return item_list;
    }
    
    
    
    /**
     * Adds a new item at the top level.
     *
     * @param item
     * @throws IllegalArgumentException if the item does not belong to the primary item class.
     */
    public void addTopItem(TOP item)
    {
        if (!primary_item_class.isInstance(item))
            throw new IllegalArgumentException("Top level items must be of class "+primary_item_class);
        addItem(item, item_browser_root);
    }
    
    public int getTopItemCount()
    {
    	return item_browser_root.getChildCount();
    	
    }
    
    public TOP getTopItem(int idx)
    {
    	DefaultMutableTreeNode tree_node = (DefaultMutableTreeNode) item_browser_root.getChildAt(idx);
    	return primary_item_class.cast(tree_node.getUserObject());
    }
    
    /** 
     * Adds a new item as the child of the selected primary item.
     * @param item
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
    private void addItem(JComponent item, DefaultMutableTreeNode parent)
    {
        DefaultMutableTreeNode view_node = new DefaultMutableTreeNode(item);
        item_browser_model.insertNodeInto(view_node,parent,parent.getChildCount());
        TreePath view_path = new TreePath(view_node.getPath());

        item_browser.scrollPathToVisible(view_path);
        item_browser.setSelectionPath(view_path);
    }
    
    /**
     * Selects a given item in the browser. 
     * 
     * @param item an item that was added to the browser before
     * @return true if the item was found, or false if not
     */
    public boolean selectItem(JComponent item)
    {
        Enumeration<TreeNode> item_enum = item_browser_root.postorderEnumeration();
        
        while (item_enum.hasMoreElements())
        {
            TreeNode visited_item = item_enum.nextElement();
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
        item_browser.setSelectionPath(null);
        //System.out.println("#*B.sI not found "+item);
        return false;
    }
    
    
    public void clearSelection()
    {
    	item_browser.clearSelection();
    }
    
    
    /**
     * Simple base class for TOP (JPanel). 
     * @deprecated
     */
    public static class PrimaryItem extends JPanel {}    
    
    private void removeNode(DefaultMutableTreeNode node)
    {
        boolean doRemove = true;
        java.util.Enumeration<?> enum_subtree = node.depthFirstEnumeration();
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
            String msg ; 
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
            showSelectedItem();
        }
    }

    
    /**
     * The message displayed when there are no available items.
     */
    private static final String NO_AVAILABLE_ITEMS = "<em>No available items.</em>";
    
    /**
     * The message displayed when there is no selection.
     */
    private static final String SELECT_AN_ITEM = "<em>Select an item.</em>";
    
    /**
     * Called whenever item selection changes
     */
    private void showSelectedItem()
    {
        JComponent selected = getSelectedItem();
        int div = content.getDividerLocation(); // we will reset so that the divider doesn't change 
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
            content.setRightComponent(selected);
            content.setDividerLocation(div);
        }        
    }    
    
    /**
     * Called when a mouse popup trigger is fired on a tree node
     */
    private void popupViewMenu(MouseEvent e)
    {
    	class PopupActionsOnNode implements ActionListener
        {
            private final DefaultMutableTreeNode node;
            private PopupActionsOnNode(DefaultMutableTreeNode node)
            {
                this.node=node;
            }
            
            @Override
            public void actionPerformed(ActionEvent E)
            {
                JMenuItem src = (JMenuItem) E.getSource();
                JComponent view = (JComponent)node.getUserObject();
                
                String menu_text = src.getText();
                if (REMOVE_NODE_TEXT.equals(menu_text))
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
                            ((ExportableData)view).saveData(export_file);
                        } catch (IOException exc)
                        {
                            AppFrame app = Session.getApp(Browser.this);
                            app.getExceptionHandler().handle(exc, "I/O error", "Error while exporting data view.");
                        }
                    }
                } else if (SAVE_NODE_TEXT.equals(menu_text))
                {
                    try {
                        SavableData<?> Sview = (SavableData<?>)view;
                        Sview.saveData(Sview.getDataFile().getFile());
                    } catch (IOException exc)
                    {
                    	AppFrame app = Session.getApp(Browser.this);
                        app.getExceptionHandler().handle(exc, "I/O error", "Error while saving data.");
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
                            SavableData<?> Sview = (SavableData<?>)view;
                            File save_file = new File(directory, file_name); 
                            Sview.saveData(save_file);
                            Browser.this.repaint(); // since name may have changed
                        } catch (IOException exc)
                        {
                        	AppFrame app = Session.getApp(Browser.this);
                            app.getExceptionHandler().handle(exc, "I/O error", "Error while saving data.");
                        }
                    }
                }
            }
        }    	
    	
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
                
                if (view instanceof SavableData)
                {
                    SavableData<?> Sview = (SavableData<?>)view;
                    DataFile<?> data_file = Sview.getDataFile();
                    if (data_file.isDirty())
                    {
                        tlabel.setText("* "+view.toString());
                    } else 
                    {
                        tlabel.setText("  "+view.toString());
                    }
                    JMenuItem save = new JMenuItem(SAVE_NODE_TEXT);
                    save.setEnabled(data_file.getFile().getParent()!=null && data_file.isDirty());
                    save.addActionListener(popup_actions);
                    JMenuItem save_as = new JMenuItem(SAVE_AS_NODE_TEXT);
                    save_as.addActionListener(popup_actions);

                    popup_options.addSeparator();
                    popup_options.add(save);
                    popup_options.add(save_as);
                } else if (view instanceof ExportableData)
                {
                    JMenuItem export = new JMenuItem(EXPORT_NODE_TEXT);
                    export.addActionListener(popup_actions);
                    popup_options.addSeparator();
                    popup_options.add(export);
                }

                popup_options.show(item_browser,e.getX(),e.getY());
            }
        }
    }
    
    /**
     * Adds a change listener that responds to changes of 
     * node selection or structure change. 
     * 
     * @param listener
     */
    public void addSelectionChangeListener(ChangeListener listener)
    {
        class ChangeAdapter implements TreeSelectionListener, TreeModelListener
        {
            private final ChangeListener listener;
            
            ChangeAdapter(ChangeListener L)
            {
                listener = L;
            }

            private void structureChanged()
            {
                JComponent src = getSelectedPrimaryItem();
                if (src == null)
                    src = Browser.this;
                
                listener.stateChanged(new ChangeEvent(src));
            }
            
            @Override
            public void treeNodesChanged(TreeModelEvent e)
            {
                // nothing
            }
            
            @Override
            public void treeNodesInserted(TreeModelEvent e)
            {
                structureChanged();
            }
            @Override
            public void treeNodesRemoved(TreeModelEvent e)
            {
                structureChanged();
            }
            @Override
            public void treeStructureChanged(TreeModelEvent e)
            {
                // nothing
            }
            
            @Override
            public void valueChanged(TreeSelectionEvent treeSelectionEvent) 
            {
                TOP daddy = getSelectedPrimaryItem();
                if (daddy != null)
                    listener.stateChanged(new ChangeEvent(daddy));
            }    	
        }
        ChangeAdapter pisl = new ChangeAdapter(listener);
        item_browser.addTreeSelectionListener(pisl);
        item_browser_model.addTreeModelListener(pisl);
    }
}
