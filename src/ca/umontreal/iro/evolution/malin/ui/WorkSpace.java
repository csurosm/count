/*
 * WorkSpace.java
 *
 * Created on November 13, 2007, 1:34 PM
 */

package ca.umontreal.iro.evolution.malin.ui;


import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JScrollPane;
import javax.swing.JOptionPane;
import javax.swing.JEditorPane;

import javax.swing.event.ChangeListener;

import java.awt.Component;
import java.awt.GridLayout;

import java.io.File;

import ca.umontreal.iro.evolution.TreeNode;
import java.awt.Color;

/**
 * Class for the workspace associated with a session. Each session has a unique 
 * phylogeny, but may have multiple data and rate sets. 
 *
 * @author  csuros
 */

public class WorkSpace extends JPanel
{
    
    /**
     * Only used by extending classes
     */
    protected WorkSpace(Dealer D)
    {
        super(new GridLayout(1,1));
        initWorkArea();
        this.dealer = D;
    }
    
    private Dealer dealer;
    
    public Dealer getDealer()
    {
        return dealer;
    }
    
    public void exceptionCaught(Exception E, String title, String more_info)
    {
        dealer.exceptionCaught(E, title, more_info);
    }
    
    /**
     * Normal construction of a new workspace: with a phylogeny and the associated file
     */
    public WorkSpace(Dealer D, TreeNode root, File path_to_tree_file)
    {
        this(D);
        setPhylogeny(root, path_to_tree_file);
    }
    
    protected JScrollPane main_panel_scroll;

    /**
     * Used by extending classes to set the phylogeny after instantiation
     *
     * Here, it adds a Tree panel tab for the main tree, and sets its appearance.
     */
    protected void setPhylogeny(TreeNode root, File path_to_tree_file)
    {
        work_area.removeAll();
        
        tree_panel = new TreePanel(root,false);
        tree_panel.setBackground(java.awt.Color.WHITE);
        tree_panel.setAssociatedFile(path_to_tree_file);
        tree_panel.setToolTipText("Main phylogeny ("+path_to_tree_file.getName()+")");
        main_panel_scroll = new JScrollPane(tree_panel);
        
        addToWorkArea(TREE_PANEL, main_panel_scroll);

    }

    protected void checkPhylogeny()
    {
        // check if there are any problems
        TreeNode[] dft = getPhylogeny().getTraversal().getDFT();
        StringBuffer problems = null;
        for (int node_idx=0; node_idx<dft.length; node_idx++)
        {
            TreeNode N = dft[node_idx];
            if (N.isLeaf())
            {
                String leaf_name = N.getName();
                if (leaf_name.indexOf(' ')>=0)
                {
                    if (problems==null) problems = new StringBuffer("<ul>");
                    problems.append("<li>Leaf name `"+leaf_name+"' has a space: in Newick tree files, underscore " +
                            "within a node name is interpreted as a space --- verify if that is consistent " +
                            "with the table column names, where underscore is left as is.</li>\n");
                }
            } else if (N.getNumChildren()==1)
            {
                if (problems==null) problems = new StringBuffer("<ul>");
                problems.append("<li>Inner node `"+N.newickName()+"' has only "+N.getNumChildren()+" child. " +
                        "This may cause problems " +
                        "in certain analyses, or even a program error. " +
                        "The input tree is supposed to have bi- or multifurcations everywhere, " +
                        "except at the leaves.</li>\n");
            }
            if (N.isRoot() && N.getNumChildren()==2)
            {
                if (problems==null) problems = new StringBuffer("<ul>");
                problems.append("<li>The root has only two children. While it is not necessarily a problem, " +
                        "it may cause ambiguities in certain inference methods, " +
                        "particularly with reversible probabilistic models " +
                        "(Felsenstein's <q>pulley principle</q>).</li>\n");
            }
        }
        if (problems != null)
        {
            problems.append("</ul>");

            JEditorPane problems_pane = new JEditorPane("text/html", "<h1>Possible problems with your tree</h1>"+problems.toString());
            problems_pane.setEditable(false);
            problems_pane.setBackground(new java.awt.Color(255, 255, 192));
            JScrollPane problems_scroll = new JScrollPane(problems_pane);
            problems_scroll.setMaximumSize(new java.awt.Dimension(500,400));
            problems_scroll.setPreferredSize(problems_scroll.getMaximumSize());
            JOptionPane.showMessageDialog(tree_panel,
                    problems_scroll,
                    "Is your tree file correct?",
                    JOptionPane.WARNING_MESSAGE
                    );
        }
    }
    
    protected static final String TREE_PANEL = "Tree"; 
    
    protected JTabbedPane work_area;
    protected TreePanel tree_panel;
    
    protected void addToWorkArea(String name, java.awt.Component C)
    {
        work_area.addTab(name, C);
        //work_area.invalidate();
    }
    
    /**
     * Selects the element that was previously added by its name
     */
    protected void selectWorkAreaElement(String name)
    {
        work_area.setSelectedIndex(work_area.indexOfTab(name));
    }
    
    public String getSelectedWorkAreaElementName()
    {
        return work_area.getTitleAt(work_area.getSelectedIndex());
    }

    protected void initWorkArea()
    {
        work_area = new JTabbedPane();
        //work_area.setPreferredSize(new java.awt.Dimension(900,900));
        add(work_area);
    }
    
    /** 
     * Adds a listener for changes of tab selection in the work area
     * 
     * @param l a ChangeListener 
     */
            
    public void addChangeListener(ChangeListener l)    
    {
        work_area.addChangeListener(l);
        //System.out.println("#*WS.aCL "+l);
    }
            
    public TreeNode getPhylogeny()
    {
        return tree_panel.getRoot();
    }
    
    
    public File getAssociatedFile()
    {
        return tree_panel.getAssociatedFile();
    }
    
    /**
     * Finds the lowest WorkSpace ancestor for a component
     */
    public static WorkSpace getWorkSpace(Component C)
    {
        while (C!=null && !(C instanceof WorkSpace))
        {
            //System.out.println("#*WS.gWS C "+C+"\tprn "+C.getParent());
            C = C.getParent();
        }
        return (WorkSpace) C;
    }
    
    
    /**
     * Finds the lowest Dealer ancestor for a component
     */
    public static Dealer getDealer(Component C)
    {
        while (C != null)
        {
            if (C instanceof WorkSpace)
            {
                WorkSpace ws = (WorkSpace) C;
                return ws.getDealer();
            } 
            C = C.getParent();
        }
        // now if we got here, we got a problem.
        return null;
    }
    
}
