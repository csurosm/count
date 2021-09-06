
package ca.umontreal.iro.evolution.malin.ui.count;

/**
 * A Swing component for showing the results of ancestral 
 * reconstruction by Dollo parsimony.
 * 
 * @author csuros
 */


import java.io.IOException;
import java.io.File;
import java.io.PrintStream;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics2D;

import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;

import ca.umontreal.iro.evolution.genecontent.NodeWithRates;
import ca.umontreal.iro.evolution.genecontent.OccurrenceTable;
import ca.umontreal.iro.evolution.genecontent.PhyleticProfile;
import ca.umontreal.iro.evolution.genecontent.TreeWithRates;

import ca.umontreal.iro.evolution.malin.DataFile;

import ca.umontreal.iro.evolution.malin.ui.EmbellishedTreePanel;

public class DolloDisplay extends AncestralReconstructionPane 
{
    /**
     * Instantiates a new Dollo parsimony browser. 
     * 
     * @param main_tree the underlying phylogeny
     * @param data_file the occurrence table structure on which the ancestral reconstruction will be carried out
     */
    public DolloDisplay(TreeWithRates main_tree, DataFile<OccurrenceTable> data_file, WorkSpaceCount W)
    {
        super(main_tree, data_file, W);
        init();
    }
    
    /**
     * Index of the node at which a family first appears according to Dollo
     */
    protected int[] family_first_at_node;
            
    @Override
    protected void initReconstructionVariables()
    {
        int num_nodes = tree_nodes.length;
        int num_families = families.length;
                
        reconstructed_family_count = new double[num_families][num_nodes];
        for (int j=0; j<reconstructed_family_count.length ; j++)
            reconstructed_family_count[j] = null;
        family_first_at_node = new int[num_families];
    }
    
    @Override
    protected double getReconstructedFamilyGain(int family_idx, int node_idx)
    {
        int parent_idx = main_tree.getParentIndex(node_idx);
        if (reconstructed_family_count[family_idx][node_idx]>0. && reconstructed_family_count[family_idx][parent_idx]==0.0)
            return 1.0;
        else
            return 0.0;
    }
   
    @Override
    protected double getReconstructedFamilyLoss(int family_idx, int node_idx)
    {
        int parent_idx = main_tree.getParentIndex(node_idx);
        if (reconstructed_family_count[family_idx][node_idx]==0. && reconstructed_family_count[family_idx][parent_idx]>0.0)
            return 1.0;
        else
            return 0.0;
    }
    
    
    @Override
    protected AncestralReconstructionPane.ComputingTask newComputingTask()
    {
        return new DolloComputingTask();
    }
    
    protected class DolloComputingTask extends AncestralReconstructionPane.ComputingTask
    {
        @Override
        protected void computeAncestralReconstruction(int family_idx)
        {
            int num_nodes = tree_nodes.length;

            reconstructed_family_count[family_idx] = new double[num_nodes];
            OccurrenceTable data_table = data_file.getData();
            PhyleticProfile profile = data_table.getProfile(family_idx);
            int[] rec = profile.computeDolloParsimony(main_tree);
            for (int node_idx=0; node_idx<num_nodes; node_idx++)
                reconstructed_family_count[family_idx][node_idx] = (double)rec[node_idx];

            family_first_at_node[family_idx] = -1;
            for (int node_idx=0; node_idx<num_nodes; node_idx++)
            {
                NodeWithRates N = tree_nodes[node_idx];
                if (!N.isRoot())
                {
                    int parent_idx = main_tree.getParentIndex(node_idx);
                    if (rec[node_idx]>0 && rec[parent_idx]==0)
                        family_first_at_node[family_idx] = node_idx;
                } else if (rec[node_idx]>0)
                    family_first_at_node[family_idx] = node_idx;
            }
        }
    }
    
    @Override
    protected double getReconstructedFamilyMulti(int family_idx, int node_idx)
    {
        return 0.0;
        //NodeWithRates N = tree_nodes[node_idx];
        //if (N.isLeaf())
        //{
        //    if (this.family_size[family_idx][node_idx]>1)
        //        return 1.0;
        //    else
        //        return 0.0;
        //} else
        //    return 0.0;
    }
    
    @Override
    protected HistoryTreePanel createTreePanel()
    {
        return new DolloTreePane();
    }

    @Override
    protected OccurrenceTableModel createFamilyTableModel()
    {
        return new FamilyTableModel();
    }
        
    @Override
    protected AncestralReconstructionPane.LineageTableModel createLineageTableModel()
    {
        return new LineageTableModel();
    }

    @Override
    protected void initComponents()
    {
        super.initComponents();
        family_table.setDefaultRenderer(FamilySizeChange.class, new SizeChangeRenderer());
    }    

    @Override
    public String toString()
    {
        return "Dollo parsimony";
    }
    
    @Override
    public void saveData(File f) throws IOException
    {
        PrintStream PS = new PrintStream(f);
        PS.println(WorkSpaceCount.getDealerCount(this).getStandardHeader(getClass()));
        
        OccurrenceTableModel table_model = family_table_scroll.getModel();
        int num_columns = table_model.getColumnCount();
        
        // header
        PS.print("Family");
        int num_nodes = main_tree.getNumNodes();
        for (int node_idx=0; node_idx<num_nodes;node_idx++)
        {
            PS.print("\t"+getShortNodeName(node_idx));
        }
        PS.println();
        
        int num_families = families.length;
        for (int family_idx=0; family_idx<num_families; family_idx++)
        {
            PS.print(families[family_idx]);
            for (int node_idx=0; node_idx<num_nodes;node_idx++)
            {
                PS.print("\t");
                NodeWithRates N = tree_nodes[node_idx];
                if (N.isLeaf())
                    PS.print(family_size[family_idx][node_idx]);
                else
                    PS.print((int)reconstructed_family_count[family_idx][node_idx]);
            }
            PS.println();
        }
    }    
    
    private class FamilySizeChange implements Comparable<FamilySizeChange>
    {
        private int value; // 0, 1 or -1
        private FamilySizeChange(int node_value, int parent_value)
        {
            this.value=node_value+parent_value*2;
        }
        public int compareTo(FamilySizeChange x){return (new Integer(value)).compareTo(x.value);}
        @Override
        public String toString()
        {
            if (value==0) return "..";
            else if (value==1) return "gain";
            else if (value==2) return "loss";
            else if (value==3) return "1:1";
            else return "?";
        }
        @Override 
        public boolean equals(Object o)
        {
            boolean retval = (o instanceof FamilySizeChange)
                    && (((FamilySizeChange)o).value == value);
            return retval;
        }
    }
    
    private class SizeChangeRenderer extends DefaultTableCellRenderer
    {
        @Override
        public Component getTableCellRendererComponent(JTable table,
                                               Object value,
                                               boolean isSelected,
                                               boolean hasFocus,
                                               int row,
                                               int column)        
        {
            Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
            if (value == null)
            {
                setToolTipText("");
            } else
            {
                if (value instanceof FamilySizeChange) // it'd better be...
                {
                    FamilySizeChange C = (FamilySizeChange)value;
                    if (C.value == 3) 
                        setToolTipText("family conserved along the lineage");
                    else if (C.value == 0)
                        setToolTipText("family absent along the lineage");
                    else if (C.value==1)
                        setToolTipText("family originating here");
                    else if (C.value==2)
                        setToolTipText("lineage-specific family loss");
                    else 
                        setToolTipText("[impossible FamilySizeChange value "+C.value+"]");    
                    
                } else
                {
                    setToolTipText("[not a FamilySizeChange] "+value);
                }
            }
            return comp;
        }
    }    
    
    
    
    protected class DolloTreePane extends AncestralReconstructionPane.HistoryTreePanel
    {
        public DolloTreePane()
        {
            super();
        }
        
        protected DolloTreePane(NodeWithRates root, EmbellishedTreePanel.LayoutStyle layout)
        {
            super(root, layout);
        }

        @Override
        protected void plotIndividualPresence(Graphics2D g2, int x, int y, int width, int height, int family_idx, int node_idx)
        {
            if (reconstructed_family_count[family_idx]!=null)
            {
                Color C = g2.getColor();

                double rec_c = getReconstructedFamilyCount(family_idx,node_idx);
                
                int ht = (int)(rec_c*height+0.5);
                if (ht>0)
                {
                    if (node_idx < family_size[family_idx].length && family_size[family_idx][node_idx]>1)
                        g2.setColor(LookAndFeel.MULTI_PRESENCE_COLOR);
                    else
                        g2.setColor(LookAndFeel.SINGLE_PRESENCE_COLOR);
                    g2.fillRect(x, y-height/2, width, ht);
                }

                if (family_first_at_node[family_idx]==node_idx)
                {
                    g2.setColor(LookAndFeel.GAIN_COLOR);
                    g2.setStroke(medium_stroke);
                } else if (!main_tree.getNode(node_idx).isRoot() && getReconstructedFamilyLoss(family_idx,node_idx)==1.0)
                {
                    g2.setColor(LookAndFeel.LOSS_COLOR);
                    g2.setStroke(medium_stroke);
                }
                else
                {
                    g2.setColor(LookAndFeel.MULTI_PRESENCE_COLOR);
                    g2.setStroke(thin_stroke);
                }
                    
                g2.drawRect(x-1, y-height/2-1, width+1, 1+height);
                g2.setColor(C);
            }
        }
   
        
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
            int num_nodes = main_tree.getNumNodes();
            int num_edges = main_tree.getNumEdges();
            int num_leaves = main_tree.getNumLeaves();
            return super.getColumnCount()
                    + (num_nodes-num_leaves)
                    +num_edges;
        }
        
        @Override
        public Object getValueAt(int row_idx, int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            int num_leaves = main_tree.getNumLeaves();
            if (column_idx<num_usual_columns)
                return super.getValueAt(row_idx, column_idx);
            else 
                column_idx = column_idx-num_usual_columns+num_leaves;

            int num_nodes = main_tree.getNumNodes();
            if (column_idx<num_nodes)
            {
                if (main_tree.getNode(column_idx).isLeaf())
                    return new Integer(family_size[row_idx][column_idx]);
                else
                    return new Integer((int)reconstructed_family_count[row_idx][column_idx]);
            } else
                column_idx -= num_nodes;
            // now it's edge index
            int r = (int)(getReconstructedFamilyGain(row_idx, column_idx)-getReconstructedFamilyLoss(row_idx, column_idx));
            int parent_idx = main_tree.getParentIndex(column_idx);
            return new FamilySizeChange((int)getReconstructedFamilyCount(row_idx, column_idx), (int)getReconstructedFamilyCount(row_idx, parent_idx));
        }
        
        @Override
        public String getColumnName(int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            int num_leaves = main_tree.getNumLeaves();
            if (column_idx<num_usual_columns)
                return super.getColumnName(column_idx);
            else 
                column_idx = column_idx-num_usual_columns+num_leaves;

            int num_nodes = main_tree.getNumNodes();
            
            if (column_idx<num_nodes)
                return getShortNodeName(column_idx);
            else
                column_idx -= num_nodes;
            
            return "\u0394"+getShortNodeName(column_idx); // Delta
        }
        
        @Override 
        public String getCellToolTip(int row_idx, int column_idx)
        {
            Object value = getValueAt(row_idx,column_idx);
            int num_usual_columns = super.getColumnCount();
            int num_leaves = main_tree.getNumLeaves();
            if (column_idx<num_usual_columns)
                return super.getCellToolTip(row_idx,column_idx);
            else 
                column_idx = column_idx-num_usual_columns+num_leaves;

            int num_nodes = main_tree.getNumNodes();
            if (column_idx<num_nodes) 
            {
                int cnt = ((Integer)value).intValue();
                String presence_info = (cnt==0?"was not present":"was present");
                return "Family "+families[row_idx]+" "+presence_info+" at node "+getLongNodeName(column_idx)
                        +" by Dollo parsimony reconstruction";
            } else 
                column_idx -= num_nodes;
            FamilySizeChange C = (FamilySizeChange)value;
            String change_info = "Family "+families[row_idx];
            if (C.value == 3)
                change_info += " was conserved along the lineage leading to node "+getLongNodeName(column_idx);
            else if (C.value == 0)
                change_info += " was absent along the lineage leading to node "+getLongNodeName(column_idx);
            else if (C.value==1)
                change_info += " appeared first at node "+getLongNodeName(column_idx);
            else if (C.value==2)
                change_info += " was lost in the lineage leading to node "+getLongNodeName(column_idx);
            
            return change_info; // Delta
        }
        
        @Override
        public Class getColumnClass(int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            int num_leaves = main_tree.getNumLeaves();
            if (column_idx<num_usual_columns)
                return super.getColumnClass(column_idx);
            else 
                column_idx = column_idx-num_usual_columns+num_leaves;

            int num_nodes = main_tree.getNumNodes();

            if (column_idx < 2+num_nodes)
                return Integer.class;
            else 
                return FamilySizeChange.class;
        }
        
        @Override
        public String getColumnHeaderToolTip(int column_idx)        
        {
            int num_usual_columns = super.getColumnCount();
            int num_leaves = main_tree.getNumLeaves();
            if (column_idx<num_usual_columns)
                return super.getColumnHeaderToolTip(column_idx);
            else 
                column_idx = column_idx-num_usual_columns+num_leaves;

            int num_nodes = main_tree.getNumNodes();
            
            if (column_idx<num_nodes)
                return "Family presence/absence (1/0) by Dollo parsimony at node "+getLongNodeName(column_idx);
            else 
                column_idx -= num_nodes;
            
            return "Event by Dollo parsimony in the lineage leading to "+getLongNodeName(column_idx);
        }
        
    }

    
    protected class LineageTableModel extends AncestralReconstructionPane.LineageTableModel
    {
        /**
         * We have no expansion / contraction columns
         * 
         * @return 5 (node, present, gain, loss)
         */
        @Override
        public int getColumnCount()
        {
            return 4;
        }

        @Override
        public String getColumnName(int column_idx)
        {
            if (column_idx<2)
                return super.getColumnName(column_idx);
            else
                return super.getColumnName(1+column_idx);
        }
        
        @Override
        public String getHeaderToolTip(int column_idx)
        {
            if (column_idx<2)
                return super.getHeaderToolTip(column_idx);
            else
                return super.getHeaderToolTip(1+column_idx);
        }

        @Override
        public Object getValueAt(int row_idx, int column_idx)
        {
            if (column_idx==0)
                return getLongNodeName(row_idx);
            if (column_idx==LINEAGE_TABLE_COLUMN_FAMILY_COUNT)
                return new Integer((int)selected_family_count[row_idx]);
            if (tree_nodes[row_idx].isRoot())
                return null;
            if (1+column_idx==LINEAGE_TABLE_COLUMN_FAMILY_GAINS)
                return new Integer((int)selected_family_gain[row_idx]);
            if (1+column_idx==LINEAGE_TABLE_COLUMN_FAMILY_LOSSES)
                return new Integer((int)selected_family_loss[row_idx]);
            // should never get here
            return null;
        }   
        
        @Override
        public Class getColumnClass(int column_idx)
        {
            if (column_idx == 0) 
                return String.class; 
            else
                return Integer.class;
        }
        
    }
    

}
