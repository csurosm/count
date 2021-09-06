
package ca.umontreal.iro.evolution.malin.ui.count;

/**
 * A Swing component for showing the results of ancestral 
 * reconstruction by Wagner parsimony.
 * 
 * @author csuros
 */


import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import ca.umontreal.iro.evolution.genecontent.OccurrenceTable;
import ca.umontreal.iro.evolution.genecontent.PhyleticProfile;
import ca.umontreal.iro.evolution.genecontent.TreeWithRates;

import ca.umontreal.iro.evolution.malin.DataFile;

public class WagnerDisplay extends AncestralReconstructionPane implements ChangeListener
{
    private static final double DEFAULT_GAIN_PENALTY = 1.0;
    
    public WagnerDisplay(TreeWithRates main_tree, DataFile<OccurrenceTable> data_file, double gain_penalty, WorkSpaceCount W)
    {
        super(main_tree, data_file, W);
        this.gain_penalty = gain_penalty;
        init();
        
    }
    public WagnerDisplay(TreeWithRates main_tree, DataFile<OccurrenceTable> data_file, WorkSpaceCount W)
    {
        this(main_tree, data_file, DEFAULT_GAIN_PENALTY, W);
    }
    
    
    /**
     * Gain penalty for asymmetric Wagner parsimony.
     */
    private double gain_penalty;
    
    /**
     * Reconstruction by Wagner parsimony
     */
    private int[][] optimal_family_sizes;
    
    private JSpinner gain_penalty_spinner;

    /**
     * Short description of this instance
     * 
     * @return a String that can be used in the browser
     */
    @Override
    public String toString()
    {
        return "Wagner parsimony";
    }
    
    
    public double getGainPenalty()
    {
        return gain_penalty;
    }
    
    @Override
    protected void initReconstructionVariables()
    {
        int num_nodes = tree_nodes.length;
        int num_families = families.length;
                
        reconstructed_family_count = new double[num_families][num_nodes];
        for (int j=0; j<reconstructed_family_count.length ; j++)
            reconstructed_family_count[j] = null;
        
        optimal_family_sizes = new int[num_families][];
    }
    
    @Override
    protected FamilyTableModel createFamilyTableModel()
    {
        return new FamilyTableModel(); // our version
    }

    @Override
    protected double getReconstructedFamilyMulti(int family_idx, int node_idx)
    {
        if (optimal_family_sizes[family_idx][node_idx]>1)
            return 1.0;
        else
            return 0.0;
    }
    
    @Override
    protected double getReconstructedFamilyGain(int family_idx, int node_idx)
    {
        int parent_idx = main_tree.getParentIndex(node_idx);
        if (optimal_family_sizes[family_idx][node_idx]>0 && optimal_family_sizes[family_idx][parent_idx]==0)
            return 1.0;
        else
            return 0.0;
    }
    
    @Override
    protected double getReconstructedFamilyLoss(int family_idx, int node_idx)
    {
        int parent_idx = main_tree.getParentIndex(node_idx);
        if (optimal_family_sizes[family_idx][node_idx]==0 && optimal_family_sizes[family_idx][parent_idx]>0)
            return 1.0;
        else
            return 0.0;
    }
    
    
    @Override
    protected double getReconstructedFamilyExpansion(int family_idx, int node_idx)
    {
        int parent_idx = main_tree.getParentIndex(node_idx);
        if (optimal_family_sizes[family_idx][node_idx]>1 && optimal_family_sizes[family_idx][parent_idx]==1)
            return 1.0;
        else
            return 0.0;
    }
    
    @Override
    protected double getReconstructedFamilyReduction(int family_idx, int node_idx)
    {
        int parent_idx = main_tree.getParentIndex(node_idx);
        if (optimal_family_sizes[family_idx][node_idx]==1 && optimal_family_sizes[family_idx][parent_idx]>1)
            return 1.0;
        else
            return 0.0;
    }

    @Override
    protected void createBottomBarLeft(Box bb)
    {
        super.createBottomBarLeft(bb);
        bb.add(createGainPenaltyPanel());
    }
    
    private JPanel createGainPenaltyPanel()
    {
        JPanel gain_penalty_panel = new JPanel(new GridLayout(1,0,20,5));
        gain_penalty_panel.setBackground(null);
        gain_penalty_panel.setBorder(BorderFactory.createEtchedBorder());
        gain_penalty_panel.setFont(new Font("Serif", Font.PLAIN, LookAndFeel.TREE_PANEL_FONT_SIZE));

        JLabel gain_penalty_info = new JLabel("Gain penalty");
        gain_penalty_info.setFont(gain_penalty_panel.getFont().deriveFont(Font.BOLD));
        gain_penalty_panel.add(gain_penalty_info);
        
        SpinnerNumberModel gain_penalty_model = new SpinnerNumberModel(gain_penalty, 0.0, Double.POSITIVE_INFINITY, 0.1);
        gain_penalty_spinner = new JSpinner(gain_penalty_model);
        gain_penalty_spinner.setFont(gain_penalty_panel.getFont());
        gain_penalty_spinner.setBackground(gain_penalty_panel.getBackground());
        gain_penalty_spinner.setToolTipText("Current penalty for asymmetric Wagner parsimony " +
                "(loss penalty equals 1)");
        final JFormattedTextField gps_text = ((JSpinner.DefaultEditor)gain_penalty_spinner.getEditor()).getTextField();
        
        gain_penalty_spinner.addChangeListener(this);
        gain_penalty_spinner.addChangeListener(new ChangeListener()
            {
                public void stateChanged(ChangeEvent ignored)
                {
                    gps_text.transferFocus(); // or else the cursor keeps blinking in there
                }
            });
        gps_text.setFont(gain_penalty_panel.getFont());
        gps_text.setBackground(LookAndFeel.GAIN_COLOR.brighter().brighter());
        gain_penalty_panel.add(gain_penalty_spinner);
        gain_penalty_panel.setMaximumSize(new Dimension(160,30));
        
        return gain_penalty_panel;
    }
    
    @Override
    protected void initComputeAll()
    {
        gain_penalty_spinner.setEnabled(false);
        //System.out.println("#*WD.iCA init ");
        super.initComputeAll();
    }
    
    @Override
    protected void finishComputeAll()
    {
        super.finishComputeAll();
        gain_penalty_spinner.setEnabled(true);
    }
    
    @Override
    protected AncestralReconstructionPane.ComputingTask newComputingTask()
    {
        return new WagnerComputingTask();
    }
    
    protected class WagnerComputingTask extends AncestralReconstructionPane.ComputingTask
    {
        @Override
        protected void computeAncestralReconstruction(int family_idx)
        {
            int num_nodes = tree_nodes.length;

            reconstructed_family_count[family_idx] = new double[num_nodes];

            OccurrenceTable data_table = data_file.getData();
            PhyleticProfile profile = data_table.getProfile(family_idx);
            int[] rec = profile.computeWagnerParsimony(main_tree, gain_penalty);
            optimal_family_sizes[family_idx] = rec;
            for (int node_idx=0; node_idx<num_nodes; node_idx++)
            {
                if (rec == null )
                    System.out.println("#**WD:cAR now rec is null [family = "+family_idx+", node="+node_idx+"]");
                if (reconstructed_family_count==null)
                    System.out.println("#**WD:cAR now rfc is null [family = "+family_idx+", node="+node_idx+"]");
                if (reconstructed_family_count[family_idx]==null)
                    System.out.println("#**WD:cAR now rfc[fi] is null [family = "+family_idx+", node="+node_idx+"]");
                reconstructed_family_count[family_idx][node_idx] = (rec[node_idx]==0?0.0:1.0);
            }
        }
        
    }

    
    /**
     * Called when gain penalty spinner's value changes
     * @param e
     */
    public void stateChanged(ChangeEvent e)
    {
        if (e.getSource() instanceof JSpinner)
        {
            JSpinner spinner = (JSpinner)e.getSource();
            double new_penalty = ((Number)spinner.getValue()).doubleValue();
            gain_penalty=new_penalty;
            computeAll();
        }
    }
    
    @Override
    public void saveData(File f) throws IOException
    {
        PrintStream PS = new PrintStream(f);
        DealerCount dc = WorkSpaceCount.getDealerCount(this);
        PS.println(dc.getStandardHeader(getClass()));
        PS.println(dc.getStandardHeader("Gain penalty: "+gain_penalty));
        
        OccurrenceTableModel table_model = family_table_scroll.getModel();
        int num_columns = table_model.getColumnCount();
        
        // header
        PS.print("# Family");
        PS.print("\tGains\tLosses\tExpansions\tContractions");

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
            double n_gains = 0.0;
            double n_losses = 0.0;
            double n_expansions = 0.0;
            double n_contractions = 0.0;

            for (int node_idx=0; node_idx<num_nodes-1;node_idx++)
            {
                n_gains += getReconstructedFamilyGain(family_idx, node_idx);
                n_losses += getReconstructedFamilyLoss(family_idx, node_idx);
                n_expansions += getReconstructedFamilyExpansion(family_idx, node_idx);
                n_contractions += getReconstructedFamilyReduction(family_idx, node_idx);
            }
            PS.print("\t"+n_gains+"\t"+n_losses+"\t"+n_expansions+"\t"+n_contractions);
            for (int node_idx=0; node_idx<num_nodes;node_idx++)
            {
                PS.print("\t"+optimal_family_sizes[family_idx][node_idx]);
            }
            PS.println();
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
            int num_usual_columns = super.getColumnCount();
            int num_leaves = main_tree.getNumLeaves();
            int num_nodes  = main_tree.getNumNodes();
            return num_usual_columns
                    +(num_nodes-num_leaves);
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
            // column_index is a node index now
            if (optimal_family_sizes[row_idx]==null)
                return new Integer(0);
            else
                return new Integer(optimal_family_sizes[row_idx][column_idx]);
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

            return getShortNodeName(column_idx);    
        }

        @Override
        public Class getColumnClass(int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            if (column_idx < num_usual_columns)
                return super.getColumnClass(column_idx);
            else
                return Integer.class;
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

            return "Ancestral family size inferred by Wagner parsimony at node "+getLongNodeName(column_idx);
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

            int v = ((Integer)value).intValue();
            
            String size_info = (v==0?"no":Integer.toString(v));
            return "Family "+families[row_idx]+" had "+size_info+" member"+(v==1?"":"s")+" at node "+getLongNodeName(column_idx)+" according to Wagner parsimony";
        }
    
    }
}
