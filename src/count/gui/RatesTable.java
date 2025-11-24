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


import javax.swing.table.AbstractTableModel;

import count.ds.IndexedTree;
import count.ds.TreeTraversal;
import count.gui.kit.RoundedDouble;
import count.gui.kit.TableScroll;
import count.model.Likelihood;
import count.model.TreeWithRates;


/**
*
* A scrollable table showing the rate model's parameters. 
* 
* @author Mikl&oacute;s Cs&#369;r&ouml;s 
*/

public class RatesTable extends TableScroll<RatesTable.Model>
{
	public RatesTable(TreeWithRates rates)
	{
		super(new Model(rates), 2);
	}

//	private void initSetup()
//    {
//        this.setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
//        this.setColumnSelectionAllowed(false);
//        this.setRowSelectionAllowed(true);
//        this.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
//        this.setFont(new Font("Serif",Font.PLAIN,AnnotatedTablePanel.TABLE_FONT_SIZE));
//        this.createDefaultColumnsFromModel();
//        this.getColumnModel().getColumn(1).setPreferredWidth(180); //  column should be wide (node name)
//        this.setDefaultRenderer(RoundedDouble.class, new RoundedDouble.Renderer());
//    }
	
    @Override
    protected int getPreferredColumnWidth(int idx)
    {
        if (idx == 1 ) // node name
            return 180;
        else
            return super.getPreferredColumnWidth(idx);
    }

    @Override
    public String getRowName(int node)
    {
    	return this.getModel().getTree().getIdent(node);
    }

    @Override
    protected String getHeaderToolTip(int column_idx)
    {
        String tt = getModel().getColumnDescription(column_idx);

        return tt +
                "; click to sort rows, drag to rearrange columns";
    }
    
    @Override
    protected void doubleClickInDataTable(int displayed_row_idx, int displayed_column_idx)
    {
    	// do nothing
    }
    
    /**
     * 
	 * Extra attribute columns per node: 2 if depth and height; 1 if depth only; 0 if none 
     */
    private static final int NODE_ATTRIBUTE_COUNT =0 ; //= 2; 
	
    /**
     * Our TableModel that tracks the updates to the underlying model parameters.
     *
     */
	public static class Model extends AbstractTableModel
	{
        private final TreeWithRates rates_model;
        private final Likelihood factory;
        
        public Model(TreeWithRates rates_model)
        {
            this.rates_model = rates_model;
            this.factory = new Likelihood(rates_model);
        }
    
        private IndexedTree getTree(){ return rates_model.getTree();}
        
        @Override
        public void fireTableDataChanged()
        {
        	factory.computeParameters();
        	super.fireTableDataChanged();
        }
        
        @Override
        public int getRowCount()
        {
            return getTree().getNumNodes();
        }

        @Override
        public int getColumnCount()
        {
            return 1 // node index
            	+1   // node name
            	+NODE_ATTRIBUTE_COUNT
            	+3	 // rate parameters
            	+1	 // absence 
            	+3   // GLD parameters
            	+3	// survival parameters
            		;
        }

        @Override
        public Object getValueAt(int node, int column)
        {
        	int this_col = 0;
        	if (column == this_col)
        		return Integer.valueOf(node);
        	++this_col;
        	if (column == this_col)
        	{
        		if (getTree().isLeaf(node))
        			return getTree().getName(node);
        		else 
        			return getTree().getIdent(node);
        	}
        	
        	
        	if (0<NODE_ATTRIBUTE_COUNT)
        	{
            	++this_col;
            	if (column == this_col)
            		return Integer.valueOf(TreeTraversal.depth(getTree(), node));
            	if (1<NODE_ATTRIBUTE_COUNT)
            	{
	            	++this_col;
	            	if (column == this_col)
	            		return Integer.valueOf(TreeTraversal.height(getTree(), node));
            	}
        	}
        	double len = rates_model.getEdgeLength(node);
        	++this_col;
        	if (column == this_col)
        		return new RoundedDouble(rates_model.getLossRate(node)*len);
        	++this_col;
        	if (column == this_col)
        		return new RoundedDouble(rates_model.getDuplicationRate(node)/rates_model.getLossRate(node));
        	++this_col;
        	if (column == this_col)
        		return new RoundedDouble(rates_model.getGainRate(node));
        	++this_col;
        	if (column == this_col)
        		return new RoundedDouble(factory.getExtinction(node));
        	++this_col;
        	if (column == this_col)
        		return new RoundedDouble(rates_model.getLossParameter(node));
        	++this_col;
        	if (column == this_col)
        		return new RoundedDouble(rates_model.getDuplicationParameter(node));
        	++this_col;
        	if (column == this_col)
        		return new RoundedDouble(rates_model.getLinearGainParameter(node));
        	++this_col;
        	if (column == this_col)
        		return new RoundedDouble(factory.getLossParameter(node));
        	++this_col;
        	if (column == this_col)
        		return new RoundedDouble(factory.getDuplicationParameter(node));
        	++this_col;
        	if (column == this_col)
        		return new RoundedDouble(factory.getLinearGainParameter(node)); // factory.getUniversalGainParameter(node));
        	
        	
        	throw new IllegalArgumentException(".getValueAt --- column index must be 0.."+this_col+" [got "+column+"]");        	
        }

        @Override
        public String getColumnName(int column)
        {
        	int this_col = 0;
        	if (column == this_col)
        		return "Order";
        	++this_col;
        	if (column == this_col)
                return "Node";
        	
        	if (0<NODE_ATTRIBUTE_COUNT)
        	{
            	++this_col;
            	if (column == this_col)
            		return "Depth";
            	if (1<NODE_ATTRIBUTE_COUNT)
            	{
	            	++this_col;
	            	if (column == this_col)
	            		return "Height";
            	}
        	}

        	++this_col;
            if (column==this_col)
                return "Length";
        	++this_col;
            if (column==this_col)
                return "Duplication rate";
        	++this_col;
            if (column==this_col)
                return "Gain rate";
        	++this_col;
            if (column==this_col)
                return "Extinction";
        	++this_col;
            if (column==this_col)
                return "p (loss)";
        	++this_col;
            if (column==this_col)
                return "q (duplication)";
        	++this_col;
            if (column==this_col)
                return "r (gain)";
        	++this_col;
            if (column==this_col)
                return "p~ (loss survival)";
        	++this_col;
            if (column==this_col)
                return "q~ (duplication survival)";
        	++this_col;
            if (column==this_col)
                return "r~ (gain survival)";
            
            throw new IllegalArgumentException(".getColumnName --- column index must be 0.."+this_col+" [got "+column+"]");
        }

        @Override
        public Class<?> getColumnClass(int column)
        {
        	int this_col = 0;
        	if (column == this_col)
        		return Integer.class;
        	++this_col;
        	if (column == this_col)
                return String.class;
        	if (0<NODE_ATTRIBUTE_COUNT)
        	{
            	++this_col;
            	if (column == this_col)
            		return Integer.class;
            	if (1<NODE_ATTRIBUTE_COUNT)
            	{
	            	++this_col;
	            	if (column == this_col)
	            		return Integer.class;
            	}
        	}
        	return RoundedDouble.class;
        }	
        
        String getColumnDescription(int column)
        {
        	int this_col = 0;
        	if (column == this_col)
        		return "Node order: leaves first, then by postfix traversal";
        	++this_col;
        	if (column == this_col)
                return "Node name";
        	if (0<NODE_ATTRIBUTE_COUNT)
        	{
            	++this_col;
            	if (column == this_col)
            		return "Node depth (root=0)";
            	if (1<NODE_ATTRIBUTE_COUNT)
            	{
	            	++this_col;
	            	if (column == this_col)
	            		return "Node height (leaf=0)";
            	}
        	}
        	++this_col;
            if (column==this_col)
                return "Edge length/total loss rate (μt)";
        	++this_col;
            if (column==this_col)
                return "Relative duplication rate (λ/μ)";
        	++this_col;
            if (column==this_col)
                return "Gain rate (γ if λ=0, κ if 0<λ)";
        	++this_col;
            if (column==this_col)
                return "Probability that a copy has no descendants at the leaves within the subtree";
        	++this_col;
            if (column==this_col)
                return "Loss parameter: probability that a copy from the parent is lost here";
        	++this_col;
            if (column==this_col)
                return "Duplication parameter for the geometric distribution of descendants from 1 parental copy (inparalogs)";
        	++this_col;
            if (column==this_col)
                return "Linear gain parameter for copies gained from outside (xenologs)";
        	++this_col;
            if (column==this_col)
                return "Survival loss parameter for conserved (non-extinct) copies";
        	++this_col;
            if (column==this_col)
                return "Survival duplication parameter for conserved (non-extinct) copies";
        	++this_col;
            if (column==this_col)
                return "Survival gain parameter for conserved (non-extinct copies)";
        	
            throw new IllegalArgumentException(".getColumnDescription --- column index must be 0.."+this_col+" [got "+column+"]");
       }
	}

}
