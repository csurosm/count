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

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import count.ds.IndexedTree;
import count.gui.kit.DrawString;
import count.gui.kit.IndexedPoint;
import count.gui.kit.RoundedDouble;
import count.gui.kit.TableScroll;

/**
 * 
 * 
 * 
 * Tree panel used in {@link HistoryView}
 * 
 * 
 *
 */
class HistoryTreePanel extends TreePanel
{
	protected static final int MAX_INDIVIDUAL_PLOT = 8;
	
	protected HistoryTreePanel(TableScroll<HistoryModel> table_scroll)
	{
		super(table_scroll.getModel().getTree(), 
				TreePanel.LayoutStyle.NODE_TABLE, 
				ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		
		this.table_scroll = table_scroll;
		this.node_totals = new ArrayList<>();
		this.node_charts = new ArrayList<>();
		
        setRootStem(0); 
        setBoundingBoxSeparation(2);
        
        Insets padding = getPadding();
        padding.left = 2*getNormalFontSize();
        
        initComponents();
//        table_scroll.getDataTable().getSelectionModel().addListSelectionListener(this);
	}
	
	private final TableScroll<HistoryModel> table_scroll;
	private JCheckBox same_scale;
	private SelectionTotals selection_totals;
	private TableScroll<SelectionTotals> totals_scroll;
	
	private final List<NodeData> node_totals;
	private final List<Chart> node_charts;
	
	
	private int chart_width = 0;
	private int chart_offset = 0;
    private int min_chart_width = 80;
	
    /**
     * Convenience method to access the tree through {@link #table_scroll}.
     * @return
     */
	private IndexedTree getTree() { return table_scroll.getModel().getTree();}

	private void initComponents()
	{
		same_scale = new JCheckBox("Same scale");
		same_scale.setSelected(true);
		same_scale.addChangeListener(e->
				{
					setValidBoundingBoxes(false);
					repaint();
				});

		selection_totals = new SelectionTotals();
	}
	
	public TableScroll<? extends TableModel> asTableScroll()
	{
		if (totals_scroll == null)
		{
			totals_scroll = new TableScroll<>(selection_totals, 2)
					{
			            @Override
			            protected int getPreferredColumnWidth(int idx)
			            {
			                if (idx==0) // order
			                	return 45;
			                else if (idx == 1 ) // node name
			                    return 180;
			                else
			                    return 60;
			            }
				
					};
			// row selection in totals_table and node selection in tree display is synchronized
					
		    TableRowSorter<SelectionTotals> sorter = new TableRowSorter<>(selection_totals);
			totals_scroll.setRowSorter(sorter);
			totals_scroll.synchronizeModelSelection(this.getSelectionModel());
			
		}
		return totals_scroll;
	}
	
	public class NodeData extends ArrayList<HistoryModel.Column<?>>
	{
		private NodeData(String header)
		{
			super(getTree().getNumNodes());
			this.header=header;
		}
		
		private NodeData(String header, List<HistoryModel.Column<?>> lineage_columns)
		{
			super(lineage_columns);
			this.header = header;
		}
		private double[] cachedTotals=null;
		private int[] last_selected_rows=null;
		private final String header;
		
		private double[] offset_values=null;

		public double getSelectionTotal(int node)
		{
			int[] selected_rows = table_scroll.getSelectedModelRows();
//			System.out.println("#**HTP.ND "+header+"\t"+node+"\tsel "+Arrays.toString(selected_rows));
			if (!Arrays.equals(selected_rows, last_selected_rows))
			{
				cachedTotals = new double[getTree().getNumNodes()];
				for (int row:selected_rows)
				{
					for (int i=0; i<cachedTotals.length; i++)
					{
						HistoryModel.Column<?> col = get(i);
						cachedTotals[i] += col.getValue(row).doubleValue();
					}
				}
//				System.out.println("#**HTP.ND "+header+"\t"+node+"\t"+Arrays.toString(cachedTotals));
				if (offset_values!=null && selected_rows.length>0 )
				{
					for (int i=0; i<cachedTotals.length; i++)
					{
						cachedTotals[i] += selected_rows.length* offset_values[i];
					}
				}
//				System.out.println("#**HTP.ND "+header+"\t"+node+"\t"+Arrays.toString(cachedTotals));
			}
			last_selected_rows = selected_rows;
			return cachedTotals[node];
		}
		
		public double getSelectionMax()
		{
			double max = 0.0;
			int num_nodes = getTree().getNumNodes();
			for (int node=0; node<num_nodes; node++)
				max = Double.max(max, getSelectionTotal(node));
			return max;
		}
		
		public void setCorrection(double[] offset)
		{
			this.offset_values = offset;
		}
		
		void clearCache()
		{
			this.last_selected_rows = null;
		}
	}
	
	public NodeData newColumns(String header, List<HistoryModel.Column<?>> lineage_columns)
	{
		NodeData col = new NodeData(header, lineage_columns);
		return col;
	}
	
	private JCheckBox addChart(Chart C)
	{
		node_charts.add(C);
		C.is_wanted.addChangeListener(selection_totals);
		setValidBoundingBoxes(false);
		C.is_wanted.setSelected(true);
		return C.is_wanted;
	}
	
	public JCheckBox showNodeStatistics(String name, NodeData count)
	{
		node_totals.add(count);
		return addChart(new Chart(name, count));
	}
	
	public JCheckBox showNodeStatistics(String name, NodeData primary, NodeData secondary)
	{
		node_totals.add(primary);
		node_totals.add(secondary);
		return addChart(new ChartSubset(name, primary, secondary));
	}
	
	public JCheckBox showChangeStatistics(String name, NodeData increase, NodeData decrease)
	{
		node_totals.add(increase);
		node_totals.add(decrease);
		return addChart(new ChartChange(name, increase, decrease));
	}
	
	public List<JCheckBox> getChartControls()
	{
		List<JCheckBox> checks = new ArrayList<>();
		for (Chart C: node_charts)
			checks.add(C.is_wanted);
		
		return checks;
	}
	
	public JCheckBox getScalingControl()
	{
		return same_scale;
	}
	
	private class SelectionTotals extends AbstractTableModel 
		implements ListSelectionListener // ListSelectionEvents (in main table) fire  table data changes
					, ChangeListener     // ChangeEvents (from checkboxes) fire table structure changes
					, TableModelListener // TableModelEvents (in main table) fire table data changes 
	{
		SelectionTotals()
		{
			super();
			// we listen to the main table
			table_scroll.getDataTable().getSelectionModel().addListSelectionListener(this);
			table_scroll.getModel().addTableModelListener(this);
		}
		
		@Override
		public int getRowCount() 
		{
			return getTree().getNumNodes();
		}

		@Override
		public int getColumnCount() 
		{
			return  1
					+1 // node name 
					+node_totals.size();
			
		}

		@Override
		public Object getValueAt(int node, int columnIndex)
		{
			if (columnIndex==0)
			{
				return (node);
			} 
			--columnIndex;
			if (columnIndex==0)
			{
				
				// node name
				return getTree().getIdent(node);
			} 
			--columnIndex;

			double tot = node_totals.get(columnIndex).getSelectionTotal(node);
			return new RoundedDouble(tot);
		}
		
        /**
         * Specify column classes to use the correct comparator in sorting
         * 
         * @param column_idx index of the column
         * @return what class the colum belongs to
         */
        @Override
        public Class<?> getColumnClass(int column_idx)
        {
            if (column_idx == 0) return Integer.class;
            --column_idx;
            if (column_idx==0) return String.class; 
            --column_idx;
            return RoundedDouble.class;
        }
		
        @Override
        public String getColumnName(int column_idx)
        {
            if (column_idx==0) return "Order";
            --column_idx; 
            if (column_idx==0) return "Node";
            --column_idx;
            
        	return node_totals.get(column_idx).header;
        }
        
        /**
         * Selection changes in main table: recompute values. 
         * Attached at instantiation.
         */
		@Override
		public void valueChanged(ListSelectionEvent e) 
		{
			if (!e.getValueIsAdjusting())
			{
//				System.out.println("#**HTP.ST.vC "+e.getFirstIndex()+".."+e.getLastIndex());
				HistoryTreePanel.this.setValidBoundingBoxes(false);
				this.fireTableRowsUpdated(0,getRowCount()-1);
				HistoryTreePanel.this.repaint();
			}
		}

		/**
		 * Our checkboxes influence tree display.
		 * Attached to {@link Chart#is_wanted} by {@link HistoryTreePanel#addChart(Chart) }
		 */
		@Override
		public void stateChanged(ChangeEvent e) 
		{
//			System.out.println("#**HTP.ST.sC "+e);
			HistoryTreePanel.this.setValidBoundingBoxes(false);
			HistoryTreePanel.this.repaint();
		}

		/**
		 * Values changing in main table: recompute if selected rows.
		 */
		@Override
		public void tableChanged(TableModelEvent e) 
		{
			boolean data_changed=false;
			if (e.getType()==TableModelEvent.UPDATE)
			{
				int first_row = e.getFirstRow();
				int last_row = e.getLastRow();
				data_changed = (last_row == Integer.MAX_VALUE);
				
				for (int row=first_row; !data_changed && row<=last_row; row++)
					if (row != TableModelEvent.HEADER_ROW)
						data_changed = table_scroll.isSelectedRow(row);

				if (data_changed)
				{
					HistoryTreePanel.this.setValidBoundingBoxes(false);
					for (NodeData column: node_totals)
						column.clearCache();
					
					this.fireTableDataChanged();
					HistoryTreePanel.this.repaint();
//					System.out.println("#**HTP.ST.tC "+e+"\tupd "+first_row+".."+last_row+"\trepaint");
				}
				else 
				{
//					System.out.println("#**HTP.ST.tC "+e+"\tupd "+first_row+".."+last_row+"\tchg "+data_changed);
				}

			} else
			{
//				System.out.println("#**HTP.ST.tC "+e+"\ttype "+e.getType()
//						+e.getFirstRow()+".."+e.getLastRow()
//						+" ignore");				
			}
		}
	}
	
	
    @Override
    protected void calculateBoundingBoxes(Graphics g)
    {
        if (g==null)
        {
        	super.calculateBoundingBoxes(g);
        	return;
        }
        
        Graphics2D g2 = (Graphics2D)g.create();
        double r = getNodeRadius();
        IndexedTree tree = getTree();
        int label_font_size = this.getLabelFontSize();
        Font label_font = new Font("Serif", Font.PLAIN, label_font_size);
        g2.setFont(label_font);

        double node_label_width = 0.0;
        for (int node=0; node<tree.getNumNodes(); node++)
        {
        	DisplayedNode N = getNode(node);
            Rectangle2D R = N.setNodeBoundingBox(g2);
            node_label_width = Double.max(node_label_width, R.getMaxX());
        }
        
        chart_offset = 1+(int) node_label_width;
        chart_width = 0;
        
        int num_charts = 0;
        for (Chart C: node_charts)
        {
        	if (C.isWanted()) num_charts++;
        }
        
        if (num_charts>0)
        {
	        int width = getWidth();
	        double hgt = getTreeHeight();
	        double tree_display_width = r*hgt;
	        
	        int sep = getTreeBoundingBoxSeparation();
	        chart_width = (int)((width
	        			-tree_display_width
	        			-node_label_width	        			
	        			)/(sep+num_charts));
	        chart_width = Integer.max(chart_width, min_chart_width);

	        DisplayedNode N0 = getNode(0);
	        int yoffset = N0.getOffsetY();
	        
    		for (int node=0; node<tree.getNumNodes(); node++)
    		{
    			DisplayedNode N = getNode(node);
    			Rectangle2D R = N.getNodeBoundingBox();
    	        int chart_pos = chart_offset;
		        for (Chart C: node_charts)
		        {
		        	if (C.isWanted()) 
		        	{
		        		chart_pos += sep;
		        		Rectangle2D chartR = new Rectangle2D.Double(chart_pos, yoffset, chart_width, label_font_size);
		        		R.add(chartR);
		        		chart_pos += chart_width;
		        	}
		        }
		        N.setNodeBoundingBox(R);
    		}
        }
    }
    
    @Override
    protected void paintNodeLabels(Graphics g)
    {
    	super.paintNodeLabels(g);
    	
    	Graphics2D g2 = (Graphics2D) g.create();
    	Color old_color = g2.getColor();

    	int chart_x0 = (int)getNode(0).getX()+chart_offset;
    	int chart_x = chart_x0;
        int sep = getTreeBoundingBoxSeparation();
        Color focus_color = new Color(128,128,0,64);

    	int num_charts = 0;
    	for (Chart C:node_charts) if (C.isWanted()) ++num_charts;
    	
    	if (num_charts>0)
    	{
	    	double common_scaling = Double.POSITIVE_INFINITY;
	    	if (same_scale.isSelected())
	    	{
		        for (Chart C: node_charts)
		        {
		        	if (C.isWanted()) 
		        	{
		        		double Cscale = C.getDefaultScaleMultiplier(g2);
		        		common_scaling = Double.min(common_scaling, Cscale);
		        	}
		        }
		        if (common_scaling<=0.0)
		        	common_scaling= 1.0;
	    	}
	    	
	    	
	    	
	    	for (Chart C: node_charts)
	    	{
	    		if (C.isWanted())
	    		{
	    			chart_x += sep;
	    			double s = same_scale.isSelected()?common_scaling:C.getDefaultScaleMultiplier(g2);
	    			
	    			C.paintChart(g2, s, chart_x);
	    			chart_x += chart_width;
	    		}
	    	}
	    	
	    	g2.setColor(focus_color);
	    	int[] selected_nodes =   getSelectionModel().getSelectedIndices();
	    	int label_font_size = getLabelFontSize();
	    	int total_chart_width = num_charts*(sep+chart_width);
	    	int focus_x = (int)getNode(0).getX();
	    	
	    	for (int node: selected_nodes)
	    	{
	    		DisplayedNode N = getNode(node); 
				int Ny = (int)(N.getY()+N.getOffsetY())-label_font_size; 
				
				g2.fillRoundRect(focus_x+N.getOffsetX(), Ny-1, chart_offset+total_chart_width+3, label_font_size+3, 3, 3);
	    	}
    	}
    }
    
	private class Chart
	{
		Chart(String name, NodeData data)
		{ 
			this.primary = data;
			this.name = name;
			is_wanted = new JCheckBox(name);
			is_wanted.setSelected(false);
		}
		final NodeData primary;
		final String name;
		final JCheckBox is_wanted;
		
		boolean isWanted() { return is_wanted.isSelected();}
		
		protected double getDefaultScaleMultiplier(Graphics2D g2)
		{
			double max = primary.getSelectionMax();
			Rectangle2D labelR = DrawString.getBoundingBoxForRotatedString(g2, RoundedDouble.toString(max), 0, 0, 0.0, 0f);
			return (chart_width-labelR.getWidth())/max;
		}
		
		protected void paintChart(Graphics2D g2, double scaling, int chart_x)
		{
			IndexedTree phylo = getTree();
			Color old_color = g2.getColor();
	        Font old_font = g2.getFont();

	        int label_font_size = getLabelFontSize();
	        Font label_font = new Font("Serif", Font.PLAIN, label_font_size);
	        g2.setFont(label_font);
	        
			int maxw=0;
			int maxy=0;
			for (int node=0; node<phylo.getNumNodes(); node++)
			{
				DisplayedNode N = getNode(node);
				int Ny = (int)(N.getY()+N.getOffsetY())-label_font_size; // align with node label
				maxy = Integer.max(maxy, Ny+label_font_size);

				double value = primary.getSelectionTotal(node);
				
				if (Double.isFinite(value) && value!=0.0)
				{
					int valw = (int)(value*scaling);
					maxw = Integer.max(maxw, valw);
					g2.setColor(N.getIcon(true).getFillColor());
					g2.fillRect(chart_x, Ny, valw, label_font_size);
					g2.setColor(Color.BLACK);
					DrawString.drawLeft(g2, RoundedDouble.toString(value), chart_x+valw+1, Ny+label_font_size);
				}
			}
			// X and Y axis
			g2.setColor(Color.BLACK);
			g2.drawLine(chart_x, 2*label_font_size, chart_x, maxy);
			g2.drawLine(chart_x, 2*label_font_size, chart_x+maxw, 2*label_font_size);
			// title
			g2.setFont(label_font.deriveFont(Font.BOLD));
			g2.drawString(name, chart_x, label_font_size-1);
			g2.setColor(old_color);
			g2.setFont(old_font);
		}
	}
	
	private class ChartSubset extends Chart
	{
		final NodeData secondary;
		ChartSubset(String name, NodeData primary, NodeData secondary)
		{
			super(name, primary);
			this.secondary = secondary;
		}
		@Override
		protected void paintChart(Graphics2D g2, double scaling, int chart_x)
		{
			IndexedTree phylo = getTree();
			Color old_color = g2.getColor();
	        int label_font_size = getLabelFontSize();
	        Font label_font = new Font("Serif", Font.PLAIN, label_font_size);
	        Font old_font = g2.getFont();
	        g2.setFont(label_font);
	        
			int maxw=0;
			int maxy=0;
			for (int node=0; node<phylo.getNumNodes(); node++)
			{
				DisplayedNode N = getNode(node);
				int Ny = (int)(N.getY()+N.getOffsetY()-label_font_size); // align with node label
				maxy = Integer.max(maxy, Ny+label_font_size);

				double value = primary.getSelectionTotal(node);
				if (Double.isFinite(value) && value!=0.0)
				{
					int valw = (int)(value*scaling);
					maxw = Integer.max(maxw, valw);
					
					double val2 = secondary.getSelectionTotal(node);
					int val2w = (int)(val2*scaling);
					
					
					g2.setColor(N.getIcon(true).getFillColor());
					g2.fillRect(chart_x, Ny+label_font_size-label_font_size/2, valw, label_font_size/2);
					g2.fillRect(chart_x, Ny, val2w, label_font_size/2);
					g2.setColor(Color.BLACK);
					DrawString.drawLeft(g2, RoundedDouble.toString(value), chart_x+valw+1, Ny+label_font_size);
				}
			}
			g2.setColor(Color.BLACK);
			g2.drawLine(chart_x, 2*label_font_size, chart_x, maxy);
			g2.drawLine(chart_x, 2*label_font_size, chart_x+maxw, 2*label_font_size);
			// title
			g2.setFont(label_font.deriveFont(Font.BOLD));
			g2.drawString(name, chart_x, label_font_size-1);
			g2.setColor(old_color);
			g2.setFont(old_font);
			
		}
		
	}
	
	private class ChartChange extends Chart
	{
		final NodeData increase;
		ChartChange(String name, NodeData increase, NodeData decrease)
		{
			super(name, decrease);
			this.increase = increase;
		}
		@Override
		protected double getDefaultScaleMultiplier(Graphics2D g2)
		{
			double max_decrease = primary.getSelectionMax();
			double max_increase = increase.getSelectionMax();
			double max_range = max_increase+max_decrease;
			Rectangle2D increaseR = DrawString.getBoundingBoxForRotatedString(g2, RoundedDouble.toString(max_increase), 0, 0, 0.0, 0f);
			Rectangle2D decreaseR = DrawString.getBoundingBoxForRotatedString(g2, RoundedDouble.toString(max_decrease), 0, 0, 0.0, 0f);
			return (chart_width-increaseR.getWidth()-decreaseR.getWidth()-2.0)/max_range;
		}
		
		@Override
		protected void paintChart(Graphics2D g2, double scaling, int chart_x)
		{
			IndexedTree phylo = getTree();
			Color old_color = g2.getColor();
	        int label_font_size = getLabelFontSize();
	        Font label_font = new Font("Serif", Font.PLAIN, label_font_size);
	        Font old_font = g2.getFont();
	        g2.setFont(label_font);
	        
			// place the axis
			double max_decrease = primary.getSelectionMax();
			double max_increase = increase.getSelectionMax();
			double max_range = max_increase+max_decrease;
			Rectangle2D decreaseR = DrawString.getBoundingBoxForRotatedString(g2, RoundedDouble.toString(max_decrease), 0, 0, 0.0, 0f);
			Rectangle2D increaseR = DrawString.getBoundingBoxForRotatedString(g2, RoundedDouble.toString(max_increase), 0, 0, 0.0, 0f);
			double plot_width = chart_width-(decreaseR.getWidth()+increaseR.getWidth()+2.0);
			
			int axisX = (int)(chart_x + decreaseR.getWidth()+1.0+ plot_width*max_decrease/max_range);
			int h = label_font_size; // bar height
	        
	        
			int max_decreaseW=0;
			int max_increaseW=0;
			int maxY=0;
			for (int node=0; node<phylo.getNumNodes(); node++)
			{
				DisplayedNode N = getNode(node);
				
				double plus = increase.getSelectionTotal(node);
				double minus = primary.getSelectionTotal(node);
				double change = plus-minus;
				int Ny = (int)(N.getY()+N.getOffsetY()-label_font_size); // align with node label
				maxY = Integer.max(maxY, Ny+h);
				
				if (Double.isFinite(plus) && plus!=0.0)
				{
					int plusW = (int)(plus*scaling);
					g2.setColor(N.getIcon(true).getFillColor());
					g2.drawRect(axisX, Ny, plusW, h);
					g2.setColor(Color.BLACK);
					DrawString.drawLeft(g2, RoundedDouble.toString(plus), axisX+plusW+1, Ny+label_font_size);
					max_increaseW = Integer.max(max_increaseW, plusW);
				}
				if (Double.isFinite(minus) && minus!=0.0)
				{
					int minusW = (int)(minus*scaling);
					g2.setColor(N.getIcon(true).getFillColor());
					g2.drawRect(axisX-minusW, Ny, minusW, h);
					g2.setColor(Color.BLACK);
					DrawString.drawRight(g2, RoundedDouble.toString(minus), axisX-minusW-1, Ny+label_font_size);
					max_decreaseW = Integer.max(max_decreaseW, minusW);
				}

				if (Double.isFinite(change) && change!=0.0)
				{
					int changeW = (int)(change*scaling);
					g2.setColor(N.getIcon(true).getFillColor());
					if (change<0.0)
					{
						g2.fillRect(axisX+changeW, Ny, -changeW, h);
					} else
					{
						g2.fillRect(axisX, Ny, changeW, h);
					}
				}
			}
			g2.setColor(Color.BLACK);
			g2.drawLine(axisX, 2*label_font_size, axisX, maxY);
			g2.drawLine(axisX-max_decreaseW, 2*label_font_size, axisX+max_increaseW, 2*label_font_size);
			// title
			g2.setFont(label_font.deriveFont(Font.BOLD));
			g2.drawString(name, chart_x, label_font_size-1);
			g2.setColor(old_color);
			g2.setFont(old_font);
		}
	}
}
