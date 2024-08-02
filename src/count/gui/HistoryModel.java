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
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.stream.Stream;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.gui.kit.RoundedDouble;
import count.io.CommandLine;

/**
 * TableModel for ancestral reconstructions
 * with  
 * one row for each family: lists-of-columns 
 * for lineage-specific statistics, 
 * and arithmetically defined columns.
 * 
 *
 * 
 * @author csuros
 *
 */
class HistoryModel extends AnnotatedTableModel 
{
	protected HistoryModel(IndexedTree tree, AnnotatedTable data_table)
	{
		super(data_table, false);
		this.history_columns = new ArrayList<>();
		this.tree = tree;
	}
	
	private List<Column<?>> history_columns;
	private final IndexedTree tree;
	
	public IndexedTree getTree() { return tree;}
	
	/**
	 * Adds a column to the TableModel. 
	 * 
	 * @param col
	 * @return index of newly added column in the model
	 */
	public int add(Column<?> col, Color color)
	{
		col.setColor(color);
		history_columns.add(col);
		return history_columns.size()-1;
	}
	
	public void clear() {history_columns.clear();} 
	
	public void resetValues()
	{
		for (Column<?> col: history_columns)
			col.reset();
	}
	
	public void addAll(Collection<Column<?>> columns, Color color)
	{
		for (Column<?> C: columns) add(C, color);
	}
	
	public Column<RoundedDouble> newDoubleColumn(String header, String tooltip)
	{
		DoubleColumn col = new DoubleColumn(header, tooltip);
		return col;
	}
	
	public Column<Integer> newIntColumn(String header, String tooltip)
	{
		IntColumn col = new IntColumn(header, tooltip);
		return col;
	}
	
	public List<Column<Integer>> newIntColumns(String header, String tooltip, boolean is_lineage_statistic)
	{
		List<Column<Integer>> cols = newColumns(Integer.class, header, tooltip, is_lineage_statistic);
		return cols;
	}

	public List<Column<RoundedDouble>> newDoubleColumns(String header, String tooltip, boolean is_lineage_statictic)
	{
		List<Column<RoundedDouble>> cols = newColumns(RoundedDouble.class, header, tooltip, is_lineage_statictic);
		return cols;
	}
		
	public <T extends Number> Column<T> newSumColumn(List<? extends Column<T>> columns, String header, String tooltip)
	{
		SumColumn<T> sumcol = new SumColumn(header, tooltip, columns);
		return sumcol;
	}
	
	public <T extends Number> Column<T> newDifferenceColumn(Column<T> left, Column<T> right, String header, String tooltip)
	{
		DifferenceColumn<T> diff = new DifferenceColumn(header, tooltip, left, right, null);
		return diff;
	}
	
	public <T extends Number> Column<T> newNonnegativeDifferenceColumn(Column<T> left, Column<T> right, String header, String tooltip)
	{
		DifferenceColumn<T> diff = new DifferenceColumn(header, tooltip, left, right, 0.0);
		return diff;
	}
	
	public Column<Integer> newThresholdColumn(Column<Integer> col, int min_value, String header, String tooltip)
	{
		Column<Integer> threshold = new Column<>(header, tooltip, Integer.class)
				{
					@Override
					public Integer getValue(int row)
					{
						int x = col.getValue(row);
						return x<min_value?0:1;
					}
				};
		return threshold;
	}
	
	/**
	 * 0/1 column for a change biggerC&gt;small;smallerC==small. 
	 * 
	 * @param biggerC
	 * @param smallerC
	 * @param small
	 * @param header
	 * @param tooltip
	 * @return
	 */
	public Column<Integer> newChangeColumn(Column<Integer> biggerC, Column<Integer> smallerC, int small, String header, String tooltip)
	{
		Column<Integer> change = new Column<>(header, tooltip, Integer.class)
				{
					@Override
					public Integer getValue(int row)
					{
						int x = biggerC==null?0:biggerC.getValue(row);
						if (x<=small) return 0;
						int y = smallerC==null?0:smallerC.getValue(row);
						return y==small?1:0;
					}
				};
		return change;
	}
	
	
	protected String getToolTip(String tooltip_info, int node, boolean edges_only)
	{
		StringBuilder tsb = new StringBuilder(tooltip_info);
		if (edges_only)
		{
			// lineage statistic
			tsb.append(" in the lineage ");
			if (tree.isRoot(node))
			{
				tsb.append("leading to root");
			} else
			{
				tsb.append(tree.getIdent(tree.getParent(node)))
				.append("→")
				.append(tree.getIdent(node));
			}
		} else
		{
			// node statistic
			tsb.append(" at node ")
			.append(tree.getIdent(node));
		}
		return tsb.toString();
	}
	
	protected String getHeader(String header, int node)
	{
		String hdr;
		if (tree.isLeaf(node))
			hdr = tree.getName(node)+header;
		else
			hdr = Integer.toString(node)+header;
		
		
		return hdr;
	}
	
	private <T extends Number> 
		List<Column<T>> newColumns(
				Class<T> cell_class,
				String header, String tooltip,
				boolean edges_only)
	{
		int num_cols = tree.getNumNodes();
		assert (tree.isRoot(num_cols-1));
		List<Column<T>> cols = new ArrayList<>(num_cols);

		int node = 0;
		int max_node = tree.getNumNodes();
		while (node<max_node) // leaves-first postorder
		{
			String hdr=getHeader(header, node);
			String ttp=getToolTip(tooltip,node, edges_only);
			
			Column<T> C;
			if (cell_class == RoundedDouble.class)
				C = (Column<T>) new DoubleColumn(hdr, ttp); // warning for unchecked type conversion but logic correct
			else
			{
				assert (cell_class == Integer.class);
				C = (Column<T>) new IntColumn(hdr, ttp); // warning for unchecked type conversion but logic correct 
			}
			cols.add(C);
			node++;
		}
		return cols;
	}
	
	public abstract class Column<T extends Number>
	{
		private Column(String header, String tooltip, Class<T> cell_class)
		{
			this(header, tooltip, cell_class, null);
		}
		private Column(String header, String tooltip, Class<T> cell_class, Color color)
		{
			this.header = header;
			this.header_tooltip = tooltip;
			this.cell_class = cell_class;
			this.color = color;
		}
		protected String header;
		protected final String header_tooltip;
		protected final Class<T> cell_class;
		protected Color color;
		
		public abstract T getValue(int row_idx);
		public void setValue(int row_idx, Number value)
		{
			// no effect by default
		}
		public void reset()
		{
			// no effect by default
		}
		
		public void setCorrection(Number offset)
		{
			// no effect by default
		}
		
		public String getHeader()
		{
			return header;
		}
		
		public double getMinimum()
		{
			double min = Double.POSITIVE_INFINITY;
			for (int row=0; row<getRowCount(); row++)
				min = Double.min(min, getValue(row).doubleValue());
			return min;
		}
		public double getMaximum()
		{
			double max = Double.NEGATIVE_INFINITY;
			for (int row=0; row<getRowCount(); row++)
				max = Double.max(max, getValue(row).doubleValue());
			return max;
		}
		public void setColor(Color color)
		{
			this.color = color;
		}
		public Color getColor()
		{
			return color;
		}
	}

	private class DoubleColumn extends Column<RoundedDouble>
	{
		private DoubleColumn(String header, String tooltip)
		{
			super(header,tooltip, RoundedDouble.class);
			this.cells = new double[getRowCount()];
			reset();
		}
		private final double[] cells;
		private double offset=0.0;
		
		
		@Override
		public RoundedDouble getValue(int row_idx) 
		{
			double y = cells[row_idx]+offset;
			return new RoundedDouble(y);
		} 
		@Override
		public void setValue(int row_idx, Number x) { cells[row_idx]=x.doubleValue();}
		@Override 
		public void reset()
		{
			Arrays.fill(cells, Double.NaN);
		}
		@Override
		public void setCorrection(Number x)
		{
			this.offset = x.doubleValue();
		}
	}
	private class IntColumn extends Column<Integer>
	{
		private IntColumn(String header, String tooltip)
		{
			super(header,tooltip, Integer.class);
			this.cells = new int[getRowCount()];
			reset();
		}
		private final int[] cells;
		@Override
		public Integer getValue(int row_idx) { return cells[row_idx];} 
		@Override
		public void setValue(int row_idx, Number x) 
		{ 
			cells[row_idx]=x.intValue();
		}
		@Override
		public void reset()
		{
			Arrays.fill(cells, -1);
		}
	}
	private class SumColumn<C extends Number> extends Column<C>
	{
		private SumColumn(String header, String tooltip, List<? extends Column<C>> column_list)
		{
			super(header, tooltip, (Class<C>) column_list.get(0).getValue(0).getClass());// warning for conversion but logically correct.
			this.column_list = column_list;
		}
		private final Collection<? extends Column<C>> column_list;
		@Override
		public C getValue(int row)
		{
			C getVal;
			if (cell_class==Integer.class)
			{
				Stream<Integer> cell_stream = column_list.stream().map(col->col.getValue(row).intValue());
				Integer sum = cell_stream.reduce(0, Integer::sum);
				getVal = cell_class.cast(sum); 
			} else
			{
				assert cell_class == RoundedDouble.class;
				Stream<Double> cell_stream = column_list.stream().map(col->col.getValue(row).doubleValue());
				Double sum = cell_stream.reduce(0.0, Double::sum);
				getVal = cell_class.cast(new RoundedDouble(sum));
			}
			return getVal;
		}
	}
	
	private class DifferenceColumn<C extends Number> extends Column<C>
	{
		DifferenceColumn(String header, String tooltip, Column<C> left, Column<C> right, Double min_value)
		{
			super(header, tooltip, 
					left==null?(Class<C>) right.getValue(0).getClass():
					(Class<C>) left.getValue(0).getClass());// warning for conversion but logically correct.			
			this.left = left;
			this.right = right;
			this.min_value = min_value;
		}
		private final Column<C> left;
		private final Column<C> right;
		private final Double min_value;

		@Override
		public C getValue(int row)
		{
			C getVal;
			if (right==null)
			{
				getVal = left.getValue(row);
			} else
			{
				if (cell_class==Integer.class)
				{
					int lval = left==null?0:left.getValue(row).intValue();
					int diff = lval-right.getValue(row).intValue();
					if (min_value!=null && diff<min_value) diff = min_value.intValue();
					getVal = cell_class.cast(Integer.valueOf(diff)); 
				} else
				{
					assert cell_class == RoundedDouble.class;
					double lval = (left==null?0.0:left.getValue(row).doubleValue());
					double diff = lval-right.getValue(row).doubleValue();
					if (min_value != null && diff<min_value) diff=min_value;
					getVal = cell_class.cast(new RoundedDouble(diff));
				}
			}
			return getVal;
		}
	}
	
    @Override
    public int getColumnCount()
    {
    	return super.getColumnCount() + history_columns.size();
    }	
    
    public int firstHistoryColumn()
    {
    	return super.getColumnCount();
    }
    
    public Column<?> getHistoryColumn(int col)
    {
    	int first = firstHistoryColumn();
    	if (col<first)    	
    		return  null;
    	else
    		return history_columns.get(col-first);
    }
    
    @Override
    public Object getValueAt(int row_idx, int column_idx)
    {
    	int first = firstHistoryColumn();
    	if (column_idx<first)
    		return super.getValueAt(row_idx, column_idx);
    	else 
    	{
    		Column<?> col = history_columns.get(column_idx-first);
    		return col.getValue(row_idx);
    	}
    }
    
    @Override
    public String getColumnName(int column_idx)
    {
    	int first = firstHistoryColumn();
    	if (column_idx<first)
    		return super.getColumnName(column_idx);
    	else 
    	{
    		Column<?> col = history_columns.get(column_idx-first);
    		return col.header;
    	}
    }
    
    /**
     * Specify column classes to use the correct comparator in sorting
     * 
     * @param column_idx index of the column
     * @return what class the column belongs to
     */
    @Override
    public Class<?> getColumnClass(int column_idx)
    {
    	int first = firstHistoryColumn();
    	if (column_idx<first)
    		return super.getColumnClass(column_idx);
    	else 
    	{
    		Column<?> col = history_columns.get(column_idx-first);
    		return col.cell_class;
    	}
    }
    
    @Override
    public String getColumnDescription(int column_idx)
    {
    	int first = firstHistoryColumn();
    	if (column_idx<first)
    		return super.getColumnDescription(column_idx);
    	else 
    	{
    		Column<?> col = history_columns.get(column_idx-first);
    		return col.header_tooltip;
    	}
    	
    }

//    @Override
//    public String getCellToolTip(int row_idx, int column_idx)
//    {
//        System.out.println("#**HM.gCTT "+row_idx+","+column_idx);
//    	int first = firstHistoryColumn();
//    	if (column_idx<first)
//    		return super.getCellToolTip(row_idx, column_idx);
//    	else
//    	{
//    		return "HistoryModel.getCellToolTip";
//    		
//    	}
//    	
//    }    
    
    public String getTextTable()
    {
    	StringBuilder sb = new StringBuilder();
    	
    	int nF = this.getRowCount();
    	int nC = this.getColumnCount();
    	int fC = firstHistoryColumn();
    	
    	// column headers
    	sb.append("Family");
    	for (int c=fC; c<nC; c++)
    		sb.append("\t").append(getColumnName(c));
    		
    	AnnotatedTable tbl = this.getTable();
    	
    	for (int f=0; f<nF; f++)
    	{
        	sb.append("\n")
        		.append(tbl.getFamilyName(f));
    		for (int c=fC; c<nC; c++)
    		{
    			Object v = getValueAt(f, c);
    			int iv = ((Number)v).intValue();
    			double dv = ((Number)v).doubleValue();
    			
    			if (iv==dv)
    			{
    				sb.append("\t").append(iv);
    			} else
    			{
    				sb.append("\t").append(dv);
    			}
    		}
    	}
    	
    	return sb.toString();
    }
    
    public void printTable(PrintStream out)
    {
    	int nF = this.getRowCount();
    	int nC = this.getColumnCount();
    	int fC = firstHistoryColumn();
    	
//    	System.out.println("#**HM.pT "+nF+"\tfamilies start\tncol "+nC+"\tfcol "+fC);
    	
    	// column headers
    	out.print("Family");
    	for (int c=fC; c<nC; c++)
    		out.printf("\t%s",getColumnName(c));
    	out.println();
    	
    	AnnotatedTable tbl = this.getTable();
    	for (int f=0; f<nF; f++)
    	{
    		out.print(tbl.getFamilyName(f));  
    		for (int c=fC; c<nC; c++)
    		{
    			out.print("\t");
    			Object v = getValueAt(f, c);
    			int iv = ((Number)v).intValue();
    			double dv = ((Number)v).doubleValue();
    			
    			if (iv==dv)
    			{
    				out.print(iv);
    			} else
    			{
    				out.print(dv);
    			}
    		}    
    		out.println();
    	}
    	
    	
//    	System.out.println("#**HM.pT "+nF+"\tfamilies done");
    }
    
    
}
