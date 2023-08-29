package count.gui.kit;
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


import java.awt.Dimension;
import java.awt.Font;
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JViewport;
import javax.swing.ListSelectionModel;
import javax.swing.RowSorter;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.RowSorterEvent;
import javax.swing.event.RowSorterListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;


/**
 * A class for a scroll pane that shows a table with 
 * frozen columns. The JScrollPane will have two JTables: a data
 * table that can be scrolled and a header table that contains the frozen columns. 
 * The two JTables have the same data model (TableModel), selection model
 * (ListeSelectionModel), and sorter (RowSorter). Scrolling, selection, and sorting  
 * is synchronized between the two tables. 
 * 
 * TODO: refine the tool tip texts {@link #getCellToolTip(int, int)}; for table panels 
 * it is redefined to use the model's getCellToolTip.  
 * 
 * @param <TMODEL> underlying table model
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 */
public class TableScroll<TMODEL extends AbstractTableModel> 
	extends JScrollPane 
	implements ScrollPaneConstants, TableModelListener
{
    /**
     * Instantiation. 
     * 
     * @param model Common model for row headers and data table.
     * @param num_header_columns number of columns frozen
     */
    public TableScroll(TMODEL model, int num_header_columns)
    {
        super();
        this.num_header_columns = num_header_columns;
        setModel(model);
        initComponents();
    }
    
    public static final int PREFERRED_COLUMN_WIDTH = 50;
    public static final int MINIMUM_COLUMN_WIDTH = 40;
    
    protected JTable data_table;
    protected JTable header_table;
    protected int num_header_columns;
    
    protected TMODEL model;
    
    
    public TMODEL getModel(){ return model;}
    
    protected void setModel(TMODEL model)
    {
        this.model = model;
        model.addTableModelListener(this);
    }
    
    /**
     * The underlying data table with the scrollable columns.
     * 
     * @return
     */
    public JTable getDataTable()
    {
        return data_table;
    }
    
    /**
     * The underlying header table with the frozen columns.
     * 
     * @return
     */
    public JTable getHeaderTable()
    {
        return header_table;
    }
    
    
    
    /**
     * Preferred width for a given column: called at
     * instantiation when the components are set up.
     * 
     * @param model_column_idx 0-based index of the column in the model
     * @return {@link #PREFERRED_COLUMN_WIDTH} by default
     */
    protected int getPreferredColumnWidth(int model_column_idx)
    {
        return PREFERRED_COLUMN_WIDTH;
    }

    /**
     * Minimum width for a given column: called at
     * instantiation when the components are set up.
     * 
     * @param model_column_idx 0-based index of the column in the model
     * @return {@link #MINIMUM_COLUMN_WIDTH} by default
     */
    protected int getMinimumColumnWidth(int model_column_idx)
    {
        return MINIMUM_COLUMN_WIDTH;
    }
    
    private TableColumnModel createTableDataColumns()
    {
        TableColumnModel table_data_columns = new DefaultTableColumnModel()
            {
                private int column_idx = 0;

                // modified to skip over the first few columns of the model,
                // which will be the row headers
                @Override
                public void addColumn(TableColumn col)
                {
                    if (column_idx<num_header_columns)
                    {
                        // nothing: do not add
                    } else
                    {
                        col.setMinWidth(getMinimumColumnWidth(column_idx));
                        col.setPreferredWidth(getPreferredColumnWidth(column_idx));
                        super.addColumn(col);
                    } 
                    column_idx++;
                }            
            };
        return table_data_columns;    
    }
    
    private TableColumnModel createTableHeaderColumns()
    {
        TableColumnModel table_header_columns = new DefaultTableColumnModel()
            {
                private int column_idx = 0;

                // only the first few model columns are actually added here
                @Override
                public void addColumn(TableColumn col)
                {
                    if (column_idx < num_header_columns)
                    {
                        col.setPreferredWidth(getPreferredColumnWidth(column_idx));
                        col.setMaxWidth(col.getPreferredWidth());
                        super.addColumn(col);
                    } else
                    {
                        // ignore others
                    }
                    column_idx++;
                }
            };        
        return table_header_columns;    
    }
    
    private void initComponents()
    {
        TableColumnModel table_data_columns = createTableDataColumns();
        TableColumnModel table_header_columns = createTableHeaderColumns();
        // set up the tables
        data_table =createTable(table_data_columns);
        data_table.setSelectionModel(new javax.swing.DefaultListSelectionModel());
        
        data_table.setSelectionMode(javax.swing.ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        data_table.setColumnSelectionAllowed(false);
        data_table.setCellSelectionEnabled(false);
        data_table.setRowSelectionAllowed(true);
        data_table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        data_table.createDefaultColumnsFromModel();
        
        // row header data_table
        header_table = createTable(table_header_columns);
        header_table.createDefaultColumnsFromModel();
        header_table.setColumnSelectionAllowed(false);
        header_table.setCellSelectionEnabled(false);
        header_table.setSelectionMode(javax.swing.ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        header_table.setRowSelectionAllowed(true);
        header_table.setSelectionModel(data_table.getSelectionModel());
        
        // add to JSpitPane
        // set up the data_table scrolling 
        setViewportView(data_table);
        
//        // scroll bars are always shown (MacOS X philosophy)
//        setVerticalScrollBarPolicy(VERTICAL_SCROLLBAR_ALWAYS);
//        setHorizontalScrollBarPolicy(HORIZONTAL_SCROLLBAR_ALWAYS);
        
        setHeaderView();

        // selection with double click
        data_table.addMouseListener(new MouseAdapter()
        {
            @Override
            public void mouseClicked(MouseEvent e)
            {
                if (e.getClickCount()==2 && SwingUtilities.isLeftMouseButton(e))
                {
                    Point click_xy = e.getPoint();
                    int displayed_row = data_table.rowAtPoint(click_xy);
                    int displayed_column = data_table.columnAtPoint(click_xy);
                    doubleClickInDataTable(displayed_row, displayed_column);
                }
            }
        });        
        
        data_table.setDefaultRenderer(RoundedDouble.class, new RoundedDouble.Renderer());
        
    } // initComponents()
    
    
    private JTable createTable(final TableColumnModel column_model)
    {
        JTable table = new JTable(model, column_model)
        {
        	
            @Override
            protected JTableHeader createDefaultTableHeader() 
            {
                return new JTableHeader(column_model) 
                {
                    @Override
                    public String getToolTipText(MouseEvent e) 
                    {
                        Point p = e.getPoint();
                        //
                        // TableColumnModel.getColumnIndexAtX is computed from column widths, not reliable

                        //int index = column_model.getColumnIndexAtX(p.x);
                        //int realIndex = column_model.getColumn(index).getModelIndex();
                        //System.out.println("#*FCT.cT#gTT idx "+index+"\treal "+realIndex+"\tcol "+;//+"\tevn "+e);
                        int col_idx = convertColumnIndexToModel(columnAtPoint(p));
                        return getHeaderToolTip(col_idx);
                    }
                };
            }
            
            @Override
            public String getToolTipText(MouseEvent event)            
            {
                Point p = event.getPoint();
                int row = convertRowIndexToModel(rowAtPoint(p));
                int col = convertColumnIndexToModel(columnAtPoint(p));   
                
                
                String origi = super.getToolTipText(event);
                String valinfo;
                String cellinfo = TableScroll.this.getCellToolTip(row, col);
                if (origi==null)
                {
                	Object v = getModel().getValueAt(row, col);
                	valinfo = "(="+v.toString()+")";
//                	Class<?> vclass = v==null?null:v.getClass();
//                	valinfo = "// no tooltip for "+v+"\tclass "+vclass;
                } else
                	valinfo = "(="+origi+")";
                	
                return cellinfo+valinfo;
            }
            @Override
            public TableCellRenderer getCellRenderer(int row, int column)
            {
            	int model_row = convertRowIndexToModel(row);
            	int model_col = convertColumnIndexToModel(column);
            	TableCellRenderer renderer = TableScroll.this.getCellRenderer(model_row, model_col);
            	if (renderer==null)
            		renderer = super.getCellRenderer(row, column);
            	return renderer;
            }
            
        };
        
        
       
        
        return table;
    }
    
    /**
     * Queried for displaying tool tips at table cells.
     * @param row row index by model
     * @param col column index by model
     * @return cell value (by default)
     */
    protected String getCellToolTip(int row, int col)
    {
        Object val = model.getValueAt(row, col);
        
        String cell_ident = model.getColumnName(col)+" at "+getRowName(row);
        return cell_ident;
    }
    
    public void resetDataColumns()
    {
        data_table.setColumnModel(createTableDataColumns());
        data_table.createDefaultColumnsFromModel();
    }
    
    @Override
    public void tableChanged(TableModelEvent e)
    {
    	if (e.getFirstRow() == TableModelEvent.HEADER_ROW)
    	{ // structure change
//    		System.out.println("#**TS.tC "+e+"\tstructure "+e.getType()
//    			+"\t"+e.getFirstRow()+".."+e.getLastRow()
//    			+"\t"+e.getColumn());
    		resetDataColumns();
    	} else
    	{
//    		System.out.println("#**TS.tC "+e+"\tignored typ "+e.getType());
    	}
    }
    
    /**
     * Called to get the tool tip text for column headers.
     * 
     * @param column_idx column index by the model
     * @return appropriate String for column header
     */
    protected String getHeaderToolTip(int column_idx)
    {
        //System.out.println("#*FCT.gHTT "+column_idx);
        if (column_idx>0 && column_idx<model.getColumnCount())
            return model.getColumnName(column_idx);
        else 
            return null;
    }
    
    
    /**
     * Creates a row header view from a JTable for the row headers. The 
     * view has synchronized scrolling for row selection with the input 
     * scroll pane.
     * 
     * @return a viewport for row headers
     */
    private JViewport setHeaderView()
    {
        JViewport row_header_view = new JViewport();
        row_header_view.setView(header_table);
        row_header_view.setPreferredSize(header_table.getMaximumSize());
        
        HeaderScrollSynchronizer sync = new HeaderScrollSynchronizer(row_header_view);
        row_header_view.addChangeListener(sync);

        setRowHeader(row_header_view);
        setCorner(UPPER_LEFT_CORNER, header_table.getTableHeader());

        return row_header_view;
    }    
    
    /**
     * Sets a common row sorter for the data and header tables: sorting  
     * will be synchronized between the two.
     * 
     * @param sorter
     */
    public void setRowSorter(RowSorter<? extends TMODEL> sorter)
    {
        SortSelectionSynchronizer sync = new SortSelectionSynchronizer();
        sorter.addRowSorterListener(sync); // synchronizer is to be notified after the JTables
        data_table.getSelectionModel().addListSelectionListener(sync); 
        
        data_table.setRowSorter(sorter);
        header_table.setRowSorter(sorter);
    }
    
    public ListSelectionModel getSelectionModel()
    {
    	assert (data_table.getSelectionModel()==header_table.getSelectionModel());
    	return data_table.getSelectionModel();
    }
    
    public void setSelectionModel(ListSelectionModel selection_model)
    {
    	data_table.setSelectionModel(selection_model);
    	header_table.setSelectionModel(selection_model);
    }
    
    public void synchronizeModelSelection(ListSelectionModel selection_model)
    {
    	ListSelectionModel our_model = getSelectionModel();
    	our_model.setSelectionMode(selection_model.getSelectionMode());
    	
    	our_model.addListSelectionListener(chg->
    		{
    			if (!chg.getValueIsAdjusting())
    			{
    				int[] selected_rows = getSelectedModelRows(); // sorted 
    				if (!Arrays.equals(selected_rows, selection_model.getSelectedIndices()))
    				{
	    				selection_model.setValueIsAdjusting(true);
	    				selection_model.clearSelection();
	    		        for (int j: selected_rows)
	    		        {
	    		            selection_model.addSelectionInterval(j, j);
	    		        }
	    		        selection_model.setValueIsAdjusting(false);
    				}
    			}
    		});
    	selection_model.addListSelectionListener(chg->
    		{
    			if (!chg.getValueIsAdjusting())
    			{
    				int[] selected_items = selection_model.getSelectedIndices();
    				if (!Arrays.equals(selected_items, getSelectedModelRows()))
    					setSelectedModelRows(selected_items);
    			}
    		});
    }
    
    
    /**
     * Row selection 
     * @return array of selected rows: row indices are by the model
     */
    public int[] getSelectedModelRows()
    {
        int[] selected_rows = data_table.getSelectedRows();
        int[] retval = new int[selected_rows.length];
        for (int i=0; i<selected_rows.length; i++)
        {
            retval[i] = data_table.convertRowIndexToModel(selected_rows[i]);
        }
        Arrays.sort(retval);
        return retval;
    }

    /**
     * Row selection
     * @return first selected row or -1 if no selection: row indices are by the model
     */
    public int getSelectedModelRow()
    {
        int selected_row_idx = data_table.getSelectedRow();
        if (selected_row_idx>=0)
            selected_row_idx = data_table.convertRowIndexToModel(selected_row_idx);
        return selected_row_idx;
    }

    /**
     * Row selection
     * 
     * @param selected_rows new set of rows selected (indexed by the model)
     */
    public void setSelectedModelRows(int[] selected_rows)
    {
        data_table.getSelectionModel().setValueIsAdjusting(true);
        data_table.clearSelection();
        for (int j=0; j<selected_rows.length; j++)
        {
            int row_idx = data_table.convertRowIndexToView(selected_rows[j]);
            //System.out.println("#**SSS.sSMR ["+j+"]\t"+selected_rows[j]+"\t-> "+row_idx);
            data_table.addRowSelectionInterval(row_idx, row_idx);
        }
        data_table.getSelectionModel().setValueIsAdjusting(false);

        // make sure at least one selected row visible
        if (selected_rows.length != 0)
            data_table.scrollRectToVisible(data_table.getCellRect(data_table.convertRowIndexToView(selected_rows[0]), data_table.convertColumnIndexToView(0), true));
    }
    
    public void setSelectedModelRow(int row)
    {
        data_table.getSelectionModel().setValueIsAdjusting(true);
        data_table.clearSelection();
    	data_table.addRowSelectionInterval(row, row);
        data_table.getSelectionModel().setValueIsAdjusting(false);
        data_table.scrollRectToVisible(data_table.getCellRect(data_table.convertRowIndexToView(row), data_table.convertColumnIndexToView(0), true));
    	
    }
    
    public boolean isSelectedRow(int row)
    {
    	int table_row = data_table.convertRowIndexToView(row);
    	return data_table.getSelectionModel().isSelectedIndex(table_row);
    }
    
    /**
     * Called (within the events thread) on a double-click with the left mouse button in the data table.
     * In default implementation, this calls {@link #selectByReference(int, int) }.
     * 
     * @param displayed_row_idx row index in the data table
     * @param displayed_column_idx column index in the data table
     */
    protected void doubleClickInDataTable(int displayed_row_idx, int displayed_column_idx)
    {
        int row_idx = data_table.convertRowIndexToModel(displayed_row_idx);
        int col = data_table.convertColumnIndexToModel(displayed_column_idx);
        if (Comparable.class.isAssignableFrom(model.getColumnClass(col)))
            selectByReference(row_idx,col);
        else
        {
            java.awt.Toolkit.getDefaultToolkit().beep();
            //System.out.println("#*FCT.dCIDT col "+col+"\t"+model.getColumnClass(col));
        }
    }
    
    /**
     * Called from within selectByReference() to identify the reference row in the GUI   
     * 
     * @param row_idx reference row
     * @return by default, just a String representation of row_idx
     */
    public String getRowName(int row_idx)
    {
        String retval = Integer.toString(row_idx);
        return retval;
    }

    /**
     * This method is called when the user double-clicks on 
     * a family table cell within a non-frozen column. 
     * Default implementation brings up a dialog and calls selectSimilarFamilies.
     */
    protected void selectByReference(int family_idx, int column_idx)
    {
        String column_name = model.getColumnName(column_idx);
        Object cell_value = model.getValueAt(family_idx, column_idx);
        Class column_class = model.getColumnClass(column_idx);

        if (cell_value instanceof Number) // including DoubleRoundedForDisplay
        {
            Number nvalue = (Number) cell_value;
            if (nvalue.intValue() != nvalue.doubleValue())
                cell_value = Double.valueOf(nvalue.doubleValue());
        }
        
        String[] choices = null;
        int default_choice_idx = 0;
        if (column_class == String.class)
        {
            choices = new String[2];
            choices[0] = "equals `"+cell_value+"'";
            choices[1] = "contains `"+cell_value+"'";
            default_choice_idx = 0;
        } else 
        {
            choices = new String[3];
            choices[0] = "is less than or equal to "+cell_value;
            choices[1] = "equals "+cell_value;
            choices[2] = "is greater than or equal to "+cell_value;
            default_choice_idx = 1;
        }
        
        String title = "Multiple selection using "+getRowName(family_idx)+" as reference";
        String description = "Select multiple families where entry in column `"+column_name+"'";

        String chosen_option = (String)JOptionPane.showInputDialog(this, description, title, JOptionPane.QUESTION_MESSAGE, new TableIcon(128,true), choices, choices[default_choice_idx]);
        
        if (column_class == String.class)
        {
            if (choices[0].equals(chosen_option))
                chosen_option = "eq";
            else if (choices[1].equals(chosen_option))
                chosen_option = "contains";
            else
                chosen_option = null;
        } else
        {
            if (choices[0].equals(chosen_option))
                chosen_option = "le";
            else if (choices[1].equals(chosen_option))
                chosen_option = "eq";
            else if (choices[2].equals(chosen_option))
                chosen_option = "ge";
            else
                chosen_option = null;
        }
        
        selectSimilarFamilies(family_idx, column_idx, chosen_option);
    }
    
    private String last_similar_families_command=null;
    
    /**
     * This method is called after the user double-clicks on 
     * a family table cell within a non-frozen column. 
     * The default behavior is to select 
     * all families that have an equal value in 
     * the specified column with the specified family. 
     * The values are obtained from the family table's model. 
     * 
     * @param family_idx the reference family 
     * @param col column in which the click was done (indexing is by table model)
     * @param command type of selection: "eq" for equality, "le" for less-than-or-equal, "ge" for more-than-or-equal, etc.
     * @return how many families are selected (0 if command is null)
     */
    public int selectSimilarFamilies(int family_idx, int col, String command)
    {
        if (command == null)
            return 0;
        
        data_table.getSelectionModel().setValueIsAdjusting(true);
        data_table.clearSelection();
        
        Comparable ref = (Comparable)model.getValueAt(family_idx, col); // warning but we compare values in the same column

        int num_selected = 0;
        
        for (int j=0; j<data_table.getRowCount(); j++)
        {
            int comp = ref.compareTo(model.getValueAt(j, col)); // this triggers a warning, but we are sure that values in the same column are comparable to each other
            boolean want_to_select = (comp == 0);
            want_to_select = want_to_select || (command.equals("le") && comp >0);
            want_to_select = want_to_select || (command.equals("ge") && comp <0);
            if (want_to_select)
            {
                int display_row = data_table.convertRowIndexToView(j);
                data_table.addRowSelectionInterval(display_row, display_row);
                num_selected++;
            }
        }
    	last_similar_families_command = model.getColumnName(col);
        if (command.equals("eq"))
        	last_similar_families_command += "=";
        else if (command.equals("le"))
        	last_similar_families_command += "\u2264";
        else if (command.equals("ge"))
        	last_similar_families_command += "\u2265";
        last_similar_families_command += model.getValueAt(family_idx, col);
    			
        data_table.getSelectionModel().setValueIsAdjusting(false);
        
        return num_selected;
    }
    
    /**
     * A JLabel with selection information 
     * that is updated when selection changes in this 
     * table.
     * 
     * @return
     */
    public JLabel createRowSelectionInfo()
    {
    	final Dimension label_dim = new Dimension(480,30);
    	final int max_listed = 8;
    	
    	class SelectionInfo extends JLabel implements ListSelectionListener
    	{
    		SelectionInfo()
    		{
    			super("...");
    	        setOpaque(false);
    	        setLabelFor(TableScroll.this);
    	        setFont(data_table.getFont().deriveFont(Font.ITALIC));
    	        setMaximumSize(label_dim); 
    	        setMinimumSize(label_dim);
    	        setPreferredSize(label_dim);
    	        getSelectionModel().addListSelectionListener(this);
    		}

			@Override
			public void valueChanged(ListSelectionEvent e) 
			{
				if (!e.getValueIsAdjusting())
				{
	                int num_selected = data_table.getSelectedRowCount();
	                StringBuilder selection_info_text = new StringBuilder();
	                if (num_selected == 0)
	                    selection_info_text.append("No rows selected.");
	                else
	                {
	                    if (num_selected==1)
	                        selection_info_text.append("One row selected");
	                    else
	                        selection_info_text.append(Integer.toString(num_selected)+
	                                (num_selected==model.getRowCount()?" (all)":""))
	                                .append(" rows selected");
						if (last_similar_families_command != null)
						{
							selection_info_text .append(" with ")
								.append(last_similar_families_command);
						}
						if (num_selected<=max_listed)
						{
							selection_info_text.append(" (");
							int[] selected_rows =  getSelectedModelRows();
							for (int i=0; i<selected_rows.length; i++)
							{
								if (i>0) selection_info_text.append(", ");
								int row = selected_rows[i];
								selection_info_text.append(getRowName(row));
							}
							selection_info_text.append(")");
						}
	                }
	                last_similar_families_command = null;		
					setText(selection_info_text.toString());
				}
			}
    	}
    	return new SelectionInfo();
    }
    
    /**
     * @param row model row
     * @param column model column
     * @return
     */
    private TableCellRenderer getCellRenderer(int row, int column)
    {
		return column_renderers.get(column);
    }
    
    private Map<Integer, TableCellRenderer> column_renderers=new HashMap<>();
    
    public void setColumnRenderer(int column, TableCellRenderer renderer)
    {
    	column_renderers.put(column, renderer);
    }
    public void removeColumnRenderer(int column)
    {
    	column_renderers.remove(column);
    }
    
    
    
    /**
     * A ChangeListener that can be attached to a row header view in a scrollable JTable. 
     * One or more columns of a table can be designated as a row view - this guy is a standard
     * listener for the row header view that ensures that
     * if multiple rows are selected in the row header,
     * then the scrollable part of the table still scrolls properly.
     * 
     * @author csuros
     */
    private class HeaderScrollSynchronizer implements ChangeListener
    {
        /**
         * Instantiates a new scroller. This scroller can then be added as a ChangeListener 
         * to the row_header_view.
         * 
         * @param row_header_view JViewport (for a JTable presumably) that gives the row headers
         */
        private HeaderScrollSynchronizer(JViewport row_header_view)
        {
            this.row_header_view = row_header_view;
        }

        private JViewport row_header_view;

        /**
         * Resets the viewport's y coordinate to match the header view's y coordinate
         * @param ignored 
         */
        @Override
        public void stateChanged(ChangeEvent ignored) 
        {
            Point Ph = row_header_view.getViewPosition();
            Point Pt = getViewport().getViewPosition();
            Point Pt_new = new Point(Pt.x,Ph.y);
            getViewport().setViewPosition(Pt_new);
        }
    }
    
    /**
     * This class recovers the selection after 
     * the displayed row order changes in a 
     * JTable. This RowSorterListener should be added 
     * to a RowSorter before it is set for a table. 
     * 
     * When two JTables use the same sorter and selection model, 
     * they both attempt to restore the selection 
     * after a row sorting events: the end result after two "corrections"
     * may be wrong (selection does change). This is the situation 
     * when frozen row headers are in one JTable and the other 
     * scrollable JTable displays the scrollable columns. 
     * 
     * The solution assumes that at sorting a 
     * RowSorterEvent with <code>Type.SORT_ORDER_CHANGED</code>
     * is generated, which is followed by another RowSorterEvent 
     * of <code>Type.SORTED</code>.
     * 
     * @author csuros
     */
    public class SortSelectionSynchronizer implements RowSorterListener, ListSelectionListener
    {
        public SortSelectionSynchronizer()
        {
            this.last_row_selection = new int[0];
            being_sorted = false;
        }

        private boolean being_sorted;

        private int[] last_row_selection;

        @Override
        public void sorterChanged(RowSorterEvent e)
        {
            if (e.getType() == RowSorterEvent.Type.SORT_ORDER_CHANGED)
                being_sorted = true;
            else
            {
                being_sorted = false;
                setSelectedModelRows(last_row_selection);
            }
        }

        @Override
        public void valueChanged(ListSelectionEvent e)
        {
            //System.out.println("#**SSS.vC select\t"+e+"\tsrc "+e.getSource());
            if (!being_sorted && !e.getValueIsAdjusting())
                last_row_selection = getSelectedModelRows();
        }

    }    
    

//    /**
//     * Random icon with row selection.
//     * 
//     * @author csuros
//     *
//     */
//    private class RowSelectionIcon implements javax.swing.Icon
//    {
//        @Override
//        public int getIconHeight()
//        {
//            return 128;
//        }
//
//        @Override
//        public int getIconWidth()
//        {
//            return 128;
//        }
//        
//        @Override
//        public void paintIcon(Component C, Graphics g, int x, int y)
//        {
//            int w = 128;
//            int h = 128;
//            int dw = 16;
//            int dh = 8;
//            Graphics gg = g.create();
//            gg.setFont(new Font("Serif",Font.PLAIN,6));
//            Random RND = new Random();
//            gg.translate(x,y);
//            gg.setColor(Color.WHITE);
//            gg.fillRect(0,0,w,h);
//            gg.setColor(Color.GRAY);
//            gg.drawRect(0,0,w,h);
//            int num_selected = 0;
//            for (int i=0; (i+1)*dh<=h; i++)
//            {
//                int row_y = dh*i;
//                boolean is_selected = (num_selected<i-1) && ((num_selected< i/6) || (RND.nextDouble()<0.333));
//                if (is_selected)
//                {
//                    num_selected++;
//                    gg.setColor(Color.BLUE);
//                    gg.fillRect(0, row_y, w, dh);
//                    
//                }
//                gg.setColor(Color.GRAY);
//                if (i!=0) 
//                    gg.drawLine(0, row_y, w, row_y);
//                for (int j=1; j*dw<w; j++)
//                    gg.drawLine(j*dw, row_y, j*dw, row_y+dh);
//            }
//            
//        }
//    }

}
