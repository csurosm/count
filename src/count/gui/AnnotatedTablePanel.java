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
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.table.TableRowSorter;

import count.ds.AnnotatedTable;
import count.ds.IntegerOrMissing;
import count.gui.kit.TableScroll;
import count.io.CommandLine;
import count.io.DataFile;
import count.io.SavableData;
import count.io.TableParser;

/**
 * Swing component (JPanel) for displaying an AnnotatedTable.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 */
public class AnnotatedTablePanel 
		// extends Browser.PrimaryItem // which is a JPanel
		extends JPanel
		implements Session.FamilySelection
					, SavableData<AnnotatedTable>
		//implements ListSelectionListener
{
	
    public static int TABLE_FONT_SIZE = 14;

    private final DataFile<AnnotatedTable> data;

    /**
     *  Initializes this display with the given OccurrenceTable file.
     * 
     *  @param T the underlying data file
     */
    public AnnotatedTablePanel(DataFile<AnnotatedTable> T)
    {
        this.data = T;
        initComponents();
    }
    
    /**
     * Table for displaying the data set;
     */
    private JTable table;
    
    /**
     * Scroll pane for holding the data set
     */
    private TableScroll<AnnotatedTableModel> table_scroll;
    
    public AnnotatedTableModel getTableModel()
    {
    	return table_scroll.getModel();
    }

    /**
     * A label about what families are selected: placed at the bottom.
     */
    private JLabel selected_rows_information;

//    /**
//     * A semaphore variable for managing the displayed
//     * information about the current family selection.
//     * Families can be selected by usual table selection,
//     * or by a double-click in a cell. In the latter case,
//     * the information can be set by the selectSimilarFamilies()
//     * method. The default information is displayed by a ListSelectionListener
//     * for the family table. This variable is set to false in
//     * order to disable the update by the selection listener
//     * when the family selection is done by a computation.
//     *
//     */
//    protected boolean update_selected_rows_information = true;

    /**
     * Returns the underlying data file.
     * 
     * @return the underlying data file
     */
    @Override
    public DataFile<AnnotatedTable> getDataFile()
    {
        return data;
    }
    
    /**
     * Initialization of components within this display.
     */
    private void initComponents()
    {
        setBackground(Color.WHITE);
        boolean show_detailed_profiles = true;
        final AnnotatedTableModel M = new AnnotatedTableModel(data.getContent(), show_detailed_profiles);
        
        class FamilyScroll extends TableScroll<AnnotatedTableModel>
        {
            private FamilyScroll()
            {
                super(M,2);
            }

            @Override
            protected int getPreferredColumnWidth(int idx)
            {
                if (idx == 1 ) // family name
                    return 120;
                else
                    return super.getPreferredColumnWidth(idx);
            }

            @Override
            public String getRowName(int family_idx)
            {
                return data.getContent().getFamilyName(family_idx);
            }

            @Override
            protected String getHeaderToolTip(int column_idx)
            {
                String tt = model.getColumnDescription(column_idx);

                return tt +
                        "; click to sort rows, drag to rearrange columns, or double-click to split the table";
            }

            @Override
            protected String getCellToolTip(int row, int col)
            {
                return model.getCellToolTip(row, col);
            }

            
            @Override 
            protected void selectByReference(int family_idx, int column_idx)        
            {
                Class column_class = model.getColumnClass(column_idx);
                if (column_class == IntegerOrMissing.class)
                {
                    Object cell_value = model.getValueAt(family_idx, column_idx);
                    IntegerOrMissing iValue = (IntegerOrMissing) cell_value;
                    if (iValue.isAmbiguous())
                        return;
                }
                super.selectByReference(family_idx, column_idx);
            }
            
//            @Override
//            public int selectSimilarFamilies(int family_idx, int col, String command)
//            {
//                if (command == null)
//                    return 0;
//                else
//                {
//                    int num_selected = super.selectSimilarFamilies(family_idx, col, command);
//                    String info = "";
//                    if (num_selected == 1)
//                    {
//                        info = "One row";
//                    } else
//                    {
//                        info = Integer.toString(num_selected)+" rows";
//                    }
//                    info += " selected with "+model.getColumnName(col);
//                    if (command.equals("eq"))
//                        info += "=";
//                    else if (command.equals("le"))
//                        info += "\u2264";
//                    else if (command.equals("ge"))
//                        info += "\u2265";
//                    info += model.getValueAt(family_idx, col);
//                    displaySelectionInfo(info);
//                    update_selected_rows_information = false;
//
//                    return num_selected;
//                }
//            }
        }

        table_scroll = new FamilyScroll();
        table = table_scroll.getDataTable();
        table.setFont(new Font("Serif",Font.PLAIN,TABLE_FONT_SIZE));
        
        M.setDefaultRenderers(table);
//        int fp = M.firstProfileColumn();
//        for (int leaf=0; leaf<data.getContent().getTaxonCount(); leaf++)
//        {
//        	int max = (int)M.getMaximumValue(fp+leaf);
//        	ColoredValueRenderer renderer =
//        			new ColoredValueRenderer(Color.WHITE, FAMILY_PRESENT_COLOR, max);
//        	table_scroll.setColumnRenderer(fp+leaf, renderer);
//        }
        
        
        // row header data_table
        JTable row_header =table_scroll.getHeaderTable();
        row_header.setFont(table.getFont());

        table_scroll.getViewport().setBackground(getBackground());
        
        // sorting
        TableRowSorter<AnnotatedTableModel> sorter = new TableRowSorter<>(M);
        table_scroll.setRowSorter(sorter);

        Font tp_font_rm = table.getFont().deriveFont(0.8f);
        Font tp_font_it = tp_font_rm.deriveFont(Font.ITALIC);

        selected_rows_information = table_scroll.createRowSelectionInfo(); //  new JLabel(":");
        selected_rows_information.setFont(tp_font_it.deriveFont(TABLE_FONT_SIZE*0.8f));
//        selected_rows_information.setOpaque(false);
//        selected_rows_information.setLabelFor(table);
//        selected_rows_information.setMaximumSize(new Dimension(520,30)); // not even a constant setup for this
//        selected_rows_information.setMinimumSize(selected_rows_information.getMaximumSize());
//        selected_rows_information.setPreferredSize(selected_rows_information.getMaximumSize());
//        
//        table.getSelectionModel().addListSelectionListener(this);

        this.setLayout(new BorderLayout());
        this.add(table_scroll, BorderLayout.CENTER);
        this.add(selected_rows_information, BorderLayout.SOUTH);

        table.setRowSelectionInterval(0,0);
    }


//    private void displaySelectionInfo(String text)
//    {
//        if (selected_rows_information != null)
//            selected_rows_information.setText(text);
//    }
//
    
//    /**
//     * We are listening to selection changes in the family table.
//     * When selection changes, the selection inormation in the bottom
//     * bar is updated.
//     *
//     * @param e a list selection event: the info is updated only if the event is not adjusting
//     */
//    @Override
//    public void valueChanged(ListSelectionEvent e)
//    {
//        if (!e.getValueIsAdjusting())
//        {
//            if (update_selected_rows_information)
//            {
//
//                int num_selected = table.getSelectedRowCount();
//
//
//                String selection_info_text = "";
//                if (num_selected == 0)
//                    selection_info_text = "No rows selected.";
//                else
//                {
//                    if (num_selected==1)
//                        selection_info_text = "One row selected";
//                    else
//                        selection_info_text = Integer.toString(num_selected)+
//                                (num_selected==table.getRowCount()?" (all)":"")
//                                +" rows selected";
//                }
//
//                displaySelectionInfo(selection_info_text);
//            }
//            else
//                update_selected_rows_information = true; // next time
//        }
//    }

    
    /**
     * Gives the name of the underlying data file (will be shown in the tree)
     * 
     * @return name of the underlying file or String for "noname" if associated file is null
     */
    @Override
    public String toString()
    {
        if (data.getFile()==null)
        {
            return "[noname]";
        } else
            return DataFile.chopFileExtension(data.getFile().getName());
    }

	@Override
	public int[] getSelectedFamilies() 
	{
		return table_scroll.getSelectedModelRows();
	}
	
	/*
	 * Savable interface
	 */
	
	@Override
	public void saveData(File f) throws IOException 
	{
        PrintStream PS = new PrintStream(f);
        PS.println(CommandLine.getStandardHeader(getClass()));
        PS.println(CommandLine.getStandardRuntimeInfo());

        PS.print(TableParser.getFormattedTable(data.getContent(), true));
        
//        JRadioButton normal_layout = new JRadioButton("Standard format (can be imported later)");
//        JRadioButton malin_layout = new JRadioButton("Transposed format (every row is the concatenated string of family sizes for an organism)");
//        ButtonGroup layout_choices = new ButtonGroup();
//        layout_choices.add(normal_layout);
//        layout_choices.add(malin_layout);
//        normal_layout.setSelected(true);
//        JPanel choice_panel = new JPanel();
//        choice_panel.add(normal_layout);
//        choice_panel.add(malin_layout);
//        JOptionPane.showConfirmDialog(this, choice_panel);
//        if (normal_layout.isSelected())
//        {
//            PS.print(TableParser.getFormattedTable(data.getContent(), true));
//        } else if (malin_layout.isSelected())
//        {
//            PS.print(TableParser.getTransposedTable(data.getContent()));
//        }
        if (PS.checkError()) // also flushes
        {
            PS.close();
            throw new IOException("Cannot write the table.");
        }
        PS.close();
	}

}
