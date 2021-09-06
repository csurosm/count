package ca.umontreal.iro.evolution.malin.ui.count;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import java.util.Vector;
import java.util.HashSet;
import java.util.Hashtable;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Font;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JViewport;
import javax.swing.ScrollPaneConstants;

//import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableRowSorter;

import ca.umontreal.iro.evolution.genecontent.OccurrenceTable;

import ca.umontreal.iro.evolution.malin.ui.Browser;
import ca.umontreal.iro.evolution.malin.ui.Dealer;
import ca.umontreal.iro.evolution.malin.ui.FrozenColumnsTable;


import ca.umontreal.iro.evolution.malin.DataFile;
import ca.umontreal.iro.evolution.malin.Saveable;


/**
 * Swing component (JPanel) for displaying an OccurrenceTable.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class FamilySizeTableDisplay 
        extends Browser.PrimaryItem // which is a JPanel
        implements Saveable, FilterableFamilies, ListSelectionListener
{
    private DataFile<OccurrenceTable> data;

    /**
     *  Initializes this display with the given OccurrenceTable file.
     * 
     *  @param T the underlying data file
     */
    public FamilySizeTableDisplay(DataFile<OccurrenceTable> T)
    {
        this.data = T;
        initComponents();
    }

    /**
     * Displayed text in menu for annotations file format 
     */
    private static final String SEPARATOR_COMMA = "Comma (',')";
    private static final String SEPARATOR_TAB = "Tab ('\t')";
    private static final String SEPARATOR_SPACE = "Space (' ')";
    private static final String SEPARATOR_WHITESPACE = "Whitespace sequence ('\\s+')";
    
    private static final String COLUMN_FAMILY = "Family";
    private static final String COLUMN_SKIPPED = "(SKIPPED)";
    private static final String COLUMN_CATEGORY = "Category";
    private static final String COLUMN_DESCRIPTION = "Description";

        
    /**
     * Table for displaying the data set;
     */
    private JTable table;
    
    /**
     * Scroll pane for holding the data set
     */
    private FrozenColumnsTable<OccurrenceTableModel> table_scroll;
    
    /**
     * Fixed row headers that are not included in the scroll
     */
    private JViewport row_header_view;

    /**
     * A label about what families are selected: placed at the bottom.
     */
    private JLabel selected_rows_information;

    /**
     * A semaphore variabel for managing the displayed
     * information about the current family selection.
     * Families can be selected by usual table selection,
     * or by a double-click in a cell. In the latter case,
     * the information can be set by the selectSimilarFamilies()
     * method. The default information is displayed by a ListSelectionListener
     * for the family table. This variable is set to false in
     * order to disable the update by the selection listener
     * when the family selection is done by a computation.
     *
     */
    protected boolean update_selected_rows_information = true;
    
    /**
     * Used when annotations are added
     */
    private AnnotationsFileDialog annotation_columns_dialog;    

    
    /**
     * Returns the underlying data file.
     * 
     * @return the underlying data file
     */
    public DataFile<OccurrenceTable> getData()
    {
        return data;
    }
    
    /**
     * Gives the name of the underlying data file
     * 
     * @return name of the underlying file or String for <q>noname</q> if associated file is null
     */
    @Override
    public String toString()
    {
        if (data.getFile()==null)
        {
            return "[noname]";
        } else
            return data.getFile().getName();
    }
    
    /**
     * Initialization of components within this display.
     */
    private void initComponents()
    {
        setBackground(Color.WHITE);
        OccurrenceTableModel M = new OccurrenceTableModel(data.getData(), true);
        table_scroll = new FamilyScroll(M);
        
        table = table_scroll.getDataTable();
        table.setFont(new Font("Serif",Font.PLAIN,LookAndFeel.TABLE_FONT_SIZE));
        M.setDefaultRenderers(table);
        
        //ColoredValueRenderer membership_renderer = new ColoredValueRenderer(LookAndFeel.ABSENCE_COLOR,LookAndFeel.SINGLE_PRESENCE_COLOR,LookAndFeel.MULTI_PRESENCE_COLOR,1,largest_membership);
        //for (int j=0; j<leaves.length; j++)
        //    table.getColumnModel().getColumn(2+j).setCellRenderer(membership_renderer);
        
        // row header data_table
        JTable row_header =table_scroll.getHeaderTable();
        row_header.setFont(table.getFont());

        table_scroll.getViewport().setBackground(getBackground());
        
        // sorting
        TableRowSorter<OccurrenceTableModel> sorter = new TableRowSorter<OccurrenceTableModel>(M);
        table_scroll.setRowSorter(sorter);
        
        Font tp_font_rm = table.getFont().deriveFont(0.8f);
        Font tp_font_it = tp_font_rm.deriveFont(Font.ITALIC);

        selected_rows_information = new JLabel(":");
        selected_rows_information.setFont(tp_font_it.deriveFont(LookAndFeel.TABLE_FONT_SIZE*0.8f));
        selected_rows_information.setOpaque(false);
        selected_rows_information.setLabelFor(table);
        selected_rows_information.setMaximumSize(new Dimension(520,30));
        selected_rows_information.setMinimumSize(selected_rows_information.getMaximumSize());
        selected_rows_information.setPreferredSize(selected_rows_information.getMaximumSize());

        table.getSelectionModel().addListSelectionListener(this);
        table.setRowSelectionInterval(0,0);

        this.setLayout(new BorderLayout());
        this.add(table_scroll, BorderLayout.CENTER);
        this.add(selected_rows_information, BorderLayout.SOUTH);

    }

    /**
     * We are listening to selection changes in the family table.
     * When selection changes, the selection inormation in the bottom
     * bar is updated.
     *
     * @param e a list selection event: the info is updated only if the event is not <q>adjusting</q>
     */
    @Override
    public void valueChanged(ListSelectionEvent e)
    {
        if (!e.getValueIsAdjusting())
        {
            if (update_selected_rows_information)
            {

                int num_selected = table.getSelectedRowCount();


                String selection_info_text = "";
                if (num_selected == 0)
                    selection_info_text = "No rows selected.";
                else
                {
                    if (num_selected==1)
                        selection_info_text = "One row selected";
                    else
                        selection_info_text = Integer.toString(num_selected)+
                                (num_selected==table.getRowCount()?" (all)":"")
                                +" rows selected";
                }

                displaySelectionInfo(selection_info_text);
            }
            else
                update_selected_rows_information = true; // next time
        }
    }

    private void displaySelectionInfo(String text)
    {
        if (selected_rows_information != null)
            selected_rows_information.setText(text);
    }

    private class FamilyScroll extends FrozenColumnsTable<OccurrenceTableModel>
    {
        private FamilyScroll(OccurrenceTableModel M)
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
            return data.getData().getFamilyName(family_idx);
        }

        @Override
        protected String getHeaderToolTip(int column_idx)
        {
            return model.getColumnHeaderToolTip(column_idx);
        }

        @Override
        protected String getCellToolTip(int row, int col)
        {
            return model.getCellToolTip(row, col);
        }

        @Override
        public int selectSimilarFamilies(int family_idx, int col, String command)
        {
            if (command == null)
                return 0;
            else
            {
                int num_selected = super.selectSimilarFamilies(family_idx, col, command);
                String info = "";
                if (num_selected == 1)
                {
                    info = "One row";
                } else
                {
                    info = Integer.toString(num_selected)+" rows";
                }
                info += " selected with "+model.getColumnName(col);
                if (command.equals("eq"))
                    info += "=";
                else if (command.equals("le"))
                    info += "\u2264";
                else if (command.equals("ge"))
                    info += "\u2265";
                info += model.getValueAt(family_idx, col);
                displaySelectionInfo(info);
                update_selected_rows_information = false;

                return num_selected;
            }
        }
    }





    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // ---------------------------------
    // --------------------------------- Filterable interface
    // ---------------------------------
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    
    /**
     * Constructs a new FamilySizeTableDisplay with the selected families.
     * 
     * @return a FamilySizeTableDisplay for the selected families only
     */
    public JComponent newDisplayWithSelectedFamilies()
    {
        int num_selected_rows = table.getSelectedRowCount();
        int[] selected_rows = table.getSelectedRows();
        if (selected_rows==null || selected_rows.length == 0)
        {
            java.awt.Toolkit.getDefaultToolkit().beep();
            return null;
        }
        
        OccurrenceTable data_table = data.getData();
        int num_families = data_table.getNumFamilies();
        boolean[] family_is_selected = new boolean[num_families];
        for (int j=0; j<selected_rows.length; j++)
        {
            int row_idx = selected_rows[j];
            int model_row_idx = table.convertRowIndexToModel(row_idx);
            family_is_selected[model_row_idx] = true;
        }
        return newDisplayWithSelectedFamilies(family_is_selected);
    }
    
    /**
     * Constructs a new FamilySizeTableDisplay with the selected families.
     * 
     * @param selected_rows indexes (by the model) of the selected families
     * @return a FamilySizeTableDisplay for the selected families only; or null if no selection
     */

    public FamilySizeTableDisplay newDisplayWithSelectedFamilies(final boolean[] family_is_selected)
    {
        OccurrenceTable filtered_table = data.getData().tableForFamilies(family_is_selected);
        File filtered_file = new File((File)null, "filt:"+data.getFile().getName());
        DataFile<OccurrenceTable> filtered_data = new DataFile<OccurrenceTable>(filtered_table,filtered_file);
        FamilySizeTableDisplay filtered_display = new FamilySizeTableDisplay(filtered_data);
        addTableModelListener(filtered_display.table_scroll); // properties are shared
        return filtered_display;
    }
    
    public FamilySizeTableDisplay createBinaryTableDisplay()
    {
        OccurrenceTable data_table = data.getData();
        OccurrenceTable binary_table = data_table.binaryTable();
        File binary_file = new File((File)null, "01:"+data.getFile().getName());
        DataFile<OccurrenceTable> binary_data = new DataFile<OccurrenceTable>(binary_table,binary_file);
        FamilySizeTableDisplay binary_display = new FamilySizeTableDisplay(binary_data);
        addTableModelListener(binary_display.table_scroll); // properties are shared
        return binary_display;
    }
    
    public void addTableModelListener(TableModelListener listener)
    {
        table_scroll.getModel().addTableModelListener(listener);
    }
    
    
    
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // ---------------------------------
    // --------------------------------- Saveable interface
    // ---------------------------------
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    
    
    /**
     *  Whether this data set is already saved in a file
     * 
     * @return the answer
     */
    public boolean hasAssociatedFile()
    {
        return data.getFile().getParent()!=null;
    }

    /**
     * Saves the underlying data_table into the associated file.
     * 
     * @throws java.io.IOException
     */
    public void saveData() throws IOException
    {
        saveData(data.getFile());
    }

    /**
     * Saves the underlying data_table into a new file. 
     * 
     * @param f the new file
     * @throws java.io.IOException
     */
    public void saveData(File f) throws IOException
    {
        PrintStream PS = new PrintStream(f);
        PS.println(WorkSpaceCount.getDealerCount(this).getStandardHeader(getClass()));
        
        OccurrenceTable T = data.getData();
        PS.print(T.getFormattedTable());
        
        PS.close();
        data.setFile(f);
        data.setDirty(false);
    }
    
    /**
     * Check the stored dirty bit for the associated file.
     * 
     * @return whether the data_table was modified (Saveable interface)
     */
    public boolean isDirty()
    {
        return data.isDirty();
    }

    /**
     * Sets the <q>dirty</q> bit for the data_table. 
     * 
     * @param dirty indicates that this data_table should be saved
     */
    public void setDirty(boolean dirty)
    {
        data.setDirty(dirty);
    }    
    
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    // ---------------------------------
    // --------------------------------- Loading annotations
    // ---------------------------------
    // ------------------------------------------------------------------------
    // ------------------------------------------------------------------------
    

    public void loadAnnotations(JFrame top_frame)
    {
        FileDialog dialog = null;
        dialog = new FileDialog(top_frame,"Open family annotations file",FileDialog.LOAD);
        dialog.setVisible(true);
       
        String file_name = dialog.getFile();
        String directory = dialog.getDirectory();
        
        Dealer dealer = WorkSpaceCount.getDealer(this); // needed for catching exceptions
        
        if (file_name != null)
        {
            try 
            {
                File selected_file = new File(directory,file_name);
                FileInputStream file_input = new FileInputStream(selected_file);
                
                BufferedReader R = new BufferedReader(new InputStreamReader(file_input));
                Vector<String> first_lines = new Vector<String>();
                do 
                {
                    String line = R.readLine();
                    if (line == null)
                        break;
                    line = line.trim();
                    if (line.startsWith("#")) // remark
                        continue;
                    first_lines.add(line);
                } while (first_lines.size()<10);
                R.close();

                annotation_columns_dialog = new AnnotationsFileDialog(top_frame, first_lines, selected_file);
                Dimension frameD = top_frame.getSize();
                annotation_columns_dialog.pack();
                annotation_columns_dialog.setBounds((int)(0.05*frameD.width),(int)(0.05*frameD.height),800,400);
                annotation_columns_dialog.setVisible(true);
            } catch (java.io.InterruptedIOException E)
            {
                // canceled
            } catch (java.io.IOException E)
            {
                dealer.exceptionCaught(E, "I/O error", "File error while reading family annotations from a file.");
            } catch (Exception E)
            {
                dealer.exceptionCaught(E, "A bug maybe?", "Error while reading family annotations from a file.");
            }
        }
    }
    
    private void annotationColumnSelectionDone(boolean is_okayed)
    {
        if (is_okayed)
        {
            Dealer dealer = WorkSpaceCount.getDealer(this); // needed for catching exceptions
            try 
            {
                File selected_file = annotation_columns_dialog.annotations_file;
                FileInputStream file_input = new FileInputStream(selected_file);
                
                BufferedReader R = new BufferedReader(new InputStreamReader(file_input));
                
                OccurrenceTable family_table = data.getData();
                Hashtable<String,Integer> family_indexes = new Hashtable<String,Integer>();
                for (int family_idx =0; family_idx<family_table.getNumFamilies(); family_idx++)
                {
                    String fname = family_table.getFamilyName(family_idx);
                    family_indexes.put(fname, new Integer(family_idx));
                }
                
                String separator = annotation_columns_dialog.getSelectedSeparator();
                AnnotationColumnsModel table_model = annotation_columns_dialog.table_model;

                int[] include_column = new int[table_model.getColumnCount()];
                for (int col_idx=1; col_idx<include_column.length; col_idx++)
                {
                    String hdr = table_model.column_headers[col_idx];
                    if (hdr.equals(COLUMN_SKIPPED))
                    {
                        include_column[col_idx]=0;
                    } else
                    {
                        include_column[col_idx]=family_table.registerProperty(hdr);
                    }
                }

                String line = null;
                do 
                {
                    line = R.readLine();
                    if (line != null)
                    {
                        line = line.trim();
                        if (line.startsWith("#")) // remark
                            continue;
                        
                        String[] fields = line.split(separator);
                        String fam = fields[0];
                        if (family_indexes.containsKey(fam))
                        {
                            int family_idx = family_indexes.get(fam);
                            for (int col_idx=1; col_idx<include_column.length; col_idx++)
                                if (include_column[col_idx]!=0)
                                {
                                    family_table.setFamilyProperty(family_idx, include_column[col_idx], fields[col_idx]);
                                }
                        } else 
                        {
                        }
                    }
                } while (line != null);
                R.close();
                OccurrenceTableModel otm = table_scroll.getModel();
                otm.fireTableStructureChanged();
                //table_scroll.resetDataColumns();
            } catch (java.io.InterruptedIOException E)
            {
                // canceled
            } catch (java.io.IOException E)
            {
                dealer.exceptionCaught(E, "I/O error", "File error while reading family annotations from a file.");
            } catch (Exception E)
            {
                dealer.exceptionCaught(E, "A bug maybe?", "Error while reading family annotations from a file.");
            }
        }
        annotation_columns_dialog.dispose();
    }
    
    private class AnnotationsFileDialog extends JDialog
                implements ActionListener
    {
        private AnnotationsFileDialog(JFrame mommy, Vector<String> first_lines, File annotations_file)
        {
            super(mommy, "Annotation columns", true);
            this.first_lines = first_lines.toArray(new String[0]);
            this.annotations_file = annotations_file;
            initComponents();
        }
        
        private String[] first_lines;
        private File annotations_file;
        private AnnotationColumnsModel table_model;
        private JTable table;
        private JComboBox separator_combo;
        
        private void initComponents()
        {
            
            HashSet<String> possible_separators = new HashSet<String>();
            possible_separators.add(SEPARATOR_WHITESPACE);
            for (int line_idx=0; line_idx<first_lines.length; line_idx++)
            {
                String[] chars = first_lines[line_idx].split("");
                for (int j=0; j<chars.length; j++)
                    if (chars[j].length()!=0)
                    {
                        boolean add_this = !Character.isWhitespace(chars[j].charAt(0));
                        if (!add_this)
                            add_this = "\t".equals(chars[j]) || " ".equals(chars[j]);
                        if (add_this)
                            possible_separators.add(chars[j]);
                    }
            }
            String[] all_possible_separators = possible_separators.toArray(new String[0]);
            java.util.Arrays.sort(all_possible_separators);
            // comma, tab, whitespace, space at the beginning
            {
                boolean has_comma = false;
                boolean has_space = false;
                boolean has_tab = false;
                for (int i=0; i<all_possible_separators.length; i++)
                    if (all_possible_separators[i].equals("\t"))
                        has_tab=true;
                    else if (all_possible_separators[i].equals(" "))
                        has_space = true;
                    else if (all_possible_separators[i].equals(","))
                        has_comma = true;
                int comma_idx = 0;
                int whitespace_idx = comma_idx+(has_comma?1:0);
                int tab_idx = whitespace_idx+1;
                int space_idx = tab_idx+(has_tab?1:0);
                for (int i=0; i<all_possible_separators.length; i++)
                {
                    String sep = all_possible_separators[i];
                    int wanted_idx = i;
                    if (sep.equals(SEPARATOR_WHITESPACE))
                        wanted_idx = whitespace_idx;
                    else if (sep.equals(","))
                    {
                        sep = SEPARATOR_COMMA;
                        wanted_idx = comma_idx;
                    }
                    else if (sep.equals("\t"))
                    {
                        sep = SEPARATOR_TAB;
                        wanted_idx = tab_idx;
                    }
                    else if (sep.equals(" "))
                    {
                        sep = SEPARATOR_SPACE;
                        wanted_idx = space_idx;
                    } 
                        
                    if (wanted_idx != i)
                    {
                        all_possible_separators[i]=all_possible_separators[wanted_idx];
                        all_possible_separators[wanted_idx]=sep;
                        i--; // check again
                    }
                }
            }
            JPanel main_panel = new JPanel();
            BoxLayout layout = new BoxLayout(main_panel, BoxLayout.PAGE_AXIS);
            main_panel.setLayout(layout);
            Box dataB = new Box(BoxLayout.LINE_AXIS);

            separator_combo = new JComboBox(all_possible_separators);
            separator_combo.setEditable(true);
            separator_combo.addActionListener(this);
            separator_combo.setMaximumSize(separator_combo.getMinimumSize());
            Box comboB = new Box(BoxLayout.PAGE_AXIS);
            
            JLabel comboL = new JLabel("Field separator:");
            comboL.setFont(new Font("Serif",Font.BOLD,LookAndFeel.TABLE_FONT_SIZE));
            comboB.add(comboL);
            comboB.add(separator_combo);
            comboB.add(Box.createVerticalGlue());
            dataB.add(comboB);

            dataB.add(Box.createHorizontalGlue());
            
            JComboBox column_headerCB = new JComboBox();
            //column_headerCB.addItem(COLUMN_FAMILY);
            column_headerCB.addItem(COLUMN_SKIPPED);
            column_headerCB.addItem(COLUMN_CATEGORY);
            column_headerCB.addItem(COLUMN_DESCRIPTION);
            column_headerCB.setEditable(true);
            final DefaultCellEditor column_header_editor = new DefaultCellEditor(column_headerCB);
            table_model = new AnnotationColumnsModel(first_lines, getSeparator(all_possible_separators[0]));
            table = new JTable(table_model)
            {
                @Override
                public boolean isCellEditable(int row_idx, int column_idx)
                {
                    int c = convertColumnIndexToModel(column_idx);
                    return (row_idx==0 && c!=0);
                }
                
                @Override
                public TableCellEditor getCellEditor(int row_idx, int col_idx)
                {
                    if (row_idx==0) return column_header_editor;
                    else return super.getCellEditor(row_idx,col_idx);
                }
            };
            table.setDefaultRenderer(String.class, new DefaultTableCellRenderer()
            {
                @Override
                public Component getTableCellRendererComponent(JTable table,
                                                               Object value,
                                                               boolean isSelected,
                                                               boolean hasFocus,
                                                               int row,
                                                               int column)
                {
                    Component C = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, (row<1?-1:row), column);
                    if (row<1)
                    {
                        C.setFont(C.getFont().deriveFont(Font.BOLD));
                        C.setBackground(LookAndFeel.SMOKY_BACKGROUND);
                    }
                    return C;
                }
            });
            table.setRowSelectionAllowed(false);
            table.setColumnSelectionAllowed(false);
            table.setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
            table.setRowHeight(column_headerCB.getPreferredSize().height);
            //table.setMaximumSize(table.getMinimumSize());
            
            JScrollPane table_scroll = new JScrollPane(table);
            table_scroll.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
            table_scroll.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
            Box tableB = new Box(BoxLayout.PAGE_AXIS);
            tableB.add(table_scroll);
            tableB.add(Box.createVerticalGlue());
            dataB.add(tableB);
            
            main_panel.add(dataB);
            main_panel.add(Box.createVerticalGlue());

            JButton ok_button = new JButton("OK");
            ok_button.setActionCommand("OK");
            ok_button.addActionListener(this);
            JButton cancel_button = new JButton("Cancel");
            cancel_button.setActionCommand("Cancel");
            cancel_button.addActionListener(this);
            
            Box button_box = new Box(BoxLayout.LINE_AXIS);
            button_box.add(Box.createHorizontalGlue());
            button_box.add(cancel_button);
            button_box.add(ok_button);
            
            main_panel.add(button_box);
            
            main_panel.add(Box.createVerticalGlue());

            add(main_panel);
        };
        
        private String getSelectedSeparator()
        {
            return getSeparator((String)separator_combo.getSelectedItem());
        }

        private String getSeparator(String menu_item_name)
        {
            if (SEPARATOR_COMMA.equals(menu_item_name))
                return ",";
            else if (SEPARATOR_WHITESPACE.equals(menu_item_name))
                return "\\s+";
            else if (SEPARATOR_TAB.equals(menu_item_name))
                return "\t";
            else if (SEPARATOR_SPACE.equals(menu_item_name))
                return " ";
            else return menu_item_name;
        }
        
        public void actionPerformed(ActionEvent E)
        {
            Object src = E.getSource();
            if (src instanceof JButton)
            {
                String cmd = E.getActionCommand();
                if ("OK".equals(cmd))
                {
                    annotationColumnSelectionDone(true);
                } else if ("Cancel".equals(cmd))
                {
                    annotationColumnSelectionDone(false);
                }
            } else if (src instanceof JComboBox)
            {
                JComboBox srcCB = (JComboBox) src;
                String separator_name = (String)srcCB.getSelectedItem();
                table_model.setSeparator(getSeparator(separator_name));
                table_model.fireTableStructureChanged();
                //table.setRowHeight(column_headerCB.getPreferredSize().height);
            } else throw new RuntimeException("Unknown ActionEvent in FamilySizeTableDisplay.AnnotationsfileDialog: "+E);
        }
    }
    
    private class AnnotationColumnsModel extends AbstractTableModel 
    {
        private AnnotationColumnsModel(String[] first_lines, String separator)
        {
            this.first_lines = first_lines;
            setSeparator(separator);
        }
        
        private String[] first_lines;
        private String separator;
        private String[][] fields;
        private String[] column_headers;
        
        private void setSeparator(String separator)
        {
            this.separator = separator;
            fields = new String[first_lines.length][];
            for (int i=0; i<first_lines.length; i++)
            {
                String line = first_lines[i];
                fields[i]=line.split(separator);
            }
            int ncol = 0;
            for (int i=0; i<fields.length;i++)
                if (fields[i].length > ncol)
                    ncol=fields[i].length;
            column_headers=new String[ncol];
            column_headers[0] = COLUMN_FAMILY;
            for (int i=1; i<ncol; i++)
                column_headers[i] = COLUMN_SKIPPED;
        }
        
        public int getColumnCount()
        {
            return column_headers.length;
        }
        
        @Override
        public String getColumnName(int column_idx)
        {
            return "C"+Integer.toString(1+column_idx);
        }

        public int getRowCount()
        {
            return 1+fields.length;
        }
        
        @Override
        public Object getValueAt(int row_idx, int col_idx)
        {
            if (row_idx==0)
            {
                return column_headers[col_idx];
            } else
                row_idx--;

            if (col_idx<fields[row_idx].length)
                return fields[row_idx][col_idx];
            else
                return null;
        }
        
        @Override 
        public Class getColumnClass(int col_idx)
        {
            return String.class;
        }
        
        @Override 
        public void setValueAt( Object value,int row_idx, int col_idx)
        {
            if (row_idx==0)
                column_headers[col_idx] = (String)value;
            else
                super.setValueAt(value, row_idx, col_idx);
        }
    }
    
    
}
