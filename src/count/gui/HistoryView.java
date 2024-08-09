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
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JEditorPane;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.SwingWorker;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.gui.kit.ColoredValueRenderer;
import count.gui.kit.RoundedDouble;
import count.gui.kit.TableScroll;
import count.io.CommandLine;
import count.io.DataFile;
import count.io.ExportableData;
import count.model.SimulatedEvolution;

import static count.gui.AnnotatedTablePanel.TABLE_FONT_SIZE;
import static count.Count.THREAD_UNIT_TASK;
import static count.Count.THREAD_PARALLELISM;
/**
 * Template component (JSplitPane) to show ancestral reconstructions.
 * Once set up, the computation for opulating the tables is launched 
 * by {@link #computeAll()}. Extending classes define {@link #computeHistoriesPrepare()}, 
 * {@link #computeFamily(int)} and {@link #computeDone()}. 
 * 
 * Data set: data in {@link #table_model} for table {@link #table_data}
 * Columnwise sums in selected rows: {@link #tree_control}
 * 
 * Event chains:
 * 1. user row selection
 * 		table:ListSelectionEvent 
 * 			. update selection
 * 2. double-click selection	
 * 		table:ListSelectionEvent  
 * 			. don't update row selection 
 * 1 or 2 
 * 		table:ListSelectionEvent
 * 			. update table values in node table; update tree 
 * 				nodetable:TableModelEvent 
 * 				
 * 
 * @author csuros
 *
 */
public class HistoryView 
		extends JSplitPane
		implements Session.FamilySelection
				, ExportableData
{
	protected final static Color FAMILY_EXPAND_COLOR = RatesTreePanel.RATES_DUPLICATION_COLOR;
	protected final static Color FAMILY_GAIN_COLOR = RatesTreePanel.RATES_GAIN_COLOR;
	protected final static Color FAMILY_LOSE_COLOR = RatesTreePanel.RATES_LOSS_COLOR;
	protected final static Color FAMILY_CONTRACT_COLOR = new Color(229,111,187); // Grape
	protected final static Color FAMILY_LL_COLOR = new Color(255,126,121); // Salmon
	protected final static Color FAMILY_PRESENT_COLOR = new Color(118, 214, 255); // Sky // new Color(128, 255, 255); // Ice

	protected HistoryView(
			IndexedTree phylo, 
			DataFile<AnnotatedTable> table_data,
			boolean has_copy_numbers, boolean has_only_integers)
	{
		super(JSplitPane.VERTICAL_SPLIT);
		this.table_data = table_data;
		this.has_copy_numbers = has_copy_numbers;
		this.has_only_integers = has_only_integers;
		this.table_model = new HistoryModel(phylo, table_data.getContent());
		
		initTableModel();
		initComponents();
		
	}
	
	static 
	{
		if (THREAD_PARALLELISM>1)
			executor = Executors.newWorkStealingPool(THREAD_PARALLELISM);
		else
			executor = null;
	}
	private static final ExecutorService executor; // only one thread pool shared by the instances 
	
	/**
	 * The data set we are working with.
	 */
	protected final DataFile<AnnotatedTable> table_data;
	/**
	 * Nonbinary or binary data
	 */
	protected final boolean has_copy_numbers;
	/**
	 * Parsimony or probability
	 */
	protected final boolean has_only_integers;
	
	/**
	 * TableModel based on our data structure.
	 */
	protected final HistoryModel table_model;

	/*
	 * GUI components
	 */
	/**
	 * Top component.
	 */
	private TableScroll<HistoryModel> table_scroll;
	/**
	 * Tree in bottom component 
	 */
	private Zoom<HistoryTreePanel> tree_control;
	
	private Box lineage_control;
	
	protected Zoom<HistoryTreePanel> getTreeControl(){ return tree_control;}
	
	protected TableScroll<HistoryModel> getTableScroll(){ return table_scroll;}
	
	protected Box getLineageControl() { return lineage_control;} 
	
	private JProgressBar computation_progress;
	protected JButton computation_cancel;
	
	/**
	 * Bottom component
	 */
	private JTabbedPane selected_rows_display;
    
	
	/*
	 * History reconstruction column codes 
	 */
	private static final String MEMBER_COUNT	= ":n";
	private static final String MEMBER_CHANGE	= ":d";
	private static final String MEMBER_BIRTH	= ":+";
	private static final String MEMBER_DEATH	= ":-";
	private static final String FAMILY_PRESENT	= ":p";
	private static final String FAMILY_MULTI	= ":m";
	private static final String FAMILY_GAIN		= ":g";
	private static final String FAMILY_LOSE		= ":l";
	private static final String FAMILY_EXPAND	= "++";
	private static final String FAMILY_CONTRACT	= "--";
	private static final String MEMBER_MAX      = ":x";
	
	/*
	 * Data structure for history reconstruction
	 */
	/**
	 * Convenience class to spare on parametric typization.
	 */
	protected class LineageColumns extends ArrayList<HistoryModel.Column<?>>
	{}
	
	public HistoryTreePanel getTreePanel()
	{
		return tree_control.getTreePanel();
	}
	
	public TableScroll<? extends TableModel> getSelectionTotals()
	{
		return tree_control.getTreePanel().asTableScroll();
	}
	
	/**
	 * Ancestral copy numbers
	 */
	protected LineageColumns table_member_count=null;
	
	protected HistoryTreePanel.NodeData lineage_member_count;
	
	protected LineageColumns table_member_max = null;
	/**
	 * Copies gained.
	 */
	protected LineageColumns table_member_birth = null;
	protected HistoryTreePanel.NodeData lineage_member_birth;
	/**
	 * Copies lost.
	 */
	protected LineageColumns table_member_death = null;
	protected HistoryTreePanel.NodeData lineage_member_death;
	/**
	 * Family with at least one member.
	 */
	protected LineageColumns table_family_present = null;
	protected HistoryTreePanel.NodeData lineage_family_present;
	
	/**
	 * Family with more than one member.
	 */
	protected LineageColumns table_family_multi = null;
	protected HistoryTreePanel.NodeData lineage_family_multi;
	
	/**
	 * Family gained here.
	 */
	protected LineageColumns table_family_gain = null;
	protected HistoryTreePanel.NodeData lineage_family_gain;

	/**
	 * Family completely lost.
	 */
	protected LineageColumns table_family_lose = null;
	protected HistoryTreePanel.NodeData lineage_family_lose;

	/**
	 * Family expanded.
	 */
	protected LineageColumns table_family_expand = null;
	protected HistoryTreePanel.NodeData lineage_family_expand;
	/**
	 * Family contracted.
	 */
	protected LineageColumns table_family_contract = null;
	protected HistoryTreePanel.NodeData lineage_family_contract;
	
	/**
	 * Profile probability
	 */
	protected HistoryModel.Column<RoundedDouble> table_family_score  = null;
	
	private void initComponents()
	{
		class FamilyScroll extends TableScroll<HistoryModel>
		{
            private FamilyScroll()
            {
//                super(table_model,2);
                super(table_model,table_model.firstHistoryColumn());
            }
            
            @Override
            protected int getPreferredColumnWidth(int idx)
            {
                if (idx == 1 ) // family name
                    return 120;
                else if (AnnotatedTable.PhyleticProfile.class == getModel().getColumnClass(idx))
                {
                	int num_nodes = table_model.getTree().getNumNodes();
                	int w;
                	if (num_nodes<100)
                		w = 2*num_nodes;
                	else
                	{
                		w = num_nodes;
                		while (w>400)
                			w = w/2;
                	}
                	return w;
                }
                else
                    return super.getPreferredColumnWidth(idx);
            }

            @Override
            public String getRowName(int family_idx)
            {
                return table_data.getContent().getFamilyName(family_idx);
            }

            @Override
            protected String getHeaderToolTip(int column_idx)
            {
                String tt = model.getColumnDescription(column_idx);

                
                return tt +
                        "; click to sort rows, drag to rearrange columns";//, or double-click to split the table";
            }
            
            @Override
            protected String getCellToolTip(int row, int col)
            {
            	String tooltip;
            	if (col<model.firstHistoryColumn())
            		tooltip = model.getCellToolTip(row, col);
            	else
            		tooltip = super.getCellToolTip(row, col);
            	return tooltip;
            }
            
		} // class FamilyScroll
//		this.setDividerLocation(240);
		this.setResizeWeight(0.25);
        this.setOneTouchExpandable(true);
        table_scroll = new FamilyScroll();
        
        JTable table = table_scroll.getDataTable();
        table.setFont(new Font("Serif",Font.PLAIN,TABLE_FONT_SIZE));
        table_model.setDefaultRenderers(table);
        
        // row header data_table
        JTable row_header =table_scroll.getHeaderTable();
        row_header.setFont(table.getFont());
        table_model.setDefaultRenderers(row_header);

        table_scroll.getViewport().setBackground(getBackground());
        
        // sorting
        TableRowSorter<HistoryModel> sorter = new TableRowSorter<>(table_model);
        table_scroll.setRowSorter(sorter);

        Font tp_font_rm = table.getFont().deriveFont(0.8f);
        Font tp_font_it = tp_font_rm.deriveFont(Font.ITALIC);

//        JLabel selected_rows_information = ;  // JLabel(":");
         
        HistoryTreePanel tree_panel = createLineageModel();
        tree_control = new Zoom<>(tree_panel);


        Box control_bar = tree_control.getControlBar();
		control_bar.removeAll();
		
		Box data_info = new Box(BoxLayout.PAGE_AXIS);
		JLabel selection_info= table_scroll.createRowSelectionInfo();
        computation_progress = new JProgressBar(0, table_data.getContent().getFamilyCount());
        computation_progress.setString("Preparing ./"+table_data.getContent().getFamilyCount());
        computation_progress.setStringPainted(true);
        computation_progress.setBorderPainted(true);
		computation_progress.setFont(selection_info.getFont());
		computation_progress.setPreferredSize(new Dimension(selection_info.getWidth(), selection_info.getFont().getSize()));
		data_info.add(computation_progress);
		data_info.add(selection_info);
		control_bar.add(data_info);
		
		computation_cancel = new JButton("Cancel");
		computation_cancel.setActionCommand( // command saves the position within the control bar 
				Integer.toString(control_bar.getComponentCount()));
		computation_cancel.setToolTipText("Stops computing the family histories");
		control_bar.add(computation_cancel);
		computation_cancel.addActionListener(clik->history_task.cancel(true));
		
		Box button_box = new Box(BoxLayout.LINE_AXIS);
		button_box.setBorder(BorderFactory.createTitledBorder("Charts"));
		
		
		List<JCheckBox> charts_wanted = tree_panel.getChartControls();
		if (has_copy_numbers)
		{
			int widx=0;
			while (widx<charts_wanted.size()) // add them in 2-1-2 pattern
			{
				if (widx+1<charts_wanted.size())
				{
					Box two_cb = new Box(BoxLayout.PAGE_AXIS);
					two_cb.add(charts_wanted.get(widx++));
					two_cb.add(charts_wanted.get(widx++));
					button_box.add(two_cb);
				}
				if (widx<charts_wanted.size())
				{
					button_box.add(charts_wanted.get(widx++));
				}
			}
		} else
		{
			for (JCheckBox cb: charts_wanted)
			{
				button_box.add(cb);
			}
		}
		control_bar.add(Box.createHorizontalGlue());
		control_bar.add(tree_panel.getScalingControl());
		control_bar.add(button_box);
		control_bar.add(Box.createHorizontalGlue());
		control_bar.add(tree_control.getSpinner());
		
        selected_rows_display = new JTabbedPane();
		selected_rows_display.addTab("Tree", tree_control);
		
		JPanel bottom_table = new JPanel(new BorderLayout());
		bottom_table.add(tree_panel.asTableScroll(), BorderLayout.CENTER);
		lineage_control = new Box(BoxLayout.LINE_AXIS);
		lineage_control.add(table_scroll.createRowSelectionInfo());
		bottom_table.add(lineage_control, BorderLayout.SOUTH);
		
//		String table_legend = "<p>The inferred counts are posterior expectations across all families."
//				+ "Copy statistics (copies, births and deaths) "
//				+ "count "
//				+ "	<b>ancestral</b> copies at each node <var>u</var>, or "
//				+ "the ancestral copy change "
//				+ "in the lineage leading to <var>u</var>, summing across selected families. "
//				+ "Ancestral copies have at least one "
//				+ "descendant ortholog copy at the leaves in the subtree of <var>u</var>."
//				+ "Copy statistics thus track the provenance of the observed copies at "
//				+ "the terminal nodes. Deaths are copies at the parent which <var>u</var> "
//				+ "does not inherit; births are anscestral copies originating at <var>u</var> "
//				+ "by duplication or gain."
//				+ "</p>"
//				+ "<p>Family statistics track the history of family repertoires, and  "
//				+ " count the families that have at least one copy "
//				+ "at the ancestor, with or without descendant orthologs. "
//				+ "(Thus, more families may be inferred than ancestral copies.) "
//				+ "Families with multiple copies are also inferred, and changes "
//				+ "to and from 0 (gain and loss), or between single- and multiple-members "
//				+ "(contraction and expansion)."
//				+"</p>"
////				+ "<p>In the family table, "
////				+ "rarity is the negative log-likelihood of the family profile, "
////				+ "as -10 log<sub>10</sub> <var>L</var>(profile). Family statistics "
////				+ "(gains, losses, expansions, contractions) "
////				+ "are summed across all nodes</p>."
//				;
//		JEditorPane table_explain = new JEditorPane("text/html", table_legend);
//        table_explain.setEditable(false);		
//		selected_rows_display.addTab("Legend", table_explain);

        
		
		selected_rows_display.addTab("Table", bottom_table);
		
		
        this.setTopComponent(table_scroll);
        this.setBottomComponent(selected_rows_display);
        
        table.setRowSelectionInterval(0,0);
        
        
	}
	
	protected JTabbedPane getSelectedRowsDisplay()
	{
		return this.selected_rows_display;
	}
	
	private HistoryTreePanel createLineageModel()
	{
		HistoryTreePanel P = new HistoryTreePanel(table_scroll);
		boolean is_binary_table = table_data.getContent().isBinaryTable();
		
		if (has_copy_numbers && !is_binary_table)
		{
			lineage_member_count = P.newColumns("Copies "+MEMBER_COUNT+"", table_member_count);
			JCheckBox count_cb = P.showNodeStatistics("Conserved copies", lineage_member_count);
			//count_cb.setForeground(FAMILY_PRESENT_COLOR);
			lineage_member_birth = P.newColumns("Births "+MEMBER_BIRTH+"", table_member_birth);
			lineage_member_death = P.newColumns("Deaths "+MEMBER_DEATH+"", table_member_death);
			P.showChangeStatistics("Copy change", lineage_member_birth, lineage_member_death);;
		}
		lineage_family_present = P.newColumns("Families "+FAMILY_PRESENT+"", table_family_present);
		if (has_copy_numbers && !is_binary_table)
		{
			lineage_family_multi = P.newColumns("Multi "+FAMILY_MULTI+"", table_family_multi);
			P.showNodeStatistics("Families", lineage_family_present, lineage_family_multi);
		} else
		{
			P.showNodeStatistics("Families", lineage_family_present);
		}
		lineage_family_gain = P.newColumns("Gains "+FAMILY_GAIN+"", table_family_gain);
		lineage_family_lose = P.newColumns("Losses "+FAMILY_LOSE+"", table_family_lose);
		P.showChangeStatistics("Loss/gain", lineage_family_gain, lineage_family_lose);
		if (has_copy_numbers && !is_binary_table)
		{
			lineage_family_expand = P.newColumns("Expansions "+FAMILY_EXPAND+"", table_family_expand);
			lineage_family_contract = P.newColumns("Contractions "+FAMILY_EXPAND+"", table_family_contract);
			P.showChangeStatistics("Contraction/expansion", lineage_family_expand, lineage_family_contract);		
		}		
		
		return P;

	}
	
	/**
	 * Creates and adds the columns to {@link #table_model} 
	 */
	private void initTableModel()
	{
		IndexedTree phylo = table_model.getTree();
		int num_nodes = phylo.getNumNodes();
		table_family_present = new LineageColumns(); // 
		table_family_gain = new LineageColumns(); // gain
		table_family_lose = new LineageColumns(); // loss
		
		boolean is_simulated = table_model.getTable() instanceof SimulatedEvolution.Table;
		
		if (has_copy_numbers)
		{
			table_member_birth = new LineageColumns();
			table_member_death = new LineageColumns();

			table_family_multi = new LineageColumns();
			table_member_count = new LineageColumns();
			table_family_expand = new LineageColumns(); // expansion 
			table_family_contract = new LineageColumns(); // contraction

			if (has_only_integers)
			{
				// numerical parsimony
				List<HistoryModel.Column<Integer>> mcL = table_model.newIntColumns(MEMBER_COUNT, "Ancestral copy number", false);
				assert (num_nodes==mcL.size());
				table_member_count.addAll(mcL);
				
				// family changes inferred
				for (int node=0; node<num_nodes; node++)
				{
					HistoryModel.Column<Integer> node_col = mcL.get(node);
					HistoryModel.Column<Integer> parent_col;
					if (phylo.isRoot(node))
						parent_col=null;
					else
						parent_col=mcL.get(phylo.getParent(node));
					table_family_gain.add(table_model.newChangeColumn(node_col, parent_col, 0, 
							table_model.getHeader(FAMILY_GAIN, node), table_model.getToolTip("Family gained", node, true)));
					table_family_lose.add(table_model.newChangeColumn(parent_col, node_col, 0, 
							table_model.getHeader(FAMILY_LOSE, node), table_model.getToolTip("Family lost", node, true)));
					table_family_expand.add(table_model.newChangeColumn(node_col, parent_col, 1, 
							table_model.getHeader(FAMILY_EXPAND, node), table_model.getToolTip("Family expanded (from size 1)", node, true)));
					table_family_contract.add(table_model.newChangeColumn(parent_col, node_col, 1, 
							table_model.getHeader(FAMILY_CONTRACT, node), table_model.getToolTip("Family contracted (to size 1)", node, true)));
					table_family_present.add(table_model.newThresholdColumn(node_col, 1, 
							table_model.getHeader(FAMILY_PRESENT, node),table_model.getToolTip("Family has at least 1 member", node, false)));
					table_family_multi.add(table_model.newThresholdColumn(node_col, 2, 
							table_model.getHeader(FAMILY_MULTI, node),table_model.getToolTip("Family has more than 1 member", node, false)));
					table_member_birth.add(table_model.newNonnegativeDifferenceColumn(node_col, parent_col,
							table_model.getHeader(MEMBER_BIRTH, node), table_model.getToolTip("Cop[y number increase", node, true)));
					table_member_death.add(table_model.newNonnegativeDifferenceColumn(parent_col, node_col,
							table_model.getHeader(MEMBER_DEATH, node), table_model.getToolTip("Cop[y number decrease", node, true)));
				}
				
			} else
			{
				// posteriors or simulation
				
				// copy count columns
				String eS = (is_simulated?"Number of ":"Expected number of ");
				table_member_count.addAll(table_model.newDoubleColumns(MEMBER_COUNT, eS+"conserved ancestral copies", false));
				table_member_birth.addAll(table_model.newDoubleColumns(MEMBER_BIRTH, eS+"gained ancestral copies", true));
				table_member_death.addAll(table_model.newDoubleColumns(MEMBER_DEATH, eS+"lost ancestral copies", true));
				
				// family change counts
				String pS = (is_simulated?"Indicator ":"Probability ");
				table_family_present.addAll(table_model.newDoubleColumns(FAMILY_PRESENT, pS+"of at least one copy", false));
				table_family_multi.addAll(table_model.newDoubleColumns(FAMILY_MULTI, pS+"of more than one copies", false));
				table_family_gain.addAll(table_model.newDoubleColumns(FAMILY_GAIN, pS+"that the family was acquired", true));
				table_family_lose.addAll(table_model.newDoubleColumns(FAMILY_LOSE, pS+"that the family was completely lost", true));
				table_family_expand.addAll(table_model.newDoubleColumns(FAMILY_EXPAND, pS+"that the family expanded (from size 1)", true));
				table_family_contract.addAll(table_model.newDoubleColumns(FAMILY_CONTRACT, pS+"that the family contracted (to size 1)", true));
			} 

		} else
		{
			// 0..1 data 
			if (has_only_integers)
			{
				// binary parsimony
				List<HistoryModel.Column<Integer>> fpL = table_model.newIntColumns(FAMILY_PRESENT, "Present", false);
				table_family_present.addAll(fpL);
				// family changes inferred
				for (int node=0; node<num_nodes; node++)
				{
					HistoryModel.Column<Integer> node_col = fpL.get(node);
					HistoryModel.Column<Integer> parent_col;
					if (phylo.isRoot(node))
						parent_col=null;
					else
						parent_col=fpL.get(phylo.getParent(node));
					table_family_gain.add(table_model.newChangeColumn(node_col, parent_col, 0, 
							table_model.getHeader(FAMILY_GAIN, node), table_model.getToolTip("Family gained", node, true)));
					table_family_lose.add(table_model.newChangeColumn(parent_col, node_col, 0, 
							table_model.getHeader(FAMILY_LOSE, node), table_model.getToolTip("Family lost", node, true)));
				}			
			} else
			{
				// posteriors presence/absence
				String pS = (is_simulated?"Presence probability":"Presence indicator");
				List<HistoryModel.Column<RoundedDouble>> fpL = table_model.newDoubleColumns(FAMILY_PRESENT,pS, false);
				table_family_present.addAll(fpL);
				String ppS = (is_simulated?"Posterior probability of ":"Indicator of ");
				table_family_gain.addAll(table_model.newDoubleColumns(FAMILY_GAIN, "gain", true));
				table_family_lose.addAll(table_model.newDoubleColumns(FAMILY_LOSE, "loss", true));
			}
		}
		if (has_only_integers)
		{
			if (has_copy_numbers)
			{
				// max in parsimony
				table_member_max = new LineageColumns();
				List<HistoryModel.Column<Integer>> mcX = table_model.newIntColumns(MEMBER_MAX, "Maximum copy number", false);			
				table_member_max.addAll(mcX);
			}
			table_family_score = table_model.newDoubleColumn("Score", "Parsimony score for history reconstruction");
		} else
		{
			if (is_simulated)
			{
				table_family_score = table_model.newDoubleColumn("Events", "Number of copy births and deaths");
			} else
			{
				table_family_score = table_model.newDoubleColumn("Rarity", "Negative log-likeihood (deciban)");
			}
		}

		HistoryModel.Column<?> tot_gain;
		HistoryModel.Column<?> tot_lose;
		
		if (has_only_integers)
		{
			tot_gain = table_model.newSumColumn(table_family_gain, "Gains", "Number of lineages in which the family was gained.");
			tot_lose = table_model.newSumColumn(table_family_lose, "Losses", "Number of lineages in which the family was lost.");
		} else
		{
			tot_gain = table_model.newSumColumn(table_family_gain, "Gains", "Expected number of lineages in which the family was gained.");
			tot_lose = table_model.newSumColumn(table_family_lose, "Losses", "Expected number of lineages in which the family was lost.");			
		}
		table_model.add(table_family_score, FAMILY_LL_COLOR);
		table_model.add(tot_gain, FAMILY_GAIN_COLOR);
		table_model.add(tot_lose, FAMILY_LOSE_COLOR);

		if (has_copy_numbers && !table_data.getContent().isBinaryTable())
		{
			HistoryModel.Column<?> tot_expand = table_model.newSumColumn(table_family_expand, "Expansions", "Number of lineages in which the family expanded.");
			HistoryModel.Column<?> tot_contract = table_model.newSumColumn(table_family_contract, "Contractions", "Number of lineages in which the family contracted");

			List<HistoryModel.Column<?>> table_member_change=new ArrayList<>(num_nodes);
			if (has_only_integers)
			{
				// parsimony: difference is thisnode.count-parent.count
				for (int node=0; node<num_nodes; node++)
				{
					HistoryModel.Column<Integer> node_col = (HistoryModel.Column<Integer>) table_member_count.get(node);
					HistoryModel.Column<Integer> parent_col;
					if (phylo.isRoot(node))
						parent_col=null;
					else 
						parent_col = (HistoryModel.Column<Integer>) table_member_count.get(phylo.getParent(node));

					HistoryModel.Column<?> diff = table_model.newDifferenceColumn(node_col, parent_col, table_model.getHeader(MEMBER_CHANGE, node), table_model.getToolTip("Net change",node,true));
					table_member_change.add(diff);
				}
			} else
			{
				// posteriors: diff is gain-loss (which is == thisnode.count-parent.count)
				for (int node=0; node<num_nodes; node++)
				{
					HistoryModel.Column<RoundedDouble> birth_col = (HistoryModel.Column<RoundedDouble>) table_member_birth.get(node);
					HistoryModel.Column<RoundedDouble> death_col = (HistoryModel.Column<RoundedDouble>) table_member_death.get(node);
					HistoryModel.Column<RoundedDouble> diff = table_model.newDifferenceColumn(birth_col, death_col, table_model.getHeader(MEMBER_CHANGE, node), table_model.getToolTip("Net change", node, true));
					table_member_change.add(diff);
				}
			}
			
			table_model.add(tot_expand, FAMILY_EXPAND_COLOR);
			table_model.add(tot_contract, FAMILY_CONTRACT_COLOR);
			
			table_model.addAll(table_member_count, FAMILY_GAIN_COLOR);
			table_model.addAll(table_member_change, null);
			
			if (has_only_integers)
			{
				table_model.addAll(table_member_max, FAMILY_GAIN_COLOR);
			}
			else
			{
				table_model.addAll(table_member_birth, FAMILY_GAIN_COLOR);
				table_model.addAll(table_member_death, FAMILY_LOSE_COLOR);
			}
			
			table_model.addAll(table_family_present, FAMILY_PRESENT_COLOR);
			table_model.addAll(table_family_multi, FAMILY_PRESENT_COLOR);

			table_model.addAll(table_family_gain, FAMILY_GAIN_COLOR);
			table_model.addAll(table_family_lose, FAMILY_LOSE_COLOR);
			table_model.addAll(table_family_expand, FAMILY_EXPAND_COLOR);
			table_model.addAll(table_family_contract, FAMILY_CONTRACT_COLOR);
			
		} else // binary
		{
//			if (!has_only_integers) // or else Dollo parsimony
			table_model.addAll(table_family_present, FAMILY_PRESENT_COLOR);
			table_model.addAll(table_family_gain, FAMILY_GAIN_COLOR);
			table_model.addAll(table_family_lose, FAMILY_LOSE_COLOR);
		}
	}
	

	
	
    public void setPhyleticProfileColoring(Color[] leaf_colors)
    {
    	table_scroll.getModel().setPhyleticProfileColoring(table_scroll.getHeaderTable(), leaf_colors);
    	table_scroll.getModel().setPhyleticProfileColoring(table_scroll.getDataTable(), leaf_colors);
    }
	
	
	
	
	protected SwingWorker<Void,Integer> history_task=null;
		
	/**
	 * Prepares the computation of the statistics for families; 
	 * executed on working thread. 
	 * 
	 */
	protected void computeHistoriesPrepare()
	{
//		computation_progress.setIndeterminate(false);
//		computation_progress.setMaximum(table_data.getContent().getFamilyCount());
//		computation_progress.setValue(0);
//		computation_progress.setForeground(new Color(1.0f,0.8f,0.1f,0.84f));
//		computation_progress.setBackground(new Color(1.0f,0.8f,0.1f,0.84f));
	}
	
	/**
	 * Computes the statistics for a family; 
	 * executed on working thread. 
	 * 
	 * @param family_idx
	 */
	protected void computeFamily(int family_idx)
	{
		// nothing here 
	}
	
	/**
	 * Called on Event thread.
	 */
	protected void computeDone()
	{
		
	}
	

	private static final boolean DEBUG_THREADS=false;
	
	/**
	 * Sets {@link #history_task} and launches execution
	 * in the background on a working thread.
	 * Called on the event thread. 
	 */
	protected void computeAll()
	{
		final int nF = table_data.getContent().getFamilyCount();
		class HistoryWorker extends SwingWorker<Void,Integer>
		{
			private int num_families_done=0;
			
			@Override
			protected Void doInBackground()
			{
				num_families_done = 0;
				computeHistoriesPrepare();
				
				
				if (executor != null && nF>THREAD_UNIT_TASK)
				{
					final BlockingQueue<Integer> families_done=new ArrayBlockingQueue<>(nF);
					
					class ComputingTask implements Runnable //Callable<int[]>
					{
						ComputingTask(int fmin, int fmax)
						{
							this.range = new int[2];
							range[0]=fmin;
							range[1]=fmax;
						}
						final int[] range;
//						@Override
//						public int[] call() 
//						{
//							run();
//							return range;
//						}
						public void run()
						{
							for (int f=range[0]; f<range[1] && !isCancelled(); f++)
							{
								if (Thread.currentThread().isInterrupted())
								{
									return;
								}
								computeFamily(f);
								boolean fits_in_queue = families_done.offer(f);
								assert (fits_in_queue);
									
							} 
						}
					} // class ComputingTask
					
					int num_tasks = Math.min(1+(nF-1)/THREAD_UNIT_TASK, 256);
					int task_size = 1+(nF-1)/num_tasks;

//					ExecutorService executor = 
//							Executors.newWorkStealingPool(THREAD_PARALLELISM); //   Executors.newCachedThreadPool();
					ExecutorCompletionService<int[]> completion = new ExecutorCompletionService<>(executor);
					List<ComputingTask> computing_tasks = new ArrayList<>();
					int fmin=0;
					while (fmin<nF)
					{
						int fmax = Integer.min(fmin+task_size,nF);
						ComputingTask T = new ComputingTask(fmin, fmax);
						computing_tasks.add(T);
						completion.submit(T, T.range);
						fmin = fmax;
					}
					try
					{
//						if (DEBUG_THREADS)
//							System.out.println("#***HV.cA/HW.dIB "+Thread.currentThread()+"\tstarting nF="+nF);
						for (int iter=0; iter<nF; iter++)
						{
							int f = families_done.take();
							publish(f);
							if (DEBUG_THREADS)
								System.out.println("#***HV.cA/HW.dIB "+Thread.currentThread()+"\tpublish "+f+"\titer "+iter+"\t/"+nF);
						}
						for (int t=0; t<computing_tasks.size(); t++)
						{
//							if (DEBUG_THREADS)
//								System.out.println("#***HV.cA/HW.dIB "+Thread.currentThread()+"\twaiting task "+t+"/"+computing_tasks.size());
							Future<int[]> this_one_done = completion.take();
							
							int[] range = this_one_done.get();
//							if (DEBUG_THREADS)
//								System.out.println("#***HV.cA/HW.dIB "+Thread.currentThread()+"\tgot "+range[0]+".."+range[1]);
							
//							for (int f=range[0]; f<range[1]; f++)
//								publish(f);
						}
					} catch (InterruptedException E)
					{
						// dunno about this 
						System.out.println("#**HV.cA/HW.dIB "+E);
						this.cancel(true);
//						throw new RuntimeException(E);
					} 
					catch (ExecutionException EE)
					{
						// or this
						System.out.println("#**HV.cA/HW.dIB "+EE);
						this.cancel(true);
//						throw new RuntimeException(EE);
					} catch (Throwable t)
					{
						System.out.println("#**HV.cA/HW.dIB "+t);
						this.cancel(true);
					}
				} else
				{
					for (int f=0; f<nF; f++)
					{
						if (Thread.currentThread().isInterrupted())
						{
							System.out.println("#**HV.cA/HW.dIB interrupted @ f="+f+"\t// "+Thread.currentThread());							
							return (Void)null;
						}
						computeFamily(f);
						publish(f);
					}
				}
				return (Void)null;
			}
			
			@Override
			public void done()
			{
				if (isCancelled())
				{
					computation_progress.setString("Canceled prematurely ");
				} else
				{
					computation_progress.setString("Done "+nF);
					computation_progress.setVisible(false);
				}
				// recolor
				for (int col=table_model.firstHistoryColumn(); col<table_model.getColumnCount(); col++)
				{
					colorHistoryColumn(col);
//					HistoryModel.Column<?> column_data = table_model.getHistoryColumn(col);
//					Color max_color = column_data.getColor();
//					if (max_color!=null)
//					{
//						double min = column_data.getMinimum();
//						if (min<1.0) min=0.0;
//						double max = Double.max(1.0,column_data.getMaximum());
//						
//						ColoredValueRenderer renderer
//						= new ColoredValueRenderer(Color.WHITE,max_color,min,max);
//						table_scroll.setColumnRenderer(col, renderer);
//					}
				}

				computation_cancel.setEnabled(false);
				//computation_cancel.setVisible(false);
				//tree_control.removeControlAt(Integer.parseInt(computation_cancel.getActionCommand()));
				computeDone();

				tree_control.repaint();
				table_scroll.repaint();
			}
			
			@Override
			protected void process(List<Integer> families)
			{
				// find range of families updated 
				int first = families.get(0);
				int last = first;
				for (int j=1; j<families.size(); j++)
				{
					int f = families.get(j);
					first = Integer.min(first, f);
					last = Integer.max(last, f);
				}
				num_families_done += families.size();
				int advance = (int)(100.0*num_families_done/nF);
				this.setProgress(advance);
				computation_progress.setString("Computing "+num_families_done+"/"+nF);
				computation_progress.setValue(num_families_done);
				table_model.fireTableRowsUpdated(first, last);
			}
		}
		history_task = new HistoryWorker();
		computation_progress.setMaximum(nF);
		computation_progress.setIndeterminate(false);
		computation_progress.setValue(0);
		computation_progress.setString("Computing 0/"+nF);
		
		history_task.execute();
		computation_cancel.setEnabled(true);
	}
	
	/**
	 * Resets min-max coloring in a column  
	 * @param col
	 */
	protected final void colorHistoryColumn(int col)
	{
		HistoryModel.Column<?> column_data = table_model.getHistoryColumn(col);
		Color max_color = column_data.getColor();
		if (max_color != null)
		{
			double min = column_data.getMinimum();
			double max = column_data.getMaximum();
			if (min==max)
			{
				table_scroll.setColumnRenderer(col, null);
			} else
			{
				if (min<1.0) min=0.0;
				max = Double.max(1.0,max);
				ColoredValueRenderer renderer
				= new ColoredValueRenderer(Color.WHITE,max_color,min,max);
				table_scroll.setColumnRenderer(col, renderer);
			}
		}
	}

	@Override
	public int[] getSelectedFamilies()
	{
		return table_scroll.getSelectedModelRows();
	}

	@Override
	public void saveData(File f) throws IOException 
	{
		// TODO: progress monitor?
		
		
        PrintStream PS = new PrintStream(f);
        PS.println(CommandLine.getStandardHeader(getClass()));
        PS.println(CommandLine.getStandardRuntimeInfo());
        
//        System.out.println("#**HV.sD "+f+"\tstart");
        
        table_model.printTable(PS);
//        PS.println(table_model.getTextTable());
        
        if (PS.checkError()) // also flushes
        {
            PS.close();
            throw new IOException("Cannot write the table.");
        }
        PS.close();
//        System.out.println("#**HV.sD "+f+"\tdone");
	}
	
	
}
