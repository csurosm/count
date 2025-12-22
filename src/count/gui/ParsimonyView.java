package count.gui;
/*
 * Copyright 2023 Mikl&oacute;s Cs&#369;r&ouml;s.
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


import count.ds.IndexedTree;

import java.awt.Color;


import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JEditorPane;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JSpinner;
import javax.swing.JTabbedPane;
import javax.swing.SpinnerNumberModel;

import count.ds.AnnotatedTable;
import count.io.DataFile;
import count.model.Parsimony;
/**
 * {@link HistoryView} with ancestral reconstruction by numerical (Wagner) parsimony.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */
public class ParsimonyView extends HistoryView
{
	private static final double SPINNER_PRECISION = 0.001;
	
	public ParsimonyView(DataFile<? extends IndexedTree> tree_data, 
			DataFile<AnnotatedTable> table_data)
	{
		super(tree_data.getContent(), table_data, true, true);
		this.tree_data = tree_data;
		this.factory = new Parsimony(tree_data.getContent(), table_data.getContent());
		initComponents();
	}
	
	private final DataFile<? extends IndexedTree> tree_data;
	
	private final Parsimony factory;
	
	private double gain_penalty = Parsimony.DEFAULT_GAIN_PENALTY;
	private double duplication_penalty = Parsimony.DEFAULT_DUPLICATION_PENALTY;
	private double loss_penalty = Parsimony.DEFAULT_LOSS_PENALTY;
	
	private JSpinner gain_spinner;
	private JSpinner duplication_spinner;
	private JSpinner loss_spinner;
	
	private void initComponents()
	{
		boolean is_binary_table = table_data.getContent().isBinaryTable();
		
		gain_spinner = createSpinner(gain_penalty,"Penalty for gain 0->1", FAMILY_GAIN_COLOR.brighter() );
		duplication_spinner = createSpinner(duplication_penalty, "Per-copy penalty factor for increase (births) 1->x", FAMILY_EXPAND_COLOR.brighter());
		
		
		loss_spinner = createSpinner(is_binary_table?1.0:loss_penalty,"Penalty for loss 1->0 (per-copy decrease/death penalty is 1.0)", FAMILY_LOSE_COLOR.brighter());

		gain_spinner.addChangeListener(e->
						{
							double pty = ((Number)gain_spinner.getValue()).doubleValue();
							if (gain_penalty != pty)
							{
								gain_penalty = pty;
								computeAll();
							}
						});
		duplication_spinner.addChangeListener(e->
						{
							double pty = ((Number)duplication_spinner.getValue()).doubleValue();
							if (duplication_penalty != pty)
							{
								duplication_penalty = pty;
								computeAll();
							}
						});
		loss_spinner.addChangeListener(e->
						{
							double pty = ((Number)loss_spinner.getValue()).doubleValue();
							if (loss_penalty != pty)
							{
								loss_penalty = pty;
								computeAll();
							}
						});
		
		int h = 12;
        JLabel gain_label = new JLabel("Gain");
        gain_label.setLabelFor(gain_spinner);
//        gain_label.setBackground(((JSpinner.DefaultEditor)gain_spinner.getEditor()).getTextField().getBackground());
//        gain_label.setOpaque(true);
        JLabel duplication_label = new JLabel("Expansion");
        duplication_label.setLabelFor(duplication_spinner);
//        duplication_label.setBackground(((JSpinner.DefaultEditor)duplication_spinner.getEditor()).getTextField().getBackground());
//        duplication_label.setOpaque(true);
        JLabel loss_label = new JLabel("Loss");
        loss_label.setLabelFor(loss_spinner);
//        loss_label.setBackground(((JSpinner.DefaultEditor)loss_spinner.getEditor()).getTextField().getBackground());
//        loss_label.setOpaque(true);

		Box control_box = new Box(BoxLayout.LINE_AXIS);
		
		control_box.add(gain_label);
		control_box.add(gain_spinner);
		if (is_binary_table)
		{
			loss_penalty =1.0;
		} else
		{
			control_box.add(duplication_label);
			control_box.add(duplication_spinner);
			control_box.add(loss_label);
			control_box.add(loss_spinner);
		}
		
		control_box.setBorder(BorderFactory.createTitledBorder("Scoring"));//  (   BorderFactory.createEtchedBorder());
		
//		control_box.setMinimumSize(new Dimension(control_box.getComponentCount()*48,h));
		computation_cancel.setVisible(false);
        getTreeControl().addControlAt(2,control_box);
        
        JEditorPane table_explain = createLegend();
        JTabbedPane selected_rows_display = this.getSelectedRowsDisplay();
		selected_rows_display.addTab("Legend", table_explain);
	}
	
	private JSpinner createSpinner(double initial_value, String tooltip_text, Color background_color)
	{
        SpinnerNumberModel penalty_model = new SpinnerNumberModel(initial_value, 0.0, Double.POSITIVE_INFINITY, SPINNER_PRECISION);
        JSpinner spinner = new JSpinner(penalty_model);
        spinner.setToolTipText(tooltip_text);
        final JFormattedTextField gps_text = ((JSpinner.DefaultEditor)spinner.getEditor()).getTextField();
        spinner.addChangeListener(e->gps_text.transferFocus()); // or else the cursor keeps blinking in there
//        gps_text.setBackground(background_color);
        
        
        ((JSpinner.DefaultEditor)spinner.getEditor()).getTextField().setColumns(5);
        
        return spinner;
	}
	
	
	private JEditorPane createLegend() {
		String table_legend = "<p>The <b>Tree</b> panel is a visual display of "
				+ "the reconstruction at each node for the selected families in the top table."
				+ "The <b>Table</b> panel gives the same values; you can sort by columns to see distributions "
				+ "across lineages. Node selections are synchronized between <b>Tree</b> and <b>Table</b>."
				+ "The numerical parsimony "
				+ "reconstruction penalizes each copy change <var>p</var>→<var>c</var> from parent with <var>p</var> copies "
				+ "to a child with <var>c</var> copies; we minimize the score of a reconstruction which is the edge-wise penalty sum on phylogeny with the "
				+ "inferred ancestral copies. The penalty function <var>P</var>(<var>p</var>,<var>c</var>) is a generalization of Wagner parsimony: "
				+ "sum of <i>step</i> penalties for each increase <var>i</var>→(<var>i</var>+1) or each decrease <var>i</var>→(<var>i</var>-1)"
				+ "from <var>p</var> to <var>c</var>. Steps from 0 or to 0 are <b>gain</b> and <b>loss</b>; "
				+ "otherwise it is an <b>expansion</b>.</p> "
				+ "<ul>"
				+ "<li><var>P</var>(<var>p</var>, <var>p</var>)=0</li>"
				+ "<li><var>P</var>(0, <var>c</var>) = <b>gain</b>+<b>expansion</b>*(<var>c</var>-1) for 1&le;<var>c</var></li>"
				+ "<li><var>P</var>(<var>p</var>, <var>p</var>+<var>inc</var>) = <b>expansion</b>*(<var>inc</var>) for 1&le;<var>p</var>, 1&le;<var>inc</var></li>"
				+ "<li><var>P</var>(<var>p</var>, 0) = <b>loss</b>+(<var>p</var>-1) for 1&le;<var>p</var></li>"
				+ "<li><var>P</var>(<var>p</var>, <var>p</var>-<var>dec</var>) = <var>dec</var> for 1&le;<var>p</var>, 1&le;<var>dec</var></li>"
				+"</ul>"
				+"<p>Use the spinners to set the penalties, or enter a value. "
				+ "Use a small deviation  from integer values (eg 2.0001 instead of 2 for gain)"
				+ " to resolve the ties gently between optimal solutions. (The table has columns for minimum-copy and maximum-copy "
				+ "reconstructions that are differnt if there are ties.)<p>"
				+ "<p>The following events are counted: copies, birth/death (copy change), family presence (at least 1 copy), family "
				+ "loss and gain (copy change to and from 0), and family contraction and expansion (copy change to and from 1).</p>"
				;
		
		JEditorPane createLegend = new JEditorPane("text/html", table_legend);
		createLegend.setEditable(false);		
		
		return createLegend;
	}
	
	
	public double getLossPenalty() {return loss_penalty;}
	public double getGainPenalty() {return gain_penalty;}
	public double getDuplicationPenalty() { return duplication_penalty;}
	
	public void setPenalties(double gain_pty, double loss_pty, double dup_pty)
	{
		this.gain_penalty = gain_pty;
		this.loss_penalty = loss_pty;
		this.duplication_penalty = dup_pty;
	}
	
	
	
	/*
	 * Profiles
	 */
	
	
	private Parsimony.Profile[] profiles=null;
	
	@Override
	public void computeAll()
	{
		
		if (profiles==null)
		{
			int nF = table_data.getContent().getFamilyCount();
			profiles = new Parsimony.Profile[nF];
		}
		
		gain_spinner.setEnabled(false);
		duplication_spinner.setEnabled(false);
		loss_spinner.setEnabled(false);
		
		
		super.computeAll();
	}
	
	@Override
	protected void computeDone()
	{
		if (!history_task.isCancelled())
		{
			gain_spinner.setEnabled(true);
			duplication_spinner.setEnabled(true);
			loss_spinner.setEnabled(true);
			
			double score = column_family_score.getSum();
			int gains = (int) column_family_gains.getSum();
			int losses = (int) column_family_losses.getSum();
			
			StringBuilder info = new StringBuilder("<html><p>");
			info.append("Score <b>").append(score).append("</b>");
			info.append("; ").append(gains).append(" gains")
			.append(", ").append(losses).append(" losses")
			.append(" (").append(gains+losses).append(" events).");
			info.append("</p></html>");
			this.setScoringInfo(info.toString());
			
		}
		colorHistoryColumn(table_model.firstHistoryColumn());
		
	}
	
	@Override
	protected void computeHistoriesPrepare()
	{
		super.computeHistoriesPrepare();
		factory.setPenalties(gain_penalty, loss_penalty, duplication_penalty);
		int nF = table_data.getContent().getFamilyCount();
		for (int f=0; f<nF; f++)
		{
			if (profiles[f]==null)
			{
				profiles[f] = factory.getProfile(f);
			}
		}
	}
	
	
	@Override 
	protected void computeFamily(int f)
	{
//		if (f<20 || f%1000==0)
//			System.out.println("#**PV.cF "+f);

		Parsimony.Profile P = profiles[f];
		int[] min_copies = P.computeSankoff(false);
		column_family_score.setValue(f, P.getSankoffScore());

		int[] max_copies = P.computeSankoff(true);
		
		IndexedTree phylo = tree_data.getContent();

		int num_nodes = phylo.getNumNodes();
		for (int node=0; node<num_nodes; node++)
		{
			table_member_count.get(node).setValue(f, min_copies[node]);
			table_member_max.get(node).setValue(f, max_copies[node]);
		}
	}
	
	
	
    @Override
    public String toString()
    {
        return "Pars @ "+DataFile.chopCommonFileExtension(tree_data.getFile().getName());
    }
	
	
	
}
