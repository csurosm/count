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
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JSpinner;
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
	}
	
	private JSpinner createSpinner(double initial_value, String tooltip_text, Color background_color)
	{
        SpinnerNumberModel penalty_model = new SpinnerNumberModel(initial_value, 0.0, Double.POSITIVE_INFINITY, 0.1);
        JSpinner spinner = new JSpinner(penalty_model);
        spinner.setToolTipText(tooltip_text);
        final JFormattedTextField gps_text = ((JSpinner.DefaultEditor)spinner.getEditor()).getTextField();
        spinner.addChangeListener(e->gps_text.transferFocus()); // or else the cursor keeps blinking in there
//        gps_text.setBackground(background_color);
        
        
        ((JSpinner.DefaultEditor)spinner.getEditor()).getTextField().setColumns(3);
        
        return spinner;
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
