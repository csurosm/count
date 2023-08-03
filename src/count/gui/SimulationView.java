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

import java.awt.Component;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.event.TableModelEvent;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.io.DataFile;
import count.model.MixedRateModel;
import count.model.Posteriors.FamilyEvent;
import count.model.SimulatedEvolution;
import count.model.SimulatedEvolution.Table.ObservedProfile;

/**
 * {@link HistoryView} for a simulated data set. T
 * The simulation includes the ancestral and observed copy numbers. The 
 * table for the observed copy numbers can be retrieved by {@link #asTable()}. 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */
public class SimulationView extends HistoryView
{
	public SimulationView(DataFile<? extends MixedRateModel> rates_data, int num_rows, int min_observed)
	{
		super(rates_data.getContent().getClassModel(0).getTree(),
				new DataFile<AnnotatedTable>(SimulatedEvolution.table(rates_data.getContent(), num_rows, min_observed),
						new File((File)null, "Sim @ "+DataFile.chopFileExtension(rates_data.getFile().getName())+"*"+num_rows)),
				true, false);
		this.rates_data = rates_data;
		// let's recover the table from super's instantiation
		this.sim_table = (SimulatedEvolution.Table)table_data.getContent();
		initDataStructures();
		initComponents();
	}

	public SimulationView(DataFile<? extends MixedRateModel> rates_data, long rnd_seed, int num_rows, int min_observed)
	{
		super(rates_data.getContent().getClassModel(0).getTree(),
				new DataFile<AnnotatedTable>(SimulatedEvolution.table(rates_data.getContent(), rnd_seed, num_rows, min_observed),
						new File((File)null, "Sim @ "+DataFile.chopFileExtension(rates_data.getFile().getName())+"*"+num_rows)),
				true, false);
		this.rates_data = rates_data;
		// let's recover the table from super's instantiation
		this.sim_table = (SimulatedEvolution.Table)table_data.getContent();
		initDataStructures();
		initComponents();
	}
	
	
	
	/**
	 * Input rate model
	 */
	private final DataFile<? extends MixedRateModel> rates_data;
	/**
	 * Same as the content of table_data. 
	 */
	private final SimulatedEvolution.Table sim_table; 
	
	private boolean want_unobserved=true;
	
	private List<HistoryModel.Column<Integer>> table_family_class;
	
	private JCheckBox obs_cb ;
	
	private AnnotatedTablePanel table_panel=null;

	
	private void initDataStructures()
	{
		table_family_class = new ArrayList<>();
		for (int ci=0; ci<sim_table.getNumClasses(); ci++)
		{
			String class_name = "Class "+Integer.toString(ci);
			HistoryModel.Column<Integer> cnt = table_model.newIntColumn(class_name, "Families in "+class_name.toLowerCase());
			table_model.add(cnt, FAMILY_PRESENT_COLOR);
			table_family_class.add(cnt);
		}
	}
	
	private void initComponents()
	{
        getTreeControl().addControlAt(2, createObservedCopiesControl());
	}
	
	private Component createObservedCopiesControl()
	{
		int minobs = sim_table.getMinimumObserved();
//		JLabel minobs_info = new JLabel(minobs+"≤copies observed");
//		Font info_font = getTableScroll().getDataTable().getFont().deriveFont(Font.ITALIC);
//		minobs_info.setFont(info_font);
//		return minobs_info;
		obs_cb = new JCheckBox("With unobserved (<"+minobs+")");
		obs_cb.setSelected(want_unobserved);
		
		obs_cb.addItemListener(e->
					{
						want_unobserved = obs_cb.isSelected();
						computeAll();
					});
		
		obs_cb.setVisible(minobs>0);
		
		return obs_cb;
	}
	
	/**
	 * The random seed used to initialize the 
	 * {@link java.util.Random} generator at instantiation. 
	 */
	public long getRandomSeed() { return sim_table.getRandomSeed();}
	
	
	public int getMinimumObserved() { return sim_table.getMinimumObserved();}
	
	public AnnotatedTablePanel asTable()
	{
		if (table_panel==null)
		{
			table_panel = new AnnotatedTablePanel(table_data);
			table_model.addTableModelListener(e->
				{
					if (e.getType()==e.UPDATE)
					{
						AnnotatedTableModel M = table_panel.getTableModel();
						M.fireTableRowsUpdated(e.getFirstRow(), e.getLastRow());
					}
				});
		}
		return table_panel;
	}
	
	@Override
	protected  void computeAll()
	{
		obs_cb.setEnabled(false);
		super.computeAll();		
	}

	@Override
	protected void computeDone()
	{
//		computation_cancel.setVisible(false);
		super.computeDone(); // disables computation_cancel button
		obs_cb.setEnabled(true);
	}
	
	@Override
	protected void computeFamily(int f)
	{
		ObservedProfile P = sim_table.getObservedProfile(f);
		table_family_score.setValue(f, P.getHistoryEventCount(want_unobserved));
		
		IndexedTree phylo = table_model.getTree() ;
		int num_nodes = phylo.getNumNodes();
		for (int node=0; node<num_nodes; node++)
		{
			double N = P.getNodeSurvivalCount(node, want_unobserved);
			double S = P.getEdgeSurvivalCount(node, want_unobserved);
			
			table_member_count.get(node).setValue(f, N);
			table_member_birth.get(node).setValue(f, N-S);
			double death = 0.0;
			
			if (!phylo.isRoot(node))
			{
				int parent = phylo.getParent(node);
				death = P.getNodeSurvivalCount(parent, want_unobserved)-S;
			}
			table_member_death.get(node).setValue(f, death);
			
			table_family_present.get(node).setValue(f, P.getFamilyPresent(node, want_unobserved));
			table_family_multi.get(node).setValue(f, P.getFamilyMulti(node, want_unobserved));
			int[] events = P.getFamilyEvents(node, want_unobserved);
			table_family_gain.get(node).setValue(f, events[FamilyEvent.GAIN.ordinal()]);
			table_family_lose.get(node).setValue(f, events[FamilyEvent.LOSS.ordinal()]);
			table_family_expand.get(node).setValue(f, events[FamilyEvent.EXPAND.ordinal()]);
			table_family_contract.get(node).setValue(f, events[FamilyEvent.CONTRACT.ordinal()]);
			
			for (int c=0;c<sim_table.getNumClasses(); c++)
			{
				HistoryModel.Column<Integer> class_cnt = table_family_class.get(c);
				class_cnt.setValue(f, P.getFamilyClass(c, want_unobserved));
			}
		}
	}
	
	public String toString()
	{
		return DataFile.chopFileExtension(table_data.getFile().getName());
	}
}
