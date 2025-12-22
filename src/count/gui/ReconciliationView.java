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

import javax.swing.JCheckBox;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.Reconciliation;
import count.gui.HistoryModel.LineageColumns;
import count.gui.kit.RoundedDouble;
import count.io.DataFile;

/**
 * 
 * Display for reconciliation statistics.
 */
public class ReconciliationView extends HistoryView {
	public ReconciliationView(DataFile<Reconciliation<?>> rec_data) {
		super(rec_data.getContent().getTree()
				, new DataFile<AnnotatedTable>(rec_data.getContent().toTable(),
						rec_data.getFile())
				, true, false, rec_data.getContent().hasIntegers());
		this.rec_data = rec_data;
		this.rec = rec_data.getContent();
		if (rec.hasIntegers())
			throw new UnsupportedOperationException("Not yet implemented for non-probabilistic reconc.");
	}
	
	private final DataFile<Reconciliation<?>> rec_data;
	private final Reconciliation<?> rec;
	private AnnotatedTablePanel table_panel=null;
	
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
	
	private final static String MEMBER_TRANSFER = "+>";
	private final static String MEMBER_DUPLICATION = "+|";
	
	private HistoryModel.Column<? extends Number> column_family_transfer;
	private HistoryModel.Column<? extends Number> column_family_duplication;
	

	
	@Override
	protected void initTableModel(){
		IndexedTree phylo = table_model.getTree();
		int num_nodes = phylo.getNumNodes();
		table_family_present = new LineageColumns(FAMILY_PRESENT); // 
		table_family_gain = new LineageColumns(FAMILY_GAIN); // gain
		table_family_lose = new LineageColumns(FAMILY_LOSE); // loss

		table_member_count = new LineageColumns(MEMBER_COUNT);
		table_member_birth = new LineageColumns(MEMBER_BIRTH);
		table_member_death = new LineageColumns(MEMBER_DEATH);
		
		table_family_expand = new LineageColumns(MEMBER_TRANSFER); // repurposed expansion 
		table_family_contract = new LineageColumns(MEMBER_DUPLICATION); // repurposed contraction

		// copy count columns
		String eS = (isSimulated?"Number of ":"Expected number of ");
		table_member_count.addAll(table_model.newDoubleColumns(MEMBER_COUNT, eS+"copies", false));
		table_member_birth.addAll(table_model.newDoubleColumns(MEMBER_BIRTH, eS+"gained copies", true));
		table_member_death.addAll(table_model.newDoubleColumns(MEMBER_DEATH, eS+"lost copies", true));
		
		table_family_expand.addAll(table_model.newDoubleColumns(MEMBER_TRANSFER, eS+"transfers+donorless originations", true));
		table_family_contract.addAll(table_model.newDoubleColumns(MEMBER_DUPLICATION, eS+"duplications", true));
		
		// family change counts
		String pS = (isSimulated?"Indicator ":"Probability ");
		table_family_present.addAll(table_model.newDoubleColumns(FAMILY_PRESENT, pS+"of at least one copy", false));
		table_family_gain.addAll(table_model.newDoubleColumns(FAMILY_GAIN, pS+"that the family was acquired", true));
		table_family_lose.addAll(table_model.newDoubleColumns(FAMILY_LOSE, pS+"that the family was completely lost", true));
	
		column_family_score = table_model.newDoubleColumn("Events", "Number of copy births and deaths");

		column_family_gains = table_model.newSumColumn(table_family_gain, HEADER_GAINS, (isSimulated?"Resolved":"Expected")+" number of lineages in which the family was gained.");
		column_family_losses = table_model.newSumColumn(table_family_lose, HEADER_LOSSES, (isSimulated?"Rewsolved":"Expected")+" number of lineages in which the family was lost.");			
	
		table_model.add(column_family_score, FAMILY_LL_COLOR);
		table_model.add(column_family_gains, FAMILY_GAIN_COLOR);
		table_model.add(column_family_losses, FAMILY_LOSE_COLOR);

		column_family_transfer = table_model.newSumColumn(table_family_expand, "Transfers", "Number of copy transfers and originations.");
		column_family_duplication = table_model.newSumColumn(table_family_contract, "Duplications", "Number of duplications");
		
		
		LineageColumns table_member_change=new LineageColumns(MEMBER_CHANGE); //<>(num_nodes);
		if (has_only_integers)
		{
			// resolved reconciliations : difference is this node.count-parent.count
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
		
		table_model.add(column_family_transfer, FAMILY_GAIN_COLOR);
		table_model.add(column_family_duplication, FAMILY_EXPAND_COLOR);
				
		table_model.addAll(table_member_count, FAMILY_GAIN_COLOR);
		table_model.addAll(table_member_change, null);
		
		table_model.addAll(table_member_birth, FAMILY_GAIN_COLOR);
		table_model.addAll(table_member_death, FAMILY_LOSE_COLOR);
		
		table_model.addAll(table_family_present, FAMILY_PRESENT_COLOR);
		
		table_model.addAll(table_family_gain, FAMILY_GAIN_COLOR);
		table_model.addAll(table_family_lose, FAMILY_LOSE_COLOR);
		table_model.addAll(table_family_expand, FAMILY_GAIN_COLOR);
		table_model.addAll(table_family_contract, FAMILY_EXPAND_COLOR);
	}	
	
	@Override
	protected HistoryTreePanel createLineageModel(HistoryTreePanel P) {


		lineage_member_count = P.newColumns("Copies "+MEMBER_COUNT+"", table_member_count);
		JCheckBox count_cb = P.showNodeStatistics("Copies", lineage_member_count);
		//count_cb.setForeground(FAMILY_PRESENT_COLOR);
		lineage_member_birth = P.newColumns("Births "+MEMBER_BIRTH+"", table_member_birth);
		lineage_member_death = P.newColumns("Deaths "+MEMBER_DEATH+"", table_member_death);
		P.showChangeStatistics("Change", lineage_member_birth, lineage_member_death);;

		lineage_family_expand = P.newColumns("Transfers "+MEMBER_TRANSFER+"", table_family_expand);
		lineage_family_contract = P.newColumns("Duplications "+MEMBER_DUPLICATION+"", table_family_contract);
		P.showNodeStatistics("Transfers", lineage_family_expand);
		P.showNodeStatistics("Duplications", lineage_family_contract);
		
		lineage_family_present = P.newColumns("Families "+FAMILY_PRESENT+"", table_family_present);
		P.showNodeStatistics("Families", lineage_family_present);
		lineage_family_gain = P.newColumns("Gains "+FAMILY_GAIN+"", table_family_gain);
		lineage_family_lose = P.newColumns("Losses "+FAMILY_LOSE+"", table_family_lose);
		P.showChangeStatistics("Loss/gain", lineage_family_gain, lineage_family_lose);
		
		return P;
	}
	
	
	@Override
	protected void computeFamily(int f){
		Reconciliation<? extends Number>.HistoryProfile hp = rec.getFamily(f);
		if (rec.hasIntegers()) {
		} else {
			IndexedTree phylo = table_model.getTree();
			int num_nodes = phylo.getNumNodes();
			
			double nevents = 0.0;
			for (int u=0; u<num_nodes; u++){
				double ncopy = hp.dNcopy(u);
				double ndup = hp.dNdup(u);
				double ntra = hp.dNtransfer(u);
				double nloss = hp.dNloss(u);
				double present = hp.dPresent(u);
				double scopy = hp.dScopy(u);
				double birth = ndup+ntra;
				table_member_birth.get(u).setValue(f, birth);
				table_member_death.get(u).setValue(f, nloss);
				table_member_count.get(u).setValue(f, ncopy);
				table_family_present.get(u).setValue(f, present);
				
				// upper bound for multi-copy family
				double multi = Double.min(ncopy-present, 1.0);
				//table_family_multi.get(u).setValue(f, multi);
				
				// scopy == ncopy-birth
				// upper bound for inheritance
				double pinherit = Double.min(ncopy-birth, present);
				double pgain = present-pinherit;
				double pp = phylo.isRoot(u)?0.0:hp.dPresent(phylo.getParent(u));
				double ploss = Double.max(0.0, pp-pinherit);
				table_family_gain.get(u).setValue(f, pgain);
				table_family_lose.get(u).setValue(f, ploss);
				table_family_expand.get(u).setValue(f, ntra);
				table_family_contract.get(u).setValue(f, ndup);
				
				nevents += ndup+ntra+nloss;
			}
			column_family_score.setValue(f, nevents);
		}
	}
	
	@Override
	protected void computeDone() {
		double nevents = column_family_score.getSum();
		double ntra = column_family_transfer.getSum();
		double ndup = column_family_duplication.getSum();
		double nloss= nevents-ntra-ndup;
		
		StringBuilder info = new StringBuilder("<html><p>Total events ");
		info.append("<b>").append(String.format("%.1f", nevents)).append("</b>");
		info.append("; ")
		.append(String.format("%.1f", ntra)).append(" transfers+originations, ")
		.append(String.format("%.1f", ndup)).append(" duplications, ")
		.append(String.format("%.1f", nloss)).append(" losses");
		info.append(".</p></html>");
		setScoringInfo(info.toString());
	}
	
	
	@Override
	public String toString(){
		return "Reconciliations"; 
	}
}
