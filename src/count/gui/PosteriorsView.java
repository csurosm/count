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
import java.awt.Component;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JEditorPane;
import javax.swing.JRadioButton;
import javax.swing.JTabbedPane;
import javax.swing.event.TableModelEvent;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.gui.kit.ColoredValueRenderer;
import count.gui.kit.RoundedDouble;
import count.gui.kit.TableScroll;
import count.io.DataFile;
import count.matek.Logarithms;
import count.model.MixedRateModel;
import count.model.MixedRatePosteriors;

import count.model.Posteriors.FamilyEvent;

/**
 * {@link HistoryView} with ancestral reconstruction by posterior probabilities.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */

public class PosteriorsView extends HistoryView 
{
	/** 
	 * Phred scale for profile rarity.
	 */
	private static final double EVIDENCE_SCALE = -10.0/Math.log(10.0);
	
	private static final boolean FAMILY_PRESENCE_BY_SURVIVAL = true;
	
	public PosteriorsView(DataFile<? extends MixedRateModel> rates_data, DataFile<AnnotatedTable> table_data)
	{
		super(rates_data.getContent().getClassModel(0).getTree(),
				table_data,
				true, false);
		this.rates_data = rates_data;
		this.posteriors = new MixedRatePosteriors(rates_data.getContent(), table_data.getContent());
		initDataStructures();
		initComponents();
	}
	
	private final DataFile<? extends MixedRateModel> rates_data;
	private final MixedRatePosteriors posteriors;
	
	private List<HistoryModel.Column<RoundedDouble>> table_family_class;
	
	private List<MinCopies> min_copiesRB;
	
//	private IndexedTree getTree()
//	{
//		return table_model.getTree();
//	}

	public int getMinimumObserved()
	{
		return posteriors.getMinimumObserved();
	}
	
	
	private void initDataStructures()
	{
		table_family_class = new ArrayList<>();

		int num_classes = posteriors.getClassCount();
		
		if (1<num_classes)
		{
			for (int c=0; c<num_classes; c++)
			{
				String class_name = "Class "+Integer.toString(c);
				HistoryModel.Column<RoundedDouble> probs = table_model.newDoubleColumn(class_name, "Posterior probability for "+class_name.toLowerCase());
				int class_column  = table_model.add(probs, FAMILY_PRESENT_COLOR);
	//			System.out.println("#**PV.iDS class "+c+"\tcol "+class_column
	//					+", after "+table_model.getColumnName(table_model.firstHistoryColumn()+class_column-1));
				
				table_family_class.add(probs);
			}
			table_model.fireTableStructureChanged();
		}
	}	
	
	private class MinCopies extends JRadioButton
    {
    	MinCopies(int n)
    	{
    		super(Integer.toString(n));
    		if (n==0)
    			this.setToolTipText("All, even empty profiles are observed; no correction to family statistics");
    		else
    			this.setToolTipText("Profiles with at least "
    					+Integer.toString(n)
    					+" member"
    					+(n>1?"s":"")
    					+"are observed; correction for unobserved empty"
    					+(n>1?" and singleton":"")
    					+"  profiles in family statistics");
    		this.addActionListener(click->
    		{
    			TableScroll<HistoryModel> table_scroll = (TableScroll<HistoryModel>) getTopComponent();			
    			int[] selected_rows = table_scroll.getSelectedModelRows();
    			posteriors.setMinimumObserved(n);
    			setUnobservedCorrections();
    			Color max_color = table_family_score.getColor();
    			double min = table_family_score.getMinimum();
    			double max = table_family_score.getMaximum();
				ColoredValueRenderer renderer
				= new ColoredValueRenderer(Color.WHITE,max_color,min,max);
				table_scroll.setColumnRenderer(table_model.firstHistoryColumn(), renderer);
    			table_scroll.setSelectedModelRows(selected_rows);
    			table_scroll.repaint();
    		});
    	}
    }
	
	/**
	 * A container with radio buttons for min.copies 0,1,2.
	 * @return
	 */
	private Component createMincopiesControl()
	{
        Box minCopiesB = new Box(BoxLayout.LINE_AXIS);
        minCopiesB.setBorder(BorderFactory.createTitledBorder("Family correction"));
        
        int data_min_copies = table_data.getContent().minCopies();
        MinCopies minCopies0 = new MinCopies(0);
        MinCopies minCopies1 = new MinCopies(1);
        MinCopies minCopies2 = new MinCopies(2);
        
        min_copiesRB = new ArrayList<>();
        min_copiesRB.add(minCopies0);
        min_copiesRB.add(minCopies1);
        min_copiesRB.add(minCopies2);
        
        
        ButtonGroup minCopiesG = new ButtonGroup();
        minCopiesG.add(minCopies0);
        minCopiesG.add(minCopies1);
        minCopiesG.add(minCopies2);
        minCopies1.setEnabled(data_min_copies>0);
        minCopies2.setEnabled(data_min_copies>1);

        minCopiesB.add(minCopies0);
        minCopiesB.add(minCopies1);
        minCopiesB.add(minCopies2);
        if (data_min_copies==0)
        	minCopies0.setSelected(true);
        else if (data_min_copies==1)
        	minCopies1.setSelected(true);
        else
        {
        	assert (data_min_copies>=2);
        	minCopies2.setSelected(true);
        }
        
        minCopies0.addActionListener(click->
		{
//			System.out.println("#**PP.cMC/AL0 "+click);
			TableScroll<HistoryModel> table_scroll = (TableScroll<HistoryModel>) getTopComponent();			
			int[] selected_rows = table_scroll.getSelectedModelRows();
			posteriors.setMinimumObserved(0);
			setUnobservedCorrections();
//			table_model.fireTableDataChanged();
			table_scroll.setSelectedModelRows(selected_rows);
		});
        minCopies1.addActionListener(click->
		{
//			System.out.println("#**PP.cMC/AL1 "+click);
			TableScroll<HistoryModel> table_scroll = (TableScroll<HistoryModel>) getTopComponent();			
			int[] selected_rows = table_scroll.getSelectedModelRows();
			posteriors.setMinimumObserved(1);
			setUnobservedCorrections();
//			table_model.fireTableDataChanged();
			table_scroll.setSelectedModelRows(selected_rows);
		});
        minCopies2.addActionListener(click->
		{
			TableScroll<HistoryModel> table_scroll = (TableScroll<HistoryModel>) getTopComponent();			
			int[] selected_rows = table_scroll.getSelectedModelRows();
			posteriors.setMinimumObserved(2);
			setUnobservedCorrections();
//			table_model.fireTableDataChanged();
			table_scroll.setSelectedModelRows(selected_rows);
		});
        return minCopiesB;
	}
	
	private void initComponents()
	{
        getTreeControl().addControlAt(2, createMincopiesControl());
        getLineageControl().add(createMincopiesControl()); // sets min_copiesRB too
        
        JTabbedPane selected_rows_display = this.getSelectedRowsDisplay();
		String table_legend = "<p>The inferred counts are posterior expectations across selected families. "
				+ "Copy statistics (copies, births and deaths) "
				+ "count "
				+ "	<b>ancestral</b> copies at each node <var>u</var>, or "
				+ "the ancestral copy change "
				+ "in the lineage leading to <var>u</var>. "
				+ "Ancestral copies have at least one "
				+ "descendant ortholog copy at the leaves in the subtree of <var>u</var>."
				+ "Copy statistics thus track the provenance of the observed copies at "
				+ "the terminal nodes. Deaths are copies at the parent which <var>u</var> "
				+ "does not inherit; births are ancestral copies originating at <var>u</var> "
				+ "by duplication or gain."
				+ "</p>";
		if (FAMILY_PRESENCE_BY_SURVIVAL)
		{
			table_legend = table_legend
					+ "<p>Family statistics track the history of family repertoires, and  "
					+ "count the families that have at least one ancestral copy "
					+ "with descendant ortholog(s). ";
		} else
		{ 
			table_legend = table_legend
				+ "<p>Family statistics track the history of family repertoires, and  "
				+ "count the families that have at least one copy "
				+ "at the ancestor, with or without descendant orthologs. "
				+ "(Thus, more families may be inferred than ancestral copies.) ";
		}
		table_legend = table_legend
				+ "Families with multiple copies are also inferred, and changes "
				+ "to and from 0 (gain and loss), or between single- and multiple-members "
				+ "(contraction and expansion)."
				+"</p>"
				+ "<p><b>Family correction</b> adds statistics for the implied unobserved families with "
				+ "fewer copy numbers than selected by the radio buttons.</p>"
//				+ "<p>In the family table, "
//				+ "rarity is the negative log-likelihood of the family profile, "
//				+ "as -10 log<sub>10</sub> <var>L</var>(profile). Family statistics "
//				+ "(gains, losses, expansions, contractions) "
//				+ "are summed across all nodes</p>."
				;
		JEditorPane table_explain = new JEditorPane("text/html", table_legend);
        table_explain.setEditable(false);		
		selected_rows_display.addTab("Legend", table_explain);
        
	}
	
//	@Override
//	protected  void computeAll()
//	{
////		System.out.println("#**PV.cA start");
//		initDataStructures();
////		System.out.println("#**PV.cA init");
//		initComponents();
////		System.out.println("#**PV.cA comp");
//		
//		super.computeAll();		
//	}
	
	/**
	 * Call before {@link #computeAll()}. 
	 * 
	 * @param absolute
	 * @param relative
	 */
	public void setCalculcationWidthThresholds(int absolute, double relative)
	{
		posteriors.setCalculationWidthThresholds(absolute, relative);
	}

//	public void computeAll(int absolute, double relative)
//	{
//		posteriors.setCalculationWidthThresholds(absolute, relative);
//		this.computeAll();
//	}
	
	
	public int getCalculationWidthAbsolute()
	{
		return posteriors.getCalculationWidthAbsolute();
	}
	
	public double getCalculationWidthRelative()
	{
		return posteriors.getCalculationWidthRelative();
	}
	
	
	private void setUnobservedCorrections()
	{
		
		int min_copies = posteriors.getMinimumObserved();
//		System.out.println("#**PV.sUC "+min_copies);
		
		
		MixedRatePosteriors.Profile[] unobserved;
		MixedRateModel rates_model = posteriors.getRateModel();
		IndexedTree phylo = table_model.getTree();
		MixedRatePosteriors empty_posteriors=null;
		MixedRatePosteriors singleton_posteriors=null;		
		if (min_copies==0)
		{
			unobserved = new MixedRatePosteriors.Profile[0];
		}
		else
		{
			empty_posteriors = new MixedRatePosteriors(rates_model, 
					ProfileTable.emptyProfile(phylo));	
			empty_posteriors.computeClasses();
			if (min_copies == 1)
			{
				unobserved = new MixedRatePosteriors.Profile[1];
				unobserved[0] = empty_posteriors.getProfile(0);
				
			} else
			{
				singleton_posteriors = new MixedRatePosteriors(rates_model, 
						ProfileTable.singletonTable(phylo));
				singleton_posteriors.computeClasses();
				unobserved = new MixedRatePosteriors.Profile[1+phylo.getNumLeaves()];
				unobserved[0] = empty_posteriors.getProfile(0);
				for (int leaf=0; leaf<phylo.getNumLeaves(); leaf++)
				{
					unobserved[1+leaf] = singleton_posteriors.getProfile(leaf);
				}
			}
		}
		
		
		double L0 = Double.NEGATIVE_INFINITY; // log(0)
		for (int i=0; i<unobserved.length; i++)
		{
			L0 = Logarithms.add(L0, unobserved[i].getLL());
		}
		double p_unobs = Math.exp(L0);
		double p_obs = -Math.expm1(L0); // 1-exp(L0)
		double Lobs = Logarithms.logToLogComplement(L0);
		
		table_family_score.setCorrection(-EVIDENCE_SCALE*Lobs);
		table_model.fireTableChanged(new TableModelEvent(table_model, 
					0, table_model.getRowCount()-1,
					table_model.firstHistoryColumn()));

		double[] factor = new double[unobserved.length];
		double[] log_factor = new double[unobserved.length];
		for (int i=0; i<unobserved.length; i++)
		{
			double Ldiff = unobserved[i].getLL()-Lobs;
			log_factor[i] = Ldiff;
			factor[i] = Math.exp(Ldiff);
			
//			System.out.println("#**PP.sUC factor "+i+"\t"+factor[i]+"\tL0 "+L0+"\tpo "+p_obs+"\tLo "+Lobs+"\tLd "+Ldiff);
		}	
		int num_nodes = phylo.getNumNodes();
		double[] delta_member_count = new double[num_nodes];
		double[] delta_member_birth;
		double[] delta_member_death;
		double[] delta_family_present;
		double[] delta_family_multi;
		double[] delta_family_gain;
		double[] delta_family_lose;
		double[] delta_family_expand;
		double[] delta_family_contract;
		if (unobserved.length == 0)
		{
			delta_member_birth
			= delta_member_death
			= delta_family_present
			= delta_family_multi
			= delta_family_gain
			= delta_family_lose
			= delta_family_expand
			= delta_family_contract
			= delta_member_count; // ==0 everywhere
		} else
		{
			delta_member_birth = new double[num_nodes];
			delta_member_death = new double[num_nodes];
			delta_family_present = new double[num_nodes];
			delta_family_multi = new double[num_nodes];
			delta_family_gain = new double[num_nodes];
			delta_family_lose = new double[num_nodes];
			delta_family_expand = new double[num_nodes];
			delta_family_contract = new double[num_nodes];
			
			
			for (int node=0; node<num_nodes; node++)
			{
				for (int i=0; i<unobserved.length; i++)
				{
					MixedRatePosteriors.Profile P = unobserved[i];
					
					if (FAMILY_PRESENCE_BY_SURVIVAL)
					{
						double[] log_pn = P.getLogNodePosteriors(node);
						double log_present = Double.NEGATIVE_INFINITY;
						double log_mean = Double.NEGATIVE_INFINITY;
						
						int n=log_pn.length-1;
						while (1<n)
						{
							log_present = Logarithms.add(log_present, log_pn[n]);
							log_mean = Logarithms.add(log_mean, log_pn[n]+Math.log(n));
							--n;
						}
						double log_multi = log_present;
						if (0<n)
						{
							log_present = Logarithms.add(log_present, log_pn[n]);
							log_mean = Logarithms.add(log_mean, log_pn[n]);
							--n;
						}
						
						delta_member_count[node] += Math.exp(log_factor[i]+log_mean);
						double birth = P.getBirthMean(node);
						delta_member_birth[node] += factor[i]*birth;
						if (!phylo.isRoot(node))
						{
							double death = P.getDeathMean(node);
							delta_member_death[node] += factor[i]*(death);
						}
						delta_family_present[node] +=  Math.exp(log_factor[i]+log_present);
						delta_family_multi[node] += Math.exp(log_factor[i]+log_multi);
						delta_family_gain[node] += Math.exp(log_factor[i]+P.getFamilyLogSurvivalEvent(node, FamilyEvent.GAIN));
						delta_family_lose[node] += Math.exp(log_factor[i]+P.getFamilyLogSurvivalEvent(node, FamilyEvent.LOSS));
						delta_family_expand[node] += Math.exp(log_factor[i]+P.getFamilyLogSurvivalEvent(node, FamilyEvent.EXPAND));
						delta_family_contract[node] += Math.exp(log_factor[i]+P.getFamilyLogSurvivalEvent(node, FamilyEvent.CONTRACT));
					} else
					{
						double N = P.getNodeMean(node);
						double S = P.getEdgeMean(node);
						delta_member_count[node] += factor[i]*N;
						double birth = P.getBirthMean(node);
						delta_member_birth[node] += factor[i]*birth;
						if (!phylo.isRoot(node))
						{
							double death = P.getDeathMean(node);
							delta_member_death[node] += factor[i]*(death);
						}
						double[] pN = P.getNodeAncestor(node);
						delta_family_present[node] += factor[i] * (1.0-pN[0]);
						delta_family_multi[node] += factor[i]*(1.0-pN[0]-pN[1]);
						delta_family_gain[node] += factor[i]*P.getFamilyEvent(node, FamilyEvent.GAIN);
						delta_family_lose[node] += factor[i]*P.getFamilyEvent(node, FamilyEvent.LOSS);
						delta_family_expand[node] += factor[i]*P.getFamilyEvent(node, FamilyEvent.EXPAND);
						delta_family_contract[node] += factor[i]*P.getFamilyEvent(node, FamilyEvent.CONTRACT);
					}
				}
			}
		}
		lineage_member_count.setCorrection(delta_member_count);
		lineage_member_birth.setCorrection(delta_member_birth);
		lineage_member_death.setCorrection(delta_member_death);
		lineage_family_present.setCorrection(delta_family_present);
		lineage_family_multi.setCorrection(delta_family_multi);
		lineage_family_gain.setCorrection(delta_family_gain);
		lineage_family_lose.setCorrection(delta_family_lose);
		lineage_family_expand.setCorrection(delta_family_expand);
		lineage_family_contract.setCorrection(delta_family_contract);
		
		getSelectionTotals().getModel().fireTableDataChanged();
		
	}
	
	@Override
	public void computeAll()
	{
		for (MinCopies RB: min_copiesRB)
		{
			RB.setEnabled(false);
		}
		super.computeAll();
	}
	
	
	@Override
	protected void computeHistoriesPrepare()
	{
		//System.out.println("#**PV.cHP start ");
		super.computeHistoriesPrepare();
		posteriors.computeClasses();
		setUnobservedCorrections();
		//System.out.println("#**PV.cHP done ");
	}
	
	@Override
	protected void computeDone()
	{
		int data_min_copies = table_data.getContent().minCopies();
		for (int m=0; m<=2; m++)
		{
			MinCopies RB = min_copiesRB.get(m);
			RB.setEnabled(m<=data_min_copies);
		}
		
		computation_cancel.setVisible(false);
		
		if (getTableScroll().getDataTable().getSelectedRowCount()==0)
				getTableScroll().getDataTable().selectAll();
	}
	
	@Override
	protected void computeFamily(int f)
	{
		MixedRatePosteriors.Profile P = posteriors.getProfile(f);
		
//		if (f==3337) // DEBUG
//			System.out.println("#**PV.cF "+f);
		
		double LL = P.getLL();
		table_family_score.setValue(f, EVIDENCE_SCALE*LL);

		IndexedTree phylo = table_model.getTree();
		int num_nodes = phylo.getNumNodes();
		for (int node=0; node<num_nodes; node++)
		{
			double birth = P.getBirthMean(node);
			table_member_birth.get(node).setValue(f, birth);
			double death;
			
			if (phylo.isRoot(node))
			{
				death = 0.0;
			} else 
			{
				death = P.getDeathMean(node);
			}
			table_member_death.get(node).setValue(f, death);
			
			
			if (FAMILY_PRESENCE_BY_SURVIVAL)
			{
				double[] log_pn = P.getLogNodePosteriors(node);
				double log_present = Double.NEGATIVE_INFINITY;
				double log_mean = Double.NEGATIVE_INFINITY;
				
				int n=log_pn.length-1;
				while (1<n)
				{
					log_present = Logarithms.add(log_present, log_pn[n]);
					log_mean = Logarithms.add(log_mean, log_pn[n]+Math.log(n));
					--n;
				}
				double log_multi = log_present;
				if (0<n)
				{
					log_present = Logarithms.add(log_present, log_pn[n]);
					log_mean = Logarithms.add(log_mean, log_pn[n]);
					--n;
				}
				table_member_count.get(node).setValue(f, Math.exp(log_mean));
				
				table_family_present.get(node).setValue(f, Math.exp(log_present));
				table_family_multi.get(node).setValue(f, Math.exp(log_multi));
				
				table_family_gain.get(node).setValue(f, Math.exp(P.getFamilyLogSurvivalEvent(node, FamilyEvent.GAIN)));
				table_family_lose.get(node).setValue(f, Math.exp(P.getFamilyLogSurvivalEvent(node, FamilyEvent.LOSS)));;
				table_family_expand.get(node).setValue(f, Math.exp(P.getFamilyLogSurvivalEvent(node, FamilyEvent.EXPAND)));
				table_family_contract.get(node).setValue(f, Math.exp(P.getFamilyLogSurvivalEvent(node, FamilyEvent.CONTRACT)));
				
			} else
			{
				double N = P.getNodeMean(node);
				double S = P.getEdgeMean(node);
				
				table_member_count.get(node).setValue(f, N);

				double[] pn = P.getNodeAncestor(node);
				assert (pn.length>=2);
	
				double sum1;
				double sum2=0.0;
				for (int i=2; i<pn.length; i++)
					sum2+=pn[i];
				sum1 = sum2+pn[1];
				
				double present = Double.max(1.0-pn[0],sum1);
				double multi = Double.max(1.0-pn[0]-pn[1], sum2);
				table_family_present.get(node).setValue(f, present);
				table_family_multi.get(node).setValue(f, multi);
				
				table_family_gain.get(node).setValue(f, P.getFamilyEvent(node, FamilyEvent.GAIN));
				table_family_lose.get(node).setValue(f, P.getFamilyEvent(node, FamilyEvent.LOSS));;
				table_family_expand.get(node).setValue(f, P.getFamilyEvent(node, FamilyEvent.EXPAND));
				table_family_contract.get(node).setValue(f,P.getFamilyEvent(node, FamilyEvent.CONTRACT));
			}
		} // for node

		
		if (posteriors.getClassCount()>1)
		{
			assert (posteriors.getClassCount() == table_family_class.size());
			for (int c=0; c<table_family_class.size(); c++)
			{
				HistoryModel.Column<RoundedDouble> C = table_family_class.get(c);
				double pc = P.getClassProbability(c);
				C.setValue(f, pc);
			}
		}
	}
	
    @Override
    public String toString()
    {
        return "Posteriors @ "+
        			DataFile.chopFileExtension(rates_data.getFile().getName());
    }
	
}
