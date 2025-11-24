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
import javax.swing.JCheckBox;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.io.DataFile;
import count.model.Parsimony;

/**
 * {@link HistoryView} with ancestral reconstruction by Dollo parsimony.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */
public class DolloView extends HistoryView
{
	public DolloView(DataFile<? extends IndexedTree> tree_data, 
			DataFile<AnnotatedTable> table_data)
	{
		super(tree_data.getContent(), table_data, false, true);
		this.tree_data = tree_data;
		this.factory = new Parsimony(tree_data.getContent(), table_data.getContent());
		initComponents();
	}
	
	private final DataFile<? extends IndexedTree> tree_data;
	
	private final Parsimony factory;

	private boolean root_surely_present = false;
	private JCheckBox root_cb;
	
	public boolean presentAtRoot() { return root_surely_present;}
	
	private void initComponents()
	{
		root_cb = new JCheckBox("Present at root");
		root_cb.addChangeListener(chg->
							{
								if (root_surely_present != root_cb.isSelected())
								{
									root_surely_present = root_cb.isSelected();
									computeAll();
								}
							});
        getTreeControl().addControlAt(2,root_cb);
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
		
		root_cb.setEnabled(false);
		super.computeAll();
	}
	@Override
	protected void computeDone()
	{
		root_cb.setEnabled(true);
	}
	@Override
	protected void computeHistoriesPrepare()
	{
		super.computeHistoriesPrepare();
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
//			System.out.println("#**DV.cF "+f);
//
		Parsimony.Profile P = profiles[f];
		int[] copies = P.computeDollo(root_surely_present);
		column_family_score.setValue(f, P.getDolloScore());
		
		IndexedTree phylo = tree_data.getContent();
		int num_nodes = phylo.getNumNodes();
		for (int node=0; node<num_nodes; node++)
		{
			table_family_present.get(node).setValue(f, copies[node]);
		}
	}

    @Override
    public String toString()
    {
        return "Dollo @ "+DataFile.chopCommonFileExtension(tree_data.getFile().getName());
    }
	
}
