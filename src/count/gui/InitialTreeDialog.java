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


import java.awt.Dimension;
import java.io.File;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Random;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;

import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.ListSelectionModel;

import count.ds.AnnotatedTable;
import count.ds.Phylogeny;
import count.ds.UPGMA;
import count.io.DataFile;

import static java.awt.Component.LEFT_ALIGNMENT;

import java.awt.BorderLayout;
import java.awt.Color;

/**
 * 
 * A small JDialog to build simple trees.
 * Possibilities: star tree (all leaves are the root's children), 
 * random binary tree over the leaves, or 
 * an UPGMA clustering hierarchy constructed from dissimilarities 
 * of copy numbers between genomes.
 *
 */
public class InitialTreeDialog extends JDialog
{
	private static final int GAP = 12;
	private static final int SPACE = 2*GAP;
	
	private final AppFrame app;
	private final String[] leaf_names;
	private final DataFile<AnnotatedTable> table_data;
	
	/**
	 * Instantiation. 
	 * 
	 * @param app enclosing frame
	 * @param title dialog title 
	 * @param table_data used for getting leaf names
	 */
	public InitialTreeDialog(AppFrame app, String title, DataFile<AnnotatedTable> table_data)
	{
		super(app, title, true); // modal dialog
		this.app = app;
		this.table_data = table_data;
		this.leaf_names = table_data.getContent().getTaxonNames();
		
		initComponents();
	}
	
	private JButton okayB;
	private JButton cancelB;

	private JRadioButton starB;
	private JRadioButton randomB;
	private JRadioButton upgmaB;
	
	private JFormattedTextField random_seedT;
	
	private JComboBox<UPGMA.Dissimilarity> upgma_distanceCB;
	private JComboBox<UPGMA.ClusteringPolicy> upgma_methodCB;
	
	/**
	 * Set by okay or cancel
	 */
	private boolean is_approved = false;
	
	public boolean isApproved() { return is_approved;}
	
	public DataFile<Phylogeny> getTreeData()
	{
		DataFile<Phylogeny> showTree = null;
		if (isApproved())
		{
			Phylogeny tree = buildTree();
			File tree_file = new File((String)null, table_data.getFile().getName()+".tre");
			showTree = new DataFile<>(tree, tree_file);
		}
		return showTree;
	}
	
	
//	public boolean isRandomSelected() { return is_approved && randomB.isSelected();}
//	public boolean isStarSelected() { return is_approved && starB.isSelected();}
//	public boolean isUPGMASelected() { return is_approved && upgmaB.isSelected();}
//	
	private long getRandomSeed() 
	{
		long seed = ((Number) random_seedT.getValue()).longValue();
		return seed;
	}
	
//	public UPGMA.Dissimilarity getUPGMADissimilarity()
//	{
//		return upgma_distanceL.getSelectedValue();
//	}
//	
//	public UPGMA.ClusteringPolicy getUPGMAClusteringPolicy()
//	{
//		return upgma_methodL.getSelectedValue();
//	}
	
	private void initComponents()
	{
		Box panel = new Box(BoxLayout.PAGE_AXIS);
		panel.setBackground(Color.WHITE);
		panel.setOpaque(true);
		
		panel.setAlignmentX(LEFT_ALIGNMENT);
		JComponent star_box = createStarBox();
		JComponent rnd_box = createRandomBox();
		JComponent upgma_box = createUPGMABox();
		panel.add(star_box);
		panel.add(rnd_box);
		panel.add(upgma_box);
		
		ButtonGroup method_buttons = new ButtonGroup();
        method_buttons.add(starB);
        method_buttons.add(randomB);
        method_buttons.add(upgmaB);
        upgmaB.setSelected(true);
        
        setLayout(new BorderLayout());
        add(panel, BorderLayout.CENTER);
        add(createButtonBox(),BorderLayout.PAGE_END);
        
        Dimension frameD = app.getSize();
        pack();
        setBounds((int)(0.5*frameD.width),(int)(0.5*frameD.height),getWidth(),getHeight());

		getRootPane().setDefaultButton(okayB);
	
	}
	
	private void synchronizeFieldStates()
	{
		random_seedT.setEnabled(randomB.isSelected());
    	random_seedT.setEditable(randomB.isSelected());		        	
    	upgma_distanceCB.setEnabled(upgmaB.isSelected());
    	
    	upgma_methodCB.setEnabled(upgmaB.isSelected());
	}
	
	private JComponent createStarBox()
	{
		starB = new JRadioButton("Star tree");
		starB.addActionListener(click->synchronizeFieldStates());

		Box star_box = new Box(BoxLayout.X_AXIS);
//        star_box.setAlignmentX(LEFT_ALIGNMENT); 
        star_box.add(starB);
        star_box.add(Box.createHorizontalGlue());
        
        star_box.setBorder(BorderFactory.createEtchedBorder());
        
        return star_box;
	}
	
	private JComponent createRandomBox()
	{
        randomB = new JRadioButton("Random tree: ");
        randomB.addActionListener(click->synchronizeFieldStates());

        Box random_box = new Box(BoxLayout.X_AXIS);
//        random_box.setAlignmentX(LEFT_ALIGNMENT);
        random_box.add(randomB);
        random_box.add(Box.createHorizontalStrut(SPACE));
        
        
        NumberFormat roundF = NumberFormat.getIntegerInstance();
        random_seedT = new JFormattedTextField(roundF);
        random_seedT.setColumns(22); // large enough for long values
        JLabel random_seedL = new JLabel("seed");
        random_seedL.setLabelFor(random_seedT);
        Random RND = new Random();
        random_seedT.setValue(RND.nextLong());
        
        random_box.add(random_seedL);
        random_box.add(Box.createHorizontalStrut(GAP));
        random_box.add(random_seedT);
        
        random_seedT.setEditable(true);
        
        random_box.add(Box.createHorizontalGlue());
        random_box.setBorder(BorderFactory.createEtchedBorder());
        return random_box;
	}
	
	private JComponent createUPGMABox()
	{
        upgmaB = new JRadioButton("UPGMA or NJ tree:");
        upgmaB.addActionListener(click->synchronizeFieldStates());

        
        
        Box upgma_box = new Box(BoxLayout.X_AXIS);
//        upgma_box.setAlignmentX(LEFT_ALIGNMENT);
        upgma_box.add(upgmaB);
        upgma_box.add(Box.createHorizontalStrut(SPACE));
        
        upgma_distanceCB = new JComboBox<>(UPGMA.Dissimilarity.values());
        JLabel distL = new JLabel("Dissimilarity statistic");
        distL.setLabelFor(upgma_distanceCB);
        upgma_box.add(distL);
        upgma_box.add(Box.createHorizontalStrut(GAP));
        upgma_box.add(upgma_distanceCB);
        upgma_distanceCB.setSelectedIndex(0);
       
        upgma_box.add(Box.createHorizontalStrut(SPACE));

        upgma_methodCB = new JComboBox<>(UPGMA.ClusteringPolicy.values());
        JLabel clusteringL = new JLabel("Algorithm");
        clusteringL.setLabelFor(upgma_methodCB);

        upgma_box.add(clusteringL);
        upgma_box.add(Box.createHorizontalStrut(GAP));
        upgma_box.add(upgma_methodCB);
        upgma_methodCB.setSelectedIndex(0);
        
        upgma_box.add(Box.createHorizontalGlue());
        upgma_box.setBorder(BorderFactory.createEtchedBorder());
        return upgma_box;
	}
	
	private JComponent createButtonBox()
	{
        okayB = new JButton("OK");
        okayB.addActionListener(click->{
        	is_approved = true; dispose();
        });
        cancelB = new JButton("Cancel");
        cancelB.addActionListener(click->{
        	is_approved = false; dispose();
        });
        Box button_box = new Box(BoxLayout.LINE_AXIS);
        button_box.add(Box.createHorizontalGlue());
        button_box.add(cancelB);
        button_box.add(Box.createRigidArea(new Dimension(10,0)));
        button_box.add(okayB);
        button_box.setBorder(BorderFactory.createEmptyBorder(0, 10, 10, 10));
        return button_box;
	}
	

	private Phylogeny buildTree()
	{
		Phylogeny buildTree = null;
		if (starB.isSelected())
		{
			buildTree = Phylogeny.starTree(leaf_names);
//			System.out.println("#**ITD.bT star "+Arrays.toString(leaf_names));
		} else if (randomB.isSelected())
		{
			long seed = getRandomSeed();
			Random RND;
			if (seed==0L)
				RND = new Random();
			else
				RND = new Random(seed);
			
//			System.out.println("#**ITD.bT random "+seed);
			buildTree = Phylogeny.randomTree(leaf_names, RND, true);
		} else if (upgmaB.isSelected())
		{
			UPGMA.Dissimilarity dist = upgma_distanceCB.getItemAt(upgma_distanceCB.getSelectedIndex());
			UPGMA.ClusteringPolicy method = upgma_methodCB.getItemAt(upgma_methodCB.getSelectedIndex());
			
//			System.out.println("#**ITD.bT upgma dist "+dist+"\tmethod "+method);
			
			if (dist != null && method != null)
			{
				buildTree = UPGMA.buildTree(table_data.getContent(), dist, method);
			}
		}
		return buildTree;
	}
	
}
