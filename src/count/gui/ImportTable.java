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

import javax.swing.border.Border;

import static count.io.CommandLine.OPT_OUTPUT;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.MouseListener;
import java.awt.event.MouseEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.IOException;

import java.text.NumberFormat;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JEditorPane;
import javax.swing.SwingWorker;

import count.ds.AnnotatedTable;
import count.ds.Phylogeny;
//import count.gui.kit.FancyFileDialog;
import count.gui.kit.ShmancyFileDialog;
import count.io.CommandLine;
import count.io.DataFile;
import count.io.TableParser;

/**
 * 
 * Elements:
 * genome.tab chooser
 * COG.csv chooser
 * COG columns 1 2 3 
 * ---
 * p2o chooser 
 * mcloutput chooser
 * 
 * @author csuros
 *
 */
public class ImportTable extends JPanel
{
	public ImportTable(AppFrame app)
	{
		super(new BorderLayout());
		this.app = app;
		initComponents();
	}
	
	private final AppFrame app;
	
	private JRadioButton cogRB;
	private JRadioButton mclRB;
	
	/*
	 * COG interface
	 */
	private ShmancyFileDialog genomeD;
	private ShmancyFileDialog cogD;
	private ShmancyFileDialog p2oD;
	private ShmancyFileDialog mclD;
	
	private JFormattedTextField col_genomeT;
	private JFormattedTextField col_geneT;
	private JFormattedTextField col_familyT;
	private JFormattedTextField col_memberT;
	
	private JCheckBox uniqueCB;
	
	
	private JRadioButton member0RB;
	private JRadioButton member1RB;
	private JRadioButton member2RB;
	private JRadioButton member3RB;
	
	//
	//
	
	Box cogB;
	Box mclB;
	
	Border mclBorder=null;
	Border cogBorder=null;
	
	/**
	 * Formatted text box size for long explanations.
	 * 
	 */
	private static final Dimension EXPLAIN_SIZE = new Dimension(480,180);
	
	private void initComponents()
	{
		this.setAlignmentX(LEFT_ALIGNMENT);
		Box formatsB = new Box(BoxLayout.PAGE_AXIS);
		formatsB.setAlignmentX(LEFT_ALIGNMENT);
		ButtonGroup formats = new ButtonGroup();
		cogB = createCOGBox();
		formatsB.add(cogB);
		formats.add(cogRB);
		
		mclB = createMCLBox();
		formatsB.add(mclB);
		formats.add(mclRB);
		this.add(formatsB, BorderLayout.CENTER);
		this.add(cogRB, BorderLayout.NORTH);
		this.add(mclRB, BorderLayout.SOUTH);
		cogRB.addItemListener(e->colorBorders());
		mclRB.addItemListener(e->colorBorders());
		
		mclB.addMouseListener(new MouseListener() {
			@Override
			public void mouseClicked(MouseEvent e) {
				mclRB.setSelected(true);
			}
			@Override
			public void mousePressed(MouseEvent e) {}
			@Override
			public void mouseReleased(MouseEvent e) {}
			@Override
			public void mouseEntered(MouseEvent e) {}
			@Override
			public void mouseExited(MouseEvent e) {}
		});
		cogB.addMouseListener(new MouseListener() {
			@Override
			public void mouseClicked(MouseEvent e) {
				cogRB.setSelected(true);
			}
			@Override
			public void mousePressed(MouseEvent e) {}
			@Override
			public void mouseReleased(MouseEvent e) {}
			@Override
			public void mouseEntered(MouseEvent e) {}
			@Override
			public void mouseExited(MouseEvent e) {}
		});
		
		mclRB.setSelected(true);
		cogRB.setSelected(true);
	}
	
	private Box createCOGBox()
	{
		Box createCOGBox = new Box(BoxLayout.PAGE_AXIS)
		{
			@Override
			public void setEnabled(boolean e)
			{
				genomeD.setEnabled(e);
				cogD.setEnabled(e);
				uniqueCB.setEnabled(e);
				member0RB.setEnabled(e);
				member1RB.setEnabled(e);
				member2RB.setEnabled(e);
				member3RB.setEnabled(e);
				col_genomeT.setEnabled(e);
				col_geneT.setEnabled(e);
				col_familyT.setEnabled(e);
				col_memberT.setEnabled(e);
			}
		};
		createCOGBox.setAlignmentX(LEFT_ALIGNMENT);
		cogRB =  new JRadioButton("Membership data");
//		createCOGBox.add(cogRB);
		
		String cog_format_description = "<em>The input file is a comma- or tab-separated text file, "
				+ "where each line describes one membership relation: genome, gene, family/domain, membership "
				+ "(e.g., <code>COG.csv</code> files). The membership class is determined by the "
				+ "alignment footprint along the protein and the family profile: "
				+ "0=most (≥75%) of profile and most ((≥75%) of protein, "
				+ "1=most of profile and part of protein, "
				+ "2=most of protein and part of profile, "
				+ "3=part of protein and part of profile. "
				+ "An optional genome id-mapping file"
				+ "can be specified, which lists in its first two columns the genome "
				+ "name used in the input file, and the corresponding taxon identifier "
				+ "to be used in the phylogeny. Leave it empty if "
				+ "the input file's genome identifiers can be used directly."
				+ "For a very generic usage, you can set the columns for the gene, genome, "
				+ "family, and membership-class values (0, 1, 6, and 7 [or 8] in COG datasets). "
				+ "If the family column is missing, or family name is empty, then "
				+ "the gene will be assumed to belong to a new orphan family. "
				+ "If the membership column is missing, then all memberships are included."
				+ "</em>";
		JEditorPane cog_format_explain = new JEditorPane("text/html", cog_format_description);
		cog_format_explain.setEditable(false);
		cog_format_explain.setPreferredSize(EXPLAIN_SIZE);
		createCOGBox.add(new JScrollPane(cog_format_explain));
		
		genomeD = new ShmancyFileDialog(app, "Open genome name-to-id map", "Genome name-to-id file");
		cogD = new ShmancyFileDialog(app, "Open membership table", "(COG) membership file");
		createCOGBox.add(genomeD.getEmbeddableChooser());
		createCOGBox.add(cogD.getEmbeddableChooser());
		
		Box uniqueB = new Box(BoxLayout.LINE_AXIS);
		String unique_description = "\u261c <em>If unchecked, all annotated domains are counted within a gene; "
				+ "otherwise only the first one.</em>";
		JEditorPane unique_explain = new JEditorPane("text/html", unique_description);
		unique_explain.setOpaque(false);
		unique_explain.setPreferredSize(new Dimension(30,20));
		unique_explain.setEditable(false);
				
		uniqueCB = new JCheckBox("Keep only the first annotation per gene");
		uniqueCB.setSelected(false);
		uniqueB.add(uniqueCB);
		uniqueB.add(unique_explain);
		createCOGBox.add(uniqueB);
		
		Box memberB = new Box(BoxLayout.LINE_AXIS);
		String member_description = "\u261c <em>Maximum membership class by alignment footprint along protein and family profile</em>";
		JEditorPane member_explain = new JEditorPane("text/html", member_description);
		member_explain.setOpaque(false);
		member_explain.setPreferredSize(new Dimension(30,20));
		member_explain.setEditable(false);
		
		member1RB = new JRadioButton("1");
		member0RB = new JRadioButton("0");
		member2RB = new JRadioButton("2");
		member3RB = new JRadioButton("3");

		ButtonGroup memberBG = new ButtonGroup();
		memberBG.add(member1RB);
		memberBG.add(member0RB);
		memberBG.add(member2RB);
		memberBG.add(member3RB);
		
		member1RB.setSelected(true);
		
		JLabel memberL = new JLabel("Membership");
		memberB.add(memberL);
		memberB.add(member0RB);
		memberB.add(member1RB);
		memberB.add(member2RB);
		memberB.add(member3RB);
		memberB.add(member_explain);
		createCOGBox.add(memberB);
		
		Box colB = new Box(BoxLayout.LINE_AXIS);
		
//		Box col_geneB = new Box(BoxLayout.LINE_AXIS);
		col_geneT = new JFormattedTextField(NumberFormat.getIntegerInstance());
		col_geneT.setColumns(8);
		col_geneT.setValue(0);
		JLabel col_geneL = new JLabel("Gene column");
		col_geneL.setLabelFor(col_geneT);
		colB.add(col_geneL);
		colB.add(col_geneT);
		colB.add(Box.createHorizontalGlue());

//		Box col_genomeB = new Box(BoxLayout.LINE_AXIS);
		col_genomeT = new JFormattedTextField(NumberFormat.getIntegerInstance());
		col_genomeT.setColumns(8);
		col_genomeT.setValue(1); // default COG.csv
		JLabel col_genomeL = new JLabel("Genome column");
		col_genomeL.setLabelFor(col_genomeT);
		colB.add(col_genomeL);
		colB.add(col_genomeT);
		colB.add(Box.createHorizontalGlue());
		
		
		
//		Box col_familyB = new Box(BoxLayout.LINE_AXIS);
		col_familyT = new JFormattedTextField(NumberFormat.getIntegerInstance());
		col_familyT.setColumns(8);
		col_familyT.setValue(6); // COG.csv default
		JLabel col_familyL = new JLabel("Family column");
		col_familyL.setLabelFor(col_familyT);
		colB.add(col_familyL);
		colB.add(col_familyT);
		
		colB.add(Box.createHorizontalGlue());
		col_memberT = new JFormattedTextField(NumberFormat.getIntegerInstance());
		col_memberT.setColumns(8);
		col_memberT.setValue(7); // COG.csv default
		JLabel col_memberL = new JLabel("Membership column");
		col_memberL.setLabelFor(col_memberT);
		colB.add(col_memberL);
		colB.add(col_memberT);
		
		createCOGBox.add(colB);
		
		createCOGBox.setBorder(BorderFactory.createTitledBorder("Import membership table (COG-style data)"));
		
		return createCOGBox;
	}
	
	private void colorBorders()
	{
		if (cogBorder==null)
		{
			cogBorder = cogB.getBorder();
			mclBorder = mclB.getBorder();
		}
		Box selected;
		Box unselected;
		Border selected_border;
		Border unselected_border;
		if (cogRB.isSelected())
		{
			selected = cogB;
			selected_border = cogBorder;
			unselected = mclB;
			unselected_border = mclBorder;
		} else
		{
			selected = mclB;
			selected_border = mclBorder;
			unselected = cogB;
			unselected_border = cogBorder;
		}
		selected.setBorder(BorderFactory.createCompoundBorder(selected_border, BorderFactory.createLineBorder(Color.black, 3)));
		unselected.setBorder(unselected_border);
		selected.setEnabled(true); // on a Box?
		unselected.setEnabled(false);
	}
	
	
	private Box createMCLBox()
	{
		Box createMCLBox = new Box(BoxLayout.PAGE_AXIS)
		{
				@Override
				public void setEnabled(boolean e)
				{
					p2oD.setEnabled(e);
					mclD.setEnabled(e);
				}
		};		
				
		createMCLBox.setAlignmentX(LEFT_ALIGNMENT);
		createMCLBox.setAlignmentY(TOP_ALIGNMENT);
		mclRB = new JRadioButton("Clustering data");
//		createMCLBox.add(mclRB);
		
		String mcl_format_description = "<em>The input file is a list of clusters (-families), "
				+ "where each line is a tab- (or comma-) separated list of cluster members (gene identifiers). "
				+ "For example, "
				+ "such is the MCL output file generated by Stijn van Dongen's mcl software."
				+ "The auxiliary "
				+ "gene-to-genome mapping file (e.g.<code>p2o.csv</code>) is mandatory: "
				+ "each line gives a gene id and taxon id in the first two columns."
				+ "Families will be named by consecutive indexes (<code>ogNNNN</code> for multi-genome "
				+ "families, <code>lseNNNN</code> for lineage-specific expansions)."
				+ "</em>";
		
		JEditorPane mcl_format_explain = new JEditorPane("text/html", mcl_format_description);
		mcl_format_explain.setEditable(false);
		mcl_format_explain.setPreferredSize(EXPLAIN_SIZE);
		createMCLBox.add(new JScrollPane(mcl_format_explain));
		p2oD = new ShmancyFileDialog(app, "Open gene-to-genome map", "Gene-to-genome file");
		mclD = new ShmancyFileDialog(app, "Open clustering output", "Clusters file");
		createMCLBox.add(p2oD.getEmbeddableChooser());
		createMCLBox.add(mclD.getEmbeddableChooser());
		
		createMCLBox.setBorder(BorderFactory.createTitledBorder("Import clustering table (MCL-style data)"));
//		createMCLBox.addMouseListener(e->mclRB.setSelected(true));
		
		return createMCLBox;
	}
	
	public SwingWorker<DataFile<AnnotatedTable>, Void> readTableTask(String[] taxon_names)
	{
		SwingWorker<DataFile<AnnotatedTable>, Void> read_task; 
		if (cogRB.isSelected())
		{
			int col_genome = ((Number)col_genomeT.getValue()).intValue();
			int col_gene = ((Number)col_geneT.getValue()).intValue();
			int col_family = ((Number)col_familyT.getValue()).intValue();
			int col_membership = ((Number)col_memberT.getValue()).intValue();
			boolean unique = uniqueCB.isSelected();
			int max_member;
			if (member1RB.isSelected())
				max_member=1;
			else if (member2RB.isSelected())
				max_member=2;
			else if (member3RB.isSelected())
				max_member = 3;
			else
			{
				assert (member0RB.isSelected());
				max_member = 0;
			}
			File table_file = new File(cogD.getDirectory(), cogD.getFile()+".txt");
			
			read_task = cogD.createTask(cogs->{ // only called with non-null argument
				BufferedReader genomeR = genomeD.getBufferedReader();
				Map<String,String> genome_map = null;
				if (genomeR != null)
				{
					genome_map = TableParser.readTwoColumnMap(genomeR);
					genomeR.close();
				}
				BufferedReader tableR  = new BufferedReader(new InputStreamReader(cogs));
//				System.out.println("#**IT.rTT start reading "+Thread.currentThread());
				AnnotatedTable taboola = TableParser.readMembershipData(tableR, col_gene, col_genome, col_family, col_membership, max_member, genome_map, taxon_names, unique);
//				System.out.println("#**IT.rTT done reading nfam "+taboola.getFamilyCount()+"\t"+ Thread.currentThread());
				return new DataFile<>(taboola, table_file);
			});
		} else
		{
			assert mclRB.isSelected();
			File table_file = new File(mclD.getDirectory(), mclD.getFile()+".txt");
			read_task = mclD.createTask(mcl->{
				BufferedReader geneR = p2oD.getBufferedReader();
				if (geneR==null) return null;
				Map<String,String> gene_map = TableParser.readTwoColumnMap(geneR);
				geneR.close();
				BufferedReader tableR = new BufferedReader(new InputStreamReader(mcl));
				AnnotatedTable taboola = TableParser.readClusteringData(tableR, gene_map, taxon_names);
				return new DataFile<>(taboola, table_file);
			});
		}
		return read_task;
	}
	
//	public DataFile<AnnotatedTable> readTable(String[] taxon_names) throws IOException
//	{
//		DataFile<AnnotatedTable> readTable;
//		
//		
//		SwingWorker<DataFile<AnnotatedTable>, Void> read_task; 
//		
//		if (cogRB.isSelected())
//		{
//			// TODO: add membership class parsing
//			System.out.println("#**IT.rT COGs "+cogD.getFile()+"\tgenomes "+genomeD.getFile());
//			int col_genome = ((Number)col_genomeT.getValue()).intValue();
//			int col_gene = ((Number)col_geneT.getValue()).intValue();
//			int col_family = ((Number)col_familyT.getValue()).intValue();
//			boolean unique = uniqueCB.isSelected();
//			File table_file = new File(cogD.getDirectory(), cogD.getFile()+".txt");
//
//			BufferedReader genomeR = genomeD.getBufferedReader();
//			BufferedReader tableR  = cogD.getBufferedReader();
//			if (tableR == null)
//				readTable = null;
//			else
//			{
//				Map<String,String> genome_map = null;
//				if (genomeR != null)
//				{
//					genome_map = TableParser.readTwoColumnMap(genomeR);
//					genomeR.close();
//				}
//				AnnotatedTable taboola = TableParser.readMembershipData(tableR, col_gene, col_genome, col_family, genome_map, taxon_names, unique);
//				tableR.close();
//				
//				readTable = new DataFile<>(taboola, table_file);
//			}
//		} else // mcl
//		{
//			assert mclRB.isSelected();
//			BufferedReader geneR = p2oD.getBufferedReader();
//			BufferedReader tableR = mclD.getBufferedReader();
//			Map<String,String> gene_map = null;
//			if (tableR == null || geneR == null)
//			{
//				readTable = null;
//			} else
//			{
//				gene_map = TableParser.readTwoColumnMap(geneR);
//				geneR.close();
//				AnnotatedTable taboola = TableParser.readClusteringData(tableR, gene_map, taxon_names);
//				tableR.close();
//				
//				File table_file = new File(mclD.getDirectory(), mclD.getFile()+".txt");
//				readTable = new DataFile<>(taboola, table_file);
//			}
//		}
//		return readTable;
//	}
	
}
