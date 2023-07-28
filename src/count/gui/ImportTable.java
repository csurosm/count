package count.gui;

import javax.swing.border.Border;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ItemEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JEditorPane;
import javax.swing.JTextField;

import count.ds.AnnotatedTable;
import count.gui.kit.FancyFileDialog;
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
	
	
	private FancyFileDialog genomeD;
	private FancyFileDialog cogD;
	private FancyFileDialog p2oD;
	private FancyFileDialog mclD;
	
	private JFormattedTextField col_genomeT;
	private JFormattedTextField col_geneT;
	private JFormattedTextField col_familyT;
	
	private JCheckBox uniqueCB;
	
	//
	//
	
	Box cogB;
	Box mclB;
	
	Border mclBorder=null;
	Border cogBorder=null;
	
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
				col_genomeT.setEnabled(e);
				col_geneT.setEnabled(e);
				col_familyT.setEnabled(e);
			}
		};
		createCOGBox.setAlignmentX(LEFT_ALIGNMENT);
		cogRB =  new JRadioButton("Membership data");
//		createCOGBox.add(cogRB);
		
		String cog_format_description = "<em>The input file is a comma- or tab-separated text file, "
				+ "where each line describes one membership relation: genome, gene, family/domain "
				+ "(e.g., <code>COG.csv</code> files). An optional genome id-mapping file"
				+ "can be specified, which lists in its first two columns the genome "
				+ "name used in the input file, and the corresponding taxon identifier "
				+ "to be used in the phylogeny. Leave it empty if "
				+ "the input file's genome identifiers can be used directly."
				+ "For a very generic usage, you can set the columns for the gene, genome and "
				+ "family values (0, 1, and 6 in COG datasets)."
				+ "</em>";
		JEditorPane cog_format_explain = new JEditorPane("text/html", cog_format_description);
		cog_format_explain.setEditable(false);
		cog_format_explain.setPreferredSize(new Dimension(40,100));
		createCOGBox.add(cog_format_explain);
		
		genomeD = new FancyFileDialog(app, "Open genome name-to-id map", "Genome name-to-id file");
		cogD = new FancyFileDialog(app, "Open membership table", "(COG) membership file");
		createCOGBox.add(genomeD.getEmbeddableChooser());
		createCOGBox.add(cogD.getEmbeddableChooser());
		
		uniqueCB = new JCheckBox("Keep only the first annotation per gene");
		uniqueCB.setSelected(false);
		createCOGBox.add(uniqueCB);
		
		Box colB = new Box(BoxLayout.LINE_AXIS);
		
		Box col_geneB = new Box(BoxLayout.LINE_AXIS);
		col_geneT = new JFormattedTextField(NumberFormat.getIntegerInstance());
		col_geneT.setColumns(8);
		col_geneT.setValue(0);
		JLabel col_geneL = new JLabel("Gene column");
		col_geneL.setLabelFor(col_geneT);
		col_geneB.add(col_geneL);
		col_geneB.add(col_geneT);
		colB.add(col_geneB);

		Box col_genomeB = new Box(BoxLayout.LINE_AXIS);
		col_genomeT = new JFormattedTextField(NumberFormat.getIntegerInstance());
		col_genomeT.setColumns(8);
		col_genomeT.setValue(1); // default COG.csv
		JLabel col_genomeL = new JLabel("Genome column");
		col_genomeL.setLabelFor(col_genomeT);
		col_genomeB.add(col_genomeL);
		col_genomeB.add(col_genomeT);
		colB.add(col_genomeB);
		
		
		
		Box col_familyB = new Box(BoxLayout.LINE_AXIS);
		col_familyT = new JFormattedTextField(NumberFormat.getIntegerInstance());
		col_familyT.setColumns(8);
		col_familyT.setValue(6); // COG.csv default
		JLabel col_familyL = new JLabel("Family column");
		col_familyL.setLabelFor(col_familyT);
		col_familyB.add(col_familyL);
		col_familyB.add(col_familyT);
		colB.add(col_familyB);
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
				+ "Families will be named by consecutive indexes (<code>ogXXXX</code> for multi-genome "
				+ "families, <code>lseXXX</code> for lineage-specific expansions)."
				+ "</em>";
		
		JEditorPane mcl_format_explain = new JEditorPane("text/html", mcl_format_description);
		mcl_format_explain.setEditable(false);
		mcl_format_explain.setPreferredSize(new Dimension(40,100));
		createMCLBox.add(mcl_format_explain);
		p2oD = new FancyFileDialog(app, "Open gene-to-genome map", "Gene-to-genome file");
		mclD = new FancyFileDialog(app, "Open clustering output", "Clusters file");
		createMCLBox.add(p2oD.getEmbeddableChooser());
		createMCLBox.add(mclD.getEmbeddableChooser());
		
		createMCLBox.setBorder(BorderFactory.createTitledBorder("Import clustering table (MCL-style data)"));
//		createMCLBox.addMouseListener(e->mclRB.setSelected(true));
		
		return createMCLBox;
	}
	
	public DataFile<AnnotatedTable> readTable(String[] taxon_names) throws IOException
	{
		DataFile<AnnotatedTable> readTable;
		if (cogRB.isSelected())
		{
			BufferedReader genomeR = genomeD.getBufferedReader();
			BufferedReader tableR  = cogD.getBufferedReader();
			if (tableR == null)
				readTable = null;
			else
			{
				Map<String,String> genome_map = null;
				if (genomeR != null)
				{
					genome_map = TableParser.readTwoColumnMap(genomeR);
					genomeR.close();
				}
				int col_genome = ((Number)col_genomeT.getValue()).intValue();
				int col_gene = ((Number)col_geneT.getValue()).intValue();
				int col_family = ((Number)col_familyT.getValue()).intValue();
				boolean unique = uniqueCB.isSelected();
				
				AnnotatedTable taboola = TableParser.readMembershipData(tableR, col_gene, col_genome, col_family, genome_map, taxon_names, unique);
				tableR.close();
				
				File table_file = new File(cogD.getDirectory(), cogD.getFile()+".txt");
				readTable = new DataFile<>(taboola, table_file);
			}
		} else // mcl
		{
			assert mclRB.isSelected();
			BufferedReader geneR = p2oD.getBufferedReader();
			BufferedReader tableR = mclD.getBufferedReader();
			Map<String,String> gene_map = null;
			if (tableR == null || geneR == null)
			{
				readTable = null;
			} else
			{
				gene_map = TableParser.readTwoColumnMap(geneR);
				AnnotatedTable taboola = TableParser.readClusteringData(tableR, gene_map, taxon_names);
				tableR.close();
				
				File table_file = new File(mclD.getDirectory(), mclD.getFile()+".txt");
				readTable = new DataFile<>(taboola, table_file);
			}
		}
		return readTable;
	}
}
