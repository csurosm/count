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
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.table.TableRowSorter;

import count.ds.IndexedTree;
import count.ds.Phylogeny;
import count.gui.kit.DiscreteGammaPlot;
import count.io.CommandLine;
import count.io.DataFile;
import count.io.NewickParser;
import count.io.RateVariationParser;
import count.io.SavableData;
import count.matek.DiscreteGamma;
import count.model.GammaInvariant;
import count.model.TreeWithRates;

import static count.gui.AnnotatedTablePanel.TABLE_FONT_SIZE;
import static count.gui.RatesTreePanel.RATES_EDGELENGTH_COLOR;
import static count.gui.RatesTreePanel.RATES_GAIN_COLOR;
import static count.gui.RatesTreePanel.RATES_LOSS_COLOR;
import static count.gui.RatesTreePanel.RATES_DUPLICATION_COLOR;;
/**
*
* A component with split panes to show {@link count.model.GammaInvariant} instance. 
* 
* @author Mikl&oacute;s Cs&#369;r&ouml;s 
*/
public class RateVariationPanel extends JPanel // Browser.PrimaryItem
	implements SavableData<GammaInvariant>
{
	public RateVariationPanel(DataFile<GammaInvariant> model_data)
	{
		super(); // vanilla instantiation
		this.model_data = model_data;
		initComponents();
	}
	
	private final DataFile<GammaInvariant> model_data;

    public static final int RATE_VARIATION_PLOT_WIDTH = 320;
    public static final int RATE_VARIATION_PLOT_HEIGHT = 60;
    public static final int RATE_VARIATION_PADDING = 5;
	
    private Zoom<RatesTreePanel> tree_control;
    private RatesTable rates_table;
    
    @Override // Savable interface
    public DataFile<GammaInvariant> getDataFile(){ return model_data;} 
    protected Zoom<RatesTreePanel> getRatesPanel(){ return tree_control;}
    protected RatesTable getRatesTable() { return rates_table;}
    
	private void initComponents()
	{
        JSplitPane main_pane = new JSplitPane();
        // horizontal split: tree at bottom
        this.setBorder(null);
        setBackground(Color.WHITE);
        
        main_pane.setBackground(getBackground());
        main_pane.setBorder(null);
        main_pane.setDividerLocation(300+getInsets().top);
        main_pane.setResizeWeight(0.5);
        main_pane.setOrientation(JSplitPane.VERTICAL_SPLIT);
        main_pane.setOneTouchExpandable(true);
        
        TreeWithRates rates = model_data.getContent().getBaseModel();
  
        
        RatesTreePanel tree_panel = new RatesTreePanel(new DataFile<>(rates));
        
//        selection_model.addListSelectionListener(e->
//        	{
//                if (!e.getValueIsAdjusting())
//                {
//                    int row_idx = selection_model.getLeadSelectionIndex();
//                    if (row_idx != -1)
//                        rates_table.scrollRectToVisible(rates_table.getCellRect(row_idx, 0, true));
//                }
//        	});  
        
        // add control for zooming and parameter selection
        tree_control = new Zoom<>(tree_panel);
        Box control_bar = tree_control.getControlBar();
        control_bar.removeAll();
        control_bar.add(tree_panel.createParameterChooser());
        control_bar.add(Box.createHorizontalGlue());
        control_bar.add(tree_control.getSpinner());
        
        rates_table = new RatesTable(rates);
        rates_table.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        TableRowSorter<RatesTable.Model> sorter = new TableRowSorter<>(rates_table.getModel());
        rates_table.setRowSorter(sorter);
        rates_table.synchronizeModelSelection(tree_panel.getSelectionModel());
//        JScrollPane table_scroll = new JScrollPane(rates_table);
//        table_scroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
//        table_scroll.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
//        table_scroll.getViewport().setBackground(getBackground());
        
        JPanel rate_multiplier_panel = new RateMultiplierPanel();
        rate_multiplier_panel.setBackground(getBackground());
        JScrollPane rate_variation_scroll = new JScrollPane(rate_multiplier_panel);
        rate_variation_scroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        rate_variation_scroll.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        
        JSplitPane top_pane = new JSplitPane();
        top_pane.setBorder(null);
        top_pane.setOrientation(JSplitPane.HORIZONTAL_SPLIT);
        top_pane.setLeftComponent(rates_table);
        top_pane.setRightComponent(rate_variation_scroll);
        top_pane.setResizeWeight(0.1);
        top_pane.setOneTouchExpandable(true);
        
        main_pane.setTopComponent(top_pane);

        
        main_pane.setBottomComponent(tree_control);
        
        BorderLayout layout = new BorderLayout();
        this.setLayout(layout);
        
        this.add(main_pane, BorderLayout.CENTER);
    }	

	protected JComponent createScoringBox()
	{
		JButton score_button = new JButton("Score");
		return score_button;
	}
	
	
	
    private class RateMultiplierPanel extends JPanel
    {
        private final Font normal_font = new Font("Serif",Font.PLAIN,TABLE_FONT_SIZE);
        private final Font title_font  = new Font("Serif", Font.BOLD, TABLE_FONT_SIZE*6/5);
        private final Font subtitle_font = normal_font.deriveFont(Font.BOLD);
        private final Font it_font = normal_font.deriveFont(Font.ITALIC);
        private final Font plot_font = new Font("Serif", Font.PLAIN, TABLE_FONT_SIZE*4/5);
        
        RateMultiplierPanel()
        {
        	super();
        }
        
        private void initComponents()
        {
            setBackground(Color.WHITE);
        }
        
        @Override 
        public Dimension getPreferredSize()
        {
            int w = RATE_VARIATION_PLOT_WIDTH*3/2 + 2*RATE_VARIATION_PADDING;
            int h = TABLE_FONT_SIZE*6/5+2*RATE_VARIATION_PADDING+4*(RATE_VARIATION_PADDING+RATE_VARIATION_PLOT_HEIGHT);
            return new Dimension(w,h);
        }
        
        @Override 
        public Dimension getMinimumSize()
        {
            return getPreferredSize();
        }
        
        @Override
        protected void paintComponent(Graphics g)
        {
            GammaInvariant rates = model_data.getContent();

            Graphics2D g2 = (Graphics2D) g.create();
            
            int current_y = title_font.getSize()+RATE_VARIATION_PADDING;
            g2.setFont(title_font);
            g2.drawString("Rate variation across families", RATE_VARIATION_PADDING, current_y);
            
            int subtitle_length = g2.getFontMetrics(subtitle_font).stringWidth("Duplication rate:"); // longest subtitle 
            current_y+=RATE_VARIATION_PADDING;
            this.paintRateMultipliers(g2, rates.getLengthAlpha(), rates.getNumLengthGammaCategories(), 0.0, current_y, RATES_EDGELENGTH_COLOR, "Edge length:", subtitle_length, "");
            current_y += RATE_VARIATION_PLOT_HEIGHT;
            current_y += RATE_VARIATION_PADDING;
            this.paintRateMultipliers(g2, rates.getLossAlpha(), rates.getNumLossGammaCategories(), rates.getLossForbidden(), current_y, RATES_LOSS_COLOR, "Loss rates:", subtitle_length, "P(no-loss)=");
            current_y += RATE_VARIATION_PLOT_HEIGHT;
            current_y += RATE_VARIATION_PADDING;
            this.paintRateMultipliers(g2, rates.getDuplicationAlpha(), rates.getNumDuplicationGammaCategories(), rates.getDuplicationForbidden(), current_y, RATES_DUPLICATION_COLOR, "Duplication rates:", subtitle_length, "P(no-duplication)=");
            current_y += RATE_VARIATION_PLOT_HEIGHT;
            current_y += RATE_VARIATION_PADDING;
            this.paintRateMultipliers(g2, rates.getGainAlpha(), rates.getNumGainGammaCategories(), rates.getGainForbidden(), current_y, RATES_GAIN_COLOR, "Gain rates:", subtitle_length, "P(no-gain)=");
        } // paintComponent
        
        private void paintRateMultipliers(Graphics2D g2, double alpha, int n, double p0, int current_y, Color color, String subtitle, int subtitle_length, String forbidden_label)
        {
            g2.setFont(subtitle_font);
            g2.setColor(color);
            g2.drawString(subtitle, RATE_VARIATION_PADDING, current_y+subtitle_font.getSize());
            if (n==1)
            {
                int txt_x = RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING;
                int txt_y = current_y+subtitle_font.getSize();
                if (p0==0.0)
                {
                    g2.setFont(it_font);
                    g2.drawString("no variation", txt_x, txt_y);
                } else
                {
                    g2.setFont(plot_font);
                    g2.drawString(forbidden_label+Double.toString(p0), txt_x, txt_y);
                }
            } else
            {
                DiscreteGammaPlot P_duplication = new DiscreteGammaPlot(new DiscreteGamma(alpha), n, color);
                g2.setFont(plot_font);
                P_duplication.paintIcon(this, g2, RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING, current_y);
                if (p0!=0.0)
                {
                    int txt_x = RATE_VARIATION_PADDING+subtitle_length+RATE_VARIATION_PADDING+RATE_VARIATION_PLOT_WIDTH+RATE_VARIATION_PADDING;
                    int txt_y = current_y+subtitle_font.getSize();
                    g2.setFont(plot_font);
                    g2.drawString(forbidden_label+Double.toString(p0), txt_x, txt_y);
                }
            }            
        }
        
        
    }	

    /**
     * Returns the associated file name (for displaying in the Browser).
     */
    @Override
    public String toString()
    {
        return DataFile.chopFileExtension(model_data.getFile().getName());
    }   
    
    
    
    @Override
    public void saveData(File f) throws IOException
    {
        PrintStream PS = new PrintStream(f);
        PS.println(CommandLine.getStandardHeader(getClass()));
        PS.println(CommandLine.getStandardRuntimeInfo());
        GammaInvariant rate_model = model_data.getContent();
        String model_description = RateVariationParser.printRates(rate_model);
        IndexedTree tree = rate_model.getBaseModel().getTree();
        if (tree instanceof Phylogeny)
        {
        	String tree_description = NewickParser.printTree((Phylogeny) tree);
        	PS.println("# tree "+tree_description);
        }
        PS.println(model_description);

        if (PS.checkError()) // also flushes
        {
            PS.close();
            throw new IOException("Cannot write the table.");
        }
        PS.close();
        model_data.setFile(f);
        model_data.setDirty(false);
    }

    /**
     * Number of free parameters in the underlying 
     * GammaInvariant model. 
     * 
     * @return
     */
    public int getModelParameterCount()
    {
    	int parameter_count = 0;
    	GammaInvariant rates_model = getDataFile().getContent();
    	TreeWithRates base_model = rates_model.getBaseModel();
    	IndexedTree phylo = base_model.getTree() ;
    	boolean has_gain = false;
    	boolean has_duplication = false;
    	for (int node=0; node<phylo.getNumNodes(); node++)
    	{
    		double len = base_model.getEdgeLength(node);
    		assert (len != 0.0); // children with 0-length edges are fused into the parent already within the phylogeny
    		if (!Double.isInfinite(len))
    		{
    			++parameter_count; // edge length
    		}
    		if (base_model.getDuplicationRate(node)!=0.0)
    		{
    			++parameter_count;
    			has_duplication = has_duplication || !phylo.isRoot(node);
    		}
    		if (base_model.getGainRate(node)!=0.0)
    		{
    			++parameter_count;
    			has_gain = has_gain || !phylo.isRoot(node);
    		}
    	}
    	int cat_length = rates_model.getNumLengthGammaCategories();
    	int cat_duplication = has_duplication
			    			?rates_model.getNumDuplicationGammaCategories()
							:1;
    	int cat_gain = has_gain
					?rates_model.getNumDuplicationGammaCategories()
					:1;
    	
    	if (cat_duplication>1) 
    		++parameter_count; // alpha
    	if (cat_gain>1)
    		++parameter_count; // alpha
    	if (cat_length>1)
    		++parameter_count; // alpha
    	return parameter_count;
    }
    
    
    
}
