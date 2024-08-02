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
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.text.DecimalFormat;
import java.text.Format;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import javax.swing.AbstractButton;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.InputVerifier;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;

import count.Count;
import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.gui.kit.CheckSelectAll;
import count.gui.kit.InputVerifiers;
import count.io.CommandLine;
import count.io.DataFile;
import count.io.ModelBundle;
import count.matek.DiscreteDistribution;
import count.matek.NegativeBinomial;
import count.matek.PointDistribution;
import count.matek.Poisson;
import count.matek.ShiftedGeometric;
import count.model.GammaInvariant;
import count.model.ML;
import count.model.MLDistribution;
import count.model.MLRateVariation;
import count.model.MixedRateModel;
import count.model.RateVariationModel;
import count.model.StraightEM;
import count.model.TreeWithRates;
import count.model.old.DirectEM;
import count.model.old.MLGamma;

/**
 * A modal dialog to select starting model and optimization parameters. 
 * 
 */
public class ModelSelectionDialog extends JDialog 
{
	
	private static final boolean DEVELOPMENT_CODE = false;
	
	private static final int OPT_ROUNDS = 256;
	private static final double OPT_EPS = 1.0/(1<<20);
	private static final int OPT_NUM_FAMILIES = 100;
	
	private static final boolean WANT_GENERAL_VARIATION = false; // if true, dup- and gain-rate discrete-Gamma variation is enabled
	/**
	 * Initializes the starting models.  
	 * 
	 * @param app used to figure out the selected rate model
	 * @param title
	 */
	public ModelSelectionDialog(AppFrame app, String title)
	{
		super(app, title, true); // modal dialog
		this.app = app;
		initDataStructures();
	}
	
	private final AppFrame app;
	
    private List<MixedRateModel> models;
    private List<JRadioButton> modelRB;
    private int selected_model = -1;
    private boolean want_command_line = false;
    
    /**
     * Button for starting the computations.
     */
    private JButton startB;
    /**
     * Button for canceling
     */
    private JButton cancelB;
    /**
     * Button for command-builder interface
     */
    private JButton cliB;
    
    /**
     * Radio buttons for model type
     */
    private JRadioButton model_gldRB;
    private JRadioButton model_dlRB;
    private JRadioButton model_glRB;
    private JRadioButton model_plRB;
    
    /**
     * Checkboxes for lineage-specific variation
     */
    private JCheckBox uniform_duplicationCB;
    private JCheckBox uniform_gainCB;
    private JCheckBox uniform_lossCB;
    
    /**
     * Radio buttons for distribution at root
    */
    private Map<Class<? extends DiscreteDistribution>, JRadioButton> root_distributionB;

    /**
     * Gamma categories
     */
    private ParameterField num_catT;
    private ParameterField gamma_gainT;
    private ParameterField gamma_lossT;
    private ParameterField gamma_duplicationT;
    
    private JComponent gamma_variationB;

    /**
     * Checkboxes for forbidden categories
     */
    private JCheckBox forbidden_gainCB;
    private JCheckBox forbidden_duplicationCB;

    
    /**
     * Fields for convergence criteria
     */
    private JFormattedTextField epsT;
    private JFormattedTextField roundT;
    private JFormattedTextField rnd_seedT;

    private JRadioButton minCopies0;
    private JRadioButton minCopies1;
    private JRadioButton minCopies2;

    /**
     * Approximate likelihood calculations
     */
    private JCheckBox approximateCB;
    private JFormattedTextField approximate_factorT;
    private JFormattedTextField approximate_thresholdT;
    private JCheckBox approximate_autoCB;
    
    /**
     * Optimization algorithm
     */
    private JRadioButton opt_emRB;
    private JRadioButton opt_ratesRB;
    private JRadioButton opt_distributionRB;
    
    private JFormattedTextField threadsT;
    
    
    /**
     * Fields for gamma parameters
     */
    private ParameterField alpha_lengthT;
    private ParameterField alpha_gainT;
    private ParameterField alpha_lossT;
    private ParameterField alpha_duplicationT;
    
    private List<ParameterField> cat_mod_lenT;
    private List<ParameterField> cat_mod_dupT;
    private List<ParameterField> cat_probT;

    /**
     * Proportion of no-duplication and no-gain families
     */
    private ParameterField forbidden_duplicationT;
    private ParameterField forbidden_gainT;
    
    /**
     * Root distribution info in parameters tab
     */
    private JFormattedTextField root_aT;
    private JFormattedTextField root_bT;
    private JFormattedTextField root_meanT;

    private JCheckBox fixed_rootCB;
    
    private JLabel root_distributionL;
    private JLabel root_aL;
    private JLabel root_bL;
    private JLabel root_meanL;

    /**
     * Lineage-specific parameters
     */
    private LineageSpecificParameters lineage_lengthP;
    private LineageSpecificParameters lineage_duplicationP;
    private LineageSpecificParameters lineage_gainP;
    private LineageSpecificParameters lineage_lossP;
    
    
    private JTabbedPane main_tabbed_pane;
    
    private long init_rnd_seed;
    
    private void initDataStructures()
    {
        this.models = new ArrayList<>();
        this.modelRB = new ArrayList<>();
        this.init_rnd_seed = (new java.util.Random()).nextLong();

       
        IndexedTree main_tree;
        MixedRateModel alt_model=null;
        DataFile<MixedRateModel> alt_file = null; 
        
        Session sesh = app.getActiveSession(); 
        ModelBundle.Entry rates_entry = sesh.getModelBrowser().getSelectedPrimaryEntry(BundleTree.RATES);
        if (rates_entry==null)
        {
        	main_tree = sesh.getModelBrowser().getSelectedPrimaryEntry(BundleTree.TREES).getTreeData().getContent();
        } else
        {
        	alt_file = rates_entry.getRatesData();
            MixedRateModel original_model = alt_file.getContent();
            
            main_tree = original_model.getBaseModel().getTree(); // same tree 
            if (original_model instanceof GammaInvariant)
            {
            	alt_model =((GammaInvariant) original_model).copy();
            } else 
            {
            	assert (original_model instanceof RateVariationModel);
            	alt_model = ((RateVariationModel) original_model).copy();
            }
        }
        
        TreeWithRates null_rates = new TreeWithRates(main_tree, new Random(init_rnd_seed));
//        GammaInvariant null_model = new GammaInvariant(null_rates, 1,1,1,1);
        RateVariationModel null_model = new RateVariationModel(null_rates);
        null_model.initConstantRates();
        addModel(null_model,"Random initial model");
        if (alt_model != null)
            addModel(alt_model,alt_file.getFile().getName());
        
        cat_mod_lenT = new ArrayList<>();
        cat_mod_dupT = new ArrayList<>();
        cat_probT = new ArrayList<>();
    }
    
    private void addModel(MixedRateModel model, String model_name)
    {
        models.add(model);
        JRadioButton rb = new JRadioButton(model_name); // initially unselected
        //rb.setToolTipText(model.getBriefModelDescription());
        modelRB.add(rb);
    }
    
    public RateOptimizationPanel showOptimization()
    {
    	return this.showOptimization(OPT_ROUNDS, OPT_EPS);
    }
    
    /**
     * Modal entrance for letting the user enter data.
     * 
     * @return
     */
    private MixedRateModel askStartingModel()
    {
    	initListeners();
        
        approximateCB.setSelected(true);
        approximate_thresholdT.setValue(12);
        approximate_factorT.setValue(3.0);
        
        // radio button initial choices
        modelRB.get(modelRB.size()-1).doClick();
        opt_ratesRB.doClick();
        
        
        model_dlRB.setEnabled(DEVELOPMENT_CODE);
        model_glRB.setEnabled(DEVELOPMENT_CODE);
        model_plRB.setEnabled(DEVELOPMENT_CODE);
        
        uniform_gainCB.setEnabled(DEVELOPMENT_CODE);
        uniform_duplicationCB.setEnabled(DEVELOPMENT_CODE);
        
        getRootPane().setDefaultButton(startB);
        
        Dimension frameD = app.getSize();

        
        pack();
        setBounds((int)(0.05*frameD.width),(int)(0.05*frameD.height),(int)(0.9*frameD.width),(int)(0.9*frameD.height));
        setVisible(true);

        if (selected_model==-1)
    		return  null;
        
        MixedRateModel starting_model = models.get(selected_model);
        return starting_model;
    }
    
    public RateOptimizationPanel showOptimization(int num_steps, double epsilon)
    {
    	initComponents(true, false);
        roundT.setValue(num_steps);
        epsT.setValue(epsilon);
        
        
        rnd_seedT.addPropertyChangeListener("value", update->updateRandomInitialModel());
        
//        approximateCB.setSelected(false);
//        approximate_thresholdT.setValue(12);
//        approximate_factorT.setValue(3.0);
//        modelRB.get(modelRB.size()-1).doClick();
//        
//        Dimension frameD = app.getSize();
//
//        pack();
//        setBounds((int)(0.05*frameD.width),(int)(0.05*frameD.height),(int)(0.9*frameD.width),(int)(0.9*frameD.height));
//        setVisible(true);
//
//        if (selected_model==-1)
//    		return  null;
        
        MixedRateModel starting_model = askStartingModel();
        if (starting_model==null) return null;
        // starting model has the parameter values copied into already
        
        // name the rate variation model
        File optimized_model_file = null;
        DataFile<AnnotatedTable> data_file = app.getActiveSession().getSelectedData();

        if (selected_model==0) // default null model
        {
            optimized_model_file = new File((File)null, "opt@"+data_file.getFile().getName());
        }
        else
        {
// BUNDLE
        	DataFile<MixedRateModel> rate_file = app.getActiveSession().getModelBrowser().getSelectedRatesEntry().getRatesData();
// BUNDLE
            optimized_model_file = new File((File)null, "opt."+rate_file.getFile().getName());
        }
        DataFile<MixedRateModel> optimized_model = new DataFile<>(starting_model,optimized_model_file);        
        
        // setup optimization
        ML opt;
//        int num_threads = ((Number)threadsT.getValue()).intValue();
//        Count.THREAD_PARALLELISM = num_threads;
        
        TreeWithRates rates = starting_model.getBaseModel();
        IndexedTree phylo = rates.getTree();

        if (opt_ratesRB.isSelected())
        {
        	if (starting_model instanceof GammaInvariant)
        	{
	        	MLGamma opt_gamma = new MLGamma((GammaInvariant)starting_model, data_file.getContent());
		    
		        boolean has_duplication = rates.hasDuplication();
		        boolean has_gain = rates.hasGain();
		        
		        if (has_gain)
		        {
		        	if (alpha_gainT.fixedCB.isSelected())
		        		opt_gamma.fixGainAlpha();
		        	if (!forbidden_gainCB.isSelected() || forbidden_gainT.fixedCB.isSelected())
		        		opt_gamma.fixForbiddenGain();
		        	
		        	if (uniform_gainCB.isSelected() && !lineage_gainP.fixed_parameters.isSelected())
		        	{
		        		opt_gamma.setUniformGain(true);
		        	} else
		        	{
			        	for (int node=0; node<lineage_gainP.numFields(); node++)
			        	{
			        		if (lineage_gainP.getField(node).fixedCB.isSelected())
			        			opt_gamma.fixGain(node);
			        	}
		        	}
		        }
		        if (has_duplication)
		        {
		        	if (alpha_duplicationT.fixedCB.isSelected())
		        		opt_gamma.fixDuplicationAlpha();
		        	if (!forbidden_duplicationCB.isSelected() || forbidden_duplicationT.fixedCB.isSelected())
		        		opt_gamma.fixForbiddenDuplication();
		        	
		        	if (uniform_duplicationCB.isSelected() && !lineage_duplicationP.fixed_parameters.isSelected())
		        	{
		        		opt_gamma.setUniformDuplication(true);
		        	} else
		        	{
			        	for (int node=0; node<lineage_duplicationP.numFields(); node++)
			        	{
			        		if (lineage_duplicationP.getField(node).fixedCB.isSelected())
			        			opt_gamma.fixDuplication(node);
			        	}
		        	}
		        }
		    	// loss
		    	if (alpha_lossT.fixedCB.isSelected())
		    		opt_gamma.fixLossAlpha();
		    	for (int node=0; node<lineage_lossP.numFields(); node++)
		    	{
		    		if (lineage_lossP.getField(node).fixedCB.isSelected())
		    			opt_gamma.fixLoss(node);
		    	}
		    	// length
		    	if (alpha_lengthT.fixedCB.isSelected())
		    		opt_gamma.fixLengthAlpha();
		    	for (int node=0; node<lineage_lengthP.numFields(); node++)
		    		if (lineage_lengthP.getField(node).fixedCB.isSelected())
		    			opt_gamma.fixLength(node);
		    	// root
		    	if (uniform_duplicationCB.isSelected() && uniform_gainCB.isSelected())
		    		opt_gamma.setStationaryRoot(has_gain || has_duplication);
		    	if (fixed_rootCB.isSelected())
		    	{
		    		int root = phylo.getRoot();
		    		opt_gamma.fixGain(root);
		    		opt_gamma.fixLoss(root);
		    		opt_gamma.fixDuplication(root);
		    		opt_gamma.fixLength(root);
		    	}
		    	
		    	opt = opt_gamma;
        	} else
        	{
        		
            	// starting_model is RateVariationModel
        		assert (starting_model instanceof RateVariationModel);
        		
        		MLRateVariation opt_rv = new MLRateVariation((RateVariationModel)starting_model, data_file.getContent());
		    	for (int node=0; node<lineage_lossP.numFields(); node++)
		    	{
	    			opt_rv.fixLoss(node, lineage_lossP.getField(node).fixedCB.isSelected());
		    	}
		    	// duplication
		    	for (int node=0; node<lineage_duplicationP.numFields(); node++)
		    	{
	    			opt_rv.fixDuplication(node, lineage_duplicationP.getField(node).fixedCB.isSelected());
		    	}
		    	// loss
		    	for (int node=0; node<lineage_gainP.numFields(); node++)
		    	{
	    			opt_rv.fixGain(node, lineage_gainP.getField(node).fixedCB.isSelected());
		    	}
            	// root
		    	if (fixed_rootCB.isSelected())
		    	{
		    		opt_rv.fixNodeParameters(starting_model.getBaseModel().getTree().getRoot(), true);
		    	}
		    	
		    	
		    	opt = opt_rv;
        	}
        } else if (opt_distributionRB.isSelected())
        {
        	MLDistribution opt_distr = new MLDistribution(rates, data_file.getContent());
        	for (int node=0; node<lineage_gainP.numFields(); node++)
        	{
        		opt_distr.fixGain(node, lineage_gainP.getField(node).fixedCB.isSelected());
        	}
        	for (int node=0; node<lineage_duplicationP.numFields(); node++)
        	{
    			opt_distr.fixDuplication(node,lineage_duplicationP.getField(node).fixedCB.isSelected());
        	}
    		opt_distr.fixNodeParameters(phylo.getRoot(), fixed_rootCB.isSelected());
        	opt = opt_distr;
        } else
        {
        	assert opt_emRB.isSelected();
        	StraightEM opt_em = new StraightEM(rates, data_file.getContent());
        	for (int node=0; node<lineage_gainP.numFields(); node++)
        	{
        		opt_em.fixNodeParameters(node, lineage_gainP.getField(node).fixedCB.isSelected() ||lineage_duplicationP.getField(node).fixedCB.isSelected() );
        	}
    		opt_em.fixNodeParameters(phylo.getRoot(), fixed_rootCB.isSelected());
        	
        	opt = opt_em;
        }
        
        int min_observed_copies;
    	if (minCopies0.isSelected())
    	{
    		min_observed_copies = 0;
    	}
    	else if (minCopies1.isSelected())
    		min_observed_copies=1;
    	else
    	{
    		assert (minCopies2.isSelected());
    		min_observed_copies=2;
    	}
		opt.setMinimumObservedCopies(min_observed_copies);
    	
        if (approximateCB.isSelected())
        {
            double relative = ((Number)approximate_factorT.getValue()).doubleValue();
            int absolute = ((Number)approximate_thresholdT.getValue()).intValue();
        	boolean auto_truncate = approximate_autoCB.isSelected();
            
            
            if (opt instanceof MLRateVariation)
            {
            	MLRateVariation opt_rv = (MLRateVariation)opt;
            	
        		opt_rv.setWantAutoTruncation(auto_truncate);
            	if (!auto_truncate)
                	opt_rv.setCalculationWidth(absolute, relative);
            } else
            {
            	opt.setCalculationWidth(absolute, relative);
            }
        }
    	
        num_steps = ((Number)roundT.getValue()).intValue();
        epsilon = ((Number)epsT.getValue()).doubleValue();
        
        //System.out.println("#**MSD.sO epsi "+epsilon+"\titer "+num_steps);

        
        if (want_command_line)
        {
        	CommandBuilderDialog cli_dialog = new CommandBuilderDialog(app, opt.getClass(), false, true);

            int nthreads = ((Number)threadsT.getValue()).intValue();
        	cli_dialog.addOption(CommandLine.OPT_THREADS, Integer.toString(nthreads), "Number of CPU threads used in optimization");
        	
        	
        	cli_dialog.addOption(CommandLine.OPT_EPS,Double.toString(epsilon), "Convergence threshold for numerical optimization");
        	cli_dialog.addOption(CommandLine.OPT_ROUNDS, Integer.toString(num_steps), "Iteration threshold for numerical optimization (number of likelihood or gradient calculations)");

        	if (init_rnd_seed != 0L)
        	{
        		cli_dialog.addOption(CommandLine.OPT_RND, Long.toString(init_rnd_seed), "Seed for initializing pseudorandom number generator");
        	}
        	
        	
        	if (approximateCB.isSelected())
        	{
                double relative = ((Number)approximate_factorT.getValue()).doubleValue();
                int absolute = ((Number)approximate_thresholdT.getValue()).intValue();
                boolean auto_truncate = approximate_autoCB.isSelected();
                
                if (opt instanceof MLRateVariation && auto_truncate)
                {
            		cli_dialog.addOption(CommandLine.OPT_TRUNCATE, "auto", "Automatic adjustment of truncated computation");
                } else
                {                
                	cli_dialog.addOption(CommandLine.OPT_TRUNCATE, Integer.toString(absolute)+","+Double.toString(relative), "Truncation parameters for maximum assumed ancestral copies");
                }
        	}
        	cli_dialog.addOption(CommandLine.OPT_MINCOPY, Integer.toString(min_observed_copies), "Minimum observed copies in a family");
        	
        	if (opt instanceof MLGamma)
        	{
        		cli_dialog.addOption(CommandLine.OPT_MODEL_LENGTH_CATEGORIES, Integer.toString(((Number)num_catT.getValue()).intValue()), "Discrete Gamma variation categories for edge length");
        		if (uniform_duplicationCB.isSelected())
        		{
            		cli_dialog.addOption(CommandLine.OPT_MODEL_UNIFORM_DUPLICATION, "true", "Same duplication rate on every edge");
        		}
        		if (uniform_gainCB.isSelected())
        		{
        			cli_dialog.addOption(CommandLine.OPT_MODEL_UNIFORM_GAIN, "true", "Same gain rate on every edge");
        		}
        	}
        	//if (opt instanceof MLGamma || opt instanceof MLDistribution)
        	{
        		count.matek.DiscreteDistribution root_prior = rates.getRootDistribution();
        		cli_dialog.addOption(CommandLine.OPT_MODEL_ROOT_PRIOR, 
        				CommandLine.encodeDistributionOption(root_prior), "Root prior distribution");
        	}
        	
        	if (opt instanceof MLRateVariation)
        	{
        		MixedRateModel sel = models.get(selected_model);
        		int ncat = ((Number) num_catT.getValue()).intValue();
        		if (sel.getNumClasses() != ncat)
        		{
        			cli_dialog.addOption(CommandLine.OPT_MODEL_LENGTH_CATEGORIES, Integer.toString(ncat), "Categories with length variation");
        			cli_dialog.addOption(CommandLine.OPT_MODEL_DUPLICATION_CATEGORIES, Integer.toString(ncat), "Categories with duplication variation");
        		}
        	}
        	
            Session sesh = app.getActiveSession(); 
        	ModelBundle.Entry tree_entry = sesh.getModelBrowser().getSelectedPrimaryEntry(BundleTree.TREES);        	
        	cli_dialog.setTreeData(tree_entry.getTreeData());

        	cli_dialog.setTableData(data_file);
        	
        	
        	if (selected_model==0)
        	{
        		cli_dialog.setRatesData(null);
        	} else
        	{
        		ModelBundle.Entry rates_entry = sesh.getModelBrowser().getSelectedPrimaryEntry(BundleTree.RATES);
        		cli_dialog.setRatesData(rates_entry.getRatesData());
        	}
        	
        	
        	
        	cli_dialog.setVisible(true);

        	return null;
        } else
        {
	    	RateOptimizationPanel opt_panel = new RateOptimizationPanel(optimized_model, opt);
	    	opt_panel.setDescendantModel(selected_model != 0);
	    	// optimization parameters
	        
	    	opt_panel.launchOptimization(app, epsilon, num_steps);
	        return opt_panel;
        }
    }
    
    public PosteriorsView showPosteriors()
    {
    	initComponents(false, false);
    	MixedRateModel rates_model = askStartingModel();
    	if (rates_model==null) return null;
    	
    	File model_file;
        if (selected_model==0) // default null model
        {
            model_file = new File((File)null, "default");
        }
        else
        {
// BUNDLE
//            DataFile<GammaInvariant> rate_file = app.getActiveSession().getRatesBrowser().getSelectedPrimaryItem().getDataFile();
            DataFile<MixedRateModel> rate_file = app.getActiveSession().getModelBrowser().getSelectedRatesEntry().getRatesData();
// BUNDLE
            model_file = new File((File)null, rate_file.getFile().getName());
        }
        DataFile<MixedRateModel> data_model = new DataFile<>(rates_model,model_file);        
        DataFile<AnnotatedTable> data_table = app.getActiveSession().getSelectedData();
    	PosteriorsView P = new PosteriorsView(data_model, data_table);
    	if (approximateCB.isSelected())
    	{
            double relative = ((Number)approximate_factorT.getValue()).doubleValue();
            int absolute = ((Number)approximate_thresholdT.getValue()).intValue();
            P.setCalculcationWidthThresholds(absolute, relative);
    	} 
    	return P;
    }
    
    public SimulationView showSimulation()
    {
    	initComponents(false, true);
    	MixedRateModel rates_model = askStartingModel();
    	if (rates_model==null)
    		return null;
    	int num_families = ((Number)roundT.getValue()).intValue();
    	int min_copies=0;
    	if (minCopies1.isSelected())
    		min_copies = 1;
    	else if (minCopies2.isSelected())
    		min_copies = 2;
    	
    	long rnd_seed = ((Number)rnd_seedT.getValue()).longValue(); // ((Number)epsT.getValue()).longValue();
    	
    	File model_file;
        if (selected_model==0) // default null model
        {
            model_file = new File((File)null, "default");
        }
        else
        {
// BUNDLE
//            DataFile<GammaInvariant> rate_file = app.getActiveSession().getRatesBrowser().getSelectedPrimaryItem().getDataFile();
            DataFile<MixedRateModel> rate_file = app.getActiveSession().getModelBrowser().getSelectedRatesEntry().getRatesData();
// BUNDLE
            model_file = new File((File)null, rate_file.getFile().getName());
        }
        DataFile<MixedRateModel> data_model = new DataFile<>(rates_model,model_file);        
    	SimulationView S = new SimulationView(data_model, rnd_seed, num_families, min_copies);
    	
//    	S.computeAll();
    	return S;
    }
    
    
    private void initComponents(boolean want_optimization, boolean want_simulation)
    {
        // model structure tab
        Box structureB = new Box(BoxLayout.PAGE_AXIS);
        JComponent initial_modelB = createInitialModelBox();
        modelRB.get(0).setVisible(want_optimization);

        JComponent model_typeB = createModelTypeBox();
        
        JComponent lineage_variationB = createLineageVariationBox();
        JComponent root_priorB = createRootPriorBox();
        JComponent family_variationB = createFamilyVariationBox();
        
        // model parameters tab
        Box parametersB = new Box(BoxLayout.PAGE_AXIS);
        
        if (want_optimization)
        {
        	structureB.add(initial_modelB);
        	structureB.add(model_typeB);
            structureB.add(lineage_variationB);
            structureB.add(root_priorB);
            structureB.add(family_variationB);
            structureB.add(createAlgorithmBox(want_optimization, want_simulation));
        }
        else
        {
        	parametersB.add(createAlgorithmBox(want_optimization, want_simulation));
        	parametersB.add(initial_modelB);
        	parametersB.add(family_variationB);
        }
        
        // parametersB.add(
        gamma_variationB = createGammaFamilyVariationParameterBox()
        //)
        ;
        
        
        parametersB.add(createRootPriorParameterBox());
        parametersB.add(createLineageSpecificParameterBox());

        JTabbedPane tabs = new JTabbedPane();
        if (want_optimization)
        {
	        tabs.addTab("Model type & algorithm", structureB);
	    	tabs.addTab("Model parameters", new JScrollPane(parametersB));
        } else
        {
	    	tabs.addTab("Algorithm & model parameters", new JScrollPane(parametersB));
        }
        // last tab : 
        
        // TODO 
        tabs.addTab("Rate variation", createRateVariationModelParameterBox(0));
        
        this.main_tabbed_pane = tabs;

        setLayout(new BorderLayout());
        add(tabs, BorderLayout.CENTER);
        add(createButtonPanel(),BorderLayout.PAGE_END);
        if (!want_optimization)
        	cliB.setVisible(false);
        
    }
    
    private void initListeners()
    {
    	// start & cancel
        cancelB.addActionListener(click -> 
        {
            selected_model = -1;
            setVisible(false);
            dispose();
        });

        startB.addActionListener(click ->
        {
            copyFieldValuesIntoModel();
            want_command_line = false;
            setVisible(false);
            dispose();
        });
        cliB.addActionListener(click -> 
        {
        	copyFieldValuesIntoModel();
        	want_command_line = true;
            setVisible(false);
        	dispose();
        });
        
        // initial model listeners
        for (int idx=0; idx<modelRB.size(); ++idx)
        {
        	final int select_idx = idx;
            JRadioButton rb = modelRB.get(idx);
            rb.addItemListener(item_event ->
            {
            	if (item_event.getStateChange() == ItemEvent.SELECTED)
                    setSelectedModel(select_idx);
            });
        }
        
        // synchronize model structure
        JRadioButton[] structure_buttons 
        = {model_gldRB, model_dlRB, model_glRB, model_plRB };
        for (JRadioButton rb : structure_buttons)
		{
			rb.addActionListener(click->synchronizeModelStructureStates());
		}
        
        for (JRadioButton b: root_distributionB.values())
            b.addActionListener(click->synchronizeRootPriorStates());
        
        root_aT.addPropertyChangeListener("value", update->synchronizeRootPriorMean());
        root_bT.addPropertyChangeListener("value", update->synchronizeRootPriorMean());
        
        num_catT.addPropertyChangeListener("value", update->synchronizeOptimizationOptions());
        gamma_duplicationT.addPropertyChangeListener("value", update->synchronizeOptimizationOptions());
        gamma_gainT.addPropertyChangeListener("value", update->synchronizeOptimizationOptions());

        // variation across families

        alpha_gainT.listenTo(gamma_gainT, 1, "gain rates do not vary across families");
        alpha_duplicationT.listenTo(gamma_duplicationT, 1, "duplication rates do not vary across families");
        alpha_lossT.listenTo(gamma_lossT, 1, "loss rates do not vary across families");
        alpha_lengthT.listenTo(num_catT, 1, "edge lengths do not vary across families");


        forbidden_duplicationT.listenTo(forbidden_duplicationCB, true, "there is no forbidden-duplication category");
        forbidden_gainT.listenTo(forbidden_gainCB, true, "there is no forbidden-gain category");        
        
        // uniform rates 
        uniform_duplicationCB.addItemListener(e->lineage_duplicationP.disableVariationIfSelected(uniform_duplicationCB));
        uniform_gainCB.addItemListener(e->lineage_gainP.disableVariationIfSelected(uniform_gainCB));
        uniform_duplicationCB.addItemListener(e->
	        {
	        	if (uniform_duplicationCB.isSelected() && uniform_gainCB.isSelected())
	        		synchronizeModelStructureStates();
	        	else
	        	{
	        		synchronizeRootPriorStates();
	        		synchronizeOptimizationOptions();
	        	}
	        });
        uniform_gainCB.addItemListener(e->
        {
        	if (uniform_duplicationCB.isSelected() && uniform_gainCB.isSelected())
        		synchronizeModelStructureStates();
        	else
        	{
        		synchronizeRootPriorStates();
        		synchronizeOptimizationOptions();
        	}
        });
        
        forbidden_gainCB.addItemListener(change->synchronizeOptimizationOptions());
        forbidden_duplicationCB.addItemListener(change->synchronizeOptimizationOptions());
        
        
        // if both are selected, set root prior, and root params; disable root prior selection

//        lineage_duplicationP.disableVariationIfSelected(uniform_duplicationCB, false);
//        lineage_gainP.disableVariationIfSelected(uniform_gainCB, false);
//
//        void listenToEqualizer(JCheckBox cb, LineageSpecificParameters lineage_lengthP)
//        {
//        	cb.addItemListener(e->
//        		{
//        			synchronizeTo(cb, false); 
//        			if (cb.isSelected()) 
//        				equalizeValues(lineage_lengthP);
//        		});
//        }
        
        num_catT.addPropertyChangeListener("value", update->this.setRateVariationModelFields());
    
    }
    
    private boolean hasRateVariation()
    {
    	
    	boolean varLength = num_catT.getValue()!=null && 1<num_catT.intValue();
    	boolean varGain = (gamma_gainT.getValue()!=null && 1<gamma_gainT.intValue()) || forbidden_gainCB.isSelected() ;
    	boolean varDuplication = (gamma_duplicationT.getValue()!=null && 1<gamma_duplicationT.intValue()) || forbidden_duplicationCB.isSelected();
    	boolean varLoss = (gamma_lossT.getValue()!=null && 1<gamma_lossT.intValue());
    	
       	boolean hasRateVariation = varLength || varGain || varDuplication || varLoss;
        
       	return hasRateVariation;
    	
    }
    
//    private boolean hasUniformRates()
//    {
//    	boolean uniformGain = uniform_gainCB.isSelected() && uniform_gainCB.isEnabled(); // disabled if model has no gain
//    	boolean uniformDuplication = uniform_duplicationCB.isSelected() && uniform_duplicationCB.isEnabled(); // disabled if model has no duplication
//    	
//    	boolean hasUniformRates = uniformGain || uniformDuplication;
//    	return hasUniformRates;
//    }
    
    
    /**
     * Copies the UI field values into the selected model (entry in {@link #models}),
     * and initializes the likelihood optimization.
     */
    private void copyFieldValuesIntoModel()
    {
    	if (selected_model!=-1)
    	{
    		assert (modelRB.get(selected_model).isSelected());
    		MixedRateModel model=models.get(selected_model);
    		TreeWithRates base_rates = model.getBaseModel();
    		for (int node=0; node<lineage_lengthP.numFields(); node++)
    		{
    			double len = lineage_lengthP.getField(node).doubleValue();
    			double grate = lineage_gainP.getField(node).doubleValue();
    			double lrate = lineage_lossP.getField(node).doubleValue();
    			double drate = lineage_duplicationP.getField(node).doubleValue();
    			
    			base_rates.setRates(node, len, grate, lrate, drate);
//    			base_rates.setEdgeLength(node, len);
//    			base_rates.setGainRate(node, grate);
//    			base_rates.setLossRate(node, lrate);
//    			base_rates.setDuplicationRate(node, drate);
    		}
    		// set root distribution
    		Class<? extends DiscreteDistribution> prior_class = getSelectedRootPriorClass();
    		
    		DiscreteDistribution root_prior;
    		
    		//int root = base_rates.getTree().getRoot();
    		assert (prior_class != null);
    		if (NegativeBinomial.class.equals(prior_class))
    		{
                double κ = ((Number)root_aT.getValue()).doubleValue();
                double q = ((Number)root_bT.getValue()).doubleValue();
                root_prior = new NegativeBinomial(κ,q);
                //base_rates.setParameters(root, κ, 1.0, q);   // loss prob=1              
    		} else if (Poisson.class.equals(prior_class))
    		{
    			double r = ((Number)root_aT.getValue()).doubleValue();
    			root_prior = new Poisson(r);
	        	//base_rates.setParameters(root, r, 1.0, 0.0); // loss prob=1 
    		} else if (ShiftedGeometric.class.equals(prior_class))
    		{
                double p0 = ((Number)root_aT.getValue()).doubleValue();
                double q = ((Number)root_bT.getValue()).doubleValue();
                root_prior = new ShiftedGeometric(p0, q);
	        	//base_rates.setParameters(root, 0.0, p0, q);
    		} else 
    		{
    			assert (PointDistribution.class.equals(prior_class)); // no other choice
                double p0 = ((Number)root_aT.getValue()).doubleValue();
                root_prior = new PointDistribution(p0);
	        	//base_rates.setParameters(root, 0.0, p0, 0.0);
    		}
    		base_rates.setRootDistribution(root_prior);
    		
    		// 
    		
    		if (model instanceof GammaInvariant)
    		{
    			GammaInvariant gamma_model = (GammaInvariant) model;
	            // rate variation across families
    			gamma_model.setGainForbidden(forbidden_gainT.doubleValue());
    			gamma_model.setDuplicationForbidden(forbidden_duplicationT.doubleValue());
	
	            int cat_length = num_catT.intValue();
	            int cat_gain = gamma_gainT.intValue();
	            int cat_loss = gamma_lossT.intValue();
	            int cat_duplication = gamma_duplicationT.intValue();
	            gamma_model.setClasses(cat_gain, cat_loss, cat_duplication, cat_length);
	
	            gamma_model.setLengthAlpha(alpha_lengthT.doubleValue());
	            gamma_model.setGainAlpha(alpha_gainT.doubleValue());
	            gamma_model.setLossAlpha(alpha_lossT.doubleValue());
	            gamma_model.setDuplicationAlpha(alpha_duplicationT.doubleValue());
    		} else
    		{
    			assert (model instanceof RateVariationModel);
    			
    			RateVariationModel variation_model = (RateVariationModel) model;
    			
    			int ncat = num_catT.intValue();
				double[] cat_prob = new double[ncat];
				double[] cat_modlen = new double[ncat];
				double[] cat_moddup = new double[ncat];
				
				for (int k=0; k<ncat; k++)
				{
					cat_prob[k] = cat_probT.get(k).doubleValue();
					cat_modlen[k] = cat_mod_lenT.get(k).doubleValue();
					cat_moddup[k] = cat_mod_dupT.get(k).doubleValue();
				}
    			
    			if (ncat == model.getNumClasses())
    			{
    				// reset existing category parameters 
    				for (int k=0; k<ncat; k++)
    				{
    					RateVariationModel.Category C = variation_model.getCategory(k);
    					C.setLogCatProbability(Math.log(cat_prob[k]));
    					C.setModifiers(cat_modlen[k], cat_moddup[k]);
    				}
    			} else
    			{
    				// reinitialize 
    				variation_model.initCategories(cat_prob, cat_modlen, cat_moddup);
    			}
    		}
    	}
    }
    
    private Class<? extends DiscreteDistribution> getSelectedRootPriorClass()
    {
    	for (Map.Entry<Class<? extends DiscreteDistribution>, JRadioButton> E: root_distributionB.entrySet())
    	{
    		if (E.getValue().isSelected())
    		{
    			//System.out.println("#**MSD.gSRPC "+E+"\t"+Arrays.toString(Thread.currentThread().getStackTrace()));
    			return E.getKey();
    		}
    	}
    	return null; 
    }
    
    
    /**
     * Panel at the bottom of the dialog with Start and Cancel buttons
     *
     * @return a container with these buttons, and initializes {@link #startB} and {@link #cancelB}
     */
    private JComponent createButtonPanel()
    {
        startB = new JButton("OK");
        cancelB = new JButton("Cancel");
        cliB = new JButton("Command-line interface");
        Box button_box = new Box(BoxLayout.LINE_AXIS);
        button_box.add(Box.createHorizontalGlue());
        button_box.add(cancelB);
        button_box.add(Box.createRigidArea(new Dimension(10,0)));
        button_box.add(cliB);
        button_box.add(Box.createRigidArea(new Dimension(10,0)));
        
        button_box.add(startB);
        
        
        
        button_box.setBorder(BorderFactory.createEmptyBorder(0, 10, 10, 10));
        return button_box;
    }
    
    /**
     * Panel with radio buttons for possible model choices.
     * Puts the created buttons into a ButtonGroup.
     * All radio buttons are initially unselected.
     *
     * @return a container
     */
    private JComponent createInitialModelBox()
    {
        JPanel starting_modelP = new JPanel();
        BoxLayout starting_modelL = new BoxLayout(starting_modelP, BoxLayout.LINE_AXIS);
        starting_modelP.setLayout(starting_modelL);
        starting_modelP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Starting model",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));


        ButtonGroup starting_modelG = new ButtonGroup();
        for (JRadioButton rb: modelRB)
        {
            starting_modelG.add(rb);
            starting_modelP.add(rb);
        }
        starting_modelP.add(Box.createHorizontalGlue());
        starting_modelP.setMaximumSize(new Dimension(Integer.MAX_VALUE,getPreferredSize().height));
        
        return new JScrollPane(starting_modelP);
    }


    /**
     * Sets up the model type radio buttons
     * ({@link #model_gldRB}, {@link #model_dlRB},{@link #model_glRB}, {@link #model_plRB}),
     * and puts them in a ButtonGroup.
     * All radio buttons are initially unselected.
     *
     * @return a container with the buttons
     */
    private JComponent createModelTypeBox()
    {
        model_gldRB = new JRadioButton("Gain-loss-duplication"); // (Cs\u0171r\u00f6s & Mikl\u00f3s)");
        model_gldRB.setActionCommand("GLD");

        model_dlRB = new JRadioButton("Duplication-loss");
        model_dlRB.setActionCommand("DL");

        //dlEB = new JRadioButton("Duplication-loss with equal rates (Hahn et al.)");
        //dlEB.setActionCommand("DL=");
        //dlEB.addActionListener(this);
        model_glRB = new JRadioButton("Gain-loss");
        model_glRB.setActionCommand("GL");

        model_plRB = new JRadioButton("Pure loss");
        model_plRB.setActionCommand("PL");

        ButtonGroup model_typeG = new ButtonGroup();
        model_typeG.add(model_gldRB);
        model_typeG.add(model_dlRB);
        model_typeG.add(model_glRB);
        model_typeG.add(model_plRB);

        JPanel model_typeP = new JPanel();
        model_typeP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Model type",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
        {
            GridBagLayout layout = new GridBagLayout();
            GridBagConstraints gc = new GridBagConstraints();
            model_typeP.setLayout(layout);

            gc.fill = GridBagConstraints.HORIZONTAL;
            gc.gridx = 0;
            gc.gridy = 0;
            gc.gridheight=1;
            gc.gridwidth=1;
            gc.anchor = GridBagConstraints.WEST;
            gc.weightx = 0.1;
            gc.weighty = 0.1;

            layout.setConstraints(model_gldRB,gc);
            model_typeP.add(model_gldRB);

            gc.gridy++;
            layout.setConstraints(model_dlRB,gc);
            model_typeP.add(model_dlRB);

            gc.gridx++;
            layout.setConstraints(model_glRB,gc);
            model_typeP.add(model_glRB);

            gc.gridx=0;
            gc.gridy++;

            //gc.gridx++;
            layout.setConstraints(model_plRB,gc);
            model_typeP.add(model_plRB);

            gc.gridheight=gc.gridy+1;
            gc.gridy=0;
            gc.gridx++;
            gc.weightx=1.0;
            gc.fill=GridBagConstraints.BOTH;
            Component fillerX = Box.createHorizontalGlue();
            layout.setConstraints(fillerX,gc);
            model_typeP.add(fillerX);
        }
        return model_typeP;
    }
    
    /**
     * Panel with uniform rate selection check boxes. Initializes the check boxes
     * {@link #uniform_duplicationCB}, {@link #uniform_gainCB}, {@link #uniform_lossCB}.
     * all check boxes are initially unselected.
     *
     * @return a container
     */
    private JComponent createLineageVariationBox()
    {
        uniform_duplicationCB = new JCheckBox("Same duplication-loss ratio in all lineages");
        uniform_gainCB = new JCheckBox("Same gain rate in all lineages");
        uniform_lossCB = new JCheckBox("Same loss rate in all lineages (fixed edge lengths)");

        JPanel lineage_variationP = new JPanel();
        lineage_variationP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Lineage-specific variation",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
        {
            GridBagLayout layout = new GridBagLayout();
            GridBagConstraints gc = new GridBagConstraints();
            lineage_variationP.setLayout(layout);

            gc.fill = GridBagConstraints.HORIZONTAL;
            gc.gridx = 0;
            gc.gridy = 0;
            gc.gridheight=1;
            gc.gridwidth=1;
            gc.anchor = GridBagConstraints.WEST;
            gc.weightx = 0.1;
            gc.weighty = 0.1;

            layout.setConstraints(uniform_gainCB,gc);
            lineage_variationP.add(uniform_gainCB);

            gc.gridx++;
            layout.setConstraints(uniform_duplicationCB,gc);
            lineage_variationP.add(uniform_duplicationCB);

            gc.gridx++;
            layout.setConstraints(uniform_lossCB,gc);
            lineage_variationP.add(uniform_lossCB);

            gc.gridheight=gc.gridy+1;
            gc.gridx++;
            gc.weightx=1.0;
            gc.fill=GridBagConstraints.BOTH;
            Component fillerX = Box.createHorizontalGlue();
            layout.setConstraints(fillerX,gc);
            lineage_variationP.add(fillerX);
        }
        return lineage_variationP;
    }
    
    /**
     * Panel with radio buttons for root prior distribution. Initializes the
     * radio buttons {@link #root_distributionB}.
     * All radio buttons are initially unselected.
     *
     * @return a container
     */
    private JComponent createRootPriorBox()
    {
        root_distributionB = new HashMap<>();

        JPanel rootP = new JPanel();
        rootP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Family size distribution at root",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
        GridBagLayout layout = new GridBagLayout();
        GridBagConstraints gc = new GridBagConstraints();
        rootP.setLayout(layout);
        gc.fill = GridBagConstraints.HORIZONTAL;
        gc.gridy = 0;
        gc.gridx = 0;
        gc.gridheight=1;
        gc.gridwidth=1;
        gc.anchor = GridBagConstraints.WEST;
        gc.weightx = 0.1;
        gc.weighty = 0.1;

        ButtonGroup rootG = new ButtonGroup();

        JRadioButton polyaB = new JRadioButton("Negative binomial (Pólya)");
        root_distributionB.put(NegativeBinomial.class, polyaB);
        layout.setConstraints(polyaB, gc);
        rootG.add(polyaB);
        rootP.add(polyaB);
        ++gc.gridx;
        
        JRadioButton poissonB = new JRadioButton("Poisson");
        root_distributionB.put(Poisson.class, poissonB);
        layout.setConstraints(poissonB, gc);
        rootG.add(poissonB);
        rootP.add(poissonB);
        ++gc.gridx;
        
        JRadioButton geomB = new JRadioButton("Shifted geometric");
        root_distributionB.put(ShiftedGeometric.class, geomB);
        layout.setConstraints(geomB, gc);
        rootG.add(geomB);
        rootP.add(geomB);
        ++gc.gridx;
        
        JRadioButton pointB = new JRadioButton("Bernoulli");
        root_distributionB.put(PointDistribution.class, pointB);
        layout.setConstraints(pointB, gc);
        rootG.add(pointB);
        rootP.add(pointB);
        ++gc.gridx;

//        JRadioButton stationaryB = new JRadioButton("Stationary");
//        root_distributionB.put(DiscreteDistribution.class, stationaryB);
//        layout.setConstraints(stationaryB, gc);
//        rootG.add(stationaryB);
//        rootP.add(stationaryB);
//        ++gc.gridx;

        gc.gridheight=gc.gridy+1;
        gc.gridx++;
        gc.weightx=1.0;
        gc.fill=GridBagConstraints.BOTH;
        Component fillerX = Box.createHorizontalGlue();
        layout.setConstraints(fillerX,gc);
        rootP.add(fillerX);

        return rootP;
    }
    
    /**
     * Panel with fields and check boxes for selecting rate variation across families.
     * Initializes {@link #gamma_duplicationT}, {@link #gamma_gainT}, {@link #num_catT}
     * and {@link #gamma_lossT}, {@link #forbidden_duplicationCB}, {@link #forbidden_gainCB}.
     * Initally, all check boxes are unselected, and all text field values are unset.
     *
     * Loss variation is disabled for now (not editable, not enabled, not visible).
     *
     * @return a container
     */
    private JComponent createFamilyVariationBox()
    {
        num_catT = new ParameterField(NumberFormat.getIntegerInstance());
        num_catT.setInputVerifier(InputVerifiers.getPositiveInstance());
        num_catT.setColumns(2);
        num_catT.setToolTipText("Number of categories  (1 means variation disabled)");

        num_catT.setEditable(DEVELOPMENT_CODE);
        //num_catT.setEnabled(DEVELOPMENT_CODE);
        
        
        //num_catT.setFocusLostBehavior(JFormattedTextField.COMMIT_OR_REVERT);
        gamma_gainT = new ParameterField(NumberFormat.getIntegerInstance());
        gamma_gainT.setInputVerifier(InputVerifiers.getPositiveInstance());
        gamma_gainT.setColumns(2);
        gamma_gainT.setToolTipText("Number of discrete Gamma categories for gain rate (1 means variation disabled)");

        gamma_duplicationT = new ParameterField(NumberFormat.getIntegerInstance());
        gamma_duplicationT.setInputVerifier(InputVerifiers.getPositiveInstance());
        gamma_duplicationT.setColumns(2);
        gamma_duplicationT.setToolTipText("Number of discrete Gamma categories for duplication rate (1 means variation disabled)");

        gamma_lossT = new ParameterField(NumberFormat.getIntegerInstance());
        gamma_lossT.setInputVerifier(InputVerifiers.getPositiveInstance());
        gamma_lossT.setColumns(2);
        gamma_lossT.setToolTipText("Number of discrete Gamma categories for loss rate (1 means variation disabled)");

        JLabel gamma_lengthL = new JLabel("");
        JLabel gamma_gainL = new JLabel("Gamma categories");
        JLabel gamma_duplicationL = new JLabel("Gamma categories");
        JLabel gamma_lossL = new JLabel("Gamma categories");

        JLabel variation_lengthL = new JLabel("Number of categories");
        JLabel variation_lossL = new JLabel("Loss rate");
        JLabel variation_gainL = new JLabel("Gain rate");
        JLabel variation_duplicationL = new JLabel("Duplication rate");

        gamma_lossT.setEditable(false);
        gamma_lossL.setVisible(false);
        gamma_lossT.setVisible(false);
        variation_lossL.setVisible(false);


        forbidden_gainCB = new JCheckBox("No-gain category");
        forbidden_duplicationCB = new JCheckBox("No-duplication category");
        
        
        gamma_gainT.setVisible(WANT_GENERAL_VARIATION);
        gamma_gainL.setVisible(gamma_gainT.isVisible());
        gamma_duplicationT.setVisible(WANT_GENERAL_VARIATION);
        gamma_duplicationL.setVisible(gamma_duplicationT.isVisible());
        forbidden_gainCB.setVisible(WANT_GENERAL_VARIATION);
        forbidden_duplicationCB.setVisible(WANT_GENERAL_VARIATION);
        variation_gainL.setVisible(WANT_GENERAL_VARIATION);
        variation_duplicationL.setVisible(WANT_GENERAL_VARIATION);

        JPanel family_variationP = new JPanel();
        family_variationP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Rate variation across families",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
        {
            GridBagLayout layout = new GridBagLayout();
            GridBagConstraints gc = new GridBagConstraints();
            family_variationP.setLayout(layout);

            gc.fill = GridBagConstraints.HORIZONTAL;
            gc.gridy = 0;
            gc.gridx = 0;
            gc.gridheight=1;
            gc.gridwidth=1;
            gc.anchor = GridBagConstraints.WEST;
            gc.weightx = 0.1;
            gc.weighty = 0.1;

            gc.fill=GridBagConstraints.NONE;
            layout.setConstraints(variation_lengthL,gc);
            family_variationP.add(variation_lengthL);
            gc.gridx++;
            gc.anchor = GridBagConstraints.EAST;
            gc.fill=GridBagConstraints.NONE;
            layout.setConstraints(num_catT, gc);
            family_variationP.add(num_catT);
            gc.gridx++;
            gc.anchor = GridBagConstraints.WEST;
            layout.setConstraints(gamma_lengthL,gc);
            family_variationP.add(gamma_lengthL);


            gc.gridx=0;
            gc.gridy++;
            gc.fill = GridBagConstraints.NONE;
            layout.setConstraints(variation_lossL, gc);
            family_variationP.add(variation_lossL);
            gc.gridx++;
            gc.anchor = GridBagConstraints.EAST;
            gc.fill = GridBagConstraints.NONE;
            layout.setConstraints(gamma_lossT, gc);
            family_variationP.add(gamma_lossT);
            gc.gridx++;
            gc.anchor = GridBagConstraints.WEST;
            layout.setConstraints(gamma_lossL,gc);
            family_variationP.add(gamma_lossL);

            gc.gridx=0;
            gc.gridy++;
            gc.fill = GridBagConstraints.NONE;
            layout.setConstraints(variation_gainL, gc);
            family_variationP.add(variation_gainL);
            gc.gridx++;
            gc.anchor = GridBagConstraints.EAST;
            gc.fill = GridBagConstraints.NONE;
            layout.setConstraints(gamma_gainT, gc);
            family_variationP.add(gamma_gainT);
            gc.gridx++;
            gc.anchor = GridBagConstraints.WEST;
            layout.setConstraints(gamma_gainL,gc);
            family_variationP.add(gamma_gainL);

            gc.gridx++;
            layout.setConstraints(forbidden_gainCB,gc);
            family_variationP.add(forbidden_gainCB);


            //gc.gridx++;
            //layout.setConstraints(forbidden_lossCB,gc);
            //family_variationP.add(forbidden_lossCB);

            gc.gridx=0;
            gc.gridy++;
            gc.fill=GridBagConstraints.NONE;
            layout.setConstraints(variation_duplicationL, gc);
            family_variationP.add(variation_duplicationL);
            gc.gridx++;
            gc.anchor = GridBagConstraints.EAST;
            gc.fill = GridBagConstraints.NONE;
            layout.setConstraints(gamma_duplicationT, gc);
            family_variationP.add(gamma_duplicationT);
            gc.gridx++;
            gc.anchor = GridBagConstraints.WEST;
            layout.setConstraints(gamma_duplicationL,gc);
            family_variationP.add(gamma_duplicationL);

            gc.gridx++;
            layout.setConstraints(forbidden_duplicationCB,gc);
            family_variationP.add(forbidden_duplicationCB);

            gc.gridheight=gc.gridy+1;
            gc.gridx++;
            gc.gridy=0;
            gc.weightx=1.0;
            gc.fill=GridBagConstraints.BOTH;
            Component fillerX = Box.createHorizontalGlue();
            layout.setConstraints(fillerX,gc);
            family_variationP.add(fillerX);
        }
        return family_variationP;
    }
    /**
     * Panel with fields for convergence criteria.
     * Initializes {@link #epsT} and {@link #roundT}.
     * Initially, field values are unset.
     *
     * @return a container
     */
    private JComponent createAlgorithmBox(boolean want_optimization, boolean want_simulation)
    {
// BUNDLE
//		AnnotatedTablePanel selected_table = app.getActiveSession().getDataBrowser().getSelectedPrimaryItem();
//		DataFile<AnnotatedTable> data_file=(selected_table==null?null:selected_table.getDataFile());
    	ModelBundle.Entry selected_table = app.getActiveSession().getDataBrowser().getSelectedTableEntry();
    	DataFile<AnnotatedTable> data_file = selected_table==null?null:selected_table.getTableData();
// BUNDLE

    	JLabel epsL, rnd_seedL;
		NumberFormat rndF =  NumberFormat.getIntegerInstance();
        rnd_seedT = new JFormattedTextField(rndF);
        rnd_seedT.setColumns(16);
        rnd_seedL = new JLabel("Random seed");   
        rnd_seedL.setLabelFor(rnd_seedT);
    	
    	if (want_simulation)
    	{
    		NumberFormat epsF =  NumberFormat.getIntegerInstance();
            epsT = new JFormattedTextField(epsF);
            //epsT.addPropertyChangeListener("value",this);
            epsT.setColumns(16);
            epsL = new JLabel("Random seed");
    	} else
    	{
    		NumberFormat epsF =  new DecimalFormat("0.############E0");
            epsF.setParseIntegerOnly(false);
            epsT = new JFormattedTextField(epsF);
            //epsT.addPropertyChangeListener("value",this);
            epsT.setInputVerifier(InputVerifiers.getNonnegativeInstance());
            epsT.setColumns(16);
            epsL = new JLabel("Convergence threshold on the likelihood");
    	}
        epsL.setLabelFor(epsT);
        
        
        

        NumberFormat roundF = NumberFormat.getIntegerInstance();
        roundT = new JFormattedTextField(roundF);
        //roundT.addPropertyChangeListener("value",this);
        roundT.setInputVerifier(InputVerifiers.getNonnegativeInstance());
        roundT.setColumns(10);
        JLabel roundL = new JLabel(
        		want_simulation?"Simulated families":
        		"Maximum number of optimization rounds");
        roundL.setLabelFor(roundT);
        
        rnd_seedT.setValue(init_rnd_seed);
        if (want_simulation)
        {
        	if (data_file == null)
        		roundT.setValue(OPT_NUM_FAMILIES);
        	else
        		roundT.setValue(data_file.getContent().getFamilyCount());
        	epsT.setValue(init_rnd_seed);
        }

        minCopies0 = new JRadioButton("0");
        minCopies1 = new JRadioButton("1");
        minCopies2 = new JRadioButton("2");
        
        int min_copies;
        if (want_simulation)
        	min_copies = 2;
        else
    	{
    		min_copies = data_file.getContent().minCopies();
    	}
        ButtonGroup minCopiesG = new ButtonGroup();
        minCopiesG.add(minCopies0);
        minCopiesG.add(minCopies1);
        minCopiesG.add(minCopies2);
        minCopies1.setEnabled(min_copies>0);
        minCopies2.setEnabled(min_copies>1);
        JLabel minCopiesL = new JLabel("Family observation threshold on minimum copies");
        //= new JLabel("Minimum total number of copies for likelihood correction");
        Box minCopiesB = new Box(BoxLayout.LINE_AXIS);
        minCopiesB.add(minCopiesL);
        minCopiesB.add(minCopies0);
        minCopiesB.add(minCopies1);
        minCopiesB.add(minCopies2);
        if (min_copies==0)
        	minCopies0.setSelected(true);
        else if (min_copies==1 || want_simulation)
        	minCopies1.setSelected(true);
        else
        {
        	assert (min_copies>=2);
        	minCopies2.setSelected(true);
        }
        
//        NumberFormat bracketF = new DecimalFormat("0.############");
//        bracketF.setParseIntegerOnly(false);
//        bracketT = new JFormattedTextField(bracketF);
//        bracketT.setInputVerifier(InputVerifiers.getPositiveInstance());
//        bracketT.setColumns(16);
//        JLabel bracketL = new JLabel("Bracketing for single-parameter optimization");
//        bracketL.setLabelFor(bracketT);

        NumberFormat approximate_likelihood_factorF = new DecimalFormat("0.0###########");
        approximate_likelihood_factorF.setParseIntegerOnly(false);
        approximate_factorT = new JFormattedTextField(approximate_likelihood_factorF);
        approximate_factorT.setColumns(16);
        JLabel approximate_factorL = new JLabel("Truncation: relative threshold");
        approximate_factorL.setLabelFor(approximate_factorT);

        NumberFormat approximate_likelihood_thresholdF = NumberFormat.getIntegerInstance();
        approximate_thresholdT = new JFormattedTextField(approximate_likelihood_thresholdF);
        approximate_thresholdT.setInputVerifier(InputVerifiers.getNonnegativeInstance());
        approximate_thresholdT.setColumns(10);
        JLabel approximate_thresholdL = new JLabel("Truncation: absolute threshold");
        approximate_thresholdL.setLabelFor(approximate_thresholdT);

        approximate_autoCB = new JCheckBox("Automatic");
        approximate_autoCB.setToolTipText("Truncation is set automatically for desired precision (faster computation)");
        approximate_autoCB.addChangeListener(e->
	        {
	        	boolean approxOK = approximateCB.isSelected();
	        	boolean autoOK = approxOK && approximate_autoCB.isSelected();
                approximate_thresholdT.setEnabled(approxOK && !autoOK);
                approximate_thresholdT.setEditable(approxOK && !autoOK);
                approximate_factorT.setEnabled(approxOK && !autoOK);
                approximate_factorT.setEditable(approxOK && !autoOK);
	        	
	        });
        
        approximateCB = new JCheckBox("Truncated likelihood computation");
        approximateCB.addActionListener(e->
            {
                boolean approxOK = approximateCB.isSelected();
	        	boolean autoOK = approxOK && approximate_autoCB.isSelected();
                approximate_thresholdT.setEnabled(approxOK && !autoOK);
                approximate_thresholdT.setEditable(approxOK && !autoOK);
                approximate_factorT.setEnabled(approxOK && !autoOK);
                approximate_factorT.setEditable(approxOK && !autoOK);
                approximate_autoCB.setEnabled(approxOK);
            });
        
        opt_ratesRB = new JRadioButton("BFGS");
        opt_ratesRB.setSelected(true);
        opt_ratesRB.setToolTipText("Slower calculations, but more precise.");
        opt_emRB = new JRadioButton("Expectation-maximization");
        opt_emRB.setToolTipText("Faster");
        opt_distributionRB = new JRadioButton("BFGS, one class");
        opt_distributionRB.setToolTipText("Slower convergence, tracking the gradient for getting to maximum-likelihood");
        ButtonGroup optG = new ButtonGroup();
        optG.add(opt_ratesRB);
        optG.add(opt_distributionRB);
        optG.add(opt_emRB);
        
        
        
        opt_ratesRB.addChangeListener(e->
        	{
        		boolean want_rates = opt_ratesRB.isSelected();
        		approximate_autoCB.setEnabled(want_rates);
        		if (!want_rates) approximate_autoCB.setSelected(false);
        	});
        
        NumberFormat threadsF = NumberFormat.getIntegerInstance();
        threadsT = new JFormattedTextField(threadsF);
        threadsT.setInputVerifier(InputVerifiers.getPositiveInstance());            

        threadsT.setColumns(6);
        threadsT.setToolTipText("Number of threads used in calculations for command-line execution");
        threadsT.setValue(Count.THREAD_PARALLELISM);
        JLabel threadsL = new JLabel("Multithreading");
        threadsL.setLabelFor(threadsT);
        Box threadsB = new Box(BoxLayout.LINE_AXIS);
        threadsB.add(threadsL);
        threadsB.add(Box.createHorizontalStrut(12));
        threadsB.add(threadsT);

        JPanel convergenceP = new JPanel();
        String panel_title;
        if (want_simulation)
        {
        	panel_title = "Simulation parameters";
        } else
        {
        	panel_title = want_optimization?"Optimization parameters":"Computation parameters";
        }
        
        
        String optimization_instructions = "<p><em>There are multiple numerical optimization methods implemented for "
        		+ "likelihood maximization. BFGS (Davidon-Fletcher-Powell)_uses the gradient to find "
        		+ "a local maximum, and may work over the rate parameters, or over distribution "
        		+ "parameters if only one rate class is used. (The latter is faster and more robust, but does not "
        		+ "generalize to multiple rate classes. Alternatively, Expectation-Maximization (EM) is possible if there is only "
        		+ "one class, which is the fastest but less robust numerically. "
        		+ "For BFGS optimization, <b>truncation</b> sets a threshold on maximum "
        		+ "ancestral copy number by multiplying the maximum copy number at the "
        		+ "descendant leaves with the <b>relative</b> value, and choosing the maximum "
        		+ "between the product and the given absolute threshold. "
        		+ "Without truncation, the likelihoods and gradients are computed exactly "
        		+ "to numerical precision. With auto-truncation, the best setting is computed automatically for the desired precision."
        		+ "The <b>convergence</b> threshold for BFGS is on "
        		+ "the relative value of the partial derivatives."
        		+ "For EM,  <b>truncation</b> is mandatory and defines the "
        		+ "ancestral copy number threshold as <var>m</var>+<var>a</var>+<var>r</var>×sqrt(1+<var>m</var>) "
        		+ "where <var>m</var> is the maximum copy number at the leaves, "
        		+ "and <var>r</var>,<var>a</var> are the specified relative and absolute "
        		+ "values. The <b>convergence</b> threshold for EM is on the relative change of the "
        		+ "log-likelihood. "
        		+ "<b>Multithreading</b> speeds up the computations; default value is "
        		+ "three-quarter of your available vCPUs."
        		+ "</em></p>"
        		+ "<p><em>Best practices: "
        		+ "(1) Start with a random model + expectation-maximization at limited truncation (e.g., relative 4.0 and absolute 4), a few (200-500) iterations should suffice; "
        		+ "(2) refine using BFGS with auto-truncation and stricter convergence (in the order of the square root of machine precision, say 1e-9) and a few thousand iterations."
        		+ "(3) If the optimization finishes with the maximum number of iterations, then the model can converge more, and you can relaunch BFGS."
        		+ "Otherwise, the model is converged to a (local) optimum. "
        		+  "You can repeat the procedure "
        		+ "with different random starting models if you suspect that there are better solutions. "
        		+ "(4) You can assess model fit by simulating a table with the given model, and comparing profile distributions against the input and the posterior reconstruction."
        		+"(5) You can validate the inference procedure by launching model optimization on a simulated table as input, in order to see how model parameters "
        		+ "are recovered, and how posterior reconstruction relates to the true history in the simulation."
        		+ "</p>";
        
        JEditorPane optimization_explain = new JEditorPane("text/html", optimization_instructions);
        optimization_explain.setEditable(false);
        
        convergenceP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),panel_title,javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
        {
            GridBagLayout layout = new GridBagLayout();
            GridBagConstraints gc = new GridBagConstraints();
            convergenceP.setLayout(layout);

            gc.fill = GridBagConstraints.HORIZONTAL;
            gc.gridy = 0;
            gc.gridx = 0;
            gc.gridheight=1;
            gc.gridwidth=1;
            gc.anchor = GridBagConstraints.WEST;
            gc.weightx = 0.1;
            gc.weighty = 0.1;
            layout.setConstraints(threadsL, gc);
            convergenceP.add(threadsL);
            gc.gridx++;
            gc.fill = GridBagConstraints.HORIZONTAL;
            layout.setConstraints(threadsT, gc);
            convergenceP.add(threadsT);
            
            threadsT.setVisible(want_optimization);
            threadsL.setVisible(threadsT.isVisible());
            
            ++gc.gridx;
            gc.fill = GridBagConstraints.HORIZONTAL;
            
            if (want_optimization)
            {
            	layout.setConstraints(rnd_seedL, gc);
            	convergenceP.add(rnd_seedL);
                ++gc.gridx;
                layout.setConstraints(rnd_seedT, gc);
                convergenceP.add(rnd_seedT);
                ++gc.gridx;
            }
            
            
            
            gc.gridy++;

            
            gc.gridx = 2;
            layout.setConstraints(roundL,gc);
            convergenceP.add(roundL);
            roundL.setVisible(want_optimization || want_simulation);

            ++ gc.gridy; //gc.gridy=1;
            if (want_simulation)
            {
	            layout.setConstraints(rnd_seedL,gc);
	            convergenceP.add(rnd_seedL);
            } else 
            {
	            layout.setConstraints(epsL,gc);
	            convergenceP.add(epsL);
	            epsL.setVisible(want_optimization || want_simulation);
            }
        
            

            gc.gridx++;
            --gc.gridy; //gc.gridy=0;
            gc.anchor = GridBagConstraints.WEST;
//            gc.fill=GridBagConstraints.NONE;
//            gc.weightx=0.0;
            layout.setConstraints(roundT,gc);
            convergenceP.add(roundT);
            roundT.setVisible(roundL.isVisible());

            ++gc.gridy; //gc.gridy=1;
            if (want_simulation)
            {
                layout.setConstraints(rnd_seedT,gc);
                convergenceP.add(rnd_seedT);
            } else
            {
                layout.setConstraints(epsT,gc);
                convergenceP.add(epsT);
                epsT.setVisible(epsL.isVisible());
            }
            

            ++gc.gridy; //gc.gridy=2;
            gc.gridx--;
            gc.gridwidth = 2;
            layout.setConstraints(minCopiesB, gc);
            convergenceP.add(minCopiesB);
            minCopiesB.setVisible(want_optimization || want_simulation);
//            layout.setConstraints(minCopiesB, gc);
//            convergenceP.add(minCopiesB);

            --gc.gridy;
            --gc.gridy;
            // gc.gridy=0;
            
            gc.gridx=0;
            //++gc.gridx;
            gc.gridwidth = 1;
            gc.gridheight = 1;
            gc.fill = GridBagConstraints.NONE;
            gc.weightx = 0.1;
            gc.weighty = 0.1;
            layout.setConstraints(approximateCB, gc);
            convergenceP.add(approximateCB);
            approximateCB.setVisible(!want_simulation);
            
            ++gc.gridx;
            layout.setConstraints(approximate_autoCB, gc);
            convergenceP.add(approximate_autoCB);
            approximate_autoCB.setVisible(want_optimization);
            --gc.gridx;
            
//            gc.gridx++;
//            layout.setConstraints(threadsB, gc);
//            convergenceP.add(threadsB);
//            --gc.gridx;
//            threadsB.setVisible(want_optimization);
            
            gc.gridwidth=1;
            gc.weightx=0.1;
            ++gc.gridy; // gc.gridy=1
            gc.fill = GridBagConstraints.HORIZONTAL;
            layout.setConstraints(approximate_factorL, gc);
            convergenceP.add(approximate_factorL);
            approximate_factorL.setVisible(!want_simulation);

            ++gc.gridx;
//            gc.weightx=0.0;
//            gc.fill=GridBagConstraints.NONE;
            layout.setConstraints(approximate_factorT,gc);
            convergenceP.add(approximate_factorT);
            approximate_factorT.setVisible(approximate_factorL.isVisible());

            --gc.gridx;
            ++gc.gridy; // gc.gridy=2
            gc.weightx=0.1;
            gc.fill = GridBagConstraints.HORIZONTAL;
            layout.setConstraints(approximate_thresholdL, gc);
            convergenceP.add(approximate_thresholdL);
            approximate_thresholdL.setVisible(approximate_factorL.isVisible());

            ++gc.gridx;
//            gc.weightx=0.0;
//            gc.fill=GridBagConstraints.NONE;
            layout.setConstraints(approximate_thresholdT, gc);
            convergenceP.add(approximate_thresholdT);
            approximate_thresholdT.setVisible(approximate_thresholdL.isVisible());

            gc.gridx = 4;
            --gc.gridy; 
            --gc.gridy; // gc.gridy = 0;

            gc.fill = GridBagConstraints.HORIZONTAL;
            gc.gridheight=1;
            gc.gridwidth=1;
            gc.anchor = GridBagConstraints.WEST;
            gc.weightx = 0.1;
            gc.weighty = 0.1;
            
            --gc.gridy;
            JLabel algoL = new JLabel("Numerical optimization algorithm");
            layout.setConstraints(algoL, gc);
            convergenceP.add(algoL);
            algoL.setVisible(want_optimization);
            ++gc.gridy;
            
            
            layout.setConstraints(opt_ratesRB, gc);
            convergenceP.add(opt_ratesRB);
            gc.gridy++; // gc.gridy=1
            layout.setConstraints(opt_distributionRB, gc);
            convergenceP.add(opt_distributionRB);
            gc.gridy++; // gc.gridy=2
            layout.setConstraints(opt_emRB,  gc);
            convergenceP.add(opt_emRB);
            
            opt_ratesRB.setVisible(want_optimization);
            opt_distributionRB.setVisible(DEVELOPMENT_CODE && opt_ratesRB.isVisible());
            //opt_distributionRB.setEnabled(DEVELOPMENT_CODE); 
            opt_emRB.setVisible(opt_ratesRB.isVisible());

            // last column 
            int last_column = 5; //want_optimization?5:4;
            
            gc.gridx = 0;
            gc.gridy += 2; // gc.gridy=4
            gc.gridwidth=last_column;
            gc.fill = GridBagConstraints.BOTH;
            layout.setConstraints(optimization_explain, gc);
            convergenceP.add(optimization_explain);
            optimization_explain.setVisible(want_optimization);
            
            
//            gc.gridx=last_column;
//            gc.gridy=0;
//            gc.gridheight=4;
//            gc.gridwidth=1;
//            gc.weightx=1.0;
//            gc.fill=GridBagConstraints.BOTH;
//            {
//                Component fillerX = Box.createHorizontalGlue();
//                layout.setConstraints(fillerX,gc);
//                convergenceP.add(fillerX);
//            }
            
            
//            Component fillerY = Box.createVerticalGlue();
//            ++gc.gridy;
//            gc.gridwidth=1;
//            gc.gridheight--;
//            gc.weighty=1.0;
//            layout.setConstraints(fillerY, gc);
//            convergenceP.add(fillerY);
            
        }
        return convergenceP;
    }
    
    /**
     * Sets the text field values in {@link #cat_probT}, {@link #cat_mod_lenT} and {@link #cat_mod_dupT}
     * 
     *
     */
    private void setRateVariationModelFields()
    {
    	// TODO update with number of categories
    	// TODO handle legacy of GammaInvariant 
    	MixedRateModel sel = models.get(selected_model);
		// System.out.println("#**MSD.sRVMF select "+selected_model+"\t"+sel.getClass());
    	if (sel instanceof GammaInvariant)
    	{
    		//  nothing to do 
    		
    	}
    	if (sel instanceof RateVariationModel)
    	{
        	RateVariationModel model = (RateVariationModel)sel;
    	
        	int model_ncat = model.getNumClasses();
	    	int ncat = cat_probT.size();
	    	
	    	int wanted_ncat = ((Number) num_catT.getValue()).intValue();
	    	if (wanted_ncat != ncat)
	    	{
	    		JComponent varB = this.createRateVariationModelParameterBox(wanted_ncat);
	        	main_tabbed_pane.setComponentAt(main_tabbed_pane.getTabCount()-1, varB);	    		
	        	ncat = wanted_ncat;
	    	}
	
	    	// System.out.println("#**MSD.sRVMF ncat "+ncat+"\tmodel "+model_ncat);
	    	
	    	double sum_p = 0.0;
	    	Random RND = null;
	    	
	
	    	double scale_model;
	    	if (ncat<model_ncat)
	    	{
	    		// model_ncat-ncat dropped categories
	    		for (int k=0; k<ncat; k++)
	    		{
	    			double p = model.getClassProbability(k);
	    			sum_p += p;
	    		}
	    		scale_model = 1.0/sum_p;
	    	} else 
	    	{
	        	//    ncat-model_ncat new categories
	    		scale_model = (model_ncat+0.0)/ncat;
	    	}
	    	
	    	for (int k=0; k<ncat; k++)
	    	{
	    		double modlen;
	    		double moddup;
	    		double p;
	    		if (k<model_ncat)
	    		{
	    			RateVariationModel.Category C = model.getCategory(k);
	    			modlen = C.getModLength();
	    			moddup = C.getModDuplication();
	    			p = scale_model * Math.exp(C.getLogCatProbability());
	    		} else
	    		{
	    			if (k==0)
	    			{
	    				modlen = 0.0;
	    				moddup = 0.0;
	    			} else 
	    			{
		    			if (RND==null) RND =new Random();
		    			modlen = Math.log(-Math.log(RND.nextDouble())); // log of exponential(1)
		    			moddup = Math.log(-Math.log(RND.nextDouble()));
	    			}
	    			p = (1.0-sum_p)/(ncat-k);
	    		}
	    		cat_probT.get(k).setValue(p);
	    		cat_mod_lenT.get(k).setValue(modlen);
	    		cat_mod_dupT.get(k).setValue(moddup);
	    	}
    	}
    }
    
    
    private JComponent createRateVariationModelParameterBox(int ncat)
    {
    	cat_mod_lenT.clear();
    	cat_mod_dupT.clear();
    	cat_probT.clear();
    	

    	NumberFormat modF = new DecimalFormat("0.############");
		modF.setParseIntegerOnly(false);    		
		InputVerifier probV = new InputVerifiers.RangeVerifier(0.0,1.0);
		
		//Box catB = new Box(BoxLayout.Y_AXIS);
		Box catB = new Box(BoxLayout.Y_AXIS);
		Box probB = new Box(BoxLayout.Y_AXIS);
		Box modlenB = new Box(BoxLayout.Y_AXIS);
		Box moddupB = new Box(BoxLayout.Y_AXIS);
		
		JLabel hdrL = new JLabel("Category");
		JLabel probL = new JLabel("Probability");
		JLabel modlenL = new JLabel("Length modifier");
		JLabel moddupL = new JLabel("Duplication modifier");
		
//		catB.add(hdrL);
		probB.add(probL);
		modlenB.add(modlenL);
		moddupB.add(moddupL);
		
    	for (int k=0; k<ncat; k++)
    	{

    		ParameterField modlenT = new ParameterField(modF);
    		ParameterField moddupT = new ParameterField(modF);
    		ParameterField probT = new ParameterField(modF);
    		probT.setInputVerifier(probV);
    		
    		modlenT.setColumns(16);
    		moddupT.setColumns(16);
    		probT.setColumns(12);
    		
    		modlenT.parameter_text = "Length modifier in category "+k;
    		moddupT.parameter_text = "Duplication modifier in category "+k;
    		probT.parameter_text = "Probability of category "+k;    		
    		probT.setEditable(false);
    		
    		
    		modlenT.fixedCB.setEnabled(false);
    		moddupT.fixedCB.setEnabled(false);
    		probT.fixedCB.setEnabled(false);
    		
    		if (k==0)
    		{
    			modlenT.freezeParameter();
    			moddupT.freezeParameter();
    			modlenT.disabled_reason =  moddupT.disabled_reason = "Base rate model";
    		}
    		
    		cat_mod_lenT.add(modlenT);
    		cat_mod_dupT.add(moddupT);
    		cat_probT.add(probT);
    		
    		JLabel catL = new JLabel(Integer.toString(k));

//    		catB.add(catL);
    		probB.add(probT);
    		modlenB.add(modlenT);
    		moddupB.add(moddupT);
    	}
    	
    	probB.add(Box.createVerticalGlue());
    	modlenB.add(Box.createVerticalGlue());
    	moddupB.add(Box.createVerticalGlue());
    	catB.add(Box.createVerticalGlue());
    	
    	Box pbox = new Box(BoxLayout.X_AXIS);
    	pbox.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Logistic shift categories for variation across families ",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
    	
    	pbox.add(catB);
    	pbox.add(probB);
    	pbox.add(modlenB);
    	pbox.add(moddupB);
    	pbox.add(probB);
    	pbox.add(modlenB);
    	pbox.add(moddupB);
    	
    	return pbox;
    	
    }
    
    /**
     * Panel for rate variation (across families) parameters.
     * Initialzes {@link #alpha_lengthT}, {@link #alpha_gainT}, {@link #alpha_duplicationT},
     * {@link #alpha_lossT}, {@link #forbidden_gainT}, {@link #forbidden_duplicationT}.
     *
     * Initially field values are unset, and check boxes are unselected.
     *
     * @return a container
     */
    private JComponent createGammaFamilyVariationParameterBox()
    {
        NumberFormat alphaF = new DecimalFormat("0.############");
        alphaF.setParseIntegerOnly(false);
        alpha_lengthT = new ParameterField(alphaF,1.0,"Fixed");
        alpha_lengthT.setInputVerifier(InputVerifiers.getPositiveInstance());
        alpha_lengthT.setColumns(16);
        alpha_lengthT.parameter_text="Gamma distribution shape parameter for length variation across families";

        alpha_gainT = new ParameterField(alphaF, 1.0, "Fixed");
        alpha_gainT.setInputVerifier(InputVerifiers.getPositiveInstance());
        alpha_gainT.setColumns(16);
        alpha_gainT.parameter_text="Gamma distribution shape parameter for gain rate variation across families";
        alpha_duplicationT = new ParameterField(alphaF, 1.0, "Fixed");
        alpha_duplicationT.setInputVerifier(InputVerifiers.getPositiveInstance());
        alpha_duplicationT.setColumns(16);
        alpha_duplicationT.parameter_text="Gamma distribution shape parameter for duplication rate variation across families";
        alpha_lossT = new ParameterField(alphaF, 1.0, "Fixed");
        alpha_lossT.setInputVerifier(InputVerifiers.getPositiveInstance());
        alpha_lossT.setColumns(16);
        alpha_lossT.parameter_text="Gamma distribution shape parameter for loss rate variation across families";

        JLabel variation_lengthL = new JLabel("Edge length");
        JLabel variation_gainL = new JLabel("Gain rate");
        JLabel variation_duplicationL = new JLabel("Duplication rate");
        JLabel variation_lossL = new JLabel("Loss rate");

        JLabel alpha_lengthL = new JLabel("Gamma shape parameter (\u03b1)");
        JLabel alpha_gainL = new JLabel("Gamma shape parameter (\u03b1)");
        JLabel alpha_duplicationL = new JLabel("Gamma shape parameter (\u03b1)");
        JLabel alpha_lossL = new JLabel("Gamma shape parameter (\u03b1)");

        JLabel forbidden_duplicationL = new JLabel("Proportion of families with no duplications");
        JLabel forbidden_gainL = new JLabel("Proportion of families with no gains");

        NumberFormat forbiddenF = new DecimalFormat("0.############");
        forbiddenF.setParseIntegerOnly(false);
        forbidden_gainT = new ParameterField(forbiddenF,0.0,"Fixed");
        forbidden_gainT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
        forbidden_gainT.setColumns(16);
        forbidden_gainT.parameter_text="Prior probability for forbidden-gain category";

        forbidden_duplicationT = new ParameterField(forbiddenF, 0.0, "Fixed");
        forbidden_duplicationT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
        forbidden_duplicationT.setColumns(16);
        forbidden_duplicationT.parameter_text="Prior probability for forbidden-duplication category";

        forbidden_gainT.freezeParameter();
        forbidden_duplicationT.freezeParameter();
        
        
        alpha_gainT.setVisible(WANT_GENERAL_VARIATION);
        variation_gainL.setVisible(alpha_gainT.isVisible());
        alpha_gainL.setVisible(alpha_gainT.isVisible());
        alpha_duplicationT.setVisible(WANT_GENERAL_VARIATION);
        alpha_duplicationL.setVisible(alpha_gainT.isVisible());
        variation_duplicationL.setVisible(WANT_GENERAL_VARIATION);
        
        forbidden_gainT.setVisible(WANT_GENERAL_VARIATION);
        forbidden_gainL.setVisible(forbidden_gainT.isVisible());
        forbidden_duplicationT.setVisible(WANT_GENERAL_VARIATION);        
        forbidden_duplicationL.setVisible(forbidden_duplicationT.isVisible());
        
        JPanel rate_variationP = new JPanel();
        rate_variationP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Rate variation across families",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
        {
            GridBagLayout layout = new GridBagLayout();
            GridBagConstraints gc = new GridBagConstraints();
            rate_variationP.setLayout(layout);

            gc.fill = GridBagConstraints.HORIZONTAL;
            gc.gridy = 0;
            gc.gridx = 0;
            gc.gridheight=1;
            gc.gridwidth=1;
            gc.anchor = GridBagConstraints.WEST;
            gc.weightx = 0.1;
            gc.weighty = 0.1;

            layout.setConstraints(variation_lengthL,gc);
            rate_variationP.add(variation_lengthL);
            gc.gridx++;
            layout.setConstraints(alpha_lengthL, gc);
            rate_variationP.add(alpha_lengthL);
            gc.gridx++;
            gc.fill=GridBagConstraints.NONE;
            layout.setConstraints(alpha_lengthT, gc);
            rate_variationP.add(alpha_lengthT);
            gc.gridx++;
            layout.setConstraints(alpha_lengthT.fixedCB, gc);
            rate_variationP.add(alpha_lengthT.fixedCB);
            gc.gridx++;

            gc.gridx=0;
            gc.gridy++;
            gc.fill = GridBagConstraints.HORIZONTAL;
            layout.setConstraints(variation_lossL,gc);
            rate_variationP.add(variation_lossL);
            gc.gridx++;
            layout.setConstraints(alpha_lossL, gc);
            rate_variationP.add(alpha_lossL);
            gc.gridx++;
            gc.fill=GridBagConstraints.NONE;
            layout.setConstraints(alpha_lossT, gc);
            rate_variationP.add(alpha_lossT);
            gc.gridx++;
            layout.setConstraints(alpha_lossT.fixedCB, gc);
            rate_variationP.add(alpha_lossT.fixedCB);
            gc.gridx++;

            variation_lossL.setVisible(false);
            alpha_lossL.setVisible(false);
            alpha_lossT.setVisible(false);

            gc.gridx=0;
            gc.gridy++;
            gc.fill = GridBagConstraints.HORIZONTAL;
            layout.setConstraints(variation_gainL,gc);
            rate_variationP.add(variation_gainL);
            gc.gridx++;
            layout.setConstraints(alpha_gainL, gc);
            rate_variationP.add(alpha_gainL);
            gc.gridx++;
            gc.fill=GridBagConstraints.NONE;
            layout.setConstraints(alpha_gainT, gc);
            rate_variationP.add(alpha_gainT);
            gc.gridx++;
            layout.setConstraints(alpha_gainT.fixedCB, gc);
            rate_variationP.add(alpha_gainT.fixedCB);
            gc.gridx++;
            gc.fill = GridBagConstraints.HORIZONTAL;
            layout.setConstraints(forbidden_gainL, gc);
            rate_variationP.add(forbidden_gainL);
            gc.gridx++;
            gc.fill = GridBagConstraints.NONE;
            layout.setConstraints(forbidden_gainT,gc);
            rate_variationP.add(forbidden_gainT);
            gc.gridx++;
            layout.setConstraints(forbidden_gainT.fixedCB, gc);
            rate_variationP.add(forbidden_gainT.fixedCB);

            gc.gridx=0;
            gc.gridy++;
            gc.fill = GridBagConstraints.HORIZONTAL;
            layout.setConstraints(variation_duplicationL,gc);
            rate_variationP.add(variation_duplicationL);
            gc.gridx++;
            layout.setConstraints(alpha_duplicationL, gc);
            rate_variationP.add(alpha_duplicationL);
            gc.gridx++;
            gc.fill=GridBagConstraints.NONE;
            layout.setConstraints(alpha_duplicationT, gc);
            rate_variationP.add(alpha_duplicationT);
            gc.gridx++;
            layout.setConstraints(alpha_duplicationT.fixedCB, gc);
            rate_variationP.add(alpha_duplicationT.fixedCB);
            gc.gridx++;
            gc.fill = GridBagConstraints.HORIZONTAL;
            layout.setConstraints(forbidden_duplicationL, gc);
            rate_variationP.add(forbidden_duplicationL);
            gc.gridx++;
            gc.fill = GridBagConstraints.NONE;
            layout.setConstraints(forbidden_duplicationT,gc);
            rate_variationP.add(forbidden_duplicationT);
            gc.gridx++;
            layout.setConstraints(forbidden_duplicationT.fixedCB, gc);
            rate_variationP.add(forbidden_duplicationT.fixedCB);

            gc.gridheight=gc.gridy+1;
            gc.gridx++;
            gc.gridy=0;
            gc.weightx=1.0;
            gc.fill=GridBagConstraints.BOTH;
            Component fillerX = Box.createHorizontalGlue();
            layout.setConstraints(fillerX,gc);
            rate_variationP.add(fillerX);
        }


        return rate_variationP;
    }
    
    /**
     * Panel for root prior parameters
     * Initializes {@link #root_aT}, {@link #root_bT}, {@link #fixed_rootCB},
     * {@link #root_distributionL},
     * {@link #root_aL} and {@link #root_bL}.
     *
     * Label texts are expected to be set when root prior is set.
     *
     * @return a container
     */
    private JComponent createRootPriorParameterBox()
    {
        NumberFormat root_aF = new DecimalFormat("0.############");
        root_aF.setParseIntegerOnly(false);
        root_aT = new JFormattedTextField(root_aF);
        root_aT.setInputVerifier(new InputVerifiers.NonnegativeInputVerifier());
        root_aT.setColumns(16);
        root_aT.setEditable(false);

        NumberFormat root_bF = new DecimalFormat("0.############");
        root_bF.setParseIntegerOnly(false);
        root_bT = new JFormattedTextField(root_bF);
        root_bT.setInputVerifier(new InputVerifiers.NonnegativeInputVerifier());
        root_bT.setColumns(16);
        root_bT.setEditable(false);

        //NumberFormat root_meanF = new DecimalFormat("0.############");
        //root_meanF.setParseIntegerOnly(false);
        root_meanT = new JFormattedTextField(new Object());
        root_meanT.setColumns(16);
        root_meanT.setEditable(false);
        root_meanT.setToolTipText("Mean of the selected distribution");
        root_meanT.setValue(null);

//        root_meanT.addPropertyChangeListener("value", new PropertyChangeListener() {
//
//            @Override
//            public void propertyChange(PropertyChangeEvent pce)
//            {
//                System.out.println("#*MSP.MD.cRPPB,mean update "+root_meanT.getValue()+"\t// "+pce);
//            }
//        });

        root_distributionL = new JLabel("");
        root_aL = new JLabel("par0");
        root_bL = new JLabel("par1");
        root_meanL = new JLabel("Mean");
        
        fixed_rootCB = new JCheckBox("Fixed");

        JPanel root_distributionP = new JPanel();
        root_distributionP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Family size distribution at root",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
        {
            GridBagLayout layout = new GridBagLayout();
            GridBagConstraints gc = new GridBagConstraints();
            root_distributionP.setLayout(layout);

            gc.fill = GridBagConstraints.NONE;
            gc.gridy = 0;
            gc.gridx = 0;
            gc.gridheight=1;
            gc.gridwidth=1;
            gc.anchor = GridBagConstraints.WEST;
            gc.weightx = 0.1;
            gc.weighty = 0.1;

            layout.setConstraints(root_distributionL,gc);
            root_distributionP.add(root_distributionL);
            gc.gridx++;
            layout.setConstraints(root_aL, gc);
            root_distributionP.add(root_aL);
            gc.gridx++;
            layout.setConstraints(root_aT, gc);
            root_distributionP.add(root_aT);
            gc.gridx++;
            layout.setConstraints(root_bL, gc);
            root_distributionP.add(root_bL);
            gc.gridx++;
            layout.setConstraints(root_bT, gc);
            root_distributionP.add(root_bT);

            gc.gridx++;
            layout.setConstraints(fixed_rootCB, gc);
            root_distributionP.add(fixed_rootCB);

            gc.gridx++;
            layout.setConstraints(root_meanL, gc);
            root_distributionP.add(root_meanL);

            gc.gridx++;
            layout.setConstraints(root_meanT, gc);
            root_distributionP.add(root_meanT);


            gc.gridheight=gc.gridy+1;
            gc.gridx++;
            gc.gridy=0;
            gc.weightx=1.0;
            gc.fill=GridBagConstraints.BOTH;
            Component fillerX = Box.createHorizontalGlue();
            layout.setConstraints(fillerX,gc);
            root_distributionP.add(fillerX);
        }

        return root_distributionP;
    }
    
    
    
    private JComponent createLineageSpecificParameterBox()
    {
    	JPanel parametersP = new JPanel();
    	
    	TreeWithRates null_model = models.get(0).getBaseModel();
    	IndexedTree tree = null_model.getTree(); // same tree as models.get(1), if alt.model exists
    	//int num_nodes = tree.getNumNodes();
    	//assert (tree.isRoot(num_nodes-1)); // root rates are not shown
    	
        lineage_lengthP = new LineageSpecificParameters(tree, InputVerifiers.getPositiveInstance());
        lineage_lengthP.parameter_text="Length";
        lineage_duplicationP = new LineageSpecificParameters(tree, InputVerifiers.getNonnegativeInstance()); // new InputVerifiers.RangeVerifier(ML.MIN_DUPLICATION_RATE, ML.MAX_DUPLICATION_RATE));
        lineage_duplicationP.parameter_text="Duplication rate";
        lineage_gainP = new LineageSpecificParameters(tree, InputVerifiers.getNonnegativeInstance()); // new InputVerifiers.RangeVerifier(ML.MIN_TRANSFER_RATE, ML.MAX_TRANSFER_RATE));
        lineage_gainP.parameter_text="Gain rate";
        lineage_lossP = new LineageSpecificParameters(tree, InputVerifiers.getPositiveInstance());
        lineage_lossP.parameter_text="Loss rate";
        
        // initial state corresponds to homogeneous loss check box unselected, gld model
        lineage_lossP.disabled_reason="rates and lengths are scaled by loss rate.";
        lineage_lossP.freezeAll();
        
        // row headers
        JLabel[] node_nameL = new JLabel[lineage_lengthP.numFields()];
        for (int node=0; node<node_nameL.length; ++node)
            node_nameL[node] = new JLabel(tree.getIdent(node));
        JLabel all_fixedL = new JLabel("all edges");
        all_fixedL.setFont(all_fixedL.getFont().deriveFont(Font.ITALIC));

        // column headers
        JLabel nodeL = new JLabel("Edge to");
        nodeL.setFont(nodeL.getFont().deriveFont(Font.BOLD));

        JLabel lengthL = new JLabel("Length");
        lengthL.setFont(lengthL.getFont().deriveFont(Font.BOLD));
        JLabel fixed_lengthL = new JLabel("Fixed");
        fixed_lengthL.setFont(fixed_lengthL.getFont().deriveFont(Font.BOLD));

        JLabel gainL = new JLabel("Gain rate");
        gainL.setFont(gainL.getFont().deriveFont(Font.BOLD));
        JLabel fixed_gainL = new JLabel("Fixed");
        fixed_gainL.setFont(fixed_gainL.getFont().deriveFont(Font.BOLD));

        JLabel duplicationL = new JLabel("Duplication rate");
        duplicationL.setFont(duplicationL.getFont().deriveFont(Font.BOLD));
        JLabel fixed_duplicationL = new JLabel("Fixed");
        fixed_duplicationL.setFont(fixed_duplicationL.getFont().deriveFont(Font.BOLD));

        JLabel lossL = new JLabel("Loss rate");
        lossL.setFont(lossL.getFont().deriveFont(Font.BOLD));
        JLabel fixed_lossL = new JLabel("Fixed");
        fixed_lossL.setFont(fixed_lossL.getFont().deriveFont(Font.BOLD));

        parametersP.setBorder(BorderFactory.createTitledBorder(
        		BorderFactory.createEtchedBorder(),
        		"Lineage-specific rates",
        		javax.swing.border.TitledBorder.LEFT,
        		javax.swing.border.TitledBorder.TOP,
        		getFont().deriveFont(Font.BOLD)));
        
        {
            GridBagLayout layout = new GridBagLayout();
            GridBagConstraints gc = new GridBagConstraints();
            parametersP.setLayout(layout);

            gc.fill = GridBagConstraints.HORIZONTAL;
            gc.anchor = GridBagConstraints.CENTER;
            gc.gridy = 0;
            gc.gridx = 0;
            gc.gridwidth=1;
            gc.gridheight=1;
            gc.weightx = 0.1;
            gc.weighty = 0.1;

            layout.setConstraints(nodeL, gc);
            parametersP.add(nodeL);
            gc.gridx++;
            layout.setConstraints(lengthL, gc);
            parametersP.add(lengthL);
            gc.gridx++;
            layout.setConstraints(fixed_lengthL, gc);
            parametersP.add(fixed_lengthL);
            gc.gridx++;
            layout.setConstraints(gainL, gc);
            parametersP.add(gainL);
            gc.gridx++;
            layout.setConstraints(fixed_gainL, gc);
            parametersP.add(fixed_gainL);
            gc.gridx++;
            layout.setConstraints(duplicationL, gc);
            parametersP.add(duplicationL);
            gc.gridx++;
            layout.setConstraints(fixed_duplicationL, gc);
            parametersP.add(fixed_duplicationL);
            gc.gridx++;
            layout.setConstraints(lossL, gc);
            parametersP.add(lossL);
            gc.gridx++;
            layout.setConstraints(fixed_lossL, gc);
            parametersP.add(fixed_lossL);

            gc.gridx=0;
            gc.gridy++;
            gc.anchor = GridBagConstraints.WEST;
            gc.fill = GridBagConstraints.HORIZONTAL;

            layout.setConstraints(all_fixedL, gc);
            parametersP.add(all_fixedL);
            gc.gridx++;

            gc.gridx++;
//            layout.setConstraints(all_fixed_lengthCB, gc);
//            parametersP.add(all_fixed_lengthCB);
            layout.setConstraints(lineage_lengthP.fixed_parameters, gc);
            parametersP.add(lineage_lengthP.fixed_parameters);
            gc.gridx++;

            gc.gridx++;
//            layout.setConstraints(all_fixed_gainCB, gc);
//            parametersP.add(all_fixed_gainCB);
            layout.setConstraints(lineage_gainP.fixed_parameters, gc);
            parametersP.add(lineage_gainP.fixed_parameters);
            gc.gridx++;

            gc.gridx++;
//            layout.setConstraints(all_fixed_duplicationCB, gc);
//            parametersP.add(all_fixed_duplicationCB);
            layout.setConstraints(lineage_duplicationP.fixed_parameters, gc);
            parametersP.add(lineage_duplicationP.fixed_parameters);
            gc.gridx++;

            gc.gridx++;
//            layout.setConstraints(all_fixed_lossCB, gc);
//            parametersP.add(all_fixed_lossCB);
            layout.setConstraints(lineage_lossP.fixed_parameters, gc);
            parametersP.add(lineage_lossP.fixed_parameters);
            gc.gridx++;

            for (int node=0; node<node_nameL.length; node++)
            {
                gc.gridx=0;
                gc.gridy++;

                gc.anchor = GridBagConstraints.WEST;
                gc.fill = GridBagConstraints.HORIZONTAL;
                layout.setConstraints(node_nameL[node], gc);
                parametersP.add(node_nameL[node]);
                gc.gridx++;
                gc.fill = GridBagConstraints.NONE;
                layout.setConstraints(lineage_lengthP.getField(node), gc);
                parametersP.add(lineage_lengthP.getField(node));
                gc.gridx++;
                gc.fill = GridBagConstraints.HORIZONTAL;
                layout.setConstraints(lineage_lengthP.getCheckBox(node), gc);
                parametersP.add(lineage_lengthP.getCheckBox(node));
                gc.gridx++;
                gc.fill = GridBagConstraints.NONE;
                layout.setConstraints(lineage_gainP.getField(node), gc);
                parametersP.add(lineage_gainP.getField(node));
                gc.gridx++;
                gc.fill = GridBagConstraints.HORIZONTAL;
                layout.setConstraints(lineage_gainP.getCheckBox(node), gc);
                parametersP.add(lineage_gainP.getCheckBox(node));
                gc.gridx++;
                gc.fill = GridBagConstraints.NONE;
                layout.setConstraints(lineage_duplicationP.getField(node), gc);
                parametersP.add(lineage_duplicationP.getField(node));
                gc.gridx++;
                gc.fill = GridBagConstraints.HORIZONTAL;
                layout.setConstraints(lineage_duplicationP.getCheckBox(node), gc);
                parametersP.add(lineage_duplicationP.getCheckBox(node));
                gc.gridx++;
                gc.fill = GridBagConstraints.NONE;
                layout.setConstraints(lineage_lossP.getField(node), gc);
                parametersP.add(lineage_lossP.getField(node));
                gc.gridx++;
                gc.fill = GridBagConstraints.HORIZONTAL;
                layout.setConstraints(lineage_lossP.getCheckBox(node), gc);
                parametersP.add(lineage_lossP.getCheckBox(node));
            } // for node_idx
            gc.gridheight=gc.gridy+1;
            gc.gridx++;
            gc.gridy=0;
            gc.weightx=1.0;
            gc.fill=GridBagConstraints.BOTH;
            Component fillerX = Box.createHorizontalGlue();
            layout.setConstraints(fillerX,gc);
            parametersP.add(fillerX);
        }
    	
    	return parametersP;
    }
    
        
    /**
     * Called when the selected initial model changes in the GUI
     * ({@link #modelRB} radio buttons).
     * Performs model type and root prior selection by
     * programmatic clicks, and
     * sets other GUI component states accordingly
     * to the newly selected model.
     *
     * @param idx the newly selected model
     */
    private void setSelectedModel(int idx)
    {
        selected_model = idx;
        
        // programmatic clicking of model structure selection
        MixedRateModel model  = models.get(selected_model);
        TreeWithRates base_model = model.getBaseModel();
        if (base_model.hasDuplication())
        {
            if (base_model.hasGain())
                model_gldRB.doClick();
            else
                model_dlRB.doClick();
        } else
        {
            if (base_model.hasGain())
                model_glRB.doClick();
            else
                model_plRB.doClick();
        } 
        
        // programmatic click for root distribution
        int root = base_model.getTree().getRoot();
        DiscreteDistribution root_distr = base_model.getRootDistribution();
        Class<? extends DiscreteDistribution> root_distr_class = root_distr.getClass();
        root_distributionB.get(root_distr_class).doClick();

        // fill root prior params display
        double[] params = root_distr.getParameters();
        root_aT.setValue(params[0]);
        if (params.length>1)
            root_bT.setValue(params[1]);

        
        if (model instanceof GammaInvariant)
        {
        	GammaInvariant gamma_model = (GammaInvariant)model;
	        
	        // fill rate variation parameters
	        num_catT.setValue(gamma_model.getNumLengthGammaCategories());
	        gamma_duplicationT.setValue(WANT_GENERAL_VARIATION?gamma_model.getNumDuplicationGammaCategories():1);
	        gamma_gainT.setValue(WANT_GENERAL_VARIATION?gamma_model.getNumGainGammaCategories():1);
	        gamma_lossT.setValue(WANT_GENERAL_VARIATION?gamma_model.getNumLossGammaCategories():1);
	
	        if (gamma_model.getDuplicationForbidden()>0.0 && WANT_GENERAL_VARIATION)
	        	forbidden_duplicationCB.setSelected(true);
	        if (gamma_model.getGainForbidden()>0.0 && WANT_GENERAL_VARIATION)
	            forbidden_gainCB.setSelected(true);
	
	        alpha_lengthT.setValue(gamma_model.getLengthAlpha());
	        alpha_duplicationT.setValue(gamma_model.getDuplicationAlpha());
	        alpha_gainT.setValue(gamma_model.getGainAlpha());
	        alpha_lossT.setValue(gamma_model.getLossAlpha());
	
	        forbidden_duplicationT.setValue(WANT_GENERAL_VARIATION?gamma_model.getDuplicationForbidden():0.0);
	        forbidden_gainT.setValue(WANT_GENERAL_VARIATION?gamma_model.getGainForbidden():0.0);
	        main_tabbed_pane.setComponentAt(main_tabbed_pane.getTabCount()-1, gamma_variationB);
	        num_catT.setValue(model.getNumClasses());
        }
        else
        {
        	// RateVariationModel
        	assert (model instanceof RateVariationModel);

        	JComponent varB = this.createRateVariationModelParameterBox(model.getNumClasses());
        	main_tabbed_pane.setComponentAt(main_tabbed_pane.getTabCount()-1, varB);
        	
        	num_catT.setValue(model.getNumClasses());
            //System.out.println("#*MSD.sSM "+idx+"\tncat "+model.getNumClasses());
        }
    	this.setRateVariationModelFields(); // not enough by property change listener 
        
        // fill lineage specific parameters
        for (int node=0; node<lineage_lengthP.numFields();node++ )
        {
        	// TODO how they get to be negative?
            double len = base_model.getEdgeLength(node);
            lineage_lengthP.getField(node).setValue(len);
            double gr = base_model.getGainRate(node);
            lineage_gainP.getField(node).setValue(gr);
            double ls = base_model.getLossRate(node);
            lineage_lossP.getField(node).setValue(ls);
            double dr = base_model.getDuplicationRate(node);
            lineage_duplicationP.getField(node).setValue(dr);
        }
        
        // lineage-specific variation
        uniform_duplicationCB.setSelected(false);
        uniform_gainCB.setSelected(false);
        uniform_lossCB.setSelected(false);
        uniform_lossCB.setVisible(false);
        
        synchronizeOptimizationOptions();
    }
  
    
    private void updateRandomInitialModel()
    {
    	if (selected_model==0)
    	{
    		MixedRateModel model  = models.get(selected_model);
            TreeWithRates rates = model.getBaseModel();
            long init_seed = ((Number)rnd_seedT.getValue()).longValue();
            if (init_seed != this.init_rnd_seed)
            {
            	if (init_seed==0L)
            		rates.setRandom(null);
            	else
            		rates.setRandom(new Random(init_seed));
            	this.init_rnd_seed = init_seed;
            	int num_nodes =  rates.getTree().getNumNodes();
            	for (int u=0; u<num_nodes; u++)
            		rates.initNodeParameters(u);
            	
            	this.setSelectedModel(selected_model);
            }
    	}
    }
    
    
    
    private void synchronizeOptimizationOptions()
    {
        if (hasRateVariation() )
        {
        	opt_ratesRB.setSelected(true);
        	opt_emRB.setEnabled(false);
        	opt_distributionRB.setEnabled(false);
        } else
        {
//        	if (opt_ratesRB.isSelected())
//        	{
//        		opt_distributionRB.setSelected(true);
//        	}
        	opt_emRB.setEnabled(true);
        	opt_distributionRB.setEnabled(true);
        }
        
    	
    }
    
    /**
     * Sets the states for fields and checkboxes related to gain
     * and duplication rate in model. Field values may be reset
     * if a rate feature is removed from the model.
     *
     * Used by listeners in {@link #initListeners() }.
     */
    private void synchronizeModelStructureStates()
    {
        boolean has_gain = model_gldRB.isSelected() || model_glRB.isSelected();
        boolean has_duplication = model_gldRB.isSelected() || model_dlRB.isSelected();

        TreeWithRates base_model = models.get(selected_model).getBaseModel();
        
        if (has_gain)
        {
            gamma_gainT.setEditable(true);
            Number g = (Number)gamma_gainT.getValue();
            if (g!=null && g.intValue()!=1)
                alpha_gainT.enableParameter();
            forbidden_gainCB.setEnabled(true);
            lineage_gainP.enableAll();
            uniform_gainCB.setEnabled(true);
            if (base_model.hasGain())
            {
            	for (int node=0; node<lineage_gainP.numFields(); node++)
            		lineage_gainP.getField(node).setValue(base_model.getGainRate(node));
            } else
            {
            	lineage_gainP.setValueAll(TreeWithRates.DEFAULT_GAIN_RATE);
            }
        } else
        {
            String reason = "model does not allow gains";
            alpha_gainT.freezeParameter();
            gamma_gainT.setEditable(false);
            gamma_gainT.setValue(1);
            forbidden_gainCB.setEnabled(false);
            forbidden_gainCB.setSelected(false);
            lineage_gainP.setDisabledReason(reason);
            lineage_gainP.freezeAll();
            lineage_gainP.resetAll(); // to 0
            uniform_gainCB.setEnabled(false);
            uniform_gainCB.setSelected(true);
        }

        if (has_duplication)
        {
            gamma_duplicationT.setEditable(true);
            Number g = (Number)gamma_duplicationT.getValue();
            if (g!=null && g.intValue()!=1)
                alpha_duplicationT.enableParameter();
            forbidden_duplicationCB.setEnabled(true);
            lineage_duplicationP.enableAll();
            uniform_duplicationCB.setEnabled(true);
            if (base_model.hasDuplication())
            {
            	for (int node=0; node<lineage_duplicationP.numFields(); node++)
            		lineage_duplicationP.getField(node).setValue(base_model.getDuplicationRate(node));
            } else
            {
            	lineage_duplicationP.setValueAll(TreeWithRates.DEFAULT_DUPLICATION_RATE);
            }
        } else
        {
            alpha_duplicationT.freezeParameter();
            gamma_duplicationT.setEditable(false);
            gamma_duplicationT.setValue(1);
            forbidden_duplicationCB.setEnabled(false);
            forbidden_duplicationCB.setSelected(false);
            lineage_duplicationP.setDisabledReason("model does not allow duplications");
            lineage_duplicationP.freezeAll();
            lineage_duplicationP.resetAll(); // to 0
            uniform_duplicationCB.setEnabled(false);
            uniform_duplicationCB.setSelected(true);
        }
        lineage_duplicationP.disableVariationIfSelected(uniform_duplicationCB);
        lineage_gainP.disableVariationIfSelected(uniform_gainCB);
        
        boolean homogeneous_loss = uniform_lossCB.isSelected(); // not implemented
        if (homogeneous_loss)
        {
            lineage_lossP.enableAll();
            lineage_lengthP.setDisabledReason("model is homogeneous");
            lineage_lengthP.freezeAll();
        } else
        {
            lineage_lossP.freezeAll();
            lineage_lossP.setDisabledReason("rates and lengths are scaled by loss rate in a heterogeneous model");
            lineage_lengthP.enableAll();
        }
        lineage_lossP.disableVariationIfSelected(uniform_lossCB);
        
        // selecting root prior for the selected model
        Class<? extends DiscreteDistribution> model_prior = null;
        if (model_gldRB.isSelected())
        	model_prior = NegativeBinomial.class;
        else if (model_glRB.isSelected())
        	model_prior = Poisson.class;
        else if (model_dlRB.isSelected())
        	model_prior = ShiftedGeometric.class;
        else if (model_plRB.isSelected())
        	model_prior = PointDistribution.class;
        assert (model_prior != null); // since no other implementations

        root_distributionB.get(model_prior).setSelected(true);  // instead of doClick() which would trigger synchronizeRootPriorStates()
        
        
        
        
		synchronizeRootPriorStates();
		
//		synchronizeOptimizationOptions();
		
    }

    /**
     * Updates the displayed value of the distribution mean ({@link #root_meanT})
     * based on the values in the GUI.
     */
    private void synchronizeRootPriorMean()
    {
        Class<? extends DiscreteDistribution> distribution_class = getSelectedRootPriorClass();
        Number a = (Number)root_aT.getValue();
        if (a==null)
        {
            root_meanT.setValue(null);
            return;
        }
        if (distribution_class==Poisson.class)
        {
            // expected value is lambda
            root_meanT.setValue(a);
        } else if (distribution_class == NegativeBinomial.class)
        {
            // expected value = q t/(1-q) ; approaches Poisson with same mean
            Number b = (Number) root_bT.getValue();
            if (b==null)
                root_meanT.setValue(null);
            else
            {
                double theta = a.doubleValue();
                double q = b.doubleValue();

                //System.out.println("#*MSP.MD.sRPM NegBin t="+theta+";\tq="+q);

                if ((1.-q) == 0.0)
                {
                    root_meanT.setValue("\u221e");
                } else
                {
                    double mean = theta * q/(1.-q);
                    root_meanT.setValue(mean);
                }
            }
        } else if (distribution_class == PointDistribution.class)
        {
            double p0 = a.doubleValue();
            root_meanT.setValue(1.-p0);
        } else if (distribution_class == ShiftedGeometric.class)
        {
            // mean is (1-p0)/(1-q)
            Number b = (Number) root_bT.getValue();
            if (b==null)
                root_meanT.setValue(null);
            else
            {
                double p0 = a.doubleValue();
                double q = b.doubleValue();
                if ((1.-q)==0.0)
                {
                    if ((1.-p0)==0.0)
                    {
                        root_meanT.setValue(0.0);
                    } else
                    {
                        root_meanT.setValue("\u221e");
                    }
                } else
                {
                    double mean = (1.-p0)/(1.-q);
                    root_meanT.setValue(mean);
                }
            }
        } else // never 
            throw new IllegalArgumentException("synchronizeRootPriorMean: distribution_class must be Poisson, NegativeBinomial, ShiftedGeometric or PointDistribution");
    }
    
    /**
     * Sets the state for displayed fields and labels for a given root prior distribution.
     * Field values are not set if they are null.
     * Field values are set (to 0.1) only if the new distribution class has
     * invalid parameter in the field.
     *
     * Used by listeners, called when
     * selection of selected model changes in the GUI.
     *
     */
    private void synchronizeRootPriorStates()
    {
        //                  a   b   p_0     p_1             mean
        // NegativeBinomial t   q   (1-q)^t t*(1-q)^t*q     t*q/(1-q)
        // ShiftedGeometric p0  q   p0      (1-p0)*(1-q)    (1-p0)/(1-q)
        // PoissonDistrib   lm  .   e^(-lm) e^{-lm)*lm      lm
        // PointDistributi  p0  .   p0      1-p0            1-p0


        // Pnt -> Geom: p0-> p0; 0->q || a'=a; b'=0
        // Pnt -> Poisson: 1-p0 -> lm || a'=mean; b'=null
        // Pnt -> NB: 0.1 -> t, ??->q by same mean || a'=0.1; b'=mean/(mean+a')
        //
        // Geom -> NB: (1-p0)/q -> t; q->q || a'=
        // Geom -> Pnt: p0 -> p0, 0->q || a'=a; b'=0
        // Geom -> Poisson : 1-p0/1-q -> lm; null->b || a'=mean; b'=null
        //
        // NB -> Geom: (1-q)^t -> p0; ??->q by same mean
        // NB -> Pnt: (1-q)^t -> p0, 0->b
        // NB -> Poisson: mean -> lm

        // Poisson -> Point: e-lm -> p0; 0->b
        // Poisson -> NB: 0.1->t, ??->q by same mean
        // Poisson -> Geom: e^-lm -> p0; 1-(1-e^-lm)/lm -> q
        //
        // p0 = 1-(1-q)*mean

        Class<?> distribution_class = getSelectedRootPriorClass();
        
        
        Number m = "\u221e".equals(root_meanT.getValue())
        		?Double.POSITIVE_INFINITY
        		:(Number)root_meanT.getValue();

        boolean stationary_root = (uniform_duplicationCB.isSelected() && uniform_gainCB.isSelected());

        		
        if (distribution_class==Poisson.class)
        {
            root_distributionL.setText("Poisson");
            root_aL.setText("r"); // lambda
            root_aT.setInputVerifier(InputVerifiers.getPositiveInstance());
            root_aT.setToolTipText("Poisson parameter r (expected value)");
            root_bL.setText("nosuch");
            root_bT.setVisible(false);
            root_bL.setVisible(false);

            // expected value is lambda
            if (stationary_root)
            {
            	double r = lineage_gainP.getField(0).doubleValue();
            	root_aT.setValue(r);
            	root_meanT.setValue(r);
            	m = r;
            } else
            {
            	root_aT.setValue(m);
            }
            root_bT.setValue(null);
        	root_aT.setEditable(!stationary_root);
        } else if (distribution_class==NegativeBinomial.class)
        {
            root_distributionL.setText("Negative binomial");
            root_aL.setText("κ"); // kappa
            root_aT.setInputVerifier(InputVerifiers.getPositiveInstance());
            root_aT.setToolTipText("stopping-time parameter κ for negative binomial");
            root_bL.setText("q");
            root_bT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
            root_bT.setToolTipText("sucess-probability parameter q for negative binomial");
            root_bT.setVisible(true);
            root_bL.setVisible(true);
            
            if (stationary_root)
            {
            	double κ = lineage_gainP.getField(0).doubleValue();
            	double q = lineage_duplicationP.getField(0).doubleValue();
            	root_aT.setValue(κ);
            	root_bT.setValue(q);
            	m = κ*q/(1.-q);
            	//root_meanT.setValue(m); // reset by value property change event on root_aT and root_bT
            } else
	            if (m!=null)
	            {
	                // expected value = q t/(1-q) ; approaches Poisson with same mean
	                // check if we have a value
	                Number b = (Number)root_bT.getValue();
	                double mean = m.doubleValue();
	
	                if (b==null || b.doubleValue()==0.0)
	                {
	                    double def_theta = 0.1;
	                    root_aT.setValue(def_theta);
	                } else
	                {
	                    // from shifted geometric
	                    double q = b.doubleValue();
	                    root_aT.setValue(mean*(1.-q)/q);
	                }
	
	                Number t = (Number)root_aT.getValue();
	                root_bT.setValue(mean/(mean+t.doubleValue()));
	
	            }
            root_aT.setEditable(!stationary_root);
            root_bT.setEditable(!stationary_root);
            
        } else if (distribution_class == PointDistribution.class)
        {
            // expected values 1-p0
            root_distributionL.setText("Point");
            root_aL.setText("p0");
            root_aT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
            root_aT.setToolTipText("failure-probability parameter p0 for Bernoulli");
            root_bT.setVisible(false);
            root_bL.setVisible(false);


            if (m!=null)
            {
                Number b = (Number)root_bT.getValue();
                if (b==null) // was poisson
                {
                    root_aT.setValue(Math.exp(-m.doubleValue()));
                } else
                {
                    double q = b.doubleValue();
                    double mq_1 = m.doubleValue()*(1.-q);
                    if (mq_1>=1.0)
                    {
                        root_aT.setValue(Math.exp(-mq_1));
                    } else
                    {
                        root_aT.setValue(1.0-mq_1);
                    }
                }
            }
            root_bT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
            root_bT.setValue(0.0);
            root_aT.setEditable(true);
            root_bT.setEditable(true);
            
        } else if (distribution_class == ShiftedGeometric.class) 
        {
            root_distributionL.setText("Shifted geometric");
            root_aL.setText("p0");
            root_aT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
            root_aT.setToolTipText("failure-probability parameter p0 for shifted geometric");
            root_bL.setText("q");
            root_bT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
            root_bT.setToolTipText("success probability parameter q for shifted geometric");
            root_bT.setVisible(true);
            root_bL.setVisible(true);

            if (stationary_root)
            {
            	double q = lineage_duplicationP.getField(0).doubleValue();
                double mq_1 = m.doubleValue()*(1.-q);
                if (mq_1>=1.0)
                {
                    root_aT.setValue(Math.exp(-mq_1));
                } else
                {
                    root_aT.setValue(1.0-mq_1);
                }
                root_bT.setValue(q);
            } else
	            if (m!=null)
	            {
	                Number b = (Number) root_bT.getValue();
	                if (b==null) // was Poisson: same mean, same p0
	                {
	                    double mean = m.doubleValue();
	                    double p0 = Math.exp(-mean);
	                    root_aT.setValue(p0);
	                    double q = 1.0+Math.expm1(-mean)/mean;
	                    root_bT.setValue(q);
	                } else
	                {
	                    double q = b.doubleValue();
	                    double mq_1 = m.doubleValue()*(1.-q);
	                    if (mq_1>=1.0)
	                    {
	                        root_aT.setValue(Math.exp(-mq_1));
	                    } else
	                    {
	                        root_aT.setValue(1.0-mq_1);
	                    }
	
	                    Number p0 = (Number)root_aT.getValue();
	                    q = 1.0-(1.0-p0.doubleValue())/m.doubleValue();
	                    root_bT.setValue(q);
	                }
	            }
            root_aT.setEditable(true);
            root_bT.setEditable(!stationary_root);
        }
        else // never
            throw new IllegalArgumentException("synchronizeRootPriorStates: distribution_class must be Poisson, NegativeBinomial, ShiftedGeometric or PointDistribution");

        Number a = (Number)root_aT.getValue();

        if (a!=null)
        {
            if (!root_aT.getInputVerifier().verify(root_aT))
                root_aT.setValue(0.1);


        }
        if (root_bT.isVisible())
        {
            Number b = (Number)root_bT.getValue();
            if (b!=null)
            {
                if (!root_bT.getInputVerifier().verify(root_bT))
                    root_bT.setValue(0.1);
            }
        }

        //synchronizeRootPriorMean();
    }

    private class LineageSpecificParameters
    {
    	LineageSpecificParameters(IndexedTree tree, InputVerifier input_verifier)
    	{
            initFields(tree, input_verifier);
    	}

        private CheckSelectAll<JCheckBox> fixed_parameters;
        private ParameterField[] lineage_parameters;
        private boolean variation_ok = true;
        private String parameter_text;
        private String[] lineage_name;
        private String disabled_reason;
    	
        private void initFields(IndexedTree tree, InputVerifier verifier)
        {
            int num_nodes = tree.getNumNodes();
            assert (tree.isRoot(num_nodes-1)); // root rates are not displayed 
            lineage_parameters = new ParameterField[num_nodes-1]; 
            fixed_parameters = new CheckSelectAll<JCheckBox>("")
            {
                @Override
                public String getToolTipText(java.awt.event.MouseEvent e)
                {
                    StringBuffer master_tooltip = new StringBuffer(parameter_text);
                    if (fixed_parameters.isSelected())
                    {
                        master_tooltip.append(" is kept constant during optimization everywhere");
                    } else
                    {
                        master_tooltip.append(" may be optimized");
                    }

                    if (!fixed_parameters.isEnabled())
                    {
                        master_tooltip.append("; disabled because ");
                        master_tooltip.append(disabled_reason);
                    }
                    return master_tooltip.toString();
                }
            };
            fixed_parameters.setToolTipText("Enables/disables optimization"); // need to set, so that Java calls getToolTipText
            NumberFormat fieldF = new DecimalFormat("0.############");
            fieldF.setParseIntegerOnly(false);
            lineage_name = new String[lineage_parameters.length];
            

            for (int node=0; node<lineage_parameters.length; ++node)
            {
                ParameterField f = new ParameterField(fieldF, 0.0, "");
                f.setInputVerifier(verifier);
                f.setColumns(16);
                lineage_parameters[node] = f;
                fixed_parameters.addCheckBox(f.fixedCB);
                lineage_name[node] = tree.getIdent(node);
            }
            fixed_parameters.addItemListener(state_changed_event ->
                {
                    if (!variation_ok && ItemEvent.DESELECTED == state_changed_event.getStateChange())
                        fixed_parameters.copySelectedToAll();
                });
            lineage_parameters[0].addPropertyChangeListener("value",
            		property_change_event -> 
                {
                    if (!variation_ok)
                    {
                        Object val = property_change_event.getNewValue();
                        for (int i=1; i<lineage_parameters.length; ++i)
                            lineage_parameters[i].setValue(val);
                    }
                });
        }
    

        void freezeAll()
        {
            for (ParameterField f: lineage_parameters)
                f.freezeParameter();
            fixed_parameters.setSelected(true);
            fixed_parameters.setEnabled(false);
        }
    
        ParameterField getField(int node) { return lineage_parameters[node];}
        JCheckBox getCheckBox(int node) { return lineage_parameters[node].fixedCB;}
        int numFields() { return lineage_parameters.length;}
        void enableAll()
        {
            fixed_parameters.setEnabled(true);
            if (variation_ok)
            {
                for (ParameterField f: lineage_parameters)
                    f.enableParameter();
            } else
            {
                lineage_parameters[0].setEditable(true);
                lineage_parameters[0].fixedCB.setEnabled(false);
                for (int i=1; i<lineage_parameters.length; ++i)
                    lineage_parameters[i].freezeParameter();
                fixed_parameters.copySelectedToAll();
            }
        }
        void setDisabledReason(String disabled_reason)
        {
            this.disabled_reason = disabled_reason;
            for (ParameterField f: lineage_parameters)
                f.disabled_reason = disabled_reason;
        }
        /**
         * Resets all parameter field values.
         */
        void resetAll()
        {
            for (ParameterField f:lineage_parameters)
                f.resetValue();
        }
        /**
         * Enables or disables lineage-specific variation
         * by a checkbox state; equalizes rates if checkbox is set
         *
         * @param cb checkbox giving the selected state used
         */
        void disableVariationIfSelected(JCheckBox cb)
        {
            boolean b = !cb.isSelected();
            setVariationEnabled(b);
            if (!b)
                for (int i=1; i<lineage_parameters.length; ++i)
                    lineage_parameters[i].disabled_reason = "lineage-specific rates are not allowed";
			if (cb.isSelected()) 
				equalizeValues(lineage_lengthP);
        }
        
        /**
         * Enables/disables lineage-specific variation.
         * Call after {@link #freezeAll()} or {@link #enableAll() }.
         *
         * @param variation_enabled
         */
        void setVariationEnabled(boolean variation_enabled)
        {
            this.variation_ok = variation_enabled;
            if (fixed_parameters.isEnabled())
                enableAll();
        }
        
        void setValueAll(Object o)
        {
            if (variation_ok)
            {
                for (ParameterField f: lineage_parameters)
                    f.setValue(o);
            } else
                lineage_parameters[0].setValue(o); // copied to others by property change listener
        }
        
//        void listenToEqualizer(JCheckBox cb, LineageSpecificParameters lineage_lengthP)
//        {
//        	cb.addItemListener(e->
//        		{
//        			disableVariationIfSelected(cb, false); 
//        			if (cb.isSelected()) 
//        				equalizeValues(lineage_lengthP);
//        		});
//        }
        
        public void equalizeValues(LineageSpecificParameters weighting)
        {
            double total_rate = 0.0;
            double total_weight = 0.0;
            if (weighting == null)
            {
                for (int i=0; i<lineage_parameters.length; ++i)
                {
                    double wt = 1.0;
                    total_rate+=wt * ((Number)lineage_parameters[i].getValue()).doubleValue();
                    total_weight += wt;
                }
            } else
            {
                for (int i=0; i<lineage_parameters.length; ++i)
                {
                    double wt = ((Number)weighting.lineage_parameters[i].getValue()).doubleValue();
                    total_rate+=wt * ((Number)lineage_parameters[i].getValue()).doubleValue();
                    total_weight += wt;
                }

            }

            double avg_rate = total_rate/total_weight;
            setValueAll(avg_rate);
        }
    }
    
    
    
    private static class ParameterField extends JFormattedTextField
    {
        /**
         * Instantiation with 0 default value.
         *
         * @param format format for this field
         */
        public ParameterField(Format format)
        {
            this(format, 0);
        }


        /**
         * Instantiation with no check box text.
         *
         * @param format format for the parameter field
         * @param default_value default value for the field (used by {@link #resetValue() })
         */
        public ParameterField(Format format, Object default_value)
        {
            this(format, default_value, "");
        }

        /**
         * Instantiation.
         *
         * @param format field format
         * @param check_box_text text for associated check box
         * @param default_value default value for parameter (when disabled)
         */
        public ParameterField(Format format, Object default_value, String check_box_text)
        {
            super(format);
            fixedCB=new JCheckBox(check_box_text)
            {
                @Override
                public String getToolTipText(MouseEvent ignored)
                {
                    return getCheckBoxToolTip();
                }
            };
            fixedCB.setToolTipText("");
            this.default_value = default_value;
            parameter_text = "Rate parameter";
            setToolTipText("");
        }

        JCheckBox fixedCB;
        private Object default_value;
        String disabled_reason;
        String parameter_text;

        double doubleValue()
        {
            return ((Number)getValue()).doubleValue();
        }

        int intValue()
        {
            return ((Number)getValue()).intValue();
        }

        private void appendDisabledText(StringBuffer sb)
        {
            sb.append("; disabled");
            if (disabled_reason!=null)
            {
                sb.append(" because ");
                sb.append(disabled_reason);
            }
        }

        private String getCheckBoxToolTip()
        {
            StringBuffer checkbox_tooltip = new StringBuffer(parameter_text);
            if (fixedCB.isSelected())
                checkbox_tooltip.append(" is kept constant during optimization");
            else
                checkbox_tooltip.append(" will be optimized");
            if (!fixedCB.isEnabled())
                appendDisabledText(checkbox_tooltip);
            return checkbox_tooltip.toString();
        }

        @Override
        public String getToolTipText(MouseEvent ignored)
        {
            StringBuffer field_tooltip = new StringBuffer(parameter_text);
            if (!isEditable())
                appendDisabledText(field_tooltip);
            return field_tooltip.toString();
        }

        @Override
        public void setVisible(boolean b)
        {
            super.setVisible(b);
            fixedCB.setVisible(b);
        }


        /**
         * Enables this parameter. The checkbox state is unchanged,
         * but the text field is set to editable and
         * both the field and the check box are enabled.
         *
         */
        void enableParameter()
        {
            fixedCB.setEnabled(true);
            setEditable(true);
            //System.out.println("#*"+getClass().getName()+".enable "+parameter_text+" fixed "+fixedCB.isSelected());
        }


        /**
         * Disables user contol over this parameter.
         * The checkbox state is set to selected,
         * the checkbox is disabled, and the field is set to uneditable.
         *
         */
        void freezeParameter()
        {
            fixedCB.setEnabled(false);
            fixedCB.setSelected(true);
            setEditable(false);
        }

        /**
         * Resets to default value.
         */
        void resetValue()
        {
           setValue(default_value);
        }
        
        @Override
        public Dimension getMaximumSize()
        {
        	return this.getPreferredSize();
        }
        

        /**
         * Attaches a listener for a field category that may enable/disable
         * this parameter. Parameter value is not changed when
         * disabling/enabling.
         *
         * @param controlTF controlling field
         * @param disable_for_value not null
         * @param disabled_reason explanation for why this field disables the parameter
         */
        void listenTo(JFormattedTextField controlTF, final Number disable_for_value, final String disabled_reason)
        {
            controlTF.addPropertyChangeListener("value", update->
                    {
                        Number newval = (Number)update.getNewValue(); // equals does not work because the field may have Integer or Long value change
                        boolean got_disabled = (newval.intValue() == disable_for_value.intValue());
                        if (got_disabled)
                        { // no more editing
                            freezeParameter();
                            ParameterField.this.disabled_reason = disabled_reason;
                        } else
                        {
                            enableParameter();
                        }
                        fixedCB.setSelected(got_disabled);
                        //System.out.println("#*"+getClass().getName()+".pC "+parameter_text+"\t"+update.getPropertyName()+" = "+update.getNewValue()+", dis "+got_disabled);
                    });

        }

        /**
         * Attaches a listener to a checkbox or other selectable that
         * enables/disables this parameter.
         * Disabled parameter is set to default value.
         *
         * @param controlSEL controlling check box
         * @param is_enabler whether SELECTED events enable(true) or disable
         * @param disabled_reason explanation for why this button disables the parameter
         */
        void listenTo(AbstractButton controlSEL, boolean is_enabler, final String disabled_reason)
        {
            if (is_enabler)
            {
                controlSEL.addItemListener(e->
                        {
                            if (ItemEvent.SELECTED == e.getStateChange())
                                enableParameter();
                            else if (ItemEvent.DESELECTED == e.getStateChange())
                            {
                                freezeParameter();
                                setValue(default_value);
                                ParameterField.this.disabled_reason =  disabled_reason;
                            }
                        });
            } else
            {
                controlSEL.addItemListener(e->
                        {
                            if (ItemEvent.DESELECTED == e.getStateChange())
                                enableParameter();
                            else if (ItemEvent.SELECTED == e.getStateChange())
                            {
                                freezeParameter();
                                setValue(default_value);
                                ParameterField.this.disabled_reason = disabled_reason;
                            }
                        });
            }
        }
    } 
}
