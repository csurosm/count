
package ca.umontreal.iro.evolution.malin.ui.count;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import java.io.File;
import java.util.HashSet;
import java.util.Set;
import java.util.Map;
import java.util.HashMap;


import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.Timer;


import java.text.NumberFormat;
import java.text.DecimalFormat;

import ca.umontreal.iro.matek.DiscreteDistribution;
import ca.umontreal.iro.matek.NegativeBinomial;
import ca.umontreal.iro.matek.PointDistribution;
import ca.umontreal.iro.matek.Poisson;

import ca.umontreal.iro.evolution.TreeNode;
import ca.umontreal.iro.evolution.genecontent.ML;
import ca.umontreal.iro.evolution.genecontent.NodeWithRates;
import ca.umontreal.iro.evolution.genecontent.OccurrenceTable;
import ca.umontreal.iro.evolution.genecontent.RateVariation;
import ca.umontreal.iro.evolution.genecontent.HomogeneousRateVariation;
import ca.umontreal.iro.evolution.genecontent.StableComputation;
import ca.umontreal.iro.evolution.genecontent.TreeWithRates;

import ca.umontreal.iro.evolution.malin.DataFile;

import ca.umontreal.iro.evolution.malin.ui.InputVerifiers;
import ca.umontreal.iro.evolution.malin.ui.CheckSelectAll;

/**
 * 
 * GUI for rate model optimization
 *
 * @author csuros
 */
public class RateOptimizationUI implements ActionListener
{
    /**
     * Instantiation of the interface: nothing happens, only 
     * local variables are initialized.
     * 
     * @param main_tree the session's tree: will be used for default null model
     * @param data_file table of family sizes
     * @param rate_file initiali rate model (may be null)
     */
    public RateOptimizationUI(TreeWithRates main_tree, DataFile<OccurrenceTable> data_file, DataFile<RateVariation> rate_file)
    {
        this.main_tree = main_tree;
        this.data_file = data_file;
        this.rate_file = rate_file;
        this.optimized_model_display = null;
    }
    
    private TreeWithRates main_tree;
    
    private boolean optimization_is_done;
    
    /**
     * Set at instantiation
     */
    private DataFile<OccurrenceTable> data_file;
    /**
     * Set at instantiation (may be null)
     */
    private DataFile<RateVariation> rate_file;
    /**
     * Dialog for the optimization parameters: set with optimizeModel()
     */
    private OptimizationParameters parameters;
    
    /**
     * optimized_model_display is used to show the progress of the optimization, and the resulting model
     */
    private ModelDisplay optimized_model_display;
    
    private JProgressBar optimization_progress;
    
    private JLabel optimization_stage;
    
    private JFormattedTextField optimization_delta;
    
    private Thread optimization_thread;
    
    private Timer optimization_timer;
    
    private JButton cancel_button ;
    
    private JButton snapshot_button;
    
    /**
     * Brings up a dialog for setting the optimizaion parameters. If cancelled, returns null.
     * Otherwise it launches the optimization in a background thread, 
     * and constructs a RateModelDisplay for the optimized model that tracks the 
     * optimization progress, which is returned. The returned display will fire 
     * a property change (<q>done</q> from false to true) when the optimization
     * finishes eithe rproperly or by cancelation. The underlying rate variation model 
     * is hidden until the optimization is done (see getRateModel()).  
     * 
     * @param dealer parent DealerCount application
     * @return null if canceled, or a GUI element for the optimization
     */
    public ModelDisplay optimizeModel(DealerCount dealer)
    {
        ModelSelectionPane parP = new ModelSelectionPane(dealer);
        parP.addNullModel(main_tree);
        if (rate_file!=null)
            parP.addModel(rate_file);
        parP.showModelOptimizationDialog(100, 0.01);


//        JFrame frame = dealer.getTopFrame();
//        parameters = new OptimizationParameters(dealer);
//        Dimension frameD = frame.getSize();
//
//        parameters.pack();
//        parameters.setBounds((int)(0.05*frameD.width),(int)(0.05*frameD.height),(int)(0.9*frameD.width),(int)(0.9*frameD.height));
//
//        parameters.setVisible(true);
        // now a modal dialog appears and we get back 
        // with the execution here when the dialog is disposed of


        DataFile<RateVariation> optimized_model = parP.getOptimizedModelFile();
        if (optimized_model == null)
            optimized_model_display = null;
        else
        {
            ML optimization = parP.getOptimization();
            optimized_model_display = new ModelDisplay(optimized_model,optimization,!parP.isDefaultModelSelected());
        }
        if (optimized_model_display != null)
            launchOptimizationProcess(dealer);
            
        return optimized_model_display; // which was set through the dialog
    }
    
    
    /**
     * This is where the actual optimization process starts in the background
     */
    private void launchOptimizationProcess(DealerCount dealer)
    {
        optimization_progress.setMaximum(optimized_model_display.optimization.getMaxRounds());
        optimization_progress.setIndeterminate(true);
        cancel_button.addActionListener(this);
        snapshot_button.addActionListener(this);

        optimization_is_done = false;
        
        optimization_thread = new Thread(new Runnable() // cannot use SwingWorker because optimize() can be used from command-;line programs too
        {
            @Override
            public void run()
            {
                ((StableComputation) optimized_model_display.optimization.getComputationModel()).allocateDataStructures();
                optimized_model_display.optimization.optimize();
            }
        });
        dealer.makeThreadExceptionHandler(optimization_thread); // for strange errors...
        optimization_thread.start();
        optimization_timer = new Timer(100,this); // start in 0.1 seconds
        optimization_timer.start();        
    }
    
    private void optimizationDone()
    {
        optimization_is_done = true;
        optimized_model_display.firePropertyChange("done", false, true);
    }
            
    @Override
    public void actionPerformed(ActionEvent E)
    {
        Object src = E.getSource();
        if (src == cancel_button)
        {
            // progress update
            optimization_thread.interrupt();
            optimization_timer.stop();
            cancel_button.setEnabled(false);
            snapshot_button.setEnabled(false);
            optimized_model_display.updateDisplay(true);

            optimizationDone();
        } else if (src == snapshot_button)
        {
            optimized_model_display.snapshot();
        }
        else if (src == optimization_timer)
        {
            // progress update
            ML optimization = optimized_model_display.optimization;
            StableComputation SC = (StableComputation) optimization.getComputationModel();
            int n_alloc = SC.getNumProfilesAllocated();
            int n_profiles =data_file.getData().getNumFamilies() ;
            if (n_alloc == n_profiles)
            {
                optimization_progress.setIndeterminate(false);
                double eps = optimization.getOptimizationCurrentLikelihoodDrop();
                double lik = optimization.getOptimizationCurrentLogLikelihood();
                String stage = optimization.getOptimizationStage();
                int round = optimization.getOptimizationRound();
                if (optimization_progress.getValue()<round)
                    optimization_progress.setValue(round);

                int max_rounds = optimization.getMaxRounds();
                optimization_progress.setString("Round "+round+" of "+max_rounds);
                if (round == max_rounds)
                {
                    // done
                    ((Timer)src).stop(); // optimization_timer
                    optimization_delta.setEnabled(false);
                    //optimization_progress.setString(null);
                    optimization_stage.setText(" Optimization done. LL="+lik);
                    optimization_progress.setEnabled(false);
                    cancel_button.setEnabled(false);
                    cancel_button.setVisible(false);
                    snapshot_button.setEnabled(false);
                    snapshot_button.setVisible(false);
                    optimized_model_display.updateDisplay(true);

                    optimizationDone();
                } else
                {
                    optimization_delta.setValue(new Double(eps));
                    optimization_stage.setText(" "+stage+". LL="+lik);
                    optimized_model_display.updateDisplay(false);
                }
            }
            else
            { // still allocating memory
                optimization_stage.setText("Allocating memory: "+n_alloc+"/"+n_profiles);
            }
        }
    }
    


    /**
     * Local class for setting the optimization parameters through a JDialog.
     *
     * States for the graphical components:
     * <ul>
     * <li> {@link #uniform_duplicationCB} enabled if model has duplication; disabled and set if model has no duplication
     * </ul>
     *
     *
     */
    private class OptimizationParameters extends JDialog implements ActionListener, ItemListener, PropertyChangeListener
    {
        /**
         * Initializes a modal dialog.
         * 
         * @param frame parent JFrame
         */
        private OptimizationParameters(DealerCount dealer)
        {
            super(dealer.getTopFrame(),"Rate optimization", true);
            this.dealer = dealer;
            initModels();
            initComponents();
            initListeners();
            initModelParameters(rate_file != null);
//            setInitialFeatureComponentStates(rate_file != null);
//            addUniformRateListeners();
        }
        private DealerCount dealer;

        private RateVariation null_model;
        private RateVariation alternative_model;
        
        private JTabbedPane tabs;
        private JPanel button_pane;
        private JButton start_button;
        private JButton cancel_button;
        
        private JPanel model_type_panel;
        private JPanel model_parameters_panel;
        //private JPanel rate_variation_panel;
        //private JPanel advanced_options_panel;
        
        
        /**
         * Radio buttons for starting model
         */
        private JRadioButton null_modelRB;
        private JRadioButton alternative_modelRB;
        
        
        /**
         * Radio buttons for model type
         */
        private ModelTypeRB model_gldB;
        private ModelTypeRB model_dlB;
        private ModelTypeRB model_glB;
        private ModelTypeRB model_plB;
        
        /**
         * Checkboxes for lineage-specific variation
         */
        private JCheckBox uniform_duplicationCB;
        private JCheckBox uniform_gainCB;
        private JCheckBox homogeneous_lossCB;
        
        /**
         * Radio buttons for distribution at root
        */
        private JRadioButton root_poissonB;
        private JRadioButton root_negbinB;
        private JRadioButton root_pointB;
        private JRadioButton root_stationaryB;
        
       
        /**
         * Gamma categories
         */
        private JFormattedTextField gamma_lengthT;
        private JFormattedTextField gamma_gainT;
        private JFormattedTextField gamma_lossT;
        private JFormattedTextField gamma_duplicationT;
        
        /**
         * Checkboxes for forbidden categories
         */
        //private JCheckBox forbidden_lossCB;
        private JCheckBox forbidden_gainCB;
        private JCheckBox forbidden_duplicationCB;
        
        /**
         * Fields for convergence criteria
         */
        private JFormattedTextField epsT;
        private JFormattedTextField roundT;
        
        /**
         * Fields for gamma parameters
         */
        private JFormattedTextField alpha_lengthT;
        private JFormattedTextField alpha_gainT;
        private JFormattedTextField alpha_lossT;
        private JFormattedTextField alpha_duplicationT;
        
        /**
         * Unoptimized gamma parameters
         */
        private CategoryParametersCB fixed_alpha_lengthCB;
        private CategoryParametersCB fixed_alpha_gainCB;
        private CategoryParametersCB fixed_alpha_duplicationCB;
        private CategoryParametersCB fixed_alpha_lossCB;
        
        /**
         * Proportion of no-duplication and no-gain families
         */
        private JFormattedTextField forbidden_duplicationT;
        private JFormattedTextField forbidden_gainT;
        
        private CategoryParametersCB fixed_forbidden_gainCB;
        private CategoryParametersCB fixed_forbidden_duplicationCB;
        
        /**
         * Root distribution
         */
        private JFormattedTextField root_aT;
        private JFormattedTextField root_bT;
        
        private JCheckBox fixed_root_CB;
        
        private JLabel root_distributionL;
        private JLabel root_aL;
        private JLabel root_bL;
        /**
         * Fields for lineage-specific rates
         */
        private JFormattedTextField[] edge_lengthT;
        private JFormattedTextField[] loss_rateT;
        private JFormattedTextField[] gain_rateT;
        private JFormattedTextField[] duplication_rateT;
        
        private JCheckBox[] fixed_lengthCB;
        private JCheckBox[] fixed_lossCB;
        private JCheckBox[] fixed_gainCB;
        private JCheckBox[] fixed_duplicationCB;
        
//        private CheckSelectAll all_fixed_lengthCB;
//        private CheckSelectAll all_fixed_gainCB;
//        private CheckSelectAll all_fixed_duplicationCB;
//        private CheckSelectAll all_fixed_lossCB;

        private LineageSpecificParametersCB lineage_lossP;
        private LineageSpecificParametersCB lineage_gainP;
        private LineageSpecificParametersCB lineage_duplicationP;
        private LineageSpecificParametersCB lineage_lengthP;

        /**
         * Sets up {@link #null_model} and {@link #alternative_model}.
         */
        private void initModels()
        {
            {
                TreeWithRates copy_tree = new TreeWithRates(NodeWithRates.copyTree(main_tree.getRoot()));
                null_model = new RateVariation(copy_tree, new Poisson(0.1), 1,1,1,1);
            }
            if (rate_file==null)
                alternative_model = null;
            else
            {
                RateVariation original_model = rate_file.getData();
                TreeWithRates copy_tree = new TreeWithRates(NodeWithRates.copyTree(original_model.getMainTree().getRoot()));
                alternative_model = original_model.sameModelForDifferentTree(copy_tree);
            }
        }

        /**
         * Sets up the Swing components for model architecture and
         * parameters, as well as the
         * tabbed layout and the start/cancel buttons.
         */
        private void initComponents()
        {
            tabs = new JTabbedPane();
            button_pane = new JPanel();
            button_pane.setLayout(new BoxLayout(button_pane,BoxLayout.LINE_AXIS));
            cancel_button = new JButton("Cancel");
            cancel_button.addActionListener(this);
            
            start_button = new JButton("Perform optimization");
            start_button.addActionListener(this);
            
            
            button_pane.add(Box.createHorizontalGlue());
            button_pane.add(cancel_button);
            button_pane.add(Box.createRigidArea(new Dimension(10,0)));
            button_pane.add(start_button);
            button_pane.setBorder(BorderFactory.createEmptyBorder(0, 10, 10, 10));
            
            model_parameters_panel = createModelParametersPanel();
            model_type_panel = createModelTypePanel();
            //rate_variation_panel = createRateVariationPanel();
            //advanced_options_panel = createAdvancedOptionsPanel();

            tabs.addTab("Model type", model_type_panel);
            tabs.addTab("Model parameters",new JScrollPane(model_parameters_panel));
            //tabs.addTab("Rate variation", rate_variation_panel);
            //tabs.addTab("Advanced options", advanced_options_panel);
            
            setLayout(new BorderLayout());
            add(tabs, BorderLayout.CENTER);
            add(button_pane,BorderLayout.PAGE_END);
        }
        

        private void initListeners()
        {
            fixed_alpha_lengthCB.listen(gamma_lengthT);
            fixed_alpha_lossCB.listen(gamma_lossT);
            fixed_alpha_duplicationCB.listen(gamma_duplicationT);
            fixed_alpha_gainCB.listen(gamma_gainT);
            fixed_forbidden_duplicationCB.listen(forbidden_duplicationCB);
            fixed_forbidden_gainCB.listen(forbidden_gainCB);

            lineage_lossP.listenToUniformParameters(homogeneous_lossCB);

            lineage_duplicationP.listenToUniformParameters(uniform_duplicationCB);
            lineage_duplicationP.listenToModelArchitecture(model_gldB, true);
            lineage_duplicationP.listenToModelArchitecture(model_plB, false);
            lineage_duplicationP.listenToModelArchitecture(model_glB, false);
            lineage_duplicationP.listenToModelArchitecture(model_dlB, true);

            lineage_gainP.listenToUniformParameters(uniform_gainCB);
            lineage_gainP.listenToModelArchitecture(model_gldB, true);
            lineage_gainP.listenToModelArchitecture(model_plB, false);
            lineage_gainP.listenToModelArchitecture(model_glB, true);
            lineage_gainP.listenToModelArchitecture(model_dlB, false);

            model_gldB.linkCheckBox(forbidden_gainCB, true, false);
            model_gldB.linkCheckBox(uniform_gainCB, true, true);
            model_gldB.linkField(gamma_gainT, true);
            model_gldB.linkCheckBox(forbidden_duplicationCB, true, false);
            model_gldB.linkCheckBox(uniform_duplicationCB, true, true);
            model_gldB.linkField(gamma_duplicationT, true);

            model_glB.linkCheckBox(forbidden_gainCB, true, false);
            model_glB.linkCheckBox(uniform_gainCB, true, true);
            model_glB.linkField(gamma_gainT, true);
            model_glB.linkCheckBox(forbidden_duplicationCB, false, false);
            model_glB.linkCheckBox(uniform_duplicationCB, false, true);
            model_glB.linkField(gamma_duplicationT, false);

            model_dlB.linkCheckBox(forbidden_gainCB, false, false);
            model_dlB.linkCheckBox(uniform_gainCB, false, true);
            model_dlB.linkField(gamma_gainT, false);
            model_dlB.linkCheckBox(forbidden_duplicationCB, true, false);
            model_dlB.linkCheckBox(uniform_duplicationCB, true, true);
            model_dlB.linkField(gamma_duplicationT, true);

            model_plB.linkCheckBox(forbidden_gainCB, false, false);
            model_plB.linkCheckBox(uniform_gainCB, false, true);
            model_plB.linkField(gamma_gainT, false);
            model_plB.linkCheckBox(forbidden_duplicationCB, false, false);
            model_plB.linkCheckBox(uniform_duplicationCB, false, true);
            model_plB.linkField(gamma_duplicationT, false);
        }

        private void initModelParameters(boolean is_alternative_model)
        {
            RateVariation selected_model = (is_alternative_model?alternative_model:null_model);
            boolean has_duplication = false;
            boolean has_gain = false;

            boolean same_duplication_loss = true;
            double common_duplication_rate = 0.;
            double common_gain_rate = 0.;
            double common_loss_rate = 0.;
            NodeWithRates[] nodes = selected_model.getMainTree().getDFT();

            for (int node_idx=0; node_idx<nodes.length ; node_idx++)
            {
                NodeWithRates N = nodes[node_idx];
                if (!N.isRoot())
                {
                    double len = N.getLength();
                    double drate = N.getDuplicationRate();
                    double grate = N.getTransferRate();
                    double lrate = N.getLossRate();

                    edge_lengthT[node_idx].setValue(len);
                    duplication_rateT[node_idx].setValue(drate);
                    gain_rateT[node_idx].setValue(grate);
                    loss_rateT[node_idx].setValue(lrate);

                    has_duplication = has_duplication || (drate*len > 0.);
                    has_gain = has_gain || (grate*len > 0.);
                    if (node_idx==0)
                    {
                        common_duplication_rate = drate;
                        common_gain_rate = grate;
                        common_loss_rate = lrate;
                    } else
                    {
                        if (!Double.isNaN(common_duplication_rate) && common_duplication_rate != drate)
                            common_duplication_rate = Double.NaN;
                        if (!Double.isNaN(common_gain_rate) && common_gain_rate != grate)
                            common_gain_rate = Double.NaN;
                        if (!Double.isNaN(common_loss_rate) && common_loss_rate != lrate)
                            common_loss_rate = Double.NaN;
                    }
                    same_duplication_loss = same_duplication_loss && (lrate == drate);
                }
            }

            if (has_gain)
            {
                if (has_duplication)
                    model_gldB.setSelected(true);
                else
                    model_glB.setSelected(true);
            } else
            {
                if (has_duplication)
                {
//                    if (same_duplication_loss)
//                        dlEB.setSelected(true);
//                    else
                        model_dlB.setSelected(true);
                } else
                    model_plB.setSelected(true);
            }

            uniform_gainCB.setSelected(has_gain && !Double.isNaN(common_gain_rate));
            uniform_duplicationCB.setSelected(has_duplication && !Double.isNaN(common_duplication_rate));
            homogeneous_lossCB.setSelected(!Double.isNaN(common_loss_rate));

            DiscreteDistribution distr = selected_model.getRootPrior();
            if (distr instanceof Poisson)
            {
                root_poissonB.setSelected(true);
                setRootDistributionComponentStates(Poisson.class);
                root_aT.setValue(new Double(distr.getParameters()[0]));
            } else if (distr instanceof NegativeBinomial)
            {
                root_negbinB.setSelected(true);
                setRootDistributionComponentStates(NegativeBinomial.class);
                root_aT.setValue(new Double(distr.getParameters()[0]));
                root_bT.setValue(new Double(distr.getParameters()[1]));
            } else if (distr instanceof PointDistribution)
            {
                root_pointB.setSelected(true);
                setRootDistributionComponentStates(PointDistribution.class);
                root_aT.setValue(new Double(distr.getParameters()[0]));
            }
            // this can be optimized
            fixed_root_CB.setSelected(false);

            // set rate variation fields
            gamma_lengthT.setValue(new Integer(selected_model.getNumEdgeLengthGammaCategories()));
            gamma_gainT.setValue(new Integer(selected_model.getNumTransferRateGammaCategories()));
            gamma_lossT.setValue(new Integer(selected_model.getNumLossRateGammaCategories()));
            gamma_duplicationT.setValue(new Integer(selected_model.getNumDuplicationRateGammaCategories()));

            alpha_lengthT.setValue(new Double(selected_model.getEdgeLengthAlpha()));
            alpha_gainT.setValue(new Double(selected_model.getTransferRateAlpha()));
            alpha_lossT.setValue(new Double(selected_model.getLossRateAlpha()));
            alpha_duplicationT.setValue(new Double(selected_model.getDuplicationRateAlpha()));

            // these can be optimized
//            fixed_alpha_lengthCB.setSelected(false);
//            fixed_alpha_gainCB.setSelected(false);
//            fixed_alpha_lossCB.setSelected(false);
//            fixed_alpha_duplicationCB.setSelected(false);

//            setGammaCategoriesComponentStates();

            forbidden_gainT.setValue(new Double(selected_model.getTransferForbiddenProportion()));
            forbidden_duplicationT.setValue(new Double(selected_model.getDuplicationForbiddenProportion()));
            forbidden_gainCB.setSelected(selected_model.getTransferForbiddenProportion()!=0.0);
            forbidden_duplicationCB.setSelected(selected_model.getDuplicationForbiddenProportion()!=0.0);
            fixed_forbidden_gainCB.setSelected(!forbidden_gainCB.isSelected()); // default state is fixed if 0
            fixed_forbidden_duplicationCB.setSelected(!forbidden_duplicationCB.isSelected()); // default state is fixed if 0
        }


//        /**
//         * When a different initial model is selected,
//         * buttons and text fields are updated accordingly.
//         *
//         * @param is_alternative_model whether default null model or alternative mode is selected
//         */
//        private void setInitialFeatureComponentStates(boolean is_alternative_model)
//        {
//            RateVariation selected_model = (is_alternative_model?alternative_model:null_model);
//
//            boolean has_duplication = false;
//            boolean has_gain = false;
//
//            if (has_gain)
//            {
//                if (has_duplication)
//                    model_gldB.setSelected(true);
//                else
//                    model_glB.setSelected(true);
//            } else
//            {
//                if (has_duplication)
//                {
////                    if (same_duplication_loss)
////                        dlEB.setSelected(true);
////                    else
//                        model_dlB.setSelected(true);
//                } else
//                    model_plB.setSelected(true);
//            }
//
//            boolean same_duplication_loss = true;
//            double common_duplication_rate = 0.;
//            double common_gain_rate = 0.;
//            NodeWithRates[] nodes = selected_model.getMainTree().getDFT();
//
//            for (int node_idx=0; node_idx<nodes.length ; node_idx++)
//            {
//                NodeWithRates N = nodes[node_idx];
//                if (!N.isRoot())
//                {
//                    double len = N.getLength();
//                    double drate = N.getDuplicationRate();
//                    double grate = N.getTransferRate();
//                    double lrate = N.getLossRate();
//
//                    edge_lengthT[node_idx].setValue(len);
//                    duplication_rateT[node_idx].setValue(drate);
//                    gain_rateT[node_idx].setValue(grate);
//                    loss_rateT[node_idx].setValue(lrate);
//
//                    has_duplication = has_duplication || (drate*len > 0.);
//                    has_gain = has_gain || (grate*len > 0.);
//                    if (node_idx==0)
//                    {
//                        common_duplication_rate = drate;
//                        common_gain_rate = grate;
//                    } else
//                    {
//                        if (!Double.isNaN(common_duplication_rate) && common_duplication_rate != drate)
//                            common_duplication_rate = Double.NaN;
//                        if (!Double.isNaN(common_gain_rate) && common_gain_rate != grate)
//                            common_gain_rate = Double.NaN;
//                    }
//                    same_duplication_loss = same_duplication_loss && (lrate == drate);
//                }
//            }
//
//            uniform_gainCB.setSelected(has_gain && !Double.isNaN(common_gain_rate));
//            uniform_duplicationCB.setSelected(has_duplication && !Double.isNaN(common_duplication_rate));
//            homogeneous_lossCB.setSelected(false);
//
//            //System.out.println("#*ROUI.sIFB "+is_alternative_model+"\tdup "+has_duplication+"\tgn "+has_gain+"\tsame_dl "+same_duplication_loss+"\tcdup "+common_duplication_rate+"\tgn "+common_gain_rate);
//
//
////            adjustEditableGainGadgets();
////            adjustEditableDuplicationGadgets();
//
//            // set the root distribution fields
//            DiscreteDistribution distr = selected_model.getRootPrior();
//            if (distr instanceof Poisson)
//            {
//                root_poissonB.setSelected(true);
//                setRootDistributionComponentStates(Poisson.class);
//                root_aT.setValue(new Double(distr.getParameters()[0]));
//            } else if (distr instanceof NegativeBinomial)
//            {
//                root_negbinB.setSelected(true);
//                setRootDistributionComponentStates(NegativeBinomial.class);
//                root_aT.setValue(new Double(distr.getParameters()[0]));
//                root_bT.setValue(new Double(distr.getParameters()[1]));
//            } else if (distr instanceof PointDistribution)
//            {
//                root_pointB.setSelected(true);
//                setRootDistributionComponentStates(PointDistribution.class);
//                root_aT.setValue(new Double(distr.getParameters()[0]));
//            }
//            // this can be optimized
//            fixed_root_CB.setSelected(false);
//
//            // set rate variation fields
//            gamma_lengthT.setValue(new Integer(selected_model.getNumEdgeLengthGammaCategories()));
//            gamma_gainT.setValue(new Integer(selected_model.getNumTransferRateGammaCategories()));
//            gamma_lossT.setValue(new Integer(selected_model.getNumLossRateGammaCategories()));
//            gamma_duplicationT.setValue(new Integer(selected_model.getNumDuplicationRateGammaCategories()));
//
//            alpha_lengthT.setValue(new Double(selected_model.getEdgeLengthAlpha()));
//            alpha_gainT.setValue(new Double(selected_model.getTransferRateAlpha()));
//            alpha_lossT.setValue(new Double(selected_model.getLossRateAlpha()));
//            alpha_duplicationT.setValue(new Double(selected_model.getDuplicationRateAlpha()));
//
//            // these can be optimized
////            fixed_alpha_lengthCB.setSelected(false);
////            fixed_alpha_gainCB.setSelected(false);
////            fixed_alpha_lossCB.setSelected(false);
////            fixed_alpha_duplicationCB.setSelected(false);
//
////            setGammaCategoriesComponentStates();
//
//            forbidden_gainT.setValue(new Double(selected_model.getTransferForbiddenProportion()));
//            forbidden_duplicationT.setValue(new Double(selected_model.getDuplicationForbiddenProportion()));
//            forbidden_gainCB.setSelected(selected_model.getTransferForbiddenProportion()!=0.0);
//            forbidden_duplicationCB.setSelected(selected_model.getDuplicationForbiddenProportion()!=0.0);
//            fixed_forbidden_gainCB.setSelected(!forbidden_gainCB.isSelected()); // default state is fixed if 0
//            fixed_forbidden_duplicationCB.setSelected(!forbidden_duplicationCB.isSelected()); // default state is fixed if 0
//
////            setForbiddenRateComponentStates();
//          for (int node_idx=0; node_idx<nodes.length; node_idx++)
//            {
//                NodeWithRates N = nodes[node_idx];
//                if (!N.isRoot())
//                {
//                    edge_lengthT[node_idx].setValue(nodes[node_idx].getLength());
//                    duplication_rateT[node_idx].setValue(nodes[node_idx].getDuplicationRate());
//                    gain_rateT[node_idx].setValue(nodes[node_idx].getTransferRate());
//                    loss_rateT[node_idx].setValue(nodes[node_idx].getLossRate());
////                    fixed_lengthCB[node_idx].setSelected(false);
////                    fixed_duplicationCB[node_idx].setSelected(false);
////                    fixed_gainCB[node_idx].setSelected(false);
////                    fixed_lossCB[node_idx].setSelected(true);
//                }
//            }
//
//            setModelFeatureButtons(has_duplication, has_gain);
//        }


        
//        /**
//         * Called when the value of a gamma_xxxT text field (number of discrete categories) changes
//         */
//        private void setGammaCategoriesComponentStates()
//        {
//            int cat_length = ((Number)gamma_lengthT.getValue()).intValue();
//            int cat_gain = ((Number)gamma_gainT.getValue()).intValue();
//            int cat_loss = ((Number)gamma_lossT.getValue()).intValue();
//            int cat_duplication = ((Number)gamma_duplicationT.getValue()).intValue();
//
//            alpha_lengthT.setEditable(cat_length>1);
//            alpha_gainT.setEditable(cat_gain>1);
//            alpha_lossT.setEditable(cat_loss>1);
//            alpha_duplicationT.setEditable(cat_duplication>1);
//
//            if (cat_length==1)
//                alpha_lengthCB.setSelected(true);
//            if (cat_gain==1)
//                alpha_gainCB.setSelected(true);
//            if (cat_duplication==1)
//                alpha_duplicationCB.setSelected(true);
//            if (cat_loss==1)
//                alpha_lossCB.setSelected(true);
//
//            alpha_lengthT.setEnabled(cat_length>1);
//            alpha_lengthCB.setEnabled(alpha_lengthT.isEnabled());
//            alpha_gainT.setEnabled(cat_gain>1);
//            alpha_gainCB.setEnabled(alpha_gainT.isEnabled());
//            alpha_lossT.setEnabled(cat_loss>1);
//            alpha_lossCB.setEnabled(alpha_lossT.isEnabled());
//            alpha_duplicationT.setEnabled(cat_duplication>1);
//            alpha_duplicationCB.setEnabled(alpha_duplicationT.isEnabled());
//        }
        
//        /**
//         * Called after user changes the forbidden rate checkboxes
//         *
//         */
//        private void setForbiddenRateComponentStates()
//        {
//            boolean gain_maybe_forbidden = forbidden_gainCB.isSelected();
//            forbidden_gainT.setEnabled(gain_maybe_forbidden);
//            fixed_forbidden_gainCB.setEnabled(gain_maybe_forbidden);
//            boolean duplication_maybe_forbidden = forbidden_duplicationCB.isSelected();
//            forbidden_duplicationT.setEnabled(duplication_maybe_forbidden);
//            fixed_forbidden_duplicationCB.setEnabled(duplication_maybe_forbidden);
//        }
        
        /**
         * Updates the fields when root distribution changes.
         * 
         * @param distribution_class one of Poisson, NegativeBinomial or PointDistribution
         */
        private void setRootDistributionComponentStates(Class distribution_class)
        {
            //System.out.println("#*ROUI.OP.sRDCS "+distribution_class);
            if (distribution_class==Poisson.class)
            {
                root_distributionL.setText("Poisson");
                root_aL.setText("\u03bb");
                root_aT.setValue(new Double(0.1));
                root_aT.setInputVerifier(InputVerifiers.getPositiveInstance());
                root_bL.setText("nosuch");
                root_bT.setVisible(false);
                root_bL.setVisible(false);
            } else if (distribution_class==NegativeBinomial.class)
            {
                root_distributionL.setText("Negative binomial");
                root_aL.setText("\u03b8");
                root_aT.setValue(new Double(0.1));
                root_aT.setInputVerifier(InputVerifiers.getPositiveInstance());
                root_bL.setText("q");
                root_bT.setValue(new Double(0.9));
                root_bT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
                root_bT.setVisible(true);
                root_bL.setVisible(true);
            } else if (distribution_class == PointDistribution.class)
            {
                root_distributionL.setText("Point");
                root_aL.setText("p(0)");
                root_aT.setValue(new Double(0.9));
                root_aT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
                root_bT.setVisible(false);
                root_bL.setVisible(false);
            } else
                throw new IllegalArgumentException("setRootdistributionComponentStates: distribuiton_class must be Poisson, NegativeBinomial or PointDistribution");

        }
        
//        /**
//         * Called when setting uniform duplication/gain rates: equalizes the lineage-specific values.
//         */
//        private void setUniformRateComponents(boolean adjust_for_gain, boolean adjust_for_duplication)
//        {
//            if (adjust_for_gain)
//            {
//                boolean lineage_specific_gain_ok = (gldB.isSelected() || glB.isSelected()) && !uniform_gainCB.isSelected();
//                for (int i=1; i<gain_rateT.length; i++)
//                {
//                    gain_rateT[i].setEditable(lineage_specific_gain_ok);
//                    fixed_gainCB[i].setEnabled(lineage_specific_gain_ok);
//                }
//            }
//            if (adjust_for_duplication)
//            {
//                boolean lineage_specific_duplication_ok = (gldB.isSelected() || dlB.isSelected()) && !uniform_duplicationCB.isSelected();
//                for (int i=1; i<duplication_rateT.length; i++)
//                {
//                    duplication_rateT[i].setEditable(lineage_specific_duplication_ok);
//                    fixed_duplicationCB[i].setEnabled(lineage_specific_duplication_ok);
//                }
//            }
//        }
//

//        private void adjustEditableGainGadgets()
//        {
//            boolean lineage_specific_gain_ok = (gldB.isSelected() || glB.isSelected()) && !uniform_gainCB.isSelected();
//            for (int i=0; i<gain_rateT.length; i++)
//            {
//                gain_rateT[i].setEditable(i==0 || lineage_specific_gain_ok);
//                fixed_gainCB[i].setEnabled(lineage_specific_gain_ok);
//            }
//            if (!lineage_specific_gain_ok)
//                all_fixed_gainCB.copySelectedToAll();
//        }
//
//        private void adjustEditableDuplicationGadgets()
//        {
//            boolean lineage_specific_duplication_ok = (gldB.isSelected() || dlB.isSelected()) && !uniform_duplicationCB.isSelected();
//            for (int i=0; i<duplication_rateT.length; i++)
//            {
//                duplication_rateT[i].setEditable(i==0 || lineage_specific_duplication_ok);
//                fixed_duplicationCB[i].setEnabled(lineage_specific_duplication_ok);
//            }
//            if (!lineage_specific_duplication_ok)
//                all_fixed_duplicationCB.copySelectedToAll();
//        }
//
//        private void adjustEditableLengthGadgets()
//        {
//            boolean length_editing_ok = !homogeneous_lossCB.isSelected();
//            for (int i=0; i<edge_lengthT.length; i++)
//            {
//                edge_lengthT[i].setEditable(length_editing_ok);
//                fixed_lengthCB[i].setEnabled(length_editing_ok);
//            }
//        }
//
//        private void adjustEditableLossGadgets()
//        {
//            boolean loss_editing_ok = homogeneous_lossCB.isSelected();
//            loss_rateT[0].setEditable(loss_editing_ok);
//        }
        
//        private void equalizeGainRates()
//        {
//            double sum = 0.0;
//            for (int i=0; i<gain_rateT.length; i++)
//            {
//                double x = ((Number)gain_rateT[i].getValue()).doubleValue();
//                sum += x;
//            }
//            Double x = new Double(sum/gain_rateT.length);
//            gain_rateT[0].setValue(x); // copied to all others by the listeners
////            fixed_gainCB[0].setSelected(false);
//        }
//
//        private void equalizeDuplicationRates()
//        {
//            double sum = 0.0;
//            for (int i=0; i<duplication_rateT.length; i++)
//            {
//                double x = ((Number)duplication_rateT[i].getValue()).doubleValue();
//                sum += x;
//            }
//            Double x = new Double(sum/duplication_rateT.length);
//            duplication_rateT[0].setValue(x);
////            fixed_duplicationCB[0].setSelected(false);
//        }



//        /**
//         * Called when lineage-specific value or "fixedness" changes
//         */
//        private void copyUniformRateComponents()
//        {
//            if (uniform_duplicationCB.isSelected())
//            {
//                Object X = duplication_rateT[0].getValue();
//                boolean fixed = fixed_duplicationCB[0].isSelected();
//                for (int i=1; i<duplication_rateT.length; i++)
//                {
//                    duplication_rateT[i].setValue(X);
//                    fixed_duplicationCB[i].setSelected(fixed);
//                }
//            }
//            if (uniform_gainCB.isSelected())
//            {
//                Object X = gain_rateT[0].getValue();
//                boolean fixed = fixed_gainCB[0].isSelected();
//                for (int i=1; i<gain_rateT.length; i++)
//                {
//                    gain_rateT[i].setValue(X);
//                    fixed_gainCB[i].setSelected(fixed);
//                }
//            }
//            if (homogeneous_lossCB.isSelected())
//            {
//                Object X = loss_rateT[0].getValue();
//                boolean fixed = fixed_lossCB[0].isSelected();
//                for (int i=1; i<loss_rateT.length; ++i)
//                {
//                    loss_rateT[i].setValue(X);
//                    fixed_lossCB[i].setSelected(fixed);
//                }
//            }
//        }
        
//        /**
//         * When a different model type is selected, the
//         * buttons are enabled/disabled accordingly.
//         *
//         * @param duplication_ok whether there is duplication in the model
//         * @param gain_ok whether there is gain in the model
//         */
//        private void setModelFeatureButtons(boolean duplication_ok, boolean gain_ok)
//        {
//            lineage_duplicationP.setEnabled(duplication_ok);
//            if (!duplication_ok)
//            {
//                lineage_duplicationP.disableParameters();
//                uniform_duplicationCB.setSelected(true);
//                forbidden_duplicationCB.setSelected(false);
//                gamma_duplicationT.setValue(new Integer(1));
//            }
//            uniform_duplicationCB.setEnabled(duplication_ok);
//            gamma_duplicationT.setEnabled(duplication_ok);
//            forbidden_duplicationCB.setEnabled(duplication_ok);
//
//            lineage_gainP.setEnabled(gain_ok);
//            if (!gain_ok)
//            {
//                lineage_gainP.disableParameters();
//                uniform_gainCB.setSelected(true);
//                forbidden_gainCB.setSelected(false);
//                gamma_gainT.setValue(new Integer(1));
//            }
//            uniform_gainCB.setEnabled(gain_ok);
//            gamma_gainT.setEnabled(gain_ok);
//            forbidden_gainCB.setEnabled(gain_ok);
//
//            alpha_duplicationT.setEnabled(duplication_ok);
////            alpha_duplicationCB.setEnabled(duplication_ok);
//            forbidden_duplicationT.setEnabled(duplication_ok);
//            alpha_gainT.setEnabled(gain_ok);
////            alpha_gainCB.setEnabled(gain_ok);
//            forbidden_gainT.setEnabled(gain_ok);
//
////            setGammaCategoriesComponentStates();
////            setForbiddenRateComponentStates();
//
////            all_fixed_duplicationCB.setSelected(!duplication_ok);
////            all_fixed_duplicationCB.setEnabled(duplication_ok);
//
////            for (int i=0; i<duplication_rateT.length; i++)
////            {
////                if (!duplication_ok)
////                {
////                    duplication_rateT[i].setValue(new Double(0.));
////                }
////                duplication_rateT[i].setEnabled(duplication_ok);
////                fixed_duplicationCB[i].setSelected(!duplication_ok);
////                fixed_duplicationCB[i].setEnabled(duplication_ok);
////            }
////
//////            all_fixed_gainCB.setSelected(!gain_ok);
//////            all_fixed_gainCB.setEnabled(gain_ok);
////
////            for (int i=0; i<gain_rateT.length; i++)
////            {
////                if (!gain_ok)
////                {
////                    gain_rateT[i].setValue(new Double(0.0));
////                }
////                gain_rateT[i].setEnabled(gain_ok);
////                fixed_gainCB[i].setSelected(!gain_ok);
////                fixed_gainCB[i].setEnabled(gain_ok);
////            }
//        }
        
        private JPanel createModelTypePanel()
        {
            model_gldB = new ModelTypeRB("Gain-loss-duplication (Cs\u0171r\u00f6s & Mikl\u00f3s)");
            model_gldB.setActionCommand("GLD");
//            gldB.addActionListener(this);
            model_dlB = new ModelTypeRB("Duplication-loss");
            model_dlB.setActionCommand("DL");
//            dlB.addActionListener(this);
            //dlEB = new JRadioButton("Duplication-loss with equal rates (Hahn et al.)");
            //dlEB.setActionCommand("DL=");
            //dlEB.addActionListener(this);
            model_glB = new ModelTypeRB("Gain-loss");
            model_glB.setActionCommand("GL");
//            glB.addActionListener(this);
            model_plB = new ModelTypeRB("Pure loss");
            model_plB.setActionCommand("PL");
//            plB.addActionListener(this);
            
            ButtonGroup model_typeG = new ButtonGroup();
            model_typeG.add(model_gldB);
//            model_gldB.setSelected(true);
            model_typeG.add(model_dlB);
            //model_typeG.add(dlEB);
            model_typeG.add(model_glB);
            model_typeG.add(model_plB);
            
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
                
                layout.setConstraints(model_gldB,gc);
                model_typeP.add(model_gldB);
                
                gc.gridy++;
                layout.setConstraints(model_dlB,gc);
                model_typeP.add(model_dlB);
                
                gc.gridx++;
                layout.setConstraints(model_glB,gc);
                model_typeP.add(model_glB);
                //layout.setConstraints(dlEB,gc);
                //model_typeP.add(dlEB);

                gc.gridx=0;
                gc.gridy++;
                //layout.setConstraints(glB,gc);
                //model_typeP.add(glB);
                
                //gc.gridx++;
                layout.setConstraints(model_plB,gc);
                model_typeP.add(model_plB);
                
                gc.gridheight=gc.gridy+1;
                gc.gridy=0;
                gc.gridx++;
                gc.weightx=1.0;
                gc.fill=GridBagConstraints.BOTH;
                Component fillerX = Box.createHorizontalGlue();
                layout.setConstraints(fillerX,gc);
                model_typeP.add(fillerX);
            }


            ButtonGroup starting_modelG = new ButtonGroup();
            null_modelRB = new JRadioButton("Default null model");
            starting_modelG.add(null_modelRB);
            null_modelRB.setActionCommand("null-model");
            null_modelRB.addActionListener(this);
            null_modelRB.setSelected(rate_file==null);
            if (rate_file != null)
            {
                alternative_modelRB = new JRadioButton(rate_file.getFile().getName());
                starting_modelG.add(alternative_modelRB);
            } else
            {
                alternative_modelRB = new JRadioButton("[no alternative model]");
                alternative_modelRB.setVisible(false);
                starting_modelG.add(alternative_modelRB);
            }
            alternative_modelRB.setActionCommand("alternative-model");
            alternative_modelRB.addActionListener(this);
            alternative_modelRB.setSelected(rate_file != null);
            
            
            JPanel starting_modelP = new JPanel();
            starting_modelP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Starting model",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
            {
                GridBagLayout layout = new GridBagLayout();
                GridBagConstraints gc = new GridBagConstraints();
                starting_modelP.setLayout(layout);
            
                gc.fill = GridBagConstraints.HORIZONTAL; 
                gc.gridx = 0;
                gc.gridy = 0;
                gc.gridheight=1;
                gc.gridwidth=1;
                gc.anchor = GridBagConstraints.WEST;
                gc.weightx = 0.1;
                gc.weighty = 0.1;
                
                layout.setConstraints(null_modelRB,gc);
                starting_modelP.add(null_modelRB);
                
                gc.gridx++;
                layout.setConstraints(alternative_modelRB,gc);
                starting_modelP.add(alternative_modelRB);
                                
                gc.gridheight=gc.gridy+1;
                gc.gridy=0;
                gc.gridx++;
                gc.weightx=1.0;
                gc.fill=GridBagConstraints.BOTH;
                Component fillerX = Box.createHorizontalGlue();
                layout.setConstraints(fillerX,gc);
                starting_modelP.add(fillerX);
            }
            
            uniform_duplicationCB = new JCheckBox("Same duplication-loss ratio in all lineages");
//            uniform_duplicationCB.setSelected(false);
            uniform_gainCB = new JCheckBox("Same gain-loss ratio in all lineages");
//            uniform_gainCB.setSelected(false);
            homogeneous_lossCB = new JCheckBox("Same loss rate in all lineages");
//            homogeneous_lossCB.setSelected(false);
            
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
                layout.setConstraints(homogeneous_lossCB,gc);
                lineage_variationP.add(homogeneous_lossCB);

                gc.gridheight=gc.gridy+1;
                gc.gridx++;
                gc.weightx=1.0;
                gc.fill=GridBagConstraints.BOTH;
                Component fillerX = Box.createHorizontalGlue();
                layout.setConstraints(fillerX,gc);
                lineage_variationP.add(fillerX);
            }
                
            root_poissonB = new JRadioButton("Poisson");
            root_poissonB.setActionCommand("root-Poisson");
            root_poissonB.addActionListener(this);
            root_negbinB = new JRadioButton("Negative binomial (P\u00f3lya)");
            root_negbinB.setActionCommand("root-negbin");
            root_negbinB.addActionListener(this);
            root_pointB = new JRadioButton("Bernoulli");
            root_pointB.setActionCommand("root-point");
            root_pointB.addActionListener(this);
            root_stationaryB = new JRadioButton("Stationary");
            root_stationaryB.setActionCommand("root-stationary");
            root_stationaryB.addActionListener(this);
            root_stationaryB.setEnabled(false);
            
            ButtonGroup rootG = new ButtonGroup();
            rootG.add(root_poissonB);
            rootG.add(root_negbinB);
            rootG.add(root_pointB);
            rootG.add(root_stationaryB);

//            root_poissonB.setSelected(true);

            JPanel rootP = new JPanel();
            rootP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Family size distribution at root",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
            {
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
                layout.setConstraints(root_poissonB,gc);
                rootP.add(root_poissonB);
                
                gc.gridx++;
                layout.setConstraints(root_negbinB, gc);
                rootP.add(root_negbinB);
                
                gc.gridx++;
                layout.setConstraints(root_pointB, gc);
                rootP.add(root_pointB);

                gc.gridx+=2;
                layout.setConstraints(root_stationaryB, gc);
                rootP.add(root_stationaryB);
                
                gc.gridheight=gc.gridy+1;
                gc.gridx++;
                gc.weightx=1.0;
                gc.fill=GridBagConstraints.BOTH;
                Component fillerX = Box.createHorizontalGlue();
                layout.setConstraints(fillerX,gc);
                rootP.add(fillerX);
            }
                        
            gamma_lengthT = new JFormattedTextField(NumberFormat.getIntegerInstance());
            gamma_lengthT.setInputVerifier(InputVerifiers.getPositiveInstance());
            gamma_lengthT.setValue(1);
            gamma_lengthT.setColumns(2);
            gamma_lengthT.addPropertyChangeListener("value", this);
            //gamma_lengthT.setFocusLostBehavior(JFormattedTextField.COMMIT_OR_REVERT);
            gamma_gainT = new JFormattedTextField(NumberFormat.getIntegerInstance());
            gamma_gainT.setInputVerifier(InputVerifiers.getPositiveInstance());
            gamma_gainT.setValue(1);
            gamma_gainT.setColumns(2);
            gamma_gainT.addPropertyChangeListener("value", this);
            gamma_duplicationT = new JFormattedTextField(NumberFormat.getIntegerInstance());
            gamma_duplicationT.setInputVerifier(InputVerifiers.getPositiveInstance());
            gamma_duplicationT.setValue(1);
            gamma_duplicationT.setColumns(2);
            gamma_duplicationT.addPropertyChangeListener("value", this);
            gamma_lossT = new JFormattedTextField(NumberFormat.getIntegerInstance());
            gamma_lossT.setInputVerifier(InputVerifiers.getPositiveInstance());
            gamma_lossT.setValue(1);
            gamma_lossT.setColumns(2);
            gamma_lossT.addPropertyChangeListener("value", this);
            
            JLabel gamma_lengthL = new JLabel("Gamma categories");
            JLabel gamma_gainL = new JLabel("Gamma categories");
            JLabel gamma_duplicationL = new JLabel("Gamma categories");
            JLabel gamma_lossL = new JLabel("Gamma categories");

            
            JLabel variation_lengthL = new JLabel("Edge length");
            JLabel variation_lossL = new JLabel("Loss rate");
            JLabel variation_gainL = new JLabel("Gain rate");
            JLabel variation_duplicationL = new JLabel("Duplication rate");

            gamma_lossT.setEditable(false);
            gamma_lossL.setVisible(false);
            gamma_lossT.setVisible(false);
            variation_lossL.setVisible(false);

            
            //forbidden_lossCB = new JCheckBox("No-loss category");
            //forbidden_lossCB.setSelected(false);
            forbidden_gainCB = new JCheckBox("No-gain category");
//            forbidden_gainCB.setSelected(false);
//            forbidden_gainCB.addItemListener(this);
            forbidden_duplicationCB = new JCheckBox("No-duplication category");
//            forbidden_duplicationCB.setSelected(false);
//            forbidden_duplicationCB.addItemListener(this);
            
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
                layout.setConstraints(gamma_lengthT, gc);
                family_variationP.add(gamma_lengthT);
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
            
            JLabel epsL = new JLabel("Convergence threshold on the likelihood");
            JLabel roundL = new JLabel("Maximum number of optimization rounds");
            NumberFormat epsF = new DecimalFormat("0.############");
            epsF.setParseIntegerOnly(false);
            epsT = new JFormattedTextField(epsF);
            //epsT.addPropertyChangeListener("value",this);
            epsT.setInputVerifier(InputVerifiers.getPositiveInstance());
            epsT.setValue(0.1);
            epsT.addPropertyChangeListener("value", this);
            epsT.setColumns(16);
            
            NumberFormat roundF = NumberFormat.getIntegerInstance();
            roundT = new JFormattedTextField(roundF);
            //roundT.addPropertyChangeListener("value",this);
            roundT.setInputVerifier(InputVerifiers.getPositiveInstance());
            roundT.setValue(100);
            roundT.addPropertyChangeListener("value", this);
            roundT.setColumns(10);
            JPanel convergenceP = new JPanel();
            convergenceP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Convergence criteria",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
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
                layout.setConstraints(roundL,gc);
                convergenceP.add(roundL);

                gc.gridy=1;
                layout.setConstraints(epsL,gc);
                convergenceP.add(epsL);

                gc.gridx=1;
                gc.gridy=0;
                gc.anchor = GridBagConstraints.WEST;
                gc.fill=GridBagConstraints.NONE;
                gc.weightx=0.0;
                layout.setConstraints(roundT,gc);
                convergenceP.add(roundT);

                gc.gridy=1;
                layout.setConstraints(epsT,gc);
                convergenceP.add(epsT);

                gc.gridx=2;
                gc.gridy=0;
                gc.gridheight=2;
                gc.weightx=1.0;
                gc.fill=GridBagConstraints.BOTH;
                Component fillerX = Box.createHorizontalGlue();
                layout.setConstraints(fillerX,gc);
                convergenceP.add(fillerX);          
            }
            
            JPanel main = new JPanel();
            main.setLayout(new BoxLayout(main, BoxLayout.PAGE_AXIS));
            main.add(starting_modelP);
            main.add(model_typeP);
            main.add(rootP);
            main.add(lineage_variationP);
            main.add(family_variationP);
            main.add(convergenceP);
            
            return main;
        }

        private JPanel createModelParametersPanel()
        {
            NumberFormat alpha_lengthF = new DecimalFormat("0.############");
            alpha_lengthF.setParseIntegerOnly(false);
            alpha_lengthT = new JFormattedTextField(alpha_lengthF);
            alpha_lengthT.setInputVerifier(InputVerifiers.getPositiveInstance());
            alpha_lengthT.addPropertyChangeListener("value", this);
            alpha_lengthT.setColumns(16);
            NumberFormat alpha_gainF = new DecimalFormat("0.############");
            alpha_gainF.setParseIntegerOnly(false);
            alpha_gainT = new JFormattedTextField(alpha_gainF);
            alpha_gainT.setInputVerifier(InputVerifiers.getPositiveInstance());
            alpha_gainT.addPropertyChangeListener("value", this);
            alpha_gainT.setColumns(16);
            NumberFormat alpha_duplicationF = new DecimalFormat("0.############");
            alpha_duplicationF.setParseIntegerOnly(false);
            alpha_duplicationT = new JFormattedTextField(alpha_duplicationF);
            alpha_duplicationT.setInputVerifier(InputVerifiers.getPositiveInstance());
            alpha_duplicationT.addPropertyChangeListener("value", this);
            alpha_duplicationT.setColumns(16);
            NumberFormat alpha_lossF = new DecimalFormat("0.############");
            alpha_lossF.setParseIntegerOnly(false);
            alpha_lossT = new JFormattedTextField(alpha_lossF);
            alpha_lossT.setInputVerifier(InputVerifiers.getPositiveInstance());
            alpha_lossT.addPropertyChangeListener("value", this);
            alpha_lossT.setColumns(16);
            
            JLabel variation_lengthL = new JLabel("Edge length");
            JLabel variation_gainL = new JLabel("Gain rate");
            JLabel variation_duplicationL = new JLabel("Duplication rate");
            JLabel variation_lossL = new JLabel("Loss rate");
            
            JLabel alpha_lengthL = new JLabel("Gamma shape parameter (\u03b1)");
            JLabel alpha_gainL = new JLabel("Gamma shape parameter (\u03b1)");
            JLabel alpha_duplicationL = new JLabel("Gamma shape parameter (\u03b1)");
            JLabel alpha_lossL = new JLabel("Gamma shape parameter (\u03b1)");
            
            fixed_alpha_lengthCB = new CategoryParametersCB(alpha_lengthT, 1.0);
            fixed_alpha_duplicationCB = new CategoryParametersCB(alpha_duplicationT, 1.0);
            fixed_alpha_gainCB = new CategoryParametersCB(alpha_gainT, 1.0);
            fixed_alpha_lossCB = new CategoryParametersCB(alpha_lossT, 1.0);
            
            JLabel forbidden_duplicationL = new JLabel("Proportion of families with no duplications");
            JLabel forbidden_gainL = new JLabel("Proportion of families with no gains");
            
            NumberFormat forbidden_gainF = new DecimalFormat("0.############");
            forbidden_gainF.setParseIntegerOnly(false);
            forbidden_gainT = new JFormattedTextField(forbidden_gainF);
            forbidden_gainT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
            forbidden_gainT.setColumns(16);
            forbidden_gainT.setValue(new Double(0.0));
            forbidden_gainT.addPropertyChangeListener("value", this);
            
            NumberFormat forbidden_duplicationF = new DecimalFormat("0.############");
            forbidden_duplicationF.setParseIntegerOnly(false);
            forbidden_duplicationT = new JFormattedTextField(forbidden_duplicationF);
            forbidden_duplicationT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
            forbidden_duplicationT.setColumns(16);
            forbidden_duplicationT.setValue(new Double(0.0));
            forbidden_duplicationT.addPropertyChangeListener("value", this);
            
            fixed_forbidden_gainCB = new CategoryParametersCB(forbidden_gainT, 0.0);
            fixed_forbidden_duplicationCB = new CategoryParametersCB(forbidden_duplicationT, 0.0);

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
                layout.setConstraints(fixed_alpha_lengthCB, gc);
                rate_variationP.add(fixed_alpha_lengthCB);
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
                layout.setConstraints(fixed_alpha_lossCB, gc);
                rate_variationP.add(fixed_alpha_lossCB);
                gc.gridx++;
                
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
                layout.setConstraints(fixed_alpha_gainCB, gc);
                rate_variationP.add(fixed_alpha_gainCB);
                gc.gridx++;
                gc.fill = GridBagConstraints.HORIZONTAL;
                layout.setConstraints(forbidden_gainL, gc);
                rate_variationP.add(forbidden_gainL);
                gc.gridx++;
                gc.fill = GridBagConstraints.NONE;
                layout.setConstraints(forbidden_gainT,gc);
                rate_variationP.add(forbidden_gainT);
                gc.gridx++;
                layout.setConstraints(fixed_forbidden_gainCB, gc);
                rate_variationP.add(fixed_forbidden_gainCB);
                
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
                layout.setConstraints(fixed_alpha_duplicationCB, gc);
                rate_variationP.add(fixed_alpha_duplicationCB);
                gc.gridx++;
                gc.fill = GridBagConstraints.HORIZONTAL;
                layout.setConstraints(forbidden_duplicationL, gc);
                rate_variationP.add(forbidden_duplicationL);
                gc.gridx++;
                gc.fill = GridBagConstraints.NONE;
                layout.setConstraints(forbidden_duplicationT,gc);
                rate_variationP.add(forbidden_duplicationT);
                gc.gridx++;
                layout.setConstraints(fixed_forbidden_duplicationCB, gc);
                rate_variationP.add(fixed_forbidden_duplicationCB);
                
                gc.gridheight=gc.gridy+1;
                gc.gridx++;
                gc.gridy=0;
                gc.weightx=1.0;
                gc.fill=GridBagConstraints.BOTH;
                Component fillerX = Box.createHorizontalGlue();
                layout.setConstraints(fillerX,gc);
                rate_variationP.add(fillerX);          
            }
            
            
            NumberFormat root_aF = new DecimalFormat("0.############");
            root_aF.setParseIntegerOnly(false);
            root_aT = new JFormattedTextField(root_aF);
            //root_aT.setInputVerifier(new InputVerifiers.NonnegativeInputVerifier());
            root_aT.setColumns(16);
            root_aT.addPropertyChangeListener("value",this);
            //root_aT.setValue(new Double(0.0));
            
            NumberFormat root_bF = new DecimalFormat("0.############");
            root_bF.setParseIntegerOnly(false);
            root_bT = new JFormattedTextField(root_bF);
            //root_bT.setInputVerifier(new InputVerifiers.NonnegativeInputVerifier());
            root_bT.setColumns(16);
            root_bT.addPropertyChangeListener("value",this);
            
            root_distributionL = new JLabel("Negative binomial");
            root_aL = new JLabel("par0");
            root_bL = new JLabel("par1");
            
            fixed_root_CB = new JCheckBox("Fixed");
//            fixed_root_CB.setSelected(false);
            
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
                layout.setConstraints(fixed_root_CB, gc);
                root_distributionP.add(fixed_root_CB);
                
                gc.gridheight=gc.gridy+1;
                gc.gridx++;
                gc.gridy=0;
                gc.weightx=1.0;
                gc.fill=GridBagConstraints.BOTH;
                Component fillerX = Box.createHorizontalGlue();
                layout.setConstraints(fillerX,gc);
                root_distributionP.add(fillerX);          
            }
            
            NodeWithRates[] nodes = main_tree.getDFT();
            edge_lengthT = new JFormattedTextField[nodes.length-1];
            gain_rateT = new JFormattedTextField[nodes.length-1];
            duplication_rateT = new JFormattedTextField[nodes.length-1];
            loss_rateT = new JFormattedTextField[nodes.length-1];
            fixed_lengthCB = new JCheckBox[nodes.length-1];
            fixed_lossCB = new JCheckBox[nodes.length-1];
            fixed_gainCB = new JCheckBox[nodes.length-1];
            fixed_duplicationCB = new JCheckBox[nodes.length-1];
            JLabel[] node_nameL = new JLabel[nodes.length-1];

            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                NodeWithRates N = nodes[node_idx];
                if (!N.isRoot())
                {
                    node_nameL[node_idx] = new JLabel(LookAndFeel.getLongNodeName(main_tree, node_idx));
                    NumberFormat lengthF = new DecimalFormat("0.############");
                    lengthF.setParseIntegerOnly(false);
                    edge_lengthT[node_idx] = new JFormattedTextField(lengthF);
                    edge_lengthT[node_idx].setInputVerifier(InputVerifiers.getPositiveInstance());
                    edge_lengthT[node_idx].setColumns(16);
                    NumberFormat lossF = new DecimalFormat("0.############");
                    lossF.setParseIntegerOnly(false);
                    loss_rateT[node_idx] = new JFormattedTextField(lossF);
                    loss_rateT[node_idx].setInputVerifier(InputVerifiers.getPositiveInstance());
                    loss_rateT[node_idx].setColumns(16);
                    loss_rateT[node_idx].setEditable(false);
                    NumberFormat gainF = new DecimalFormat("0.############");
                    gainF.setParseIntegerOnly(false);
                    gain_rateT[node_idx] = new JFormattedTextField(gainF);
                    gain_rateT[node_idx].setInputVerifier(new InputVerifiers.RangeVerifier(ML.MIN_TRANSFER_RATE, ML.MAX_TRANSFER_RATE));
                    gain_rateT[node_idx].setColumns(16);
                    NumberFormat duplicationF = new DecimalFormat("0.############");
                    duplicationF.setParseIntegerOnly(false);
                    duplication_rateT[node_idx] = new JFormattedTextField(duplicationF);
                    duplication_rateT[node_idx].setInputVerifier(new InputVerifiers.RangeVerifier(ML.MIN_DUPLICATION_RATE, ML.MAX_DUPLICATION_RATE));
                    duplication_rateT[node_idx].setColumns(16);
                    fixed_lengthCB[node_idx]=new JCheckBox("");
                    fixed_lossCB[node_idx]=new JCheckBox("");
                    fixed_lossCB[node_idx].setSelected(true);
                    fixed_lossCB[node_idx].setEnabled(false);
                    fixed_duplicationCB[node_idx]=new JCheckBox("");
                    fixed_gainCB[node_idx]=new JCheckBox("");
                }
            }
            
            JLabel all_fixedL = new JLabel("all edges");
            all_fixedL.setFont(all_fixedL.getFont().deriveFont(Font.ITALIC));

//            all_fixed_lengthCB = new CheckSelectAll("");
//            all_fixed_gainCB = new CheckSelectAll("");
//            all_fixed_duplicationCB = new CheckSelectAll("");
//            all_fixed_lossCB = new CheckSelectAll("");
//
//            for (JCheckBox cb: fixed_lengthCB)
//                all_fixed_lengthCB.addCheckBox(cb);
//            for (JCheckBox cb: fixed_gainCB)
//                all_fixed_gainCB.addCheckBox(cb);
//            for (JCheckBox cb: fixed_lossCB)
//                all_fixed_lossCB.addCheckBox(cb);
//            for (JCheckBox cb: fixed_duplicationCB)
//                all_fixed_duplicationCB.addCheckBox(cb);

            lineage_lossP = new LineageSpecificParametersCB();
            lineage_lossP.addCheckBoxes(fixed_lossCB);
            lineage_lossP.addFields(loss_rateT);

            lineage_duplicationP = new LineageSpecificParametersCB();
            lineage_duplicationP.addCheckBoxes(fixed_duplicationCB);
            lineage_duplicationP.addFields(duplication_rateT);

            lineage_gainP = new LineageSpecificParametersCB();
            lineage_gainP.addCheckBoxes(fixed_gainCB);
            lineage_duplicationP.addFields(gain_rateT);

            lineage_lengthP = new LineageSpecificParametersCB();
            lineage_lengthP.addCheckBoxes(fixed_lengthCB);
            lineage_lengthP.addFields(edge_lengthT);

            lineage_lossP.setEnabled(false);

//            homogeneous_lossCB.addItemListener(lineage_lossP);
//            homogeneous_lossCB.addItemListener(this);
//            uniform_duplicationCB.addItemListener(lineage_duplicationP);
//            uniform_gainCB.addItemListener(lineage_gainP);

            
//            ActionListener set_all_fixed_length = new ActionListener(){
//                @Override
//                public void actionPerformed(ActionEvent E)
//                {
//                    boolean all_fixed = true;
//                    for (int i=0; i<fixed_lengthCB.length && all_fixed; i++)
//                        all_fixed = fixed_lengthCB[i].isSelected();
//                    all_fixed_lengthCB.setSelected(all_fixed);
//                }
//            };
//
//            for (int i=0; i<fixed_lengthCB.length; i++)
//                fixed_lengthCB[i].addActionListener(set_all_fixed_length);
//            all_fixed_lengthCB.addActionListener(new ActionListener()
//                {
//                    @Override
//                    public void actionPerformed(ActionEvent E)
//                    {
//                        boolean all_fixed = all_fixed_lengthCB.isSelected();
//                        for (int i=0; i<fixed_lengthCB.length; i++)
//                            if (fixed_lengthCB[i].isSelected() != all_fixed)
//                                fixed_lengthCB[i].setSelected(all_fixed);
//                    }
//            });
//
//            ActionListener set_all_fixed_gain = new ActionListener(){
//                @Override
//                public void actionPerformed(ActionEvent E)
//                {
//                    boolean all_fixed = true;
//                    for (int i=0; i<fixed_gainCB.length && all_fixed; i++)
//                        all_fixed = fixed_gainCB[i].isSelected();
//                    all_fixed_gainCB.setSelected(all_fixed);
//                }
//            };
//
//            for (int i=0; i<fixed_gainCB.length; i++)
//                fixed_gainCB[i].addActionListener(set_all_fixed_gain);
//            all_fixed_gainCB.addActionListener(new ActionListener()
//                {
//                    @Override
//                    public void actionPerformed(ActionEvent E)
//                    {
//                        boolean all_fixed = all_fixed_gainCB.isSelected();
//                        for (int i=0; i<fixed_gainCB.length; i++)
//                            if (fixed_gainCB[i].isSelected() != all_fixed)
//                                fixed_gainCB[i].setSelected(all_fixed);
//                    }
//            });
//
//            ActionListener set_all_fixed_duplication = new ActionListener(){
//                @Override
//                public void actionPerformed(ActionEvent E)
//                {
//                    boolean all_fixed = true;
//                    for (int i=0; i<fixed_duplicationCB.length && all_fixed; i++)
//                        all_fixed = fixed_duplicationCB[i].isSelected();
//                    all_fixed_duplicationCB.setSelected(all_fixed);
//                }
//            };
//
//            for (int i=0; i<fixed_duplicationCB.length; i++)
//                fixed_duplicationCB[i].addActionListener(set_all_fixed_duplication);
//            all_fixed_duplicationCB.addActionListener(new ActionListener()
//                {
//                    @Override
//                    public void actionPerformed(ActionEvent E)
//                    {
//                        boolean all_fixed = all_fixed_duplicationCB.isSelected();
//                        for (int i=0; i<fixed_duplicationCB.length; i++)
//                            if (fixed_duplicationCB[i].isSelected() != all_fixed)
//                                fixed_duplicationCB[i].setSelected(all_fixed);
//                    }
//            });
//

//            all_fixed_lossCB.setEnabled(false);
//            all_fixed_lossCB.setSelected(true);
//            all_fixed_lossCB.addActionListener(new ActionListener()
//            {
//                @Override
//                public void actionPerformed(ActionEvent E)
//                {
//                    boolean all_fixed = all_fixed_lossCB.isSelected();
//                    for (int i=0; i<fixed_lossCB.length; ++i)
//                        if (fixed_lossCB[i].isSelected() != all_fixed)
//                            fixed_lossCB[i].setSelected(all_fixed);
//                }
//            });
            
            // headers
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
            
            JPanel parametersP = new JPanel();
            parametersP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Lineage-specific rates",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
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
//                layout.setConstraints(all_fixed_lengthCB, gc);
//                parametersP.add(all_fixed_lengthCB);
                layout.setConstraints(lineage_lengthP, gc);
                parametersP.add(lineage_lengthP);
                gc.gridx++;
                
                gc.gridx++;
//                layout.setConstraints(all_fixed_gainCB, gc);
//                parametersP.add(all_fixed_gainCB);
                layout.setConstraints(lineage_gainP, gc);
                parametersP.add(lineage_gainP);
                gc.gridx++;
                
                gc.gridx++;
//                layout.setConstraints(all_fixed_duplicationCB, gc);
//                parametersP.add(all_fixed_duplicationCB);
                layout.setConstraints(lineage_duplicationP, gc);
                parametersP.add(lineage_duplicationP);
                gc.gridx++;

                gc.gridx++;
//                layout.setConstraints(all_fixed_lossCB, gc);
//                parametersP.add(all_fixed_lossCB);
                layout.setConstraints(lineage_lossP, gc);
                parametersP.add(lineage_lossP);
                gc.gridx++;
                
                for (int node_idx=0; node_idx<nodes.length; node_idx++)
                {
                    NodeWithRates N = nodes[node_idx];
                    if (!N.isRoot())
                    {
                        gc.gridx=0; 
                        gc.gridy++;

                        gc.anchor = GridBagConstraints.WEST;
                        gc.fill = GridBagConstraints.HORIZONTAL;
                        layout.setConstraints(node_nameL[node_idx], gc);
                        parametersP.add(node_nameL[node_idx]);
                        gc.gridx++;
                        gc.fill = GridBagConstraints.NONE;
                        layout.setConstraints(edge_lengthT[node_idx], gc);
                        parametersP.add(edge_lengthT[node_idx]);
                        gc.gridx++;
                        gc.fill = GridBagConstraints.HORIZONTAL;
                        layout.setConstraints(fixed_lengthCB[node_idx], gc);
                        parametersP.add(fixed_lengthCB[node_idx]);
                        gc.gridx++;
                        gc.fill = GridBagConstraints.NONE;
                        layout.setConstraints(gain_rateT[node_idx], gc);
                        parametersP.add(gain_rateT[node_idx]);
                        gc.gridx++;
                        gc.fill = GridBagConstraints.HORIZONTAL;
                        layout.setConstraints(fixed_gainCB[node_idx], gc);
                        parametersP.add(fixed_gainCB[node_idx]);
                        gc.gridx++;
                        gc.fill = GridBagConstraints.NONE;
                        layout.setConstraints(duplication_rateT[node_idx], gc);
                        parametersP.add(duplication_rateT[node_idx]);
                        gc.gridx++;
                        gc.fill = GridBagConstraints.HORIZONTAL;
                        layout.setConstraints(fixed_duplicationCB[node_idx], gc);
                        parametersP.add(fixed_duplicationCB[node_idx]);
                        gc.gridx++;
                        gc.fill = GridBagConstraints.NONE;
                        layout.setConstraints(loss_rateT[node_idx], gc);
                        parametersP.add(loss_rateT[node_idx]);
                        gc.gridx++;
                        gc.fill = GridBagConstraints.HORIZONTAL;
                        layout.setConstraints(fixed_lossCB[node_idx], gc);
                        parametersP.add(fixed_lossCB[node_idx]);
                    }
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

            JPanel main = new JPanel();
            main.setLayout(new BoxLayout(main, BoxLayout.PAGE_AXIS));
            main.add(root_distributionP);
            main.add(rate_variationP);
            main.add(parametersP);
            
            return main;
            
        }


        /**
         * A utility class that enables/disables graphics components
         * for lineage-specific parameters. Lineage-specific parameter
         * components are JFormattedTextFields for values, and
         * accompanying "fixed" check boxes, as well as a master
         * "Fix all" check box. When lineage-specific variation is
         * enabled, all fields can be edited, and all components are enabled.
         * When lineage-specific variation is not allowed,
         * fields are disabled, except for one at index [0], and
         * lineage-specific check boxes are disabled,
         * and take the value of "all fixed."
         *
         * An instance of this class is an "all fixed" check box.
         *
         */
        private class LineageSpecificParametersCB extends CheckSelectAll<JCheckBox>
        {
            LineageSpecificParametersCB(){super(""); init();}

            private void init()
            {
                this.parameter_fields = new JFormattedTextField[1];
                parameter_fields[0] = new JFormattedTextField(); // dummy field
                addItemListener(new ItemListener()
                    {
                        @Override
                        public void itemStateChanged(ItemEvent event)
                        {
                            if (!variation_ok && ItemEvent.DESELECTED == event.getStateChange()) // super's listener will copy SELECTED events
                                copySelectedToAll();
                        }
                });
            }

            void addFields(JFormattedTextField[] parameter_fields)
            {
                this.parameter_fields = parameter_fields;
                listenToField0();
            }

            private void listenToField0()
            {
                parameter_fields[0].addPropertyChangeListener("value", new PropertyChangeListener()
                    {
                        @Override
                        public void propertyChange(PropertyChangeEvent pce)
                        {
//                            System.out.println("#*ROUI.OP.LSP.lF0.pC "+variation_ok+"\t"+pce);
                            if (!variation_ok)
                            {
                                Object x = pce.getNewValue();
                                copyToFields(x);
                            }
                        }
                });
            }

            /**
             * Copies the value to all fields except at index 0
             * @param x value to be copied
             */
            private void copyToFields(Object x)
            {
                for (int field_idx=1; field_idx<parameter_fields.length; ++field_idx)
                    parameter_fields[field_idx].setValue(x);

            }

            private JFormattedTextField[] parameter_fields;
            private boolean variation_ok;

            /**
             * Enables/disables lineage-specific variation
             *
             * @param enabled
             */
            void setVariationEnabled(boolean enabled)
            {
                this.variation_ok = enabled;
                for (int field_idx=1; field_idx<parameter_fields.length; ++field_idx)
                {
                    JFormattedTextField t = parameter_fields[field_idx];
                    t.setEditable(enabled);
                }
                if (!variation_ok)
                    equalizeValues();
                setEnabledAll(variation_ok);
            }

            private void equalizeValues()
            {
                if (parameter_fields.length>0)
                {
                    double sum = 0.0;
                    for (JFormattedTextField ftf:parameter_fields)
                    {
                        double v = ((Number)ftf.getValue()).doubleValue();
                        sum += v;
                    }
                    double avg = sum / parameter_fields.length;
                    parameter_fields[0].setValue(avg);
                }
            }


            /**
             * When enabled, fields and checkboxes are active by the 
             * variation rules. When disabled, 
             * this check box stays at selected, and 
             * none of the fields/checkboxes can change.
             * 
             * @param enabled
             */
            @Override
            public void setEnabled(boolean enabled)
            {
                if (parameter_fields.length>0) parameter_fields[0].setEditable(enabled);
                if (!enabled)
                    setSelected(true);
                setEnabledAll(enabled && variation_ok);
                super.setEnabled(enabled);
            }

            /**
             * Attaches a selection listener to
             * a "uniform parameters" (lineage-specific parameters are the same)
             * check box. When the component is
             * selected, variation is enabled, otherwise it's
             * disabled.
             *
             * @param selection_component "uniform parameters" selection
             */
            void listenToUniformParameters(java.awt.ItemSelectable selection_component)
            {
                selection_component.addItemListener(
                    new ItemListener()
                    {
                        @Override
                        public void itemStateChanged(ItemEvent e)
                        {
                            if (ItemEvent.SELECTED == e.getStateChange())
                                setVariationEnabled(false);
                            else
                            if (ItemEvent.DESELECTED == e.getStateChange())
                                setVariationEnabled(true);
                        }
                    }
                );
            }

            /**
             * Attaches a listener to changes in model architecture.
             * The model architecture may allow for these types of lineage-specific
             * parameters: the selection triggers enabling the associated components (lineage
             * specific check boxes and fields) and this guy ("all fixed check box").
             * Model architecture may not use thse types of parameters:
             * selection triggers a call to {@link #disableParameters() }.
             *
             * @param selection_component model architecture selection component.
             * @param enable_when_selected whether selection of the model architecture component triggers enabling or disabling these associated components
             */
            void listenToModelArchitecture(java.awt.ItemSelectable selection_component, boolean enable_when_selected)
            {
                if (enable_when_selected)
                {
                    selection_component.addItemListener(
                        new ItemListener()
                        {
                            @Override
                            public void itemStateChanged(ItemEvent e)
                            {
                                if (ItemEvent.SELECTED == e.getStateChange())
                                    setEnabled(true);
                            }
                        }
                    );
                } else
                {
                    selection_component.addItemListener(
                        new ItemListener()
                        {
                            @Override
                            public void itemStateChanged(ItemEvent e)
                            {
                                if (ItemEvent.SELECTED == e.getStateChange())
                                    disableParameters();
                            }
                        }
                    );
                }
            }

            /**
             * Disables parameters of this kind. Text field entries are set to 0.
             */
            void disableParameters()
            {
                setEnabled(false);
                if (parameter_fields.length>0)
                {
                    parameter_fields[0].setValue(0.0);
                    if (variation_ok) copyToFields(0.0);
                }
            }
        }


        private class CategoryParametersCB extends JCheckBox
        {
            /**
             * Instantiation.
             *
             * @param valueTF field for entering and displaying parameter value
             * @param default_value value for parameter when disabled
             */
            CategoryParametersCB(JFormattedTextField valueTF, Object default_value)
            {
                super("Fixed");
                this.valueTF = valueTF;
                this.default_value = default_value;
            }

            private JFormattedTextField valueTF;
            private Object default_value;

            void listen(JFormattedTextField num_categoriesTF)
            {
                num_categoriesTF.addPropertyChangeListener("value",
                    new PropertyChangeListener()
                    {
                        @Override
                        public void propertyChange(PropertyChangeEvent update)
                        {
                            int num_categories = ((Number)update.getNewValue()).intValue();
                            if (num_categories==1)
                            { // no variation
                                disableCategories(false);
                            } else
                            {
                                enableCategories(false);
                            }
                        }
                });
            }

            void enableCategories(boolean set_default_value)
            {
                setSelected(false);
                setEnabled(true);
                valueTF.setEditable(true);
                valueTF.setEnabled(true);
                if (set_default_value)
                   valueTF.setValue(default_value);
            }

            void disableCategories(boolean set_default_value)
            {
                setSelected(true);
                setEnabled(false);
                valueTF.setEditable(false);
                valueTF.setEnabled(false);
                if (set_default_value)
                   valueTF.setValue(default_value);
            }

            void listen(JCheckBox variation_allowed)
            {
                variation_allowed.addItemListener(new ItemListener()
                    {
                        @Override
                        public void itemStateChanged(ItemEvent e)
                        {
                            if (ItemEvent.SELECTED == e.getStateChange())
                                enableCategories(false);
                            else if (ItemEvent.DESELECTED == e.getStateChange())
                                disableCategories(true);
                        }
                });
            }
        }


        private class ModelTypeRB extends JRadioButton
        {
            ModelTypeRB(String text)
            {
                super(text);
                init();
            }
            
            private Map<JCheckBox,Boolean> dependent_boxes;
            private Map<JCheckBox,Boolean> dependent_boxes_set;
            private Map<JFormattedTextField,Boolean> dependent_fields;
            
            private void init()
            {
                dependent_boxes = new HashMap<JCheckBox, Boolean>();
                dependent_boxes_set = new HashMap<JCheckBox,Boolean>();
                dependent_fields = new HashMap<JFormattedTextField, Boolean>();
                this.addItemListener(new ItemListener()
                    {
                        @Override
                        public void itemStateChanged(ItemEvent e)
                        {
                            if (ItemEvent.SELECTED == e.getStateChange())
                            {
                                for (JCheckBox cb: dependent_boxes.keySet())
                                {
                                    boolean gets_enabled = dependent_boxes.get(cb);
                                    if (!gets_enabled)
                                        cb.setSelected(dependent_boxes_set.get(cb));
                                    cb.setEnabled(gets_enabled);
                                }
                                for (JFormattedTextField f: dependent_fields.keySet())
                                {
                                    boolean gets_enabled = dependent_fields.get(f);
//                                    f.setEnabled(true);
                                    f.setEditable(gets_enabled);
                                    if (!gets_enabled) f.setValue(1);
//                                    f.setEnabled(gets_enabled);
                                }
                            }
                        }
                });
            }
            
            
            void linkCheckBox(JCheckBox cb, boolean gets_enabled, boolean set_to)
            {
                dependent_boxes.put(cb, gets_enabled);
                dependent_boxes_set.put(cb, set_to);
            }

            /**
             * Linked fields are set to 1 when disabled.
             *
             * @param f text field
             * @param gets_enabled whether selection should enable it
             */
            void linkField(JFormattedTextField f, boolean gets_enabled)
            {
                dependent_fields.put(f, gets_enabled);
            }
        }
        
        
//        private void addModelParameterListeners()
//        {
//            lineage_lossP.listenToUniformParameters(homogeneous_lossCB);
//
//            lineage_duplicationP.listenToUniformParameters(uniform_duplicationCB);
//            lineage_duplicationP.listenToModelArchitecture(gldB, true);
//            lineage_duplicationP.listenToModelArchitecture(plB, false);
//            lineage_duplicationP.listenToModelArchitecture(glB, false);
//            lineage_duplicationP.listenToModelArchitecture(dlB, true);
//
//            lineage_gainP.listenToUniformParameters(uniform_gainCB);
//            lineage_gainP.listenToModelArchitecture(gldB, true);
//            lineage_gainP.listenToModelArchitecture(plB, false);
//            lineage_gainP.listenToModelArchitecture(glB, true);
//            lineage_gainP.listenToModelArchitecture(dlB, false);
//        }


//        private void addUniformRateListeners()
//        {
//            homogeneous_lossCB.addItemListener(lineage_lossP);
//            homogeneous_lossCB.addItemListener(this);
//            uniform_duplicationCB.addItemListener(lineage_duplicationP);
//            uniform_gainCB.addItemListener(lineage_gainP);
//        }
//            uniform_gainCB.addItemListener(this);
//            uniform_duplicationCB.addItemListener(this);
//            homogeneous_lossCB.addItemListener(this);
//
//            gain_rateT[0].addPropertyChangeListener("value", this);
//            duplication_rateT[0].addPropertyChangeListener("value", this);
//            loss_rateT[0].addPropertyChangeListener("value",this);
//
//            fixed_gainCB[0].addItemListener(this);
//            fixed_duplicationCB[0].addItemListener(this);
////            fixed_lossCB[0].addItemListener(this);
////            fixed_lengthCB[0].addItemListener(this);
//        }
        
        //private JPanel createRateVariationPanel()
        //{
        //    return new JPanel();
        //}
        
        //private JPanel createAdvancedOptionsPanel()
        //{
        //    return new JPanel();
        //}
        
        
        private void copyFieldValuesIntoModel()
        {
            boolean is_alternative_model = !null_modelRB.isSelected();
            RateVariation selected_model = (is_alternative_model?alternative_model:null_model);

            TreeWithRates selected_tree = selected_model.getMainTree();
            
            // update by modified values in GUI
            
            // root distribution
            DiscreteDistribution distr = null;
            if (root_poissonB.isSelected())
            {
                double lambda = ((Number)root_aT.getValue()).doubleValue();
                distr = new Poisson(lambda);
            } else if (root_negbinB.isSelected())
            {
                double r = ((Number)root_aT.getValue()).doubleValue();
                double q = ((Number)root_bT.getValue()).doubleValue();
                distr = new NegativeBinomial(r,q);
            } else if (root_pointB.isSelected())
            {
                double p0 = ((Number)root_aT.getValue()).doubleValue();
                distr = new PointDistribution(p0);
            }
            selected_model.setRootPrior(distr);
            
            // rate variation across families
            double forbidden_gain = ((Number)forbidden_gainT.getValue()).doubleValue();
            double forbidden_duplication = ((Number)forbidden_duplicationT.getValue()).doubleValue();
            selected_model.setInvariantFractions(forbidden_duplication, selected_model.getLossForbiddenProportion(), forbidden_gain);
            int cat_length = ((Number)gamma_lengthT.getValue()).intValue();
            int cat_gain = ((Number)gamma_gainT.getValue()).intValue();
            int cat_loss = ((Number)gamma_lossT.getValue()).intValue();
            int cat_duplication = ((Number)gamma_duplicationT.getValue()).intValue();
            selected_model.setNumberOfDiscreteCategories(cat_duplication, cat_loss, cat_gain, cat_length);
            double alpha_gain = ((Number)alpha_gainT.getValue()).doubleValue();
            double alpha_length = ((Number)alpha_lengthT.getValue()).doubleValue();
            double alpha_loss = ((Number)alpha_lossT.getValue()).doubleValue();
            double alpha_duplication = ((Number)alpha_duplicationT.getValue()).doubleValue();
            selected_model.setEdgeLengthAlpha(alpha_length);
            selected_model.setTransferRateAlpha(alpha_gain);
            selected_model.setLossRateAlpha(alpha_loss);
            selected_model.setDuplicationRateAlpha(alpha_duplication);
            
            // lineage-specific rates
            NodeWithRates[] nodes = selected_tree.getDFT();
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                NodeWithRates N = nodes[node_idx];
                if (!N.isRoot())
                {
                    double len = ((Number)edge_lengthT[node_idx].getValue()).doubleValue();
                    double ls = ((Number)loss_rateT[node_idx].getValue()).doubleValue();
                    double gn = ((Number)gain_rateT[node_idx].getValue()).doubleValue();
                    double dp = ((Number)duplication_rateT[node_idx].getValue()).doubleValue();
                    selected_model.setEdgeLength(node_idx, len);
                    selected_model.setNodeLossRate(node_idx, ls);
                    selected_model.setNodeTransferRate(node_idx, gn);
                    selected_model.setNodeDuplicationRate(node_idx, dp);
                }
            }
            
            if (homogeneous_lossCB.isSelected())
            {
                TreeNode original_root = dealer.getActiveWorkSpaceCount().getOriginalNode(main_tree.getNumNodes()-1);
                selected_model = new HomogeneousRateVariation(original_root,selected_model);
            }
            
            // set up the optimization
            StableComputation SC = new StableComputation();
            SC.initWithoutAllocations(data_file.getData(),selected_model);
            ML optimization = new ML(SC);
            // convergence
            int num_steps = ((Number)roundT.getValue()).intValue();
            double eps = ((Number)epsT.getValue()).doubleValue();
            optimization.setOptimizationParameters(num_steps, eps);
            
            boolean uniform_duplication = uniform_duplicationCB.isSelected();
            boolean uniform_gain = uniform_gainCB.isSelected();
            optimization.setUniformEdgeParameters(uniform_duplication, uniform_gain, false);
            
            boolean forbidden_gain_ok = !fixed_forbidden_gainCB.isSelected();
            boolean forbidden_duplication_ok = !fixed_forbidden_duplicationCB.isSelected();
            optimization.setForbiddenCategories(forbidden_gain_ok, forbidden_duplication_ok);
            
            //System.out.println("#*ROUI.OP.cFVIM forbidden dup"+selected_model.getDuplicationForbiddenProportion()+"/"+forbidden_duplication_ok
            //        +"\tgain "+selected_model.getTransferForbiddenProportion()+"/"+forbidden_gain_ok);
            
            boolean fixed_alpha_length = (cat_length==1 || fixed_alpha_lengthCB.isSelected());
            boolean fixed_alpha_gain = (cat_gain == 1 || fixed_alpha_gainCB.isSelected());
            boolean fixed_alpha_duplication = (cat_duplication == 1 || fixed_alpha_duplicationCB.isSelected());
            boolean fixed_alpha_loss = (cat_loss == 1 || fixed_alpha_lossCB.isSelected());
            optimization.setOptimizationFixedVariation(fixed_alpha_length, fixed_alpha_duplication, fixed_alpha_loss, fixed_alpha_gain);
            
            boolean fixed_root = fixed_root_CB.isSelected();
            boolean fixed_length[] = null;
            boolean fixed_duplication[] = null;
            boolean fixed_gain[] = null;
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                NodeWithRates N = nodes[node_idx];
                if (!N.isRoot())
                {
                    if (fixed_lengthCB[node_idx].isSelected())
                    {
                        if (fixed_length == null)
                            fixed_length = new boolean[nodes.length-1]; // initialized with all false
                        fixed_length[node_idx]=true;
                    }
                    if (fixed_gainCB[node_idx].isSelected())
                    {
                        if (fixed_gain == null)
                            fixed_gain = new boolean[nodes.length-1];
                        fixed_gain[node_idx]=true;
                    }
                    if (fixed_duplicationCB[node_idx].isSelected())
                    {
                        if (fixed_duplication == null)
                            fixed_duplication = new boolean[nodes.length-1];
                        fixed_duplication[node_idx]=true;
                    }
                }
            }

            optimization.setOptimizationFixedParameters(fixed_root, fixed_duplication, fixed_gain, fixed_length);
            
            // name the rate variation model
            File optimized_model_file = null;
            if (is_alternative_model)
                optimized_model_file = new File((File)null, "opt:"+rate_file.getFile().getName());
            else
                optimized_model_file = new File((File)null, "opt@"+data_file.getFile().getName());
            DataFile<RateVariation> optimized_model = new DataFile<RateVariation>(selected_model,optimized_model_file);
            
            optimized_model_display = new ModelDisplay(optimized_model,optimization,is_alternative_model);
        }
        
        /**
         * Capturing button clicks.
         * @param E the action event generated by pressing one of our buttons
         */
        @Override
        public void actionPerformed(ActionEvent E)
        {
            Object src = E.getSource();
            if (src == cancel_button)
            {
                dispose();
            } else if (src == start_button)
            {
                // set up the rate model display for the optimized model
                copyFieldValuesIntoModel();
                dispose();
            } else 
            {
                String command = E.getActionCommand();
                if ("alternative-model".equals(command))
                {
                    initModelParameters(true);
                } else if ("null-model".equals(command))
                {
                    initModelParameters(false);
                } else if ("root-Poisson".equals(command))
                {
                    setRootDistributionComponentStates(Poisson.class);
                } else if ("root-negbin".equals(command))
                {
                    setRootDistributionComponentStates(NegativeBinomial.class);
                } else if ("root-point".equals(command))
                {
                    setRootDistributionComponentStates(PointDistribution.class);
//                } else if ("GLD".equals(command))
//                {
//                    setModelFeatureButtons(true, true);
//                } else if ("DL".equals(command))
//                {
//                    setModelFeatureButtons(true, false);
//                } else if ("DL=".equals(command))
//                {
//                    setModelFeatureButtons(true, false);
//                } else if ("GL".equals(command))
//                {
//                    setModelFeatureButtons(false, true);
//                } else if ("PL".equals(command))
//                {
//                    setModelFeatureButtons(false, false);
                } else 
                {
                    //System.out.println("#*ROUI.OP.aP unhandled ActionEvent '"+command+"'\t// "+E);
                }
            }
        }    
        
        /**
         * We are listening to changes in text fields' values
         * 
         * @param E the associated property change event
         */
        @Override
        public void propertyChange(PropertyChangeEvent E)
        {
//            Object source = E.getSource();
//            if (source == gamma_lengthT || source == gamma_gainT || source == gamma_lossT || source == gamma_duplicationT)
//            {
//                if (source == gamma_lengthT)
//                    alpha_lengthCB.setSelected(false);
//                else if (source == gamma_gainT)
//                    alpha_gainCB.setSelected(false);
//                else if (source == gamma_lossT)
//                    alpha_lossCB.setSelected(false);
//                else if (source == gamma_duplicationT)
//                    alpha_duplicationCB.setSelected(false);
//                setGammaCategoriesComponentStates();
//
////            }
////            else if (source == gain_rateT[0] || source == duplication_rateT[0] || source==loss_rateT[0])
////            {
////                copyUniformRateComponents();
//            } else
//            {
//                //System.out.println("#*ROUI.OP.pC unhandled property change "+E);
//            }
//            //Component comp = (Component)source;
//            //if (comp.isFocusOwner() )
//            //    comp.transferFocusUpCycle();
        }
        
        @Override
        public void itemStateChanged(ItemEvent E)
        {
            Object source = E.getSource();
//            if (source == forbidden_gainCB || source==forbidden_duplicationCB)
//            {
//                if (E.getStateChange() == ItemEvent.DESELECTED)
//                {
//                    if (source == forbidden_gainCB)
//                    {
//                        forbidden_gainT.setValue(new Double (0.0));
//                        fixed_forbidden_gainCB.setSelected(true);
//                    } else
//                    {
//                        forbidden_duplicationT.setValue(new Double(0.0));
//                        fixed_forbidden_duplicationCB.setSelected(true);
//                    }
//                } else if (E.getStateChange() == ItemEvent.SELECTED)
//                {
//                    if (source == forbidden_gainCB)
//                        fixed_forbidden_gainCB.setSelected(false);
//                    else
//                        fixed_forbidden_duplicationCB.setSelected(false);
//                }
//                setForbiddenRateComponentStates();
////            } else if (source == uniform_gainCB)
////            {
////                if (ItemEvent.SELECTED == E.getStateChange())
////                    lineage_gainP.setVariationEnabled(false);
////                else if (ItemEvent.DESELECTED == E.getStateChange())
////                    lineage_gainP.setVariationEnabled(true);
////
//////                adjustEditableGainGadgets();
//////                if (E.getStateChange() == ItemEvent.SELECTED)
//////                    equalizeGainRates();
////            } else if (source == uniform_duplicationCB)
////            {
////                if (ItemEvent.SELECTED == E.getStateChange())
////                    lineage_duplicationP.setVariationEnabled(false);
////                else if (ItemEvent.DESELECTED == E.getStateChange())
////                    lineage_duplicationP.setVariationEnabled(true);
//////                adjustEditableDuplicationGadgets();
//////                if (E.getStateChange() == ItemEvent.SELECTED)
//////                    equalizeDuplicationRates();
//////            } else if (source == fixed_duplicationCB[0] || source == fixed_gainCB[0] || source == fixed_lossCB[0])
//////            {
//////                copyUniformRateComponents();
//            } else
                if (source == homogeneous_lossCB)
            {
                if (ItemEvent.SELECTED==E.getStateChange())
                {
                    copyEdgeLengths();
                    lineage_lengthP.setSelected(true);
                    lineage_lengthP.setEnabledAll(false);
                } else if (ItemEvent.DESELECTED==E.getStateChange())
                {
                    copyEdgeLengths();
                    lineage_lengthP.setVariationEnabled(true);
                    lineage_lengthP.setEnabledAll(true);
                }
            }
        }

        private void copyEdgeLengths()
        {
            if (homogeneous_lossCB.isSelected())
            {
                WorkSpaceCount ws = dealer.getActiveWorkSpaceCount();
                double original_total_length = 0.0;
                double total_length = 0.0;
                double total_loss = 0.0;
                for (int edge_idx=0; edge_idx<main_tree.getNumEdges(); ++edge_idx)
                {
                    double len = ((Number)edge_lengthT[edge_idx].getValue()).doubleValue();
                    total_length += len;
                    double original_length = ws.getOriginalEdgeLength(edge_idx);
                    edge_lengthT[edge_idx].setInputVerifier(null);
                    edge_lengthT[edge_idx].setValue(original_length);
                    original_total_length += original_length;
                    double ls = len*((Number)loss_rateT[edge_idx].getValue()).doubleValue();
                    total_loss += ls;
                }
                double common_loss_rate = total_loss / original_total_length;
                for (int edge_idx=0; edge_idx<main_tree.getNumEdges(); ++edge_idx)
                {
                    // now loss rate is allowed to vary
                    loss_rateT[edge_idx].setInputVerifier(null);
                    loss_rateT[edge_idx].setValue(common_loss_rate);
//                    fixed_lossCB[edge_idx].setEnabled(edge_idx==0);
//                    loss_rateT[edge_idx].setEditable(edge_idx==0);
//                    // but edge length must be fixed
//                    edge_lengthT[edge_idx].setEditable(false);
                }
            } else // heterogeneous model
            {
                for (int edge_idx=0; edge_idx<main_tree.getNumEdges(); ++edge_idx)
                {
                    edge_lengthT[edge_idx].setInputVerifier(new InputVerifiers.RangeVerifier(ML.MIN_EDGE_LENGTH, ML.MAX_EDGE_LENGTH));
                    edge_lengthT[edge_idx].setValue(main_tree.getNode(edge_idx).getLength());
                    loss_rateT[edge_idx].setInputVerifier(new InputVerifiers.RangeVerifier(ML.MIN_LOSS_RATE, ML.MAX_LOSS_RATE));
                    loss_rateT[edge_idx].setValue(main_tree.getNode(edge_idx).getLossRate());
                    fixed_lossCB[edge_idx].setEnabled(false);
                    loss_rateT[edge_idx].setEditable(false);
                }
            }
        }
    }
    
    
    /**
     * 
     */
    public class ModelDisplay extends RateModelDisplay
    {
        private ModelDisplay(DataFile<RateVariation>rate_data, ML optimization, boolean is_derived_model)
        {
            super(rate_data);
            this.optimization = optimization;
            this.is_derived_model = is_derived_model;
        }
        
        private ML optimization;
        private boolean is_derived_model;
        
        /**
         * Whether the optimization starts with the default null model, or 
         * the model selected in the Rates browser. 
         * 
         * @return whether this is a derived model
         */
        public boolean isDerivedModel()
        {
            return is_derived_model;
        }
        
        public boolean isOptimizationDone()
        {
            return optimization_is_done;
        }
        
        @Override
        protected void createBottomBarLeft(Box bb)
        {
            super.createBottomBarLeft(bb);
            initProgressTrackingComponents(bb);
        }
        
        /**
         * Returns the underlying rate variation model
         * when the optimization is done. 
         * 
         * @return null while the optimization is running, or the final optimized model afterwards
         */
        @Override
        public DataFile<RateVariation> getRateModel()
        {
            if (optimization_is_done)
                return super.getRateModel();   
            else
                return null;
        }
        

        private void snapshot()
        {
            int round = optimization.getOptimizationRound();
            DataFile<RateVariation> current_file = super.getRateModel();
            RateVariation current_rates =  current_file.getData();
            TreeWithRates copy_tree = new TreeWithRates(NodeWithRates.copyTree(current_rates.getMainTree().getRoot()));
            RateVariation copy_rates= current_rates.sameModelForDifferentTree(copy_tree);

            File snapshot_file = new File((File)null,current_file.getFile().getName()+"."+Integer.toString(round));

            WorkSpaceCount ws = (WorkSpaceCount) WorkSpaceCount.getWorkSpace(this);
            ws.addRates(copy_rates, snapshot_file, false);
        }
        
        private void initProgressTrackingComponents(Box bb)
        {
            bb.add(Box.createHorizontalGlue());
            Font tp_font_rm = new Font("Serif", Font.PLAIN, LookAndFeel.TREE_PANEL_FONT_SIZE);

            optimization_progress = new JProgressBar();
            optimization_progress.setIndeterminate(false);
            optimization_progress.setStringPainted(true);
            optimization_progress.setString("");
            //optimization_progress.setPreferredSize(new Dimension(200,30));
            //optimization_progress.setMinimumSize(optimization_progress.getPreferredSize());
            optimization_progress.setFont(tp_font_rm);
            //selected_stats_progress.setEnabled(false);
            optimization_progress.setBorderPainted(false);

            bb.add(optimization_progress);
            
            optimization_stage = new JLabel("");
            optimization_stage.setFont(tp_font_rm);
            optimization_stage.setPreferredSize(new Dimension(400,30));
            optimization_stage.setMinimumSize(optimization_stage.getPreferredSize());
            optimization_stage.setToolTipText("Current optimization stage and log-likelihood (LL)");
            bb.add(optimization_stage);
            
            NumberFormat deltaF = new DecimalFormat("0.######");
            deltaF.setParseIntegerOnly(false);
            optimization_delta = new JFormattedTextField(deltaF);
            optimization_delta.setColumns(9);
            optimization_delta.setMaximumSize(optimization_delta.getPreferredSize());
            optimization_delta.setEditable(false);
            //optimization_delta.setBackground(bb.getBackground());
            optimization_delta.setFont(tp_font_rm);
            
            optimization_delta.setToolTipText("Decrease of the log-likelihood in the previous round");
            
            JLabel delta_label = new JLabel("\u0394LL:");
            delta_label.setFont(tp_font_rm);
            bb.add(delta_label);
            bb.add(optimization_delta);
            
            snapshot_button = new JButton("Snapshot");
            snapshot_button.setToolTipText("Takes a snapshot of the current model");
            bb.add(snapshot_button);
            
            cancel_button = new JButton("Stop");
            cancel_button.setToolTipText("Stops the optimization");
            bb.add(cancel_button);
            
            
            bb.add(Box.createHorizontalGlue());
        }
        
        private void updateDisplay(boolean recompute_layout)
        {
            if (recompute_layout)
                tree_panel.recomputeTreeLayout();
            else
                tree_zoom.setValidBoundingBoxes(false);
            repaint();
        }
        
    }



    private static void testButtonEvents()
    {
        JFrame appF = new JFrame("test");
        appF.setSize(600,400);
        appF.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        JPanel appP = new JPanel();
        BoxLayout appL = new BoxLayout(appP, BoxLayout.PAGE_AXIS);
        appP.setLayout(appL);

        class BM extends javax.swing.JToggleButton.ToggleButtonModel
        {
            BM(){super();}
            @Override
            public void setArmed(boolean b)
            {
                System.out.println("BM/"+getActionCommand()+"\tarmed "+b);
                super.setArmed(b);
            }
            @Override
            public void setPressed(boolean b)
            {
                System.out.println("BM/"+getActionCommand()+"\tpressed "+b);
                super.setPressed(b);
            }
            @Override
            public void setSelected(boolean b)
            {
                System.out.println("BM/"+getActionCommand()+"\tselected "+b);
                super.setSelected(b);
            }
        }

        Box tfB = new Box(BoxLayout.LINE_AXIS);

        class TF extends JFormattedTextField
        {
            TF(java.text.Format f){super(f);}
            @Override
            public void setValue(Object o)
            {
                System.out.println("TF\tsetValue\t"+o);
                super.setValue(o);
            }
        }

        NumberFormat egyF = new DecimalFormat("0.############");
        egyF.setParseIntegerOnly(false);

        final JFormattedTextField egyTF = new TF(egyF);
        egyTF.setColumns(10);
        egyTF.setInputVerifier(new InputVerifiers.RangeVerifier(0, 1000)
        {
            @Override
            public boolean shouldYieldFocus(javax.swing.JComponent c)
            {
                JFormattedTextField f = (JFormattedTextField) c;
                boolean b = super.shouldYieldFocus(c);
                System.out.println("IV\tYield "+b+"\t"+c);
                return b;
            }

            @Override
            public boolean verify(javax.swing.JComponent c)
            {
                boolean b = super.verify(c);
                System.out.println("IV\tVerify "+b+"\t"+c);
                return b;
            }
        });
        egyTF.setFocusLostBehavior(JFormattedTextField.COMMIT_OR_REVERT);
        egyTF.addPropertyChangeListener("value", new PropertyChangeListener()
            {
                @Override
                public void propertyChange(PropertyChangeEvent e)
                {
                    JFormattedTextField t = (JFormattedTextField)e.getSource();
                    System.out.println("TF\tProperty\tval "+t.getValue()+"\t"+e);
                }

            });
        egyTF.setValue(100);


        tfB.add(egyTF);
        final JCheckBox enabledCB = new JCheckBox("enable", true);
        final JCheckBox editableCB = new JCheckBox("editable", true);

        enabledCB.addItemListener(new ItemListener()
        {
            @Override
            public void itemStateChanged(ItemEvent e)
            {
                egyTF.setEnabled(enabledCB.isSelected());
            }
        });
        editableCB.addItemListener(new ItemListener()
        {
            @Override
            public void itemStateChanged(ItemEvent e)
            {
                egyTF.setEditable(editableCB.isSelected());
            }
        });
        tfB.add(enabledCB);
        tfB.add(editableCB);

        appP.add(tfB);

        ActionListener AL = new ActionListener()
            {
                @Override
                public void actionPerformed(ActionEvent e)
                {
                    javax.swing.AbstractButton cb = (javax.swing.AbstractButton) e.getSource();
                    System.out.println(cb.getText()+"\tAction\t"+e+"\n");
                }
            };

        Box rbB = new Box(BoxLayout.LINE_AXIS);

        class RB extends JRadioButton
        {
            RB(String txt){super(txt);}
            RB(String txt, boolean b){super(txt, b);}
            @Override
            public void setSelected(boolean b)
            {
                System.out.println(getText()+"\tsetSelected\told "+isSelected()+"\tnew "+b);
                super.setSelected(b);
            }
        }

        final JRadioButton egyRB = new RB("R1", true);
        egyRB.setActionCommand(egyRB.getText());
        final JRadioButton ketRB = new RB("R2");
        ketRB.setActionCommand(ketRB.getText());

        ItemListener IL = new ItemListener()
            {
                @Override
                public void itemStateChanged(ItemEvent e)
                {
                    javax.swing.AbstractButton b = (javax.swing.AbstractButton) e.getSource();
                    System.out.println(b.getText()+"\tItem\tsel "+b.isSelected()+"\t"+e);
                    if ("R2".equals(b.getText()))
                        egyTF.setValue(ketRB.isSelected()?-1:1);
                }
            };
        egyRB.setModel(new BM());
        egyRB.getModel().setActionCommand("m"+egyRB.getText());
        ketRB.setModel(new BM());
        ketRB.getModel().setActionCommand("m"+ketRB.getText());

        ButtonGroup BG = new ButtonGroup();
        BG.add(egyRB);
        BG.add(ketRB);

        egyRB.addActionListener(AL);
        ketRB.addActionListener(AL);

        egyRB.addItemListener(IL);
        ketRB.addItemListener(IL);

        rbB.add(egyRB);
        rbB.add(ketRB);
        appP.add(rbB);

        Box cbB = new Box(BoxLayout.LINE_AXIS);
        class CB extends JCheckBox
        {
            CB(String text){super(text);}
            CB(String text, boolean b){super(text,b);}
            @Override
            public void setSelected(boolean b)
            {
                System.out.println(getText()+"\tsetSelected\told "+isSelected()+"\tnew "+b);
                super.setSelected(b);
            }
        }


        final JCheckBox egyCB = new CB("C1", true);
        egyCB.setActionCommand(egyCB.getText());
        final JCheckBox ketCB = new CB("C2");
        ketCB.setActionCommand(ketCB.getText());

        egyCB.setModel(new BM());
        egyCB.getModel().setActionCommand("m"+egyCB.getText());
        ketCB.setModel(new BM());
        ketCB.getModel().setActionCommand("m"+ketCB.getText());

        egyCB.addActionListener(AL);
        ketCB.addActionListener(AL);

        ItemListener cbIL = new ItemListener()
            {
                @Override
                public void itemStateChanged(ItemEvent e)
                {
                    javax.swing.AbstractButton cb = (javax.swing.AbstractButton) e.getSource();
                    System.out.println(cb.getText()+"\tItem\tsel "+cb.isSelected()+"\t"+e);
                    if ("C1".equals(cb.getText()))
                    {
                        ketCB.setSelected(cb.isSelected());
                        egyRB.setSelected(cb.isSelected());
                    }
                }
            };
        egyCB.addItemListener(cbIL);
        ketCB.addItemListener(cbIL);

//        egyCB.setSelected(false);
//        egyCB.setSelected(true);
//        ketCB.setSelected(true);
//        ketCB.setSelected(true);
//        ketCB.setSelected(false);

        cbB.add(egyCB);
        cbB.add(ketCB);

        ketCB.setEnabled(false);
        appP.add(cbB);

        appF.getContentPane().add(appP);

        appF.setVisible(true);

        System.out.println("# -------- appF launched.");
    }

    public static void main(String[] args)
    {
        if (args.length>0)
        {
            if ("gombok".equals(args[0]))
            {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {testButtonEvents();} });
            }
        }
    }
}
