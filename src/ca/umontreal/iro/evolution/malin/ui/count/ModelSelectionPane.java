package ca.umontreal.iro.evolution.malin.ui.count;

import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import java.io.File;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;

import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeEvent;

import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.text.Format;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.InputVerifier;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;

import ca.umontreal.iro.matek.DiscreteDistribution;
import ca.umontreal.iro.matek.NegativeBinomial;
import ca.umontreal.iro.matek.PointDistribution;
import ca.umontreal.iro.matek.Poisson;
import ca.umontreal.iro.matek.ShiftedGeometric;

import ca.umontreal.iro.evolution.TreeNode;

import ca.umontreal.iro.evolution.malin.DataFile;

import ca.umontreal.iro.evolution.malin.ui.CheckSelectAll;
import ca.umontreal.iro.evolution.malin.ui.InputVerifiers;

import ca.umontreal.iro.evolution.genecontent.ApproximateComputation;
import ca.umontreal.iro.evolution.genecontent.HomogeneousRateVariation;
import ca.umontreal.iro.evolution.genecontent.ML;
import ca.umontreal.iro.evolution.genecontent.NodeWithRates;
import ca.umontreal.iro.evolution.genecontent.OccurrenceTable;
import ca.umontreal.iro.evolution.genecontent.RateVariation;
import ca.umontreal.iro.evolution.genecontent.StableComputation;
import ca.umontreal.iro.evolution.genecontent.TreeWithRates;


/**
 * Class for selecting model parameters.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class ModelSelectionPane
{
    public ModelSelectionPane(DealerCount dealer)
    {
        this.dealer = dealer;
        initDataStructures();
    }

    private DealerCount dealer;
    private ArrayList<RateVariation> models;
    private ArrayList<JRadioButton> modelRB;
    private int selected_model;

    private ML optimization;

    private void initDataStructures()
    {
        this.models = new ArrayList<RateVariation>();
        this.modelRB = new ArrayList<JRadioButton>();
        this.selected_model = -1;
    }

    public void addModel(RateVariation model, String model_name)
    {
        models.add(model);
        JRadioButton rb = new JRadioButton(model_name); // initally unselected
        rb.setToolTipText(model.getBriefModelDescription());
        modelRB.add(rb);
    }

    /**
     * Adds a copy of the null model
     *
     * @param main_tree
     */
    public void addNullModel(TreeWithRates main_tree)
    {
        TreeWithRates copy_tree = new TreeWithRates(NodeWithRates.copyTree(main_tree.getRoot()));
        RateVariation null_model = new RateVariation(copy_tree, new Poisson(0.1), 1,1,1,1);
        addModel(null_model,"Default null model");
    }

    /**
     * Adds a copy of the model from a rate variation file
     *
     * @param rate_file
     */
    public void addModel(DataFile<RateVariation> rate_file)
    {
        RateVariation original_model = rate_file.getData();
        TreeWithRates copy_tree = new TreeWithRates(NodeWithRates.copyTree(original_model.getMainTree().getRoot()));
        RateVariation model = original_model.sameModelForDifferentTree(copy_tree);
        addModel(model, rate_file.getFile().getName());
    }

    public RateVariation showModelOptimizationDialog(int num_steps, double epsilon)
    {
        JFrame frame = dealer.getTopFrame();
        ModelDialog parameters = new ModelDialog("Rate optimization", num_steps>0);
        parameters.setOptimizationParameters(num_steps, epsilon);
        Dimension frameD = frame.getSize();

        parameters.pack();
        parameters.setBounds((int)(0.05*frameD.width),(int)(0.05*frameD.height),(int)(0.9*frameD.width),(int)(0.9*frameD.height));

        parameters.setVisible(true);

        // modal dialog

        return getSelectedModel();
    }

    public RateVariation getSelectedModel()
    {
        if (selected_model==-1)
            return null;
        else
            return models.get(selected_model);
    }

    public boolean isDefaultModelSelected()
    {
        return (selected_model==0);
    }

    public ML getOptimization()
    {
        return optimization;
    }

    public DataFile<RateVariation> getOptimizedModelFile()
    {
        RateVariation selected_rates = getSelectedModel();
        if (selected_rates == null)
            return null;
        // name the rate variation model
        File optimized_model_file = null;

        if (selected_model==0)
        {
            DataFile<RateVariation> rate_file = dealer.getActiveWorkSpaceCount().getSelectedRateModel();
            optimized_model_file = new File((File)null, "opt:"+rate_file.getFile().getName());
        }
        else
        {
            DataFile<OccurrenceTable> data_file = dealer.getActiveWorkSpaceCount().getSelectedDataset();
            optimized_model_file = new File((File)null, "opt@"+data_file.getFile().getName());
        }
        DataFile<RateVariation> optimized_model = new DataFile<RateVariation>(selected_rates,optimized_model_file);
        return optimized_model;
    }

    private class ModelDialog extends JDialog
    {
        ModelDialog(String title, boolean want_convergence)
        {
            super(dealer.getTopFrame(), title, true); // modal dialog
            initComponents(want_convergence);
            initListeners();

            modelRB.get(modelRB.size()-1).doClick();
//
//            setSelectedModelIndex(0);
        }


        private JButton start_button;
        private JButton cancel_button;

        /**
         * Radio buttons for model type
         */
        private JRadioButton model_gldRB;
        private JRadioButton model_dlRB;
        private JRadioButton model_glRB;
        private JRadioButton model_plRB;

        /**
         * Combo box for selecting grouping
         */
        private GroupingCombo family_groupingCM;

        /**
         * Checkboxes for lineage-specific variation
         */
        private JCheckBox uniform_duplicationCB;
        private JCheckBox uniform_gainCB;
        private JCheckBox uniform_lossCB;

        /**
         * Radio buttons for distribution at root
        */
        private ArrayList<DistributionSelectionButton> root_distributionB;

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
        private JCheckBox forbidden_gainCB;
        private JCheckBox forbidden_duplicationCB;

        /**
         * Fields for convergence criteria
         */
        private JFormattedTextField epsT;
        private JFormattedTextField roundT;

        /**
         * Approximate likelihood calculations
         */
        private JCheckBox approximateCB;
        private JFormattedTextField approximate_factorT;
        private JFormattedTextField approximate_thresholdT;

        /**
         * Fields for gamma parameters
         */
        private ParameterField alpha_lengthT;
        private ParameterField alpha_gainT;
        private ParameterField alpha_lossT;
        private ParameterField alpha_duplicationT;

        /**
         * Proportion of no-duplication and no-gain families
         */
        private ParameterField forbidden_duplicationT;
        private ParameterField forbidden_gainT;

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
         * Lineage-specific parameters
         */
        private LineageSpecificParameters lineage_lengthP;
        private LineageSpecificParameters lineage_duplicationP;
        private LineageSpecificParameters lineage_gainP;
        private LineageSpecificParameters lineage_lossP;


        private Class getRootPriorClass()
        {
            Class distribution_class = null;
            for (DistributionSelectionButton b: root_distributionB)
                if (b.isSelected())
                {
                    distribution_class = b.getDistributionClass();
                    break;
                }
            return distribution_class;
        }

        private DiscreteDistribution getRootPriorDistribution()
        {
            Class c = getRootPriorClass();
            DiscreteDistribution distr = null;
            if (c==NegativeBinomial.class)
            {
                double r = ((Number)root_aT.getValue()).doubleValue();
                double q = ((Number)root_bT.getValue()).doubleValue();
                distr = new NegativeBinomial(r,q);
            } else if (c==Poisson.class)
            {
                double lambda = ((Number)root_aT.getValue()).doubleValue();
                distr = new Poisson(lambda);
            } else if (c==PointDistribution.class)
            {
                double p0 = ((Number)root_aT.getValue()).doubleValue();
                distr = new PointDistribution(p0);
            } else if (c==ShiftedGeometric.class)
            {
                double p0 = ((Number)root_aT.getValue()).doubleValue();
                double q = ((Number)root_bT.getValue()).doubleValue();
                distr = new ca.umontreal.iro.matek.ShiftedGeometric(p0, q);
            }
            return distr;
        }

        /**
         * Copies the UI field values into the selected model (entry in {@link #models}),
         * and initializes the likelihood optimization ({@link #optmization}).
         */
        private void copyFieldValuesIntoModel()
        {
            selected_model = -1;
            for (int model_idx=0;model_idx<models.size() && selected_model==-1; ++model_idx)
            {
                if (modelRB.get(model_idx).isSelected())
                    selected_model = model_idx;
            }
            if (selected_model >= 0)
            {
                RateVariation model = models.get(selected_model);

                // copy values into model
                boolean has_gain = model_gldRB.isSelected() || model_glRB.isSelected();
                model.setGainAllowed(has_gain);

                boolean has_duplication = model_gldRB.isSelected() || model_dlRB.isSelected();
                model.setDuplicationAllowed(has_duplication);

                // root distribution
                model.setRootPrior(getRootPriorDistribution());

                // rate variation across families
                double forbidden_gain = ((Number)forbidden_gainT.getValue()).doubleValue();
                double forbidden_duplication = ((Number)forbidden_duplicationT.getValue()).doubleValue();
                model.setInvariantFractions(forbidden_duplication, model.getLossForbiddenProportion(), forbidden_gain);

                int cat_length = ((Number)gamma_lengthT.getValue()).intValue();
                int cat_gain = ((Number)gamma_gainT.getValue()).intValue();
                int cat_loss = ((Number)gamma_lossT.getValue()).intValue();
                int cat_duplication = ((Number)gamma_duplicationT.getValue()).intValue();
                model.setNumberOfDiscreteCategories(cat_duplication, cat_loss, cat_gain, cat_length);

                double alpha_length = ((Number)alpha_lengthT.getValue()).doubleValue();
                model.setEdgeLengthAlpha(alpha_length);

                double alpha_gain = ((Number)alpha_gainT.getValue()).doubleValue();
                model.setTransferRateAlpha(alpha_gain);

                double alpha_loss = ((Number)alpha_lossT.getValue()).doubleValue();
                model.setLossRateAlpha(alpha_loss);

                double alpha_duplication = ((Number)alpha_duplicationT.getValue()).doubleValue();
                model.setDuplicationRateAlpha(alpha_duplication);

                // lineage-specific rates
                NodeWithRates[] nodes = model.getMainTree().getDFT();
                for (int node_idx=0; node_idx<nodes.length; node_idx++)
                {
                    NodeWithRates N = nodes[node_idx];
                    if (!N.isRoot())
                    {
                        double len = lineage_lengthP.getField(node_idx).doubleValue();
                        model.setEdgeLength(node_idx, len);

                        double ls = lineage_lossP.getField(node_idx).doubleValue();
                        model.setNodeLossRate(node_idx, ls);

                        double gn = lineage_gainP.getField(node_idx).doubleValue();
                        model.setNodeTransferRate(node_idx, gn);

                        double dp = lineage_duplicationP.getField(node_idx).doubleValue();
                        model.setNodeDuplicationRate(node_idx, dp);
                    }
                }

                if (uniform_lossCB.isSelected())
                {
                    TreeNode original_root = dealer.getActiveWorkSpaceCount().getOriginalNode(nodes.length-1);
                    model = new HomogeneousRateVariation(original_root,model);
                }


                // set up the optimization
                DataFile<OccurrenceTable> data_file = dealer.getActiveWorkSpaceCount().getSelectedDataset();

                // likelihood computing
                if (approximateCB.isSelected())
                {
                    ApproximateComputation AC = new ApproximateComputation();
                    double factor = ((Number)approximate_factorT.getValue()).doubleValue();
                    int threshold = ((Number)approximate_thresholdT.getValue()).intValue();
                    AC.setCutoffFactor(factor);
                    AC.setCutoffThreshold(threshold);
                    AC.initWithoutAllocations(data_file.getData(), model);
                    optimization = new ML(AC);
                } else
                {
                    StableComputation SC = new StableComputation();
                    SC.initWithoutAllocations(data_file.getData(),model);
                    optimization = new ML(SC);
                }

                // convergence
                int num_steps = ((Number)roundT.getValue()).intValue();
                double eps = ((Number)epsT.getValue()).doubleValue();
                optimization.setOptimizationParameters(num_steps, eps);

                // uniform rates
                boolean uniform_duplication = uniform_duplicationCB.isSelected();
                boolean uniform_gain = uniform_gainCB.isSelected();
                optimization.setUniformEdgeParameters(uniform_duplication, uniform_gain, false);

                //
                boolean forbidden_gain_ok = !forbidden_gainT.isFixed();
                boolean forbidden_duplication_ok = !forbidden_duplicationT.isFixed();
                optimization.setForbiddenCategories(forbidden_gain_ok, forbidden_duplication_ok);
//
//                //System.out.println("#*ROUI.OP.cFVIM forbidden dup"+selected_model.getDuplicationForbiddenProportion()+"/"+forbidden_duplication_ok
//                //        +"\tgain "+selected_model.getTransferForbiddenProportion()+"/"+forbidden_gain_ok);
//
                boolean fixed_alpha_length = (cat_length==1 || alpha_lengthT.isFixed());
                boolean fixed_alpha_gain = (cat_gain == 1 || alpha_gainT.isFixed());
                boolean fixed_alpha_duplication = (cat_duplication == 1 || alpha_duplicationT.isFixed());
                boolean fixed_alpha_loss = (cat_loss == 1 || alpha_lossT.isFixed());
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
                        if (lineage_lengthP.getField(node_idx).isFixed())
                        {
                            if (fixed_length == null)
                                fixed_length = new boolean[nodes.length-1]; // initialized with all false
                            fixed_length[node_idx]=true;
                        }
                        if (lineage_gainP.getField(node_idx).isFixed())
                        {
                            if (fixed_gain == null)
                                fixed_gain = new boolean[nodes.length-1];
                            fixed_gain[node_idx]=true;
                        }
                        if (lineage_duplicationP.getField(node_idx).isFixed())
                        {
                            if (fixed_duplication == null)
                                fixed_duplication = new boolean[nodes.length-1];
                            fixed_duplication[node_idx]=true;
                        }
                    }
                }

                optimization.setOptimizationFixedParameters(fixed_root, fixed_duplication, fixed_gain, fixed_length);
            }
        }

        private void setSelectedModelIndex(int idx)
        {
            System.out.println("#*MSP.MD.sSMI "+idx);
            selected_model = idx;
            clickModelType();
            clickRootPrior();
            fillRootPriorParameters();
            selectFamilyVariation();
            fillFamilyVariationParameters();
            selectLineageVariation();
            fillLineageSpecificParameters();
        }

//        private int getSelectedModelIndex()
//        {
//            int selected_model = -1;
//            for (int model_idx=0;model_idx<models.size() && selected_model==-1; ++model_idx)
//            {
//                if (modelRB.get(model_idx).isSelected())
//                    selected_model = model_idx;
//            }
//            return selected_model;
//        }

        /**
         * Initializes the graphics components.
         */
        private void initComponents(boolean want_convergence)
        {
            // two tabs: one for model architecture, and another for model parameters
            JTabbedPane tabs = new JTabbedPane();
            tabs.addTab("Model type", createModelStructurePanel(want_convergence));
            tabs.addTab("Model parameters", new JScrollPane(createModelParametersPanel()));

            setLayout(new BorderLayout());
            add(tabs, BorderLayout.CENTER);
            add(createButtonPanel(),BorderLayout.PAGE_END);
        }


        private void initListeners()
        {
            attachButtonListeners();
            attachInitialModelListeners();
            attachModelStructureListeners();
            attachRootPriorListeners();
            attachFamilyVariationListeners();
            attachLineageVariationListeners();
        }

        private JComponent createModelStructurePanel(boolean want_optimization)
        {
            Box structureB = new Box(BoxLayout.PAGE_AXIS);
            if (models.size()!=0)
                structureB.add(createInitialModelBox());

            structureB.add(createModelTypeBox());
            structureB.add(createGroupingBox());

            structureB.add(createLineageVariationBox());
            structureB.add(createRootPriorBox());
            structureB.add(createFamilyVariationBox());
            if (want_optimization)
                structureB.add(createOptimizationBox());

            return structureB;
        }

        private JComponent createModelParametersPanel()
        {
            Box parametersB = new Box(BoxLayout.PAGE_AXIS);

            parametersB.add(createFamilyVariationParameterBox());
            parametersB.add(createRootPriorParameterBox());
            parametersB.add(createLineageSpecificParameterBox());

            return parametersB;
        }

        /**
         * Panel at the bottom of the dialog with Start and Cancel buttons
         *
         * @return a container with these buttons, and initializes {@link start_button} and {@link cancel_button}
         */
        private JComponent createButtonPanel()
        {
            start_button = new JButton("OK");
            cancel_button = new JButton("Cancel");
            Box button_box = new Box(BoxLayout.LINE_AXIS);
            button_box.add(Box.createHorizontalGlue());
            button_box.add(cancel_button);
            button_box.add(Box.createRigidArea(new Dimension(10,0)));
            button_box.add(start_button);
            button_box.setBorder(BorderFactory.createEmptyBorder(0, 10, 10, 10));
            return button_box;
        }

        void setStartButtonTrext(String txt)
        {
            start_button.setText(txt);
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

//            starting_modelP.setBackground(java.awt.Color.RED);

            ButtonGroup starting_modelG = new ButtonGroup();
            for (JRadioButton rb: modelRB)
            {
                starting_modelG.add(rb);
                starting_modelP.add(rb);
            }
            starting_modelP.add(Box.createHorizontalGlue());
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
//            gldB.addActionListener(this);
            model_dlRB = new JRadioButton("Duplication-loss");
            model_dlRB.setActionCommand("DL");
//            dlB.addActionListener(this);
            //dlEB = new JRadioButton("Duplication-loss with equal rates (Hahn et al.)");
            //dlEB.setActionCommand("DL=");
            //dlEB.addActionListener(this);
            model_glRB = new JRadioButton("Gain-loss");
            model_glRB.setActionCommand("GL");
//            glB.addActionListener(this);
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
            uniform_gainCB = new JCheckBox("Same gain-loss ratio in all lineages");
            uniform_lossCB = new JCheckBox("Same loss rate in all lineages");

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
         * radio buttons {@link #root_pointB}, {@link #root_negbinB}, {@link #root_poissonB} and {@link #root_stationaryB}.
         * All radio buttons are initially unselected.
         *
         * @returna container
         */
        private JComponent createRootPriorBox()
        {
            root_distributionB = new ArrayList<DistributionSelectionButton>();
            root_distributionB.add(new DistributionSelectionButton("Poisson", Poisson.class));
            root_distributionB.add(new DistributionSelectionButton("Negative binomial (P\u00f3lya)", NegativeBinomial.class));
            root_distributionB.add(new DistributionSelectionButton("Bernoulli",PointDistribution.class));


            DistributionSelectionButton root_stationaryB = new DistributionSelectionButton("Stationary",DiscreteDistribution.class);
            root_stationaryB.setEnabled(false);
            root_distributionB.add(root_stationaryB);

            ButtonGroup rootG = new ButtonGroup();
            for (DistributionSelectionButton b: root_distributionB)
                rootG.add(b);

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

                for (DistributionSelectionButton b: root_distributionB)
                {
                    layout.setConstraints(b, gc);
                    rootP.add(b);
                    ++gc.gridx;
                }
                gc.gridheight=gc.gridy+1;
                gc.gridx++;
                gc.weightx=1.0;
                gc.fill=GridBagConstraints.BOTH;
                Component fillerX = Box.createHorizontalGlue();
                layout.setConstraints(fillerX,gc);
                rootP.add(fillerX);
            }
            return rootP;
        }


        /**
         * Panel with fields and check boxes for selecting rate variation across families.
         * Initializes {@link #gamma_duplicationT}, {@link #gamma_gainT}, {@link #gamma_lengthT}
         * and {@link #gamma_lossT}, {@link #forbidden_duplicationCB}, {@link #forbidden_gainCB}.
         * Initally, all check boxes are unselected, and all text field values are unset.
         *
         * Loss variation is disabled for now (not editable, not enabled, not visible).
         *
         * @return a container
         */
        private JComponent createFamilyVariationBox()
        {
            gamma_lengthT = new JFormattedTextField(NumberFormat.getIntegerInstance());
            gamma_lengthT.setInputVerifier(InputVerifiers.getPositiveInstance());
            gamma_lengthT.setColumns(2);
            gamma_lengthT.setToolTipText("Number of discrete Gamma categories for edge length (1 means variation disabled)");

            //gamma_lengthT.setFocusLostBehavior(JFormattedTextField.COMMIT_OR_REVERT);
            gamma_gainT = new JFormattedTextField(NumberFormat.getIntegerInstance());
            gamma_gainT.setInputVerifier(InputVerifiers.getPositiveInstance());
            gamma_gainT.setColumns(2);
            gamma_gainT.setToolTipText("Number of discrete Gamma categories for gain rate (1 means variation disabled)");

            gamma_duplicationT = new JFormattedTextField(NumberFormat.getIntegerInstance());
            gamma_duplicationT.setInputVerifier(InputVerifiers.getPositiveInstance());
            gamma_duplicationT.setColumns(2);
            gamma_duplicationT.setToolTipText("Number of discrete Gamma categories for duplication rate (1 means variation disabled)");

            gamma_lossT = new JFormattedTextField(NumberFormat.getIntegerInstance());
            gamma_lossT.setInputVerifier(InputVerifiers.getPositiveInstance());
            gamma_lossT.setColumns(2);
            gamma_lossT.setToolTipText("Number of discrete Gamma categories for loss rate (1 means variation disabled)");

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


            forbidden_gainCB = new JCheckBox("No-gain category");
            forbidden_duplicationCB = new JCheckBox("No-duplication category");

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
            return family_variationP;
        }


        /**
         * Panel with fields for convergence criteria.
         * Initializes {@link #epsT} and {@link #roundT}.
         * Initially, field values are unset.
         *
         * @return a container
         */
        private JComponent createOptimizationBox()
        {
            JLabel epsL = new JLabel("Convergence threshold on the likelihood");
            JLabel roundL = new JLabel("Maximum number of optimization rounds");
            NumberFormat epsF = new DecimalFormat("0.############");
            epsF.setParseIntegerOnly(false);
            epsT = new JFormattedTextField(epsF);
            //epsT.addPropertyChangeListener("value",this);
            epsT.setInputVerifier(InputVerifiers.getNonnegativeInstance());
            epsT.setColumns(16);

            NumberFormat roundF = NumberFormat.getIntegerInstance();
            roundT = new JFormattedTextField(roundF);
            //roundT.addPropertyChangeListener("value",this);
            roundT.setInputVerifier(InputVerifiers.getNonnegativeInstance());
            roundT.setColumns(10);

            NumberFormat approximate_likelihood_factorF = new DecimalFormat("0.0###########");
            approximate_likelihood_factorF.setParseIntegerOnly(false);
            approximate_factorT = new JFormattedTextField(approximate_likelihood_factorF);
            approximate_factorT.setColumns(16);
            JLabel approximate_factorL = new JLabel("Truncation: factor applied to max. occurrence");
            approximate_factorL.setLabelFor(approximate_factorT);

            NumberFormat approximate_likelihood_thresholdF = NumberFormat.getIntegerInstance();
            approximate_thresholdT = new JFormattedTextField(approximate_likelihood_thresholdF);
            approximate_thresholdT.setInputVerifier(InputVerifiers.getNonnegativeInstance());
            approximate_thresholdT.setColumns(10);
            JLabel approximate_thresholdL = new JLabel("Truncation: minimum cutoff threshold");
            approximate_thresholdL.setLabelFor(approximate_thresholdT);

            approximateCB = new JCheckBox("Truncated likelihood computation");
            approximateCB.addActionListener(new ActionListener()
            {
                @Override
                public void actionPerformed(ActionEvent e)
                {
                    boolean approxOK = approximateCB.isSelected();
                    approximate_thresholdT.setEnabled(approxOK);
                    approximate_thresholdT.setEditable(approxOK);
                    approximate_factorT.setEnabled(approxOK);
                    approximate_factorT.setEditable(approxOK);
                }
            });
            approximateCB.setSelected(true);



            JPanel convergenceP = new JPanel();

            convergenceP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Optimization parameters",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
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
                gc.gridheight=3;
                gc.weightx=1.0;
                gc.fill=GridBagConstraints.BOTH;
                {
                    Component fillerX = Box.createHorizontalGlue();
                    layout.setConstraints(fillerX,gc);
                    convergenceP.add(fillerX);
                }

                ++gc.gridx;
                gc.gridy=0;
                gc.gridwidth = 2;
                gc.gridheight = 1;
                gc.fill = GridBagConstraints.NONE;
                gc.weightx = 0.2;
                layout.setConstraints(approximateCB, gc);
                convergenceP.add(approximateCB);

                gc.gridwidth=1;
                gc.weightx=0.1;
                ++gc.gridy;
                gc.fill = GridBagConstraints.HORIZONTAL;
                layout.setConstraints(approximate_factorL, gc);
                convergenceP.add(approximate_factorL);

                ++gc.gridx;
                gc.weightx=0.0;
                gc.fill=GridBagConstraints.NONE;
                layout.setConstraints(approximate_factorT,gc);
                convergenceP.add(approximate_factorT);

                --gc.gridx;
                ++gc.gridy;
                gc.weightx=0.1;
                gc.fill = GridBagConstraints.HORIZONTAL;
                layout.setConstraints(approximate_thresholdL, gc);
                convergenceP.add(approximate_thresholdL);

                ++gc.gridx;
                gc.weightx=0.0;
                gc.fill=GridBagConstraints.NONE;
                layout.setConstraints(approximate_thresholdT, gc);
                convergenceP.add(approximate_thresholdT);


                ++gc.gridx;
                gc.gridy=0;
                gc.gridheight=3;
                gc.weightx=1.0;
                gc.fill=GridBagConstraints.BOTH;
                {
                    Component fillerX = Box.createHorizontalGlue();
                    layout.setConstraints(fillerX,gc);
                    convergenceP.add(fillerX);
                }
            }
            return convergenceP;
        }



        private JComponent createGroupingBox()
        {
            OccurrenceTable table = dealer.getActiveWorkSpaceCount().getSelectedDataset().getData();

            int num_families = table.getNumFamilies();

            String[] data_props = table.getKnownProperties();

            String[] combo_choices = new String[1+data_props.length];
            String[] combo_descriptions = new String[combo_choices.length];
            int[] num_groups = new int[combo_choices.length];

            for (int prop_idx=-1; prop_idx<data_props.length; ++prop_idx)
            {
                int combo_idx = prop_idx+1;
                num_groups[combo_idx]=0;
                Map<String, Integer> property_counts = new HashMap<String, Integer>();
                if (prop_idx==-1)
                {

                    combo_choices[combo_idx] = "[no grouping]";
                    num_groups[combo_idx] = 1;
                } else
                {
                    combo_choices[combo_idx]=data_props[prop_idx];

                    for (int family_idx=0; family_idx<num_families; ++family_idx)
                    {
                        String p = table.getFamilyProperty(family_idx, prop_idx);
                        if (property_counts.containsKey(p))
                        {
                            int cnt = property_counts.get(p);
                            property_counts.put(p, 1+cnt);
                        } else
                        {
                            property_counts.put(p, 1);
                            num_groups[combo_idx]++;
                        }
                    }
                }

                if (num_groups[combo_idx]==num_families)
                {
                    combo_descriptions[combo_idx] = "(every family has a different set of parameters)";
                } else if (num_groups[combo_idx] == 1)
                {
                    combo_descriptions[combo_idx]="("+(prop_idx==-1?"":"grouping by "+data_props[prop_idx]+": ")
                            +"same model applies to every family)";
                } else
                {
                    int min_size = num_families;
                    int max_size = 0;
                    for (String val: property_counts.keySet())
                    {
                        int cnt = property_counts.get(val);
                        if (cnt<min_size) min_size = cnt;
                        if (cnt>max_size) max_size = cnt;
                    }
                    combo_descriptions[combo_idx] = "(families grouped by "+data_props[prop_idx]
                            +": "+num_groups[combo_idx]+" different models for groups of size"
                            +(min_size<max_size
                                ?"s "+min_size+" ... "+max_size
                                :" "+min_size+" each")
                            +")";
                }
            } // props

            family_groupingCM = new GroupingCombo(combo_choices,combo_descriptions,num_groups);
            JLabel family_groupingL = family_groupingCM.getLabel();

            JPanel family_groupingP = new JPanel();
            family_groupingP.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(),"Family grouping",javax.swing.border.TitledBorder.LEFT,javax.swing.border.TitledBorder.TOP,getFont().deriveFont(Font.BOLD)));
            {
                GridBagLayout layout = new GridBagLayout();
                GridBagConstraints gc = new GridBagConstraints();
                family_groupingP.setLayout(layout);

                gc.fill = GridBagConstraints.NONE;
                gc.gridx = 0;
                gc.gridy = 0;
                gc.gridheight=1;
                gc.gridwidth=1;
                gc.anchor = GridBagConstraints.WEST;
                gc.weightx = 0.1;
                gc.weighty = 0.1;

                layout.setConstraints(family_groupingCM, gc);
                family_groupingP.add(family_groupingCM);

                gc.gridy++;
                layout.setConstraints(family_groupingL, gc);
                family_groupingP.add(family_groupingL);
            }

            return family_groupingP;
        }
        
        private class GroupingCombo extends JComboBox
        {
            GroupingCombo(String[] choices, String[] desc, int[] num_groups)
            {
                super(choices);
                this.desc = desc;
                this.num_groups = num_groups; 
                
                initComponents();
                initListener();
            }
            private String[] desc;
            private int[] num_groups;
            
            private JLabel grouping_infoL;
            private void initComponents()
            {
                this.setEditable(false);
                
                grouping_infoL = new JLabel(desc[0]);
                Font infoF = grouping_infoL.getFont();
                grouping_infoL.setFont(infoF.deriveFont(Font.ITALIC).deriveFont(0.8f*infoF.getSize()));
            }
            
            private void initListener()
            {
                this.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {
                        int idx = getSelectedIndex();
                        if (idx<0) // shouldn't happen
                        {
                            grouping_infoL.setText("");
                        } else
                        {
                            grouping_infoL.setText(desc[idx]);
                        }
                    }
                });
                
            }
            
            int getSelectedGroupCount()
            {
                int idx = getSelectedIndex();
                return (idx<0?0:num_groups[idx]);
            }
            
            JLabel getLabel()
            {
                return grouping_infoL;
            }
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
        private JComponent createFamilyVariationParameterBox()
        {
            NumberFormat alphaF = new DecimalFormat("0.############");
            alphaF.setParseIntegerOnly(false);
            alpha_lengthT = new ParameterField(alphaF,1.0,"Fixed");
            alpha_lengthT.setInputVerifier(InputVerifiers.getPositiveInstance());
            alpha_lengthT.setColumns(16);
            alpha_lengthT.setParameterName("Gamma distribution shape parameter for length variation across families");

            alpha_gainT = new ParameterField(alphaF, 1.0, "Fixed");
            alpha_gainT.setInputVerifier(InputVerifiers.getPositiveInstance());
            alpha_gainT.setColumns(16);
            alpha_gainT.setParameterName("Gamma distribution shape parameter for gain rate variation across families");
            alpha_duplicationT = new ParameterField(alphaF, 1.0, "Fixed");
            alpha_duplicationT.setInputVerifier(InputVerifiers.getPositiveInstance());
            alpha_duplicationT.setColumns(16);
            alpha_duplicationT.setParameterName("Gamma distribution shape parameter for duplication rate variation across families");
            alpha_lossT = new ParameterField(alphaF, 1.0, "Fixed");
            alpha_lossT.setInputVerifier(InputVerifiers.getPositiveInstance());
            alpha_lossT.setColumns(16);
            alpha_lossT.setParameterName("Gamma distribution shape parameter for loss rate variation across families");

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
            forbidden_gainT.setParameterName("Prior probability for forbidden-gain category");

            forbidden_duplicationT = new ParameterField(forbiddenF, 0.0, "Fixed");
            forbidden_duplicationT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
            forbidden_duplicationT.setColumns(16);
            forbidden_duplicationT.setParameterName("Prior probability for forbidden-duplication category");

            forbidden_gainT.freezeParameter();
            forbidden_duplicationT.freezeParameter();
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
                layout.setConstraints(alpha_lengthT.getCheckBox(), gc);
                rate_variationP.add(alpha_lengthT.getCheckBox());
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
                layout.setConstraints(alpha_lossT.getCheckBox(), gc);
                rate_variationP.add(alpha_lossT.getCheckBox());
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
                layout.setConstraints(alpha_gainT.getCheckBox(), gc);
                rate_variationP.add(alpha_gainT.getCheckBox());
                gc.gridx++;
                gc.fill = GridBagConstraints.HORIZONTAL;
                layout.setConstraints(forbidden_gainL, gc);
                rate_variationP.add(forbidden_gainL);
                gc.gridx++;
                gc.fill = GridBagConstraints.NONE;
                layout.setConstraints(forbidden_gainT,gc);
                rate_variationP.add(forbidden_gainT);
                gc.gridx++;
                layout.setConstraints(forbidden_gainT.getCheckBox(), gc);
                rate_variationP.add(forbidden_gainT.getCheckBox());

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
                layout.setConstraints(alpha_duplicationT.getCheckBox(), gc);
                rate_variationP.add(alpha_duplicationT.getCheckBox());
                gc.gridx++;
                gc.fill = GridBagConstraints.HORIZONTAL;
                layout.setConstraints(forbidden_duplicationL, gc);
                rate_variationP.add(forbidden_duplicationL);
                gc.gridx++;
                gc.fill = GridBagConstraints.NONE;
                layout.setConstraints(forbidden_duplicationT,gc);
                rate_variationP.add(forbidden_duplicationT);
                gc.gridx++;
                layout.setConstraints(forbidden_duplicationT.getCheckBox(), gc);
                rate_variationP.add(forbidden_duplicationT.getCheckBox());

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
         * Initializes {@link #root_aT}, {@link #root_bT}, {@link #fixed_root_CB},
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

            NumberFormat root_bF = new DecimalFormat("0.############");
            root_bF.setParseIntegerOnly(false);
            root_bT = new JFormattedTextField(root_bF);
            root_bT.setInputVerifier(new InputVerifiers.NonnegativeInputVerifier());
            root_bT.setColumns(16);

            root_distributionL = new JLabel("");
            root_aL = new JLabel("par0");
            root_bL = new JLabel("par1");

            fixed_root_CB = new JCheckBox("Fixed");

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

            return root_distributionP;
        }

        private JComponent createLineageSpecificParameterBox()
        {
            TreeWithRates tree = models.get(0).getMainTree();
            NodeWithRates[] nodes = tree.getDFT();
            int num_edges = tree.getNumEdges();

            lineage_lengthP = new LineageSpecificParameters(tree, InputVerifiers.getPositiveInstance());
            lineage_lengthP.setParameterName("Length");
            lineage_duplicationP = new LineageSpecificParameters(tree, InputVerifiers.getNonnegativeInstance()); // new InputVerifiers.RangeVerifier(ML.MIN_DUPLICATION_RATE, ML.MAX_DUPLICATION_RATE));
            lineage_duplicationP.setParameterName("Duplication-loss ratio");
            lineage_gainP = new LineageSpecificParameters(tree, InputVerifiers.getNonnegativeInstance()); // new InputVerifiers.RangeVerifier(ML.MIN_TRANSFER_RATE, ML.MAX_TRANSFER_RATE));
            lineage_gainP.setParameterName("Gain-loss ratio");
            lineage_lossP = new LineageSpecificParameters(tree, InputVerifiers.getPositiveInstance());
            lineage_lossP.setParameterName("Loss rate");

            // initial state corresponds to homogeneous loss check box unselected, gld model
            lineage_lossP.setDisabledReason("rates and lengths are scaled by loss rate in a heterogeneous model");
            lineage_lossP.freezeAll();

            // row headers
            JLabel[] node_nameL = new JLabel[num_edges];
            for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
                node_nameL[edge_idx] = new JLabel(LookAndFeel.getLongNodeName(tree, edge_idx));
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
                layout.setConstraints(lineage_lengthP.getMasterCheckBox(), gc);
                parametersP.add(lineage_lengthP.getMasterCheckBox());
                gc.gridx++;

                gc.gridx++;
//                layout.setConstraints(all_fixed_gainCB, gc);
//                parametersP.add(all_fixed_gainCB);
                layout.setConstraints(lineage_gainP.getMasterCheckBox(), gc);
                parametersP.add(lineage_gainP.getMasterCheckBox());
                gc.gridx++;

                gc.gridx++;
//                layout.setConstraints(all_fixed_duplicationCB, gc);
//                parametersP.add(all_fixed_duplicationCB);
                layout.setConstraints(lineage_duplicationP.getMasterCheckBox(), gc);
                parametersP.add(lineage_duplicationP.getMasterCheckBox());
                gc.gridx++;

                gc.gridx++;
//                layout.setConstraints(all_fixed_lossCB, gc);
//                parametersP.add(all_fixed_lossCB);
                layout.setConstraints(lineage_lossP.getMasterCheckBox(), gc);
                parametersP.add(lineage_lossP.getMasterCheckBox());
                gc.gridx++;

                for (int node_idx=0; node_idx<num_edges; node_idx++)
                {
                    NodeWithRates N = nodes[node_idx];
                    assert (!N.isRoot());

                    gc.gridx=0;
                    gc.gridy++;

                    gc.anchor = GridBagConstraints.WEST;
                    gc.fill = GridBagConstraints.HORIZONTAL;
                    layout.setConstraints(node_nameL[node_idx], gc);
                    parametersP.add(node_nameL[node_idx]);
                    gc.gridx++;
                    gc.fill = GridBagConstraints.NONE;
                    layout.setConstraints(lineage_lengthP.getField(node_idx), gc);
                    parametersP.add(lineage_lengthP.getField(node_idx));
                    gc.gridx++;
                    gc.fill = GridBagConstraints.HORIZONTAL;
                    layout.setConstraints(lineage_lengthP.getCheckBox(node_idx), gc);
                    parametersP.add(lineage_lengthP.getCheckBox(node_idx));
                    gc.gridx++;
                    gc.fill = GridBagConstraints.NONE;
                    layout.setConstraints(lineage_gainP.getField(node_idx), gc);
                    parametersP.add(lineage_gainP.getField(node_idx));
                    gc.gridx++;
                    gc.fill = GridBagConstraints.HORIZONTAL;
                    layout.setConstraints(lineage_gainP.getCheckBox(node_idx), gc);
                    parametersP.add(lineage_gainP.getCheckBox(node_idx));
                    gc.gridx++;
                    gc.fill = GridBagConstraints.NONE;
                    layout.setConstraints(lineage_duplicationP.getField(node_idx), gc);
                    parametersP.add(lineage_duplicationP.getField(node_idx));
                    gc.gridx++;
                    gc.fill = GridBagConstraints.HORIZONTAL;
                    layout.setConstraints(lineage_duplicationP.getCheckBox(node_idx), gc);
                    parametersP.add(lineage_duplicationP.getCheckBox(node_idx));
                    gc.gridx++;
                    gc.fill = GridBagConstraints.NONE;
                    layout.setConstraints(lineage_lossP.getField(node_idx), gc);
                    parametersP.add(lineage_lossP.getField(node_idx));
                    gc.gridx++;
                    gc.fill = GridBagConstraints.HORIZONTAL;
                    layout.setConstraints(lineage_lossP.getCheckBox(node_idx), gc);
                    parametersP.add(lineage_lossP.getCheckBox(node_idx));
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

        private class DistributionSelectionButton extends JRadioButton
        {
            DistributionSelectionButton(String txt, Class distribution_class)
            {
                super(txt);
                this.distribution_class=distribution_class;
            }
            private Class distribution_class;
            Class getDistributionClass()
            {
                return distribution_class;
            }
        }


        private class LineageSpecificParameters
        {

            LineageSpecificParameters(TreeWithRates tree, InputVerifier verifier)
            {
                initFields(tree, verifier);
            }
            private CheckSelectAll<JCheckBox> fixed_parameters;
            private ParameterField[] lineage_parameters;
            private boolean variation_ok = true;
            private String parameter_text;
            private String[] lineage_name;
            private String disabled_reason;

            void setParameterName(String txt)
            {
                this.parameter_text = txt;
                for (int j=0; j<lineage_parameters.length; ++j)
                {
                    ParameterField f = lineage_parameters[j];
                    f.setParameterName(parameter_text+" for lineage "+lineage_name[j]);
                }
            }

            void setDisabledReason(String disabled_reason)
            {
                this.disabled_reason = disabled_reason;
                for (ParameterField f: lineage_parameters)
                {
                    f.setDisabledReason(disabled_reason);
                }
            }

            private void initFields(TreeWithRates tree, InputVerifier verifier)
            {
                int num_edges = tree.getNumEdges();
                lineage_parameters = new ParameterField[num_edges];
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
                for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
                {
                    ParameterField f = new ParameterField(fieldF, 0.0, "");
                    f.setInputVerifier(verifier);
                    f.setColumns(16);
                    lineage_parameters[edge_idx] = f;
                    fixed_parameters.addCheckBox(f.getCheckBox());
                }
                fixed_parameters.addItemListener(new ItemListener()
                {
                    @Override
                    public void itemStateChanged(ItemEvent e)
                    {
                        if (!variation_ok && ItemEvent.DESELECTED == e.getStateChange())
                            fixed_parameters.copySelectedToAll();
                    }
                });
                lineage_parameters[0].addPropertyChangeListener("value",new PropertyChangeListener()
                {
                    @Override
                    public void propertyChange(PropertyChangeEvent e)
                    {
                        if (!variation_ok)
                        {
                            Object val = e.getNewValue();
                            for (int i=1; i<lineage_parameters.length; ++i)
                                lineage_parameters[i].setValue(val);
                        }
                    }
                });
                lineage_name = new String[num_edges];
                for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
                {
                    lineage_name[edge_idx] = LookAndFeel.getShortNodeName(tree, edge_idx);
                }
            }

            /**
             * The parameter field associated with a lineage
             *
             * @param edge_idx index by {@link TreeWithRates}
             * @return corresponding parameter field
             */
            ParameterField getField(int edge_idx)
            {
                return lineage_parameters[edge_idx];
            }

            /**
             * The master checkbox for "fixed"
             *
             * @return checkbox
             */
            CheckSelectAll<JCheckBox> getMasterCheckBox()
            {
                return fixed_parameters;
            }

            /**
             * The "Fixed" checkbox associated with a lineage.
             *
             * @param edge_idx lineage index
             * @return associated fixed checkbox
             */
            JCheckBox getCheckBox(int edge_idx)
            {
                return lineage_parameters[edge_idx].getCheckBox();
            }

            /**
             * Resets all paramtere field values.
             */
            void resetAll()
            {
                for (ParameterField f:lineage_parameters)
                    f.resetValue();
            }

            void freezeAll()
            {
                for (ParameterField f: lineage_parameters)
                    f.freezeParameter();
                fixed_parameters.setSelected(true);
                fixed_parameters.setEnabled(false);
            }

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
                    lineage_parameters[0].getCheckBox().setEnabled(false);
                    for (int i=1; i<lineage_parameters.length; ++i)
                        lineage_parameters[i].freezeParameter();
                    fixed_parameters.copySelectedToAll();
                }
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

            /**
             * Enables/disables lineage-specific variation.
             * Call after {@ink #freezeAll()} or {@link #enableAll() }.
             *
             * @param variation_enabled
             */
            void setVariationEnabled(boolean variation_enabled)
            {
                this.variation_ok = variation_enabled;
                if (fixed_parameters.isEnabled())
                    enableAll();
            }

            /**
             * Enables or disables lineage-specific variation
             * by a checkbox state.
             *
             * @param cb checkbox giving the selected state used
             * @param enable_when_selected variation is enabled if selected stated equals this parameter, otherwise disable
             */
            void synchronizeTo(JCheckBox cb, boolean enable_when_selected)
            {
                boolean b = cb.isSelected() == enable_when_selected;
                setVariationEnabled(b);
                if (!b)
                    for (int i=1; i<lineage_parameters.length; ++i)
                        lineage_parameters[i].setDisabledReason("lineage-specific rates are not allowed");
            }

            void equalizeValues(LineageSpecificParameters weighting)
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

        /**
         * Sets the state for displayed fields and labels for a given root prior distribution.
         * Field values are not set if they are null.
         * Field values are set (to 0.1) only if the new distribution class has
         * invalid parameter in the field.
         *
         * Used by listeners of {@link #attachRootPriorListeners() }.
         *
         */
        private void synchronizeRootPriorStates()
        {

            Class distribution_class = getRootPriorClass();
            if (distribution_class==Poisson.class)
            {
                root_distributionL.setText("Poisson");
                root_aL.setText("\u03bb");
                root_aT.setInputVerifier(InputVerifiers.getPositiveInstance());
                root_aT.setToolTipText("Poisson parameter (expected value)");
                root_bL.setText("nosuch");
                root_bT.setVisible(false);
                root_bL.setVisible(false);
            } else if (distribution_class==NegativeBinomial.class)
            {
                root_distributionL.setText("Negative binomial");
                root_aL.setText("\u03b8");
                root_aT.setInputVerifier(InputVerifiers.getPositiveInstance());
                root_aT.setToolTipText("Stopping-time parameter for negative binomial");
                root_bL.setText("q");
                root_bT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
                root_bT.setToolTipText("Sucess-probability parameter for negative binomial");
                root_bT.setVisible(true);
                root_bL.setVisible(true);
            } else if (distribution_class == PointDistribution.class)
            {
                root_distributionL.setText("Point");
                root_aL.setText("p0");
                root_aT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
                root_aT.setToolTipText("Failure-probability parameter for Bernoulli");
                root_bT.setVisible(false);
                root_bL.setVisible(false);
            } else if (distribution_class == ShiftedGeometric.class)
            {
                root_distributionL.setText("Shifted geometric");
                root_aL.setText("p0");
                root_aT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
                root_aT.setToolTipText("Failure-probability parameter for shfted geometric");
                root_bL.setText("q");
                root_bT.setInputVerifier(new InputVerifiers.RangeVerifier(0.0, 1.0));
                root_bT.setToolTipText("Success probability parameter for shifted geometric");
                root_bT.setVisible(true);
                root_bL.setVisible(true);
            }
            else
                throw new IllegalArgumentException("setRootdistributionComponentStates: distribution_class must be Poisson, NegativeBinomial or PointDistribution");

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
        }

        /**
         * Sets the states for fields and checkboxes related to gain
         * and duplication rate in model. Field values may be reset
         * if a rate feature is rmeoved from the model.
         *
         * Used by listeners of {@link #attachModelStructureListeners() }.
         */
        private void synchronizeModelStructureStates()
        {
            boolean has_gain = model_gldRB.isSelected() || model_glRB.isSelected();
            boolean has_duplication = model_gldRB.isSelected() || model_dlRB.isSelected();
//            System.out.println("#*MSP.MD.sMSS gain "+has_gain+"\tdup "+has_duplication);

            if (has_gain)
            {
                gamma_gainT.setEditable(true);
                Number g = (Number)gamma_gainT.getValue();
                if (g!=null && g.intValue()!=1)
                    alpha_gainT.enableParameter();
                forbidden_gainCB.setEnabled(true);
                lineage_gainP.enableAll();
                uniform_gainCB.setEnabled(true);
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
                lineage_gainP.resetAll();
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
            } else
            {
                alpha_duplicationT.freezeParameter();
                gamma_duplicationT.setEditable(false);
                gamma_duplicationT.setValue(1);
                forbidden_duplicationCB.setEnabled(false);
                forbidden_duplicationCB.setSelected(false);
                lineage_duplicationP.setDisabledReason("model does not allow duplications");
                lineage_duplicationP.freezeAll();
                lineage_duplicationP.resetAll();
                uniform_duplicationCB.setEnabled(false);
                uniform_duplicationCB.setSelected(true);
            }
        }

        private void synchronizeLineageVariationState()
        {
//            System.out.println("#*MSP.MD.sLVS");
            lineage_duplicationP.synchronizeTo(uniform_duplicationCB, false);
            lineage_gainP.synchronizeTo(uniform_gainCB, false);
            boolean homogeneous_loss = uniform_lossCB.isSelected();
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
            lineage_lossP.synchronizeTo(uniform_lossCB, false);
        }

        /**
         * Attaches ActionListeners reacting to
         * user's changing the prior distribution
         * ({@link #root_distributionB}).
         */
        private void attachRootPriorListeners()
        {
            ActionListener root_listener = new ActionListener()
            {
                @Override
                public void actionPerformed(ActionEvent e)
                {
                    synchronizeRootPriorStates();
                }
            };
            for (DistributionSelectionButton b: root_distributionB)
                b.addActionListener(root_listener);
        }

        /**
         * Attaches ActionListeners to {@link #cancel_button} and {@link #start_button}.
         */
        private void attachButtonListeners()
        {
            cancel_button.addActionListener(new ActionListener()
            {
                @Override
                public void actionPerformed(ActionEvent e)
                {
                    selected_model = -1;
                    dispose();
                }
            });

            start_button.addActionListener(new ActionListener()
            {
                @Override
                public void actionPerformed(ActionEvent e)
                {
                    copyFieldValuesIntoModel();
                    dispose();
                }
            });
        }

        private void attachInitialModelListeners()
        {
            class InitialModelListener implements ItemListener
            {
                InitialModelListener(int idx)
                {
                    this.index = idx;
                }
                private int index;
                @Override
                public void itemStateChanged(ItemEvent e)
                {
                    if (e.getStateChange() == ItemEvent.SELECTED)
                    {
                        setSelectedModelIndex(index);
                    }
                }
            }

            for (int idx=0; idx<modelRB.size(); ++idx)
            {
                JRadioButton rb = modelRB.get(idx);
                rb.addItemListener(new InitialModelListener(idx));
            }
        }

        /**
         * Attaches ActionListeners to the model structure buttons
         * {@link #model_gldRB}, {@link #model_dlRB}, {@link #model_glRB}
         * and {@link #model_plRB}.
         */
        private void attachModelStructureListeners()
        {
            ActionListener gain_loss_listener = new ActionListener()
            {
                @Override
                public void actionPerformed(ActionEvent e)
                {
                    synchronizeModelStructureStates();
                    synchronizeLineageVariationState();
                }
            };
            model_gldRB.addActionListener(gain_loss_listener);
            model_dlRB.addActionListener(gain_loss_listener);
            model_glRB.addActionListener(gain_loss_listener);
            model_plRB.addActionListener(gain_loss_listener);
        }

        /**
         * Attaches the listeners to parameters related to
         * rate variation across families:
         * {@link #alpha_lengthT}, {@link #alpha_duplicationT},
         * {@link #alpha_gainT}, {@link #alpha_lossT},
         * {@link #forbidden_gainT} and {@link #forbidden_gainT}.
         */
        private void attachFamilyVariationListeners()
        {
            alpha_gainT.listenTo(gamma_gainT, 1, "gain rates do not vary across families");
            alpha_duplicationT.listenTo(gamma_duplicationT, 1, "duplication rates do not vary across families");
            alpha_lossT.listenTo(gamma_lossT, 1, "loss rates do not vary across families");
            alpha_lengthT.listenTo(gamma_lengthT, 1, "edge lengths do not vary across families");

            forbidden_duplicationT.listenTo(forbidden_duplicationCB, true, "there is no forbidden-duplication category");
            forbidden_gainT.listenTo(forbidden_gainCB, true, "there is no forbidden-gain category");
        }


        /**
         * Attaches listeners to disable/enable lineage-specific variations
         * of loss, duplication and gain.
         */
        private void attachLineageVariationListeners()
        {
            class RateEqualizer implements ItemListener
            {
                RateEqualizer(LineageSpecificParameters params)
                {
                    this.params = params;
                }
                private LineageSpecificParameters params;
                @Override
                public void itemStateChanged(ItemEvent e)
                {
                    JCheckBox uniform_ratesCB = (JCheckBox) e.getSource();
                    params.synchronizeTo(uniform_ratesCB, false);
                    if (!uniform_ratesCB.isSelected())
                        params.equalizeValues(lineage_lengthP);
                }
            }

            uniform_duplicationCB.addItemListener(new RateEqualizer(lineage_duplicationP));
            uniform_gainCB.addItemListener(new RateEqualizer(lineage_gainP));
            ItemListener length_selector = new ItemListener()
            {
                @Override
                public void itemStateChanged(ItemEvent e)
                {
                    if (ItemEvent.SELECTED == e.getStateChange()
                            || ItemEvent.DESELECTED == e.getStateChange())
                    {
                        RateVariation model = models.get(selected_model);
                        TreeWithRates tree = model.getMainTree();
                        int num_edges = tree.getNumEdges();
                        boolean homogeneous_rates = ((JCheckBox)e.getSource()).isSelected();
                        if (homogeneous_rates)
                        {
                            lineage_lengthP.setDisabledReason("model is homogeneous");
                            lineage_lengthP.freezeAll();
                            lineage_lossP.enableAll();
                            lineage_lossP.synchronizeTo(uniform_lossCB, false);

                            WorkSpaceCount ws = dealer.getActiveWorkSpaceCount();
                            double original_total_length = 0.0;
                            double total_loss = 0.0;
                            for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
                            {
                                double len = ((Number)lineage_lengthP.getField(edge_idx).getValue()).doubleValue();
                                double original_length = ws.getOriginalEdgeLength(edge_idx);
                                lineage_lengthP.getField(edge_idx).setValue(original_length);
                                original_total_length += original_length;
                                double ls = len*((Number)lineage_lossP.getField(edge_idx).getValue()).doubleValue();
                                total_loss += ls;
                            }
                            double common_loss_rate = total_loss / original_total_length;
                            lineage_lossP.setValueAll(common_loss_rate);
                        } else
                        {
                            lineage_lengthP.enableAll();
                            lineage_lengthP.getMasterCheckBox().setSelected(false);
                            lineage_lengthP.getMasterCheckBox().copySelectedToAll(); // can change now
                            lineage_lossP.setDisabledReason("rates and lengths are scaled by loss rate in a heterogeneous model");
                            lineage_lossP.freezeAll();
                            lineage_lossP.synchronizeTo(uniform_lossCB, false);

                            for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
                            {
                                NodeWithRates N = tree.getNode(edge_idx);
                                double len = ((Number)lineage_lengthP.getField(edge_idx).getValue()).doubleValue();
                                double new_len = N.getLength();
                                lineage_lengthP.getField(edge_idx).setValue(new_len);
                                double ls = len*((Number)lineage_lossP.getField(edge_idx).getValue()).doubleValue();
                                double new_ls = ls / new_len;
                                lineage_lossP.getField(edge_idx).setValue(new_ls);
                            }
                        }
                    }
                }
            };
            uniform_lossCB.addItemListener(length_selector);
        }

        /**
         * Performs a programmating click for selecting
         * the model radio button that corresponds to the selected model.
         */
        private void clickModelType()
        {
            RateVariation model  = models.get(selected_model);
            if (model.hasDuplication())
            {
                if (model.hasGain())
                    model_gldRB.doClick();
                else
                    model_dlRB.doClick();
            } else
            {
                if (model.hasGain())
                    model_glRB.doClick();
                else
                    model_plRB.doClick();
            }
        }

        /**
         * Sets the selected state of the uniform rate checkboxes
         * according to the selected model.
         */

        private void selectLineageVariation()
        {
            RateVariation model  = models.get(selected_model);
            uniform_duplicationCB.setSelected(!model.hasLineageSpecificDuplication());
            uniform_gainCB.setSelected(!model.hasLineageSpecificGain());
            uniform_lossCB.setSelected(!model.hasLineageSpecificLoss());
        }

        private void clickRootPrior()
        {
            RateVariation model  = models.get(selected_model);
            Class root_distribution = model.getRootPrior().getClass();
            for (DistributionSelectionButton b:root_distributionB)
            {
                if (root_distribution.equals(b.getDistributionClass()))
                {
                    b.doClick();
                    break;
                }
            }
        }

        private void selectFamilyVariation()
        {
            RateVariation model  = models.get(selected_model);
            gamma_lengthT.setValue(model.getNumEdgeLengthGammaCategories());
            gamma_duplicationT.setValue(model.getNumDuplicationRateGammaCategories());
            gamma_gainT.setValue(model.getNumTransferRateGammaCategories());
            gamma_lossT.setValue(model.getNumLossRateGammaCategories());

            if (model.getDuplicationForbiddenProportion()>0.0)
                forbidden_duplicationCB.setSelected(true);
            if (model.getTransferForbiddenProportion()>0.0)
                forbidden_gainCB.setSelected(true);
        }

        private void fillRootPriorParameters()
        {
            RateVariation model  = models.get(selected_model);
            DiscreteDistribution root_prior = model.getRootPrior();
            double[] params = root_prior.getParameters();
            root_aT.setValue(params[0]);
            if (params.length>1)
                root_bT.setValue(params[1]);
        }

        private void fillFamilyVariationParameters()
        {
            RateVariation model  = models.get(selected_model);
            alpha_lengthT.setValue(model.getEdgeLengthAlpha());
            alpha_duplicationT.setValue(model.getDuplicationRateAlpha());
            alpha_gainT.setValue(model.getTransferRateAlpha());
            alpha_lossT.setValue(model.getLossRateAlpha());

            forbidden_duplicationT.setValue(model.getDuplicationForbiddenProportion());
            forbidden_gainT.setValue(model.getTransferForbiddenProportion());
        }

        private void fillLineageSpecificParameters()
        {
//            System.out.println("#*MSP.MD.fLSP "+selected_model);
            RateVariation model  = models.get(selected_model);
            TreeWithRates tree = model.getMainTree();
            int num_edges = tree.getNumEdges();
            boolean homogeneous_model = !model.hasLineageSpecificLoss();
            if (homogeneous_model)
            {
                WorkSpaceCount ws = dealer.getActiveWorkSpaceCount();
                double original_total_length = 0.0;
                double total_length = 0.0;
                double total_loss = 0.0;
                for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
                {
                    NodeWithRates N = tree.getNode(edge_idx);
                    assert (!N.isRoot());
                    double len = N.getLength();
                    total_length += len;
                    double original_length = ws.getOriginalEdgeLength(edge_idx);
                    lineage_lengthP.getField(edge_idx).setValue(original_length);
                    original_total_length += original_length;
                    double ls = len*N.getLossRate();;
                    total_loss += ls;
                }
                double common_loss_rate = total_loss / original_total_length;
                lineage_lossP.setValueAll(common_loss_rate);
            } else
            {
                assert (!uniform_lossCB.isSelected());
                for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
                {
                    NodeWithRates N = tree.getNode(edge_idx);
                    assert (!N.isRoot());

                    double len = N.getLength();
                    lineage_lengthP.getField(edge_idx).setValue(len);
                    double ls = N.getLossRate();
                    lineage_lossP.getField(edge_idx).setValue(ls);
                }
            }

            for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
            {
                NodeWithRates N = tree.getNode(edge_idx);
                assert (!N.isRoot());
                double gr = N.getTransferRate();
                double dr = N.getDuplicationRate();
                lineage_gainP.getField(edge_idx).setValue(gr);
                lineage_duplicationP.getField(edge_idx).setValue(dr);
            }
        }

        private void setOptimizationParameters(int num_rounds, double epsilon)
        {
            roundT.setValue(num_rounds);
            epsT.setValue(epsilon);
            approximateCB.setSelected(true);
            approximate_thresholdT.setValue(ApproximateComputation.DEFAULT_CUTOFF_THRESHOLD);
            approximate_factorT.setValue(ApproximateComputation.DEFAULT_CUTOFF_FACTOR);
        }
    }



    /**
     * This a class for a text field
     * containing a model parameter
     * and an associated "fixed" check box.
     *
     * Tooltip text is set to parameter name, other info
     * (enabled, fixed) is attached automatically both for the field and
     * the checkbox.
     *
     */
    private class ParameterField extends JFormattedTextField
    {
        /**
         * Instantiation with 0 default value.
         *
         * @param format format for this field
         */
        ParameterField(Format format)
        {
            this(format, 0);
        }


        /**
         * Instantiation with no check box text.
         *
         * @param format format for the parameter field
         * @param default_value default value for the field (used by {@link #resetValue() })
         */
        ParameterField(Format format, Object default_value)
        {
            this(format, default_value, "");
        }

        /**
         * Instantiation.
         *
         * @param format field format
         * @param check_box_txt text for associated check box
         * @param default_value default value for parameter (when disabled)
         */
        ParameterField(Format format, Object default_value, String check_box_text)
        {
            super(format);
            initCheckBox(check_box_text);
            setDefaultValue(default_value);
            setParameterName("Rate parameter");
            setToolTipText("");
        }

        private JCheckBox fixedCB;
        private Object default_value;
        private String disabled_reason;
        private String parameter_text;

        public double doubleValue()
        {
            return ((Number)getValue()).doubleValue();
        }

        public double intValue()
        {
            return ((Number)getValue()).intValue();
        }

        /**
         * Sets the name of this parameter.
         *
         * @param txt parameter name, not null
         */
        void setParameterName(String txt)
        {
            this.parameter_text = txt;
        }

        /**
         * Sets the explanation for why the field is disabled.
         *
         * @param txt reason for disabling this parameter
         */
        void setDisabledReason(String txt)
        {
            this.disabled_reason = txt;
        }

        private void initCheckBox(String txt)
        {
            fixedCB=new JCheckBox(txt)
            {
                @Override
                public String getToolTipText(MouseEvent ignored)
                {
                    return getCheckBoxToolTip();
                }
            };
            fixedCB.setToolTipText("");
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


        void setDefaultValue(Object value)
        {
            this.default_value = value;
        }

        /**
         * Resets to default value.
         */
        void resetValue()
        {
           setValue(default_value);
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
         * The associated check box, set at instantiation.
         *
         * @return a check box
         */
        JCheckBox getCheckBox()
        {
            return fixedCB;
        }

        boolean isFixed()
        {
            return fixedCB.isSelected();

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
        void listenTo(JFormattedTextField controlTF, final Object disable_for_value, final String disabled_reason)
        {
            controlTF.addPropertyChangeListener("value",
                new PropertyChangeListener()
                {
                    @Override
                    public void propertyChange(PropertyChangeEvent update)
                    {
                        if (update.getNewValue().equals(disable_for_value))
                        { // no more editing
                            freezeParameter();
                            setDisabledReason(disabled_reason);
                        } else
                        {
                            enableParameter();
                        }
                    }
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
        void listenTo(javax.swing.AbstractButton controlSEL, boolean is_enabler, final String disabled_reason)
        {
            if (is_enabler)
            {
                controlSEL.addItemListener(new ItemListener()
                    {
                        @Override
                        public void itemStateChanged(ItemEvent e)
                        {
                            if (ItemEvent.SELECTED == e.getStateChange())
                                enableParameter();
                            else if (ItemEvent.DESELECTED == e.getStateChange())
                            {
                                freezeParameter();
                                resetValue();
                                setDisabledReason(disabled_reason);
                            }
                        }
                });
            } else
            {
                controlSEL.addItemListener(new ItemListener()
                    {
                        @Override
                        public void itemStateChanged(ItemEvent e)
                        {
                            if (ItemEvent.DESELECTED == e.getStateChange())
                                enableParameter();
                            else if (ItemEvent.SELECTED == e.getStateChange())
                            {
                                freezeParameter();
                                resetValue();
                                setDisabledReason(disabled_reason);
                            }
                        }
                });
            }
        }
    }


}
