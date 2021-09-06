package ca.umontreal.iro.evolution.malin.ui.count;

/**
 *
 * @author csuros
 */

import java.beans.PropertyChangeEvent;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import ca.umontreal.iro.banality.StringSplit;

import ca.umontreal.iro.evolution.genecontent.OccurrenceTable;
import ca.umontreal.iro.evolution.genecontent.Posteriors;
import ca.umontreal.iro.evolution.genecontent.RateVariation;
import ca.umontreal.iro.evolution.genecontent.TreeWithRates;
        
import ca.umontreal.iro.evolution.malin.DataFile;
import ca.umontreal.iro.evolution.malin.DoubleRoundedForDisplay;

import ca.umontreal.iro.evolution.malin.ui.ZoomableTreePanel;

public class PosteriorDisplay extends AncestralReconstructionPane 
{
    public PosteriorDisplay(DataFile<OccurrenceTable> data_file, DataFile<RateVariation> rate_file, WorkSpaceCount W)
    {
        super(rate_file.getData().getMainTree(), data_file, W);
        //Verbose.setVerbose(true);
        this.rate_file = rate_file;
        init();
    }
    
    protected PosteriorDisplay( TreeWithRates main_tree, DataFile<OccurrenceTable> data_file,WorkSpaceCount W)
    {
        super(main_tree, data_file, W);
        this.rate_file = null;
        init();
    }
    
    
    private DataFile<RateVariation> rate_file;
    private Posteriors posteriors;
    
    public DataFile<RateVariation> getRateModel()
    {
        return rate_file;
    }
    
    /**
     * Probability for more than one representative in the family at the nodes
     */
    private double[][] reconstructed_family_multi;

    /**
     * Probability for a first member appearing at the nodes
     */
    private double[][] reconstructed_family_gain;
    
    /**
     * Probability for a first member appearing at the nodes
     */
    private double[][] reconstructed_family_loss;
    
    /**
     * Probability for family size expansion (but no gain) at the nodes
     */
    private double[][] reconstructed_family_expansion;

    /**
     * Probability for family size reduction (but no loss) at the nodes
     */
    private double[][] reconstructed_family_reduction;

    /**
     * Statistics for all-0 profile: return value from Posteriors.getAllPosteriors() 
     * 
     * Array values are indexed as [e][n] where n is node 
     * and e is the ordinal value for a Posteriors.PosteriorType.
     */
    private double[][] reconstructed_family_absent;
    
    private double[][] posterior_family_class;
    
    private double[] posterior_absent_class;

    /**
     * Probability for an all-0 profile
     */
    private double probability_absent;
    
    /**
     * Whether the displayed totals in the lineage table include absent families
     */
    protected JCheckBox lineage_totals_include_absent;

    @Override
    public String toString()
    {
        return "Posteriors @ "+rate_file.getFile().getName();
    }
            
    @Override
    protected void init()
    {
        super.init();
    }
    
    @Override
    protected void initGeneralVariables()
    {
        super.initGeneralVariables();
        posteriors = null;
    }
    
    /**
     * Here we check if rate file is set yet, and do the computation only if it is
     */
    @Override
    protected void computeAll()
    {
        if (rate_file != null)
            super.computeAll();
    }
    
    /**
     * Superclass's method disabled because  
     * there is nothing to initialize here: 
     * initComputeAll() will reallocate the reconstructed)[] arrays.
     */
    @Override
    protected void initReconstructionVariables()
    {
        int num_nodes = tree_nodes.length;
        int num_edges = tree_nodes.length-1;
        int num_families = families.length;
                
        reconstructed_family_count = new double[num_families][num_nodes];
        reconstructed_family_multi = new double[num_families][num_nodes];
        for (int j=0; j<reconstructed_family_count.length ; j++)
            reconstructed_family_count[j] = null;
        
        posterior_family_class = new double[families.length][];
    }
    
    @Override
    protected double getReconstructedFamilyMulti(int family_idx, int node_idx)
    {
        if (reconstructed_family_multi[family_idx]!=null)
            return reconstructed_family_multi[family_idx][node_idx];
        else
            return Double.NaN;
    }
    
    
    @Override
    protected double getReconstructedFamilyGain(int family_idx, int node_idx)
    {
        if (reconstructed_family_gain[family_idx]!=null)
            return reconstructed_family_gain[family_idx][node_idx];
        else
            return Double.NaN;
    }

    /**
     * Total number of gains across all lineages
     * 
     * @param family_idx index of the family
     * @return expected number of lineages in which the family was gained
     */
    protected double getReconstructedFamilyGain(int family_idx)
    {
        int num_nodes = tree_nodes.length;
        double g = 0.0;
        for (int node_idx=0; node_idx<num_nodes; node_idx++)
        {
            g+= getReconstructedFamilyGain(family_idx, node_idx);
        }
        return g;
    }
    
    
    @Override
    protected double getReconstructedFamilyLoss(int family_idx, int node_idx)
    {
        if (reconstructed_family_loss[family_idx]!=null)
            return reconstructed_family_loss[family_idx][node_idx];
        else
            return Double.NaN;
    }


    /**
     * Total number of losses across all lineages
     * 
     * @param family_idx index of the family
     * @return expected number of lineages in which the family was lost
     */
    protected double getReconstructedFamilyLoss(int family_idx)
    {
        int num_nodes = tree_nodes.length;
        double l = 0.0;
        for (int node_idx=0; node_idx<num_nodes; node_idx++)
        {
            l += getReconstructedFamilyLoss(family_idx, node_idx);
        }
        return l;
    }
    
    @Override
    protected double getReconstructedFamilyExpansion(int family_idx, int node_idx)
    {
        if (reconstructed_family_expansion[family_idx]!=null)
            return reconstructed_family_expansion[family_idx][node_idx];
        else
            return Double.NaN;
    }

    /**
     * Total number of expansions across all lineages
     * 
     * @param family_idx index of the family
     * @return expected number of lineages in which the family expanded
     */
    protected double getReconstructedFamilyExpansion(int family_idx)
    {
        int num_nodes = tree_nodes.length;
        double e = 0.0;
        for (int node_idx=0; node_idx<num_nodes; node_idx++)
        {
            e += getReconstructedFamilyExpansion(family_idx, node_idx);
        }
        return e;
    }
    
    @Override
    protected double getReconstructedFamilyReduction(int family_idx, int node_idx)
    {
        if (reconstructed_family_reduction[family_idx]!=null)
            return reconstructed_family_reduction[family_idx][node_idx];
        else
            return Double.NaN;
    }


    /**
     * Total number of reductions across all lineages
     * 
     * @param family_idx index of the family
     * @return expected number of lineages in which the family contracted
     */
    protected double getReconstructedFamilyReduction(int family_idx)
    {
        int num_nodes = tree_nodes.length;
        double r = 0.0;
        for (int node_idx=0; node_idx<num_nodes; node_idx++)
        {
            r += getReconstructedFamilyReduction(family_idx, node_idx);
        }
        return r;
    }
    
    /**
     * Total number of arrivals (gains + presence at root) across all lineages
     * 
     * @param family_idx index of the family
     * @return expected number of times the family newly appeared
     */
    protected double getReconstructedFamilyArrivals(int family_idx)
    {
        int num_nodes = tree_nodes.length;
        int root_idx = num_nodes-1;
        return getReconstructedFamilyGain(family_idx)+this.getReconstructedFamilyCount(family_idx, root_idx);
    }
    
    
    /**
     * Modifies super's method to include absent family information
     * if necessary.
     * 
     * @param selected_rows rows in which the summing is carried out
     */
    @Override
    protected void sumInRows(int[] selected_rows)
    {
        super.sumInRows(selected_rows);
        //System.out.println("#*PD.sIR "+selected_rows.length);
        if (lineage_totals_include_absent.isSelected() && reconstructed_family_absent != null)
        {
            double exp_missing = probability_absent * selected_rows.length/(1.-probability_absent);
            //System.out.println("#*PD.sIR "+selected_rows.length+" expected missing "+exp_missing);
            for (int node_idx =0; node_idx<tree_nodes.length; node_idx++)
            {
                selected_family_count[node_idx] += exp_missing * reconstructed_family_absent[Posteriors.PosteriorType.PRESENT.ordinal()][node_idx];
                selected_family_multi[node_idx] += exp_missing * reconstructed_family_absent[Posteriors.PosteriorType.MULTI.ordinal()][node_idx];
            }
            for (int edge_idx=0; edge_idx<tree_nodes.length-1; edge_idx++)
            {
                selected_family_gain[edge_idx] += exp_missing * reconstructed_family_absent[Posteriors.PosteriorType.GAIN.ordinal()][edge_idx];
                selected_family_loss[edge_idx] += exp_missing * reconstructed_family_absent[Posteriors.PosteriorType.LOSS.ordinal()][edge_idx];
                selected_family_expansion[edge_idx] += exp_missing * reconstructed_family_absent[Posteriors.PosteriorType.EXPANSION.ordinal()][edge_idx];
                selected_family_reduction[edge_idx] += exp_missing * reconstructed_family_absent[Posteriors.PosteriorType.REDUCTION.ordinal()][edge_idx];
            }
        }
    }    

    @Override
    protected void initComputeAll()
    {
        super.initComputeAll();
        initAllReconstructionVariables();
    }

    /**
     * Allocates reconstructed_family_XXX with XXX=count, muli, gain, loss, expansion and reduction
     */
    private void initAllReconstructionVariables()
    {
        int n = data_file.getData().getNumFamilies();
        reconstructed_family_count = new double[n][];
        reconstructed_family_multi = new double[n][];
        reconstructed_family_gain = new double[n][];
        reconstructed_family_loss = new double[n][];
        reconstructed_family_expansion = new double[n][];
        reconstructed_family_reduction = new double[n][];
    }
        
    @Override 
    protected void finishComputeAll()
    {
        super.finishComputeAll();
        posteriors = null;
    }
    
    @Override
    protected AncestralReconstructionPane.ComputingTask newComputingTask()
    {
        return new PosteriorComputingTask();
    }
    

    @Override
    public void propertyChange(PropertyChangeEvent evt)
    {
        
        String prop = evt.getPropertyName();
        //System.out.println("#*PD.pC "+prop+"\t// evt "+evt);
        if ("progress".equals(prop))
        {
            int pct = (Integer) evt.getNewValue();
            computation_progress.setIndeterminate(false);
            //System.out.println("#*PD.pC "+prop+" percentage "+pct);
            
            computation_progress.setValue(pct);
        } else if ("compute:class".equals(prop) || "compute:post".equals(prop))
        {
            String phase = ("compute:class".equals(prop)?"\u2190":"\u2192");
            String progress_message = phase + " "+evt.getNewValue();
            //System.out.println("#*PD.pC "+prop+"\t"+progress_message);
            computation_progress.setString(progress_message);
        }
    }
    
    protected class PosteriorComputingTask extends ComputingTask
    {
        /**
         * We estimate the remaining time by assuming a quadratic behavior in 
         * the family size. These values are computed by prepareComputation()
         */
        private long   total_work; // used to estimate remaining time
        private long[] cumulative_work;
        
        /**
         * Initializes the <var>posteriors</var> variable 
         */
        @Override
        protected void prepareComputation()
        {
//            System.out.println("#*PD.PCT.pC ");
            // let's estimate how long this will take 
            cumulative_work = new long[families.length]; 
            int num_leaves = main_tree.getNumLeaves();
            for (int family_idx=0; family_idx<families.length; family_idx++)
            {
                long s = member_count[family_idx]+num_leaves;
                cumulative_work[family_idx] = s*s;
            }
            for (int family_idx=1; family_idx<families.length; family_idx++)
                cumulative_work[family_idx]+=cumulative_work[family_idx-1];
            total_work = cumulative_work[families.length-1];
//            System.out.println("#*PD.PCT.pC total "+total_work);
            OccurrenceTable original_table =  data_file.getData();
//            System.out.println("#*PD.PCT.pC table class "+original_table.getClass().getCanonicalName()+"\twith "+original_table.getNumFamilies());
            // need to make a copy because Posteriors() will insert an ABSENT family
            OccurrenceTable filtered_table = original_table.filterByMaximumSize(Integer.MAX_VALUE, 0);
//            System.out.println("#*PD.PCT.pC filtered class "+filtered_table.getClass()+"\twith "+filtered_table.getNumFamilies()+" families");
//            for (int pidx=0; pidx<filtered_table.getNumFamilies(); pidx++)
//            {
//                System.out.println("#*PD.PCT.pC pidx "+pidx+"\t"+filtered_table.getProfile(pidx).getPatternString());
//            }

            posteriors = new Posteriors(filtered_table, rate_file.getData(), 20)
            {
                /**
                 * Overridden to display progree information.
                 */
                @Override
                protected void setSavedPosteriorClass(int pidx)
                {
                    super.setSavedPosteriorClass(pidx);
                    pidx--;
                    if (pidx>=0) // ABSENT
                    {
                        long work_done = total_work-(pidx==0?0:cumulative_work[pidx-1]);
                        int pct = (int)((50.0*work_done)/total_work);
//                        System.out.println("#*PD.PCT.sSPC "+pidx+"\t"+pct+"\t"+work_done+"/"+total_work);
                        setProgress(pct);
                        firePropertyChange("compute:class",null,families[pidx]+" ("+Integer.toString(families.length-pidx)+"/"+Integer.toString(families.length)+")");
                    }
                }
            };
            //System.out.println("#*PD.PCT.pC posteriors init");

            reconstructed_family_absent = posteriors.getAllPosteriors(0);
            posterior_absent_class = posteriors.getSavedPosteriorClass(0);
            probability_absent = posteriors.getAbsentProfileProbability(1);
            super.prepareComputation();
            //System.out.println("#*PD.PCT.pC done.");
        }
        
        @Override
        protected void computeAncestralReconstruction(int family_idx)
        {    
            //System.out.println("#**PD.cAR "+family_idx+"/"+data_file.getData().getFamilyName(family_idx));

            double[][] post = posteriors.getAllPosteriors(family_idx+1); // in posteriors, ABSENT was inserted with index 0 
            reconstructed_family_count[family_idx] = post[Posteriors.PosteriorType.PRESENT.ordinal()];
            reconstructed_family_multi[family_idx] = post[Posteriors.PosteriorType.MULTI.ordinal()];
            reconstructed_family_gain[family_idx] = post[Posteriors.PosteriorType.GAIN.ordinal()];
            reconstructed_family_loss[family_idx] = post[Posteriors.PosteriorType.LOSS.ordinal()];
            reconstructed_family_expansion[family_idx] = post[Posteriors.PosteriorType.EXPANSION.ordinal()];
            reconstructed_family_reduction[family_idx] = post[Posteriors.PosteriorType.REDUCTION.ordinal()];
            
            posterior_family_class[family_idx] = posteriors.getSavedPosteriorClass(family_idx+1);
            //System.out.println("#**PD.cAR "+family_idx+"\tdone");
        }
        
        @Override 
        protected void reportFamilyDone(int family_idx)
        {
            int pct = Math.min(99, 50+(int)((50.0*cumulative_work[family_idx])/total_work));
            //System.out.println("#*PD.rFD "+family_idx+" pct "+pct);
            setProgress(pct);
            firePropertyChange("compute:post", null, families[family_idx]+" ("+Integer.toString(1+family_idx)+"/"+Integer.toString(families.length)+")");
        }
    }
    
    
    /**
     * Exportable interface 
     * 
     * @param f File into which the output is written
     * 
     * @throws java.io.IOException when something goes wrong 
     */
    @Override
    public void saveData(File f) throws IOException
    {
        JFrame top_frame = WorkSpaceCount.getDealerCount(this).getTopFrame();
        ColumnSelectionDialog dialog = new ColumnSelectionDialog(top_frame, f);
        
        Dimension frameD = top_frame.getSize();
        dialog.pack();
        dialog.setBounds(frameD.width/20,frameD.height/20,Math.min(frameD.width*9/10,1000),240);
        dialog.setVisible(true);
        
        // modal dialog -> will come back here when selection is done
        if (dialog.do_save)
        {
            PrintStream PS = new PrintStream(dialog.file);
            PS.println(WorkSpaceCount.getDealerCount(this).getStandardHeader(getClass()));
            PS.println("#| Rate file: "+rate_file.getFile());
            saveColumns(
                    PS,
                    dialog.include_absent.isSelected(),
                    dialog.include_profile.isSelected(), 
                    dialog.include_annotations.isSelected(),
                    dialog.include_class_membership.isSelected(),
                    dialog.include_presence.isSelected(),
                    dialog.include_multi.isSelected(),
                    dialog.include_family_totals.isSelected(),
                    dialog.include_gain.isSelected(),
                    dialog.include_loss.isSelected(),
                    dialog.include_expansion.isSelected(),
                    dialog.include_reduction.isSelected());
            PS.close();
        }

    }    
    
    public void saveColumns(PrintStream PS,
            boolean show_absent,
            boolean save_profile,
            boolean save_annotations,
            boolean save_class_memberships,
            boolean save_node_presence,
            boolean save_node_multi,
            boolean save_totals,
            boolean save_gains,
            boolean save_losses,
            boolean save_expansions,
            boolean save_reductions)
    {
        
        //System.out.println("#*PD.sC absent "+show_absent
        //        +" profile "+save_profile
        //        +" ann "+save_annotations
        //        +" class "+save_class_memberships
        //        +" pres "+save_node_presence
        //        +" multi "+save_node_multi
        //        +" tot "+save_totals
        //        +" gain "+save_gains
        //        +" loss "+save_losses
        //        +" exp "+save_expansions
        //        +" red "+save_reductions);
        
        OccurrenceTable table = data_file.getData();
        RateVariation   rates = rate_file.getData(); 
        
        int num_leaves = main_tree.getNumLeaves();
        int num_nodes = main_tree.getNumNodes();
        int num_edges = main_tree.getNumEdges();
        
        // header line 
        PS.print("# ");
        PS.print("Family\tPattern");
        if (save_annotations)
        {
            int num_props = table.getKnownPropertiesCount();
            for (int prop_idx=1; prop_idx<num_props; prop_idx++) // idx==0 for family name
                PS.print("\t"+table.getPropertyName(prop_idx));
        }
        if (save_class_memberships)
        {
            int num_classes = rates.getNumClasses();
            for (int class_idx=0; class_idx<num_classes; class_idx++)
                if (rates.getClassProbability(class_idx)!=0.0)
                {
                    int cat_dup = rates.getDuplicationRateCategory(class_idx);
                    int cat_loss = rates.getLossRateCategory(class_idx);
                    int cat_trans = rates.getTransferRateCategory(class_idx);
                    int cat_edge = rates.getEdgeLengthCategory(class_idx);
                    String name = "C"+Integer.toString(class_idx)+"/e"+cat_edge+",d"+cat_dup+",l"+cat_loss+",t"+cat_trans;
                    PS.print("\t"+name);
                }
        }
        if (save_totals)
        {
            PS.print("\tArrivals" +
                    "\tGains" +
                    "\tLosses" +
                    "\tExpansions" +
                    "\tContractions");
        }
        
        if (save_profile)
        {
            for (int leaf_idx=0; leaf_idx<num_leaves; leaf_idx++)
                PS.print("\t"+getShortNodeName(leaf_idx));
        }
        if (save_node_presence)
        {
            for (int node_idx=0; node_idx<num_nodes; node_idx++)
                PS.print("\t"+getShortNodeName(node_idx)+":present");
        }
        if (save_node_multi)
        {
            for (int node_idx=0; node_idx<num_nodes; node_idx++)
                PS.print("\t"+getShortNodeName(node_idx)+":multi");
        }
        if (save_gains)
        {
            for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                PS.print("\t"+getShortNodeName(edge_idx)+":gain");
        }
        if (save_losses)
        {
            for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                PS.print("\t"+getShortNodeName(edge_idx)+":loss");
        }
        if (save_expansions)
        {
            for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                PS.print("\t"+getShortNodeName(edge_idx)+":expansion");
        }
        if (save_reductions)
        {
            for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                PS.print("\t"+getShortNodeName(edge_idx)+":contraction");
        }
        PS.println();

        int num_families = families.length;
        
        if (show_absent)
        {
            double exp_missing = probability_absent * num_families/(1.-probability_absent);
            PS.print("ABSENT p0="+probability_absent+"\t"+exp_missing);

            //char[] all0 = new char[num_leaves];
            //java.util.Arrays.fill(all0, '0');
            //PS.print(new String(all0));
            if (save_annotations)
            {
                int num_props = table.getKnownPropertiesCount();
                for (int prop_idx=1; prop_idx<num_props; prop_idx++) // idx==0 for family name
                    PS.print("\tn/a");
            }
            if (save_class_memberships)
            {
                int num_classes = rates.getNumClasses();
                for (int class_idx=0; class_idx<num_classes; class_idx++)
                    if (rates.getClassProbability(class_idx)!=0.0)
                    {
                        PS.print("\t"+posterior_absent_class[class_idx]);
                    }
            }
            
            if (save_totals)
            {
                
                double[] total = new double[reconstructed_family_absent.length];
                for (int type=0; type<total.length; type++)
                {
                    total[type]=0.0;
                    for (int idx=0; idx<reconstructed_family_absent[type].length; idx++)
                        total[type]+=reconstructed_family_absent[type][idx];
                }
                double root_presence = reconstructed_family_absent[Posteriors.PosteriorType.PRESENT.ordinal()][main_tree.getNumNodes()-1];
                double count_arrivals = total[Posteriors.PosteriorType.GAIN.ordinal()]+root_presence;
                PS.print("\t"+count_arrivals);
                PS.print("\t"+total[Posteriors.PosteriorType.GAIN.ordinal()]);
                PS.print("\t"+total[Posteriors.PosteriorType.LOSS.ordinal()]);
                PS.print("\t"+total[Posteriors.PosteriorType.EXPANSION.ordinal()]);
                PS.print("\t"+total[Posteriors.PosteriorType.REDUCTION.ordinal()]);
            }
            if (save_profile)
            {
                for (int leaf_idx=0; leaf_idx<num_leaves; leaf_idx++)
                    PS.print("\t0");
            }
            if (save_node_presence)
            {
                for (int node_idx=0; node_idx<num_nodes; node_idx++)
                    PS.print("\t"+reconstructed_family_absent[Posteriors.PosteriorType.PRESENT.ordinal()][node_idx]);
            }
            if (save_node_multi)
            {
                for (int node_idx=0; node_idx<num_nodes; node_idx++)
                    PS.print("\t"+reconstructed_family_absent[Posteriors.PosteriorType.MULTI.ordinal()][node_idx]);
            }
            if (save_gains)
            {
                for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                    PS.print("\t"+reconstructed_family_absent[Posteriors.PosteriorType.GAIN.ordinal()][edge_idx]);
            }
            if (save_losses)
            {
                for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                    PS.print("\t"+reconstructed_family_absent[Posteriors.PosteriorType.LOSS.ordinal()][edge_idx]);
            }
            if (save_expansions)
            {
                for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                    PS.print("\t"+reconstructed_family_absent[Posteriors.PosteriorType.EXPANSION.ordinal()][edge_idx]);
            }
            if (save_reductions)
            {
                for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                    PS.print("\t"+reconstructed_family_absent[Posteriors.PosteriorType.REDUCTION.ordinal()][edge_idx]);
            }
            PS.println();
            
        }
        
        for (int family_idx=0; family_idx<num_families; family_idx++)
        {
            int[] family_sizes = table.getSizes(family_idx);

            PS.print(families[family_idx]);
            PS.print("\t"+table.getProfile(family_idx).getPatternString());
            if (save_annotations)
            {
                int num_props = table.getKnownPropertiesCount();
                for (int prop_idx=1; prop_idx<num_props; prop_idx++) // idx==0 for family name
                    PS.print("\t"+table.getFamilyProperty(family_idx, prop_idx));
            }
            if (save_class_memberships)
            {
                int num_classes = rates.getNumClasses();
                double[] pc = posterior_family_class[family_idx];
                if (pc==null)
                {
                    pc = new double[num_classes];
                    java.util.Arrays.fill(pc, Double.NaN);
                }
                for (int class_idx=0; class_idx<num_classes; class_idx++)
                    if (rates.getClassProbability(class_idx)!=0.0)
                    {
                        PS.print("\t"+pc[class_idx]);
                    }
            }
            if (save_totals)
            {
                PS.print("\t"+this.getReconstructedFamilyArrivals(family_idx));
                PS.print("\t"+this.getReconstructedFamilyGain(family_idx));
                PS.print("\t"+this.getReconstructedFamilyLoss(family_idx));
                PS.print("\t"+this.getReconstructedFamilyExpansion(family_idx));
                PS.print("\t"+this.getReconstructedFamilyReduction(family_idx));
            }
            if (save_profile)
            {
                for (int leaf_idx=0; leaf_idx<num_leaves; leaf_idx++)
                    PS.print("\t"+family_sizes[leaf_idx]);
            }
            if (save_node_presence)
            {
                for (int node_idx=0; node_idx<num_nodes; node_idx++)
                    PS.print("\t"+getReconstructedFamilyCount(family_idx, node_idx));
            }
            if (save_node_multi)
            {
                for (int node_idx=0; node_idx<num_nodes; node_idx++)
                    PS.print("\t"+this.getReconstructedFamilyMulti(family_idx, node_idx));
            }
            if (save_gains)
            {
                for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                    PS.print("\t"+this.getReconstructedFamilyGain(family_idx, edge_idx));
            }
            if (save_losses)
            {
                for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                    PS.print("\t"+this.getReconstructedFamilyLoss(family_idx, edge_idx));
            }
            if (save_expansions)
            {
                for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                    PS.print("\t"+this.getReconstructedFamilyExpansion(family_idx,edge_idx));
            }
            if (save_reductions)
            {
                for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                    PS.print("\t"+this.getReconstructedFamilyReduction(family_idx, edge_idx));
            }
            PS.println();
        }
        
    }

    protected void savePosteriors(PrintStream out)
    {
        saveColumns(out, true, false, false, true, true, true, false, true, true, true, true);        
    }
    
    /**
     * Sets up the reconstruction variables
     * 
     * @param reader input (format is defined by savePosteriors)
     * @param rate_file here you must set the rate file
     */
    protected void loadPosteriors(Reader reader, DataFile<RateVariation> rate_file) throws IOException
    {
        this.rate_file = rate_file;
        initAllReconstructionVariables();
        
        // how many non-0 classes?
        RateVariation rates = rate_file.getData();
        
        int num_all_categories = rates.getNumClasses();
        double[] prob_category = new double[num_all_categories];
        for (int class_idx=0; class_idx<num_all_categories; class_idx++)
            prob_category[class_idx] = rates.getClassProbability(class_idx);
        
        BufferedReader BR=new BufferedReader(reader);
        String line=null;
        int family_idx=0;
        do 
        {
            line = BR.readLine();
            if (line == null || line.trim().length()==0 || line.startsWith("#"))
                continue;
            String[] fields=StringSplit.splitAtTabs(line);        
            String family_name = fields[0];
            
            // [0] is name
            // []1 is pattern
            // [2..nc] are categories
            double[] posterior_input_class = new double[num_all_categories];  
            int col_idx = 2;
            for (int class_idx=0; class_idx<num_all_categories; class_idx++)
              if (prob_category[class_idx]!=0.0)
              {
                  posterior_input_class[class_idx] = Double.parseDouble(fields[col_idx]);
                  col_idx++;
              }
            // :present columns
            int num_nodes = main_tree.getNumNodes();
            double[] reconstructed_input_count = new double[num_nodes];
            for (int node_idx=0; node_idx<num_nodes; node_idx++)
            {
              reconstructed_input_count[node_idx] = Double.parseDouble(fields[col_idx]);
              col_idx++;
            }
            // :multi columns
            double[] reconstructed_input_multi = new double[num_nodes];
            for (int node_idx=0; node_idx<num_nodes; node_idx++)
            {
              reconstructed_input_multi[node_idx] = Double.parseDouble(fields[col_idx]);
              col_idx++;
            }
            // :gain columns
            int num_edges = main_tree.getNumEdges();
            double[] reconstructed_input_gain = new double[num_nodes];
            for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
            {
              reconstructed_input_gain[edge_idx] = Double.parseDouble(fields[col_idx]);
              col_idx++;
            }
            // :loss columns
            double[] reconstructed_input_loss = new double[num_nodes];
            for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
            {
              reconstructed_input_loss[edge_idx] = Double.parseDouble(fields[col_idx]);
              col_idx++;
            }
            // :expansion columns
            double[] reconstructed_input_expansion = new double[num_nodes];
            for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
            {
              reconstructed_input_expansion[edge_idx] = Double.parseDouble(fields[col_idx]);
              col_idx++;
            }
            // :contraction columns
            double[] reconstructed_input_reduction = new double[num_nodes];
            for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
            {
              reconstructed_input_reduction[edge_idx] = Double.parseDouble(fields[col_idx]);
              col_idx++;
            }            
            
            if (family_name.startsWith("ABSENT"))
            {
                int i = family_name.indexOf("p0=");
                probability_absent = Double.parseDouble(family_name.substring(i+3));
                posterior_absent_class = posterior_input_class;
                reconstructed_family_absent = new double[6][];
                reconstructed_family_absent[Posteriors.PosteriorType.PRESENT.ordinal()] = reconstructed_input_count;
                reconstructed_family_absent[Posteriors.PosteriorType.MULTI.ordinal()] = reconstructed_input_multi;
                reconstructed_family_absent[Posteriors.PosteriorType.GAIN.ordinal()] = reconstructed_input_gain;
                reconstructed_family_absent[Posteriors.PosteriorType.LOSS.ordinal()] = reconstructed_input_loss;
                reconstructed_family_absent[Posteriors.PosteriorType.EXPANSION.ordinal()] = reconstructed_input_expansion;
                reconstructed_family_absent[Posteriors.PosteriorType.REDUCTION.ordinal()] = reconstructed_input_reduction;
            } else
            { // regular family
                posterior_family_class[family_idx] = posterior_input_class;
                reconstructed_family_count[family_idx] = reconstructed_input_count;
                reconstructed_family_multi[family_idx] = reconstructed_input_multi;
                reconstructed_family_gain[family_idx] = reconstructed_input_gain;
                reconstructed_family_loss[family_idx] = reconstructed_input_loss;
                reconstructed_family_expansion[family_idx] = reconstructed_input_expansion;
                reconstructed_family_reduction[family_idx] = reconstructed_input_reduction;

                family_idx++;
            }
            
        } while (line != null);
        BR.close();
        finishComputeAll();
    }
    
    

    
    /**
     * Class used to select the columns that need to be written in a file
     */ 
    private class ColumnSelectionDialog extends JDialog implements ActionListener
    {
        private ColumnSelectionDialog(JFrame mommy, File file)
        {
            super(mommy, "Column selection for exporting posterior reconstruction", true);
            this.file = file;
            this.do_save = false;
            initComponents();
        }
        private JCheckBox include_profile;
        private JCheckBox include_annotations;

        private JCheckBox include_class_membership;
        private JCheckBox include_presence;
        private JCheckBox include_multi;

        private JCheckBox include_family_totals;
        
        private JCheckBox include_gain;
        private JCheckBox include_loss;
        private JCheckBox include_expansion;
        private JCheckBox include_reduction;
        
        private JCheckBox include_absent;
        
        private File file;
        
        private boolean do_save;
                
        private void initComponents()
        {
            JPanel main_panel = new JPanel();
            BoxLayout main_layout = new BoxLayout(main_panel, BoxLayout.PAGE_AXIS);
            main_panel.setLayout(main_layout);
            
            JPanel column_selection_panel = new JPanel();
            GridLayout column_selection_layout = new GridLayout(0,4);
            column_selection_panel.setLayout(column_selection_layout);
            
            JLabel file_info = new JLabel("(exporting into "+file.getName()+")");
            file_info.setFont(new Font("Serif",Font.ITALIC,LookAndFeel.TABLE_FONT_SIZE*6/5));
            JLabel selection_info = new JLabel("   Select the columns to include:");
            selection_info.setFont(new Font("Serif",Font.BOLD,LookAndFeel.TABLE_FONT_SIZE*6/5));
            
            include_absent = new JCheckBox("Include absent families");
            include_absent.setSelected(lineage_totals_include_absent.isSelected());

            column_selection_panel.add(selection_info);
            column_selection_panel.add(file_info);
            column_selection_panel.add(Box.createHorizontalGlue());
            column_selection_panel.add(include_absent);
            
            include_profile = new JCheckBox("Family profile");
            include_profile.setSelected(true);
            
            include_annotations = new JCheckBox("Annotations");
            include_annotations.setSelected(false);
            include_annotations.setEnabled(data_file.getData().getKnownPropertiesCount()>1);
            
            include_class_membership = new JCheckBox("Rate categories");
            int nc = rate_file.getData().getNumClasses();
            for (int i=0; i<rate_file.getData().getNumClasses(); i++)
                if (rate_file.getData().getClassProbability(i)==0.0) nc--;
            include_class_membership.setSelected(nc>1);
            include_class_membership.setEnabled(nc>1);
            
            column_selection_panel.add(include_profile);
            column_selection_panel.add(Box.createHorizontalGlue());
            column_selection_panel.add(include_annotations);
            column_selection_panel.add(include_class_membership);
            
            include_presence = new JCheckBox("Presence at ancestral nodes");
            include_presence.setSelected(true);
            
            include_multi = new JCheckBox("Multiple homologs at ancestors");
            include_multi.setSelected(true);
            
            include_family_totals = new JCheckBox("Totals across lineages");
            include_family_totals.setSelected(true);

            column_selection_panel.add(include_presence);
            column_selection_panel.add(include_multi);
            column_selection_panel.add(Box.createHorizontalGlue());
            column_selection_panel.add(include_family_totals);
            
            include_gain = new JCheckBox("Lineage-specific family gains");
            include_gain.setSelected(true);
            include_loss = new JCheckBox("Lineage-specific famiy losses");
            include_loss.setSelected(true);
            include_expansion = new JCheckBox("Lineage-specific expansions");
            include_expansion.setSelected(true);
            include_reduction = new JCheckBox("Lineage-specific contractions");
            include_reduction.setSelected(true);
            
            column_selection_panel.add(include_gain);
            column_selection_panel.add(include_loss);
            column_selection_panel.add(include_expansion);
            column_selection_panel.add(include_reduction);
            
            main_panel.add(column_selection_panel);
            main_panel.add(Box.createVerticalGlue());
            
            
            JPanel button_panel = new JPanel();
            BoxLayout button_layout = new BoxLayout(button_panel, BoxLayout.LINE_AXIS);
            button_panel.setLayout(button_layout);
            
            JButton cancel_button = new JButton("Cancel");
            cancel_button.setActionCommand("cancel");
            cancel_button.addActionListener(this);
            JButton ok_button = new JButton("Save");
            ok_button.setActionCommand("save");
            ok_button.addActionListener(this);
            
            button_panel.add(Box.createHorizontalGlue());
            button_panel.add(cancel_button);
            button_panel.add(ok_button);
            
            main_panel.add(button_panel);
            
            add(main_panel);
        }
        
        
        
        public void actionPerformed(ActionEvent e)
        {
            String cmd = e.getActionCommand();
            if ("cancel".equals(cmd))
            {
                do_save = false;
                dispose();
            } else if ("save".equals(cmd))
            {
                do_save = true;
                dispose();
            }
        }
    }
    
    /**
     * Includes a checkbox for absent families.
     * 
     * @param bb Box for the bottom bar
     */
    @Override
    protected void createBottomBarLeft(Box bb)
    {
        super.createBottomBarLeft(bb);
        lineage_totals_include_absent = new JCheckBox("Include absent families");
        lineage_totals_include_absent.setFont(new Font("Serif", Font.BOLD, LookAndFeel.TREE_PANEL_FONT_SIZE));
        lineage_totals_include_absent.setSelected(false);
        lineage_totals_include_absent.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e)
            {
                //System.out.println("#*PD.cBBL.aP "+e);
                last_row_selection = null; // so that recomputing is enforced
                recomputeLineageTableForCurrentSelection();            
            }
        });
        bb.add(lineage_totals_include_absent);
        bb.add(Box.createHorizontalGlue());
    }
    
    @Override
    protected ZoomableTreePanel createZoomableTreePanel()
    {
        return new PosteriorZoomablePane();
    }

    protected class PosteriorZoomablePane extends ARZoomablePane
    {
        @Override
        protected void initBottomBarElements(Box bb)
        {
            computation_progress.setPreferredSize(new Dimension(400,computation_progress.getPreferredSize().height));
            super.initBottomBarElements(bb);
        }
        
    }
    
    @Override
    protected FamilyTableModel createFamilyTableModel()
    {
        return new FamilyTableModel();
    }
    
    
    
    protected class FamilyTableModel extends OccurrenceTableModel
    {     
        protected FamilyTableModel()
        {
            super(data_file.getData(), false);
        }

        @Override
        public int getColumnCount()
        {
            int num_nodes = main_tree.getNumNodes();
            int num_edges = main_tree.getNumEdges();
            int num_leaves = main_tree.getNumLeaves();
            int num_usual_columns = super.getColumnCount();
            return num_usual_columns
                    +5 // arrivals, gins, losses, expansions, contractions
                    +2*(num_nodes-num_leaves) //presence & multi
                    +4*num_edges;
        }
        
        @Override
        public Object getValueAt(int row_idx, int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            if (column_idx < num_usual_columns)
                return super.getValueAt(row_idx, column_idx);
            else
                column_idx-=num_usual_columns;

            // sums across lineages
            if (column_idx == 0)
                return new DoubleRoundedForDisplay(getReconstructedFamilyArrivals(row_idx));
            else 
                column_idx--;
            
            if (column_idx == 0)
                return new DoubleRoundedForDisplay(getReconstructedFamilyGain(row_idx));
            else 
                column_idx--;
            
            if (column_idx == 0)
                return new DoubleRoundedForDisplay(getReconstructedFamilyLoss(row_idx));
            else 
                column_idx--;
            
            if (column_idx == 0)
                return new DoubleRoundedForDisplay(getReconstructedFamilyExpansion(row_idx));
            else 
                column_idx--;

            if (column_idx == 0)
                return new DoubleRoundedForDisplay(getReconstructedFamilyReduction(row_idx));
            else column_idx--;
            
            int num_nodes = main_tree.getNumNodes();
            int num_edges = main_tree.getNumEdges();
            int num_leaves = main_tree.getNumLeaves();
            if (column_idx < num_nodes-num_leaves)
                return new DoubleRoundedForDisplay(getReconstructedFamilyCount(row_idx,column_idx+num_leaves));
            else 
                column_idx -= num_nodes-num_leaves;
            
            if (column_idx < num_nodes-num_leaves)
                return new DoubleRoundedForDisplay(getReconstructedFamilyMulti(row_idx,column_idx+num_leaves));
            else 
                column_idx -= num_nodes-num_leaves;
            if (column_idx < num_edges)
                return new DoubleRoundedForDisplay(getReconstructedFamilyGain(row_idx,column_idx));
            else 
                column_idx -= num_edges;
            if (column_idx < num_edges)
                return new DoubleRoundedForDisplay(getReconstructedFamilyLoss(row_idx,column_idx));
            else 
                column_idx -= num_edges;
            if (column_idx < num_edges)
                return new DoubleRoundedForDisplay(getReconstructedFamilyExpansion(row_idx,column_idx));
            else 
                column_idx -= num_edges;
            if (column_idx < num_edges)
                return new DoubleRoundedForDisplay(getReconstructedFamilyReduction(row_idx,column_idx));
            else 
                return "[This column does not exist.]";
        }
        
        @Override
        public String getColumnName(int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            if (column_idx < num_usual_columns)
                return super.getColumnName(column_idx);
            else
                column_idx-=num_usual_columns;
            
            if (column_idx == 0)
                return "Arrivals";
            else 
                column_idx--;
            
            if (column_idx == 0)
                return "Gains";
            else 
                column_idx--;
            
            if (column_idx == 0)
                return "Losses";
            else 
                column_idx--;
            
            if (column_idx == 0)
                return "Expansions";
            else 
                column_idx--;

            if (column_idx == 0)
                return "Contractions";
            else 
                column_idx--;
            
            int num_nodes = main_tree.getNumNodes();
            int num_edges = main_tree.getNumEdges();
            int num_leaves = main_tree.getNumLeaves();
            if (column_idx < num_nodes-num_leaves)
                return getShortNodeName(column_idx+num_leaves);
            else
                column_idx -= num_nodes-num_leaves;
            if (column_idx < num_nodes-num_leaves)
                return getShortNodeName(column_idx+num_leaves)+":m";
            else 
                column_idx -= num_nodes-num_leaves;

            if (column_idx < num_edges)
                return getShortNodeName(column_idx)+":g";
            else 
                column_idx -= num_edges;
            if (column_idx < num_edges)
                return getShortNodeName(column_idx)+":l";
            else 
                column_idx -= num_edges;
            if (column_idx < num_edges)
                return getShortNodeName(column_idx)+"++";
            else 
                column_idx -= num_edges;
            if (column_idx < num_edges)
                return getShortNodeName(column_idx)+"--";
            else 
                return "[This column does not exist.]";
        }
        
        @Override
        public String getColumnHeaderToolTip(int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            if (column_idx < num_usual_columns)
                return super.getColumnHeaderToolTip(column_idx);
            else
                column_idx-=num_usual_columns;
            
            if (column_idx == 0)
                return "Number of times the family appeared (presence at root + lineage-specific gains)";
            else 
                column_idx--;
            
            if (column_idx == 0)
                return "Number of lineages in which the family was gained.";
            else 
                column_idx--;
            
            if (column_idx == 0)
                return "Number of lineages in which the family was lost.";
            else 
                column_idx--;
            
            if (column_idx == 0)
                return "Number of lineages in which the famiy underwent an expansion (from size 1).";
            else 
                column_idx--;

            if (column_idx == 0)
                return "Number of lineages in which the family underwent contractions (to size 1).";
            else 
                column_idx--;
            
            int num_nodes = main_tree.getNumNodes();
            int num_edges = main_tree.getNumEdges();
            int num_leaves = main_tree.getNumLeaves();
            
            if (column_idx < num_nodes-num_leaves)
                return "Probability of at least one homolog present at node "+getLongNodeName(column_idx+num_leaves)+" in the family history";
            else
                column_idx -= num_nodes-num_leaves;
            if (column_idx < num_nodes-num_leaves)
                return "Probability that more than one homologs are present at node "+getLongNodeName(column_idx+num_leaves)+" in the family history";
            else 
                column_idx -= num_nodes-num_leaves;

            String edge_info = getLongNodeName(main_tree.getParentIndex(column_idx % num_edges))+"-"+getLongNodeName(column_idx % num_edges);

            if (column_idx < num_edges)
                return "Probability that the family was acquired in the lineage "+edge_info;
            else 
                column_idx -= num_edges;
            if (column_idx < num_edges)
                return "Probability that the family was completely lost in the lineage "+edge_info;
            else 
                column_idx -= num_edges;
            if (column_idx < num_edges)
                return "Probability that the family expanded (from size 1) in the lineage "+edge_info;
            else 
                column_idx -= num_edges;
            if (column_idx < num_edges)
                return "Probability that the family contracted (to size 1) in the lineage "+edge_info;
            else 
                return "[This column does not exist.]";            
        }
        
        @Override 
        public String getCellToolTip(int row_idx, int column_idx)
        {
            Object value = getValueAt(row_idx, column_idx);

            int num_usual_columns = super.getColumnCount();
            if (column_idx < num_usual_columns)
                return super.getCellToolTip(row_idx,column_idx);
            else
                column_idx-=num_usual_columns;
            
            int num_leaves = main_tree.getNumLeaves();
            int num_nodes = main_tree.getNumNodes();
            int num_edges = main_tree.getNumEdges();
            
            DoubleRoundedForDisplay dvalue = (DoubleRoundedForDisplay)value;
            if (Double.isNaN(dvalue.doubleValue()))
                return "not available (yet)";
            
            String family_info = "Family "+families[row_idx];
 
            if (column_idx == 0)
                return family_info+" newly appeared (presence at root + lineage-specific gains) an expected number of "+dvalue.doubleValue()+" times";
            else 
                column_idx--;
            
            if (column_idx == 0)
                return family_info +" was newly gained (lineage-specific gains) an expected number of "+dvalue.doubleValue()+" lineages";
            else 
                column_idx--;
            
            if (column_idx == 0)
                return family_info+" was completely lost in an expected number of "+dvalue.doubleValue()+" lineages";
            else 
                column_idx--;
            
            if (column_idx == 0)
                return family_info +" underwent an expansion (from size 1) in an expected number of "+dvalue.doubleValue()+" lineages";
            else 
                column_idx--;

            if (column_idx == 0)
                return family_info +" underwent contractions (to size 1) in an expected number of "+dvalue.doubleValue()+" lineages";
            else 
                column_idx--;
            
            
            if (column_idx < num_nodes-num_leaves)
                return family_info+" had at least one member present at node "+getLongNodeName(column_idx+num_leaves)+" with probability "+dvalue.doubleValue();
            else
                column_idx -= num_nodes-num_leaves;
            if (column_idx < num_nodes-num_leaves)
                return family_info+" had more than one member present at node "+getLongNodeName(column_idx+num_leaves)+" with probability "+dvalue.doubleValue();
            else 
                column_idx -= num_nodes-num_leaves;

            String edge_info = getLongNodeName(main_tree.getParentIndex(column_idx % num_edges))+"-"+getLongNodeName(column_idx % num_edges);
            
            if (column_idx < num_edges)
                return family_info+" was newly acquired in the lineage "+edge_info+" with probability "+dvalue.doubleValue();
            else 
                column_idx -= num_edges;
            if (column_idx < num_edges)
                return family_info +" was completely lost in the lineage "+edge_info+" with probability "+dvalue.doubleValue();
            else 
                column_idx -= num_edges;
            if (column_idx < num_edges)
                return family_info +" expanded (from size 1) in the lineage "+edge_info+" with probability "+dvalue.doubleValue();
            else 
                column_idx -= num_edges;
            if (column_idx < num_edges)
                return family_info +" contracted (to size 1) in the lineage "+edge_info+" with probability "+dvalue.doubleValue();
            else 
                return "[This column does not exist.]";            
        }
        
        
        @Override
        public Class getColumnClass(int column_idx)
        {
            int num_usual_columns = super.getColumnCount();
            if (column_idx < num_usual_columns)
                return super.getColumnClass(column_idx);
            else
                return DoubleRoundedForDisplay.class;
        }
        
    }

}
