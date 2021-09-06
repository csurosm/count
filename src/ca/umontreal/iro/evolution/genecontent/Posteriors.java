package ca.umontreal.iro.evolution.genecontent;

import java.util.Hashtable;

import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.Verbose;

import ca.umontreal.iro.matek.DiscreteDistribution;

/**
 *
 * This class is used to compute posterior probabilities for presence
 * at ancestral nodes and events on the edges.
 * The computation is optimized for sequential queries
 * about the profiles.
 *
 * @author csuros
 */
public class Posteriors extends BasicExecutable
{
    private static boolean LINEAGE_TOTALS_ONLY=false;
    
    /**
     * Initializes the class for a data set and rate model
     * 
     * @param table data set of phylogenetic profiles
     * @param model rate variation model used in computing the posteriors 
     */
    public Posteriors(OccurrenceTable table, RateVariation model)
    {
        this(table, model, 100);
    }
    public Posteriors(OccurrenceTable table, RateVariation model, int segment_size)
    {
        setSegmentSize(segment_size);
        init(table, model);
    }

    /**
     * Used by convenience methods: data structures will be initialized later
     */
    private Posteriors(){setSegmentSize(100);}
    
    private OccurrenceTable main_table;
    private RateVariation model;
    
    private StableComputation SC;
    
    private OccurrenceTable segment_table[];
    
    //private int family_segment_idx[];
    
    private int segment_size;
    private void setSegmentSize(int s)
    {
        this.segment_size = s;
    }
    
    private static final String ABSENT_PROFILE_NAME = "ABSENT";
    
    /**
     * Initializes the data structures.
     * Also adds a profile with all-0 sizes.
     * 
     * @param current_table occurrence current_table for phylogenetic profiles
     * @param model rate variation model
     */
    private void init(OccurrenceTable table, RateVariation model)
    {
//        System.out.println("#*P.i model ==========\n"+model.tableRates());
//        for (int c=0; c<model.getNumClasses(); c++)
//            if (model.getClassProbability(c)!=0.0)
//            {
//                TreeWithRates class_tree = model.getRateTree(c);
//                System.out.println("#*P.i ct "+c+"\t"+class_tree+" [g "+class_tree.hasGain()+" d "+class_tree.hasDuplication()+"]");
//            }
        this.main_table = table;
        this.model = model;
        {
            int num_leaves = model.getMainTree().getNumLeaves();
            int [] all_zero_profile = new int[num_leaves];
            main_table.addProfile(all_zero_profile, ABSENT_PROFILE_NAME);
        }
        
        int num_families = table.getNumFamilies();
        int num_family_segments = num_families/segment_size;
        if (num_family_segments*segment_size<num_families)
            num_family_segments++;
        segment_table = new OccurrenceTable[num_family_segments];
//        System.out.println("#*P.init segment "+num_family_segments);
        for (int segment_idx=0; segment_idx<num_family_segments; segment_idx++)
        {
            int start_idx = segment_idx * segment_size;
            int end_idx = Math.min(start_idx+segment_size ,num_families);
            segment_table[segment_idx] = table.tableForFamilies(start_idx, end_idx-1);
        }
        active_segment_idx = -1;
        
        initSavedPosteriorClass();
    }
    
    private int active_segment_idx=-99;
    private int active_subtree_segment_idx = -99;
    
    private double[][][] upper_likelihood0;
    private double[][][] upper_likelihood1;
    private double[][][] side_likelihood0;
    private double[][][] side_likelihood1;

    private void initSegment(int family_segment_idx, boolean need_subtree_likelihoods)
    {
        if (active_segment_idx == family_segment_idx
                && (!need_subtree_likelihoods || active_subtree_segment_idx == family_segment_idx))
            return;
        //System.out.println("#*P.iS "+family_segment_idx);
        Verbose.message("P.iS "+family_segment_idx);
        
        OccurrenceTable current_table = segment_table[family_segment_idx];
        SC = new StableComputation(current_table, model);
        if (need_subtree_likelihoods)
        {
            TreeWithRates main_tree = model.getMainTree();
            int num_nodes = main_tree.getNumNodes();
            Hashtable <String, Integer> main_node_indices = new Hashtable<String,Integer>();
            for (int i=0; i<main_tree.getNumLeaves(); i++)
            {
                NodeWithRates leaf = main_tree.getNode(i);
                main_node_indices.put(leaf.getName(), new Integer(i));
            }
            //SC0 = new StableComputation[num_nodes];
            //SC1 = new StableComputation[num_nodes];
            //SC_edge = new StableComputation[num_nodes];

            int num_families = current_table.getNumFamilies();
            int num_classes = model.getNumClasses();

            upper_likelihood0 = new double[num_nodes][num_families][num_classes];
            upper_likelihood1 = new double[num_nodes][num_families][num_classes];
            side_likelihood0  = new double[num_nodes][num_families][num_classes];
            side_likelihood1  = new double[num_nodes][num_families][num_classes];

            for (int node_idx=0; node_idx<num_nodes; node_idx++)
            {
                NodeWithRates N = main_tree.getNode(node_idx);
                if (!N.isLeaf() && !N.isRoot())
                {
                    // compute pruned tree
                    NodeWithRates pruned_root = NodeWithRates.copyTree(main_tree.getRoot());
                    TreeWithRates pruned_tree = new TreeWithRates(pruned_root);
                    NodeWithRates pruned_node = pruned_tree.getNode(node_idx);
                    pruned_node.clearChildren(); // nodes are in the ame order
                    pruned_tree = new TreeWithRates(pruned_root);
                    RateVariation pruned_tree_rates = model.sameModelForDifferentTree(pruned_tree);
                    // compute occurrence tables for this guy
                    NodeWithRates[] pruned_leaves = pruned_tree.getLeaves();
                    int[] leaf_index_map = new int[pruned_leaves.length];
                    int pruned_node_idx = -1;
                    for (int leaf_idx=0; leaf_idx<pruned_leaves.length; leaf_idx++)
                    {
                        NodeWithRates leaf = pruned_leaves[leaf_idx];
                        if (leaf == pruned_node)
                        {
                            leaf_index_map[leaf_idx]=-1;
                            pruned_node_idx = leaf_idx;
                        }
                        else
                            leaf_index_map[leaf_idx] = main_node_indices.get(leaf.getName()).intValue();
                    }
                    { // set T0
                        OccurrenceTable T0 = new OccurrenceTable(pruned_leaves);
                        int[][] pruned_values=new int[num_families][pruned_leaves.length];
                        for (int family_idx=0; family_idx<num_families; family_idx++)
                        {
                            PhyleticProfile PP = current_table.getProfile(family_idx);
                            for (int leaf_idx=0; leaf_idx<pruned_leaves.length; leaf_idx++)
                            {
                                int original_idx = leaf_index_map[leaf_idx];
                                if (original_idx==-1)
                                    pruned_values[family_idx][leaf_idx] = 0;
                                else
                                    pruned_values[family_idx][leaf_idx]=PP.get(original_idx);
                            }
                        }
                        T0.setTable(pruned_values);

                        StableComputation SC0 = new StableComputation(T0, pruned_tree_rates);
                        for (int pridx=0; pridx<num_families; pridx++)
                            for (int cidx=0; cidx<num_classes; cidx++)
                                if (model.getClassProbability(cidx)!=0.0)
                                    upper_likelihood0[node_idx][pridx][cidx] = SC0.getLikelihood(pridx, cidx);
                    } 
                    {
                        // set T1
                        OccurrenceTable T1 = new OccurrenceTable(pruned_leaves);
                        int[][] pruned_values=new int[num_families][pruned_leaves.length];
                        for (int family_idx=0; family_idx<num_families; family_idx++)
                        {
                            PhyleticProfile PP = current_table.getProfile(family_idx);
                            for (int leaf_idx=0; leaf_idx<pruned_leaves.length; leaf_idx++)
                            {
                                int original_idx = leaf_index_map[leaf_idx];
                                if (original_idx==-1)
                                    pruned_values[family_idx][leaf_idx] = 1;
                                else
                                    pruned_values[family_idx][leaf_idx]=PP.get(original_idx);
                            }
                        }
                        T1.setTable(pruned_values);
                        StableComputation SC1 = new StableComputation(T1, pruned_tree_rates);
                        for (int pridx=0; pridx<num_families; pridx++)
                            for (int cidx=0; cidx<num_classes; cidx++)
                                if (model.getClassProbability(cidx)!=0.0)
                                    upper_likelihood1[node_idx][pridx][cidx] = SC1.getLikelihood(pridx, cidx);
                    }
                    //SC0[node_idx] = null;
                    //SC1[node_idx] = null;
                    //SC0[node_idx] = new StableComputation(T0,pruned_tree_rates);
                    //SC1[node_idx] = new StableComputation(T1,pruned_tree_rates);
                } // for each inner node
                if (!N.isRoot())
                {
                    int parent_idx = main_tree.getParentIndex(node_idx);
                    NodeWithRates pruned_root = NodeWithRates.copyTree(main_tree.getRoot());
                    TreeWithRates pruned_tree = new TreeWithRates(pruned_root);
                    NodeWithRates parent = pruned_tree.getNode(parent_idx);

                    int nc = parent.getNumChildren();
                    parent.clearChildren();
                    for (int ci=0; ci<nc; ci++)
                    {
                        int child_idx = pruned_tree.getChildIndex(parent_idx, ci);
                        if (child_idx != node_idx)
                            parent.addChild(pruned_tree.getNode(child_idx));
                    }

                    parent.clearParent();
                    pruned_tree = new TreeWithRates(parent);

                    RateVariation pruned_tree_rates = model.sameModelForDifferentTree(pruned_tree);

                    NodeWithRates[] pruned_leaves = pruned_tree.getLeaves();
                    OccurrenceTable T = new OccurrenceTable(pruned_leaves);
                    int[] leaf_index_map = new int[pruned_leaves.length];
                    for (int leaf_idx=0; leaf_idx<leaf_index_map.length; leaf_idx++)
                    {
                        NodeWithRates leaf = pruned_leaves[leaf_idx];
                        leaf_index_map[leaf_idx] = main_node_indices.get(leaf.getName()).intValue();
                    }
                    {
                        int[][] pruned_values = new int[num_families][pruned_leaves.length];
                        for (int family_idx=0; family_idx<num_families; family_idx++)
                        {
                            PhyleticProfile PP = current_table.getProfile(family_idx);
                            for (int leaf_idx=0; leaf_idx<pruned_leaves.length; leaf_idx++)
                            {
                                int original_idx = leaf_index_map[leaf_idx];
                                pruned_values[family_idx][leaf_idx]=PP.get(original_idx);
                            }
                        }
                        T.setTable(pruned_values);
                        //System.out.println("# ------------------------ pruned "+node_idx);
                        //System.out.println(T.getFormattedTable());
                    }
                    //SC_edge[node_idx] = null;
                    StableComputation SC_edge = new StableComputation(T,pruned_tree_rates);
                    for (int pridx=0; pridx<num_families; pridx++)
                        for (int cidx=0; cidx<num_classes; cidx++)
                            if (model.getClassProbability(cidx)!=0.0)
                            {
                                double[] lik = SC_edge.getRootLikelihoods(pridx, cidx);
                                int pruned_root_idx = SC_edge.getModel().getMainTree().getNumNodes()-1;
                                double D = SC_edge.getExtinctionProbability(pruned_root_idx, cidx);
                                side_likelihood0[node_idx][pridx][cidx] = lik[0];
                                side_likelihood1[node_idx][pridx][cidx] = lik[0]*D;
                                if (lik.length>1)
                                    side_likelihood1[node_idx][pridx][cidx] += lik[1]*(1.-D);
                            }

                    //Verbose.message("P.ini SC_edge "+node_idx+"/"+N.newickName()+"\t"+SC_edge[node_idx].getModel().getMainTree().getRoot().newickTree());
                }
            } // for each node
            active_subtree_segment_idx = family_segment_idx;           
        } // if upper_ and side_likelihoods need to be set
        else
        {
            upper_likelihood0 = upper_likelihood1 = side_likelihood0 = side_likelihood1 = null;
            active_subtree_segment_idx = -99;
        }
        active_segment_idx = family_segment_idx;
    }
    
    private int getSegmentProfileIndex(int family_idx, boolean need_subtree_likelihoods)
    {
        int segment_idx = family_idx / segment_size;
        int pidx = family_idx - segment_idx * segment_size;
        initSegment(segment_idx, need_subtree_likelihoods);
        return pidx;
    }
    
    /**
     * Computes the posterior distribution for rate variation categories
     * 
     * @param family_idx family index in the occurrence current_table
     * @return array of posterior probabilities for each rate variation class 
     */
    public double[] getPosteriorClassDistribution(int family_idx)
    {
        int pidx = getSegmentProfileIndex(family_idx, false);
        int num_classes = model.getNumClasses();
        double[] p_class = new double[num_classes];
        double L_total = 0.0;
        for (int cidx=0; cidx<num_classes; cidx++)
        {
            double pc = model.getClassProbability(cidx);
            if (pc==0.0)
                p_class[cidx]=0.0;
            else
            {
                double L = SC.getLikelihood(pidx,cidx);
                p_class[cidx] = pc*L;
                //Verbose.message("P.gPCD "+pidx+"\tclass "+cidx+"\t"+L+"\tpc "+pc);
                L_total += p_class[cidx]; 
            }
        }
        for (int cidx=0;cidx<num_classes; cidx++)
            p_class[cidx] /= L_total;
        return p_class;
    }
    
    /**
     * Posterior probabilites for classes are saved here
     */
    private double[][] saved_posterior_class;

    /**
     * Sets up the saved_posterior_class array
     */
    protected void initSavedPosteriorClass()
    {
        int num_families = main_table.getNumFamilies();
        saved_posterior_class = new double[num_families][];
        for (int pidx=num_families-1; pidx>=0; pidx--)
            setSavedPosteriorClass(pidx);
    }
    
    /**
     * Called from within initSavedPosteriorClass, first when
     * the posterior class probabilities need to be computed. 
     * 
     * @param pidx family index (goes down to 0 when called from the iteration)
     */
    protected void setSavedPosteriorClass(int pidx)
    {
        saved_posterior_class[pidx]=getPosteriorClassDistribution(pidx);
    }
    
    public double[] getSavedPosteriorClass(int pidx)
    {
        return saved_posterior_class[pidx];
    }
    
    
    /**
     * Computes the posterior distribution for presence at a node. 
     * 
     * @param node_idx index of the node in the tree
     * @param family_idx index of the family
     * @return 2-element array for posterior probabilties of 0 and 1 family aizes 
     */
    public double[] getPosteriorCountDistribution(int node_idx, int family_idx)
    {
//        System.out.println("#*P.gPCD "+node_idx+"\tfam "+family_idx);
        if (saved_posterior_class==null)
            initSavedPosteriorClass();
        
        int pridx = getSegmentProfileIndex(family_idx, true);

        double[] prob = new double[2];

        TreeWithRates main_tree = model.getMainTree();
        NodeWithRates N = main_tree.getNode(node_idx);
        if (N.isLeaf())
        {
            int cnt = main_table.getProfile(family_idx).get(node_idx);
            if (cnt==0)
                prob[0]=1.0;
            else if (cnt==1)
                prob[1]=1.0;
        } else 
        {
            //double L_total = SC.getLikelihood(pidx);
            int num_classes = model.getNumClasses();
            double[] u0; 
            double[] u1; 
            if (N.isRoot())
            {
                u0 = new double[num_classes];
                u1 = new double[num_classes];
                for (int cidx=0; cidx<num_classes; cidx++)
                    if (model.getClassProbability(cidx)!=0.0)
                    {
                        DiscreteDistribution rootD = model.getRootPrior();
                        double[] p = rootD.getDistribution(1);
                        u0[cidx]=p[0];
                        u1[cidx]=p[1];
                    }
            } else
            {
                u0 = upper_likelihood0[node_idx][pridx];
                u1 = upper_likelihood1[node_idx][pridx];
                //for (int cidx=0; cidx<num_classes; cidx++)
                //{
                //    if (model.getClassProbability(cidx)!=0.0)
                //    {
                //        u0[cidx] = SC0[node_idx].getLikelihood(pridx, cidx);
                //        u1[cidx] = SC1[node_idx].getLikelihood(pridx, cidx);
                //    }
                //        
                //}
            }
            double p0=0.0;
            double p1=0.0;
            for (int cidx=0; cidx<num_classes; cidx++)
                if (model.getClassProbability(cidx)!=0.0)
                {
                    double pc = saved_posterior_class[family_idx][cidx];
                    if (pc != 0.0)
                    {
                        double[] p = SC.getNodeLikelihoods(node_idx, pridx, cidx);
                        double p0c = u0[cidx]*p[0];
                        double D = SC.getExtinctionProbability(node_idx, cidx);
                        double lower1 = D*p[0];
                        if (p.length>1)
                            lower1 += (1.-D)*p[1];
                        double p1c = u1[cidx]*lower1;
                        double Lc = SC.getLikelihood(pridx, cidx);
                        p0c /= Lc;
                        p1c /= Lc;
                        p0 += p0c*pc;
                        p1 += p1c*pc;
//                        System.out.println("#*P.gPCD "+family_idx+"\t"+node_idx+"/"+N.newickName()
//                                +"\tc "+cidx+"\tp0c "+p0c+"\tp1c "+p1c+"\tpc "+pc+"\tL "+Lc
//                                +"\t u0 "+u0[cidx]+", u1 "+u1[cidx]
//                                +"\t l0 "+p[0]+", l1 "+lower1
//                                +"\t"+(p0c+p1c)
//                                +(p0c+p1c>1.0?"\t// %%%%%%%%%%%%%%%%%% "+node_idx:"")
//                                +"\tp[]="+java.util.Arrays.toString(p));
                    }
                }
            prob[0]=p0;
            prob[1]=p1;
        }
        return prob;
    }
    
    public enum EventType {STAYS0, GAIN, LOSS, CONSERVED, XGD};
    
    /**
     * Computes the posterior distribution for events on an edge. 
     * 
     * @param node_idx index of the node in the tree (child node for the edge)
     * @param family_idx family index
     * @return 5-element array for posterior probabilities of events in the order of EventType
     */
    public double[] getPosteriorEventDistribution(int node_idx, int family_idx)
    {
        if (saved_posterior_class==null)
            initSavedPosteriorClass();
        
        int pridx = getSegmentProfileIndex(family_idx, true);

        TreeWithRates main_tree = model.getMainTree();
        NodeWithRates N = main_tree.getNode(node_idx);

        int num_classes = model.getNumClasses();
        double[] u0 = new double[num_classes];
        double[] u1 = new double[num_classes];
        
        int parent_idx = main_tree.getParentIndex(node_idx);

        if (main_tree.getNode(parent_idx).isRoot())
        {
            for (int cidx=0; cidx<num_classes; cidx++)
                if (model.getClassProbability(cidx)!=0.0)
                {
                    DiscreteDistribution rootD = model.getRootPrior();
                    double[] p = rootD.getDistribution(1);
                    u0[cidx]=p[0];
                    u1[cidx]=p[1];
                }
        } else
        {
            u0 = upper_likelihood0[parent_idx][pridx];
            u1 = upper_likelihood1[parent_idx][pridx];
            //for (int cidx=0; cidx<num_classes; cidx++)
            //{
            //    if (model.getClassProbability(cidx)!=0.0)
            //    {
            //        u0[cidx] = SC0[parent_idx].getLikelihood(pridx, cidx);
            //        u1[cidx] = SC1[parent_idx].getLikelihood(pridx, cidx);
            //    }
            //}
        }
        
        double[] s0 = side_likelihood0[node_idx][pridx]; 
        double[] s1 = side_likelihood1[node_idx][pridx];
        
        //for (int cidx=0; cidx<num_classes; cidx++)
        //{
        //    if (model.getClassProbability(cidx)!=0.0)
        //    {
        //        double[] lik = SC_edge[node_idx].getRootLikelihoods(pridx, cidx);
        //        int pruned_root_idx = SC_edge[node_idx].getModel().getMainTree().getNumNodes()-1;
        //        double D = SC_edge[node_idx].getExtinctionProbability(pruned_root_idx, cidx);
        //        
        //        s0[cidx] = lik[0];
        //        s1[cidx] = lik[0]*D;
        //        if (lik.length>1)
        //            s1[cidx]+=lik[1]*(1.-D);
                
                
                //if (false && Verbose.isVerbose())
                //{
                //    Verbose.message("P.go s\t"+node_idx+"/"+N.newickName()+"\tclass "+cidx+"\ts0 "+lik[0]+"\ts1 "+s1[cidx]+"\tD "+D+"\t//"+(lik.length>1?lik[1]:0.0));
                //    TreeWithRates Et = SC_edge[node_idx].getModel().getMainTree();
                //    Verbose.message("P.go s\t"+node_idx+"/"+N.newickName()+"\tclass "+cidx+"\ttree "+Et.getRoot().newickTree());
                //    for (int ni=0; ni<Et.getNumNodes(); ni++)
                //    {
                //        double[] nl = SC_edge[node_idx].getNodeLikelihoods(ni, pridx, cidx);
                //        Verbose.message("P.go s\t"+node_idx+"/"+N.newickName()+"\tclass "+cidx+"\tnode "+ni+"\tl0 "+nl[0]+"\tl1 "+(nl.length>1?nl[1]:-1.)+"\t// "+Et.getNode(ni));
                //    }
                //    
                //}
        //   }
        //}
        
        double[] l0 = new double[num_classes];
        double[] l1 = new double[num_classes];
        for (int cidx=0; cidx<num_classes; cidx++)
            if (model.getClassProbability(cidx)!=0.0)
            {
                double[] lik = SC.getNodeLikelihoods(node_idx, pridx, cidx);
                l0[cidx] = lik[0];
                double D = SC.getExtinctionProbability(node_idx, cidx);
                l1[cidx] = D*lik[0];
                if (lik.length>1)
                    l1[cidx]+=(1.-D)*lik[1];
                Verbose.message("P.gPED "+node_idx+"/"+N.newickName()+"\tclass "+cidx+"\tlik0 "+lik[0]+"\tl1 "+l1[cidx]+"\tD "+D);
            }
        
        double[] edge_events = new double[EventType.values().length];
        for (int cidx=0; cidx<num_classes; cidx++)
            if (model.getClassProbability(cidx)!=0.0)
            {
                double pc = saved_posterior_class[family_idx][cidx];
                if (pc != 0.0)
                {
                    double Lc = SC.getLikelihood(pridx,cidx);
                    TreeWithRates class_tree = model.getRateTree(cidx);

                    NodeWithRates E = class_tree.getNode(node_idx);
                    double[] h = E.getTransferDistribution().getDistribution(2);
                    double[] g = E.getDuplicationDistribution().getDistribution(2);
                    double p00 = u0[cidx]*s0[cidx]*h[0]*l0[cidx]/Lc;
                    double p01 = u0[cidx]*s0[cidx]*h[1]*l1[cidx]/Lc;
                    double p10 = u1[cidx]*s1[cidx]*h[0]*g[0]*l0[cidx]/Lc;
                    double p11 = u1[cidx]*s1[cidx]*h[0]*g[1]*l1[cidx]/Lc; // xenologous gene displacement counted separately
                    double pxgd = u1[cidx]*s1[cidx]*h[1]*g[0]*l1[cidx]/Lc;

                    edge_events[EventType.STAYS0.ordinal()]    += pc*p00;
                    edge_events[EventType.GAIN.ordinal()]      += pc*p01;
                    edge_events[EventType.LOSS.ordinal()]      += pc*p10;
                    edge_events[EventType.CONSERVED.ordinal()] += pc*p11;
                    edge_events[EventType.XGD.ordinal()]       += pc*pxgd;

    //                System.out.println("#*P.gPED "+family_idx+"/"+pridx+"\t"+node_idx+"/"+N.newickName()+"\tclass "+cidx
    //                        +"\t-- "+p00
    //                        +"\tgain  "+p01
    //                        +"\tloss "+p10
    //                        +"\t++ "+p11
    //                        +"\txgd "+pxgd
    //                        +"\tsum "+(p00+p01+p10+p11+pxgd)
    //                        +"\t// h "+h[0]+", "+h[1]+", "+h[2]
    //                        +"\t g "+g[0]+", "+g[1]+", "+g[2]
    //                        +"\t u "+u0[cidx]+", "+u1[cidx]
    //                        +"\t l "+l0[cidx]+", "+l1[cidx]
    //                        +"\t s "+s0[cidx]+", "+s1[cidx]
    //                        +"\t Lc "+Lc
    //                        +"\t pc "+pc
    //                        +"\tct "+class_tree+" [g "+class_tree.hasGain()+" d "+class_tree.hasDuplication()
    //                        +"\tE "+E);
                }
            }
        
        return edge_events;
    }
    
    public enum PosteriorType {PRESENT, MULTI, GAIN, LOSS, EXPANSION, REDUCTION};
    
    /**
     * Computes all the posteriors probabilities: presence, multi-member-presence, 
     * gain, loss, expansion, reduction. 
     * 
     * @param family_idx index of the family
     * @return array of expected numbers for events: [e][n] has the expectation for event 
     *     <var>e</var> (in the order of PosteriorType) at node <var>n</var>. 
     */
    public double[][] getAllPosteriors(int family_idx)
    {
//        System.out.println("#*P.gAP "+family_idx);
        TreeWithRates main_tree = model.getMainTree();
        int num_nodes = main_tree.getNumNodes();
        double[][] retval = new double[PosteriorType.values().length][num_nodes];
        double[][] post = new double[num_nodes][];
        double[][] edge = new double[num_nodes-1][];
        
        double family_occurrence_count = 1.0;
        if (ABSENT_PROFILE_NAME.equals(main_table.getFamilyName(family_idx)))
        {
            double absent_profile_prob1 = SC.getAbsentProfileProbability(1);
            double absent_profile_prob = SC.getAbsentProfileProbability(MIN_PRESENT_LINEAGES);
            int num_families = main_table.getNumFamilies()-1;
            double exp_missing = absent_profile_prob * num_families/(1.-absent_profile_prob);
            family_occurrence_count = exp_missing;
        }

        for (int node_idx=0; node_idx<num_nodes; node_idx++)
        {
//            System.out.println("#*P.gAP "+family_idx+"\t: "+node_idx+"/"+main_tree.getNode(node_idx).newickName()+"\tpost/edge (num_nodes="+num_nodes+")");
            post[node_idx] = getPosteriorCountDistribution(node_idx,family_idx);
//            System.out.println("#**P.gAP "+family_idx+"\t: "+node_idx+"/"+main_tree.getNode(node_idx).newickName()+"\tpost done.");
            
            if (!main_tree.getNode(node_idx).isRoot())
            {
                edge[node_idx]=getPosteriorEventDistribution(node_idx,family_idx);
//                System.out.println("#**P.gAP "+family_idx+"\t: "+node_idx+"/"+main_tree.getNode(node_idx).newickName()+"\tedge done.");
            }
        }
//        System.out.println("#**P.gAP "+family_idx+"\t done with post/edge");
        
        for (int node_idx=0; node_idx<num_nodes; node_idx++)
        {
            double [] p = post[node_idx];
            //System.out.print("\t"+p[0]);
            retval[PosteriorType.PRESENT.ordinal()][node_idx] = 1.-p[0];
            retval[PosteriorType.MULTI.ordinal()][node_idx] = 1.-p[0]-p[1];

            if (!main_tree.getNode(node_idx).isRoot())
            {
                double[] e = edge[node_idx];
                int parent_idx = main_tree.getParentIndex(node_idx);
                double gain = post[parent_idx][0]-e[EventType.STAYS0.ordinal()];
                if (Math.abs(gain)<1e-15) gain=0.0;
                double loss = p[0]-e[EventType.STAYS0.ordinal()];
                if (Math.abs(loss)<1e-15) loss=0.0;
                retval[PosteriorType.GAIN.ordinal()][node_idx] = gain;
                retval[PosteriorType.LOSS.ordinal()][node_idx] = loss;
//                System.out.println("#**P.gAP "+family_idx+"\t: "+node_idx+"/"+main_tree.getNode(node_idx).newickName()+"\tp0 "+p[0]+"\tp1 "+p[1]+"\te "+e[EventType.CONSERVED.ordinal()]+"\tloss "+loss);

                double expansion = post[parent_idx][1]-e[EventType.CONSERVED.ordinal()]-e[EventType.LOSS.ordinal()]-e[EventType.XGD.ordinal()];
                if (Math.abs(expansion)<1e-15) expansion = 0.0;
                double reduction = p[1]-e[EventType.CONSERVED.ordinal()]-e[EventType.GAIN.ordinal()]-e[EventType.XGD.ordinal()];
                if (Math.abs(reduction)<1e-15) reduction = 0.0;
                
                retval[PosteriorType.EXPANSION.ordinal()][node_idx] = expansion;
                retval[PosteriorType.REDUCTION.ordinal()][node_idx] = reduction;
            }        
        }
//        System.out.println("#**P.gAP "+family_idx+"\t done with retval");

        return retval;
        
    }
    
    public double getAbsentProfileProbability(int min_lineages_present)
    {
        return SC.getAbsentProfileProbability(min_lineages_present);
    }
            
    /**
     * Called at the command-line invokation.
     * 
     * @param args arguments to main, after switches are stripped (tree, current_table, and rates)
     * @throws java.lang.Exception if input files do not have the right format
     */
    private void go(String[] args) throws Exception
    {
        if (args.length != 3)
        {
            System.err.println("Call as java "+getClass().getCanonicalName()+" tree table rates");
            System.exit(2008);
        }
        reportLaunch(args);
        reportOtherArguments("Max. paralogs "+MAX_PARALOGS+", min. lineages "+MIN_PRESENT_LINEAGES+"\ttotals "+LINEAGE_TOTALS_ONLY);
        
        String tree_file = args[0];
        String table_file = args[1];
        String rates_file = args[2];

        TreeWithRates input_tree = new TreeWithRates(ca.umontreal.iro.evolution.Parser.readNewick(new java.io.FileReader(tree_file)));
        OccurrenceTable input_table = new OccurrenceTable(input_tree.getLeaves());
        input_table.readTable(new java.io.FileReader(table_file));
        OccurrenceTable filtered_table = input_table.filterByMaximumSize(MAX_PARALOGS,MIN_PRESENT_LINEAGES);
        int num_families = filtered_table.getNumFamilies();

        RateVariation input_model = RateVariation.read(new java.io.FileReader(rates_file),input_tree);

        init(filtered_table, input_model);
        
        NodeWithRates[] nodes = input_tree.getDFT();
        System.out.print("Family\tPattern");
        for (int class_idx=0; class_idx<model.getNumClasses(); class_idx++)
            if (model.getClassProbability(class_idx)!=0.0)
            {
                int cat_dup = model.getDuplicationRateCategory(class_idx);
                int cat_loss = model.getLossRateCategory(class_idx);
                int cat_trans = model.getTransferRateCategory(class_idx);
                int cat_edge = model.getEdgeLengthCategory(class_idx);
                String name = "C"+Integer.toString(class_idx)+"/e"+cat_edge+",d"+cat_dup+",l"+cat_loss+",t"+cat_trans;
                System.out.print("\t"+name);
            }
        for (int node_idx=0; node_idx<nodes.length; node_idx++)
        {
            String name = nodes[node_idx].newickName();
            //System.out.print("\t"+name+":0");
            System.out.print("\t"+name+":1");
            System.out.print("\t"+name+":m");
            if (!nodes[node_idx].isRoot())
            {
                System.out.print("\t"+name+":gain");
                System.out.print("\t"+name+":loss");
                System.out.print("\t"+name+":expansion");
                System.out.print("\t"+name+":reduction");
            }
        }
        System.out.println();

        double[] total_present_1 = new double[nodes.length];
        double[] total_present_m = new double[nodes.length];
        double[] total_gain = new double[nodes.length];
        double[] total_loss = new double[nodes.length];
        double[] total_expansion = new double[nodes.length];
        double[] total_reduction = new double[nodes.length];
        
        double[] total_class = new double[model.getNumClasses()];
        
        double absent_profile_prob1 = SC.getAbsentProfileProbability(1);
        double absent_profile_prob = SC.getAbsentProfileProbability(MIN_PRESENT_LINEAGES);
        double exp_missing = absent_profile_prob * num_families/(1.-absent_profile_prob);

        for (int pidx=0; pidx<filtered_table.getNumFamilies(); pidx++)
        {
            double[][] post = new double[nodes.length][];
            double[][] edge = new double[nodes.length-1][];
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                post[node_idx] = getPosteriorCountDistribution(node_idx,pidx);
                if (!input_tree.getNode(node_idx).isRoot())
                {
                    edge[node_idx]=getPosteriorEventDistribution(node_idx,pidx);
                    double[] e = edge[node_idx];
                    Verbose.message("P.go "+pidx+"\t"+node_idx+"/"+input_tree.getNode(node_idx).newickName()
                            +"\tg "+e[EventType.GAIN.ordinal()]
                            +"\tl "+e[EventType.LOSS.ordinal()]
                            +"\tc "+e[EventType.CONSERVED.ordinal()]
                            +"\tz "+e[EventType.STAYS0.ordinal()]
                            +"\tx "+e[EventType.XGD.ordinal()]);
                }
            }
            
            String family_name = filtered_table.getFamilyName(pidx);
            double family_occurrence_count = (ABSENT_PROFILE_NAME.equals(family_name)?exp_missing:1.0);
            
            if (!LINEAGE_TOTALS_ONLY)
            {
                System.out.print(family_name);
                System.out.print("\t"+filtered_table.getProfile(pidx).getPatternString());
            }
            double[] pc = saved_posterior_class[pidx];
            for (int cidx=0; cidx<model.getNumClasses(); cidx++)
                if (model.getClassProbability(cidx)!=0.0)
                {
                    total_class[cidx] += family_occurrence_count * pc[cidx];
                    if (!LINEAGE_TOTALS_ONLY)
                        System.out.print("\t"+pc[cidx]);
                }
            
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                double [] p = post[node_idx];
                //System.out.print("\t"+p[0]);
                if (!LINEAGE_TOTALS_ONLY)
                {
                    System.out.print("\t"+p[1]);
                    System.out.print("\t"+(1.-p[0]-p[1]));
                }

                total_present_1[node_idx] += family_occurrence_count*p[1];
                total_present_m[node_idx] += family_occurrence_count*(1.-p[0]-p[1]);
                
                if (!nodes[node_idx].isRoot())
                {
                    double[] e = edge[node_idx];
                    int parent_idx = input_tree.getParentIndex(node_idx);
                    double gain = post[parent_idx][0]-e[EventType.STAYS0.ordinal()];
                    if (Math.abs(gain)<1e-15) gain=0.0;
                    double loss = p[0]-e[EventType.STAYS0.ordinal()];
                    if (Math.abs(loss)<1e-15) loss=0.0;
                    if (!LINEAGE_TOTALS_ONLY)
                    {
                        System.out.print("\t"+gain);
                        System.out.print("\t"+loss);
                    }
                    total_gain[node_idx] += family_occurrence_count*gain;
                    total_loss[node_idx] += family_occurrence_count*loss;
                            
                    double expansion = post[parent_idx][1]-e[EventType.CONSERVED.ordinal()]-e[EventType.LOSS.ordinal()]-e[EventType.XGD.ordinal()];
                    if (Math.abs(expansion)<1e-15) expansion = 0.0;
                    double reduction = p[1]-e[EventType.CONSERVED.ordinal()]-e[EventType.GAIN.ordinal()]-e[EventType.XGD.ordinal()];
                    if (Math.abs(reduction)<1e-15) reduction = 0.0;
                    if (!LINEAGE_TOTALS_ONLY)
                    {
                        System.out.print("\t"+expansion);
                        System.out.print("\t"+reduction);
                    }
                    total_expansion[node_idx] += family_occurrence_count*expansion;
                    total_reduction[node_idx] += family_occurrence_count*reduction;
                }
            }
            
            if (!LINEAGE_TOTALS_ONLY)
                System.out.println();
        } // for all profiles
        
        if (LINEAGE_TOTALS_ONLY)
        {
            System.out.print("TOTALS");
            System.out.print("\t"+exp_missing);
            for (int cidx=0; cidx<model.getNumClasses(); cidx++)
                if (model.getClassProbability(cidx)!=0.0)
                {
                    System.out.print("\t"+total_class[cidx]);
                }
            
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                System.out.print("\t"+total_present_1[node_idx]);
                System.out.print("\t"+total_present_m[node_idx]);
                if (!nodes[node_idx].isRoot())
                {
                    System.out.print("\t"+total_gain[node_idx]);
                    System.out.print("\t"+total_loss[node_idx]);
                    System.out.print("\t"+total_expansion[node_idx]);
                    System.out.print("\t"+total_reduction[node_idx]);
                }
            }
            System.out.println();
        } else
        {
            System.out.println("#ABSENT Profiles not in the data set (presence in 0 or 1 lineage)");
            System.out.println("#ABSENT\tall-0\t"+absent_profile_prob1);
            if (MIN_PRESENT_LINEAGES>1)
            {
                double[] single_lineage_probs = SC.getSingleLineageProbability();
                for (int leaf_idx=0; leaf_idx < input_tree.getNumLeaves(); leaf_idx++)
                {
                    System.out.println("#ABSENT\t"+leaf_idx+"/"+nodes[leaf_idx].newickName()+"\t"+single_lineage_probs[leaf_idx]);
                }
            }
            System.out.println("#ABSENT expected missing "+exp_missing);
        }
        
        
    }
    
    private static int MAX_PARALOGS = 150;
    private static int MIN_PRESENT_LINEAGES = 1;
    
    /**
     * Outputs the posterior probabilities in a current_table format to standard output
     * 
     * @param args Command-line arguments
     */
    public static void main(String[] args)
    {
        Verbose.setVerbose(false);
                
        Posteriors O = new Posteriors(); //InclusionExclusionComputation());
        
        int num_switches = 0;
        try {
            while (args.length>2*num_switches && args[2*num_switches].startsWith("-"))
            {
                String arg_switch = args[2*num_switches].substring(1);
                if (arg_switch.equals("h"))
                    O.go(new String[0]); // will throw an Exception
                if (args.length==2*num_switches+1)
                    throw new IllegalArgumentException("Missing argument for switch "+args[2*num_switches]);
                String arg_value = args[2*num_switches+1];
                if (arg_switch.equals("v"))
                {
                    Verbose.setVerbose(arg_value.equals("true"));
                } else if (arg_switch.equals("max_paralogs"))
                {
                    MAX_PARALOGS = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("min_lineages"))
                {
                    MIN_PRESENT_LINEAGES = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("lineage_totals"))
                {
                    LINEAGE_TOTALS_ONLY = "true".equals(arg_value);
                }
                else 
                    throw new IllegalArgumentException("Switch not recognized: '"+args[2*num_switches]+"'");
                    
                num_switches++;
            }
            
            String[] rest=new String[args.length-2*num_switches];
            for (int j=0; j<rest.length; j++)
                rest[j]=args[2*num_switches+j];
            O.go(rest);
        } catch (Exception E)
        {
            die(E);
        }        
    }

}
