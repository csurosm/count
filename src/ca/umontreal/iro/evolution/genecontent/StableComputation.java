
package ca.umontreal.iro.evolution.genecontent;

import java.util.Arrays;

/**
 * Stable likelihood computations.
 *
 * @author csuros
 */

import java.util.Hashtable;

import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.Verbose;
import ca.umontreal.iro.banality.Functions;


import ca.umontreal.iro.matek.DiscreteDistribution;
import ca.umontreal.iro.matek.NegativeBinomial;
import ca.umontreal.iro.matek.PointDistribution;
import ca.umontreal.iro.matek.Poisson;
import ca.umontreal.iro.matek.ShiftedGeometric;

public class StableComputation extends BasicExecutable implements LikelihoodComputation
{
	private static final boolean DEBUG22 = true;
	
    /**
     * Instantiation with allocations for the data structures
     *
     * @param table family size table
     * @param model rate variation model
     */
    public StableComputation(OccurrenceTable table, RateVariation model)
    {
        init(table, model);
    }

    private static final int MAX_PRECOMPUTED_TRANSITIONS = 80;

    /**
     * Instantiation with no memory allocations.
     * Use {@link #initWithoutAllocations(ca.umontreal.iro.evolution.genecontent.OccurrenceTable, ca.umontreal.iro.evolution.genecontent.RateVariation) initWithoutAllocations}
     * to specify the table and the rate model, and
     * allocate the data structures with {@link #allocateDataStructures() allocateDataStructures} before
     * doing the computations.
     */
    public StableComputation(){}
    
    private OccurrenceTable table;
    private RateVariation model;

    @Override
    public OccurrenceTable getTable(){return table;}

    @Override
    public RateVariation getModel(){return model;}
    
//    private long timeP=0L; // timing computations within .Probabilities
//    private long timeLC=0L; // timing computations within .LowerConditionals

    @Override
    public void init(OccurrenceTable table, RateVariation model)
    {
        //Verbose.setVerbose(true);
        initWithoutAllocations(table, model);
        allocateDataStructures();
    }

    /**
     * Registers the family size table and the rate model, but allocates no
     * memory for the data structures.
     *
     * @param table family size table
     * @param model rate model
     */
    public void initWithoutAllocations(OccurrenceTable table, RateVariation model)
    {
        this.table = table;
        this.model = model;
        this.profiles = table.getProfiles();
        this.main_tree = model.getMainTree();
//        System.out.println("#*SC.iWA table "+table.getClass().getName()+" with "+table.getNumFamilies()+" families");
    }

    /**
     * Allocates memory and sets up the data structures.
     */
    public void allocateDataStructures()
    {
        //Verbose.message("SC.i "+Integer.toString(profiles.length)+" families");
        initClassStructures();
        //Verbose.message("SC.i classes done.");
        
        { // initialize strage for memoization of n-choose-m
            int Mmax = 0;
            int root_idx = main_tree.getNumNodes()-1;

            for (int profile_idx=0; profile_idx<profiles.length; profile_idx++)
            {
                int m = family_sizes[profile_idx][root_idx];
                if (m>Mmax)
                    Mmax = m;
            }
            n_choose_m = new double[Mmax+1][];
            for (int i=0; i<=Mmax; i++)
                n_choose_m[i] = new double[i+1];
        }
        
        initClassProbabilities();
        initProfileParameters();
    }
    
    
    private double[][] n_choose_m;
    
    /**
     * Computes n choose m using a memoized implementation (computing Pascal's triangle) 
     * 
     * @param n first argument
     * @param m second argument
     * @return n choose m: in how many ways can you pick m balls from n without replacement
     */
    private double getNchooseM(int n, int m)
    {
        //Verbose.message("SC.gNcM "+n+"\t"+m);
        if (n_choose_m[n][m]==0)
        {
            if (m==n || m==0)
                n_choose_m[n][m]=1.0;
            else if (m==1 || m==n-1)
                n_choose_m[n][m]=(double)n;
            else if (m>n-m)
                n_choose_m[n][m] = getNchooseM(n,n-m);
            else 
            {
                double a = getNchooseM(n-1,m-1);
                double b = getNchooseM(n-1,m);
                n_choose_m[n][m] = a+b;
            }
        }
        return n_choose_m[n][m];
    }

    private PhyleticProfile[] profiles;
    
    private TreeWithRates main_tree;

    /**
     * Computes various types of
     * family size aggregates
     * for subtrees and subgraphs.
     */
    private void initClassStructures()
    {
        computeFamilySizes();
        computeBestChildOrdering();
        computeCombinedSizes();
        computeExtremeSizes();
    }
    
    /** 
     * Family size sums at each subtree.
     * Indexing: family_sizes[family_index][node_index]
     */
    private int[][] family_sizes;
    
    private int[] max_family_sizes;

    /**
     * sets up the {@link #family_sizes family_sizes} array
     */
    private void computeFamilySizes()
    {
        family_sizes = new int[profiles.length][];
        max_family_sizes = new int[main_tree.getNumNodes()];
        
        for (int family_idx = 0; family_idx<profiles.length; family_idx++)
        {
            int[] this_size = profiles[family_idx].computeSubtreeSizes(main_tree);
            family_sizes[family_idx] = this_size;
            for (int node_idx=0; node_idx<main_tree.getNumNodes(); node_idx++)
                if (this_size[node_idx]>max_family_sizes[node_idx])
                    max_family_sizes[node_idx]=this_size[node_idx];
        }
        
            
    }

    /**
     * Gives an optimal ordering of children in the computations.
     * Indexing: best_child_ordering[node_index][processing_step_index]
     */
    private int[][] best_child_ordering;
    /**
     * Inverse of best_child_ordering.
     * Indexing: best_child_ordering[node_index][child_position]
     */
    private int[][] best_child_ordering_inverse;
    
    /**
     * Gives the optimal ordering of child node. 
     * 
     * @param node_idx parent node
     * @param step_index processing order: 0,1,...,numchildren-1
     * @return which child should be processed at the i-th step
     */
    public int getChildForStep(int node_idx, int step_index)
    {
        return best_child_ordering[node_idx][step_index];
    }
    
    /**
     * Inverse of getChildForStep(). 
     * 
     * @param node_idx parent node
     * @param child_idx position of the child among the parent's children: 0,1,...,numchildren-1
     * @return in which step this child should be processed child should be processed at the i-th step
     */
    public int getStepForChild(int node_idx, int child_idx)
    {
        return best_child_ordering_inverse[node_idx][child_idx];
    }
    
    /**
     * calculates how much child order matters in calculating the 
     * likelihood with the computationally stable method
     *
     *
     * For child sizes n1 n2 n3 ..., the computation time is proportional to 
     * n1 n2^2+(n1+n2) n3^2+(n1+n2+n3)n4^2 ...
     */
    private void computeBestChildOrdering()
    {
        NodeWithRates[] tree_nodes = main_tree.getDFT();
        
        best_child_ordering = new int[tree_nodes.length][];
        best_child_ordering_inverse = new int[tree_nodes.length][];
        
//        long total_best = 0L; // used in debug messages
//        long total_worst = 0L; // used in debug messages
//        long total_avg = 0L; // used in debug messages
//
        for (int node_idx=0; node_idx<tree_nodes.length; node_idx++)
        {
            NodeWithRates N = tree_nodes[node_idx];
            int nc = N.getNumChildren();
//            if (true || nc <=1) // let's not complicate the matter too much: different child orderings do not give significantly different running times
//            {
                best_child_ordering[node_idx] = new int[nc];
                best_child_ordering_inverse[node_idx] = new int[nc];
                for (int i=0; i<nc; i++)
                {
                    int j = best_child_ordering[node_idx][i]=i;
                    best_child_ordering_inverse[node_idx][j]=i;
                }
//            } else
//            {
//
//                int num_permutations = (int)Functions.factorial(nc);
//                Verbose.message("SC.cBCO ---------------------------- node "+node_idx+"/"+N.newickName()+"\t"+nc+" children,\t"+num_permutations+" permutations");
//                int [] permutation_scores = new int[num_permutations];
//                PermutationEnumeration PE=new PermutationEnumeration(nc);
//                int[] child_order=new int[nc];
//                for (int family_idx=0; family_idx<profiles.length; family_idx++)
//                {
//                    PE.reset();
//                    for (int permutation_idx=0; PE.hasMoreElements(); permutation_idx++)
//                    {
//                        int this_score = 0;
//                        PE.nextElement(child_order);
//                        int first_child_idx = main_tree.getChildIndex(node_idx,child_order[0]);
//                        int subsum = family_sizes[family_idx][first_child_idx];
//                        for (int ci=1; ci<nc; ci++)
//                        {
//                            int child_idx = main_tree.getChildIndex(node_idx,child_order[ci]);
//                            int s = family_sizes[family_idx][child_idx];
//                            this_score += (subsum+1) * (s+1);// * (s+1);
//                            subsum += s;
//                        }
//                        permutation_scores[permutation_idx] += this_score;
//                    }
//                }
//                int best_permutation_idx = -1;
//                int worst_permutation_idx = -1; // used in debug messages
//                long sum_all = 0L; // used in debug messages
//                for (int pi=0; pi<num_permutations; pi++)
//                {
//                    if (best_permutation_idx==-1 || permutation_scores[pi]<permutation_scores[best_permutation_idx])
//                    {
//                        best_permutation_idx = pi;
//                        best_child_ordering[node_idx] = new int[nc];
//                        for (int i=0; i<nc; i++)
//                            best_child_ordering[node_idx][i]=child_order[i];
//                    }
//                    if (worst_permutation_idx==-1 || permutation_scores[pi]>permutation_scores[worst_permutation_idx])
//                        worst_permutation_idx = pi;
//                    sum_all += permutation_scores[pi];
//                }
//                // statistics for debug messages
//                long avg_score=sum_all / num_permutations;
//                double ratio_bw = ((int)(permutation_scores[best_permutation_idx]/((double)permutation_scores[worst_permutation_idx])*100+0.5))/100.0;
//                double ratio_ba = ((int)(permutation_scores[best_permutation_idx]/((double)avg_score)*100+0.5))/100.0;
//                total_best += permutation_scores[best_permutation_idx];
//                total_worst += permutation_scores[worst_permutation_idx];
//                total_avg += (int)avg_score;
//
//                Verbose.message("SC.cBCO scores best "+permutation_scores[best_permutation_idx]
//                        +"\tb/w "+ratio_bw
//                        +"\tb/avg "+ratio_ba
//                        +"\tworst "+permutation_scores[worst_permutation_idx]
//                        +"\tavg "+avg_score);
//                if (Verbose.isVerbose())
//                {
//                    PE.reset();
//                    for (int permutation_idx=0; PE.hasMoreElements(); permutation_idx++)
//                    {
//                        PE.nextElement(child_order);
//                        StringBuffer msg = new StringBuffer();
//                        msg.append("SC.cBCO permutation "+(1+permutation_idx));
//                        if (permutation_idx==best_permutation_idx)
//                            msg.append('+');
//                        if (permutation_idx==worst_permutation_idx)
//                            msg.append('-');
//                        msg.append("\t");
//                        for (int ci=0; ci<nc; ci++)
//                        {
//                            msg.append(Integer.toString(child_order[ci]));
//                        }
//                        msg.append("\t"+permutation_scores[permutation_idx]);
//                        if (permutation_idx==best_permutation_idx)
//                        {
//                            msg.append("\t");
//                            for (int ci=0; ci<nc; ci++)
//                            {
//                                msg.append(' ');
//                                msg.append(main_tree.getNode(main_tree.getChildIndex(node_idx, child_order[ci])).newickName());
//                            }
//                        }
//                        Verbose.message(msg.toString());
//                    }
//                }
//
//                // make sure the inverse is also set up
//
//                best_child_ordering_inverse[node_idx] = new int[nc];
//                for (int j=0; j<nc; j++)
//                {
//                    int idx = best_child_ordering[node_idx][j];
//                    best_child_ordering_inverse[node_idx][idx] = j;
//                }
//            } // if at least 2 children
        } // for node
//        //if (Verbose.isVerbose())
//        //{
//        //    double ratio_bw = ((int)(total_best/((double)total_worst)*100+0.5))/100.0;
//        //    double ratio_ba = ((int)(total_best/((double)total_avg)*100+0.5))/100.0;
//        //    Verbose.message("SC.cBCO ---------------------------- TOTAL");
//        //    Verbose.message("SC.cBCO scores best "+total_best+"\tb/w "+ratio_bw+"\tb/avg "+ratio_ba+"\tworst "+total_worst+"\tavg "+total_avg);
//        //}
    }
    
    /**
     * Sum of family sizes below a node and its left siblings: <q>left</q> is defined in 
     * terms of best child ordering.
     * Indexing: combined_sizes[family_index][node_index].
     */
    private int[][] combined_sizes;

    /**
     * Maximum value of combined_sizes at each node.
     * Indexing: max_combined_sizes[node_idx]; not computed for root.
     */
    private int[] max_combined_sizes;
    
    /**
     * Sets up the combined_sizes[] array
     */
    private void computeCombinedSizes()
    {
        NodeWithRates[] tree_nodes = main_tree.getDFT();
        combined_sizes = new int[profiles.length][tree_nodes.length];
        max_combined_sizes = new int[tree_nodes.length-1];
        for (int family_idx=0; family_idx<profiles.length; family_idx++)
        {
            for (int node_idx=0; node_idx<tree_nodes.length; node_idx++)
            {
                NodeWithRates N = tree_nodes[node_idx];
                if (N.isRoot()) // before we forget to deal with the root
                    combined_sizes[family_idx][node_idx] = family_sizes[family_idx][node_idx];
                if (!N.isLeaf())
                {
                    int num_children = N.getNumChildren();
                    int cumulative = 0;
                    for (int ci=0; ci<num_children; ci++) // take N's children in computed order
                    {
                        int child_node_idx = main_tree.getChildIndex(node_idx, getChildForStep(node_idx,ci));
                        cumulative += family_sizes[family_idx][child_node_idx];
                        combined_sizes[family_idx][child_node_idx]=cumulative;
                        if (cumulative>max_combined_sizes[child_node_idx])
                            max_combined_sizes[child_node_idx]=cumulative;
                    }
                } // for leaves
            } // for each node
        } // for each family
    }
    
    /**
     * Extreme sizes at every node. For each node u, extreme_sizes[u] is an array with entries at 
     * 0...m where m is the maximum combined size at its left sibling. If u has no left sibling, 
     * m is set at 0. extreme_sizes[u][i] gives the maximum combined size j at node u 
     * when its left sibling has combined size i or less. 
     */
    private int extreme_sizes[][];
    
    /**
     * Sets up the extreme_sizes array.
     */
    private void computeExtremeSizes()
    {
        NodeWithRates[] tree_nodes = main_tree.getDFT();
        extreme_sizes = new int[tree_nodes.length-1][];
        for (int parent_idx=0; parent_idx<tree_nodes.length; parent_idx++)
        {
            NodeWithRates P = main_tree.getNode(parent_idx);
            if (!P.isLeaf())
            {
                for (int ci=0; ci<P.getNumChildren(); ci++)
                {
                    int node_idx = main_tree.getChildIndex(parent_idx, getChildForStep(parent_idx,ci));
                    if (ci==0)
                    {
                        extreme_sizes[node_idx]=new int[1];
                        extreme_sizes[node_idx][0] = max_combined_sizes[node_idx];
                    } else
                    {
                        int left_sibling_idx = main_tree.getChildIndex(parent_idx, getChildForStep(parent_idx,ci-1));
                        extreme_sizes[node_idx]=new int[max_combined_sizes[left_sibling_idx]+1];
                        for (int family_idx=0; family_idx<profiles.length; family_idx++)
                        {
                            int size_at_node = combined_sizes[family_idx][node_idx];
                            int size_at_left_sibling = combined_sizes[family_idx][left_sibling_idx];
                            if (extreme_sizes[node_idx][size_at_left_sibling]<size_at_node)
                                extreme_sizes[node_idx][size_at_left_sibling]=size_at_node;
                        } // for each family
                    }
                    
                    // make sure it's monotone
                    for (int lsize=extreme_sizes[node_idx].length-2; lsize>=0; lsize--)
                    {
                        int max_above = extreme_sizes[node_idx][lsize+1];
                        if (max_above>extreme_sizes[node_idx][lsize])
                            extreme_sizes[node_idx][lsize]=max_above;
                    }
                    
                    //if (Verbose.isVerbose())
                    //{
                    //    Verbose.message("SC.cES node "+node_idx+"/"+main_tree.getNode(node_idx).newickName()+"\tM1="+(extreme_sizes[node_idx].length-1));
                    //    for (int t=0; t<extreme_sizes[node_idx].length; t++)
                    //        Verbose.message("SC.cES "+node_idx+"\t\tt "+t+"\tmax "+extreme_sizes[node_idx][t]);
                    //}
                } // every child
                
            } // not a leaf
        } // for each node
    }
    
    private Probabilities[] class_probabilities;

    /**
     * Initializes data structures for computing
     * probabilities defined by the rate model
     * only (and not the profiles).
     */
    private void initClassProbabilities()
    {
        int num_classes = model.getNumClasses();
        class_probabilities = new Probabilities[num_classes]; 
        for (int cidx=0; cidx<num_classes; cidx++)
        {
            TreeWithRates rate_tree = model.getRateTree(cidx);
//            System.out.println("#*SC.iCP class "+cidx
//                    +"\tlen* "+model.getEdgeLengthCategories()[model.getEdgeLengthCategory(cidx)]
//                    +"\tloss* "+model.getLossRateCategories()[model.getLossRateCategory(cidx)]
//                    +"\tdup* "+model.getDuplicationRateCategories()[model.getDuplicationRateCategory(cidx)]
//                    +"\tgain* "+model.getTransferRateCategories()[model.getTransferRateCategory(cidx)]
//                    );
//            for (int node_idx=0; node_idx<rate_tree.getNumEdges(); node_idx++)
//            {
//                NodeWithRates N = rate_tree.getNode(node_idx);
//                NodeWithRates N0 = main_tree.getNode(node_idx);
//                System.out.println("#*SC.iCP node "+N+"\tN0 "+N0);
//            }
            class_probabilities[cidx]=new Probabilities(rate_tree);
        }
    }
    
    private LowerConditionals[][] profile_parameters;
    private int[] first_profile_occurrences;
    

    /**
     * Initializes the data structures supporting
     * conditional likelihoods for each profile
     * (profile_parameters and first_profile_occurrences).
     */
    private void initProfileParameters()
    {
        
        initProfileOccurrences();
        initLowerconditionals();
    }

    private void initLowerconditionals()
    {
        int num_classes = model.getNumClasses();
        int num_profiles = profiles.length;

        profile_parameters = new LowerConditionals[num_profiles][num_classes];
        for (int pidx=0; pidx<profiles.length; pidx++)
        {
            int first_occurrence = first_profile_occurrences[pidx];
            {
                if (first_occurrence == pidx)
                    for (int cidx=0; cidx<num_classes; cidx++)
                        profile_parameters[pidx][cidx] = newLowerConditionals(profiles[pidx],class_probabilities[cidx]);
                else
                        profile_parameters[pidx] = profile_parameters[first_occurrence];
            }
            setNumProfilesAllocated(pidx+1);
        }
    }

    protected LowerConditionals newLowerConditionals(PhyleticProfile profile, Probabilities probs)
    {
        return new LowerConditionals(profile,probs);
    }
    
    private void initProfileOccurrences()
    {
        int num_profiles = profiles.length;
        
        first_profile_occurrences = new int[num_profiles];
        {
            Hashtable<String,Integer> first_pattern_occurrence = new Hashtable<String,Integer>();
            int num_unique_profiles=0;
            for (int pidx=0;pidx<num_profiles; pidx++)
            {
                String pattern = profiles[pidx].getPatternString();
                int pos = pidx;
                if (first_pattern_occurrence.containsKey(pattern))
                {
                    pos = first_pattern_occurrence.get(pattern).intValue();
                    //Verbose.message("SC.iPP repeated profile "+pattern+"\t"+pidx+"\t"+pos);
                }
                else
                {
                    first_pattern_occurrence.put(pattern, new Integer(pos));
                    num_unique_profiles++;
                }
                first_profile_occurrences[pidx]=pos;
            } // pidx
            //Verbose.message("SC.iPP unique profiles "+num_unique_profiles);
        }
    }


    private int num_profiles_allocated=0;
    public synchronized int getNumProfilesAllocated()
    {
        return num_profiles_allocated;
    }

    private synchronized void setNumProfilesAllocated(int x)
    {
        num_profiles_allocated = x;
    }
    
    public void recomputeCategorySupport(int cidx)
    {
        if (model.getClassProbability(cidx)!=0.0)
            class_probabilities[cidx].recomputeAll();
    }
    
    public void recomputeCategorySupport(int cidx, int node_idx)
    {
        if (model.getClassProbability(cidx)!=0.0)
           class_probabilities[cidx].recomputeOnPathToRoot(node_idx);
    }
    
    public DiscreteDistribution getRootPrior(int cidx)
    {
        return class_probabilities[cidx].survivingRootPrior(model.getRootPrior());
    }
    
    public double[] getRootLikelihoods(int pidx, int cidx)
    {
        int root_idx  = main_tree.getNumNodes()-1; // last node in postorder traversal
        return getNodeLikelihoods(root_idx,pidx,cidx);
    }
    
    public double[] getNodeLikelihoods(int node_idx, int pidx, int cidx)
    {
        LowerConditionals LC = profile_parameters[pidx][cidx];
        double[] L = LC.getLikelihoods(node_idx);
        return L;
    }
    
    
    public void recomputeProfileSupport(int pidx, int cidx, int node_idx)
    {
        if (first_profile_occurrences[pidx]==pidx)
        {
            if (model.getClassProbability(cidx)!=0.0)
            {
                if (!main_tree.getNode(node_idx).isRoot())
                {
                    //Verbose.message("SC.RPS profile "+pidx+"/"+profiles[pidx].getPatternString()+"\tnode "+node_idx+"\tclass "+cidx+" ---------------------------------------");
                    profile_parameters[pidx][cidx].recomputeOnPathToRoot(main_tree.getParentIndex(node_idx));
                }
            }
        }
    }
    
    public void recomputeProfileSupport(int pidx,int cidx)
    {
        if (first_profile_occurrences[pidx]==pidx)
        {
            if (model.getClassProbability(cidx)!=0.0)
            {
                //Verbose.message("SC.RPS profile "+pidx+"/"+profiles[pidx].getPatternString()+"\tclass "+cidx+" ---------------------------------------\n\n\n\n\n");
                profile_parameters[pidx][cidx].recomputeAll();
            }
        }
    }
    
    public double[] getSingleLineageProbability()
    {
        int num_classes = model.getNumClasses();
        int num_leaves = main_tree.getNumLeaves(); 
        double[] p1 = new double[num_leaves];
        for (int cidx=0; cidx<num_classes; cidx++)
        {
            if (model.getClassProbability(cidx)!=0.)
            {
                DiscreteDistribution class_root_prior = getRootPrior(cidx);
                double[] p1c = class_probabilities[cidx].getSingleLineagePresence(class_root_prior);
                double pc = model.getClassProbability(cidx);
                for (int leaf_idx=0; leaf_idx<num_leaves; leaf_idx++)
                    p1[leaf_idx]+=pc*p1c[leaf_idx];
            }
        }
        return p1;
    }
    
//    public double[] getAbsentProfileClassPosteriors(int min_lineages)
//    {
//        if (min_lineages == 0)
//            return null;
//        int num_classes = model.getNumClasses();
//        if (min_lineages == 1)
//        {
//
//        }
//
//        return null;
//    }
    
    public double getAbsentProfileProbability(int min_lineages)
    {
        if (min_lineages == 0)
            return 0.0;
        if (min_lineages>2)
            throw new IllegalArgumentException("min_lineages must be 0,1 or 2");
        int num_classes = model.getNumClasses();
        double p0 = 0.0;
        int root_idx = main_tree.getNumNodes()-1;
        DiscreteDistribution main_prior = model.getRootPrior();
        
        if (min_lineages == 1)
        {
            for (int cidx=0; cidx<num_classes; cidx++)
            {
                if (model.getClassProbability(cidx)!=0.)
                {
                    double p_class = getRootPrior(cidx).getDistribution(0)[0];
                    double z = class_probabilities[cidx].getAllAbsentProbability(); //LC0.getLikelihood(root_idx,0);
                    p0 += model.getClassProbability(cidx)*p_class*z;
                }
            }
        
        } else if (min_lineages == 2)
        {
            for (int cidx=0; cidx<num_classes; cidx++)
            {
                if (model.getClassProbability(cidx)!=0.0)
                {
                    double z = class_probabilities[cidx].getAbsentOrSingle(main_prior);
                    p0 += z*model.getClassProbability(cidx);
                    //Verbose.message("SC.gAPP "+cidx+" z "+z+"\tpc "+model.getClassProbability(cidx)+"\tp0 "+p0);
                }
            }
        }
        
        return p0;
    }
    
    /**
     * Computes the likelihood of a profile. (No correction is employed for mising profiles.)
     *  
     * @param pidx profile index
     * @return probability of the profile
     */
    public double getLikelihood(int pidx)
    {
        int num_classes = model.getNumClasses();
        double p = 0.0;
        for (int cidx=0;cidx<num_classes;cidx++)
            if (model.getClassProbability(cidx)!=0.0)
            {
                double pc = getLikelihood(pidx,cidx);
                double prior = model.getClassProbability(cidx);
                System.out.println("#*SC.gL "+pidx+"\tc "+cidx+"\tpc "+pc+"\tpri "+prior+"\tp "+p);
                p += pc * prior;
            }
        return p;
    }
    
    public double getLikelihood(int pidx, int cidx)
    {
        double[] L = getRootLikelihoods(pidx,cidx);
        DiscreteDistribution R = getRootPrior(cidx);
        double[] pr = R.getDistribution(L.length-1);
        double pc = 0.0;
        for (int m=0; m<L.length; m++)
        {
//            System.out.println("#*SC.gL "+pidx+"\tc "+cidx+"\tm "+m+"\tpr "+pr[m]+"\tL "+L[m]+"\tpc "+pc);
            pc += pr[m]*L[m];
        }
//        System.out.println("#*SC.gL "+pidx+"\tc "+cidx+"\tL[] "+L.length+"\tpc "+pc+"\tL[]="+java.util.Arrays.toString(L));
        return pc;
    }
    
    public double getExtinctionProbability(int node_idx, int cidx)
    {
        return class_probabilities[cidx].getExtinction(node_idx);
    }
        
    /**
     * Class for computing and recomputing survival probabilities.
     */
    protected class Probabilities
    {
        /**
         * Instantiates the class. 
         * 
         * @param rate_tree a tree with the same topology as the one used for the enclosing ConditionedSurvival instance: here, we use the rates too 
         */
        private Probabilities(TreeWithRates rate_tree)
        {
            this.rate_tree = rate_tree;
            initDataStructures();
            recomputeAll();
        }
        
        /**
         * the rates used in this instance
         */
        private TreeWithRates rate_tree;
        
        public TreeWithRates getTree()
        {
            return rate_tree;
        }
        
        /**
         * This is probability G(0) for each edge, i.e., the probability that 
         * an individual at the parent has no surviving descendants 
         * within the subtree rooted at the child node.
         * Indexing: dup0[u] is G(0) for the edge leading to u.
         */
        private double[] dup0;
        
        /**
         * Extinction probabilities at the nodes.
         * Indexing: extinction_probability[u] is the oprobability that an 
         * individual at u has no descendants at any of the leaves in the subtree rooted at u.
         */
        private double[] extinction_probability;
        
        /**
         * Product of G(0) over left children. 
         * Indexing: partial_extinction_probability[u][i] is the 
         *    produt of G(0) over children 0..i of u (children indexed in processing order) 
         */
        private double[][] partial_extinction_probability;
        
        private double[][] partial_extinction_powers;
        private double[][] partial_survival_powers;
        
        
        /**
         * Transition probabilities on each edge
         * Indexing: trans[u] is the transitions on the edge leading to node u
         */
        private EdgeTransitions[] trans;        
        
        /**
         * intializes the arrays used here
         */
        private void initDataStructures()
        {
            NodeWithRates[] tree_nodes = rate_tree.getDFT();
            trans = new EdgeTransitions[tree_nodes.length-1]; // no edge from the root up
            for (int node_idx=0; node_idx<trans.length; node_idx++)
                trans[node_idx]=new EdgeTransitions(node_idx);
            extinction_probability = new double[tree_nodes.length];
            dup0 = new double[tree_nodes.length-1];
            partial_extinction_probability = new double[tree_nodes.length][];
            for (int node_idx=0; node_idx<tree_nodes.length; node_idx++)
            {
                NodeWithRates N = rate_tree.getNode(node_idx);
                if (!N.isLeaf())
                    partial_extinction_probability[node_idx]=new double[N.getNumChildren()];
            }
            partial_extinction_powers = new double[tree_nodes.length][];
            partial_survival_powers = new double[tree_nodes.length][];
            for (int parent_idx=0; parent_idx<tree_nodes.length; parent_idx++)
            {
                NodeWithRates P = rate_tree.getNode(parent_idx);
                if (!P.isLeaf())
                {
                    for (int ci=0; ci<P.getNumChildren(); ci++)
                    {
                        int child_idx = rate_tree.getChildIndex(parent_idx,getChildForStep(parent_idx,ci));
                    }
                }
            }
        }
        
        /**
         * recomputes the transition probabilities for a branch.

         * @param node_idx child node of the branch
         */
        public void recomputeAtNode(int node_idx)
        {
//            long T0 = System.currentTimeMillis();
            NodeWithRates N = rate_tree.getNode(node_idx);
            // compute the extinction probability
            if (N.isLeaf())
                extinction_probability[node_idx]=0.0;
            else
            {
                double p=1.0;
                for (int ci=0; ci<N.getNumChildren(); ci++)
                {
                    int child_idx = rate_tree.getChildIndex(node_idx, getChildForStep(node_idx,ci));
                    p*=dup0[child_idx];
                    //partial_extinction_probability[node_idx][ci]=p;
                }
                extinction_probability[node_idx]=p;
//                System.out.println("#*SC.P.rAN "+node_idx+"\textinct "+extinction_probability[node_idx]);
            }
            
            if (!N.isRoot())
            {
                // compute duplication distribution at this guy
                DiscreteDistribution G = N.getDuplicationDistribution(extinction_probability[node_idx]);
                dup0[node_idx] = G.getDistribution(0)[0];
                
//                System.out.println("#*SC.P.rAN "+node_idx+"\tdup0 "+dup0[node_idx]);

                EdgeTransitions E = trans[node_idx];
                E.recompute(N, extinction_probability[node_idx]);
                recomputePartialsAtNode(node_idx);
            } // not root
//            timeP += System.currentTimeMillis()-T0;
        }



        /**
         * Updates the partial_extinction and partial_survival values at all the right siblings too
         * 
         * @param node_idx
         */
        private void recomputePartialsAtNode(int node_idx)
        {
            int parent_idx = rate_tree.getParentIndex(node_idx);
            int num_children = rate_tree.getNode(parent_idx).getNumChildren();
            int child_order = getStepForChild(parent_idx,rate_tree.getNode(node_idx).getIndexAtParent());
            for (;child_order<num_children;child_order++)
            {
                int child_position = getChildForStep(parent_idx, child_order);
                node_idx = rate_tree.getChildIndex(parent_idx, child_position);
                if (child_order==0)
                    partial_extinction_probability[parent_idx][child_order] = dup0[node_idx];
                else 
                    partial_extinction_probability[parent_idx][child_order]
                                = dup0[node_idx]*partial_extinction_probability[parent_idx][child_order-1];

                // compute partial_extinction_powers and partial_survival_powers
                if (child_order!=0)
                {
                    double D2 = partial_extinction_probability[parent_idx][child_order];
                    int M = max_combined_sizes[node_idx];
                    partial_survival_powers[node_idx] = Functions.powers(1.-D2, M);
//                    Verbose.message("SC.P.rAN psur "+node_idx+"/"+main_tree.getNode(node_idx).newickName()+" ("+parent_idx+":"+child_order+")\tD2 "+D2+"\tM "+M+"\t[] "+java.util.Arrays.toString(partial_survival_powers[node_idx]));
                }
                if (child_order != rate_tree.getNode(parent_idx).getNumChildren()-1)
                {
                    double D1 = partial_extinction_probability[parent_idx][child_order];
                    int M = max_combined_sizes[rate_tree.getChildIndex(parent_idx, getChildForStep(parent_idx,child_order+1))];
                    partial_extinction_powers[node_idx] = Functions.powers(D1, M);
//                    Verbose.message("SC.P.rAN pext "+node_idx+"/"+main_tree.getNode(node_idx).newickName()+" ("+parent_idx+":"+child_order+")\tD1 "+D1+"\tM "+M+"\t[] "+java.util.Arrays.toString(partial_extinction_powers[node_idx]));
                }
            }
        }
        
        /**
         * This is the array of (1-D)^k where D is the partial exinction probability 
         * in the lineages below a node and the node's left siblings (<q>left</q> defined
         * in terms of optimal child ordering).
         * 
         * @param node_idx Index of the child node on the edge 
         * @return array of [(1-D)^k]
         */
        public double[] getPartialSurvivalPowers(int node_idx)
        {
            return partial_survival_powers[node_idx];
        }
        
        /**
         * This is the array of D^k where D is the partial exinction probability 
         * in the lineages below a node and the node's left siblings (<q>left</q> defined
         * in terms of optimal child ordering).
         * 
         * @param node_idx Index of the child node on the edge 
         * @return array of [D^k]
         */
        public double[] getPartialExtinctionPowers(int node_idx)
        {
            return partial_extinction_powers[node_idx];
        }
        
        /**
         * Recomputes the transition probabilities on a path to the root (useful if branch parameters are updated)
         *  
         * @param node_idx child node of the lowest branch
         */
        public void recomputeOnPathToRoot(int node_idx)
        {
            recomputeAtNode(node_idx);
            if (!main_tree.getNode(node_idx).isRoot())
                recomputeOnPathToRoot(main_tree.getParentIndex(node_idx));
        }

        /**
         * Recomputes the transition probabilities on all edges.
         */
        public void recomputeAll()
        {
            for (int node_idx=0; node_idx<rate_tree.getNumNodes(); node_idx++)
            {
                NodeWithRates N = main_tree.getNode(node_idx);
                if (!N.isLeaf())
                {
                    for (int ci = 0; ci<N.getNumChildren(); ci++)
                    {
                        int child_idx = main_tree.getChildIndex(node_idx,getChildForStep(node_idx,ci));
                        recomputeAtNode(child_idx);
                    }
                }
            }
            recomputeAtNode(main_tree.getNumNodes()-1); // root
        }
        
//        /**
//         * Access to transition probabilities.
//         * <var>p</var>(<var>m</var>,<var>s</var>,<var>t</var>) is the probability that
//         * <var>Y</var>+<b>sum</b><sub><var>i</var>=1</sub><sup><var>s</var>+<var>t</var></sup><var>X</var><sub>i</sub>
//         * =<var>m</var> and <var>X</var><sub>i</sub>&gt;0
//         * for all <var>i</var>=1,...,<var>s</var>, where
//         * <var>Y</var> has the distribution of surviving xenologs <var>H</var>(<var>.</var>)
//         * and <var>X</var><sub>i</sub> are independent having identical distributions
//         * with surviving inparalogs <var>G</var>(.) for an individual.
//         *
//         * @param node_idx index for the child node
//         * @param m the population size at the child node; m &gt; s-1 is assumed
//         * @param s number of surviving lineages towards the child node
//         * @param t number of additional lineages towards the child node
//         * @return p(m,s,t)
//         */
//        public double getSurvival(int node_idx, int m, int s, int t)
//        {
//            //if (t!=0)
//            //{
//            //    throw new IllegalArgumentException("StableComputation.Probablities.getSurvival must have t=0: node "+node_idx+"\ts "+s+"\tt "+t+"\tm "+m);
//            //}
//            return trans[node_idx].get(m, s, t);
//        }

        public void computeTransitionsAtNode(int node_idx, int i_min, int length, double[] previous_W, double[] current_W)
        {
            trans[node_idx].recursion(rate_tree.getNode(node_idx), extinction_probability[node_idx], i_min, previous_W, current_W, length);
        }

        public double[] getSurvival(int node_idx, int s)
        {
            if (s<=trans[node_idx].getMaxAvailable_s())
                return trans[node_idx].get(s);
            else
                return null;
        }

        /**
         * Returns the probability for 0 surviving inparalogs on an edge, i.e., G(0).
         * 
         * @param edge_idx is the index of the child node on this edge
         * @return probability of complete extinction for a single individual's offsprings along this edge, or below its subtree
         */
        public double getLineageExtinction(int edge_idx)
        {
            return dup0[edge_idx];
        }
        
        
        /**
         * Returns the extinction probability for a subtree
         * 
         * @param node_idx index of the subtree root
         * @return probability that an individual at the node has no descendants at any of the leaves below
         */
        public double getExtinction(int node_idx)
        {
            return extinction_probability[node_idx];
        }
        
        public double getPartialExtinction(int node_idx, int child_step)
        {
            return partial_extinction_probability[node_idx][child_step];
        }
        
          /**
         * Computes the distribution for non-extinct individuals 
         * from the root.
         *
         * @param root_prior original prior distribution at root: must be one of NegativeBinomial, Poisson, PointDistribution, or ShiftedGeometric
         * @return null if unrecognized distribution
         */
        public DiscreteDistribution survivingRootPrior(DiscreteDistribution root_prior)
        {
            int root_idx = rate_tree.getNumNodes()-1;
            double D = extinction_probability[root_idx];
            return survivingRootPrior(root_prior,D);
        }
            
        /**
         * Computes the distribution for non-extinct individuals 
         * from the root.
         *
         * @param root_prior original prior distribution at root: must be one of NegativeBinomial, Poisson, PointDistribution, or ShiftedGeometric
         * @param D extinction probability
         * @return null if unrecognized distribution
         */
        public DiscreteDistribution survivingRootPrior(DiscreteDistribution root_prior, double D)
        {
            double[] root_params = root_prior.getParameters();
            DiscreteDistribution corrected_prior = null;

            if (root_prior instanceof Poisson)
            {
                double r = root_params[0];
                double class_r = r*(1.-D);
                corrected_prior = new Poisson(class_r);
            } else if (root_prior instanceof ShiftedGeometric)
            {
                double p = root_params[0];
                double q = root_params[1];
                double u = 1.-q*D;
                double class_p = (p*(1.-D)+(1.-q)*D)/u;
                double class_q = q*(1.-D)/u;
                corrected_prior=new ShiftedGeometric(class_p,class_q);
            } else if (root_prior instanceof NegativeBinomial)
            {
                double t = root_params[0];
                double q = root_params[1];
                double u = 1.-q*D;
                double class_q = q*(1.-D)/u;
                corrected_prior=new NegativeBinomial(t,class_q);
            } else if (root_prior instanceof PointDistribution)
            {
                double p = root_params[0];
                double class_p = p+(1.-p)*D;
                corrected_prior=new PointDistribution(class_p);
            }
            return corrected_prior;
        }
      
        /**
         * Computes the probability of all-0: it's simply the product of H_y(0) over all edges xy
         * 
         * @return probability that there is no surviving descendant at any leaf
         */
        public double getAllAbsentProbability()
        {
            double p_absent = 1.0;
            int num_nodes = rate_tree.getNumNodes();
            for (int node_idx=0; node_idx<num_nodes; node_idx++)
            {
                NodeWithRates N = rate_tree.getNode(node_idx);
                if (!N.isRoot()) // i.e., node_idx<num_nodes-1
                {
                    EdgeTransitions E = trans[node_idx];
                    p_absent *= E.get(0,0);
                }
            }
           return p_absent;
        }

        
        /**
         * Probability that a family is represented in at most one lineage
         * 
         * @param root_prior distribution at the root
         * @return probability
         */
        public double getAbsentOrSingle(DiscreteDistribution root_prior)
        {
            int num_leaves = rate_tree.getNumLeaves();
            double[] present_at = getSingleLineagePresence(root_prior);
            double psum = 0.0;
            for (int leaf_idx=0; leaf_idx<num_leaves; leaf_idx++)
            {
                psum += present_at[leaf_idx];
            }
            double p_absent = getAllAbsentProbability()*survivingRootPrior(root_prior).getDistribution(0)[0];
            double pcorr = psum+ p_absent;
            //Verbose.message("SC.P.gAOS sum "+psum+"\tabs "+p_absent+"\tdiff "+pcorr+"\tp1 "+p1);
            
            return pcorr;
        }
        
        private double[] getSingleLineagePresence(DiscreteDistribution root_prior)
        {

            int num_nodes = rate_tree.getNumNodes();
            int num_edges = num_nodes-1;
            int root_idx = num_edges;
            double [] outside_extinct = new double[num_nodes];
            double[] H0 = new double[num_edges];
            double[] G0 = new double[num_edges];
            System.arraycopy(extinction_probability,0,outside_extinct,0,num_nodes);
            System.arraycopy(dup0,0,G0,0,num_edges);
            for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
            {
                NodeWithRates N = rate_tree.getNode(edge_idx);
//                System.out.println("#*SC.P.gSLP node "+edge_idx+"/"+N.newickName()+"\trp "+root_prior+"\tN "+N);
                DiscreteDistribution Hy = N.getTransferDistribution(extinction_probability[edge_idx]);
//                System.out.println("#*SC.P.gSLP node "+edge_idx+"\tp0 "+extinction_probability[edge_idx]+"\tH "+Hy);
//                System.out.flush();

                H0[edge_idx] = Hy.getDistribution(0)[0];
            }
            int num_leaves = rate_tree.getNumLeaves();
            double[] present_at = new double[num_leaves];
            for (int leaf_idx=0; leaf_idx<num_leaves; leaf_idx++)
            {
                NodeWithRates Leaf = rate_tree.getNode(leaf_idx);
                int leaf_height=0;
                for (int n=leaf_idx; n!=num_nodes-1; n=rate_tree.getParentIndex(n))
                    leaf_height++;
                int[] path = new int[leaf_height+1];
                leaf_height=0;
                for (int n=leaf_idx; n!=num_nodes-1; n=rate_tree.getParentIndex(n))
                {
                    path[leaf_height]=n;
                    //Verbose.message("SPF.R.gAS "+leaf_idx+" path "+leaf_height+"\t"+rate_tree.getNode(n));
                    leaf_height++;
                }
                path[path.length-1]=num_nodes-1; // root

                // recompute extinction probability along the path to the root
                for (int path_idx=1; path_idx<path.length; path_idx++)
                {
                    int ancestor_idx = path[path_idx];
                    NodeWithRates A = rate_tree.getNode(ancestor_idx);
                    outside_extinct[ancestor_idx]=1.0;
                    for (int ci=0; ci<A.getNumChildren(); ci++)
                    {
                        int child_idx = rate_tree.getChildIndex(ancestor_idx,ci);
                        if (child_idx != leaf_idx) // only fails with the immediate parent
                            outside_extinct[ancestor_idx]*=G0[child_idx];
                    }
                    if (!A.isRoot())
                    {
                        double D = outside_extinct[ancestor_idx];
                        double g0 = A.getDuplicationDistribution(D).getDistribution(0)[0];
                        double h0 = A.getTransferDistribution(D).getDistribution(0)[0];
                        G0[ancestor_idx] = g0;
                        H0[ancestor_idx] = h0;
                    }
                }
                double product_h0 = 1.0;
                for (int edge_idx=0; edge_idx<num_edges; edge_idx++)
                    if (edge_idx != leaf_idx)
                        product_h0 *= H0[edge_idx];
                double p0 = survivingRootPrior(root_prior,outside_extinct[root_idx]).getDistribution(0)[0];
                present_at[leaf_idx] = product_h0*p0;
                //Verbose.message("SPF.R.gAOS "+leaf_idx+"/"+rate_tree.getNode(leaf_idx).getName()+"\th "+product_h0+"\tp0 "+p0);

                // reset values
                for (int path_idx=1; path_idx<path.length; path_idx++)
                {
                    int ancestor_idx = path[path_idx];
                    NodeWithRates A = rate_tree.getNode(ancestor_idx);
                    outside_extinct[ancestor_idx]=1.0;
                    for (int ci=0; ci<A.getNumChildren(); ci++)
                    {
                        int child_idx = rate_tree.getChildIndex(ancestor_idx,ci);
                        outside_extinct[ancestor_idx]*=G0[child_idx];
                    }
                    if (!A.isRoot())
                    {
                        double D = outside_extinct[ancestor_idx];
                        double g0 = A.getDuplicationDistribution(D).getDistribution(0)[0];
                        double h0 = A.getTransferDistribution(D).getDistribution(0)[0];
                        G0[ancestor_idx] = g0;
                        H0[ancestor_idx] = h0;
                    }
                }
            }
            double p_absent = getAllAbsentProbability()*survivingRootPrior(root_prior).getDistribution(0)[0];
            for (int leaf_idx=0; leaf_idx<num_leaves; leaf_idx++)
            {
                present_at[leaf_idx]-=p_absent;
            }
            return present_at;
        }        
        
    }
    



    
    /**
     * Controlled access to p(m,s) transition probabilities.
     * The probabilities are not computed here, this
     * class is only for storage.
     */
    private class EdgeTransitions
    {

        /**
         * Initializes the data structure for a specific edge
         * 
         * @param node_idx index of the child node the the edge leads to
         */
        private EdgeTransitions(int node_idx)
        {
            this.node_idx = node_idx;
            initArray();
        }
        
        private int node_idx;
        private double[][] sp;
        
        /**
         * Allocates the appropriate space for the sp array.
         */
        private void initArray()
        {
            // we only need p(m,s,0) where s = 0...M, m=s...M, and M is the sum of sizes below this node
            int max_size = max_combined_sizes[node_idx];
            sp = new double[Math.min(max_size, MAX_PRECOMPUTED_TRANSITIONS)+1][];
            for (int s=0; s<sp.length; s++)
            {
                int full_size = max_size-s+1;
                int precomputed_size = full_size;
                if (s>0 && full_size > MAX_PRECOMPUTED_TRANSITIONS)
                    precomputed_size = MAX_PRECOMPUTED_TRANSITIONS;
                sp[s]=new double[precomputed_size];
            }
        }

        /**
         * Access to transition probabilities:
         * <var>w</var>*(<var>m</var>|<var>s</var>) is the probability
         * that
         * <var>Y</var>+&sum;<sub><var>i</var>=1,...,<var>s</var></sub><var>X</var><sub>i</sub>
         * =<var>m</var> and <var>X</var><sub>i</sub>&gt;0
         * for all <var>i</var>=1,...,<var>s</var>, where 
         * <var>Y</var> has the distribution of surviving xenologs <var>H</var>(<var>.</var>)
         * and <var>X</var><sub>i</sub> are independent having identical distributions 
         * with surviving inparalogs <var>G</var>(.) for an individual.
         * 
         * @param m the population size at the child node; m &gt; s-1 is assumed
         * @param s number of surviving lineages towards the child node
         * @return <var>w</var>*(<var>m</var>|<var>s</var>)
         */
        public final double get(int m, int s)
        {
            return sp[s][m-s];
        }
        
        /**
         * Access to transition probabilities:
         * returns an array <var>W</var>[]
         * where <var>W</var>[<var>i</var>] = <var>w<var>*(<var>i</var>+<var>s</var> | <var>s</var>).
         *
         * @see #get(int, int)
         *
         * @param s number of surviving lineages towards the child node
         * @return <var>W</var>[] of transition probabilities
         */
        public final double[] get(int s)
        {
            return sp[s];
        }
        
        /**
         * Sets <var>w*</var>(<var>m</var> | <var>s</var>).
         * 
         * @param m the population size at the child node; m &ge; s is assumed
         * @param s number of surviving lineages towards the child node
         * @param w the probability <var>w*</var>(<var>m</var> | <var>s</var>)
         */
        public final void set(int m, int s, double w)
        {
            sp[s][m-s]=w;
        }


        /**
         * Sets all <var>w*</var>( ... | <var>s</var>).
         *
         *
         * @param s number of surviving lineages towards the child node
         * @param W array <var>W</var>[<var>i</var>] = <var>w<var>*(<var>i</var>+<var>s</var> | <var>s</var>)
         */
        public final void set(int s, double[] W)
        {
            System.arraycopy(W, 0, sp[s], 0, sp[s].length);
        }
        
        /**
         * Maximum value of s for which p(*,s) is available.
         * 
         * @return maximum of s
         */
        public final int getMaxAvailable_s()
        {
            return sp.length-1;
        }
        
        /**
         * Maximum value of m for which <var>w</var>*(<var>m</var> | <var>s</var>) is precomputed.
         * 
         * @param s number of surviving lineages towards the child node
         * @return maximum of m
         */
        public final int getMaxAvailable_m(int s)
        {
            return s+sp[s].length-1;
        }

        private void recursion(double dup1, double q, int pos, double[] previous, double[] current, int len)
        {
            System.arraycopy(previous, pos, current, pos, len-pos);
            for (int i=pos; i<len; i++)
                current[i] *= dup1;
            if (pos==0) pos++;
            for (int i=pos; i<len; i++)
                current[i] += current[i-1]*q;
        }

        /**
         * Recursion for computing the
         * array <var>W</var>[<var>i</var>] = <var>w</var>*(<var>i</var>+<var>s</var> | <var>s</var>)
         *
         * @param N child node on edge
         * @param p_extinct probability of extinction at N
         * @param pos first <var>i</var> for which <var>w</var>*(<var>i</var>+<var>s</var> | <var>s</var>) is computed
         * @param previous array <var>W</var> for <var>s</var>-1
         * @param current here is where the results are placed
         * @param len computation is requested for <var>W</var>[pos..len-1]
         */
        private void recursion(NodeWithRates N, double p_extinct, int pos, double[] previous, double[] current, int len)
        {
            DiscreteDistribution G = N.getDuplicationDistribution(p_extinct);
            double[] params = G.getParameters();
            double p1 = 1.0-params[0];
            double dup1 = p1;
            double q = 0.0;
            if (params.length == 2) // shifted geometric
            {
                q = params[1];
                dup1 *= (1.-q);
            }
            recursion(dup1, q, pos, previous, current, len);
        }

        public void recompute(NodeWithRates N, double p_extinct)
        {
            double p,q,p1;
            {
                // compute duplication distribution at this guy
                DiscreteDistribution G = N.getDuplicationDistribution(p_extinct);
                // G might be a point distribution if duplication rate is 0.
                if (G instanceof ShiftedGeometric)
                {
                    double[] params = G.getParameters();
                    p = params[0];
                    q = params[1];
                } else
                {
                    // duplication rate is 0...
                    p = G.getParameters()[0];
                    q = 0.0;
                }
                p1 = 1.0-p;
            }

            {
                // compute p(m,0,0)
//                System.out.println("#*SC.ET.r pext "+p_extinct+"\tN "+N);
//                if (N.getLossRate()==0.0)
//                {
//                    NodeWithRates N0 = main_tree.getNode(N.getId());
//                    System.out.println("#*SC.ET.r N0 "+N0);
//                }
                DiscreteDistribution H = N.getTransferDistribution(p_extinct);
//                System.out.println("#*SC.ET.r pext "+p_extinct+"\tH "+H+"\t// "+N);
//                System.out.flush();

                int mmax = getMaxAvailable_m(0);
                double Hm[] = H.getDistribution(mmax);
                System.arraycopy(Hm, 0, sp[0], 0, Hm.length);
            }
            double dup1=p1*(1.-q);

            // compute p(m,s) for all m>=s and s>0
            for (int s=1; s<sp.length; s++)
            {
                recursion(dup1, q, 0, sp[s-1], sp[s], sp[s].length);
            } // for s
        }
    }    
    
    protected class LowerConditionals
    {
        protected LowerConditionals(PhyleticProfile profile, Probabilities survival)
        {
            this.survival = survival;
            this.profile = profile;
            rate_tree = survival.getTree();
            initDataStructures();
            recomputeAll();
        }
        
        protected Probabilities survival;
        protected PhyleticProfile profile;
        protected TreeWithRates rate_tree;

        /**
         * Conditional likelihoods for each edge: survival is considered 
         * in the subtree below this edge and its siblings to the left.
         */
        protected double[][] edge_likelihood;
        
        /**
         * Conditional likelihoods at nodes: for a non-leaf node, 
         * node_likelihood[u] is the same as the edge_likelihood on
         * the edge leading to u's rightmost child.
         */
        protected double[][] node_likelihood;
        
        //private int[] nonzero_family_size;
        
        public double[] getLikelihoods(int node_idx)
        {
//            System.out.println("#*SC.gL prof "+profile.getPatternString()+"\tnode "+node_idx+"\t"+java.util.Arrays.toString(node_likelihood[node_idx]));

            return node_likelihood[node_idx];
        }
        
        protected void initDataStructures()
        {
            NodeWithRates[] nodes = rate_tree.getDFT();
            edge_likelihood = new double[nodes.length-1][];
            node_likelihood = new double[nodes.length][];
            //nonzero_family_size = new int[nodes.length];
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                NodeWithRates N = rate_tree.getNode(node_idx);
                if (N.isLeaf())
                {
                    int n = profile.get(node_idx);
                    node_likelihood[node_idx] = new double[n+1];
                } else
                {
                    int num_children = N.getNumChildren();
                    int sum_m = 0;
                    for (int ci=0; ci<num_children; ci++)
                    {
                        int child_position = getChildForStep(node_idx,ci);
                        int child_idx = rate_tree.getChildIndex(node_idx, child_position);
                        int m = node_likelihood[child_idx].length-1;
                        sum_m += m;
                        
                        if (ci==0)
                            edge_likelihood[child_idx]=new double[sum_m+1];
                        else
                            edge_likelihood[child_idx]=new double[sum_m+1];
                        
                        
                    }
                    int last_child_idx = rate_tree.getChildIndex(node_idx, getChildForStep(node_idx,num_children-1));
                    node_likelihood[node_idx] = edge_likelihood[last_child_idx];
                }
            }
            //if (Verbose.isVerbose())
            //{
            //    int[] s = profile.computeSubtreeSizes(rate_tree);
            //    for (int node_idx=0; node_idx<nodes.length; node_idx++)
            //    {
            //        Verbose.message("SC.LC.iDS node "+node_idx+"/"+rate_tree.getNode(node_idx).newickName()
            //                +"\tsize "+s[node_idx]+"\tnl[] "+node_likelihood[node_idx].length);
            //        for (int ci=0; ci<rate_tree.getNode(node_idx).getNumChildren(); ci++)
            //        {
            //            int child_idx = rate_tree.getChildIndex(node_idx, getChildForStep(node_idx,ci));
            //            Verbose.message("SC.LC.iDS edge "+ci+"/"+rate_tree.getNode(child_idx).newickName()
            //                    +"\tsize "+s[child_idx]+"\tel[] "+edge_likelihood[child_idx].length);
            //        }
            //    }
                
            //}
            
        }
        
        private void recomputeAtNode(int node_idx)
        {
//            long T0 = System.currentTimeMillis();
            NodeWithRates N = rate_tree.getNode(node_idx);
            if (N.isLeaf())
            {
                java.util.Arrays.fill(node_likelihood[node_idx],0.0);
                int v = profile.get(node_idx);
                node_likelihood[node_idx][v] = 1.0;
                
                
                //nonzero_family_size[node_idx]=node_likelihood[node_idx].length-1;
                //if (Verbose.isVerbose())
                //{
                //    Verbose.message("SC.LC.rAN leaf "+node_idx+"/"+N.newickName()+" ================\t"+N);
                //    Verbose.message("SC.LC.rAN "+profile.getPatternString()+"\t"+node_idx+"\tL["+v+"]="+node_likelihood[node_idx][v]);
                //}
            } else
            {
                // figure out the processing order
                int num_children = N.getNumChildren();
                int[] child_indexes = new int[num_children];
                for (int ci=0; ci<num_children; ci++)
                {
                    int child_position = getChildForStep(node_idx,ci);
                    int child_idx = rate_tree.getChildIndex(node_idx, child_position);
                    child_indexes[ci]=child_idx;
                }
//                double extinct=0.0;
                int M1;
                int M2;
                
                if (num_children == 1)
                {
                    int node1 = child_indexes[0];
                    double extinct = survival.getLineageExtinction(node1);
                    M1 = node_likelihood[node1].length-1; //nonzero_family_size[node1];
                    double[] D1pow = Functions.powers(1.-extinct, M1);
                    double[] w_current = new double[M1+1];
                    double[] w_previous = new double[M1+1];
                    for (int s=0; s<=M1; s++)
                    {
                        int desired_length = M1-s+1;
                        { // compute w*( | )
                            double[] p_array = survival.getSurvival(node1, s);
                            if (p_array == null || p_array.length<desired_length)
                            {
                                int i_min = 0;
                                if (p_array != null)
                                {
                                    i_min = p_array.length;
                                    System.arraycopy(p_array,0,w_current,0,i_min);
                                }
                                survival.computeTransitionsAtNode(node1, i_min, desired_length, w_previous, w_current);

                            } else
                                System.arraycopy(p_array,0,w_current,0,desired_length);
                        }

                        double L1 = 0.0;
                        for (int m=s; m<=M1; m++)
                        {
                            double p = w_current[m-s]; // survival.getSurvival(node1, m, t, 0);
                            double Lm = node_likelihood[node1][m];
                            L1 += p*Lm;
                            //Verbose.message("SC.LC.rAN left "+node1+"\tt "+t+"\tm "+m+"\tp "+p+"\tLm "+Lm+"\tL1 "+L1);
                        }
                        edge_likelihood[node1][s]=L1/D1pow[s];

                        if (s!=M1)
                            System.arraycopy(w_current,0,w_previous,0,desired_length);
                    }
                } else 
                {   // compute edge likelihood for second child
                    int node1 = child_indexes[0];
                    int node2 = child_indexes[1];
                    
//                    extinct = survival.getLineageExtinction(node1)*survival.getLineageExtinction(node2);
//                    System.out.println("#*SC.LC.rAN node "+node_idx+"\tc1 "+node1+"\tc2 "+node2+"\te "+extinct+"\te1 "+survival.getExtinction(node1)+"\te2 "+survival.getExtinction(node2));
                    
                    M1 = node_likelihood[node1].length-1; //nonzero_family_size[node1];
                    M2 = node_likelihood[node2].length-1;//nonzero_family_size[node2];
                    for (int n=0; n<edge_likelihood[node2].length; n++)
                        edge_likelihood[node2][n]=0.0;
                    // powers of (1-extinct)
                    double[] D2pow = survival.getPartialSurvivalPowers(node2); // Functions.powers(1.-extinct, M1+M2);
                    // powers of G0
                    double[] Gpow = survival.getPartialExtinctionPowers(node1); // Functions.powers(survival.getLineageExtinction(node1), M2);
                    
                    double[] L2ts_previous=new double[M2+1]; // for t-1
                    double[] L2ts_current = new double[M2+1]; // for t

                    double[] w1_current = new double[M1+1];
                    double[] w1_previous = new double[M1+1];
                    for (int t=0; t<=M1; t++)
                    {
                        double L1 = 0.0;
                        {
                            int desired_length = M1-t+1;
                            {
                                double[] p_array = survival.getSurvival(node1, t);
//                                System.out.println("#*SC.LC.rAN node1 "+node1+"/"+main_tree.getNode(node1).newickName()
//                                        +"\tdesired "+desired_length
//                                        +"\tsurvival["+0+"] "
//                                        +(p_array==null?"null":java.util.Arrays.toString(p_array)));

                                if (p_array == null || p_array.length<desired_length)
                                {
                                    int i_min = 0;
                                    if (p_array != null)
                                    {
                                        i_min = p_array.length;
                                        System.arraycopy(p_array,0,w1_current,0,i_min);
                                        {
                                            boolean b = false;
                                            for (int i=0; i<p_array.length && !b; ++i)
                                                b = Double.isNaN(p_array[i]);
                                            if (b)
                                            {
                                                System.out.println("#SC.LC.rAN NUMERROR1 numc 2\tsurvival node1 "+node1+"/"+main_tree.getNode(node1).newickName()+"\tt "+t+"\tp[] "+java.util.Arrays.toString(p_array));
                                                System.out.flush();
                                            }
                                            assert (!b);
                                        }
                                    }
                                    if (t!=0)
                                        survival.computeTransitionsAtNode(node1, i_min, desired_length, w1_previous, w1_current);
                                    else
                                        throw new RuntimeException("Cannot compute recursion for s=0");

                                } else
                                    System.arraycopy(p_array,0,w1_current,0,desired_length);

//                                if (t!=0 && desired_length>2)
//                                    survival.computeTransitionsAtNode(node1,2, desired_length, w1_previous, w1_current);

                            }

                            for (int m=t; m<=M1; m++)
                            {
                                double p = w1_current[m-t]; // survival.getSurvival(node1, m, t, 0);
                                double Lm = node_likelihood[node1][m];
                                if (Double.isNaN(p))
                                {
                                    System.out.println("#SC.LC.rAN NUMERROR2 numc 2\tnode1 "+node1+"/"+main_tree.getNode(node1).newickName()
                                            +"\tt "+t+"\tm "+m+"\tw1[] "
                                            +java.util.Arrays.toString(w1_current)
                                            +"\tprev "+java.util.Arrays.toString(w1_previous));
                                    System.out.flush();

                                }
                                assert (!Double.isNaN(p));
                                assert (!Double.isNaN(Lm));
                                
                                
                                L1 += p*Lm;
//                                if (DEBUG22)
//                                {
//                                	System.out.println("#*SC.LC.rAN L1 "+node1+"\tell "+m+"\tt "+t+"\ts "+(m-t)+"\tw* "+p+"\tL(ell) "+Lm);
//                                }
                            }
                            if (t!=M1)
                                System.arraycopy(w1_current, 0, w1_previous, 0, desired_length);
                        }
                        
//                        if (DEBUG22)
//                        {
//                        	L1 *= Math.exp(t*Math.log(1.0-Gpow[1]));
//                        }
                        
//                        System.out.println("SC.LC.rAN left "+node1+"\tt "+"\tL1 "+L1+"\tLL1 "+Math.log(L1));
                        
                    	if (t==0)
                        {
                            double[] w2_current = new double[M2+1];
                            double[] w2_previous = new double[M2+1];
                            for (int s=0; s<=M2; s++)
                            {
                                int desired_length = M2-s+1;
                                {
                                    double[] p_array2 = survival.getSurvival(node2, s);
                                    if (p_array2 == null || p_array2.length<desired_length)
                                    {
                                        int i_min = 0;
                                        if (p_array2 != null)
                                        {
                                            i_min = p_array2.length;
                                            System.arraycopy(p_array2,0,w2_current,0,i_min);
                                        }
                                        if (s!=0)
                                            survival.computeTransitionsAtNode(node2, i_min, desired_length, w2_previous, w2_current);
                                        else
                                            throw new RuntimeException("Cannot compute recursion for s=0");
                                    } else
                                        System.arraycopy(p_array2,0,w2_current,0,desired_length);

//                                    if (s!=0 && desired_length>2)
//                                        survival.computeTransitionsAtNode(node2, 2, desired_length, w2_previous, w2_current);
    
                                }
                                double L2 = 0.0;
                                for (int m=s; m<=M2; m++)
                                {
                                    double p=w2_current[m-s];
                                    double Lm = node_likelihood[node2][m];
                                    //Verbose.message("SC.LC.rAN right "+node2+"\tt "+t+"\ts "+s+"\tm "+m+"\tp "+p+"\tLm "+Lm+"\tL2 "+L2);                                
                                    assert (!Double.isNaN(p));
                                    assert (!Double.isNaN(Lm));
                                    L2 += p*Lm;
                                }
                                L2ts_current[s] = L2;
                                if (s!=M2)
                                    System.arraycopy(w2_current, 0, w2_previous, 0, desired_length);
                            }
                        } else
                        {
                            double G0 =  survival.getLineageExtinction(node2);
                            assert (!Double.isNaN(G0));
                            for (int s=0; s<M2; s++)
                            {
                                //double L2 = 0.0;
                                //double[] p_array2 = survival.getSurvival(node2, s, t);
                                //for (int m=s; m<=M2; m++)
                                //{
                                //    double p=p_array2[m-s];
                                //    double Lm = node_likelihood[node2][m];
                                //    //Verbose.message("SC.LC.rAN right "+node2+"\tt "+t+"\ts "+s+"\tm "+m+"\tp "+p+"\tLm "+Lm+"\tL2 "+L2);                                
                                //    L2 += p*Lm;
                                //}
                                L2ts_current[s] = L2ts_previous[s+1] + G0*L2ts_previous[s];
                            }
                            L2ts_current[M2] = G0*L2ts_previous[M2];
                        } // t>0

                        double previous_multiplier = 0.0;
                        int smax = Math.min(M2,edge_likelihood[node2].length-1-t);
                        for (int s=0; s<=smax; s++)
                        {
                            double L2 = L2ts_current[s];
                            double ts_choose_s = getNchooseM(t+s,s);
                            double f = ts_choose_s*(Gpow[s]/D2pow[t+s]);
                            
                            
                            if (s != 0 && (Double.isNaN(f) || Double.isInfinite(f)))
                            {
                                // (t+s)!/(t!*s!) * a^s/b^(t+s) = ((t+s-1)!/(t!*(s-1)!) * a^(s-1)/b^(t+s-1)) * ((t+s)/s)*(a/b)
                                f = previous_multiplier / (s+0.) * Gpow[1]/D2pow[1] * (t+s+0.);
                            }
                            double term = f*L1*L2;

                            if (Double.isNaN(term) || Double.isInfinite(term) || term<0.0)
                            {
//                                System.out.println("#*SC.LC.rAN NUMERROR3 node "+node_idx+"/"+N.getTaxonName()+"\tt "+t+"\ts "+s
//                                        +"\tterm "+term+"\ttsc "+ts_choose_s+"\tD2^ "+D2pow[t+s]+"\tG "+Gpow[s]
//                                        +"\tL1 "+L1+"\tL2 "+L2
//                                        +"\tf "+f+"\tprev "+previous_multiplier
//                                        +"\tL2ts "+java.util.Arrays.toString(L2ts_current));
//                                System.out.print("#*SC.LC.rAN NUMERROR3 L*G {");
//                                for (int u=0; u<smax; u++)
//                                {
//                                    if (u!=0)
//                                        System.out.print(", ");
//                                    System.out.print(L2ts_current[u]*Gpow[u]);
//                                }
//                                System.out.println();
//                                System.out.flush();

                                if (L1==0.0 || L2==0.0)
                                    term = 0.0;
                                else
                                { // try with logarithms
                                    double logf = Functions.bicoln(t+s,s)+s*Math.log(Gpow[1])-(t+s)*Math.log(D2pow[1]);
                                    if (DEBUG22) f=Math.exp(logf);
                                    	
                                    double log1 = Math.log(L1);
                                    double log2 = Math.log(L2);
                                    term = Math.exp(logf+log1+log2);

                                    if (Double.isNaN(term) || Double.isInfinite(term) || term<0.0)
                                    {
                                        System.out.println("#*SC.LC.rAN NUMERROR3 node "+node_idx+"/"+N.getTaxonName()+"\tt "+t+"\ts "+s
                                                +"\tterm "+term+"\ttsc "+ts_choose_s+"\tD2^ "+D2pow[t+s]+"\tG "+Gpow[s]
                                                +"\tL1 "+L1+"\tL2 "+L2
                                                +"\tf "+f+"\tprev "+previous_multiplier
                                                +"\tL2ts "+java.util.Arrays.toString(L2ts_current));
                                        System.out.print("#*SC.LC.rAN NUMERROR3 L*G {");
                                        for (int u=0; u<smax; u++)
                                        {
                                            if (u!=0)
                                                System.out.print(", ");
                                            System.out.print(L2ts_current[u]*Gpow[u]);
                                        }
                                        System.out.println();
                                        System.out.flush();
                                    }


                                }
                            }

                            assert (!Double.isNaN(term));
                            
                            double e = edge_likelihood[node2][s+t]+term;
                            assert (!Double.isNaN(e));
                            
//                            if (DEBUG22)
//                            {
//                            	System.out.println("#*SC.LC.rAN "+N.getId()+"+"+node2+"\tell "+(s+t)+"\ts "+s+"\tt "+t+"\t"+term
//                            		+"\tL1/Ct "+L1+"\tL2/Ds "+L2+"\tf "+f);
//                            }
                            
                            edge_likelihood[node2][s+t] = e;
                            L2ts_previous[s] = L2ts_current[s];
                            previous_multiplier = f;
                        }
                    }
                    M1+=M2;
                    //while(M1>0 && edge_likelihood[node2][M1]==0.0)
                    //    M1--;
                }
                
                // compute edge likelihoods for other children
                for (int ci=2; ci<num_children; ci++)
                {
                    int node1 = child_indexes[ci-1];
                    int node2 = child_indexes[ci];
                    
                    M2 = node_likelihood[node2].length-1;//nonzero_family_size[node2];
                    //Verbose.message("SC.rAN "+node_idx+" child "+ci+" M1 "+M1+", M2 "+M2);
                    for (int n=0; n<edge_likelihood[node2].length; n++)
                        edge_likelihood[node2][n]=0.0;

//                    double previous_extinct = extinct;
//                    extinct *= survival.getLineageExtinction(node2);
                    
                    // compute powers of (1.-extinct)
                    double[] D2pow = survival.getPartialSurvivalPowers(node2);// Functions.powers(1.-extinct, M1+M2);
                    double[] D1pow = survival.getPartialExtinctionPowers(node1); // Functions.powers(previous_extinct, M2);
                    double[] D1_1pow = survival.getPartialSurvivalPowers(node1); // Functions.powers(1.-previous_extinct, M1);

                    double[] L2ts_previous=new double[M2+1]; // for t-1
                    double[] L2ts_current = new double[M2+1]; // for t

                    for (int t=0; t<=M1; t++)
                    {
                        double L1 = D1_1pow[t]*edge_likelihood[node1][t];
                        
                        if (t==0)
                        {
                            double[] w2_current = new double[M2+1];
                            double[] w2_previous = new double[M2+1];
                            for (int s=0; s<=M2; s++)
                            {
                                int desired_length = M2-s+1;
                                double L2 = 0.0;
                                {
                                    double[] p_array2 = survival.getSurvival(node2, s);
                                    if (p_array2 == null || p_array2.length<desired_length)
                                    {
                                        int i_min = 0;
                                        if (p_array2 != null)
                                        {
                                            i_min = p_array2.length;
                                            System.arraycopy(p_array2,0,w2_current,0,i_min);
                                        }
                                        if (s!=0)
                                            survival.computeTransitionsAtNode(node2, i_min, desired_length, w2_previous, w2_current);
                                        else
                                            throw new RuntimeException("Cannot compute recursion for s=0");
                                    } else
                                        System.arraycopy(p_array2,0,w2_current,0,desired_length);
//                                    if (s!=0 && desired_length>2)
//                                        survival.computeTransitionsAtNode(node2, 2, desired_length, w2_previous, w2_current);
                                }
                                for (int m=s; m<=M2; m++)
                                {
                                    double p=w2_current[m-s];
                                    double Lm = node_likelihood[node2][m];
                                    //Verbose.message("SC.LC.rAN right "+node2+"\tt "+t+"\ts "+s+"\tm "+m+"\tp "+p+"\tLm "+Lm+"\tL2 "+L2);                                
                                    L2 += p*Lm;
                                }
                                L2ts_current[s] = L2;
                                if (s!=M2)
                                    System.arraycopy(w2_current, 0, w2_previous, 0, desired_length);
                            }
                        } else
                        {
                            double G0 = survival.getLineageExtinction(node2);
                            for (int s=0; s<M2; s++)
                            {
                                
                                L2ts_current[s] = L2ts_previous[s+1] + G0*L2ts_previous[s];
                            }
                            L2ts_current[M2]=G0*L2ts_previous[M2];
                        } // t>0

                        double previous_multiplier = 0.0;
                        int smax = Math.min(M2,edge_likelihood[node2].length-1-t);
                        for (int s=0; s<=smax; s++)
                        {
                            double L2 = L2ts_current[s]; 
                            
                            //double L2 = 0.0;
                            //double[] p_array = survival.getSurvival(node2, s, t);
                            //for (int m=s; m<=M2; m++)
                            //{
                            //    double p = p_array[m-s];
                            //    double Lm = node_likelihood[node2][m];
                            //    L2 += p*Lm;
                            //}
                            double ts_choose_s = getNchooseM(t+s,s); 
                            double f = ts_choose_s/D2pow[t+s]*D1pow[s];
                            if (s != 0 && (Double.isNaN(f) || Double.isInfinite(f)))
                            {
//                                double fff = f;
                                f = previous_multiplier / (s+0.) * D1pow[1]/D2pow[1] * (t+s+0.);
//                                System.out.println("#*SC.LC.rAN RECOMP "+node_idx+"\tf "+fff+"\t-> "+f+"\tprev "+previous_multiplier+"\tt "+t+"\ts "+s+"\tD1 "+D1pow[1]+"\tD2 "+D2pow[1]);
                            }
                            double term = f*L1*L2;
                                                        
                            edge_likelihood[node2][s+t]+=term;
                            L2ts_previous[s] = L2ts_current[s];
                            previous_multiplier = f;
                            //Verbose.message("SC.LC.rAN edge["+ci+"] "+node2+"/"+rate_tree.getNode(node2).newickName()+"\tt "+t+"\ts "+s+"\tterm "+term+"\t// ts "+ts_choose_s+"\tL1 "+L1+"\tL2 "+L2+"\tex "+extinct+"\tD2 "+(1.-D2pow[1]));
                        } // for s
                    } // for t
                    M1+=M2;
                    //while(M1>0 && edge_likelihood[node2][M1]==0.0)
                    //    M1--;
                    
                } // for children
                
                // node_likelihood[node_idx] is now properly set, since it 
                // is the same as edge_likelihood[last_child]
            } // for non-leaf node
            
//            timeLC += System.currentTimeMillis()-T0;
//            if (DEBUG22)
//            {
//            	System.out.println("#*SC.LC.rAN node "+N.getId()+"\t"+Arrays.toString(node_likelihood[node_idx])+"\t// "+N);
//            }
            
            
        }
        
        public void recomputeOnPathToRoot(int node_idx)
        {
            recomputeAtNode(node_idx);
            if (!rate_tree.getNode(node_idx).isRoot())
                recomputeOnPathToRoot(rate_tree.getParentIndex(node_idx));
        }

        public void recomputeAll()
        {
            for (int node_idx=0; node_idx<rate_tree.getNumNodes(); node_idx++)
            {
                NodeWithRates N = rate_tree.getNode(node_idx);
                if (!N.isLeaf())
                {
                    for (int ci = 0; ci<N.getNumChildren(); ci++)
                    {
                        int child_idx = rate_tree.getChildIndex(node_idx,getChildForStep(node_idx,ci));
                        recomputeAtNode(child_idx);
                    }
                }
            }
            recomputeAtNode(rate_tree.getNumNodes()-1); // root
        }        
    }
    
    public String reportTime(int num_calls)
    {
//        StringBuffer sb = new StringBuffer();
//        double tp = (timeP/1000.0)/num_calls;
//        sb.append(".Probabilities.rAN "+tp+" s");
//        double tlc = (timeLC/1000.0)/num_calls;
//        sb.append(",\t.LowerConditionals.rAN "+tlc+" s");
//        return sb.toString();
        return "[StableComputation.reportTime(): Time profiling disabled.]";
    }
    
    private void go (String[] args) throws Exception
    {
        reportLaunch(args);
        reportOtherArguments("Maximum number of paralogs: "+ML.MAX_PARALOGS);
        reportOtherArguments("Minimum number of families: "+ML.MIN_PRESENT_LINEAGES);
        ML.testComputation(new StableComputation(), args);
    }
    
    public static void main(String[] args) throws Exception
    {
        Verbose.setVerbose(false);
                
        StableComputation O = new StableComputation();
        
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
                    ML.MAX_PARALOGS = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("min_lineages"))
                {
                    ML.MIN_PRESENT_LINEAGES = Integer.parseInt(arg_value);
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
