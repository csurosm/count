
package ca.umontreal.iro.evolution.genecontent;


import java.util.Hashtable;

import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.Verbose;

import ca.umontreal.iro.matek.DiscreteDistribution;
import ca.umontreal.iro.matek.Poisson;
import ca.umontreal.iro.matek.NegativeBinomial;
import ca.umontreal.iro.matek.FunctionMinimization;
import ca.umontreal.iro.matek.OneParameterFunction;
import ca.umontreal.iro.matek.DerivableMultiParameterFunction;

import ca.umontreal.iro.evolution.genecontent.HomogeneousRateVariation;

import ca.umontreal.iro.evolution.DiscreteGamma;

/**
 * General likelihood maximization framework for any LikelihoodComputation method.
 *
 * @author csuros
 */
public class ML extends BasicExecutable
{
    /**
     * rate range for loss and gain
     */
    public static double MIN_DUPLICATION_RATE = 1e-5;
    public static double MIN_TRANSFER_RATE = 1e-5;
    public static double MIN_EDGE_LENGTH = 1e-5;
    public static double MAX_DUPLICATION_RATE = 99.;
    public static double MAX_TRANSFER_RATE = 99.;
    public static double MAX_EDGE_LENGTH = 9.;
    public static double MIN_LOSS_RATE = 1e-5;
    public static double MAX_LOSS_RATE = 1e-5;


    /**
     * convergence criterion for Brent's line minimization on edge parameters
     */
    public static double EDGEWISE_BRACKET = 1e-2;

    /**
     * how often are edge rates optimized individually?
     */
    public static int EDGEWISE_OPTIMIZATION_FREQUENCY = 20;
    
    /**
     * Minimum number of lineages in which families are present in the input table.
     * Must be 0,1, or 2.
     */
    public static int MIN_PRESENT_LINEAGES = 1;
    
    /**
     * Maximum total number of paralogs in the input tables
     */
    public static int MAX_PARALOGS = Integer.MAX_VALUE;
    
    /*
     * Number of duplication categories in the input model.
     * -1 means that the input model's specification is used.
     */
    public static int CATEGORY_DUPLICATION = -1;
    
    /*
     * Number of transfer categories in the input model.
     * -1 means that the input model's specification is used.
     */
    public static int CATEGORY_TRANSFER = -1;
    
    /*
     * Number of length categories in the input model.
     * -1 means that the input model's specification is used.
     */
    public static int CATEGORY_EDGE_LENGTH = -1;
    
    /*
     * Number of loss categories in the input model.
     * -1 means that the input model's specification is used.
     */
    public static int CATEGORY_LOSS = -1;
    
    /**
     * Alpha parameter for the underlying discrete gamma model of duplication categories.
     * Double.NaN means that the input model's specification is used.
     */
    public static double ALPHA_DUPLICATION = Double.NaN;
    
    /**
     * Alpha parameter for the underlying discrete gamma model of transfer categories.
     * Double.NaN means that the input model's specification is used.
     */
    public static double ALPHA_TRANSFER = Double.NaN;
    
    /**
     * Alpha parameter for the underlying discrete gamma model of edge length categories.
     * Double.NaN means that the input model's specification is used.
     */
    public static double ALPHA_EDGE_LENGTH = Double.NaN;
    
    /**
     * Alpha parameter for the underlying discrete gamma model of loss categories.
     * Double.NaN means that the input model's specification is used.
     */
    public static double ALPHA_LOSS = Double.NaN;
    
    /**
     * Whether uniform duplication rates are assumed
     */
    public static boolean UNIFORM_DUPLICATION = true;
    
    /**
     * Whether uniform transfer rates are assumed
     */
    public static boolean UNIFORM_TRANSFER = true;
    
    /**
     * Whether uniform edge lengths are assumed
     */
    public static boolean UNIFORM_EDGE_LENGTH = false;
    
    /**
     * Whether there is a <q>forbidden transfer</q> rate category
     */
    public static boolean FORBIDDEN_TRANSFER = false;
    
    /**
     * Whether there is a <q>forbidden duplication</q> rate category
     */
    public static boolean FORBIDDEN_DUPLICATION = false;
    
    
    
    /**
     * Maximum number of rounds in the optimization
     */
    public static int OPTIMIZATION_ROUNDS=100;
    
    /**
     * Convergence criterion for the optimization
     */
    public static double OPTIMIZATION_LL_DELTA=1e-3;

    /**
     * Initializes the class with a LikelihoodComputation model.
     * 
     * @param comp likelihood computation framework
     */
    public ML(LikelihoodComputation comp) 
    {
        init(comp);
    }
    
    /**
     * An argument-less constructor, used in test code.
     * 
     */
    protected ML(){}    
    
    /**
     * Called to initialize the data structures for a given table-model pair
     * and computation model.
     * 
     * @param comp computation model
     */
    private void init(LikelihoodComputation comp)
    {
        this.computation_model= comp;
        init();
    }

    public LikelihoodComputation getComputationModel()
    {
        return computation_model;
    }

    /**
     * Should be called when the computation model is reinitalized.
     */
    private void init()
    {
        this.table = computation_model.getTable();
        this.model = computation_model.getModel(); 
        if (table != null && model !=null)
            initDataStructures();
    }

    /**
     * The computation model implementing the algorithms for manipulating 
     * likelihoods.
     */
    private LikelihoodComputation computation_model;
    
    /**
     * The underlying rate variation model
     */
    private RateVariation model;
    /**
     * The underlying occurrence table
     */
    protected OccurrenceTable table;
    
    /**
     * Main tree (for topology) from the rate variation model
     */
    protected TreeWithRates main_tree;
    
    /**
     * Profiles of the occurence table.
     * Indexing: profiles[pidx]
     */
    protected PhyleticProfile[] profiles;
    
    /**
     * Array of corrected toor priors (one distribution for each rate category)
     * Indexing: corrected_root_priors[cidx]
     */
    protected DiscreteDistribution[] corrected_root_priors;
    
    /**
     * Largest number of paralogs across all families
     */
    protected int maximum_profile_sum;
    
    /**
     * Called to set up the supporting data structures after 
     * changes to the rate variation model or the occurrence table
     */
    protected void initDataStructures()
    {
        uniform_loss = (model instanceof HomogeneousRateVariation);
        main_tree = model.getMainTree();
//        for (int node_idx=0; node_idx<main_tree.getNumNodes(); node_idx++)
//        {
//            System.out.println("#*ML.iDS node "+node_idx+"/"+main_tree.getNode(node_idx).newickName()+"\t"+main_tree.getNode(node_idx));
//        }

        
        profiles = table.getProfiles();
        
        int num_classes = model.getNumClasses();
        corrected_root_priors = new DiscreteDistribution[num_classes];
        
        maximum_profile_sum=0;
        for (int pidx=0; pidx<profiles.length; pidx++)
        {
            int[] pattern = profiles[pidx].getProfile();
            int sum=0; 
            for (int j=0; j<pattern.length; j++)
                sum += pattern[j];
            if (sum>maximum_profile_sum)
                maximum_profile_sum = sum;
        }
        
        //computeSharedProfiles();
    }
    
    /**
     * index of first occurrence for profiles within subtrees: used by computeSharedProfiles
     */
    private int[][] first_occurrence; // indexed as [profile index][node index]
    
    /**
     * experimental code
     * @deprecated
     */
    private void computeSharedProfiles()
    {
        int num_leaves = main_tree.getNumLeaves();
        NodeWithRates[] nodes = main_tree.getDFT();
        first_occurrence = new int[profiles.length][nodes.length];

        // do leaves first
        for (int node_idx=0; node_idx<num_leaves; node_idx++)
        {
            // find maximum at this leaf
            int nmax=0;
            for (int pidx=0; pidx<profiles.length; pidx++)
            {
                int n = profiles[pidx].getProfile()[node_idx];
                if (n>nmax)
                    nmax=n;
            }
            // stre first occurrence of each size
            int[] occurrence = new int[nmax+1];
            for (int n=0;n<occurrence.length;n++) occurrence[n]=-1; // init

            for (int pidx=0; pidx<profiles.length; pidx++)
            {
                int n = profiles[pidx].getProfile()[node_idx];
                if (occurrence[n]==-1) // never seen before
                    occurrence[n]=pidx;
                
                first_occurrence[pidx][node_idx] = occurrence[n];
            }
            
            //System.out.print("#**SPF.cC leaf "+node_idx+"/"+nodes[node_idx].getName());
            //for (int n=0; n<=nmax; n++)
            //    if (occurrence[n]!=-1)
            //        System.out.print("\t"+n+":"+occurrence[n]);
            //System.out.println();
        }
        
        // go through other nodes
        
        for (int node_idx=num_leaves; node_idx<nodes.length; node_idx++)
        {
            Hashtable<String,Integer> mapper = new Hashtable<String,Integer>();
            NodeWithRates N = nodes[node_idx];
            int num_children = N.getNumChildren();
            for (int pidx=0; pidx<profiles.length; pidx++)
            {
                StringBuffer code = new StringBuffer();
                for (int ci=0; ci<num_children; ci++)
                {
                    int occ = first_occurrence[pidx][main_tree.getChildIndex(node_idx,ci)];
                    code.append(Integer.toString(occ));
                    code.append(',');
                }
                String key = code.toString();
                int occ0 = pidx;
                if (mapper.containsKey(key))
                {
                    occ0 = mapper.get(key).intValue();
                } else
                {
                    mapper.put(key,occ0);
                }
                first_occurrence[pidx][node_idx]=occ0;
            }
        }

        int[][] subtree_sizes = new int[nodes.length][profiles.length];
        for (int profile_idx = 0; profile_idx<profiles.length; profile_idx++)
        {
            int[] s = profiles[profile_idx].computeSubtreeSizes(main_tree);
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
                subtree_sizes[node_idx][profile_idx]=s[node_idx];
        }
        
        int total_sq = 0;
        int total_compressed_sq=0;
        for (int node_idx=num_leaves; node_idx<nodes.length; node_idx++)
        {
            int num_different_patterns = 0;
            int[] ndiff = new int[10000];
            int sum_compressed_sq = 0;
            int sum_sq = 0;
            for (int pidx=0; pidx<profiles.length; pidx++)
            {
                int s = subtree_sizes[node_idx][pidx];
                int sq = (s+1)*(s+1);
                if (first_occurrence[pidx][node_idx]==pidx)
                {
                    num_different_patterns++;
                    ndiff[s]++;
                    sum_compressed_sq += sq;
                }
                sum_sq += sq;
            }
            double frac = sum_compressed_sq/((double)sum_sq);
            Verbose.message("ML.cSP node "+node_idx+"/"+nodes[node_idx].getTaxonName()+"\tpatterns "+num_different_patterns+"/"+(profiles.length+1)+"\tsquares "+sum_compressed_sq+"/"+sum_sq+"\t"+frac);
            total_compressed_sq += sum_compressed_sq;
            total_sq += sum_sq;
        }
        double frac = total_compressed_sq / ((double)total_sq);
        Verbose.message("ML.cSP total squares "+total_compressed_sq+"/"+total_sq+"\t"+frac);
    }        

    /**
     * Recomputes the root priors for surviving populations 
     */
    private void setCorrectedRootPriors()
    {
        int num_classes = model.getNumClasses();
        for (int cidx=0; cidx<num_classes; cidx++)
            corrected_root_priors[cidx] = computation_model.getRootPrior(cidx);
    }
    
    private long timeGLL=0L;
    private int callsGLL=0;
    
    /**
     * Computes the log-likelihood of the data set.
     *
     * @return the natural logarithm of the likelihood (of course, a negative number). 
     */
    public double getLogLikelihood()
    {
//        long T0 = System.currentTimeMillis();
        int num_classes = model.getNumClasses();

        double[][] p_root = new double[num_classes][];
        for (int cidx=0; cidx<num_classes; cidx++)
            p_root[cidx] = corrected_root_priors[cidx].getDistribution(maximum_profile_sum);
        
        double ll = 0.0;

        for (int pidx=0; pidx<profiles.length; pidx++ )
        {
            double pL = 0.0;
            for (int cidx=0; cidx<num_classes; cidx++)
                if (model.getClassProbability(cidx)!=0.)
                {
                    double[] L = computation_model.getRootLikelihoods(pidx, cidx);
                    double pc=0.0;
                    for (int i=0; i<L.length; i++)
                    {
                        double w = L[i]*p_root[cidx][i];
                        pc += w;
//                        Verbose.message("ML.gLL profile "+pidx+"\t"+cidx+"\t"+i+"\t"+L[i]+"\t"+p_root[cidx][i]+"\t"+w);
                    }
                    System.out.println("#*ML.gLL profile "+pidx+"\tclass "+cidx+"\tpc "+pc+"\tprob(clas) "+model.getClassProbability(cidx)+"\tpL "+pL);
                    pL += pc*model.getClassProbability(cidx);
                }
            double logpL = 0.0;

            if (pL == 0.0)
            {
                // cannot take logarithm...
                logpL = Math.log(Double.MIN_VALUE)-1;
            } else
                logpL = Math.log(pL);
            Verbose.message("ML.gLL profile "+pidx+"\t"+profiles[pidx].getPatternString()+"\tlogpL "+logpL+"\tll "+ll);
            ll += logpL;
        }
        if (Verbose.isVerbose())
        {
        }
        double p0 = computation_model.getAbsentProfileProbability(MIN_PRESENT_LINEAGES);
//        Verbose.message("ML.gLL profile PRESENT\t"+PhyleticProfile.getPatternString(new int[main_tree.getNumLeaves()])+"\tlogpL "+Math.log1p(-p0)+"\tll "+Double.toString(-profiles.length*Math.log1p(-p0)));
        
        double corr_ll = ll-profiles.length*Math.log(1.-p0);
        
        //Verbose.message("ML.gLL log-likelihood "+ll+" corrected "+corr_ll+"\t(p0 "+p0+")");
        
//        long dT = System.currentTimeMillis()-T0;
//        timeGLL += dT;
//        callsGLL ++;
//        if (callsGLL % 100==0)
//        {
//            double avg = (timeGLL/1000.0)/callsGLL;
//            String x = "";
//            if (computation_model instanceof StableComputation)
//                x = ";\t"+((StableComputation)computation_model).reportTime(callsGLL);
//            Verbose.message("ML.gLL time "+callsGLL+"\tavg "+avg+" seconds"+x);
//        }
        return corr_ll;
    }
    
    
    protected int optimize_maxsteps=100;
    protected double optimize_ll_delta=1e-2;

    /**
     * Sets parameters for stopping the iterations in optimizeRates().
     * @param maxsteps maximum number of steps when rates are not already set (call of optimizeRates(false))
     * @param delta minimum change in ln-likelihood in an iteration step for continuing
     */
    public void setOptimizationParameters(int maxsteps, double delta)
    {
        this.optimize_maxsteps = maxsteps;
        this.optimize_ll_delta = delta;
    }
    
    public int getMaxRounds()
    {
        return optimize_maxsteps;
    }
    
    public double getConvergenceDelta()
    {
        return optimize_ll_delta;
    }
            
    
    /**
     * Here you can fix the edge gain/loss rates (so that they are not 
     * changed in calls to optimizeEdgeRates(int)
     */
    private boolean[] fixed_edge_duplication=null;
    private boolean[] fixed_edge_transfer=null;
    private boolean[] fixed_edge_length=null;
    private boolean fixed_root_prior=false;
    private boolean uniform_duplication = true;//false;
    private boolean uniform_transfer = true;//false;
    private boolean uniform_edge_length = false;//false;
    private boolean uniform_loss = false;
    
    
    private boolean fixed_duplication_gamma = false;
    private boolean fixed_length_gamma = false;
    private boolean fixed_transfer_gamma = false;
    private boolean fixed_loss_gamma = false;

    /**
     * Sets whether the same parameters are applied on all edges
     * @param uniform_duplication all duplication rates are the same on all edges
     * @param uniform_transfer all transfer rates are the same on all edges
     * @param uniform_edge_length all edge lengths are the same 
     */
    public void setUniformEdgeParameters(
            boolean uniform_duplication,
            boolean uniform_transfer, 
            boolean uniform_edge_length)
    {
        this.uniform_duplication = uniform_duplication;
        this.uniform_transfer = uniform_transfer;
        this.uniform_edge_length = uniform_edge_length;
    }


    /**
     * Excludes certain parameters from the optimization. 
     *
     * @param fixed_root_prior should be true iff the prior distribution is not to be optimized
     * @param fixed_duplication_rates fdr[node_idx] should be true iff the duplication rate on the branch leading to node_idx is not to be optimized
     * @param fixed_transfer_rates ftr[node_idx] should be true iff the transfer rate on the branch leading to node_idx is not to be optimized
     * @param fixed_lengths fl[node_idx] should be true iff the edge length on the branch leading to node_idx is not to be optimized
     */
    public void setOptimizationFixedParameters(boolean fixed_root_prior, boolean[] fixed_duplication_rates, boolean[] fixed_transfer_rates, boolean[] fixed_lengths)
    {
        this.fixed_edge_duplication = fixed_duplication_rates;
        this.fixed_edge_transfer = fixed_transfer_rates;
        this.fixed_edge_length = fixed_lengths;
        this.fixed_root_prior = fixed_root_prior;
    }
    
    /**
     * Whether gain (transfer) rates are potentially non-zero
     * 
     * @return true if either the underlying tree has gain rates, or they are optimized
     */
    public boolean hasGain()
    {
        boolean has_gain = main_tree.hasGain();
        if (!has_gain && fixed_edge_transfer != null)
        {
            for (int i=0; i<fixed_edge_transfer.length && !has_gain; i++)
                has_gain = has_gain || !fixed_edge_transfer[i];
        }
        return has_gain;
    }

    /**
     * Whether duplication rates are potentially non-zero
     * 
     * @return true if either the underlying tree has duplication rates, or they are optimized
     */
    public boolean hasDuplication()
    {
        boolean has_duplication = main_tree.hasDuplication();
        if (!has_duplication && fixed_edge_duplication != null)
        {
            for (int i=0; i<fixed_edge_duplication.length && !has_duplication; i++)
                has_duplication = has_duplication || !fixed_edge_duplication[i];
        }
        return has_duplication;
    }
    
    
    /**
     * Excudes rate variation parameters from the optimization.
     * 
     * @param fixed_length_gamma whether shape parameter for edge length is excluded (only consulted if there is at least 1 category)
     * @param fixed_duplication_gamma whether shape parameter for duplication rate is excluded (only consulted if there is at least 1 category)
     * @param fixed_loss_gamma whether shape parameter for loss rate is excluded (only consulted if there is at least 1 category)
     * @param fixed_transfer_gamma whether shape parameter for transfer (gain)  rate is excluded (only consulted if there is at least 1 category)
     */
    public void setOptimizationFixedVariation(boolean fixed_length_gamma, boolean fixed_duplication_gamma, boolean fixed_loss_gamma, boolean fixed_transfer_gamma)
    {
        this.fixed_length_gamma = fixed_length_gamma;
        this.fixed_loss_gamma = fixed_loss_gamma;
        this.fixed_duplication_gamma = fixed_duplication_gamma;
        this.fixed_transfer_gamma = fixed_transfer_gamma;
    }
    
    
    private boolean forbidden_transfer = false;
    private boolean forbidden_duplication = false;
    
    /**
     * Exludes the optimization of forbidden proportions
     * 
     * @param forbidden_transfer whether proportion of families with no transfer is optimized
     * @param forbidden_duplication whether proportion of families with no duplication is optimized
     */
    public void setForbiddenCategories(boolean forbidden_transfer, boolean forbidden_duplication)
    {
        this.forbidden_transfer = forbidden_transfer;
        this.forbidden_duplication = forbidden_duplication;
    }
    
    private int optimization_round;
    private String optimization_stage;
    private double current_loglikelihood_value;
    private double current_likelihood_drop;
    
    /**
     * Returns the optimization round; useful for progress indicators
     * 
     * @return the optimization round currently executed
     * 
     */
    public synchronized int getOptimizationRound()
    {
        return optimization_round;
    }

    /**
     * Returns a descriptive message about the optimization stage; useful for progress indicators
     * @return 
     */
    public synchronized String getOptimizationStage()
    {
        return optimization_stage;
    }
    
    /**
     * Current value of the log-likelihood; useful for progress indicators
     * @return current value of the the log-likelihood
     */
    public synchronized double getOptimizationCurrentLogLikelihood()
    {
        return current_loglikelihood_value;
    }
    
    /**
     * Current value of the change in the log-likelihood; useful for progress indicators
     * @return current value of the the change in log-likelihood
     */
    public synchronized double getOptimizationCurrentLikelihoodDrop()
    {
        return current_likelihood_drop;
    }
    
    /**
     * Counts the number of true entries in a boolean array
     */
    private static int countTrue(boolean[] array)
    {
        int nt = 0;
        if (array != null)
        {
            for (int i=0; i<array.length; i++)
                if (array[i])
                    nt++;
        }
        return nt;
    }
    
    /** 
     * Computes the posterior distribution for number of surviving members at each node.
     * The computed array gives the probabilities P[u][m] that at node u, there are 
     * m individuals that have descendants within u's subtree, where the probability is conditioned 
     * on the observed profile. 
     * 
     * @param profile_idx index of the profile
     * @return array of distributions for this profile
     */
    //public abstract double[][] getPosteriorSurvivingCount(int profile_idx);
    
    /**
     * Computes the posterior probability for belonging in the classes
     * 
     * @param profile_idx index of the profile
     * @return array of posterior probabilities: one entry for each rate category 
     */
    public double[] getPosteriorClassProbabilities(int profile_idx)
    {
        PhyleticProfile PP = profiles[profile_idx];
        int root_idx  = main_tree.getNumNodes()-1; // last node in postorder traversal

        int num_classes = model.getNumClasses();
        double[] p_class = new double[num_classes];
        double p_normalize = 0.0;
        
        for (int cidx=0; cidx<num_classes; cidx++)
        {
            double[] L = computation_model.getRootLikelihoods(profile_idx,cidx);
            double[] p_root = corrected_root_priors[cidx].getDistribution(L.length);
            double pc=0.0;
            for (int i=0; i<L.length; i++)
            {
                double w = L[i]*p_root[i];
                pc += w;
            }
            p_normalize += pc;
            p_class[cidx]=pc; // class prior distribution is uniform ...
        }
        for (int cidx=0; cidx<num_classes; cidx++)
            p_class[cidx] /= p_normalize;
        
        return p_class;
    }
    
    private void recomputeAll()
    {
        for (int cidx=0; cidx<model.getNumClasses(); cidx++)
            computation_model.recomputeCategorySupport(cidx);
        
        setCorrectedRootPriors();
        
        for (int pidx=0; pidx<profiles.length; pidx++ )
            for (int cidx=0; cidx<model.getNumClasses(); cidx++)
                computation_model.recomputeProfileSupport(pidx,cidx);
    }

    private void recomputeOnPathToRoot(int node_idx)
    {
        for (int cidx=0; cidx<model.getNumClasses(); cidx++)
            computation_model.recomputeCategorySupport(cidx, node_idx);
        
        setCorrectedRootPriors();

        for (int pidx=0; pidx<profiles.length; pidx++)
            for (int cidx=0; cidx<model.getNumClasses(); cidx++)
                computation_model.recomputeProfileSupport(pidx, cidx, node_idx);
    }
    
    
    /**
     * @return the negative log-likelihood, or NaN if interrupted (Thread.interrupt)
     */
    public double optimize()
    {
        recomputeAll();

        int num_nodes = main_tree.getNumNodes();        
        int max_rounds = optimize_maxsteps;
        double opt_ll=-getLogLikelihood();
        
        current_loglikelihood_value = -opt_ll;
        optimization_stage = "Optimization start";
        optimization_round = 0;
        current_likelihood_drop = 0.;
        
        double previous_delta = Double.POSITIVE_INFINITY;
        
        for (int round=0; round<max_rounds; round++)
        {
            optimization_round = round+1;

            double delta = 0.0;
            
            // optimize edge parameters
            boolean do_edgewise_optimization 
                    = (previous_delta>100.0 || round%EDGEWISE_OPTIMIZATION_FREQUENCY==0); // get good starting values using line minimization - it's safer because derivatives are not used;
            if (!do_edgewise_optimization)
            {
                optimization_stage = "Setting branch parameters by multidimensional optimization";

                double new_ll = optimizeEdgeParameters(opt_ll);
                //
                // check interrupt here
                //
                if (Thread.interrupted())
                {
                    return Double.NaN;
                }

                if (Double.isNaN(new_ll))
                {
                    Verbose.message("ML.oEP dfpmin unfeasible");
                    //recomputeAll();
                    do_edgewise_optimization = true;
                    //System.out.println("======================\n\n\n"+ model.tableRates());
                    //System.exit(666);
                } else
                {
                    Verbose.message("ML.o edges \t"+opt_ll+" ->\t"+new_ll);
                    delta += opt_ll-new_ll;
                    opt_ll = new_ll;

                    current_loglikelihood_value = -opt_ll;
                }
            }

            
            // uniform rates on the edges
            if (uniform_duplication || uniform_transfer || uniform_edge_length || uniform_loss)
            {
                optimization_stage = "Optimizing uniform edge parameters";
                double new_ll = optimizeUniformRates(opt_ll);
                delta += opt_ll - new_ll;
                opt_ll = new_ll;
                
                current_loglikelihood_value = -opt_ll;
            } 
            
            if (do_edgewise_optimization)
            {
                for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // not at root
                {
                    NodeWithRates N = main_tree.getNode(node_idx);
                    optimization_stage = "Optimizing rates on branch "+N.newickName();
                    double new_ll=optimizeEdgeParameters(node_idx);

                    //
                    // check interrupt here
                    //
                    if (Thread.interrupted())
                    {
                        return Double.NaN;
                    }

                    if (new_ll != 0.0)
                    {
                        Verbose.message("ML.o edge "+node_idx+"/"+N.getTaxonName()+"\t\t"+opt_ll+" ->\t"+new_ll);
                        delta += opt_ll-new_ll;
                        opt_ll=new_ll;

                        current_loglikelihood_value = -opt_ll;
                    }
                }
            }
            
            if (forbidden_transfer || forbidden_duplication)
            {
                optimization_stage = "Optimizing forbidden rate categories";
                double new_ll = optimizeForbiddenCategories(opt_ll);
                delta += opt_ll - new_ll;
                opt_ll = new_ll;
                
                current_loglikelihood_value = -opt_ll;
            }
            
            //
            // check interrupt here
            //
            if (Thread.interrupted())
            {
                return Double.NaN;
            }

            // optimize root prior
            if (!fixed_root_prior)
            {
                optimization_stage = "Optimizing prior size distribution at root.";
                double new_ll = optimizeRootPrior();
                Verbose.message("ML.o root \t"+opt_ll+" ->\t"+new_ll);
                delta += opt_ll-new_ll;
                opt_ll=new_ll;
                
                current_loglikelihood_value = -opt_ll;
            }
            
            //
            // check interrupt here
            //
            if (Thread.interrupted())
            {
                return Double.NaN;
            }

            // optimize classes
            if (model.getNumClasses()>1)
            {
                optimization_stage = "Optimizing rate variation parameters";
                double new_ll = optimizeRateVariation(opt_ll);
                delta += opt_ll - new_ll;
                opt_ll = new_ll;
                
                current_loglikelihood_value = -opt_ll;
            }

            //
            // check interrupt here
            //
            if (Thread.interrupted())
            {
                return Double.NaN;
            }
            
            
            current_likelihood_drop = delta;
            
            Verbose.message("ML.o round "+(1+round)+" ll "+opt_ll+" delta "+delta);
            if (delta<optimize_ll_delta && previous_delta<optimize_ll_delta)
                    break;
            previous_delta = delta;
            
            if (Double.isNaN(opt_ll))
                throw new LikelihoodOptimizationException("Numerical error in ML.optimize()");
        }
        
        optimization_round = max_rounds;
        optimization_stage = "Optimization done";
        
        double retval = -getLogLikelihood();
        return retval;
    }
    
    /**
     * Multidimensional optimization for all edge parameters
     * @param opt_ll starting value of the log-likelihood
     * @return value of log-likelihood after optimizing all edge parameters
     */
    protected double optimizeEdgeParameters(double opt_ll)
    {
        OptimizeAllEdges O = new OptimizeAllEdges();
        // so that steps are not too big: the starting point may have a very steep gradient
        // In FunctionMinimization, dfpmin uses line search 
        // in some direction: initial values are DFP_STMAX*scale,
        // where scale is either the Euclidean length of the parameter vector 
        // (in our case, all parameters are rates and lengths, 
        // and their logarithm is used in the vector), or the number of 
        // parameters optimized together.
        FunctionMinimization.DFP_STPMX = 2.0; 
        FunctionMinimization.DFP_ITMAX = 2*main_tree.getNumNodes();
        try 
        {
            double[] x = O.getCurrentValue(true);
            double[] starting_values = O.getCurrentValue(false);
            double new_ll = Double.NaN;
            try 
            {
                new_ll=FunctionMinimization.dfpmin(x, optimize_ll_delta, O);
            } catch (LikelihoodOptimizationException E)
            {
                Verbose.message("ML.oEP caught "+E);
                if (O.best_x != null && O.best_ll < opt_ll-optimize_ll_delta)
                {
                    Verbose.message("ML.oEP using "+O.best_ll);
                    O.set(O.best_x, true);
                    return O.best_ll;
                } else
                {
                    O.set(starting_values, false);
                    Verbose.message("ML.oEP resetting original values");
                    return Double.NaN;
                }
            }
            O.set(x, true);
            Verbose.message("ML.OEP dfpmin OK "+new_ll+"\twas "+opt_ll+"\tbest "+O.best_ll+"\tcalls "+O.num_eval_calls
                    +(O.best_ll<new_ll?" \tSTOP AT INFERIOR":""));
            if (opt_ll<=new_ll)
            {
                if (O.best_x != null && O.best_ll < opt_ll-optimize_ll_delta)
                {
                    Verbose.message("ML.oEP using "+O.best_ll);
                    O.set(O.best_x, true);
                    return O.best_ll;
                } else
                    return Double.NaN;
            }
            opt_ll=new_ll;
        } catch (FunctionMinimization.OptimizationException E)
        {
            Verbose.message("ML.oEP   **** dfpmin interrupted ");
            opt_ll = -getLogLikelihood();
        }
        return opt_ll;
    }
    
    /**
     * Maximizes the likelihood in function of rates and length for a particular edge
     * @param node_idx child node for the edge
     * @return value of the log-likelihood after optimization
     */
    protected double optimizeEdgeParameters(int node_idx)
    {
        double opt_ll=0.;
        boolean optimize_duplication = !uniform_duplication && (fixed_edge_duplication ==null || !fixed_edge_duplication[node_idx]);
        boolean optimize_transfer = !uniform_transfer && (fixed_edge_transfer == null || !fixed_edge_transfer[node_idx]);
        boolean optimize_length = !uniform_edge_length && (fixed_edge_length == null || !fixed_edge_length[node_idx]);
        NodeWithRates N =main_tree.getNode(node_idx);
        double old_length = N.getLength();
        double old_duplication = N.getDuplicationRate();
        double old_transfer = N.getTransferRate();
        
        if (optimize_length)
        {
            double len = N.getLength();
            OptimizeEdgeLength O = new OptimizeEdgeLength(node_idx);
            double l0 = getLogLikelihood();
            O.set(len);
            double l1 = getLogLikelihood();
            double[] opt = FunctionMinimization.brent(MIN_EDGE_LENGTH, len, MAX_EDGE_LENGTH, O, EDGEWISE_BRACKET);
            //Verbose.message("ML.oEP "+node_idx+"/"+N.getTaxonName()+" length opt_ll "+opt[1]+"\tlen -> "+opt[0]+"\t["+MIN_EDGE_LENGTH+".."+MAX_EDGE_LENGTH+"]\tl0 "+l0+"\tl1 "+l1);       
            O.set(opt[0]);
            opt_ll = opt[1];
        }
        if (optimize_duplication)
        {
            double dup = N.getDuplicationRate();
            OptimizeEdgeDuplication O = new OptimizeEdgeDuplication(node_idx);
            double[] opt = FunctionMinimization.brent(MIN_DUPLICATION_RATE, dup, MAX_DUPLICATION_RATE, O, EDGEWISE_BRACKET);
            //Verbose.message("ML.oEP "+node_idx+"/"+N.getTaxonName()+" duplication opt_ll "+opt[1]);       
            O.set(opt[0]);
            opt_ll = opt[1];
            
        }
        if (optimize_transfer)
        {
            double tra = N.getTransferRate();
            OptimizeEdgeTransfer O = new OptimizeEdgeTransfer(node_idx);
            double[] opt = FunctionMinimization.brent(MIN_TRANSFER_RATE, tra, MAX_TRANSFER_RATE, O, EDGEWISE_BRACKET);
            //Verbose.message("ML.oEP "+node_idx+"/"+N.getTaxonName()+" transfer opt_ll "+opt[1]);       
            O.set(opt[0]);
            opt_ll = opt[1];
            
        }
        
        if (Verbose.isVerbose() && (optimize_length || optimize_duplication || optimize_transfer))
        {
            double new_length = N.getLength();
            double new_duplication = N.getDuplicationRate();
            double new_transfer = N.getTransferRate();


            double lpct = ((int)(1000.0*(new_length-old_length)/old_length+0.5))/10.0;
            double dpct = ((int)(1000.0*(new_duplication-old_duplication)/old_duplication+0.5))/10.0;
            double tpct = ((int)(1000.0*(new_transfer-old_transfer)/old_transfer+0.5))/10.0;

            Verbose.message("ML.oEP "+node_idx+"/"+N.getTaxonName()+"\t(l"+lpct+"%, d"+dpct+"%, t"+tpct+"%)" 
                    + "\t(l "+old_length+", d "+old_duplication+", t "+old_transfer+") " 
                    + "->\t("+new_length+", "+new_duplication+", "+new_transfer+")");

        }
        return opt_ll;
    }
    
    protected double optimizeRootPrior()
    {
        double opt_ll = -getLogLikelihood();
        DiscreteDistribution D = model.getRootPrior();
        int num_params = D.getNumParameters();
        
        for (int round=0; round<(num_params==1?1:num_params*2); round++)
        {
            int pidx = round % num_params;
            double original_value = D.getParameters()[pidx];
            double pmin = 0.0;
            double pmax = 1.0;        
            if (D instanceof Poisson)
            {
                pmin=0.001;
                pmax=20.0;
            } else if (D instanceof NegativeBinomial && pidx==0) 
            {
                pmin=0.001;
                pmax=100.0;
            }
            OptimizeRootPrior O = new OptimizeRootPrior(pidx);
            double[] opt = FunctionMinimization.brent(pmin,original_value,pmax,O,1e-5 );
            O.set(opt[0]);
            Verbose.message("ML.oRP root "+D.getClass().getSimpleName()+"["+pidx+"]/"+round+"\t"+opt_ll+" ->\t"+opt[1]+"\t"+original_value+" ->\t"+opt[0]);
            opt_ll = opt[1];
        }
        return opt_ll;
    }
    
    protected double optimizeUniformRates(double opt_ll)
    {
        
        NodeWithRates N0 =main_tree.getNode(0);
        double old_length = N0.getLength();
        double old_duplication = N0.getDuplicationRate();
        double old_transfer = N0.getTransferRate();
        String info_actions="";
                
        if (uniform_edge_length)
        {
            double len = N0.getLength();
            OptimizeEdgeLength O = new OptimizeEdgeLength(-1);
            double[] opt = FunctionMinimization.brent(MIN_EDGE_LENGTH, len, MAX_EDGE_LENGTH, O, EDGEWISE_BRACKET);
            O.set(opt[0]);
            info_actions = "\tlen "+len+" ->\t"+opt[0]+"\t change "+(opt_ll-opt[1]);
            opt_ll = opt[1];
        }
        if (uniform_duplication)
        {
            double dup = N0.getDuplicationRate();
            OptimizeEdgeDuplication O = new OptimizeEdgeDuplication(-1);
            double[] opt = FunctionMinimization.brent(MIN_DUPLICATION_RATE, dup, MAX_DUPLICATION_RATE, O, EDGEWISE_BRACKET);
            O.set(opt[0]);
            info_actions = info_actions+"\tdup "+dup+" ->\t"+opt[0]+"\t change "+(opt_ll-opt[1]);
            opt_ll = opt[1];
        }
        if (uniform_transfer)
        {
            double tra = N0.getTransferRate();
            OptimizeEdgeTransfer O = new OptimizeEdgeTransfer(-1);
            double[] opt = FunctionMinimization.brent(MIN_TRANSFER_RATE, tra, MAX_TRANSFER_RATE, O, EDGEWISE_BRACKET);
            O.set(opt[0]);
            info_actions = info_actions+"\ttra "+tra+" ->\t"+opt[0]+"\t change "+(opt_ll-opt[1]);
            opt_ll = opt[1];
        }
        if (uniform_loss)
        {
            double ls = N0.getLossRate();
            OptimizeEdgeLoss O = new OptimizeEdgeLoss();
            double[] opt = FunctionMinimization.brent(MIN_LOSS_RATE, ls, MAX_LOSS_RATE, O, EDGEWISE_BRACKET);
            info_actions = "\tloss "+ls+" ->\t"+opt[0]+"\t change "+(opt_ll-opt[1]);
            O.set(opt[0]);
            opt_ll = opt[1];
        }
        
        if (Verbose.isVerbose())
        {
            double new_length = N0.getLength();
            double new_duplication = N0.getDuplicationRate();
            double new_transfer = N0.getTransferRate();

//            double lpct = ((int)(1000.0*(new_length-old_length)/old_length+0.5))/10.0;
//            double dpct = ((int)(1000.0*(new_duplication-old_duplication)/old_duplication+0.5))/10.0;
//            double tpct = ((int)(1000.0*(new_transfer-old_transfer)/old_transfer+0.5))/10.0;

            Verbose.message("ML.oUR"+info_actions);
        }
        return opt_ll;
    }
    
    protected double optimizeForbiddenCategories(double opt_ll)
    {
        if (forbidden_duplication)
        {
            double p0 = model.getDuplicationForbiddenProportion();
            optimization_stage = "Optimizing proportion of families with forbidden duplications";
            OptimizeForbiddenDuplication O = new OptimizeForbiddenDuplication();
            double[] opt = FunctionMinimization.brent(0.0, p0, 1.0, O, 1e-5);
            O.set(opt[0]);
            Verbose.message("ML.oFC dup  "+opt_ll+" ->\t"+opt[1]+"\t"+p0+" ->\t"+opt[0]);
            opt_ll = opt[1];
        }
        if (forbidden_transfer)
        {
            double p0 = model.getDuplicationForbiddenProportion();
            optimization_stage = "Optimizing proportion of families with forbidden gains";
            OptimizeForbiddenTransfer O = new OptimizeForbiddenTransfer();
            double[] opt = FunctionMinimization.brent(0.0, p0, 1.0, O, 1e-5);
            O.set(opt[0]);
            Verbose.message("ML.oFC trans "+opt_ll+" ->\t"+opt[1]+"\t"+p0+" ->\t"+opt[0]);
            opt_ll = opt[1];
        }
        
        return opt_ll;
    }
    
    protected double optimizeRateVariation(double opt_ll)
    {
        if (model.getNumDuplicationRateGammaCategories() > 1 && !fixed_duplication_gamma)
        {
            double a0 = model.getDuplicationRateAlpha();
            optimization_stage = "Optimizing variation of duplication rates";
            OptimizeDuplicationRateGamma O = new OptimizeDuplicationRateGamma();
            double[] opt = FunctionMinimization.brent(0.0, a0, DiscreteGamma.RECOMMENDED_MAXIMUM_ALPHA, O, 1e-5);
            O.set(opt[0]);
            Verbose.message("ML.oRV dup "+opt_ll+" ->\t"+opt[1]+"\t"+a0+" ->\t"+opt[0]+"\tk= "+model.getNumDuplicationRateGammaCategories());
            opt_ll = opt[1];
        }
        if (model.getNumTransferRateGammaCategories() > 1 && !fixed_transfer_gamma)
        {
            double a0 = model.getTransferRateAlpha();
            optimization_stage = "Optimizing variation of transfer rates";
            OptimizeTransferRateGamma O = new OptimizeTransferRateGamma();
            double[] opt = FunctionMinimization.brent(0.0, a0, DiscreteGamma.RECOMMENDED_MAXIMUM_ALPHA, O, 1e-5);
            O.set(opt[0]);
            Verbose.message("ML.oRV tran "+opt_ll+" ->\t"+opt[1]+"\t"+a0+" ->\t"+opt[0]+"\tk= "+model.getNumTransferRateGammaCategories());
            opt_ll = opt[1];
        }
        if (model.getNumLossRateGammaCategories() > 1 && !fixed_loss_gamma)
        {
            double a0 = model.getLossRateAlpha();
            optimization_stage = "Optimizing variation of loss rates";
            OptimizeLossRateGamma O = new OptimizeLossRateGamma();
            double[] opt = FunctionMinimization.brent(0.0, a0, DiscreteGamma.RECOMMENDED_MAXIMUM_ALPHA, O, 1e-5);
            O.set(opt[0]);
            Verbose.message("ML.oRV loss "+opt_ll+" ->\t"+opt[1]+"\t"+a0+" ->\t"+opt[0]+"\tk= "+model.getNumLossRateGammaCategories());
            opt_ll = opt[1];
        }
        if (model.getNumEdgeLengthGammaCategories() > 1 && !fixed_length_gamma)
        {
            double a0 = model.getEdgeLengthAlpha();
            optimization_stage = "Optimizing variation of edge lengths";
            OptimizeEdgeLengthGamma O = new OptimizeEdgeLengthGamma();
            double[] opt = FunctionMinimization.brent(0.0, a0, DiscreteGamma.RECOMMENDED_MAXIMUM_ALPHA, O, 1e-5);
            O.set(opt[0]);
            Verbose.message("ML.oRV length "+opt_ll+" ->\t"+opt[1]+"\t"+a0+" ->\t"+opt[0]);
            opt_ll = opt[1];
        }
        return opt_ll;
    }

    class OptimizeEdgeLoss implements OneParameterFunction
    {
        /**
         * Initializes a new object for line minimization of uniform loss rate.
         */
        OptimizeEdgeLoss()
        {
        }
        private int num_calls=0;

        @Override
        public double eval(double dup)
        {
            set(dup);
            double ll = getLogLikelihood();
            num_calls++;
            return -ll;
        }

        public int getNumCalls()
        {
            return num_calls;
        }

        public void set(double rate)
        {
            assert (model instanceof HomogeneousRateVariation);
            model.setNodeLossRate(0, rate);
            recomputeAll();
        }

    }
           
    class OptimizeEdgeDuplication implements OneParameterFunction 
    {
        /**
         * Initializes a new object for line minimization of duplication rate. 
         *
         * @param node_idx index for the bottom node on the branch, or -1 if the same rate is imposed on all branches
         */ 
        OptimizeEdgeDuplication(int node_idx)
        {
            this.node_idx = node_idx;
        }
        private int node_idx;
        private int num_calls=0;

        @Override
        public double eval(double dup) 
        {
            set(dup);
            double ll = getLogLikelihood();
            num_calls++;
            return -ll;
        }
        
        public int getNumCalls()
        {
            return num_calls;
        }
        
        public void set(double rate)
        {
            if (node_idx != -1)
            {
                model.setNodeDuplicationRate(node_idx, rate);
                recomputeOnPathToRoot(node_idx);
                //double l1=getLogLikelihood();
                //recomputeAll();
                //double l2 = getLogLikelihood();
                //Verbose.message("ML.OED.s ["+node_idx+"] "+rate+"\tl1 "+l1+"\tl2 "+l2);
            } else
            {
                for (int j=0; j<main_tree.getNumNodes()-1; j++)
                    if (fixed_edge_duplication == null || !fixed_edge_duplication[j] )                    
                        model.setNodeDuplicationRate(j, rate);
                recomputeAll();
            }
        }
        
    }
         
    class OptimizeEdgeTransfer implements OneParameterFunction 
    {
        /**
         * Initializes a new object for line minimization of transfer rate. 
         *
         * @param node_idx index for the bottom node on the branch, or -1 if the same rate is imposed on all branches
         */ 
        OptimizeEdgeTransfer(int node_idx)
        {
            this.node_idx = node_idx;
        }
        private int node_idx;
        private int num_calls=0;
        
        @Override
        public double eval(double dup) 
        {
            set(dup);
            double ll = getLogLikelihood();
            num_calls++;
            return -ll;
        }
        
        public int getNumCalls()
        {
            return num_calls;
        }
        
        public void set(double rate)
        {
            if (node_idx != -1)
            {
                model.setNodeTransferRate(node_idx, rate);
                recomputeOnPathToRoot(node_idx);
            } else
            {
                for (int j=0; j<main_tree.getNumNodes()-1; j++)
                    if (fixed_edge_transfer == null || !fixed_edge_transfer[j] )
                        model.setNodeTransferRate(j, rate);
                recomputeAll();
            }
        }
    }
    
    class OptimizeEdgeLength implements OneParameterFunction 
    {
        /**
         * Initializes a new object for line minimization of edge length. 
         *
         * @param node_idx index for the bottom node on the branch, or -1 if the same rate is imposed on all branches
         */ 
        OptimizeEdgeLength(int node_idx)
        {
            this.node_idx = node_idx;
        }
        private int node_idx;
        private int num_calls=0;
        
        @Override
        public double eval(double length) 
        {
            double l0 = getLogLikelihood();
            set(length);
            double ll = getLogLikelihood();
            num_calls++;
            //Verbose.message("ML.OEL.e "+node_idx+"\tlength "+length+"\tll "+ll+"\tcalls "+num_calls+"\t// was "+l0);
            return -ll;
        }
        
        public int getNumCalls()
        {
            return num_calls;
        }
        
        public void set(double rate)
        {
            if (node_idx != -1)
            {
                model.setEdgeLength(node_idx, rate);
                recomputeOnPathToRoot(node_idx);
            } else
            {
                for (int j=0; j<main_tree.getNumNodes()-1; j++)
                    if (fixed_edge_length == null || !fixed_edge_length[j] )
                        model.setEdgeLength(j, rate);
                recomputeAll();
            }
        }
    }

    class OptimizeDuplicationRateGamma implements OneParameterFunction
    {
        public void set (double alpha)
        {
            model.setDuplicationRateAlpha(alpha);
            recomputeAll();
        }

        @Override
        public double eval(double alpha)
        {
            set(alpha);
            double ll = getLogLikelihood();
            return -ll;
        }
    }
    
    class OptimizeLossRateGamma implements OneParameterFunction
    {
        public void set (double alpha)
        {
            model.setLossRateAlpha(alpha);
            recomputeAll();
        }
        
        @Override
        public double eval(double alpha)
        {
            set(alpha);
            double ll = getLogLikelihood();
            return -ll;
        }
    }

    
    class OptimizeTransferRateGamma implements OneParameterFunction
    {
        public void set (double alpha)
        {
            model.setTransferRateAlpha(alpha);
            recomputeAll();
        }
        
        @Override
        public double eval(double alpha)
        {
            set(alpha);
            double ll = getLogLikelihood();
            return -ll;
        }
    }
    
    class OptimizeEdgeLengthGamma implements OneParameterFunction
    {
        public void set (double alpha)
        {
            model.setEdgeLengthAlpha(alpha);
            recomputeAll();
        }
        
        @Override
        public double eval(double alpha)
        {
            set(alpha);
            double ll = getLogLikelihood();
            return -ll;
        }
    }
    
    

    private static final int PARAMETER_INDEX_DUPLICATION=0;
    private static final int PARAMETER_INDEX_TRANSFER=1;
    private static final int PARAMETER_INDEX_LENGTH=2;
    
    class OptimizeAllEdges implements DerivableMultiParameterFunction
    {
        
        OptimizeAllEdges()
        {
            initMapping();
            //double l0 = getLogLikelihood();
            //set(getCurrentValue());
            //recomputeAll();
            //double l1 = getLogLikelihood();
            //Verbose.message("#**ML.OAE() "+l0+"\t"+l1);
        }
        
        private void initMapping()
        {
            int num_edges = main_tree.getNumNodes()-1;
            int num_fixed_duplication=countTrue(fixed_edge_duplication);
            if (uniform_duplication)
                num_fixed_duplication = num_edges;
            int num_fixed_transfer=countTrue(fixed_edge_transfer);
            if (uniform_transfer)
                num_fixed_transfer = num_edges;
            int num_fixed_length=countTrue(fixed_edge_length);
            if (uniform_edge_length)
                num_fixed_length = num_edges;
            int num_parameters = num_edges*3;
            parameter_map = new int[num_parameters-(num_fixed_duplication+num_fixed_transfer+num_fixed_length)];
            last_x_set = new double[parameter_map.length];
            best_x = null;
            int parameter_idx=0;
            for (int i=0; i<num_edges; i++)
                if (!uniform_duplication && (num_fixed_duplication==0 || !fixed_edge_duplication[i]))
                {
                    parameter_map[parameter_idx] = i*3+PARAMETER_INDEX_DUPLICATION;
                    parameter_idx++;
                }
            for (int i=0; i<num_edges; i++)
                if (!uniform_transfer && (num_fixed_transfer==0 || !fixed_edge_transfer[i]))
                {
                    parameter_map[parameter_idx] = i*3+PARAMETER_INDEX_TRANSFER;
                    parameter_idx++;
                }
            for (int i=0; i<num_edges; i++)
                if (!uniform_edge_length && (num_fixed_length==0 || !fixed_edge_length[i]))
                {
                    parameter_map[parameter_idx] = i*3+PARAMETER_INDEX_LENGTH;
                    parameter_idx++;
                }
        }
        
        private int[] parameter_map;
        
        /**
         * A vector with all non-fixed duplication and transfer rates and edge lengths that can be used as
         * starting point for the optimization
         */
        public double[] getCurrentValue(boolean do_transformation)
        {
            double[] x =new double[parameter_map.length];
            for (int param_idx=0; param_idx<parameter_map.length; param_idx++)
            {
                int node_idx = parameter_map[param_idx] / 3;
                int parameter_type = parameter_map[param_idx] % 3;
                NodeWithRates N=main_tree.getNode(node_idx);
                double v = 0.0;
                if (parameter_type == PARAMETER_INDEX_DUPLICATION)
                {
                    v = N.getDuplicationRate();
                }
                else if (parameter_type == PARAMETER_INDEX_TRANSFER)
                    v = N.getTransferRate();
                else 
                    v =N.getLength();
                if (do_transformation)
                    v = Math.log(v);
                //Verbose.message("ML.OAE.gCV ["+param_idx+"] "+v+"\t"+node_idx+":"+parameter_type);
                x[param_idx] =v;
            }

            //System.arraycopy(x, 0, last_x_set, 0, x.length);
            
            return x;
        }
        
        private double[] last_x_set;
        
        public void set(double[] x, boolean do_transformation)
        {
            for (int param_idx=0; param_idx<parameter_map.length; param_idx++)
            {
                double value_set = 999.999;
                int node_idx = parameter_map[param_idx] / 3;
                int parameter_type = parameter_map[param_idx] % 3;
                if (parameter_type == PARAMETER_INDEX_DUPLICATION)
                    value_set = setNodeDuplicationRate(node_idx, x[param_idx], do_transformation);
                else if (parameter_type == PARAMETER_INDEX_TRANSFER)
                    value_set = setNodeTransferRate(node_idx,x[param_idx], do_transformation);
                else
                    value_set = setEdgeLength(node_idx,x[param_idx], do_transformation);
                //Verbose.message("ML.OAE.s ["+param_idx+"] "+x[param_idx]+"\t"+node_idx+":"+parameter_type+"\t"+value_set+"\t"+main_tree.getNode(node_idx));
            }
            
            recomputeAll();
                        
            System.arraycopy(x, 0, last_x_set, 0, x.length);
        }
        
        private double setNodeDuplicationRate(int node_idx, double param_x, boolean do_transformation)
        {
            if (do_transformation)
                param_x = Math.exp(param_x);
            //Verbose.message("ML.OAE.sNDR ["+node_idx+"] "+main_tree.getNode(node_idx).getDuplicationRate()+" ->\t"+v);
            model.setNodeDuplicationRate(node_idx,param_x);
            return param_x;
        }   
        
        private double setNodeTransferRate(int node_idx, double param_x, boolean do_transformation)
        {
            if (do_transformation)
                param_x = Math.exp(param_x);
            //Verbose.message("ML.OAE.sNTR ["+node_idx+"] "+main_tree.getNode(node_idx).getTransferRate()+" ->\t"+v);
            model.setNodeTransferRate(node_idx,param_x);
            return param_x;
        }
        
        private double setEdgeLength(int node_idx, double param_x, boolean do_transformation)
        {
            if (do_transformation)
                param_x = Math.exp(param_x);
            //Verbose.message("ML.OAE.sEL ["+node_idx+"] "+main_tree.getNode(node_idx).getLength()+" ->\t"+v);
            model.setEdgeLength(node_idx,param_x);
            return param_x;
        }
                
        private static final double EPS=1e-7; // cubic root of machine precision for numerical derivation
        @Override
        public double[] dfunc(double[] x) 
        {
            for (int i=0; i<x.length; i++)
                if (x[i] != last_x_set[i])
                {
                    set(x, true);
                    break;
                }
            
            double y0=-getLogLikelihood();
    
            double[] retval = new double[x.length]; 
            for (int param_idx=0; param_idx<x.length; param_idx++)
            {
                int node_idx = parameter_map[param_idx] / 3;
                int parameter_type = parameter_map[param_idx] % 3;
                double param_x = x[param_idx];
                double h = Math.abs(param_x*EPS);
                {
                    double  tmp = param_x+h;
                    h = (node_idx==node_idx?tmp-param_x:h); // obfuscated code to avoid compiler optimization
                }
                if (parameter_type == PARAMETER_INDEX_DUPLICATION)
                {
                    double value_set = setNodeDuplicationRate(node_idx,param_x+h,true);
                    recomputeOnPathToRoot(node_idx);
                    double y2 = -getLogLikelihood();
                    retval[param_idx] = (y2-y0)/h;
                    //Verbose.message("ML.OAE.df ["+param_idx+"] "+param_x+"\t"+(param_x+h)+"\t"+node_idx+":"+parameter_type+"\t"+value_set+"\tlik "+y2+"\tdf "+retval[param_idx]);
                    setNodeDuplicationRate(node_idx,param_x,true);
                    recomputeOnPathToRoot(node_idx);
                } else if (parameter_type == PARAMETER_INDEX_TRANSFER)
                {
                    double value_set = setNodeTransferRate(node_idx,param_x+h,true);
                    recomputeOnPathToRoot(node_idx);
                    double y2 = -getLogLikelihood();
                    retval[param_idx] = (y2-y0)/h;
                    //Verbose.message("ML.OAE.df ["+param_idx+"] "+param_x+"\t"+(param_x+h)+"\t"+node_idx+":"+parameter_type+"\t"+value_set+"\tlik "+y2+"\tdf "+retval[param_idx]);
                    setNodeTransferRate(node_idx,param_x,true);
                    recomputeOnPathToRoot(node_idx);
                } else
                {
                    double value_set = setEdgeLength(node_idx,param_x+h,true);
                    recomputeOnPathToRoot(node_idx);
                    double y2 = -getLogLikelihood();
                    retval[param_idx] = (y2-y0)/h;
                    //Verbose.message("ML.OAE.df ["+param_idx+"] "+param_x+"\t"+(param_x+h)+"\t"+node_idx+":"+parameter_type+"\t"+value_set+"\tlik "+y2+"\tdf "+retval[param_idx]);
                    setEdgeLength(node_idx,param_x,true);
                    recomputeOnPathToRoot(node_idx);
                }
                //Verbose.message("ML.OAE.df\t"+node_idx+"\tdg "+retval[2*node_idx]+"\tdl "+retval[2*node_idx+1]+"\t("+Math.exp(gain)+","+Math.exp(loss)+")");
            }
            
            return retval;
        }
        
        @Override
        public double eval(double[] x) 
        {
            set(x,true);
            double ll = -getLogLikelihood();
            num_eval_calls++;
            Verbose.message("ML.OAE.e "+num_eval_calls+"\t"+ll);
            
            if (Double.isNaN(ll))
            {
                //for (int node_idx=0; node_idx<main_tree.getNumNodes()-1; node_idx++)
                //{
                //    Verbose.message("*****ML.OAE.e bad params node "+node_idx+"\t"+main_tree.getNode(node_idx));
                //}
                throw new LikelihoodOptimizationException("Bad parameter values in ML.OptimizeAllEdges.eval - likelihood is NaN");
            } else
            {
                if (ll < best_ll)
                {
                    if (best_x == null)
                        best_x = new double[x.length];
                    System.arraycopy(x,0,best_x,0,x.length);
                    best_ll = ll;
                }
            }
            return ll;
        }
        
        private double[] best_x;
        private double best_ll=Double.POSITIVE_INFINITY;
        
        private int num_eval_calls=0;
    }    
    
    class OptimizeRootPrior implements OneParameterFunction
    {
        OptimizeRootPrior(int param_idx)
        {
            this.param_idx =param_idx;
        }
        
        private int param_idx;
        
        public void set(double param)
        {
            model.getRootPrior().setParameter(param_idx, param);
            //setCorrectedRootPriors();
            recomputeAll();
        }
        
        @Override
        public double eval(double x)
        {
            set(x);
            double ll = getLogLikelihood();
            return -ll;
        }
    }
    
    class OptimizeForbiddenTransfer implements OneParameterFunction
    {
        public void set (double p)
        {
            model.setInvariantFractions(model.getDuplicationForbiddenProportion(), model.getLossForbiddenProportion(), p);
            recomputeAll();
        }
        
        @Override
        public double eval(double alpha)
        {
            set(alpha);
            double ll = getLogLikelihood();
            return -ll;
        }
    }
    
    class OptimizeForbiddenDuplication implements OneParameterFunction
    {
        public void set (double p)
        {
            model.setInvariantFractions(p, model.getLossForbiddenProportion(), model.getTransferForbiddenProportion());
            recomputeAll();
        }
        
        @Override
        public double eval(double alpha)
        {
            set(alpha);
            double ll = getLogLikelihood();
            return -ll;
        }
    }
    
    
    public class LikelihoodOptimizationException extends RuntimeException
    {
        private LikelihoodOptimizationException(String s)
        {
            super(s);
        }
    }
    
    private void go(String[] args) throws Exception
    {
        if (args.length < 2 || args.length>3)
        {
            System.err.println("Call as java "+this.getClass().getCanonicalName()+" [options] phylogeny table [rates]");
            System.err.println("| Options:");
            System.err.println("| \t-v");
            System.err.println("| \t-duplication_k -length_k -transfer_k -loss_k");
            System.err.println("| \t-duplication_a -length_a -transfer_a -loss_a");
            System.err.println("| \t-uniform_duplication -uniform_length -uniform_transfer");
            System.err.println("| \t-max_paralogs -min_lineages");
            System.err.println("| \t-opt_rounds -opt_delta");
            
            System.exit(2008);
        }
        reportLaunch(args);
        reportOtherArguments("Max. paralogs "+MAX_PARALOGS+", min. lineages "+MIN_PRESENT_LINEAGES);
        
        String tree_file = args[0];
        String table_file = args[1];

        TreeWithRates input_tree = new TreeWithRates(ca.umontreal.iro.evolution.Parser.readNewick(new java.io.FileReader(tree_file)));
        OccurrenceTable input_table = new OccurrenceTable(input_tree.getLeaves());
        input_table.readTable(new java.io.FileReader(table_file));
        OccurrenceTable filtered_table = input_table.filterByMaximumSize(MAX_PARALOGS,MIN_PRESENT_LINEAGES);

        RateVariation input_model = null;
        if (args.length == 3)
        {
            String rates_file = args[2];
            input_model = RateVariation.read(new java.io.FileReader(rates_file),input_tree);
        } else
        {
            input_model = new RateVariation(input_tree, new Poisson(0.1), 1,1,1,1);
        }
        // check model switches
        {
            // number of discrete gamma categories
            int ndup = input_model.getNumDuplicationRateGammaCategories();
            if (CATEGORY_DUPLICATION != -1)
            {
                ndup = CATEGORY_DUPLICATION;
                reportOtherArguments("Duplication rate Gamma categories sought: "+ndup);
            }
                
            int nloss = input_model.getNumLossRateGammaCategories();
            if (CATEGORY_LOSS != -1)
            {
                nloss = CATEGORY_LOSS;
                reportOtherArguments("Loss rate Gamma categories sought: "+nloss);
            }

            int ntra = input_model.getNumTransferRateGammaCategories();
            if (CATEGORY_TRANSFER != -1)
            {
                ntra = CATEGORY_TRANSFER;
                reportOtherArguments("Transfer rate Gamma categories sought: "+ntra);
            }
            int nlen = input_model.getNumEdgeLengthGammaCategories();
            if (CATEGORY_EDGE_LENGTH != -1)
            {
                nlen = CATEGORY_EDGE_LENGTH;
                reportOtherArguments("Edge length Gamma categories sought: "+nlen);
            }
            input_model.setNumberOfDiscreteCategories(ndup,nloss,ntra,nlen);
            boolean fdup = (input_model.getDuplicationForbiddenProportion()!=0.0 || FORBIDDEN_DUPLICATION);
            boolean ftrans = (input_model.getTransferForbiddenProportion()!=0.0 || FORBIDDEN_TRANSFER);

            // alpha parameter for gamma distribution
            if (!Double.isNaN(ALPHA_DUPLICATION))
                input_model.setDuplicationRateAlpha(ALPHA_DUPLICATION);
            if (!Double.isNaN(ALPHA_LOSS))
                input_model.setLossRateAlpha(ALPHA_LOSS);
            if (!Double.isNaN(ALPHA_TRANSFER))
                input_model.setTransferRateAlpha(ALPHA_TRANSFER);
            if (!Double.isNaN(ALPHA_EDGE_LENGTH))
                input_model.setEdgeLengthAlpha(ALPHA_EDGE_LENGTH);
            
            reportOtherArguments("Uniform rates: "
                    +"duplication - "+(UNIFORM_DUPLICATION?"yes":"no")
                    +", transfer - "+(UNIFORM_TRANSFER?"yes":"no")
                    +", length - "+(UNIFORM_EDGE_LENGTH?"yes":"no"));
            reportOtherArguments("0-rate transfer category: "+(ftrans?"yes":"no")
                    +", 0-rate duplication category: "+(fdup?"yes":"no"));
            // check if uniform rates
            setUniformEdgeParameters(UNIFORM_DUPLICATION, UNIFORM_TRANSFER, UNIFORM_EDGE_LENGTH);
            setForbiddenCategories(ftrans,fdup);
        }

        
        // optimization parameters
        this.optimize_maxsteps = OPTIMIZATION_ROUNDS;
        this.optimize_ll_delta = OPTIMIZATION_LL_DELTA;
        reportOtherArguments("Maximum number of rounds in the optimization: "+OPTIMIZATION_ROUNDS);
        reportOtherArguments("Optimization convergence delta: "+OPTIMIZATION_LL_DELTA);
        
        computation_model.init(filtered_table,input_model);
        init();
        
        double ll = optimize();
        
        System.out.println("# final likelihood "+ll+"\t(p0="+computation_model.getAbsentProfileProbability(MIN_PRESENT_LINEAGES)+")");
        System.out.println("# "+input_tree.getRoot().newickTree());
        System.out.println(input_model.tableRates());
    }
    
    public static void testComputation(LikelihoodComputation LC, String[] args) throws Exception
    {
        if (args.length!=3)
        {
            System.err.println("Call as $0 phylogeny table rates");
            System.exit(2008);
        }
        String tree_file = args[0];
        String table_file = args[1];
        String rates_file = args[2];

        TreeWithRates input_tree = new TreeWithRates(ca.umontreal.iro.evolution.Parser.readNewick(new java.io.FileReader(tree_file)));
        OccurrenceTable input_table = new OccurrenceTable(input_tree.getLeaves());
        input_table.readTable(new java.io.FileReader(table_file));
        OccurrenceTable filtered_table = input_table.filterByMaximumSize(MAX_PARALOGS,MIN_PRESENT_LINEAGES);
        System.out.println("# Filtered table: "+filtered_table.getNumFamilies()+" families.");
//        System.out.println("unconstrained log-likelihood "+filtered_table.getUnconstrainedLogLikelihood());

        RateVariation input_model = RateVariation.read(new java.io.FileReader(rates_file),input_tree);        
        LC.init(filtered_table, input_model);
        
        ML O = new ML(LC);
        O.recomputeAll();
        
        double L = O.getLogLikelihood();
        System.out.println("log-likelihood "+L+"\tp0 "+LC.getAbsentProfileProbability(MIN_PRESENT_LINEAGES));
    }
    
    public static void main(String[] args) 
    {
        Verbose.setVerbose(false);
                
        ML O = new ML(new StableComputation()); //InclusionExclusionComputation());
        
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
                } else if (arg_switch.equals("duplication_k"))
                {
                    CATEGORY_DUPLICATION = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("transfer_k") || arg_switch.equals("gain_k"))
                {
                    CATEGORY_TRANSFER = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("length_k"))
                {
                    CATEGORY_EDGE_LENGTH = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("loss_k"))
                {
                    CATEGORY_LOSS = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("duplication_a"))
                {
                    ALPHA_DUPLICATION = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("transfer_a") || arg_switch.equals("gain_a"))
                {
                    ALPHA_TRANSFER = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("length_a"))
                {
                    ALPHA_EDGE_LENGTH = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("loss_a"))
                {
                    ALPHA_LOSS = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("uniform_duplication"))
                {
                    UNIFORM_DUPLICATION = "true".equals(arg_value);
                } else if (arg_switch.equals("uniform_transfer") || arg_switch.equals("uniform_gain"))
                {
                    UNIFORM_TRANSFER = "true".equals(arg_value);
                } else if (arg_switch.equals("uniform_length"))
                {
                    UNIFORM_EDGE_LENGTH = "true".equals(arg_value);
                } else if (arg_switch.equals("max_paralogs"))
                {
                    MAX_PARALOGS = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("min_lineages"))
                {
                    MIN_PRESENT_LINEAGES = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("opt_rounds"))
                {
                    OPTIMIZATION_ROUNDS = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("opt_delta"))
                {
                    OPTIMIZATION_LL_DELTA = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("forbidden_transfer") || arg_switch.equals("forbidden_gain"))
                {
                    FORBIDDEN_TRANSFER = "true".equals(arg_value);
                } else if (arg_switch.equals("forbidden_duplication"))
                {
                    FORBIDDEN_DUPLICATION = "true".equals(arg_value);
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
