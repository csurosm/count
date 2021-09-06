/*
 * SurvivalProbabilityFactory.java
 *
 * Created on April 21, 2008, 11:42 AM
 *
 */

package ca.umontreal.iro.evolution.genecontent;

/**  
 *
 * Facilitates the aggregate computations over a set of profiles
 *
 * @author csuros
 */

import ca.umontreal.iro.banality.Verbose;

import ca.umontreal.iro.matek.DiscreteDistribution;
import ca.umontreal.iro.matek.ShiftedGeometric;
import ca.umontreal.iro.matek.Poisson;
import ca.umontreal.iro.matek.NegativeBinomial;
import ca.umontreal.iro.matek.PointDistribution;

public class SurvivalProbabilityFactory 
{
    
    /**
     * Initializes the factory with the given tree and profile set
     *
     * @param tree the common tree structure TreeWithRates used in the calculations
     * @param profile_array entry [i][j] is the size of the j-th family in the i-th profile 
     */
    public SurvivalProbabilityFactory(TreeWithRates tree, OccurrenceTable T) 
    {
        this.tree = tree;
        this.profiles = T.getProfiles();
        computeExtremeSizes();
    }
    
    public SurvivalProbabilityFactory(TreeWithRates tree, PhyleticProfile[] profiles)
    {
        this.tree = tree;
        this.profiles = profiles;
        computeExtremeSizes();
    }
    
    private TreeWithRates tree;
    
    private PhyleticProfile[] profiles;

    /**
     * Returns the main tree used at instantiation
     */
    public TreeWithRates getTree() 
    {
        return tree;
    }
    
    /**
     * Returns the array of profiles for this factory
     */
    public PhyleticProfile[] getProfiles()
    {
        return profiles;
    }
        
    private int[] max_survivals[];
    
    /**
     * Computes extreme sizes for survivals. These are the last n,m entries for survival probabilities
     */
    private void computeExtremeSizes()
    {
        //computeSharedProfiles();
        
        NodeWithRates[] nodes = tree.getDFT();
        max_survivals = new int[nodes.length-1][]; // last entry is the root at index nodes.length-1
        // compute maximum sizes at each non-leaf node
        int[] max_sizes = new int[nodes.length];
        int[][] subtree_sizes = new int[nodes.length][profiles.length];
        Verbose.message("SPF.cES nodes "+nodes.length+" profiles "+profiles.length);
        
        for (int profile_idx = 0; profile_idx<profiles.length; profile_idx++)
        {
            int[] s = profiles[profile_idx].computeSubtreeSizes(tree);
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
                subtree_sizes[node_idx][profile_idx]=s[node_idx];
        }
        
        for (int node_idx = 0; node_idx<nodes.length; node_idx++)
        {
            NodeWithRates N = nodes[node_idx];
            if (!N.isLeaf())
            {
                int[] ss = subtree_sizes[node_idx];
                int n_max = ss[0];
                for (int profile_idx = 1; profile_idx<profiles.length; profile_idx++)
                {
                    int node_size = ss[profile_idx];
                    if (n_max<node_size)
                        n_max = node_size;
                }
                max_sizes[node_idx] = n_max;
                //Verbose.message("#**PPF.cES "+node_idx+"/"+N.getTaxonName()+"\t"+n_max+"\t// "+N);
            } 
        }
        
        // compute the maximum subtree sizes at each node for each parental size
        for (int node_idx = 0; node_idx<nodes.length-1; node_idx++)
        {
            int parent_idx = tree.getParentIndex(node_idx);
            int n_max = max_sizes[parent_idx];
            max_survivals[node_idx] = new int[n_max+1];
            int[] node_ss = subtree_sizes[node_idx];
            int[] parent_ss = subtree_sizes[parent_idx];
            
            for (int profile_idx = 0; profile_idx<profiles.length; profile_idx++)
            {
                int parent_size = parent_ss[profile_idx];
                int node_size = node_ss[profile_idx];
                if (max_survivals[node_idx][parent_size]<node_size)
                    max_survivals[node_idx][parent_size]=node_size;
            }
        }
        
        // make sure the max_survivals array is monotone 
        for (int node_idx = 0; node_idx<nodes.length-1; node_idx++)
        {
            int n_max = max_survivals[node_idx].length-1;
            for (int n=n_max-1; n>=0; n--)
            {
                int m_previous = max_survivals[node_idx][n+1];
                int m_current = max_survivals[node_idx][n];
                if (m_current<m_previous)
                    max_survivals[node_idx][n] = m_previous;
            }
            // experimental ...
            for (int n=0; n<n_max; n++)
            {
                max_survivals[node_idx][n] += 10;
            }
            
        }
    }
    
    /**
     * Returns the array of maximum m values for m|n on the edge.
     */
    public int[] getExtremeSizes(int node_idx)
    {
        return max_survivals[node_idx];
    }
    
    /**
     *
     * @param rate_tree tree with the same structure (but possibly different branch parameters) as the tree used at initialization  
     *
     */ 
    public Rates getSurvivalProbability(TreeWithRates rate_tree)
    {
        return new Rates(rate_tree);
    }
    
    public class Rates
    {
        private Rates(TreeWithRates rate_tree)
        {
            this.rate_tree = rate_tree;
            initDataStructures();
            recomputeAll();
        }
        
        private TreeWithRates rate_tree;
        
        public void recomputeAtNode(int node_idx)
        {
            NodeWithRates N = rate_tree.getNode(node_idx);

            // compute extinction probability for this guy
            if (N.isLeaf())
            {
                extinction_probability[node_idx] = 0.0;
            } else
            {
                int num_children = N.getNumChildren();
                extinction_probability[node_idx] = 1.0;
                for (int ci=0;ci<num_children;ci++)
                {
                    int child_idx = tree.getChildIndex(node_idx,ci);
                    extinction_probability[node_idx] *= dup0[child_idx];
                }
            }
            if (!N.isRoot())
            {
                double p,q,pq1,p1;
                {
                    // compute duplication distribution at this guy
                    DiscreteDistribution G = N.getDuplicationDistribution(extinction_probability[node_idx]);
                    dup0[node_idx] = G.getDistribution(0)[0];

                    // G might be a point distribution if duplication rate is 0.
                    if (G instanceof ShiftedGeometric)
                    {
                        double[] params = G.getParameters();
                        p = params[0];
                        q = params[1];
                        pq1 = 1.0-p-q;
                    } else 
                    {
                        // duplication rate is 0...
                        p = G.getParameters()[0];
                        q = 0.0;
                        pq1 = 1.-p;
                    }
                    p1 = 1.0-p;
                }

                double[][] sp = survival_probability[node_idx];
                // survival probabilities for n==0
                {
                    DiscreteDistribution H = N.getTransferDistribution(extinction_probability[node_idx]);

                    int mmax0 = survival_probability[node_idx][0].length-1;
                    double Hm[] = H.getDistribution(mmax0);
                    for (int m=0; m<mmax0; m++)
                    {
                        sp[0][m] = Hm[m];
                        //if (Double.isNaN(Hm[m]))
                        //{
                        //    double[] params = H.getParameters();
                        //    Verbose.message("SPF.R.rAN H["+m+"] is NaN\t"+params[0]+"\t"+params[params.length-1]);
                        //}
                    }
                }

                int nmax = sp.length-1;
                for (int n=1; n<=nmax; n++)
                {
                    // m==0
                    int mmax = sp[n].length-1;
                    //Verbose.message("SPF.R.re "+node_idx+"\tn "+n+" m 0");
                    sp[n][0] = sp[n-1][0]*p;
                    if (mmax !=0)
                    {
                        sp[n][1] = sp[n-1][1]*p + sp[n-1][0]*p1;
                        for (int m=2; m<=mmax; m++)
                        {
                            sp[n][m] = pq1*sp[n-1][m-1]+p*sp[n-1][m];
                            if (q!=0.)
                                sp[n][m] += q*sp[n][m-1];
                        }
                    }
                }
                //for (int n=0; n<survival_probability[node_idx].length; n++)
                //    for (int m=0; m<survival_probability[node_idx][n].length; m++)
                //        if (Double.isNaN(survival_probability[node_idx][n][m]))
                //        {
                //            Verbose.message("SPF.R.rAN sp "+node_idx+"/"+rate_tree.getNode(node_idx)+"\t"+m+"|"+n+" NaN");
                //            Verbose.message("SPF.R.rAN "+node_idx+"/"+rate_tree.getNode(node_idx)+"\textinction "+extinction_probability[node_idx]);
                //            System.exit(666);
                //        }

            } // for non-root node
        } 
        
        public void recomputeOnPathToRoot(int node_idx)
        {
            recomputeAtNode(node_idx);
            if (!tree.getNode(node_idx).isRoot())
                recomputeOnPathToRoot(tree.getParentIndex(node_idx));
        }

        /**
         * Recomputes the survival and extinction probabilities (must be called when rate tree changes)
         */
        public void recomputeAll()
        {
            for (int node_idx=0; node_idx<tree.getNumNodes(); node_idx++)
                recomputeAtNode(node_idx);
        }
        
        /**
         * Returns the extinction probability for a subtree
         * 
         * @param node_idx index of the subtree root
         */
        public double getExtinction(int node_idx)
        {
            return extinction_probability[node_idx];
        }
        
        /**
         * Probability that a family is represented in at most one lineage
         *
         */
        public double getAbsentOrSingle(DiscreteDistribution root_prior)
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
                DiscreteDistribution Hy = N.getTransferDistribution(extinction_probability[edge_idx]);
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
            
            
            double psum = 0.0;
            for (int leaf_idx=0; leaf_idx<num_leaves; leaf_idx++)
            {
                psum += present_at[leaf_idx];
            }
            double p_absent = getAllAbsentProbability()*survivingRootPrior(root_prior).getDistribution(0)[0];
            double pcorr = psum- (num_leaves-1.0)*p_absent;
            double p1 = psum - num_leaves*p_absent;
            //Verbose.message("SPF.R.gAOS sum "+psum+"\tabs "+p_absent+"\tdiff "+pcorr+"\tp1 "+p1);
            
            return pcorr;
        }
        
        /**
         * Computes the probability of all-0: it's simply the product of H_y(0) over all edges xy
         */
        public double getAllAbsentProbability()
        {
            double p_absent = 1.0;
            int num_nodes = rate_tree.getNumNodes();
            for (int node_idx=0; node_idx<num_nodes; node_idx++)
            {
                NodeWithRates N = rate_tree.getNode(node_idx);
                if (!N.isRoot()) // i.e., node_idx<num_nodes-1
                    p_absent *= survival_probability[node_idx][0][0];
            }
           return p_absent;
        }

        
        /**
         * Survival probabilities p(m|n) on an edge
         *
         * @
         */
        public double[] getSurvival(int node_idx, int n)
        {
            return survival_probability[node_idx][n];
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

        private double[] extinction_probability;
        private double[][][] survival_probability;
        private double[] dup0;
        
        private void initDataStructures()
        {
            NodeWithRates[] nodes = tree.getDFT();
            extinction_probability = new double[nodes.length]; // none for the root
            survival_probability = new double[nodes.length-1][][]; //none for the root
            //Verbose.message("SPF.R.iDS sp[0:"+survival_probability.length+"][][]");
            for (int node_idx=0; node_idx<nodes.length-1; node_idx++)
            {
                int n_max = max_survivals[node_idx].length-1;
                survival_probability[node_idx] =  new double[n_max+1][];
                //Verbose.message("SPF.R.iDS sp["+node_idx+"][0:"+survival_probability[node_idx].length+"][]");
                for (int n=0; n<=n_max; n++)
                {
                    int m_max = max_survivals[node_idx][n];
                    survival_probability[node_idx][n] = new double[m_max+1];
                    //Verbose.message("SPF.R.iDS sp["+node_idx+"]["+n+"][0:"+survival_probability[node_idx][n].length+"]");
                }
            }
            dup0 = new double[nodes.length-1];
        }
        
        public TreeWithRates getTree()
        {
            return rate_tree;
        }
    }
}
