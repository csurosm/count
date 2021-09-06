package ca.umontreal.iro.evolution.genecontent;

/**
 *
 * @author csuros
 */

import ca.umontreal.iro.banality.Verbose;

import ca.umontreal.iro.matek.LowerTriangularMatrix;
import ca.umontreal.iro.matek.DiscreteDistribution;
import ca.umontreal.iro.matek.NegativeBinomial;
import ca.umontreal.iro.matek.PointDistribution;
import ca.umontreal.iro.matek.Poisson;
import ca.umontreal.iro.matek.ShiftedGeometric;

public class InclusionExclusionComputation implements LikelihoodComputation
{
    public InclusionExclusionComputation(OccurrenceTable table, RateVariation model)
    {
        init(table, model);
    }
    
    public InclusionExclusionComputation()
    {}
    
    private OccurrenceTable table;
    private RateVariation model;
    
    public OccurrenceTable getTable(){return table;}
    public RateVariation getModel(){return model;}

    public void init(OccurrenceTable table, RateVariation model)
    {
        this.table = table;
        this.model = model;
        this.profiles = table.getProfiles();
        this.main_tree = model.getMainTree();
        computeExtremeSizes();
        initClassParameters();
        initProfileParameters();
    }
    
    private PhyleticProfile[] profiles;
    
    private TreeWithRates main_tree;
    
    private int[] max_survivals[];
    
    private Rates[] class_parameters;
    
    private void initClassParameters()
    {
        int num_classes = model.getNumClasses();
        class_parameters = new Rates[num_classes];
        for (int cidx=0;cidx<num_classes; cidx++)
            class_parameters[cidx]=new Rates(model.getRateTree(cidx));
    }
    
    /**
     * conditional likelihoods for each rate category and profile.
     * Indexing: profile_parameters[profile index][category index]
     */
    private LowerConditionals[][] profile_parameters;
    
    private void initProfileParameters()
    {
        int num_classes = model.getNumClasses();
        int num_profiles = profiles.length;
        profile_parameters = new LowerConditionals[num_profiles][num_classes];
        for (int pidx=0; pidx<profiles.length; pidx++)
            for (int cidx=0; cidx<num_classes; cidx++)
                profile_parameters[pidx][cidx] = new LowerConditionals(profiles[pidx],class_parameters[cidx]);
    }
    
    public void recomputeCategorySupport(int cidx)
    {
        class_parameters[cidx].recomputeAll();
    }
    
    public void recomputeCategorySupport(int cidx, int node_idx)
    {
        class_parameters[cidx].recomputeAtNode(node_idx);
    }

    public DiscreteDistribution getRootPrior(int cidx)
    {
        return class_parameters[cidx].survivingRootPrior(model.getRootPrior());
    }
    
    public double[] getRootLikelihoods(int pidx, int cidx)
    {
        int root_idx  = main_tree.getNumNodes()-1; // last node in postorder traversal
        return getNodeLikelihoods(root_idx, pidx, cidx);
    }

    public double[] getNodeLikelihoods(int node_idx, int pidx, int cidx)
    {
        LowerConditionals LC = profile_parameters[pidx][cidx];
        double[] L = LC.getLikelihoods(node_idx);
        return L;
    }

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
                if (model.getClassProbability(cidx)!=0.)
                {
                    double p_class = getRootPrior(cidx).getDistribution(0)[0];
                    double z = class_parameters[cidx].getAllAbsentProbability(); //LC0.getLikelihood(root_idx,0);
                    p0 += model.getClassProbability(cidx)*p_class*z;
                }
        
        } else if (min_lineages == 2)
        {
            for (int cidx=0; cidx<num_classes; cidx++)
                if (model.getClassProbability(cidx)!=0.)
                {
                    double z = class_parameters[cidx].getAbsentOrSingle(main_prior);
                    p0 += model.getClassProbability(cidx)*z;
                }
        }
        //Verbose.message("ML.gAEP p(absent) "+p0);
        
        return p0;
    }
    
    
    public void recomputeProfileSupport(int pidx, int cidx, int node_idx)
    {
        if (!main_tree.getNode(node_idx).isRoot())
        {
            profile_parameters[pidx][cidx].recomputeOnPathToRoot(main_tree.getParentIndex(node_idx));
        }
    }
    
    public void recomputeProfileSupport(int pidx,int cidx)
    {
        profile_parameters[pidx][cidx].recomputeAll();
    }
    
    
    /**
     * Computes extreme sizes for survivals. These are the last n,m entries for survival probabilities
     */
    private void computeExtremeSizes()
    {
        //computeSharedProfiles();
        
        NodeWithRates[] nodes = main_tree.getDFT();
        max_survivals = new int[nodes.length-1][]; // last entry is the root at index nodes.length-1
        // compute maximum sizes at each non-leaf node
        int[] max_sizes = new int[nodes.length];
        int[][] subtree_sizes = new int[nodes.length][profiles.length];
        Verbose.message("SPF.cES nodes "+nodes.length+" profiles "+profiles.length);
        
        for (int profile_idx = 0; profile_idx<profiles.length; profile_idx++)
        {
            int[] s = profiles[profile_idx].computeSubtreeSizes(main_tree);
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
            int parent_idx = main_tree.getParentIndex(node_idx);
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
     * Machine precision for double
     */
    private static final double EPS=1e-14;
    
    /**
     * Computes the binomial inversion for a series of values. 
     * Given the array a[0..], the formula b[n] = sum_{i=0}^n (n choose i) p^{n-1} (1-p)^i a[i] 
     * defines the entries of the array b[0..]. 
     * Returns the array a.
     *
     * @param b an array of at least nmax+1 elements
     * @param p probability parameter in the inversion formula
     * @param nmax maximum index needed
     * @return array of length nmax+1
     */
    public static double[] binomialInversion(double[] b, double p, int nmax)
    {
        double[] a = new double[nmax+1];
        if (p == 0.0)
        {
            System.arraycopy(b,0,a,0,nmax+1);
        } else if (p==1.0)
        {
            java.util.Arrays.fill(a,b[0]);
        } else
        {
            a[0] = b[0];    
            double d[] = new double[nmax+1]; d[0]=1.0;
            double pn[] = new double[nmax+1]; pn[0]=1.0;
            double q = (1.-p) / p;
            for (int n=1; n<=nmax; n++)
            {
                if (n%2 == 1)
                {
                    d[n] = d[n-1]*(1-p);
                    pn[n] = pn[n-1]*p;
                }
                else
                {
                    d[n] = d[n/2]*d[n/2];
                    pn[n] = pn[n/2]*pn[n/2];
                }
                //Verbose.message("PP.bI power "+n+"\td "+d[n]+"\t"+Math.pow(1.-p,n)+"\tpn "+pn[n]+"\t"+Math.pow(p,n));
                
                
                if (true)
                {
                    double f = pn[n];
                    double subtract = 0.0;
                    for (int i=0; i<n; i++)
                    {
                        subtract += a[i] * f;
                        if (i!=n-1)
                            f *= q*(n-i)/(i+1.0);
                        //Verbose.message("PP.bI\t"+n+"\ti "+i+"\t"+subtract+"\t"+f);
                    }

                    if (b[n] != 0.0) 
                    { // check for roundoff problems
                        double h = 1.0-subtract/b[n];
                        if (h<EPS)
                        {
                            // they are too close
                            a[n] = 0.0;
                            //Verbose.message("IEC.bI\t"+n+"\tb "+b[n]+"\tsub "+subtract+"\t= "+a[n]+"\tdiv "+d[n]+"\t// ++++ CLOSE\th "+h);
                        } else 
                        {
                            a[n] = (b[n]-subtract)/d[n];
                            //Verbose.message("IEC.bI\t"+n+"\tb "+b[n]+"\tsub "+subtract+"\t= "+a[n]+"\tdiv "+d[n]+"\th "+h);
                        }
                    }
                } 
            } // for n
            double[][] factors = new double[nmax+1][nmax+1];
            for (int n=0; n<=nmax; n++)
            {
                double x = Math.pow(p,n);
                for (int i=0; i<=n; i++)
                {
                    factors[n][i] = x;
                    if (i!=n)
                    {
                        x *= q*(n-i)/(i+1.0);
                    }
                }
            }
            double [] bb = new double[nmax+1];
            System.arraycopy(b,0,bb,0,nmax+1);
            for (int n=0; n<=nmax; n++)
            {
                for (int i=0; i<n; i++)
                    factors[n][i] /= factors[n][n];
                bb[n] /= factors[n][n];
                //Verbose.message("PP.bI "+n+"\tnormalize "+factors[n][n]+"\t"+bb[n]);
                factors[n][n]=1.0;
            }
            LowerTriangularMatrix M = new LowerTriangularMatrix(factors);
            double[] x = M.solve(bb);
            double delta = Double.POSITIVE_INFINITY;
            int rr = 0;
            for (rr=0; rr<10 && delta>1e-6; rr++)
            {
                delta = M.improveSolution(x,bb);
                //Verbose.message("IEC.bI iteration "+rr+"\t"+delta);
            }
            //for (int i=0; i<x.length; i++)
            //{
            //    Verbose.message("IEC.bI "+i+"\tx "+x[i]+"\ta "+a[i]+"\tb "+bb[i]+"\t// "+b[i]+"\t// delta "+delta+"\t"+rr);
            //}
            a = x;
        }
        //System.exit(99);
        return a;
    }
    
    
    
    /**
     * Inner class for storing probabilities for population change: one 
     * instance for each rate category.
     */
    private class Rates
    {
        private Rates(TreeWithRates rate_tree)
        {
            this.rate_tree = rate_tree;
            initDataStructures();
            //recomputeAll();
        }
        
        private TreeWithRates rate_tree;
        
        public TreeWithRates getTree()
        {
            return rate_tree;
        }
        
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
                    int child_idx = main_tree.getChildIndex(node_idx,ci);
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
                        sp[n][1] = sp[n-1][1]*p + sp[n-1][0]*p1*(1.-q);
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
        
        private double[] extinction_probability;
        private double[][][] survival_probability;
        private double[] dup0;
        
        private void initDataStructures()
        {
            NodeWithRates[] nodes = main_tree.getDFT();
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
        
        public void recomputeOnPathToRoot(int node_idx)
        {
            recomputeAtNode(node_idx);
            if (!main_tree.getNode(node_idx).isRoot())
                recomputeOnPathToRoot(main_tree.getParentIndex(node_idx));
        }

        /**
         * Recomputes the survival and extinction probabilities (must be called when rate tree changes)
         */
        public void recomputeAll()
        {
            for (int node_idx=0; node_idx<main_tree.getNumNodes(); node_idx++)
                recomputeAtNode(node_idx);
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

        /**
         * Returns the extinction probability for a subtree
         * 
         * @param node_idx index of the subtree root
         * @return probability that an individual at this node has no descendants at the leaves within the node's subtree
         */
        public double getExtinction(int node_idx)
        {
            return extinction_probability[node_idx];
        }
        
        /**
         * Survival probabilities p(m|n) on an edge
         *
         * @param node_idx lower node on the edge
         * @param n population size at parent
         * @return distribution of surviving population size at the child
         */
        public double[] getSurvival(int node_idx, int n)
        {
            return survival_probability[node_idx][n];
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
        
    }    
    
    /**
     * Inner class for computing likelihoods for a profile
     */
    private class LowerConditionals
    {
        private LowerConditionals(PhyleticProfile profile, Rates survival)
        {
            this.survival = survival;
            this.profile = profile;
            rate_tree = survival.getTree();
            initDataStructures();
            //recomputeAll();
        }
        
        private PhyleticProfile profile;
        private Rates survival;
        private TreeWithRates rate_tree;
        private int[] subtree_sizes;
        
        private void initDataStructures()
        {
            NodeWithRates[] nodes = rate_tree.getDFT();
            likelihood = new double[nodes.length][];
            uncorrected_likelihood = new double[nodes.length][];
            subtree_sizes = profile.computeSubtreeSizes(rate_tree);
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                likelihood[node_idx] = new double[subtree_sizes[node_idx]+1];
                uncorrected_likelihood[node_idx] = new double[subtree_sizes[node_idx]+1];
            }
        }

        private double[][] likelihood;
        private double[][] uncorrected_likelihood;
 
        public double getLikelihood(int node_idx, int n)
        {
            return likelihood[node_idx][n];
        }
        
        public double[] getLikelihoods(int node_idx)
        {
            return likelihood[node_idx];
        }
        
        public double[] getUncorrectedLikelihoods(int node_idx)
        {
            return uncorrected_likelihood[node_idx];
        }
        
        
        public void recomputeAtNode(int node_idx)
        {
            NodeWithRates N = rate_tree.getNode(node_idx);
            double[] piciL = uncorrected_likelihood[node_idx];
            java.util.Arrays.fill(piciL,0.0);
            
            if (N.isLeaf())
            {
                java.util.Arrays.fill(likelihood[node_idx],0.0);
                int v = profile.get(node_idx);
                piciL[v]=likelihood[node_idx][v] = 1.0;
                if (Verbose.isVerbose())
                {
                    Verbose.message("IEC.LC.rAN leaf "+node_idx+"/"+N.newickName()+" ================\t"+N);
                    Verbose.message("SC.LC.rAN "+profile.getPatternString()+"\t"+node_idx+"\tL["+v+"]="+likelihood[node_idx][v]);
                }
            } else
            {
                double extinct = survival.getExtinction(node_idx);
                int num_children = N.getNumChildren();

                Verbose.message("IEC.LC.rAN "+profile.getPatternString()+"\t"+node_idx+"/"+rate_tree.getNode(node_idx)+"\textinct "+extinct);
                for (int n=0; n<piciL.length; n++)
                {
                    piciL[n] = 1.0;
                    for (int ci=0; ci<num_children; ci++)
                    {
                        double Ledge = 0.0;
                        int child_idx = rate_tree.getChildIndex(node_idx,ci);
                        double[] Lc = likelihood[child_idx]; 
                        double[] pmn = survival.getSurvival(child_idx,n);
                        int Mmax = subtree_sizes[child_idx];
                        for (int m=0; m<=Mmax; m++)
                        {
                            double z = Lc[m] * pmn[m];
                            //Verbose.message("IEC.LC.rAN "+node_idx+"\tn "+n+"\tci "+ci+"\tm "+m+"\tLc "+Lc[m]+"\tpmn "+pmn[m]+"\t// "+rate_tree.getNode(child_idx));
                            Ledge += z;
                        }
                        //Verbose.message("IEC.LC.rAN "+node_idx+"\tn "+n+"\tci "+ci+"\tLedge "+Ledge);
                        piciL[n] *= Ledge;
                    }
                    //Verbose.message("IEC.LC.rAN "+node_idx+"\tn "+n+"\tpiciL "+piciL[n]);
                    
                }
                likelihood[node_idx] = binomialInversion(piciL, extinct, piciL.length-1);  
                for (int n=0; n<piciL.length; n++)
                {
                    if (Double.isInfinite(likelihood[node_idx][n]) || Double.isNaN(likelihood[node_idx][n]))
                    {
                        for (int j=0; j<piciL.length; j++)
                        {
                            System.out.println("#### NUMERICAL_ERROR IEC.LC.rAN "+node_idx+"/"+rate_tree.getNode(node_idx).getTaxonName()+"\t["+j+"]\tL "+likelihood[node_idx][j]+"\tl "+piciL[j]);
                        }
                        System.exit(2008);
                    }
                }

                if (Verbose.isVerbose())
                {
                    StringBuffer distr = new StringBuffer();
                    double[] nG = N.getDuplicationDistribution(extinct).getDistribution(3);
                    double[] nH = N.getTransferDistribution(extinct).getDistribution(3);
                    distr.append("G[]={");
                    for (int i=0; i<nG.length; i++)
                    {
                        if (i!=0) distr.append(",");
                        distr.append(nG[i]);
                    }
                    distr.append("}  H[]={");
                    for (int i=0; i<nH.length; i++)
                    {
                        if (i!=0) distr.append(", ");
                        distr.append(nH[i]);
                    }
                    
                    Verbose.message("IEC.LC.rAN "+node_idx+"/"+N.newickName()+" extinct "+extinct+" ================\t"+distr);
                    for (int m=0; m<likelihood[node_idx].length; m++)
                    {
                        Verbose.message("IEC.LC.rAN "+node_idx+"\tL["+m+"]="+likelihood[node_idx][m]+"\tpici "+piciL[m]);
                    }
                }
                
                //    L[n] = piciL[n];
                //    if (n!=0)
                //    {
                //        double corr = 0.0;
                //        double factor = Math.pow(extinct,n);
                //        double DD1 = (1.-extinct)/extinct;
                //        for (int i=0; i<n; i++)
                //        {
                //            corr += L[i]*factor;
                //            if (i!=n-1)
                //            {
                //                factor *= DD1*(n-i)/(i+1.0);
                //            }
                //            //Verbose.message("PP.LLC.rAN "+node_idx+"/"+N.getTaxonName()+"\t"+n+"\ti "+i+"\t"+corr+"\t"+factor);
                //        }
                //        // check for numerical error
                //        if (L[n]!=0.)
                //        {
                //            double h = 1.0-corr/L[n];
                //            if (h<EPS)
                //            {
                //                // they are too close
                //                //Verbose.message("PP.LC.rec n "+n+" piciL "+piciL+"\tcorr "+corr);
                //                L[n] = 0.0;
                //           } else 
                //            {
                //                double Dn1 = Math.pow(1.0-extinct,n);
                //                L[n] = (L[n]-corr)/Dn1;
                //                //Verbose.message("PP.LC.rAN "+node_idx+"/"+N.getTaxonName()+"\t"+n+"\tpL "+piciL[n]+"\tcorr "+corr+"\t= "+L[n]+"\t"+Dn1);
                //            }
                //        }
                //    } // n != 0
                //} // for n
            } // N is not a leaf
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
                recomputeAtNode(node_idx);
        }        
        
        
    }
    
    
    public static void main(String[] args) throws Exception
    {
        Verbose.setVerbose(true);
        ML.testComputation(new InclusionExclusionComputation(), args);
    }
}
