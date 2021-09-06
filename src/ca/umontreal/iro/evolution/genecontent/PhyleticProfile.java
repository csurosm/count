/*
 * PhyleticProfile.java
 *
 * Created on April 19, 2008, 12:00 AM
 *
 */

package ca.umontreal.iro.evolution.genecontent;

/**
 * Computations over profiles (size distribution across tree leaves)
 *
 * @author csuros
 */

import java.util.Arrays;
import java.util.Vector;

import ca.umontreal.iro.matek.DiscreteDistribution;
import ca.umontreal.iro.matek.LowerTriangularMatrix;

import ca.umontreal.iro.banality.Verbose;

public class PhyleticProfile 
{

    /**
     * Instantiaton method for extending classes. setProfile() will have to be called 
     * separately.
     */
    protected PhyleticProfile(){}

    /**
     * Sets the profile here.
     * 
     * @param profile the profile: size ditribution at the leaves
     */
    protected void setProfile(int[] profile){this.profile=profile;}
    
    /**
     * Instantiates the class
     *
     * @param profile array of sizes across the terminal taxa
     */
    public PhyleticProfile(int[] profile) 
    {
        setProfile(profile);
    }
    
    private int[] profile;

    /**
     * Profile with which this instance was created
     * 
     * @return the array with which this instance was created
     */
    public int[] getProfile()
    {
        return profile;
    }
    
    /**
     * Profile entry (family size) for a particular node.
     * 
     * @param leaf_idx index of a leaf (0..length of profile-1)
     * @return profile value at leaf_idx
     */
    public int get(int leaf_idx)
    {
        return profile[leaf_idx];
    }
    
    /**
     * Reduces the extended profile into a binary profile.
     * 
     * @return only 0-1 [absence/presence] profile: larger family sizes are replaced by 1
     */
    public PhyleticProfile binaryProfile()
    {
        int[] ap = new int[profile.length];
        for (int i=0; i<profile.length; i++)
            ap[i] = (profile[i]==0?0:1);
        return new PhyleticProfile(ap);
    }

    /**
     * Number of zero entries in the profile
     * @return how many times there is a 0 in the profile
     */
    public int getAbsenceCount()
    {
        int num_zeros=0;
        for (int i=0; i<profile.length; i++)
            if (profile[i]==0)
                num_zeros++;
        return num_zeros;
    }

    /**
     * Sum of the family sizes in the profile
     * @return sum of sizes
     */
    public int sumSizes()
    {
        int sum = 0;
        for (int i=0; i<profile.length; i++)
            sum += Math.max(0,profile[i]);
        return sum;
    }
    
    /**
     * Computes the sum of sizes seen in each subtree
     *
     * @param tree phylogeny for the terminal taxa: rates are ignored, but the fixed postorder traversal is used 
     * 
     * @return subtree sizes in the same order as the nodes in the tree traversal
     */
    public int[] computeSubtreeSizes(TreeWithRates tree)
    {
        NodeWithRates[] nodes = tree.getDFT();
        int[] subtree_sizes = new int[nodes.length];
        for (int node_idx=0; node_idx<nodes.length; node_idx++)
        {
            NodeWithRates N = nodes[node_idx];
            if (N.isLeaf())
            {
                subtree_sizes[node_idx] = profile[node_idx];
            } else
            {
                subtree_sizes[node_idx] = 0;
                for (int cidx=0; cidx<N.getNumChildren(); cidx++)
                {
                    NodeWithRates C = (NodeWithRates)N.getChild(cidx);
                    int Ci = C.getId();
                    subtree_sizes[node_idx]+=subtree_sizes[Ci];
                }
            }
            // experimental ...
            // subtree_sizes[node_idx]+=10;
        }
        return subtree_sizes;
    }

    /**
     * String representation of the pattern: numbers with at least two 
     * digits are enclosed in parentheses, one-digit numbers are simply listed.
     * Example: <tt>017(12)4</tt>.
     * 
     * @return a String representation of the extended phyletic pattern
     */
    public String getPatternString()
    {
        return getPatternString(profile);
    }
    
    /**
     * String representation of the pattern: numbers with at least two 
     * digits are enclosed in parentheses, one-digit numbers are simply listed.
     * Example: <tt>017(12)4</tt>.
     * 
     * @param profile array of sizes 
     * @return a String representation of the extended phyletic pattern
     */
    public static String getPatternString(int[] profile)
    {
        StringBuffer sb = new StringBuffer();
        for (int i=0; i<profile.length; i++)
        {
            if (profile[i]<0)
                sb.append('?');
            else if (profile[i]<10)
                sb.append(Integer.toString(profile[i]));
            else
            {
                sb.append('(');
                sb.append(Integer.toString(profile[i]));
                sb.append(')');
            }
        }
        return sb.toString();
    }

    /**
     * Computes the ancestral reconstruction by asymmetric Wagner parsimony for this profile. 
     * For details, see Mikl&oacute;s Cs&#369;r&ouml;s <q>Ancestral reconstruction by asymmetric 
     * Wagner parsimony over continuous characters and squared parsimony over distributions</q>, 
     * at <em>RECOMB Workshop on Comparative Genomics</em> (Paris, France, October 2008), 
     * Springer Lecture Notes in Computer Science, volume 5267.
     * 
     * @param tree phylogeny for the terminal taxa: rates are not used, only the fixed postorder traversal
     * @param gain_penalty penalty for gain, relative to loss (i.e., loss prenalty is 1.0)
     * @return array of values for all the nodes of the tree (in the fixed traversal order), ancestral values minimize asymmetric Wagner parsimony
     */
    public int[] computeWagnerParsimony(TreeWithRates tree, double gain_penalty)
    {
        //Verbose.message("PP.cWP "+this.getPatternString());
        NodeWithRates[] nodes = tree.getDFT();
        double[][] subtree_slope = new double[nodes.length][];
        int[][] subtree_breakpoint = new int[nodes.length][];
        double[] subtree_shift = new double[nodes.length];
        double[][] stem_slope = new double[nodes.length][];
        int[][] stem_breakpoint = new int[nodes.length][]; 
        double[] stem_shift = new double[nodes.length];
        int[] stem_left_bp = new int[nodes.length];
        int[] stem_right_bp = new int[nodes.length];
        
        for (int node_idx=0; node_idx<nodes.length; node_idx++)
        {
            NodeWithRates N = nodes[node_idx];
            if (N.isLeaf())
            {
                // compute stem weight function
                if (profile[node_idx]<0)
                {
                    stem_slope[node_idx]=new double[1];
                    stem_breakpoint[node_idx]=new int[1];
                    stem_slope[node_idx][0]=0.0;
//                    stem_slope[node_idx][1]=0.0;
//                    stem_breakpoint[node_idx][1]=0;
                    stem_shift[node_idx]=0.0;
                } else
                {
                    stem_slope[node_idx]=new double[2];
                    stem_breakpoint[node_idx]=new int[2];
                    stem_slope[node_idx][0]=-gain_penalty;
                    stem_slope[node_idx][1]=1.0;

                    stem_breakpoint[node_idx][1]=profile[node_idx];
                    stem_shift[node_idx]=gain_penalty*profile[node_idx];
                }
                
                //Verbose.message("PP.cWP leaf "+node_idx+"/"+N.getTaxonName()+"\t"+profile[node_idx]);
            } else
            {
                // compute subtree weight functions
                int num_children = N.getNumChildren();
                Vector<Double> slopeV = new Vector<Double>();
                Vector<Integer> breakpointV = new Vector<Integer>();
                subtree_shift[node_idx]=0.0;
                
                for (int ci=0; ci<num_children; ci++)
                {
                    int child_idx = tree.getChildIndex(node_idx,ci);
                    sumPiecewiseLinear(slopeV, breakpointV, stem_slope[child_idx], stem_breakpoint[child_idx]);
                    subtree_shift[node_idx]+=stem_shift[child_idx];
                }
                int k = slopeV.size();
                subtree_slope[node_idx]=new double[k];
                subtree_breakpoint[node_idx]=new int[k];
                for (int i=0; i<k; i++)
                {
                    subtree_slope[node_idx][i] = slopeV.get(i).doubleValue();
                    subtree_breakpoint[node_idx][i] = breakpointV.get(i).intValue();
                }
                //if (Verbose.isVerbose())
                //{
                //   for (int i=0; i<subtree_slope[node_idx].length; i++)
                //   {
                //       Verbose.message("PP.cWP subtree "+node_idx+"/"+N.getTaxonName()+"\t"+i
                //               +"\t"+subtree_slope[node_idx][i]
                //               +"\t"+subtree_breakpoint[node_idx][i]);
                //   }
                //}
                if (!N.isRoot())
                {
                    // compute stem weight function
                    int i_left = 0;
                    for (int i=1; i<k; i++)
                        if (subtree_slope[node_idx][i]>=-gain_penalty)
                        {
                            i_left=i;
                            break;
                        }
                    
                    // what is the shift here?
                    double phi = subtree_shift[node_idx];
                    for (int i=1; i<=i_left; i++)
                    {
                        if (i==1)
                            phi += subtree_slope[node_idx][0]*subtree_breakpoint[node_idx][1];
                        else
                            phi += subtree_slope[node_idx][i-1]*(subtree_breakpoint[node_idx][i]-subtree_breakpoint[node_idx][i-1]);
                    }
                    stem_shift[node_idx] = phi + gain_penalty * subtree_breakpoint[node_idx][i_left];
                    
                    int i_right = k-1;
                    while (i_right>=i_left && subtree_slope[node_idx][i_right]>=1.0)
                        i_right--;
                    i_right++;
                    //Verbose.message("PP.cWP stem @ "+i_left+".."+i_right+" ["+k+"] shift "+stem_shift[node_idx]);
                    stem_left_bp[node_idx]=i_left;
                    stem_right_bp[node_idx]=i_right;
                    
                    stem_slope[node_idx]=new double[i_right-i_left+2];
                    stem_breakpoint[node_idx]=new int[i_right-i_left+2];
                    stem_slope[node_idx][0]=-gain_penalty;
                    for (int i=i_left; i<=i_right; i++)
                    {
                        stem_slope[node_idx][i-i_left+1]=subtree_slope[node_idx][i];
                        stem_breakpoint[node_idx][i-i_left+1]=subtree_breakpoint[node_idx][i];
                    }
                    stem_slope[node_idx][i_right-i_left+1]=1.0;
                    //if (Verbose.isVerbose())
                    //{
                    //   for (int i=0; i<stem_slope[node_idx].length; i++)
                    //   {
                    //       Verbose.message("PP.cWP stem "+node_idx+"/"+N.getTaxonName()+"\t"+i+"\t"+stem_slope[node_idx][i]+"\t"+stem_breakpoint[node_idx][i]);
                    //   }
                    //}

                }
            }
        } // for all nodes
        int[] retval = new int[nodes.length];
        
        // find minimum at root
        int imin=1;
        while (subtree_slope[nodes.length-1][imin]<0) imin++;
        retval[nodes.length-1] = subtree_breakpoint[nodes.length-1][imin];
        NodeWithRates root=tree.getRoot();
        String root_name = (root.getName()==null?"root":root.getName());
        //Verbose.message("PP.cWP solution "+root_name+"\t"+retval[nodes.length-1]+"\t// "+getPatternString());        
        
        for (int node_idx=nodes.length-2; node_idx>=0; node_idx--)
        {
            NodeWithRates N = nodes[node_idx];
            int parent_idx = tree.getParentIndex(node_idx);
            int y = retval[parent_idx];
            
            if (N.isLeaf())
            {
                if (profile[node_idx]<0)
                {
                    retval[node_idx] = retval[parent_idx];
                } else
                {
                    retval[node_idx]=profile[node_idx];
                }
            }
            else
            {
                int x0 = subtree_breakpoint[node_idx][stem_left_bp[node_idx]];
                int x1 = subtree_breakpoint[node_idx][stem_right_bp[node_idx]];
                if (y<x0)
                    retval[node_idx]=x0;
                else if (y>x1)
                    retval[node_idx]=x1;
                else
                    retval[node_idx]=y;
                //Verbose.message("PP.cWP solution "+node_idx+"/"+N.getTaxonName()+"\t"+retval[node_idx]+"\tprn "+y+"\t["+x0+", "+x1+"]");
            }
            
        }
        
        if (false)
        { // check by regular Sankoff
            int max_value=0;
            for (int i=0; i<profile.length; i++)
                if (profile[i]>max_value)
                    max_value=profile[i];
            double[][] score = new double[nodes.length][max_value+1];
            int[][] best_value = new int[nodes.length][max_value+1];
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                NodeWithRates N = nodes[node_idx];
                if (N.isLeaf())
                {
                    Arrays.fill(score[node_idx],Double.POSITIVE_INFINITY);
                    score[node_idx][profile[node_idx]]=0.0;
                } else
                {
                    int num_children = N.getNumChildren();
                    for (int v=0; v<=max_value; v++)
                    {
                        double s = 0.0;
                        for (int ci=0; ci<num_children; ci++)
                        {
                            int child_idx = tree.getChildIndex(node_idx,ci);
                            int best_x = 0;
                            double best_score = Double.POSITIVE_INFINITY;
                            for (int x=0; x<=max_value; x++)
                            {
                                double pty = (x<v?(v-x):(x-v)*gain_penalty);
                                double stem_score = pty + score[child_idx][x];
                                //Verbose.message("PP.cWP Sankoff "+child_idx+"/"+nodes[child_idx].getTaxonName()+"\tchange "+v+"\t->"+x+"\t"+stem_score);
                                
                                if (stem_score<best_score)
                                {
                                    best_x = x;
                                    best_score=stem_score;
                                }
                            }
                            best_value[child_idx][v]=best_x;
                            s+=best_score;
                            //Verbose.message("PP.cWP Sankoff "+child_idx+"/"+nodes[child_idx].getTaxonName()+"\tbest "+v+"\t->"+best_x+"\t"+best_score);
                        }
                        score[node_idx][v]=s;
                        //Verbose.message("PP.cWP Sankoff "+node_idx+"/"+nodes[node_idx].getTaxonName()+"\t"+v+"\t"+s);
                    }
                }
            } // for node_idx
            // backtrack
            int root_best=0;
            for (int v=1; v<=max_value; v++)
                if (score[nodes.length-1][v]<score[nodes.length-1][root_best])
                    root_best=v;
            int[] sankoff = new int[nodes.length];
            sankoff[nodes.length-1]=root_best;
            for (int node_idx=nodes.length-2; node_idx>=0; node_idx--)
            {
                int parent_idx = tree.getParentIndex(node_idx);
                sankoff[node_idx] = best_value[node_idx][sankoff[parent_idx]];
            }
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                Verbose.message("PP.cWP Sankoff "+node_idx+"/"+nodes[node_idx].getTaxonName()
                    +"\t"+retval[node_idx]+"\t"+sankoff[node_idx]
                    +"\t"+(retval[node_idx]!=sankoff[node_idx]?"XXXXXXXXXX":"")+"\t// "+score[node_idx][retval[node_idx]]+"\t"+score[node_idx][sankoff[node_idx]]);
            }
        }
        
        return retval;
    }
    
    private void sumPiecewiseLinear(Vector<Double> slopes, Vector<Integer> breakpoints, double[] slopes2, int[] breakpoints2)
    {
        if (slopes.size()==0)
        { // first call
            for (int i=0; i<slopes2.length; i++)
            {
                slopes.add(new Double(slopes2[i]));
                breakpoints.add(new Integer(breakpoints2[i]));
            }
        } else
        {
            int n1 = slopes.size();
            double[] slopes1 = new double[n1];
            int[] breakpoints1 = new int[n1];
            for (int i=0; i<n1; i++)
            {
                slopes1[i] = slopes.get(i).doubleValue();
                breakpoints1[i] = breakpoints.get(i).intValue();
            }
            slopes.clear();
            breakpoints.clear();

            slopes.add(new Double(slopes1[0]+slopes2[0]));
            breakpoints.add(new Integer(0)); // dummy placeholder
            int n2 = slopes2.length;
            
            int i1=1; 
            int i2=1;
            while (i1<n1 || i2<n2)
            {
                if (i1==n1 || (i2<n2 && breakpoints2[i2]<breakpoints1[i1]))
                {
                    // add breakpoint2
                    int x = breakpoints2[i2];
                    double a = slopes1[i1-1]+slopes2[i2];
                    breakpoints.add(new Integer(x));
                    slopes.add(new Double(a));
                    i2++;
                } else if (i2==n2 || (i1<n1 && breakpoints1[i1]<breakpoints2[i2]))
                {
                    // add breakpoint1
                    int x = breakpoints1[i1];
                    double a = slopes1[i1]+slopes2[i2-1];
                    breakpoints.add(new Integer(x));
                    slopes.add(new Double(a));
                    i1++;
                } else 
                { // equality
                    int x = breakpoints1[i1];
                    double a = slopes1[i1]+slopes2[i2];
                    breakpoints.add(new Integer(x));
                    slopes.add(new Double(a));
                    i1++;
                    i2++;
                }
            }
            //if (Verbose.isVerbose())
            //{
            //    for (int i=0; i<n1; i++)
            //        Verbose.message("PL.sPL AAAA "+slopes1[i]+"\t"+breakpoints1[i]);
            //    for (int i=0; i<n2; i++)
            //        Verbose.message("PP.sPL BBBB "+slopes2[i]+"\t"+breakpoints2[i]);
            //    for (int i=0; i<slopes.size(); i++)
            //        Verbose.message("PP.sPL cccc "+slopes.get(i)+"\t"+breakpoints.get(i));
            //}
            
        }
    }
    
    
    /**
     * Computes most parsimonious history for this profile using Dollo criterion (every family is gained only once). 
     * 
     * @param tree phylogeny for the terminal taxa: only the fixed node traversal is used 
     * 
     * @return array with 0 and 1 for family presence/absence: one entry for each node
     */ 
    public int[] computeDolloParsimony(TreeWithRates tree)
    {
        int num_nodes = tree.getNumNodes();
        // whether at least one leaf in the subtree has this family
        boolean[] some_present = new boolean[num_nodes];
        for (int node_idx=0; node_idx<num_nodes; node_idx++)
        {
            NodeWithRates N = tree.getNode(node_idx);
            if (N.isLeaf())
            {
                some_present[node_idx]=(profile[node_idx]>0);// includes missing data
            } else 
            {
                int num_children = N.getNumChildren();
                for (int i=0; i<num_children; i++)
                {
                    int child_node_idx = tree.getChildIndex(node_idx, i);
                    if (some_present[node_idx]=some_present[child_node_idx]) // not == but = 
                        break;
                }
            }
        } // for node
        
        // do the dynamic programming
        int[] score0 = new int[num_nodes];
        int[] score1 = new int[num_nodes];
        boolean[] loss = new boolean[num_nodes]; // whether the 1->0 is better than 1->1 on this edge 
        boolean[] gain = new boolean[num_nodes];
        // bottom-up
        
        //
        // Integer.MAX_VALUES is used for 'infinity'
        for (int node_idx=0; node_idx<num_nodes; node_idx++)
        {
            NodeWithRates N = tree.getNode(node_idx);
            if (N.isLeaf())
            {
                score0[node_idx]=score1[node_idx]=0;
                if (profile[node_idx]==0)
                {
                    score1[node_idx]=Integer.MAX_VALUE;
                }
                else if (profile[node_idx]>0)
                {
                    score0[node_idx]=Integer.MAX_VALUE;
                }
            } else
            {
                int num_children = N.getNumChildren();
                // count how many subtrees have "some_present"
                int num_present = 0;
                for (int i=0; i<num_children; i++)
                    if (some_present[tree.getChildIndex(node_idx,i)])
                        num_present++;
                if (num_present>1)
                {
                    // clear-cut case, here we must have "present"
                    score0[node_idx]=Integer.MAX_VALUE;
                    score1[node_idx]=0;
                    for (int i=0; i<num_children; i++)
                    {
                        int child_node_idx = tree.getChildIndex(node_idx,i);
                        if (some_present[child_node_idx])
                        {
                            // here we cannot switch to '0' because then we would have another 0->1 in the subtree
                            score1[node_idx]+=score1[child_node_idx];
                            loss[child_node_idx]=false;
                        } else
                        {
                            // here we could have 0 or 1 but the optimum will surely put 0 at the child who will be all-0
                            score1[node_idx]++;
                            loss[child_node_idx]=true;
                        }
                    } // for all children
                }  else if (num_present==1)
                {
                    // now only one guy has some_present, we can be 0 or 1
                    score0[node_idx]=score1[node_idx]=0;
                    for (int i=0; i<num_children; i++)
                    {
                        int child_node_idx = tree.getChildIndex(node_idx,i);
                        if (some_present[child_node_idx])
                        {
                            // 1->0 is not allowed here
                            score1[node_idx]+=score1[child_node_idx];
                            loss[child_node_idx]=false;
                            // but 0->1 and 0->0 are both possible
                            if (score0[child_node_idx]-1>score1[child_node_idx]) // important to use '-1' on left instead of '+1' on right to avoid performing Integer.MAXINT+1
                            {
                                score0[child_node_idx] += score1[child_node_idx]+1;
                                gain[child_node_idx]=true;
                            } else
                            {
                                score0[node_idx] += score0[child_node_idx];
                                gain[child_node_idx]=false;
                            }
                        } else
                        {
                            // nobody in this subtree, so we better label it with '0'
                            score1[node_idx]++;
                            loss[child_node_idx]=true;
                            // score0 does not change
                            gain[child_node_idx]=false;
                        }
                    } // children
                } else 
                {
                    // now everybody is absent below us
                    score0[node_idx]=0;
                    score1[node_idx]=num_children;
                    for (int i=0; i<num_children; i++)
                    {
                        int child_node_idx = tree.getChildIndex(node_idx,i);
                        loss[child_node_idx]=true;
                        gain[child_node_idx]=false;
                    }
                }
            } // non-leaf
        } // for all nodes 
        
        // find best labels top-down fashion
        int label[] = new int[num_nodes];
        
        // root label
        if (score0[num_nodes-1]>score1[num_nodes-1]) // egality does not happen with a trifurcation
            label[num_nodes-1]=1;
        else 
            label[num_nodes-1]=0;
        for (int node_idx=num_nodes-2; node_idx>=0; node_idx--)
        {
            int parent_node_idx = tree.getParentIndex(node_idx);
            boolean present 
                = (label[parent_node_idx]==0 && gain[node_idx])
                || (label[parent_node_idx] == 1 && !loss[node_idx]);
            if (present)
                label[node_idx]=1;
            else 
                label[node_idx]=0;
        }
        
        return label;
    }    
    
    
    public LowerConditionals getLowerConditional(SurvivalProbabilityFactory.Rates survival)
    {
        return new LowerConditionals(survival);
    }
    
    public UpperConditionals  getUpperConditionals(LowerConditionals LC, DiscreteDistribution root_prior)
    {
        return new UpperConditionals(LC,root_prior);
    }
        
    /**
     * Machine precision for double
     */
    private static final double EPS=1e-14;
    /**
     * Computes the binomial inversion for a series of values. 
     * Given the array a[0..], the formula b[n] = sum_{i=0}^n (n choose i) p^{n-1} (1-)^i a[i] 
     * defines the entries of the array b[0..]. 
     * Returns the array a.
     *
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
                            Verbose.message("PP.bI\t"+n+"\tb "+b[n]+"\tsub "+subtract+"\t= "+a[n]+"\tdiv "+d[n]+"\t// ++++ CLOSE\th "+h);
                        } else 
                        {
                            a[n] = (b[n]-subtract)/d[n];
                            Verbose.message("PP.bI\t"+n+"\tb "+b[n]+"\tsub "+subtract+"\t= "+a[n]+"\tdiv "+d[n]+"\th "+h);
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
                Verbose.message("PP.bI iteration "+rr+"\t"+delta);
            }
            for (int i=0; i<x.length; i++)
            {
                Verbose.message("PP.bI "+i+"\tx "+x[i]+"\ta "+a[i]+"\tb "+bb[i]+"\t// "+b[i]+"\t// delta "+delta+"\t"+rr);
            }
            a = x;
        }
        //System.exit(99);
        return a;
    }
    
    
    public class LowerConditionals 
    {
        private LowerConditionals(SurvivalProbabilityFactory.Rates survival)
        {
            this.survival = survival;
            rate_tree = survival.getTree();
            initDataStructures();
            recomputeAll();
        }
        
        private SurvivalProbabilityFactory.Rates survival;
        private TreeWithRates rate_tree;
        private int[] subtree_sizes;
        
        private void initDataStructures()
        {
            NodeWithRates[] nodes = rate_tree.getDFT();
            likelihood = new double[nodes.length][];
            uncorrected_likelihood = new double[nodes.length][];
            subtree_sizes = computeSubtreeSizes(rate_tree);
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                likelihood[node_idx] = new double[subtree_sizes[node_idx]+1];
                uncorrected_likelihood[node_idx] = new double[subtree_sizes[node_idx]+1];
            }
        }

        private double[][] likelihood;
        private double[][] uncorrected_likelihood;
        
        
        public void recomputeAtNode(int node_idx)
        {
            NodeWithRates N = rate_tree.getNode(node_idx);
            double[] piciL = uncorrected_likelihood[node_idx];
            java.util.Arrays.fill(piciL,0.0);
            
            if (N.isLeaf())
            {
                java.util.Arrays.fill(likelihood[node_idx],0.0);
                int v = profile[node_idx];
                piciL[v]=likelihood[node_idx][v] = 1.0;
            } else
            {
                double extinct = survival.getExtinction(node_idx);
                int num_children = N.getNumChildren();

                Verbose.message("PP.LC.rAN "+getPatternString()+"\t"+node_idx+"/"+rate_tree.getNode(node_idx));
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
                            //Verbose.message("PP.LC.rAN "+node_idx+"\tn "+n+"\tci "+ci+"\tm "+m+"\tLc "+Lc[m]+"\tpmn "+pmn[m]+"\t// "+rate_tree.getNode(child_idx));
                            Ledge += z;
                        }
                        //Verbose.message("PP.LC.rAN "+node_idx+"\tn "+n+"\tci "+ci+"\tLedge "+Ledge);
                        piciL[n] *= Ledge;
                    }
                    //Verbose.message("PP.LC.rAN "+node_idx+"\tn "+n+"\tpiciL "+piciL[n]);
                    
                }
                likelihood[node_idx] = binomialInversion(piciL, extinct, piciL.length-1);  
                for (int n=0; n<piciL.length; n++)
                {
                    if (Double.isInfinite(likelihood[node_idx][n]) || Double.isNaN(likelihood[node_idx][n]))
                    {
                        for (int j=0; j<piciL.length; j++)
                        {
                            System.out.println("#### NUMERICAL_ERROR PP.LC.rAN "+node_idx+"/"+rate_tree.getNode(node_idx).getTaxonName()+"\t["+j+"]\tL "+likelihood[node_idx][j]+"\tl "+piciL[j]);
                        }
                        System.exit(2008);
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
    }
    
    public class UpperConditionals 
    {
        private UpperConditionals(LowerConditionals LC, DiscreteDistribution root_prior)
        {
            this.LC =LC;
            this.root_prior = root_prior;
            initDataStructures();
        }
        
        private LowerConditionals LC;
        private DiscreteDistribution root_prior;
        
        double[][] likelihood;
        
        private void initDataStructures()
        {
            NodeWithRates[] nodes = LC.rate_tree.getDFT();
            likelihood = new double[nodes.length][];
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                likelihood[node_idx] = new double[LC.subtree_sizes[node_idx]+1];
            }
        }
        
        public void recomputeAtNode(int node_idx)
        {
            TreeWithRates rate_tree = LC.rate_tree;
            NodeWithRates N = rate_tree.getNode(node_idx);
            int M_y = likelihood[node_idx].length-1;
            if (N.isRoot())
            {
                likelihood[node_idx] = LC.survival.survivingRootPrior(root_prior).getDistribution(M_y);
            } else
            { // non-root node
                int parent_idx = rate_tree.getParentIndex(node_idx);
                NodeWithRates Parent = rate_tree.getNode(parent_idx);
                int num_children = Parent.getNumChildren();
                int M_x = LC.subtree_sizes[parent_idx];
                double [] lex = new double[M_x+1];
                double [][] Lchild = new double[num_children][];
                for (int ci=0; ci<num_children; ci++)
                {
                    int child_idx = rate_tree.getChildIndex(parent_idx,ci);
                    if (child_idx != node_idx)
                        Lchild[ci] = LC.getLikelihoods(child_idx);
                }
                for (int n=0; n<=M_x; n++)
                {
                    lex[n] = 1.0;
                    for (int ci=0; ci<num_children; ci++)
                    {
                        int child_idx = rate_tree.getChildIndex(parent_idx,ci);
                        if (child_idx != node_idx)
                        {
                            double[] pmn = LC.survival.getSurvival(child_idx,n);
                            double Ledge=0.0;
                            for (int m=0; m<Lchild[ci].length; m++)
                            {
                                double z=pmn[m]*Lchild[ci][m];
                                Ledge += z;
                            }
                            lex[n] *= Ledge;
                        }
                    }
                }
                double[][] w = new double[M_y+1][M_x+1];
                for (int n=0; n<=M_x; n++)
                {
                    double[] pmn = LC.survival.getSurvival(node_idx,n);
                    for (int m=0; m<=M_y; m++)
                        w[m][n] = pmn[m];
                }
                double[][] W=new double[M_y+1][];
                double extinct = LC.survival.getExtinction(node_idx);
                for (int m=0; m<=M_y; m++)
                    W[m]=binomialInversion(w[m],extinct,M_x);
                double[] parentL = likelihood[parent_idx];
                for (int m=0; m<=M_y; m++)
                {
                    double Lm = 0.0;
                    for (int n=0; n<=M_x; n++)
                    {
                        double z = W[m][n]*parentL[n];
                        Lm += z;
                    }
                    likelihood[node_idx][m] = Lm;
                }
            }
        }
        
        public void recomputeAll()
        {
            int num_nodes = LC.rate_tree.getNumNodes();
            for (int node_idx = num_nodes-1; node_idx >=0; node_idx--)
                recomputeAtNode(node_idx);
        }
        
        public double[] getLikelihoods(int node_idx)
        {
            return likelihood[node_idx];
        }

        public PhyleticProfile.LowerConditionals getLC() {
            return LC;
        }
        
        public double[] getPosteriors(int node_idx)
        {
            double[] p = new double[likelihood[node_idx].length];
            double[] L = LC.getLikelihoods(node_idx);
            double sum = 0.0;
            double[] U = likelihood[node_idx];
            for (int m=0; m<p.length; m++)
            {
                p[m] = L[m]*U[m];
                sum += p[m];
            }
            for (int m=0; m<p.length; m++)
                p[m] /= sum;
                    
            return p;
        }
    }
    
}
