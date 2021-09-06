/*
 * Likelihood.java
 *
 * Created on September 9, 2005, 12:37 AM
 */

package ca.umontreal.iro.evolution.genecontent;

import java.util.Hashtable;
//import java.math.BigDecimal;
import ca.umontreal.iro.banality.Verbose;

/**
 *
 * @author  csuros
 * @deprecated
 */
public class Likelihood {
    private static final double EPS=1e-14;
    
    /** Creates a new instance of Likelihood */
    public Likelihood(NodeWithRates root, String[] gene_names, int[][] copies) {
        this.root=root;
        this.gene_names=gene_names;
        this.copies=copies;
        init();
    }
    
    private NodeWithRates root;
    private String[] gene_names;
    private int[][] copies;
    
    /**
     * Private index variables for fast traversal of the tree.
     */
    private NodeWithRates[] nodes;
    private Hashtable node_indices;
    private int[][] idx_children;
    private int[] idx_parent;
    private NodeWithRates[] leaves;
    private Hashtable leaf_indices;
    /**
     * Sets up the private index variables.
     */
    private void init(){
        // integer indexes for nodes 
        nodes = root.getDFT();
        node_indices =  new Hashtable();
        for (int idx=0; idx<nodes.length; idx++)
            node_indices.put(nodes[idx], new Integer(idx));
        // integer indexes for children & parent
        idx_children = new int[nodes.length][];
        idx_parent = new int[nodes.length];
        for (int idx=0; idx<nodes.length; idx++){
            NodeWithRates N=nodes[idx];
            int n_child = N.getNumChildren();
            idx_children[idx] = new int[n_child];
            for (int child_idx=0; child_idx<n_child; child_idx++)
                idx_children[idx][child_idx]=((Integer)node_indices.get(N.getChild(child_idx))).intValue();
            if (N.isRoot()){
                idx_parent[idx]=-1;
            } else {
                idx_parent[idx]= ((Integer)node_indices.get(N.getParent())).intValue();
            }
        }
        // integer indexes for leaves
        leaves = root.getLeaves();
        leaf_indices =  new Hashtable();
        for (int idx=0; idx<leaves.length; idx++)
            leaf_indices.put(leaves[idx], new Integer(idx));
        
        calculateSubtreeSizes();
        extinction_prob = new double[nodes.length];
        newParameters();
    }
    
    /**
     * Extinction probabilities for a gene at a node.
     */
    private double[] extinction_prob;
    /**
     * Allocates and fills up auxiliary data structures..
     * Must be called when rate or branch length parameters change (not for every family though).
     */
    void newParameters(){
        for (int node_idx = 0; node_idx<nodes.length; node_idx++){
            NodeWithRates N=nodes[node_idx];
            if (N.isLeaf())
                extinction_prob[node_idx]=0.0;
            else {
                int num_children = N.getNumChildren();
                double p=1.0;
                for (int child_idx=0; child_idx<num_children;child_idx++){
                    NodeWithRates C =(NodeWithRates) N.getChild(child_idx);
                    double beta   = C.getBeta();
                    double mu     = C.getLossRate();
                    double lambda = C.getDuplicationRate();
                    double D      = extinction_prob[idx_children[node_idx][child_idx]];
                    double mb = mu*beta;
                    double lb = lambda*beta;
                    double y = (1.0-mb)*(1.-lb)*D/(1.0-lb*D);
                    p *= (mb+y);
                }
                extinction_prob[node_idx]=p;
                //Verbose.message("L.nP node "+N.getTaxonName()+" D="+p);
            }
        }
    }
    
    /**
     * Number of copies in the subtree below this node: [i][j] is for family i and node j
     */
    private int[][] copies_in_subtree; 
    /**
     * Allocates and fills up copies_in_subtree[] array
     */
    private void calculateSubtreeSizes(){
        copies_in_subtree=new int[copies.length][nodes.length];
        for (int family_idx=0; family_idx<copies.length; family_idx++){
            for (int node_idx = 0; node_idx<nodes.length; node_idx++){
                NodeWithRates N=nodes[node_idx];
                if (N.isLeaf()){
                    int leaf_idx =  ((Integer)leaf_indices.get(N)).intValue();
                    copies_in_subtree[family_idx][node_idx]=copies[family_idx][leaf_idx];
                }
                else {
                    int num_children = N.getNumChildren();
                    int M=0;
                    for (int child_idx=0; child_idx<num_children;child_idx++){
                        M += copies_in_subtree[family_idx][idx_children[node_idx][child_idx]];
                    }
                    copies_in_subtree[family_idx][node_idx]=M;
                }
            } // for nodes
            //Verbose.message("L.cSS family @"+family_idx+" "+gene_names[family_idx]+" size "+copies_in_subtree[family_idx][nodes.length-1]);
            
        } // for families
    }
    
    /**
     * @return array of survival probabilities for given family.
     *     p[y][m][n] is the probability on the branch leading to y, for n copies at parent and m surviving copies at y.
     */
    private double[][][] getSurvivalProbabilities(int family_idx){
        //Verbose.message("L.gSP "+family_idx);
        double[][][] p = new double[nodes.length][][];
        
        for (int node_idx = 0; node_idx<nodes.length-1; node_idx++){
            //Verbose.message("L.gSP ... "+family_idx+" "+node_idx);
            NodeWithRates N=nodes[node_idx];
            int M_parent=copies_in_subtree[family_idx][idx_parent[node_idx]];
            int M = copies_in_subtree[family_idx][node_idx];
            //Verbose.message("L.gSP .... "+family_idx+" "+node_idx+" "+M_parent+" "+M);
            double[][] P = new double[M_parent+1][M+1];
            //Verbose.message("L.gSP ..... "+family_idx+" "+node_idx);
            double beta = N.getBeta();
            double lambda = N.getDuplicationRate();
            double mu = N.getLossRate();
            double kappa = N.getTransferRate();
            double D = extinction_prob[node_idx];
            
            double theta = kappa/lambda; 
            double lb = lambda*beta;
            double mb = mu*beta;
            double Dlb = D*lb;
            double Dlb_1 = 1.0-Dlb;
            
            double x = (1.0-D)* lb/Dlb_1;
            double H0 = Math.pow((1.0-lb)/Dlb_1,theta);
            double G0 = 1.0-(1.0-mb)*(1.0-D)/Dlb_1;
            double G1 = (1.0-mb)*(1.0-lb)*(1.0-D)/(Dlb_1*Dlb_1);
            double Gcorr = G1-x*G0;

            /*
            {
                double[] G=new double[40];
                G[0]=G0;
                G[1]=G1;
                for (int j=2; j<40; j++){
                    G[j]=G[j-1]*x;
                }
                double sum=0.;
                for (int j=0; j<40; j++){
                    sum += G[j];
                    Verbose.message("L.gSP G["+j+"] "+G[j]+"\t"+sum);
                }
            }
            
            {
                double[] H=new double[40];
                H[0]=H0;
                for (int m=1; m<40; m++){
                    H[m]=H[m-1]*((theta+m-1.)/((double)m))*x;
                }
                double sum=0.;
                for (int j=0; j<40;j++){
                    sum += H[j];
                    Verbose.message("L.gSP H["+j+"] "+H[j]+"\t"+sum);
                }
            }
            */

            
            P[0][0] = H0;
            for (int m=1; m<=M; m++)
                P[0][m] = P[0][m-1]*(theta+m-1.)/((double)m)*x;
            for (int n=1; n<=M_parent; n++){
                P[n][0] = P[n-1][0] * G0;
                if (M>0) P[n][1] = G0*P[n-1][1]+G1*P[n-1][0];
                for (int m=2; m<=M; m++){
                    P[n][m] = x*P[n][m-1]+G0*P[n-1][m]+Gcorr*P[n-1][m-1];
                }
            }
            
            
            /*
            { // print data
                for (int n=0; n<P.length; n++){
                    double sum = 0.0;
                    for (int m=0; m<P[n].length; m++){
                        sum += P[n][m];
                        Verbose.message("L.gSP "+N.getTaxonName()+"\tp("+m+" | "+n+") "+ P[n][m]+"\t"+sum);
                    }
                }
            }
             */
             
            //Verbose.message("L.gSP "+N.getTaxonName()+" done.");
            p[node_idx]=P;
        }
        
        
        
        return p;
    }
    
    
    /**
     * @return array of conditional likelihoods for given family. L[m] is for root's m surviving copies.
     * 
     */
    private double[] getConditionalLikelihoods(int family_idx){
        //Verbose.message("L.gCL "+family_idx);
        double[][][] p = getSurvivalProbabilities(family_idx);
        //Verbose.message("L.gCL ... "+family_idx);

        double[][] L = new double[nodes.length][]; //BigDecimal[][] BD_L=new BigDecimal[nodes.length][];
        
        for (int node_idx=0; node_idx<nodes.length; node_idx++){
            NodeWithRates N=nodes[node_idx];
            int M = copies_in_subtree[family_idx][node_idx];            
            L[node_idx]=new double[M+1]; //BD_L[node_idx]=new BigDecimal[M+1];
            
            if (N.isLeaf()){
                for (int m=0; m<=M; m++){
                    L[node_idx][m]=0.0; //BD_L[node_idx][m]=new BigDecimal(0.0);
                }
                L[node_idx][M]=1.0; //BD_L[node_idx][M]= new BigDecimal(1.0);
                
            } else {
                int num_children = N.getNumChildren();
                for (int n=0; n<=M; n++){
                    L[node_idx][n]=1.0; //BD_L[node_idx][n]=new BigDecimal(1.0);
                    L[node_idx][n]=1.0; //BD_L[node_idx][n]=new BigDecimal(1.0);
                    for (int ci=0; ci<num_children; ci++){
                        int child_idx = idx_children[node_idx][ci];
                        int Mc = copies_in_subtree[family_idx][child_idx];
                        double ll = 0.0; //BigDecimal BD_ll = new BigDecimal(0.0);
                        for (int m=0; m<=Mc; m++){
                            double q = p[child_idx][n][m]*L[child_idx][m];
                            ll += q;
                            //BigDecimal BD_q = (new BigDecimal(p[child_idx][n][m])).multiply(BD_L[child_idx][m]); 
                            //BD_ll = BD_ll.add(BD_q);
                        }
                        L[node_idx][n] *= ll;
                        //BD_L[node_idx][n] = BD_L[node_idx][n].multiply(BD_ll);
                    }
                } // for n
                /*
                {
                    double numerical_error = (BD_L[node_idx][0].doubleValue()==0.?L[node_idx][0]:Math.abs((L[node_idx][0]-BD_L[node_idx][0].doubleValue())/BD_L[node_idx][0].doubleValue()));
                    Verbose.message("L.gCL "+node_idx
                        +" L[0]\t "+L[node_idx][0]
                            +"; error "+numerical_error+" BD_L "+BD_L[node_idx][0].doubleValue());
                }
                */
                
                //Verbose.message("L.cCL ... "+N.getTaxonName()+" L[0] "+L[node_idx][0]);

                double D=extinction_prob[node_idx];
                double x = (1.0-D)/D;
                
                for (int n=1; n<=M; n++){
                    //Verbose.message("L.cCL ... "+N.getTaxonName()+" L["+n+"] "+L[node_idx][n]);
                    double y=Math.pow(D,n); //BigDecimal BD_y = new BigDecimal(y);
                    double sub_L = 0.0; //BigDecimal BD_sub_L=new BigDecimal(0.0);
                    for (int i=0; i<n; i++){
                        if (i>0){
                            double z = x*(n-i+1.0)/((double)i);
                            y = y*z;
                            //BD_y = BD_y.multiply(new BigDecimal(z));
                        }
                        //Verbose.message("L.cCL ... "+N.getTaxonName()+" L["+n+"] "+L[node_idx][n]+"\ty "+y+"\tx "+x);                    
                        sub_L += y*L[node_idx][i];
                        //BD_sub_L = BD_sub_L.add(BD_L[node_idx][i].multiply(BD_y));
                    }
                    
                    
                    
                    /*
                    double diff = L[node_idx][n]-sub_L; //BigDecimal BD_diff = BD_L[node_idx][n].subtract(BD_sub_L);
                    double diff2 = BD_diff.doubleValue();
                    double numerical_error = (diff2==0.?diff:Math.abs((diff-diff2)/diff2));
                    Verbose.message("L.gCL "+node_idx
                        +" ell["+n+"]\t "+(sub_L>L[node_idx][n]?"**** ":"     ")+L[node_idx][n]+" sub "+sub_L+" diff "+diff
                            +"; error "+numerical_error+" BD_L "+BD_L[node_idx][n].doubleValue()+" BD_sub "+BD_sub_L.doubleValue()+" BD_diff "+BD_diff.doubleValue());
                    BD_L[node_idx][n] = BD_diff;
                    */
                    /*
                    if (sub_L>L[node_idx][n] || sub_L < L[node_idx][n]-1.0){
                        L[node_idx][n]=0.0;
                    } else 
                        L[node_idx][n] -= sub_L;
                     */
                    
                    if (L[node_idx][n] != 0.){
                        double h = 1.0-sub_L/L[node_idx][n];
                        if (h<EPS){ // they are too close: L[node_idx][n] is less than (1+EPS)*sub_L --- does not make sense to take the difference.
                            L[node_idx][n]=0.0;
                        } else {
                            L[node_idx][n] -= sub_L;
                        }
                    }
                        
                } // for n
                
            } // node is not leaf
            
            
            /*
            {
                for (int n=0; n<=M;n++){
                    Verbose.message("L.gCL "+N.getTaxonName()+" L["+n+"] "+L[node_idx][n]);
                }
            }*/
            
            
        } // node_idx
        return L[nodes.length-1];
    }
    
    /**
     * @return the complete log-likelihood for the gene family
     */
    public double getLikelihood(int gene_family){
        //Verbose.message("L.gL "+gene_family);
        double[] L = getConditionalLikelihoods(gene_family);
        //Verbose.message("L.gL ... "+gene_family);

        double lambda=root.getDuplicationRate();
        double mu = root.getLossRate();
        double kappa = root.getTransferRate();
        double D = extinction_prob[nodes.length-1];
        double theta = kappa/lambda;
        double ml = lambda/mu;
        
        double f = Math.pow((1.-ml)/(1.-ml*D),theta);
        double x = (1.0-D)*ml/(1.-ml*D); 
        
        //Verbose.message("L.gL ... D "+D+" ml "+ml);
        
        double Lroot = f*L[0];
        double[] factors = new double[L.length];
        factors[0]=f;
        for (int j=1; j<L.length; j++){
            f *= x*((theta+j-1.0)/((double)j));
            Lroot += f*L[j];
            //Verbose.message("L.gL ... Lroot "+Lroot+" f "+f+" Ln "+L[j]);
            factors[j]=f;
        }
        /*
        {
            for (int j=0; j<L.length; j++){
                Verbose.message("L.gL "+gene_family+"/"+gene_names[gene_family]+" L["+j+"] "+L[j]);
            }
        }*/
        
        if (Lroot <= 0. || Lroot>1.0){
            Verbose.message("L.gL "+gene_family+"/"+gene_names[gene_family]+" numerror: Lroot "+Lroot+" x="+x);
            for (int j=0; j<L.length; j++){
                System.out.println("#**L.gL "+gene_family+" numerror L["+j+"] "+L[j]+"\t* "+factors[j]);
            }
        }
        
        return Lroot;
    }
}
