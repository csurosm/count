/*
 * NodeWithRates.java
 *
 * TreeNode Object with gain,loss, and duplication rates.
 *
 * Created on September 7, 2005, 12:02 PM
 */

package ca.umontreal.iro.evolution.genecontent;

import java.util.Vector;

import ca.umontreal.iro.banality.Verbose;

import ca.umontreal.iro.matek.DiscreteDistribution;
import ca.umontreal.iro.matek.Poisson;
import ca.umontreal.iro.matek.PointDistribution;
import ca.umontreal.iro.matek.NegativeBinomial;
import ca.umontreal.iro.matek.ShiftedGeometric;

import ca.umontreal.iro.evolution.TreeNode;


/**
 * A tree node with loss, duplication and gain/transfer rates.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class NodeWithRates extends TreeNode 
{
    public static final double DEFAULT_EDGE_LENGTH=0.1;
    public static final double DEFAULT_LOSS_RATE = 1.0;
    public static final double DEFAULT_DUPLICATION_RATE = 0.5;
    public static final double DEFAULT_GAIN_RATE = 0.4;

    private static final boolean SHOW_SIZE_DISTRIBUTIONS_IN_PARAMSTRING = false;

    /** Creates a new instance of NodeWithRates
     * 
     */
    public NodeWithRates() 
    {
        super();
        setDefaultParams();
    }
    /**
     * sets default node parameters: 1.0 for loss, 0.5 for duplication, 0.4 for transfer, and length 0.1
     */
    public void setDefaultParams()
    {
        loss_rate=DEFAULT_LOSS_RATE;
        duplication_rate=DEFAULT_DUPLICATION_RATE;
        transfer_rate=DEFAULT_GAIN_RATE;
        setLength(DEFAULT_EDGE_LENGTH);
    }
    
    private double duplication_rate; // lambda
    private double loss_rate; // mu
    private double transfer_rate; // kappa 
    
    public double getLossRate(){return loss_rate;}
    public double getDuplicationRate(){return duplication_rate;}
    public double getTransferRate(){return transfer_rate;}
    public void setLossRate(double rate){loss_rate=rate;}
    public void setDuplicationRate(double rate){duplication_rate=rate;}
    public void setTransferRate(double rate){transfer_rate=rate;}
        
    /**
     * Takes a tree already built with TreeNodes and returns and identical tree with 
     * NodeWithRates.
     * Either the default parameters are used for rates (in case the tree nodes are
     * not NodeWitRates objects), or the same parameters are copied.
     */
    public static NodeWithRates copyTree(TreeNode root)
    {
        NodeWithRates new_node=null;
        try 
        {
            new_node=copyTree(root, null);
        } catch (java.io.IOException E) 
        {
            // no such thing
        }
        return new_node;
    }
    
    /**
     * Takes a tree already built with TreeNodes and returns and identical tree with 
     * NodeWithRates. Also sets the gain/loss/transfer rates based on a previously saved 
     * tableRates() output read through a BufferedReader. 
     *
     * @param R a Reader for the string previously produced by tableRates or null.
     * If R is null, then either the default parameters are used for rates (in case the tree nodes are
     * not NodeWitRates objects), or the same parameters are copied. 
     * 
     */
    public static NodeWithRates copyTree(TreeNode old_node, java.io.BufferedReader R) throws java.io.IOException 
    {
        NodeWithRates new_node=new NodeWithRates();
        new_node.setLength(old_node.getLength());        
        new_node.setName(old_node.getName());
        new_node.setId(old_node.getId());
        int nc=old_node.getNumChildren();
        for (int ci=0; ci<nc; ci++)
        {
            NodeWithRates N=copyTree(old_node.getChild(ci), R);
            new_node.addChild(N);
        }
        if (R != null && !old_node.isRoot())
        {
            String line=null;
            do line=R.readLine(); while (line != null &&  (line.startsWith("#") || line.matches("^\\s*$")));
            if (line == null)
                throw new java.io.IOException("Premature end of rates table at node "+new_node);
            //System.out.println("#**NWR.tT line "+line);
            
            String[] fields=line.split("\\s+"); 
            double t=Double.parseDouble(fields[0]);
            double g=Double.parseDouble(fields[1]);
            double l=Double.parseDouble(fields[2]);
            double h=Double.parseDouble(fields[3]);
            new_node.duplication_rate=g;
            new_node.loss_rate=l;
            new_node.setLength(t);
            new_node.transfer_rate=h;
        } else if (R == null)
        {
            if (old_node instanceof NodeWithRates)
            {
                new_node.duplication_rate = ((NodeWithRates)old_node).duplication_rate;
                new_node.loss_rate = ((NodeWithRates)old_node).loss_rate;
                new_node.transfer_rate = ((NodeWithRates)old_node).transfer_rate;
                new_node.setLength(old_node.getLength());
                
            }
            else    
                new_node.setDefaultParams();
        }
        //Verbose.message("NWR.tT "+old_node+" -> "+new_node);
        return new_node;
    }

    /** 
     * Puts the nodes in an array suitable for depth-first traversal of the subtree rooted at this guy.
     */
    public NodeWithRates[] getDFT()
    {
        return getDFT(this);
    }
    
    /** 
     * Puts the nodes in an array suitable for depth-first traversal of the subtree rooted dft_root.
     */
    public static NodeWithRates[] getDFT(NodeWithRates dft_root)
    {
        Vector<NodeWithRates> dft_vector=new Vector<NodeWithRates>();
        collectNodes(dft_root, dft_vector, false);
        return (NodeWithRates[]) dft_vector.toArray(new NodeWithRates[0]);
    }

    
    /** 
     * @return array of leaves of the subtree rooted at this guy.
     */
    public NodeWithRates[] getLeaves()
    {
        return getLeaves(this);
    }
    
    /** 
     * @return array of leaves of the subtree rooted at this guy.
     */
    public static NodeWithRates[] getLeaves(NodeWithRates subtree_root)
    {
        Vector<NodeWithRates> leaf_vector = new Vector<NodeWithRates>();
        collectNodes(subtree_root, leaf_vector, true);
        return (NodeWithRates[]) leaf_vector.toArray(new NodeWithRates[0]);
    }
    
    /**
     * Helper recursive function for filling up the vector of nodes for depth-first traversal
     */
    private static void collectNodes(NodeWithRates N, Vector<NodeWithRates> V, boolean leaf_only)
    {
        for (int child_idx=0; child_idx<N.getNumChildren(); child_idx++)
        {
            NodeWithRates Child = (NodeWithRates)N.getChild(child_idx);
            collectNodes(Child,V,leaf_only);
        }
        if (!leaf_only || N.isLeaf())
            V.add(N);
    }
    
    
    /**
     * @return a tab-separated string for edge gain/loss probabilities in a depth-first traversal order
     */
    public String tableRates()
    {
        StringBuffer sb=new StringBuffer("# length (t)\tduplication (lambda)\tloss (mu)\ttransfer (kappa)\t// node info\n");
        tableRates(sb);
        return sb.toString();
    }
    
    private void tableRates(StringBuffer sb)
    {
        if (!isExternal())
        {
            for (int i=0; i<this.getNumChildren(); i++)
                ((NodeWithRates)getChild(i)).tableRates(sb);
        }
        if (!isRoot())  
            sb.append(getLength()+"\t"+duplication_rate+"\t"+loss_rate+"\t"+transfer_rate+"\t// "+toString()+"\n");        
    }
    
    /**
     * Sets the parameters for this edge using a line from the table produced by tableRates()
     */
    public void setEdgeParametersFromTable(String line)
    {
        String[] fields=line.split("\\s+");
        setLength(Double.parseDouble(fields[0]));
        setDuplicationRate(Double.parseDouble(fields[1]));
        setLossRate(Double.parseDouble(fields[2]));
        setTransferRate(Double.parseDouble(fields[3]));
    }
    
    /**
     * A safe way to get the name: either the true name from the newick file 
     * or a machine-constructed identifier.
     */
    public String getTaxonName()
    {
        if (getName()==null){
            return shortDesc();
        } else
            return getName();
    }
    
    
    /**
     * Distribution for xenologs
     */
    public DiscreteDistribution getTransferDistribution()
    {
        return getTransferDistribution(0.0);
    }
    
    public DiscreteDistribution getTransferDistribution(double extinction_probability)
    {
        DiscreteDistribution DD = null;

        if (transfer_rate == 0.0 || getLength()==0.0)
            DD = new PointDistribution(1.0);
        else
        {

            if (duplication_rate == 0.0)
            {
                // Math.expm1(x) = e^x-1
                // 1-e^(-y) = -(e^(-y)-1) = -Math.expm1(-y);
                double lm = (loss_rate>0.0
                        ?(transfer_rate*(-Math.expm1(-loss_rate*getLength()))/loss_rate)
                        :(transfer_rate*getLength()));
//                assert (lm>0.0);
                lm *= (1.0-extinction_probability);
                DD = (lm==0.0?new PointDistribution(1.0):new Poisson(lm));
            } else
            {
                double t = transfer_rate / duplication_rate;
                double q=0.;
                if (duplication_rate == loss_rate)
                {
                    double dr = duplication_rate*getLength();
                    q = dr/(1.+dr);
                } else 
                {
                    double b = getBeta();
                    q = duplication_rate * b;
                }
                // if (extinction_probability != 0.0) // formula OK with ext_prob=0.0
                q = q*(1.-extinction_probability)/(1.-q*extinction_probability);
                if (q<0.0) q = 0.0;
                if (q>1.0) q = 1.0;

                DD = (q==0.0)?(new PointDistribution(1.0)):(new NegativeBinomial(t,q));
            }
        }
        
        return DD;
    }
    
    public DiscreteDistribution getDuplicationDistribution()
    {
        return getDuplicationDistribution(0.0);
    }
    
    public DiscreteDistribution getDuplicationDistribution(double extinction_probability)
    {
        DiscreteDistribution DD = null;
        if (duplication_rate == 0.0)
        {
            double ml =loss_rate * getLength();
            double p1 = Math.exp(-ml);
            if (extinction_probability != 0.0)
                p1 *= 1.0-extinction_probability;
            // avoid numerical roundoff errors
            if (p1<0.0) p1=0.;
            if (p1>1.0) p1 = 1.0;
            DD = new PointDistribution(1.-p1);
        } else
        {
            double p=0.0, q=0.0;
            if (duplication_rate == loss_rate)
            {
                double dl = duplication_rate * getLength();
                p = q = dl/(1.+dl);
            } else
            {
                double b= getBeta();
                p = loss_rate * b;
                q = duplication_rate * b;
            }
            double d = 1.-q*extinction_probability;
            double pmod =  (p*(1.-extinction_probability)+(1.0-q)*extinction_probability)/d;
            double qmod = q*(1-extinction_probability)/d;
//            if (Double.isNaN(pmod) || Double.isNaN(qmod))
//            {
//                System.out.println("#*NWR.gDR NUMERROR1 p "+p+"\tq "+q+"\text "+extinction_probability
//                        +"\tb "+getBeta()
//                        +"\td "+d+"\tpmod "+pmod+"\tqmod "+qmod+"\tnode "+toString());
//                System.out.flush();
//                System.exit(999);
//            } else
            {
                p = pmod;
                q = qmod;
            }

            // avoid numerical roundoff errors
            if (p<0.0) p = 0.0;
            if (p>1.0) p = 1.0;
            if (q<0.0) q = 0.0;
            if (q>1.0) q = 1.0;

            assert (!Double.isNaN(p));
            assert (!Double.isNaN(q));

            DD = new ShiftedGeometric(p,q);
        }
        
        return DD;
    }
    /**
     * @return a vector of h_t(i): probabilities for i=0..n copies that 
     * trace back to a horizontal tranfer event on this branch. 
     */
    //public double[] getBasicProbabilityTransfer(int n)
    //{
    //    double[] retval = new double[n+1];
    //    double theta = transfer_rate/gain_rate;
    //    double beta = getBeta();
    //    double L =gain_rate*beta;
    //    retval[0]=Math.pow(1.0-L, theta);
    //    double sum=retval[0];
    //    for (int j=1; j<=n; j++){
    //        retval[j] = retval[j-1]*((theta+j-1.0)/((double)j))*L;
    //        sum += retval[j];
    //        Verbose.message("NWR.gBPT "+this.getTaxonName()+"\tj "+j+"\t"+retval[j]+"\t"+sum);
    //    }
    //    
    //    return retval;
    //}
            
     /**
     * @return a vector of g_t(i): probabilities for i=0..n copies that 
     * trace back to a single ancestral gene on this branch. 
     */
   //public double[] getBasicProbabilityDuplication(int n){
   //     double[] retval = new double[n+1];
   //     double beta = getBeta();
   //     double mb = loss_rate*beta;
   //     double lb = gain_rate*beta;
   //     retval[0]=mb;
   //     if (n>0){
   //         retval[1]=(1.0-mb)*(1.0-lb);
   //         for (int j=2; j<=n; j++)
   //             retval[j] = retval[j-1]*lb;
   //     }
   //     return retval;
   // }
    
    /**
     * Computes <var>beta</var>(\var>t</var>) for the
     * distribution formulas.
     *
     * @return beta(t), that is, (1-e^{-(\mu-\lambda)t})/(\mu-\lambda e^{-(\mu-\lambda)t}).
     */
    public double getBeta()
    {
        double d=loss_rate-duplication_rate;
        double minus_dl = -d*getLength();
        double E = Math.exp(minus_dl);
        double y = loss_rate-duplication_rate*E;
        if (minus_dl<1.0) // this is the expected behavior: loss rate should be larger than duplication rate
        {
            double x = -Math.expm1(minus_dl); //1.0-E;
            return x/y;
        } else // duplication rate is larger than loss_rate + 1.0!
        {
            // need to be careful with E here:
            // minus_dl is a large positive number,
            // E is a large positive number, possibly infinity
            // y is a large negative number, possible negative inifnity
            // d is negative

            // x/y may become NaN (infinity / infinity)

            // (1-exp(-dl)) / (mu-lambda exp^(-dl))
            // = 1/lambda * ((mu-lambda*E) + lambda-mu )/ (mu-lambda * E)
            // = 1/lambda * (1+(lambda-mu)/(mu-lambda*E))
            return (1.0 -d / y)/duplication_rate;
        }

    }
    
    protected String paramString()
    {
        StringBuffer sb=new StringBuffer(super.paramString());
        sb.append("; dup ");
        sb.append(duplication_rate);
        sb.append(", loss ");
        sb.append(loss_rate);
        sb.append(", transfer ");
        sb.append(transfer_rate);
        //sb.append(", p11 ");
        //sb.append(getProbability(true,true));
        if (SHOW_SIZE_DISTRIBUTIONS_IN_PARAMSTRING)
        {
            int m = 3;
            sb.append(", g[]={");
            double [] dup = getDuplicationDistribution().getDistribution(m);
            for (int i=0; i<=m; i++)
            {
                if (i!=0)
                    sb.append(',');
                sb.append(dup[i]);
            }
            sb.append('}');

            sb.append(", h[]={");
            double [] h = getTransferDistribution().getDistribution(3);
            for (int i=0; i<=m; i++)
            {
                if (i!=0)
                    sb.append(',');
                sb.append(h[i]);
            }
            sb.append('}');
        }
        return sb.toString();
    }

  
}