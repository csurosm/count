/*
 * ShiftedGeometric.java
 *
 * Created on March 20, 2008, 10:59 AM
 *
 */

package ca.umontreal.iro.matek;

/**
 * 
 * Shifted geometric distribution. 
 * Specifically, the point mass function [<var>p</var><sub>i</sub>]
 * is defined as follows.
 * 
 * <var>p</var><sub>0</sub>=<var>p</var>; 
 * <var>p</var><sub>1</sub>=(1-<var>p</var>)(1-<var>q</var>); 
 * <var>p</var><sub>k</sub>=(1-<var>p</var>)(1-<var>q</var>)<var>q</var><sup><var>k</var>-1</sup>.
 * 
 * @author csuros
 */
public class ShiftedGeometric implements DiscreteDistribution
{

    /*
     *
     * @param p distribution parameter: probability of 0
     * @param q distribution parameter: failure probability
     */
    public ShiftedGeometric(double p, double q) 
    {
        this.p = p;
        this.q = q;
    }
    
    private double p;
    private double q;
    
    /**
     * Computes point mass function [<var>p</var><sub>i</sub>] for a shifted geometric distribution.
     * Reminder: <var>p</var><sub>k</sub>=(1-<var>p</var>)(1-<var>q</var>)<var>q</var><sup><var>k</var>-1</sup>.
     *
     * @param n limit on how far the probabilities are computed
     * @return [<var>p</var><sub>0</sub> <var>p</var><sub>1</sub> ... <var>p</var><sub>n</sub>]
     */
    public double[] getDistribution(int n)
    {
        double[] d = new double[n+1];
        d[0] = p;
        //System.out.println("#*SG.gD n "+0+"\t"+d[0]+"\t// p "+p+"\tq "+q);
        if (n>=1)
        {
            d[1]=(1-p)*(1-q);
            //System.out.println("#*SG.gD n "+1+"\t"+d[1]+"\t// p "+p+"\tq "+q);
        }
        if (q!=0.0)
        {
            for (int i=2; i<=n; i++)
            {
                d[i] = d[i-1]*q;
                //System.out.println("#*SG.gD n "+n+"\t"+d[n]+"\t// p "+p+"\tq "+q);
            }
        }
        return d;
    }    
    
    /**
     * Returns p and q for this distribution
     *
     * @return {p,q}
     */
    public double[] getParameters()
    {
        double[] parm = new double[2];
        parm[0] = p;
        parm[1] = q;
        return parm;
    }
    
    public void setParameter(int param_idx, double value)
    {
        if (param_idx==0)
            p = value;
        else if (param_idx==1)
            q = value;
        else 
            throw new IllegalArgumentException("First argument to setParameter must be 0 or 1.");
    }
    
    public final int getNumParameters(){return 2;}
    
    
}
