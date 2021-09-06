/*
 * PointDistribution.java
 *
 * Created on March 20, 2008, 12:36 PM
 *
 */

package ca.umontreal.iro.matek;

/**
 * Point distribution: p_0=p; p_1=1-p.
 *
 * @author csuros
 */
public class PointDistribution implements DiscreteDistribution
{
    
    /** Creates a new instance of PointDistribution.
     * @param p probability of X=0
     */
    public PointDistribution(double p) 
    {
        set(p);
    }
    
    private double p;

    /**
     * Probability of X=k for k=0,1,..., n.
     * @param n last index for probabilities
     * @return array of probabilities
     */
    @Override
    public double[] getDistribution(int n)
    {
        double[] d= new double[n+1];
        d[0] = p;
        if (n>=1) d[1]=1.-p;
        return d;
    }

    /**
     * Number of parameters for this distribution.
     *
     * @return 1
     */
    @Override
    public final int getNumParameters(){return 1;}
    
    /**
     * Returns the probability of x==0
     * 
     * @return p
     */
    @Override
    public double[] getParameters()
    {
        double[] parm = new double[1];
        parm[0]=p;
        return parm;
    }
    
    /**
     * Sets the distribution's parameter p (i.e., probability of X==0)
     * @param p probability of X=0
     */
    public void set(double p)
    {
        this.p = p;
    }
    
    /**
     * Sets the single parameter of this distribution
     *
     * @param must_be0 must be 0
     * @param value probability of X==0, a value between 0 and 1 (inclusive)
     */
    @Override
    public void setParameter(int must_be0, double value)
    {
        if (must_be0==0)
            set(value);
        else
            throw new IllegalArgumentException("First argument to setParameter must be 0");
    }
    
    
}
