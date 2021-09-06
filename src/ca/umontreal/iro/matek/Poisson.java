/*
 * Poisson.java
 *
 * Created on March 20, 2008, 12:19 PM
 */

package ca.umontreal.iro.matek;

/**
 * Negative binomial distribution. 
 * Specifically, the point mass function [<var>p</var><sub>i</sub>]
 * is defined as follows.
 * 
 * <var>p</var><sub>k</sub>=e<sup>-<var>r</var></sup> <var>r</var><sup><var>k</var></sup>/<var>k</var>!.
 * 
 * @author csuros
 */
public class Poisson implements DiscreteDistribution
{

    /**
     * Poisson distribution.
     * @param r distribution parameter
     */
    public Poisson(double r) 
    {
        set(r);
    }
    
    private double r;
    
    /**
     * Computes point mass function [<var>p</var><sub>i</sub>] for a Poisson distribution.
     * Reminder: <var>p</var><sub>k</sub>=e<sup>-<var>r</var></sup> <var>r</var><sup><var>k</var></sup>/<var>k</var>!.
     */
    public double[] getDistribution(int n)
    {
        //System.out.println("#**Poisson  r="+r);
        double[] d = new double[n+1];
        d[0] = Math.exp(-r);
        //System.out.println("#**Poisson[0]\t"+d[0]);
        for (int j=1; j<=n; j++)
        {
            double f = r/j;
            d[j] = d[j-1]*f;
            //System.out.println("#**Poisson["+j+"]\t"+d[j]);
        }
        return d;
    }
    
    /**
     * @return {r}
     */
    public double[] getParameters()
    {
        double[] parm = new double[1];
        parm[0]=r;
        return parm;
    }
    
    public void setParameter(int parameter_idx, double r)
    {
        if (parameter_idx == 0)
            set(r);
        else
            throw new IllegalArgumentException("First argument to setParameter must be 0");
    }
    
    /**
     * Sets the single parameter of the Poisson distribution
     */
    public void set(double r)
    {
        this.r= r;
    }
    
    public final int getNumParameters(){return 1;}
    
}
