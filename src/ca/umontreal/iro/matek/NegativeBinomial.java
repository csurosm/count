/*
 * NegativeBinomial.java
 *
 * Created on March 20, 2008, 12:10 PM
 *
 */

package ca.umontreal.iro.matek;

/**
 * Negative binomial distribution. 
 * Specifically, the point mass function [<var>p</var><sub>i</sub>]
 * is defined as follows.
 * 
 * <var>p</var><sub>k</sub>=Binom(<var>t</var>+<var>n</var>-1,<var>n</var>)(1-<var>q</var>)<sup>t</sup><var>q</var><sup><var>k</var></sup>(1-<var>p</var>)(1-<var>q</var>)<var>q</var><sup><var>k</var>-1</sup>.
 *
 * @author csuros
 */
public class NegativeBinomial implements DiscreteDistribution
{
    
    /**
     * @param t distribution parameter.
     * @param q distribution parameter 
     */
    public NegativeBinomial(double t, double q) 
    {
        this.t = t;
        this.q = q;
    }
    
    private double t;
    private double q;
    

    /**
     * Computes point mass function [<var>p</var><sub>i</sub>] for a negative binomial distribution.
     * Reminder: <var>p</var><sub>k</sub>=Binom(<var>t</var>+<var>n</var>-1,<var>n</var>)(1-<var>q</var>)<sup>t</sup><var>q</var><sup><var>k</var></sup>(1-<var>p</var>)(1-<var>q</var>)<var>q</var><sup><var>k</var>-1</sup>.
     */
    @Override
    public double[] getDistribution(int n)
    {
        double[] d = new double[n+1];
        if (t>2e9) // largest integer is about 2.147e9
        {
            // use logarithms in the calculation
            
            // (1-q)^t = e^{-tq) when t is large
            d[0]=-q*t;
            double logq = Math.log(q);
            for (int j=1; j<=n; j++)
            {
                double f = Math.log((t+j-1.0)/((double)j));
                d[j] = d[j-1]+f+logq; 
                //System.out.println("#*NB.gD logd["+j+"]= "+d[j]+"\t"+d[j-1]+"\tf "+f+"\tq "+q+"\tt "+t+"\t"+(1.-q));
            }
            for (int j=0; j<d.length; j++)
            {
                d[j] = Math.exp(d[j]);
            }
        }
        else
        {
            double q1 = 1.-q;
            d[0] = Math.pow(q1,t);
            for (int j=1; j<=n; j++)
            {
                double f = (t+j-1.0)/((double)j);
                d[j] = d[j-1]*f*q; 
                //System.out.println("#*NB.gD logd["+j+"]= "+d[j]+"\t"+d[j-1]+"\tf "+f+"\tq "+q+"\tt "+t+"\t"+(1.-q));
            }
        }
        
        //for (int j=0; j<d.length; j++)
        //{
        //    System.out.println("#*NB.gD d["+j+"]= "+d[j]);
        //}
        
        return d;
    }
    
    /**
     * @return {t,q}
     */
    @Override
    public double[] getParameters()
    {
        double[] parms = new double[2];
        parms[0]=t;
        parms[1]=q;
        
        return parms;
    }
    
    @Override
    public void setParameter(int param_idx, double value)
    {
        if (param_idx==0)
            this.t = value;
        else if (param_idx==1)
            this.q = value;
        else
            throw new IllegalArgumentException("First argument to setParameter must be 0 or 1");
    }
    
    @Override
    public final int getNumParameters(){return 2;}

    @Override
    public String toString()
    {
        return getClass().getSimpleName()+"[t "+t+"; q "+q+"]";
    }
}
