/*
 * DiscreteDistribution.java
 *
 * Created on March 20, 2008, 10:41 AM
 *
 */

package ca.umontreal.iro.matek;

/**
 * Common interface for working with probabilistic distributions.
 *
 * @author csuros
 */
public interface DiscreteDistribution 
{
    /**
     * Returns the first few elements of the probability mass function
     *
     * @param n limit on how far the probabilities are computed
     * @return array of <var>n</var>+1 elements: [<var>p</var><sub>0</sub> <var>p</var><sub>1</sub> ... <var>p</var><sub>n</sub>]
     */
    
    public double[] getDistribution(int n);
    
    /**
     * Returns the distribution's parameters: depends on the implementation what they mean 
     */
    public double[] getParameters();
    
    
    /**
     *  Sets one of the distribution parameters
     *
     * @param parameter_idx index of the parameter (parameters are indexed the same order as with getParameters())
     * @param value New value of the parameter
     */
    public void setParameter(int parameter_idx, double value);
    

    /**
     * Returns the number of parameters for defining the distribution: same as the length for the array returned by getParameters()
     */
    public int getNumParameters();
    
}
