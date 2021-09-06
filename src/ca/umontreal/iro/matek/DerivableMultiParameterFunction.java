/*
 * DerivableMultiParameterFunction.java
 *
 * Created on November 3, 2006, 1:28 PM
 */

package ca.umontreal.iro.matek;

/**
 * A function for which the partial derivatives are computed. 
 * An interface used by function minimization algorithms using gradient.
 *
 * @author  csuros
 */
public interface DerivableMultiParameterFunction extends MultiParameterFunction 
{
    /**
     * Partial derivatives at a point.
     * @param x vector of point coordinates
     * @return vector of partial derivatives at the point
     */
    public double[] dfunc(double[] x);
}
