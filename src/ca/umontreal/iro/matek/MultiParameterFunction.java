/*
 * MultiparameterFunction.java
 *
 * Created on October 28, 2004, 12:40 PM
 */

package ca.umontreal.iro.matek;

/**
 *
 * @author  csuros
 */
public interface MultiParameterFunction {
    
    /**
     * Function value at a given point
     * 
     * @param param vector of point coordinates
     * @return function value at the point
     */
    public double eval(double[] param);
    
}
