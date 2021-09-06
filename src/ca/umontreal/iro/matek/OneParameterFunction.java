/*
 * OneParameterFunction.java
 *
 * Created on October 26, 2004, 8:22 PM
 */

package ca.umontreal.iro.matek;

/**
 * A common interface to functions with one parameter: used by FunctionMinimization class.
 *
 * @author  csuros
 */
public interface OneParameterFunction {
    /**
     * Calculates the value of a function <var>f</var> at a point.
     * @param param <var>x</var> coordinate
     * @return function value <var>f</var>(<var>x</var>)
     */
    public double eval(double param);
    
}
