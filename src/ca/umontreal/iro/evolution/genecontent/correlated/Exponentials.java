package ca.umontreal.iro.evolution.genecontent.correlated;

import java.util.ArrayList;

/**
 * Helper class for representing functions that are 
 * linear combinations of exponentials written 
 * in the form <var>f</var>(<var>t</var>)=<var>a</var><sub>1</sub> e<sup>-<var>u</var><sub>1</sub> <var>t</var></sup>
 * + <var>a</var><sub>2</sub> e<sup>-<var>u</var><sub>2</sub> <var>t</var></sup>
 * + ... + <var>a</var><sub><var>n</var></sub> e<sup>-<var>u</var><sub><var>n</var></sub> <var>t</var></sup>.<br>
 * It is assumed that <var>u</var><sub><var>i</var></sub> are nonnegative. <br>
 * The class gives some limited support for terms of the form 
 * <var>a</var> <var>t</var> e<sup>-<var>u</var> <var>t</var></sup> that can 
 * arise as the result of combine() and are properly handled by eval(). 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s and Louis Philippe B. Bouvrette
 */
public class Exponentials
{
    /**
     * Creates a new function.
     * 
     * @param coefficients array of coefficients <var>a</var>
     * @param exponents array of exponents <var>u</var> (all non-negative)
     */
    public Exponentials(double[] coefficients, double[] exponents)
    {
	if (coefficients.length != exponents.length)
	    throw new IllegalArgumentException("Coefficients and exponents must have the same length");
	this.coeff = new double[coefficients.length];
	System.arraycopy(coefficients, 0, coeff, 0, coefficients.length);
	this.expon = new double[exponents.length];
	System.arraycopy(exponents, 0, expon, 0, exponents.length);
    }
    
    /**
     * Creates a new constant function (no dependence on <var>t</var>)
     * 
     * @param constant the constant value 
     */
    public Exponentials(double constant)
    {
	this.coeff = new double[1];
	this.expon = new double[1];
	coeff[0] = constant;
	expon[0] = 0.0;
    }
    
    /**
     * Instanciation without parameters: used by internal methods
     * 
     * @param num_terms how many terms there will be
     */
    private Exponentials(int num_terms)
    {
	this.coeff = new double[num_terms];
	this.expon = new double[num_terms];
    }
    
    private double[] coeff;
    private double[] expon;
    
    /**
     * Evaluates the function for a given <var>t</var>
     * 
     * @param t where the function value is computed
     * @return the function value <var>f</var>(<var>t</var>)
     */
    public double eval(double t)
    {
	double result = 0.0;
	for (int i=0; i<coeff.length; i++)
	    if (expon[i] != 0.0)
                {
                    if (expon[i]>0.0)
                        result += coeff[i]*Math.exp(-expon[i]*t);
                    else // negative exponent is used for t*exp(-u t)
                        result += coeff[i]*t*Math.exp(expon[i]*t);
                } else
		result += coeff[i];
	return result;
    }
    
    /**
     * Computes a sort of a convolution, defined as 
     * the integral of mu * exp(-mu s)*f(s)*g(t-s) ds where 
     * s goes from 0 to t.    
     * 
     * @param mu the mu parameter of the convolution  
     * @param g the second function of the convolution
     * @return the convolution result
     */
    public Exponentials convolution(double mu, Exponentials g)
    {
	ArrayList<Double> result_coeff = new ArrayList<Double>();
	ArrayList<Double> result_expon = new ArrayList<Double>();
        
	for (int i=0; i<this.coeff.length; i++)
            {
                for (int j=0; j<g.coeff.length; j++)
		    {
			if (mu + expon[i] != g.expon[j])
			    {
				double c = mu*coeff[i]*g.coeff[j]/(mu+expon[i]-g.expon[j]);
				result_coeff.add(c);
				result_expon.add(g.expon[j]);
				result_coeff.add(-c);
				result_expon.add(mu+expon[i]);
			    } else
			    {
				double c = mu*coeff[i]*g.coeff[j];
				result_coeff.add(c);
				result_expon.add(-g.expon[j]);
			    }
		    }
            }
	
	Exponentials result = new Exponentials(result_coeff.size());
	for (int i=0; i<result_coeff.size(); i++)
            {
                result.coeff[i] = result_coeff.get(i);
                result.expon[i] = result_expon.get(i);
            }
	return result;
    }
    
    /**
     * Adds two Exponentials
     * 
     * @param x the other term in the addition
     * @return result (<code>this</code>+<var>x</var>)
     */
    public Exponentials add(Exponentials x)
    {
	Exponentials result = new Exponentials(coeff.length+x.coeff.length);
	for (int i=0; i<coeff.length; i++)
            {
                result.coeff[i] = coeff[i];
                result.expon[i] = expon[i];
            }
	for (int j=0; j<x.coeff.length; j++)
            {
                int k = coeff.length+j;
                result.coeff[k] = x.coeff[j];
                result.expon[k] = x.expon[j];
            }
	return result;
    }
    
    
    /**
     * Subtracts two Exponentials
     * 
     * @param x the term to be subtracted
     * @return result (<code>this</code>-<var>x</var>)
     */
    public Exponentials subtract(Exponentials x)
    {
	Exponentials result = new Exponentials(coeff.length+x.coeff.length);
	for (int i=0; i<coeff.length; i++)
            {
                result.coeff[i] = coeff[i];
                result.expon[i] = expon[i];
            }
	for (int j=0; j<x.coeff.length; j++)
            {
                int k = coeff.length+j;
                result.coeff[k] = -x.coeff[j];
                result.expon[k] = x.expon[j];
            }
	return result;
    }
        
    /**
     * Multiplies two Exponentials
     * 
     * @param x the other factor in the product
     * @return result (<code>this</code>*<var>x</var>)
     */
    public Exponentials multiply(Exponentials x)
    {
	Exponentials result = new Exponentials(coeff.length*x.coeff.length);
	int k=0;
	for (int i=0; i<coeff.length; i++)
	    for (int j=0; j<x.coeff.length; j++)
                {
                    result.coeff[k] = coeff[i]*x.coeff[j];
                    result.expon[k] = expon[i]+x.expon[j];
                    k++;
                }
	return result;
    }
    
    public String toString()
    {
        StringBuffer sb = new StringBuffer("Exp[");
        for (int i=0; i<coeff.length; i++)
        {
            if (i!=0)
                sb.append("+ ");
            sb.append(Double.toString(coeff[i]));
            sb.append("*e^");
            sb.append(Double.toString(expon[i]));
            sb.append("t");
        }
        sb.append("]");
        return sb.toString();
    }
}

