/*
 * ApproximateDerivation.java
 *
 * Created on November 3, 2006, 9:35 PM
 */

package ca.umontreal.iro.matek;

/**
 *
 * @author  csuros
 */
public class ApproximateDerivation implements DerivableMultiParameterFunction
{
    
    /** Creates a new instance of ApproximateDerivation */
    public ApproximateDerivation(MultiParameterFunction func) 
    {
    }
    
    private MultiParameterFunction func;
    
    private static final double EPS=1e-7; // cubic root of machine precision
    
    public double[] dfunc(double[] x) 
    {
        double[] retval = new double[x.length]; 
        for (int i=0; i<x.length; i++)
        {
            double[] xh=new double[x.length];
            double h = x[i]*EPS;
            for (int j=0; j<x.length; j++)
                xh[j]=x[j];
            xh[i]=x[i]+h;
            double y2=func.eval(xh);
            xh[i]=x[i]-h;
            double y1=func.eval(xh);
            retval[i]=(y2-y1)/(2.0*h);
        }
            
        return retval;
    }
    
    public double eval(double[] param) 
    {
        return func.eval(param);
    }
    
}
