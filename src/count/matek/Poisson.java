/*
 * Copyright 2021 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package count.matek;

/**
 * Poisson distribution. 
 * Specifically, the point mass function [<var>p</var><sub>i</sub>]
 * is defined as follows.
 * 
 * <var>p</var><sub>k</sub>=e<sup>-<var>r</var></sup> <var>r</var><sup><var>k</var></sup>/<var>k</var>!.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 * 
 * @since March 20, 2008, 12:19 PM
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
    @Override
    public double[] getPointMassFunction(int n)
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
    @Override
    public double[] getParameters()
    {
        double[] parm = new double[1];
        parm[0]=r;
        return parm;
    }
    
    @Override
    public void setParameter(int parameter_idx, double r)
    {
        if (parameter_idx == 0)
            set(r);
        else
            throw new IllegalArgumentException("First argument to setParameter must be 0");
    }
    
    /**
     * Sets the single parameter of the Poisson distribution
     * @param r Poisson parameter (= mean)
     */
    public final void set(double r)
    {
        this.r= r;
    }
    
    @Override
    public final int getNumParameters(){return 1;}
    
    @Override
    public String toString()
    {
        return getClass().getSimpleName()+"[r "+r+"]";
    }
    
}
