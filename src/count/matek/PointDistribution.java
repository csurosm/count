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
 * Point (Bernoulli) distribution: <var>p</var>(0)=<var>p</var>; <var>p</var>(1)=1-p.
 *
 * @since March 20, 2008, 12:36 PM
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class PointDistribution implements DiscreteDistribution
{
    /**
     * X=1 with probability 1.
     */
    public PointDistribution(){ this(0.0);}
    
    /** 
     * Creates a new instance of PointDistribution.
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
        // d[2..] are all 0.0 if those were requested
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
    public final void set(double p)
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
