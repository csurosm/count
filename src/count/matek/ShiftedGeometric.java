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
 * 
 * Shifted geometric distribution. 
 * Specifically, the point mass function [<var>p</var><sub>i</sub>]
 * is defined as follows.
 * 
 * <var>p</var><sub>0</sub>=<var>p</var>; 
 * <var>p</var><sub>1</sub>=(1-<var>p</var>)(1-<var>q</var>); 
 * <var>p</var><sub>k</sub>=(1-<var>p</var>)(1-<var>q</var>)<var>q</var><sup><var>k</var>-1</sup>.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 * 
 * @since March 20, 2008, 10:59 AM
 */
public class ShiftedGeometric implements DiscreteDistribution
{
    /**
     *
     * @param p distribution parameter: probability of 0
     * @param q distribution parameter: failure probability
     */
    public ShiftedGeometric(double p, double q) 
    {
        this.p = p;
        this.q = q;
    }
    
    private double p;
    private double q;
    
    /**
     * Computes point mass function [<var>p</var><sub>i</sub>] for a shifted geometric distribution.
     * Reminder: <var>p</var><sub>k</sub>=(1-<var>p</var>)(1-<var>q</var>)<var>q</var><sup><var>k</var>-1</sup>.
     *
     * @param n limit on how far the probabilities are computed
     * @return [<var>p</var><sub>0</sub> <var>p</var><sub>1</sub> ... <var>p</var><sub>n</sub>]
     */
    @Override
    public double[] getDistribution(int n)
    {
        double[] d = new double[n+1];
        d[0] = p;
        //System.out.println("#*SG.gD n "+0+"\t"+d[0]+"\t// p "+p+"\tq "+q);
        if (n>=1)
        {
            d[1]=(1-p)*(1-q);
            //System.out.println("#*SG.gD n "+1+"\t"+d[1]+"\t// p "+p+"\tq "+q);
        }
        if (q!=0.0)
        {
            for (int i=2; i<=n; i++)
            {
                d[i] = d[i-1]*q;
                //System.out.println("#*SG.gD n "+n+"\t"+d[n]+"\t// p "+p+"\tq "+q);
            }
        }
        return d;
    }    
    
    /**
     * Returns p and q for this distribution
     *
     * @return {p,q}
     */
    @Override
    public double[] getParameters()
    {
        double[] parm = new double[2];
        parm[0] = p;
        parm[1] = q;
        return parm;
    }
    
    public final static int LOSS_IDX=0;
    public final static int DUPLICATION_IDX=1;
    
    @Override
    public void setParameter(int param_idx, double value)
    {
        if (param_idx==0)
            p = value;
        else if (param_idx==1)
            q = value;
        else 
            throw new IllegalArgumentException("First argument to setParameter must be 0 or 1.");
    }
    
    @Override
    public final int getNumParameters(){return 2;}
    
    @Override
    public String toString()
    {
        return getClass().getSimpleName()+"[p "+p+"; q "+q+"]";
    }
    
}