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

import java.util.function.DoubleFunction;

/**
 * Ziheng Yang's discrete gamma approximation for rates-across-sites variation
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 * @since March 21, 2008, 9:23 AM
 */
public final class DiscreteGamma 
{
    /**
     * Initializes the regular gamma distribution with mean 1
     * The parameter should be less than about 100.0 to avoid numerical errors.
     * 
     * @param alpha shape parameter
     */
    public DiscreteGamma(double alpha)
    {
        setAlpha(alpha);
    }
    
    public static final double RECOMMENDED_MAXIMUM_ALPHA = 100.0;
    
    /**
     * Initializes with alpha=1.0
     */
    public DiscreteGamma()
    {
        this(1.0);
    }
    
    private double alpha;
    
    public double getAlpha()
    {
        return alpha;
    }
    
    public void setAlpha(double alpha)
    {
        this.alpha = alpha;
    }
    
    /**
     * Computes the partition of the Gamma distribution into 
     * n categories. The partition is defined by the percentage points for
     * 1/n, 2/n, ..., 1-1/n.
     *
     * @param n number of partitions: must be at least 1
     * @return array of n-1 elements corresponding to the percentage points i/n i=1,...,n-1
     */
    public double[] getPartition(int n)
    {
        double[] x = new double[n-1];
        // find 1/n percentage point
        
        // define brackets
        double y1 =1.0/n;
        PercentagePoint PP1 = new PercentagePoint(y1);
        // lower bracket endpoint
        double a1=0.1;
        while (PP1.apply(a1)>0.)
        {
            a1 *= 0.5;
        }
        double b1 = 2.0*a1;
        while (PP1.apply(b1)<0.)
            b1 *= 2.0;
        //System.out.println("#**DG.gP /1 "+a1+"\t"+b1);
        x[0]=FunctionMinimization.zbrent(PP1,a1,b1,(b1-a1)*1e-6);
        for (int j=1; j<x.length; j++)
        {
            double yj = (1.0+j)/n;
            PercentagePoint PPj = new PercentagePoint(yj);
            double b = x[j-1]*2.0;
            while (PPj.apply(b)<0.)
                b *= 2.0;
            x[j] = FunctionMinimization.zbrent(PPj,x[j-1],b,(b-x[j-1])*1e-6);
        }
        return x;
    }
        
    /**
     * Takes the partition of the Gamma distribution into 
     * n categories. The partition is defined by the percentage points for
     * 1/n, 2/n, ..., 1-1/n. Computes the mean value within each interval 
     * [i/n, (i+1)/n] for i=0..n-1
     *
     * @param n number of partitions: supposed to be at least 1
     * @return array of n elements corresponding to the means within the partitions
     */
    public double[] getPartitionMeans(int n)
    {
        return getPartitionMeans(n, null);
    }
    
    /**
     * Takes the partition of the Gamma distribution into 
     * n categories. The partition is defined by the percentage points for
     * 1/n, 2/n, ..., 1-1/n. Computes the mean value within each interval 
     * [i/n, (i+1)/n] for i=0..n-1
     *
     * @param n number of partitions: supposed to be at least 1
     * @param partition if not null, then this must be an array of length at least n-1, where the partition boundaries are saved (result of getPartition)
     * @return array of n elements corresponding to the means within the partitions
     */
    public double[] getPartitionMeans(int n, double[] partition)
    {
    
        double[] x = getPartition(n);
        double[] m = new double[n];
        
        double prev_z = 0.0;
        for (int j=0;j<n; j++)
        {
            double current_z = 1.0;
            if (j!=n-1)
            {
                current_z = Functions.gammp(alpha+1.0,alpha*x[j]);
            }
            m[j] = n*(current_z-prev_z);
            prev_z = current_z;
        }
        
        if (partition != null)
            System.arraycopy(x,0,partition,0,n-1);
        return m;
    }
    
    
    
    /**
     * Computes the point at which the cdf F(x)=u
     */
    public double inverseGammaCDF(double u)
    {
        PercentagePoint PP = new PercentagePoint(u);
        // lower bracket endpoint
        double a1=0.1;
        while (PP.apply(a1)>0.)
        {
            a1 *= 0.5;
        }
        double b1 = 2.0*a1;
        while (PP.apply(b1)<0.)
            b1 *= 2.0;
        //System.out.println("#**DG.gP /1 "+a1+"\t"+b1);
        double x =FunctionMinimization.zbrent(PP,a1,b1,(b1-a1)*1e-6);
        return x;
    }
    
    private class PercentagePoint implements DoubleFunction<Double>
    {
        private PercentagePoint(double cutoff)
        {
            this.cutoff = cutoff;
        }
        
        private double cutoff;
        
        @Override
        public Double apply(double x)
        {
            double P = Functions.gammp(alpha,x*alpha);
            //System.out.println("#**DG.PP alpha "+alpha+"\tx "+x+"\tP "+P);
            
            return P-cutoff;
        }
    }
    
    private void mainmain(String[] args)
    {
        if (args.length != 2)
        {
            throw new IllegalArgumentException("Call as "+getClass().getName()+" alpha n // alpha: parameter for Gamma distribution; n: number of discrete categories");
        }
        double a = Double.parseDouble(args[0]);
        int n = Integer.parseInt(args[1]);
        
        setAlpha(a);
        double[] x = getPartition(n);
        for (int i=0; i<x.length; i++)
            System.out.println("# Threshold x["+i+"]\t"+x[i]);
        double[] m = getPartitionMeans(n);
        for (int i=0; i<m.length; i++)
            System.out.println("# Mean value in category "+i+"\t"+m[i]);
        
        System.out.println("#Plotting: y(multiplier)\tx-mid");
        System.out.println("0.0\t0.0");
        for (int j=0; j<m.length; j++)
        {
            double x0 = (j==0?0.0:x[j-1]);
            double x1 = (j==m.length-1?x0:x[j]);
            double xm = (x0+x1)/2.0;
            System.out.println(m[j]+"\t"+xm);
        }        
    }
    
    public static void main(String[] args)
    {
        DiscreteGamma O = new DiscreteGamma();
        O.mainmain(args);

    }
}
