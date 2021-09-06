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
 * Arithmetics on log-scale. 
 * 
 * @author csuros
 *
 */
public class Logarithms 
{

	/**
	 * Adds two values. 
	 * @param x
	 * @param y
	 * @return ln(e<sup><var>x</var></sup>+e<sup><var>y</var></sup>)
	 */
	public static double add(double x, double y)
	{
//		++cnt_op_add;
        // ln (e^x + e^y) = x + ln (1+e^{y-x})
        if (y>=x)
        {
    		if (Double.isInfinite(x))
    			return y; // covers 0+0 
    		else
    			return y + Math.log1p(Math.exp(x-y));
        } else
        {
			return x + Math.log1p(Math.exp(y-x));
        }		
	}
	
//	private static int cnt_op_add = 0;

	/**
	 * Sums multiple log-scale values. 
	 * 
	 * @param x array ; unchanged
	 * @param n prefix over which summation is carried out 
	 * @return ln(exp(<var>x</var>[0]+...+<var>x</var>[n-1]) (negative infinity if n=0)
	 */
	public static double sum(double[] x, int n)
	{
		if (n==0) return Double.NEGATIVE_INFINITY;
		
		double sum = x[0];
		for (int i=1; i<n; i++)
			sum = add(sum, x[i]);
		return sum;
	}
	
	private static final double LOG_EPS = Math.log(Functions.EPS);

	/**
	 * Sums multiple log-scale values 
	 * @param x input array; values will be rearranged 
	 * @param n
	 * @return ln(exp(<var>x</var>[0]+...+<var>x</var>[n-1]) (negative infinity if n=0)
	 */
	public static double fastSum(double[] x, int n)
	{
		if (n==0) return Double.NEGATIVE_INFINITY;
		// heapify
        for (int i=n/2-1; i>=0; i--)
        {
            sink(x, x[i], i, n);
        }
        double sum=Double.NEGATIVE_INFINITY;
        int m=n;
        while (m>0)
        {
        	--m;
        	double v = x[0];
        	
        	if (sum == Double.NEGATIVE_INFINITY)
        		sum = v;
        	else 
        	{
            	double d = v-sum;

            	double a = Math.log1p(Math.exp(d));

            	if (d<LOG_EPS)
            		break;
            	
	        	sum += a;
        	}
        	double y = x[m];
        	x[m]=v; // saving the element, but we could just forget about it ... 
        	sink(x, y, 0, m);
        }
        return sum;
        
	}
	
    /**
     * Sink with max-heap.
     * 
     * @param H array with its prefix [0..n-1] in maxheap order
     * @param x value to be placed at H[i] or below
     * @param i index for placement
     * @param n number of elements on the heap
     */
    private static void sink(double[] H, double x, int i, int n)
    {
        assert (n<=H.length); 
        
        int c=2*i+1; // heap arithmetics: 0->(1,2) 1->(3,4), 2->(5,6), ...
        while (c<n)
        {
            double vc = H[c]; 
            int c2 = c+1;
            if (c2<n)
            {
                double v2 = H[c2]; if (v2 > vc){ c=c2; vc = v2;}
            }
            if (vc<=x) break;
            H[i] = vc; i=c; c= 2*i+1; 
        }
        H[i] = x;
    }
    
    
    public static void main(String[] args)
    {
    	int N = 20;
    	if (args.length>0)
    		N = Integer.parseInt(args[0]);
    	
    	double[] x = new double[N];
    	double[] y = new double[N];
    	double[] z = new double[N];
    	
    	java.util.Random RND = new java.util.Random(2021);
    	
    	double sum_x = 0.0;
    	double sum_z = 0.0;
    	
    	for (int i=0; i<N; i++)
    	{
    		double r = RND.nextDouble();
    		x[i] = 1.0/r;
    		y[i] = Math.log(x[i]);
    		z[i] = Math.exp(x[i]);
    		sum_x += x[i];
    		sum_z += z[i];
    		System.out.println("# "+i+"\t"+x[i]+"\t"+y[i]+"\t"+z[i]);
    	}
    	
    	double log_sumx1 = Logarithms.sum(y, N); 
//    	int cnt_opx1 = cnt_op_add; cnt_op_add=0;
    	double log_sumx2 = Logarithms.fastSum(y, N); 
//    	int cnt_opx2 = cnt_op_add; cnt_op_add=0;
    	double log_sumz1 = Logarithms.sum(x, N); 
//    	int cnt_opz1 = cnt_op_add; cnt_op_add = 0;
    	double log_sumz2 = Logarithms.fastSum(x, N); 
//    	int cnt_opz2 = cnt_op_add; cnt_op_add=0;
    	
    	double log_sumx = Math.log(sum_x);
    	double log_sumz = Math.log(sum_z);
    			
    	double dx1 = (log_sumx1-log_sumx)/log_sumx;
    	double dx2 = (log_sumx2-log_sumx)/log_sumx;
    	System.out.printf("Sum x:\t%.12f\t%.12f\t(%.2f)\t%.12f\t(%.2f)\n", 
    			log_sumx, log_sumx1, dx1, log_sumx2, dx2);
    	double dz1 = (log_sumz1-log_sumz)/log_sumz;
    	double dz2 = (log_sumz2-log_sumz)/log_sumz;
    	System.out.printf("Sum z:\t%.12f\t%.12f\t(%.2f)\t%.12f\t(%.2f)\n", 
    			log_sumz, log_sumz1, dz1, log_sumz2, dz2);
    	
    	System.out.println("add0: "+Logarithms.add(Double.NEGATIVE_INFINITY, y[0])+"\t"+y[0]);
    }
    
	
	

}
