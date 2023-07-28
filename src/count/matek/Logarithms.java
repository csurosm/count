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

import java.util.Arrays;

/**
 * Arithmetics on log-scale. 
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
	
//	private static int cnt_op_add = 0; // no op counting bc multithreading

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
		
		double sum;
		if (n<=FAST_CUTOFF)
		{
			sum = x[0];
			for (int i=1; i<n; i++)
				sum = add(sum, x[i]);
		} else
		{
			sum = fastSum(x.clone(), n);
//			System.out.println("#**L.sum  n "+n+"\tsum "+sum+"\tfsum "+fsum);
		}
		return sum;
	}
	
	private static final int FAST_CUTOFF = 6; // 6 gives 33%, 14 gives 16% speedup on 59-leaf tree //Integer.MAX_VALUE; // 
	
	private static final double LOG_EPS = 53.0*Math.log(0.5); //  machine epsilon half

	/**
	 * Sums multiple log-scale values, from largest to smallest; terminating when remaining small elements will 
	 * make no contribution at double precision. 
	 * 
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
        int m=n; // number of summing terms on a max-heap
        while (m>0)
        {
        	double v = x[0]; // peek: max remaining element
        	
        	if (sum == Double.NEGATIVE_INFINITY) // == Math.log(0.0)
        		sum = v;
        	else 
        	{ 
            	double d = v-sum; // log-ratio of new term and partial sum
            	double logm = m<logn.length?logn[m]:logn(m); // synchronized function is called only if m is too large
            	double tail_bound = d+logm; // maximum contribution of remaining elements
            	if (tail_bound<LOG_EPS) // would not change sum
            	{
//        			System.out.println("#**L.fastsum  m "+m+"/"+n+"\tsum "+sum+"\td "+d);
            		break;
            	}

            	double a = Math.log1p(Math.exp(d));

            	
	        	sum += a;
        	}
        	--m;
        	// deletion of term at x[0]: decrease heap size and sink from position 0
        	double y = x[m];
        	x[m]=v; // saving the element in the position after the heap, but we could just forget about it; at the end the x[] suffix is sorted, and the prefix is heap-ordered small terms 
        	sink(x, y, 0, m);
        }
        return sum;
        
	}
	
	/**
	 * Requesting a log-value that is not precalculated: extend the {@link #logn} array, and fill it in. 
	 * Synchronized, so thread-blocking. Likely, never called bc logn[] is instantiated with a large capacity already.  
	 * 
	 * @param n
	 * @return
	 */
	private static synchronized double logn(int n)
	{
		int len = logn.length;
		if (n<len) return logn[n];
		
		
		int cap = len;
		do
		{
			cap += cap%3==0?cap/3:cap/2;
		} while (cap <= n);
		
		logn = Arrays.copyOf(logn, cap);
		for (int i=len; len<cap; len++)
			logn[i]=Math.log(i);
		return logn[n];
	}
	
	/**
	 * Precalculated values of log(0)=NEGATIVE_INFINITY, log(1)=0, log(2), log(3), ...
	 */
	private static double[] logn;
	static 
	{
		final int maxn = 1<<16; // power of 2 ; it is unlikely we need to sum 64k terms bc the term count is bounded by the sum of copy numbers in a family
		logn=new double[maxn];
		for (int n=0; n<logn.length; n++)
		{
			logn[n]=Math.log(n);
		}
//		System.out.println("#**Logarithms.fastSum cutoff "+FAST_CUTOFF);
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
    
    
    /**
     * Test code.
     * 
     * @param args
     */
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
