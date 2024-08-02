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
 * Arithmetics on log- and logistic scales. 
 * 
 */
public class Logarithms 
{

	public static double LN2 = Math.log(2.0);
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
	
	/**
	 * Calculates <var>x</var> = logit(<var>p</var>) from <var>p</var>.
	 * 
	 * @param p between 0.0 and 1.0 (inclusively)
	 * @return logit(<var>p</var>)
	 */
	public static double toLogit(double p)
	{
		double x = Math.log(p)-Math.log1p(-p);
		return x;
	}
	
	/**
	 * Calculates ln <var>p</var> from <var>x</var> = logit(<var>p</var>)= ln(<var>p</var>)-ln(1-<var>p</var>) 
	 * 
	 * @param x
	 * @return
	 */
	public static double logitToLogValue(double x)
	{
		double logp;
		if (0.0<=x)
		{
			logp = -Math.log1p(Math.exp(-x));
		} else
		{
			logp = x-Math.log1p(Math.exp(x));
		}
		return logp;
	}
	
	/**
	 * Calculates ln(1-<var>p</var>) from <var>x</var> = logit(<var>p</var>)= ln(<var>p</var>)-ln(1-<var>p</var>) 
	 * 
	 * @param x
	 * @return
	 */
	public static double logitToLogComplement(double x)
	{
		double log1_p;
		if (0.0<=x)
		{
			log1_p = -x-Math.log1p(Math.exp(-x));
		} else
		{
			log1_p = -Math.log1p(Math.exp(x));
		}
		return log1_p;
	}
	
	/**
	 * Calculates ln(-ln(1-<var>p</var>)) from <var>x</var> = logit(<var>p</var>)= ln(<var>p</var>)-ln(1-<var>p</var>) 
	 * 
	 * @param x
	 * @return
	 */
	public static double logitToLogLogComplement(double x)
	{
		double log1_p = logitToLogComplement(x);
		double loglog1_p = Math.log(-log1_p);
		if (loglog1_p==Double.NEGATIVE_INFINITY)
		{
			// p is so small that ln(1-p) = -p ; or p==0.0
			double log_p = logitToLogValue(x);
			loglog1_p = log_p;
		}
		return loglog1_p;
	}
	
	
	
	/**
	 * Multiplication  on logistic scale: 
	 * calculates logit(ab) = ln (ab/(1-ab)) when we know x=logit(a)
	 * and y=logit(b).
	 * 
	 * @param x
	 * @param y
	 * @return result of the multiplication on logistic scale
	 */
	public static double mulLogit(double x, double y)
	{
		double z;
		
		if (x<=y)
		{
			if (0<=y) // x<=y, x<=x+y
			{
				if (Double.isInfinite(y)) return x;
				z = x-Math.log1p(Math.exp(x-y)+Math.exp(-y));
			} else // x and y are negative: x+y <= x<= y<0
			{
				z = (x+y)-Math.log1p(Math.exp(x)+Math.exp(y)); // ok if y=-infty or x=-infty
			}
		} else
		{
			if (0<=x) // y<x, y<= x+y
			{
				if (Double.isInfinite(x)) return y;
				z = y-Math.log1p(Math.exp(y-x)+Math.exp(-x));
			} else // x and y are negative: x+y < y < x <0
			{
				z = (x+y)-Math.log1p(Math.exp(x)+Math.exp(y)); // ok if y=-infty or x=-infty
			}
		}
		
		return z; 
	}
	
	/**
	 * Calculates ln(1+exp(x))-ln(1+exp(y)); i.e. 
	 * the log-ratio of the complement parameters ln((1-q)/(1-p))
	 * with x=logit(p) and y=logit(q)
	 * 
	 * 
	 * @param x
	 * @param y
	 * @return negative value if y&le;x
	 */
	public static double logLogitComplementRatio(double x, double y)
	{
		
		double log_ratio;
		if (x<y)
		{
			log_ratio = -logLogitComplementRatio(y,x);
		} else if (x==y)
		{
			log_ratio = 0.0;
		}
		else // y<x 
		{
			if (0<=x)
			{
				if (0<=y)
				{
					double f = Math.exp(-y);
					double g = -Math.expm1(y-x);
					double h = 1.0+Math.exp(-x);
					log_ratio = (x-y)-Math.log1p(f*g/h);
					
//					System.out.println("#**L.lLCR x "+x+"\ty "+y+"\tf "+f+"\tg "+g+"\th "+h+"\tlr "+log_ratio);
					
				} else if (-x<=y) // y<0, 0<=x+y
				{
					double f = Math.exp(y);
					double g = -Math.expm1(-x-y);
					double h = 1.0+Math.exp(-x);
					log_ratio = x-Math.log1p(f*g/h);
				} else // y<-x; x+y<0
				{
					double f = Math.exp(-x);
					double g = -Math.expm1(x+y);
					double h = 1.0+Math.exp(y);
					log_ratio = x+Math.log1p(f*g/h);
				}
			} else // y<=x<0
			{
				double f = Math.exp(x);
				double g = -Math.expm1(y-x);
				double h = 1.0+Math.exp(y);
				log_ratio = Math.log1p(f*g/h);
			}
		}
		return log_ratio;
	}
	
	/**
	 * Calculates (exp(x)-exp(y))/(exp(x)+exp(x+y)); i.e.
	 * the logarithm of the complement of parameter ratios ln(1-q/p)
	 * with x=logit(p) and y=logit(q); y&le;x
	 * 
	 * @param x larger
	 * @param y smaller 
	 */
	public static double logLogitRatioComplement(double x, double y)
	{
		assert (y<=x);
		double log_complement;
		if (x==y)
		{
			log_complement = Double.NEGATIVE_INFINITY;
		} else
		{
			if (0<=y)
			{
				log_complement = -y + Math.log(-Math.expm1(y-x))-Math.log1p(Math.exp(-y));
			} else
			{
				log_complement = Math.log(-Math.expm1(y-x))-Math.log1p(Math.exp(y));
			}
		}
		return log_complement;
	}
	
	/**
	 * Calculates (exp(y)+exp(x+y))/(exp(x)-exp(y)); i.e., logit(q/p)
	 * with x=logit(p) and y=logit(q); assuming y&le;x (q&le;p).
	 * 
	 * @param x
	 * @param y
	 * @return
	 */
	public static double logitParameterRatio(double x, double y)
	{
		assert (y<=x);
		double logitR;
		double log_denom = Math.log(-Math.expm1(y-x));
		if (0<=x)
		{
			logitR = y+Math.log1p(Math.exp(-x))-log_denom; 
		} else
		{
			logitR = y-x+Math.log1p(Math.exp(x))-log_denom;
		}
		return logitR;
	}
	
	/**
	 * Calculates ln(1-exp(-logp))
	 * @param log_p log-probability (must be non-positive)
	 * @return
	 */
	public static double logToLogComplement(double log_p)
	{
		assert(log_p<=0.0);
		double log1_p;
		if (-log_p < LN2) // p>1/2 -logp < ln 2
		{
			// p > 1/2; 1-p < 1/2 
			// ln(1-exp(x)) = ln (-(exp(x)-1))
			double one_minus_p =  -Math.expm1(log_p);
			log1_p = Math.log(one_minus_p);
			// if log_p = 0 & p=1, then Math.log(-0.0)=-infty 
		} else
		{
			// p<= 1/2
			// ln(1-exp(x)) = ln(1+(-exp(x)))
			double p = Math.exp(log_p);
			log1_p = Math.log1p(-p);
			// if log_p = -infty & p=0, then Math.log1p(0.0)=0.0
		}
		return log1_p;
	}
	
	/**
	 * Calculates logit (exp(logp))
	 * @param log_p
	 * @return
	 */
	public static double logToLogit(double log_p)
	{
		double log1_p = logToLogComplement(log_p);
		return log_p-log1_p;
	}
	
	
	/*
	 * 
	 * Log-difference methods. A 2-element array is 
	 * used to store [a,b] representing the difference exp(a)-exp(b).
	 */
	
	private static final int LDIFF_POS = 0;
	private static final int LDIFF_NEG = 1;
	
	
	/**
	 * 
	 * Creates a 2-element array that is  
	 * used to store (<var>a</var>,<var>b</var>) 
	 * representing the difference exp(<var>a</var>)-exp(<var>b</var>).
	 * @param a logarithm of positive term
	 * @param b logarithm of negative term
	 * @return
	 */
	public static double[] ldiff(double a, double b)
	{
		double[] ld = new double[2];
		ld[LDIFF_POS] = a;
		ld[LDIFF_NEG] = b;
		return ld;
	}
	
	/**
	 * Log-difference with value 0
	 */
	public static double[] ldiff()
	{
		return ldiff(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
	}
	
	/**
	 * Evaluates the difference, and sets one of the terms to log(0).
	 * Calculates the balance <var>d</var>=log(|exp(<var>a</var>)-exp(<var>b</var>)|),
	 * and sets the entries to (<var>d</var>, -infty) if <var>b</var>&le;<var>a</var>
	 * (non-negative balance),
	 * or to (-infty, <var>d</var>) if <var>a</var>&lt;<var>b</var> (negative balance).
	 * 
	 * 
	 * @param ldiff updated to log-difference from 0 for the same result, 
	 */
	public static void ldiffBalance(double[] ldiff)
	{
		ldiffBalance(ldiff, ldiff);
	}
	
	
	/**
	 * Shifts the log-difference and sets one of the terms to log(0).
	 * Calculates the balance <var>d</var>=log(|exp(<var>a</var>)-exp(<var>b</var>)|),
	 * and sets the result entries to (<var>d</var>, -infty) if <var>b</var>&le;<var>a</var>
	 * (non-negative balance),
	 * or to (-infty, <var>d</var>) if <var>a</var>&lt;<var>b</var> (negative balance).
	 * 
	 * 
	 * @param ldiff to be evaluated
	 * @param result where the resul should be put; if null, a new array is created
	 * @return result with updated entries if non-null, or a new array if input result is null
	 */
	private static double[] ldiffBalance(double[] ldiff, double[] result)
	{
		
		double a = ldiff[LDIFF_POS];
		double b = ldiff[LDIFF_NEG];
		
		boolean is_positive = (b<=a);
		if (!is_positive)
		{
			// exchange a and b 
			double t = a; a=b; b=t;
		}
		double d = b-a;
		
//		if (!(d<=0.0)) // DEBUG
//		{
//			System.out.println("#**L.ldB (a,b) "+a+"\t"+b+"\t"+Arrays.toString(ldiff)+"\td="+d+"\tispos "+is_positive);
//		}
		assert (d<=0.0); // because we exchanged a and b otherwise 
		
		// log(exp(a)-exp(b)) = a + log(1-exp(b-a)) with b<=a 
		double log_balance = a + logToLogComplement(d);
		
		
		if (result == null) result = new double[2];

		if (is_positive)
		{
			result[LDIFF_POS] = log_balance;
			result[LDIFF_NEG] = Double.NEGATIVE_INFINITY;
		} else
		{
			result[LDIFF_POS] = Double.NEGATIVE_INFINITY;
			result[LDIFF_NEG] = log_balance;
		}
		return result;
	}
	
	/**
	 * Evaluates the log-difference (<var>a</var>, <var>b</var>)
	 * as exp(<var>a</var>)-exp(<var>b</var>). 
	 * 
	 * @param ldiff
	 * @return
	 */
	public static double ldiffValue(double[] ldiff)
	{
		
		if (ldiff[LDIFF_POS]!=Double.NEGATIVE_INFINITY
				&& ldiff[LDIFF_NEG]!=Double.NEGATIVE_INFINITY)
		{
			// not balanced
			ldiff = ldiffBalance(ldiff, null); // replace with 0-offset difference
		}
		double ldValue;
		if (ldiff[LDIFF_POS]==Double.NEGATIVE_INFINITY)
		{ // negative balance
			ldValue = -Math.exp(ldiff[LDIFF_NEG]);
		} else
		{
			assert (ldiff[LDIFF_NEG]==Double.NEGATIVE_INFINITY) ;
			ldValue = Math.exp(ldiff[LDIFF_POS]);
		}
		return ldValue;
	}
	
	public static double ldiffLogValue(double[] ldiff)
	{
		if (ldiff[LDIFF_POS]!=Double.NEGATIVE_INFINITY
				&& ldiff[LDIFF_NEG]!=Double.NEGATIVE_INFINITY)
		{
			// not balanced
			ldiff = ldiffBalance(ldiff, null); // replace with 0-offset difference
		}
		double ldiffLogValue = ldiff[LDIFF_POS]; // negative infty if non-positive balance 
		return ldiffLogValue;
	}

	/**
	 * 
	 * Calculates a+b*c on a log-difference scale 
	 * 
	 * @param ld_a log-difference for a
	 * @param log_b log of b 
	 * @param ld_c log-difference for c 
	 * @param result where to put the result; if null, a new array is created
	 * @return result on log-difference scale 
	 */
	public static double[] ldiffAddMultiply(double[] ld_a, double log_b, double[] ld_c, double[] result)
	{
		if (result==null) result = ldiff();
		if (ldiffIsZero(ld_c) || log_b == Double.NEGATIVE_INFINITY)
		{
			result[LDIFF_POS]  = ld_a[LDIFF_POS];
			result[LDIFF_NEG]  = ld_a[LDIFF_NEG];
		} else
		{
			result[LDIFF_POS] =  add(ld_a[LDIFF_POS], log_b+ld_c[LDIFF_POS]);
			result[LDIFF_NEG] =  add(ld_a[LDIFF_NEG], log_b+ld_c[LDIFF_NEG]);
		}
		return result;
	}
	
	public static boolean ldiffIsZero(double[] ldiff)
	{
		return ldiff[LDIFF_POS]==ldiff[LDIFF_NEG]
				;
	}
	
	public static boolean ldiffIsNegative(double[] ldiff)
	{
		return ldiff[LDIFF_POS]<ldiff[LDIFF_NEG];
	}
	
	public static boolean ldiffIsPositive(double[] ldiff)
	{
		return ldiff[LDIFF_NEG]<ldiff[LDIFF_POS];
	}
	
	
	/**
	 * 
	 * Calculates a-b*c on a log-difference scale 
	 * 
	 * @param ld_a log-difference for a
	 * @param log_b log of b 
	 * @param ld_c log-difference for c 
	 * @param result where to put the result; if null, a new array is created
	 * @return result on log-difference scale 
	 */
	public static double[] ldiffSubtractMultiply(double[] ld_a, double log_b, double[] ld_c, double[] result)
	{
		if (result==null) result = ldiff();
		if (ldiffIsZero(ld_c))
		{
			result[LDIFF_POS]  = ld_a[LDIFF_POS];
			result[LDIFF_NEG]  = ld_a[LDIFF_NEG];
		} else
		{
			double lcpos = ld_c[LDIFF_POS];
			double lcneg = ld_c[LDIFF_NEG];
			result[LDIFF_POS] =  add(ld_a[LDIFF_POS], log_b+lcneg);
			result[LDIFF_NEG] =  add(ld_a[LDIFF_NEG], log_b+lcpos);
		}
		return result;
	}
	/**
	 * Calculates a*ax-b*bx on a log-difference scale 
	 * 
	 * @param ld_a log-difference for a 
	 * @param log_ax log of multiplier ax
	 * @param ld_b log-difference for b
	 * @param log_bx log of multiplier bx
	 * @param result where to put the result; if null, a new array is created
	 * @return result on log-difference scale 
	 */
	public static double[] ldiffSubtractMultiply(double log_ax, double[] ld_a, double log_bx, double[] ld_b, double[] result)
	{
		double[] ld_aax = ldiffMultiply(ld_a, log_ax, null);
		result = ldiffSubtractMultiply(ld_aax, log_bx, ld_b, result==null?ld_aax:result);
		return result;
	}	
	
	/**
	 * Calculates a*b on a log-difference scale.
	 *  
	 * @param ld_a log-difference for a 
	 * @param log_b log of multiplier b
	 * @param result
	 * @param result where to put the result; if null, a new array is created
	 * @return result on log-difference scale 
	 */
	public static double[] ldiffMultiply(double[] ld_a, double log_b, double[] result)
	{
		if (result==null) result = ldiff();
		double rpos,rneg;
		if (ldiffIsZero(ld_a))
		{
			rpos = ld_a[LDIFF_POS];
			rneg = ld_a[LDIFF_NEG];
		} else
		{
			rpos = ld_a[LDIFF_POS]+log_b;
			rneg = ld_a[LDIFF_NEG]+log_b;
		}
		
		if (Double.isNaN(rpos) || Double.isNaN(rneg)) // DEBUG
		{
			System.out.println("#**L.ldM lda "+Arrays.toString(ld_a)+"\tlogb "+log_b+"\trpos "+rpos+"\trneg "+rneg);
		}
		result[LDIFF_POS] = rpos ;
		result[LDIFF_NEG] = rneg ;
		
		assert (!Double.isNaN(result[LDIFF_POS]));
		assert (!Double.isNaN(result[LDIFF_NEG]));
		
		return result;
	}
	
	public static double[] ldiffInverse(double[] ld_a, double[] result)
	{
		if (result==null) result = ldiff();
		double ldpos = ld_a[LDIFF_POS];
		double ldneg = ld_a[LDIFF_NEG];
		
		result[LDIFF_POS] = ldneg;
		result[LDIFF_NEG] = ldpos;
		return result;
	}
	
	
	
	
	/*
	 * 
	 * Addition on logarithmic scale
	 * 
	 */

	
	
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
		if (n<=0) return Double.NEGATIVE_INFINITY;
		
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
	
	
	/* Summing in an array on logarithmic scale:  by sorting the elements */
	
	private static final int FAST_CUTOFF = Integer.MAX_VALUE;//   6; // 6 gives 33%, 14 gives 16% speedup on 59-leaf tree //Integer.MAX_VALUE; // 
	
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
    	
    	java.util.Random RND = new java.util.Random(2021);

    	String test = "add";
    	
    	if (args.length>1)
    		test = args[1];
    	
    	if ("add".equals(test))
    	{
    		
	    	double[] x = new double[N];
	    	double[] y = new double[N];
	    	double[] z = new double[N];
	    	
	    	
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
    	else if ("logit".equals(test))
    	{
    		double[] x = new double[N];
    		double[] y = new double[N];
    		double[] p = new double[N];
    		double[] q = new double[N];
	    	for (int i=0; i<N; i++)
	    	{
	    		p[i] = RND.nextDouble();
	    		double t = RND.nextDouble();
	    		if (p[i]<t)
	    		{
	    			q[i]=p[i];
	    			p[i]=t;
	    		} else
	    		{
	    			q[i] = t;
	    		}
	    		x[i] = Math.log(p[i]/(1.0-p[i]));
	    		y[i] = Math.log(q[i]/(1.0-q[i]));
	    	}
	    	
	    	System.out.println("#logitValues\tti\tx[i]\ty[i]\terr1\terr3\tp[i]\tq[i]\tln((1-q)/(1-p)):xy\t==:pq\tlogit(q/p):xy\t==:pq");
	    	for (int i=0; i<N; i++)
	    	{
	    		double lr1 = logLogitComplementRatio(x[i],y[i]);
	    		double lr1d = Math.log1p(-q[i])-Math.log1p(-p[i]);
	    		double diff1 = (lr1-lr1d)/Math.abs(lr1d);
	    		
	    		double lr2 = logLogitRatioComplement(x[i],y[i]);
	    		double lr2d = Math.log1p(-q[i]/p[i]);
	    		double diff2 = (lr2-lr2d)/Math.abs(lr2d);
	    		
	    		double lr3 = logitParameterRatio(x[i],y[i]);
	    		double lr3d = Math.log(q[i]/p[i])-Math.log(1.0-q[i]/p[i]);
	    		double diff3 = (lr3-lr3d)/Math.abs(lr3d);
	    		
	    		System.out.printf("%d\t%.12f\t%.12f\t%.3g\t%.3g\tp %f\tq %f\tlr1 %f/%f\tlr2 %f/%f\t(x+y) %f\n"
	    				, i, x[i], y[i]
	    				, diff1
	    				//, diff2 
	    				, diff3 
	    				, p[i], q[i]
	    				, lr1, lr1d
	    				//, lr2, lr2d
	    				, lr3, lr3d
	    				, x[i]+y[i]
	    				);
	    		
	    	}
    		
    	} else if ("ldiff".equals(test))
    	{
    		double[] A = new double[N];
    		double[] B = new double[N];
    		
    		
    		double[] a = new double[N];
    		double[] b = new double[N];
    		
    		double[][] ldiff = new double[N][];
    		double eps =Math.ulp(1.0);
    		
    		double d = 0.0;
    		
    		System.out.println("#ldiffValue\ti\ta[i]\tb[i]\tldval\texpected\tinverse");
	    	for (int i=0; i<N; i++)
	    	{
	    		A[i] = RND.nextDouble();
	    		B[i] = RND.nextDouble();
	    		a[i] = Math.log(A[i]);
	    		b[i] = Math.log(B[i]);
	    		ldiff[i] = ldiff(a[i], b[i]);
	    		double A_B =  A[i]-B[i];
	    		double ldval = ldiffValue(ldiff[i]);
	    		System.out.printf("#ldiffValue\t%d\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n",
	    				i, a[i], b[i]
	    						, ldval
	    						, A_B
	    						, ldiffValue(ldiffInverse(ldiff[i], null))
	    						);
	    		d += Math.abs((ldval-A_B)/A_B);
	    	}
	    	d/=N;
	    	System.out.printf("ldiffValue error %g\t(%.12f ulp)\n",d,d/eps);

    		System.out.println("#ldiffMul\ti\tc\tldval\texpected");
    		d=0.0;
	    	for (int i=0; i<N; i++ )
	    	{
	    		double c = RND.nextDouble();
	    		double[] ldr = ldiffMultiply(ldiff[i],Math.log(c),null);
	    		double ldval = ldiffValue(ldr);
	    		double r = c*(A[i]-B[i]);
	    		System.out.printf("#ldiffMul\t%d\t%.12f\t%.12f\t%.12f\n"
									, i
									, c
									, ldval, r
				    				);
	    		d = d+Math.abs((ldval-r)/r);
	    	}	    	
	    	d/=N;
	    	System.out.printf("ldiffMul error %g\t(%.12f ulp)\n",d,d/eps);
	    	
    		System.out.println("#ldiffAdd/Sub\ti\tc\tldval\texpected");
    		d=0.0;
	    	for (int i=0; i<N; i+=2 )
	    	{
	    		double [] lda = ldiff[i];
	    		double [] ldb = ldiff[i+1];
	    		double c = 2.0*(RND.nextDouble()-0.5); // between 1 and -1 
	    		double[] ldr;
	    		String method;
	    		if (c<0.0)
	    		{
	    			ldr = ldiffSubtractMultiply(lda, Math.log(-c), ldb, null);
	    			method = "ldiffSub"; // should be via Logarithms.getClass().getMethod(...).toString
	    		} else 
	    		{
	    			ldr = ldiffAddMultiply(lda, Math.log(c), ldb, null);
	    			method = "ldiffAdd";
	    		}
	    		double r = (A[i]-B[i])+c*(A[i+1]-B[i+1]);
	    		
	    		double ldval = ldiffValue(ldr);
	    		System.out.printf("#ldiff%s\t%d\t%.12f\t%.12f\t%.12f\n"
	    							, method, i
	    							, c
	    							, ldval, r
	    				);
	    		d += Math.abs((ldval-r)/r);
	    	}
	    	d/=0.5*N;
	    	System.out.printf("ldiffAdd/Sub error %g\t(%.12f ulp)\n",d,d/eps);
	    	
    		System.out.println("#ldiffBMul\ti\tc\tldval\texpected");
    		d=0.0;
	    	for (int i=0; i<N; i++ )
	    	{
	    		double c = RND.nextDouble();
	    		ldiffBalance(ldiff[i]);
	    		double[] ldr = ldiffMultiply(ldiff[i],Math.log(c),null);
	    		double ldval = ldiffValue(ldr);
	    		double r = c*(A[i]-B[i]);
	    		System.out.printf("#ldiffBMul\t%d\t%.12f\t%.12f\t%.12f\n"
									, i
									, c
									, ldval, r
				    				);
	    		d = d+Math.abs((ldval-r)/r);
	    	}	    	
	    	d/=N;
	    	System.out.printf("ldiffBalancedMultiply error %g\t(%.12f ulp)\n",d,d/eps);
	    	
    		System.out.println("#ldiffBAdd/Sub\ti\tc\tldval\texpected");
    		d=0.0;
	    	for (int i=0; i<N; i+=2 )
	    	{
	    		double [] lda = ldiff[i];
	    		double [] ldb = ldiff[i+1];
	    		double c = 2.0*(RND.nextDouble()-0.5); // between 1 and -1 
	    		double[] ldr;
	    		String method;
	    		if (c<0.0)
	    		{
	    			ldr = ldiffSubtractMultiply(lda, Math.log(-c), ldb, null);
	    			method = "ldiffBSub"; // should be via Logarithms.getClass().getMethod(...).toString
	    		} else 
	    		{
	    			ldr = ldiffAddMultiply(lda, Math.log(c), ldb, null);
	    			method = "ldiffBAdd";
	    		}
	    		double r = (A[i]-B[i])+c*(A[i+1]-B[i+1]);
	    		
	    		double ldval = ldiffValue(ldr);
	    		System.out.printf("#ldiff%s\t%d\t%.12f\t%.12f\t%.12f\n"
	    							, method, i
	    							, c
	    							, ldval, r
	    				);
	    		d += Math.abs((ldval-r)/r);
	    	}
	    	d/=0.5*N;
	    	System.out.printf("ldiffBalancedAdd/Sub error %g\t(%.12f ulp)\n",d,d/eps);
    	} else
    	{
    		throw new IllegalArgumentException("Use as "+Logarithms.class.toString()+" [N [test]]; with test {\"add\",\"logit\",\"ldiff\"} ");
    	}
    }
    
	
	

}
