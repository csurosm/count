package count.machine;

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
//			++cnt_op_add;
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
}
