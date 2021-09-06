/*
 * Functions.java
 *
 * Created on November 7, 2003, 2:58 PM
 */

package ca.umontreal.iro.banality;

/**
 * Some functions from Numerical Recipes.
 *
 * @author  csuros
 */
public class Functions {
    
    /**
    * static methods only.
    */
    private Functions(){}


    public static final double EPS=4e-16;

    /**
    * Logarithm of the Gamma function.
    *
    * (based on NR ch.6.1)
    *
    * @param xx positive value
    * @return ln(Gamma(x))
    */
    public static final double gammln(double xx){
        double x,y,tmp,ser;
        int j;
        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*Math.log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += gamma_cof[j]/++y;
        return -tmp+Math.log(2.5066282746310005*ser/x);
    }

    /**
    * coefficients for calculating gammln()
    */
    private static final double gamma_cof[]
        ={76.18009172947146,-86.50532032941677,
          24.01409824083091,-1.231739572450155,
          0.1208650973866179e-2,-0.5395239384953e-5};
      

    private static int factorial_ntop=0;
    private static final int MAX_DIRECTLY_COMPUTED_FACTORIAL = 32;
    private static double factorial_value[]=new double[MAX_DIRECTLY_COMPUTED_FACTORIAL+1];
    /**
     * (based on NR ch. 6.1)
     *
     * @param n non-ngative integer
     * @return the value n! as a floating-point number.
     */
    public static final double factorial(int n){
        int j;
        if (n < 0) throw new IllegalArgumentException("Negative factorial in routine factorial");
        if (n > MAX_DIRECTLY_COMPUTED_FACTORIAL) return Math.exp(gammln(n+1.0));
        if (factorial_ntop<1){
            factorial_value[0]=1.;
            factorial_ntop=1;
        }
        for (;factorial_ntop<=n; factorial_ntop++)
            factorial_value[factorial_ntop]=factorial_value[factorial_ntop-1]*factorial_ntop;
        return factorial_value[n];
    }

    /**
     * (based on NR ch. 6.1)
     *
     * @param n integer
     * @return ln(n!)
     */
    public static final double factln(int n)
    {
        if (n>MAX_DIRECTLY_COMPUTED_FACTORIAL)
        {
            return gammln(n+1.0);
        } else
        {
            return Math.log(factorial(n));
        }
    }

    /**
     * Computes ln (n choose k)
     * @param n
     * @param k
     * @return ln (n!/(k!*(n-k)!))
     */
    public static final double bicoln(int n, int k)
    {
        return factln(n)-factln(k)-factln(n-k);
    }
    
    /**
     * Computes the powers of some number.
     * 
     * @param base base for the exponentiation
     * @param max_exponent maximum value of k for which base^k is computed (inclusive)
     * @return array of [1.0, base, base^2, ..., base^max_exponent]
     */
    public static final double[] powers(double base, int max_exponent)
    {
        double[] A = new double[max_exponent+1];
        A[0] = 1.0;
        for (int i=1; i<=max_exponent;i++)
            if (i%2==0)
                A[i] = A[i/2]*A[i/2];
            else
                A[i]=A[i-1]*base;
        return A;
    }
    

    /**
    * Used by betai: evaluates continued fraction for incomplete beta function
    * by modifed Lentz's method.
    *
    * (NR ch. 6.4)
    */
    public static final double betacf(double a, double b, double x){
        int m,m2;
        double aa,c,d,del,h,qab,qam,qap;
        qab=a+b;
        qap=a+1.0;
        qam=a-1.0;
        c=1.0; // First step of Lentz's method.
        d=1.0-qab*x/qap;
        if (Math.abs(d) < FPMIN) d=FPMIN;
        d=1.0/d;
        h=d;
        for (m=1;m<=MAXIT;m++) {
            m2=2*m;
            aa=m*(b-m)*x/((qam+m2)*(a+m2));
            d=1.0+aa*d; // One step (the even one) of the recurrence.
            if (Math.abs(d) < FPMIN) d=FPMIN;
            c=1.0+aa/c;
            if (Math.abs(c) < FPMIN) c=FPMIN;
            d=1.0/d;
            h *= d*c;
            aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
            d= 1.0+aa*d; // Next step of the recurrence (the odd one).
            if (Math.abs(d) < FPMIN) d=FPMIN;
            c=1.0+aa/c;
            if (Math.abs(c) < FPMIN) c=FPMIN;
            d=1.0/d;
            del=d*c;
            h *= del;
            if (Math.abs(del-1.0) < EPS) break; // Are we done?
        }
        if (m > MAXIT)
          throw new ArithmeticException("a or b too big, or MAXIT too small in betacf");
        return h;
    }

    private static final double FPMIN=Double.MIN_VALUE/EPS;
    private static final int MAXIT=100;

    /**
    * Incomplete beta function.
    *
    * @param x must be between 0 and 1
    * (NR ch. 6.4)
    */
    public static final double betai(double a, double b, double x){
        double bt;
        if (x < 0.0 || x > 1.)
          throw new IllegalArgumentException("Bad x in routine betai.");
        if (x == 0.0 || x == 1.0) bt=0.0;
        else
          bt=Math.exp(gammln(a+b)-gammln(a)-gammln(b)+a*Math.log(x)+b*Math.log(1.0-x));
        if (x < (a+1.0)/(a+b+2.0))
          return bt*betacf(a,b,x)/a;
        else
          return 1.0-bt*betacf(b,a,1.0-x)/b;
    }

    /**
    * @return integral(-t,t) Student_df(x)
    */
    public static final double Student_t(int degrees_of_freedom, double t){
    return betai(degrees_of_freedom/2., .5, degrees_of_freedom/(degrees_of_freedom+t*t));
    }

    /**
    * Complementary error function erfc(x).
    *
    * \frac{2}{\sqrt{\pi}} \int_{x}&{\infty} e^(-t^2) dt
    *
    * (NR ch. 6.2)
    */
    public static final double erfcc(double x){
        double t,z=x,ans;
        if (x<0.)
          z=-x;
        t=1.0/(1.0+0.5*z);
        ans=t*Math.exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
        t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
        t*(-0.82215223+t*0.17087277)))))))));
        return x >= 0.0 ? ans : 2.0-ans;
    }
    
    public static final double log_erfcc(double x)
    {
        if (x<0.)
            throw new IllegalArgumentException("Functions.log_erfcc: argument must be non-negative.");

        double t=1.0/(1.0+0.5*x);
        double ans=Math.log(t)+(-x*x-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
         t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
         t*(-0.82215223+t*0.17087277)))))))));
      
        return ans;
    }
    
    /**
     * Cumulative distribution function for the standard normal distribution
     */
    public static final double normal_cdf(double z){
        if (z>0)
            return 1.0-0.5*erfcc(z/SQRT2);
        else 
            return 0.5*erfcc(-z/SQRT2);
    }
    
    public static final double log_normal_cdf(double z)
    {
        if (z>=0.)
            throw new IllegalArgumentException("Functions.log_normal_cdf: argument should be negative.");
        return Math.log(0.5)+log_erfcc(-z/SQRT2);
    }
    
    
    public final static double SQRT2=Math.sqrt(2.);
    
    
    /**
     * (NR ch. 6.2)
     *
     * @return incomplete gamma function Q(a,x)=1-P(a,x).
     */
    public static final double gammq(double a, double x){
        if (x < 0.0 || a <= 0.0) throw new IllegalArgumentException("Invalid arguments in routine gammq");
        if (x < (a+1.0)) { //Use the series representation
            return 1.0-gser(a,x);
        } else { // Use the continued fraction representation.
            return gcf(a,x);
        }
    }
    
    /**
     * (NR ch. 6.2)
     *
     * @return incomplete gamma function P(a,x)
     */
    public static final double gammp(double a, double x)
    {
        if (x < 0.0 || a <= 0.0) throw new IllegalArgumentException("Negative x or non-positive a in routine gammp [a="+a+", x="+x+"]");
        if (x<a+1.0)
            return gser(a,x);
        else
            return 1.0-gcf(a,x);
    }
    
    /**
     * NR ch. 6.2
     *
     * @return the incomplete gamma function P(a, x) evaluated by its series representation.
     */ 
    public static final double gser(double a, double x){
        double gln=gammln(a);
        if (x <= 0.0) {
            if (x < 0.0) throw new IllegalArgumentException("x less than 0 in routine gser");
            return 0.0;
        } else {
            double ap=a;
            double del=1.0/a;
            double sum=del;
            for (int n=1;n<=MAXIT;n++) {
                ++ap;
                del *= x/ap;
                sum += del;
                if (Math.abs(del) < Math.abs(sum)*EPS) {
                    return sum*Math.exp(-x+a*Math.log(x)-(gln));
                }
            }
            throw new ArithmeticException("a too large, ITMAX too small in routine gser");
        }
    }
    
    /**
     * (NR. ch 6.2)
     *
     * @return the incomplete gamma function Q(a, x) evaluated by its continued fraction representation
     */
    public static final double gcf(double a, double x){
        double gln=gammln(a);
        double b=x+1.0-a; // Set up for evaluating continued fraction by modified Lentz's method with b0=0.
        double c=1.0/FPMIN;
        double d=1.0/b;
        double h=d;
        int i;
        for (i=1;i<=MAXIT;i++) { // Iterate to convergence.
            double an = -i*(i-a);
            b += 2.0;
            d=an*d+b;
            if (Math.abs(d) < FPMIN) d=FPMIN;
            c=b+an/c;
            if (Math.abs(c) < FPMIN) c=FPMIN;
            d=1.0/d;
            double del=d*c;
            h *= del;
            if (Math.abs(del-1.0) < EPS) break;
        }
        if (i > MAXIT) throw new ArithmeticException("a too large, ITMAX too small in gcf");
        return Math.exp(-x+a*Math.log(x)-(gln))*h; // Put factors in front.
    }
    
    /**
     * @return the value of \sum_i=0^{k-1} Pr(X=i) where X is distributed by Poisson(lambda).
     */
    public static final double Poisson_cumulative(double lambda, int k){
        return gammq(k,lambda);
    }
    
    /**
     * @return the value of \sum_i=k^\infty Pr(X=i) where X is distributed by Poisson(lambda).
     */
    public static final double Poisson_tail(double lambda, int k){
        if (lambda < 0.0 || k <= 0.0) throw new IllegalArgumentException("Invalid arguments in routine Poisson_tail");
        if (lambda < (k+1.0)) { //Use the series representation
            return gser(k,lambda);
        } else { // Use the continued fraction representation.
            return 1.0-gcf(k,lambda);
        }        
    }
    
    /**
     * Probability that chi-square random variable exceeds a particular value
     * @param nu number of degrees of freedom
     * @param chi_square threshold
     */
    public static final double Chi_square_tail(double nu, double chi_square)
    {
        return gammq(0.5*nu, 0.5*chi_square);
    }
    
    /**
     * Computes the chi-squared test statistic.
     * Based on Numerical Recipes 14.3
     *
     * @param data binned data
     * @param expected expected values in the bins by the model
     * @param num_constraints number of constraints that determine the degrees of freedom: 
     *     if expected is computed from data, then +1, plus number of parameters fitted
     * @return nothing; sets the test_statistic --- chi-squure statistic and P-value (entries 0 and 1)
     */ 
    public static final void Chi_square_test(double[] data, double[] expected, int num_constraints, double[] test_statistic)
    {
        int n = data.length;
        int df= n-num_constraints;
        double chisq=0.0;
        for (int i=0; i<n; i++)
        {
            double d = (data[i]-expected[i]);
            chisq += d*d/expected[i];
        }
        double p = Chi_square_tail(df, chisq);
        test_statistic[0] = chisq;
        test_statistic[1] = p;
    }
    
    
}
