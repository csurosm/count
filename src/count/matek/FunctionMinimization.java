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
import java.util.List;
import java.util.function.DoubleFunction;
import java.util.function.Function;



/**
*
* This class implements a number of numerical optimization methods 
* for function minimization and root finding in one or multiple dimensions.
* 
* @author Mikl&oacute;s Cs&#369;r&ouml;s
* @since October 26, 2004, 8:22 PM
*/
public class FunctionMinimization 
{
	/**
	 * Messages on standard output / console.
	 */
	private static final boolean REPORT_UNUSUAL = false; 

    private FunctionMinimization(){}
    
 
    
    /**
     * Machine precision: smallest value with
     * 1.0+EPS &gt; 1.
     */
    
    
    public static double EPS = 1.0/(1L<<52); // 52 bits: 1/2^52 = 2.22 10^-16
    
    /**
     * Convergence criterion on x values ({@link #dfpmin(double[], double, int, Function, Function, List)})
     */
    public static double DFP_TOLX = (4.*EPS); 

    /**
     * Convergence criterion on delta x ({@link #lnsrch(double[], double, double[], double[], double[], double[], double, Function)})
     */
    public static double LNSRCH_TOLX=DFP_TOLX; // 1.0e-7;
    
    
    /**
     * Maximum number of iterations in {@link #zbrent(DoubleFunction, double, double, double)}
     * and {@link #brent(double, double, double, DoubleFunction, double)}
     */
    private static final int BRENT_ITMAX=100;


    private static double sign(double a, double b)
    {
        if (b>=0.0)
            return Math.abs(a);
        else
            return -Math.abs(a);
    }
    
    private static double SQR(double x) { return x*x;}
    
    /**
     * Length/magnitude of the vector p[0..n-1]
     * @param p
     * @param n
     * @return
     */
    public static double euclideanNorm(double[] p, int n)
    {
    	assert (n<=p.length);
    	double sum = 0.0;
    	for (int i=0; i<n; i++)
    	{
    		sum +=p[i]*p[i];
    	}
    	double euclideanNorm = Math.sqrt(sum);
    	return euclideanNorm;
    }
    /**
     * Vector length/magnitude. 
     * 
     * @param p
     * @return
     */
    public static double euclideanNorm(double[] p)
    {
    	return euclideanNorm(p, p.length);
    }
    
    
    /**
     * Root finding: using Brent's method, find the root of a function func known to lie between x1 and x2. 
     * The root will be refined until its accuracy is tol. Specifically, the 
     * iteration stops when the root is bracketed within an interval of length tol.
     *
     * Based on NR 9.3.
     *
     * @param func function to be minimized
     * @param x1 left bracket for x value
     * @param x2 right bracket for x value
     * @param tol bracketing accuracy
     * 
     * @return x value that zeros the function within the bracket
     */ 
    public static double zbrent(DoubleFunction<Double> func, double x1, double x2, double tol) 
    { 
        double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2; 
        double fa=func.apply(a), fb=func.apply(b),fc,p,q,r,s,tol1,xm; 
        if((fa>0.0&&fb> 0.0)||(fa< 0.0&&fb<0.0)) 
            throw new OptimizationException("Root must be bracketed in zbrent"); 
        fc=fb; 
        for(int iter=0;iter<BRENT_ITMAX;iter++)
        { 
            if((fb>0.0 && fc>0.0) || (fb<0.0 && fc< 0.0))
            { 
                c=a; // Rename a, b,c and adjust bounding interval D
                fc=fa; 
                e=d=b-a; 
            } 
            if(Math.abs(fc)<Math.abs(fb))
            { 
                a=b; 
                b=c; 
                c=a; 
                fa=fb; 
                fb=fc; 
                fc=fa; 
            } 
            tol1=2.0*EPS*Math.abs(b)+0.5*tol; // Convergence check. 
            xm=0.5*(c-b); 
            if(Math.abs(xm)<=tol1 || fb==0.0)
                return b; 
            if(Math.abs(e)>=tol1 && Math.abs(fa)> Math.abs(fb))
            { 
                s=fb/fa; // Attempt inverse quadratic interpolation. 
                if(a==c)
                { 
                    p=2.0*xm*s; 
                    q=1.0-s; 
                } else
                { 
                    q=fa/fc; 
                    r=fb/fc; 
                    p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0)); 
                    q=(q-1.0)*(r-1.0)*(s-1.0); 
                } 
                if(p>0.0)
                    q=-q; // Check whether in bounds. 
                p=Math.abs(p); 
                min1=3.0*xm*q-Math.abs(tol1*q); 
                min2=Math.abs(e*q); 
                if(2.0*p<(min1<min2? min1:min2))
                { 
                    e=d; // Accept interpolation. 
                    d=p/q; 
                } else
                { 
                    d=xm; // Interpolation failed, use bisection. 
                    e=d; 
                } 
            } else { // Bounds decreasing too slowly, use bisection. 
                d=xm; 
                e=d;
            } 
            a=b; // Move last best guess to a. 
            fa=fb; 
            if(Math.abs(d)>tol1) // Evaluate new trial root. 
                b +=d; 
            else 
                b +=sign(tol1,xm); 
            fb=func.apply(b); 
        } // iter
        throw new OptimizationException("Maximum number of iterations exceeded in zbrent"); 
    }
    
    /**
     * The golden ratio.
     */
    private static final double R=0.5*(Math.sqrt(5.0)-1.0); // 0.618
    /**
     * Complement golden ratio 1.0-{@link #R}
     */
    private static final double C=(1.0-R);


    /**
     * Golden Section search
     * 
     * Based on Numerical Recipes 10.1.
     *
     * Given a function func, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
     * between ax and cx, and func(bx) is less than both func(ax) and func(cx)), this routine performs a
     * golden section search for the minimum, isolating it to a fractional precision of about tol. 
     *
     * @param ax bracket interval endpoint
     * @param bx bracket midpoint 
     * @param cx bracket interval endpoint
     * @param func the function to be minimized
     * @param tol relative length for the smaller bracket length for stopping the bracket shrinking (set to sqrt(precision))
     * @return [x,y] where y is minimum function value, and x is its abscissa
     */
    public static double[] golden(double ax, double bx, double cx, DoubleFunction<Double> func, double tol){
        double x0=ax; // At any given time we will keep track of four points, x0,x1,x2,x3. 
        double x3=cx;
        double x1, x2;
        if (Math.abs(cx-bx) > Math.abs(bx-ax)) { // Make x0 to x1 the smaller segment, and fill in the new point to be tried.
            x1=bx;
            x2=bx+C*(cx-bx); 
        } else {
            x2=bx;
            x1=bx-C*(bx-ax);
        }
        // The initial function evaluations. Note that we never need to evaluate the function
        // at the original endpoints.
        double f1=func.apply(x1); 
        double f2=func.apply(x2);
        
        while (Math.abs(x3-x0) > tol*(Math.abs(x1)+Math.abs(x2))) {
            if (f2 < f1) { // One possible outcome,
                x0=x1; x1=x2; x2=R*x1+C*x3;
                f1=f2; f2=func.apply(x2);
            } else {
                x3=x2; x2=x1; x1=R*x2+C*x0;
                f2=f1; f1=func.apply(x1);
            }
        } 
        // Back to see if we are done.
        double[] retval=new double[2];
        if (f1 < f2) { // We are done. Output the best of the two current values. 
            retval[0]=x1;
            retval[1]=f1;
        } else {
            retval[0]=x2;
            retval[1]=f2;
        }
        return retval;
    }
            
    
    /** Used in {@link #brent(double, double, double, DoubleFunction, double)}, {@link #ZEPS} 
     * is a small number that protects against trying to achieve fractional accuracy for a minimum that
     * happens to be exactly zero.
     */
    private static final double ZEPS=1.0e-10;    
    /**
     * Brent's method of function minimization
     *
     * 
     * Based on Numerical Recipes 10.2.
     *
     * Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
     * between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
     * the minimum to a fractional precision of about tol using Brent's method. 
     * 
     * @param ax bracket interval endpoint
     * @param bx bracket midpoint 
     * @param cx bracket interval endpoint
     * @param func the function to be minimized
     * @param tol relative length for the smaller bracket length for stopping the bracket shrinking (set to sqrt(precision))
     * @return [x,y] where y is minimum function value, and x is its abscissa
     */

    public static double[] brent(double ax, double bx, double cx, DoubleFunction<Double> func, double tol){
        double x,w,v,fx,fw,fv,u,fu;
        double e=0.0;
        double a=(ax < cx ? ax : cx); // a and b must be in ascending order, but input abscissas need not be. 
        double b=(ax > cx ? ax : cx);
        x=w=v=bx; 
        fw=fv=fx=func.apply(x);

        double d=0.0;
        for (int iter=1;iter<=BRENT_ITMAX;iter++) { // Main program loop.
            double xm=0.5*(a+b);
            double tol1;//=tol*Math.abs(x)+ZEPS;
            double tol2=2.0*(tol1=tol*Math.abs(x)+ZEPS);
            double threshold=tol2-0.5*(b-a);
            //System.out.println("#**B.brent "+iter+" x "+x+", xm "+xm+"; tol1 "+tol1+", tol2 "+tol2+"; a "+a+", b "+b+"; thresh "+threshold+", xdiff "+Math.abs(x-xm));
            if (Math.abs(x-xm) <= threshold) { // Test for done here.
                double[] retval=new double[2];
                retval[0]=x;
                retval[1]=fx;
                return retval;
            }
            if (Math.abs(e) > tol1) { // Construct a trial parabolic fit.
                double r=(x-w)*(fx-fv);
                double q=(x-v)*(fx-fw);
                double p=(x-v)*q-(x-w)*r;
                q=2.0*(q-r);
                if (q > 0.0) p = -p;
                q=Math.abs(q);
                double etemp=e;
                e=d;
                if (Math.abs(p) >= Math.abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                    d=C*(e=(x >= xm ? a-x : b-x));
                    // The above conditions determine the acceptability of the parabolic fit. Here we
                    // take the golden section step into the larger of the two segments.
                else {
                    d=p/q;//  Take the parabolic step.
                    u=x+d;
                    if (u-a < tol2 || b-u < tol2)
                    d=sign(tol1,xm-x);
                }
            } else {
                e=(x >= xm ? a-x : b-x);
                d=C*e;
            }
            u=(Math.abs(d) >= tol1 ? x+d : x+sign(tol1,d));
            fu=func.apply(u);
            // This is the one function evaluation per iteration.
            if (fu <= fx) { // Now decide what to do with our function evaluation. 
                //System.out.println("#**B.brent "+iter+" fu<=fx u "+u+", fu "+fu+"; fx "+fx+"; v "+v+", w "+w);
                if (u >= x) a=x; else b=x;
                v=w; w=x; x=u;
                fv=fw; fw=fx; fx=fu;
            } else {
                //System.out.println("#**B.brent "+iter+" fu>>fx u "+u+", fu "+fu+"; fx "+fx+"; v "+v+", w "+w);
                if (u < x) a=u; else b=u;
                if (fu <= fw || w == x) {
                    v=w; w=u;
                    fv=fw; fw=fu;
                } else if (fu <= fv || v == x || v == w) {
                    v=u;
                    fv=fu;
                }
            } // Done with housekeeping. Back for another iteration. 
        }
        throw new OptimizationException("Too many iterations in brent");
    }
    
    /**
     * Brent's method of function minimization
     *
     * 
     * Based on Numerical Recipes 10.2.
     *    Given a function f and its derivative function df,
     *     and given a bracketing triplet of abscissas ax, bx, cx 
     *     [such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)], 
     *     this routine isolates the minimum to a fractional precision of about tol 
     *     using a modification of Brent’s method that uses derivatives. 
     *     
     * @param ax bracket interval endpoint
     * @param bx bracket midpoint 
     * @param cx bracket interval endpoint
     * @param f the function to be minimized
     * @param df first derivative of f
     * @param tol relative length for the smaller bracket length for stopping the bracket shrinking (set to sqrt(precision))
     * @return [x,y] where y is minimum function value, and x is its abscissa
     */     

    public static double[] dbrent(double ax, double bx, double cx, DoubleFunction<Double> f, DoubleFunction<Double> df, double tol){

    	double a,b;
    	double d=0.0, e=0.0;
    	double du,dv,dw,dx,fu,fv,fw,fx;
    	double olde;
    	double u,u1,u2,v,w,x;
    	
    	
    	a=(ax < cx ? ax : cx); 
    	b=(ax > cx ? ax : cx); 
    	x=w=v=bx; 
    	fw=fv=fx=f.apply(x); // All our housekeeping chores are doubled by the necessity of 
    	dw=dv=dx=df.apply(x); // moving derivative values around as well as function values.
    	
    	for (int iter=1;iter<=BRENT_ITMAX;iter++){ 
    		double xm=0.5*(a+b);
	        double tol1;//=tol*Math.abs(x)+ZEPS;
	        double tol2=2.0*(tol1=tol*Math.abs(x)+ZEPS);
    	
            double threshold=tol2-0.5*(b-a);
            //System.out.println("#**B.brent "+iter+" x "+x+", xm "+xm+"; tol1 "+tol1+", tol2 "+tol2+"; a "+a+", b "+b+"; thresh "+threshold+", xdiff "+Math.abs(x-xm));
            if (Math.abs(x-xm) <= threshold) { // Test for done here.
                double[] retval=new double[2];
                retval[0]=x;
                retval[1]=fx;
                return retval;
            }
            if (Math.abs(e) > tol1) {
        		double d1=2.0*(b-a); // Initialize these d’s to an out-of-bracket value.
        		double d2=d1;
    			if (dw!= dx) d1=(w-x)*dx/(dx-dw); // Secant method with one point.
				if (dv!= dx) d2=(v-x)*dx/(dx-dv); // And the other.
				// Which of these two estimates of d shall we take?
				// We will insist that they be within the bracket, and on the 
				// side pointed to by the derivative at x:
				u1=x+d1;
				u2=x+d2;
				boolean ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0; 
				boolean ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
				// Movement on the step before last. 
				olde=e;
				e=d;
				if (ok1 || ok2) {
					// Take only an acceptable d,
					if (ok1 && ok2)  // and if both are acceptable, then take 
						d=(Math.abs(d1) < Math.abs(d2) ? d1 : d2); // the smallest one.
					else if (ok1) d=d1;
					else d=d2; // ok2
					if (Math.abs(d) <= Math.abs(0.5*olde)) { 
						u=x+d;
						if (u-a < tol2 || b-u < tol2) d=sign(tol1,xm-x);
					} else {
						// Bisect, not golden section. 
						d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
						// Decide which segment by the sign of the derivative.
					}
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x)); }
            } else {
            	d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
            }
            if (Math.abs(d) >= tol1) {
            	u=x+d;
            	fu=f.apply(u); 
            } else {
            	u=x+sign(tol1,d); fu=f.apply(u);
            	// If the minimum step in the downhill direction takes us uphill, then we are done.    	
            	if (fu > fx) {
                    double[] retval=new double[2];
                    retval[0]=x;
                    retval[1]=fx;
                    return retval;
            	}
            }
            du=df.apply(u);
            //   Now all the housekeeping, sigh.
            if (fu <= fx) {
            	if (u >= x) a=x; else b=x; 
            	v=w; fv=fw; dv=dw;
            	// MOV3(v,fv,dv, w,fw,dw)
            	w=x; fw=fx; dw=dx;
            	//MOV3(w,fw,dw, x,fx,dx) 
            	x=u; fx=fu; dx=du;
            	// MOV3(x,fx,dx, u,fu,du)
            } else {
            	if (u < x) a=u; else b=u; 
            	if (fu <= fw || w == x) {
            		v=w; fv=fw; dv=dw;
            		// MOV3(v,fv,dv, w,fw,dw)
            		w=u; fw=fu; dw=du;
            		// MOV3(w,fw,dw, u,fu,du) 
            	} else if (fu < fv || v == x || v==w) {
            		v=u; fv=fu; dv=du;	
        			// MOV3(v,fv,dv, u,fu,du)
            	}	 
            }
    	} // for iter
        throw new OptimizationException("Too many iterations in brent");
    }

    /**
     * GOLD is the default ratio by which successive intervals are magnified in mbrak()
     */
    private static final double GOLD=0.5*(Math.sqrt(5.0)+1.0);
    /**
     * GLIMIT is the maximum magnification allowed for a parabolic-fit step of mbrak().
     */
    public static double GLIMIT=100.0;
    private static final double MBRAK_TINY=1.0e-20;

    /**
     * Initial bracketing routine; given starting guesses ax and bx, brackets the minimum.
     * 
     * Based on Numerical Recipes 10.1.
     *
     * @param ax interval [left] endpoint
     * @param bx interval [right] endpoint
     * @param func function for which minimization is needed
     * @return [a,b, c, func(a), func(b), func(c)] with midpoint a&lt;c&lt;b or a&gt;b&gt;c such that f(b)&lt; min(f(a), f(c))
     */
    public static double[] mnbrak(double ax, double bx, DoubleFunction<Double> func){
        double[] retval=new double[6];

        // Given a function func, and given distinct initial points ax and bx, this routine searches in
        // the downhill direction (defined by the function as evaluated at the initial points) and returns
        // new points ax, bx, cx that bracket a minimum of the function. 
        double fa=func.apply(ax);
        double fb=func.apply(bx);
        if (fb > fa) { // Switch roles of a and b so that we can go downhill in the direction from a to b. 
            double tmp=fa;
            fa=fb; fb=tmp;
            tmp=ax;
            ax=bx; bx=tmp;
        }
        double cx=(bx)+GOLD*(bx-ax); // First guess for c.
        double fc=func.apply(cx);
        while (fb > fc) { // Keep returning here until we bracket.
            // Compute u by parabolic extrapolation from a; b; c. 
            // TINY is used to prevent any possible division by zero.        
            double r=(bx-ax)*(fb-fc); 
            double q=(bx-cx)*(fb-fa);

            double u=bx-((bx-cx)*q-(bx-ax)*r)/
                (2.0*sign(Math.max(Math.abs(q-r),MBRAK_TINY),q-r));
            double ulim=(bx)+GLIMIT*(cx-bx); 
            double fu;
            // We won't go farther than this. Test various possibilities:
            if ((bx-u)*(u-cx) > 0.0) { // Parabolic u is between b and c: try it.
                //System.out.println("#**FM.mbrak case 1");
                fu=func.apply(u);
                if (fu < fc) { // Got a minimum between b and c.
                    retval[0]=bx; retval[3]=fb;
                    retval[1]=u;  retval[4]=fu;
                    retval[2]=cx; retval[5]=fc;
                    return retval;
                } else if (fu > fb) { // Got a minimum between between a and u.
                    retval[0]=ax; retval[3]=fa;
                    retval[1]=bx; retval[4]=fb;
                    retval[2]=u;  retval[5]=fu;
                    return retval;
                }
                u=(cx)+GOLD*(cx-bx); // Parabolic fit was no use. Use default magnification. 
                fu=func.apply(u);
            } else if ((cx-u)*(u-ulim) > 0.0) { // Parabolic  fit is between c and its allowed limit. 
                //System.out.println("#**FM.mbrak case 2");
                fu=func.apply(u);
                if (fu < fc) {
                    bx=cx; cx=u; u=cx+GOLD*(cx-bx);
                    fb=fc; fc=fu; fu=func.apply(u);
                }
            } else if ((u-ulim)*(ulim-cx) >= 0.0) { // Limit parabolic u to maximum allowed value. 
                //System.out.println("#**FM.mbrak case 3");
                u=ulim;
                fu=func.apply(u);
            } else {  // Reject parabolic u, use default magnification.
                //System.out.println("#**FM.mbrak case 4");
                u=(cx)+GOLD*(cx-bx);
                fu=func.apply(u);
            }
            ax=bx; bx=cx; cx=u;
            fa=fb; fb=fc; fc=fu; // Eliminate oldest point and continue.
        }
        retval[0]=ax; retval[3]=fa;
        retval[1]=bx; retval[4]=fb;
        retval[2]=cx;  retval[5]=fc;
        return retval;
    }


    
    
    /**
     * Ensures sufficient decrease in function value. 
     */
    public static double LNSRCH_ALF=1.0e-4; 

    /**
     * Given a point xold[1..n], and a direction p[1..n], finds a new point x[1..n] 
     * along the direction p from 
     * xold where the function func has decreased <q>sufficiently.</q>
     *
     * @param xold starting point
     * @param fold function value at starting point
     * @param g gradient at starting point
     * @param p direction for line search, usually the Newton direction
     * @param x new point (returned)
     * @param f new function value (returned)
     * @param stpmax limits the length of the steps so that you do not try to 
     *   evaluate the function in regions where it is undefined or subject to overflow
     * @param func the function
     *
     * @return false on a normal exit. It is true when x is too close to xold. 
     *   In a minimization algorithm, this usually signals convergence and can 
     *   be ignored. However, in a zero-finding algorithm the calling program should check whether the 
     *   convergence is spurious.
     */ 
    public static boolean lnsrch(double []xold, double fold, double[] g, double[] p, double[] x, 
        double f[], double stpmax, Function<double[], Double> func) 
    {
    	
    	double myminf = fold;
    	double[] myminx = xold;
    	
    	
        //inti; 
        //floata,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp, 
        //test,tmplam; 
        //*check=0; 
        double sum=0.0;
        
        double max_p = 0.0;
        int n = xold.length; // vector length
        for(int i=0;i<n;i++) 
            sum+=p[i]*p[i]; 
        sum=Math.sqrt(sum); 
        if(sum>stpmax) 
            for(int i=0;i<n;i++)
                p[i]*=stpmax/sum; // Scale if attempted step is too big. 
        double slope = 0.0;
        for (int i=0;i<n;i++) 
            slope+=g[i]*p[i]; 

        if(slope>=0.0)
        {
        	if (REPORT_UNUSUAL)
        		System.out.println("#**FM.lnsrch slope "+slope+"\tg "+Arrays.toString(g)+"\tp "+Arrays.toString(p));
        	
        	System.arraycopy(xold, 0, x, 0, xold.length);
        	f[0] = fold;
        	
        	return true;
        	
            // throw new OptimizationException("Roundoff problem in lnsrch."); 
        }
        // compute lambda_min
        double test=0.0;
        for (int i=0; i<n; i++)
        { 
            double temp=Math.abs(p[i])/Math.max(Math.abs(xold[i]),1.0); 
            if(temp>test)test=temp; 
        } 
        double alamin=LNSRCH_TOLX/test; 
        double alam=1.0; // Always try full Newton step first. 
        double alam2=Double.NaN;
        double f2=Double.NaN;
        
//        System.out.println("#**FM.lnsrch sum "+sum+"\tscale (if <1 rescaled direction) "+(stpmax/sum)+"\ttest "+test+"\talamin "+alamin+"\tstpmax "+stpmax);        
        // avoid boundary with function value infinity
        	
        
        
        for(int iter=0; true; iter++)
        { // Start of iteration loop. 
//        	do
//        	{
	            for(int i=0;i<n;i++)
	                x[i]=xold[i]+alam*p[i]; 
            
            
//            try
//            {
//            	f[0]=1e99;
            
	            f[0]=func.apply(x); 
//	            System.out.println("#**FM.lnsrch iter "+iter+"\tf "+f[0]+"\talam "+alam+"\t|p| "+euclideanNorm(p));
//	            if (iter % 5==0)
//	            {
//	            	for (int i=0; i<n; i++)
//	            	{
//	            		System.out.println("#**FM.lnsrch iter "+iter+"\talam "+alam+"\ti "+i+"\tx "+x[i]+"\tp "+p[i]+"\txo "+xold[i]);
//	            	}
////	            	System.out.println("#**FM.lnsrch iter "+iter+"\talam "+alam+"\tp "+Arrays.toString(p)+"\t; x "+Arrays.toString(x));
//	            }
//	            
////	            if (!Double.isFinite(f[0]) && iter==0)
////	            {
////	            	alam *= 0.5; // too big
////	            	System.out.println("#**FM.lnsrch reducing alam "+alam);
////	            }
//        	} while (iter==0 && !Double.isFinite(f[0]) && alam>=alamin);
            	
        	if (f[0]<myminf)
        	{
        		myminf = f[0];
        		myminx = x.clone();
        	}
            
            double tmplam=Double.NaN;
            
            if(alam<alamin)
            { // Convergence on Delta x. For zero finding, the calling program should verify the convergence. 
                if (myminf<f[0])
                {
                	if (REPORT_UNUSUAL)
                		System.out.println("#**FM.lnsrch skip/+\t"+f[0]+"->"+myminf);
                    for(int i=0;i<n;i++)
                    {
                        x[i]=myminx[i]; 
                    }
                	f[0]=myminf;
                } else
                {
	                for(int i=0;i<n;i++)
	                {
	                    x[i]=xold[i]; 
	                }
                }
                return true;
            }
            else if(f[0]<=fold+LNSRCH_ALF*alam*slope)
            {
                {
                    if (myminf<f[0])
                    {
                    	if (REPORT_UNUSUAL)
                    		System.out.println("#**FM.lnsrch skip/-\t"+f[0]+"->"+myminf);
                        for(int i=0;i<n;i++)
                        {
                            x[i]=myminx[i]; 
                        }
                    	f[0]=myminf;
                    } 
                	return false ; // Sufficient function decrease. 
                }
            }
            else
            { // Backtrack.
                if(alam==1.0) 
                    tmplam=-slope/(2.0*(f[0]-fold-slope)); // First time. 
                else
                { // Subsequent backtracks. 
                    double rhs1=f[0]-fold-alam*slope; 
                    double rhs2=f2-fold-alam2*slope; 
                    double a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2); 
                    double b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2); 
                    if(a==0.0)
                        tmplam=-slope/(2.0*b); 
                    else
                    { 
                        double disc=b*b-3.0*a*slope; 
                        if(disc<0.0)
                            tmplam=0.5*alam; 
                        else if(b<=0.0)
                            tmplam=(-b+Math.sqrt(disc))/(3.0*a); 
                        else tmplam=-slope/(b+Math.sqrt(disc)); 
                    } 
                    if(tmplam>0.5*alam) 
                        tmplam=0.5*alam; // lambda  \le 0.5 \lambda_1
                } 
            } 
            alam2=alam; 
            f2=f[0]; 
            alam=Math.max(tmplam,0.1*alam); // \lambda \ge 0.1 \lambda_1 
        } // Try again. 
    } 
    

    /**
     * Maximum allowed number of iterations (dfpmin) 
     */
    public static int DFP_ITMAX=200;


    /**
     * Scaled maximum step length allowed in line searches (dfpmin).
     */
    public static double DFP_STPMX=100.0; //5.0; // Numerical Recipes: 100.0;  
    
    /**
     * Given a starting point, the Broyden-Fletcher-Goldfarb-Shanno
     * variant of Davidon-Fletcher-Powell minimization is performed
     * on a function, using 
     * its gradient.  
     * The routine lnsrch is called to perform approximate line minimizations. 
     *
     * Based on Numerical Recipes ch. 10.7
     *
     * @param p starting point, also the location of the minimum on return
     * @param gtol the convergence requirement on zeroing the gradient
     * @param func the function to be minimized
     * 
     * @return the value of the minimum, and p is updated to show its location
    */
    public static double dfpmin(double[] p, double gtol, Function<double[], Double> func, Function<double[], double[]>  gradient)
    {
    	return dfpmin(p, gtol, DFP_ITMAX, func, gradient, null);

    }
    
    public static boolean DFP_BFGS_UPDATE = true; // if false, original Davidon-Fletcher-Powell update
    private static final boolean CG_POLAK_RIBIERE = true; // if false, 
    
    private static final boolean CHATTY_DFP = false;
    
    /**
     * Given a starting point, the Broyden-Fletcher-Goldfarb-Shanno
     * variant of Davidon-Fletcher-Powell minimization is performed
     * on a function, using 
     * its gradient.  
     * The routine lnsrch is called to perform approximate line minimizations. 
     *
     * Based on Numerical Recipes ch. 10.7
     *
     * @param p starting point, also the location of the minimum on return
     * @param gtol the convergence requirement on zeroing the gradient
     * @param func the function to be minimized
     * @param iterations tracks the intermediate function values: one element added  per iteration
     * 
     * @return the value of the minimum, and p is updated to show its location
    */
    public static double dfpmin(double[] p, double gtol, int dfp_itmax, Function<double[], Double> func, Function<double[], double[]>  gradient, List<Double> iterations)
    {
        int n=p.length;
        
        double[] fret = new double[1];

        double fp=func.apply(p); // Calculate starting function value and gradient
        double[] g=gradient.apply(p);
        
        if (iterations!=null) iterations.add(fp);
        
        double[][] hessin=new double[n][n];
        double[] xi=new double[n]; 
        double sum=0.0;
        double max_p = 0.0; // to avoid extreme steps with skewed gradients
        for(int i=0;i<n;i++)
        { // and initialize the inverse Hessian to the unit matrix.
            for(int j=0;j<n;j++)
                hessin[i][j]=0.0; 
            hessin[i][i]=1.0; 
            xi[i]=-g[i]; // Initial line direction. 
            {
            	double c = p[i]*p[i];
//            	System.out.println("#**FM.dfpmin i "+i+"\tpi "+p[i]+"\tp^2i "+c+"\tgi "+g[i]);
            }
            sum+=p[i]*p[i]; 
            max_p = Double.max(max_p, Math.abs(p[i]));
        }
        double stpmax= DFP_STPMX*Math.max(Math.sqrt(sum),(double)n); // versions 23.xx < 23.0820 DFP_STPMX*Math.max(Math.sqrt(sum), 1.0);  // Numerical Recipes: DFP_STPMX*Math.max(Math.sqrt(sum),(double)n);   
//        // DEBUG
//        double stpmaxmax =  DFP_STPMX*Math.max(max_p, 1.0); 
//        System.out.println("#*FM.dfpmin stpmax "+stpmax+"\t(maxmax "+stpmaxmax+")\tsum "+sum+"\tn "+n+"\titmax "+dfp_itmax+"\tNRstpmax "+DFP_STPMX*Math.max(Math.sqrt(sum),(double)n)
//        		+"\tsqrsum "+Math.sqrt(sum)+"\tmaxp "+max_p);
        
        //stpmax = stpmaxmax;
        
        double[] pnew=new double[n]; 
        double[] dg=new double[n];
        double[] hdg=new double[n]; 
        int its;
        for(its=1;its<=dfp_itmax;its++)
        { // Main loop over the iterations. 
        	
//        	try
//        	{
        		lnsrch(p,fp,g,xi,pnew,fret,stpmax,func); 
//        	} catch (Throwable t)
//        	{
//        		System.out.println("#**FM.dfpmin numerr its "+its+"\tstpmax "+stpmax+"\tsum "+sum);
//        		throw t;
//        	}

            // The new function evaluation occurs in lnsrch; save the function value in fp for the 
            // next linesearch. It is usually safe to ignore the return value. 
            fp=fret[0]; 
            if (iterations!=null) iterations.add(fp);
            
            for(int i=0;i<n;i++)
            { 
                xi[i]=pnew[i]-p[i]; // Update the line direction, 
                p[i]=pnew[i]; // and the current point. 
            } 
            double test=0.0; // Test for convergence on Delta x. 
            for(int i=0;i<n;i++)
            { 
                double temp=Math.abs(xi[i])/Math.max(Math.abs(p[i]),1.0); 
                if(temp>test)
                    test=temp; 
            } 
            if(test<DFP_TOLX)
            { 
            	if (CHATTY_DFP)
            		System.out.println("#**FM.dfp "+its+"\tdeltax small\t"+test+"\tdone/delta"+"\t// "+DFP_TOLX); // DEBUG
                return fret[0]; 
            } 
            double test_deltax = test;
            
            for(int i=0;i<n;i++)
                dg[i]=g[i]; // Save the old gradient, 
            g=gradient.apply(p); // and get the new gradient. 
            test=0.0; // Test for convergence on zero gradient. 
            double den= Math.max(Math.abs(fret[0]),1.0);    // Numerical recipes: Math.max(fret[0],1.0); 
            for(int i=0;i<n;i++)
            { 
                double temp=Math.abs(g[i])*Math.max(Math.abs(p[i]),1.0)/den; 
                if(temp>test)
                    test=temp; 
            } 
            if(test<gtol)
            { 
            	if (CHATTY_DFP)
            		System.out.println("#**FM.dfp "+its+"\tgrad "+test+"\tdone/gradient"); // DEBUG
                return fret[0]; 
            } 
            else {
            	
//            	System.out.println("#**FM.dfp "+its+"\tgrad "+test+"\tdeltax "+test_deltax+"\tgtol,tolx "+gtol+","+DFP_TOLX); 
            }
            
            for(int i=0;i<n;i++)
                dg[i]=g[i]-dg[i]; // Compute difference of gradients, 
            for(int i=0;i<n;i++)
            { // and difference times current matrix. 
                hdg[i]=0.0; 
                for(int j=0;j<n;j++)
                    hdg[i]+=hessin[i][j]*dg[j]; 
            } 
            double fac=0., fad=0., fae=0., sumdg=0., sumxi=0.; // Calculate dot products for the denominators.
            for(int i=0;i<n;i++)
            { 
                fac+=dg[i]*xi[i]; 
                fae+=dg[i]*hdg[i]; 
                sumdg+=dg[i]*dg[i]; //      not Math.sqrt(dg[i]); 
                sumxi+=xi[i]*xi[i]; //      not Math.sqrt(xi[i]); 
            } 
            //
            // DFP updating formula (10.7.8)
            //		H = H + xi**xi/fac - hdg**hdg/fae
            //
            // vector for BFGS
            //		u = xi / fac - hdg / fae
            // BFGS updating formula (10.7.9)	
            //		H = H + xi**xi/fac - hdg**hdg/fac - fae u**u
            //
            // DFP: 
            //		H[i][j] = H[i][j] + xi[i] * xi[j] / fac + hdg[i]*hdg[j]/fae  
            // BFGS:
            //		H[i][j] = H[i][j] + xi[i]*xi[j]/fac - hdg[i]*hdg[j]/fae + fae*(xi[i]/fac - hdg[i]/fae)*(xi[j]/fac-hdg[j]/fae)
            //				= H[i][j] + xi[i]*xi[j]/fac 
            if(fac>Math.sqrt(EPS*sumdg*sumxi))
            { // Skip update if fac not sufficiently positive.
                fac=1.0/fac; 
                fad=1.0/fae; 
                
                if (DFP_BFGS_UPDATE)
                {
	                
	                // The vector that makes BFGS different from DFP: 
	                for(int i=0;i<n;i++)
	                    dg[i]=fac*xi[i]-fad*hdg[i]; // == u from (10.7.10)
	                for(int i=0;i<n;i++)
	                { // The BFGS updating formula: 
	                    for(int j=i;j<n;j++)
	                    { 
	                        hessin[i][j]+=fac*xi[i]*xi[j] -fad*hdg[i]*hdg[j]+  fae*dg[i]*dg[j]; 
	                        hessin[j][i]=hessin[i][j]; 
	                    } 
	                } 
                } else
                {
	                for(int i=0;i<n;i++)
	                { // The DFP updating formula: 
	                    for(int j=i;j<n;j++)
	                    { 
	                        hessin[i][j]+=fac*xi[i]*xi[j] -fad*hdg[i]*hdg[j]; 
	                        hessin[j][i]=hessin[i][j]; 
	                    } 
	                } 
                }
            } 
            for(int i=0;i<n;i++)
            { // Now calculate the next direction to go, 
                xi[i]=0.0; 
                for(int j=0;j<n;j++)
                    xi[i]-=hessin[i][j]*g[j]; 
            } 
            
            if (Thread.currentThread().isInterrupted()) // keeps interrupt status
            {
                return Double.NaN;
            }
            
        } // and go back for another iteration. 
        if (iterations == null) 
        	if (REPORT_UNUSUAL)
        		System.out.println("#**FM.dfpmin: Too many iterations in dfpmin"); 
        
    	if (CHATTY_DFP)
    		System.out.println("#**FM.dfp "+its+"\tdone/iter"); // DEBUG

        
        return fret[0];
    }    
    /**
     * calls powell with unit directions
     */
    public static double powell(double[] p, double ftol, Function<double[], Double> func)
    {
//        int n=p.length;
//        double [][] xi=new double[n][n];
//        for (int i=0; i<n; i++){
//            Arrays.fill(xi[i],0.0);
//            xi[i][i]=1.0;
//        }
        return powell(p, ftol, POWELL_ITMAX, func, null);
    }
    
    /**
     * calls powell with unit directions
     */
    public static double powell(double[] p, double ftol, int powell_itmax, Function<double[], Double> func, List<Double> iterations)
    {
        int n=p.length;
        double [][] xi=new double[n][n];
        for (int i=0; i<n; i++){
            Arrays.fill(xi[i],0.0);
            xi[i][i]=1.0;
        }
        return powell(p, xi, ftol, powell_itmax, func, iterations);
    }
    
    private final static double POWELL_TINY=1.0e-25; // A small number.
    public static int POWELL_ITMAX=200; // Maximum allowed iterations.
    /**
     * Minimization of a function func of n variables. Input consists of an initial starting point
     * p[1..n]; an initial matrix xi[1..n][1..n], whose rows contain the initial set of directions
     * (usually the n unit vectors); and ftol, the fractional tolerance in the function value
     * such that failure to decrease by more than this amount on one iteration signals doneness. 
     * On output, p is set to the best point found, xi is the then-current direction set.
     * The routine linmin is used.
     * 
     * Based on Numerical Recipes ch. 10.5.
     *
     * @param p initial point
     * @param xi initial set of directions
     * @param ftol tolerance in function value change
     * @param func function to be minimized
     * @return minimum function value
     */ 
    public static double powell(double p[], double[][] xi, double ftol, Function<double[], Double> func)
    {
    	return powell(p, xi, ftol, POWELL_ITMAX, func, null);
    }    
    
    /**
     * Minimization of a function func of n variables. Input consists of an initial starting point
     * p[1..n]; an initial matrix xi[1..n][1..n], whose rows contain the initial set of directions
     * (usually the n unit vectors); and ftol, the fractional tolerance in the function value
     * such that failure to decrease by more than this amount on one iteration signals doneness. 
     * On output, p is set to the best point found, xi is the then-current direction set.
     * The routine linmin is used.
     * 
     * Based on Numerical Recipes ch. 10.5.
     *
     * @param p initial point
     * @param xi initial set of directions
     * @param ftol tolerance in function value change
     * @param powell_itmax maximum number of iterations for outer loop
     * @param func function to be minimized
     * @param iterations if not null, intermediate function values are added to the list at each iteration
     * @return minimum function value
     */ 
    public static double powell(double p[], double[][] xi, double ftol, int powell_itmax, Function<double[], Double> func, List<Double> iterations)
    {
        int n=p.length;
        double[] pt=new double[n];
        double[] ptt=new double[n];
        double[] xit=new double[n];
        double fret=func.apply(p);
        double fptt;
        System.arraycopy(p, 0, pt, 0, n);
        //for (int j=0;j<n;j++) pt[j]=p[j];  // Save the initial point.
        for (int iter=1;iter <=powell_itmax;iter++) 
        {    	
            double fp=fret;
            int ibig=0;    
            double del=0.0; // Will be the biggest function decrease.
            for (int i=0;i<n;i++) {  // In each iteration, loop over all directions in the set.
                System.arraycopy(xi[i], 0, xit, 0, n);
                //for (int j=0;j<n;j++) xit[j]=xi[i][j]; // Copy the direction,
                fptt=fret;
                fret=linmin(p,xit,func); // minimize along it,
                if (fptt-(fret) > del) { // and record it if it is the largest decrease so far. 
                    del=fptt-(fret);
                    ibig=i;
                }
                // DEBUG
                System.out.println("#**FM.powell iter:dir "+iter+":"+i+"\tf "+fret+"\tdel "+(fptt-fret));
            }
            
            // DEBUG
            System.out.println("#**FM.powell iter "+iter+"/"+powell_itmax+"\tf "+fret+"\tdel "+del+"\t@ibig "+ibig+"/"+n);
            
            if (iterations!=null) iterations.add(fret); // record this step
            
            if (2.0*(fp-(fret)) <= ftol*(Math.abs(fp)+Math.abs(fret))+POWELL_TINY) {
                return fret;
            }
            for (int j=0;j<n;j++) { // Construct the extrapolated point and the average direction moved. Save the old starting point.
                ptt[j]=2.0*p[j]-pt[j];
                xit[j]=p[j]-pt[j];
                pt[j]=p[j];
            }
            fptt=func.apply(ptt); // Function value at extrapolated point.
            if (fptt < fp) {
                double t=2.0*(fp-2.0*fret+fptt)*SQR(fp-(fret)-del)-del*SQR(fp-fptt); // not 2.0*(fp-2.0*fret+fptt)*Math.sqrt(fp-(fret)-del)-del*Math.sqrt(fp-fptt);
                if (t < 0.0) {
                    fret=linmin(p,xit,func); // Move to the minimum of the new direction, and save the new direction. 
                    System.arraycopy(xi[n-1], 0, xi[ibig], 0, n);
                    System.arraycopy(xit, 0, xi[n-1], 0, n);
                    //for (int j=0;j<n;j++) {
                    //    xi[ibig][j]=xi[n-1][j];
                    //    xi[n-1][j]=xit[j];
                    //}
                }
            }
            if (Thread.currentThread().isInterrupted()) // keeps interrupt status
            {
                return Double.NaN;
            }
            
//            System.out.println("#*FM.powell iter "+iter+"\tfp "+fp+"\tfret "+fret);
        } // Back for another iteration.
        if (iterations == null)
        	if (REPORT_UNUSUAL)
        		System.out.println("#**FM.powell Too many iterations in powell");

    	return fret;
    }
    
    public static int FRPRMN_ITMAX = 200;
    
    public static double frprmn(double p[], double ftol, Function<double[], Double> func, Function<double[], double[]>  dfunc)
    {
    	return frprmn(p, ftol, FRPRMN_ITMAX, func, dfunc, null);
    }
    
    
    /**
     * Conjugate gradient method for minimization. 
     * Given a starting point p[], 
     * Fletcher-Reeves-Polak-Ribiere minimization is performed on a function func, 
     * using its gradient as calculated by a routine dfunc. 
     * The convergence tolerance on the function value is input as ftol. 
     * Returned quantities are p (the location of the minimum), 
     * fret (the minimum value of the function). 
     * 
     * Based on NumericaL Recipes 10.6.
     * 
     * @param p starting point
     * @param ftol tolerance on function value 
     * @param frprmn_itmax maximum number of iterations
     * @param func function to be minimized
     * @param dfunc gradient
     * @param iterations records the successive minimums in the iterations (null ok)
     * @return minimum function value
     */
    public static double frprmn(double p[], double ftol, int frprmn_itmax, Function<double[], Double> func, Function<double[], double[]>  dfunc, List<Double> iterations)
    {
    	// Initializations.
    	int n = p.length;
    	double[] g=new double[n]; 
    	double[] h=new double[n]; 
    	double fp=func.apply(p); 
    	double[] xi = dfunc.apply(p);
    	for (int j=0;j<n;j++) 
    	{ 
    		g[j] = -xi[j]; xi[j]=h[j]=g[j]; 
    	}
    	double fret = fp;
    	for  (int its=1;its<=frprmn_itmax;its++) {
			// Loop over iterations.
			
			fret = dlinmin(p, xi, func, dfunc); // fret = linmin(p,xi,func);
            if (iterations!=null) iterations.add(fret); // record this step
			// Next statement is the normal return: 
			if (2.0*Math.abs(fret-fp) <= ftol*(Math.abs(fret)+Math.abs(fp)+EPS)) 
			{
				return fret;
			}
			fp = fret;
			xi = dfunc.apply(p);
			double dgg=0.0;
			double gg=0.0;
			for (int j=0;j<n;j++) 
			{
				gg += g[j]*g[j];
				if (CG_POLAK_RIBIERE)
				{
					dgg += (xi[j]+g[j])*xi[j];
				} else
				{
					// Fletcher-Reeves
					dgg += xi[j]*xi[j]; 
				}
			}
			if (gg == 0.0) { // unlikely case of gradient == 0
				return fret;
			}
			double gam=dgg/gg;
			for (int j=0;j<n;j++) {
				g[j] = -xi[j];
				xi[j]=h[j]=g[j]+gam*h[j]; 
			}
		} // iterations
        if (iterations == null)
        	if (REPORT_UNUSUAL)
        		System.out.println("#**FM.frprmn Too many iterations in frprmn");
        return fret;
    }
    
    
    
    /**
     * initial bracketing values for line minimization within powell()
     */
    public static double POWELL_BRACKET_A=0.0;
    public static double POWELL_BRACKET_B=1.0;
    public static final double POWELL_TOL=2e-4;
    
    /**
     * Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and
     * resets p to where the function func(p) takes on a minimum along the direction xi from p,
     * and replaces xi by the actual vector displacement that p was moved. 
     * @return the value of func at the returned location p. 
     * This is actually all accomplished by calling the routines mnbrak and brent.
     * 
     * Based on Numerical Recipes 10.5. 
     */
    public static double linmin(double[] p, double[] xi, Function<double[], Double> func)
    {
        double ax=POWELL_BRACKET_A;
        double xx=POWELL_BRACKET_B;
        class PowellDirection implements DoubleFunction<Double>
        {
            PowellDirection()
            {
                xt=new double[p.length];
            }
            
            private final double[] xt; // reused 
            
            @Override
            public Double apply(double x)
            {
                for (int i=0; i<p.length; i++)
                {
                    xt[i]=p[i]+x*xi[i];
                    //System.out.println("FM.PD.e "+x+"\t"+i+"\tp "+p[i]+"\txi "+xi[i]+"\txt "+xt[i]);
                }
                return func.apply(xt);
            }
        }
            
        DoubleFunction<Double> f1dim=new PowellDirection();
        double[] brak=mnbrak(ax, xx, f1dim);
        double[] z=brent(brak[0], brak[1], brak[2], f1dim, POWELL_TOL);
        double fret=z[1];
        double xmin=z[0];
        int n=p.length;
        for (int j=0; j<n; j++)
        {
            xi[j] *= xmin;
            p[j] += xi[j];
        }
        return fret;
    }
    
    /**
     * Given an n-dimensional point p[] and an n-dimensional direction xi[], 
     * moves and resets p to where the function func(p) takes on a minimum along the 
     * direction xi from p, and replaces xi by the actual vector displacement that p was moved. 
     * Also returns as fret the value of func at the returned location p. 
     * This is actually all accomplished by calling the routines mnbrak and dbrent.
     * 
     * Based on Numerical Recipes 10.6. 
     * 
     * @param p
     * @param xi
     * @param func
     * @param dfunc
     * @return
     */
    public static double dlinmin(double[] p, double[] xi, Function<double[], Double> func, Function<double[], double[]> dfunc)
    {
        double ax=POWELL_BRACKET_A;
        double xx=POWELL_BRACKET_B;
        final int n=p.length;
        class F implements DoubleFunction<Double>
        {
            F()
            {
                xt=new double[n];
            }
            
            private final double[] xt; // reused 
            
            @Override
            public Double apply(double x)
            {
                for (int i=0; i<n; i++)
                {
                    xt[i]=p[i]+x*xi[i];
//                    System.out.println("FM.PD.e "+x+"\t"+i+"\tp "+p[i]+"\txi "+xi[i]+"\txt "+xt[i]);
                }
                return func.apply(xt);
            }
        }    	
        
        class DF implements DoubleFunction<Double>
        {
        	DF()
        	{
                xt=new double[n];
        	}
        	private final double[] xt; // reused 
            @Override
            public Double apply(double x)
            {
                for (int i=0; i<n; i++)
                {
                    xt[i]=p[i]+x*xi[i];
//                    System.out.println("FM.PD.e "+x+"\t"+i+"\tp "+p[i]+"\txi "+xi[i]+"\txt "+xt[i]);
                }
                double[] df = dfunc.apply(xt);
                double dfdx = 0.0;
                for (int j=0; j<n; j++)
                {
                	dfdx += df[j]*xi[j];
                }
                return dfdx;
            }
        }
        DoubleFunction<Double> f1dim=new F();
        DoubleFunction<Double> df1dim = new DF();
        double[] brak=mnbrak(ax, xx, f1dim);
        double[] z=dbrent(brak[0], brak[1], brak[2], f1dim, df1dim, POWELL_TOL);      
        double fret=z[1];
        double xmin=z[0];
        
        for (int j=0; j<n; j++)
        {
            xi[j] *= xmin;
            p[j] += xi[j];
        }
        return fret;
    }
    
    
    public static final double ARMIJO_C = 0.9; //0.5; //Math.exp(Math.log(0.5)/2.0);
    public static final double ARMIJO_MUL = 0.5; //Math.exp(Math.log(0.5)/2.0);
    public static final double ARMIJO_TMIN = 32.0;
    
    public static double armijo(double[] x, double fx, double[] df, Function<double[], Double> func)
    {
    	double[] x0 = new double[x.length];
    	for (int i=0; i<x0.length; i++) x0[i]=x[i];
    	double f0 = fx;
    	
    	assert (x.length == df.length);
    	
    	double[] p = new double[df.length]; // unit direction opposite to df0
    	double df_sumsq = 0.0;
    	for (int i=0; i<df.length; i++)
    	{
    		double d = df[i];
    		p[i] = -d;
    		df_sumsq += d*d;
    	}
    	double df_length = Math.sqrt(df_sumsq); // euclidean norm
    	for (int i=0; i<p.length; i++)
    	{
    		p[i] /= df_length;
    	}
    	double m = -df_length; // dot product between gradient and p 
    	
    	double t = -ARMIJO_C*m;
    	if (t>ARMIJO_TMIN) t=ARMIJO_TMIN;

    	int iter = 0;
    	double alpha_step = 16.0/ARMIJO_MUL;
    	double delta_f;
    	double threshold;
    	do {
    		alpha_step *= ARMIJO_MUL;
    		for (int i=0; i<x.length; i++ )
    		{
    			x[i] = x0[i] + alpha_step*p[i];
    		}
    		fx = func.apply(x);
    		
    		delta_f = f0-fx;
    		threshold = alpha_step * t;

    		iter++;
//    		System.out.println("#**FM.arm iter "+iter+"\tdelta "+delta_f+"\tthreshold "+threshold+"\tdf "+df_length);
    	} while (delta_f < threshold  && delta_f != 0.0);
    	return fx;
    }
    
    public static boolean GD_ARMIJO = true; // if false, then line minimization
    
    public static double gradientDescent(double[] x, double gtol, int gd_itmax, Function<double[], Double> func, Function<double[], double[]>  gradient, List<Double> iterations)
    {
    	double[] x0 = new double[x.length];
    	for (int i=0; i<x0.length; i++)
    		x0[i]=x[i];
    	double f0 = func.apply(x0);
		double[] df0 = gradient.apply(x0);
		
    	for (int iter=0; iter<gd_itmax; iter++)
    	{
    		double f1;
    		if (GD_ARMIJO)
    		{
    			f1 = armijo(x, f0, df0, func); // sets x 
    		} else
    		{ 
    			double df_sumsq=0.0;
    			for (int j=0; j<df0.length; j++) 
    			{
    				double d = df0[j];
    				df0[j]=-d;
    	    		df_sumsq += d*d;
    			}
    			double df_len = Math.sqrt(df_sumsq);
    			for (int j=0; j<df0.length; j++)
    			{
    				df0[j]/=df_len;
    			}
    			f1 = dlinmin(x, df0, func, gradient); // sets x 
    			//f1 = linmin(x, df0, func);
    		}
            if (iterations!=null) iterations.add(f1);
            
            double[] df = gradient.apply(x);
            
            double dlen = euclideanNorm(df);
            double fdiff = f0-f1;
            double rel_change = Math.abs(fdiff/f0);
            double rel_dlen = dlen/f1;
            
//            System.out.println("#**FM.gD iter "+iter+"\tf1 "+f1+"\tfdiff "+fdiff+"\tdlen "+dlen+"\trellen "+rel_dlen+"\trelchg "+rel_change);
            
            if (dlen < Math.abs(gtol*f1) || rel_change<gtol)
            	return f1;
            
            f0 = f1;
            df0 = df;
    	}
    	return f0; 
    }
    
    
    
    
    
    public static double GRADIENT_DIFF_EPS = 1.0/(1<<17); // around cubic root of machine precision  
    /**
     * Straightforward one-point estimation of partial derivatives; no guarantees 
     * for accuracy. (Don't do numerical derivation ...) 
     * 
     * @param func function for which gradient is calculated
     * @param x at this point 
     * @return vector of partial derivatives; same length as x 
     */
    public static double[] numericalGradient(Function<double[],Double> func, double[] x)
    {
		double[] D = new double[x.length]; // return value
		// we leave the function call at x for the last, so that
		// return state yields the function evaluation at x 
		for (int p=0; p<x.length; p++)
		{
			double θ = x[p];
			double h = Math.abs(θ*GRADIENT_DIFF_EPS);
			x[p] = θ+h;
			double fd = func.apply(x);
			D[p]=fd;
			
			x[p] = θ;
		}
		double fx = func.apply(x);
		for (int p=0; p<x.length; p++)
		{
			double θ = x[p];
			double h = Math.abs(θ*GRADIENT_DIFF_EPS);
			double fd = D[p];
			double delta = fd-fx;
			double dfdθ = delta/h;
			D[p] = dfdθ;
		}
		return D;
    }

    /**
     * Straightforward two-point estimation of partial derivatives; no guarantees 
     * for accuracy. (Don't do numerical derivation ...) 
     * 
     * @param func function for which gradient is calculated
     * @param x at this point 
     * @return vector of partial derivatives; same length as x 
     */
    public static double[] numericalGradientTwoPoint(Function<double[],Double> func, double[] x)
    {
		double[] D = new double[x.length]; // return value
		
		for (int p=0; p<x.length; p++)
		{
			double θ = x[p];
			double h = Math.abs(θ*GRADIENT_DIFF_EPS);
			x[p] = θ+h;
			double fd = func.apply(x);
			D[p]=fd;
			
			x[p] = θ;
		}
		for (int p=0; p<x.length; p++)
		{
			double θ = x[p];
			double h = Math.abs(θ*GRADIENT_DIFF_EPS);
			x[p] = θ-h;
			double fx = func.apply(x);
			double fd = D[p];
			double delta = fd-fx;
			double dfdθ = delta/(2.0*h);
			D[p] = dfdθ;
			
			x[p] = θ;
		}
		func.apply(x);
		return D;
    }    
    
//    
    /**
     * Exception thrown when too many iteration in one of the routines.
     */
    public static class OptimizationException extends ArithmeticException
    {
        private OptimizationException(String message){
            super(message);
        }
    }


    
}
