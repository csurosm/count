/*
 * FunctionMinimization.java
 *
 * Created on October 26, 2004, 8:22 PM
 */

package ca.umontreal.iro.matek;

import java.util.Arrays;

/**
 *
 * This class implements a number of numerical optimization methods 
 * for function minimization in one or multiple dimensions.
 * 
 * @author  csuros
 */
public class FunctionMinimization {
    

    private FunctionMinimization(){}
    
    /**
     * GOLD is the default ratio by which successive intervals are magnified in mbrak()
     */
    private static final double GOLD=0.5*(Math.sqrt(5.0)+1.0);
    /**
     * GLIMIT is the maximum magnification allowed for a parabolic-fit step of mbrak().
     */
    public static double GLIMIT=100.0;
    private static final double MBRAK_TINY=1.0e-20;
    
    private static final double sign(double a, double b){
        if (b>=0.0)
            return Math.abs(a);
        else
            return -Math.abs(a);
    }
    

    /**
     * Using Brent's method, find the root of a function func known to lie between x1 and x2. 
     * The root will be refined until its accuracy is tol. Specifically, the 
     * iteration stops when the root is bracketed within an interval of length tol.
     *
     * Based on NR 9.3.
     *
     * @param tol bracketing accuracy
     */ 
    public static double zbrent(OneParameterFunction func,double x1, double x2, double tol) 
    { 
        double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2; 
        double fa=func.eval(a), fb=func.eval(b),fc,p,q,r,s,tol1,xm; 
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
            fb=func.eval(b); 
        } // iter
        throw new OptimizationException("Maximum number of iterations exceeded in zbrent"); 
    }
    
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
    public static double[] mnbrak(double ax, double bx, OneParameterFunction func){
        double[] retval=new double[6];

        // Given a function func, and given distinct initial points ax and bx, this routine searches in
        // the downhill direction (defined by the function as evaluated at the initial points) and returns
        // new points ax, bx, cx that bracket a minimum of the function. 
        double fa=func.eval(ax);
        double fb=func.eval(bx);
        if (fb > fa) { // Switch roles of a and b so that we can go downhill in the direction from a to b. 
            double tmp=fa;
            fa=fb; fb=tmp;
            tmp=ax;
            ax=bx; bx=tmp;
        }
        double cx=(bx)+GOLD*(bx-ax); // First guess for c.
        double fc=func.eval(cx);
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
                fu=func.eval(u);
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
                fu=func.eval(u);
            } else if ((cx-u)*(u-ulim) > 0.0) { // Parabolic  fit is between c and its allowed limit. 
                //System.out.println("#**FM.mbrak case 2");
                fu=func.eval(u);
                if (fu < fc) {
                    bx=cx; cx=u; u=cx+GOLD*(cx-bx);
                    fb=fc; fc=fu; fu=func.eval(u);
                }
            } else if ((u-ulim)*(ulim-cx) >= 0.0) { // Limit parabolic u to maximum allowed value. 
                //System.out.println("#**FM.mbrak case 3");
                u=ulim;
                fu=func.eval(u);
            } else {  // Reject parabolic u, use default magnification.
                //System.out.println("#**FM.mbrak case 4");
                u=(cx)+GOLD*(cx-bx);
                fu=func.eval(u);
            }
            ax=bx; bx=cx; cx=u;
            fa=fb; fb=fc; fc=fu; // Eliminate oldest point and continue.
        }
        retval[0]=ax; retval[3]=fa;
        retval[1]=bx; retval[4]=fb;
        retval[2]=cx;  retval[5]=fc;
        return retval;
    }
    
    private static final double R=0.5*(Math.sqrt(5.0)-1.0); // The golden ratios.
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
    public static double[] golden(double ax, double bx, double cx, OneParameterFunction func, double tol){
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
        double f1=func.eval(x1); 
        double f2=func.eval(x2);
        while (Math.abs(x3-x0) > tol*(Math.abs(x1)+Math.abs(x2))) {
            if (f2 < f1) { // One possible outcome,
                x0=x1; x1=x2; x2=R*x1+C*x3;
                f1=f2; f2=func.eval(x2);
            } else {
                x3=x2; x2=x1; x1=R*x2+C*x0;
                f2=f1; f1=func.eval(x1);
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
    
    private static final int BRENT_ITMAX=200;
    private static final double ZEPS=1.0e-10;
    // Here ITMAX is the maximum allowed number of iterations; 
    // ZEPS is a small number that protects against trying to achieve fractional accuracy for a minimum that
    // happens to be exactly zero.
    
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

    public static double[] brent(double ax, double bx, double cx, OneParameterFunction func, double tol){
        double x,w,v,fx,fw,fv,u,fu;
        double e=0.0;
        double a=(ax < cx ? ax : cx); // a and b must be in ascending order, but input abscissas need not be. 
        double b=(ax > cx ? ax : cx);
        x=w=v=bx; 
        fw=fv=fx=func.eval(x);

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
            fu=func.eval(u);
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
     * calls powell with unit directions
     */
    public static double powell(double[] p, double ftol, MultiParameterFunction func){
        int n=p.length;
        double [][] xi=new double[n][n];
        for (int i=0; i<n; i++){
            Arrays.fill(xi[i],0.0);
            xi[i][i]=1.0;
        }
        return powell(p, xi, ftol, func);
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
    public static double powell(double p[], double[][] xi, double ftol, MultiParameterFunction func){
        int n=p.length;
        double[] pt=new double[n];
        double[] ptt=new double[n];
        double[] xit=new double[n];
        double fret=func.eval(p);
        double fptt;
        System.arraycopy(p, 0, pt, 0, n);
        //for (int j=0;j<n;j++) pt[j]=p[j];  // Save the initial point.
        for (int iter=1;iter <=POWELL_ITMAX;++(iter)) {
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
            }
            if (2.0*(fp-(fret)) <= ftol*(Math.abs(fp)+Math.abs(fret))+POWELL_TINY) {
                return fret;
            }
            for (int j=0;j<n;j++) { // Construct the extrapolated point and the average direction moved. Save the old starting point.
                ptt[j]=2.0*p[j]-pt[j];
                xit[j]=p[j]-pt[j];
                pt[j]=p[j];
            }
            fptt=func.eval(ptt); // Function value at extrapolated point.
            if (fptt < fp) {
                double t=2.0*(fp-2.0*fret+fptt)*Math.sqrt(fp-(fret)-del)-del*Math.sqrt(fp-fptt);
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
        } // Back for another iteration.
        throw new OptimizationException("Too many iterations in powell");
    }
    
    /**
     * initial bracketing values for line minimization within powell()
     */
    public static double POWELL_BRACKET_A=0.0;
    public static double POWELL_BRACKET_B=1.0;
    private static final double POWELL_TOL=2e-4;
    
    /**
     * Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and
     * resets p to where the function func(p) takes on a minimum along the direction xi from p,
     * and replaces xi by the actual vector displacement that p was moved. 
     * @return the value of func at the returned location p. 
     * This is actually all accomplished by calling the routines mnbrak and brent.
     */
    private static double linmin(double[] p, double[] xi, MultiParameterFunction func){
        double ax=POWELL_BRACKET_A;
        double xx=POWELL_BRACKET_B;
        OneParameterFunction f1dim=new PowellDirection(p, xi, func);
        double[] brak=mnbrak(ax, xx, f1dim);
        double[] z=brent(brak[0], brak[1], brak[2], f1dim, POWELL_TOL);
        double fret=z[1];
        double xmin=z[0];
        int n=p.length;
        for (int j=0; j<n; j++){
            xi[j] *= xmin;
            p[j] += xi[j];
        }
        return fret;
    }
    
    static class PowellDirection implements OneParameterFunction 
    {
        PowellDirection(double[] p, double[] xi, MultiParameterFunction func){
            this.p=p;
            this.xi=xi;
            this.func=func;
            n=p.length;
            xt=new double[n];
        }
        
        private double[] p; // initial point
        private double[] xi; // direction
        private MultiParameterFunction func;
        private double[] xt;
        private int n;
        
        public double eval(double x){
            for (int i=0; i<n; i++)
            {
                xt[i]=p[i]+x*xi[i];
                //System.out.println("FM.PD.e "+x+"\t"+i+"\tp "+p[i]+"\txi "+xi[i]+"\txt "+xt[i]);
            }
            return func.eval(xt);
        }
    }
    

    /**
     * Ensures sufficient decrease in function value. 
     */
    public static double LNSRCH_ALF=1.0e-4; 
    /**
     * Convergence criterion on delta x. 
     */
    public static double LNSRCH_TOLX=1.0e-7;

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
        double f[], double stpmax, MultiParameterFunction func) 
    {
        //inti; 
        //floata,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp, 
        //test,tmplam; 
        //*check=0; 
        double sum=0.0;
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
            throw new OptimizationException("Roundoff problem in lnsrch."); 
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
        while(true)
        { // Start of iteration loop. 
            for(int i=0;i<n;i++)
                x[i]=xold[i]+alam*p[i]; 
            f[0]=func.eval(x); 
            
            double tmplam=Double.NaN;
            
            if(alam<alamin)
            { // Convergence on Delta x. For zero finding, the calling program should verify the convergence. 
                for(int i=0;i<n;i++)
                    x[i]=xold[i]; 
                return true;
            }
            else if(f[0]<=fold+LNSRCH_ALF*alam*slope)
                return false ; // Sufficient function decrease. 
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
     * Machine precision.
     */
    public static double EPS = 1e-22; 
    /**
     * Convergence criterion on x values (dfpmin)
     */
    public static double DFP_TOLX = (4.*EPS); 
    /**
     * Scaled maximum step length allowed in line searches (dfpmin).
     */
    public static double DFP_STPMX=100.0;  
            
    /*
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
    public static double dfpmin(double[] p, double gtol, DerivableMultiParameterFunction func)
    {
        int n=p.length;
        
        double[] fret = new double[1];

        double fp=func.eval(p); // Calculate starting function value and gradient
        //System.out.println("#**FM.dfpmin "+fp);
        double[] g=func.dfunc(p);
        //for (int i=0; i<g.length; i++)
        //    System.out.println("#**FM.dfpmin g"+i+"\t"+g[i]);
        
        double[][] hessin=new double[n][n];
        double[] xi=new double[n]; 
        double sum=0.0;
        for(int i=0;i<n;i++)
        { // and initialize the inverse Hessian to the unit matrix.
            for(int j=0;j<n;j++)
                hessin[i][j]=0.0; 
            hessin[i][i]=1.0; 
            xi[i]=-g[i]; // Initial line direction. 
            sum+=p[i]*p[i]; 
        }
        double stpmax=DFP_STPMX*Math.max(Math.sqrt(sum),(double)n);    
        
        double[] pnew=new double[n]; 
        double[] dg=new double[n];
        double[] hdg=new double[n]; 
        for(int its=1;its<=DFP_ITMAX;its++)
        { // Main loop over the iterations. 
            lnsrch(p,fp,g,xi,pnew,fret,stpmax,func); 
            
            // The new function evaluation occurs in lnsrch; save the function value in fp for the 
            // next linesearch. It is usually safe to ignore the return value. 
            fp=fret[0]; 
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
                return fret[0]; 
            } 
            for(int i=0;i<n;i++)
                dg[i]=g[i]; // Save the old gradient, 
            g=func.dfunc(p); // and get the new gradient. 
            test=0.0; // Test for convergence on zero gradient. 
            double den=Math.max(fret[0],1.0); 
            for(int i=0;i<n;i++)
            { 
                double temp=Math.abs(g[i])*Math.max(Math.abs(p[i]),1.0)/den; 
                if(temp>test)
                    test=temp; 
            } 
            if(test<gtol)
            { 
                return fret[0]; 
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
                sumdg+=Math.sqrt(dg[i]); 
                sumxi+=Math.sqrt(xi[i]); 
            } 
            if(fac>Math.sqrt(EPS*sumdg*sumxi))
            { // Skip update if fac not sufficiently positive.
                fac=1.0/fac; 
                fad=1.0/fae; 
                // The vector that makes BFGS different from DFP: 
                for(int i=0;i<n;i++)
                    dg[i]=fac*xi[i]-fad*hdg[i]; 
                for(int i=0;i<n;i++)
                { // The BFGS updating formula: 
                    for(int j=i;j<n;j++)
                    { 
                        hessin[i][j]+=fac*xi[i]*xi[j] -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j]; 
                        hessin[j][i]=hessin[i][j]; 
                    } 
                } 
            } 
            for(int i=0;i<n;i++)
            { // Now calculate the next direction to go, 
                xi[i]=0.0; 
                for(int j=0;j<n;j++)
                    xi[i]-=hessin[i][j]*g[j]; 
            } 
        } // and go back for another iteration. 
        throw new FunctionMinimization.OptimizationException("Too many iterations in dfpmin"); 
    }
    /**
     * Exception throws when too many iteration in one of the routines.
     */
    public static class OptimizationException extends RuntimeException 
    {
        private OptimizationException(String message){
            super(message);
        }
    }
}
