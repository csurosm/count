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
import static count.matek.Functions.EPS;
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

    private FunctionMinimization(){}
    
 
    /**
     * Maximum number of iterations in {@link #zbrent(DoubleFunction, double, double, double)}
     */
    private static final int BRENT_ITMAX=200;

    private static double sign(double a, double b)
    {
        if (b>=0.0)
            return Math.abs(a);
        else
            return -Math.abs(a);
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
     * @return x value that minimizes the function within the bracket
     */ 
    public static double zbrent(DoubleFunction<Double> func,double x1, double x2, double tol) 
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
     * Exception thrown when too many iteration in one of the routines.
     */
    public static class OptimizationException extends ArithmeticException
    {
        private OptimizationException(String message){
            super(message);
        }
    }


    
}
