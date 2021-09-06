/*
 * LowerTriangularMatrix.java
 *
 * Created on September 16, 2005, 11:00 AM
 */

package ca.umontreal.iro.matek;

/**
 *
 * @author  csuros
 */
public class LowerTriangularMatrix {
    
    /** Creates a new instance of LowerTriangularMatrix.
     * @param A is a lower triangular matrix such that A[i][j]=0 if j>i and A[i][i]=1.0. 
     * (Actually, A[i][j] is never accessed with j>=i so the row vectors can have different size.)
     */
    public LowerTriangularMatrix(double[][] A) {
        this.A=A;
        this.n=A.length;
    }
    
    private double[][] A;
    private int n;
    
    /**
     * Solves the linear equation Ax=b for x by forward substitution.
     * @param b is a vector of length identical to the number of rows in A.
     * @return the solution x
     */
    public double[] solve(double[] b){
        double[] x=new double[n];
        x[0]=b[0];
        for (int i=1; i<n; i++){
            x[i]=b[i];
            for (int j=0; j<i; j++)
                x[i]-=x[j]*A[i][j];
        }
        return x;
    }
    
    /**
     * Calculates y=Ax.
     * @param x is a vactor of length identical to the number of rows in A.
     * @return the solution y.
     */
    public double[] multiply(double[] x){
        double[] y=new double[n];
        for (int i=0; i<n; i++){
            y[i]=x[i];
            for (int j=0; j<i; j++)
                y[i]+=A[i][j]*x[j];
        }
        return y;
    }
    
    /**
     * Iterative improvement of a solution for Ax=b.
     * The idea is to set A(x+e)=b+d where e and d are error terms.
     * (x+e) is a current solution that we want to improve, and find 
     * the solution for x with Ae=d.
     * Write Ae = A(x+e)-b and solve for e.
     * The idea is from Numerical Recipes 2.5.
     * @param x is a previous solution
     * @param b is the right hand side of the original linear equation
     * @return maximum relative change to the x[i] values ;updates x by subtracting e.
     */ 
    public double improveSolution(double[] x, double[] b){
        double[] f=multiply(x);
        for (int j=0; j<n; j++)
            f[j]-=b[j];
        // now Ae=f
        double[] e=solve(f);
        double max_chg = 0.0;
        for (int j=0; j<n; j++){
            double rel_error = 0.0;
            if (x[j]!=0.0) rel_error=Math.abs(e[j]/x[j]);
            //System.out.println("#**LTM.iS e["+j+"] "+e[j]+" x "+x[j]+" rel "+rel_error);
            max_chg = Math.max(rel_error,max_chg);
            x[j]-=e[j];
        }
        return max_chg;
    }
}
