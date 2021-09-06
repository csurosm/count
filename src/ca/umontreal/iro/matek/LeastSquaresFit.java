/*
 * LeastSquaresFit.java
 *
 * Created on August 8, 2006, 1:53 PM
 */

package ca.umontreal.iro.matek;

/**
 * Programs for computing least-squares fit to data points.
 *
 * @author  csuros
 */

import ca.umontreal.iro.banality.Verbose;

public class LeastSquaresFit 
{
    
    /** Not used. **/
    private LeastSquaresFit() {}
    
    private static int mfit;
    private static double ochisq;
    private static double[] atry;
    private static double[] beta;
    private static double[] da;
    private static double[][] oneda;
    
    
    /**
     * Levenberg-Marquardt method, attempting to reduce the value chi-square of 
     * a fit between a set of data 
     * points x[0..ndata-1], y[0..ndata-1] with individual standard deviations 
     * sig[0..ndata-1], 
     * and a non-linear function dependent on ma coefficients a[0..ma-1]. 
     * The input array ia[0..ma-1] 
     * indicates by true entries those components of a that should be fitted for, and by false 
     * entries those components that should be held fixed at their input values. 
     * The program returns current best-fit values for the parameters a[0..ma-1], and chi-square. 
     * The arrays covar[0..ma-1][0..ma-1], alpha[0..ma-1][0..ma-1] are used as working space during most 
     * iterations. Supply a routine funcs(x,a,dyda) that evaluates the fitting function, 
     * and its derivatives dyda[0..ma-1] with respect to the fitting parameters a at x. On 
     * the first call provide an initial guess for the parameters a, and set alamda<0 for initialization 
     * (which then sets alamda=.001). If a step succeeds chisq becomes smaller and alamda de- 
     * creases by a factor of 10. If a step fails alamda grows by a factor of 10. You must call this    
     * routine repeatedly until convergence is achieved. Then, make one final call with alamda=0, so 
     * that covar[0..ma-1][0..ma-1] returns the covariance matrix, and alpha the curvature matrix. 
     * (Parameters held fixed will return zero covariances.) 
     *
     * Based on Numerical Recipes 15.5
     *
     * @param alamda only entry [0] is used
     *
     * @return chi-square
     *
     */   
    public static double mrqmin(double x[], double y[],double sig[], double[] a, boolean[] ia, 
	double[][] covar, double[][] alpha, Model F, double[] alamda)
    {  
        // intj,k,l;intj,k,l; 
        int ma = ia.length;
        
        double chisq;
        
        if(alamda[0]<0.0)
        { // Initialization. 
            atry=new double[ma]; 
            beta=new double[ma]; 
            da=new double[ma]; 
            mfit=0;
            for(int j=0;j<ma;j++) 
                if(ia[j])mfit++; 
            oneda = new double[mfit][1];
            alamda[0]=0.001; 
            ochisq=chisq = mrqcof(x,y,sig,a,ia,alpha,beta,F); 
            for(int j=0;j<ma;j++)
                atry[j]=a[j]; 
        } 
        for(int j=0;j<mfit;j++)
        { 
            // Alter linearized fitting matrix, by augmenting diagonal elements.
            for(int k=0;k<mfit;k++)
                covar[j][k]=alpha[j][k]; 
            covar[j][j]=alpha[j][j]*(1.0+(alamda[0])); 
            oneda[j][0]=beta[j]; 
            //Verbose.message("LSF.mrqmin descent "+j+" "+beta[j]+"/("+alpha[j][j]+"*"+alamda[0]+")");
        } 
        if (false)
        {
            StringBuffer sb=new StringBuffer();
            sb.append("covar = {");
            for (int i=0; i<covar.length; i++)
            {
                sb.append("{");
                for (int j=0; j<covar[i].length; j++)
                {
                    if (j!=0)
                        sb.append(", ");
                    sb.append(covar[i][j]);
                }
                sb.append("}");
            }
            sb.append("} "); 
            sb.append(" oneda={");
            for (int j=0; j<oneda.length; j++)
            {
                if (j!=0)
                    sb.append(", ");
                sb.append(oneda[j][0]);
            }
            sb.append("}");
            Verbose.message("LSF.mrqmin solve "+sb);
        }
        GaussJordan.gaussj(covar,oneda); // Matrix solution. 
        if (false)
        {
            StringBuffer sb=new StringBuffer();
            sb.append("covar = {");
            for (int i=0; i<mfit; i++)
            {
                sb.append("{");
                for (int j=0; j<mfit; j++)
                {
                    if (j!=0)
                        sb.append(", ");
                    sb.append(covar[i][j]);
                }
                sb.append("}");
            }
            sb.append("} "); 
            sb.append(" oneda={");
            for (int j=0; j<mfit; j++)
            {
                if (j!=0)
                    sb.append(", ");
                sb.append(oneda[j][0]);
            }
            sb.append("}");
            Verbose.message("LSF.mrqmin gaussj "+sb);
        }
        for(int j=0;j<mfit;j++)
            da[j]=oneda[j][0]; 
        if(alamda[0]==0.0)
        { // Once converged, evaluate covariance matrix. 
            covsrt(covar,ia); 
            covsrt(alpha,ia); // Spread out alpha to its full size too. 
            return ochisq; 
        } 
        for(int j=0,l=0;l<ma;l++) // Did the trial succeed? 
            if(ia[l])
                atry[l]=a[l]+da[j++]; 
        //Verbose.message("LSF.mrqmin try "+atry[0]+"\t"+atry[1]+"\t"+atry[2]+"\t"+atry[3]);
        chisq=mrqcof(x,y,sig,atry,ia,covar,da,F); 
        if(chisq<ochisq){ // Success, accept the new solution. 
            //Verbose.message("LSF.mrqmin better "+atry[0]+"\t"+atry[1]+"\t"+atry[2]+"\t"+atry[3]+"\t"+ochisq+" > "+chisq);
            alamda[0]*=0.1; 
            ochisq=chisq; 
            for(int j=0;j<mfit;j++)
            { 
                for(int k=0;k<mfit;k++)
                    alpha[j][k]=covar[j][k]; 
                beta[j]=da[j]; 
            } 
            for(int l=0;l<ma;l++)
                a[l]=atry[l]; 
        } else
        { // Failure, increase alamda and return. 
            //Verbose.message("LSF.mrqmin fail "+atry[0]+"\t"+atry[1]+"\t"+atry[2]+"\t"+atry[3]+"\t"+ochisq+" < "+chisq);
            
            alamda[0]*=10.0; 
            chisq=ochisq; 
        } 
        
        return chisq;
    }
    
    /**
     * Used by mrqmin to evaluate the linearized fitting matrix alpha, and vector beta, 
     * and calculate chi-square. 
     *
     * @return chi-square
     */
    private static double mrqcof(double x[], double y[], double sig[], double a[], boolean ia[], 
        double[][] alpha, double beta[], Model F)
    {
        int ndata = x.length; // number of data points
        int ma = ia.length; // number of parameters
        double[] dyda=new double[ma]; // where partial derivatives are stored
        
        for(int j=0;j<mfit;j++)
        { 
            // Initialize alpha, beta. 
            for(int k=0;k<mfit;k++)alpha[j][k]=0.0; 
            beta[j]=0.0; 
        } 
        //
        // beta[j] equals
        // \sum (y_i-y(x;a))/sigma_i^2 * dyda[j]
        //
        
        double chisq=0.0; 
        for(int i=0;i<ndata;i++)
        { // Summation loop over all data. 
            double ymod = F.funcs(x[i], a, dyda);
            // now ymod is the predicted value
            double sig2i=1.0/(sig[i]*sig[i]); 
            double dy=y[i]-ymod; 
            int j=0;
            for(int l=0;l<ma;l++)
            { 
                // for all parameters
                if(ia[l])
                { 
                    double wt=dyda[l]*sig2i; 
                    int k=0;
                    for(int m=0;m<=l;m++) 
                        if(ia[m])
                            alpha[j][k++]+=wt*dyda[m];
                    beta[j]+=dy*wt; 
                    j++;
                } 
            } 
            // And find chi-square. 
            chisq+=dy*dy*sig2i; 
            //Verbose.message("LSF.mrqcof "+i+"\t"+x[i]+"\t"+y[i]+"\t"+ymod+"\tdy "+dy
            //    +"\t"+dyda[0]+"\t"+dyda[1]+"\t"+dyda[2]+"\t"+dyda[3]);
        } 
        for(int j=1;j<mfit;j++) // Fill in the symmetric side. 
            for(int k=0;k<j;k++)
                alpha[k][j]=alpha[j][k]; 
        
        if (Double.isInfinite(alpha[0][0]))
        {
            for(int i=0;i<ndata;i++)
            { // Summation loop over all data. 
                double ymod = F.funcs(x[i], a, dyda);
                // now ymod is the predicted value
                double sig2i=1.0/(sig[i]*sig[i]); 
                double dy=y[i]-ymod; 
                int j=0;
                for(int l=0;l<ma;l++)
                { 
                    // for all parameters
                    if(ia[l])
                    { 
                        double wt=dyda[l]*sig2i; 
                        int k=0;
                        for(int m=0;m<=l;m++) 
                            if(ia[m])
                                alpha[j][k++]+=wt*dyda[m];
                        beta[j]+=dy*wt; 
                        j++;
                    } 
                } 
                Verbose.message("LSF.mrqcof  ???? "+i+"\t"+x[i]+"\t"+y[i]+"\t"+ymod+"\tdy "+dy
                    +"\t"+dyda[0]+"\t"+dyda[1]+"\t"+dyda[2]+"\t"+dyda[3]+"\t"+sig2i+" "+sig[i]);
                
            } 
            
        }
        
        
        return chisq;
    }
    
    
    /**
     * Expand in storage the covariance matrix covar, so as to take into account parameters that are 
     * being held fixed. (For the latter, return zero covariances.) 
     *
     * Based on Numerical Recipes 15.4
     * 
     * @param covar ma*ma matrix covar[0..ma-1][0..ma-1]
     * @param ia vector of length ma ia[0..ma-1] (to denote which variables are not frozen)
     * @param mfit compute from index mfit+1 on
     *
     */
    private static void covsrt(double[][] covar, boolean[] ia)
    {
        int ma = covar.length;
        for(int i=mfit+1;i<ma;i++) 
            for(int j=0;j<i;j++) covar[i][j]=covar[j][i]=0.0; 
        int k=mfit; 
        for(int j=ma-1;j>=0;j--)
        { 
            if(ia[j])
            { 
                for(int i=0;i<ma;i++)
                {
                    double tmp=covar[i][k];
                    covar[i][k]=covar[i][j];
                    covar[i][j]=tmp;
                }
                for(int i=0;i<ma;i++)
                {
                    double tmp = covar[k][i];
                    covar[k][i]=covar[j][i];
                    covar[j][i]=tmp;
                }
                k--; 
            } 
        } 
    }
    
    public static abstract class Model
    {
        /**
         * Function that computes the derivatives for the fit
         * by Levenberg-Marquardt.
         * 
         * @param x value of x at which the function needs to be evaluated
         * @param a parameter vector
         * @param dyda vector of length ma dyda[0..ma-1] which are the partial derivatives at x
         * @return function value y
         */
        public abstract double funcs(double x, double[] a, double[] dyda);
    }
}
