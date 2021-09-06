/*
 * PowerLaw.java
 *
 * Created on December 14, 2005, 11:14 AM
 */

package ca.umontreal.iro.matek;

import ca.umontreal.iro.banality.Verbose;

/**
 * (a+n)^{-b}
 *
 * A class for computing power law distriobutions 
 * @author  csuros
 */
public class PowerLaw {
    
    /** Creates a new instance of PowerLaw */
    public PowerLaw(double a, double b) {
        shift = a;
        exponent = -b;
        computeSum();
    }
    
    
    public PowerLaw(double b){
        this(0.,b);
    }
    
    public void setShift(double a){
        shift = a;
        computeSum();
    }
    
    public double getShift(){ return shift;}
    
    public void setExponent(double b){
        exponent  = -b;
        computeSum();
    }
    
    public double getExponent(){ return -exponent;}
    
    /**
     * This is the sum of (n+a)^-b
     */
    public double getScale(){ return scaler;}
    
    private double shift;
    private double exponent;
    private double scaler;
    private double log_scaler;
    
    private static double EPS = 1e-12;
    
    /**
     * Approximates the Riemann-Zeta sum n^{-s} using Euler-Maclaurin
     *
     * @param s a positive value grater than 1
     */
    public static double Zeta(double s){
        return ShiftedZeta(0.,s);
    }
    
    /**
     * Approximates the sum (n+x)^{-s} using Euler-Maclaurin (for x=0 this is the Riemann-zeta)
     *
     * @param s a positive value grater than 1
     * @param x a value greater than -1
     */
    public static double ShiftedZeta(double x, double s){
        double sum=0.;
        int N=1;
        double s_=-s;
        for (; ; N++){
            double v = Math.pow(x+N,s_);
            sum += v;
            if (v/sum<EPS || N==1000000){
                break;
            }
        }
        // Euler-Maclaurin
        double corr=Math.pow(N+x, 1.-s)/(s-1.0);
        corr += Math.pow(N+x, s_)/2.0;
        sum += corr;
        //Verbose.message("PL.SZ s="+s+" x="+x+" done @"+N+" "+sum+" corr "+corr);
        return sum;
    }

    private void computeSum(){
        scaler=ShiftedZeta(shift,-exponent);
        log_scaler = Math.log(scaler);
    }
    
    public double valueAt(int k){
        return Math.pow(shift+k,exponent)/scaler;
    }
    
    public double logValueAt(int k){
        return exponent*Math.log(k+shift)-log_scaler;
    }
    
    public void FitToDistribution(int[] frequencies, boolean fit_shift){
        ShiftParameter SP=new ShiftParameter(frequencies);
        ExponentParameter EP=new ExponentParameter(frequencies);
        
        double[] opt = new double[2];
        for (int rep=0; rep<5; rep++){
            opt = FunctionMinimization.brent(1.0, 2.0, 5.0, EP, 1e-4);
            setExponent(opt[0]);
            if (fit_shift){
                opt = FunctionMinimization.brent(0.0, 1.0, 4.0, SP, 1e-4);
                setShift(opt[0]);
            }
        }
        
        Verbose.message("PL.FTD shift "+shift+" b "+-exponent);
    }
    
    
    private class ShiftParameter implements OneParameterFunction {
        private int[] distribution;
        
        private ShiftParameter(int[] d){
            this.distribution = d;
        }
        
        /**
         * @return negative log-likelihood
         */
        public double eval(double param) {
            setShift(param);
            
            double ll=0.;
            for (int i=0; i<distribution.length; i++)
                if (distribution[i]>0){
                    double v = logValueAt(i+1)*distribution[i];
                    ll -= v;
                }
            //Verbose.message("PL.SP.e "+param+" ll "+ll); 
            
            return ll;
        }
    }
    

    private class ExponentParameter implements OneParameterFunction {
        private int[] distribution;
        
        private ExponentParameter(int[] d){
            this.distribution = d;
        }
        
        /**
         * @return negative log-likelihood
         */
        public double eval(double param) {
            setExponent(param);
            
            double ll=0.;
            for (int i=0; i<distribution.length; i++)
                if (distribution[i]>0){
                    double v = logValueAt(i+1)*distribution[i];
                    ll -= v;
                }
            //Verbose.message("PL.EP.e "+param+" ll "+ll); 
            
            return ll;
        }
    }

    public static void main(String[] args){
        Verbose.setVerbose(true);
        PowerLaw PL = new PowerLaw(Double.parseDouble(args[0]),Double.parseDouble(args[1]));
        System.out.println("10 "+PL.valueAt(10));
        System.out.println("200 "+PL.logValueAt(200));
    }
    
}
