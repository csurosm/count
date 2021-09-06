/*
 * ParetoLogNormal.java
 *
 * Created on February 23, 2006, 11:23 PM
 */

package ca.umontreal.iro.matek;

/**
 * Class for modeling one-sided (right) Pareto-lognormal distribution
 *
 * @author  csuros
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Vector;
import java.util.Random;

import ca.umontreal.iro.banality.Verbose;
import ca.umontreal.iro.banality.Functions;
import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.StringSplit;

public class ParetoLogNormal extends BasicExecutable {
    private static final int MAX_PRECOMPUTED=3000;
    
    /** 
     * @param nu expectation for the lognormal
     * @param tau standard deviation for the lognormal
     * @param alpha parameter for the Pareto
     */
    public ParetoLogNormal(double alpha, double nu, double tau) {
        this.nu = nu;
        this.tau = tau;
        this.alpha = alpha;
        
        init();
    }
    
    private ParetoLogNormal(){
        setNu(1.2);
        setTau(0.8);
        setAlpha (2.0);
    }
    
    private double nu;
    private double tau;
    private double alpha;
    private double A;
    private double aA;
    private double nat2; // nu+alpha*t^2
    
    private double[] precomputed;
    
    private void setNu(double nu){
        this.nu = nu;
        init(0);
    }
    
    private void setTau(double tau){
        this.tau = tau;
        init(0);
    }
    
    private void setAlpha(double alpha){
        this.alpha=alpha;
        init(0);
    }
    
    private void init(){
        init(MAX_PRECOMPUTED);
    }
    
    /**
     * Sets up additional private variables and precomputes 
     * some values for the density.
     */
    private void init(int num_precomputed){
        A=Math.exp(alpha*nu+alpha*alpha*tau*tau/2.0);
        aA = A*alpha;
        logaA = (alpha*nu+alpha*alpha*tau*tau/2.0)+Math.log(alpha);
        nat2 = nu+alpha*tau*tau;
        precomputed = new double[num_precomputed];
        for (int i=0; i<precomputed.length; i++){
            precomputed[i]=getDensity(i+1.0);
        }
        //System.out.println("nu "+nu+" tau "+tau+" alpha "+alpha+" A "+A+" nat2 "+nat2);
    }

    /**
     * Computes the value of the density function at x (never uses precomputed values)
     * @param x a positive value
     */ 
    public double getDensity(double x){
        double y = (Math.log(x)-nat2)/tau;
        double phi = Functions.normal_cdf(y);
        double xa = Math.pow(x, alpha+1.);
        double f = aA*phi/xa;
        //Verbose.message("PLN.gLD x "+x+" `y "+y+" phi "+phi+" xa "+xa+" aA "+aA+" f "+f);
        
        return f;
    }
    
    public double getLogDensity(double x){
        double y = (Math.log(x)-nu)/tau-alpha*tau;
        double phi = Functions.normal_cdf(y);
        if (phi==0.0){
            phi = Double.MIN_VALUE;
            //Verbose.message("PLN.gLD ZERO x "+x+" y "+y+" phi "+phi+" nu "+nu+" tau "+tau+" "+" a "+alpha+" erfc "+Functions.erfcc(y/Functions.SQRT2));
        }
        double logf = logaA+Math.log(phi)-(alpha+1)*Math.log(x);

        return logf;
        
    }
    
    private double logaA;
    
    /**
     * Computes the value of the density function at x (uses precomputed values for small x)
     *
     * @param x a positive integer
     */
    public double getDensity(int x){
        if (x>precomputed.length)
            return getDensity((double)x);
        else
            return precomputed[x-1];
    }
    
    /**
     * Computes the value of the cumulative distribution function at x.
     */
    public double getCDF(double x)
    {
        double logx = Math.log(x);
        double z1 = (logx-nu)/tau;
        double z2 = (logx-nat2)/tau;
        double xpow_alpha = Math.pow(x, -alpha);
        double term1 = Functions.normal_cdf(z1);
        double term2 = Functions.normal_cdf(z2)*A*xpow_alpha;
        double retval = term1-term2;
        Verbose.message("PLN.gCDF nu "+nu+" tau "+tau+" alpha "+alpha+" x "+x+" ln "+logx+" F "+retval+" z1 "+z1+" z2 "+z2+" xa "+xpow_alpha+" t1 "+term1+" t2 "+term2+" d "+getDensity(x));
        return retval;
    }
    
    /**
     * Computes the probability that the random variable with this distribution 
     * takes the value between x-0.5 and x+0.5.
     * 
     * @param x a non-negative integer value
     */
    public double getMass(int x)
    {
        double f2 = getCDF(x+0.5);
        double f1 = 0.0;
        if (x>0)
            f1 = getCDF(x-0.5);
        double df = f2-f1;
        Verbose.message("PLN.gM "+x+" f1 "+f1+" f2 "+f2+" df "+df);
        return df;
    }
    
    /**
     * Computes the expected value of the distribution.
     */
    public double getMean(){
        double E = alpha*Math.exp(nu+tau*tau/2.0)/(alpha-1.0);
        return E;
    }
    
    /**
     * Sets the distribution parameters so that they fit the given 
     * sample, using likelihood maximization.
     */
    public void FitToDistribution(double[] frequencies){
        
        double m = 0.0; 
        double lm = 0.0;
        { 
            double s=0.; 
            for (int i = 0; i<frequencies.length; i++) {
                m+= (i+1.)*frequencies[i]; 
                lm+=Math.log(i+1.)*frequencies[i];
                s+=frequencies[i];
            }
            m /= s;
            lm /= s;
        }
        
        for (int i_alpha=0; i_alpha<26; i_alpha++)
        {
            setAlpha(1.0001+0.25*i_alpha);
            for (int i_nu=0; i_nu<=40; i_nu++)
            {
                setNu(lm-2.0+0.1*i_nu);
                for (int i_tau=1; i_tau<26; i_tau++)
                {
                    setTau(i_tau*0.25);
                    double ll = getLogLikelihood(frequencies);
                    System.out.println(alpha+"\t"+nu+"\t"+tau+"\t"+ll);
                }
            }
        }
        
        
        /*
        NuParameter NP = new NuParameter(frequencies, lm);
        TauParameter TP = new TauParameter(frequencies);
        AlphaParameter AP = new AlphaParameter(frequencies, lm);
        
        Verbose.message("PLN.FtD 0 nu "+nu+" t "+tau+" a "+alpha+" mean "+getMean()+" dist_avg "+m+" avg_log "+lm);

        Random RND=new Random();
        for (int rnd=0; rnd<10; rnd++)
        {
            double a = Math.max(1.0,1.0/lm)+RND.nextInt(50)/10.0;
            setAlpha(a);
            setNu(lm-1./a);
            setTau(RND.nextInt(20)/10.0);

            //setNu(lm-0.25);
            //setTau(0.7);
            //setAlpha(4.0);

            Verbose.message("PLN.FtD RND"+rnd+" nu "+nu+" t "+tau+" a "+alpha+" mean "+getMean()+" dist_avg "+m+" avg_log "+lm);

            for (int rep=0; rep<5; rep++)
            {

                {
                    double[] opt = FunctionMinimization.brent(lm-1.0,nu,lm,NP,1e-4);
                    NP.set(opt[0]);
                }
                {
                    double[] opt = FunctionMinimization.brent(0.1, tau, 10.0, TP, 1e-4);
                    TP.set(opt[0]);
                }
                {
                    double[] opt = FunctionMinimization.brent(Math.max(1./lm,1.0)+1.0001, alpha, 5.0, AP, 1e-4);
                    AP.set(opt[0]);
                }
                Verbose.message("PLN.FtD "+(1+rep)+" nu "+nu+" t "+tau+" a "+alpha+" ll "+getLogLikelihood(frequencies)+" mean "+getMean()+" dist_avg "+m);

                //for (int i=0; i<100; i++){
                //    double ld = getLogDensity(i+1.0);
                //    double d = getDensity(i+1.0);
                //    Verbose.message("PLN.FtD "+(i+1)+" cnt "+frequencies[i]+" dens "+d+" logd "+ld);
                //}
                //for (int i=0; i<100; i++){
                //    double ld = getLogDensity(i+1000.0);
                //    double d = getDensity(i+1000.0);
                //    Verbose.message("PLN.FtD "+(i+1000)+" cnt "+frequencies[i+999]+" dens "+d+" logd "+ld);
                //}
                //if (true)
                //    System.exit(1);
            }
        }
         */
               
        init();
    }
    
    
    /**
     * Class used in Maximum Likelihood optimization for setting parameter values.
     */
    private class NuParameter implements OneParameterFunction {
        private double[] distribution;
        private double y;
        private NuParameter(double[] d, double avg_log){
            this.distribution = d;
            this.y = avg_log;
        }
        
        public void set(double nu){
            setNu(nu);
            //setAlpha(1.0/(y-nu));
        }
        
        public double eval(double nu) {
            set(nu);
            return getLogLikelihood(distribution);
        }
    }

    /**
     * Class used in Maximum Likelihood optimization for setting parameter values.
     */
    private class TauParameter implements OneParameterFunction {
        private double[] distribution;
        private TauParameter(double[] d){
            this.distribution = d;
        }
        
        public void set(double tau){
            setTau(tau);
        }
        
        public double eval(double tau) {
            set(tau);
            double ll = 
                getLogLikelihood(distribution);
            Verbose.message("PLN.TP e "+tau+" "+ll+" nu "+nu+" alpha "+alpha);
            return ll;
        }
    }
    
    /**
     * Class used in Maximum Likelihood optimization for setting parameter values.
     */
    private class AlphaParameter implements OneParameterFunction {
        private double[] distribution;
        private double y;
        private AlphaParameter(double[] d, double avg_log){
            this.distribution = d;
            this.y = avg_log;
        }
        
        public void set(double alpha){
            setAlpha(alpha);
            //setNu(y-1.0/alpha);
        }
        
        public double eval(double alpha) {
            set(alpha);
            double ll = 
                getLogLikelihood(distribution);
            Verbose.message("PLN.AP.e "+alpha+" "+ll+" nu "+nu+" tau "+tau);
            return ll;
        }
    }
    

    /**
     * Computes the log-likelihood of a distribution for the given parameters.
     * Does not use precomputed values.
     */
    private double getLogLikelihood(double[] distribution)
    {
        double ll=0.;
        //Verbose.message("PLN.gLL nu "+nu+" tau "+tau+" alpha "+alpha);
        for (int i=1; i<distribution.length; i++) //0.9*distribution.length; i++) // i==0 is not considered in the computations
        {
            if (distribution[i]>0)
            {
                double p = getMass(i);
                // double d = 0.0;
                //if (df>0)
                //    d=Math.log(df);
                //else if (i!=0) 
                //    d=getLogDensity(i);
                double d = Math.log(p);
                double f = d*distribution[i];
                //if (alpha == 7.2501 && tau==3.0 && nu<0)
                //Verbose.message("PLN.gLL freq "+i+" p "+p+" logp "+d+" f "+f+" nkeys "+distribution[i]);
                ll -= f;
            }
        }
        return ll;
    }
    
    private void go(String[] args) throws java.io.IOException {
        this.reportLaunch(args);
        String file = args[0];
        
        Vector Vfreq=new Vector();
        Vector Vcount=new Vector();
        
        BufferedReader BR = new BufferedReader(new FileReader(file));
        String line = null;
        int max_freq = 0;
        double sum=0.0;
        double mean = 0.0;
        do {
            line = BR.readLine();
            if (line != null && line.trim().length() != 0 && !line.startsWith("#") ){
                String[] col = StringSplit.splitAtTabs(line);
                int freq = Integer.parseInt(col[0]);
                double count = Double.parseDouble(col[1]);
                Vfreq.add(new Integer(freq));
                Vcount.add(new Double(count));
                if (freq>max_freq) max_freq = freq;
                sum += count;
                mean += count*freq;
                
            }
        } while (line != null);
        BR.close();
        
        mean /= sum;
        
        int cutoff = max_freq; //(int)Math.min(max_freq*0.9+0.5, 4.0*mean);
        double[] distribution = new double[cutoff];
        for (int i=0; i<Vfreq.size(); i++){
            int freq = ((Integer)Vfreq.get(i)).intValue();
            double count = ((Double)Vcount.get(i)).doubleValue();
            if (freq<distribution.length ){
                distribution[freq]=count;
                //Verbose.message("PLN.g distr "+freq+" "+count);
            }
        }
            
        
        FitToDistribution(distribution);
        
        System.out.println("Fitted parameters");
        System.out.println("nu (center of lognormal) "+nu);
        System.out.println("tau (standard deviation of lognormal) "+tau);
        System.out.println("alpha (Pareto parameter) "+alpha);
    }

    /**
     * Fits the parameters to a distribution.
     * Call with file name as parameter. File syntax: '#' for remarks; 2 columns: frequency, nkeys separated by tabs.
     */
    public static void main(String args[])  {
        Verbose.setVerbose(true);
        
        //if (true)
        //{
        //    ParetoLogNormal D = new ParetoLogNormal(2.9,0.3,0.6);
        //    System.out.println("cdf[2] "+D.getCDF(2.0));
        //    System.out.println("dens[3] "+D.getDensity(3.0));
        //    System.out.println("mass[4] "+D.getMass(4));
        //}
        //else
        //{
        
        int num_switches = 0;
        try {
            while (args.length>2*num_switches && args[2*num_switches].startsWith("-")){
                String arg_switch = args[2*num_switches].substring(1);
                if (arg_switch.equals("h"))
                    (new ParetoLogNormal()).go(new String[0]); // will throw an Exception
                if (args.length==2*num_switches+1)
                    throw new IllegalArgumentException("Missing argument for switch "+args[2*num_switches]);
                String arg_value = args[2*num_switches+1];
                if (arg_switch.equals("v"))
                    Verbose.setVerbose(Boolean.getBoolean(arg_value));
                else
                    throw new IllegalArgumentException("Switch not recognized: "+args[2*num_switches]);
                    
                num_switches++;
            }
            
            String[] rest=new String[args.length-2*num_switches];
            for (int j=0; j<rest.length; j++)
                rest[j]=args[2*num_switches+j];
            (new ParetoLogNormal()).go(rest);
        }
         catch (Exception E){
             die(E);
         }
        //}
    }
    
    
    
}
