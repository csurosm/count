/*
 * DoubleParetoLogNormal.java
 *
 * Created on August 4, 2006, 12:20 PM
 */

package ca.umontreal.iro.matek;

/**
 * Class for modeling two-sided Pareto-lognormal distribution
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

public class DoubleParetoLogNormal extends BasicExecutable 
{
    
    private static String EPS_OUTPUT = null;
    private static boolean ML_ONLY = false;
    private static int FIT_MIN = 0;
    private static int FIT_MAX = Integer.MAX_VALUE;
    private static boolean RIGHT_TAIL=true;
    private static boolean LEFT_TAIL=true;
    
    /** 
     * Creates an object representing a double Pareto-lognormal distribution
     * with the given parameters.
     *
     * @param nu center parameter for the lognormal
     * @param tau standard deviation parameter or the lognormal
     * @param alpha slope parameter for the Pareto on the right-hand side
     * @param beta slope parameter for the Pareto on the left-hand side
     */
    public DoubleParetoLogNormal(double alpha, double beta, double nu, double tau) 
    {
        this.nu = nu;
        this.tau = tau;
        this.alpha = alpha;
        this.beta = beta;
        
        init();
    }
    
    public DoubleParetoLogNormal()
    {
        this(2.0, // alpha
             1.1, // beta
             1.2, // nu
             0.8 // tau
             );
    }

    /**
     * Pareto slope for right tail
     */
    private double alpha;
    
    /**
     * Pareto slope for left tail
     */
    private double beta;
    
    /**
     * Center for lognormal
     */
    private double nu;
    
    /**
     * Width for lognormal
     */
    private double tau;
    
    /**
     * Sets the nu parameter (center for lognormal)
     */
    public void setNu(double nu)
    {
        this.nu = nu;
        init();
    }
    
    public double getNu()
    {
        return nu;
    }
    
    /**
     * Sets the tau parameter (sd for lognormal)
     */
    public void setTau(double tau)
    {
        this.tau = tau;
        init();
    }
    
    public double getTau()
    {
        return tau;
    }
    

    /**
     * Sets the alpha parameter (Pareto slope on the right)
     */
    public void setAlpha(double alpha)
    {
        this.alpha=alpha;
        init();
    }
    
    public double getAlpha(){return alpha;}

    
    /**
     * Sets the beta parameter (Pareto slope on the left)
     */
    public void setBeta(double beta)
    {
        this.beta = beta;
        init();
    }
    
    public double getBeta(){return beta;}
    
    /**
     * A(alpha, nu, tau)
     */
    private double A_alpha;
    
    /**
     * A(-beta, nu, tau)
     */
    private double A_beta;

    /**
     * Sets up additional private variables.
     */
    private void init()
    {
        if (Double.isInfinite(alpha))
        {
            A_alpha = Double.POSITIVE_INFINITY;
            alpha_times_beta_div_alpha_plus_beta = beta;
        }
        else
            A_alpha = funcA(alpha);
        if (Double.isInfinite(beta))
        {
            A_beta=Double.POSITIVE_INFINITY;
            alpha_times_beta_div_alpha_plus_beta = alpha;
        } 
        else
            A_beta = funcA(-beta);
        if (!Double.isInfinite(alpha) && !Double.isInfinite(beta))
        {
            alpha_times_beta_div_alpha_plus_beta = alpha*beta/(alpha+beta);
        }
        nu_div_tau_plus_alpha_tau = nu/tau + alpha*tau;
        nu_div_tau_minus_beta_tau = nu/tau - beta*tau;
        
        //System.out.println("nu "+nu+" tau "+tau+" alpha "+alpha+" A "+A+" nat2 "+nat2);
    }
    
    /**
     * This is the A(theta,nu,tau) expression defined by Eq. (9) in the Reed-Jorgensen paper.
     */
    private double funcA(double theta)
    {
        double d1 = theta*nu;
        double d2 = theta*tau;
        
        return Math.exp(d1+d2*d2*0.5);
    }
    
    private double alpha_times_beta_div_alpha_plus_beta;
    private double nu_div_tau_plus_alpha_tau;
    private double nu_div_tau_minus_beta_tau;

    /**
     * Computes the value of the density function at x. 
     * Works also if alpha or beta (but not both) are infinite.
     *
     * @param x a positive value
     */ 
    public double getDensity(double x)
    {
        double y = Math.log(x)/tau;
        double t1=0.;
        double t2=0.;
        if (!Double.isInfinite(alpha))
        {
            double z1 = y - nu_div_tau_plus_alpha_tau;
            double F1 = Functions.normal_cdf(z1);
            double xa = Math.pow(x,-alpha-1.0);
            t1 = A_alpha * xa * F1;
        }
        if (!Double.isInfinite(beta))
        {
            double z2 = y -nu_div_tau_minus_beta_tau;
            double F2 = 1.0-Functions.normal_cdf(z2);
            double xb = Math.pow(x,beta-1.0);
            t2 = A_beta * xb * F2;
        }
        return alpha_times_beta_div_alpha_plus_beta*(t1+t2);
    }
    
    /**
     * Computes the value of the cumulative distribution function at x.
     * Works also if alpha, beta, or both are infinite.
     *
     * @param x a positive value (if negative or 0, returns 0.0)
     */
    public double getCDF(double x)
    {
        if (x<=0.)
            return 0.;
        double y = Math.log(x)/tau;
        double t1=0., t2=0.;
        if (!Double.isInfinite(alpha))
        {
            double z1 = y - nu_div_tau_plus_alpha_tau;
            double F1 = Functions.normal_cdf(z1);
            double xa = Math.pow(x,-alpha);
            t1 = A_alpha * xa * F1;
            if (Double.isInfinite(A_alpha) && z1<0.)
            {
                double logF = Functions.log_normal_cdf(z1);
                double logA = alpha*(nu+alpha*tau*tau*0.5);
                double logxa = -Math.log(x)*alpha;
                double logt = logA+logxa+logF;
                t1 = Math.exp(logt);
            }
        }
        if (!Double.isInfinite(beta))
        {
            double z2 = y -nu_div_tau_minus_beta_tau;
            double F2 = Functions.normal_cdf(-z2);
            double xb = Math.pow(x,beta);
            t2 = A_beta * xb * F2;
            if (Double.isInfinite(A_beta) && z2>0.)
            {
                double logF = Functions.log_normal_cdf(-z2);
                double logB = -beta*(nu-beta*tau*tau*0.5);
                double logxb = Math.log(x)*beta;
                double logt = logF+logB+logxb;
                t2 = Math.exp(logt);
            }
        }
        
        double z0 = y - nu/tau;
        double F0 = Functions.normal_cdf(z0);
        
        double retval = F0;
        if (Double.isInfinite(alpha) && Double.isInfinite(beta))
            return retval;
        if (Double.isInfinite(alpha))
            retval += t2;
        else if (Double.isInfinite(beta))
            retval -= t1;
        else 
            retval += (t2*alpha-t1*beta)/(alpha+beta);

        //Verbose.message("DPln.gCDF "+x+" t1 "+t1+" t2 "+t2+" F0 "+F0+" retval "+retval);

        return retval;
    }
    
    /**
     * Computes the probability that the random variable with this distribution 
     * takes the value between x-0.5 and x+0.5.
     * Woks also if alpha or beta or both are infinite.
     * 
     * @param x a non-negative integer value
     */
    public double getMass(double x)
    {
        double f2 = getCDF(x+0.5);//x+0.5;
        double f1 = 0.0;
        if (x>0.5)
            f1 = getCDF(x-0.5); // x-0.5
        double df = f2-f1;
        //Verbose.message("PLN.gM "+x+" f1 "+f1+" f2 "+f2+" df "+df);
        return df;
    }
        
    /**
     * Computes the expected value of the distribution.
     * Works also if alpha or beta or both are infinite.
     */
    public double getMean(){
        if (alpha<=1.0)
            throw new ArithmeticException("Mean does not exist when alpha<=1");
        
        double a = 1.0;
        if (!Double.isInfinite(alpha))
            a *= alpha/(alpha-1.0);
        if (!Double.isInfinite(beta))
            a *= beta/(beta+1.0);
        double b = Math.exp(nu+tau*tau/2.0);
        
        return a*b;
    }

    
    /**
     * Computes partial derivatives of the CDF
     *
     * @param x derivatives evaluated at this x value
     * @param dyda an array of length four to store the values
     */
    public void partialDerivativesCDF(double x, double[] dyda)
    {
        double y = Math.log(x)/tau;
        double U=0., V=0.;
        double z1 = y - nu_div_tau_plus_alpha_tau;
        double z2 = y - nu_div_tau_minus_beta_tau;
        if (!Double.isInfinite(alpha))
        {
            double F1 = Functions.normal_cdf(z1);
            double xa = Math.pow(x,-alpha);
            U = A_alpha * xa * F1;

            if (z1<0. && (Double.isNaN(U) || Double.isInfinite(U)))
            {
                double logA = alpha*(nu+alpha*tau*tau*0.5);
                double logF = Functions.log_normal_cdf(z1);
                double logxa = Math.log(x)*(-alpha);
                double logU = logA+logF+logxa;
                double Uu = Math.exp(logU);
                //Verbose.message("DPln.pDCDF "+x+" U="+U+" ... log "+logU+"("+logA+"+"+logF+"+"+logxa+") U' "+Uu);
                U=Uu;
            }
        }
        if (!Double.isInfinite(beta))
        {
            double F2 = Functions.normal_cdf(-z2);
            double xb = Math.pow(x,beta);
            V = A_beta * xb * F2;

            if (z2>0. && (Double.isNaN(V) || Double.isInfinite(V)))
            {
                double logB = -beta*(nu-beta*tau*tau*0.5);
                double logF = Functions.log_normal_cdf(-z2);
                double logxb = Math.log(x)*(beta);
                double logV = logB+logF+logxb;
                double Vv = Math.exp(logV);
                //Verbose.message("DPln.pDCDF "+x+" U="+U+" ... log "+logU+"("+logA+"+"+logF+"+"+logxa+") U' "+Uu);
                V=Vv;
            }
        }

        
        double z0 = y-nu/tau; // y = ln(x)/tau
        double phi = Math.exp(-z0*z0/2.0)/Math.sqrt(2.*Math.PI);
        
        if (!Double.isInfinite(alpha) && !Double.isInfinite(beta))
        {
            double alpha_plus_beta = alpha+beta;
            double alpha_plus_beta_squared = alpha_plus_beta * alpha_plus_beta;
            double UbVa = (U*beta-V*alpha)/alpha_plus_beta_squared; 
            dyda[PARAMETER_ALPHA]
                = UbVa + ((phi + U * z1) * beta * tau + V) / alpha_plus_beta;
            dyda[PARAMETER_BETA]
                = UbVa + ((-phi + V * z2) * alpha * tau - U) / (alpha_plus_beta);
            dyda[PARAMETER_NU]
                = - alpha_times_beta_div_alpha_plus_beta*(U+V);
            dyda[PARAMETER_TAU]
                = alpha_times_beta_div_alpha_plus_beta * tau * (-alpha * U + beta * V);
        } else if (Double.isInfinite(alpha) && !Double.isInfinite(beta))
        {
            dyda[PARAMETER_ALPHA] = 0.0;
            dyda[PARAMETER_BETA] = (-phi+z2*V)*tau;
            dyda[PARAMETER_NU] = -beta * V;
            dyda[PARAMETER_TAU] = beta*(beta*tau*V-phi);
        } else if (Double.isInfinite(beta) && !Double.isInfinite(alpha))
        {
            dyda[PARAMETER_ALPHA] = (phi+z1*U)*tau;
            dyda[PARAMETER_BETA] = 0.0;
            dyda[PARAMETER_NU] = -alpha * U;
            dyda[PARAMETER_TAU] = -alpha*(alpha*tau*U-phi);
        } else {
            // simple lognormal...
            dyda[PARAMETER_ALPHA] = 0.;
            dyda[PARAMETER_BETA] = 0.;
            dyda[PARAMETER_NU] = -phi/tau;
            dyda[PARAMETER_TAU] = -phi*z0/tau;
        }
        
        /*
        if (Double.isNaN(dyda[PARAMETER_ALPHA])
            || Double.isNaN(dyda[PARAMETER_BETA])
            || Double.isNaN(dyda[PARAMETER_NU])
            || Double.isNaN(dyda[PARAMETER_TAU]))
        {
            System.out.println("DPLN.pDCDF ??? x="+x
                +" da="+dyda[PARAMETER_ALPHA]
                +" db="+dyda[PARAMETER_BETA]
                +" dn="+dyda[PARAMETER_NU]
                +" dt="+dyda[PARAMETER_TAU]
                +"; z "+z0+" "+z1+" "+z2
                +", F "+phi+" "+F1+" "+F2
                +", xpow "+xa+" "+xb
                +" U "+U+" V "+V+" UbVa "+UbVa
                +" A "+A_alpha+" "+A_beta);
        }
        */
        // Verbose.message("DPln.pD x "+x+" U "+U+" V "+V+" z1 "+z1+" F1 "+F1+" z2 "+z2+" F2 "+F2+" z0 "+z0+" phi "+phi);
    }
    
    public static final int PARAMETER_ALPHA=0;
    public static final int PARAMETER_BETA=1;
    public static final int PARAMETER_NU=2;
    public static final int PARAMETER_TAU=3;
    
    public String toString()
    {
        StringBuffer sb = new StringBuffer(getClass().getName());
        sb.append('[');
        sb.append("alpha=");
        sb.append(alpha);
        sb.append(", beta=");
        sb.append(beta);
        sb.append(", nu=");
        sb.append(nu);
        sb.append(", tau=");
        sb.append(tau);
        if (alpha>1.0)
        {
            sb.append("; mean=");
            sb.append(getMean());
        } else
        {
            sb.append("; no mean");
        }
        sb.append(']');
        return sb.toString();
            
    }
    
    private boolean LIST_LL=false;

    
    public double getEntropy(int max)
    {
        
        double Ent = 0.0;
        for (int i=2; i<max; i++)
        {
            double p = getMass(i);
            if (p != 0.)
            Ent -= p * Math.log(p);
        }
        return Ent/Math.log(2.0);
    }
    
   /**
     * Computes the log-likelihood of a distribution for the given parameters.
     * Does not use precomputed values.
     */
    private double getLogLikelihood(double[] distribution)
    {
        double ll=0.;
        //Verbose.message("PLN.gLL nu "+nu+" tau "+tau+" alpha "+alpha);
        {
            double p = getCDF(1.5);
            if (p!=0.)
            {
                double d = Math.log(p);
                double f = d * distribution[0]+distribution[1];
                ll = -f;
                if (LIST_LL)
                    Verbose.message("PLN.gLL freq 0+1 p\t"+p+" logp\t"+(-d)+" f\t"+(-f)+" nkeys\t"+(distribution[0]+distribution[1]));
            }
        }
        for (int i=distribution.length-1; i>1; i--) //0.9*distribution.length; i++) // i==0 is not considered in the computations
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
                ll -= f;
                if (LIST_LL)
                    Verbose.message("PLN.gLL freq "+i+" p\t"+p+" logp\t"+(-d)+" f\t"+(-f)+" nkeys\t"+distribution[i]+" ll \t"+ll);
            }
        }
        return ll;
    }
    
    public void FitToDistribution(int[] frequencies)
    {
        double[] f=new double[frequencies.length];
        for (int j=0; j<frequencies.length; j++)
        {
            //System.out.println(j+"\t"+frequencies[j]+"\t"+f[j]);
            f[j]=frequencies[j];
        }
        FitToDistribution(f);
    }
    
    /**
     * Sets the distribution parameters so that they fit the given 
     * sample, using likelihood maximization.
     */
    public void FitToDistribution(double[] frequencies){
        int fit_min = FIT_MIN;
        int fit_max = Math.min(FIT_MAX, frequencies.length-1); 
        
        double m = 0.0; 
        double lm = 0.0;
        double nkeys=0.;
        nkeys=frequencies[0]; 
        for (int i = 1; i<frequencies.length; i++) 
        {
            m+= i*frequencies[i]; 
            lm+=Math.log(i)*frequencies[i];
            nkeys+=frequencies[i];
        }
        // lm is average of logs, m is average of all
        Verbose.message("DPLN.FtD nkeys "+nkeys+" total length "+m);
        m /= nkeys;
        lm /= nkeys;

        
        //if (Double.isInfinite(alpha) && Double.isInfinite(beta))
        //{
        //   // simple log-normal
        //    setNu((lm*nkeys)/(nkeys-frequencies[0]));
        //    double V=0.;
        //    for (int i=1; i<frequencies.length; i++)
        //    {
        //        double d = Math.log(i)-nu;
        //        V += d*d;
        //    }
        //    setTau(Math.sqrt(V/(nkeys-frequencies[0])));
        //    Verbose.message("DPLN.FtD lognormal "+nu+" "+tau
        //        +" ("+lm+", "+Math.sqrt(V*(nkeys-frequencies[0])/nkeys)+")");
        //} else 
        {
            Random RND=new Random();
            ParameterAlpha PA = new ParameterAlpha(frequencies);
            ParameterBeta  PB = new ParameterBeta(frequencies);
            ParameterNu    PN = new ParameterNu(frequencies);
            ParameterTau   PT = new ParameterTau(frequencies);
            double[] best=new double[5];
            best[0] = Double.NaN;
            for (int rnd=0; rnd<20; rnd++)
            {
                if (!Double.isInfinite(alpha))
                {
                    setAlpha(1.1+rnd*0.1);
                    setNu(lm-1./alpha);
                } else
                    setNu(lm-rnd*0.1);
                setTau(0.1+RND.nextInt(10)/10.0);        
                if (!Double.isInfinite(beta))
                    setBeta(0.1+RND.nextInt(20)/10.0);
                try 
                {
                    for (int rep=0; rep<4; rep++)
                    {
                        if (!Double.isInfinite(beta))
                        {
                            double[] opt = FunctionMinimization.brent(Math.max(1./lm,1.0)+1.0001, beta, 5.0, PB, 1e-7);
                            PB.set(opt[0]);
                        }
                        {
                            double[] opt = FunctionMinimization.brent(lm-1.0,nu,lm,PN,1e-7);
                            PN.set(opt[0]);
                        }
                        {
                            double[] opt = FunctionMinimization.brent(0.1, tau, 10.0, PT, 1e-7);
                            PT.set(opt[0]);
                        }
                    }
                    if (!Double.isInfinite(alpha))
                    {
                        double[] opt = FunctionMinimization.brent(Math.max(1./lm,1.0)+1.0001, alpha, 5.0, PA, 1e-7);
                        PA.set(opt[0]);
                    }
                } catch (FunctionMinimization.OptimizationException E)
                {
                    // cannot use brent
                    System.err.println("DPln.FTD optimization fails in round "+(rnd+1));
                }
                //LIST_LL=true;
                double ll = getLogLikelihood(frequencies);
                //LIST_LL=false;
                if (Double.isNaN(best[0]) || ll<best[0])
                {
                    best[0] = ll;
                    best[1] = alpha;
                    best[2] = beta;
                    best[3] = nu;
                    best[4] = tau;
                }

                Verbose.message("DPln.FTD "+(1+rnd)+"\t"+ll+" "+(ll==best[0]?"*":".")+"\t"+alpha+"\t"+beta+"\t"+nu+"\t"+tau+"\t"+getMean()+"\t"+m);
            }
            setAlpha(best[1]);
            setBeta(best[2]);
            setNu(best[3]);
            setTau(best[4]);
            
            //{
            //    Verbose.message("DPln.FTD entropy "+getEntropy(frequencies.length)+" bits per key");
            //    double nbits = 0.;
            //    for (int i=1; i<frequencies.length; i++)
            //    {
            //        double p = this.getMass(i);
            //        if (p != 0.)
            //        {
            //            nbits -= frequencies[i]*Math.log(p);
            //        }
            //    }
            //    nbits /= Math.log(2.0);
            //    nbits /= Math.log(nkeys)/Math.log(4.0);
            //    System.out.println("DPln.FTD compressed sequence "+nbits+" bits ");
            //}
        }
        
        double[] dfreq = new double[frequencies.length];
        for (int i=0; i<frequencies.length; i++)
            dfreq[i] = frequencies[i];
        try {
            double[] chisq_ML = new double[2];
            {
                double[] expected = new double[frequencies.length];
                for (int i=0; i<expected.length; i++)
                    expected[i] = getMass(i)*nkeys;
                expected[0] = dfreq[0];
                Functions.Chi_square_test(dfreq, expected, 6, chisq_ML);
            }

            Verbose.message("DPln.FTD fit ML "+this+" chisq="+chisq_ML[0]+" df="+(frequencies.length-5)+" P="+chisq_ML[1]);
        } catch (ArithmeticException E)
        {
            Verbose.message("DPln.FTD fit ML "+this+" chisq test fails "+E);
        }

        //for (int fitrep=0; fitrep<3; fitrep++)
        if (!ML_ONLY)
        {
            //Verbose.message("DPln.FTD ########################### "+fitrep);
            double[] datax = new double[fit_max-fit_min+1];
            double[] datay = new double[fit_max-fit_min+1];
            dfreq = new double[fit_max-fit_min+1];
            //{
            //    datax[0]=1.0;
            //    datay[0]=(frequencies[0]+frequencies[1])/nkeys;
            //}
            for (int i=0; i<datax.length ; i++)
            {
                datax[i]=fit_min+i;
                datay[i]=frequencies[fit_min+i]/nkeys;
                dfreq[i] = frequencies[i+fit_min];
            }

            double[] alamda = new double[1];
            alamda[0]=-1.;
            
            double[] params = new double[4];
            params[PARAMETER_ALPHA] = encode1(alpha);
            params[PARAMETER_BETA] = encode0(beta);
            params[PARAMETER_NU] = nu;
            params[PARAMETER_TAU] = encode0(tau);
            
            DPlnModel M = new DPlnModel();
            
            double[] sig = new double[datax.length];
            for (int i=0; i<sig.length; i++)
            {
                double p = getMass(datax[i]);
                sig[i] = Math.sqrt(p*(1.-p)/nkeys);
                if (sig[i]==0.)
                    sig[i] = Math.sqrt(1./nkeys);
                
            }
 
            boolean[] ia = new boolean[4];
            ia[PARAMETER_ALPHA] = !Double.isInfinite(alpha);
            ia[PARAMETER_BETA] = !Double.isInfinite(beta);
            ia[PARAMETER_NU]=true;
            ia[PARAMETER_TAU]=true;
            int num_fitted_parameters = 2;
            if (ia[PARAMETER_ALPHA])
                num_fitted_parameters++;
            if (ia[PARAMETER_BETA])
                num_fitted_parameters++;
            
            double[][] covar = new double[num_fitted_parameters][num_fitted_parameters];
            double[][] fit_alpha = new double[num_fitted_parameters][num_fitted_parameters];

            double oldchi=Double.NaN;
            for (int rep=0; rep<100; rep++)
            {
                double chi=LeastSquaresFit.mrqmin(datax, datay,sig,params,ia,covar,fit_alpha,M,alamda);
                Verbose.message("DPln.FTD chi "+chi+" alamda "+alamda[0]
                    +"\t"+decode1(params[PARAMETER_ALPHA])
                    +"\t"+decode0(params[PARAMETER_BETA])
                    +"\t"+params[PARAMETER_NU]
                    +"\t"+decode0(params[PARAMETER_TAU]));   
                if (!Double.isNaN(oldchi))
                {
                    double diff = oldchi-chi;
                    if (diff>0 && diff<0.01)
                        break;
                }
                oldchi=chi;
            }
            if (!Double.isNaN(oldchi))
            {
                alamda[0]=0.0;
                double chi=oldchi;//LeastSquaresFit.mrqmin(datax, datay,sig,params,ia,covar,alpha,M,alamda);
                setAlpha(decode1(params[PARAMETER_ALPHA]));
                setBeta(decode1(params[PARAMETER_BETA]));
                setNu(params[PARAMETER_NU]);
                setTau(decode0(params[PARAMETER_TAU]));
                
                try {
                    Verbose.message("DPln.FTD least-squares fit chi="+chi+" df "+(datax.length-4.)+" P "+Functions.Chi_square_tail(datax.length-4.,chi));
                } catch (ArithmeticException E)
                {
                    Verbose.message("DPln.FTD least-squares fit chisq test fails "+E);
                }
                {
                    double[] chisq_Fit=new double[2];
                    double[] expected = new double[datay.length];
                    for (int i=0; i<datay.length; i++)
                        expected[i] = getMass(datax[i])*nkeys;
                    try 
                    {
                        Functions.Chi_square_test(dfreq, expected, 6, chisq_Fit);
                        Verbose.message("DPln.FTD fit Chi "+this+" chisq="+chisq_Fit[0]+" df="+(datax.length-5)+" P="+chisq_Fit[1]);
                    
                    } catch (ArithmeticException E)
                    {
                        Verbose.message("DPln.FTD fit Chi "+this+" chisq test fails chisq="+Chi_square(dfreq, expected)+"df="+(datax.length-5)+"; "+E);
                    }
                }
                
            }
            
        }

        
    }
    
    private static double Chi_square(double[] data, double[] expected)
    {
        double chisq=0.0;
        for (int i=0; i<data.length; i++)
        {
            double d = (data[i]-expected[i]);
            chisq += d*d/expected[i];
        }
        return chisq;
    }
    
    //
    //
    //
    //
    //
    //
    // --------------------------------------------------------------------------------
    //
    //
    //
    //
    //
    //
    
    
    private void go(String[] args) throws java.io.IOException {
        this.reportLaunch(args);
        String file = args[0];
        reportOtherArguments("Fit min "+FIT_MIN);
        reportOtherArguments("Fit max "+FIT_MAX);
        
        Vector Vfreq=new Vector();
        Vector Vcount=new Vector();
        
        BufferedReader BR = new BufferedReader(new FileReader(file));
        String line = null;
        int max_freq = 0;
        double nkeys=0.0;
        double mean = 0.0;
        do {
            line = BR.readLine();
            if (line != null && line.trim().length() != 0 && !line.startsWith("#") ){
                String[] col = StringSplit.splitAtTabs(line);
                int freq = Integer.parseInt(col[0].trim());
                double count = Double.parseDouble(col[1].trim());
                Vfreq.add(new Integer(freq));
                Vcount.add(new Double(count));
                if (freq>max_freq) max_freq = freq;
                nkeys += count;
                mean += count*freq;
                
            }
        } while (line != null);
        BR.close();
        
        mean /= nkeys;
        
        int cutoff = max_freq; //(int)Math.min(max_freq*0.9+0.5, 4.0*mean);
        double[] distribution = new double[cutoff+1];
        for (int i=0; i<Vfreq.size(); i++){
            int freq = ((Integer)Vfreq.get(i)).intValue();
            double count = ((Double)Vcount.get(i)).doubleValue();
            if (freq<distribution.length ){
                distribution[freq]=count;
                //Verbose.message("PLN.g distr "+freq+" "+count);
            }
        }
        
        if (!LEFT_TAIL)
            beta = Double.POSITIVE_INFINITY;
        if (!RIGHT_TAIL)
            alpha = Double.POSITIVE_INFINITY;
            
        FitToDistribution(distribution);
        
        if (EPS_OUTPUT != null)
        {
            double empirical_avg_frequency=0.;

            System.out.println("% data and fit");
            System.out.println("[");
            for (int i=0; i<Math.min(distribution.length,65530/2); i++)
                if (distribution[i]!=0.){
                    System.out.println("   "+(i==0?0.5:i)+"\t"+distribution[i]);
                    empirical_avg_frequency += i*distribution[i];
                }
            System.out.println("] /data.observation-"+EPS_OUTPUT+" exch def");
            System.out.println("[");
            for (int i=0; i<Math.min(distribution.length,65530/2); i++){
                double Ekeys = getMass(i)*nkeys;
                double diff = Ekeys-distribution[i];
                System.out.println("   "+(i==0?0.5:i)+"  \t"+Ekeys+"\t% "+((int)(100.0*(diff*diff)/Ekeys)));
                if (Math.log(i+0.1) > nu && Ekeys<1e-2)
                    break;
            }
            System.out.println("] /data.fitted-"+EPS_OUTPUT+" exch def");
            System.out.println(
                "/param.alpha-"+EPS_OUTPUT+" "+(Double.isInfinite(alpha)?9e9:alpha)+" def\n"
                +"/param.beta-"+EPS_OUTPUT+" "+(Double.isInfinite(beta)?9e9:beta)+" def\n"
                +"/param.nu-"+EPS_OUTPUT+" "+nu+" def\n"
                +"/param.tau-"+EPS_OUTPUT+" "+tau+" def");
            System.out.println("/nkeys-"+EPS_OUTPUT+" "+nkeys+" def");
            
            System.out.println("% average key frequencies: empirical "+empirical_avg_frequency/nkeys+" DPL "+getMean());

            //double true_alpha = alpha;
            //
            //setAlpha(Double.POSITIVE_INFINITY);
            //System.out.println("[");
            //for (int i=0; i<distribution.length; i++){
            //    double Ekeys = getMass(i)*nkeys;
            //    System.out.println("   "+(i==0?0.5:i)+"  \t"+Ekeys);
            //    if (Math.log(i+0.1) > nu && Ekeys<1e-2)
            //        break;
            //}
            //System.out.println("] /data.fitted-"+EPS_OUTPUT+"-noright exch def");
            
            //double true_beta = beta;
            //
            //setBeta(Double.POSITIVE_INFINITY);
            //System.out.println("[");
            //for (int i=0; i<distribution.length; i++){
            //    double Ekeys = getMass(i)*nkeys;
            //    System.out.println("   "+(i==0?0.5:i)+"  \t"+Ekeys);
            //    if (Math.log(i+0.1) > nu && Ekeys<1e-2)
            //        break;
            //}
            //System.out.println("] /data.fitted-"+EPS_OUTPUT+"-ln exch def");
            //setAlpha(true_alpha);
            //setBeta(true_beta);
        }

        Verbose.message("DPln.go Fitted parameters");
        Verbose.message("DPln.go alpha (right-hand Pareto parameter) "+alpha);
        Verbose.message("DPln.go beta (left-hand Pareto parameter) "+beta);
        Verbose.message("DPln.go nu (center of lognormal) "+nu);
        Verbose.message("DPln.go tau (standard deviation of lognormal) "+tau);

    }
    

    /**
     * Fits the parameters to a distribution.
     * Call with file name as parameter. File syntax: '#' for remarks; 2 columns: frequency, nkeys separated by tabs.
     */
    public static void main(String args[])  {
        Verbose.setVerbose(true);
        //if (true)
        //{
        //    DoubleParetoLogNormal D = new DoubleParetoLogNormal(1.7,Double.POSITIVE_INFINITY,0.8,0.6);
        //    double F = D.getCDF(3.);
        //    System.out.println(F);
        //    double[] dyda = new double[4];
        //    D.partialDerivativesCDF(3.0, dyda);
        //    System.out.println(dyda[PARAMETER_ALPHA]+"\t"+dyda[PARAMETER_NU]+"\t"+dyda[PARAMETER_TAU]);
        //    return;
        //}
        //if (true)
        //{
        //    DoubleParetoLogNormal D = new DoubleParetoLogNormal(3.0,9.0,-0.5,0.5);
        //    double K = (1L<<32)*(1.-D.getCDF(40));
        //    System.out.println(K);
        //    return;
        //}
                
        DoubleParetoLogNormal O = new DoubleParetoLogNormal();
        int num_switches = 0;
        try {
            while (args.length>2*num_switches && args[2*num_switches].startsWith("-")){
                String arg_switch = args[2*num_switches].substring(1);
                if (arg_switch.equals("h"))
                    O.go(new String[0]); // will throw an Exception
                if (args.length==2*num_switches+1)
                    throw new IllegalArgumentException("Missing argument for switch "+args[2*num_switches]);
                String arg_value = args[2*num_switches+1];
                if (arg_switch.equals("v"))
                    Verbose.setVerbose(arg_value.equals("true"));
                else if (arg_switch.equals("fitmin"))
                    FIT_MIN = Integer.parseInt(arg_value);
                else if (arg_switch.equals("fitmax"))
                    FIT_MAX = Integer.parseInt(arg_value);
                else if (arg_switch.equals("mlonly"))
                    ML_ONLY = arg_value.equals("true");
                else if (arg_switch.equals("eps"))
                    EPS_OUTPUT = new String(arg_value);
                else if (arg_switch.equals("right"))
                    RIGHT_TAIL = arg_value.equals("true");
                else if (arg_switch.equals("left"))
                    LEFT_TAIL = arg_value.equals("true");
                else
                    throw new IllegalArgumentException("Switch not recognized: '"+args[2*num_switches]+"'");
                    
                num_switches++;
            }
            
            String[] rest=new String[args.length-2*num_switches];
            for (int j=0; j<rest.length; j++)
                rest[j]=args[2*num_switches+j];
            O.go(rest);
        }
         catch (Exception E){
             die(E);
         }
    }
    
    //
    //
    // ------------------------------------- FITTING A DISTRIBUTION
    //
    //
    //
    
    
    private class ParameterAlpha implements OneParameterFunction
    {
        private ParameterAlpha(double[] d)
        {
            this.distribution=d;
        }
        
        private double[] distribution;
        public void set(double a)
        {
            setAlpha(a);
        }
        
        public double eval(double a)
        {
            set(a);
            double ll = 
                getLogLikelihood(distribution);
            return ll;
        }
    }
    
    private class ParameterBeta implements OneParameterFunction
    {
        private ParameterBeta(double[] d)
        {
            this.distribution=d;
        }
        
        private double[] distribution;
        public void set(double b)
        {
            setBeta(b);
        }
        
        public double eval(double b)
        {
            set(b);
            double ll = 
                getLogLikelihood(distribution);
            return ll;
        }
    }
    
    private class ParameterNu implements OneParameterFunction
    {
        private ParameterNu(double[] d)
        {
            this.distribution=d;
        }
        
        private double[] distribution;
        public void set(double x)
        {
            setNu(x);
        }
        
        public double eval(double x)
        {
            set(x);
            double ll = 
                getLogLikelihood(distribution);
            return ll;
        }
    }
    
    private class ParameterTau implements OneParameterFunction
    {
        private ParameterTau(double[] d)
        {
            this.distribution=d;
        }
        
        private double[] distribution;
        public void set(double x)
        {
            setTau(x);
        }
        
        public double eval(double x)
        {
            set(x);
            double ll = 
                getLogLikelihood(distribution);
            return ll;
        }
    }
    
    
    /**
     * Class used in Maximum Likelihood optimization for setting parameter values.
     */
    private class Parameters implements MultiParameterFunction 
    {

        private Parameters(double[] d, double avg_log)
        {
            this.distribution = d;
            this.lm = avg_log;
        }
        private double[] distribution;
        private double lm;
        
        public void set(double[] param)
        {
            setAlpha(to_constrained(param[0],1.0,8.));
            setBeta(to_constrained(param[1],0.1,8.0));
            setNu(to_constrained(param[2],lm-2.,lm+0.5));
            setTau(to_constrained(param[3],0.05,3.));
        }
        
        public void optimizeParameters()
        {
            FunctionMinimization.powell(init_values(),1e-4,this);
        }
        
        public double[] init_values()
        {
            double[] param = new double[4];
            param[0] = to_unconstrained(alpha,1.0,8.0);
            param[1] = to_unconstrained(beta,0.1,8.0);
            param[2] = to_unconstrained(nu,lm-2.,lm+0.5);
            param[3] = to_unconstrained(tau,0.05,3.);
            return param;
        }
        
        public double eval(double[] param) {
            set(param);
            double ll = 
                getLogLikelihood(distribution);
            Verbose.message("DPlN.P.e "+alpha+"\t"+beta+"\t"+nu+"\t"+tau+"\t"+ll);
            return ll;
        }
        
        private double to_unconstrained(double x, double xmin, double xmax)
        {
            double xp = (x-xmin)/(xmax-x);
            double y = Math.log(xp);
            return y;
        }
        
        private double to_constrained(double y, double xmin, double xmax)
        {
            double yy = Math.exp(y);
            double x = (xmax*yy+xmin)/(1.+yy);
            return x;
        }
    
    }

    
    /**
     * Encoding for parameters in the range 1..infinity
     */
    public static final double encode1(double alpha)
    { 
        return alpha;
        //return Math.log(alpha-1.);
    }

    /**
     * Decoding for parameters in the range 1..infinity
     */
    public static final double decode1(double a)
    {
        return a;
        //return 1.+Math.exp(a);
    }

    /**
     * Derivative for encoded parameters in the range 1..infinity
     */
    public static final double derivative1(double a)
    {
        return 1.0;
        //return Math.exp(a);
    }
    

    /** 
     * Encoding for positive parameters
     */
    public static final double encode0(double beta)
    {
        return beta;
        //return Math.log(beta);
    }

    /**
     * Decoding for positive parameters
     */
    public static final double decode0(double b)
    {
        return b;
        //return Math.exp(b);
    }

    /**
     * Derivative for encoded positive parameters
     */
    public static final double derivative0(double b)
    {
        return 1.0;
        //return Math.exp(b);
    }
    
    
    private class DPlnModel extends LeastSquaresFit.Model
    {
        
        private double[] dyda2=null;
        private DPlnModel()
        { 
            dyda2 = new double[4];
        }
        
        public double funcs(double x, double[] a, double[] dyda) 
        {
            setAlpha(decode1(a[PARAMETER_ALPHA]));
            setBeta(decode0(a[PARAMETER_BETA]));
            setNu(a[PARAMETER_NU]);
            setTau(decode0(a[PARAMETER_TAU]));
            partialDerivativesCDF(x+0.5, dyda);
            if (x>0.5)
            {
                partialDerivativesCDF(x-0.5, dyda2);
                for (int i=0; i<4; i++)
                    dyda[i] -= dyda2[i];
            }
            
            dyda[PARAMETER_ALPHA]*=derivative1(a[PARAMETER_ALPHA]);
            dyda[PARAMETER_BETA] *=derivative0(a[PARAMETER_BETA]);
            dyda[PARAMETER_TAU]  *=derivative0(a[PARAMETER_TAU]);
            
            //for (int i=0; i<4; i++)
            //{
            //    if (Double.isInfinite(dyda[i]) || Double.isNaN(dyda[i]))
            //    {
            //        for (int j=0; j<4; j++)
            //            System.out.println("ops x="+x+" "+j+" "+dyda[j]+" "+dyda2[j]+" "+getMass(x)+" "+DoubleParetoLogNormal.this);
            //        System.exit(88);
            //    }
            //}
            
            return getMass(x);
        }
    }
    
}
