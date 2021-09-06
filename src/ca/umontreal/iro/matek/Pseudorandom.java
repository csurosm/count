/*
 * Pseudorandom.java
 *
 * Created on November 12, 2006, 9:04 PM
 */

package ca.umontreal.iro.matek;

/**
 * Collection of pseudo-random generators.
 *
 * @author  csuros
 */

import java.util.Random;

public class Pseudorandom 
{
    
    /** Creates a new instance of Pseudorandom */
    public Pseudorandom(Random RND) 
    {
        this.RND=RND;
    }
    
    public Pseudorandom()
    {
        this(new Random());
    }
    
    private Random RND;

    /**
     * A random value by Exponential(lambda). Uses the classic transformation method. 
     *
     * @param lambda the rate parameter (expected value 1/lambda)
     *
     */
    public final double nextExponential(double lambda)
    {
        return nextExponential(lambda, RND);
    }
    
    
    /**
     * A random value by Exponential(lambda). Uses the classic transformation method. 
     *
     * @param lambda the rate parameter (expected value 1/lambda)
     * @param RND random uniform[0,1] number generator 
     *
     */
    public static final double nextExponential(double lambda, Random RND)
    {
        return -Math.log(RND.nextDouble())/lambda;        
    }
    

    /**
     * Beta(a,b) over [0,1] for a<1 using Johnk's method (Devroye IX.3.5).
     */
    public final double nextBetaJohnk(double a, double b)
    {
        return nextBetaJohnk(a, b, RND);
    }
    
    /**
     * Beta(a,b) over [0,1] for a<1 using Johnk's method (Devroye IX.3.5).
     */
    public static final double nextBetaJohnk(double a, double b, Random RND)
    {
        double x,y;
        do 
        {
            double u = RND.nextDouble();
            double v = RND.nextDouble();
            x = Math.pow(u, 1./a);
            y = Math.pow(v, 1./b);
        } while (x+y>= 1.);
        double z = x/(x+y);
        return z;
    }
    
    /**
     * Geometric random variable Pr{X=k} = p*(1-p)^{k-1}, k=1,2,...
     */
    public final int nextGeometric(double p)
    {
        return nextGeometric(p, RND);
    }


    /**
     * Lognormal random variable: ln X is Normal(mu,tau)
     * @param mu expected value of the (natural) logarithm
     * @param tau standard deviation of the (natural) logarithm
     * @return a lognormally distributed random variable
     */
    public final double nextLogNormal(double mu, double tau)
    {
        return nextLogNormal(mu,tau,RND);
    }
    
    /**
     * Lognormal random variable: ln X is Normal(mu,tau) 
     * @param mu expected value of the (natural) logarithm
     * @param tau standard deviation of the (natural) logarithm
     * @param RND random number generatpr
     * @return a lognormally distributed random variable
     */
    public final double nextLogNormal(double mu, double tau, Random RND)
    {
        return Math.exp(mu+RND.nextGaussian()*tau);
    }

    /**
     * Geometric random variable Pr{X=k} = p*(1-p)^{k-1}, k=1,2,...
     */
    public static final int nextGeometric(double p, Random RND)
    {
        // Then Pr{x=k} = p*(1-p)^(k-1) and Pr{x<=k} = 1-(1-p)^k.
        // In order to generate a r.v. with this distribution,
        // we use ceiling(log U/log (1-p)), where U is uniform(0,1).
        
        double z = Math.log(1.-p);
        double U = RND.nextDouble();
        double R = Math.log(U)/z;
        int retval = (int)(R+1.); // poor man's ceiling 
        return retval;
    }
    
    /**
     * Yule or Yule-Simon distribution
     */
    public final int nextYule(double rho)
    {
        return nextYule(rho, RND);
    }
    
    
    /**
     * Uniform[0,1] random variable (using Random.nextDouble())
     */
    public final double nextUniform()
    {
        return RND.nextDouble();
    }
    

    /**
     * Uniform[0..range] integer-valued random variable (using Random.nextInt(n))
     */
    public final int nextUniform(int range)
    {
        return RND.nextInt(range);
    }
    
    /**
     * Yule or Yule-Simon distribution
     */
    public static final int nextYule(double rho, Random RND)
    {
        double w = nextExponential(rho, RND);
        double p = Math.exp(-w); // Uniform ^ {1/rho} 
        return nextGeometric(p, RND);
    }
}
