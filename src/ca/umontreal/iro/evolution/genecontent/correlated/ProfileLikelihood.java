
package ca.umontreal.iro.evolution.genecontent.correlated;

import ca.umontreal.iro.evolution.TreeNode;
import ca.umontreal.iro.matek.FunctionMinimization;
import ca.umontreal.iro.matek.OneParameterFunction;

/**
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s and Louis Philippe B. Bouvrette
 */
public abstract class ProfileLikelihood
{
        public ProfileLikelihood(TreeNode root, int[] profile, double mu, double lambda, double delta)
        {
            nodes = root.getTraversal().getDFT();
            NOGD.mapping(nodes);        // vraiment utile ??? --- necessary for proper recursion in the peeling algo
            this.profile = profile;
            this.lambda = lambda;
            this.delta = delta;
            this.mu = mu;
            // find max. rate value
            scaled_rate_max = getMaxScaledRate(nodes);
            absent_profile = new int[profile.length]; // all 0
        }

        protected double lambda ;
        protected double delta;
        protected double mu;
        protected int[] profile;
        protected TreeNode[] nodes;
        protected double scaled_rate_max;
        protected int[] absent_profile;
        
        protected static final double RATE_MAX = 10.0; // rate*edge length should not be larger than this
        protected static final double TOL = 0.0001; // relative accuracy in function minimization
        
        /**
         * Computes the maximum for any rate (loss, gain). 
         * The return value multiplied by tehe longest edge length gives the constant 
         * RATE_MAX.
         * 
         * @param nodes array of tree nodes (last entry for the root is ignored)
         * @return maximum rate value used in optimization
         */
        public static double getMaxScaledRate(TreeNode[] nodes)
        {
            double rate_max = Double.POSITIVE_INFINITY;
            for (int node_idx=0; node_idx<nodes.length-1; node_idx++) // every node except the root as the last element of nodes[]
            {
                TreeNode N = nodes[node_idx];
                double len =  N.getLength(); // length of the edge leading to N
                double this_rate_max = RATE_MAX/len;
                if (this_rate_max<rate_max)
                   rate_max = this_rate_max;
            }
            return rate_max;
        }

        public abstract double getLikelihood(int[] profile);

        final public double getLikelihood()
        {
            double p = getLikelihood(profile);
            double p0 = getLikelihood(absent_profile);
            
            //System.out.println("#*PL.gL "+p+"\t"+p0+"\t"+toString());
            
            return p/(1.-p0);        
        }
        
        final public double getNegativeLogLikelihood()
        {
            double lik = getLikelihood();
            double nloglik = -Math.log(lik); 
            return nloglik;        
        }
        
        public double optimize()
        { //... optimizeRates() dans une boucle;
            double opt_lh = getNegativeLogLikelihood();
               for (int i=0; i<100; i++){
                   double new_lh = optimizeRates();
                   //System.out.println("#*PL.o "+i+"\tlh "+new_lh+"\twas "+opt_lh+"\t"+this);
                   double zeta = opt_lh - new_lh;
                   opt_lh = new_lh;
                   
                   if (i>2 && zeta < 0.001){
                       break;
                   }
               }
               return opt_lh;
        }

        public double optimizeRates()
        {//... utiliser les classes internes
            double opt_lh = 0.0;
            {
                MuOptimization opt = new MuOptimization();
                double[] brent = FunctionMinimization.brent(0.0, mu, scaled_rate_max, opt, TOL);
                opt.set(brent[0]); 
                opt_lh = brent[1]; // [1] ? vraiment ? [0] = x ; [1] = fx ?
            }
            {
                DeltaOptimization opt = new DeltaOptimization();
                double[] brent = FunctionMinimization.brent(getDeltaMin(), delta, scaled_rate_max, opt, TOL); 
                opt.set(brent[0]);
                opt_lh = brent[1]; 
            }
            //if (false)
            {
                LambdaOptimization opt = new LambdaOptimization();
                double[] brent = FunctionMinimization.brent(0.0, lambda, getLambdaMax(), opt, TOL);
                opt.set(brent[0]);
                opt_lh = brent[1]; 
            }
           return opt_lh;
           
        }
        
        /**
         * Maximum value of lambda in the optimization
         * 
         * @return max. lambda value
         */
        protected double getLambdaMax()
        {
            return delta;
        }
        
        /**
         * Maximum value of delta in the optimization
         * 
         * @return max. delta value
         */
        protected double getDeltaMin()
        {
            return lambda;
        }
      //  ... classes internes pour l'optimisation
        
    protected abstract class RateOptimization implements OneParameterFunction
    {
        protected RateOptimization()
        {
            //scanAllValues();
        }
        
        /**
         * Used to see whether line optimization is OK: are there multiple maxima?
         */
        protected void scanAllValues()
        {
            double o_mu = mu;
            double o_lambda = lambda; 
            double o_delta = delta;
            for (double rate = 2.0; rate>1e-5; rate /= 1.01)
            {
                set(rate);
                double lik = getLikelihood();
                System.out.println("PL.RO.sAV\t"+rate+"\t"+lik+"\t"+getClass().getSimpleName()+"\t"+ProfileLikelihood.this.paramString());
            }
            delta = o_delta;
            lambda = o_lambda;
            mu = o_mu;
        }
        
        @Override
        public double eval(double rate)
        {
            set(rate);
            double loglik = getNegativeLogLikelihood();
            double lik = getLikelihood();
            //System.out.println("#*PL."+getClass().getSimpleName()+".eval\t"+rate+"\tll "+loglik+"\tlik "+lik);
            return loglik;                   
        }
        
        public abstract void set(double rate);
    }
    
        
    class MuOptimization extends RateOptimization 
    {        
        public void set(double rate)
        {
            mu = rate;
        }
    }
    class LambdaOptimization extends RateOptimization
    {
       public void set(double rate)
        {
            lambda = rate;
        }
    }
    class DeltaOptimization extends RateOptimization
    {
         public void set(double rate)
        {
            delta = rate;
        }
    }
    
    @Override
    public String toString()
    {
        return getClass().getSimpleName()+"["+paramString()+"]";
    }
    
    protected String paramString()
    {
        return "\tmu "+mu+"\tlm "+lambda+"\tdl "+delta;
        
    }
    
}
