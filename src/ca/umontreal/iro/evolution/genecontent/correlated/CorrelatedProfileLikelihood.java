
package ca.umontreal.iro.evolution.genecontent.correlated;

import ca.umontreal.iro.evolution.TreeNode;
import ca.umontreal.iro.matek.FunctionMinimization;

/**
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s and Louis Philippe B. Bouvrette
 */
public class CorrelatedProfileLikelihood extends ProfileLikelihood{
        
    public CorrelatedProfileLikelihood(TreeNode root, int[] profile, double mu, double lambda, double nu, double delta)
        {
                super(root, profile, mu, lambda, delta);
                this.nu = nu;
        }
        protected double nu;

        @Override
        public double getLikelihood(int[] profile)
        {
                return NOGD.getCorrelatedLikelihood(profile, nodes, mu, lambda, nu, delta);
        }
        
        @Override
        public double optimizeRates() // redéfinition pour optimiser nu aussi
        {   
              double opt_lh = super.optimizeRates();
              NuOptimization opt = new NuOptimization();
              double[] brent = FunctionMinimization.brent(lambda, nu, getNuMax(), opt, TOL);
              opt.set(brent[0]);
              opt_lh = brent[1];
              return opt_lh;
        }

        @Override
        protected double getLambdaMax()
        {
            return Math.min(nu,delta);
        }
        
        //protected double getDeltaMin()
        //{
        //    return nu;
        //}
        
        protected double getNuMax()
        {
            return scaled_rate_max;
        }
        
        protected class NuOptimization extends RateOptimization
        {
            @Override
            protected void scanAllValues()
            {
                double o_nu = nu;
                super.scanAllValues();
                nu = o_nu;
            }
            public void set(double rate)
            {
                nu = rate;
            }
        } 
    
        
    protected String paramString()
    {
        return super.paramString()+"\tnu "+nu;
    }
}
