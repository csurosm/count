
package ca.umontreal.iro.evolution.genecontent.correlated;

import ca.umontreal.iro.evolution.TreeNode;

/**
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s and Louis Philippe B. Bouvrette
 */
public class UncorrelatedProfileLikelihood extends ProfileLikelihood{

    public UncorrelatedProfileLikelihood(TreeNode root, int[] profile, double mu, double lambda, double delta)
    {
            super(root, profile, mu, lambda,delta);
    }

    @Override
    public double getLikelihood(int[] profile)
    {
            return NOGD.getUncorrelatedLikelihood(profile, nodes, mu, lambda, delta);
    }
    
}
