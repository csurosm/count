package ca.umontreal.iro.evolution.genecontent.correlated;

import ca.umontreal.iro.evolution.TreeNode;

/**
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 */
public class UnconditionalProfileLikelihood extends CorrelatedProfileLikelihood
{
    public UnconditionalProfileLikelihood(TreeNode root, int[] profile, double mu, double lambda, double nu, double delta)
    {
        super(root, profile, mu, lambda, nu, delta);
    }

    @Override
    public double getLikelihood(int[] profile)
    {
        return NOGD.getUnconditionalLikelihood(profile, nodes, mu, lambda, nu, delta);
    }
}

        
        