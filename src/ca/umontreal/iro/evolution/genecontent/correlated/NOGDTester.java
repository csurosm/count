
package ca.umontreal.iro.evolution.genecontent.correlated;

import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;

import ca.umontreal.iro.evolution.TreeNode;
import ca.umontreal.iro.evolution.Parser;
import ca.umontreal.iro.evolution.genecontent.SimpleOccurrenceTable;
import ca.umontreal.iro.banality.Functions;


/**
 * 
 * Testing NOGD computations
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class NOGDTester 
{
 
    private NOGDTester(){}
    
    private void testXY12(String[] args) throws IOException, Parser.ParseException
    {
        if (args.length!=5)
        {
            System.err.println("Call as java "+getClass().getName()+" tree table X Y Y'");
            System.exit(2009);
        }
        String tree_file = args[0];
        String table_file = args[1];
        int family_X = Integer.parseInt(args[2]);
        int family_Y1 = Integer.parseInt(args[3]);
        int family_Y2 = Integer.parseInt(args[4]);
        
        TreeNode root = Parser.readNewick(new BufferedReader(new FileReader(tree_file)));
        TreeNode[] leaves = root.getTraversal().getLeaves();

        SimpleOccurrenceTable tbl = new SimpleOccurrenceTable(leaves);
        tbl.readTable(new FileReader(table_file));        
        
        int[] profile = NOGD.getTripleProfile(tbl, family_X, family_Y1, family_Y2);
        System.out.print("# Profile: X="+tbl.getFamilyName(family_X)+"\tY1="+tbl.getFamilyName(family_Y1)+"\tY2="+tbl.getFamilyName(family_Y2));
        
        for (int i=0; i<profile.length; i++)
            System.out.print("\t"+leaves[i].newickName()+":"+Integer.toBinaryString(8+profile[i]).substring(1));
        System.out.println();
        
        
        double max_rate = ProfileLikelihood.getMaxScaledRate(root.getTraversal().getDFT());
        double mu0 = max_rate * 0.001;
        double lambda0 = max_rate * 0.0001;
        double nu0 = max_rate * 0.01;
        double delta0 = max_rate * 0.04;
        
        UncorrelatedProfileLikelihood P_uncorr = new UncorrelatedProfileLikelihood(root, profile, mu0, lambda0, delta0);
        double l_uncorr = P_uncorr.optimize();
        CorrelatedProfileLikelihood P_corr = new CorrelatedProfileLikelihood(root, profile, P_uncorr.mu, P_uncorr.lambda, P_uncorr.delta, P_uncorr.delta);
        double l_corr = P_corr.optimize();
        UnconditionalProfileLikelihood P_uncond = new UnconditionalProfileLikelihood(root, profile, P_corr.mu, P_corr.lambda, P_corr.nu, P_corr.delta);
        double l_uncond = P_uncond.optimize();
        //MixtureProfileLikelihood P_mixt = new MixtureProfileLikelihood(root, profile, P_corr.mu, P_corr.lambda, P_corr.nu, P_corr.delta, 0.5);
        //double l_mixt = P_mixt.optimize();
        
        double chi_uncond = (l_uncorr>l_uncond?Functions.Chi_square_tail(1, 2.0*(l_uncorr-l_uncond)):1.0);
        double chi_corr = (l_uncorr>l_corr?Functions.Chi_square_tail(1, 2.0*(l_uncorr-l_corr)):1.0);
        //double chi_mixt = (l_uncorr>l_mixt?Functions.Chi_square_tail(2, 2.0*(l_uncorr-l_mixt)):1.0);
        
        System.out.println("# Model\tLog-likelihood\tmu\tlambda\tdelta\tnu\ttheta\tchi");
        System.out.println("UU\t"+(-l_uncorr)+"\t"+P_uncorr.mu+"\t"+P_uncorr.lambda+"\t"+P_uncorr.delta);
        System.out.println("UC\t"+(-l_uncond)+"\t"+P_uncond.mu+"\t"+P_uncond.lambda+"\t"+P_uncond.delta+"\t"+P_uncond.nu+"\t0\t"+chi_uncond);
        System.out.println("CC\t"+(-l_corr)+"\t"+P_corr.mu+"\t"+P_corr.lambda+"\t"+P_corr.delta+"\t"+P_corr.nu+"\t1\t"+chi_corr);
       // System.out.println("MC\t"+(-l_mixt)+"\t"+P_mixt.mu+"\t"+P_mixt.lambda+"\t"+P_mixt.delta+"\t"+P_mixt.nu+"\t"+P_mixt.theta+"\t"+chi_mixt);
        
    }
    
    public static void main(String[] args) throws Exception
    {
        NOGDTester O = new NOGDTester();
        O.testXY12(args);
        
        
    }

    
    
}
