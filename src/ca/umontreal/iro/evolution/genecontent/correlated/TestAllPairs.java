package ca.umontreal.iro.evolution.genecontent.correlated;

import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;

import java.util.Hashtable;

import ca.umontreal.iro.evolution.TreeNode;
import ca.umontreal.iro.evolution.Parser;
import ca.umontreal.iro.evolution.genecontent.SimpleOccurrenceTable;
import ca.umontreal.iro.banality.Functions;

/**
 *
 * @author csuros
 */
public class TestAllPairs implements Runnable
{
    private TestAllPairs(){}
    
    private TreeNode root;
    private SimpleOccurrenceTable table;
    private int X_idx;
    
    @Override
    public void run()
    {
        String X_name = table.getFamilyName(X_idx);
        double max_rate = ProfileLikelihood.getMaxScaledRate(root.getTraversal().getDFT());
        int[] px = table.getSizes(X_idx);
        int presentx = 0;
        for (int i=0; i<px.length; i++)
            if (px[i]>0) presentx++;

        int min_present1 = 2*presentx/3;
        int min_present2 = 4;
        
        Hashtable<String,StringBuffer> seen = new Hashtable<String, StringBuffer>();
        
        System.out.println("# X\tY1\tY2\tcomp\tamp_uncond\tamp_corr\tl_uncorr\tl_uncond\tl_corr\tchi_uncond\tchi_corr\tpattern");
        
        int num_families = table.getNumFamilies();
        for (int Y1_idx=0; Y1_idx<num_families; Y1_idx++)
            if (Y1_idx != X_idx)
            {
                int[] p1 = table.getSizes(Y1_idx);
                int present1 = 0;
                for (int i=0; i<p1.length && present1 < min_present1; i++)
                    if (p1[i]>0) present1++;
                if (present1 >= min_present1)
                {
                    String Y1_name = table.getFamilyName(Y1_idx);
                    for (int Y2_idx=Y1_idx+1; Y2_idx<num_families; Y2_idx++)
                        if (Y2_idx != X_idx)
                        {
                            int[] p2 = table.getSizes(Y2_idx);
                            int present2 = 0;
                            for (int i=0; i<p2.length && present2 < min_present2; i++)
                                if (p2[i]>0) present2++;
                            if (present2>=min_present2)
                            {
                                String Y2_name = table.getFamilyName(Y2_idx);
                                int[] profile = NOGD.getTripleProfile(table, X_idx, Y1_idx, Y2_idx);

                                int n01 = 0;
                                int n10 = 0;
                                int n11 = 0;
                                for (int i=0;i<profile.length; i++)
                                {
                                    int t = profile[i] & 3;
                                    if (t==1) n01++;
                                    else if (t==2) n10++;
                                    else if (t==3) n11++;
                                } 
                                double comp = (n01+n10-n11)/(n01+n10+n11+0.0);
                                if (comp>0.5)
                                {
                                    StringBuffer line_start = new StringBuffer(X_name);
                                    line_start.append('\t');
                                    line_start.append(Y1_name);
                                    line_start.append('\t');
                                    line_start.append(Y2_name);
                                    //line.append('\t');
                                    //for (int i=0; i<profile.length; i++)
                                    //{
                                    //    if (i!=0)
                                    //        line.append(".");
                                    //    line.append(Integer.toBinaryString(8+profile[i]).substring(1));
                                    //}
                                    line_start.append('\t');
                                    line_start.append(comp);
                                    
                                    char[] pattern = new char[profile.length];
                                    for (int i=0; i<pattern.length; i++)
                                    {
                                        int z = profile[i] & 3;
                                        if (z==0) pattern[i] = '0';
                                        else if (z==1) pattern[i] = '1';
                                        else if (z==2) pattern[i] = '2';
                                        else if (z==3) pattern[i] = '3';
                                    }
                                    String pattern_key = new String(pattern);
                                    
                                    if (seen.containsKey(pattern_key))
                                        System.out.println(line_start+"\t"+seen.get(pattern_key)+"\t...");
                                    else
                                    {
                                        StringBuffer line = null;

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

                                        double amp1 = 10.0*Math.log10(P_uncond.delta/P_uncond.nu);
                                        line.append('\t');
                                        line.append(amp1);

                                        double amp2 = 10.0*Math.log10(P_corr.delta/P_corr.nu);
                                        line.append('\t');
                                        line.append(amp2);

                                        double chi_uncond = (l_uncorr>l_uncond?Functions.Chi_square_tail(1, 2.0*(l_uncorr-l_uncond)):1.0);
                                        double chi_corr = (l_uncorr>l_corr?Functions.Chi_square_tail(1, 2.0*(l_uncorr-l_corr)):1.0);

                                        line.append('\t');
                                        line.append(-l_uncorr);
                                        line.append('\t');
                                        line.append(-l_uncond);
                                        line.append('\t');
                                        line.append(-l_corr);

                                        line.append('\t');
                                        line.append(chi_uncond);
                                        line.append('\t');
                                        line.append(chi_corr);
                                        
                                        line.append('\t');
                                        line.append(pattern_key);

                                        if (chi_uncond<0.01)
                                            line.append("\tUcorr **** "+Integer.toString((int)(amp1+0.5)));

                                        if (chi_corr<0.01)
                                            line.append("\tCcorr **** "+Integer.toString((int)(amp2+0.5)));


                                        seen.put(pattern_key, line);
                                        for (int i=0; i<pattern.length; i++)
                                        {
                                            int c = pattern[i];
                                            if (c=='2') pattern[i] = '1';
                                            else if (c=='1') pattern[i] = '2';
                                        }
                                        seen.put(new String(pattern), line);

                                        System.out.println(line_start+"\t"+line);
                                    }
                                }
                            }
                        }
                }
            }
    }
    
    private void run(String[] args) throws IOException, Parser.ParseException
    {
        if (args.length!=3)
        {
            System.err.println("Call as java "+getClass().getName()+" tree table X");
            System.exit(2009);
        }
        String tree_file = args[0];
        String table_file = args[1];
        X_idx = Integer.parseInt(args[2]);
        
        root = Parser.readNewick(new BufferedReader(new FileReader(tree_file)));
        TreeNode[] leaves = root.getTraversal().getLeaves();

        table = new SimpleOccurrenceTable(leaves);
        table.readTable(new FileReader(table_file));    
        
        run();
    }
    
    public static void main(String[] args) throws Exception
    {
        TestAllPairs T = new TestAllPairs();
        T.run(args);
    }

}
