package ca.umontreal.iro.evolution.genecontent.correlated;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import java.util.Hashtable;
import java.util.Random;

import ca.umontreal.iro.evolution.TreeNode;
import ca.umontreal.iro.evolution.Parser;

import ca.umontreal.iro.banality.Functions;

/**
 *
 * @author csuros
 */
public class Simulator 
{
    private Simulator(){}

    public Simulator(TreeNode root)
    {
        setRoot(root);
    }
    
    private void setRoot(TreeNode root)
    {
        RND= new Random();
        this.root = root;
        nodes = root.getTraversal().getDFT();
        Hashtable<TreeNode, Integer> node_indexes = new Hashtable<TreeNode, Integer>();
        for (int node_idx = 0; node_idx<nodes.length; node_idx++)
        {
            TreeNode N = nodes[node_idx];
            node_indexes.put(N, node_idx);
        }
        node_parent = new int[nodes.length];
        leaf_index = new int[nodes.length];
        int leaf_idx = 0;
        for (int node_idx=0; node_idx<nodes.length; node_idx++)
        {
            TreeNode N = nodes[node_idx];
            if (N.isRoot())
                node_parent[node_idx] = -1;
            else
            {
                TreeNode P = N.getParent();
                int parent_idx = node_indexes.get(P);
                node_parent[node_idx] = parent_idx;
            }
            if (N.isLeaf())
            {
                leaf_index[node_idx] = leaf_idx;
                leaf_idx++;
            } else
                leaf_index[node_idx]=-1;
        }
        
        leaves = root.getTraversal().getLeaves();
    }
    
    private Random RND;
    private TreeNode root;
    private TreeNode[] nodes;
    private TreeNode[] leaves;
    private int[] node_parent;
    private int[] leaf_index;
    
    public int[] randomCorrelatedProfile(double mu, double lambda, double nu, double delta)
    {
        return randomProfile(mu, lambda, nu, delta, true, true);
    }
    
    public int[ ] randomUnconditionalProfile(double mu, double lambda, double nu, double delta)
    {
        return randomProfile(mu, lambda, nu, delta, true, false);
    }
    
    public int[ ] randomUncorrelatedProfile(double mu, double lambda, double delta)
    {
        return randomProfile(mu, lambda, delta, delta, false, false);
    }
    
    private enum ProfileType {CORRELATED, UNCONDITIONAL, UNCORRELATED};
    
    private int[] randomProfile(double mu, double lambda, double nu, double delta, boolean is_conditional, boolean is_correlated)
    {
        int[] label = new int[nodes.length];
        double p00 = delta*nu;
        double p01 = delta*lambda;
        double p10 = p01;
        double p11 = lambda*lambda;
        {
            double scale = p00+p01+p10+p11;
            p00 /= scale;
            p01 /= scale;
            p10 /= scale;
            p11 /= scale;
        }
        int root_idx = label.length-1;
        {
            double p = RND.nextDouble();
            if (p<p00)
                label[root_idx] = 4;
            else 
            {
                p -= p00;
                if (p<p01)
                    label[root_idx] = 5;
                else
                {
                    p-= p01;
                    if (p<p10)
                        label[root_idx] = 6;
                    else 
                        label[root_idx] = 7;
                }
            }
        }
        for (int node_idx = root_idx-1; node_idx>=0; node_idx--)
        {
            TreeNode N = nodes[node_idx];
            int pidx = node_parent[node_idx];
            int parent_label = label[pidx];
            label[node_idx] = simulateBranch(parent_label, N, mu, lambda, nu, delta, is_conditional, is_correlated);
        }
        
        int[] retval = new int[leaves.length];
        for (int node_idx=0; node_idx<nodes.length; node_idx++)
        {
            TreeNode N = nodes[node_idx];
            if (N.isLeaf())
            {
                int lidx = leaf_index[node_idx];
                retval[lidx] = label[node_idx];
            }
        }

        return retval;
    }
    
    private int simulateBranch(int parent_label, TreeNode N, double mu, double lambda, double nu, double delta, boolean is_conditional, boolean is_correlated)
    {
        boolean X_is_present = (parent_label & 4) !=0;
        boolean Y1_is_present = (parent_label & 2) != 0;
        boolean Y2_is_present = (parent_label & 1) != 0;
        
        double time_left = N.getLength();
        while(time_left > 0.0)
        {
            Event E = null;
            {
                double rate = nu;
                if (Y1_is_present)
                {
                    if (Y2_is_present && is_correlated)
                        if (!is_conditional || X_is_present)
                            rate = delta;
                } else
                    rate = lambda;
                E = randomEvent(EventType.Y1, rate);
            }
            {
                Event E2 = null;
                double rate = nu;
                if (Y2_is_present)
                {
                    if (Y1_is_present && is_correlated)
                        if (!is_conditional || X_is_present)
                            rate = delta;
                } else
                    rate = lambda;
                E2 = randomEvent(EventType.Y2,rate);
                if (E2.time<E.time)
                    E = E2;
            }
            if (X_is_present)
            {
                Event Ex = randomEvent(EventType.X,mu);
                if (Ex.time < E.time)
                    E = Ex;
            }
            if (E.time < time_left)
            {
                if (E.type == EventType.Y1)
                    Y1_is_present = !Y1_is_present;
                else if (E.type == EventType.Y2)
                    Y2_is_present = !Y2_is_present;
                else if (E.type == EventType.X)
                    X_is_present = false;
                
                //System.out.println("#*S.sB "+N.newickName()+"\t"+E);
            }
            time_left -= E.time;
        }
        int child_label = 0;
        if (X_is_present)
            child_label = child_label | 4;
        if (Y1_is_present)
            child_label = child_label | 2;
        if (Y2_is_present)
            child_label = child_label | 1;
        return child_label;
    }
    
    private double randomEventTime(double rate)
    {
        return -Math.log(RND.nextDouble())/rate;
    }
    
    
    public enum EventType {X, Y1, Y2};
    
    private Event randomEvent(EventType type, double rate)
    {
        return new Event(type, randomEventTime(rate));
    }
    
    private static class Event 
    {
        
        private Event(EventType type, double time)
        {
            this.type = type;
            this.time = time;
        }
        private EventType type;
        private double time;
        
        @Override
        public String toString()
        {
            return type.name()+" @ "+time;
        }
    }


    private int sampleCorrelatedProfile(double mu, double lambda, double nu, double delta)
    {
        int[] profile = randomCorrelatedProfile(mu, lambda, nu, delta);
        boolean all_zero = true;
        char[] pattern = new char[profile.length];
        for (int i=0; i<profile.length; i++)
        {
            if (profile[i]!=0)
                all_zero = false;
            pattern[i] = Integer.toString(profile[i]).charAt(0);
        }
        System.out.println("#*S.sCP "+(new String(pattern)));
        if (all_zero)
            return 0;
        
        double max_rate = ProfileLikelihood.getMaxScaledRate(nodes);
        double mu0 = max_rate * 0.001;
        double lambda0 = max_rate * 0.0001;
        double nu0 = max_rate * 0.01;
        double delta0 = max_rate * 0.04;
        
        ProfileLikelihood P_true = new CorrelatedProfileLikelihood(root, profile, mu, lambda, nu, delta);    
        double l_true = P_true.getNegativeLogLikelihood();

        UncorrelatedProfileLikelihood P_uncorr = new UncorrelatedProfileLikelihood(root, profile, mu0, lambda0, delta0);
        double l_uncorr = P_uncorr.optimize();
        CorrelatedProfileLikelihood P_corr = new CorrelatedProfileLikelihood(root, profile, mu, lambda, nu, delta);
        double l_corr = P_corr.optimize();
        double chi_corr = (l_uncorr>l_corr?Functions.Chi_square_tail(1, 2.0*(l_uncorr-l_corr)):1.0);
        
        System.out.println(l_true+"\t"+l_uncorr+"\t"+l_corr+"\t"+chi_corr+"\t"+P_corr.mu+"\t"+P_corr.lambda+"\t"+P_corr.nu+"\t"+P_corr.delta+"\t// "+P_uncorr.mu+"\t"+P_uncorr.lambda+"\t"+P_uncorr.delta);
        return 1;
    }
        
    private int sampleUnconditionalProfile(double mu, double lambda, double nu, double delta)
    {
        int[] profile = randomUnconditionalProfile(mu, lambda, nu, delta);
        boolean all_zero = true;
        for (int i=0; i<profile.length; i++)
            if (profile[i]!=0)
            {
                all_zero = false;
                break;
            }
        if (all_zero)
            return 0;
        
        double max_rate = ProfileLikelihood.getMaxScaledRate(nodes);
        double mu0 = max_rate * 0.001;
        double lambda0 = max_rate * 0.0001;
        double nu0 = max_rate * 0.01;
        double delta0 = max_rate * 0.04;
        
        ProfileLikelihood P_true = new UnconditionalProfileLikelihood(root, profile, mu, lambda, nu, delta);    
        double l_true = P_true.getNegativeLogLikelihood();

        UncorrelatedProfileLikelihood P_uncorr = new UncorrelatedProfileLikelihood(root, profile, mu0, lambda0, delta0);
        double l_uncorr = P_uncorr.optimize();

        UnconditionalProfileLikelihood P_uncond = new UnconditionalProfileLikelihood(root, profile, P_uncorr.mu, P_uncorr.lambda, P_uncorr.delta, P_uncorr.delta);
        double l_uncond = P_uncond.optimize();
        double chi_uncond = (l_uncorr>l_uncond?Functions.Chi_square_tail(1, 2.0*(l_uncorr-l_uncond)):1.0);
        
        System.out.println(l_true+"\t"+l_uncorr+"\t"+l_uncond+"\t"+chi_uncond+"\t"+P_uncond.mu+"\t"+P_uncond.lambda+"\t"+P_uncond.nu+"\t"+P_uncond.delta+"\t// "+P_uncorr.mu+"\t"+P_uncorr.lambda+"\t"+P_uncorr.delta);
        return 1;
    }

    private int sampleUncorrelatedProfile(double mu, double lambda, double delta)
    {
        int[] profile = randomUncorrelatedProfile(mu, lambda, delta);
        boolean all_zero = true;
        char[] pattern = new char[profile.length];
        for (int i=0; i<profile.length; i++)
        {
            if (profile[i]!=0)
                all_zero = false;
            pattern[i] = Integer.toString(profile[i]).charAt(0);
        }
        System.out.println("#*S.sUncorr "+(new String(pattern)));

        if (all_zero)
            return 0;
        
        double max_rate = ProfileLikelihood.getMaxScaledRate(nodes);
        double mu0 = max_rate * 0.001;
        double lambda0 = max_rate * 0.0001;
        double nu0 = max_rate * 0.01;
        double delta0 = max_rate * 0.04;
        
        ProfileLikelihood P_true = new UncorrelatedProfileLikelihood(root, profile, mu, lambda, delta);    
        double l_true = P_true.getNegativeLogLikelihood();

        UncorrelatedProfileLikelihood P_uncorr = new UncorrelatedProfileLikelihood(root, profile, mu0, lambda0, delta0);
        double l_uncorr = P_uncorr.optimize();

        System.out.println(l_true+"\t"+l_uncorr+"\t"+P_uncorr.mu+"\t"+P_uncorr.lambda+"\t"+P_uncorr.delta);
        return 1;
    }    
    
    public void run(String[] args) throws IOException, Parser.ParseException
    {
        int num_switches = 0;
        while (args.length>2*num_switches && args[2*num_switches].startsWith("-"))
        {
            String arg_switch = args[2*num_switches].substring(1);
            if (args.length==2*num_switches+1)
                throw new IllegalArgumentException("Missing argument for switch "+args[2*num_switches]);
            String arg_value = args[2*num_switches+1];
            
            if ("type".equals(arg_switch))
            {
                if ("correlated".equals(arg_value))
                    this.profile_type = ProfileType.CORRELATED;
                else if ("unconditional".equals(arg_value))
                    this.profile_type = ProfileType.UNCONDITIONAL;
                else if ("uncorrelated".equals(arg_value))
                    this.profile_type = ProfileType.UNCORRELATED;
                else throw new IllegalArgumentException("-type must be uncorrelated, unconditional, or correlated");
            } else if ("mu".equals(arg_switch))
                mu = Double.parseDouble(arg_value);
            else if ("lambda".equals(arg_switch))
                lambda = Double.parseDouble(arg_value);
            else if ("nu".equals(arg_switch))
                nu = Double.parseDouble(arg_value);
            else if ("delta".equals(arg_switch))
                delta = Double.parseDouble(arg_value);
            else if ("rep".equals(arg_switch))
                repeats = Integer.parseInt(arg_value);

            num_switches++;
        }
        
        
        if (args.length!=2*num_switches+1)
        {
            System.err.println("Call as java "+getClass().getName()+" [switches] tree ");
            System.exit(2009);
        }
        String tree_file = args[2*num_switches];
        
        setRoot(Parser.readNewick(new BufferedReader(new FileReader(tree_file))));
        
        run();
    }
    
    public void run()
    {
        if (profile_type.equals(ProfileType.CORRELATED))
        {
            System.out.println("# Ltrue\tLuncorr\tLcorr\tchi\tmu("+mu+")\tlm("+lambda+")\tnu("+nu+")\tdl("+delta+")");
            for (int rep=0; rep<repeats; )
                rep += sampleCorrelatedProfile(mu, lambda, nu, delta);
        } else if (profile_type.equals(ProfileType.UNCONDITIONAL))
        {
            System.out.println("# Ltrue\tLuncorr\tLuncond\tchi\tmu("+mu+")\tlm("+lambda+")\tnu("+nu+")\tdl("+delta+")");
            for (int rep=0; rep<repeats; )
                rep += sampleUnconditionalProfile(mu, lambda, nu, delta);
        } else if (profile_type.equals(ProfileType.UNCORRELATED))
        {
            System.out.println("# Ltrue\tLuncorr\tmu("+mu+")\tlm("+lambda+")\tdl("+delta+")");
            for (int rep=0; rep<repeats; )
                rep += sampleUncorrelatedProfile(mu, lambda, delta);
            
        }
    } 
    
    private double mu = 0.6;
    private double lambda = 0.2;
    private double delta = 3.0;
    private double nu = 0.6;
    private ProfileType profile_type = ProfileType.CORRELATED;
    private int repeats = 10;

    public static void main(String[] args) throws Exception
    {
        Simulator S = new Simulator();
        S.run(args);
    }
}
