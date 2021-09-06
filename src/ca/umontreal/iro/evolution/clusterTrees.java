package ca.umontreal.iro.evolution;

/**
 * A utility for counting how many different topologies 
 * are in a file of trees (Newick format).
 * 
 * @author csuros
 */

import java.io.BufferedReader;

import java.util.Vector;

import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.Verbose;

public class clusterTrees extends BasicExecutable
{
    
    private void go(String[] args) throws Exception
    {
        if (args.length != 1)
        {
            System.err.println("Call as java "+getClass().getCanonicalName()+" file-of-trees");
            System.exit(2008);
        }
        
        reportLaunch(args);        
        String trees_file = args[0];
        Vector<TreeNode> trees = new Vector<TreeNode>();
        
        BufferedReader BR = new java.io.BufferedReader(new java.io.FileReader(trees_file));
        //{
        //    while (true)
        //    {
        //        TreeNode R = Parser.readNewick(BR);
        //        if (R==null)
        //            break;
        //        R = ConsensusTree.reroot(R);
        //        trees.add(R);
        //    }
        //}
        
            
        //String line = null;
        //do
        //{
        //    line = BR.readLine();
        //    if (line != null)
        //    {
        //        if (line.indexOf(';')==-1)
        //            continue;
        //        StringReader lineR = new StringReader(line);
        //        TreeNode R = Parser.readNewick(lineR);
        //        R = ConsensusTree.reroot(R);
        //        trees.add(R);
        //   }
        //} while (line != null);
        TreeNode[] trees_in_the_file = Parser.readTrees(BR);
        BR.close();

        int num_trees = trees_in_the_file.length;
        for (int tree_idx=0; tree_idx<num_trees; tree_idx++)
        {
            trees.add(ConsensusTree.reroot(trees_in_the_file[tree_idx]));
        }

        Verbose.message("cT.go Read "+num_trees+" trees");
        
        int[] canonical_topology=new int[num_trees];
        //canonical_topology[0] =0;

        Vector<Integer> canonical_trees = new Vector<Integer>();
        //canonical_trees.add(new Integer(0));
        int[] canonical_tree_count = new int[num_trees];
        
        for (int tree_idx=0; tree_idx<num_trees; tree_idx++)
        {
            TreeNode root = trees.get(tree_idx);
            canonical_topology[tree_idx] = -1;
            for (int other_idx=0; other_idx<canonical_trees.size(); other_idx++)
            {
                int other_tree_idx = canonical_trees.get(other_idx).intValue();
                TreeNode other_root = trees.get(other_tree_idx);
                int rf = ConsensusTree.RobinsonFoulds(root, other_root);
                if (rf==0)
                {
                    canonical_topology[tree_idx] = other_tree_idx;
                    break;
                }
            }
            if (canonical_topology[tree_idx]==-1)
            {
                canonical_trees.add(new Integer(tree_idx));
                canonical_topology[tree_idx]=tree_idx;
            }
            canonical_tree_count[canonical_topology[tree_idx]]++;
        }
        
        Verbose.message("cT.go Found "+canonical_trees.size()+" different topologies");
        for (int i=0; i<trees.size(); i++)
            if (canonical_topology[i]==i)
            {
                TreeNode R = trees.get(i);
                System.out.println("[tree"+(i)+" appears "+canonical_tree_count[i]+" times]\t"+R.newickTree());
            }

    }
    
    public static void main(String[] args)
    {
        Verbose.setVerbose(false);
        
        clusterTrees O = new clusterTrees(); //InclusionExclusionComputation());
        
        int num_switches = 0;
        try {
            while (args.length>2*num_switches && args[2*num_switches].startsWith("-"))
            {
                String arg_switch = args[2*num_switches].substring(1);
                if (arg_switch.equals("h"))
                    O.go(new String[0]); // will throw an Exception
                if (args.length==2*num_switches+1)
                    throw new IllegalArgumentException("Missing argument for switch "+args[2*num_switches]);
                String arg_value = args[2*num_switches+1];
                if (arg_switch.equals("v"))
                {
                    Verbose.setVerbose(arg_value.equals("true"));
                } 
                
                else 
                    throw new IllegalArgumentException("Switch not recognized: '"+args[2*num_switches]+"'");
                    
                num_switches++;
            }
            
            String[] rest=new String[args.length-2*num_switches];
            for (int j=0; j<rest.length; j++)
                rest[j]=args[2*num_switches+j];
            O.go(rest);
        } catch (Exception E)
        {
            die(E);
        }        
        
    }

}
