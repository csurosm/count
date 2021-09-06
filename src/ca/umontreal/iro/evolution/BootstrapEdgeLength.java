
package ca.umontreal.iro.evolution;

/**
 * A utility to calculate bootstrap confidence intervals for edges 
 * of a phylogeny. Bootstrap trees are listed in a file 
 * (in Newick format, one after the other), and 
 * the main phylogeny is specified in a different file. 
 * Bootstrap trees are supposed to have the same shape as
 * the main phylogeny. 
 *
 * @author csuros
 */

import java.io.BufferedReader;
import java.util.Hashtable;

import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.Verbose;


public class BootstrapEdgeLength extends BasicExecutable
{
    private void go(String[] args) throws Exception
    {
        if (args.length != 2)
        {
            System.err.println("Call as java "+getClass().getCanonicalName()+" main_tree file-of-trees");
            System.exit(2008);
        }
        
        reportLaunch(args);       
        String main_tree_file = args[0];
        
        TreeNode main_tree = null;
        {
            BufferedReader BR = new java.io.BufferedReader(new java.io.FileReader(main_tree_file));
            main_tree = Parser.readNewick(BR);
            BR.close();
        }

        TreeNode original_root = main_tree;
        Verbose.message("BEL.go original "+main_tree.newickTree());
        TreeNode[] main_nodes = main_tree.subtreeNodes();
        
        Hashtable<String,TreeNode> original_names = new Hashtable<String,TreeNode>();
        Hashtable<TreeNode, Integer> original_idx = new Hashtable<TreeNode, Integer>();
        for (int node_idx=0; node_idx<main_nodes.length; node_idx++)
        {
            TreeNode N = main_nodes[node_idx];
            if (N.isLeaf())
                original_names.put(N.getName(), N);
            original_idx.put(N, new Integer(node_idx));
        }
        
        // read bootstrap trees
        String trees_file = args[1];
        TreeNode[] trees = null;
        {
            BufferedReader BR = new java.io.BufferedReader(new java.io.FileReader(trees_file));
            trees = Parser.readTrees(BR);
            BR.close();
        }
        int num_trees = trees.length;
        Verbose.message("BEL.go Read "+num_trees+" trees");
    
    
        double[][] lengths = new double[main_nodes.length][num_trees];
        
        for (int tree_idx=0; tree_idx<num_trees; tree_idx++)
        {
            TreeNode[] nodes = trees[tree_idx].subtreeNodes();
            Hashtable<TreeNode,Integer> node_map = new Hashtable<TreeNode,Integer>();
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                TreeNode N = nodes[node_idx];
                if (N.isLeaf())
                {
                    TreeNode mainN = original_names.get(N.getName());
                    node_map.put(N, original_idx.get(mainN));
                } else
                {
                    int child_index = node_map.get(N.getChild(0));
                    TreeNode main_child = main_nodes[child_index];
                    TreeNode main_parent = main_child.getParent();
                    node_map.put(N, original_idx.get(main_parent));
                }
            }
            
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                TreeNode N = nodes[node_idx];
                int original_node_idx = node_map.get(N).intValue();
                if (!N.isRoot())
                {
                    lengths[original_node_idx][tree_idx] = N.getLength();
                }
                if (Verbose.isVerbose() && tree_idx<10)
                {
                    Verbose.message("BEL.go tree "+tree_idx+"\tnode "+node_idx+"/"+N.newickName()+"\tmain "+original_node_idx+"/"+main_nodes[original_node_idx].newickName());
                }
            }
        }
        
        // header line
        System.out.print("Node");
        for (int tree_idx=0; tree_idx<trees.length; tree_idx++)
            System.out.print("\tbs"+(1+tree_idx));
        System.out.println();
        
        for (int node_idx=0; node_idx<main_nodes.length; node_idx++)
        {
            TreeNode N = main_nodes[node_idx];
            if (!N.isRoot())
            {
                System.out.print(N.newickName());
                for (int tree_idx=0; tree_idx<trees.length; tree_idx++)
                    System.out.print("\t"+lengths[node_idx][tree_idx]);
                System.out.println();
            }
        }
        
    }
    
    public static void main(String[] args)
    {
        Verbose.setVerbose(false);
        
        BootstrapEdgeLength O = new BootstrapEdgeLength(); 
        
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
