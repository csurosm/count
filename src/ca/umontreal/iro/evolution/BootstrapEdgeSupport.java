
package ca.umontreal.iro.evolution;

/**
 * A utility to calculate bootstrap support for edges 
 * of a phylogeny. Bootstrap trees are listed in a file 
 * (in Newick format, one after the other), and 
 * the main phylogeny is specified in a different file. 
 * 
 * @author csuros
 */


import java.io.BufferedReader;

import java.util.Hashtable;

import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.Verbose;
        
public class BootstrapEdgeSupport extends BasicExecutable
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
        
        // save original parent-child relationsips
        Hashtable<TreeNode,TreeNode> original_parents = new Hashtable<TreeNode, TreeNode>();
        {
            TreeNode.Traversal traversal = main_tree.getTraversal();
            TreeNode[] node = traversal.getDFT();
            for (int node_idx=0;node_idx<node.length; node_idx++)
            {
                TreeNode N = node[node_idx];
                if (!N.isRoot())
                {
                    TreeNode P = N.getParent();
                    original_parents.put(N, P);
                }
            }
        }
        TreeNode original_root = main_tree;
        
        Verbose.message("BES.go original "+main_tree.newickTree());
        // reroot by alphabetically first leaf
        main_tree = ConsensusTree.reroot(main_tree);
        Verbose.message("BES.go original (rerooted) "+main_tree.newickTree());
        
        // read bootstrap trees
        String trees_file = args[1];
        TreeNode[] trees = null;
        {
            BufferedReader BR = new java.io.BufferedReader(new java.io.FileReader(trees_file));
            TreeNode[] original = Parser.readTrees(BR);
            trees = new TreeNode[original.length];
            for (int tree_idx=0; tree_idx<original.length; tree_idx++)
            {
                trees[tree_idx] = ConsensusTree.reroot(original[tree_idx]);
            }

            BR.close();
        }
        int num_trees = trees.length;
        Verbose.message("BES.go Read "+num_trees+" trees");
        

        TreeNode[] main_nodes = main_tree.getTraversal().getDFT();
        
        // initialize counters
        Hashtable<TreeNode, Integer> edge_support = new Hashtable<TreeNode, Integer>();
        {
            for (int node_idx=0; node_idx<main_nodes.length; node_idx++)
            {
                TreeNode N = main_nodes[node_idx];
                if (!N.isLeaf())
                    edge_support.put(N, new Integer(0));
            }
        }
                    
        ConsensusTree C = new ConsensusTree(main_tree);
        for (int tree_idx=0; tree_idx<num_trees; tree_idx++)
        {
            TreeNode other_tree = trees[tree_idx];
            TreeNode[] matching_nodes = C.consensus(other_tree);
            int match_counter = 0;
            for (int node_idx=0; node_idx<matching_nodes.length; node_idx++)
            {
                TreeNode N = matching_nodes[node_idx];
                if (!N.isLeaf())
                {
                    int cnt = edge_support.get(N).intValue();
                    edge_support.put(N, new Integer(cnt+1));
                    Verbose.message("BES.go tree "+tree_idx+"\t"+N.newickName()+"\tmatch ("+(1+cnt)+")\t"+N);
                    match_counter++;
                }
            }
            Verbose.message("BES.go tree "+tree_idx+"\tsupports "+match_counter+" edges.");
        }
        
        for (int node_idx=0; node_idx<main_nodes.length; node_idx++)
        {
            TreeNode N = main_nodes[node_idx];
            //if (!N.isLeaf())
            {
                int cnt = -1;
                TreeNode cP = N.getParent();
                TreeNode oP = original_parents.get(N);
                Verbose.message("BES.go node "+N+"\tcP "+cP+"\toP "+oP);
                if (edge_support.containsKey(N))
                    cnt = edge_support.get(N).intValue();
                TreeNode original_edge = (cP==oP?N:cP);
                if (cnt != -1 && original_edge!=null)
                    System.out.println(original_edge.newickName()+"\t"+(cP==oP)+"\t"+N.newickName()+"\t"+cnt+"\t"+N.newickSubtree());
            }
        }
    }
    
    public static void main(String[] args)
    {
        Verbose.setVerbose(false);
        
        BootstrapEdgeSupport O = new BootstrapEdgeSupport(); 
        
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
