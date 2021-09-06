/*
 * TreeWithRates.java
 *
 * Created on April 18, 2008, 3:33 PM
 *
 */

package ca.umontreal.iro.evolution.genecontent;


import ca.umontreal.iro.evolution.TreeNode;
import java.io.IOException;
import java.util.Hashtable;
import java.util.HashMap;
import java.util.Map;

/**
 * This class is useful for representing a tree that is used often in
 * depth-first traversal algorithms.
 *
 * @author  csuros
 */
public class TreeWithRates 
{
    
    /**
     * Initializes the tree with the given root.
     * The tree topology is assumed to be immutable, so that the 
     * postfix traversal of the nodes does not change.  
     * 
     * @param root the root of the tree topology
     */ 
    public TreeWithRates(NodeWithRates root) 
    {
        setTree(root);
    }
    
    /**
     * Initializes the tree with the given root.
     * The tree topology is assumed to be immutable, so that the 
     * postfix traversal of the nodes does not change.  
     * 
     * @param root the root of the tree topology
     */ 
    public TreeWithRates(TreeNode root) 
    {
        NodeWithRates my_root = NodeWithRates.copyTree(root);
        
        setTree(my_root);
    }
    
    private NodeWithRates[] depth_first_traversal;
    private NodeWithRates[] leaves;
    
    private NodeWithRates root;
    
    /**
     * (Re)sets the tree.
     * Also sets the node ids so that they are indices into the 
     * depth-first traversal array and the leaves array.
     * 
     * @param fixed_tree_root new root for the tree
     */ 
    public void setTree(NodeWithRates fixed_tree_root)
    {
        root=fixed_tree_root;
        depth_first_traversal = root.getDFT();
        node_indexes = null;

        // now reorder so that leaves are at the front
        {
            NodeWithRates[] reordered = new NodeWithRates[depth_first_traversal.length];
            int idx=0;
            for (int node_idx=0; node_idx<depth_first_traversal.length; node_idx++)
                if (depth_first_traversal[node_idx].isLeaf())
                {
                    reordered[idx]=depth_first_traversal[node_idx];
                    idx++;
                }
            for (int node_idx=0; node_idx<depth_first_traversal.length; node_idx++)
                if (!depth_first_traversal[node_idx].isLeaf())
                {
                    reordered[idx]=depth_first_traversal[node_idx];
                    idx++;
                }
            depth_first_traversal=reordered;
        }
        
        parent_idx = new int[depth_first_traversal.length-1];
        
        leaves = root.getLeaves();
        for (int i=depth_first_traversal.length-1; i>=0; i--)
        {
            NodeWithRates N = depth_first_traversal[i];
            N.setId(i);
            if (!N.isRoot())
            {
                TreeNode P = N.getParent();
                parent_idx[i] = P.getId();
            }
        }
        
        child_indexes = new int[depth_first_traversal.length][];
        for (int i=0; i<depth_first_traversal.length; i++)
        {
            NodeWithRates N = depth_first_traversal[i];
            int num_children = N.getNumChildren();
            child_indexes[i] = new int[num_children];
            for (int ci=0; ci<num_children; ci++)
                child_indexes[i][ci]=N.getChild(ci).getId();
        }
        
    }
    
    private int[] parent_idx;
    private int[][] child_indexes;

    /**
     * Returns the index of the parent in the traversal 
     * 
     * @param node_idx Index of a non-root node in the traversal 
     * @return index of the parent
     */
    public int getParentIndex(int node_idx)
    {
        return parent_idx[node_idx];
    }
    
    /**
     * Returns the index of a child in the traversal 
     * 
     * @param node_idx index of a node in the traversal 
     * @param child_position position (left,right,...) of the child at the parent, can be 0..nc-1 where nc is number of the children 
     * @return index of the child (for using with getNode())
     */
    public int getChildIndex(int node_idx, int child_position)
    {
        return child_indexes[node_idx][child_position];
    }
    
    
    /**
     * Returns an array of the tree nodes, suitable for depth-first traversal. 
     * Tree leaves are always in the front of the array, and tree nodes have an id
     * that gives their index in this array.
     *
     * @return array of tree nodes
     */
    public NodeWithRates[] getDFT()
    {
        return depth_first_traversal;
    }
    
    /**
     * Number of terminal nodes in the tree
     * @return number of treeleaves
     */
    public int getNumLeaves()
    {
        return leaves.length;
    }
    
    /**
     * Number of nodes in the tree
     * @return number of nodes
     */
    public int getNumNodes()
    {
        return depth_first_traversal.length;
    }
    
    /**
     * Number of edges in the tree: number of nodes -1
     * @return numbr of edges
     */
    public int getNumEdges()
    {
        return getNumNodes()-1;
    }
    
    /**
     * Returns an array of the tree leaves. 
     * Tree leaves have an id
     * that gives their index in this array.
     *
     * @return array of tree leaves
     */
    public NodeWithRates[] getLeaves()
    {
        return leaves;
    }

    /**
     * Returns the root of the tree.
     * @return the tree root
     */
    public NodeWithRates getRoot()
    {
        return root;
    }

    /**
     * Returns the node with the given index in the depth-first traversal array.
     * @param node_idx index of a node: 0 to numNodes()-1
     * @return the corresponding tree node
     */
    public NodeWithRates getNode(int node_idx)
    {
        return depth_first_traversal[node_idx];
    }
    
    public int getNodeIndex(NodeWithRates N)
    {
        if (node_indexes==null)
        {
            node_indexes = new Hashtable<NodeWithRates,Integer>();
            for (int i=0; i<depth_first_traversal.length; i++)
                node_indexes.put(depth_first_traversal[i], i);
        }
        return node_indexes.get(N);
    }
    
    private Hashtable<NodeWithRates,Integer> node_indexes;
    
    
    
    
    public int[] getSubtreeSizes()
    {
        int[] size = new int[depth_first_traversal.length];
        for (int node_idx=0; node_idx<depth_first_traversal.length; node_idx++)
        {
            NodeWithRates N = depth_first_traversal[node_idx];
            if (N.isLeaf())
                size[node_idx]=1;
            else
            {
                int num_children = N.getNumChildren();
                for (int ci=0; ci<num_children; ci++)
                {
                    int child_idx = getChildIndex(node_idx,ci);
                    size[node_idx]+=size[child_idx];
                }
            }
        }
        return size;
    }
    
    /**
     * @param len sequence length
     * @param alpha alphabet size
     */
    public void reportCompression(int len, int alpha)
    {
        int exponent_cutoff = 0;
        int A_k=1;
        while (A_k<len)
        {
            exponent_cutoff++;
            A_k *= alpha;
        }
        int[] A_exp = new int[exponent_cutoff];
        A_exp[0]=1;
        for (int j=1; j<A_exp.length; j++)
            A_exp[j]=A_exp[j-1]*alpha;
        
        int[] size = getSubtreeSizes();
        int total = 0;
        for (int node_idx=0; node_idx<size.length; node_idx++)
        {
            int s = 0;
            if (size[node_idx]>=exponent_cutoff)
                s = len;
            else 
                s = A_exp[size[node_idx]];
            
            total += s;
            NodeWithRates N = depth_first_traversal[node_idx];
            System.out.println("# TWR compression "+node_idx+"/"+N.getTaxonName()+"\t"+size[node_idx]+"\t"+s);
        }
        System.out.println("# TWR compression total\t"+total+"\t// nl="+depth_first_traversal.length*len);
    }
    
    /**
     * Whether there is an edge with positive duplication rate on it.
     * @return true if there is at least one edge where duplication rate is positive.
     */
    public boolean hasDuplication()
    {
        for (int node_idx=0; node_idx<depth_first_traversal.length; node_idx++)
        {
            NodeWithRates N = depth_first_traversal[node_idx];
            if (!N.isRoot())
            {
                if (N.getDuplicationRate()>0.)
                    return true;
            }
        }
        return  false;
    }

    /**
     * Whether there is an edge with positive gain (transfer) rate on it.
     * @return true if there is at least one edge where gain/transfer rate is positive.
     */
    public boolean hasGain()
    {
        for (int node_idx=0; node_idx<depth_first_traversal.length; node_idx++)
        {
            NodeWithRates N = depth_first_traversal[node_idx];
            if (!N.isRoot())
            {
                if (N.getTransferRate()>0.)
                    return true;
            }
        }
        return  false;
    }

    /**
     * Computes a mapping between nodes of trees with the same
     * topology (same leaf names, same branching structure),
     * but represented by TreeNode and TreeWithRates.
     *
     * @param original_tree_root alternative tree representation
     * @return mapping from this tree's nodes to the other tree's nodes
     */
    public Map<NodeWithRates,TreeNode> getTreeNodeMapping(TreeNode original_tree_root)
    {
        Map<NodeWithRates,TreeNode> node_mapping = new HashMap<NodeWithRates,TreeNode>();
        // all nodes
        TreeNode[] original_nodes = original_tree_root.getTraversal().getDFT();
        NodeWithRates[] new_nodes = getDFT();
        // map the leaves
        Map<String,TreeNode> original_leaves = new HashMap<String,TreeNode>();
        for (int node_idx=0; node_idx<original_nodes.length; node_idx++)
        {
            TreeNode N = original_nodes[node_idx];
            if (N.isLeaf())
                original_leaves.put(N.getName(), N);
        }
        for (int node_idx=0; node_idx<new_nodes.length; node_idx++)
        {
            NodeWithRates N = new_nodes[node_idx];
            if (N.isLeaf())
            {
                String name = N.getName();
                TreeNode oN = original_leaves.get(name);
                node_mapping.put(N, oN);
            } else
            {
                // an inner node
                // the trees have the same topology: it's enough to look at one child
                NodeWithRates C = (NodeWithRates)N.getChild(0);
                TreeNode oC = node_mapping.get(C);
                TreeNode oN = oC.getParent();
                node_mapping.put(N, oN);
            }
        }
        return node_mapping;
    }

    /**
     * Average loss rate in this tree
     *
     * @param weight_by_length whether a weighted or arithmetic average is wanted
     * @return average rate
     */
    public double averageLossRate(boolean weight_by_length)
    {
        int num_edges = getNumEdges();
        // average rates across edges to get a common homogeneous rate
        double common_loss_rate = 0.0;

        double total_edge_length = 0.0;
        for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
        {
            NodeWithRates N = getNode(edge_idx);
            double Nlen = (weight_by_length?N.getLength():1.0);
            double Nloss = Nlen*N.getLossRate();
            common_loss_rate += Nloss;
            total_edge_length += Nlen;
        }

        common_loss_rate /= total_edge_length;
        return common_loss_rate;
    }

    /**
     * Average duplication rate in this tree
     *
     * @param weight_by_length whether a weighted or arithmetic average is wanted
     * @return average rate
     */
    public double averageDuplicationRate(boolean weight_by_length)
    {
        int num_edges = getNumEdges();
        // average rates across edges to get a common homogeneous rate
        double common_rate = 0.0;

        double total_edge_length = 0.0;
        for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
        {
            NodeWithRates N = getNode(edge_idx);
            double Nlen = (weight_by_length?N.getLength():1.0);
            double Nloss = Nlen*N.getDuplicationRate();
            common_rate += Nloss;
            total_edge_length += Nlen;
        }

        common_rate /= total_edge_length;
        return common_rate;
    }

    /**
     * Average gain rate in this tree
     *
     * @param weight_by_length whether a weighted or arithmetic average is wanted
     * @return average rate
     */
    public double averageTransferRate(boolean weight_by_length)
    {
        int num_edges = getNumEdges();
        // average rates across edges to get a common homogeneous rate
        double common_rate = 0.0;

        double total_edge_length = 0.0;
        for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
        {
            NodeWithRates N = getNode(edge_idx);
            double Nlen = (weight_by_length?N.getLength():1.0);
            double Nloss = Nlen*N.getTransferRate();
            common_rate += Nloss;
            total_edge_length += Nlen;
        }

        common_rate /= total_edge_length;
        return common_rate;
    }

    
    /**
     * Test code: computes compression bound for this tree
     * Call as $0 <tree> <length> <alphabet size>
     * @param args command-line arguments
     * @throws java.io.IOException if there is a problem while reading the files 
     * @throws ca.umontreal.iro.evolution.Parser.ParseException if the tree is not in proper Newick format
     */
    public static void main(String[] args) throws java.io.IOException, ca.umontreal.iro.evolution.Parser.ParseException
    {
        if (args.length != 3)
        {
            System.err.println("Call as $0 tree sequence_length alphabet_size");
            System.exit(2008);
            
        }
        
        String tree_file = args[0];
        int seq_len = Integer.parseInt(args[1]);
        int alphabet_size = Integer.parseInt(args[2]);
        
        TreeWithRates main_tree = new TreeWithRates(ca.umontreal.iro.evolution.Parser.readNewick(new java.io.FileReader(tree_file)));
        
        main_tree.reportCompression(seq_len,alphabet_size);
    }
    
}
