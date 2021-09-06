/*
 * ConsensusTree.java
 *
 * Created on August 18, 2007, 6:16 PM
 */

package ca.umontreal.iro.evolution;

/**
 * Implementation of linear-time methods for computing strict consensus 
 * or Robinson-Foulds distance from Day's paper (J. Classification 2:7-28, 1985)
 *
 * @author  csuros
 */

import java.util.Vector;
import java.util.Hashtable;

import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.Verbose;

import ca.umontreal.iro.evolution.Parser;

public class ConsensusTree extends BasicExecutable
{
    
    /**
     * Initializes the data structures with the first tree.
     * Leaves are supposed to have id's 0..n-1
     */
    public ConsensusTree(TreeNode root) 
    {
        setRoot(root);
    }
    
    private ConsensusTree()
    {}
    
    private void setRoot(TreeNode root)
    {
        computePostorder(root);
        computeSplits();
    }
    
    /**
     * Computes the strict consensus with a second tree
     * @param root2 root of the second tree
     * @return array of matching nodes (from tree 1)
     */
    public TreeNode[] consensus(TreeNode root2)
    {
        boolean[] matches = new boolean[N];
        int[] range = new int[2];
        range[0]=N+1;
        range[1]=-1;
        checkSplits(root2, range, matches);
        
        Vector<TreeNode> M=new Vector<TreeNode>();
        for (int j=0; j<N; j++)
            if (matches[j])
                M.add(postorder[splits[j][2]].v);
        
        return M.toArray(new TreeNode[0]);
    }

    
    /**
     * Computes the size of the strict consensus between a second tree and the one used at instantiation.
     *
     * @param root2 root of the second tree
     * @return number of interor nodes of tree 1 that have a match in tree 2
     */
    public int sizeConsensus(TreeNode root2)
    {
        return consensus(root2).length;
    }
    
    /**
     * Auxiliary function for performing a postorder traversal while verifying splits.
     * @param v subtree root
     */
    private int checkSplits(TreeNode v, int[] range, boolean[] matches)
    {
        if (v.isLeaf())
        {
            int vidx = ((Integer)leafOrder.get(v.getName())).intValue();
            range[0]=range[1]=vidx;
            return 1;
        } else 
        {
            int num_leaves_in_subtree = 0;
            for (int cidx=0; cidx<v.getNumChildren(); cidx++)
            {
                int[] subtree_range = new int[2];
                subtree_range[0]=N+1;
                subtree_range[1]=-1;
                num_leaves_in_subtree += checkSplits(v.getChild(cidx), subtree_range, matches);
                if (subtree_range[0]<range[0])
                    range[0]=subtree_range[0];
                if (subtree_range[1]>range[1])
                    range[1]=subtree_range[1];
            }
            //Verbose.message("CT.cS range "+range[0]+".."+range[1]+", n="+num_leaves_in_subtree+" // "+v);
            if (num_leaves_in_subtree == range[1]-range[0]+1)
            {
                // this might be a match
                if (splits[range[0]]!=null && splits[range[0]][0]==range[0] && splits[range[0]][1]==range[1])
                {
                    matches[range[0]]=true;
                    //Verbose.message("CT.cS match @ "+range[0]+".."+range[1]);
                }
                else if (splits[range[1]]!=null && splits[range[1]][0]==range[0] && splits[range[1]][1]==range[1])
                {
                    matches[range[1]]=true;
                    //Verbose.message("CT.cS match @ "+range[0]+".."+range[1]);
                }
            }
            return num_leaves_in_subtree;
        }
    }
    
    /**
     * Mapping (String->Integer) of leaf names to their postorder index in the first tree
     * This table correspnds to the third column X[*,3] in Day's paper
     */
    private Hashtable<String,Integer> leafOrder;
    
    /**
     * Splits: entry i is a 2-element array [L,R] for L..R range of leaves in the first tree
     * Using Day's encoding [L,R] is either at splits[L] or splits[R]
     * This array corresponds to columns 1 and 2 X[*,1] X[*,2] in Day's paper 
     * The idea is to store [L,R] split at index R iff its subtree root is a rightmost child. 
     * The third element of splits[i] is the index into postorder[] for the subtree root.
     */
    private int[][] splits;

    /**
     * Fills up the leafOrder and splits
     */ 
    private void computeSplits()
    {
       leafOrder = new Hashtable<String,Integer>();
       splits = new int[N][];
       
       int leaf_index=0;
       // index into postorder
       int parcours=0;
       int R=-1;
       while (parcours<N+M)
       {
           Subtree S = postorder[parcours];
           if (S.v.isLeaf())
           {
               leafOrder.put(S.v.getName(),new Integer(leaf_index));
               //Verbose.message("CT.cS leaf "+S.v.getName()+"->"+leaf_index);
               R = leaf_index;
               leaf_index++;
               parcours++;
           } else
           {
               // find out left-hand side for the current subtree
               TreeNode Lv = postorder[parcours-S.w].v;
               int L = ((Integer) leafOrder.get(Lv.getName())).intValue();
               int loc = L;
               if (S.v.isRoot() || postorder[parcours+1].v.isLeaf())
                   loc = R;
               splits[loc]=new int[3];
               splits[loc][0]=L;
               splits[loc][1]=R;
               splits[loc][2]=parcours;
               //Verbose.message("CT.cS split "+loc+" -> ["+L+", "+R+"]");
               parcours++;
           }
           //for (int j=0; j<splits.length; j++)
           //    System.out.println("["+j+"]"
           //     +(splits[j]==null
           //         ?"null"
           //         :(splits[j][0]+","+splits[j][1])));
       }
    }

    /** 
     * Array of nodes in a postorder traversal
     */
    private Subtree[] postorder;

    /** 
     * Number of leaves
     */
    private int N;
    
    /**
     * Number of interior nodes
     */
    private int M;

    /**
     * Computes leaf counts, interior node counts, and postorder traversal
     * @param root root for the tree
     */
    private void computePostorder(TreeNode root)
    {
        N=M=0;
        Vector<Subtree> V = new Vector<Subtree>();
        computePostorder(root, V);
        postorder=V.toArray(new Subtree[0]);
    }

    /**
     * The recursive auxiliary function for initPostorder()
     * @param v subtree root
     * @param V Vector collecting nodes in postorder  traversal
     * @return size of the subtree 
     */
    private int computePostorder(TreeNode v, Vector<Subtree> V)
    {
        int w=0;
        for (int cidx=0; cidx<v.getNumChildren(); cidx++)
            w+=computePostorder(v.getChild(cidx),V);
        V.add(new Subtree(v,w));
        if (v.isLeaf())
            N++;
        else
            M++;
        return w+1;
    }
    
    /**
     * A simple class to represent (node, subtree size) pairs
     */
    private class Subtree
    {
        Subtree(TreeNode v, int w)
        {
            this.v = v;
            this.w = w;
        }
        
        TreeNode v;
        int w;
    }

    /**
     * Finds the first node in the tree
     * 
     * @param v root of the tree
     * @return the leaf node that is alphabeticaly first among v's descendants
     */
    private static TreeNode first_leaf_alphabetically(TreeNode v)
    {
        if (v.isLeaf())
            return v;
        else
        {
             TreeNode first = null;
            for (int cidx=0; cidx<v.getNumChildren(); cidx++)
            {
                TreeNode fc = first_leaf_alphabetically(v.getChild(cidx));
                if (first == null || fc.getName().compareTo(first.getName())<0)
                    first = fc;
            }
            return first;
        }
    }
    
    /**
     * Reroots a tree at the parent of the leaf that is alphabetically first
     *
     * @param root old root
     * @return new root
     */
    public static TreeNode reroot(TreeNode root)
    {
        //Verbose.message("CT.rer before "+root.newickTree());
        TreeNode first = first_leaf_alphabetically(root);
        TreeNode P = first.getParent();
        P.reroot();
        //Verbose.message("CT.rer after "+P.newickTree());
        return P;
    }
    
    private void go(String[] args) throws Exception
    {
        if (args.length !=2)
        {
            System.err.println("Call as $0 tree1 tree2");
            System.exit(2007);
        }
        String tree_file1 = args[0];
        String tree_file2 = args[1];
        TreeNode R = Parser.readNewick(new java.io.FileReader(tree_file1));
        R = reroot(R);
        setRoot(R);
        
        TreeNode R2 = Parser.readNewick(new java.io.FileReader(tree_file2));
        R2 = reroot(R2);
        
        TreeNode[] C = consensus(R2);
        for (int j=0; j<C.length; j++)
            System.out.println("Match "+(1+j)+"\t"+C[j].newickSubtree());
        
        int t1 = M; 
        int t2 = R2.numNodes()-N;
        Verbose.message("CT.go t1="+t1+" t2="+t2);
        System.out.println("Robinsons-Foulds distance "+(t1+t2-2*C.length)/2);
    }
    
    /**
     * Computes the Robinson-Foulds distance between the tree used at initialization.
     * It is asusmed that child nodes are ordered alphabetically: call reroot() to achieve that
     * 
     * @param root1 root of the first tree
     * @param root2 root of the second tree
     * @return Robinson-Foulds distance between the two trees 
     */
    public static int RobinsonFoulds(TreeNode root1, TreeNode root2)
    {
        //root1 = reroot(root1);
        ConsensusTree CT = new ConsensusTree(root1);
        //root2 = reroot(root2);
        TreeNode[] c = CT.consensus(root2);
        int t1 = CT.M;
        int t2 = root2.numNodes()-CT.N;
        int rf = (t1+t2-2*c.length)/2;
        return rf;
    }
    
    public static void main(String[] args)
    {
        Verbose.setVerbose(false);
        ConsensusTree O = new ConsensusTree();
        int num_switches = 0;
        try {
            while (args.length>2*num_switches && args[2*num_switches].startsWith("-")){
                String arg_switch = args[2*num_switches].substring(1);
                if (arg_switch.equals("h"))
                    O.go(new String[0]); // will throw an Exception
                if (args.length==2*num_switches+1)
                    throw new IllegalArgumentException("Missing argument for switch "+args[2*num_switches]);
                String arg_value = args[2*num_switches+1];
                if (arg_switch.equals("v"))
                    Verbose.setVerbose(arg_value.equals("true"));
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
