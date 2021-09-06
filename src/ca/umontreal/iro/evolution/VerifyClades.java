/*
 * MatchClades.java
 *
 * Created on May 13, 2007, 10:58 AM
 */

package ca.umontreal.iro.evolution;

/**
 * This class compares splits between a rooted and an unrooted phylogenetic tree 
 * (reference and query, respectively). 
 *
 * @author  csuros
 */


import java.util.HashSet;

import ca.umontreal.iro.banality.Verbose;
import ca.umontreal.iro.banality.BasicExecutable;
import ca.umontreal.iro.banality.StringSplit;

public class VerifyClades extends BasicExecutable
{
    
    /** Not instantiated from outside */
    private VerifyClades() 
    {
    }
    
    /**
     * This is essentially the true <q>main</q>, but needs instantiation
     */
    private void go(String[] args) throws Parser.ParseException, java.io.IOException
    {
        if (args.length!=3)
        {
            System.err.println("Call as $0 reference_tree query_tree clade_list");
            System.exit(99);
        }
        reportLaunch(args);
        
        String reference_tree_file= args[0];
        String query_tree_file   = args[1];
        String list_of_matched_clades = args[2].trim();

        TreeNode reference_root = Parser.readNewick(new java.io.FileReader(reference_tree_file));        
        TreeNode query_root     = Parser.readNewick(new java.io.FileReader(query_tree_file));
        
        TreeNode[] reference_nodes = reference_root.subtreeNodes();
        for (int j=0; j<reference_nodes.length; j++)
            reference_nodes[j].setId(j);
        TreeNode[] query_nodes = query_root.subtreeNodes();
        for (int j=0; j<query_nodes.length; j++)
            query_nodes[j].setId(j);

        String[] clade = StringSplit.splitAt(list_of_matched_clades,',');
        
        for (int clade_idx=0; clade_idx<clade.length; clade_idx++)
        {
            // processing clade[clade_idx]
            
            // find which node this is in the reference tree
            int ref_node_idx = -1;
            for (int j=0; j<reference_nodes.length; j++)
                if (clade[clade_idx].equals(reference_nodes[j].getName()))
                {
                    ref_node_idx=j;
                    break;
                }
            if (ref_node_idx==-1)
                throw new IllegalArgumentException("There is no node '"+clade[clade_idx]+"' in the reference tree.");
            
            // collect the terminal taxa under this node: it will be easier to check if a taxon belongs here or not
            // at the end, reference_clade is a set of taxon names (String) in the query clade, 
            HashSet<String> reference_clade = new HashSet<String>();
            aggregateClade(reference_nodes[ref_node_idx], reference_clade);
            int num_reference_taxa = reference_clade.size();
            int num_all_taxa = reference_root.subtreeLeaves().length;
            // compute the complement of the reference clade 
            //HashSet complement_reference_clade = new HashSet();
            //for (int j=0; j<num_all_taxa; j++)
            //    if (reference_nodes[j].isLeaf())
            //{
            //    String rname = reference_nodes[j].getName();
            //    complement_reference_clade.add(rname);
            //}
                            
            // node labeling for subtree content
            // each node in the query tree is labeled by TAXA_IN_CLADE, TAXA_IN_COMPLEMENT or TAXA_MIXED
            // TAXA_IN_CLADE: all the leaves in this subtree belong to the reference clade
            // TAXA_IN_COMPLEMENT : none of the leaves in the subtree belong to the reference clade
            // TAXA_MIXED: some but not all leaves in the subtree belong to the reference clade
            int[] query_node_label_i = new int[query_nodes.length];
            int num_present = 0;
            int num_query_leaves = query_root.subtreeLeaves().length;
            for (int qi=0; qi<query_nodes.length-1; qi++) // root is not labeled ...
            {
                if (query_nodes[qi].isLeaf())
                {
                    String qname = query_nodes[qi].getName();
                    if (reference_clade.contains(qname))
                    {
                        query_node_label_i[qi]=TAXA_IN_CLADE;
                        num_present++;
                        Verbose.message("VC.go found "+qname);
                    }
                    else 
                    {
                        query_node_label_i[qi]=TAXA_IN_COMPLEMENT;
                        Verbose.message("VC.go complement "+qname);
                    }
                } else
                {
                    int z = TAXA_UNKNOWN;
                    for (int cidx=0; cidx<query_nodes[qi].getNumChildren(); cidx++)
                    {
                        int child_label = query_node_label_i[query_nodes[qi].getChild(cidx).getId()];
                        if (cidx==0)
                            z = child_label;
                        else if (z!=child_label)
                            z = TAXA_MIXED;
                    }
                    query_node_label_i[qi]=z;
                }
            } // for qi
            if (num_present == 0 || num_present == num_query_leaves)
            {
                System.out.println(clade[clade_idx]+"\tabsent");
                continue;
            }
            
            // node labeling for the complement of the subtree
            // each node in the query tree is labeled by TAXA_IN_CLADE, TAXA_IN_COMPLEMENT or TAXA_MIXED
            // TAXA_IN_CLADE: all the leaves outside of this subtree belong to the reference clade
            // TAXA_IN_COMPLEMENT : none of the leaves outside of the subtree belong to the reference clade
            // TAXA_MIXED: some but not all leaves outside of the subtree belong to the reference clade
            int[] query_node_label_o = new int[query_nodes.length];
            for (int qi=query_nodes.length-2; qi>=0; qi--)
            {
                int z = TAXA_UNKNOWN;
                TreeNode qP = query_nodes[qi].getParent();
                if (!qP.isRoot())
                    z = query_node_label_o[qP.getId()];
                for (int cidx=0;cidx<qP.getNumChildren(); cidx++)
                {
                    int sibling_idx = qP.getChild(cidx).getId();
                    if (sibling_idx != qi)
                    {
                        if (z == TAXA_UNKNOWN)
                            z = query_node_label_i[sibling_idx];
                        else if (z != query_node_label_i[sibling_idx])
                            z = TAXA_MIXED;
                    }
                }
                query_node_label_o[qi]=z;
                Verbose.message("VC.go qnode "+z+" "+query_nodes[qi]);
            } // for qi
            
            // now verify if there is a node that is good for us
            int matched_clade_idx = -1;
            boolean matches_complement = false;
            int mismatch_evidence = -1;
            for (int qi=0; qi<query_nodes.length-1; qi++)
            {
                Verbose.message("VC.go "+qi+"\ti "+query_node_label_i[qi]+" o "+query_node_label_o[qi]+"\t"+query_nodes[qi].newickSubtree(false,false,false,false));
                boolean matches_clade = (query_node_label_i[qi]==TAXA_IN_CLADE && query_node_label_o[qi]==TAXA_IN_COMPLEMENT);
                matches_complement = (query_node_label_o[qi]==TAXA_IN_CLADE && query_node_label_i[qi]==TAXA_IN_COMPLEMENT);
                if (matches_clade || matches_complement)
                {
                    matched_clade_idx = qi;
                    break;
                } 
                if (query_node_label_i[qi] == TAXA_MIXED && query_node_label_o[qi] == TAXA_MIXED)
                {
                    mismatch_evidence = qi;
                    break;
                }
            }
            if (matched_clade_idx==-1)
                System.out.println(clade[clade_idx]+"\tmismatched\t"
                +"\t\t"+query_nodes[mismatch_evidence].newickSubtree());
            else
                System.out.println(clade[clade_idx]+"\tmatches\t"
                    +num_present+"/"+num_reference_taxa
                    +","+(num_query_leaves-num_present)+"/"+(num_all_taxa-num_reference_taxa)+"\t"
                    +(matches_complement?"complement/":"")
                    +query_nodes[matched_clade_idx].newickSubtree(false,false,false,false));
        } // for clade_idx
    }
    
    private static final int TAXA_UNKNOWN = 999;
    private static final int TAXA_IN_CLADE = 901;
    private static final int TAXA_IN_COMPLEMENT = 910;
    private static final int TAXA_MIXED = 911;
    
    /**
     * Auxiliary recursive function that constructs the set of leaves in a subtree
     */
    private void aggregateClade(TreeNode subtree_root, HashSet<String> clade_so_far)
    {
        if (subtree_root.isLeaf())
            clade_so_far.add(subtree_root.getName());
        else 
            for (int cidx=0; cidx<subtree_root.getNumChildren(); cidx++)
                aggregateClade(subtree_root.getChild(cidx), clade_so_far);
    }
    
    public static void main(String[] args)
    {
        Verbose.setVerbose(false);
                
        VerifyClades O = new VerifyClades();
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
        } 
        catch (Exception E)
        {
            die(E);
        }
    }
    
}
