package ca.umontreal.iro.evolution.genecontent;

import java.util.Map;
import java.util.Arrays;

import ca.umontreal.iro.evolution.TreeNode;
import ca.umontreal.iro.matek.DiscreteDistribution;

/**
 *
 * @since Tue 27 Sep 2011
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class HomogeneousRateVariation extends RateVariation
{
    private TreeNode scaled_root;
    private double common_edge_scaling;
    private double common_loss_rate;
    private Map<NodeWithRates,TreeNode> node_mapping;

    protected HomogeneousRateVariation(){ super();}

    public HomogeneousRateVariation(TreeNode scaled_root, TreeWithRates main_tree,
            DiscreteDistribution root_prior,
            int duplication_rate_categories,
            int loss_rate_categories,
            int transfer_rate_categories,
            int edge_length_categories)
    {
        super(main_tree, root_prior, duplication_rate_categories, loss_rate_categories, transfer_rate_categories, edge_length_categories);
        setScaledTree(scaled_root);
    }

    public HomogeneousRateVariation(TreeNode scaled_root, RateVariation model)
    {
        this(scaled_root, model.main_tree,model.getRootPrior(),
                model.getNumDuplicationRateGammaCategories(),
                model.getNumLossRateGammaCategories(),
                model.getNumTransferRateGammaCategories(),
                model.getNumEdgeLengthGammaCategories());
    }

    private void setScaledTree(TreeNode root)
    {
        scaled_root = root;

        // figure out common scaling
        TreeNode[] scaled_nodes = scaled_root.getTraversal().getDFT();
        double[] edge_lengths = new double[scaled_nodes.length-1];
        int num_edges = main_tree.getNumEdges();
        for (int node_idx=0; node_idx<num_edges; ++node_idx) // skip root in the last posiiton of the table
            edge_lengths[node_idx] = scaled_nodes[node_idx].getLength();
        Arrays.sort(edge_lengths);

        double len_max = edge_lengths[edge_lengths.length-1];
        double len_med = edge_lengths[edge_lengths.length/2];

        double range = ML.MAX_EDGE_LENGTH/ML.MIN_EDGE_LENGTH;
        assert (range>1.0);
        double len_min = len_max/range;
        for (int j=0; edge_lengths[j]<len_min; ++j) edge_lengths[j] = len_min;
        len_min = edge_lengths[0];

        double scale_max = ML.MAX_EDGE_LENGTH/len_max;
        double scale_med = NodeWithRates.DEFAULT_EDGE_LENGTH/len_med;
        common_edge_scaling = Math.min(scale_max, scale_med);
        // round to the nearest power of 10
        common_edge_scaling = Math.pow(10.0, (int)Math.log10(common_edge_scaling));

        // average rates across edges to get a common homogeneous rate
        common_loss_rate = main_tree.averageLossRate(true);
//        System.out.println("#*HRV.sST length scale "+common_edge_scaling+"\tloss "+common_loss_rate);

        node_mapping = main_tree.getTreeNodeMapping(scaled_root);

        for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
        {
            NodeWithRates N = main_tree.getNode(edge_idx);
            TreeNode original_node = node_mapping.get(N);
            double original_length = original_node.getLength();
            double new_length = original_length*common_edge_scaling;
            super.setEdgeLength(edge_idx, new_length);
            double Nlen = N.getLength();
            double Nloss = Nlen*common_loss_rate;
            double Ngain = Nlen*N.getTransferRate();
            double Ndupl = Nlen*N.getDuplicationRate();
            super.setNodeLossRate(edge_idx, Nloss/new_length);
            super.setNodeDuplicationRate(edge_idx, Ndupl/new_length);
            super.setNodeTransferRate(edge_idx, Ngain/new_length);
            super.setEdgeLength(edge_idx, new_length);
//            System.out.println("#*HRV.sST node "+edge_idx+"/"+N.newickName()+"\t"+N);
        }

        updateClassTrees();
    }

    /**
     * Sets the loss rate for all edges
     *
     * @param node_idx ignored
     * @param rate common rate for all nodes
     */
    @Override
    public void setNodeLossRate(int node_idx, double rate)
    {

        this.common_loss_rate = rate;
        int num_edges = main_tree.getNumEdges();
        for (int edge_idx=0; edge_idx<num_edges; ++edge_idx)
            super.setNodeLossRate(edge_idx, common_loss_rate);
    }


    @Override
    public HomogeneousRateVariation sameModelForDifferentTree(TreeWithRates other_tree)
    {
        HomogeneousRateVariation R = new HomogeneousRateVariation(scaled_root, other_tree, root_prior_distribution,
                    getNumDuplicationRateGammaCategories(),
                    getNumLossRateGammaCategories(),
                    getNumTransferRateGammaCategories(),
                    getNumEdgeLengthGammaCategories());
        copyVariationParameters(R);
        return R;
    }

    /**
     * Edge length cannot be reset in this model
     * @param node_idx ignored
     * @param length ignored
     */
    @Override
    public void setEdgeLength(int node_idx, double length)
    {
        // no effect
    }

    @Override
    public boolean hasLineageSpecificLoss()
    {
        return false;
    }


    
}
