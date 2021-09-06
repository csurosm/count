/*
 * RateVariation.java
 *
 * Created on April 22, 2008, 11:26 AM
 *
 */

package ca.umontreal.iro.evolution.genecontent;


import java.io.BufferedReader;
import java.io.Reader;
import java.io.IOException;

import ca.umontreal.iro.banality.Verbose;
import ca.umontreal.iro.banality.BasicExecutable;

import ca.umontreal.iro.matek.DiscreteDistribution;
import ca.umontreal.iro.matek.Poisson;
import ca.umontreal.iro.matek.ShiftedGeometric;
import ca.umontreal.iro.matek.PointDistribution;
import ca.umontreal.iro.matek.NegativeBinomial;

import ca.umontreal.iro.evolution.DiscreteGamma;

/**
 * Rate variation model: gamma + invariant rate factors
 * for edge length, loss, duplication and transfer.
 * It is probably not a good idea to use edge length variation
 * together with some rate variation.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class RateVariation extends BasicExecutable
{
    /**
     * Instantiates a rate-variation model. In such a model, 
     * a cross-product of discrete Gamma categories is used for
     * rates (duplication,loss,transfer) and branch lengths.
     * For parameters on which variation is not imposed, the number of categories should be 1.
     * @param main_tree tree for lineage-specific rates
     * @param root_prior prior size distribution at root
     * @param duplication_rate_categories number of duplication rate categories for rate variation across sites
     * @param loss_rate_categories number of loss rate categories for rate variation across sites
     * @param transfer_rate_categories number of transfer rate categories for rate variation across sites
     * @param edge_length_categories number of edge_length_categories for rate variation across sites
     */
    public RateVariation(TreeWithRates main_tree,
            DiscreteDistribution root_prior,
            int duplication_rate_categories, 
            int loss_rate_categories, 
            int transfer_rate_categories, 
            int edge_length_categories) 
    {
        this.main_tree = main_tree;
        
        this.root_prior_distribution = root_prior;
        setNumberOfDiscreteCategories(
            duplication_rate_categories, 
            loss_rate_categories, 
            transfer_rate_categories, 
            edge_length_categories);
    }

    /**
     * Produces a RateVariation model that has the same
     * rate variation across families, and
     * the same root distribution.
     *
     * @param other_tree another phylogeny
     * @return a rate variation model 
     */
    public RateVariation sameModelForDifferentTree(TreeWithRates other_tree)
    {
        RateVariation R = new RateVariation(other_tree, root_prior_distribution, 
                    getNumDuplicationRateGammaCategories(), 
                    getNumLossRateGammaCategories(),
                    getNumTransferRateGammaCategories(),
                    getNumEdgeLengthGammaCategories());
        copyVariationParameters(R);
        return R;
    }

    protected final void copyVariationParameters(RateVariation target)
    {
        target.setDuplicationForbidden(this.duplication_forbidden);
        target.setTransferForbidden(this.transfer_forbidden);
        target.setDuplicationRateAlpha(this.duplication_rate_alpha);
        target.setLossRateAlpha(this.loss_rate_alpha);
        target.setTransferRateAlpha(this.transfer_rate_alpha);
        target.setEdgeLengthAlpha(this.edge_length_alpha);
    }

    /**
     * Sets the proportions of sites where the invariant model (0 factor) applies instead of Gamma variation
     * 
     * @param duplication_forbidden proportion of families with no duplication (between 0.0 and 1.0)
     * @param loss_forbidden proportion of families with no loss (between 0.0 and 1.0)
     * @param transfer_forbidden proportion of families with no transfer (between 0.0 and 1.0)
     */
    public void setInvariantFractions(double duplication_forbidden, double loss_forbidden, double transfer_forbidden)
    {
        setDuplicationForbidden(duplication_forbidden);
        setLossForbidden(loss_forbidden);
        setTransferForbidden(transfer_forbidden);
    }
    
    private void setDuplicationForbidden(double duplication_forbidden)
    {
        this.duplication_forbidden = duplication_forbidden;
    }
    
    private void setLossForbidden(double loss_forbidden)
    {
        this.loss_forbidden = loss_forbidden;
    }
    
    private void setTransferForbidden(double transfer_forbidden)
    {
        this.transfer_forbidden = transfer_forbidden;
    }
    
    /**
     * Sets the number of categories for duplication rate, loss rate, transfer rate and edge length (in this oder of the 
     * arguments). For parameters on which variation is not imposed, the number of categories should be 1. 
     *
     * @param duplication_rate_categories a positive number, use 1 for no rate variation
     * @param loss_rate_categories a positive number, use 1 for no rate variation
     * @param transfer_rate_categories a positive number, use 1 for no rate variation
     * @param edge_length_categories a positive number, use 1 for no rate variation
     */
    public void setNumberOfDiscreteCategories(
            int duplication_rate_categories, 
            int loss_rate_categories, 
            int transfer_rate_categories, 
            int edge_length_categories)
    {
        setNumDuplicationRateCategories(duplication_rate_categories);
        setNumLossRateCategories(loss_rate_categories);
        setNumTransferRateCategories(transfer_rate_categories);
        setNumEdgeLengthCategories(edge_length_categories);
        
        initDataStructures();
    }

    protected RateVariation(){}
    
    private double[] duplication_rate_categories; 
    private double[] loss_rate_categories; 
    private double[] transfer_rate_categories; 
    private double[] edge_length_categories; 
    
    private double duplication_rate_alpha=1.0;
    private double loss_rate_alpha=1.0;
    private double transfer_rate_alpha=1.0;
    private double edge_length_alpha=1.0;
    
    private double duplication_forbidden = 0.0;
    private double loss_forbidden = 0.0;
    private double transfer_forbidden = 0.0;

    private boolean has_duplication = true;
    private boolean has_gain = true;

    /**
     * Main underlying rate tree specifying edge-specific rates and lengths.
     */
    protected TreeWithRates main_tree;

    /**
     * Rate trees for categories.
     */
    protected TreeWithRates[] rate_tree;
    
    protected DiscreteDistribution root_prior_distribution;
        
    private void initDataStructures()
    {
        int num_classes = getNumClasses();
        rate_tree = new TreeWithRates[num_classes];
        
        NodeWithRates main_root = main_tree.getRoot();
        for (int cidx=0; cidx<num_classes; cidx++)
        {
            NodeWithRates class_root = NodeWithRates.copyTree(main_root);
            TreeWithRates class_tree = new TreeWithRates(class_root);
            rate_tree[cidx] = class_tree;
        }
        updateClassTrees();
        this.setRatesAllowed();
    }
    
    protected final void updateClassTrees()
    {
        int num_classes = getNumClasses();
        NodeWithRates[] main_nodes = main_tree.getDFT();
        
        for (int cidx=0; cidx<num_classes; cidx++)
        {
            TreeWithRates class_tree = rate_tree[cidx];

            NodeWithRates[] nodes = class_tree.getDFT();
            double xd = duplication_rate_categories[getDuplicationRateCategory(cidx)];
            double xl = loss_rate_categories[getLossRateCategory(cidx)];
            double xt = transfer_rate_categories[getTransferRateCategory(cidx)];
            double xe = edge_length_categories[getEdgeLengthCategory(cidx)];

            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                NodeWithRates N = nodes[node_idx];
                if (!N.isRoot())
                {
                    NodeWithRates protoN = main_nodes[node_idx];
                    N.setDuplicationRate(xd*protoN.getDuplicationRate());
                    N.setLossRate(xl*protoN.getLossRate());
                    N.setTransferRate(xt*protoN.getTransferRate());
                    N.setLength(xe*protoN.getLength());
                }
            }
        }
    }

    public void setDuplicationAllowed(boolean b){this.has_duplication = b;}
    public void setGainAllowed(boolean b){this.has_gain = b;}
    public void setRatesAllowed()
    {
        NodeWithRates[] nodes = getMainTree().getDFT();
        boolean hd = false;
        boolean hg = false;
        for (int node_idx=0; node_idx<nodes.length; ++node_idx)
        {
            NodeWithRates N = nodes[node_idx];
            if (!N.isRoot())
            {
                double len = N.getLength();
                double drate = len*N.getDuplicationRate();
                double grate = len*N.getTransferRate();
                hd |= (drate*len>0.0);
                hg |= (grate*len>0.0);
            }
        }
        setDuplicationAllowed(hd);
        setGainAllowed(hg);
    }

    public boolean hasDuplication(){return has_duplication;}
    public boolean hasGain(){return has_gain;}

    public boolean hasLineageSpecificDuplication()
    {
        double common_dup = 0.0;
        NodeWithRates[] nodes = getMainTree().getDFT();
        for (int node_idx=0; node_idx<nodes.length && !Double.isNaN(common_dup); ++node_idx)
        {
            NodeWithRates N = nodes[node_idx];
            if (!N.isRoot())
            {
                if (node_idx==0)
                    common_dup = N.getDuplicationRate();
                else if (N.getDuplicationRate() != common_dup)
                    common_dup = Double.NaN;
            }
        }
        return Double.isNaN(common_dup);
    }

    public boolean hasLineageSpecificGain()
    {
        double common_gain = 0.0;
        NodeWithRates[] nodes = getMainTree().getDFT();
        for (int node_idx=0; node_idx<nodes.length && !Double.isNaN(common_gain); ++node_idx)
        {
            NodeWithRates N = nodes[node_idx];
            if (!N.isRoot())
            {
                if (node_idx==0)
                    common_gain = N.getTransferRate();
                else if (N.getTransferRate() != common_gain)
                    common_gain = Double.NaN;
            }
        }
        return Double.isNaN(common_gain);
    }

    public boolean hasLineageSpecificLoss()
    {
        double common_loss = 1.0;
        NodeWithRates[] nodes = getMainTree().getDFT();
        for (int node_idx=0; node_idx<nodes.length && !Double.isNaN(common_loss); ++node_idx)
        {
            NodeWithRates N = nodes[node_idx];
            if (node_idx==0)
                common_loss = N.getLossRate();
            else if (N.getLossRate() != common_loss)
                common_loss = Double.NaN;
        }
        return !(!Double.isNaN(common_loss) && (common_loss != 1.0));
    }


    /**
     * Sets the duplication rate on an edge
     * @param node_idx child node on the edge
     * @param rate new duplication rate
     */
    public void setNodeDuplicationRate(int node_idx, double rate)
    {
        NodeWithRates protoN = main_tree.getNode(node_idx);
        //Verbose.message("RV.sNDR "+node_idx+"/"+protoN.getTaxonName()+"\t"+protoN.getDuplicationRate()+" -> "+rate);
        protoN.setDuplicationRate(rate);
        for (int cidx=0; cidx<rate_tree.length; cidx++)
        {
            TreeWithRates class_tree = rate_tree[cidx];
            NodeWithRates N = class_tree.getNode(node_idx);
            double xd = duplication_rate_categories[getDuplicationRateCategory(cidx)];
            N.setDuplicationRate(xd*rate);
        }
    }
    
    /**
     * Sets the transfer/gain rate on an edge
     * @param node_idx child node on the edge
     * @param rate new transfer/gain rate
     */
    public void setNodeTransferRate(int node_idx, double rate)
    {
        NodeWithRates protoN = main_tree.getNode(node_idx);
        //Verbose.message("RV.sNTR "+node_idx+"/"+protoN.getTaxonName()+"\t"+protoN.getTransferRate()+" -> "+rate);
        protoN.setTransferRate(rate);
        for (int cidx=0; cidx<rate_tree.length; cidx++)
        {
            TreeWithRates class_tree = rate_tree[cidx];
            NodeWithRates N = class_tree.getNode(node_idx);
            double xt = transfer_rate_categories[getTransferRateCategory(cidx)];
            N.setTransferRate(xt*rate);
        }
    }
    
    /**
     * Sets the loss rate on an edge
     * @param node_idx child node on the edge
     * @param rate new loss rate
     */
    public void setNodeLossRate(int node_idx, double rate)
    {
        NodeWithRates protoN = main_tree.getNode(node_idx);
        protoN.setLossRate(rate);
        for (int cidx=0; cidx<rate_tree.length; cidx++)
        {
            TreeWithRates class_tree = rate_tree[cidx];
            NodeWithRates N = class_tree.getNode(node_idx);
            double xl = loss_rate_categories[getLossRateCategory(cidx)];
            N.setLossRate(xl*rate);
        }
    }
    
    /**
     * Sets the length an edge
     * @param node_idx child node on the edge
     * @param len new edge length
     */
    public void setEdgeLength(int node_idx, double len)
    {
        NodeWithRates protoN = main_tree.getNode(node_idx);
        //Verbose.message("RV.sEL "+node_idx+"/"+protoN.getTaxonName()+"\t"+protoN.getLength()+" -> "+len);
        protoN.setLength(len);
        for (int cidx=0; cidx<rate_tree.length; cidx++)
        {
            TreeWithRates class_tree = rate_tree[cidx];
            NodeWithRates N = class_tree.getNode(node_idx);
            double xe = edge_length_categories[getEdgeLengthCategory(cidx)];
            N.setLength(xe*len);
        }
    }
    
    private static void setArrayValues(double[] A, double alpha, int num_gamma_classes)
    {
        if (num_gamma_classes == 1)
        {
            A[0]=1.0;
        } else
        {
            DiscreteGamma G = new DiscreteGamma(alpha);

            double[] val = G.getPartitionMeans(num_gamma_classes);
            for (int i=0; i<val.length; i++)
            {
                A[i] = val[i];
                //Verbose.message("RV.sAV "+i+"/"+A.length+"\t"+A[i]+"\t// alpha="+alpha);
            }
        }
        if (num_gamma_classes != A.length)
        {
            A[A.length-1] = 0.0;
        }
    }
    
    /**
     * Sets the gamma rate parameter for duplication and updates the rate-class trees
     * @param alpha new shape parameter
     */
    public void setDuplicationRateAlpha(double alpha)
    {
        duplication_rate_alpha = alpha;
        setArrayValues(duplication_rate_categories, alpha, getNumDuplicationRateGammaCategories());
        updateClassTrees();
    }
    
    /**
     * Returns the alpha parameter for the Gamma distribution on duplication rates.
     * @return shape parameter
     */
    public double getDuplicationRateAlpha()
    {
        return duplication_rate_alpha;
    }

    /**
     * Sets the number of duplication rate categories.
     *
     * @param num_categories a positive number
     */
    private void setNumDuplicationRateCategories(int num_categories)
    {
        duplication_rate_categories = new double[num_categories+1];
        setArrayValues(duplication_rate_categories,duplication_rate_alpha,num_categories);
    }
    
    public int getNumDuplicationRateGammaCategories()
    {
        return duplication_rate_categories.length-1;
    }
    
    public double getDuplicationForbiddenProportion()
    {
        return duplication_forbidden;
    }
    
    private void initDuplicationRateCategories(int num_categories, double alpha)
    {
        duplication_rate_alpha = alpha;
        setNumDuplicationRateCategories(num_categories);
    }
    
    /**
     * Returns the multipliers for duplication rates in the rate categories
     * (computed as mean values for percentiles of the Gamma distribution)
     * @return array of multipliers
     *
     */
    public double[] getDuplicationRateCategories()
    {
        return duplication_rate_categories;
    }
        
    
    /**
     * Sets the gamma rate parameter for loss and updates the rate-class trees
     * @param alpha new shape parameter
     */
    public void setLossRateAlpha(double alpha)
    {
        loss_rate_alpha = alpha;
        setArrayValues(loss_rate_categories, alpha, getNumLossRateGammaCategories());
        updateClassTrees();
    }

    /**
     * Returns the alpha parameter for the Gamma distribution on transfer rates.
     * @return shape parameter
     */
    public double getLossRateAlpha()
    {
        return loss_rate_alpha;
    }
    
    /**
     * Sets the number of loss rate categories.
     *
     * @param num_categories a positive number
     */
    private void setNumLossRateCategories(int num_categories)
    {
        loss_rate_categories = new double[num_categories+1];
        setArrayValues(loss_rate_categories,loss_rate_alpha,num_categories);
    }
    
    public int getNumLossRateGammaCategories()
    {
        return loss_rate_categories.length-1;
    }
         
    public double getLossForbiddenProportion()
    {
        return loss_forbidden;
    }
    
    private void initLossRateCategories(int num_categories, double alpha)
    {
        loss_rate_alpha = alpha;
        setNumLossRateCategories(num_categories);
    }
    

    /**
     * Returns the multipliers for loss rates in the rate categories
     * (computed as mean values for percentiles of the Gamma distribution)
     *
     * @return array of rate multipliers
     */
    public double[] getLossRateCategories()
    {
        return loss_rate_categories;
    }
    
    
    /**
     * Sets the gamma rate parameter for transfer and updates the rate-class trees
     * @param alpha new shape parameter
     */
    public void setTransferRateAlpha(double alpha)
    {
        transfer_rate_alpha = alpha;
        setArrayValues(transfer_rate_categories, alpha, getNumTransferRateGammaCategories());
        updateClassTrees();
    }
    
    /**
     * Returns the alpha parameter for the Gamma distribution on transfer rates.
     * @return shape parameter
     */
    public double getTransferRateAlpha()
    {
        return transfer_rate_alpha;
    }

    /**
     * Returns the multipliers for transfer rates in the rate categories
     * (computed as mean values for percentiles of the Gamma distribution)
     *
     * @return array of rate multipliers
     */
    public double[] getTransferRateCategories()
    {
        return transfer_rate_categories;
    }
        
    
    /**
     * Sets the number of transfer rate categories.
     *
     * @param num_categories a positive number
     */
    private void setNumTransferRateCategories(int num_categories)
    {
        transfer_rate_categories = new double[num_categories+1];
        setArrayValues(transfer_rate_categories,transfer_rate_alpha,num_categories);
    }
    
    private void initTransferRateCategories(int num_categories, double alpha)
    {
        transfer_rate_alpha = alpha;
        setNumTransferRateCategories(num_categories);
    }
    
    public int getNumTransferRateGammaCategories()
    {
        return transfer_rate_categories.length-1;
    }
    
    public double getTransferForbiddenProportion()
    {
        return transfer_forbidden;
    }
    
    /**
     * Sets the gamma rate parameter for edge length and updates the rate-class trees
     * @param alpha new shape parameter
     */
    public void setEdgeLengthAlpha(double alpha)
    {
        edge_length_alpha=alpha;
        setArrayValues(edge_length_categories, alpha, getNumEdgeLengthGammaCategories());
        updateClassTrees();
    }
    
    /**
     * Returns the alpha parameter for the Gamma distribution on edge lengths.
     * @return shape parameter
     */
    public double getEdgeLengthAlpha()
    {
        return edge_length_alpha;
    }

    /**
     * Sets the number of edge length categories.
     *
     * @param num_categories a positive number
     */
    private void setNumEdgeLengthCategories(int num_categories)
    {
        edge_length_categories = new double[num_categories];
        setArrayValues(edge_length_categories,edge_length_alpha,num_categories);
    }
    
    private void initEdgeLengthCategories(int num_categories, double alpha)
    {
        edge_length_alpha = alpha;
        setNumEdgeLengthCategories(num_categories);
    }
    
    
    public int getNumEdgeLengthGammaCategories()
    {
        return edge_length_categories.length;
    }
    
    /**
     * Returns the multipliers for edge lengths in the rate categories
     * (computed as mean values for percentiles of the Gamma distribution)
     *
     * @return array of length multipliers
     */
    public double[] getEdgeLengthCategories()
    {
        return edge_length_categories;
    }

    /**
     * Returns the appropriate class index for these categories
     *
     * @param duplication_category category for duplication rate
     * @param loss_category category for loss rate
     * @param transfer_category category for gain rate
     * @param edge_length_category category for edge length
     * @return combined rate class index
     */
    public final int getClassIndex(int duplication_category, int loss_category, int transfer_category, int edge_length_category)
    {
        int nd = duplication_rate_categories.length;
        int nl = loss_rate_categories.length;
        int nt = transfer_rate_categories.length;
        //int ne = edge_length_categories.length;
        
        int class_idx = duplication_category + nd*(loss_category + nl*(transfer_category + nt*edge_length_category));

        return class_idx;
    }
    
    /**
     * Returns the duplication rate category for this class
     *
     * @param class_idx combined rate class index
     * @return index for duplication rate category
     */
    public final int getDuplicationRateCategory(int class_idx)
    {
        int nd = duplication_rate_categories.length;

        return class_idx % nd;
    }

    /**
     * Prior probability for a rate class
     *
     * @param class_idx combined class index
     * @return prior probability for the class
     */
    public double getClassProbability(int class_idx)
    {
        int dup_idx = getDuplicationRateCategory(class_idx);
        int dup_n= getNumDuplicationRateGammaCategories();
        double dup_p = (dup_idx == dup_n)
                       ?duplication_forbidden
                       :(1.-duplication_forbidden)/dup_n;
        int loss_idx = getLossRateCategory(class_idx);
        int loss_n = getNumLossRateGammaCategories();
        double loss_p = (loss_idx == loss_n)
                       ?loss_forbidden
                       :(1.-loss_forbidden)/loss_n;
        int transfer_idx = getTransferRateCategory(class_idx);
        int transfer_n = getNumTransferRateGammaCategories();
        double transfer_p = (transfer_idx == transfer_n)
                            ?transfer_forbidden
                            :(1.-transfer_forbidden)/transfer_n;
        int edge_n = edge_length_categories.length;
        double edge_p = 1.0/edge_n;
        return dup_p * loss_p * transfer_p * edge_p;
    }
    
    /**
     * Returns the loss rate category for this class
     *
     * @param class_idx combined rate class index
     * @return loss rate category in this class
     */
    public final int getLossRateCategory(int class_idx)
    {
        int nd = duplication_rate_categories.length;
        int nl = loss_rate_categories.length;
        
        return (class_idx / nd ) % nl;
    }
    
    /**
     * Returns the transfer rate category for this class
     *
     * @param class_idx combined rate class index
     * @return gain rate category for this class
     */
    public final int getTransferRateCategory(int class_idx)
    {
        int nd = duplication_rate_categories.length;
        int nl = loss_rate_categories.length;
        int nt = transfer_rate_categories.length;
        
        return (class_idx / (nd *nl) ) % nt;
    }
    
    /**
     * Returns the edge length category for this class
     *
     * @param class_idx combined class index
     * @return edge length category for this class
     */
    public final int getEdgeLengthCategory(int class_idx)
    {
        int nd = duplication_rate_categories.length;
        int nl = loss_rate_categories.length;
        int nt = transfer_rate_categories.length;
        
        return (class_idx / (nd*nl*nt));
    }
    
    /**
     * Returns the number of discrete classes
     *
     * @return number of classes
     */
    public final int getNumClasses()
    {
        int nd = duplication_rate_categories.length;
        int nl = loss_rate_categories.length;
        int nt = transfer_rate_categories.length;
        int ne = edge_length_categories.length;
        return nd*nl*nt*ne;
    }
    
    /**
     * Returns a class-specific tree, in which 
     * branch parameters are scaled appropriately
     *
     * @param class_idx combined class index
     * @return tree with class-specific rates and lengths for the edges
     */
    public TreeWithRates getRateTree(int class_idx)
    {
        TreeWithRates class_tree = rate_tree[class_idx];
        return class_tree;
    }
    
    public TreeWithRates getMainTree()
    {
        return main_tree;
    }
    
    public DiscreteDistribution getRootPrior()
    {
        return root_prior_distribution;
    }
    
    public void setRootPrior(DiscreteDistribution root_prior)
    {
        this.root_prior_distribution = root_prior;
    }


    public String getBriefModelDescription()
    {
        StringBuffer sb = new StringBuffer();
        if (hasDuplication())
        {
            if (hasGain())
            {
                sb.append("GLD");
            } else
                sb.append("DL");
        } else
        {
            if (hasGain())
            {
                sb.append("GL");
            } else
                sb.append("PL");
        }
        sb.append(';');

        if (hasLineageSpecificLoss())
        {
            sb.append("heterogeneous");
        } else
            sb.append("homogeneous");
        if (!hasLineageSpecificDuplication())
            sb.append(",uniform dup");
        if (!hasLineageSpecificGain())
            sb.append(",uniform gain");

        sb.append(";root:");
        DiscreteDistribution distr = getRootPrior();
        if (distr instanceof NegativeBinomial)
        {
            sb.append("NegBin");
        } else if (distr instanceof ShiftedGeometric)
            sb.append("ShGeom");
        else if (distr instanceof PointDistribution)
            sb.append("Point");
        else sb.append(distr.getClass().getSimpleName());

        StringBuffer cat_sb = new StringBuffer();
        int nlen = getNumEdgeLengthGammaCategories();
        if (nlen>1)
        {
            cat_sb.append("length:");
            cat_sb.append(nlen);
            cat_sb.append("G");
        }
        int ngain = getNumTransferRateGammaCategories();
        double pnog = getTransferForbiddenProportion();
        if (ngain>1 || pnog>0.0)
        {
            if (cat_sb.length()!=0)
                cat_sb.append(',');
            cat_sb.append("gain:");
            if (ngain>1)
            {
                cat_sb.append(ngain);
                cat_sb.append("G");
                if (pnog>0.0)
                    cat_sb.append("+Z");
            } else
                cat_sb.append("Z");
        }
        int ndup = getNumDuplicationRateGammaCategories();
        double pnod = getDuplicationForbiddenProportion();
        if (ndup>1 || pnod>0.0)
        {
            if (cat_sb.length()!=0)
                cat_sb.append(',');
            cat_sb.append("dup:");
            if (ndup>1)
            {
                cat_sb.append(ndup);
                cat_sb.append("G");
                if (pnod>0.)
                    cat_sb.append("+Z");
            } else
                cat_sb.append("Z");
        }
        int nloss = getNumLossRateGammaCategories();
        if (nloss>1)
        {
            if (cat_sb.length()!=0)
                cat_sb.append(',');
            cat_sb.append("loss:");
            cat_sb.append(nloss);
            cat_sb.append("G");
        }
        if (cat_sb.length()==0)
            cat_sb.append("no site variation");
        else
            cat_sb.append(" site categories");
        sb.append(';');
        sb.append(cat_sb);
        return sb.toString();
    }

    /**
     * The string used in the rate file to mark the rate variation parameters
     */
    static String RATE_VARIATION_PREFIX = "|variation";
    static String ROOT_PRIOR_PREFIX = "|root";
    static String MODEL_END = "|End";
    
    public String tableRates()
    {
        StringBuffer sb = new StringBuffer(main_tree.getRoot().tableRates());
        sb.append(RATE_VARIATION_PREFIX);
        sb.append("\tduplication\t");
        sb.append(getNumDuplicationRateGammaCategories());
        sb.append('\t');
        sb.append(duplication_rate_alpha);
        sb.append('\t');
        sb.append(duplication_forbidden);
        sb.append("\n");
        sb.append(RATE_VARIATION_PREFIX);
        sb.append("\tloss\t");
        sb.append(getNumLossRateGammaCategories());
        sb.append('\t');
        sb.append(loss_rate_alpha);
        sb.append('\t');
        sb.append(loss_forbidden);
        sb.append("\n");
        sb.append(RATE_VARIATION_PREFIX);
        sb.append("\ttransfer\t");
        sb.append(getNumTransferRateGammaCategories());
        sb.append('\t');
        sb.append(transfer_rate_alpha);
        sb.append('\t');
        sb.append(transfer_forbidden);
        sb.append("\n");
        sb.append(RATE_VARIATION_PREFIX);
        sb.append("\tlength\t");
        sb.append(getNumEdgeLengthGammaCategories());
        sb.append('\t');
        sb.append(edge_length_alpha);
        sb.append("\n");
        sb.append(ROOT_PRIOR_PREFIX);
        sb.append('\t');
        sb.append(root_prior_distribution.getClass().getSimpleName());
        double[] params = root_prior_distribution.getParameters();
        for (int i=0; i<params.length; i++)
        {
            sb.append('\t');
            sb.append(params[i]);
        }
        sb.append("\n");
        sb.append(MODEL_END);
        sb.append("\n");
            
        return sb.toString();
    }
    
    private void initFromFile(Reader reader) throws FileFormatException, IOException
    {
        BufferedReader BR = new BufferedReader(reader);
        NodeWithRates[] tree_node = main_tree.getRoot().getDFT(); // important to use DFT for the root and not TreeWithRates ...
        
        int node_idx= 0;
        {
            String line=null;
            do
            {
                line=BR.readLine();
                if (line != null)
                {
                    if (line.startsWith(MODEL_END))
                        break;

                    int pos=0;
                    while (pos<line.length() && Character.isWhitespace(line.charAt(pos))) pos++;
                    if (pos==line.length() || line.charAt(pos)=='#')
                        continue;
                    if (line.startsWith(RATE_VARIATION_PREFIX))
                    {
                        String variation_data = line.substring(RATE_VARIATION_PREFIX.length()+1);
                        String[] fields = variation_data.split("\\s+");
                        if (fields.length<3)
                            throw new FileFormatException("Rate variation line has bad syntax: "+line);
                            
                        int num_categories = Integer.parseInt(fields[1]);
                        double alpha = Double.parseDouble(fields[2]);
                        double zero = 0.0;
                        if (fields.length>3)
                        {
                            zero = Double.parseDouble(fields[3]);
                        }
                        if ("loss".equals(fields[0]))
                        {
                            initLossRateCategories(num_categories, alpha);
                            setLossForbidden(zero);
                        } else if ("duplication".equals(fields[0]))
                        {
                            initDuplicationRateCategories(num_categories, alpha);
                            setDuplicationForbidden(zero);
                        } else if ("transfer".equals(fields[0]))
                        {
                            initTransferRateCategories(num_categories, alpha);
                            setTransferForbidden(zero);
                        } else if ("length".equals(fields[0]))
                        {
                            initEdgeLengthCategories(num_categories, alpha);
                        } else 
                        {
                            throw new FileFormatException("Variation type '"+fields[0]+"' is not recognized in the line '"+line+"'");
                        }
                    } else if (line.startsWith(ROOT_PRIOR_PREFIX))
                    {
                        String prior_data = line.substring(ROOT_PRIOR_PREFIX.length()+1);
                        String[] fields = prior_data.split("\\s+");
                        double[] params = new double[fields.length-1];
                        for (int i=0; i<params.length; i++)
                            params[i] = Double.parseDouble(fields[i+1]);
                        if (fields[0].equals(Poisson.class.getSimpleName()))
                            root_prior_distribution = new Poisson(params[0]);
                        else if (fields[0].equals(NegativeBinomial.class.getSimpleName()))
                            root_prior_distribution = new NegativeBinomial(params[0], params[1]);
                        else if (fields[0].equals(PointDistribution.class.getSimpleName()))
                            root_prior_distribution = new PointDistribution(params[0]);
                        else if (fields[0].equals(ShiftedGeometric.class.getSimpleName()))
                            root_prior_distribution = new ShiftedGeometric(params[0], params[1]);
                        else
                            throw new FileFormatException("Root prior distribution '"+fields[0]+"' is unknown in line '"+line+"'");
                    } else
                    {
                        //System.out.println("#*RV.iFF node "+node_idx+" pos "+pos+" line `"+line+"'"+"\t"+line.charAt(pos)+"\t"+(line.charAt(pos)=='#'));
                        tree_node[node_idx].setEdgeParametersFromTable(line);
                        node_idx++;
                    }
                }     
            } while (line != null);
        }
        BR.close();
        if (root_prior_distribution == null)
            throw new FileFormatException("Root prior line (starts with '"+ROOT_PRIOR_PREFIX+"' is missing");
        if (duplication_rate_categories == null)
            throw new FileFormatException("Duplication rate variation info is missing");
        if (transfer_rate_categories == null)
            throw new FileFormatException("Transfer rate variation info is missing");
        if (loss_rate_categories == null)
            throw new FileFormatException("Loss rate variation info is missing");
        if (edge_length_categories == null)
            throw new FileFormatException("Edge length variation info is missing");

        initDataStructures();
    }

    /**
     * Reads in a RateVariation from saved rates file.
     *
     * @param R reader for saved rates (will not be closed)
     * @param main_tree underlying tree
     * @return rate model
     *
     * @throws FileFormatException if input is not in a parsable format
     * @throws IOException if something goes wrong at reading
     */
    public static RateVariation read(Reader R, TreeWithRates main_tree) throws FileFormatException, IOException
    {
        RateVariation M = new RateVariation();
        M.main_tree = main_tree;
        M.initFromFile(R);
        return M;
    }

    /**
     * Our own exception type.
     */
    public static class FileFormatException extends IOException
    {
        private FileFormatException(String msg)
        {
            super(msg);
        }
    }

    /**
     * Used with {@link #getScaledTree(EdgeLengthType) }
     */
    private enum EdgeLengthType {DUPLICATION, LOSS, TRANSFER, MAX};

    /**
     * Calculates edge lengths according to the requested scaling. When scaling by rates,
     * edge lengths are set by the total rate (length times lineage-specific rate).
     *
     * @param type type of scaling:
     * @return scaled tree (lengths and rates are rescaled)
     */
    private TreeWithRates getScaledTree(EdgeLengthType type)
    {
        NodeWithRates main_root = main_tree.getRoot();
        TreeWithRates scaled_tree = new TreeWithRates(NodeWithRates.copyTree(main_root));
        
        NodeWithRates nodes[] = scaled_tree.getDFT();
        for (int node_idx=0; node_idx<nodes.length; node_idx++)
        {
            NodeWithRates mainN = main_tree.getNode(node_idx);
            if (!mainN.isRoot())
            {
                double len = mainN.getLength();
                double xd = len*mainN.getDuplicationRate();
                double xt = len*mainN.getTransferRate();
                double xl = len*mainN.getLossRate();
                if (type == EdgeLengthType.DUPLICATION)
                    len = xd;
                else if (type == EdgeLengthType.TRANSFER)
                    len = xt;
                else if (type == EdgeLengthType.LOSS)
                    len = xl;
                else // (type == EdgeLengthType.MAX)
                {
                    len = Math.max(xd,xl);
                    len = Math.max(len,xt);
                }
                NodeWithRates scaledN = scaled_tree.getNode(node_idx);
                scaledN.setLength(len);
                if (len != 0.0)
                {
                    scaledN.setDuplicationRate(xd/len);
                    scaledN.setTransferRate(xt/len);
                    scaledN.setLossRate(xl/len);
                }
                // inner nodes in the scaled tree have no name
                // so that Phylip can plot the trees easily
                if (!scaledN.isLeaf())
                    scaledN.setName(null);
            }
        } // for all nodes
            
        return scaled_tree;
    }

    /**
     * Computes a scaled version of the main tree that defines the same 
     * joint size distribution at the nodes. 
     * Loss rate is 1.0 at every node in the scaled tree.
     * 
     * @return scaled tree
     */
    public TreeWithRates getLossTree()
    {
        return getScaledTree(EdgeLengthType.LOSS);
    }
    
    /**
     * Computes a scaled version of the main tree that defines the same 
     * joint size distribution at the nodes. 
     * Duplication rate is 1.0 at every node in the scaled tree.
     * 
     * @return scaled tree
     */
    public TreeWithRates getDuplicationTree()
    {
        return getScaledTree(EdgeLengthType.DUPLICATION);
    }
    
    /**
     * Computes a scaled version of the main tree that defines the same 
     * joint size distribution at the nodes. 
     * Ttransfer rate is 1.0 at every node in the scaled tree.
     * 
     * @return scaled tree
     */
    public TreeWithRates getTransferTree()
    {
        return getScaledTree(EdgeLengthType.TRANSFER);
    }
    
  
    /**
     * Computes a scaled version of the main tree that defines the same 
     * joint size distribution at the nodes. 
     * The largest among loss, duplication and transfer rates is 1.0 at every node in the scaled tree.
     * 
     * @return scaled tree
     */
    public TreeWithRates getMaxRateTree()
    {
        return getScaledTree(EdgeLengthType.MAX);
    }
    
    private void go(String[] args) throws Exception
    {
        reportLaunch(args);
        if (args.length < 1)
        {
            System.err.println("Call as java "+RateVariation.class.getCanonicalName()+" [switches] tree [rate file]");
            System.err.println("\tSwitches");
            System.err.println("\t-gain_factor d\tMultiplier applied to all lineage-specific gain (transfer) rates");
            System.err.println("\t-loss_factor d\tMultiplier applied to all lineage-specific loss rates");
            System.err.println("\t-duplication_factor d\tMultiplier applied to all lineage-specific duplication rates");
            System.err.println("\t-length_factor d\tMultiplier applied to all branch lengths");
            System.err.println("\t-gain_a d\tGamma distribution shape parameter for gain rate variation");
            System.err.println("\t-loss_a d\tGamma distribution shape parameter for loss rate variation");
            System.err.println("\t-duplication_a d\tGamma distribution shape parameter for duplication rate variation");
            System.err.println("\t-length_a d\tGamma distribution shape parameter for edge length variation");
            System.err.println("\t-gain_k\tNumber of discrete Gamma categories for gain rate variation (=1 for no rate variation)");
            System.err.println("\t-loss_k\tNumber of discrete Gamma categories for loss rate variation (=1 for no rate variation)");
            System.err.println("\t-duplication_k\tNumber of discrete Gamma categories for duplication rate variation (=1 for no rate variation)");
            System.err.println("\t-length_k\tNumber of discrete Gamma categories for edge length variation (=1 for no rate variation)");
            System.err.println("\t-root D,p0...\tRoot prior: D may be 'Poisson', 'NegativeBinomial', 'PointDistribution', 'ShiftedGeometric',\n\t\t\tfollowed by the parameters in comma-separated list p0,p1,...");
            System.err.println("\t-forbidden_duplication\tFraction of families with no duplications (=0.0 for disabling this category)");
            System.err.println("\t-forbidden_gain\tFraction of families with no gains (=0.0 for disabling this category)");

            System.exit(2010);
        }



        if (GAIN_FACTOR != 1.0)
        {
            reportOtherArguments("Stretch for gain rate: "+GAIN_FACTOR);
            NodeWithRates[] nodes = main_tree.getDFT();
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                NodeWithRates N = nodes[node_idx];
                if (!N.isRoot())
                {
                    setNodeTransferRate(node_idx,N.getTransferRate()*GAIN_FACTOR);
                }
            }
        }
        if (LOSS_FACTOR != 1.0)
        {
            reportOtherArguments("Stretch for loss rate: "+LOSS_FACTOR);
            NodeWithRates[] nodes = main_tree.getDFT();
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                NodeWithRates N = nodes[node_idx];
                if (!N.isRoot())
                {
                    setNodeLossRate(node_idx,N.getLossRate()*LOSS_FACTOR);
                }
            }
        }
        if (DUPLICATION_FACTOR != 1.0)
        {
            reportOtherArguments("Stretch for duplication rate: "+DUPLICATION_FACTOR);
            NodeWithRates[] nodes = main_tree.getDFT();
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                NodeWithRates N = nodes[node_idx];
                if (!N.isRoot())
                {
                    setNodeDuplicationRate(node_idx,N.getDuplicationRate()*DUPLICATION_FACTOR);
                }
            }
        }
        if (EDGE_LENGTH_FACTOR != 1.0)
        {
            reportOtherArguments("Stretch for edge length: "+EDGE_LENGTH_FACTOR);
            NodeWithRates[] nodes = main_tree.getDFT();
            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                NodeWithRates N = nodes[node_idx];
                if (!N.isRoot())
                {
                    setEdgeLength(node_idx,N.getLength()*EDGE_LENGTH_FACTOR);
                }
            }
        }
        if (!Double.isNaN(NEW_GAIN_ALPHA ))
        {
            reportOtherArguments("New gain rate variation alpha: "+NEW_GAIN_ALPHA);
            setTransferRateAlpha(NEW_GAIN_ALPHA);
        }
        if (!Double.isNaN(NEW_LOSS_ALPHA ))
        {
            reportOtherArguments("New loss rate variation alpha: "+NEW_LOSS_ALPHA);
            setLossRateAlpha(NEW_LOSS_ALPHA);
        }
        if (!Double.isNaN(NEW_DUPLICATION_ALPHA ))
        {
            reportOtherArguments("New duplication rate variation alpha: "+NEW_DUPLICATION_ALPHA);
            setDuplicationRateAlpha(NEW_DUPLICATION_ALPHA);
        }
        if (!Double.isNaN(NEW_EDGE_LENGTH_ALPHA ))
        {
            reportOtherArguments("New edge length variation alpha: "+NEW_EDGE_LENGTH_ALPHA);
            setEdgeLengthAlpha(NEW_EDGE_LENGTH_ALPHA);
        }
        if (NEW_GAIN_CATEGORIES != 0)
        {
            reportOtherArguments("New gain rate categories: "+NEW_GAIN_CATEGORIES);
            setNumTransferRateCategories(NEW_GAIN_CATEGORIES);
        }
        if (NEW_LOSS_CATEGORIES != 0)
        {
            reportOtherArguments("New loss rate categories: "+NEW_LOSS_CATEGORIES);
            setNumLossRateCategories(NEW_LOSS_CATEGORIES);
        }
        if (NEW_DUPLICATION_CATEGORIES != 0)
        {
            reportOtherArguments("New duplication rate categories: "+NEW_DUPLICATION_CATEGORIES);
            setNumDuplicationRateCategories(NEW_DUPLICATION_CATEGORIES);
        }
        if (NEW_EDGE_LENGTH_CATEGORIES != 0)
        {
            reportOtherArguments("New edge length categories: "+NEW_EDGE_LENGTH_CATEGORIES);
            setNumEdgeLengthCategories(NEW_EDGE_LENGTH_CATEGORIES);
        }
        if (NEW_FORBIDDEN_GAIN != null)
        {
            reportOtherArguments("Forbidden gain: "+NEW_FORBIDDEN_GAIN);
            this.setTransferForbidden(NEW_FORBIDDEN_GAIN);
        }
        if (NEW_FORBIDDEN_DUPLICATION != null)
        {
            reportOtherArguments("Forbidden duplication: "+NEW_FORBIDDEN_DUPLICATION);
            this.setDuplicationForbidden(NEW_FORBIDDEN_DUPLICATION);
        }


        if (NEW_ROOT_PRIOR != null)
        {
            StringBuffer distr = new StringBuffer(NEW_ROOT_PRIOR.getClass().getSimpleName());
            for (int i=0; i<NEW_ROOT_PRIOR.getNumParameters(); i++)
            {
                distr.append(',');
                distr.append(Double.toString(NEW_ROOT_PRIOR.getParameters()[i]));
            }
            reportOtherArguments("New root prior distribution: "+distr.toString());
            setRootPrior(NEW_ROOT_PRIOR);
        }
        
        if (getNumEdgeLengthGammaCategories()>1)
        {
            System.out.print("# Edge length multipliers: ");
            for (int i=0; i<edge_length_categories.length; i++)
                System.out.print("\t"+edge_length_categories[i]);
            System.out.println();
        }
        if (getNumDuplicationRateGammaCategories()>1)
        {
            System.out.print("# Duplication rate multipliers: ");
            for (int i=0; i<duplication_rate_categories.length; i++)
                System.out.print("\t"+duplication_rate_categories[i]);
            System.out.println();
        }
        if (getNumLossRateGammaCategories()>1)
        {
            System.out.print("# Loss rate multipliers: ");
            for (int i=0; i<loss_rate_categories.length; i++)
                System.out.print("\t"+loss_rate_categories[i]);
            System.out.println();
        }
        
        if (getNumTransferRateGammaCategories()>1)
        {
            System.out.print("# Transfer rate multipliers: ");
            for (int i=0; i<transfer_rate_categories.length; i++)
                System.out.print("\t"+transfer_rate_categories[i]);
            System.out.println();
        }
        
        System.out.println("# Loss tree: "+getLossTree().getRoot().newickTree());
        System.out.println("# Duplication tree: "+getDuplicationTree().getRoot().newickTree());
        System.out.println("# Gain tree: "+getTransferTree().getRoot().newickTree());
        System.out.println("# Max-rate tree: "+getMaxRateTree().getRoot().newickTree());
        
        System.out.println(tableRates());        
    }

    private static double EDGE_LENGTH_FACTOR = 1.0;
    private static double DUPLICATION_FACTOR = 1.0;
    private static double GAIN_FACTOR = 1.0;
    private static double LOSS_FACTOR = 1.0;
    private static double NEW_EDGE_LENGTH_ALPHA = Double.NaN;
    private static double NEW_DUPLICATION_ALPHA = Double.NaN;
    private static double NEW_GAIN_ALPHA = Double.NaN;
    private static double NEW_LOSS_ALPHA = Double.NaN;
    private static int NEW_EDGE_LENGTH_CATEGORIES = 0;
    private static int NEW_GAIN_CATEGORIES = 0;
    private static int NEW_DUPLICATION_CATEGORIES = 0;
    private static int NEW_LOSS_CATEGORIES = 0;
    private static DiscreteDistribution NEW_ROOT_PRIOR=null;
    private static Double NEW_FORBIDDEN_DUPLICATION = null;
    private static Double NEW_FORBIDDEN_GAIN = null;

    public static void main(String[] args) 
    {
        Verbose.setVerbose(false);

        int num_switches = 0;
        try {
            while (args.length>2*num_switches && args[2*num_switches].startsWith("-"))
            {
                String arg_switch = args[2*num_switches].substring(1);
                if (arg_switch.equals("h"))
                {
                    (new RateVariation()).go(new String[0]); // will throw an Exception
                }
                if (args.length==2*num_switches+1)
                    throw new IllegalArgumentException("Missing argument for switch "+args[2*num_switches]);
                String arg_value = args[2*num_switches+1];
                if (arg_switch.equals("v"))
                {
                    Verbose.setVerbose(arg_value.equals("true"));
                } else if (arg_switch.equals("duplication_factor"))
                {
                    DUPLICATION_FACTOR = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("length_factor"))
                {
                    EDGE_LENGTH_FACTOR = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("gain_factor"))
                {
                    GAIN_FACTOR = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("loss_factor"))
                {
                    LOSS_FACTOR = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("duplication_a"))
                {
                    NEW_DUPLICATION_ALPHA = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("gain_a"))
                {
                    NEW_GAIN_ALPHA = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("loss_a"))
                {
                    NEW_LOSS_ALPHA = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("length_a"))
                {
                    NEW_EDGE_LENGTH_ALPHA = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("gain_k"))
                {
                    NEW_GAIN_CATEGORIES = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("loss_k"))
                {
                    NEW_LOSS_CATEGORIES = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("duplication_k"))
                {
                    NEW_DUPLICATION_CATEGORIES = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("length_k"))
                {
                    NEW_EDGE_LENGTH_CATEGORIES = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("root"))
                {
                    String[] fields = arg_value.split(",");
                    double[] params = new double[fields.length-1];
                    for (int i=0; i<params.length; i++)
                        params[i] = Double.parseDouble(fields[i+1]);
                    if (fields[0].equals(Poisson.class.getSimpleName()))
                        NEW_ROOT_PRIOR = new Poisson(params[0]);
                    else if (fields[0].equals(NegativeBinomial.class.getSimpleName()))
                        NEW_ROOT_PRIOR = new NegativeBinomial(params[0], params[1]);
                    else if (fields[0].equals(PointDistribution.class.getSimpleName()))
                        NEW_ROOT_PRIOR = new PointDistribution(params[0]);
                    else if (fields[0].equals(ShiftedGeometric.class.getSimpleName()))
                        NEW_ROOT_PRIOR = new ShiftedGeometric(params[0], params[1]);
                    else
                        throw new FileFormatException("Root prior distribution '"+fields[0]+"' is unknown.");
                } else if (arg_switch.equals("forbidden_duplication"))
                {
                    NEW_FORBIDDEN_DUPLICATION = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("forbidden_gain"))
                {
                    NEW_FORBIDDEN_GAIN = Double.parseDouble(arg_value);
                }
                else
                    throw new IllegalArgumentException("Switch not recognized: '"+args[2*num_switches]+"'");

                num_switches++;
            }

            String[] rest=new String[args.length-2*num_switches];
            if (rest.length == 0)
                (new RateVariation()).go(new String[0]); // will throw an Exception


            for (int j=0; j<rest.length; j++)
                rest[j]=args[2*num_switches+j];
            String tree_file = rest[0];

            TreeWithRates main_tree = new TreeWithRates(ca.umontreal.iro.evolution.Parser.readNewick(new java.io.FileReader(tree_file)));
            String rate_file = null;
            if (rest.length>1)
                rate_file = rest[1];
            RateVariation RV = null;

            if (rate_file == null)
            {
                RV = new RateVariation(main_tree, new Poisson(1.0), 1, 3, 1, 1);
            } else
            {
                RV = RateVariation.read(new java.io.FileReader(rate_file), main_tree);
            }

            RV.go(rest);
        } catch (Exception E)
        {
            die(E);
        }

        

    }
    
}
