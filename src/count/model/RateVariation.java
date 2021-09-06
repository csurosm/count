/*
 * Copyright 2021 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package count.model;

import count.ds.IndexedTree;
import count.matek.DiscreteGamma;


/**
 * Rate variation model: gamma + invariant rate factors
 * for edge length, loss, duplication and transfer.
 *
 * <h2>Rate variation</h2>
 * Rate variation across families is by using a multiplier with
 * a discretized Gamma distribution ({@link DiscreteGamma}),
 * and possibly a separate category for a 0 multiplier.
 * Rate variation is implemented for gain, duplication, loss, and edge length.
 * (Loss rate variation often causes numerical errors though.)
 * The model's rate variation classes
 * correspond to a Cartesian product for
 * the individual rate categories.
 * Each rate variation has a Gamma shape category (methods <code>getXXXAlpha</code>),
 * and some discrete categories, each with their own multipliers, computed as mean values in percentiles
 * (methods <code>getXXXRateCategoryMultipliers()</code>).
 *
 * As an experimental feature, there are
 * possible forbidden duplication and gain categories.
 * The corresponding rate multipliers are 0 in those categories.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class RateVariation implements MixedRateModel
{
    /**
     * Instantiates a rate-variation model. In such a model, 
     * a cross-product of discrete Gamma categories is used for
     * rates (duplication,loss,transfer) and branch lengths.
     * For parameters on which variation is not imposed, the number of categories should be 1.
     * @param base_rates tree with lineage-specific parameters
     * @param root_prior prior size distribution at root
     * @param duplication_rate_categories number of duplication rate categories for rate variation across families
     * @param loss_rate_categories number of loss rate categories for rate variation across families
     * @param gain_rate_categories number of gain rate categories for rate variation across families
     * @param edge_length_categories number of edge length categories for rate variation across families
     */
    public RateVariation(RateModel.GLD base_rates,
            int duplication_rate_categories, 
            int loss_rate_categories, 
            int gain_rate_categories, 
            int edge_length_categories) 
    {
        this(base_rates);
        setClasses(
            duplication_rate_categories, 
            loss_rate_categories, 
            gain_rate_categories, 
            edge_length_categories);
    }
    
    /**
     * Instantiation for setting up without specifying the classes. 
     * @param base_rates tree with lineage-specific parameters
     */
    public RateVariation(RateModel.GLD base_rates)
    {
        this.base_rates = base_rates;
        this.main_tree = base_rates.getTree();
        // we don't know the classes yet
    } 
    
    private RateModel.GLD base_rates;
    private IndexedTree main_tree;
    private double[] gain_multipliers;
    private double[] loss_multipliers;
    private double[] duplication_multipliers;
    private double[] length_multipliers;
    
    private RateModel.GLD[] class_rates; // parametric type is lost  in array, but these are phylogenies
    
    private double gain_alpha=1.0;
    private double loss_alpha=1.0;
    private double duplication_alpha=1.0;
    private double length_alpha=1.0;
    
    private double gain_forbidden = 0.0;
    private double loss_forbidden = 0.0;
    private double duplication_forbidden = 0.0;
    
    
    public RateModel.GLD getBaseModel()
    {
    	return base_rates;
    }
    
    /**
     * Sets the number of categories for duplication rate, loss rate, transfer rate and edge length (in this oder of the 
     * arguments). For parameters on which variation is not imposed, the number of categories should be 1. 
     *
     * @param gain_rate_categories a positive number, use 1 for no rate variation
     * @param loss_rate_categories a positive number, use 1 for no rate variation
     * @param duplication_rate_categories a positive number, use 1 for no rate variation
     * @param edge_length_categories a positive number, use 1 for no rate variation
     */
    public void setClasses(int gain_rate_categories, 
	    int loss_rate_categories, 
	    int duplication_rate_categories, 
	    int edge_length_categories)

    {
    	gain_multipliers 
		= setGammaMultipliers(new double[gain_rate_categories+1],gain_alpha,gain_rate_categories);
    	loss_multipliers 
		= setGammaMultipliers(new double[loss_rate_categories+1],loss_alpha,loss_rate_categories);
    	duplication_multipliers 
		= setGammaMultipliers(new double[duplication_rate_categories+1],duplication_alpha,duplication_rate_categories);
    	length_multipliers 
    	= setGammaMultipliers(new double[edge_length_categories], length_alpha, edge_length_categories);
    	
    	class_rates = new RateModel.GLD[getNumClasses()]; 
    	setClassRates(); // adjusts the rates on the class trees 
    }

    /**
     * Sets the proportions of sites where the invariant model (0 factor) applies instead of Gamma variation
     * 
     * @param gain_forbidden proportion of families with no transfer (between 0.0 and 1.0)
     * @param loss_forbidden proportion of families with no loss (between 0.0 and 1.0)
     * @param duplication_forbidden proportion of families with no duplication (between 0.0 and 1.0)
     */
    public void setForbiddenFractions(double gain_forbidden, double loss_forbidden, double duplication_forbidden)
    {
    	boolean reset_class_rates = 
    			(gain_forbidden==0.0) != (this.gain_forbidden==0.0);
    	this.gain_forbidden = gain_forbidden;
    	reset_class_rates = reset_class_rates 
    			|| (loss_forbidden == 0.0) != (this.loss_forbidden==0.0);
    	this.loss_forbidden = loss_forbidden;
    	reset_class_rates = reset_class_rates 
    			|| (duplication_forbidden == 0.0) != (this.duplication_forbidden==0.0);
    	this.duplication_forbidden = duplication_forbidden;
    	if (reset_class_rates) 
    		setClassRates(); // only if forbidden / not-forbidden changes 
    }
    
    public double getGainForbidden() { return gain_forbidden;}
    public double getLossForbidden() { return loss_forbidden;}
    public double getDuplicationForbidden() {return duplication_forbidden;}
    
    
    
    
    private void setClassRates()
    {
    	int ng = gain_multipliers.length;
    	int nl = loss_multipliers.length;
    	int nd = duplication_multipliers.length;
    	int ne = length_multipliers.length;
    	
    	int num_nodes = main_tree.getNumNodes();
    	
    	int class_idx = 0;
    	// class_idx =  dup_idx + nd*(loss_idx + nl*(gain_idx + ng*length_idx)
		for (int length_idx=0; length_idx<ne; length_idx++)
		{
			double xe = length_multipliers[length_idx];
			int gain_idx = 0;
			while (gain_idx<ng)
			{
				boolean is_excluded = (gain_idx==ng-1 && gain_forbidden==0.0);
				double xg = gain_multipliers[gain_idx];
				++gain_idx;
	    		int loss_idx = 0;
	    		while (loss_idx<nl)
	    		{
	    			is_excluded = is_excluded || (loss_idx==nl-1 && loss_forbidden==0.0);
	    			double xl = loss_multipliers[loss_idx];
	    			++loss_idx;
			    	int dup_idx = 0;
			    	while (dup_idx<nd)
		    		{
			    		is_excluded  = is_excluded || (dup_idx==nd-1 && duplication_forbidden == 0.0);
			            double xd = duplication_multipliers[dup_idx];
			            ++dup_idx;
    					TreeWithRates R;
    					if (is_excluded)
    					{
    						R=null;
    					} else
    					{
    						R = new TreeWithRates(base_rates);
    						if (xl!= 1.0)
    						{
	    						for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // root keeps the base rates 
	    						{
	    							double yl = xl*R.getLossRate(node_idx);
	    							R.setLossRate(node_idx, yl);
	    						}
    						}
    						if (xd != 1.0)
    						{
	    						for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // root keeps the base rates 
	    						{
	    							double yd = xd*R.getDuplicationRate(node_idx);
	    							R.setDuplicationRate(node_idx, yd);
	    						}
    						}
    						for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // root keeps the base rates 
    						{
    							double yg = xg*R.getGainRate(node_idx);
    							
    							// legacy logic for multi-rate variations 
    							double yd = R.getDuplicationRate(node_idx);
    							if (yd==0.0)
    							{
    								double yl = R.getLossRate(node_idx);
    								if (yl==0.0)	
    									R.setGainRate(node_idx, yg);
    								else 
    									R.setGainRate(node_idx, yg*base_rates.getLossRate(node_idx)/yl);
    							} else
    							{
    								R.setGainRate(node_idx,  yg*base_rates.getDuplicationRate(node_idx)/yd);
    							}
    						}
    						if (xe != 1.0)
    						{
	    						for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // root keeps the base rates 
	    						{
	    							double ye = xe*R.getEdgeLength(node_idx);
	    							R.setEdgeLength(node_idx, ye);
	    						}
    						}
    					}
//    					System.out.println("#*class "+class_idx+"\t"+R);
    					class_rates[class_idx++] = R;
    				} // while dup
    			} // while loss
    		} // while gain
    	} //  for length
    }
    
    public int getNumGainGammaCategories(){return gain_multipliers.length-1;}
    public int getNumLossGammaCategories(){return loss_multipliers.length-1;}
    public int getNumDuplicationGammaCategories(){return duplication_multipliers.length-1;}
    public int getNumLengthGammaCategories(){return length_multipliers.length;}
    
    @Override 
    public int getNumClasses()
    {
    	int ng = gain_multipliers.length;
    	int nl = loss_multipliers.length;
    	int nd = duplication_multipliers.length;
    	int ne = length_multipliers.length;
    	
    	return ng*nl*nd*ne;
    }
    
    /**
     * Returns the duplication rate category for this class
     *
     * @param class_idx combined rate class index
     * @return index for duplication rate category
     */
    public final int getDuplicationRateCategory(int class_idx)
    {
        int nd = duplication_multipliers.length;

        return class_idx % nd;
    }
    
    /**
     * Returns the loss rate category for this class
     *
     * @param class_idx combined rate class index
     * @return loss rate category in this class
     */
    public final int getLossRateCategory(int class_idx)
    {
        int nl = loss_multipliers.length;
        int nd = duplication_multipliers.length;
        
        return (class_idx / nd ) % nl;
    }
    
    /**
     * Returns the gain rate category for this class
     *
     * @param class_idx combined rate class index
     * @return gain rate category for this class
     */
    public final int getGainRateCategory(int class_idx)
    {
        int ng = gain_multipliers.length;
        int nl = loss_multipliers.length;
        int nd = duplication_multipliers.length;
        
        return (class_idx / (nd *nl) ) % ng;
    }
    
    /**
     * Returns the edge length category for this class
     *
     * @param class_idx combined class index
     * @return edge length category for this class
     */
    public final int getEdgeLengthCategory(int class_idx)
    {
        int ng = gain_multipliers.length;
        int nl = loss_multipliers.length;
        int nd = duplication_multipliers.length;
        
        return (class_idx / (nd*nl*ng));
    }
    
    /**
     * A null class has 0.0 probability. 
     * 
     * @param class_idx
     * @return true if the class has zero probability 
     */
    public boolean isNullClass(int class_idx)
    {
    	return getClassModel(class_idx)==null;
    }
    
    @Override
    public double getClassProbability(int class_idx)
    {
        int nd = duplication_multipliers.length;
        int dup_idx = class_idx % nd;
        int dup_n= nd-1;
        double dup_p = (dup_idx == dup_n)
                       ?duplication_forbidden
                       :(1.-duplication_forbidden)/dup_n;
        int nl = loss_multipliers.length;
        int loss_idx = (class_idx / nd ) % nl;
        int loss_n = nl-1;
        double loss_p = (loss_idx == loss_n)
                       ?loss_forbidden
                       :(1.-loss_forbidden)/loss_n;
        int ng = gain_multipliers.length;
        int gain_idx = (class_idx / (nd *nl) ) % ng;
        int gain_n = ng-1;
        double gain_p = (gain_idx == gain_n)
                            ?gain_forbidden
                            :(1.-gain_forbidden)/gain_n;
        int edge_n = length_multipliers.length;
        double edge_p = 1.0/edge_n;
        return dup_p * loss_p * gain_p * edge_p;
    }
    
    /**
     * Returns the model for the class; null if class probability is 0
     * (ie not pertinent). 
     */
    @Override 
    public RateModel.GLD getClassModel(int class_idx)
    {
    	return class_rates[class_idx];
    }

    public double getGainAlpha() { return gain_alpha;}
    public double getLossAlpha() { return loss_alpha;}
    public double getDuplicationAlpha() { return duplication_alpha;}
    public double getLengthAlpha() { return length_alpha;}
    
    public void setAlpha(double gain_alpha, double loss_alpha, double duplication_alpha, double length_alpha)
    {
    	this.gain_alpha = gain_alpha;
    	this.loss_alpha = loss_alpha;
    	this.duplication_alpha = duplication_alpha;
    	this.length_alpha = length_alpha;
    	this.setClassRates();
    }
    
    /**
     * Sets the multipliers for Gamma+Zero variation.
     * 
     * @param A multipliers filled in with Gamma values; if there is an extra cell 0.0 is placed there (for forbidden) 
     * @param alpha
     * @param num_gamma_classes
     * @return A
     */
    private static double[] setGammaMultipliers(double[] A, double alpha, int num_gamma_classes)
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
            }
        }
        if (num_gamma_classes != A.length)
        {
            A[A.length-1] = 0.0;
        }
        return A;
    }
    

}
