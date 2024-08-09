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
 * Not used anymore.
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
 * 
 */
public class GammaInvariant implements MixedRateModel.RateMultipliers
{
    /**
     * Instantiates a rate-variation model. In such a model, 
     * a cross-product of discrete Gamma categories is used for
     * rates (duplication,loss,transfer) and branch lengths.
     * For parameters on which variation is not imposed, the number of categories should be 1.
     * @param base_rates tree with lineage-specific parameters
     * @param duplication_rate_categories number of duplication rate categories for rate variation across families
     * @param loss_rate_categories number of loss rate categories for rate variation across families
     * @param gain_rate_categories number of gain rate categories for rate variation across families
     * @param edge_length_categories number of edge length categories for rate variation across families
     */
    public GammaInvariant(TreeWithRates base_rates,
            int gain_rate_categories, 
            int loss_rate_categories, 
            int duplication_rate_categories, 
            int edge_length_categories) 
    {
        this(base_rates);
        setClasses(
            gain_rate_categories, 
            loss_rate_categories, 
            duplication_rate_categories, 
            edge_length_categories);
    }
    
    /**
     * Instantiation for setting up without specifying the classes. 
     * @param base_rates tree with lineage-specific parameters
     */
    public GammaInvariant(TreeWithRates base_rates)
    {
        this.base_rates = base_rates;
        this.main_tree = base_rates.getTree();
        // we don't know the classes yet
    } 
    
    /**
     * Null model with 1 rate class. 
     * @param main_tree
     * @return
     */
    public static GammaInvariant nullModel(IndexedTree main_tree)
    {
    	return nullModel(main_tree, null);
    }
    
    public static GammaInvariant nullModel(IndexedTree main_tree, java.util.Random RND)
    {
        TreeWithRates null_rates = new TreeWithRates(main_tree, RND);
        GammaInvariant null_model = new GammaInvariant(null_rates, 1,1,1,1);
    	return null_model;
    }
    
    private TreeWithRates base_rates;
    private IndexedTree main_tree;
    private double[] gain_multipliers;
    private double[] loss_multipliers;
    private double[] duplication_multipliers;
    private double[] length_multipliers;
    
//    private TreeWithRates[] class_rates; // parametric type is lost  in array, but these are phylogenies
    
    private double gain_alpha=1.0;
    private double loss_alpha=1.0;
    private double duplication_alpha=1.0;
    private double length_alpha=1.0;
    
    private double gain_forbidden = 0.0;
    private double loss_forbidden = 0.0;
    private double duplication_forbidden = 0.0;
    
    
    @Override
    public TreeWithRates getBaseModel()
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
    	
//    	class_rates = new TreeWithRates[getNumClasses()]; 
//    	setClassRates(); // adjusts the rates on the class trees 
    }

    public double getGainForbidden() { return gain_forbidden;}
    public double getLossForbidden() { return loss_forbidden;}
    public double getDuplicationForbidden() {return duplication_forbidden;}
    
    /**
     * Sets the proportions of sites where the invariant model (0 factor) applies instead of Gamma variation
     * 
     * @param gain_forbidden proportion of families with no transfer (between 0.0 and 1.0)
     */
    public void setGainForbidden(double gain_forbidden)
    {
    	this.gain_forbidden = gain_forbidden;
    }
    
    /**
     * Sets the proportions of sites where the invariant model (0 factor) applies instead of Gamma variation
     * 
     * @param loss_forbidden proportion of families with no loss (between 0.0 and 1.0)
     */
    public void setLossForbidden(double loss_forbidden)
    {
    	this.loss_forbidden = loss_forbidden;
    }
    
    /**
     * Sets the proportions of sites where the invariant model (0 factor) applies instead of Gamma variation
     * 
     * @param duplication_forbidden proportion of families with no duplication (between 0.0 and 1.0)
     */
    public void setDuplicationForbidden(double duplication_forbidden)
    {
    	this.duplication_forbidden = duplication_forbidden;
    }
    
//    private void setClassRates()
//    {
//    	int ng = gain_multipliers.length;
//    	int nl = loss_multipliers.length;
//    	int nd = duplication_multipliers.length;
//    	int ne = length_multipliers.length;
//    	
//    	int num_nodes = main_tree.getNumNodes();
//    	
//    	int class_idx = 0;
//    	// class_idx =  dup_idx + nd*(loss_idx + nl*(gain_idx + ng*length_idx)
//		for (int length_idx=0; length_idx<ne; length_idx++)
//		{
//			double xe = length_multipliers[length_idx];
//			int gain_idx = 0;
//			while (gain_idx<ng)
//			{
//				boolean is_excluded = (gain_idx==ng-1 && gain_forbidden==0.0);
//				double xg = gain_multipliers[gain_idx];
//				++gain_idx;
//	    		int loss_idx = 0;
//	    		while (loss_idx<nl)
//	    		{
//	    			is_excluded = is_excluded || (loss_idx==nl-1 && loss_forbidden==0.0);
//	    			double xl = loss_multipliers[loss_idx];
//	    			++loss_idx;
//			    	int dup_idx = 0;
//			    	while (dup_idx<nd)
//		    		{
//			    		is_excluded  = is_excluded || (dup_idx==nd-1 && duplication_forbidden == 0.0);
//			            double xd = duplication_multipliers[dup_idx];
//			            ++dup_idx;
//    					TreeWithRates R;
//    					if (is_excluded)
//    					{
//    						R=null;
//    					} else
//    					{
//    						R = new TreeWithRates(base_rates);
//    						if (xl!= 1.0)
//    						{
//	    						for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // root keeps the base rates 
//	    						{
//	    							double yl = xl*R.getLossRate(node_idx);
//	    							R.setLossRate(node_idx, yl);
//	    						}
//    						}
//    						if (xd != 1.0)
//    						{
//	    						for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // root keeps the base rates 
//	    						{
//	    							double yd = xd*R.getDuplicationRate(node_idx);
//	    							R.setDuplicationRate(node_idx, yd);
//	    						}
//    						}
//    						for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // root keeps the base rates 
//    						{
//    							double yg = xg*R.getGainRate(node_idx);
//    							
//    							// legacy logic for multi-rate variations 
//    							double yd = R.getDuplicationRate(node_idx);
//    							if (yd==0.0)
//    							{
//    								double yl = R.getLossRate(node_idx);
//    								if (yl==0.0)	
//    									R.setGainRate(node_idx, yg);
//    								else 
//    									R.setGainRate(node_idx, yg*base_rates.getLossRate(node_idx)/yl);
//    							} else
//    							{
//    								R.setGainRate(node_idx,  yg*base_rates.getDuplicationRate(node_idx)/yd);
//    							}
//    						}
//    						if (xe != 1.0)
//    						{
//	    						for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // root keeps the base rates 
//	    						{
//	    							double ye = xe*R.getEdgeLength(node_idx);
//	    							R.setEdgeLength(node_idx, ye);
//	    						}
//    						}
//    					}
////    					System.out.println("#*class "+class_idx+"\t"+R);
//    					class_rates[class_idx++] = R;
//    				} // while dup
//    			} // while loss
//    		} // while gain
//    	} //  for length
//    }
    
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
     * Number of classes that have positive probability. 
     * @return
     */
    public int  getNumActiveClasses()
    {
    	int nc = getNumClasses();
    	int num_active_classes = 0;
    	for (int c=0; c<nc; c++)
    	{
    		if (getClassProbability(c)>0.0)
    			num_active_classes++;
    	}
    	assert (num_active_classes>0);
    	return num_active_classes;
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
    
    public double getGainRateMultiplier(int class_idx)
    {
    	return gain_multipliers[getGainRateCategory(class_idx)];
    }
    
    public double getLossRateMultiplier(int class_idx)
    {
    	return loss_multipliers[getLossRateCategory(class_idx)];
    }
    
    public double getDuplicationRateMultiplier(int class_idx)
    {
    	return duplication_multipliers[getDuplicationRateCategory(class_idx)];
    }
    
    public double getEdgeLengthMultiplier(int class_idx)
    {
    	return length_multipliers[getEdgeLengthCategory(class_idx)];
    }
    
//    /**
//     * A null class has 0.0 probability. 
//     * 
//     * @param class_idx
//     * @return true if the class has zero probability 
//     */
//    public boolean isNullClass(int class_idx)
//    {
//    	return getClassProbability(class_idx)==0.0;
//    			
//    			//getClassModel(class_idx)==null;
//    }
    
    @Override
    public double getClassProbability(int class_idx)
    {
        int nd = duplication_multipliers.length;
        int dup_idx = getDuplicationRateCategory(class_idx); // class_idx % nd;
        int dup_n= nd-1;
        double dup_p = (dup_idx == dup_n)
                       ?duplication_forbidden
                       :(1.-duplication_forbidden)/dup_n;
        int nl = loss_multipliers.length;
        int loss_idx = getLossRateCategory(class_idx); // (class_idx / nd ) % nl;
        int loss_n = nl-1;
        double loss_p = (loss_idx == loss_n)
                       ?loss_forbidden
                       :(1.-loss_forbidden)/loss_n;
        int ng = gain_multipliers.length;
        int gain_idx = getGainRateCategory(class_idx); //(class_idx / (nd *nl) ) % ng;
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
    public TreeWithRates getClassModel(int class_idx)
    {
//    	return class_rates[class_idx];
    	
    	double xg = getGainRateMultiplier(class_idx);
    	double xl = getLossRateMultiplier(class_idx);
    	double xd = getDuplicationRateMultiplier(class_idx);
    	double xe = getEdgeLengthMultiplier(class_idx);
    	
		TreeWithRates R = new TreeWithRates(base_rates);
		int num_nodes = main_tree.getNumNodes();
		
		for (int node=0; node<num_nodes-1; node++) // not for root==num_nodes-1
		{
			double yl = xl*R.getLossRate(node);
			double yd = xd*R.getDuplicationRate(node);
			double yg = xg*R.getGainRate(node);
			double ye = xe*R.getEdgeLength(node);

			// legacy logic for gain-rate variations 
			if (yd==0.0)
			{
				if (yl!=0.0)	
				{
					yg = yg*base_rates.getLossRate(node)/yl;
				}
			} else
			{
				yg = yg*base_rates.getDuplicationRate(node)/yd;
			}
			R.setRates(node, ye, yg, yl, yd);
		}
		
//		
//		if (xl!= 1.0)
//		{
//			for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // root keeps the base rates 
//			{
//				double yl = xl*R.getLossRate(node_idx);
//				R.setLossRate(node_idx, yl);
//			}
//		}
//		if (xd != 1.0)
//		{
//			for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // root keeps the base rates 
//			{
//				double yd = xd*R.getDuplicationRate(node_idx);
//				R.setDuplicationRate(node_idx, yd);
//			}
//		}
//		for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // root keeps the base rates 
//		{
//			double yg = xg*R.getGainRate(node_idx);
//			
//			// legacy logic for multi-rate variations 
//			double yd = R.getDuplicationRate(node_idx);
//			if (yd==0.0)
//			{
//				double yl = R.getLossRate(node_idx);
//				if (yl==0.0)	
//					R.setGainRate(node_idx, yg);
//				else 
//					R.setGainRate(node_idx, yg*base_rates.getLossRate(node_idx)/yl);
//			} else
//			{
//				R.setGainRate(node_idx,  yg*base_rates.getDuplicationRate(node_idx)/yd);
//			}
//		}
//		if (xe != 1.0)
//		{
//			for (int node_idx=0; node_idx<num_nodes-1; node_idx++) // root keeps the base rates 
//			{
//				double ye = xe*R.getEdgeLength(node_idx);
//				R.setEdgeLength(node_idx, ye);
//			}
//		}
		return R;
    }

    public double getGainAlpha() { return gain_alpha;}
    public double getLossAlpha() { return loss_alpha;}
    public double getDuplicationAlpha() { return duplication_alpha;}
    public double getLengthAlpha() { return length_alpha;}
    
    public void setGainAlpha(double gain_alpha)
    {
    	this.gain_alpha = gain_alpha;
    	gain_multipliers 
		= setGammaMultipliers(gain_multipliers,gain_alpha,gain_multipliers.length-1);
    }
    
    public void setLossAlpha(double loss_alpha)
    {
    	this.loss_alpha = loss_alpha;
    	loss_multipliers 
		= setGammaMultipliers(loss_multipliers,loss_alpha,loss_multipliers.length-1);
    }
    
    public void setDuplicationAlpha(double duplication_alpha)
    {
    	this.duplication_alpha = duplication_alpha;
    	duplication_multipliers 
		= setGammaMultipliers(duplication_multipliers,duplication_alpha,duplication_multipliers.length-1);
    }
    
    public void setLengthAlpha(double length_alpha)
    {
    	this.length_alpha = length_alpha;
    	length_multipliers 
    	= setGammaMultipliers(length_multipliers, length_alpha, length_multipliers.length);
    }
    
    /**
     * Creates a copy of this model with the same underlying phylogeny. 
     * 
     * @return
     */
    public GammaInvariant copy()
    {
    	TreeWithRates base_copy = new TreeWithRates(base_rates);
    	GammaInvariant copy = new GammaInvariant(base_copy, 
    			getNumGainGammaCategories(),
    			getNumLossGammaCategories(),
    			getNumDuplicationGammaCategories(),
    			getNumLengthGammaCategories());
    	copy.setDuplicationAlpha(getDuplicationAlpha());
    	copy.setLossAlpha(getLossAlpha());
    	copy.setGainAlpha(getGainAlpha());
    	copy.setLengthAlpha(getLengthAlpha());
    	copy.setDuplicationForbidden(getDuplicationForbidden());
    	copy.setGainForbidden(getGainForbidden());
    	copy.setLossForbidden(getLossForbidden());
    	
    	return copy;
    }
   
    
    /**
     * Sets the multipliers for Gamma+Zero variation.
     * 
     * @param A multipliers filled in with Gamma values; if there is an extra cell 0.0 is placed there (for forbidden) 
     * @param alpha
     * @param num_gamma_classes (=A.length or A.length-1 for 0.0 in last cell)
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
                A[i] = val[i]; // val[i]/val[val.length-1];
            }
        }
        if (num_gamma_classes != A.length)
        {
            A[A.length-1] = 0.0;
        }
        return A;
    }
    

}
