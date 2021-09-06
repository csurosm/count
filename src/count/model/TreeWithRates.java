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

import java.util.Arrays;

import count.ds.IndexedTree;
import count.ds.Phylogeny;
import count.matek.DiscreteDistribution;
import count.matek.NegativeBinomial;
import count.matek.PointDistribution;
import count.matek.Poisson;
import count.matek.ShiftedGeometric;

public class TreeWithRates implements RateModel.GLD
{
	public TreeWithRates(IndexedTree tree)
	{
		this.tree = tree;
		int num_nodes = tree.getNumNodes();
		this.gain_rates = new double[num_nodes];
		this.loss_rates = new double[num_nodes]; 
		this.duplication_rates = new double[num_nodes];
		this.edge_lengths = new double[num_nodes];
		Arrays.fill(gain_rates, DEFAULT_GAIN_RATE);
		Arrays.fill(loss_rates, DEFAULT_LOSS_RATE);
		Arrays.fill(duplication_rates, DEFAULT_DUPLICATION_RATE);
		if (tree instanceof Phylogeny)
			scaleByPhylogeny();
		else
			Arrays.fill(edge_lengths, 1.0);
	}
		
	/**
	 * Creates a copy linked to the same phylogeny; rate parameters are copied. 
	 * @param same_phylogeny base model
	 */
	public TreeWithRates(RateModel same_phylogeny)
	{
		this.tree = same_phylogeny.getTree(); 
		int num_nodes = tree.getNumNodes();
		this.gain_rates = new double[num_nodes];
		this.loss_rates = new double[num_nodes]; 
		this.duplication_rates = new double[num_nodes];
		this.edge_lengths = new double[num_nodes];
		for (int node=0; node<num_nodes; node++)
		{
			gain_rates[node] = same_phylogeny.getGainRate(node);
			loss_rates[node] = same_phylogeny.getLossRate(node);
			duplication_rates[node] = same_phylogeny.getDuplicationRate(node);
			edge_lengths[node] = same_phylogeny.getEdgeLength(node);
		}  
	}

	
	private IndexedTree tree;
	private final double[] gain_rates;
	private final double[] loss_rates;
	private final double[] duplication_rates;
	private final double[] edge_lengths;
	
    public static final double DEFAULT_GAIN_RATE = 0.2;
    public static final double DEFAULT_LOSS_RATE = 1.0;
    public static final double DEFAULT_DUPLICATION_RATE = 0.5;
    public static final double DEFAULT_EDGE_LENGTH = 1.0;
    
    
    /**
     * Powers of ten.
     */
    private static final double[] POW10 = {1.,10.,100.,1000.,10000.,100000.,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19,1e20,1e21,1e22,1e23,1e24,1e25,1e26,1e27,1e28,1e29,1e30,1e31};

    private void scaleByPhylogeny()
	{
    	Phylogeny phylo = (Phylogeny) tree;
    	// rescale edge lengths 
    	int num_nodes = phylo.getNumNodes(); 
    	int node_idx = num_nodes;
    	assert node_idx>0;
    	double d[] = new double[node_idx];// depths in the tree
    	--node_idx;
    	d[node_idx] = 0.0 ; // root
    	double dmax = 0.0;
    	while (node_idx>0)
    	{
    		--node_idx;
    		int parent_idx = tree.getParent(node_idx);
    		double parent_depth = d[parent_idx];
    		double node_depth = parent_depth + phylo.getLength(node_idx);
    		if (node_depth > dmax) dmax = node_depth;
    		d[node_idx]  = node_depth;
    	}
    	// find a normalizing factor
        double factor = 1.0;
        int dlog = (int)Math.log10(dmax);
        if (dlog>0) 
        {
        	while (dlog>=POW10.length){factor *= 0.1;--dlog;}
        	factor *= 1.0/POW10[dlog];
        } else if (dlog<0)
        {
        	while (dlog<=-POW10.length) {factor *=10.0; ++dlog;}
        	factor *= POW10[dlog];
        }
        while (node_idx<num_nodes)
        {
        	edge_lengths[node_idx] = phylo.getLength(node_idx)*factor;
        	++node_idx;
        }
	}
	
	@Override
	public IndexedTree getTree()
	{
		return tree;
	}

	@Override
	public double getGainRate(int node_idx)
	{
		return gain_rates[node_idx];
	}

	public void setGainRate(int node_idx, double gain_rate)
	{
		gain_rates[node_idx]=gain_rate;
	}

	@Override
	public double getLossRate(int node_idx)
	{
		return loss_rates[node_idx];
	}

	public void setLossRate(int node_idx, double loss_rate)
	{
		loss_rates[node_idx]=loss_rate;
	}

	@Override 
	public double getDuplicationRate(int node_idx)
	{
		return duplication_rates[node_idx];
	}
	
	public void setDuplicationRate(int node_idx, double dup_rate)
	{
		duplication_rates[node_idx]=dup_rate;
	}
	
	
	@Override
	public double getEdgeLength(int node_idx)
	{
		return edge_lengths[node_idx];
	}
	
	public void setEdgeLength(int node_idx, double length)
	{
		edge_lengths[node_idx]=length;
	}

	
	
	/**
	 * Sets the gain, loss and duplication rates. 
	 * 
	 * @param node_idx
	 * @param gain_param
	 * @param loss_param
	 * @param dup_param
	 */
	public void setParameters(int node_idx, double gain_param, double loss_param, double dup_param)
	{
		
		double loss_rate, dup_rate;

		if (dup_param == loss_param)
		{ // includes ==0.0 
			dup_rate = loss_rate = loss_param/(1.0+loss_param);
		} else
		{		
			if (loss_param == 0.0)
			{
				loss_rate = 0.0;
				dup_rate = -Math.log(1.0-dup_param);
			} else if (dup_param==0.0)
			{
				dup_rate = 0.0;
				loss_rate = -Math.log(1.0-loss_param);
				
			} else
			{
				double dl_ratio = dup_param/loss_param;
				double delta = 1.0-dl_ratio;
				loss_rate = Math.log((1.0-dup_param)/(1.0-loss_param))/delta;
				dup_rate = loss_rate * dl_ratio;
			}
		}
		double gain_rate;
		if (gain_param==0.0)
		{
			gain_rate = 0.0;
		} else
		{
			if (dup_param == 0.0 && loss_param != 0.0)
				gain_rate = gain_param/loss_param;
			else
				gain_rate = gain_param;
		}
		double t = getEdgeLength(node_idx);
		
		setGainRate(node_idx, gain_rate);
		setLossRate(node_idx, loss_rate/t);
		setDuplicationRate(node_idx, dup_rate/t);
	}
}
