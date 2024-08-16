package count.model;
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

import count.ds.IndexedTree;

/**
 * Likelihood computations for ancestral reconstruction at a node.
 * The trick is to set extinction at 0.0 for the node where the 
 * ancestor copies are to be observed.
 *  
 * 
 * 
 * @author csuros
 *
 */
public class Ancestor extends LikelihoodParametrized
{
	public Ancestor(Likelihood factory, int node)
	{
		super(factory.rates, factory.table);
		this.ancestor = node;
		setSameCalculationWidthThresholds(factory);
		computeParameters();
	}
	
	private int ancestor=-1; 
	private int ancestor_calculation_width = 2;
	
	public void setCalculationWidth(int width)
	{
		this.ancestor_calculation_width = width;
	}
	
	@Override
	public void computeParameters()
	{
		if (ancestor>=0) // but not when super's instantiation calls and ==-1
			super.computeParameters();
	}
	
	@Override
	public void computeNodeParameters(int node)
	{
		if (node==ancestor) 
		{
			// TODO this.computeNodeParameters(node, Double.NEGATIVE_INFINITY); // logit(0.0)
			super.setNodeParameters(node, 0.0);
		} else
			super.computeNodeParameters(node);
	}
	
	@Override
	protected Likelihood.Profile getProfileLikelihood(int family_idx)
    {
    	return this.new Profile(family_idx);
    }	
	
	private class Profile extends LikelihoodParametrized.Profile
	{
		private Profile(int family_idx)
		{
			super(family_idx);
		}
		
		/**
		 * Calls super's method but 
		 * extends the calculation widths to include 
		 * copy number 1 at the ancestor.
		 */
		@Override
		public void computeLikelihoods()
		{
			IndexedTree tree = rates.getTree();
			
			int diff_width = ancestor_calculation_width-getCalculationWidth(ancestor);
			// we need [0], [1], ..., [width-1] on the root to the root
			if (diff_width>0)
			{
				int node = ancestor;
				do
				{
					int w = getCalculationWidth(node);
					setCalculationWidth(node, w+diff_width); 
					
					node = tree.getParent(node);
				} while (node >= 0);
			} else
			{
				
			}
			super.computeLikelihoods();
		}
		
		@Override 
		public void copyLikelihoods(Likelihood.Profile that)
		{
			IndexedTree tree = rates.getTree();
			int diff_width = ancestor_calculation_width-getCalculationWidth(ancestor);
			if (diff_width!=0)
			{
				setCalculationWidth(ancestor, ancestor_calculation_width);
				int node = tree.getParent(ancestor);
				while (node>=0)
				{
					setCalculationWidth(node);
					node = tree.getParent(node);
				}
			}
			super.copyLikelihoods(that);
		}
	}

}
