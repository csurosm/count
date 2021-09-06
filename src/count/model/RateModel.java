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
import count.matek.DiscreteDistribution;
import count.matek.NegativeBinomial;
import count.matek.Poisson;
import count.matek.PointDistribution;
import count.matek.ShiftedGeometric;

/**
 * Public interface used by the algorithms:  
 * gain, loss and duplication rates, and edge lengths over a 
 * tree. 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 */
public interface RateModel
{
	/**
	 * The underlying tree. 
	 * 
	 * @return
	 */
	public abstract IndexedTree getTree();
	
	/**
	 * Gain rate (<var>κ</var>) on the edge leading to a node, or for root prior. 
	 * 
	 * @param node_idx any node 
	 * @return non-negative gain rate 
	 */
	public abstract double getGainRate(int node_idx);

	/**
	 * Loss rate (<var>μ</var>) on the edge leading to a node.
	 * 
	 * @param node_idx non-leaf node
	 * @return non-negative loss rate 
	 */
	public abstract double getLossRate(int node_idx);

	/**
	 * Duplication rate (<var>λ</var>) on the edge leading to a node, or for root prior. 
	 * 
	 * @param node_idx any node 
	 * @return non-negative duplication rate 
	 */
	public abstract double getDuplicationRate(int node_idx);

	
	/**
	 * Edge length (<var>t</var>) on the edge leading to a node, or for root prior. 
	 * 
	 * @param node_idx any node 
	 * @return non-negative gain rate 
	 */
	public abstract double getEdgeLength(int node_idx);
	
	public interface GLD extends RateModel, GLDParameters
	{
		//
		// GLD parameric interface
		// 
		@Override
		public default double getGainParameter(int node_idx)
		{
			double gain_rate = getGainRate(node_idx);
			double gainParameter;
			if (getLossRate(node_idx)==0.0 || getDuplicationRate(node_idx)!=0.0 || gain_rate==0.0)
				gainParameter = gain_rate;
			else
				gainParameter = gain_rate * getLossParameter(node_idx); // gain-loss-noduplication = Poisson
			return gainParameter;
		}
		
		@Override
		public default double getLossParameter(int node_idx)
		{
			double loss_rate = getLossRate(node_idx);
			double lossParameter;
			if (loss_rate==0.0)
				lossParameter = 0.0;
			else
			{

				double dup_rate = getDuplicationRate(node_idx);
				double t = getEdgeLength(node_idx);
				if (loss_rate == dup_rate)
				{
					double mu_t = loss_rate*t;
					lossParameter = mu_t/(1.0+mu_t);
				} else
				{
					if (dup_rate == 0.0)
					{
						// p = 1-e^{-mu t}
						lossParameter = -Math.expm1(-loss_rate*t);
					} else
					{
						double d = (loss_rate-dup_rate)*t; // ok if negative
						double E = Math.exp(-d);
						double E1 = -Math.expm1(-d);
						// p = (mu-mu*E)/(mu-lambda*E)
						lossParameter = loss_rate * E1/(loss_rate - dup_rate*E);
					}
				}
			}
			return lossParameter;
		}
		
		@Override
		public default double getDuplicationParameter(int node_idx)
		{
			double dup_rate = getDuplicationRate(node_idx);
			double t = getEdgeLength(node_idx);
			double duplicationParameter;
			if (dup_rate==0.0)
				duplicationParameter = 0.0;
			else 
			{
				double loss_rate = getLossRate(node_idx);
				if (loss_rate == dup_rate)
				{
					double lambda_t = dup_rate*t;
					duplicationParameter = lambda_t/(1.0+lambda_t);
				} else
				{
					double d = (loss_rate-dup_rate)*t;
					double E = Math.exp(-d);
					double E1 = -Math.expm1(-d);
					// q = (lambda-lambda*E)/(mu-lambda*E)
					duplicationParameter = dup_rate * E1/(loss_rate - dup_rate*E);
				}
			}
			return duplicationParameter;
		}
	}
	
}
