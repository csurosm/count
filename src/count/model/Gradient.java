
 /* Copyright 2021 Mikl&oacute;s Cs&#369;r&ouml;s.
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


import count.ds.ProfileTable;
import count.ds.IndexedTree;

import static count.model.Likelihood.PARAMETER_GAIN;
import static count.model.Likelihood.PARAMETER_LOSS;
import static count.model.Likelihood.PARAMETER_DUPLICATION;


public class Gradient extends Posteriors
{
	public Gradient(RateModel.GLD rates, ProfileTable table)
	{
		super(rates, table);
	}
	
		
	double[] getSurvivalGradient()
	{
		// first, compute the posterior means and tails
		int num_nodes = factory.tree.getNumNodes();
		double[] node_means = new double[num_nodes];
		double[] edge_means = new double[num_nodes];
		double[][] node_tails = new double[num_nodes][];
		double[][] edge_tails = new double[num_nodes][];
		for (int node=0; node<num_nodes; node++)
		{
			double[] Ntail = node_tails[node] = new double[factory.max_family_size[node]];
			edge_tails[node] = new double[Ntail.length];
		}
		
		int nF = factory.table.getFamilyCount();
		for (int family_idx=0; family_idx<nF; family_idx++)
		{
			Profile post = getPosteriors(family_idx);
			double[] N = post.getNodeMeans();
			double[] S = post.getEdgeMeans();
			for (int node=0; node<num_nodes; node++)
			{
				node_means[node]+=N[node];
				edge_means[node]+=S[node];
				double[] Ncdf = post.getNodeCDF(node);
				double[] Ntail = node_tails[node];
				for (int ell=0; ell<Ncdf.length-1; ell++) // last entry is 1.0
					Ntail[ell]=1.0-Ncdf[ell];
				double[] Scdf = post.getEdgeCDF(node);
				double[] Stail = edge_tails[node];
				for (int s=0; s<Scdf.length-1; s++)
					Stail[s]=1.0-Scdf[s];
			}
		}
		// now compute the gradients 
		
		double[] dL = new double[factory.tree.getNumNodes()*3];// return value
		double logL0 = factory.getEmptyLL();
		double p_not0 = -Math.expm1(logL0);
		
		for (int node=0; node<num_nodes; node++)
		{
			double q = factory.getDuplicationParameter(node);
			if (q==0.0)
			{
				// Poisson
				double r = factory.getGainParameter(node);
				if (r!=0.0)
				{
					dL[3*node+PARAMETER_GAIN] 
							= (node_means[node] - edge_means[node])/r - nF/p_not0; 
				}
			} else
			{
				// Pólya
				double κ = factory.getGainParameter(node);
				if (κ!=0.0)
				{
					double dLdk = 0.0;
					double[] Ntail = node_tails[node];
					double[] Stail = edge_tails[node];
					for (int i=0; i<Ntail.length; i++)
					{
						dLdk += (Ntail[i]-Stail[i])/(κ + i);
					}
					double log1_q = Math.log1p(-q);
					dL[3*node+PARAMETER_GAIN] 
							= dLdk + nF * log1_q/p_not0;
				}
				dL[3*node+PARAMETER_DUPLICATION]
						= (node_means[node]-edge_means[node])/q
						- (edge_means[node]+ κ/p_not0)/(1.-q);
			}
			if (!factory.tree.isRoot(node))
			{
				double p = factory.getLossParameter(node);
				int parent = factory.tree.getParent(node);
				double epsi = factory.extinction[parent]/p;
				dL[3*node+PARAMETER_LOSS]
						= node_means[node]/p - (1.0-epsi)*edge_means[node]/(1.0-p);
			}
		}

		return dL;
	}
	
	public double[] getGradient()
	{
		double[] dL = getSurvivalGradient();
		int num_nodes = factory.tree.getNumNodes();
		double[] de = new double[num_nodes];
		
		int node = factory.tree.getRoot();
		double q = factory.rates.getDuplicationParameter(node);
		if (q==0.0)
		{
			// Poisson
			double r = factory.rates.getGainParameter(node);
			if (r!=0.0)
			{
				double dLdr =  dL[3*node+PARAMETER_GAIN];
				dL[3*node+PARAMETER_GAIN] = dLdr*(1.0-factory.extinction[node]); 
				de[node] = dLdr*(-r);
			}
		} else
		{
			// Pólya or geometric
			double dLdq = dL[3*node+PARAMETER_DUPLICATION];
			double a = 1.0-q*factory.extinction[node];
			a = a*a;
			dL[3*node+PARAMETER_DUPLICATION] = dLdq * (1.0-factory.extinction[node])/a;
			de[node] = -dLdq*(1.0-q)*q/a;
		}
		while (node>0)
		{
			--node;
			int parent = factory.tree.getParent(node);
			q = factory.rates.getDuplicationParameter(node);
			double p = factory.rates.getLossParameter(node);
			double epsi = factory.extinction[parent]/factory.getLossParameter(node);
			double dLdpe = dL[3*node+PARAMETER_LOSS]+epsi*de[parent];
			if (q==0.0)
			{
				// Poisson
				double r = factory.rates.getGainParameter(node);
				if (r!=0.0)
				{
					double dLdr =  dL[3*node+PARAMETER_GAIN];
					dL[3*node+PARAMETER_GAIN] = dLdr * (1.0-factory.extinction[node]); 
					dL[3*node+PARAMETER_LOSS] = dLdpe * (1.0-factory.extinction[node]);
					de[node] = dLdpe * (1.0-p) - dLdr * r;
				}
			} else
			{
				// Pólya or geometric
				double a = 1.0-q*factory.extinction[node];
				dL[3*node+PARAMETER_LOSS] = dLdpe * (1.0-factory.extinction[node]) / a;
				double dLdq = dL[3*node+PARAMETER_DUPLICATION];
				a = a*a;
				dL[3*node+PARAMETER_DUPLICATION] = (dLdq - dLdpe * (1.0-p) * factory.extinction[node]) 
							* (1.0-factory.extinction[node])/a;
				de[node] = (dLdpe * (1.0-p) - dLdq * q) * (1.0-q)/a;
			}
		}
		return dL;
	}
	
	public double[] getLogisticGradient()
	{
		double[] dL = getGradient();
		int num_nodes = factory.tree.getNumNodes();
		
		for (int node=0; node<num_nodes; node++)
		{
			double q = factory.rates.getDuplicationParameter(node);
			dL[3*node+PARAMETER_DUPLICATION] *= q*(1.0-q);
			if (!factory.tree.isRoot(node))
			{
				double p = factory.rates.getLossParameter(node);
				dL[3*node+PARAMETER_LOSS] *= p*(1.0-p);
			}
			dL[3*node+PARAMETER_GAIN] *= factory.rates.getGainParameter(node);
		}
		return dL;
	}
	
}
