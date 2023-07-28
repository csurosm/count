package count.model;
/*
 * Copyright 2023 Mikl&oacute;s Cs&#369;r&ouml;s.
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



import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.TreeComparator;
import count.io.CommandLine;
import count.matek.Logarithms;

/**
 * Joint evaluation of a set of rate models by bootstrapping. 
 * Add models one by one using {@link #addRates(TreeWithRates)},
 * then call {@link #getBootstrapWinners(Random, int)} and {@link #getBootstrapSupport(int[], TreeComparator)}
 * to translate into edge support for a chosen phylogeny.  
 */
public class RELL 
{
	public RELL(AnnotatedTable table)
	{
		this.table = table;
		this.models = new ArrayList<>();
		this.modelLL = new ArrayList<>();
		this.min_copies = Integer.min(2,table.minCopies());

	}
	
	private final AnnotatedTable table;
//	private final List<AnnotatedTable> projected_tables;
	private final List<TreeWithRates> models;
	private final List<double[]> modelLL;
	
//	private TreeComparator cmp;
	
	private int min_copies;
	/**
	 * Sets the observation bias : minimum number of observed copies in the families. 
	 * @param min_copies 0, 1, or 2
	 */
	public void setMinimumObservedCopies(int min_copies)
	{
		if (min_copies<0 || min_copies>2)
			throw new UnsupportedOperationException("Minimum copies must be 0,1 or 2.");
		this.min_copies = min_copies;
	}
	
	protected int getMinimumObservedCopies()
	{
		return min_copies;
	}
	
	public void clear()
	{
		models.clear();
//		cmp = null;
		modelLL.clear();
	}
	
	public void addRates(TreeWithRates rates)
	{
		models.add(rates);
		IndexedTree rates_tree = rates.getTree();
		AnnotatedTable rates_table = table.mappedToTree(rates_tree);
//		if (models.size()==1)
//		{
//			cmp = new TreeComparator(rates_tree);
//		} 
		double[] fLL = familyLL(rates, rates_table);
		modelLL.add(fLL);
	}
	
	private double getUnobservedLL(Likelihood factory)
	{
		double unobservedLL = Double.NEGATIVE_INFINITY;
		if (min_copies>0)
		{
			unobservedLL = factory.getEmptyLL();
			if (min_copies>1)
			{
				assert (min_copies==2);
				unobservedLL = Logarithms.add(unobservedLL, factory.getSingletonLL());
			}
		}
		return unobservedLL;
	}
	
	private double[] familyLL(TreeWithRates rates, AnnotatedTable rates_table)
	{
		double[] familyLL = new double[table.getFamilyCount()];
		Likelihood factory = new Likelihood(rates,rates_table);
		double log_punobs = getUnobservedLL(factory) ;
		double p_obs = -Math.expm1(log_punobs); 
		double log_pobs = Math.log(p_obs);
		
		for (int f=0; f<familyLL.length; f++)
		{
			Likelihood.Profile P = factory.getProfileLikelihood(f);
			familyLL[f] = P.getLogLikelihood()-log_pobs; // corrected likelihood
			
//			if (f<100 || f%1000==0) // DEBUG
//			{
//				System.out.println("#**RELL.fLL "+f+"\t"+familyLL[f]);
//			}
		}
		
		
		return familyLL;
	}
	
	/**
	 * Log-likelihood over the original table.
	 * 
	 * @return array of log-likelihoods per model
	 */
	public double[] getCorrectedLL()
	{
		
		double[] LL = new double[models.size()];
		for (int model_idx=0; model_idx<LL.length; model_idx++)
		{
			double[] fLL = modelLL.get(model_idx);
			for (int f=0; f<fLL.length; f++)
				LL[model_idx] += fLL[f];
		}
		return LL;
	}
	
	private static double sum(double[] LL, int[] multiplicities)
	{
		double sum = 0.0;
		for (int f=0; f<LL.length; f++)
			sum = Math.fma(multiplicities[f], LL[f], sum);
		return sum;
	}
	
	/**
	 * Log-likelihood of a random bootstrap sample across 
	 * tyhe different models. 
	 * 
	 * @param RND pseudorandom number generator
	 * @return array of log-likelihood for the added rate models 
	 */
	public double[] getBootstrapLL(Random RND)
	{
		int nF = table.getFamilyCount();
		int[] bootstrap = new int[nF]; // number of times each family is selected
		for (int f=0; f<nF; f++)
		{
			int random_f = RND.nextInt(nF);
			bootstrap[random_f]++;
		}
		double[] bootstrapLL = new double[models.size()];
		for (int model_idx=0; model_idx<models.size(); model_idx++)
		{
			double mLL = sum(modelLL.get(model_idx), bootstrap);
			bootstrapLL[model_idx] = mLL;
		}
		return bootstrapLL;
	}
	
	/**
	 * Translates winning statistics into bootstrap edge support for first model, 
	 * for the unrooted phylogeny. 
	 * 
	 * @param winners array of bootstrap winning statistics from {@link #getBootstrapWinners(Random, int)} 
	 * @param cmp the reference tree whose edges are to be annotated (use one of the tested models)
	 * @return bootstrap support: how many times each edge was present in bootstrap winners 
	 */
	public int[] getBootstrapSupport(int[] winners, TreeComparator cmp)
	{
		assert (winners.length == models.size());
		int[][] node_maps = new int[models.size()][];
		int num_nodes = cmp.getReferenceTree().getNumNodes();
//		int[] ident = new int[num_nodes];
//		for (int node=0; node<num_nodes; node++)
//			ident[node] = node;
//		node_maps[0] = ident;
		for (int m=0; m<models.size(); m++)
		{
			IndexedTree tree = models.get(m).getTree();
			node_maps[m] = cmp.map(tree).toReference();
		}
		int[] support = new int[num_nodes];
		for (int m=0; m<winners.length; m++)
			if (winners[m]!=0)
			{
				int[] ref_node = node_maps[m];
				for (int node=0; node<ref_node.length; node++)
				{
					int rnode = ref_node[node];
					if (rnode!=-1)
						support[rnode]+=winners[m];
				}
			}
		return support;
	}
	
	/**
	 * Statistics on how many times each model is best on bootstrap samples
	 * 
	 * @param RND ranomd number generator
	 * @param num_bootstrap_samples so many samples are generated randomly (with replacement from original families)
	 * @return array[m] is the number of times model m was best
	 */
	
	public int[] getBootstrapWinners(Random RND, int num_bootstrap_samples)
	{
		int[] winners = new int[models.size()];
		for (int smp = 0; smp<num_bootstrap_samples; smp++)
		{
			double[] smpLL = getBootstrapLL(RND);
			int w = -1;
			double maxLL = Double.NEGATIVE_INFINITY;
			for (int m=0; m<models.size(); m++)
			{
				if (smpLL[m]>maxLL)
				{
					maxLL= smpLL[m];
					w = m;
				}
			}
			winners[w]++;
		}
		return winners;
	}
	
	public static void main(String[] args) throws Exception
	{
		
		CommandLine cli = new CommandLine(args, RELL.class);
		AnnotatedTable table= cli.getTable();
		TreeWithRates rates = cli.getModel().getBaseModel();
		
		RELL R = new RELL(table);
		R.addRates(rates);
		
		Random RND = cli.getOptionRND(System.out);
		if (RND==null) RND = new Random(2023);
		
		double[] bLL = R.getBootstrapLL(RND);
		
		System.out.println("#**RELL.main bLL "+Arrays.toString(bLL));
		
		int[] win = R.getBootstrapWinners(RND, 20);
		System.out.println("#**RELL.main win "+Arrays.toString(win));
		
		int[] support = R.getBootstrapSupport(win, new TreeComparator(rates.getTree()));
		System.out.println("#**RELL.main sup "+Arrays.toString(support));
		
		
		
	}
}
