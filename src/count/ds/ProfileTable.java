package count.ds;
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



import java.util.Arrays;
import java.util.Comparator;


/**
 * Public interface used by the algorithms:  
 * a table of family profiles. A profile is an array of integers; negative values denote 
 * ambiguity. They are the copy numbers for a set of taxa (leaves in the phylogeny). 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 */
public interface ProfileTable 
{
	/**
	 * Number of rows (=gene families) in the table. 
	 * 
	 * @return non-negative integer. 
	 */
	public abstract int getFamilyCount();
	
	/**
	 * Number of columns in the table. 
	 * 
	 * @return non-negative integer.
	 */
	public abstract int getTaxonCount();
	
	/**
	 * Profile (a row) for a given family.  
	 * @param family_idx one of the rows 0..{@link #getFamilyCount()}-1
	 */
	public abstract int[] getFamilyProfile(int family_idx);
	
	public default String getFamilyName(int family_idx)
	{
		return Integer.toString(family_idx);
	}
	
	public default String[] getTaxonNames()
	{
		String[] names = new String[getTaxonCount()];
		for (int t=0; t<names.length; t++)
			names[t] = Integer.toString(t);
		return names;
	}
	
	/**
	 * Calculates an array of maximum family size per node (for an ancestral node, 
	 *    family size is the sum of copy numbers at the leaves in its subtree).
	 * 
	 * @param tree tree over which the summation is carried out 
	 * @return an array with as many entries as the tree nodes 
	 */
	public default int[] getMaxFamilySizes(IndexedTree tree)
	{
		int num_nodes = tree.getNumNodes();
		// maximum copy numbers
		int[] max_m = new int[num_nodes];
		int[] sum_m = new int[num_nodes]; // reused for summing in subtrees
		for (int fidx=0; fidx<getFamilyCount(); fidx++)
		{
			int[] profile = getFamilyProfile(fidx);
			for (int node=0; node<num_nodes; node++)
			{
				int n;
				if (tree.isLeaf(node))
					n = Math.max(0,profile[node]); // copy number at this node
				else // ancestors
				{
					n = 0;
					int num_children = tree.getNumChildren(node);
					for (int ci=0; ci<num_children; ci++)
					{
						int child = tree.getChild(node, ci);
						n += sum_m[child];
					}
				}
				if (n>max_m[node]) max_m[node]=n;
				sum_m[node] = n;
			}
		}	
		return max_m;
	}
	
	/**
	 * Calculates an array of maximum copy number per node (for an ancestral node, 
	 *    max copy number is the maximum of copy numbers at the leaves in its subtree).
	 * 
	 * @param tree tree over which the summation is carried out 
	 * @return an array with as many entries as the tree nodes 
	 */
	public default int[] getMaxCopies(IndexedTree tree)
	{
		int num_nodes = tree.getNumNodes();
		int num_leaves = tree.getNumLeaves();
		// maximum copy numbers
		int[] max_m = new int[num_nodes];
		for (int fidx=0; fidx<getFamilyCount(); fidx++)
		{
			int[] profile = getFamilyProfile(fidx);
			
			for (int node=0; node<num_leaves; node++)
			{
				max_m[node]=Integer.max(profile[node], max_m[node]);
			}
		}
		for (int node=num_leaves; node<num_nodes; node++)
		{
			int num_children = tree.getNumChildren(node);
			for (int ci=0; ci<num_children; ci++)
			{
				int child = tree.getChild(node, ci);
				max_m[node] = Integer.max(max_m[child],max_m[node]);
			}
		}	
		return max_m;
	}
	
	
	/**
	 * Mean copy number at the leaves.
	 * 
	 * @param geometric whether quadratic or arithmetic mean across families
	 * @return
	 */
	public default double getMeanCopies(boolean quadratic)
	{
		double sum_avg=0.0;
		double sum_avgsq = 0.0;
		
		int nF = getFamilyCount();
		for (int f=0; f<nF; f++)
		{
			int[] profile = getFamilyProfile(f);
			double sum = 0.0;
			for (int u=0; u<profile.length; u++)
			{
				sum = sum+Integer.max(0, profile[u]);
			}
			double avg = sum/profile.length;
			sum_avg += avg;
			sum_avgsq += avg*avg;
		}
		double avg_avg = sum_avg/nF;
		double avg_avgsq = sum_avgsq/nF;
		
		//System.out.println("#**PT.gMC arimean "+avg_avg+"\tgeomean "+Math.sqrt(avg_avgsq));
		
		return quadratic?Math.sqrt(avg_avgsq):avg_avg;
	}
	
	public default double getMeanMaxCopies(boolean quadratic)
	{
		double sum_max=0.0;
		double sum_maxsq = 0.0;
		
		int nF = getFamilyCount();
		for (int f=0; f<nF; f++)
		{
			int m = maxCopies(f);
			sum_max += m;
			sum_maxsq += m*m;
		}
		double avg_max = sum_max / nF;
		double avg_maxsq = sum_maxsq / nF;
		
		return quadratic?Math.sqrt(avg_maxsq):avg_max;
	}
	
	public default double getMeanMeanMaxCopies(IndexedTree tree, boolean quadratic)
	{
		if (tree==null) return Double.NaN;
		int num_nodes = tree.getNumNodes();
		int num_leaves = tree.getNumLeaves();
		double[] node_avg_max = new double[num_nodes];
		double[] node_avg_maxsq = new double[num_nodes];
		int nF = getFamilyCount();
		double[] pmax = new double[num_nodes]; // reused
		for (int f=0; f<nF; f++)
		{
			Arrays.fill(pmax, 0);
			int[] profile = getFamilyProfile(f);
			int node = 0;
			while (node<num_leaves)
			{
				double m = pmax[node] = profile[node];
				node_avg_max[node]+=m/nF;
				node_avg_maxsq[node]+=m*m/nF;
				++node;
			}
			while (node<num_nodes)
			{
				double m = pmax[node] = 0.0;
				for (int ci=0; ci<tree.getNumChildren(node); ci++)
				{
					int child = tree.getChild(node, ci);
					double cmax = pmax[child];
					m = Double.max(m, cmax);
				}
				pmax[node] = m;
				node_avg_max[node]+=m/nF;
				node_avg_maxsq[node]+=m*m/nF;
				++node;
			}
		}		
		// per node average  
		double avg_max = 0.0;
		double avg_maxsq = 0.0;
		
		for (int node=0; node<num_nodes; node++)
		{
			avg_max += node_avg_max[node]/num_nodes;
			avg_maxsq += node_avg_maxsq[node]/num_nodes;
		}
		
		return quadratic?Math.sqrt(avg_maxsq):avg_max;
	}
	
	
	public default double getMeanSumCopies(IndexedTree tree, boolean quadratic) {
		if (tree==null) return Double.NaN;
		int num_nodes = tree.getNumNodes();
		int num_leaves = tree.getNumLeaves();
		double[] node_avg_sum = new double[num_nodes];
		double[] node_avg_sumsq = new double[num_nodes];
		int nF = getFamilyCount();
		double[] psum = new double[num_nodes]; // reused
		for (int f=0; f<nF; f++)
		{
			Arrays.fill(psum, 0.0);
			int[] profile = getFamilyProfile(f);
			int node = 0;
			while (node<num_leaves)
			{
				double m = psum[node] = Integer.max(0,profile[node]);
				//if (plus1) m+= 1.0;
				node_avg_sum[node]+=m/nF;
				node_avg_sumsq[node]+=m*m/nF;
				++node;
			}
			while (node<num_nodes)
			{
				double m = psum[node] = 0.0;
				for (int ci=0; ci<tree.getNumChildren(node); ci++)
				{
					int child = tree.getChild(node, ci);
					double csum = psum[child];
					m = m+csum;
				}
				psum[node] = m;
				//if (plus1) m+= 1.0;
				node_avg_sum[node]+=m/nF;
				node_avg_sumsq[node]+=m*m/nF;
				++node;
			}
			
		}
		// per node average  
		double avg_sum = 0.0;
		double avg_sumsq = 0.0;
		
		for (int node=0; node<num_nodes; node++)
		{
			avg_sum += node_avg_sum[node]/num_nodes;
			avg_sumsq += node_avg_sumsq[node]/num_nodes;
		}
		
		return quadratic?Math.sqrt(avg_sumsq):avg_sum;
		
	}
	
	public default int minCopies()
	{
		int min=Integer.MAX_VALUE;
		
		for (int fidx=0; fidx<getFamilyCount(); fidx++)
		{
			int[] profile = getFamilyProfile(fidx);
			int sum_copies = 0;
			for (int c: profile)
			{
				if (c>=0)
					sum_copies += c;
			}
			if (sum_copies==0) return 0;
			if (sum_copies<min) min=sum_copies;
		}
		return min;
	}
	/**
	 * Maximum copy count at leaves
	 * within a family.
	 * 
	 * @param family
	 * @return
	 */
	public default int maxCopies(int family)
	{
		int[] profile = getFamilyProfile(family);
		int max = -1;
		for (int n: profile) max = Integer.max(max, n);
		return max;
	}
	
	public default boolean isBinaryTable()
	{
		int nF = getFamilyCount();
		for (int f=0; f<nF; f++)
		{
			int[] profile = getFamilyProfile(f);
			for (int n: profile)
				if (n>1) return false;
		}
		return true;
	}
	
	/**
	 * Number of copies across the leaves.
	 * 
	 * @param family
	 * @return
	 */
	public default int getMemberCount(int family)
	{
		int[] profile = getFamilyProfile(family);
		return Arrays.stream(profile).map(n->Integer.max(0,n)).sum();
//		
//		int cnt = 0;
//		for (int i: profile)
//			cnt+=Integer.max(0,i);
//		return cnt;
	}
	
	public default int getLineageCount(int family)
	{
		int[] profile = getFamilyProfile(family);
		return Arrays.stream(profile).map(n->(n>0?1:0)).sum();
		
//		int cnt = 0;
//		for (int i: profile)
//			if (i>0) cnt++;
//		return cnt;
		
	}
	
	/**
	 * Hashcode computed over the families: with columns 
	 * in lexicographic order of taxon names, so that 
	 * different trees/orderings of the same taxa give the same hashcode.
	 * 	  
	 * @return an integer value 
	 */
	public default int tableHashCode()
	{
		String[] taxa = getTaxonNames();
		Integer[] order = new Integer[taxa.length];
		for (int i=0; i<order.length; i++) order[i]=i;
		Arrays.sort(order, new Comparator<Integer>() {
			@Override
			public int compare(Integer o1, Integer o2) {
				return taxa[o1].compareTo(taxa[o2]);
			}
		});
		
		int h = 0;
		for (int fam=0; fam<getFamilyCount(); fam++)
		{
			int[] profile = getFamilyProfile(fam);
			int hf = 0;
			for (int i=0; i<order.length; i++)
				hf = 31*hf + Integer.hashCode(profile[order[i]]);
			// int hf = Arrays.hashCode(profile);
			h = 17*h+hf;
		}	
		return h;
	}
	
	public default String tableStatistics(IndexedTree tree) {
		StringBuilder sb = new StringBuilder();
		int nfam = this.getFamilyCount();
		sb.append("nfam ").append(nfam);
		int num_nodes = tree.getNumNodes();
		sb.append(", ntaxa ").append(tree.getNumLeaves()).append("/").append(num_nodes);
		sb.append("; qavg ").append(getMeanCopies(true));
		double nqmean = getMeanSumCopies(tree, true);
		double namean = getMeanSumCopies(tree,false);
		double complexity = (nqmean*nqmean + namean)*nfam*num_nodes/1e6;
		
		sb.append(", qsumavg ").append(nqmean);
		
		sb.append("; time complexity ").append(complexity).append("M");
		return sb.toString();
	}
	
	
	/**
	 * A table with as many profiles as leaves: each profile has a single copy 
	 * at a single leaf. Stores only the instantiating tree; profiles 
	 * for {@link #getFamilyProfile(int)} are created on the fly. 
 
	 * @param tree
	 * @return
	 */
	public static ProfileTable singletonTable(IndexedTree tree)
	{
		ProfileTable singletons
		 = new ProfileTable()
		 {
			@Override 
			public int getFamilyCount()
			{
				return tree.getNumLeaves();
			}
			
			@Override
			public int getTaxonCount()
			{
				return tree.getNumLeaves();
			}
			
			@Override
			public int[] getFamilyProfile(int fam)
			{
				int[] profile = new int[getTaxonCount()];
				profile[fam]=1;
				return profile;
			}
			@Override
			public String[] getTaxonNames()
			{
				return tree.getLeafNames();
			}
			@Override
			public String getFamilyName(int fam)
			{
				return "SINGLE_"+tree.getName(fam);
			}
		 };		
		 return singletons;
	}
	
	/** 
	 * A table with a single profile : all-0. 
	 * @param tree
	 * @return
	 */
	public static ProfileTable emptyProfile(IndexedTree tree)
	{
		ProfileTable empty 
		 = new ProfileTable()
		 {
			@Override 
			public int getFamilyCount()
			{
				return 1;
			}
			
			@Override
			public int getTaxonCount()
			{
				return tree.getNumLeaves();
			}
			
			@Override
			public String[] getTaxonNames()
			{
				return tree.getLeafNames();
			}
			
			@Override
			public int[] getFamilyProfile(int fam)
			{
				int[] profile = new int[getTaxonCount()];
				return profile;
			}
			
			@Override
			public String getFamilyName(int fam)
			{
				return "EMPTY";
			}
		 };	
		 return empty;

		
	}
	
	/**
	 * An empty table with no profiles. 
	 * 
	 * @param tree
	 * @return
	 */
	public static ProfileTable emptyTable(IndexedTree tree)
	{
		ProfileTable empty = new ProfileTable()
		 {
			@Override 
			public int getFamilyCount()
			{
				return 0;
			}
			
			@Override
			public int getTaxonCount()
			{
				return tree.getNumLeaves();
			}
			@Override
			public String[] getTaxonNames()
			{
				return tree.getLeafNames();
			}
			
			@Override
			public int[] getFamilyProfile(int fam)
			{
				return null;
			}
		 };	
		 return empty;
	}
	
	
}
