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


package count.ds;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import javax.xml.transform.stream.StreamSource;

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
	
	public default int tableHashCode()
	{
		int h = 0;
		for (int fam=0; fam<getFamilyCount(); fam++)
		{
			int[] profile = getFamilyProfile(fam);
			int hf = Arrays.hashCode(profile);
			h = 17*h+hf;
		}	
		return h;
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
