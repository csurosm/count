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
	 * @param family_idx one of the rows 0..{@link #getFamilyCount()-1}
	 */
	public abstract int[] getFamilyProfile(int family_idx); 
	
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
	

}
