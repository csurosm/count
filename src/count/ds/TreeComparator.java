package count.ds;
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

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Class implementing Day's linear-time 
 * algorithm for comparing evolutionary tree topologies. 
 * 
 *
 */
public class TreeComparator 
{
	public TreeComparator(IndexedTree tree)
	{
		this.ref_tree = tree;
		this.leaf_indices = new HashMap<>();
		this.ref_subtrees = new int[tree.getNumLeaves()][];
		initDataStructures();
	}
	/**
	 * Indexing for values in int[] array
	 */
	private static final int PAIRED_LEAF=0;
	private static final int SUBTREE_SIZE=1;
	private static final int SUBTREE_ROOT=2;

	private final int[][] ref_subtrees;
	private final Map<String, Integer> leaf_indices;
	
	private final IndexedTree ref_tree;
	
	public IndexedTree getReferenceTree() { return ref_tree;}
	
	
	private void initDataStructures()
	{
		for (int leaf=0; leaf<ref_tree.getNumNodes(); leaf++)
			leaf_indices.put(ref_tree.getName(leaf), leaf);	
		initSubtrees(ref_tree, ref_tree.getRoot(), ref_subtrees);
	}
	
	public NodeMap map(IndexedTree query)
	{
		return new NodeMap(query);
	}
	
	/**
	 * Whether the same unrooted topology. 
	 * 
	 * @param query
	 * @return
	 */
	
	public boolean sameTopology(IndexedTree query)
	{
		if (ref_tree.getNumNodes()!=query.getNumNodes())
			return false;
		NodeMap map = new NodeMap(query);
		for (int node: map.fromReference())
			if (node==-1) return false;
		return true;
	}
	
	public boolean sameRootedTopology(IndexedTree query)
	{
		if (ref_tree.getNumNodes()!=query.getNumNodes())
			return false;
		NodeMap map = new NodeMap(query);
		int[] ridx = map.toReference();
		for (int node: ridx)
			if (node==-1) return false;
		int qroot = query.getRoot();
		return ref_tree.getRoot() == ridx[qroot];
	}
	
	public static boolean sameTopology(boolean rooted, IndexedTree ref, IndexedTree query)
	{
		TreeComparator cmp = new TreeComparator(ref);
		boolean sameTopology = rooted?cmp.sameRootedTopology(query):cmp.sameTopology(query);
		return sameTopology;
	}

	public static boolean sameTopology(IndexedTree ref, IndexedTree query)
	{
		TreeComparator cmp = new TreeComparator(ref);
		return cmp.sameRootedTopology(query);
	}
	
	
	/**
	 * Robinson-Foulds distance.
	 * 
	 * @param query
	 * @return
	 */
	public int getRFdistance(IndexedTree query)
	{
		NodeMap map = new NodeMap(query);
		return map.getRFdistance();
	}
	
	/**
	 * Robinson-Foulds distance between the two trees. 
	 * 
	 * @param ref
	 * @param query
	 * @return
	 */
	public static int getRFdistance(IndexedTree ref, IndexedTree query)
	{
		NodeMap map = (new TreeComparator(ref)).new NodeMap(query);
		return map.getRFdistance();
	}
	
	public class NodeMap
	{
		private NodeMap(IndexedTree query_tree)
		{
			this.query_tree = query_tree;
			this.query_subtrees = new int[ref_subtrees.length][];
			this.initQuerySubtrees();
		}
		private final IndexedTree query_tree;
		private final int[][] query_subtrees;
		private void initQuerySubtrees()
		{
			initSubtrees(query_tree, query_tree.getRoot(), query_subtrees);
		}
		
		/**
		 * Query node index for each reference node.
		 * 
		 * @return -1 entry for unmapped reference nodes
		 */
		public int[] fromReference()
		{
			int[] qry_map = new int[ref_tree.getNumNodes()];
			Arrays.fill(qry_map, -1);
			for (int qleaf=0; qleaf<query_tree.getNumLeaves(); qleaf++)
			{
				String qname = query_tree.getName(qleaf);
				if (leaf_indices.containsKey(qname))
					qry_map[leaf_indices.get(qname)]=qleaf;
			}
			for (int rleaf=0; rleaf<ref_tree.getNumLeaves(); rleaf++)
			{
				int other_leaf = ref_subtrees[rleaf][PAIRED_LEAF];
				int size = ref_subtrees[rleaf][SUBTREE_SIZE];
				if (query_subtrees[rleaf]==null)
				{
					// nobody maps there 
					assert (qry_map[rleaf]==-1);
				} else
				{
					int qry_other = query_subtrees[rleaf][PAIRED_LEAF];
					int rnode = ref_subtrees[rleaf][SUBTREE_ROOT];
					if (qry_other==other_leaf && query_subtrees[rleaf][SUBTREE_SIZE]==size)
					{
						int qnode = query_subtrees[rleaf][SUBTREE_ROOT];
						qry_map[rnode]=qnode;
					} else if (query_subtrees[other_leaf]!=null)
					{
						qry_other = query_subtrees[other_leaf][PAIRED_LEAF];
						if (qry_other==rleaf && query_subtrees[other_leaf][SUBTREE_SIZE]==size)
						{
							int qnode = query_subtrees[other_leaf][SUBTREE_ROOT];
							qry_map[rnode]=qnode;
						}
					}
				}
			}
			return qry_map;
		}
		
		/**
		 * Array of reference tree indices 
		 * for each query node. 
		 * 
		 * @return -1 entries for unmapped nodes in query
		 */
		public int[] toReference()
		{
			int[] ref_map = new int[query_tree.getNumNodes()];
			Arrays.fill(ref_map, -1);
			
			for(int qleaf = 0; qleaf<query_tree.getNumLeaves(); qleaf++)
			{
				String qname = query_tree.getName(qleaf);
				if (leaf_indices.containsKey(qname))
					ref_map[qleaf]=leaf_indices.get(qname);
				else
					ref_map[qleaf]=-1;
			}
			for (int rleaf=0; rleaf<ref_tree.getNumLeaves(); rleaf++)
			{
				if (query_subtrees[rleaf]==null)
				{
					// nobody maps here
				} else
				{
					int other_leaf = query_subtrees[rleaf][PAIRED_LEAF];
					int size = query_subtrees[rleaf][SUBTREE_SIZE];
					int qnode = query_subtrees[rleaf][SUBTREE_ROOT];
					if (ref_subtrees[rleaf][PAIRED_LEAF]==other_leaf
							&& (ref_subtrees[rleaf][SUBTREE_SIZE]==size))
					{
						int rnode = ref_subtrees[rleaf][SUBTREE_ROOT];
						ref_map[qnode]=rnode;	
					} else if (ref_subtrees[other_leaf][PAIRED_LEAF]==rleaf
							&& ref_subtrees[other_leaf][SUBTREE_SIZE]==size)
					{
						int rnode = ref_subtrees[other_leaf][SUBTREE_ROOT];
						ref_map[qnode]=rnode;		
					}
				}
			}
			return ref_map;
		}
		
		/**
		 * Robinson-Foulds distance
		 * 
		 * @return
		 */
		public int getRFdistance()
		{
			int rf_dist = 0;
			
			for (int node: fromReference())
				if (node==-1) rf_dist++;
			for (int node: toReference())
				if (node==-1) rf_dist++;
			return rf_dist;
		}
	}
	
	
	private int[] initSubtrees(IndexedTree tree, int node, int[][] subs)
	{
		assert subs.length == ref_tree.getNumLeaves();
		int[] min_max_n = new int[3];
		if (tree.isLeaf(node))
		{
			String taxon_name = tree.getName(node);
			assert (leaf_indices.containsKey(taxon_name));
			int leaf_idx = leaf_indices.get(taxon_name);
			min_max_n[0]=min_max_n[1]=leaf_idx;
			min_max_n[2]=1;
			subs[leaf_idx]=new int[3];
			subs[leaf_idx][PAIRED_LEAF]=leaf_idx;
			subs[leaf_idx][SUBTREE_SIZE]=1;
			subs[leaf_idx][SUBTREE_ROOT]=node;
		} else
		{
			// init min & max 
			min_max_n[0]=ref_tree.getNumLeaves();
			min_max_n[1]=-1;
			min_max_n[2]=0;
			for (int ci=0; ci<tree.getNumChildren(node); ci++)
			{
				int child =tree.getChild(node, ci);
				int[] cminmax = initSubtrees(tree,child,subs);
				min_max_n[0] = Integer.min(min_max_n[0],cminmax[0]);
				min_max_n[1] = Integer.max(min_max_n[1],cminmax[1]);
				min_max_n[2] += cminmax[2];
			}
			int min = min_max_n[0];
			int max = min_max_n[1];
			
			subs[min][PAIRED_LEAF]=max;
			subs[max][PAIRED_LEAF]=min;
			
			subs[min][SUBTREE_SIZE]=
			subs[max][SUBTREE_SIZE]=min_max_n[2];
			subs[max][SUBTREE_ROOT]=
			subs[min][SUBTREE_ROOT]=node;
		}
		return min_max_n;
	}
	
	
	public static void main(String[] args) throws Exception
	{
		String tree_file1 = args[0];
		String tree_file2 = args[1];
    	IndexedTree tree1 = count.io.NewickParser.readTree(new java.io.FileReader(tree_file1));
		TreeComparator C = new TreeComparator(tree1);
		for (IndexedTree tree2: count.io.NewickParser.readAllTrees(new java.io.FileReader(tree_file2)))
		{
			System.out.println("#TREE "+count.io.NewickParser.printTree(tree2));
//			IndexedTree tree2 = count.io.NewickParser.readTree(new java.io.FileReader(tree_file2));
			TreeComparator.NodeMap M = C.map(tree2);
			
			int[] ref_mapping = M.fromReference();
			int num_unmapped = 0;
			for (int node=0; node<ref_mapping.length; node++)
			{
				if (ref_mapping[node]==-1)
				{
					System.out.println("#MISSING\t"+node+"\t"+tree1.toString(node));
					num_unmapped++;
				}
			}
			System.out.println("Mapping from reference:\tunmapped\t"+num_unmapped+"\t"+Arrays.toString(ref_mapping));
			
			int[] og_nodes = M.toReference();
			int num_notref = 0;
			for (int j:og_nodes) num_notref += (j==-1?1:0);
			System.out.println("Mapping to reference: \tnotref\t"+num_notref+"\t"+Arrays.toString(og_nodes));
		}
		
		
	}
}
