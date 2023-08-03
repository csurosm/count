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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.IntConsumer;

/**
*
* Traversal algorithms for an indexed tree. 
* 
*/
public class TreeTraversal 
{
	public TreeTraversal(IndexedTree tree)
	{
		this.tree = tree;
	}
	
	private final IndexedTree tree;
	
	/**
	 * Postorder traversal for the tree.
	 * 
	 * @param tree
	 * @return array of node indices as visited in postorder (parents after children)
	 */
	public static int[] postOrder(IndexedTree tree)
	{
		return (new TreeTraversal(tree)).postOrder();
	}
	
	/**
	 * Preorder traversal for the tree. 
	 * 
	 * @param tree
	 * @return array of node indices as visited in preorder (parents before children)
	 */
	public static int[] preOrder(IndexedTree tree)
	{
		return (new TreeTraversal(tree)).preOrder();
	}
	
	public int[] postOrder()
	{
		int[] postOrder = new int[tree.getNumNodes()];
		if (postOrder.length>0)
		{
			int root = postOrder.length-1;
			int n = visit(root, 0, null, postOrder);
			assert (n==postOrder.length); // everybody was visited
		}
		return postOrder;
	}
	
	public int[] preOrder()
	{
		int[] preOrder = new int[tree.getNumNodes()];
		if (preOrder.length>0)
		{
			int root = preOrder.length-1;
			int n = visit(root, 0, preOrder, null);
			assert (n==preOrder.length); // everybody was visited
		}
		return preOrder;
	}
	
	/**
	 * Reverses the order in the input array.
	 * @param order
	 * @return same array, bit elements are now in reverse order
	 */
	public static int[] reverse(int[] order)
	{
		int i=0;
		int j=order.length-1;
		while (i<j)
		{
			int oj = order[j];
			order[j]=order[i];
			order[i]=oj;
			++i;
			--j;
		}
		return order;
	}
	
	public static int[] levelOrder(IndexedTree tree)
	{
		return (new TreeTraversal(tree)).levelOrder();
	}
	
	public  int[] levelOrder()
	{
		int[] visited_nodes = new int[tree.getNumNodes()]; // used as a FIFO queue for level-order visit
		if (visited_nodes.length==0) return visited_nodes; // empty tree? 
		int next_put=0;
		int next_get=0;
		visited_nodes[next_put++] = tree.getRoot();
		while (next_get<next_put && next_put<visited_nodes.length) // stop when all leaves are enqueued
		{
			int node = visited_nodes[next_get++]; // dequeue
			if (!tree.isLeaf(node))
			{
				for (int ci=0; ci<tree.getNumChildren(node); ci++)
					visited_nodes[next_put++] = tree.getChild(node, ci); // enqueue
			}
		}
		return visited_nodes;
	}
	
	
	public static int[] indexOrder(IndexedTree tree)
	{
		return countUp(tree.getNumNodes());
	}
	
	/**
	 * Increasing sequence of naturals 0,1,...,<var>n</var>-1. 
	 * 
	 * @param n number of elements 
	 * @return
	 */
	public static int[] countUp(int n)
	{
		int[] countUp = new int[n];
		while (n>0)
		{
			--n;
			countUp[n]=n;
		}
		return countUp;
	}
	
//	public static void reverse(int[] A)
//	{
//		int i=0, j=A.length-1;
//		while (i<j)
//		{
//			int a = A[i];
//			int b = A[j];
//			A[i] = b;
//			A[j] = a;
//		}
//	}
	
	public static int[] getSubtreeNodes(IndexedTree tree, int subtree_root)
	{
		List<Integer> subtree_nodes = new ArrayList<>();
		TreeTraversal T = new TreeTraversal(tree);
		T.traverse(subtree_root, node->subtree_nodes.add(node), null);
		int[] get = new int[subtree_nodes.size()];
		for (int j=0; j<get.length; j++)
			get[j] = subtree_nodes.get(j);
		return get;
	}
	
	
	/**
	 * Tree traversal by recursion 
	 * 
	 * @param node
	 * @param preVisit called in prefix order, with node index 
	 * @param postVisit called in postfix order, with node index 
	 */
	public void traverse(int node, IntConsumer preVisit, IntConsumer postVisit)
	{
		if (preVisit != null)
			preVisit.accept(node);
		
		int num_children = tree.getNumChildren(node); 
		for (int ci=0; ci<num_children; ci++)
		{
			int child = tree.getChild(node, ci);
			traverse(child, preVisit, postVisit);
		}
		
		if (postVisit!=null)
			postVisit.accept(node);
	}
	
	
	private int visit(int node, int visit_pos, int[] preOrder, int[] postOrder)
	{
		if (preOrder!=null)
		{
			preOrder[visit_pos] = node;
		}
		int num_children = tree.getNumChildren(node); 
		for (int ci=0; ci<num_children; ci++)
		{
			int child = tree.getChild(node, ci);
			visit_pos = visit(child, visit_pos, preOrder, postOrder);
		}
		if (postOrder != null)
		{
			postOrder[visit_pos] = node;
		}
		return visit_pos+1;
	}
	
	/**
	 * Adapter for any indexed tree to respect the leaves-first order. 
	 * Our indexing is based on postorder: leaves first, and children have lower indices 
	 * than the parents. 
	 * 
	 * @author csuros
	 *
	 */
	public static class WellIndexed implements IndexedTree
	{
		/**
		 * Builds a leaves-first postOrder index for a given tree, using 
		 * a postorder traversal. 
		 * 
		 * @param tree a tree with nodes indexed 0..<var>n</var>-1 (any order)
		 */
		public WellIndexed(IndexedTree tree)
		{
			int num_nodes = tree.getNumNodes();
			this.traversal = new TreeTraversal(tree);
			this.original_index = new int[num_nodes];
			this.our_index = new int[num_nodes];
			buildIndex();
		}
		/**
		 * Instantiated with the input tree. 
		 */
		private final TreeTraversal traversal;
		/**
		 * Mapping from our index to original tree index. 
		 */
		private final int[] original_index;
		/**
		 * Inverse of {@link #original_index}: mapping from original tree index to our index.
		 */
		private final int[] our_index;
		
		/**
		 * Builds the data structure.
		 */
		private void buildIndex()
		{
			int[] postOrder = traversal.postOrder();
			int leaf_idx = traversal.tree.getNumLeaves();
			int ancestor_idx = traversal.tree.getNumNodes();
			int n = postOrder.length;
			assert (n==ancestor_idx); // by design
			while (n>0) // reverse postorder
			{
				--n;
				int node = postOrder[n];
				int node_idx;
				if (traversal.tree.isLeaf(node))
					node_idx = --leaf_idx;
				else
					node_idx = --ancestor_idx;
				original_index[node_idx] = node;
				our_index[node] = node_idx;
			}
			assert (leaf_idx==0); // getNumLeaves was correct
			assert (ancestor_idx==0);  // getnumNodes was correct
		}
		
		@Override
		public int getNumNodes()
		{
			return our_index.length;
		}
		
		@Override
		public int getNumLeaves()
		{
			return traversal.tree.getNumLeaves();
		}

		@Override
		public int getParent(int node)
		{
			int i = original_index[node];
			int pi = traversal.tree.getParent(i);
			return pi<0?pi:our_index[pi];
		}
		
		@Override 
		public int getNumChildren(int node)
		{
			int i = original_index[node];
			return traversal.tree.getNumChildren(i);
		}
		
		@Override
		public String getName(int node)
		{
			int i = original_index[node];
			return traversal.tree.getName(i);
		}
		
		@Override
		public int getChild(int node, int cidx)
		{
			int i = original_index[node];
			int c = traversal.tree.getChild(i, cidx);
			return our_index[c];
		}
	}
	
	/**
	 * Array of subtree sizes --- number of nodes within the subtree --- 
	 * rooted at each node.
	 * 
	 * @param tree
	 * @return
	 */
    public static int[] getSubtreeSizes(IndexedTree tree)
    {
        final int num_nodes = tree.getNumNodes();
        final int[] size = new int[num_nodes];
        for (int node=0; node<num_nodes; node++)
        {
            if (tree.isLeaf(node))
                size[node]=1;
            else
            {
                int s = 0;
                final int nc = tree.getNumChildren(node);
                for (int ci=0; ci<nc; ci++)
                {
                    int child = tree.getChild(node, ci);
                    s += size[child];
                }
                size[node]=s+1;
            }
        }
        return size;
    }
    
    public static int depth(IndexedTree tree, int node)
    {
    	int d = 0;
    	while (!tree.isRoot(node))
    	{
    		node = tree.getParent(node);
    		d++;
    	}
    	return d;
    }
    
    public static int[] getPathToRoot(IndexedTree tree, int node)
    {
    	int[] path = new int[depth(tree, node)+1];
    	int i=0;
    	do
    	{
    		path[i] = node;
    		node = tree.getParent(node);
    		i++;
    	} while (i<path.length);
    	return path;
    }
    
    public static int[] getHeights(IndexedTree tree)
    {
    	final int num_nodes = tree.getNumNodes();
    	final int[] height = new int[num_nodes];
    
    	// stays 0 for leaves
    	for (int node = tree.getNumLeaves(); node<num_nodes; node++)
    	{
			int h = 0;
			for (int ci=0; ci<tree.getNumChildren(node); ci++)
			{
				int child = tree.getChild(node, ci);
				h = Integer.max(h, height[child]);
			}
			height[node] = 1+h;
    	}
    	return height;
    }
    
    public static int[] getDepths(IndexedTree tree)
    {
    	final int num_nodes = tree.getNumNodes();
    	final int[] depth = new int[num_nodes];
    
    	int node = num_nodes-1;
    	assert (tree.isRoot(node)); 
    	// depth[node]=0;
    	while (node>0)
    	{
    		--node;
    		depth[node] = 1+depth[tree.getParent(node)];
    	}
    	return depth;
    	
    }
    
	
}
