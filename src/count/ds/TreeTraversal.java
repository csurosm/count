package count.ds;

import java.util.Arrays;

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
	
	
	public static int[] indexOrder(IndexedTree tree)
	{
		return countUp(tree.getNumNodes());
	}
	
	/**
	 * Increasing sequence of naturals 0,1,...,</var>n</var>-1. 
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
		public int getChild(int node, int cidx)
		{
			int i = original_index[node];
			int c = traversal.tree.getChild(i, cidx);
			return our_index[c];
		}
	}
}
