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
import java.util.HashMap;
import java.util.Map;

import count.ds.Phylogeny.Node;

/**
 * Public interface used by the algorithms:  
 * a rooted ordered tree in which every node is assigned a unique index defined by a 
 * depth-first traversal. 
 *  
 * Node indices are consecutive 0,1,...,<var>n</var>-1, where <var>n</var> is the total number of nodes. 
 * Leaves are placed at the lower indices: 0,1,...,<var>m</var>-1 are 
 * leaves, and indices <var>m</var>, <var>m</var>+1,...,<var>n</var>-1 are ancestral nodes
 * with non-null children.  
 * Parent's index is always larger than the child's index; the root has index <var>n</var>-1.  
 * Children are ordered by the node's indexing {@link Node#getChild(int) }.
 * 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public interface IndexedTree
{
    /**
     * Number of nodes in the tree (including internal nodes and external nodes / leaves). 
     * 
     * @return a nonnegative integer 
     */
    public abstract int getNumNodes();
    
    /**
     * Number of leaves / terminals / external nodes in the tree. 
     * 
     * @return a nonnegative integer
     */
    public abstract int getNumLeaves();
    
    /**
     * Index of the root node. 
     * 
     * @return {@link #getNumNodes()}-1 by default.
     */
    public default int getRoot()
    {
    	return getNumNodes()-1;
    }

    /**
     * Index of parent node. 
     * 
     * Default implementation uses {@link Node#getParent() }.
     * 
     * @param node_idx node index; node may be root
     * @return parent's index; negative value for root
     */
    public abstract int getParent(int node_idx);
    
    
    /**
     * Description for node. 
     * 
     * @param node_idx node index
     * 
     */
    public default String toString(int node_idx)
    {
    	return Integer.toString(node_idx);
    }
    
    /**
     * Node name 
     * 
     * @param node
     * @return
     */
    public String getName(int node);
    
    public static String LEAF_IDENT_PREFIX = "T";
    public static String NODE_IDENT_PREFIX = "U";

    public default String getIdent(int node)
    {
    	String ident 
    	= (isLeaf(node)?LEAF_IDENT_PREFIX:NODE_IDENT_PREFIX)+Integer.toString(node);
    	String name = getName(node);
    	if (name != null)
    		ident = ident + "[" + name +"]";
    	return ident;
    }
    
    /**
     * Whether edge lengths are set. 
     * 
     * @return false by default
     */
    public default boolean hasLength()
    {
    	return false;
    }
    
    /**
     * Length of edge leading to node.
     * @param node
     * @return
     */
    public default double getLength(int node)
    {
    	if (isRoot(node))
    		return Double.POSITIVE_INFINITY;
    	else
    		return 1.0;
    }

    /**
     * Number of children. Default implementation calls {@link Node#getNumChildren() }.
     * 
     * @param node_idx node index
     * @return number of children; 0 for leaf/external node
     */
    public abstract int getNumChildren(int node_idx);
    
    /**
     * Index of a child node. 
     * 
     * @param node_idx node index;
     * @param cidx order of the child node (0=left, 1=right for binary); 0..{@link Node#getNumChildren()}-1
     * @return index of that child node; unspecified/exception for bad index 
     */
    public abstract int getChild(int node_idx, int cidx);

    /**
     * Test for being a terminal node / leaf. Default implementation calls {@link #getNumChildren(int) }.
     * 
     * @param node_idx node index
     * @return whether the indexed node is a leaf
     */
    public default boolean isLeaf(int node_idx)
    {
        return getNumChildren(node_idx)==0;
    }
    
    /**
     * Test for root. 
     * 
     * @param node_idx
     * @return
     */
    public default boolean isRoot(int node_idx)
    {
        return getParent(node_idx)<0;
    }
    
    /**
     * Array of leaf names (to be used for ordering the columns of the input table)
     * 
     * @return a non-null array of leaf names 
     */
    public default String[] getLeafNames()    
    {
    	int num_leaves = getNumLeaves();
        String[] leaf_names = new String[num_leaves];
        for (int leaf=0;leaf<num_leaves; leaf++)
        {
            leaf_names[leaf] = getName(leaf);
        }
        return leaf_names;
    }
    
    /**
     * Shortest positive-length edge in 
     * the phylogeny (not considering root's "length").  
     * 
     * @return {@link Double#POSITIVE_INFINITY} if all edges are 0, or a finite positive number
     */
    public default double shortestEdgeLength()
    {
    	double shortest = Double.POSITIVE_INFINITY;
    	for (int node=0; node<getNumNodes(); node++)
    	{
    		if (!isRoot(node))
    		{
    			double len = getLength(node);
    			if (len>0)
    				shortest = Double.min(shortest, len);
    		}
    	}
    	return shortest;
    }
    
    /**
     * Longest finite-length edge in 
     * the phylogeny (not considering root's "length").  
     * 
     * @return 0.0 if all edges are infinite or 0.0, or else a finite positive number
     */
    public default double longestEdgeLength()
    {
    	double longest = 0.0;
    	for (int node=0; node<getNumNodes(); node++)
    	{
    		if (!isRoot(node))
    		{
    			double len = getLength(node);
    			if (len != Double.POSITIVE_INFINITY)
    				longest = Double.max(longest, len);
    		}
    	}
    	return longest;
    }
    
    public default double[] edgeLengthQuantiles()
    {
    	double[] edges = new double[getNumNodes()];
    	int num_edges = 0;
    	for (int node=0; node<edges.length; node++)
    	{
    		if (!isRoot(node))
    		{
    			double len = getLength(node);
    			if (len != 0.0 && len != Double.POSITIVE_INFINITY)
    				edges[num_edges++] = len;
    		}
    	}
    	double[] quantiles = new double[3];
    	if (num_edges > 0)
    	{
	    	Arrays.sort(edges,0,num_edges);
	    	quantiles[1] = edges[edges.length/2];
	    	quantiles[0] = edges[edges.length/4];
	    	quantiles[2] = edges[3*edges.length/4];
    	}
    	return quantiles;
    }
    
}
