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
 * @param <U> node types in this tree
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
     * @
     * 
     */
    public default String toString(int node_idx)
    {
    	return Integer.toString(node_idx);
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
    
    
//    /**
//     * Node selection by index.
//     * 
//     * @param node_idx index of node in standard traversal {@link #getNodes}; assumed to be valid 
//     * 
//     * @return tree node with the given index 
//     */
//    public abstract U getNode(int node_idx);

    
//    /**
//     * Array of nodes in the index order. 
//     * 
//     * Node with index <var>i</var> occupies cell <var>i</var> in the returned array.
//     * The indexes correspond to a depth-first traversal (parents 
//     * after children) visiting leaves first. Leaves occupy lower indices, 
//     * and the root in the array's last cell. 
//     * 
//     * @return an array of nodes (not null, but length-0 array for empty tree).
//     */
//    public abstract U[] getNodes();
    
    
//    /**
//     * Array of leaves in the index order. 
//     * 
//     * @return an array of leaf nodes
//     */
//    public default U[] getLeaves()
//    {
//        U[] leaves = Arrays.copyOf(getNodes(), getNumLeaves());
//        return leaves;
//    }
    
//    /**
//     * Number of edges in the tree. 
//     * 
//     * @return Number of nodes {@#getNumNodes()} minus 1. 
//     */
//    public default int getNumEdges()
//    {
//        return getNumNodes()-1;
//    }
    

    
//    /**
//     * Tree root.
//     * 
//     * @return null if empty tree; otherwise a node without parent
//     */
//    public default U getRoot()
//    {
//        int n = getNumNodes();
//        return n==0?null:getNode(n-1);
//    }
    
    
}
