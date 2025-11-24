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

import static count.io.CommandLine.OPT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;

import count.io.CommandLine;
import count.io.GeneralizedFileReader;
import count.io.NewickParser;

/**
 * Standard implementation for rooted tree with arbitrary 
 * arities at the nodes.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class Phylogeny implements IndexedTree
{
    /**
     * Initializes a phylogeny with a non-null root.
     */
    public Phylogeny()
    {
        this.nodes = null;
        this.has_length = false; // no lengths for now
        this.root= new Node();
    }

    /**
     * Initializes a phylogeny with the same topology as the given tree.
     */
    public Phylogeny(IndexedTree tree)
    {
        this();
        copyFromTree(tree);
    }    
    
    
    public Phylogeny(IndexedTree tree, int subtree_root) {
    	this();
    	copyFromTree(tree, subtree_root);
    }
        
    private Node root;
    /**
     * True if edge lengths are meaningful.
     */
    private boolean has_length;
    
    /**
     * Inverted index for retrieving nodes by identifier. 
     */
    private Node[] nodes;
    
    public final Node getRootNode()
    {
        return root;
    }
    
    
    @Override 
    public int getNumLeaves()
    {
        if (root==null) return  0;
        else return root.num_leaves_in_subtree;
    }
    
    @Override
    public int getNumNodes()
    {
        if (root==null) return 0;
        else return root.num_nodes_in_subtree;
    }
    
    @Override 
    public String toString(int node)
    {
    	return nodes[node].toString();
    }
    
    
    
    /**
     * Whether this phylogeny has edge lengths, or just the topology.
     * 
     * @return 
     */
    @Override
    public boolean hasLength() 
    {
        return has_length;
    }
    
    /**
     * Affirms or denies that the tree edges have meaningful lengths. 
     * @param has_length whether this tree has edges with lengths.
     */
    public void hasLength(boolean has_length)
    { 
    	for (Node N: getNodes()) N.has_length = has_length;
        this.has_length = has_length;
    }
    
    public static Phylogeny starTree(String[] taxon_names)
    {
    	Phylogeny star = new Phylogeny();
    	for (int t=0; t<taxon_names.length; t++)
    	{
    		Node leaf = new Node(t);
    		leaf.setName(taxon_names[t]);
    		star.root.addChild(leaf);
    	}
    	star.computeIndexes();
    	return star;
    }
    
    public static Phylogeny randomTree(String[] taxon_names, java.util.Random RND, boolean yule_harding)
    {
    	if (RND==null) RND = new java.util.Random();
    	int num_leaves = taxon_names.length;
    	
    	int leaf = 0;
    	Phylogeny randomTree = new Phylogeny();
    	Node R = randomTree.getRootNode();
    	
    	List<Node> tree_nodes = new ArrayList<>();
    	if (!yule_harding)
    	{
	    	// do not add parent for Yule-Harding
	    	tree_nodes.add(R);
    	}
    	
    	while (leaf<num_leaves)
    	{
    		// random permutation of names
    		int j = leaf + RND.nextInt(num_leaves-leaf);
    		String leaf_name = taxon_names[j];
    		if (j!=leaf)
    		{
    			taxon_names[j]=taxon_names[leaf];
    			taxon_names[leaf] = leaf_name;
    		}
    		Node L = new Node(); // next leaf for random tree --- Node is static class 
    		L.setName(leaf_name);
    		
    		// add leaf on a random edge 
    		if (leaf==0 || leaf==1)
    		{
    			// first and second are added as root children
    			R.addChild(L);
    		} else
    		{
        		// find random placement in rooted tree
    			int edge = RND.nextInt(tree_nodes.size()); // 
    			R = tree_nodes.get(edge); // new sibling = regraft node
    			randomTree.regraft(L, R, false);
    			
    			// add parent if uniform tree distribution is wanted
    			// do not add parent for Yule-Harding (random binary search tree topology) distribution
    			if (!yule_harding)
    				tree_nodes.add(L.getParent());
    		}
			tree_nodes.add(L);
    		
    		leaf++;
    	}
    	randomTree.computeIndexes();
    	return randomTree;
    }
    
    private static Node[] EMPTY_ARRAY = new Node[0];
    public Node[] getNodes()
    {
    	List<Node> nodes = root.listNodes(null, N->N.isLeaf());
    	nodes = root.listNodes(nodes, N->!N.isLeaf());
    	return nodes.toArray(EMPTY_ARRAY);
//        if (nodes==null) computeIndexes();
//        return Arrays.copyOf(nodes, nodes.length);
    }
    
    public Node[] getLeaves()
    {
    	List<Node> leaves = root.listNodes(null,N->N.isLeaf());
    	return leaves.toArray(EMPTY_ARRAY);
//        if (nodes==null) computeIndexes();
//        return Arrays.copyOf(nodes, getNumLeaves());
    }
    
    /**
     * Array of leaf names (to be used for ordering the columns of the input table)
     * 
     * @return a non-null array of leaf names {@link Node#getName()}
     */
    @Override
    public String[] getLeafNames()
    {
        String[] leaf_names = new String[getNumLeaves()];
        int leaf_idx=0;
        for (Node L: getLeaves())
        {
            leaf_names[leaf_idx++] = L.name;
        }
        return leaf_names;
    }
    
    @Override
    public String getName(int node)
    {
    	return getNode(node).getName();
    }
    
//    @Override 
//    public String getIdent(int node)
//    {
//    	String name = getName(node);
//    	if (isLeaf(node) && false)
//    		return name;
//    	else 
//    	{
//    		String node_ident = IndexedTree.NODE_IDENT_PREFIX+Integer.toString(node);
//    		if (name != null)
//    			return node_ident + "[" + name+ "]";
//    		else
//    			return node_ident;
//    	}
//    }
    
    
    
    /**
     * @throws ArrayIndexOutOfBoundsException if called with bad index
     */
    public Node getNode(int node_idx)
    {
        if (nodes==null) computeIndexes();
        return nodes[node_idx];
    }
    
    /**
     * @throws NullPointerException if called with the root
     * @throws ArrayIndexOutOfBoundsException if called with bad index
     */
    @Override 
    public int getParent(int node_idx)
    {
    	Node parent = getNode(node_idx).getParent();
    	if (parent == null) 
    		return -1;
    	else 
    		return parent.getIndex();
    }
    
    /**
     * @throws ArrayIndexOutOfBoundsException if called with bad index
     */
    @Override
    public int getNumChildren(int node_idx)
    {
    	return getNode(node_idx).getNumChildren();
    }
    
    /**
     * @throws ArrayIndexOutOfBoundsException if called with bad index for node or the child
     */
    @Override 
    public int getChild(int node_idx, int child_idx)
    {
    	return getNode(node_idx).getChild(child_idx).getIndex();
    }
    
    /**
     * Length of edge leading to a node, or for root prior. 
     * 
     * @param node_idx
     * @return
     */
    @Override
    public double getLength(int node_idx)
    {
    	return getNode(node_idx).getLength();
    }
    /**
     * Reroots the phylogeny at an arbitrary node.
     * 
     * @param new_root
     */
    public void reroot(Node new_root)
    {
    	new_root.reroot();
    	this.root = new_root;
    	cleanUnaryNodes();
    	computeIndexes(); // need to redo node indexing  
    }
    
    public void rerootEdge(Node child)
    {
    	assert(!child.isRoot());
    	double len = child.getLength();

    	Node parent = child.getParent();
    	assert(!parent.isRoot() || parent.getNumChildren()>2);
    	parent.reroot();

    	child.detach();
    	Node new_root = new Node();
    	new_root.has_length = parent.has_length;
    	new_root.setLength(parent.getLength());
    	new_root.addChild(child);
    	new_root.addChild(parent);
    	this.root = new_root;
    	if (child.has_length)
    	{
    		assert (Double.isFinite(len) && len != 0.0);
    		child.setLength(len/2.0);
    		parent.has_length = true;
    		parent.setLength(len/2.0);
    	} else
    	{
    		parent.has_length = false;
    	}
    	assert (parent.getLength() != Double.POSITIVE_INFINITY);
    	assert (child.getLength() != Double.POSITIVE_INFINITY);
    	cleanUnaryNodes();
    	computeIndexes();
    }
    
    
//    /**
//     * Creates a new root with two children: the current root
//     * and the other phylogeny's root. 
//     * 
//     * @param other
//     */
//    public void join(Phylogeny other)
//    {
//    	assert (this.hasLength() == other.hasLength());
//    	
//    	double len = other.hasLength()?other.root.getLength():1.0;
//    	double rlen = root.getLength();
//    	
//    	Node join = new Node(2);
//    	join.setLength(rlen);
//    	root.setLength(len);
//    	join.addChild(root);
//    	join.addChild(other.root);
//    	this.root = join;
//    	computeIndexes();
//    }
    
    /**
     * Fuses a node into its parent. 
     * 
     * @param node
     */
    public void fuseIntoParent(Node node)
    {
    	node.fuseIntoParent();
    	computeIndexes(); // need to redo node indexing  
    }
    
    /**
     * Subtree prune and regraft.
     * 
     * @param pruned_node this is the node cut with its subtree
     * @param graft_position where it gets attached
     * @param graft_into whether node should become a new child, or rather a sibling (new node inserted on edge to graft position) 
     */
    public void pruneAndRegraft(Node pruned_node, Node graft_position, boolean graft_into)
    {
    	this.regraft(pruned_node, graft_position, graft_into);
    	cleanUnaryNodes();
    	computeIndexes();
    }
    
    
    /**
     * Reverses the order of leaves (changes node indexing). 
     */
    public void reverseLeafOrder() {
    	getRootNode().reverseChildOrdering();
    	this.computeIndexes();
    }
    
    /**
     * Removes the unselected leaves, and the resulting unary inner nodes.
     * 
     * @param leaves_kept
     */
    public void filterLeaves(String[] leaves_kept)
    {
    	int nl = getNumLeaves();
    	Map<String, Node> keep = new HashMap<>();
    	
    	Node none = new Node();
    	for (String s: leaves_kept)
    		keep.put(s, none);
    	
    	int ndel = 0;
    	int nkeep = 0;
    	
    	for (int leaf=0; leaf<nl; leaf++)
    	{
    		Node L = getNode(leaf);
    		String Lname = L.getName();
    		if (keep.containsKey(Lname))
    		{
//    			System.out.println("#**P.fL keep "+Lname+"\t"+L);
    			keep.put(Lname,  L);
    			++nkeep;
    		} else
    		{
    			Node P = L.getParent();
//    			System.out.println("#**P.fL delete "+Lname+"\t"+L+"\tparent "+P);
    			L.detach();
    			++ndel;
    			ndel += P.cleanUnaryNodes();
    			
    		}
    	}
    	
    	this.cleanUnaryNodes(); // removes single-child root if necessary
    	
    	int nmiss = 0;
    	for (String s: leaves_kept)
    	{
    		if (keep.get(s)==none)
    		{
    	    	System.out.println("#**P.fL miss "+s);
    			nmiss++;
    		}
    	}
    	System.out.println("#**P.fL kept "+nkeep+"\tmissed "+nmiss+"\tdelleaves "+(nl-nkeep)+"\tdelnodes "+ndel);
    	computeIndexes();
    }
    
    
    /**
     * Removes single-child nodes (including root if necessary).
     * @return number of nodes removed
     */
    private int cleanUnaryNodes()
    {
        int num_cleaned = root.cleanUnaryNodes();
        while (root.getNumChildren()==1)
        {
            Node child = root.getChild(0);
//            System.out.println("#**P.cUN root "+root+"\tchild "+child);

            child.setLength(root.getLength());
            child.has_length = root.has_length;
            root = child;
            root.detach();
            ++num_cleaned;
        }
        return num_cleaned;
    }
    
    /**
     * Subtree prune and regraft: node is detached and placed at a new position.
     * If grafted <em>into</em>, then the pruned node becomes a new child of the graft node. 
     * If not grafted into, then the graft node is replaced by a 
     * new node, with graft node on left and pruned node on right. 
     * In the latter case, if the graft node has length, then 
     * the new node is placed at halfway. 
     * The pruned node always retains its length.
     * 
     * @param pruned_node ths node gets detached and become s rightmost child at the new placement; or if no parent there is no detachment
     * @param graft_position where it gets attached
     * @param graft_into whether node should become a new child, or rather a sibling (new node inserted on edge to graft position) 
     */
    private void regraft(Node pruned_node, Node graft_position, boolean graft_into)
    {
    	if (!pruned_node.isRoot()) pruned_node.detach();
        Node p = graft_position;
        if (!graft_into)
        {
            p =  new Node();
            p.has_length = graft_position.hasLength();
            // need to create new parent
            if (graft_position.isRoot())
            {
                p.setLength(graft_position.getLength());
//                if (!p.hasLength())
//                    graft_position.setLength(1.0);
                p.addChild(graft_position);
                Phylogeny.this.root = p;
                if (graft_position.hasLength())
                {
                	assert pruned_node.hasLength();
                	graft_position.setLength(pruned_node.getLength());
                }
            } else
            {
                if (graft_position.hasLength())
                {
                    double len = graft_position.getLength()/2.0;
                    graft_position.setLength(len);
                    p.setLength(len);
                }
                int gidx = graft_position.getIndexAtParent();
                graft_position.getParent().setChild(p, gidx);
                p.addChild(graft_position);
            }
        }
        p.addChild(pruned_node);
    }
    

    private void copyFromTree(IndexedTree tree, final int subtree_root) {
    	TreeTraversal traversal = new TreeTraversal(tree);
    	int[] og_nodes = traversal.preOrder(subtree_root);
    	Node[] our_nodes = new Node[tree.getNumNodes()];

    	for (int node:og_nodes) {
    		Node N;
    		if (node==subtree_root) {
    			N = root; // created at instantiation
    		} else {
    			N = new Node(tree.getNumChildren(node));
    			int parent = tree.getParent(node);
    			Node P = our_nodes[parent];
    			P.addChild(N);
    	        if (tree.hasLength())
    	        	N.setLength(tree.getLength(node));
    		}
    		N.setName(tree.getName(node));
			our_nodes[node] = N;
    	}
    	this.hasLength(tree.hasLength());
    	computeIndexes(); // indexes are computed automatically but why wait 
    }    
    
    private void copyFromTree(IndexedTree tree)
    {
        int num_nodes = tree.getNumNodes();
        
        Node[] copied_nodes = new Node[num_nodes];
        int nidx=0;
        while (nidx<num_nodes) // with root
        {
            Node node;
            if (tree.isRoot(nidx))
            	node = root; // we were initialized with a single root
            else 
            	node = new Node(tree.getNumChildren(nidx));
            copied_nodes[nidx]= node;
            node.setName(tree.getName(nidx));
            node.node_idx = nidx;
            nidx++;
        }        
        nidx=0;
        if (tree.hasLength())
        {
            while (nidx<num_nodes)
            {
                copied_nodes[nidx].setLength(tree.getLength(nidx));
                nidx++;
            }
            nidx=0;
        }
        while (nidx<num_nodes)
        {
            Node node = copied_nodes[nidx];
            if (!tree.isLeaf(nidx))
            {
            	int nc = tree.getNumChildren(nidx);

                for (int ci=0; ci<nc; ci++)
                {
                    int cidx = tree.getChild(nidx, ci);
                    Node child = copied_nodes[cidx];
                    node.addChild(child);
                }
            }
            nidx++;
        }
    	this.hasLength(tree.hasLength());
    }
    
    
    /**
     * Sets the {@link Node#node_idx} fields, and fills up the {@link #nodes} array. 
     */
    protected void computeIndexes()
    {
        int num_leaves = getNumLeaves();
        int num_nodes = getNumNodes();
        nodes = new Node[num_nodes];
        getRootNode().computeIndexes(nodes, 0,num_leaves);
        for (Node N: nodes)
        	N.has_length = hasLength();
    }    
    
    
    /**
     * Replaces 0-length  edges with a 
     * tiny positive length; (does not fuse 0-length inner edges 
     * into parent).
     * 
     * @return number of edges modified (lengths changed from 0)
     */
    public int fixZeroEdges()
    {
    	double shortest = shortestEdgeLength();
    	double longest = longestEdgeLength();
    	double tiny = shortest*(shortest/longest); // tiny:shortest == shortest:longest
    	int fixed = root.fixZeroEdges(tiny);
    	int cleaned = cleanUnaryNodes();
//    	if (cleaned!=0)
//    	{
//    		System.out.println("#**P.fZE "+cleaned+" unary nodes cleaned");
//    	}
    	computeIndexes();
    	return fixed;
    }
    

    /**
     * Implementation for a node in a phylogeny. 
     * 
     */
    public static class Node 
    	implements Comparable<Node> // natural order by indexing
    {
        /**
         * Parent node.
         */
        private Node momma;
        /**
         * Number of children
         */
        private int num_children;
        /**
         * Dynamically allocated array of children.
         */
        private Node[] children;
        /**
         * Length of lineage leading to this node. 
         */
        private double length;
        /**
         * Name of this taxon; may be null
         */
        private String name;
        /**
         * Number of nodes in subtree; used for computing the indexes.
         */
        private int num_nodes_in_subtree;
        /**
         * Number of leaves in subtree; used for computing the indexes.
         */
        private int num_leaves_in_subtree;

        /**
         * Node index.
         */
        private int node_idx;
        
        private boolean has_length;
                
        /**
         * Instantiation of a new node with no children. 
         * 
         * Initial attributes are appropriate for root; length is set to 1.
         */
        private Node()
        {
            this.momma = null;
            this.children = null;
            this.num_children = 0;
            this.length = 1.0;
            this.num_nodes_in_subtree = 1;
            this.num_leaves_in_subtree = 1;
        }
        
        /**
         * Instantiation of a new node with capacity for children. 
         * 
         */
        private Node(int children_capacity)
        {
            this();
            if (children_capacity != 0)
                children = new Node[children_capacity];
        }


        /**
         * Parent access
         * 
         * @return parent node; null for root
         */
        public Node getParent() 
        {
            return momma;
        }
        
        /**
         * Number of children. 
         * 
         * @return number of children; 0 for leaf
         */
        public int getNumChildren() 
        {
            return num_children;
        }

        /**
         * Children from left (cidx=0) to right. 
         * (Behavior is unspecified for negative argument or too large of an index.)  
         * 
         * @param cidx child index: 0,1,...,{@link #getParent() }-1
         * @return child node with given index
         */
        public Node getChild(int cidx) 
        {
            // assert (idx>=0 && idx<num_children);
            return children[cidx];
        }

        /**
         * Length of edge leading to this node. 
         * 
         * @return 
         */
        public double getLength() 
        {
            return length;
        }

        public void setLength(double d) 
        {
            this.length = d;
        }
        
        public boolean hasLength()
        {
        	return has_length;
        }

        /**
         * Name that was set using {@link #setName(java.lang.String)}
         * @return 
         */
        public String getName() 
        {
//        	System.out.println("#**P.N.gN "+this.getNodeIdentifier()+"\tname "+name);
            return name;
        }
        
        
        /**
         * Sets the node name.
         * 
         * @param s
         */
        public void setName(String s)
        {
//        	System.out.println("#**P.N.sN "+this+"\tname "+name);
            this.name = s;
        }

        /**
         * Index of the node in the {@link IndexedTree} order.
         * 
         * @return
         */
        public int getIndex() 
        {
            return node_idx;
        }
        
        /**
         * Natural order by {@link #getIndex()}.
         */
        @Override 
        public int compareTo(Node other)
        {
        	return Integer.compare(node_idx, other.node_idx);
        }
        
        
        /**
         * Whether this node has no children.
         * 
         * @return true if {@link #getNumChildren() } is 0 
         */
        public boolean isLeaf()
        {
            return getNumChildren()==0;
        }
        
        /**
         * Whether this node has no parent. 
         * 
         * @return true if {@link #getParent() } is null.
         */
        public boolean isRoot()
        {
            return getParent()==null;
        }
        

        /**
         * 
         * @return {@link #getNumChildren()} if not found 
         */
        public int getIndexAtParent()
        {
            Node p = getParent();
            int idx=0; 
            while (idx<p.getNumChildren() && p.getChild(idx)!=this) idx++;
            return idx;
        }
        
        /**
         * Whether this node's parent b {@link #getParent()}
         * lists this node as a child.
         * 
         * @return
         */
        public boolean isDisownedChild()
        {
        	return !isRoot() && getIndexAtParent() == getParent().getNumChildren();
        }
        
        /**
         * Creates and adds new child for this node. 
         * 
         * @return the new child node 
         */
        public Node newChild() 
        {
            Node N = new Node();
            this.addChild(N);
            return N;
        }
        
        /**
         * Non-leaf nodes are binary by default.  
         */
        private static final int DEFAULT_NUM_CHILDREN=2;
        /**
         * Adds a new child; the child's parent is set to this node,
         * 
         * @param N new child to be added
         */
        public void addChild(Node N)
        {
            num_children++;
//            // ensure initialization and capacity
//            if (children == null)
//            {
//                children = new Node[DEFAULT_NUM_CHILDREN];
//            } else if (num_children==children.length)
//            {
//                int new_capacity=children.length+(children.length%3==0?children.length/3:children.length/2);
//                // doubling every second time 0,2,3,4,6,8,12,16,...
//                children = Arrays.copyOf(children, new_capacity);
//            }
        	ensureCapacity(num_children);
            setChild(N, num_children-1);
        }        
        
        
        
        
        /**
         * Detaches this node from its parent; updates node counts. 
         */
        private void detach()
        {
        	assert (!isRoot());
        	
        	Node P = getParent();
        	int idx = getIndexAtParent();
            for (int i=idx; i+1<P.num_children; i++)
            	P.children[i] = P.children[i+1];
            P.num_children --;
            P.children[P.num_children] = null; // no dangling pointers
            this.momma = null;
            P.updateNodeCounts();
        }
        
        
        /**
         * Allocates enough cells in the {@link #children} array. 
         * 
         * @param child_capacity
         */
        private void ensureCapacity(int child_capacity)
        {
            // ensure initialization and capacity
            if (children == null)
            {
                children = new Node[DEFAULT_NUM_CHILDREN];
            }
            int clen = children.length;
            if (clen < child_capacity)
            {
            	do 
            	{
                    clen=clen+(clen%3==0?clen/3:clen/2);
                    // doubling every second time 0,2,3,4,6,8,12,16,...
            	} while (clen < child_capacity);
            	children = Arrays.copyOf(children, clen);
            }
        }
        
        /**
         * Adds a child to this node and updates the node counts.
         * 
         * Sets the pointers between child and parent; 
         * sets child edge length to 1. Parent's number of children is unaffected. 
         * 
         * @param Nchild child  node
         * @param cidx at which index should this child be set
         * @return previous child node at this position (possibly null)
         */
        private Node setChild(Node Nchild, int cidx)
        {
            Node old_node = this.children[cidx];
            this.children[cidx]=Nchild;
            Nchild.momma=this;
            Nchild.length = 1.0;
            this.updateNodeCounts();
            return old_node;
        }

        /**
         * Updates the node counts at this node and all its ancestors. 
         */
        private void updateNodeCounts()
        {
        	updateNodeCounts(0, this);
        }
        
        private void updateNodeCounts(int rec_depth, Node N)
        {
        	
        	if (rec_depth>0 && N==this) // 
        		throw new IllegalStateException("Infinite loop because ancestors form a cycle of length "+rec_depth+" for "+N);
        	
            this.num_leaves_in_subtree = isLeaf()?1:0; // count this guy if it is a leaf
            this.num_nodes_in_subtree  = 1; // counting the node itself
            for (int ci=0; ci<getNumChildren(); ci++)
            {
                Node C = this.getChild(ci);
                this.num_leaves_in_subtree += C.num_leaves_in_subtree;
                this.num_nodes_in_subtree += C.num_nodes_in_subtree;
            }
            if (!isRoot())
                getParent().updateNodeCounts(rec_depth+1, N);
        	
        }
        
        /**
         * Recursive computation of node indexes.
         * 
         * @param nodes
         * @param next_leaf_idx
         * @param next_intl_idx
         */
        private void computeIndexes(Node[] nodes, int next_leaf_idx, int next_intl_idx)
        {            
            if (isLeaf())
            {
                nodes[next_leaf_idx]=this;
//                System.out.println("#**P.N.cI "+this.getNodeIdentifier()+"->"+next_leaf_idx);
                this.node_idx = next_leaf_idx;
            } else
            {
                for (int cidx=0; cidx<getNumChildren();cidx++)
                {
                    Node C = this.getChild(cidx);
                    C.computeIndexes(nodes, next_leaf_idx, next_intl_idx);
                    next_leaf_idx += C.num_leaves_in_subtree;
                    next_intl_idx += C.num_nodes_in_subtree-C.num_leaves_in_subtree;
                }
                nodes[next_intl_idx]=this;
//                System.out.println("#**P.N.cI "+this.getNodeIdentifier()+"->"+next_intl_idx);
                this.node_idx=next_intl_idx;
            }
        }
        
        /**
         * Lists the nodes in a postorder traversal of the subtree.
         * listNodes(null, null) gives all nodes in the subtree (including this node).
         * 
         * @param node_list may be null; if not, it is used for adding nodes in postorder
         * @param node_selector condition for including a node; if null, all nodes are selected
         * @return list of nodes (never null) in postorder traversal (same list as first argument, if not null)
         */
        public List<Node> listNodes(List<Node> node_list, Predicate<Node> node_selector)
        {
        	if (node_list == null)
        	{
        		node_list = new ArrayList<>();
        	}
        	for (int ci=0; ci<getNumChildren(); ci++)
        	{
        		Node child = getChild(ci);
        		node_list = child.listNodes(node_list, node_selector);
        	}
        	if (node_selector == null || node_selector.test(this))
        		node_list.add(this);
        	return node_list;
        }
        
        /**
         * Lists the ancestors up to the root (excluding this root)
         * 
         * @param node_list may be null; if not, it is used to add the noes 
         * @return list of nodes (never null) from parent up to the root
         */
        public List<Node> listAncestors(List<Node> node_list)
        {
        	if (node_list==null)
        	{
        		node_list = isRoot()?Collections.EMPTY_LIST:new ArrayList<>();
        	}
        	if (!isRoot())
        	{
            	node_list.add(momma);
        		node_list = momma.listAncestors(node_list);
        	}
        	return node_list;
        }
        
        /**
         * List of edges ending undirected paths of given distance range.
         * 
         * @param from a child or parent with recursion on other edges; or else (null or this) recursion on all edges 
         * @param distance nonnegative undirected distance, maximum (inclusive)
         * @param mindist nonnegative undirected distance, minimum (inclusive) 
         * @param neighbor_list must not be null;  neighbors are added on the given list 
         * @return neighbor_list 
         */
        private List<Node> neighbors(Node from, int distance, int mindist, List<Node> neighbor_list)
        {
//    		System.out.println("#**P.N.nghbr "+this+"\tfrom "+from+"\td "+distance+"\tmin "+mindist);
        	assert (neighbor_list != null);
        	
        	if (from == getParent()) // truly from parent, or starting search at root with from=null
        	{
            	if (mindist == 0)
            	{
            		//if (neighbor_list == null) neighbor_list = new ArrayList<>();
            		neighbor_list.add(this);
            	} else
            		mindist--;
            	if (distance>0)
            	{
	            	// continue downwards
	    			for (int ci=0; ci<getNumChildren(); ci++)
	    			{
	    				Node C = getChild(ci);
	    				neighbor_list = C.neighbors(this, distance-1, mindist, neighbor_list);
	    			}
            	}
        	} else // arriving from a child, or from == this
        	{
        		if (mindist <= 1)
        		{
        			assert (distance>=1);
            		//if (neighbor_list == null) neighbor_list = new ArrayList<>();
        			neighbor_list.add(this); // bottom of upward edge ending a path
        		}
        		if (mindist>0)
        			--mindist;
        		if (distance>0)
        		{
        			Node P = getParent();
        			if (distance>1 && P!=null)
        			{
        				// continue upwards
            			neighbor_list = P.neighbors(this, distance-1, mindist, neighbor_list);
        			}
        			// continue downwards but not where we came from
        			for (int ci=0; ci<getNumChildren(); ci++)
        			{
        				Node C = getChild(ci);
        				if (C!=from)
        				{
        					neighbor_list = C.neighbors(this, distance-1, mindist, neighbor_list);
        				}
        			}
        		}
        	}
        	//if (neighbor_list==null) neighbor_list = new ArrayList<>(); //Collections.EMPTY_LIST;
        	return neighbor_list;
        }
        
//        public List<Node> neighbors(int min_distance, int max_distance)
//        {
//        	return neighbors(this, max_distance, min_distance, null);
//        }
//        
        /**
         * List of nodes outside of subtree at a given minimum and maxmimum distance. 
         * 
         * @param min_distance (
         * @param max_distance
         * @return
         */
        public List<Node> graftNeighbors(int min_distance, int max_distance)
        {
        	List<Node> graft_neighbors;
        	if (isRoot())
        		graft_neighbors = Collections.EMPTY_LIST;
        	else
        		graft_neighbors = getParent().neighbors(this, max_distance-1, min_distance-1, new ArrayList<>());
        	
//        	System.out.println("#**P.N.gN "+this+"\tmin "+min_distance+"\tmax "+max_distance+"\tsize "+graft_neighbors.size());
        	return graft_neighbors;
        }
        
        
        /**
         * Reroots the tree by this node using recursion.
         */
        private void reroot()
        {
        	if (isRoot()) return; // already rooted here
            double len = getLength();
        	Node P = getParent();
        	P.reroot();
        	double rlen = P.getLength();
        	this.detach();
        	
        	this.addChild(P);
        	P.setLength(len);

        	this.setLength(rlen);
        }
        
        
        /**
         * Fuses the node into its parent (all its children become children of the original parent).
         */
        private void fuseIntoParent()
        {
        	if (isLeaf())
        		throw new IllegalArgumentException("Thou shalt not fuse a leaf into its parent.");
        	if (isRoot())
        		throw new IllegalArgumentException("Thou shalt not fuse the root into its (null) parent");

        	Node P = getParent();
        	int idx = getIndexAtParent();
        	
        	
        	int fused_nc = P.num_children-1+this.num_children;
        	P.ensureCapacity(fused_nc);
//        	System.out.println("#**P.N.fIP start "+this+"\tinto "+P+"\t// cap "+P.children.length); 
            { // shift old children to make place for the new
                int old_idx=P.num_children-1;
                int new_idx=fused_nc-1;
                
                while (old_idx>idx)
                {
                	Node child = P.children[old_idx];
//                	System.out.println("#**P.N.fIP shift "+child+"\t"+old_idx+"\t-> "+new_idx); 
                	P.children[new_idx] = child;
                    old_idx--;
                    new_idx--;
                }
            }
            double len = (hasLength()?getLength():0.0);
            for (int ci=num_children-1; ci>=0; ci--)
            {
                Node child = getChild(ci);
                child.setLength(child.getLength()+len);
                P.children[idx+ci]=child;
                child.momma = P;

//            	System.out.println("#**P.N.fIP add "+child+"\tat "+(idx+ci)); 
                
            }
            
            P.num_children = fused_nc;
//        	System.out.println("#**P.N.fIP children "+Arrays.toString(P.children)); 
            
            
            // update node statistics 
            do 
            {
            	P.num_nodes_in_subtree--;
//            	System.out.println("#**P.N.fIP sizeupdate "+P); 
            	P = P.getParent();
            } while (P!=null);
            
            // detach ourselves from the tree 
            this.momma = null;
            this.children = null;
            this.num_children = 0;
            this.num_leaves_in_subtree = this.num_nodes_in_subtree  = 1;
        }
        
        /**
         * Fuses nodes with a single child into their parents using postfix traversal of this subtree.
         * @return number of nodes cut
         */
        private int cleanUnaryNodes()
        {
            int nc = getNumChildren();
            int num_cleaned = 0;
            
            for (int j=0; j<nc; j++) // skipped when nc==0
                num_cleaned+=getChild(j).cleanUnaryNodes();

            if (nc==1 && !isRoot())
            {
                Node child = getChild(0);
                double ln = this.getLength();
                double lc = this.getLength();
                double len = (hasLength())?(ln+lc):1.0;
                Node parent = getParent();
                int idx = getIndexAtParent();
                // this.momma = null; // keep parent pointer in case we need to find ncestors
                parent.setChild(child, idx); // updates node counts at parent and all ancestors 
                child.setLength(len);
                ++num_cleaned;
            } 
            return num_cleaned;
        }
        
//        /**
//         * Subtree prune and regraft: this node is detached and placed at a new posiiton.
//         * 
//         * @param graft_position where it gets attached
//         * @param graft_into whether node should become a new child, or rather a sibling (new node inserted on edge to graft position) 
//         */
//        private void regraft(Node graft_position, boolean graft_into)
//        {
//        	detach();
//            Node p = graft_position;
//            if (!graft_into)
//            {
//                p =  new Node();
//                // need to create new parent
//                if (graft_position.isRoot())
//                {
//                    p.setLength(graft_position.getLength());
//                    if (!hasLength())
//                        graft_position.setLength(1.0);
//                    p.addChild(graft_position);
//                    Phylogeny.this.root = p;
//                } else
//                {
//                    if (hasLength())
//                    {
//                        double len = graft_position.getLength()/2.0;
//                        graft_position.setLength(len);
//                        p.setLength(len);
//                    }
//                    int gidx = graft_position.getIndexAtParent();
//                    graft_position.getParent().setChild(p, gidx);
//                    p.addChild(graft_position);
//                }
//            }
//            p.addChild(this);
//        }
        
        private int fixZeroEdges(double tiny_length_for_leaves)
        {
        	int fixed=0;
        	if (isLeaf())
        	{
        		double len = getLength();
        		if (len==0.0)
        		{
        			setLength(tiny_length_for_leaves);
//        			System.out.println("#**P.M.fZE leaf "+this);
        			++fixed;
        		}
        	} else
        	{
        		// save children bc they will get fused here, messing with the child indexing
        		int nc = getNumChildren();
        		Node[] children = new Node[nc];
        		for (int ci=0; ci<nc; ci++)
        			children[ci]=getChild(ci);
        		for (Node C: children)
        			fixed += C.fixZeroEdges(tiny_length_for_leaves);
        		if (!isRoot())
        		{
        			double len = getLength();
        			if (len==0.0)
        			{
//        				fuseIntoParent();
//            			System.out.println("#**P.M.fZE node "+this);
        				setLength(tiny_length_for_leaves);
        				fixed++;
        			}
        		}
        	}
        	return fixed;
        }
        
        
        private void reverseChildOrdering() {
        	int nc = getNumChildren();
        	for (int c=0; c<nc; c++) {
        		Node child = getChild(c);
        		child.reverseChildOrdering();
        	}
        	int i=0, j=nc-1; 
        	while (i<j) {
        		Node ci = children[i];
        		Node cj = children[j];
        		children[i]=cj;
        		children[j]=ci;
        		++i;
        		--j;
        	}
        }
        
        
        /**
         * A short identifying string for this node: leaf or ancestral + index.  
         * 
         * @return 
         */
        public String getNodeIdentifier()
        {
            return (isLeaf()?LEAF_IDENT_PREFIX:NODE_IDENT_PREFIX)+Integer.toString(this.node_idx);
        }
        
        /**
         * Node identifier + name (if set).
         * 
         * @return 
         */
        public String getFullName()
        {
            String nnn = getName();
            return getNodeIdentifier()+(nnn==null?"":"/"+nnn);
        }
        
        /**
         * Debug info on node.
         * 
         * @return info String with attribute values
         */ 
        protected String paramString()
        {
            StringBuilder sb=new StringBuilder();
            sb.append(getFullName());
            sb.append(" len ");
            sb.append(length);
            sb.append(" prnt ");
            if (getParent() != null)
                sb.append(getParent().getFullName());
            else
                sb.append('-');
            sb.append(" chld {");
            for (int i=0; i<num_children; i++)
            {
                if (i>0)
                    sb.append(", ");
                sb.append(getChild(i).getFullName());
            }
            sb.append("}");
            sb.append(", nl ").append(num_leaves_in_subtree);
            sb.append(", nn ").append(num_nodes_in_subtree);
            return sb.toString();
        }        
        
        /**
         * Debug info on node.
         * 
         * @return a String for log messages (class and {@link #paramString()})
         */
        @Override
        public String toString()
        {
            String simple_name = getClass().getSimpleName();
            String class_id_name = "";
            for (int i=0; i<simple_name.length(); i++)
            {
                char c = simple_name.charAt(i);
                if (Character.isUpperCase(c))
                    class_id_name = class_id_name + c;
            }
            if ("".equals(class_id_name))
                class_id_name = simple_name.substring(0,Math.min(simple_name.length(),5));
            StringBuilder sb=new StringBuilder(class_id_name);
            sb.append("[");
            sb.append(paramString());
            sb.append("]");
            return sb.toString();
        }
        
        public String prettySubtree()
        {
        	StringBuffer list_sb = this.prettySubtree(null, 0);
        	return list_sb.toString();
        }
        
        private StringBuffer prettySubtree(StringBuffer list_sb, int depth)
        {
        	
        	if (list_sb==null)
        		list_sb = new StringBuffer();
        	String space = "  ";
        	for (int tab=0; tab<depth; tab++)
        		list_sb.append(space);
        	list_sb.append(toString()).append("\n");
        	for (int ci=0; ci<getNumChildren(); ci++)
        	{
        		Node child = getChild(ci);
        		list_sb = child.prettySubtree(list_sb, depth+1);
        	}
        	return  list_sb;
        }


        
     }
    
    private int[] mapLCA(Phylogeny that)
    {
    	Map<String,Integer> leaf_indices=new HashMap<>();
    	for (int leaf=0; leaf<this.getNumLeaves(); leaf++)
    	{
    		leaf_indices.put(this.getName(leaf), leaf);
    	}
    	int[] mapLCA = new int[that.getNumNodes()];
    	for (int that_node=0; that_node<mapLCA.length; that_node++)
    	{
    		if (that.isLeaf(that_node))
    		{
    			Integer this_leaf = leaf_indices.get(that.getName(that_node));
    			if (this_leaf == null)
    				throw new IllegalArgumentException("Phylogeny argument must contain leaves with names known here");
    			else
    				mapLCA[that_node] = this_leaf;
    		} else
    		{
    			int this_node=-1;
    			for (int ci=0; ci<that.getNumChildren(that_node); ci++)
    			{
    				int that_child = that.getChild(that_node, ci);
    				int this_child = mapLCA[that_child];
    				if (ci==0) this_node = this_child;
    				else this_node = this.getLCA(this_node, this_child);
    			}
    			mapLCA[that_node] = this_node;
    		}
    	}
    	return mapLCA;
    }
    
    /**
     * Test code: outputs table of LCA mapping from a tree to downsampled phylogenies.
     * 
     * @param args list of tree files 
     * @throws Exception
     */
	public static void main(String[] args) throws Exception
	{
		Class<?> our_class = java.lang.invoke.MethodHandles.lookup().lookupClass();
	
		CommandLine cli = new CommandLine(args, our_class, 1);

        PrintStream out = System.out; 
    	String out_file = cli.getOptionValue(OPT_OUTPUT);
    	if (out_file!=null)
    	{
    		out = new PrintStream(out_file);
    	}	    
		out.println(CommandLine.getStandardHeader(our_class));
	    out.println(CommandLine.getStandardRuntimeInfo(our_class, args));
    
	    Map<String,Phylogeny> input_trees = new HashMap<>();
	    Phylogeny main_tree = cli.getTree();
	    assert (main_tree != null);
	    
	    if (main_tree == null)
	    	throw new IllegalArgumentException("Call with at least 1 tree");
	    
		for (int ti=0; ti<cli.getExtraArgumentCount(); ti++)
		{
			String treefile = cli.getExtraArgument(ti);
			Phylogeny phylo = NewickParser.readTree(GeneralizedFileReader.guessReaderForInput(treefile));
			input_trees.put(treefile,phylo);
		}
	    
		List<String> tree_order = new ArrayList<>(input_trees.keySet());
		Collections.sort(tree_order, new java.util.Comparator<>() {
			@Override
			public int compare(String o1, String o2) {
				Phylogeny t1 = input_trees.get(o1);
				Phylogeny t2 = input_trees.get(o2);
				return Integer.compare(t2.getNumNodes(), t1.getNumNodes()); // descending
			}
		});
		
		String main_treefile = cli.getTreeData().getFile().toString();
	    input_trees.put(main_treefile, main_tree);
	    tree_order.add(0, main_treefile);
	    
	    
		// header
		out.println("# Node mapping from first tree");
		for (int ti=0; ti<tree_order.size(); ti++)
		{
			String treefile = tree_order.get(ti);
			if (0<ti) out.print("\t");
			out.print(treefile);
		}
		out.println();
//		
//		final Phylogeny main_tree = input_trees.get(tree_order.get(0));
//		TreeComparator TC = new TreeComparator(main_tree);
		List<int[]> node_maps = new ArrayList<>(); // in tree order
		for (int ti=0; ti<tree_order.size(); ti++)
		{
			Phylogeny tree = input_trees.get(tree_order.get(ti));
			int[] toref = main_tree.mapLCA(tree);
			int[] invtoref = new int[main_tree.getNumNodes()];
			Arrays.fill(invtoref, -1);
			for (int node=0; node<tree.getNumNodes(); node++)
			{
				int rnode = toref[node];
				if (rnode!=-1)
					invtoref[rnode]=node;
			}
			
			node_maps.add(invtoref);
		}
		
		for(int u = 0; u<main_tree.getNumNodes(); u++)
		{
			for (int ti=0; ti<tree_order.size(); ti++)
			{
				int[] ref = node_maps.get(ti);
				if (0<ti) out.print("\t");
				out.print(ref[u]);
			}				
			out.println();
		}
	}
    
}
