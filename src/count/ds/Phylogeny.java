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

/**
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
        this(new Node());
    }
    /**
     * Creates a phylogeny with the given root.
     * 
     * @param root phylogeny root
     */
    public Phylogeny (Node  root)
    {
        this.root = root;
        this.nodes = null;
        this.has_length = false; // no lengths for now
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
        this.has_length = has_length;
    }
    
    public Node[] getNodes()
    {
        if (nodes==null) computeIndexes();
        return Arrays.copyOf(nodes, nodes.length);
    }
    
    public Node[] getLeaves()
    {
        if (nodes==null) computeIndexes();
        return Arrays.copyOf(nodes, getNumLeaves());
    }
    
    /**
     * Array of leaf names (to be used for ordering the columns of the input table)
     * 
     * @return a non-null array of leaf names {@link Node#getName()}
     */
    public String[] getLeafNames()
    {
        String[] leaf_names = new String[getNumLeaves()];
        for (Node L: getLeaves())
        {
            leaf_names[L.node_idx] = L.name;
        }
        return leaf_names;
    }
    
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
    public double getLength(int node_idx)
    {
    	return getNode(node_idx).getLength();
    }
    
    /**
     * Sets the {@link Node#node_idx} fields, and fills up the {@link #nodes} array. 
     */
    private void computeIndexes()
    {
        int num_leaves = getNumLeaves();
        int num_nodes = getNumNodes();
        nodes = new Node[num_nodes];
        getRootNode().computeIndexes(nodes, 0,num_leaves);
    }    

    /**
     * Standard implementation for a node in a phylogeny. 
     * 
     */
    public static class Node
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

        /**
         * Name that was set using {@link #setName(java.lang.String)}
         * @return 
         */
        public String getName() 
        {
            return name;
        }
        
        /**
         * Sets the node name.
         * 
         * @param s
         */
        public void setName(String s)
        {
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
            // ensure initialization and capacity
            if (children == null)
            {
                children = new Node[DEFAULT_NUM_CHILDREN];
            } else if (num_children==children.length)
            {
                int new_capacity=children.length+(children.length%3==0?children.length/3:children.length/2);
                // doubling every second time 0,2,3,4,6,8,12,16,...
                children = Arrays.copyOf(children, new_capacity);
            }
            num_children++;
            setChild(N, num_children-1);
        }        

        /**
         * Connects a child to this node.
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
            this.num_leaves_in_subtree = isLeaf()?1:0; // count this guy if it is a leaf
            this.num_nodes_in_subtree  = 1; // counting the node itself
            for (int ci=0; ci<getNumChildren(); ci++)
            {
                Node C = this.getChild(ci);
                this.num_leaves_in_subtree += C.num_leaves_in_subtree;
                this.num_nodes_in_subtree += C.num_nodes_in_subtree;
            }
            if (!isRoot())
                getParent().updateNodeCounts();
        }
        
        private void computeIndexes(Node[] nodes, int next_leaf_idx, int next_intl_idx)
        {            
            if (isLeaf())
            {
                nodes[next_leaf_idx]=this;
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
                this.node_idx=next_intl_idx;
            }
        }

        /**
         * A short identifying string for this node: leaf or ancestral + index.  
         * 
         * @return 
         */
        public String getNodeIdentifier()
        {
            return (isLeaf()?"L":"N")+Integer.toString(this.node_idx);
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
        
     }
    
}
