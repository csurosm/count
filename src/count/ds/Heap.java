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

import java.io.PrintStream;

/**
 * A binary heap structure (min-heap).
 *
 * @author  csuros
 */

import java.util.AbstractCollection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * A class implementing (min-)priority queue and directory (search for element).  
 * 
 * @author csuros
 *
 * @param <T> element type
 */
public class Heap<T extends Comparable<? super T>> extends AbstractCollection<T>
{
	
    private List<Node> H;
    private final Comparator<? super T> C;
    private final Map<T, Node> items;
    
    /** 
     * A new heap with the given comparator
     * 
     * @param C a Comparator to be used for the min-heap ordering 
     */
    public Heap(Comparator <? super T> C) 
    {
        H = new ArrayList<>();
        this.items = new HashMap<T, Node>();
        this.C = C;

        H.add(new Node()); // 1-based heap
    }
    
    /**
     * A new heap with a class-specific Comparator (using the compareTo method of the item class)
     * 
     */
    public Heap()
    {
        H = new ArrayList<>();
        this.items = new HashMap();
        this.C = new Comparator<T>()
            {
                public int compare(T o1, T o2)
                {
                    return o1.compareTo(o2);
                }
            };
        H.add(new Node());
    }
    
    /**
	 * Adds a new element to the heap
	 * 
	 * @param O the new element
	 * @return false if element is already on heap (not added)
	 */
    @Override
    public boolean add(T O)
	{
    	if (O==null)
    		throw new NullPointerException();    		
    	if (items.containsKey(O))
    		return false;
		Node node = new Node(O);
		int i = H.size();
		H.add(node);
		int j=node.swim(i);

		items.put(O, node);
		return true;
	}
    
    /**
     * Finds the element in the heap that equals the argument. 
     * 
     * @param element
     * @return
     */
    public T get(T element)
    {
    	Node N = items.get(element);
    	if (N==null)
    		return null;
    	else
    		return N.payload;
    }
    
    @Override
    public boolean contains(Object o)
    {
    	return items.containsKey(o);
    }
    
    @Override
    public void clear()
    {
    	Node head = H.get(0);
    	H.clear();
    	H.add(head);
    	items.clear();
    }
    
    @Override
    public boolean isEmpty()
    {
    	return H.size()==1;
    }
    
    @Override
    public int size()
    {
    	assert (items.size()==H.size()-1);
    	return items.size();
    }
    
    /**
     * Deletes and returns an element with minimum value. 
     * 
     * @return
     */
    public T deleteLeast()
	{
    	if (isEmpty())
    		throw new java.util.EmptyStackException(); // "Heap is empty --- deleteMin is not possible");
    	Node head = H.get(1);
    	T min = head.payload;
    	items.remove(min);
    	int lasti = H.size()-1;
    	Node last = H.get(lasti);
    	H.remove(lasti);
    	if (last != head)
    		last.sink(1);
    	return min;
	}
    
    /**
     * Returns but does not delete element with minimum value. 
     * 
     * @return
     */
    public T peek()
    {
    	if (isEmpty())
    		return null;
    	return H.get(1).payload;
    }
    
//    /**
//     * 
//     * 
//     * @param heap_pos must be between 1 and {@link #size()} [inclusively]; 1 is tp (least element)
//     * @return null if heap_pos=0
//     */
//    public T peek(int heap_pos)
//    {
////    	if (heap_pos<1) // not 
////    	{
////    		throw new IndexOutOfBoundsException();
////    	}
//    	return H.get(heap_pos).payload;
//    }
    
    /**
     * Call if priority changes for an element already in the heap. 
     * 
     * @param value
     */
    public void updateOrder(T value)
    {
    	if (!items.containsKey(value))
    		return ; //-1;
    	Node node = items.get(value);
    	int i = node.heap_pos;
    	Node parent = H.get(i/2);
    	
    	if (node.compareTo(parent)<0)
    		i=node.swim(i);
    	else
    		i=node.sink(i);
//    	return i;
    }
    
    @Override
    public boolean remove(Object value)
    {
    	if (!items.containsKey(value))
    		return false;
    	Node node = items.get(value);
//    	{ // DEBUG
//    		System.out.println("#**H.remove "+value+"\tnode "+node);
//    	}
    	
    	
    	
    	items.remove(value);
    	int i = node.heap_pos;
    	
    	int lasti = H.size()-1;
    	Node last = H.get(lasti);
    	H.remove(lasti);
    	
    	if (i<lasti)
    	{
	    	if (node.compareTo(last)<0)
	    	{
	    		last.sink(i);
	    	} else
	    	{
	    		last.swim(i);
	    	}
    	}
    	return true;
    }
	
    /**
     * Implemented as a series of {@link #remove(Object)} 
     * operations instead of using the iterator's remove.
     * 
     * @param elements elements to be deleted 
     * @return true if the heap changed 
     * 
     */
    @Override
    public boolean removeAll(Collection<?> elements)
    {
    	boolean removeAll = false;
    	for (Object o: elements)
    	{
    		boolean contains = contains(o);
    		boolean remove = remove(o);
//    		{ // DEBUG
//	    		System.out.println("#**H.rALL "+o+"\tcontains "+contains+"/"+contains(o)+"\tremove "+remove);
//    		}
    		removeAll = removeAll || remove;
    	}
    	return removeAll;
    }
    
    
    
    private Node get(int i)
    {
    	if (i<H.size())
    		return H.get(i);
    	else 
    		return null;
    }

	private Node minChild(int i)
	{
		assert (i!=0);
		int c = 2*i;
		Node child = get(c);
		if (child != null)
		{
    		Node child2 = get(c+1);
			if (child.compareTo(child2)>0) // works also if child2==null
				child = child2;
		}
		return child;
	}
    
	/**
	 * Cell content for binary heap, tracking heap position. 
	 * 
	 * @author csuros
	 *
	 */
    private class Node implements Comparable<Node>
    {
    	int heap_pos;
    	final T payload;
    	
    	/**
    	 * Sentinel node for position 0: with payload=null, 
    	 * placed before everybody else. 
    	 */
    	Node()
    	{
    		this(null);
    		this.heap_pos=0;
    	}
    	
    	Node(T payload)
    	{
    		this.payload = payload;
    	}
    	
    	boolean isRoot()
    	{
    		return payload==null;
    	}
    	
    	@Override
    	public int hashCode()
    	{
    		if (payload==null) return this.hashCode();
    		else return payload.hashCode();
    	}
    	
    	@Override
    	public boolean equals(Object o)
    	{
    		if (o==null) return false;
    		if (o instanceof Heap.Node)
    		{
	    		Node other = (Node)o; // unchecked conversion but payload type is immaterial
	    		if (this.isRoot())
	    			return other.isRoot();
	    		else
	    			return payload.equals(other.payload);
    		} else 
    			return super.equals(o);
    	}
    	
    	@Override
    	public String toString()
    	{
    		StringBuilder sb = new StringBuilder("[");
    		sb.append(payload).append("@").append(heap_pos).append("]");
    		return sb.toString();
    	}
    	
    	private void placeAt(int i)
    	{
    		this.heap_pos = i;
    		H.set(i, this);
    	}
    	
    	@Override
    	public int compareTo(Node other)
    	{
    		if (this.isRoot() || other==null)
    			return -1;
    		else
    			return other.isRoot()?+1:C.compare(this.payload, other.payload);
    	}
    	
    	/**
    	 * Placement of a node at position i, and percolating 
    	 * it up towards the root.   
    	 * 
    	 * @param i
    	 * @return
    	 */
    	int swim(int i)
    	{
    		int p = i/2;
    		Node parent = get(p);
    		while (this.compareTo(parent)<0) //C.compare(H.get(p),v)>0)
    		{
    			parent.placeAt(i);
    			i = p;
    			p = i/2;
    			parent = get(p);
    		}
    		assert (i>0);
    		this.placeAt(i);
    		return i;
    	}    
    	int sink(int i)
    	{
    		Node child = minChild(i);
    		while (child!=null && child.compareTo(this)<0)
    		{
    			int ci = child.heap_pos;
    			child.placeAt(i);
    			i = ci;
    			child = minChild(i);
    		}
    		this.placeAt(i);
    		return i;
    	}
    }
    
    public static void main(String[] args)
    {
    	PrintStream out = System.out;
    	int n = 9;
    	if (args.length>0)
    		n = Integer.parseInt(args[0]);
    	// generate n random numbers
    	double[] x = new double[n];
    	long seed = 2022;
    	Random RND = new Random(seed);
    	for (int i=0; i<n; i++)
    	{
    		double r = RND.nextDouble();
    		x[i]=r;
    	}
    	// add them on a heap
    	Heap<Double> heap = new Heap<>();
    	for (int i=0; i<n; i++)
    	{
    		double r = x[i];
    		heap.add(r);
    		out.println("add\t"+i+"\t"+r);
    	}
    	
//    	out.println("# heap: "+heap.H.toString());
    	
    	int head = 3;
    	int i=0; 
    	while (i<head)
    	{
    		double r = heap.deleteLeast();
    		out.println("dmin\t"+i+"\t"+r);
    		i++;
    	}
//    	out.println("# heap: "+heap.H.toString());

    	int del = 2*head;
    	for (int j=0; j<del; j++)
    	{
    		double r = x[j];
    		boolean d = heap.remove(r);
    		out.println("del\t"+i+"\t"+r+"\t"+d);
    		if (d) i++;
    	}
//    	out.println("# heap: "+heap.H.toString());

    	while (!heap.isEmpty())
    	{
    		double r = heap.deleteLeast();
    		out.println("dmin\t"+i+"\t"+r);
    		i++;
    	}
    }

	@Override
	public Iterator<T> iterator() 
	{
		final Iterator<Node> node_iterator = H.iterator();
		node_iterator.next(); // skip head sentinel
		return new Iterator<T>()
				{

					@Override
					public boolean hasNext() 
					{
						return node_iterator.hasNext();
					}

					@Override
					public T next() 
					{
						Node next = node_iterator.next();
						return next==null?null:next.payload;
					}
					
					@Override
					public void remove()
					{
						throw new UnsupportedOperationException();
					}
					
				};
	}
}
