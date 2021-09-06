/*
 * Heap.java
 *
 * Created on November 17, 2007, 12:15 PM
 */

package ca.umontreal.iro.banality;

/**
 * A binary heap structure (min-heap).
 *
 * @author  csuros
 */

import java.util.Comparator;
import java.util.Vector;

public class Heap<ItemClass extends Comparable<? super ItemClass>>
{
    
    /** 
     * A new heap with the given comparator
     * 
     * @param C a Comparator to be used for the heap-ordering 
     */
    public Heap(Comparator <? super ItemClass> C) 
    {
        H = new Vector<ItemClass>();
        this.C = C;
        H.add(null);
    }
    
    /**
     * A new heap with a class-specific Comparator (using the compareTo method of the item class)
     * 
     */
    public Heap()
    {
        H = new Vector<ItemClass>();
        this.C = new Comparator<ItemClass>()
            {
                public int compare(ItemClass o1, ItemClass o2)
                {
                    return o1.compareTo(o2);
                }
            };
            H.add(null);
    }
    
    private Vector<ItemClass> H;
    private Comparator<? super ItemClass> C;

    /**
     * Adds a new element to the heap
     * 
     * @param O the new element
     */
    public void add(ItemClass O)
    {
        H.add(O);
        swim(O, H.size()-1);
        //System.out.print(listHeap("add: "+O));
    }
    
    private void swim(ItemClass v, int pos)
    {
        int p = (pos%2==0?pos/2:(pos-1)/2); // ceil of (pos-1)/2
        while (p!=0 && C.compare(H.get(p),v)>0)
        {
            H.set(pos,H.get(p));
            pos = p;
            p = (pos%2==0?pos/2:(pos-1)/2);
        }
        H.set(pos,v);
    }
    
    /**
     * Removes the smallest element from the heap. 
     *  
     * @return the smallest element within the heap 
     * @throws java.util.EmptyStackException if the heap is empty
     */
    public ItemClass deleteMin()
    {
        if (H.size()==1)
            throw new java.util.EmptyStackException(); // "Heap is empty --- deleteMin is not possible");
        ItemClass r = H.get(1);
        if (H.size()>2)
        {
            ItemClass v= H.get(H.size()-1);
            H.remove(H.size()-1);
            sink(v,1);
        } else 
            H.remove(1);
        //System.out.print(listHeap("deleteMin: "+r));
        return r;
    }

    /**
     * Number of elements on the heap.
     * @return heap size
     */
    public int size()
    {
        return H.size()-1;
    }

    /**
     * Verifies if the heap is empty.
     * 
     * @return true of there are no elements on the heap
     */
    public boolean isEmpty()
    {
        return H.size()==1;
    }
    
    private void sink(ItemClass v, int pos)
    {
        int c = minChild(pos);
        while (c!=0 && C.compare(H.get(c),v)<0)
        {
            H.set(pos, H.get(c));
            pos = c;
            c = minChild(pos);
        }
        H.set(pos,v);
    }
    
    private int minChild(int i)
    {
        int c1 = i*2;
        int c2 = i*2+1;
        if (c2<H.size())
        {
            if (C.compare(H.get(c1),H.get(c2))<0)
            {
                return c1;
            } else
                return c2;
        } else if (c1<H.size())
        {
            return c1;
        } else
            return 0;
    }
    
    private String listHeap(String reason)
    {
        StringBuffer sb = new StringBuffer();
        sb.append("| "+reason+"\t heap of size "+size());
        for (int i=1; i<=size(); i++)
        {
            sb.append("\n| "+i+"\t"+H.get(i));
        }
        sb.append("\n");
        return sb.toString();
    }
}
