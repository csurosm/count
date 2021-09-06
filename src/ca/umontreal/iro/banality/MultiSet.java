/*
 * MultiSet.java
 *
 * Created on August 3, 2003, 6:03 PM
 */

package ca.umontreal.iro.banality;

import java.util.Collection;
import java.util.Hashtable;
import java.util.Iterator;

/**
 * A multiset implementation based on an underlying Hashtable. 
 * When a new elemenbt is added, it is used as a key in the table to a counter 
 * object, which is maintained as items are deleted. 
 *
 * @author  csuros
 */
public class MultiSet implements Collection, java.io.Serializable {

    private Hashtable elements;
    private int totalSize=0;
    
    private static final int DEFAULT_CAPACITY = 16;
    private static final float DEFAULT_LOAD_FACTOR = 0.75f;
    
    /** Constructs a new, empty multiset with default initial capacity (16) and load factor (0.75)
     */
    public MultiSet() {
        elements=new Hashtable(DEFAULT_CAPACITY, DEFAULT_LOAD_FACTOR);
    }
    
    /** Constructs a new, empty multiset with the given initial capacity, and the default load factor
     * (0.75).
     */
    public MultiSet(int initialCapacity){
        elements=new Hashtable(initialCapacity, DEFAULT_LOAD_FACTOR);
    }
    
    /** Constructs a new, empty multiset with the given initial capacity and load factor.0
     */
    public MultiSet(int initialCapacity, float loadFactor){
        elements=new Hashtable(initialCapacity, loadFactor);
    }

    /**
     * Constructs a new multiset with the given elements.
     */
    public MultiSet(Collection C){
        elements=new Hashtable();
        addAll(C);
    }

    /** 
     * @return how many times the given object appears in the set (0 if the set does not contain it)
     */
    public int getMultiplicity(Object obj){
        if (elements.containsKey(obj)){
            Counter c=(Counter)elements.get(obj);
            return c.getValue();            
        } else 
            return 0;
    }
    
    /**
     * Increases the multiplicity of the given object.
     * @return true if the object was not in the set before
     */
    public boolean add(Object obj) {
        totalSize++;
        if (elements.containsKey(obj)){
            Counter c=(Counter)elements.get(obj);
            c.increment();
            return false;
        } else {
            Counter c=new Counter();
            elements.put(obj,c);
            return true; 
        }
    }
    
    /**
     * Increases the multiplicity of the given object.
     * @param multiplicity how many times the objetc needs to be added. It should be non-negative; 
     *     if equals 0, then nothing changes but the return value may be useful for tests. 
     *
     * @return true if the object was not in the set before
     */
    public boolean add(Object obj, int multiplicity) {
        totalSize+=multiplicity;
        if (elements.containsKey(obj)){
            Counter c=(Counter)elements.get(obj);
            c.increment(multiplicity);
            return false;
        } else if (multiplicity!=0){
            Counter c=new Counter(multiplicity);
            elements.put(obj,c);
            return true; 
        } else 
            return true;
    }
    
    public boolean addAll(Collection collection) {
        Iterator I=collection.iterator();
        boolean changed=false;
        while (I.hasNext()){
            changed |= add(I.next());
        }
        return changed;
    }
    
    public void clear() {
        elements.clear();
        totalSize=0;
    }
    
    public boolean contains(Object obj) {
        return elements.containsKey(obj);
    }
    
    public boolean containsAll(java.util.Collection collection) {
        Iterator I=collection.iterator();
        while (I.hasNext()){
            if (!contains(I.next()))
                return false;
        }
        return true;
    }
    
    public boolean isEmpty() {
        return elements.isEmpty();
    }

    /**
     * The iterator goes through the elements: every element is listed only once,
     * regardless of its multiplicity
     */
    public Iterator iterator() {
        return elements.keySet().iterator();
    }

    /**
     * Decreases the given object's multiplicity.
     */
    public boolean remove(Object obj) {
        if (elements.containsKey(obj)){
            Counter c=(Counter) elements.get(obj);
            if (c.decrement())
                elements.remove(obj);
            totalSize--;
            return true;
        } else
            return false;
    }
    
    public boolean removeAll(Collection collection) {
        Iterator I=collection.iterator();
        boolean changed=false;
        while (I.hasNext())
            changed |= remove(I.next());
        return changed;            
    }

    /**
     * Removes all the appearances of the given object.
     */
    public boolean purge(Object obj){
        if (elements.contains(obj)){
            Counter c=(Counter) elements.get(obj);
            totalSize -= c.getValue();
            elements.remove(obj);
            return true;
        } else 
            return false;
        
    }
    
    /**
     * Not supported.
     */
    public boolean retainAll(java.util.Collection collection) {
        throw new UnsupportedOperationException();
    }
    
    /**
     * @return number of elements with at least one occurrence.
     */
    public int numElements() {
        return elements.size();
    }
    
    /**
     * @return total number of elements, respecting their multiplicity
     */ 
    public int size(){
        return totalSize;
    }
    /**
     * An array containing all the values in this set (every element is includee only 
     * once, regardless of multiplicity).
     */ 
    public Object[] toArray() {
        return elements.keySet().toArray();
    }
    
    public Object[] toArray(Object[] obj) {
        return elements.keySet().toArray(obj);
    }
    
    public Collection elementsWithMinMultiplicity(int min_multiplicity){
        return new MultiplicityThreshold(min_multiplicity);
    }
    
    private static class Counter implements java.io.Serializable {
        private int cnt;
        private Counter(){
            this(1);
        }
        
        private Counter(int starting_value){
            cnt=starting_value;
        }
        private void increment(){ cnt++;}
        /**
         * @return true if counter gets to 0
         */
        private boolean decrement(){ 
            if (cnt>=0) cnt--;
            return (cnt==0);
        }
        
        private void increment(int inc){cnt+=inc;}
        
        private int getValue(){return cnt;}
    }
       
    private class MultiplicityThreshold extends java.util.AbstractCollection implements java.io.Serializable {
        private int threshold;
        private int num_abundant_elements;
        private MultiplicityThreshold(int threshold){
            this.threshold=threshold;
            num_abundant_elements=0;
            Iterator I=elements.values().iterator();
            while (I.hasNext()){
                Counter c=(Counter) I.next();
                if (c.getValue()>=threshold)
                    num_abundant_elements++;
            }
        }
        
        public boolean contains(Object obj){
            return getMultiplicity(obj)>=threshold;
        }
        
        public Iterator iterator() {
            return new AbundantElements(threshold);
        }        
        
        public int size() {
            return num_abundant_elements;
        }
        
    }
    
    private class AbundantElements implements Iterator, java.io.Serializable {

        private int min_multiplicity;
        private Iterator I;
        private Object nextGuy;
        private AbundantElements(int min_multiplicity){
            this.min_multiplicity=min_multiplicity;
            I=iterator();
            turn();
        }
        
        public boolean hasNext() {
            return (nextGuy != null);
        }
        
        public Object next() {
            if (nextGuy==null)
                throw new java.util.NoSuchElementException();
            Object retval=nextGuy;
            turn();
            return retval;            
        }
        
        private void turn(){
            while (I.hasNext()){
                nextGuy = I.next();
                Counter cnt=(Counter)elements.get(nextGuy);
                if (cnt.getValue()>=min_multiplicity)
                    return;
            }
            nextGuy=null;
        }
        
        /**
         * Not supported.
         */
        public void remove() {
            throw new UnsupportedOperationException();
        }
        
    }
    
}


