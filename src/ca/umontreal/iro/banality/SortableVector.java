/*
 * SortableVector.java
 *
 * Created on August 3, 2003, 2:04 AM
 */

package ca.umontreal.iro.banality;

/**
 * An extension of the Vector class with the possibility of sorting the
 * Vector elements in place.
 *
 * @author  csuros
 */
public class SortableVector extends java.util.Vector {
    
    public SortableVector() {super();}
    public SortableVector(java.util.Collection C){super(C);}
    public SortableVector(int initialCapacity){super(initialCapacity);}
    public SortableVector(int initialCapacity, int capacityIncrement){super(initialCapacity, capacityIncrement);}
    
    /**
     * Sorts the elements in place using java.util.Arrays.sort().
     */
    public synchronized void sort(){
        java.util.Arrays.sort(elementData,0,elementCount);
    }
    
    /**
     * Sorts the elements in place using java.util.Arrays.sort(Comparator).
     */
    public synchronized void sort(java.util.Comparator C){
        java.util.Arrays.sort(elementData,0,elementCount,C);
    }
    
    /**
     * Maximum number of elements listed by toString()
     */
    public static int MAX_ELEMENTS_LISTED=2;
    
    
    /**
     * Lists at most MAX_ELEMENTS_LISTED elements (default is 2).
     */
    public String toString(){
        return toString(MAX_ELEMENTS_LISTED);
    }

    /**
     * Gives the class name, vector size, and list of the first elements in the vector.
     *
     * @param max_elements_listed lists at most the given number of elements
     */
    public String toString(int max_elements_listed){
        StringBuffer sb=new StringBuffer(getClass().getName());
        sb.append('[');
        sb.append(elementCount);
        sb.append(" {");
        for (int i=0; i<Math.min(elementCount, max_elements_listed); i++){
            if (i>0) sb.append(", ");            
            sb.append(elementData[i]);
        }
        if (elementCount>max_elements_listed)
            sb.append(",...");
        sb.append("}]");
        return sb.toString();
    }
    
}
