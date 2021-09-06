/*
 * PartitionEnumerator.java
 *
 * Created on April 2, 2002, 11:46 PM
 */

package ca.umontreal.iro.combinatorics;
import java.util.Enumeration;
import java.util.NoSuchElementException;

/**
 *
 * creates an enumeration of ordered integer partitions with a fixed sum and length,
 * i.e. vectors of type x<sub>1</sub>,x<sub>2</sub>,...,x<sub>n</sub> with
 * 0 &lt;= x<sub>i</sub>&lt;= M, sum(x<sub>i</sub>)=N
 * 
 * n is argument <var>length</var> of the constructor
 * M is argument <var>maxValue</var> of the constructor
 * N is argument <var>sumValue</var> of the constructor  
 *
 * @author  miki
 */
public final class OrderedPartitionEnumeration implements IntArrayEnumeration {

    private int length;
    private int sumValue;
    private int maxValue;

    private int firstElement;
    private OrderedPartitionEnumeration otherElements;
    
    public OrderedPartitionEnumeration(int length, int sumValue, int maxValue) {
        this.sumValue=sumValue;
        this.maxValue=Math.min(maxValue,sumValue);
        this.length=length;
        reset();
    }
    
    private boolean setOtherElementsEnumeration(){
        this.otherElements
            =new OrderedPartitionEnumeration(
                    this.length-1,
                    this.sumValue-firstElement,
                    this.maxValue);
        return otherElements.hasMoreElements();
    }
    
    public boolean hasMoreElements(){
        return (firstElement <= maxValue);
    }

    public Object nextElement() throws NoSuchElementException {
        return this.nextElement(new int[length],0);
    }
        
    public Object nextElement(int[] result) throws NoSuchElementException {
        return this.nextElement(result, 0);
    }
    
    public Object nextElement(int[] result, int pos) throws NoSuchElementException {
        if (length > 1){
            otherElements.nextElement(result, pos+1);
        }
        result[pos]=firstElement;
        
        if (length==1 || !otherElements.hasMoreElements()){
            do {
                firstElement++;
            } while (firstElement<=maxValue 
                    && length>1 && !setOtherElementsEnumeration());
        }
        
        return result;
    }
    
    public void reset() {
        if (length==1)
            firstElement=sumValue; // so that nextElement gives fE=sV
        else
            firstElement=0;
        if (this.length>1)
            while (firstElement <= maxValue && !setOtherElementsEnumeration()) 
                firstElement++;
    }    
    
    
}

    
    

