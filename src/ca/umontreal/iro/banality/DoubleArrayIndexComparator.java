/*
 * DoubleArrayIndexComparator.java
 *
 * Created on February 14, 2004, 12:04 PM
 */

package ca.umontreal.iro.banality;

import java.util.Arrays;

/**
 *
 * A class for sorting indexes into a double[] array.
 *
 * @author  csuros
 */
public class DoubleArrayIndexComparator implements java.util.Comparator {
        
    private double[] values;

    public DoubleArrayIndexComparator(double[] values){
        this.values=values;
    }

    /**
     * Permits the use of java.util.Arrays.sort(Object[], Comparator c) in conjunction with an Integer array,
     * in decreasing order of values[]
     *
     * @param obj1,obj2 are Number objects (.intValue() will be applied to them), represent 
     * indexes into the array of values with which this comparator was instantiated.
     */
    public int compare(Object obj1, Object obj2) {
        int i1=((Number)obj1).intValue();
        int i2=((Number)obj2).intValue();
        return -Double.compare(values[i1],values[i2]);
    }

    /**
     * @return an Integer array of size array_length, entry 0 is 0, entry 1 is 1 and so on
     */
    public static Integer[] indexes(int array_length){
        Integer[] I=new Integer[array_length];
        for (int i=0; i<array_length; i++)
            I[i]=new Integer(i);
        return I;
    }

    /**
     * @return an index array for the double[] we were initialized with
     */
    public Integer[] indexes(){
        return indexes(values.length);
    }
    
    public int[] decreasingOrder(){
        Integer[] I=indexes();
        java.util.Arrays.sort(I,this);
        int[] retval=new int[I.length];
        for (int i=0; i<I.length; i++)
            retval[i]=I[i].intValue();
        return retval;
    }
    
    
    public int[] increasingOrder(){
        Integer[] I=indexes();
        java.util.Arrays.sort(I,this);
        int n=values.length-1;
        int[] retval=new int[n+1];
        for (int i=0; i<=n; i++)
            retval[i]=I[n-i].intValue();
        return retval;
        
    }

}
    
