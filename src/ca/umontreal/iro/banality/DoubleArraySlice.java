/*
 * DoubleArraySlice.java
 *
 * Created on March 30, 2004, 9:01 AM
 */

package ca.umontreal.iro.banality;

/**
 * This class is useful in cases when a double[] array A is used,
 * its elements are entered consecutively and only a few recent elements 
 * are accessed after setting aany element. In particular, after using 
 * A[i] as an lvalue, only A[i], A[i-1], ..., A[i-m+1] can be used as rvalues.  
 *
 * @author  csuros
 */
public class DoubleArraySlice {
    
    /**
     * Sets the slice size (m).  
     */
    public DoubleArraySlice(int slice) {
        data = new double[slice];
        this.slice_size=slice;
        current_idx=-1;
        slice_idx=-1;
    }
    
    private int slice_size;
    private int current_idx;
    private double[] data;
    private int slice_idx;
    
    /*
     * @return index of the last element set
     */    
    public int getCurrentIndex(){
        return current_idx;
    }
    
    public void add(double x){
        slice_idx = (slice_idx+1) % slice_size;
        data[slice_idx]=x;
        current_idx++;
    }
    
    /**
     * Same as add(x) but checks if index is OK [it must be equal to getCurrentIndex()].
     */
    public void add(int idx, double x){
        if (idx != current_idx+1)
            throw new java.lang.ArrayIndexOutOfBoundsException("Can only set element ["+(current_idx+1)+"] and not ["+idx+"]");
        add(x);
    }
    
    /**
     * @param idx a value in the range c-m+1..c where c is the current index and m is the slice size
     * @return array element at the current index.
     */
    public double get(int idx){
        if (idx>current_idx || idx<=current_idx-slice_size || idx<0){
            throw new java.lang.ArrayIndexOutOfBoundsException("Accessible array indexes are "+(current_idx-slice_size+1)+".."+current_idx+" [and not "+idx+"]");
        }
        int i=(slice_idx+slice_size-(current_idx-idx)) % slice_size;
        return data[i];
    }
    
    /**
     * tests
     */
    public static void main(String[] args){
        DoubleArraySlice A=new DoubleArraySlice(5);
        for (int i=0; i<12; i++)
            A.add((double)i);
        System.out.println("A[10]="+A.get(10));
        A.add(12.0);
        System.out.println("A[11]="+A.get(11));
        A.add(13.0);
        System.out.println("A[13]="+A.get(13));
        System.out.println("A[3]="+A.get(3));
        
        
    }
    
}
