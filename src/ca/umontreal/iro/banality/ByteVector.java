/*
 * ByteVector.java
 *
 * Created on January 24, 2004, 8:57 PM
 */

package ca.umontreal.iro.banality;

/**
 * A class implementing the funcionality of the Vector class, for byte elements.
 * 
 * @author  csuros
 */
public class ByteVector  {
    
    /** Creates a new instance of ByteVector */
    public ByteVector() {
        this(10);
    }
    
    public ByteVector(int initialCapacity){
        init(initialCapacity);
    }
    
    protected byte[] elementData;
    protected int elementCount;
    protected int capacity;
    
    public byte[] getElementData(){
        return elementData;
    }
    
    private void init(int capacity){
        this.capacity=capacity;
        elementData= new byte[capacity];
        elementCount=0;
    }
    
    /**
     *
     * Returns the component at the specified index.
     *
     * @param index an index into this vector.
     * @return the byte at the specified index.
     * @throws ArrayIndexOutOfBoundsException if the index is negative or not less than the 
     * current size of this ByteVector object is given.
     */
    public byte get(int index){
        if (index<elementCount)
            return elementData[index];
        else 
            throw new ArrayIndexOutOfBoundsException("ByteVector of size "+elementCount+"; attempted access to element #"+index);        
    }
    
    /**
     * Sets the size to 0.
     */
    public void clear(){
        this.elementCount=0;
    }
    
    public boolean isEmpty(){
        return (elementCount==0);
    }
    
    public int size(){
        return elementCount;
    }
    
    public int getCapacity(){return capacity;}

    /**
     * Replaces the element at the specified position with the given value.
     *
     * @param index index of element to replace
     * @param value new value at that position
     * @return the old value before replacement
     * @throws ArrayIndexOutOfBoundsException index out of range
     *
     */
    public byte set(int index, byte value){
        if (index<elementCount){
            byte old=elementData[index];
            elementData[index]=value;
            return old;
        } else
            throw new ArrayIndexOutOfBoundsException("ByteVector of size "+elementCount+"; attempted access to element #"+index);        
    }
    
    public void ensureCapacity(int minCapacity){
        if (capacity<minCapacity){
            try {
            //System.out.println("#**BV.eC current "+capacity+" new "+minCapacity);
                byte[] newElementData=new byte[minCapacity];
                System.arraycopy(elementData, 0, newElementData, 0, elementCount);
                elementData=newElementData;
                capacity=minCapacity;
            } catch (java.lang.OutOfMemoryError E){
                // not that easy to extend
                System.out.println("#**BV.lM low memory ... need "+minCapacity+ " more bytes (currently using "+ capacity+")");
                throw new java.lang.OutOfMemoryError();
            }
        }
    }
    
    /**
     * Adds a new element at the end of the vector.
     */
    public void add(byte value){
        if (elementCount == capacity)
            ensureCapacity(2*capacity);
        elementData[elementCount]=value;
        elementCount++;
    }
    
    /**
     * Adds an array of values to the end of the vector
     */
    public void addAll(byte[] values){
        int n=values.length;
        ensureCapacity(elementCount+n);
        System.arraycopy(values, 0, elementData, elementCount, n);
        elementCount += n;
    }
    
    /**
     * Adds a slice of an array to the end of the vector.
     * @param offset index of first element in values which needs to be added
     * @param length number of elements to be added from values array
     */
    public void addAll(byte[] values, int offset, int length){
        ensureCapacity(elementCount+length);
        System.arraycopy(values, offset, elementData, elementCount, length);
        elementCount += length;
    }
    
    /**
     * Computes a String from the element data: starting with the i-th element, 
     * length elements are included in total. 
     * @param i starting position for substring returned (0-based)
     * @param length lenght of substring returned
     */
    public String getSubString(int i, int length){
        byte[] chars=new byte[length];
        System.arraycopy(elementData, i, chars,0, length);
        return new String(chars);
    }
    
    /**
     * @return a String built from the bytes in elementData
     */
    public String toString(){
        return getSubString(0,elementCount);
    }
    
}
