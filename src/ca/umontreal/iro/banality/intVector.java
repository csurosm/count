/*
 * IntegerVector.java
 *
 * Created on May 30, 2005, 8:50 PM
 */

package ca.umontreal.iro.banality;

/**
 * Vector-like class for storing int values.
 *
 * @author  csuros
 */
public class intVector {
    
    /**
     * Creates a new intVector of initial capacity 10.
     */
    public intVector() {
        this(10);
    }
    
    public intVector(int initialCapacity){
        init(initialCapacity);
    }

    /**
     * Underlying table for storing the elements.
     */
    protected int[] elementData;
    /**
     * Current size of this vector
     */
    protected int elementCount;
    /**
     * Current capacity of this vector.
     */
    protected int capacity;
    
    /**
     * Returns the array of elements used by this object; 
     * it may be larger than the vector size.
     * Typically, toArray() is a better choice.
     */  
    public int[] getElementData(){
        return elementData;
    }
    
    /**
     * Initializes (size 0) the vector with the given capacity.
     */
    private void init(int capacity){
        this.capacity=capacity;
        elementData= new int[capacity];
        elementCount=0;
    }
    
    /**
     *
     * Returns the component at the specified index.
     *
     * @param index an index into this vector.
     * @return the int at the specified index.
     * @throws ArrayIndexOutOfBoundsException if the index is negative or not less than the 
     *        current size of this intVector object is given.
     */
    public int get(int index){
        if (index<elementCount)
            return elementData[index];
        else 
            throw new ArrayIndexOutOfBoundsException("intVector of size "+elementCount+"; attempted access to element #"+index);        
    }
    
    /**
     * Sets the size to 0.
     */
    public void clear(){
        this.elementCount=0;
    }
    
    /**
     * Sets all elements to 0 (whether added or not added yet).
     */
    public void erase(){
        java.util.Arrays.fill(elementData, 0);
    }
    
    /**
     * Tests whether the vector has size 0.
     */
    public boolean isEmpty(){
        return (elementCount==0);
    }
    
    /**
     * Returns the size of the vector.
     */
    public int size(){
        return elementCount;
    }
    

    /**
     * Gives the current capacity of this vector.
     */
    public int getCapacity(){return capacity;}

    /**
     * Replaces the element at the specified position with the given value (the element must exist already)
     *
     * @param index index of element to replace
     * @param value new value at that position
     * @return the old value before replacement
     * @throws ArrayIndexOutOfBoundsException index out of range
     *
     */
    public int set(int index, int value){
        if (index<elementCount){
            int old=elementData[index];
            elementData[index]=value;
            return old;
        } else
            throw new ArrayIndexOutOfBoundsException("intVector of size "+elementCount+"; attempted access to element #"+index);        
    }
    
    /**
     * Enlarges the capacity if necessary.
     */
    public void ensureCapacity(int minCapacity){
        if (capacity<minCapacity){
            //System.out.println("#**BV.eC current "+capacity+" new "+minCapacity);
            int[] newElementData=new int[minCapacity];
            System.arraycopy(elementData, 0, newElementData, 0, elementCount);
            elementData=newElementData;
            capacity=minCapacity;
        }
    }
    
    /**
     * Increments the value at the given position: if it does not exist yet, it is added with value 1.
     */
    public void increment(int index){
        if (index>=elementCount)
            add(index,1);
        else 
            elementData[index]++;
    }

    /**
     * Adds a new element at the end of the vector.
     */
    public void add(int value){
        add(elementCount, value);
    }
    
    /**
     * Adds a new element at the given index. If the element exists already (ie index&lt;elementCount), 
     * then the behavior is 
     * the same as with set(), otherwise, the underlying array is enlarged if necessary to ensure 
     * that the element can be added.
     */
    public void add(int index, int value){
        if (index<elementCount)
            set(index,value);
        else{
            if (index >= capacity)
                ensureCapacity(Math.max(2*capacity,index+1));
            elementData[index]=value;
            elementCount=index+1;
        }
    }
    
    /**
     * Deletes the element at the given position.
     * Elements after the given position are shifted down by one.
     *
     * @return the deleted element
     */
    public int remove(int index){
        int retval = elementData[index];
        System.arraycopy(elementData, index+1, elementData, index, elementCount-index-1);
        elementCount--;
        return retval;
    }
    
    /**
     * Inserts the element at the given position.
     * Elements at and after the given position are shifted up by one.
     */
    public void insert(int index, int value){
        ensureCapacity(elementCount+1);
        System.arraycopy(elementData, index, elementData, index+1, elementCount-index);
        elementData[index]=value;
        elementCount++;
    }
    
    /**
     * Adds an array of values to the end of the vector.
     */
    public void addAll(int[] values){
        int n=values.length;
        ensureCapacity(elementCount+n);
        System.arraycopy(values, 0, elementData, elementCount, n);
    }
    
    /**
     * Adds an array of values at the given offset; vector is enlarged if necessary.
     */
    public void addAll(int[] values, int offset, int length){
        ensureCapacity(elementCount+length);
        System.arraycopy(values, offset, elementData, elementCount, length);
    }
    
    /**
     * @return an int[] array of size [exactly] equal to the number of elements.
     */
    public int[] toArray(){
        int[] retval = new int[elementCount];
        System.arraycopy(elementData, 0, retval, 0, elementCount);
        return retval;
    }
    
    /**
     * @return an int[] array for elements of index offset..elementCount-1.
     */
    public int[] toArray(int offset)
    {
        int len = elementCount-offset;
        return toArray(offset, len);
    }
    
    /**
     * @return an int[] array for elements of index offset..offset+length-1.
     */
    public int[] toArray(int offset, int length)
    {
        int[] retval = new int[length];
        System.arraycopy(elementData,offset,retval,0,length);
        return retval;
    }
       

    /**
     * Returns a comma-delimited list of elements
     */
    public String elementListString()
    {
        StringBuffer sb = new StringBuffer();
        for (int i=0; i<elementCount; i++)
        {
            if (i!=0)
                sb.append(",");
            sb.append(Integer.toString(get(i)));
        }
        return sb.toString();
    }
}
