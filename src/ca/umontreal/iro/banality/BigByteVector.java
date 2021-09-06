/*
 * BigByteVector.java
 *
 * Created on June 9, 2005, 11:27 AM
 */

package ca.umontreal.iro.banality;

/**
 *
 * A class implementing the funcionality of the Vector class, for byte elements.
 * The memory is allocated in chunks of 8Mb, unless the initial capacity is less than half of that
 * (in which case that capacity is kept as long as possible, but at the first necessary capacity increase 
 * the behavior changes to allocation by chunks).
 * If the capacity needs to be increased, it will be increased by a multiple of the chunk size.
 * The implementation has the advantage that there is almost no extra memory allocated when the 
 * space needs to be increased. 
 *
 * @author  [Miklos] : csuros @@@ iro.umontreal.ca
 */
public class BigByteVector {
    
    /** Creates a new instance of BigByteVector */
    public BigByteVector() {
        this(BLOCK_SIZE);
    }
    
    public BigByteVector(int initialCapacity){
        init(initialCapacity);
    }
    
    protected static final int BLOCK_BITS=23; // 20=1Mb 21=2Mb 22=4Mb 23=8Mb 24=16Mb 25=32Mb 26=64Mb 27=128Mb 
    protected static final int BLOCK_SIZE=1<<BLOCK_BITS;
    protected static final int BLOCK_OFFSET=BLOCK_SIZE-1;// for masking lower bits
    
    protected byte[][] elementData;
    
    protected int capacity;
    protected int elementCount;
    
    /**
     * Allocates space in elementData, initializes capacity and elementCount variables.
     */
    private void init(int initialCapacity){
        if (initialCapacity < BLOCK_SIZE/2){
            elementData=new byte[1][initialCapacity];
            capacity = initialCapacity;
        } else {
            int nblocks = (initialCapacity >> BLOCK_BITS)+1;
            elementData = new byte[nblocks][BLOCK_SIZE];
            capacity = nblocks<<BLOCK_BITS;
        }
        elementCount=0;
    }
    
    public byte[] toArray(){
        return toArray(new byte[elementCount]);
    }
    
    /**
     * Returns an array containing all of the elements in this vector in the correct order.
     * If the vector fits in the specified array, it is returned therein. Otherwise, a new array is allocated 
     * with the size of this vector.
     * If the vector fits in the specified array with room to spare (i.e., the array has more elements than the vector), 
     * the element in the array immediately following the end of the vector is set to the null byte. 
     * This is useful in determining the length of the vector only if the caller knows that the vector does not
     * contain any null bytes.
     */
    public byte[] toArray(byte[] arr){
        if (arr.length<elementCount)
            return toArray(new byte[elementCount]);
        
        for (int ncopied=0; ncopied<elementCount; ){
            int hi=ncopied >> BLOCK_BITS;
            int in_this_block=Math.min(elementCount-ncopied, BLOCK_SIZE);
            System.arraycopy(elementData[hi],0,arr,ncopied,in_this_block);
            ncopied+=in_this_block;
        }
        return arr;
    }

    /**
     * Computes a slice of the vector and returns it as an array. Does not check whether
     * the indexes surpass the element count (for speed).
     *
     * @param pos starting index for the slice
     * @param length length of the slice
     * @return a byte[] array of given length for the elements pos..pos+length-1
     */
    public byte[] slice(int pos, int length){
        return slice(pos, length, new byte[length]);
    }
    
    /**
     * Computes a slice of the vector and returns it as an array. Does not check whether
     * the indexes surpass the element count (for speed). The elements of the slice are put 
     * in the array given as a parameter: possible elements beyond the length are not changed.
     *
     * @param pos starting index for the slice
     * @param length length of the slice
     * @param destination where the elements are copied (with indexes 0 to length-1)
     * @return destination
     */
    public byte[] slice(int pos, int length, byte[] destination){
        int offset = pos & BLOCK_OFFSET;
        int hi = pos >> BLOCK_BITS;
        for (int ncopied=0; ncopied<length; ){
            int in_this_block = Math.min(BLOCK_SIZE-offset,length-ncopied);
            System.arraycopy(elementData[hi], offset, destination, ncopied, in_this_block);
            offset=0;
            hi++;
            ncopied+=in_this_block;
        }
        return destination;
    }
    
    /**
     * Enlarges the underlying data structure when necessary.
     *
     * @param minCapacity minimum capacity after enlarging (it is rounded up to the nearest multiple of BLOCK_SIZE if necessary). 
     */  
    public void ensureCapacity(int minCapacity){
        if (minCapacity > capacity && capacity < BLOCK_SIZE/2){ // was a tiny vector
            byte[] newData=new byte[BLOCK_SIZE];
            System.arraycopy(elementData[0],0,newData,0,elementCount);
            elementData[0]=null;
            elementData[0]=newData;
            double d=((int)(0.5+capacity/((1<<19)+0.)))*0.1;
            //System.out.println("#**BBV.eC capacity increase "+d+"Mb -> "+(1<<(BLOCK_BITS-20))+"Mb ");
            capacity = BLOCK_SIZE;
        }
        if (minCapacity>capacity){
            int oblocks = elementData.length;
            int nblocks = (minCapacity>>BLOCK_BITS)+1;
            //System.out.println("#**BBV.eC capacity increase "+(oblocks<<(BLOCK_BITS-20))+"Mb -> "+(nblocks<<(BLOCK_BITS-20))+"Mb ");
            byte[][] newData = new byte[nblocks][];
            System.arraycopy(elementData, 0, newData, 0, oblocks);
            elementData = newData;
            for (int i=oblocks; i<nblocks; i++)
                elementData[i]=new byte[BLOCK_SIZE];
            capacity = nblocks<<BLOCK_BITS;
        }
    }
    
    /**
     * Access to individual bytes in the data structure. For the sake of 
     * speed, it does not verify that index is within the allowable range 
     * (less than size()): the return value is a 0 byte 
     * in that case.
     * 
     * @return the element at the given position
     */
    public byte get(int index){
        int hi = index>>BLOCK_BITS;
        int lo = index & BLOCK_OFFSET;
        
        return elementData[hi][lo];
    }
    
     /* Access to indivudal bytes in the data structure. For the sake of 
     * speed, it does not verify that index is within the allowable range 
     * (less than size()): the return value is a 0 byte 
     * in that case.
     * 
     * @return the old element at the given position
     */
    public byte set(int index, byte element){
        int hi = index>>BLOCK_BITS;
        int lo = index & BLOCK_OFFSET;
        byte was = elementData[hi][lo];
        elementData[hi][lo]=element;
        return was;
    }
    
    /**
     * @return number of elements addded so far
     */
    public int size(){return elementCount;}
    
    /**
     * Sets the size to 0. Does not release memory!
     */
    public void clear(){
        this.elementCount=0;
    }
    
    public boolean isEmpty(){
        return (elementCount==0);
    }
    
    /**
     * Adds an array of values to the end of the vector.
     */
    public void addAll(byte[] values){
        addAll(values, 0, values.length);
    }
    
    /**
     * Adds a slice of an array to the end of the vector.
     * @param offset index of first element in values which needs to be added
     * @param length number of elements to be added from values array
     */
    public void addAll(byte[] values, int offset, int length){
        ensureCapacity(elementCount+length);
        
        for (int to_go=length; to_go>0;){
            int lo=elementCount & BLOCK_OFFSET;
            int hi=elementCount >> BLOCK_BITS;
		
            int ncopied = Math.min(to_go, BLOCK_SIZE-lo);
            //System.out.println("#**BBV.aA hi "+hi+"/"+elementData.length+" lo "+lo+" ncop "+ncopied+"/"+BLOCK_SIZE);
            System.arraycopy(values, offset, elementData[hi], lo, ncopied);
            elementCount += ncopied;
            to_go -= ncopied;
            offset += ncopied;
        }
    }
    
    /**
     * Adds a new element at the end of the vector.
     */
    public void add(byte element){
        //System.out.println("#**BBV.a "+element+" ["+elementCount+"]");
        ensureCapacity(elementCount+1);
        elementCount++;
        set(elementCount-1,element);
    }

    
    private void test(String[] args) throws Exception {
        int n=1<<24;
        int dn=8200;
        for (int i=0; i<n+dn; i++)
            add((byte)(i&63));
        System.out.println("#** Size "+size());
        
        byte[] a=new byte[n+dn];
        for (int i=0; i<n+dn; i++)
            a[i]=(byte)(i&63);
        
        addAll(a, dn, n);
        
        for (int i=0; i<2*n+dn; i++)
            if ((i&63) != get(i))
                System.out.println("#** i "+i+" emt "+get(i)+" "+(i&63));
            else if (i % 100000==0)
                System.out.println("#** i"+i);
        
    }
    
    /**
     * Runs debug tests.
     */
    public static void main(String[] args) throws Exception {
        (new BigByteVector()).test(args);
    }
    
//    public final Cursor getCursor(int starting_position){
//        return new Cursor(starting_position);
//    }
//    
//    public final Cursor getCursor(){
//        return new Cursor(0);
//    }
//    
//    /**
//     * A class for rapidly iterating through the values of this vector. The caller needs to keep track of the 
//     * current position: nothing is checked explicitly.
//     */
//    public final class Cursor {
//        private int current_hi;
//        private int current_lo;
//        private byte E[];
//        
//        private Cursor(int starting_pos){
//            setPosition(starting_pos);
//        }
//        
//        /**
//         * Moves the cursor forward
//         */
//        public final void stepForward(){
//            current_lo++;
//            if (current_lo==BLOCK_SIZE){
//                current_hi++;
//                E = elementData[current_hi];
//                current_lo=0;
//            }
//        }
//        
//        /**
//         * Moves the cursor backward
//         */
//        public final void stepBackward(){
//            current_lo--;
//            if (current_lo==0){
//                current_hi--;
//                current_lo=BLOCK_OFFSET;
//                E=elementData[current_hi];
//            }
//        }
//        
//        /**
//         * @return value at the current position of the cursor
//         */
//        public final byte get(){
//            return E[current_lo];
//        }
//        
//        /**
//         * Sets the value at the cursor
//         */
//        public final void set(byte b){
//            E[current_lo]=b;
//        }
//        
//
//        /**
//         * Sets the position of the cursor.
//         */
//        public final void setPosition(int current_pos){
//            this.current_hi = current_pos>>BLOCK_BITS;
//            this.current_lo = current_pos & BLOCK_OFFSET;
//            this.E = elementData[current_hi];
//        }
//        
//        /**
//         * @return the current position of the cursor
//         */
//        public final int getPosition(){return (current_hi<<BLOCK_BITS)+current_lo;}
//    }
}
