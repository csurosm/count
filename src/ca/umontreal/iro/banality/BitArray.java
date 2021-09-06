/*
 * BitVector.java
 *
 * Created on June 8, 2005, 3:15 PM
 */

package ca.umontreal.iro.banality;

/**
 *
 * Simple implementation of a bit array that uses 1+o(1) bits per entry.
 * java.util.BitSet provides the same functionality but this is not implementation-dependent, 
 * and uses int values for indexing the positions. (Internally, a byte[] array is used that
 * is indexed with int values, and thus the bit array may be of length up to 8 times MAX_INT =2^31-1.)
 *
 * @author  csuros
 */
public class BitArray implements SetOfLongs {
    
    /** Creates a new instance of BitArray of given length: initialized to all-false bits */
    public BitArray(long length) {
        long max = 7L+8L*Integer.MAX_VALUE;
        if (length > max)
            throw new IllegalArgumentException("Bit array too large: wanted ["+length+"], max value "+(max));
        this.length=length;
        init();
    }
    
    private long length;
    
    /**
     * internal representation using bytes
     */
    private byte[] data;
    
    /**
     * Allocates the necessary space.
     */
    private void init(){
        
        int tru_length = (int)((length+7L)/8L);
        //System.out.println("#**BA.i allocate "+length+" tru "+tru_length);
        data = new byte[tru_length];
    }
    
    /**
     * @return the bit value at the given position
     */
    public final boolean get(long pos){
        int tru_pos = (int)(pos >>> 3);
        int offset = (int)(pos & 7L);
        return ((data[tru_pos]& bit_in[offset])!=0);
    }
    
    public boolean get(int pos){
        return get(toLong(pos));
    }
    
    /**
     * Sets the bit in the given position to the value
     */
    public void set(long pos, boolean value){
        if (value)
            set(pos);
        else
            clear(pos);
    }
    
    public void set(int pos, boolean value){
        set(toLong(pos),value);
    }
    
    private final long toLong(int pos){
        long retval=0L;
        if (pos<0) retval= (1L<<32)+(long)pos;
        else retval= (long)pos;
        //if (retval<0)
        //    System.err.println("#**BA.tL "+retval+" <- "+pos);
        return retval;
    }
    
    
    /** 
     * Sets the bit in the given position to true
     * @return whether it was set before
     */
    public final boolean getset(long pos){
        int tru_pos = (int)(pos >>> 3);
        int offset = (int)(pos & 7L);
        boolean is_there = ((data[tru_pos]& bit_in[offset])!=0);
        if (is_there)
            return true;
        else {
            data[tru_pos] |= bit_in[offset];
            return false;
        }
    }
    
    /** 
     *Sets the bit in the given position to true
     */
    public final void set(long pos){
        int tru_pos = (int)(pos >>> 3);
        int offset = (int)(pos & 7L);
        //if (offset<0 || tru_pos<0)
        //    System.err.println("#**BA.set pos "+pos+" offset "+offset+" tru "+tru_pos);
        data[tru_pos] |= bit_in[offset];
    }
    
    /** 
     * Sets the bit in the given position to true. 
     * Failsafe: negative ints are converted to positive longs.
     */
    public void set(int pos){
        set(toLong(pos));
    }
    
    /**
     * Sets the bit in the given position to false
     */
    public void clear(long pos){
        int tru_pos = (int)(pos >>> 3);
        int offset = (int)(pos & 7L);
        data[tru_pos] &= bit_out[offset];
    }
    
    public void clear(int pos){
        clear(toLong(pos));
    }
    
    /**
     * @return number of bits set to true in this array.
     */
    public int getCardinality(){
        int cnt=0;
        int max_pos = data.length;
        for (int tru_pos=0; tru_pos<max_pos; tru_pos++){
            cnt += getWeight(data[tru_pos]);
        }
        return cnt;
    }
    
    /**
     * Calculates the number of bits set in a byte
     */
    public static final int getWeight(byte b){
        if (b>=0)
            return bit_weight[b];
        else 
            return bit_weight[256+b];
    }
    
    private static final int[] bit_in = {1,2,4,8,16,32,64,128};
    private static final int[] bit_out = {254, 253, 251, 247, 239, 223, 191, 127};
    private static final int[] bit_weight 
        = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
            1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
            1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
            2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
            1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
            2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
            2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
            3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
    };
    
    
    /**
     * @return  number of bits that can be used 
     */
    public long size() {
        return length;
    }
    
    /**
     * Tests
     */
    public static void main(String[] args){
        int num_arrays = 10;
        BitArray[] A=new BitArray[num_arrays];
        byte[][] zed=new byte[num_arrays][];
        for (int i=0; i<num_arrays; i++){
            System.out.println("Allocate "+i+"/"+num_arrays+" zlen "+zed.length+" "+zed);
            A[i]=new BitArray(1<<24L);
            //zed[i] = new byte[1<<21];
        }
        java.util.Random RND = new java.util.Random();
        for (int j=0; j<10; j++){
            int arr_idx = RND.nextInt(num_arrays);
            long val = RND.nextInt() >>> 8; 
            System.out.println(arr_idx+"/"+val+" "+A[arr_idx].getset(val));
        }
    }
}