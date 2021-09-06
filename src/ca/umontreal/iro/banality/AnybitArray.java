
package ca.umontreal.iro.banality;


import java.io.IOException;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.OutputStream;
import java.io.InputStream;
import java.io.ObjectOutputStream;
import java.io.ObjectInputStream;


/**
 * A class for packing b-bit integers into long arrays.
 * The structure is based on a 2D array, so that
 * many elements can be stored together
 * (more than the alllowed size for arrays in Java).
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class AnybitArray implements java.io.Serializable
{
    private static final long serialVersionUID = 20111025L;

    private static final int PAGE_SIZE = 1<<20; // this many long's need 8 Mbytes

    private static final boolean DEBUG_MESSAGES = true;

    private int[]  bit_offset;
    private long[] low_bits;
    private long[] high_bits;
    private int[]  num_low_bits;
    private long[] low_zeros;
    private long[] high_zeros;

    private int word_size;

    private long size;

    private long[][] pages;
    private int page_size;

    /**
     * Dummy instantiation: need to set word size and array size explicitly (in this order)
     */
    protected AnybitArray(){};

    /**
     * Instantiation of an array.
     *
     * @param size array length
     * @param word_size number of bits per entry.
     * @throws IllegalArgumentException if the total number of bits size*word_size is not divisible by 64 
     */
    public AnybitArray(long size, int word_size)
    {
        long num_bits = size*word_size;
        if ((num_bits & 63L) != 0L)
            throw new IllegalArgumentException("Size*word_size must be a multiple of 64 [size="+size+", word="+word_size+" requested]");
//        if (num_bits > (1L<<37))
//            throw new IllegalArgumentException("Maximum size is "+(1L<<37)/word_size+" for word size "+word_size+" [size="+size+"requested]");
        setWordSize(word_size);
        setArraySize(size);
    }


    /**
     * Sets page and word sizes; calculates bitmasks for shifting at long boundaries
     *
     * @param word_size bits per entry
     */
    protected void setWordSize(int word_size)
    {
        this.word_size = word_size;
        this.page_size = PAGE_SIZE *word_size;

        bit_offset = new int[64];
        low_bits = new long[64];
        high_bits = new long[64];
        num_low_bits = new int[64];
        low_zeros = new long[64];
        high_zeros = new long[64];
        long[] bits = new long[word_size+1];
        bits[0] = 0;
        for (int i=1; i<=word_size; i++)
            bits[i] = (bits[i-1] << 1) | 1L; // 2^i-1
        for (int i=0; i<64; i++)
            bit_offset[i] = (i*word_size) & 63; // mod 64
        long all_one = -1L;
        for (int i=0; i<64; i++)
        {
            int next_idx = (i+1) & 63;
            if (bit_offset[i]<bit_offset[next_idx]) // this word does not overlap boundaries
                num_low_bits[i] = word_size;
            else // this word does overlap boundaries
                num_low_bits[i] = 64-bit_offset[i];

            low_bits[i] = bits[num_low_bits[i]] << bit_offset[i]; // or bits[word_size] <<< bit_offset, as upper bits will "fall off" on the left
            high_bits[i] = bits[word_size-num_low_bits[i]];

            low_zeros[i] = all_one ^ low_bits[i];
            high_zeros[i] = all_one ^ high_bits[i];
        }
        if (DEBUG_MESSAGES)
        {
            for (int i=0; i<64; i++)
            {
                System.out.println("#*AA.sWS("+word_size+") "
                        +i
                        +"\toffset "+bit_offset[i]
                        +"\tnlo "+num_low_bits[i]
                        +"\tlo "+toBinaryString(low_bits[i])
                        +"\thi "+toBinaryString(high_bits[i])
                        +"\tlz "+toBinaryString(low_zeros[i])
                        +"\thz "+toBinaryString(high_zeros[i]));
            }
        }
    }

    /**
     * Allocates the pages for a given size and word length. Array entries are initialized to zero.
     * Must be called after {@link #setWordSize(int)}. 
     *
     * @param size array size (number of entries)
     * @throws IllegalArgumentException if the total number of bits size*word_size is not divisible by 64
     */
    protected void setArraySize(long size)
    {
        pages = null;
        this.size = size;
        long num_bits = size*word_size;
        if ((num_bits & 63L) != 0L)
            throw new IllegalArgumentException("Size*word_size must be a multiple of 64 [size="+size+", word="+word_size+" requested]");
        long data_size = num_bits >>> 6;
        long num_pages = data_size / page_size;
        if (num_pages * page_size < data_size)
            num_pages ++;

        if (DEBUG_MESSAGES)
        {
            double mem = (data_size * 8.0)/ (1<<20);
            System.out.println("#*AA.i allocating "+data_size+" longs\tmem "+mem+"M\tnum_bits "+num_bits+"\tsize "+size+"\tnpages "+num_pages+"\tpsize "+page_size);
        }
        if (num_pages > Integer.MAX_VALUE)
            throw new IllegalArgumentException("Cannot allocate so many pages ["+num_pages+"] requested");

        pages = new long[(int)num_pages][];
        long num_need_to_allocate = data_size;
        {
            int page_idx = 0;
            while (num_need_to_allocate>page_size)
            {
                pages[page_idx++] = new long[page_size];
                num_need_to_allocate -= page_size;
            }
            if (num_need_to_allocate>0L)
                pages[page_idx] = new long[(int)num_need_to_allocate];
        }
    }

    /**
     * Resets the array length.
     *
     * @param new_size new number of entries; must be larger than current size
     */
    protected void enlargeArray(long new_size)
    {
        assert (new_size>=size);
        long num_bits = new_size*word_size;
        if ((num_bits & 63L) != 0L)
            throw new IllegalArgumentException("Size*word_size must be a multiple of 64 [size="+new_size+", word="+word_size+" requested]");
        long data_size = num_bits >>> 6;
        long num_pages = (data_size+page_size-1) / page_size; // rounded to ceiling
//        if (num_pages * page_size < data_size)
//            num_pages ++;
        long[][] new_pages = new long[(int)num_pages][];

        int num_full_pages = pages.length-1;
        System.arraycopy(pages, 0, new_pages, 0, num_full_pages);
        long num_need_to_allocate = data_size - num_full_pages*page_size;
        assert (num_need_to_allocate>0L);
        {
            long[] extended_last_page = new long[Math.min((int)num_need_to_allocate, page_size)];
            assert (extended_last_page.length>=pages[num_full_pages].length);
            System.arraycopy(pages[num_full_pages], 0, extended_last_page, 0, pages[num_full_pages].length);
            num_need_to_allocate -= extended_last_page.length;
            new_pages[num_full_pages]=extended_last_page;
        }
        pages = new_pages;
        int page_idx = pages.length;
        while (num_need_to_allocate>page_size)
        {
            pages[page_idx++] = new long[page_size];
            num_need_to_allocate -= page_size;
        }
        if (num_need_to_allocate>0L)
            pages[page_idx] = new long[(int)num_need_to_allocate];
        size = new_size;
    }



    /**
     * Number of entries stored in this structure.
     *
     * @return array size
     */
    public long size(){ return this.size;}

    /**
     * Number of bits used per entry.
     *
     * @return word length
     */
    public int getWordSize(){ return this.word_size;}

    /**
     * Entry at a given index
     *
     * @param idx array position index (0..{@link #size()}-1)
     *
     * @return then value stored there; only the lower bits store information (number of bits given by {@link #getWordSize() }
     */
    public final long get(final long idx)
    {
        long full_idx = (idx*word_size) >>> 6;
        int page_idx = (int)(full_idx / page_size);
        long[] data = pages[page_idx];
        int data_idx = (int)(full_idx % page_size);
        long data_val=data[data_idx];

        int idx_mod = (int)(idx & 63);
        long lo = (data_val & low_bits[idx_mod]) >>> bit_offset[idx_mod]; //same as  >>> (idx*word_size); only lower 64 bits count

        if (num_low_bits[idx_mod]==word_size) // fits into one word
            return lo;
        else
        {
            // this would be OK even without the test for num_low_bits (hi=0 when ==word_size),
            // bit it's faster this way and we don't need to check data_idx==data.length-1 (page length is a multiple of word length)
            long hi = data[data_idx+1] & high_bits[idx_mod];
            return (hi << num_low_bits[idx_mod]) | lo;
        }
    }


    /**
     * Sets the entry at a given index
     *
     * @param idx array position index (0..{@link #size()}-1)
     * @param value the value to be stored at that index: only the lower bits are used (number of bits given by {@link #getWordSize() }
     */
    public final void set(final long idx, long value)
    {
        long full_idx = (idx*word_size) >>> 6;
        int page_idx = (int)(full_idx / page_size);
        long[] data = pages[page_idx];
        int data_idx = (int)(full_idx % page_size);
        long data_val=data[data_idx];

        int idx_mod = (int)(idx & 63);
        value = value & low_bits[0];
        long lo = value << bit_offset[idx_mod];
        data[data_idx] = (data_val & low_zeros[idx_mod]) | lo;

        if (num_low_bits[idx_mod]!=word_size)
        {
            // set high bits too
            long hi = value >>> num_low_bits[idx_mod];
            data[data_idx+1] = (data[data_idx+1] & high_zeros[idx_mod]) | hi;
//            System.out.println("#*AA.set2 "+idx+"\tval "+value+"\tmod "+idx_mod+"\tnlo "+num_low_bits[idx_mod]+"<<\tdata_idx "+data_idx);
        } else
        {
//            System.out.println("#*AA.set1 "+idx+"\tval "+value+"\tmod "+idx_mod+"\tnlo "+num_low_bits[idx_mod]+"==\tdata_idx "+data_idx);
        }

    }

    /**
     * Fills the array with a given value
     *
     * @param value the value to be stored at every entry: only the lower bits are used (number of bits given by {@link #getWordSize() }
     */
    public final void fill(final long value)
    {
        for (int idx=0; idx<Math.min(64,size); idx++)
            set(idx, value);
        // the first 64 words fit into the first word_size long entries in data[] array
        // data[] is periodic: data[j]=data[j-word_size]
        {
            long[] data = pages[0];
            for (int t=0; t<word_size; t++)
                for (int j=t; j<data.length; j+=word_size)
                {
                    data[j]=data[t];
                }
        }
        for (int page_idx=1; page_idx<pages.length; page_idx++)
            System.arraycopy(pages[0],0,pages[page_idx],0,pages[page_idx].length);
    }

    /**
     * String representation of 64-bit values: like {@link Long#toBinaryString(long) } but padded with zeroes for high bits
     *
     * @param x a long value
     * @return a string of exactly 64 0/1 characters; first character is most significant bit
     */
    public static String toBinaryString(long x)
    {
        String prefix = "0000000000000000000000000000000000000000000000000000000000000000";
        String retval = prefix+Long.toBinaryString(x);
        int l = retval.length();
        return retval.substring(l-64,l);
    }

    /**
     * Used in serialization
     *
     * @param out
     * @throws IOException
     */
    protected void writeObject(ObjectOutputStream out) throws IOException
    {
        out.writeInt(word_size);
        out.writeLong(size);
        for (int page_idx=0; page_idx<pages.length; page_idx++)
        {
            long[] data = pages[page_idx];
            for (int i=0; i<data.length; i++)
                out.writeLong(data[i]);
    //            System.out.println("#*AA.sA wrote "+i);
        }
    }

    /**
     * Used in serialization
     *
     * @param in
     * @throws IOException
     * @throws ClassNotFoundException
     */
    protected void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException
    {
        setWordSize(in.readInt());
        setArraySize(in.readLong());
        for (int page_idx=0; page_idx<pages.length; page_idx++)
        {
            long[] data = pages[page_idx];
            for (int i=0; i<data.length; i++)
            {
    //            System.out.println("#*AA.sA reading "+i);
                data[i] = in.readLong();
            }
        }
    }

    /**
     * Saves the array in a file (uses GZIP compression)
     *
     * @param file
     * @throws IOException
     */
    public final void saveArray(String file) throws IOException
    {
        OutputStream O = new GZIPOutputStream(new java.io.BufferedOutputStream(new java.io.FileOutputStream(file)));
        ObjectOutputStream D = new ObjectOutputStream(O);
        D.writeUTF(getClass().getCanonicalName());
        D.writeLong(serialVersionUID);
        writeObject(D);
        D.flush();
        D.close();
        O.close();
    }


    /**
     * Loads the array from a file (uses GZIP compression)
     *
     * @param file
     * @throws IOException
     */
    public void loadArray(String file) throws IOException
    {
        InputStream I = new GZIPInputStream(new java.io.FileInputStream(file));
        ObjectInputStream D = new ObjectInputStream(new java.io.BufferedInputStream(I));
        String file_class = D.readUTF();
        if (!getClass().getName().equals(file_class))
            throw new IllegalArgumentException("Bad header format: missing or wrong class identifier [got "+file_class+"]");
        long other_version = D.readLong();
        if (other_version>serialVersionUID || other_version<0L)
            throw new IllegalArgumentException("Bad header format: incompatible version id [got "
                    +other_version+", need "+serialVersionUID+" or lower]");
        setWordSize(D.readInt());
        setArraySize(D.readLong());
        for (int page_idx=0; page_idx<pages.length; page_idx++)
        {
            long[] data = pages[page_idx];
            for (int i=0; i<data.length; i++)
            {
    //            System.out.println("#*AA.sA reading "+i);
                data[i] = D.readLong();
            }
        }
        D.close();
        I.close();
    }

    public static class ExtendableArray extends AnybitArray
    {
        private long true_size;
        protected ExtendableArray()
        {
            super();
            true_size = 0L;
        }
        public ExtendableArray(long capacity, int word_size)
        {
            super(capacity, word_size);
            if (((capacity*word_size)&127L)!=0L)
                throw new IllegalArgumentException("Capacity*word_size must be divisible by 128"); // for expansion
            true_size = 0L;
        }
        public void add(long element)
        {
            if (true_size == super.size())
            {
                long expand = (true_size%3L==0L?true_size/3L:true_size/2L); // doubled at every two expansions
                enlargeArray(true_size+expand);
            }
            set(true_size, element);
            ++true_size;
        }
        @Override
        public long size()
        {
            return true_size;
        }

        @Override
        protected void writeObject(ObjectOutputStream out) throws IOException
        {
            super.writeObject(out);
            out.writeLong(true_size);
        }

        @Override
        protected void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException
        {
            super.readObject(in);
            true_size = in.readLong();
        }
    }


    public static void main(String[] args) throws Exception
    {
        int array_size = 4096;
        if (args.length ==0 )
            throw new IllegalArgumentException("Call as $0 file [word_size [array_size]] --- saves and loads the file storing an array of "+array_size+" words");

        String file = args[0];

        int word_size = 12;
        if (args.length > 1 )
            word_size = Integer.parseInt(args[1]);
        if (args.length > 2)
            array_size = Integer.parseInt(args[2]);

        AnybitArray AA = new AnybitArray(array_size, word_size);
        long[] A = new long[array_size];
//        long fill_val  = ((2009L) << (64-word_size))>>>(64-word_size);
//        AA.fill(fill_val);
//        for (int i=0; i<1024; i++)
//            System.out.println("#*AA.main fill "+i+"\t"+AA.get(i)+"\t// should be "+fill_val);

        for (int i=0; i<array_size; i++)
        {
            long val  = (((long)i) << (64-word_size))>>>(64-word_size);
            AA.set(i, val);
            A[i] = val;
        }
//        for (int i=0; i<64; i++)
//            System.out.println("#*AA.main set "+i+"\t"+AA.get(i));

        java.util.Random RND = new java.util.Random(1969);

        System.out.println("Testing array with "+word_size+"-bit integers, length "+array_size);

        for (int i=0; i<100; i++)
        {
            int idx = RND.nextInt(32);
            int v = RND.nextInt(1<<word_size);
            A[idx]=v;
            AA.set(idx,v);
        }
        int num_setok =0;
        for (int i=0; i<64; i++)
        {
            String check_msg = "";
            long stored_val = AA.get(i);
            boolean stored_ok = (stored_val==A[i]);
            if (stored_ok)
            {
                check_msg = "value OK";
                ++num_setok;
            } else
                check_msg = " ****** DIFFERENT *****\tshould be "+A[i];
            System.out.println("#*AA.main get "+i+"\t"+stored_val+"\t// "+check_msg);
        }
        for (int i=64; i<array_size; ++i)
        {
            long stored_val = AA.get(i);
            boolean stored_ok = (AA.get(i)==A[i]);
            if (stored_ok)
                ++num_setok;
        }
        System.out.println("------- Store/retrieve test: "+num_setok+" correct out of "+array_size);

        AA.saveArray(file);

        AA = new AnybitArray();
        AA.loadArray(file);

        int num_loadok=0;
        for (int i=0; i<64; i++)
        {
            String check_msg = "";
            long stored_val = AA.get(i);
            boolean stored_ok = (stored_val==A[i]);
            if (stored_ok)
            {
                check_msg = "value OK";
                ++num_loadok;
            } else
                check_msg = " ****** DIFFERENT *****\tshould be "+A[i];
            System.out.println("#*AA.main load "+i+"\t"+stored_val+"\t// "+check_msg);
        }
        for (int i=64; i<array_size; ++i)
        {
            long stored_val = AA.get(i);
            boolean stored_ok = (stored_val==A[i]);
            if (stored_ok)
                ++num_loadok;
        }
        System.out.println("------- Save/load test: "+num_loadok+" correct out of "+array_size);

    }


}
