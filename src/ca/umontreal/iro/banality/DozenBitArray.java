package ca.umontreal.iro.banality;

import java.io.IOException;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.BufferedReader;
import java.io.PrintStream;

/**
 * 
 * This class implements an array of 12-bit integers using about <var>n</var>*1.5 bytes
 * for <var>n</var> elements. 
 * The implementation puts a certain number of elements into a <em>page</em>, 
 * which is just an integer array of some specified size (default setting: 768k). 
 * A page contains <var>ps</var>*8/3 elements, where <var>ps</var> 
 * is the size for the integer array supporting the page. The number of pages is 
 * <var>n</var>/(<var>ps</var>*8/3), which yields a negligible overhead in memory usage
 * (12 bytes per page). For example, when <var>n</var>=2<sup>30</sup>, the default setting
 * uses 512 pages, meaning an overhead of 6kbytes over the necessary 1.5Gbytes.
 * 
 * @author csuros
 */

public class DozenBitArray 
{
    private static final double VERSION_ID = 10.0619;
    private static final int PAGE_SIZE = 3<<18; // this many int's (768k) need 3 Mbytes: can store 2M int12 values
    
    /**
     * Initializes an array to store a given number of elements
     * 
     * @param num_elements how many int12 values will be used. 
     */
    public DozenBitArray(long num_elements)
    {
        this(num_elements, PAGE_SIZE);
    }

    
    /**
     * Initializes an array to store a given number of elements
     * 
     * @param num_elements how many int12 values will be used. 
     * @param page_size engineering parameter for storage: this many (32-bit) int values are allocated in a block
     */
    public DozenBitArray(long num_elements, int page_size)
    {
        this.num_elements = num_elements;
        this.page_size = page_size;
        if (page_size % 3 != 0)
            throw new IllegalArgumentException("Page size must be a multiple of three (got "+page_size+" instead)");
        init();
    }
    
    private int page_size;
    private long num_elements;
    private long page_size_int12;
    private int[][] data;

    //private int[] debug_array;
    
    /**
     * Number of 12-bit integers in this array. 
     * 
     * @return array size
     */
    public long size()
    {
        return num_elements;
    }
    
    /**
     * Initializes the data array after num_elements and page_size are set
     */
    private void init()
    {
        long num_entries_int = (num_elements*12L)/32L;
        if (32L*num_entries_int < 12L*num_elements)
            num_entries_int++;
        long num_pages = (num_entries_int / page_size);
        if (num_pages*page_size<num_entries_int)
            num_pages++;
        if (num_pages > Integer.MAX_VALUE)
            throw new IllegalArgumentException("Cannot allocate storage in DozenBitArray: increase page size or decrease the array size (now need "+num_pages+" pages)");
        
        data = new int[(int)num_pages][page_size];
        page_size_int12 = (8L*(long)page_size)/3L; 
        
        //debug_array = new int[(int)num_elements]; 
    }
    
    private static final int BITS_UPPER_4  = 0xF0000000; 
    private static final int BITS_LOWER_28 = 0x0FFFFFFF;
    private static final int BITS_UPPER_8  = 0xFF000000;
    private static final int BITS_LOWER_24 = 0x00FFFFFF;
    
    private static final int BITS_LOWER_4  = 0x0000000F;
    private static final int BITS_UPPER_28 = 0xFFFFFFF0;
    private static final int BITS_LOWER_8  = 0x000000FF;
    private static final int BITS_UPPER_24 = 0xFFFFFF00;

    private static final int BITS_DOZEN_0  = 0xFFF00000;
    private static final int BITS_DOZEN_C0 = 0x000FFFFF;
    private static final int BITS_DOZEN_1  = 0x0FFF0000;
    private static final int BITS_DOZEN_C1 = 0xF000FFFF;
    private static final int BITS_DOZEN_2  = 0x00FFF000;
    private static final int BITS_DOZEN_C2 = 0xFF000FFF;
    private static final int BITS_DOZEN_3  = 0x000FFF00;
    private static final int BITS_DOZEN_C3 = 0xFFF000FF;
    private static final int BITS_DOZEN_4  = 0x0000FFF0;
    private static final int BITS_DOZEN_C4 = 0xFFFF000F;
    private static final int BITS_DOZEN_5  = 0x00000FFF;
    private static final int BITS_DOZEN_C5 = 0xFFFFF000;
    
    
    /**
     * Access to an element of the array
     * @param index element index
     * @return a 12-bit integer value (0..4095) stored at the given index
     */
    public final int get(long index)
    {
        int page_idx = (int)(index / page_size_int12);
        int page_offset = (int)(index % page_size_int12);
        
        int[] intval = data[page_idx];
        
        int int_idx = (page_offset * 3) / 8;
        int halfbyte_offset = (page_offset*3) % 8;
        
        int retval=-1;
        
        //System.out.println("#**DBA.get("+index+")\tpage "+page_idx+"+"+page_offset+"\tint "+int_idx+"+"+halfbyte_offset);

        if (halfbyte_offset == 0)
        {
            // top 3 halfbytes
            retval = intval[int_idx] & BITS_DOZEN_0;
            retval = retval >>> 20;
        } else if (halfbyte_offset == 1)
        {
            retval = intval[int_idx] & BITS_DOZEN_1;
            retval = retval >>> 16;
        } else if (halfbyte_offset == 2)
        {
            retval = intval[int_idx] & BITS_DOZEN_2;
            retval = retval >>> 12;
        } else if (halfbyte_offset == 3)
        {
            retval = intval[int_idx] & BITS_DOZEN_3;
            retval = retval >>> 8;
        } else if (halfbyte_offset == 4)
        {
            retval = intval[int_idx] & BITS_DOZEN_4;
            retval = retval >>> 4;
        } else if (halfbyte_offset == 5)
        {
            retval = intval[int_idx] & BITS_DOZEN_5;
        } else if (halfbyte_offset == 6)
        {
            // lower 8 bits at this index and 4 upper bits from next 
            int a = intval[int_idx] & BITS_LOWER_8;
            int b = (intval[int_idx+1] & BITS_UPPER_4) >>> 28;
            retval = (a<<4) | b;
        } else if (halfbyte_offset == 7)
        {
            // lower 4 bits at this index and 8 upper bits in next int
            int a = intval[int_idx] & BITS_LOWER_4;
            int b = (intval[int_idx+1] & BITS_UPPER_8) >>> 24;
            retval = (a<<8) | b;
        }
        
        //{
        //   int r = debug_array[(int)index];;
        //    if (r!=retval)
        //    {
        //        System.out.println("Error!: should have "+r+"/"+Integer.toHexString(r)+"\tgot "+retval+"/"+Integer.toHexString(retval));
        //        System.out.println("Half-byte offset "+halfbyte_offset);
        //        System.out.println("int "+Integer.toHexString(intval[int_idx])+", "+Integer.toHexString(intval[int_idx+1]));
        //    }
        //}
        
        return retval;
    }
    
    public final void set(long index, int x)
    {
        int page_idx = (int)(index / page_size_int12);
        int page_offset = (int)(index % page_size_int12);
        
        int[] intval = data[page_idx];
        
        int int_idx = (page_offset * 3) / 8;
        int halfbyte_offset = (page_offset*3) % 8;
        
        
        if (halfbyte_offset == 0)
        {
            // top 3 halfbytes
            intval[int_idx] = intval[int_idx] & BITS_DOZEN_C0;
            intval[int_idx] = intval[int_idx] | (x<<20);
        } else if (halfbyte_offset == 1)
        {
            intval[int_idx] = intval[int_idx] & BITS_DOZEN_C1;
            intval[int_idx] = intval[int_idx] | (x<<16);
        } else if (halfbyte_offset == 2)
        {
            intval[int_idx] = intval[int_idx] & BITS_DOZEN_C2;
            intval[int_idx] = intval[int_idx] | (x<<12);
        } else if (halfbyte_offset == 3)
        {
            intval[int_idx] = intval[int_idx] & BITS_DOZEN_C3;
            intval[int_idx] = intval[int_idx] | (x<<8);
        } else if (halfbyte_offset == 4)
        {
            intval[int_idx] = intval[int_idx] & BITS_DOZEN_C4;
            intval[int_idx] = intval[int_idx] | (x<<4);
        } else if (halfbyte_offset == 5)
        {
            intval[int_idx] = intval[int_idx] & BITS_DOZEN_C5;
            intval[int_idx] = intval[int_idx] | x;
        } else if (halfbyte_offset == 6)
        {
            // lower 8 bits of this index and 4 upper bits at next int entry
            int xu = x >>> 4;
            int xl = x & BITS_LOWER_4;
            intval[int_idx] = intval[int_idx] & BITS_UPPER_24;
            intval[int_idx] = intval[int_idx] | xu;
            int next_idx = int_idx+1;
            intval[next_idx] = intval[next_idx] & BITS_LOWER_28;
            intval[next_idx] = intval[next_idx] | (xl << 28);
        } else if (halfbyte_offset == 7)
        {
            // lower 4 bits at this index and 8 upper bits in next int
            int xu = x >>> 8;
            int xl = x & BITS_LOWER_8;
            intval[int_idx] = intval[int_idx] & BITS_UPPER_28;
            intval[int_idx] = intval[int_idx] | xu;
            int next_idx = int_idx+1;
            intval[next_idx] = intval[next_idx] & BITS_LOWER_24;
            intval[next_idx] = intval[next_idx] | (xl << 24);
        }        

        //System.out.println("#**DBA.set("+index+")"+"\t<- "+x+"\tpage "+page_idx+"+"+page_offset+"\tint "+int_idx+"+"+halfbyte_offset+"\tcontent "+Integer.toHexString(intval[int_idx]));
    
        
        //debug_array[(int)index]=x;
    }
    
    /**
     * Fills up the array with a given value: every entry is set to this
     * 
     * @param x the value
     */
    public final void fill(int x)
    {
        int e0 = x;
        e0 = (x<<12) | x;
        int x_top_8 = x>>>4;
        e0 = (e0<<8) | x_top_8;

        int x_bottom_4 = x & BITS_LOWER_4;
        int e1 = x_bottom_4;
        e1 = (e1<<12) | x;
        e1 = (e1<<12) | x;
        int x_top_4 = x>>>8;
        e1 = (e1<<4) | x_top_4;
        
        int x_bottom_8 = x & BITS_LOWER_8;
        int e2 = x_bottom_8;
        e2 = (e2<<12) | x;
        e2 = (e2<<12) | x;
        
        //System.out.println("#**DBA.f "+x+"/"+Integer.toHexString(x)
        //        +"\t"+Integer.toHexString(e0)
        //        +", "+Integer.toHexString(e1)
        //        +", "+Integer.toHexString(e2));
        
        for (int pidx=0; pidx<data.length; pidx++)
        {
            int[] A = data[pidx];
            for (int j=0; j<A.length; j+=3)
                A[j]=e0;
            for (int j=1; j<A.length; j+=3)
                A[j]=e1;
            for (int j=2; j<A.length; j+=3)
                A[j]=e2;
            
        }
    }

    public void saveArray(String file) throws IOException
    {
        PrintStream P = new PrintStream(new GZIPOutputStream(new java.io.BufferedOutputStream(new java.io.FileOutputStream(file))));

        writeArray(P);

        P.flush();
        P.close();
    }

    /**
     * Writes the attay content into a text file.
     *
     * @param P stream to which the output is written.
     */
    public void writeArray(PrintStream P)
    {
        P.println("# "+getClass().getCanonicalName()+"\tV"+VERSION_ID+"\t"+num_elements+"\t"+page_size);
        HalfbyteArray.writePages(data,P);
    }

    public void readArray(BufferedReader R) throws IOException
    {
        String header = R.readLine();
        if (header == null || !header.startsWith("# "))
        {
            throw new IOException("Truncated file with no header");
        }
        String[] header_fields = header.split("\t");
        if (header_fields.length <4
                || !header_fields[0].equals("# "+getClass().getName())
                || !header_fields[1].startsWith("V"))
        {
            throw new IllegalArgumentException("Bad header format");
        }

        double file_version = Double.parseDouble(header_fields[1].substring(1));

        if (file_version > VERSION_ID)
        {
            throw new IllegalArgumentException("File version ("+file_version+") is higher than this version ("+VERSION_ID+")");
        }
        num_elements = Long.parseLong(header_fields[2]);
        page_size = Integer.parseInt(header_fields[3]);
        init();

        int num_pages = data.length;
        data = null;

//        System.out.println("#*DBA.rA numpages "+num_pages+"\tpagesize "+page_size+"\tnumemts "+num_elements);

        data = HalfbyteArray.readPages(R, num_pages, page_size);
    }

    public void loadArray(String file) throws IOException
    {
        BufferedReader R
            = new BufferedReader(new java.io.InputStreamReader(new GZIPInputStream(new java.io.FileInputStream(file))));
        readArray(R);
        R.close();
    }
    
    /**
     * Test code.
     * 
     * @param args
     */
    public static void main(String[] args) throws IOException
    {

        if (args.length != 1)
            throw new IllegalArgumentException("Call as $0 file ");
        DozenBitArray DBA = new DozenBitArray(1L<<30);
        int[] A = new int[4096];
        DBA.fill(2009);
//        for (int i=0; i<64; i++)
//            System.out.println(i+"\t"+DBA.get(i));

        for (int i=0; i<4096; i++)
        {
            DBA.set(i, i);
            A[i] = i;
        }
        //for (int i=0; i<4096; i++)
        //    System.out.println(i+"\t"+DBA.get(i));
        
        java.util.Random RND = new java.util.Random(1969);
        
        for (int i=0; i<100; i++)
        {
            int idx = RND.nextInt(32);
            int v = RND.nextInt(4096);
            A[idx]=v;
            DBA.set(idx,v);
        }
        for (int i=0; i<64; i++)
            System.out.println("#*DBA.main saved "+i+"\t"+DBA.get(i)+"\t// "+A[i]);
        
        DBA.saveArray(args[0]);
        DBA.loadArray(args[0]);
        for (int i=0; i<64; i++)
            System.out.println("#*DBA.main loaded "+i+"\t"+DBA.get(i)+"\t// "+A[i]);

    }
    
}
