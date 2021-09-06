
package ca.umontreal.iro.banality;

/**
 * 
 * This class implements an array of 4-bit integers using about <var>n</var>*0.5 bytes
 * for <var>n</var> elements. 
 * The implementation puts a certain number of elements into a <em>page</em>, 
 * which is just an int array of some specified size (default setting: 1M). 
 * 
 * @author csuros
 */

import java.io.PrintStream;
import java.io.BufferedReader;

public class HalfbyteArray 
{
    private static final String VERSION_ID = "1";
    
    private static final int PAGE_SIZE = 1<<20; // 1M integers = 4Mbytes = 8M halfbytes
    
    /**
     * Initializes an array to store a given number of elements
     * 
     * @param num_elements how many half-bytes will be used. 
     */
    public HalfbyteArray(long num_elements)
    {
        this(num_elements, PAGE_SIZE);
    }

    
    /**
     * Initializes an array to store a given number of elements
     * 
     * @param num_elements how many half-bytes will be used. 
     * @param page_size engineering parameter for storage: this many int values are allocated in a block
     */
    public HalfbyteArray(long num_elements, int page_size)
    {
        this.num_elements = num_elements;
        this.page_size = page_size;
        init();
    }
    
    private int page_size;
    private long num_elements;
    private int[][] data;
    private long page_size_halfbytes;

    /**
     * Number of half-bytes in this array. 
     * 
     * @return array size
     */
    public long size()
    {
        return num_elements;
    }
    
    /**
     * Initializes the data array
     */
    private void init()
    {
        long num_entries_int = (num_elements+7L)/8L ;
        long num_pages = (num_entries_int / page_size);
        if (num_pages*page_size<num_entries_int)
            num_pages++;
        if (num_pages > Integer.MAX_VALUE)
            throw new IllegalArgumentException("Cannot allocate storage in HalfbyteArray: increase page size or decrease the array size (now need "+num_pages+" pages)");
        
        data = new int[(int)num_pages][page_size];
        page_size_halfbytes = 8L*page_size;
    }
        

    /**
     * Access to an element of the array
     * @param index element index
     * @return a 4-bit integer value (0..15) stored at the given index
     */
    public final int get(long index)
    {
        int page_idx = (int)(index / page_size_halfbytes);
        int page_offset = (int)(index % page_size_halfbytes);
        
        int int_idx = page_offset / 8;
        int int_val = data[page_idx][int_idx];

        int int_offset = page_offset % 8;
        if (int_offset>0)
        {
            int bit_offset = 4*int_offset;
            int_val = int_val >>> bit_offset;
        }
        return int_val & 0xF;
    }
    
    public final void set(long index, int x)
    {
        int page_idx = (int)(index / page_size_halfbytes);
        long page_offset = (int)(index % page_size_halfbytes);
        
        int[] dpage = data[page_idx];
        int int_idx = (int)(page_offset / 8L);
        int int_val = dpage[int_idx];

        int int_offset = (int)(page_offset % 8L);
        int bit_offset = 4*int_offset;
        int hb_selector = 0xF<<bit_offset;
        int mask = 0xFFFFFFFF ^ hb_selector;
        int_val = (int_val & mask) | (x<<bit_offset);
        dpage[int_idx] = int_val;
    }
    
    /**
     * Writes the attay content into a text file.
     *
     * @param P stream to which the output is written.
     */
    public void writeArray(PrintStream P)
    {
        P.println("# "+getClass().getCanonicalName()+"\tV"+VERSION_ID+"\t"+num_elements+"\t"+page_size);
        writePages(data,P);
    }
    
    public static void writePages(int[][] pages, PrintStream P) 
    {
        int num_pages = pages.length;
        for (int pidx=0; pidx<num_pages; pidx++)
        {
            int[] page = pages[pidx];
            int elements_in_line = 0;
            for (int i=0; i<page.length; i++)
            {
                if (elements_in_line==8)
                {
                    P.println();
                    elements_in_line = 0;
                }
                if (elements_in_line != 0)
                    P.print('\t');
                P.print(Integer.toString(page[i],16));
                elements_in_line++;
            }
            P.println();
        }
        P.flush();
//        System.out.println("#*HA.wP numpages "+num_pages+"\tpsize "+pages[0].length);
    }
    
    public static int[][] readPages(BufferedReader R, int num_pages, int psize) throws java.io.IOException
    {
        int[][] x = new int[num_pages][psize];
        long num_elements_read = 0L;
        long num_elements_needed = ((long)num_pages)*((long)psize);
        int current_page_idx = 0;
        int current_page_offset = 0;
        do
        {
            String line = R.readLine();
            if (line == null) // end-of-read
                throw new java.io.IOException("Truncated data: have "+num_elements_read+" of "+num_elements_needed+" ("+num_pages+"*"+psize+")");
            if (line.startsWith("#")) // comments
                continue;
            String[] values = line.split("\\s+");
            for (int i=0; i<values.length; i++)
            {
//                System.out.println("#*HA.rP offset "+current_page_offset+"\tidx "+current_page_idx+"\tfield "+i+"\tval "+values[i]);
                int v = Integer.parseInt(values[i],16);
                x[current_page_idx][current_page_offset]=v;
                current_page_offset++;
                if (current_page_offset==psize)
                {
                    current_page_offset=0;
                    current_page_idx++;
                }
            }
            num_elements_read+=values.length;
        } while (num_elements_read<num_elements_needed);
        return x;
    }
    
    /**
     * Test code.
     * 
     * @param args
     */
    public static void main(String[] args) throws Exception
    {
        int N = 129;
        HalfbyteArray HBA = new HalfbyteArray(N,16);
        int[] A = new int[N];

        for (int i=0; i<N; i++)
        {
            int v = i & 0xF;
            HBA.set(i, v);
            A[i] = v;
        }
        
        java.util.Random RND = new java.util.Random();
        
        for (int i=0; i<49; i++)
        {
            int idx = RND.nextInt(32);
            int v = RND.nextInt(16);
            A[idx]=v;
            HBA.set(idx,v);
        }
        for (int i=0; i<N; i++)
            System.out.println(i+"\t"+HBA.get(i)+"\t// "+A[i]);
        
        HBA.writeArray(System.out);
    }

}
