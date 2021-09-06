package ca.umontreal.iro.banality;

/**
 * This class is for storing integers of a fixed number of bits. 
 * A word array is an array of w-bit words. For a size N array, 
 * Nw+o(Nw) bits are used.   
 * 
 * @author csuros
 */
public class WordArray 
{
    private static final int PAGE_UNIT = 1<<18; // 0.25M: so many int's [4 bytes] take 1 Mbyte
    
    
    private int word_length;
    private int num_words;
    
    private int[][] data;
    
    private static int[] bit_position_mask=null;
    private static int[][] bit_interval_mask=null;
    
    public WordArray(int word_length, int num_words)
    {
        this.word_length = word_length;
        this.num_words = num_words;
        init();
    }
    
    private void init()
    {
        if (bit_position_mask==null)
        {
            bit_position_mask = new int[32];
            for (int i=0; i<32; i++)
                bit_position_mask[i] = 1<<i;
        }

        int page_size_ints = word_length * PAGE_UNIT; // this many ints in one page
        int page_size = PAGE_UNIT*32; // this many words in one page
        int num_pages = num_words/page_size;
        if (num_pages*page_size<num_words) num_pages++; 
        data = new int[num_pages][page_size_ints];
    }

    
    public static int getWord(int[] array, int word_length, int word_idx)
    {
        int bit_idx = word_length*word_idx; // position for first bit 
        int j = bit_idx / 32; // position of word within the array
        int u = array[j];
        int upper_length = bit_idx-32*j;
        if (upper_length <= word_length) // OK, all of the word fits within this guy
        {
            int w = u>>>(32-word_length);
            return w;
        } else // word overlaps between two int entries
        {
            u = u>>>(32-upper_length);
            int lower_length = word_length-upper_length;
            u = u<<lower_length;
            
        }
        return 0;
    }

}
