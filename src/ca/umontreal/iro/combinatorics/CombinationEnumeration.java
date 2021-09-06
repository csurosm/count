/*
 * CombinationEnumeration.java
 *
 * Created on April 3, 2002, 1:46 PM
 */

package ca.umontreal.iro.combinatorics;
import java.util.NoSuchElementException;

/**
 * combination of given size with values from the set 0...<var>maxValue<var>-1
 * 
 * @author  miki
 */
public class CombinationEnumeration implements IntArrayEnumeration {
    
    private int length;
    private int maxValue;
    private int[] element;
    private boolean thereIsMore;
    
    /** Creates a new instance of CombinationEnumeration 
     * @param length number of elements chosen
     * @param maxValue the set 0..<var>maxValue</var>-1 is used for the selection of elements
     */
    public CombinationEnumeration(int length, int maxValue) {
        this.length=length;
        this.maxValue=maxValue;
        if (length>maxValue)
            throw new IllegalArgumentException("Length ["+length+"] cannot be larger than maxValue ["+maxValue+"]");
        element=new int[length];
        reset();
    }
    
    public void reset(){
        //System.out.println("#**CE.rs > "+toString());
        for(int i=0; i<length; i++)
            element[i]=length-i-1;
        thereIsMore = (length<=maxValue);
        //System.out.println("#**CE.rs <"+toString());
    }
    
    private void turn(){
        int pos=0;
        while(pos<length && element[pos]==maxValue-1-pos) pos++;
        if (pos<length){
            element[pos]++;
            for (int i=pos-1; i>=0 ; i--)
                element[i]=element[i+1]+1;
        } else 
            thereIsMore=false;
    }

    public boolean hasMoreElements(){return thereIsMore;}
    /**
     * @return int[] array of length given to the constructor: the elements are in decreasing order
     */
    public Object nextElement() throws NoSuchElementException{
        return nextElement(new int[length]);
    }
    
    public Object nextElement(int[] result) throws NoSuchElementException {
        if (!thereIsMore) throw new NoSuchElementException();
        if (length==0){
            turn();
            return result;
        }
        if (length>2) System.arraycopy(element,0,result,0,length);
        else {
            result[0]=element[0];
            if (length>1) result[1]=element[1];
        }
        turn();
        return result;
    }
    
    protected String paramString(){
        StringBuffer sb=new StringBuffer();
        sb.append("len=");
        sb.append(length);
        sb.append(" max=");
        sb.append(maxValue);
        sb.append(" more=");
        sb.append(thereIsMore);
        sb.append(" emt={");
        for (int i=0; i<length; i++){
            if (i>0) sb.append(',');
            sb.append(element[i]);
        }
       sb.append('}');
       return sb.toString();
    }
    
    public String toString(){
        StringBuffer sb=new StringBuffer(this.getClass().getName());
        sb.append('[');        
        sb.append(paramString());
        sb.append(']');
        return sb.toString();
    }
    
    /**
     * @return a String of style "{a, b, }" for the array
     */
    public static String toString(int[] array){
        StringBuffer sb=new StringBuffer("{");
        for (int i=0; i<array.length; i++){
            if (i>0)
                sb.append(",");
            sb.append(array[i]);
        }
        sb.append('}');
        return sb.toString();
    }

    /**
     * @param k number of elements chosen
     * @param n base set size
     * @return number of combinations k choose n
     */ 
    public static long getNumCombinations(int k, int n){
        if (k>n/2) k=n-k;
        long c=1;
        for (int i=1; i<=k; i++){
            c = (c*(n-i+1))/i;
        }
        return c;
    }
    
    private static java.util.Random my_own_RND=null; // to keep between successive calls
    /**
     * @param k number of elements chosen 
     * @param n size of the int array from which they are chosen (range 0..n-1, inclusive)
     * @return a (k choose n) combination picked uniform randomly
     */
    public static int[] getRandomCombination(int k, int n){
        if (my_own_RND==null)
            my_own_RND = new java.util.Random();
        return getRandomCombination(k,n,my_own_RND);
    }
    
    /**
     * @param k number of elements chosen 
     * @param n size of the int array from which they are chosen (range 0..n-1, inclusive)
     * @param RND the random number generator to use
     *
     * @return a (k choose n) combination picked uniform randomly
     */
    public static int[] getRandomCombination(int k, int n, java.util.Random RND){
        if (k>n)
            throw new IllegalArgumentException("Number of elements chosen must be between 0 and "+n+" [want "+k+"]");
        int[] chosen = new int[k];
        int num_chosen=0;
        for (int pos=0; pos<n; pos++){
            double p=(k-num_chosen)/(n-pos+0.);
            if (RND.nextDouble()<p){
                chosen[num_chosen]=pos;
                num_chosen++;
            } 
        }
        return chosen;
    }
    
    /**
     * Tests
     */
    public static void main(String[] args){
        System.out.println("(30 choose 20)="+getNumCombinations(20,30));
        for (int i=0; i<3; i++){
            int [] a=getRandomCombination(10,20);
            System.out.println("random(10 choose 20) "+toString(a));
        }
        
    }
}