/*
 * PermutationEnumeration.java
 *
 * Created on September 8, 2008, 3:18 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package ca.umontreal.iro.combinatorics;

/**
 *
 * Enumerates permutations of a multiset in lexicographically increasing order.
 * A direct implementation of Algorithm L from Knuth's TAOCP (Section 7.2.1.2 in Pre-fascicle 2B). 
 *
 * @author csuros
 */
public class PermutationEnumeration implements IntArrayEnumeration
{
    
    /**
     * Initializes an enumeration with an empty set of elements.
     */
    protected PermutationEnumeration() 
    {
        this.elements = new int[0];
    }
    
    /**
     * Where the elements are stored.
     */
    protected int[] elements;
    
    /**
     * Initializes an enumeration for the elements 0..n-1. Must be a positive number.
     */
    public PermutationEnumeration(int n)
    {
        this();
        if (n<0)            
            throw new IllegalArgumentException("Number of elements must be nonnegative at the instantiation of PermutationEnumeration(int).");
        int[] values = new int[n];
        for (int i=0; i<n; i++)
            values[i] = i;
        setElements(values);
        reset();
    }
    
    /**
     * Initializes an enumeration for a multiset of elements.
     */
    public PermutationEnumeration(int[] elements)
    {
        if (elements == null)
            throw new IllegalArgumentException("PermutationEnumeration cannot be instantiated with a null array.");
        this.elements = elements;
        reset();
    }

    /**
     * Sets the elements for the enumeration; does not call reset()
     */
    protected void setElements(int[] elements)
    {
        this.elements = elements;
    }

    /**
     * Next permutation for display
     */
    private int[] next_permutation;
    
    /**
     * Whether there are more permutations in te lexicographic enumeration.
     */
    public boolean hasMoreElements() 
    {
        return (elements.length > 1 && nextIndexToIncrease()>=0) || (elements.length==1 && just_reset);
    }

    /**
     * Finds the next element in the enumeration.
     * @param result the next element will be copied into this array too (must have equal length with the element set used at initialization)
     * @returns an int[] array with identical content as result: the next element in the enumeration
     */
    public Object nextElement(int[] result) throws java.util.NoSuchElementException
    {
        if (!just_reset)
        {
            int j=nextIndexToIncrease();
            if (j<0)
                throw new java.util.NoSuchElementException();
            int l=next_permutation.length-1;
            while (next_permutation[j]>=next_permutation[l])
                l--;
            {
                int x = next_permutation[j];
                next_permutation[j] = next_permutation[l];
                next_permutation[l] = x;
            }
            int k = j+1;
            int m = next_permutation.length-1;
            while (k<m)
            {
                int x = next_permutation[k];
                next_permutation[k] = next_permutation[m];
                next_permutation[m] = x;
                m--;
                k++;
            }
        }
        just_reset = false;
        System.arraycopy(next_permutation,0,result,0,next_permutation.length);
        return result;
    }

    /**
     * Calculates the next element in the lexicographic enumeration.
     * @return an int[] array of equal size as the array used at instantiation
     */
    public Object nextElement() throws java.util.NoSuchElementException
    {
        return nextElement(new int[next_permutation.length]);
    }

    /**
     * Resets the enumeration, so that the next element displayed will be lexicographically the first (i.e., elements sorted in increasing order)
     */
    public void reset()
    {
        next_permutation = new int[elements.length];
        System.arraycopy(elements,0,next_permutation,0,elements.length);
        java.util.Arrays.sort(next_permutation);
        //System.out.print("#PE.res");
        //for (int j=0; j<next_permutation.length; j++)
        //    System.out.print(" "+next_permutation[j]);
        //System.out.println();
        just_reset = true;
    }
    
    private boolean just_reset;
    
    /**
     * Finds the next index that will change
     *
     * @return -1 if done with the enumeration or -2 if array of length 0
     */
    private int nextIndexToIncrease()
    {
        int j=next_permutation.length-2;
        while (j>=0 && next_permutation[j]>=next_permutation[j+1])
            j--;
        return j;
    }

    /**
     * Test code: enumerates the permutations for the arguments
     */
    public static void main(String[] args)
    {
        PermutationEnumeration PE = null;
        if (args.length==1 && Integer.parseInt(args[0])<0)
        {
            int size = - Integer.parseInt(args[0]);
            PE = new PermutationEnumeration(size);
        } else
        {
            int[] emt = new int[args.length];
            for (int i=0; i<args.length; i++)
                emt[i] = Integer.parseInt(args[i]);
            PE = new PermutationEnumeration(emt);
        }
        while (PE.hasMoreElements())
        {
            int[] res = (int[])PE.nextElement();
            for (int i=0; i<res.length; i++)
            {
                if (i!=0)
                    System.out.print(' ');
                System.out.print(res[i]);
            }
            System.out.println();
        }
        PE.reset();
        while (PE.hasMoreElements())
        {
            int[] res = (int[])PE.nextElement();
            System.out.print("#BIS ");
            for (int i=0; i<res.length; i++)
            {
                if (i!=0)
                    System.out.print(' ');
                System.out.print(res[i]);
            }
            System.out.println();
        }
        

    }
}
