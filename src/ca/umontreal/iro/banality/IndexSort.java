/*
 * IndexSort.java
 *
 * Created on February 14, 2004, 1:57 PM
 */

package ca.umontreal.iro.banality;

import java.util.Random;
import java.util.List;

/**
 *
 * @author  csuros
 */
public class IndexSort {
    
    /** Creates a new instance of IndexSort */
    private IndexSort() {}
    
    /**
     * Performs QuickSort on the elements of an array. Only indexes are sorted, the 
     * array itself is untouched. Memory requirements are n+2*log_2(n)+O(1) int variables,
     * where n is the length of the array argument. 
     *
     * Implementation based on NR 8.4
     *
     * @return index of elements in ascending order
     */
    public static int[] quickSort(int[] v){
        int n=v.length;
        if (n==0)
            return new int[0];
        else if (n==1)
        {
            int[] retval = new int[1];
            retval[0] = 0;
            return retval;
        }
        // to store subarrays that need to be sorted instead of a recursion
        int jstack=0;
        int[] istack=new int[2*log2(n)+2];
        // this will be the return value
        int[] index=new int[n];
        for (int j=0; j<n; j++) index[j]=j;

        int ir=n-1;
        int l=0;
        
        while (true){
            if (ir-l < SMALL_ARRAY) {
                // insertion sort
                for (int j=l+1;j<=ir;j++) {
                    int indxt=index[j];
                    int a=v[indxt];
                    int i;
                    for (i=j-1;i>=l;i--) {
                        if (v[index[i]] <= a) break;
                        index[i+1]=index[i];
                    }
                    index[i+1]=indxt;
                }
                if (jstack == 0) break;
                ir=istack[jstack--];
                l=istack[jstack--];
            } else {
                int k=(l+ir) >> 1;
                {int temp=index[k]; index[k]=index[l+1]; index[l+1]=temp;}
                if (v[index[l]] > v[index[ir]]) {
                    int temp=index[l]; index[l]=index[ir]; index[ir]=temp;
                }
                if (v[index[l+1]] > v[index[ir]]) {
                    int temp=index[l+1]; index[l+1]=index[ir]; index[ir]=temp;
                }
                if (v[index[l]] > v[index[l+1]]) {
                    int temp=index[l]; index[l]=index[l+1]; index[l+1]=temp;
                }
                int i=l+1;
                int j=ir;
                int indxt=index[l+1];
                
                int a=v[indxt];
                while(true) {
                    do i++; while (v[index[i]]<a);
                    do j--; while (v[index[j]]>a);
                    if (j < i) break;
                    int temp=index[i]; index[i]=index[j]; index[j]=temp;
                }
                index[l+1]=index[j];
                index[j]=indxt;
                jstack += 2;
                if (ir-i+1 >= j-l) {
                    istack[jstack]=ir;
                    istack[jstack-1]=i;
                    ir=j-1;
                } else {
                    istack[jstack]=j-1;
                    istack[jstack-1]=l;
                    l=i;
                }                
            }
        }
        return index;
    }

    /**
     * Performs QuickSort on the elements of an array. Only indexes are sorted, the 
     * array itself is untouched. Memory requirements are n+2*log_2(n)+O(1) int variables,
     * where n is the length of the array argument. 
     *
     * Implementation based on NR 8.4
     *
     * @return index of elements in ascending order
     */
    public static int[] quickSort(double[] v){
        int n=v.length;
        if (n==0)
            return new int[0];
        else if (n==1)
        {
            int[] retval = new int[1];
            retval[0] = 0;
            return retval;
        }
        // to store subarrays that need to be sorted instead of a recursion
        int jstack=0;
        int[] istack=new int[2*log2(n)+2];
        // this will be the return value
        int[] index=new int[n];
        for (int j=0; j<n; j++) index[j]=j;

        int ir=n-1;
        int l=0;
        
        while (true){
            if (ir-l < SMALL_ARRAY) {
                // insertion sort
                for (int j=l+1;j<=ir;j++) {
                    int indxt=index[j];
                    double a=v[indxt];
                    int i;
                    for (i=j-1;i>=l;i--) {
                        if (v[index[i]] <= a) break;
                        index[i+1]=index[i];
                    }
                    index[i+1]=indxt;
                }
                if (jstack == 0) break;
                ir=istack[jstack--];
                l=istack[jstack--];
            } else {
                int k=(l+ir) >> 1;
                {int temp=index[k]; index[k]=index[l+1]; index[l+1]=temp;}
                if (v[index[l]] > v[index[ir]]) {
                    int temp=index[l]; index[l]=index[ir]; index[ir]=temp;
                }
                if (v[index[l+1]] > v[index[ir]]) {
                    int temp=index[l+1]; index[l+1]=index[ir]; index[ir]=temp;
                }
                if (v[index[l]] > v[index[l+1]]) {
                    int temp=index[l]; index[l]=index[l+1]; index[l+1]=temp;
                }
                int i=l+1;
                int j=ir;
                int indxt=index[l+1];
                
                double a=v[indxt];
                while(true) {
                    do i++; while (v[index[i]]<a);
                    do j--; while (v[index[j]]>a);
                    if (j < i) break;
                    int temp=index[i]; index[i]=index[j]; index[j]=temp;
                }
                index[l+1]=index[j];
                index[j]=indxt;
                jstack += 2;
                if (ir-i+1 >= j-l) {
                    istack[jstack]=ir;
                    istack[jstack-1]=i;
                    ir=j-1;
                } else {
                    istack[jstack]=j-1;
                    istack[jstack-1]=l;
                    l=i;
                }                
            }
        }
        return index;
    }

    /**
     * Performs QuickSort on the elements of an array. Only indexes are sorted, the
     * array itself is untouched. Memory requirements are n+2*log_2(n)+O(1) int variables,
     * where n is the length of the array argument.
     *
     * Implementation based on NR 8.4
     *
     * @return index of elements in ascending order
     */
    public static int[] quickSort(String[] v){
        int n=v.length;
        if (n==0)
            return new int[0];
        else if (n==1)
        {
            int[] retval = new int[1];
            retval[0] = 0;
            return retval;
        }
        // to store subarrays that need to be sorted instead of a recursion
        int jstack=0;
        int[] istack=new int[2*log2(n)+2];
        // this will be the return value
        int[] index=new int[n];
        for (int j=0; j<n; j++) index[j]=j;

        int ir=n-1;
        int l=0;

        while (true){
            if (ir-l < SMALL_ARRAY) {
                // insertion sort
                for (int j=l+1;j<=ir;j++) {
                    int indxt=index[j];
                    String a=v[indxt];
                    int i;
                    for (i=j-1;i>=l;i--) {
                        if (v[index[i]].compareTo(a) <= 0) break;
                        index[i+1]=index[i];
                    }
                    index[i+1]=indxt;
                }
                if (jstack == 0) break;
                ir=istack[jstack--];
                l=istack[jstack--];
            } else {
                int k=(l+ir) >> 1;
                {int temp=index[k]; index[k]=index[l+1]; index[l+1]=temp;}
                if (v[index[l]].compareTo(v[index[ir]])>0) {
                    int temp=index[l]; index[l]=index[ir]; index[ir]=temp;
                }
                if (v[index[l+1]].compareTo(v[index[ir]])>0) {
                    int temp=index[l+1]; index[l+1]=index[ir]; index[ir]=temp;
                }
                if (v[index[l]].compareTo(v[index[l+1]])>0) {
                    int temp=index[l]; index[l]=index[l+1]; index[l+1]=temp;
                }
                int i=l+1;
                int j=ir;
                int indxt=index[l+1];

                String a=v[indxt];
                while(true) {
                    do i++; while (v[index[i]].compareTo(a)<0);
                    do j--; while (v[index[j]].compareTo(a)>0);
                    if (j < i) break;
                    int temp=index[i]; index[i]=index[j]; index[j]=temp;
                }
                index[l+1]=index[j];
                index[j]=indxt;
                jstack += 2;
                if (ir-i+1 >= j-l) {
                    istack[jstack]=ir;
                    istack[jstack-1]=i;
                    ir=j-1;
                } else {
                    istack[jstack]=j-1;
                    istack[jstack-1]=l;
                    l=i;
                }
            }
        }
        return index;
    }


    /**
     * Performs QuickSort on the elements of an array. Only indexes are sorted, the
     * array itself is untouched. Memory requirements are n+2*log_2(n)+O(1) int variables,
     * where n is the length of the array argument.
     *
     * Implementation based on NR 8.4
     *
     * @return index of elements in ascending order
     */
    public static int[] quickSort(long[] v){
        int n=v.length;
        if (n==0)
            return new int[0];
        else if (n==1)
        {
            int[] retval = new int[1];
            retval[0] = 0;
            return retval;
        }
        // to store subarrays that need to be sorted instead of a recursion
        int jstack=0;
        int[] istack=new int[2*log2(n)+2];
        // this will be the return value
        int[] index=new int[n];
        for (int j=0; j<n; j++) index[j]=j;

        int ir=n-1;
        int l=0;

        while (true){
            if (ir-l < SMALL_ARRAY) {
                // insertion sort
                for (int j=l+1;j<=ir;j++) {
                    int indxt=index[j];
                    long a=v[indxt];
                    int i;
                    for (i=j-1;i>=l;i--) {
                        if (v[index[i]] <= a) break;
                        index[i+1]=index[i];
                    }
                    index[i+1]=indxt;
                }
                if (jstack == 0) break;
                ir=istack[jstack--];
                l=istack[jstack--];
            } else {
                int k=(l+ir) >> 1;
                {int temp=index[k]; index[k]=index[l+1]; index[l+1]=temp;}
                if (v[index[l]] > v[index[ir]]) {
                    int temp=index[l]; index[l]=index[ir]; index[ir]=temp;
                }
                if (v[index[l+1]] > v[index[ir]]) {
                    int temp=index[l+1]; index[l+1]=index[ir]; index[ir]=temp;
                }
                if (v[index[l]] > v[index[l+1]]) {
                    int temp=index[l]; index[l]=index[l+1]; index[l+1]=temp;
                }
                int i=l+1;
                int j=ir;
                int indxt=index[l+1];

                long a=v[indxt];
                while(true) {
                    do i++; while (v[index[i]]<a);
                    do j--; while (v[index[j]]>a);
                    if (j < i) break;
                    int temp=index[i]; index[i]=index[j]; index[j]=temp;
                }
                index[l+1]=index[j];
                index[j]=indxt;
                jstack += 2;
                if (ir-i+1 >= j-l) {
                    istack[jstack]=ir;
                    istack[jstack-1]=i;
                    ir=j-1;
                } else {
                    istack[jstack]=j-1;
                    istack[jstack-1]=l;
                    l=i;
                }
            }
        }
        return index;
    }

    /**
     * Performs QuickSort on the elements of an array. Only indexes are sorted, the
     * list itself is untouched. Memory requirements are n+2*log_2(n)+O(1) int variables,
     * where n is the length of the list argument.
     *
     * Implementation based on NR 8.4
     *
     * @param <T> Type of the elements in the list
     * @param L list of elements
     * @return index of elements in ascending order
     */
    public static <T extends Comparable<? super T>> int[] quickSort(List<T> L)
    {
        int n=L.size();
        if (n==0)
            return new int[0];
        else if (n==1)
        {
            int[] retval = new int[1];
            retval[0] = 0;
            return retval;
        }
        // to store subarrays that need to be sorted instead of a recursion
        int jstack=0;
        int[] istack=new int[2*log2(n)+2];
        // this will be the return value
        int[] index=new int[n];
        for (int j=0; j<n; j++) index[j]=j;

        int ir=n-1;
        int l=0;

        while (true){
            if (ir-l < SMALL_ARRAY) {
                // insertion sort
                for (int j=l+1;j<=ir;j++) {
                    int indxt=index[j];
                    T a = L.get(indxt);
                    int i;
                    for (i=j-1;i>=l;i--) {
                        if (L.get(index[i]).compareTo(a) <= 0) break;
                        index[i+1]=index[i];
                    }
                    index[i+1]=indxt;
                }
                if (jstack == 0) break;
                ir=istack[jstack--];
                l=istack[jstack--];
            } else {
                int k=(l+ir) >> 1;
                {int temp=index[k]; index[k]=index[l+1]; index[l+1]=temp;}
                if (L.get(index[l]).compareTo(L.get(index[ir]))>0) {
                    int temp=index[l]; index[l]=index[ir]; index[ir]=temp;
                }
                if (L.get(index[l+1]).compareTo(L.get(index[ir]))>0) {
                    int temp=index[l+1]; index[l+1]=index[ir]; index[ir]=temp;
                }
                if (L.get(index[l]).compareTo(L.get(index[l+1]))>0) {
                    int temp=index[l]; index[l]=index[l+1]; index[l+1]=temp;
                }
                int i=l+1;
                int j=ir;
                int indxt=index[l+1];

                T a=L.get(indxt);
                while(true) {
                    do i++; while (L.get(index[i]).compareTo(a)<0);
                    do j--; while (L.get(index[j]).compareTo(a)>0);
                    if (j < i) break;
                    int temp=index[i]; index[i]=index[j]; index[j]=temp;
                }
                index[l+1]=index[j];
                index[j]=indxt;
                jstack += 2;
                if (ir-i+1 >= j-l) {
                    istack[jstack]=ir;
                    istack[jstack-1]=i;
                    ir=j-1;
                } else {
                    istack[jstack]=j-1;
                    istack[jstack-1]=l;
                    l=i;
                }
            }
        }
        return index;
    }

    
    /**
     * Computes the permutated version of an array. 
     * 
     * @param x input array
     * @param indexes order of elements in result
     * @return a permutated y copy of x: y[i] = x[indexes[i]]
     */
    public static double[] permute(double[] x, int[] indexes)
    {
        double[] retval = new double[x.length];
        for (int i=0; i<indexes.length; i++)
            retval[i] = x[indexes[i]];
        return retval;
    }
    
    
    /**
     * Size of a <q>small</q> array: for arrays with this size or less, insertion sort is carried out.
     */
    private static final int SMALL_ARRAY=7;
    
    /**
     * Calculates log of n base 2. Running time is O(log n).
     *
     * @param n a positive integer
     * @return floor of log2(n)
     * @exception IllegalArgumentException if n is not positive
     */
    public static int log2(int n){
        if (n<1)
            throw new IllegalArgumentException("Argument of log2 must be positive.");
        int l=-1;
        while (n>0){
            l++;
            n >>= 1;
        }
        return l;
    }
    
    /**
     * tests quicksort() on a small random array.
     */
    public static void main(String[] args){
        Random RND=new Random();
        int size=2+RND.nextInt(20);
        double[] array=new double[size];
        for (int i=0; i<size; i++)
            array[i]=RND.nextDouble();
        int[] order=quickSort(array);
        boolean[] seen=new boolean[size];
        for (int i=0; i<size; i++){
            int idx=order[i];
            System.out.println("#** i "+idx+" v "+array[idx]);
            seen[idx]=true;            
        }
        for (int i=0; i<size; i++)
            if (!seen[i])
                System.out.println(" not seen "+i);
    }
    
}
