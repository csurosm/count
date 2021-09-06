/*
 * IntArrayEnumeration.java
 *
 * Created on April 3, 2002, 1:39 PM
 */

package ca.umontreal.iro.combinatorics;

/**
 *
 * @author  miki
 */
public interface IntArrayEnumeration extends java.util.Enumeration {
    public Object nextElement(int[] result) throws java.util.NoSuchElementException;
    public void reset();
}
