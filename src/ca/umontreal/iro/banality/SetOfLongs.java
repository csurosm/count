/*
 * LongSet.java
 *
 * Created on June 13, 2005, 2:55 PM
 */

package ca.umontreal.iro.banality;

/**
 *
 * @author  csuros
 */
public interface SetOfLongs {
    
    /**
     * @return whether the given value is in the set
     */
    public boolean get(long value);
    
    /**
     * Sets the presence of value to true or false (add to or remove from the set, respectively)
     */
    public void set(long value, boolean present);
    
    /** 
     * Adds the value to the set.
     */
    public void set(long value);
    
    /**
     * Removes value from the set.
     */
    public void clear(long value);
    
    /**
     * @return number of elements in the set.
     */
    public int getCardinality();
    
        
    /** 
     * Adds the value to the set.
     *
     * @return whether it was there before already
     */
    public boolean getset(long value);
    
    
    
    
}
