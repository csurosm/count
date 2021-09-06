/*
 * Removable.java
 *
 * Created on January 27, 2008, 9:40 PM
 *
 */

package ca.umontreal.iro.evolution.malin;

/**
 * Interface for anything that can be removed.
 *
 * @author csuros
 */
public interface Removable 
{
    /**
     * Prepares the removal of this object (useful e.g., if there are some associated threads that need to be stopped)
     *
     * @return whether the operation was successful
     */
    public boolean remove();
    
}
