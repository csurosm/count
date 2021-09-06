/*
 * Saveable.java
 *
 * Created on December 28, 2007, 3:39 PM
 *
 */

package ca.umontreal.iro.evolution.malin;

/**
 *
 * @author csuros
 */
import java.io.File;
import java.io.IOException;

public interface Saveable extends Exportable
{
    /**
     * Whether an associated file is set already for this Object (i.e., <q>Save</q> can be used)
     */
    public boolean hasAssociatedFile();

    /**
     * Saves the data with the current associated file
     */
    public void saveData() throws IOException;
    
    
    /**
     * Whether a save is necessary.
     */
    public boolean isDirty();
    
    /**
     * Sets whether save is necessary.
     */
    public void setDirty(boolean dirty);
    
}
