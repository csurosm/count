/*
 * Exportable.java
 *
 * Created on December 29, 2007, 3:54 PM
 *
 */

package ca.umontreal.iro.evolution.malin;

import java.io.File;
import java.io.IOException;

/**
 *
 * Interface for anything that can be saved.
 *
 * @author csuros
 */
public interface Exportable 
{
    /**
     * Saves the data in a new file.
     */
    public void saveData(File f) throws IOException;    
    
}
