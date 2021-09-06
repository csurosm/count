/*
 * DataFile.java
 *
 * Created on December 10, 2007, 9:10 PM
 *
 */

package ca.umontreal.iro.evolution.malin;

/**
 * A class for storing different types of Objects (data) with an associated file. 
 * 
 * @author csuros
 */

import java.io.File;

public class DataFile<Datatype> 
{
    
    /** Creates a new instance of DataFile */
    public DataFile(Datatype data, File associated_file) 
    {
        this.data = data;
        this.associated_file = associated_file;
        setDirty(associated_file == null || associated_file.getParent()==null);
    }

    /**
     * Initialization with null file. 
     */
    public DataFile(Datatype data) 
    {
        this(data, null);
    }
    
    
    /**
     * Returns the file associated with this data set
     */
    public File getFile()
    {
        return associated_file;
    }

    /**
     * Returns the data set
     */
    public Datatype getData()
    {
        return data;
    }
    

    /**
     * Sets the associated file
     */
    public void setFile(File F)
    {
        this.associated_file = F;
    }
    

    /**
     * Sets the data
     */
    public void setData(Datatype D)
    {
        this.data = D;
    }
    
    /**
     * Whether the data set was modified since the last save
     */
    public boolean isDirty()
    {
        return dirty;
    }
    
    /**
     * Sets the <q>dirty</q> bit (data set is not saved)
     */ 
    public void setDirty(boolean dirty)
    {
        this.dirty=dirty;
    }
    
    private Datatype data;
    private File associated_file;
    
    private boolean dirty;
}
