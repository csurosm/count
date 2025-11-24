package count.io;
/*
 * Copyright 2023 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


import java.io.File;
import java.util.Random;

/**
 * A class for storing different types of Objects (data) with an associated file.
 *
 * @param <DATA> whatever data is associated with the file
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 * @since December 10, 2007, 9:10 PM
 */
public final class DataFile<DATA> 
{
    
    /** Creates a new instance of DataFile */
    public DataFile(DATA data, File associated_file) 
    {
        this.data = data;
        this.associated_file = associated_file;
        setDirty(associated_file == null || associated_file.getParent()==null);
    }

    /**
     * Initialization with null file. 
     *
     * @param data the object stored here
     */
    public DataFile(DATA data) 
    {
        this(data, null);
    }
    
    
    /**
     * Returns the file associated with this data set
     * @return the associated file (may be null)
     */
    public File getFile()
    {
        return associated_file;
    }

    /**
     * Returns the data set
     */
    public DATA getContent()
    {
        return data;
    }
    
    /**
     * Same as {@link #getContent()} but permits null argument. 
     * 
     * @param <DataType>
     * @param data
     * @return
     */
    public static <DataType> DataType getContent(DataFile<DataType> data)
    {
    	return data==null?null:data.getContent();
    }

    /**
     * Sets the associated file
     */
    public void setFile(File F)
    {
        this.associated_file = F;
        setDirty(associated_file == null || associated_file.getParent()==null);        
    }
    

    /**
     * Sets the data
     */
    public void setData(DATA D)
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
     * Sets the dirty bit (data set is not saved)
     */ 
    public void setDirty(boolean dirty)
    {
        this.dirty=dirty;
    }
    
    public String getNote() {return notes;}
    public void setNote(String notes) {this.notes=notes;}
    
    private DATA data;
    private File associated_file;
    
    private boolean dirty;
    
    private String notes=null;
    
    public static String chopDirectory(String s)
    {
    	if (s==null) return null;
        int last_dir = s.lastIndexOf(File.separatorChar);
        s = s.substring(last_dir+1); // ok if not found (-1+1=0)
        return s;
    }
    
    /**
     * Cuts of common file extensions from the file name
     * 
     * @param filename
     * @return
     */
    public static String chopCommonFileExtension(String filename) {
    	String[] commonSuffixes = {"gz","txt","log","tre"};
    	
        if (filename==null) return null;
        int last_dir = filename.lastIndexOf(File.separatorChar);
        filename = filename.substring(last_dir+1); // ok if not found (-1+1=0)
        int last_dot = filename.lastIndexOf('.');
        while (last_dot != -1) {
        	String ext=filename.substring(last_dot+1);
        	boolean common = false;
        	for (String suf: commonSuffixes) 
        		if (common = suf.equals(ext)) break;
        	if (common) {
            	filename = filename.substring(0, last_dot);
        		last_dot = filename.lastIndexOf('.');
        	} else break;
        }

        return filename;    	
    }
    
    public static String chopFileExtension(String s)
    {
        if (s==null) return null;
        int last_dir = s.lastIndexOf(File.separatorChar);
        s = s.substring(last_dir+1); // ok if not found (-1+1=0)
        int last_dot = s.lastIndexOf('.');
        if (last_dot != -1)
        {
        	String ext=s.substring(last_dot+1);
        	s = s.substring(0, last_dot);
        	if ("gz".equals(ext))
        	{
        		last_dot = s.lastIndexOf('.');
        		s = s.substring(0,last_dot);
        	}
        }
        return s;
//        
//        int first_dot = s.indexOf('.', last_dir+1); // +1 ok if not found
//        return s.substring(last_dir+1, first_dot==-1?s.length():first_dot);
    }
    
    public static String[] random_identifiers 
    = { "Alfa", 
        "Bravo", 
        "Charlie", 
        "Delta", 
        "Echo", 
        "Fox", 
        "Golf", 
        "Hotel", 
        "India", 
        "Juliet", 
        "Kilo", 
        "Lima", 
        "Mike", 
        "Nancy", 
        "Oscar", 
        "Papa", 
        "Quebec", 
        "Romeo", 
        "Sierra", 
        "Tango", 
        "Uncle", 
        "Victor", 
        "Whisky", 
        "Xray", 
        "Yankee", 
        "Zulu"};
//    = {"Anna",
//        "Bob",
//        "Chelsea",
//        "Doug",
//        "Ella",
//        "Frank",
//        "Gina",
//        "Hugo",
//        "Ivy",
//        "Jules",
//        "Kate",
//        "Louis",
//        "Mimi",
//        "Nimrod",
//        "Olga",
//        "Peter",
//        "Queen",
//        "Roland",
//        "Susan",
//        "Timmy",
//        "Ulrike",
//        "Vince",
//        "William",
//        "Xena",
//        "Yannick"
//    };
	public static String anyIdentifier()
	{
		Random RND = new Random();
		return random_identifiers[RND.nextInt(random_identifiers.length)];
	}
	
	
	public static String createIdentifier(int idx)
	{
		if (idx<random_identifiers.length)
		    return random_identifiers[idx];
		else 
		    return "Z"+Integer.toString(idx);
	}
	
	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder("DF[");
		sb.append(getFile()).append(dirty?"/dirty":"/clean");
		sb.append("(").append(data==null?data:data.getClass().getSimpleName()).append(")");
		sb.append("]");
		return sb.toString();
	}
    
    
}
