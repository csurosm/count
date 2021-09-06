/*
 * SimpleObject.java
 *
 * Created on April 10, 2003, 10:03 AM
 */

package ca.umontreal.iro.banality;

/**
 * Defines the standard toString: class name+paramString.
 *
 * @author  csuros
 */
public class SimpleObject {

    /**
     * Called by toString().
     */
    protected String paramString(){
        return "";
    }

    /**
     * Standard toString() procedure: uses class name and paramString()
     */
    public String toString(){
        StringBuffer sb=new StringBuffer(getClass().getName());
        sb.append('[');
        sb.append(paramString());
        sb.append(']');
        return sb.toString();
    }

}
