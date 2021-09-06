/*
 * MultiIntArrayEnumeration.java
 *
 * Created on April 3, 2002, 2:08 PM
 */

package ca.umontreal.iro.combinatorics;
import java.util.NoSuchElementException;

/**
 *
 * @author  miki
 */
public class MultiIntArrayEnumeration implements java.util.Enumeration {
    
    private IntArrayEnumeration[] elements;
    
    private int length;
    
    private boolean thereIsMore;
    
    private int[][] values;
    
    private boolean firstCall;
    
    /** Creates a new instance of MultiIntArrayEnumeration */
    public MultiIntArrayEnumeration(IntArrayEnumeration[] elements) {
        this.length=elements.length;
        this.elements=elements;
        values=new int[length][];
        
        firstCall=true;

        thereIsMore=true;
        checkIfThereIsMore();
    }
    
    public boolean hasMoreElements() {return thereIsMore;}
    
    public Object nextElement() throws java.util.NoSuchElementException {
        if (!thereIsMore)
                throw new NoSuchElementException();
        turn();
        return values;
    }
    
    private void checkIfThereIsMore() {
        for (int i=0; i<length; i++)
            if (elements[i] != null && !elements[i].hasMoreElements()){
                thereIsMore=false;
                break;
            }
    }
    
    private void turn(){
        if (thereIsMore)
            if (firstCall){            
                for (int i=0; i<length; i++)
                    if (elements[i]==null)
                        values[i]=new int[0];
                    else
                        values[i]=(int[]) elements[i].nextElement();
                firstCall=false;
            } else {
                int pos=0;
                while (pos<length && (elements[pos]==null || !elements[pos].hasMoreElements())) pos++;
                if (pos==length){
                    thereIsMore=false;
                } else {
                    //System.out.println("#**MIAE.t pos="+pos+" e="+enum[pos]);
                    elements[pos].nextElement(values[pos]);
                    for (int i=pos-1; i>=0; i--)
                        if (elements[i]!=null){
                            elements[i].reset();
                            elements[i].nextElement(values[i]);
                        }
                    
                }
            }
        
    }
    
    protected String paramString(){
        StringBuffer sb=new StringBuffer();
        sb.append("len=");
        sb.append(length);
        sb.append(" enum={");
        for (int i=0; i<length; i++){
            if (i>0) sb.append(',');
            sb.append(elements[i].toString());
        }
        return sb.toString();
    }
    
    public final String toString(){
        StringBuffer sb=new StringBuffer(this.getClass().getName());
        sb.append('[');
        sb.append(paramString());
        sb.append(']');
        return sb.toString();
    }
    
}
