/*
 * StringCollector.java
 *
 * Created on May 31, 2005, 2:26 PM
 */

package ca.umontreal.iro.banality;

/**
 *
 * A simple class that makes sure that identical strings are stored only once.
 *
 * @author  csuros
 */

import java.util.Hashtable;

public class StringCollector {
    
    /** Creates a new instance of StringCollector */
    public StringCollector() {
        eddigi=new Hashtable();
    }
    
    private Hashtable eddigi;
    
    /**
     * @return a String that is identical to S.
     */
    public String string(String S){
        if (eddigi.contains(S)){
            return (String) eddigi.get(S);
        } else {
            String S2=new String(S);
            eddigi.put(S2,S2);
            //System.out.println("#**SC.s add "+S2);
            return S2;
        }
    }
}
