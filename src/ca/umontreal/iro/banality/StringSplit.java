/*
 * StringSplit.java
 *
 * Created on June 1, 2005, 1:08 PM
 */

package ca.umontreal.iro.banality;

/**
 *
 * @author  csuros
 */
public final class StringSplit {
    
    /** Creates a new instance of StringSplit */
    private StringSplit() {}
    
    /**
     * Splits a String at TABs (\t) : faster than using String.split()
     */
    public static final String[] splitAtTabs(String S){
        return splitAt(S, '\t');
    }
    
    /**
     * Splits a String at sep : faster than using String.split()
     */
    public static final String[] splitAt(String S, char sep){
        intVector tab_positions = new intVector();
        int pos = S.indexOf(sep);
        while (pos != -1){
            tab_positions.add(pos);
            pos = S.indexOf(sep,pos+1);
        }
        int num_tabs = tab_positions.size();
        String[] fields = new String[num_tabs+1];
        int current_pos=0;
        for (int i=0; i<num_tabs; i++){
            int next_tab = tab_positions.get(i);
            fields[i]=S.substring(current_pos,next_tab);
            current_pos=next_tab+1;
        }
        fields[num_tabs]=S.substring(current_pos);
        return fields;
        
    }
    
    
}
