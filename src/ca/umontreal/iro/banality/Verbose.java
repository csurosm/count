/*
 * Verbose.java
 *
 * Created on June 24, 2005, 12:12 AM
 */

package ca.umontreal.iro.banality;

/**
 *
 * @author  csuros
 */
public final class Verbose {
    private Verbose(){}
    private static boolean verbose = false;
    
    public static void setVerbose(boolean v){verbose=v;}
    public static boolean isVerbose(){return verbose;}
    
    public static void message(String msg){
        if (verbose)
            System.out.println(prefix+msg);
    }
    
    private static String prefix = "#**";
    public static void setPrefix(String s){prefix=s;}
    public static String getPrefix(){return prefix;}
    
}
