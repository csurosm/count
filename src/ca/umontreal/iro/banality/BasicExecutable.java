/*
 * BasicExecutable.java
 *
 * Created on June 25, 2005, 8:56 PM
 */

package ca.umontreal.iro.banality;

/**
 * A simple class that provides the reporting of 
 * calling parameters etc, and handling of Exceptions. 
 *
 * @author  csuros
 */
public class BasicExecutable {
    
    private static final String REPORT_PREFIX="#| ";

    protected String executableInfo()
    {
        return getClass().getName();
    }

    public void reportLaunch(String[] args){
        java.util.Date now = java.util.Calendar.getInstance().getTime();
        java.util.Properties Props=System.getProperties();
        String system = Props.getProperty("os.name", "[unknown OS]")+" "+Props.getProperty("os.version","[Unknown version]")+" "+Props.getProperty("os.arch","[Unknown architecture]");
        String java = Props.getProperty("java.vm.name","[Unknown VM]")+" "+Props.getProperty("java.vm.version","[Unknown version]")+" ("+Props.getProperty("java.runtime.version","[Unknwon runtime]")+") "+Props.getProperty("java.vm.info","")+", "+Props.getProperty("java.vm.vendor","[Uknown vendor]");
        System.out.println(REPORT_PREFIX+executableInfo()+" "+now);
        System.out.println(REPORT_PREFIX+"System: "+system);
        System.out.println(REPORT_PREFIX+"Java: "+java);
        System.out.println(REPORT_PREFIX+"Current directory: "+Props.getProperty("user.dir","[unknown]"));
        System.out.print(REPORT_PREFIX+"Arguments:");
        for (int i=0; i<args.length; i++) System.out.print(" "+args[i]);
        System.out.println();
    }
    
    public void reportOtherArguments(String message){
        System.out.println(REPORT_PREFIX+message);
    }
    
    
    public static void die(Exception E){
        System.out.println("Awfully sorry - something went wrong.\nIf this is surprising, then write to Miklos Csuros [csuros"+"@"+"iro.umontreal.ca] - copy the error message below.");
        E.printStackTrace();
        System.exit(69);
    }
}
