package count.machine;


/**
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 */
public class App 
{
    private App(){}
    
    public static final String TITLE = "Evolution Machine";
    public static final String VERSION = "I";
    
    public static String getFullPackageName(){ return TITLE+" "+VERSION;}
    
    public static final String HDR_PREFIX = "#| ";
    
    public static String getStandardHeader(Class<?> C)
    {
        return getStandardHeader(getFullPackageName()+"::"+C.getName());
    }

    public static String getStandardHeader(String info)
    {
        return HDR_PREFIX+info;
    }
    
    public static String getStandardRuntimeInfo()
    {
        java.util.Date now = java.util.Calendar.getInstance().getTime();
        
        java.util.Properties Props=System.getProperties();

        String system = Props.getProperty("os.name", "[unknown OS]")+" "+Props.getProperty("os.version","[Unknown version]")+" "+Props.getProperty("os.arch","[Unknown architecture]");
        String java = Props.getProperty("java.vm.name","[Unknown VM]")+" "+Props.getProperty("java.vm.version","[Unknown version]")+" ("+Props.getProperty("java.runtime.version","[Unknwon runtime]")+") "+Props.getProperty("java.vm.info","")+", "+Props.getProperty("java.vm.vendor","[Uknown vendor]");
        String cwd = Props.getProperty("user.dir","[unknown]");

        StringBuilder message = new StringBuilder();
        message.append(getStandardHeader("System: ")).append(system);
        message.append("\n").append(getStandardHeader("Java engine:")).append(java);
        message.append("\n").append(getStandardHeader("Date: ")).append(now);
        message.append("\n").append(getStandardHeader("Current directory: ")).append(cwd);
        return message.toString();
    }
    
}