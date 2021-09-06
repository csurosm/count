/*
 * Copyright 2021 Mikl&oacute;s Cs&#369;r&ouml;s.
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
package count.util;

/**
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 */
public class Executable 
{
    private Executable(){}
    
    public static final String TITLE = "Count";
    public static final String VERSION = "XXII 21.0727(dev)";
    
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