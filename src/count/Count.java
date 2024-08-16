package count;

/*
 * Copyright 2021-2023 Mikl&oacute;s Cs&#369;r&ouml;s.
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


import java.io.PrintStream;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import java.util.concurrent.ForkJoinPool;

import count.gui.AppFrame;
import count.gui.Session;
import count.io.CommandLine;


/**
 * GUI application launcher (Swing).
 * 
 * @author Miklós Csűrös
 *
 */
public class Count 
{
	/**
	 * Official app name.
	 */
    public static final String APP_TITLE = "Count XXIV";
    /**
     * Version string starting with a float value (integer part=major version); 
     * set to the date during development.
     */
    public static final String APP_VERSION = "24.0815::BundleTree.Node.addTree";
    
    /**
     * Name of the application.
     * 
     * @return app title
     */
    private final static String getAppTitle()
    {
    	return APP_TITLE;
    }
    
    /**
     * Version info for application.
     * 
     * @return app version string
     */
    private final static String getAppVersion()
    {
    	return APP_VERSION;
    }
    
    /**
     * Numerical part of app version.
     * 
     * @return
     */
    public final static double getAppVersionNumber()
    {
    	int dot = -1;
    	int digits = 0;
    	String version = APP_VERSION;
    	while (digits<version.length())
    	{
    		char c = version.charAt(digits);
    		if (Character.isDigit(c))
    		{
    			++digits;
    		} else if (c=='.')
    		{
    			if (0<=dot) break; // already had a dot
    			dot = digits;
    			++digits;
    		} else
    			break;
    	}
    	String prefix = version.substring(0,digits);
    	double getAppVersionNumber = Double.parseDouble(prefix);
    	return getAppVersionNumber;
    }
    
    public final static String getAppFullName()
    { return getAppTitle()+" "+getAppVersion();}
    
    
	/**
	 * Number of threads allocated in thread pools
	 * ; set to 1 for no multithreading. 
	 * (Different classes use their own static 
	 * thread pools.)  
	 * Default value is 0.75*<var>vCPU</var> rounded up; so if <var>vCPU</var>&lt;4, 
	 * then equals <var>vCPU</var>.
	 */
	public static int THREAD_PARALLELISM = 
			(3*Runtime.getRuntime().availableProcessors()+3)/4; // 
	/**
	 * Number of families handled by one thread
	 * at a time. Set to {@link Integer#MAX_VALUE} to 
	 * calculate on a single thread. 
	 */
	public static int THREAD_UNIT_TASK = 64; //Integer.MAX_VALUE; //64; // so many families are processed by one thread in a unit task 
	
	
	/**
	 * Empty interface to indicate that a class is multithreaded
	 * 
	 */
	public static interface UsesThreadpool{}
	
	/**
	 * Instantiates a new thread pool based on the 
	 * current value of {@link #THREAD_PARALLELISM} 
	 * @return null if {@link #THREAD_PARALLELISM} is not greater than 1
	 */
	public static ForkJoinPool threadPool()
	{
		ForkJoinPool pool = null;
		if (1<THREAD_PARALLELISM)
		{
			pool = new ForkJoinPool(THREAD_PARALLELISM);	
//			System.out.println("#**C.threadPool init: "+THREAD_PARALLELISM+" threads for "+pool+"\ton thread "+Thread.currentThread());
		}
		return pool;
	}
	
	public static int unitTask(int task_count)
	{
		return unitTask(task_count, THREAD_UNIT_TASK);
	}
	
	public static int unitTask(int task_count, int max_unit_task)
	{
		if (max_unit_task<=0) throw new IllegalArgumentException("Count.unitTask: max_unit_task must be positive (got "+max_unit_task+")");
		if (max_unit_task == Integer.MAX_VALUE || THREAD_PARALLELISM==1)
			return task_count;

		int avg_load  = (task_count+THREAD_PARALLELISM-1)/THREAD_PARALLELISM;
		int unitTask = Integer.min(avg_load, max_unit_task); 
//		System.out.println("#**C.unitTask "+task_count+"/"+THREAD_PARALLELISM+" threads: avg "+avg_load+"\tmax "+max_unit_task+"\tset "+unitTask);
		return unitTask;
		// want u*	THREAD_PARALLELISM <= task_count
		// u = task_count/THREAD_PARALLELISM;
	}
	
	public static PrintStream out = System.out;
	
    /**
     * Initialize from command line
     * 
     * @param args command-line arguments
     * @throws Exception if there is a problem with the command-line arguments
     */
    public Count(String[] args) throws Exception
    {
    	this.cli = new CommandLine(args, getClass(), 0);
    }
    /**
     * Frame in which the execution is done.
     */
    private AppFrame top_frame;
	
    private final CommandLine cli;
    
    /**
     * Launches the Swing application; call from event queue
     */
    private void createAndShowGUI() 
    {
        top_frame = new AppFrame(this);
        count.gui.UncaughtExceptionHandler handler = top_frame.getExceptionHandler();
        Thread.currentThread().setUncaughtExceptionHandler(handler);
        	
        if (cli.getTree()==null) 
        {
	        top_frame.doAbout(true);
        } else
        {
        	Session sesh = new Session(top_frame, cli.getTreeData());
        	top_frame.addSession(sesh);
        	if (cli.getTableData() != null)
        	{
        		sesh.addDataSet(cli.getTableData());
        	}
        	if (cli.getMixedRateModelData()!= null)
        	{
        		sesh.addRates(cli.getMixedRateModelData(), true);
        	}
        }
        top_frame.setVisible(true);
    }
    
    
    
    /**
     * Called when user selects the Quit menu point, or programmatically. 
     * 
     * Default implementation calls System.exit(0).
     */
    public void doQuit()
    {
        //System.out.println("#**This is doQuit!");
        System.exit(0);
    }
    
    private static void setMacLook()
    {
        try
        {
            String system = System.getProperty("os.name", "[unknown OS]");

            if (system.startsWith("Mac OS"))
            {
                // MacOS X 10.4+ --- Java runtime properties; no effect on other systems 
                System.setProperty("apple.laf.useScreenMenuBar", "true");
                System.setProperty("com.apple.mrj.application.apple.menu.about.name", APP_TITLE);
                System.setProperty("apple.awt.showGrowBox","true");
                //System.setProperty("apple.awt.brushMetalLook","true");
                
                
                // System.setProperty( "apple.awt.application.name", 
                // 
                System.setProperty( "apple.awt.application.appearance", "system" );
            } 
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception E) // ClassNotFoundException, InstantiationException
        {
            System.err.println("Ooops - no Look&Feel");
            E.printStackTrace();
        }
    }
    
    public static boolean isMacOS()
    {
        String system = System.getProperty("os.name", "[unknown OS]");

        return (system.startsWith("Mac OS"));
    }
    
	public static void main(String[] args) throws Exception
	{
		setMacLook();
		
		Count app  = new Count(args);
		SwingUtilities.invokeLater(()->app.createAndShowGUI());
    }
    
}

