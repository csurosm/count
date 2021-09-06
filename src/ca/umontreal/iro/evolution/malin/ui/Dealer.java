/*
 * Dealer.java
 *
 * Created on November 13, 2007, 12:04 PM
 */

package ca.umontreal.iro.evolution.malin.ui;


//import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.InvocationHandler;
import java.lang.reflect.Method;
import java.lang.reflect.Proxy;

import java.net.URL;
import java.net.URLClassLoader;

//import javax.swing.Box;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.Icon;
import javax.swing.JCheckBox;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;


import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Toolkit;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;

import java.io.File;
import java.io.PrintWriter;
import java.io.StringWriter;

import java.util.Vector;

import ca.umontreal.iro.evolution.Parser;
import ca.umontreal.iro.evolution.TreeNode;


/**
 *  Main executable class: supposed to be extended for different applications.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public abstract class Dealer implements 
        InvocationHandler // for Mac application events
            
{
    protected static String VERSION_NUMBER = "0.0002";
    
    /**
     * Initialize with a frame.
     * 
     * @param frame root container
     */
    public Dealer(JFrame frame)
    {
        top_frame = frame;
        //frame.setUndecorated(true);
        init();
    }
    
    /**
     * Frame in which the execution is done.
     */
    protected JFrame top_frame;
        
    protected CardLayout layout;
    
    protected JPanel main_panel;
    

    protected Vector<WorkSpace> work_space_list;
    protected int active_work_space=-1;
    
    public JFrame getTopFrame()
    {
        return top_frame;
    }

    protected abstract String getExecutableTitle();

    /**
     * Version number displayed in error messages.
     * 
     * @return a string that is appended to 'v' for display
     */
    protected String getVersionNumber(){ return VERSION_NUMBER;}
    
    
    
    /**
     * Called after top-level container is set.
     */
    private void init()
    {
        layout = new CardLayout();
        main_panel=new JPanel(layout);
        main_panel.setOpaque(true);        
        top_frame.setContentPane(main_panel);
        
        work_space_list = new Vector<WorkSpace>();
        setTitle(getExecutableTitle()); // +" v"+VERSION_NUMBER);

        initApplication();

        initMenu();
        
    }
    
    public void add(java.awt.Component C, String key)
    {
        main_panel.add(C,key);
        show(key);
    }
    
    /** 
     * Brings the work space with the given key in the foreground
     * 
     * @param component_key title for the work space
     */
    public void show(String component_key)
    {
        layout.show(main_panel,component_key);
        top_frame.setTitle(component_key);
    }

    public void addWorkSpace(WorkSpace work_space)
    {
        work_space_list.add(work_space);

        String key =  work_space.getAssociatedFile().getName();//Integer.toString(work_space_list.size()); // work_space.getAssociatedFile().getName()
        //main_panel.invalidate();
        this.add(work_space, key);
        
        active_work_space = work_space_list.size()-1;  
        initMenu();
        main_panel.validate();
        //System.out.println("#**D.aWS frame "+top_frame);
        //System.out.println("#**D.aWS panel "+main_panel);
        //System.out.println("#**D.aWS layout "+layout);
        //System.out.println("#**D.aWS workspace "+work_space);
    }    
    
    public void removeWorkSpace(WorkSpace work_space)
    {
        work_space_list.remove(work_space);
        main_panel.remove(work_space);
        //layout.removeLayoutComponent(work_space);
        active_work_space = work_space_list.size()-1;  
        if (active_work_space == -1)
            setTitle(getExecutableTitle()); // +" v"+VERSION_NUMBER);
        else
        {
            String key =  getActiveWorkSpace().getAssociatedFile().getName();//Integer.toString(work_space_list.size()); // work_space.getAssociatedFile().getName()
            show(key);
        }
        
        initMenu();

        main_panel.validate();
    }
    
    /**
     * Menu bar at the top-level container
     */
    protected JMenuBar top_menu_bar;
    
    protected JMenu session_menu;
    protected JMenu data_menu;
    protected JMenu rate_menu;
    protected JMenu analysis_menu;
    
    protected JMenu window_menu;
    protected JMenu help_menu;
    
    /**
     * Initializes the top menu.
     */
    protected void initMenu()
    {
        if (session_menu==null)
            session_menu = new JMenu("Session");
        
        setSessionMenu();
        
        if (data_menu==null)
            data_menu = new JMenu("Data");
        
        setDataMenu();
        
        if (rate_menu==null)
            rate_menu = new JMenu("Rates");
        setRateMenu();
            
        if (analysis_menu==null)
            analysis_menu = new JMenu("Analysis");
        
        setAnalysisMenu();
        
        if (help_menu==null)
            help_menu = new JMenu("Help");
        
        setHelpMenu();

        if (top_menu_bar == null)
        {
            top_menu_bar = new JMenuBar();
            
            top_menu_bar.add(session_menu);
            
            top_menu_bar.add(data_menu);
            top_menu_bar.add(rate_menu);        
            top_menu_bar.add(analysis_menu);
            top_menu_bar.add(help_menu);

            top_frame.setJMenuBar(top_menu_bar);
        }

        
        //window_menu = new JMenu("Window");
        //setWindowMenu();
        //top_menu_bar.add(window_menu);
        
        
        //top_menu_bar.setHelpMenu(new JMenu("Help is on the way"));
        
        
    }

    public WorkSpace getActiveWorkSpace()
    {
        if (active_work_space == -1)
            return null;
        else
            return (WorkSpace)work_space_list.get(active_work_space);
    }
    
//    protected void setWindowMenu()
//    {
//        window_menu.removeAll();
//        window_menu.add(new JMenuItem("Close"));
//        window_menu.add(new JMenuItem("Minimize"));
//        window_menu.addSeparator();
//        if (!work_space_list.isEmpty())
//        {
//            for (int j=0; j<work_space_list.size(); j++)
//            {
//                WorkSpace work_space = (WorkSpace) work_space_list.get(j);
//                JMenuItem item = new JMenuItem(work_space.getAssociatedFile().getName());
//                item.addActionListener(new windowActivator(j));
//                window_menu.add(item);
//            }
//        }
//    }
    
    private class windowActivator implements ActionListener
    {
        private int window_idx;
        windowActivator(int window_idx)
        {
            this.window_idx=window_idx;
        }
        
        @Override
        public void actionPerformed(ActionEvent e)
        {
            WorkSpace ws = (WorkSpace) work_space_list.get(window_idx);
            String key = ws.getAssociatedFile().getName();
            //System.out.println("activate "+ key);
            show(key);//ws.getAssociatedFile().getName());
            active_work_space = window_idx;
            setSessionMenu();
            setDataMenu();
            setRateMenu();
            setAnalysisMenu();
        }
    }
    
    /**
     * Used by addSessionMenuItems to check if the menu items are all there 
     */
    protected int session_menu_selection_count=-1;
    /**
     * Called by the session menu: after this, the Quit 
     * menu item will be added on non-Mac executions.
     *  
     */
    protected void addSessionMenuItems()
    {
        if (session_menu_selection_count != work_space_list.size())
        {
            JMenuItem load_phylo = new JMenuItem("Start new session ...");
            load_phylo.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {loadPhylogeny();}
                });
            load_phylo.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, 
                Toolkit.getDefaultToolkit().getMenuShortcutKeyMask())); 

            session_menu.add(load_phylo);
        
            if (!work_space_list.isEmpty())
            {
                session_menu.addSeparator();
                for (int j=0; j<work_space_list.size(); j++)
                {
                    WorkSpace work_space = (WorkSpace) work_space_list.get(j);
                    JMenuItem item = null;
                    if (j==active_work_space) 
                        item = new JMenuItem(work_space.getAssociatedFile().getName(), new SelectedSessionIcon());
                    else 
                        item = new JMenuItem(work_space.getAssociatedFile().getName());
                    item.addActionListener(new windowActivator(j));
                    session_menu.add(item);
                }
            }
            session_menu_selection_count = work_space_list.size();
        
            session_menu.addSeparator();
            JMenuItem close_session = new JMenuItem("Close session");
            close_session.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {
                        removeWorkSpace(getActiveWorkSpace());
                    }
                });
            close_session.setEnabled(getActiveWorkSpace() != null);
            session_menu.add(close_session);        
        }
    }
    
    private class SelectedSessionIcon implements Icon
    {
        @Override
        public int getIconWidth(){return 30;}
        @Override
        public int getIconHeight(){return 20;}
        @Override
        public void paintIcon(Component c, Graphics g, int x, int y) // x y is top left 
        {
            String icon_string = "\u261e";//"\u2717";//
            Graphics myg = g.create();
            myg.setFont(myg.getFont().deriveFont(24f).deriveFont(Font.BOLD));
            myg.setColor(Color.BLACK);
            int w = myg.getFontMetrics().stringWidth(icon_string);
            myg.drawString(icon_string, x+15-w/2, y+20);
        }
    }
    
    /**
     * Sets up the items in the <q>Session</q> menu.
     */
    protected void setSessionMenu()
    {
        if (session_menu_selection_count != work_space_list.size())
        {
                
            session_menu.removeAll();
            addSessionMenuItems();

            if (Mac_application == null)
            {
                JMenuItem item = new JMenuItem("Quit");
                item.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent E)
                    {
                        doQuit();
                    }
                });
                item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Q, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
                session_menu.add(item);
            }
            
        }
    }
    
    protected abstract void setDataMenu();
    
    protected abstract void setRateMenu();
    
    protected abstract void setAnalysisMenu();
    
    /**
     * The default implementation is that an About item is added when necessary.
     */
    protected void setHelpMenu()
    {
        if (help_menu.getItemCount()==0)
        {
            if (Mac_application == null)
            {
                JMenuItem about_item = new JMenuItem(ABOUT_MENU_TEXT);
                about_item.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent E)
                    {
                        doAbout();
                    }
                });
                help_menu.add(about_item);
            }
        }
    }
        
    protected String ABOUT_MENU_TEXT = "About";
    
    /**
     * Loads the file with Newick-format phylogeny.
     * Opens a new session.
     */
    protected void loadPhylogeny() 
    {
        FileDialog dialog = null;
        dialog = new FileDialog(top_frame,"Open Newick-format phylogeny file",FileDialog.LOAD);
        dialog.setVisible(true);
       
        String file_name = dialog.getFile();
        String directory = dialog.getDirectory();
        TreeNode main_phylo = null;

        if (file_name != null)
            for (int i=0; i<work_space_list.size(); i++)
                if (file_name.equals(work_space_list.get(i).getAssociatedFile().getName()))
                {
                    exceptionCaught(new RuntimeException("Duplicate file names"), "Duplicate names", "Cannot have two open sessions with the same phylogeny file names");
                    file_name = null;
                    break;
                } 
        
        if (file_name != null)
        {
            try 
            {
                java.io.InputStreamReader R 
                = new java.io.InputStreamReader(new java.io.FileInputStream(directory+file_name));
                main_phylo = Parser.readNewick(R,false,true,false);
            } catch (java.io.InterruptedIOException E)
            {
                // load canceled
                main_phylo = null;
            } catch (java.io.IOException E)
            {
                exceptionCaught(E, "I/O error", "File error while attempting ro read Newick file "+directory+file_name+".");  
            } catch (Parser.ParseException E)
            {
                exceptionCaught(E, "Newick file parsing error", "Parsing error while reading Newick file "+directory+file_name+".");
            }
            if (main_phylo != null)
            {
                TreeNode[] dft=main_phylo.getTraversal().getDFT();
                boolean all_zero = true;
                for (int j=0; j<dft.length-1; j++)
                    if (dft[j].getLength()!=0.)
                    {
                        all_zero=false;
                        break;
                    }
                if (all_zero)
                {
                    for (int j=0; j<dft.length-1; j++) // don't change root
                        dft[j].setLength(1.0);
                }
                WorkSpace ws = newWorkSpace(main_phylo, new File(directory,file_name));
                
                String ws_title = file_name; 
                int ri = file_name.lastIndexOf('.');
                if (ri != -1)
                    ws_title = file_name.substring(0,ri);
                try {
                    setTitle(ws_title);
                    addWorkSpace(ws);
                    //ws.setVisible(true);
                    //System.out.println("#**D.lP add "+ws_title+" / "+ws);
                } catch (Exception E)
                {
                    exceptionCaught(E, "You found a programming bug!", "Error while computing layout of phylogeny");
                }
            }
            setSessionMenu();
            setDataMenu();
            setRateMenu();
            setAnalysisMenu();
            top_frame.validate();
        }
    }
    
    private class ReportUncaughtException implements Thread.UncaughtExceptionHandler
    {
        @Override
        public void uncaughtException(Thread thread, final Throwable T)
        {
            SwingUtilities.invokeLater(new Runnable()
            {
                @Override
                public void run()
                {
                    exceptionCaught(T, "A baffling error...");
                }
            });
        }
    }
    
    protected abstract WorkSpace newWorkSpace(TreeNode root, File F);

    /**
     * Turns this guy into the uncaught exception handler for the thread
     * @param T a thread
     */
    public void makeThreadExceptionHandler(Thread T)
    {
        T.setUncaughtExceptionHandler(new ReportUncaughtException());
    }
    
    /**
     * Popup dialog for generic error message. This should be called from within the Event Dispatch Thread.
     * 
     * @param E the exception that is reported here
     * @param title title for the dialog window
     * 
     */
    public void exceptionCaught(Throwable E, String title)
    {
        exceptionCaught(E, title, null);
    }
    
    /**
     * Popup dialog for generic error message. This should be called from within the Event Dispatch Thread.
     * 
     * @param E the exception that is reported here
     * @param title title for the dialog window
     * @param more_info Additional information about the error
     */
    public void exceptionCaught(Throwable E, String title, String more_info)
    {
        StringWriter SW = new StringWriter();
        E.printStackTrace(new PrintWriter(SW));
        String stack_trace = SW.toString().trim();
        
        StringBuffer short_message = new StringBuffer();
        short_message.append("<p><b>Awfully sorry - something went wrong.</b></p>");
        if (more_info != null)
            short_message.append("<p>"+more_info+"</p>");
        short_message.append("<p>If this is surprising, " +
            "then please write to Mikl&oacute;s Cs&#369;r&ouml;s [csuros"+"@"+"iro.umontreal.ca] about it:<br />" +
            "click below to display the technical information, and copy the text into your e-mail.</p>");
        
        // system info
        java.util.Date now = java.util.Calendar.getInstance().getTime();
        java.util.Properties Props=System.getProperties();
        String system = Props.getProperty("os.name", "[unknown OS]")+" "+Props.getProperty("os.version","[Unknown version]")+" "+Props.getProperty("os.arch","[Unknown architecture]");
        String java = Props.getProperty("java.vm.name","[Unknown VM]")+" "+Props.getProperty("java.vm.version","[Unknown version]")+" ("+Props.getProperty("java.runtime.version","[Unknwon runtime]")+") "+Props.getProperty("java.vm.info","")+", "+Props.getProperty("java.vm.vendor","[Uknown vendor]");
        
        StringBuffer message = new StringBuffer();
        message.append("<p>Version: "+getExecutableTitle()+" v"+getVersionNumber());
        message.append("<br>System: "+system+"<br>Java engine: "+java);
        message.append("<br>UIManager: "+javax.swing.UIManager.getLookAndFeel().getName());
        message.append("</p>");
        
        message.append("<p><i>Stack trace:</i></p>");
        message.append("<pre>");
        message.append(SW.getBuffer());
        message.append("</pre>");
        
        JEditorPane shortEP = new JEditorPane("text/html", short_message.toString());
        shortEP.setEditable(false);
        shortEP.setBackground(new Color(255,192,192)); // Banana
        shortEP.setBorder(BorderFactory.createRaisedBevelBorder());
        
        final JEditorPane EP = new JEditorPane("text/html",message.toString());
        EP.setEditable(false);
        EP.setBorder(BorderFactory.createEtchedBorder());
        
        JPanel msg_display = new JPanel();
        BoxLayout msg_layout  = new javax.swing.BoxLayout(msg_display, BoxLayout.PAGE_AXIS);
        msg_display.setLayout(msg_layout);
        final JCheckBox show_details = new JCheckBox("Show technical information");
        show_details.addActionListener(new ActionListener()
            {
                @Override
                public void actionPerformed(ActionEvent E)
                {
                    EP.setVisible(show_details.isSelected());
                }
            });
        show_details.setSelected(false); 
        EP.setVisible(false);
        msg_display.add(shortEP);
        msg_display.add(show_details);
        msg_display.add(EP);
        
        JScrollPane msg_pane = new JScrollPane(msg_display);
        msg_pane.setPreferredSize(new Dimension(900, 600));
        msg_pane.setMinimumSize(new Dimension(80, 30));   
        
        Object[] button_text = {"OK", "Quit "+getExecutableTitle()};
        
        int ans = JOptionPane.showOptionDialog(top_frame,msg_pane,title,JOptionPane.OK_CANCEL_OPTION,JOptionPane.ERROR_MESSAGE,null,button_text,button_text[0]);
        if (ans==1)
            doQuit();
    }

    /**
     * Width of the top-level container 
     * @return correct width if frame, a constant 1000 if applet
     */ 
    public int getWidth()
    {
        return top_frame.getWidth();
    }

    /**
     * Height of the top-level container 
     *
     * @return correct width if frame, a constant 700 if applet
     */ 
    public int getHeight()
    {
        return top_frame.getHeight();
    }
    
    /**
     * Sets the title for the top-level container.
     * 
     * @param title new title for the root container
     * 
     */
    public void setTitle(String title)
    {
        top_frame.setTitle(title);
    }
    
    /**
     * Sets the visibility of the top-level container.
     * 
     * @param b whether container should be visible
     */
    public void setVisible(boolean b)
    {
        top_frame.setVisible(b);
    }
    
    
    // 
    // MacOS X integration --- code based on Eirik Bjorsnos org.simplericity.macify.eawt
    //
    
    protected Object Mac_application = null;
    
    public boolean isMac()
    {
        return Mac_application != null;
    }
    
    /**
     * Initializes Mac_application
     */
    private void initApplication()
    {
       try 
       {
            final File file = new File("/System/Library/Java");
            if(file.exists()) 
            {
                ClassLoader scl = ClassLoader.getSystemClassLoader();
                Class clc = scl.getClass();
                if(URLClassLoader.class.isAssignableFrom(clc)) 
                {
                    Method addUrl  = URLClassLoader.class.getDeclaredMethod("addURL", new Class[] {URL.class});
                    addUrl.setAccessible(true);
                    addUrl.invoke(scl, new Object[] {file.toURI().toURL()});
                }
            }

            Class<?> application_class = Class.forName("com.apple.eawt.Application");
            Mac_application = application_class.getMethod("getApplication", new Class[0]).invoke(null, new Object[0]);
            Class application_listener_class = Class.forName("com.apple.eawt.ApplicationListener");
            Object listener = Proxy.newProxyInstance(getClass().getClassLoader(), 
                    new Class[]{application_listener_class},
                    this);
            callMethod("addApplicationListener",new Class[]{application_listener_class}, new Object[]{listener});
        } catch (ClassNotFoundException e) 
        {
            Mac_application = null;
        } catch (Exception e) 
        {
            exceptionCaught(e,"A bug?","There was a problem with the initialization of the application.");
        }
    }
    
    /**
     * Called when user selects the About menu point
     */
    protected void doAbout()
    {
        //System.out.println("#**D.dA ");
    }
    
    /**
     * Called when application is about to quit
     * Default implementation calls System.exit(0).
     *
     */
    protected void doQuit()
    {
        System.exit(0);
    }
    
    
    private Object callMethod(String methodname) 
    {
        return callMethod(methodname, new Class[0], new Object[0]);
    }

    private Object callMethod(String methodname, Class[] classes, Object[] arguments) 
    {
        try 
        {
            if (classes == null) 
            {
                classes = new Class[arguments.length];
                for (int i = 0; i < classes.length; i++) 
                {
                    classes[i] = arguments[i].getClass();

                }
            }
            Method method = Mac_application.getClass().getMethod(methodname, classes);
            return method.invoke(Mac_application, arguments);
        } catch (Exception E) 
        {
            exceptionCaught(E, "A bug?", "Problem when trying to invoke method "+methodname+" for "+Mac_application);
            return null;
        }
    }
    
    /**
     * Called with application events
     */
    public Object invoke(Object proxy,
                     Method apple_method,
                     Object[] args)
              throws Throwable   
    {
        try 
        {
            String method_name = apple_method.getName();
            if ("handleAbout".equals(method_name))
            {
                doAbout();
                // make sure no other window pops up ...
                Class<?> event_class = Class.forName("com.apple.eawt.ApplicationEvent");
                Method set_handled = event_class.getMethod("setHandled", new Class[]{boolean.class});
                set_handled.invoke(args[0],new Object[]{Boolean.TRUE});
            }
            else if ("handleQuit".equals(method_name))
            {
                    doQuit();
            } else
            {
                //System.out.println("#**D.invoke "+method_name);
            }
        } 
        catch (Exception T)
        {
            exceptionCaught(T, "A bug?", "An error occurred while trying to handle an Apple event");
        }
        return null;
    }


    
    
}
