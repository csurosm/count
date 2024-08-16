package count.gui;
/*
 * Copyright 2022 Mikl&oacute;s Cs&#369;r&ouml;s.
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

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.PrintWriter;
import java.io.StringWriter;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JEditorPane;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

import static count.Count.APP_TITLE;
import static count.Count.APP_VERSION;

/**
 * Class for handling exception via a popup dialog that gives the details. 
 * 
 */
public class UncaughtExceptionHandler implements Thread.UncaughtExceptionHandler
{
	
	UncaughtExceptionHandler(AppFrame app)
	{
		this.app = app;
	}
	
	private final AppFrame app;
	
    /**
     * Internal counter for number of exceptions caught. 
     */
    private int num_caught_exceptions = 0;
    
    /**
     * Cutoff for too many exceptions. 
     * When multiple exceptions are handled and the pop-up 
     * windows start to pile up (counted by {@link #num_caught_exceptions}, the executable simply quits. 
     */
    private static final int TOO_MANY_EXCEPTIONS = 10;
    

    
    private static final String anyError()
    {
    	String[] synonyms = 
            {   "puzzling", 
                    "bewildering", 
                    "perplexing", 
                    "mystifying", 
                    "baffling", 
                    "mysterious", 
                    "peculiar", 
                    "curious", 
                    "bizarre", 
                    "strange", 
                    "weird",
                    "unfathomable",
                    "abstruse",
                    "enigmatic",
                    "unanticipated",
                    "unforeseen",
                    "unexpected",
                    "startling"
        };
        java.util.Random RND = new java.util.Random();
        String error_attribute = synonyms[RND.nextInt(synonyms.length)];
        String article = "aeiou".indexOf(error_attribute.charAt(0))==-1?"A ":"An ";
        return article+error_attribute+" error...";    	
    }

    
    
    /**
     * Calls {@link #handle(Throwable, String)} with 
     * this error.
     * 
     * @param thread where the exception was thrown
     * @param T the exception that was thrown
     * 
     */
    @Override
    public void uncaughtException(Thread thread, final Throwable T)
    {
        SwingUtilities.invokeLater(new Runnable()
        {
            @Override
            public void run()
            {
//                String[] synonyms = 
//                {   "puzzling", 
//                    "bewildering", 
//                    "perplexing", 
//                    "mystifying", 
//                    "baffling", 
//                    "mysterious", 
//                    "peculiar", 
//                    "curious", 
//                    "bizarre", 
//                    "strange", 
//                    "weird",
//                    "unfathomable",
//                    "abstruse",
//                    "enigmatic",
//                    "unanticipated",
//                    "unforeseen",
//                    "unexpected",
//                    "startling"};
//                java.util.Random RND = new java.util.Random();
//                String error_attribute = synonyms[RND.nextInt(synonyms.length)];
//                String article = "aeiou".indexOf(error_attribute.charAt(0))==-1?"A ":"An ";
//                handle(T, article+error_attribute+" error...");
            	handle(T);
            }
        });
    }
    
    /**
     * Pop-up dialog for generic error message. 
     * This should be called from within the Event Dispatch Thread.
     * 
     * @param E the exception that is reported here
     * @param title title for the dialog window
     * 
     */
    public void handle(Throwable E, String title)
    {
        handle(E, title, null);
    }
    
    /**
     * Pop-up dialog for generic error message.
     */     
    public void handle(Throwable T)
    {
    	handle(T, anyError());
    }
    
    /**
     * Pop-up dialog for generic error message. 
     * This should be called from within the Event Dispatch Thread.
     * 
     * @param E the exception that is reported here
     * @param title title for the dialog window
     * @param more_info Additional information about the error
     */
    public void handle(Throwable E, String title, String more_info)
    {
        synchronized(this)
        {
            num_caught_exceptions++;
            if (num_caught_exceptions==TOO_MANY_EXCEPTIONS)
            {
                System.err.println("TOO MANY ERRORS in too many threads: "+E);
                E.printStackTrace(System.err);
                System.exit(2022);
            }
        }

        StringWriter SW = new StringWriter();
        PrintWriter PW = new PrintWriter(SW);
        E.printStackTrace(PW);
        
//        while ((E=E.getCause())!=null) // prinStackTrace inclues the cause
//        {
//        	PW.println(".. caused by ");
//        	E.printStackTrace(PW);
//        }
        
        //String stack_trace = SW.toString().trim();
        
        StringBuilder short_message = new StringBuilder();
        short_message.append("<p><b>Awfully sorry - something went astray.</b></p>");
        if (more_info != null)
            short_message.append("<p>").append(more_info).append("</p>");
        short_message.append("<p>If this is surprising, " +
            "then please write to Mikl&oacute;s Cs&#369;r&ouml;s [csurosm@gmail.com] about it:<br />" +
            "click below to display the technical information, and copy the text into your e-mail.</p>");
        
        // system info
        java.util.Date now = java.util.Calendar.getInstance().getTime();
        java.util.Properties Props=System.getProperties();
        String system = Props.getProperty("os.name", "[unknown OS]")+" "+Props.getProperty("os.version","[Unknown version]")+" "+Props.getProperty("os.arch","[Unknown architecture]");
        String java = Props.getProperty("java.vm.name","[Unknown VM]")+" "+Props.getProperty("java.vm.version","[Unknown version]")+" ("+Props.getProperty("java.runtime.version","[Unknwon runtime]")+") "+Props.getProperty("java.vm.info","")+", "+Props.getProperty("java.vm.vendor","[Uknown vendor]");
        
        StringBuilder message = new StringBuilder();
        message.append("<p>Version: ").append(APP_TITLE).append(" v").append(APP_VERSION);
        message.append("<br>System: ").append(system);
        message.append("<br>Java engine: ").append(java);
        message.append("<br>UIManager: ").append(javax.swing.UIManager.getLookAndFeel().getName());
        message.append("<br>Thread: ").append(Thread.currentThread()); // EDT unless coding mistake
        message.append("<br>Date: ").append(now);
        message.append("</p>");
        
        message.append("<p><i>Stack trace:</i></p>");
        message.append("<pre>");
        message.append(SW.getBuffer());
        message.append("</pre>");
        
        JEditorPane shortEP = new JEditorPane("text/html", short_message.toString());
        shortEP.setEditable(false);
        shortEP.setBackground(app.WARNING_COLOR);
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
                public void actionPerformed(ActionEvent ignored)
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
        
        Object[] button_text = {"So be it", "Quit "+APP_TITLE};
        
        int ans = JOptionPane.showOptionDialog(app,msg_pane,title,
                JOptionPane.OK_CANCEL_OPTION,JOptionPane.ERROR_MESSAGE,
                null,button_text,button_text[0]);

        //System.out.println("#*"+getClass().getName()+".handleException "+EP.getText()+"\n// "+Thread.currentThread());

        if (ans==1)
            app.doQuit();
    }    
}