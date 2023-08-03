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

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.KeyStroke;
import javax.swing.ProgressMonitorInputStream;
import javax.swing.SwingWorker;

import count.Count;
import count.ds.AnnotatedTable;
import count.ds.Phylogeny;
import count.gui.kit.CountActions;
import count.gui.kit.FancyFileDialog;
import count.gui.kit.StringIcon;
import count.io.DataFile;
import count.io.GeneralizedFileReader;
import count.io.ModelBundle;
import count.io.NewickParser;
import count.io.RateVariationParser;
import count.io.TableParser;
import count.model.GammaInvariant;

import static count.Count.APP_TITLE;
import static count.Count.APP_VERSION;

/**
 * Application frame; the main component of the GUI.  
 * 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */
public class AppFrame extends JFrame
{
    public static final Color WARNING_COLOR = new Color(255, 255, 192); // Banana
    public static final Color SMOKY_BACKGROUND = new Color(120,120,200,50);

    
    
    public AppFrame(Count app)
	{
		super();
		this.app = app;
        this.session_list = new ArrayList<>();
        this.active_session = -1;
        
        this.exception_handler = new UncaughtExceptionHandler(this);
        Thread.setDefaultUncaughtExceptionHandler(exception_handler);
//        Thread.currentThread().setUncaughtExceptionHandler(exception_handler);

        layout = new CardLayout();
        main_panel = new JPanel(layout);
        main_panel.setOpaque(true);  
        
// TOOLBAR	
//      this.setContentPane(main_panel);	
        this.tool_bar = createToolBar();
        JPanel panel_with_toolbar = new JPanel(new BorderLayout());
        panel_with_toolbar.add(tool_bar, BorderLayout.PAGE_START);
        panel_with_toolbar.add(main_panel, BorderLayout.CENTER);
        this.setContentPane(panel_with_toolbar);
// TOOLBAR	
        
        
        menu_bar = new MenuBar();
        this.setJMenuBar(menu_bar);
        
        setDesktop();
        
        setIconImage(getCountIcon().getImage());        

        setBounds(25,25,
            (int)Toolkit.getDefaultToolkit().getScreenSize().getWidth()-50,
            (int)Toolkit.getDefaultToolkit().getScreenSize().getHeight()-150);
        
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        addRootPanel();
        
    }
    
    private static final String ROOT_KEY = "\n__root__\n";
	
	private final Count app;

    /**
     * Main panel (content pane for top frame). 
     */
    private final JPanel main_panel;
    
    /**
     * Main panel layout.
     */
    private final CardLayout layout;
    
    private final MenuBar menu_bar;
    
    private final CountTools tool_bar;
    
    /**
     * Exception handler.
     */
    private final UncaughtExceptionHandler exception_handler;
    
    
    private final List<Session> session_list;
    private int active_session;
    
    public UncaughtExceptionHandler getExceptionHandler()
    {
    	return exception_handler;
    }
    
    /**
     * Sets the Quit and About handlers.
     */
    private void setDesktop()
    {
    	Desktop desktop = Desktop.getDesktop();
    	desktop.setQuitHandler((event,response)->doQuit());
    	desktop.setAboutHandler(event->doAbout(false));
    }
    
    public Count getApp()
    {
    	return app;
    }
    
    public Session getActiveSession()
    {
    	if (active_session == -1)
    		return null;
    	else
    		return session_list.get(active_session);
    }
    
    public int getSessionCount()
    {
    	return session_list.size();
    }
    
    /** 
     * Brings the session with the given key in the foreground
     * and sets the {@link #active_session}. 
     * 
     * @param session_idx index of the workspace to select
     */
    public void showSession (int session_idx)
    {
        if (session_idx == -1)
        {
                setTitle(APP_TITLE); // +" v"+VERSION_NUMBER);
                layout.show(main_panel, ROOT_KEY);
        }
        else
        {
            Session sesh = session_list.get(session_idx);
            String component_key = sesh.getSessionKey();
            layout.show(main_panel,component_key);
            setTitle(sesh.getSessionTitle());
        }
        this.active_session = session_idx;
    	setMenusAndTitle();
    	tool_bar.updateEnabledStates();
    }
    
    public void removeSession(Session session)
    {
    	int confirm_remove = JOptionPane.showConfirmDialog(this
    			, "Please confirm: do you really want to remove this session? The removal cannot be undone."
    			, "Removing "+session.getSessionKey()
    			, JOptionPane.YES_NO_OPTION
    			, JOptionPane.WARNING_MESSAGE);
    	if (confirm_remove == JOptionPane.YES_OPTION)		
    	{
	    	session_list.remove(session);
	    	main_panel.remove(session);
	    	showSession(session_list.size()-1);
    	}
    }
    
    /**
     * The root panel is displayed when no 
     * sessions are available: gives some basic instructions.
     * 
     */
    private void addRootPanel()
    {
    	Dimension D = this.getSize();
    	String instructions = "<h2><span style=\"color:rgb(128,0,0)\">Steps for using Count</span></h2><ol>"
//    			+"<h4>Sessions</h4>"
    			+"<li> In Count, you work with <span style=\"color:rgb(128,0,0)\"><b>session</b>s</span>, "
    			+ "defined by a set of terminal taxon names "
    			+ "and one more phylogenies over the same leaf set. "
    			+ "As a first step, "
    			+ "you can either"
    			+ "<ul>"
    			+ "<li>start a <b>new session</b> by loading a phylogeny in Newick format, which will be the main "
    			+ "tree associated with the session; or</li>"
    			+ "<li>load a table in a <b>no-tree session</b> or <b>import a table</b> of family-sizes (copy-numbers), and "
    			+ "assume a simple initial tree (star tree, random tree, or UPGMA calculated from copy number "
    			+ "profiles); or</li>"
    			+ "<li><b>load a session<</b> previously saved by Count (XML format)</li>."
    			+ "</ul>"
    			+ "</li>"
//    			+ "<h4>Data tables</h4>"
    			+ "<li>Count works with a family size <span style=\"color:rgb(128,0,0)\"><b>table</b></span>: "
    			+ "a set of genomes described "
    			+ "as a multiset of families their genes belong to. "
    			+ "Such a table is tab-delimited text file that "
    			+ "starts with a header line that "
    			+ "lists the taxon names as column headers. "
    			+ "The first column has the family names. "
    			+ "Additional columns can be loaded as family annotation columns. "
    			+ "You can either <b>load a table</b> in that format, or "
    			+ "<b>import a table</b>. from COG-style csv files or MCL clustering output. "
    			+ "You can derive further tables by <b>filtering the rows</b> of "
    			+ "an existing table (condition on cell values in a column), "
    			+ "or by converting it to a <b>binary</b> presence-absence table.</li>"
    			+ ""
//    			+ "<h4>Phylogenies</h4>"
    			+ "<li>Every session has a main phylogeny, and can include other rooted "
    			+ "trees over the same leaf set. You can either <b>load a tree</b> "
    			+ "in Newick format, or <b>build a tree</b>, or modify one by hand, using the <b>edit tree</b> "
    			+ "function. (Tres can be built with a simple UPGMA procedure from "
    			+ "profiles, or set to random or star topology. "
    			+ "Edit operations include rerooting, edge-contraction, and "
    			+ "subtree-prune-and-regraft.) </li>"
    			+ ""
//    			+ "<h4>Rate models</h4>"
    			+ "<li>Rate <span style=\"color:rgb(128,0,0)\"><b>models</b></span> equip the selected phylogeny's edges with "
    			+ "a linear birth-and-death process defined by "
    			+ "rates and length: loss, duplication and gain rates. (For convenience, "
    			+ "edge lengths and rates are scaled so that loss rate=1.0. "
    			+ "Phylogenies with multifurcations are OK, but a "
    			+ "completely resolved model uses a binary tree.) "
    			+ "You can <b>load a rate model</b> (previously saved from Count), "
    			+ "or make one by <b>model optimization</b>. "
    			+ "(Models are optimized by maximizing their likelihood). "
    			+ "You should set the observation (aka <em>ascertainment</em>) bias for your input table, depending "
    			+ "on how the homolog families were constructed: "
    			+ "whether they have a minimum of 0 (if even the families without any members are known), "
    			+ "1 (if the table covers all genes from all genomes), "
    			+ "or 2 members (if families correspond to alignments)."
    			+ "</li>"
    			+ ""
//    			+ "<h4>Ancestral reconstructions</h4>"
    			+ "<li>Given a table, you can do <span style=\"color:rgb(128,0,0)\"><b>ancestral reconstructions</b></span>. "
    			+ "by either <b>numerical parsimony</b> (a generalized Wagner penalization), "
    			+ "or by <b>Dollo parsimony</b> if the data is binary. "
    			+ "If you have a rate model, then you can also reconstruct the "
    			+ "ancestral genomes by <b>posteriors</b>, giving probabilities and "
    			+ "expectations for family memberships; the reconstruction includes correction "
    			+ "for the observation bias of minimum copies.</li>"
    			+ ""
//    			+ "<h4>Saving your data</h4>"
    			+ "<li>You can <b>save</b> phylogenies, rate models, and tables, or "
    			+ "<b>export</b> ancestral reconstructions. "
    			+ "With the exception of phylogenies, "
    			+ "all other data is in tab-delimited files. "
    			+ "You can also <b>save all</b>  all your sessions "
    			+ "if you want to come back to them later."
    			+ "(The sessions are saved in an XML-format file.)"
    			+ "</li>"
    			+ "<li><b>Tooltips</b> give further instructions and informations "
    			+ "on the diplayed elements in the GUI. "
    			+ "</ol>";
    	JEditorPane instruction_pane = new JEditorPane("text/html", instructions);
    	//instruction_pane.setMaximumSize(D);
    	JPanel instruction_container = new JPanel(null);
    	instruction_pane.setBorder(BorderFactory.createLoweredSoftBevelBorder());
    	
    	instruction_container.add(instruction_pane);
    	
//    	D.width/=2;
//    	D.height/=2;
//    	instruction_pane.setSize(D);
//    	instruction_pane.setMinimumSize(D);
//    	instruction_pane.setMaximumSize(D);
    	instruction_pane.setBounds(D.width/8,D.height/8,D.width*3/4,D.height*3/4);
    	
    	
//    	Dimension D = this.getSize();
//    	D.width = D.width/2;
//    	D.height = D.height/2;

//    	instruction_container.setSize(D);
//    	instruction_container.setMinimumSize(D);
//    	instruction_container.setMaximumSize(D);
    	main_panel.add(instruction_container,ROOT_KEY);
    }
    
    public Session addSession(Session session)
    {
//    	session.addChangeListener(e->setMenusAndTitle());
//// BUNDLE
//    	session.getModelBrowser().addSelectionChangeListener(e->setMenusAndTitle());
//// BUNDLE
//    	session.getDataBrowser().addSelectionChangeListener(e->setMenusAndTitle());
    	session.addSelectionChangeListener(e->setMenusAndTitle());
    	session.addSelectionChangeListener(tool_bar);
    	
    	
    	
    	session_list.add(session);
        String key =  session.getSessionKey(); //Integer.toString(session_list.size()); // session.getAssociatedFile().getName()
        
        
        //main_panel.invalidate();
        main_panel.add(session, key);
        showSession(session_list.size()-1);
        
        
//        setMenusAndTitle();
//        tool_bar.updateEnabledStates();
        main_panel.validate();
        return session;
    }
    
    
    
    /**
     * Gets the JPEG-based icon from the jarfile, or via http. 
     * 
     * @return null if picture file is not found
     */
    public static ImageIcon getCountIcon()
    {
        java.net.URL icon_url = ClassLoader.getSystemResource("img/count-icon.jpg");
        try 
        {
            if (icon_url == null)
            {
            	System.out.println("#**AF.gCI Getting icon image by http");
                icon_url = new java.net.URL("http://www.iro.umontreal.ca/~csuros/gene_content/images/count-icon.jpg");
            }
        
            return new ImageIcon(icon_url);
        } catch (Exception E)
        {
            return null;
        }
    }
    
    
    private JComponent getAboutComponent()
    {
    	String ABOUT_TEXT 
        = "<h1>The "+APP_TITLE+" (v"+APP_VERSION+")</h1>" +
        "<p>The Count is a software package " +
        "for the evolutionary analysis of phylogenetic profiles, and numerical characters in general, " +
        "written entirely in Java. The main page for the software is <a href=\"https://github.com/csurosm/count\">https://github.com/csurosm/count</a> " +
        "</p>" +
        "<p>Author: Mikl&oacute;s Cs&#369;r&ouml;s <a href=\"http://www.iro.umontreal.ca/~csuros/\">http://www.iro.umontreal.ca/~csuros/</a></p> " +
        "<h2>Copyright and licenses</h2>"+
        "<ul><li>Licensed under the Apache License, Version 2.0 <a href=\"http://www.apache.org/licenses/LICENSE-2.0\">http://www.apache.org/licenses/LICENSE-2.0</a>."+
        "</li>"+
        "<li> \u261c The background image (Moon and Mars) of the Count logo (batty batty bat batty bat batty bat batty bat...) " +
        "is used with permission from the copyright owner, John Harms." +
        "</li>"+
        "<li>Some numerical optimization routines were adapted from " +
        "<em>Numerical Recipes in C: The Art of Scientific Computing</em> " +
        "[W. H. Press, S. A. Teukolsky, W. V. Vetterling and B. P. Flannery; " +
        "Second edition, Cambridge University Press, 1997].</li>" +
        "<li>Some icons are from the Java Look and Feel Graphics Repository, licensed under the Oracle Binary Code License Agreement for Java SE.</li></ul>"+
        "<h2>References</h2>"+
        "<p>Please <b>cite</b> Count as "+UsersGuide.COUNT_REFERENCE+"</p>" +
        "<p>Algorithmic ideas uderlying the Count package were described in the following publications.</p>" +
        "<ul>" +
        UsersGuide.METHOD_REFERENCES+
        "</ul>" +
//        "<p>Montr&eacute;al/Budapest/Ann Arbor/Hong Kong/Nha Trang/Bengaluru/Amsterdam/Szentendre/Szeged, 2010-2023</p>";
        "<p>Montr&eacute;al, 2010–2023</p>";

    	JEditorPane EP = new JEditorPane("text/html",ABOUT_TEXT);
        EP.setEditable(false);
        EP.setBackground(this.getBackground());
        EP.setPreferredSize(new Dimension(800,600));
        JScrollPane ep_scroll = new JScrollPane(EP);
        EP.setCaretPosition(0);
//        ep_scroll.getVerticalScrollBar().setValue(0);
        return ep_scroll;
    }
    
    private void setMenusAndTitle()
    {
    	menu_bar.setMenus();
    	Session sesh = getActiveSession();
    	if (sesh == null)
    	{
            setTitle(app.APP_TITLE); // +" v"+VERSION_NUMBER);
    	} else
    	{
    		setTitle(sesh.getSessionTitle());    	
    	}
    }
    
    private CountTools createToolBar()
    {
    	final Dimension gap = new Dimension(12,0);
    	CountTools createToolBar = new CountTools(this, "");
    	JButton bouton ; // reused for tool tip texts 
    	bouton = createToolBar.addLoadSession("Load session", e->doLoadSessions());
    	bouton.setToolTipText("Load previously saved sessions (limited support for pre-2023 versions)");
    	bouton = createToolBar.addNewSession("New session", e->doStartSession());
    	bouton.setToolTipText("Load a Newick-format phylogeny for a new session");
    	bouton = createToolBar.addInitSession("No-tree session", e->doInitSession());
    	bouton.setToolTipText("Load a copy-number table for a new session, with a simple initial phylogeny (star, random or UPGMA tree)");
    	bouton = createToolBar.addSaveSession("Save all", e->doSaveSessions());
    	bouton.setToolTipText("Save all current sessions");
    	bouton = createToolBar.addCloseSession("Close session", e->removeSession(getActiveSession()));
    	bouton.setToolTipText("Remove the session without saving");
    	createToolBar.addSeparator();
    	bouton = createToolBar.addLoadTree("Load tree", e->doLoadTree());
    	bouton.setToolTipText("Load an alternative phylogeny into the current session");
    	bouton = createToolBar.addEditTree("Edit tree", e->doEditTree());
    	bouton.setToolTipText("Build a tree by hand: prune-and-regraft, contract edges, or reroot");
    	bouton = createToolBar.addBuildTree("Build tree", e->doBuildTree());
    	bouton.setToolTipText("Infer a phylogeny from the copy numbers or randomly");
    	createToolBar.addSeparator(gap);
    	
    	bouton = createToolBar.addLoadAnnotations("Load table", e->doLoadTable(true));
    	bouton.setToolTipText("Load a copy-number table, keep columns with family annotations");
    	bouton = createToolBar.addImportTable("Import table", e->doImportTableData());
    	bouton.setToolTipText("Load copy-number data from membership (COG) or clustering (MCL) data file.");
    	
    	bouton = createToolBar.addSimulation("Simulate table", e->getActiveSession().doSimulation());
    	bouton.setToolTipText("Create a copy-number table by simulating evolution with the selected rate model");
    	createToolBar.addSeparator();
    	bouton = createToolBar.addFilterRows("Filter rows", e->getActiveSession().showFilteredFamilies());
    	bouton.setToolTipText("Extract the selected families into a new table");
    	bouton = createToolBar.addBinaryTable("Binary", e->getActiveSession().showBinaryProfiles());
    	bouton.setToolTipText("Convert to presence-absence table");
    	createToolBar.addSeparator(gap);
    	
    	bouton=createToolBar.addLoadRates("Load model", e->doLoadRates());
    	bouton.setToolTipText("Load a rate model");
    	bouton = createToolBar.addOptimizeRates("Optimize model", e->getActiveSession().doOptimize());
    	bouton.setToolTipText("Set model parameters by numerical optimization");

    	createToolBar.addSeparator();
    	bouton = createToolBar.addDollo("Dollo", e->getActiveSession().doDollo());
    	bouton.setToolTipText("Ancestral reconstruction by Dollo parsimony (from binary profiles)");
    	bouton = createToolBar.addParsimony("Parsimony", e->getActiveSession().doParsimony());
    	bouton.setToolTipText("Ancestral reconstruction by numerical parsimony (generalized Wagner parsimony)");
    	bouton = createToolBar.addPosteriors("Posteriors", e->getActiveSession().doPosteriors());
    	bouton.setToolTipText("Ancestral reconstruction by posteriors (with the selected table and rate model)");
    	createToolBar.addSeparator();
    	bouton = createToolBar.addSave("Export", e->getActiveSession().doSaveData());
    	bouton.setToolTipText("Save the currently displayed item");
    	bouton = createToolBar.addRemove("Remove", e->getActiveSession().doRemoveItem());
    	bouton.setToolTipText("Delete the currently displayed item without saving");
    	return createToolBar; 
    }
    
    private class MenuBar extends JMenuBar
    {
    	final JMenu session_menu;
    	final JMenu data_menu;
    	final JMenu rate_menu;
    	final JMenu analysis_menu;      
    	final JMenu help_menu;
    	
    	JMenuItem load_session;
    	JMenuItem start_session;
    	JMenuItem init_session;
    	
    	MenuBar()
    	{
            session_menu = new JMenu("Session");
            data_menu = new JMenu("Data");
            rate_menu = new JMenu("Model");
            analysis_menu = new JMenu("Analysis");
            help_menu = new JMenu("Help");
            
            
            add(session_menu);
            add(data_menu);
            add(rate_menu);
            add(analysis_menu);
            add(help_menu);

            setMenus();
    	}
    	
    	void setMenus()
    	{
    		setSessionMenu();
    		setDataMenu();
    		setRateMenu();
    		setAnalysisMenu();
    	}
    	
    	void setSessionMenu()
    	{
    		load_session = new JMenuItem(CountActions.createLoadSession("Open previously saved session(s) ...", e->doLoadSessions()));
//            load_session = new JMenuItem("Open previously saved session(s) ...");
//            load_session.addActionListener(e->doLoadSession());
            load_session.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, 
                Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx())); 
            load_session.setToolTipText("This tooltip has no information"); // no tooltips in menus
            
            JMenuItem save_session = new JMenuItem(CountActions.createSaveAll("Save everything...", e->doSaveSessions()));
//            JMenuItem save_session = new JMenuItem("Save everything ...");
//            save_session.addActionListener(e->doSaveSession());
            save_session.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, 
                Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx())); 
//            JMenuItem load_session_url = new JMenuItem("Load session by path name or URL...");    		
            
            JMenuItem close_session = new JMenuItem(CountActions.createDelete("Close session (without saving)", e->removeSession(getActiveSession())));
            
            start_session = new JMenuItem(CountActions.createNewSession("Start new session with a known phylogeny...", e->doStartSession()));
//            start_session = new JMenuItem("Start new session ...");
//            start_session.addActionListener(e->doStartSession());
            start_session.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N, 
                Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx())); 

            init_session = new JMenuItem(CountActions.createInitSession("Start new session with a table and no phylogeny...", e->doInitSession()));
            
            
            JMenuItem nog_session = new JMenuItem("Start an EggNOG session ...");
            nog_session.addActionListener(e->doNOGSession());

            JMenuItem add_tree = new JMenuItem(CountActions.createLoadTree("Load phylogeny into active session ...", e->doLoadTree()));
//            JMenuItem add_tree = new JMenuItem("Load phylogeny into active session ...");
//            add_tree.addActionListener(e->doLoadTree());
            
            JMenuItem edit_tree = new JMenuItem(CountActions.createEditTree("Edit tree", e->doEditTree()));
//            JMenuItem edit_tree = new JMenuItem("Edit tree");
//            edit_tree.addActionListener(e->doEditTree());
            
            JMenuItem build_tree = new JMenuItem(CountActions.createBuildTree("Build tree", e->doBuildTree())); 
            
            
            session_menu.removeAll();
            
            session_menu.add(load_session);
//            session_menu.add(load_session_url);
            session_menu.add(start_session);
            session_menu.add(init_session);
//            session_menu.add(nog_session);
            
            class SessionActivator implements ActionListener
            {
            	private final int session_idx;
            	SessionActivator(int session_idx)
            	{ this.session_idx = session_idx;}
            	
                @Override
                public void actionPerformed(ActionEvent e)
                {
                	showSession(session_idx);
                }
            }
            
            if (!session_list.isEmpty())
            {
                session_menu.addSeparator();
                for (int j=0; j<session_list.size(); j++)
                {
                    Session sesh = session_list.get(j);
                    JMenuItem item = null;
                    if (j==active_session) 
                        item = new JMenuItem(sesh.getSessionKey(), StringIcon.createRightPointingFinger());
                    else 
                        item = new JMenuItem(sesh.getSessionKey());
                    item.addActionListener(new SessionActivator(j));
                    session_menu.add(item);
                }
            }
            session_menu.addSeparator();
            session_menu.add(add_tree);
            session_menu.add(edit_tree);
            session_menu.add(build_tree);
            
            session_menu.addSeparator();
            session_menu.add(close_session);
            session_menu.add(save_session);
            
            Session active = getActiveSession();
            load_session.setEnabled(session_list.isEmpty());
//            load_session_url.setEnabled(load_session.isEnabled() && false); // not implemented
            start_session.setEnabled(true);
            init_session.setEnabled(start_session.isEnabled());
            nog_session.setEnabled(start_session.isEnabled());
            add_tree.setEnabled(active != null);
            edit_tree.setEnabled(active != null);
            build_tree.setEnabled(active != null);
            save_session.setEnabled(active != null);
            close_session.setEnabled(active != null);
    	}
    	
    	void setDataMenu()
    	{
    		data_menu.removeAll();
    		
    		JMenuItem data_load_table = new JMenuItem(CountActions.createLoadTable("Open table... (keep only the columns with matching leaf names)",e->doLoadTable(false) ));
//            JMenuItem data_load_table = new JMenuItem("Open table... (keep only the columns with matching leaf names)");
//            data_load_table.addActionListener(e->doLoadTable(false));
            data_load_table.setAccelerator(
            		KeyStroke.getKeyStroke(KeyEvent.VK_T,
            				Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx()));
            data_menu.add(data_load_table);       

            JMenuItem data_load_annotated = new JMenuItem(CountActions.createLoadAnnotations("Open annotated table ... (keep nomatching columns as annotations)", e->doLoadTable(true)));
//            JMenuItem data_load_annotated = new JMenuItem("Open annotated table... (keep other columns as annotations)");
//            data_load_annotated.addActionListener(e->doLoadTable(true));
            data_menu.add(data_load_annotated);

//            JMenuItem data_load_url = new JMenuItem("Load table by path name or URL...");
//            data_load_url.addActionListener(e->doLoadTablePath());
//            data_menu.add(data_load_url);
            
            JMenuItem data_import = new JMenuItem(CountActions.createImportTable("Import family profile data from clustering or membership data (MCL and COG formats)", e->doImportTableData()));
            data_menu.add(data_import);
            
            JMenuItem data_simulation = new JMenuItem(CountActions.createSimulation("Simulate a table by rate model...", e->getActiveSession().doSimulation()));
//            JMenuItem data_simulation = new JMenuItem("Simulate a table by rate model...");
//            data_simulation.addActionListener(e->getActiveSession().doSimulation());
            data_menu.add(data_simulation);
            
            data_menu.addSeparator();
         
            JMenuItem data_filter_lines = new JMenuItem(CountActions.createFilterRows("Extract selected families into a new table", e->getActiveSession().showFilteredFamilies()));
//            JMenuItem data_filter_lines = new JMenuItem("Extract selected families into a new table");
//            data_filter_lines.addActionListener(e->getActiveSession().showFilteredFamilies());

            data_menu.add(data_filter_lines);

//            JMenuItem data_split_table = new JMenuItem("Split table by entries in a column...");
//            data_split_table.addActionListener(e->getActiveSession().showSplitTable());
//            data_menu.add(data_split_table);


            JMenuItem data_binary = new JMenuItem(CountActions.createBinaryTable("Convert to binary profiles", e->getActiveSession().showBinaryProfiles()));
//            JMenuItem data_binary = new JMenuItem("Convert to binary profiles");
//            data_binary.addActionListener(e->getActiveSession().showBinaryProfiles());
            data_menu.add(data_binary);   
            
            
            // enable/disable 
            Session sesh = getActiveSession();
 // BUNDLE
            ModelBundle.Entry table_entry=null;
 // BUNDLE
            if (sesh!=null)
            {
// BUNDLE
//            	table_panel = sesh.getDataBrowser().getSelectedPrimaryItem();
            	table_entry = sesh.getDataBrowser().getSelectedTableEntry();
// BUNDLE
            }
            data_load_table.setEnabled(sesh != null);
            data_load_annotated.setEnabled(data_load_table.isEnabled());
//            data_load_url.setEnabled(data_load_table.isEnabled());
            data_filter_lines.setEnabled(table_entry!=null && 
            		sesh.getDataBrowser().getSelectedItem() instanceof Session.FamilySelection); 
//            data_split_table.setEnabled(table_panel!=null);
            data_binary.setEnabled(table_entry!=null
            		&& !table_entry.getTableData().getContent().isBinaryTable());
            
            ModelBundle.Entry rates_entry
								= (sesh==null)
								? null
								: sesh.getModelBrowser().getSelectedRatesEntry();
            data_simulation.setEnabled(rates_entry != null);
            
            // not yet implemented
//            data_load_url.setEnabled(false);
//            data_filter_lines.setEnabled(false);
//            data_binary.setEnabled(false);
    	}
    	
    	void setRateMenu()
    	{
    		rate_menu.removeAll();
    		
    		JMenuItem rates_load = new JMenuItem(CountActions.createLoadRates("Load rates...", e->doLoadRates()));
//            JMenuItem rates_load = new JMenuItem("Load rates...");
//            rates_load.addActionListener(e->doLoadRates(null)); // null argument triggers FileDialog
            rates_load.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_R,
            		Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx()));
            rate_menu.add(rates_load);
        
            rate_menu.addSeparator();
            JMenuItem rates_optimize = new JMenuItem(CountActions.createOptimizeRates("Optimize rates...", e->getActiveSession().doOptimize()));
//            JMenuItem rates_optimize = new JMenuItem("Optimize rates...");
//            rates_optimize.addActionListener(e->getActiveSession().doOptimize());
            rate_menu.add(rates_optimize);

//            JMenuItem rates_score = new JMenuItem("Assess model fit...");
//            rate_menu.add(rates_score); 
            
            // enable/disable 
            Session sesh = getActiveSession();
// BUNDLE           AnnotatedTablePanel table_panel
            ModelBundle.Entry table_entry 
            					= (sesh==null)
            						?null
            						:sesh.getDataBrowser().getSelectedTableEntry();
            rates_load.setEnabled(sesh!=null);
            rates_optimize.setEnabled(table_entry!=null);
//            rates_score.setEnabled(false); // not implemented yet
            
    	}
    	
    	void setAnalysisMenu()
    	{
    		analysis_menu.removeAll();
    		
    		JMenuItem analysis_dollo = new JMenuItem(CountActions.createDollo("Family history by Dollo parsimony", e->getActiveSession().doDollo()));
//            JMenuItem analysis_dollo = new JMenuItem("Family history by Dollo parsimony");
//            analysis_dollo.addActionListener(e->getActiveSession().doDollo());
            analysis_menu.add(analysis_dollo);
            
            JMenuItem analysis_wagner = new JMenuItem(CountActions.createParsimony("Family history by numerical parsimony", e->getActiveSession().doParsimony()));
//            JMenuItem analysis_wagner = new JMenuItem("Family history by numerical parsimony");
//            analysis_wagner.addActionListener(e->getActiveSession().doParsimony());
            analysis_menu.add(analysis_wagner);

            JMenuItem analysis_posterior = new JMenuItem(CountActions.createPosteriors("Family history by posterior probabilities", e->getActiveSession().doPosteriors()));
//            JMenuItem analysis_posterior = new JMenuItem("Family history by posterior probabilities");
//            analysis_posterior.addActionListener(e->getActiveSession().doPosteriors());
            analysis_menu.add(analysis_posterior);

//            analysis_menu.addSeparator();
//            JMenuItem analysis_PGL = new JMenuItem("PGL: propensity for gene loss (Krylov-Wolf-Rogozin-Koonin)");
//            analysis_menu.add(analysis_PGL);
//            analysis_PGL.addActionListener(e->getActiveSession().showPGL());
            
            // enable/disable 
            Session sesh = getActiveSession();
// BUNDLE            AnnotatedTablePanel table_panel
            ModelBundle.Entry table_entry 
            					= (sesh==null)
            						?null
            						:sesh.getDataBrowser().getSelectedTableEntry();
//  BUNDLE           RateVariationPanel rates_panel
            ModelBundle.Entry rates_entry 
            					= (sesh==null)
            					? null
            					: sesh.getModelBrowser().getSelectedRatesEntry();
            analysis_dollo.setEnabled(table_entry!=null && table_entry.getTableData().getContent().isBinaryTable());
            analysis_wagner.setEnabled(table_entry != null);
            analysis_posterior.setEnabled(table_entry != null && rates_entry != null);
            
    	}
    }
    
    /**
     * Opens up a file dialog to read a phylogeny from a Newick-format file.
     * 
     * @return null if failed 
     */
    private DataFile<Phylogeny> loadPhylo()
    {
        FancyFileDialog dialog = null;
        dialog = new FancyFileDialog(this,"Open Newick-format phylogeny file",FileDialog.LOAD);
        dialog.setVisible(true);
       
        String file_name = dialog.getFile();
        String directory = dialog.getDirectory();

        Phylogeny main_phylo=null;
        
        if (file_name != null)
        {
        	String new_session_key = Session.getSessionKey(file_name);
        	for (Session session: session_list)
        		if (new_session_key.equals(session.getSessionKey()))
                {
                    exception_handler.handle(new RuntimeException("Duplicate session names"), "Duplicate names", "Cannot have two open sessions with the same phylogeny file names");
                    file_name = null;
                    break;
                } 
        }
        
        DataFile<Phylogeny> data_read = null;
        
        if (file_name != null)
        {
            try 
            {
                java.io.Reader R = dialog.getReader();
//                = new java.io.InputStreamReader(new java.io.FileInputStream(directory+file_name)); // we shun FileReader for no particular reason
                main_phylo = NewickParser.readTree(R);
            } catch (java.io.InterruptedIOException E)
            {
                // load canceled
                main_phylo = null;
            } catch (NewickParser.ParseException E)
            {
            	exception_handler.handle(E, "Newick file parsing error", "Parsing error while reading Newick file "+directory+file_name+".");
            } catch (java.io.IOException E)
            {
            	exception_handler.handle(E, "I/O error", "File error while attempting ro read Newick file "+directory+file_name+".");  
            }
            
            if (main_phylo != null)
            {
                boolean all_zero = true;
                int num_nodes = main_phylo.getNumNodes();
                for (int j=0; j<num_nodes && all_zero; j++) 
                	if (!main_phylo.isRoot(j)) // not for root
                		all_zero = (main_phylo.getLength(j)==0.);
                
                if (all_zero)
                {
                    for (int j=0; j<num_nodes; j++) 
                    	if (!main_phylo.isRoot(j))// don't change root
                    		main_phylo.getNode(j).setLength(1.0);
                }
                
                int num_fixed_edges = main_phylo.fixZeroEdges();
//            	if (num_fixed_edges>0)
//            	{
//            		System.out.println("#**AF.lP fixed "
//            			+num_fixed_edges
//            			+"\tn/l "+num_nodes+"->"+main_phylo.getNumNodes()+"/"+main_phylo.getNumLeaves());
//            	}
                
                data_read = new DataFile<>(main_phylo, new File(directory,file_name));
            }        
        }
        
        return data_read;
    }    
    
    
    /**
     * Shortcut to {@link count.Count#doQuit()}
     */
    public void doQuit() { app.doQuit();}
    
    
    /**
     * A modal dialog about Count. 
     * 
     * @param splash whether this is the starting splash screen 
     */
    public void doAbout(boolean splash)
    {
        String title = (splash?"Welcome to ":"About ")+APP_TITLE;
        if (splash) //splash)
        {
            String[] first_step = {"Just close the splash screen", 
            		menu_bar.load_session.getText(),  menu_bar.start_session.getText()};
            int selected = JOptionPane.showOptionDialog(this, getAboutComponent(), 
            		title, JOptionPane.DEFAULT_OPTION, 
            		JOptionPane.INFORMATION_MESSAGE, 
            		getCountIcon(), 
            		first_step, first_step[0]);
            if (selected == 1)
            {
                doLoadSessions();
            } else if (selected==2)
            {
                doStartSession();
            }
        } else
        {
            JOptionPane.showMessageDialog(this, getAboutComponent(), 
            		title, JOptionPane.INFORMATION_MESSAGE,
                    getCountIcon());

        }
    }
    
    
    private void doLoadSessions()
    {
        FancyFileDialog dialog = null;
        dialog = new FancyFileDialog(this,"Load sessions",FileDialog.LOAD);
        dialog.setVisible(true);
       
        final String file_name = dialog.getFile();
        final String directory = dialog.getDirectory();
        
        if (file_name != null)
        {
            final JDialog wait_until_loaded = new JDialog(this,"Loading sessions",true);
            BoxLayout load_layout = new BoxLayout(wait_until_loaded.getContentPane(),BoxLayout.PAGE_AXIS);
            wait_until_loaded.setLayout(load_layout);
            wait_until_loaded.add(Box.createVerticalGlue());
            wait_until_loaded.add(new JLabel("Reading and parsing file content: please be patient. "));
            JProgressBar progress = new JProgressBar();
            progress.setIndeterminate(true);
            Box progress_box = new Box(BoxLayout.LINE_AXIS);
            progress_box.add(Box.createHorizontalStrut(12));
            progress_box.add(progress);
            progress_box.add(Box.createHorizontalStrut(12));
            wait_until_loaded.add(progress_box);
            wait_until_loaded.add(Box.createVerticalGlue());

            
            SwingWorker<Void,Void> load_worker = new SwingWorker<Void,Void>()
            {
                @Override
                public Void doInBackground()
                {
                    try 
                    {
                        File selected_file = new File(directory,file_name);
                        BufferedReader BR = GeneralizedFileReader.guessBufferedReaderForInput(selected_file.getPath());
                        List<ModelBundle> bundle_list = ModelBundle.readBundle(BR);
                        BR.close();
                        for (ModelBundle bundle: bundle_list)
                        {
                        	Session sesh = new Session(AppFrame.this, bundle);
                        	addSession(sesh);
                        	main_panel.repaint();
                        }
                    } catch (java.io.IOException E)
                    {
                        exception_handler.handle(E, "I/O error", "File error while loading sessions.");
                    } catch(org.xml.sax.SAXException E)
                    {
                    	exception_handler.handle(E, "File format problem", "The file is not formatted correctly");
                    } catch (Exception E)
                    {
                    	exception_handler.handle(E, "A bug maybe?", "Error while loading sessions.");
                    }

                    return null;
                }
                
                @Override
                public void done()
                {
                	
                    wait_until_loaded.setVisible(false);
                }
            };
            load_worker.execute();
            wait_until_loaded.pack();
            java.awt.Dimension frameD = this.getSize();
            wait_until_loaded.setBounds(frameD.width/3, frameD.height/3, frameD.width/3, frameD.height/6);
            wait_until_loaded.setVisible(true);
        }
    	
//        throw new UnsupportedOperationException();
    }
    
    private void doSaveSessions()
    {
        FileDialog dialog = null;
        dialog = new FileDialog(this,"Save all work",FileDialog.SAVE);
        
        dialog.setVisible(true);
        
        final String file_name = dialog.getFile();
        final String directory = dialog.getDirectory();
        
        if (file_name != null)
        {
            final JDialog wait_until_saved = new JDialog(this,"Saving sessions",true);
            BoxLayout saving_layout = new BoxLayout(wait_until_saved.getContentPane(),BoxLayout.PAGE_AXIS);
            wait_until_saved.getContentPane().setLayout(saving_layout);
            wait_until_saved.add(Box.createVerticalGlue());
            wait_until_saved.add(new JLabel("Saving session data to disk: please wait ..."));
            JProgressBar progress = new JProgressBar();
            progress.setIndeterminate(true);
            Box progress_box = new Box(BoxLayout.LINE_AXIS);
            progress_box.add(Box.createHorizontalStrut(12));
            progress_box.add(progress);
            progress_box.add(Box.createHorizontalStrut(12));
            wait_until_saved.add(progress_box);
            wait_until_saved.add(Box.createVerticalGlue());

            SwingWorker<Void,Void> save_worker = new SwingWorker<Void,Void>()
            {
                @Override
                public Void doInBackground()
                {
                    try 
                    {
                        File selected_file = new File(directory,file_name);
                        PrintStream out = new PrintStream(selected_file);
                        List<ModelBundle> bundle_list = new ArrayList<>();
                        for (Session sesh: session_list)
                        	bundle_list.add(sesh.getModelBundle());
                        ModelBundle.printBundles(out, bundle_list);
                        out.close();
                    }
                    catch (java.io.IOException E)
                    {
                        exception_handler.handle(E, "I/O error", "File error while saving sessions.");
                    } catch (Exception E)
                    {
                        exception_handler.handle(E, "A bug maybe?", "Error while saving the sessions.");
                    }

                    return null;
                }

                @Override
                public void done()
                {
                    wait_until_saved.setVisible(false);
                }    
            }; 
            save_worker.execute();
            wait_until_saved.pack();
            java.awt.Dimension frameD = getSize();
            wait_until_saved.setBounds(frameD.width/3, frameD.height/3, frameD.width/3, frameD.height/6);
            wait_until_saved.setVisible(true);            
            
        }        
    }
    

    private void doStartSession()
    {
        DataFile<Phylogeny> did_load = loadPhylo();
        if (did_load != null)
            addSession(new Session(this, did_load));
    }
    
    private void doInitSession()
    {
        String dialog_title = ("Open family size table");
        FileDialog dialog  = new FancyFileDialog(this,dialog_title,FileDialog.LOAD);
        dialog.setVisible(true);
        
        String file_name = dialog.getFile();
        if (file_name != null)
        {
            DataFile<AnnotatedTable> table_data = null;
            String directory = dialog.getDirectory();
            try
            {
	        	AnnotatedTable table 
	        	= TableParser
	        	.readTable(null, 
	        			GeneralizedFileReader.guessReaderForInput(directory+file_name), 
	        			false);
	        		        	
	            checkTableEmpty(table);
	            table_data = new DataFile<>(table, new File(directory, file_name));
	        } catch (java.io.InterruptedIOException E)
	        {
	            // Canceled: nothing to do
	        } catch (java.io.IOException E)
	        {
	        	exception_handler.handle(E, "I/O error", "File error while reading occurrence table from a file.");
	        } catch (IllegalArgumentException E)
	        {
	        	exception_handler.handle(E, "Parsing error", "Cannot parse occurrence table in file.");
	        } catch (Exception E)
	        {
	        	exception_handler.handle(E, "A bug!", "Error while reading occurrence table from a file.");
	        }   
            if (table_data != null)
            	initSession(table_data);
        }
    }
    
    private void doImportTableData()
    {
		Session sesh = getActiveSession();
    	ImportTable importer = new ImportTable(this);
    	
    	int wantsto = JOptionPane.showConfirmDialog(this, importer, "Import family profile data", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, CountActions.createLoadTableIcon(CountActions.SIZE_L));
    	if (wantsto == JOptionPane.OK_OPTION)
    	{
    		String[] taxon_names;
    		if (sesh==null)
    		{
    			taxon_names =null;
    		} else
    		{
    			taxon_names = sesh.getModelBrowser().getMainPhylogeny().getLeafNames();
    		}
    		final SwingWorker<DataFile<AnnotatedTable>,Void> read_task = importer.readTableTask(taxon_names);
    		read_task.addPropertyChangeListener(event->
    		{
    			if (! read_task.isCancelled()
    					&&
    					"state".equals(event.getPropertyName())
    	                 && SwingWorker.StateValue.DONE == event.getNewValue())
    			{
    				DataFile<AnnotatedTable> table_data=null;
//    				System.out.println("#**AF.bIT task done");
    				try {table_data = read_task.get();} 
    				catch (InterruptedException ignored) {}
					catch (ExecutionException ignored) {}
    		        if (table_data != null)
    		        {
    		            checkTableEmpty(table_data.getContent());
    		        	if (sesh==null)
    		        		initSession(table_data);
    		        	else
    		        	{
    		        		// no need for mapping bc readTable used the taxon name order from the main phylo
    		        		sesh.addDataSet(table_data);
    		        	}
    		        }
    			}
    		});
    		read_task.execute();
    	}
    	
    }
    
//    private void foregroundImportTableData()
//    {
//		Session sesh = getActiveSession();
//    	ImportTable importer = new ImportTable(this);
//    	DataFile<AnnotatedTable> table_data = null;
//        try
//        { 
//        	int wantsto = JOptionPane.showConfirmDialog(this, importer, "Import family profile data", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, CountActions.createLoadTableIcon(CountActions.SIZE_L));
//        	if (wantsto == JOptionPane.OK_OPTION)
//        	{
//        		String[] taxon_names;
//        		if (sesh==null)
//        		{
//        			taxon_names =null;
//        		} else
//        		{
//        			taxon_names = sesh.getModelBrowser().getMainPhylogeny().getLeafNames();
//        		}
//        		table_data = importer.readTable(taxon_names);
//        	}
//        } catch (IOException E)
//        {
//            exception_handler.handle(E, "I/O error", "File error while importing table data.");
//        } catch (Exception E)
//        {
//            exception_handler.handle(E, "A bug maybe?", "Error while importing table data.");
//        }
//        if (table_data != null)
//        {
//            checkTableEmpty(table_data.getContent());
//        	if (sesh==null)
//        		initSession(table_data);
//        	else
//        	{
//        		// no need for mapping bc readTable used the taxon name order from the main phylo
//        		sesh.addDataSet(table_data);
//        	}
//        }
//    }
    
    private void initSession(DataFile<AnnotatedTable> table_data)
    {
        InitialTreeDialog init_dialog = new InitialTreeDialog(this, "Build initial phylogeny for "+table_data.getFile(), table_data);
        init_dialog.setVisible(true);
        // wait for answer
        DataFile<Phylogeny> tree_data = init_dialog.getTreeData();
    	if (tree_data != null)
    	{
    		Session sesh = new Session(this, tree_data);
            addSession(sesh);
    		AnnotatedTable mapped_table = table_data.getContent().mappedToTree(tree_data.getContent());
    		DataFile<AnnotatedTable> mapped_table_data = new DataFile(mapped_table, table_data.getFile());
    		sesh.addDataSet(mapped_table_data);
    	}
    }
    
    private void doBuildTree()
    {
    	Session sesh = getActiveSession();
    	DataFile<AnnotatedTable> table_data = sesh.getSelectedData(); // this is mapped to the currently selected tree 
    	
        InitialTreeDialog init_dialog = new InitialTreeDialog(this, "Build initial phylogeny", table_data);
        init_dialog.setVisible(true);
        // wait for it
        DataFile<Phylogeny> tree_data = init_dialog.getTreeData();
        if (tree_data != null)
        {
        	sesh.addTopTree(renamedTreeData(tree_data));
        }
    }
    
    private void doNOGSession()
    {
    	String taxon_id = JOptionPane.showInputDialog("Enter a NCBI Taxonomy ID");
    	if (taxon_id != null)
    	{
    		String server = "eggnog5.embl.de";
    		String URL_root = "http://"+server+"/download/eggnog_5.0/per_tax_level";
    		String URL_dir = URL_root+"/"+taxon_id;
    		String URLprefix = URL_dir+"/"+taxon_id+"_";
    		String URL_members = URLprefix+"members.tsv.gz";
    		String URL_annotations = URLprefix+"annotations.tsv.gz";
    		String URL_trees = URLprefix+"trees.tsv.gz";
    		
    		System.out.println("#**AF.dNS downloading "+URL_members);
    		
    		try
    		{
    			InputStream members_stream = GeneralizedFileReader.guessInputStreamForInput(URL_members);
	    		Reader members_reader = 
	    				new InputStreamReader(
	    						new ProgressMonitorInputStream(AppFrame.this, 
	    								"Downloading "+URL_members, 
	    								members_stream));
	    		InputStream annotations_stream = GeneralizedFileReader.guessInputStreamForInput(URL_annotations);
	    		Reader annotations_reader = new InputStreamReader(
						new ProgressMonitorInputStream(AppFrame.this, 
								"Downloading "+URL_annotations, 
								annotations_stream));
	    		InputStream trees_stream = GeneralizedFileReader.guessInputStreamForInput(URL_trees);
	    		Reader trees_reader = new InputStreamReader(
						new ProgressMonitorInputStream(AppFrame.this, 
								"Downloading "+URL_trees, 
								trees_stream));
	    		
	    		SwingWorker<AnnotatedTable, Void> download_task
	    		= new SwingWorker<>()
	    		{
	    			private List<DataFile<Phylogeny>> gene_trees = new ArrayList<>();
	    			private AnnotatedTable nog_table;

	    			@Override
	    			protected AnnotatedTable doInBackground() throws IOException
	    			{
	    				try {
		    	    		System.out.println("#**AF.dNS/SW.dIB start "+nog_table);
		    	    		nog_table = TableParser.eggNOG(members_reader, annotations_reader, trees_reader, gene_trees);
		    	    		System.out.println("#**AF.dNS/SW.dIB done "+nog_table+" trees "+gene_trees.size());
		    	    		System.out.flush();
		    	    		members_stream.close();
		    	    		annotations_stream.close();
		    	    		trees_stream.close();
	    				} catch (Throwable t)
	    				{
	    					t.printStackTrace();
	    					System.exit(999);
	    				}
	    	    		return nog_table;
	    			}
	    			
	    			@Override
	    			protected void done()
	    			{
	    	    		System.out.println("#**AF.dNS/SW.done "+nog_table+" trees "+gene_trees.size());
	    				if (!isCancelled())
	    				{
	    					if (gene_trees.isEmpty())
	    					{
	    						JOptionPane.showMessageDialog(AppFrame.this, 
	    								"No available gene trees (need exactly 1 copy at each leaf)"	
	    								, "Cannot open session", JOptionPane.ERROR_MESSAGE);;
	    					} else
	    					{
			    	    		File nog_file = new File(URL_root, taxon_id);
			    	    		DataFile<AnnotatedTable> nog_table_data = new DataFile<>(nog_table, nog_file);
			    	    		nog_table_data.setDirty(true);
			    	    		
			    	    		
			    	    		DataFile<Phylogeny> main_data;
//			    	    		if (gene_trees.isEmpty())
//			    	    		{
//			    	    			System.out.println("#**AF.dNS no trees for "+nog_table);
//			    	    			
//			    	    			Phylogeny star = Phylogeny.starTree(nog_table.getTaxonNames());
//			    	    			main_data = new DataFile<>(star, new File((File)null, taxon_id));
//			    	    		} else
		    	    			main_data = gene_trees.get(0);
		
			    	    		Session sesh = addSession(new Session(AppFrame.this, main_data));
			    	    		sesh.addDataSet(nog_table_data);
			    	    		
			    	    		for (int t=1; t<gene_trees.size(); t++)
			    	    			sesh.addTopTree(gene_trees.get(t));
			    	    		
			    	    		StringBuilder sb = new StringBuilder();
			    	    		sb.append("<h1>").append(taxon_id).append("</h1>");
			    	    		sb.append("<p><b>EggNOG dataset for ")
			    	    			.append(taxon_id)
			    	    			.append(" loaded successfully!</b></p>");
			    	    		sb.append("<h2>Copy number data</h2>");
			    	    		sb.append("<p>A data table named <tt>")
			    	    			.append(nog_table_data.getFile().getName())
			    	    			.append("</tt> was created with ")
			    	    			.append(nog_table.getFamilyCount())
			    	    			.append(" families across ")
			    	    			.append(nog_table.getTaxonCount())
			    	    			.append(" terminal taxa (leaves).</p>");
			    	    		sb.append("<h2>Gene trees</h2>");
			    	    		if (gene_trees.isEmpty())
			    	    		{
			    	    			sb.append("<p>There were no gene trees with exactly 1 gene for "
			    	    					+ "each terminal taxon. A star tree is used.</p>");
			    	    		} else
			    	    		{
				    	    		sb.append("<p>")
				    	    		.append(gene_trees.size())
				    	    		.append(" phylogenies were added, based on "
				    	    				+ "families with exactly one member "
				    	    				+ "at each terminal node.</p>");
				    	    		sb.append("<ul>");
				    	    		for (int t=0; t<gene_trees.size(); t++)
				    	    		{
				    	    			DataFile<Phylogeny> phylo_data = gene_trees.get(t);
				    	    			sb.append("<li>");
				    	    			if (t==0)
				    	    				sb.append("<b>");
				    	    			sb.append(phylo_data.getFile().getName());
				    	    			if (t==0)
				    	    			{
				    	    				sb.append("</b> (main phylogeny)");
				    	    			}
				    	    			sb.append(":");
				    	    			String note = phylo_data.getNote();
				    	    			if (note == null)
				    	    			{
				    	    				sb.append(phylo_data.getContent().getNumNodes())
				    	    					.append(" nodes, added without adjustment");
				    	    			} else if (note.startsWith("star"))
				    	    			{
				    	    				sb.append(note);
				    	    			} else
				    	    			{
				    	    				sb.append(phylo_data.getContent().getNumNodes())
			    	    					.append(" nodes, added by "
			    	    							+ "eliminating 0-length edges ")
			    	    					.append(note);
				    	    			}
				    	    			sb.append("</li>");
				    	    		}
				    	    		sb.append("</ul>");
			    	    		}
			    	    		
			    	            JEditorPane EP = new JEditorPane("text/html", sb.toString());
			    	            EP.setEditable(false);
			    	            EP.setBackground(WARNING_COLOR);
			    	            EP.setBorder(BorderFactory.createRaisedBevelBorder());
			    	    		
			    	            JScrollPane msg_pane = new JScrollPane(EP);
			    	            
			    	            
			    	            // msg_pane.setPreferredSize(new Dimension(900, 600));
			    	            msg_pane.setPreferredSize(new Dimension(900, 600));
			    	            int minw = 360;
			    	            int minh = 120;
			    	            msg_pane.setMinimumSize(new Dimension(minw, minh));
			    	            Dimension frame_size = AppFrame.this.getSize();
			    	            msg_pane.setMaximumSize(new Dimension(Integer.max(frame_size.width-100, minw),600));
			    	            
			    	            
			    	    		JOptionPane.showMessageDialog(AppFrame.this, msg_pane, "EggNOG session",
			    	    				JOptionPane.INFORMATION_MESSAGE,
			    	    				new StringIcon(36,36,"\u270c")); // victory
	    					}
	    				}
	    			}
	    		};
	    		download_task.execute();
    		} catch (Exception E)
    		{
	        	exception_handler.handle(E, "I/O error", "Error while downloading "+taxon_id+" from "+server);
    		}
    	}
    }
    
    private void doLoadTree()
    {
        DataFile<Phylogeny> did_load = loadPhylo();
        if (did_load != null)
            getActiveSession().addTopTree(did_load);
    }
    
    private void doEditTree()
    {
        TreeEditDialog edit_dialog = new TreeEditDialog(this, "Edit tree");
        
        Dimension frameD = this.getSize();

        edit_dialog.pack();
        edit_dialog.setBounds((int)(0.05*frameD.width),(int)(0.05*frameD.height),(int)(0.9*frameD.width),(int)(0.9*frameD.height));

        edit_dialog.setVisible(true);
        DataFile<Phylogeny> did_edit = edit_dialog.getEditedTree();
        if (did_edit != null)
        {
            // rename it
// BUNDLE        	
//            TreeCollection all_trees = getActiveSession().getTreesBrowser();
// BUNDLE            
//            int rename_idx=0;
//            String renamed_ident = DataFile.createIdentifier(rename_idx);
//            while (all_trees.getTree(renamed_ident)!=null) // avoid name clashes
//            {
//                rename_idx++;
//                renamed_ident = DataFile.createIdentifier(rename_idx);
//            }
//
//            File renamed = new File(did_edit.getFile().getParent(), renamed_ident);
//            did_edit = new DataFile<>(did_edit.getContent(), renamed);
            
            getActiveSession().addTree(renamedTreeData(did_edit));
        }    
    }
    
    private DataFile<Phylogeny> renamedTreeData(DataFile<Phylogeny> tree_data)
    {
    	BundleBrowser all_trees = getActiveSession().getModelBrowser();

        int rename_idx=0;
        String renamed_ident = DataFile.createIdentifier(rename_idx);
        while (all_trees.getTree(renamed_ident)!=null) // avoid name clashes
        {
            rename_idx++;
            renamed_ident = DataFile.createIdentifier(rename_idx);
        }

        File renamed = new File(tree_data.getFile().getParent(), renamed_ident);
        DataFile<Phylogeny> renamedTreeData = new DataFile<>(tree_data.getContent(), renamed);
    	return renamedTreeData;
    }
    
     private void doLoadTable(boolean with_annotation_columns)
    {
        String dialog_title = (with_annotation_columns?"Open annotated family size table":"Open family size table");
        FancyFileDialog dialog  = new FancyFileDialog(this,dialog_title,FileDialog.LOAD);
        dialog.setVisible(true);
        String file_name = dialog.getFile();
        
        if (file_name != null)
        {
        	
        	Session sesh = getActiveSession();
        	String[] terminal_names = sesh.getModelBrowser().getMainPhylogeny().getLeafNames();
        	String directory = dialog.getDirectory();
        	File table_file = new File(directory, file_name);
        	SwingWorker<DataFile<AnnotatedTable>, Void> load_task =
        			dialog.createTask(stream->
        			{
        				Reader tableR = new InputStreamReader(stream);
        				AnnotatedTable table 
        	        	= TableParser.readTable(terminal_names,tableR, with_annotation_columns);
        	            checkTableEmpty(table);
        				DataFile<AnnotatedTable> table_data = new DataFile<>(table, table_file);
        				return table_data;
        			});
        	load_task.addPropertyChangeListener(event->
    		{
    			if (! load_task.isCancelled()
    					&&
    					"state".equals(event.getPropertyName())
    	                 && SwingWorker.StateValue.DONE == event.getNewValue())   
    			{
    				DataFile<AnnotatedTable> table_data=null;
//    				System.out.println("#**AF.bIT task done");
    				try {table_data = load_task.get();} 
    				catch (InterruptedException swallowed) {}
    				catch (ExecutionException swallowed) {}
    				if (table_data!=null)
    		            sesh.addDataSet(table_data);
    			}
    		});
			load_task.execute();
        }
    }

    
//    private void foregroundLoadTable(boolean with_annotation_columns)
//    {
//        String dialog_title = (with_annotation_columns?"Open annotated family size table":"Open family size table");
//        FancyFileDialog dialog  = new FancyFileDialog(this,dialog_title,FileDialog.LOAD);
//        dialog.setVisible(true);
//        
//        String file_name = dialog.getFile();
//        
//        if (file_name != null)
//        {
//        	Session sesh = getActiveSession();
//// BUNDLE        	
////        	String[] terminal_names = sesh.getTreesBrowser().getLeafNames();
//        	String[] terminal_names = sesh.getModelBrowser().getMainPhylogeny().getLeafNames();
//// BUNDLE        	
//            String directory = dialog.getDirectory();
//            try
//            {
//	        	AnnotatedTable table 
//	        	= TableParser.readTable(terminal_names, dialog.getReader(),
////	        			GeneralizedFileReader.guessReaderForInput(directory+file_name), 
//	        			with_annotation_columns);
//	            checkTableEmpty(table);
//	            DataFile<AnnotatedTable> table_data = new DataFile<>(table, new File(directory, file_name));
//	            
//	            sesh.addDataSet(table_data);
//	        } catch (java.io.InterruptedIOException E)
//	        {
//	            // Canceled: nothing to do
//	        } catch (java.io.IOException E)
//	        {
//	        	exception_handler.handle(E, "I/O error", "File error while reading occurrence table from a file.");
//	        } catch (IllegalArgumentException E)
//	        {
//	        	exception_handler.handle(E, "Parsing error", "Cannot parse occurrence table in file.");
//	        } catch (Exception E)
//	        {
//	        	exception_handler.handle(E, "A bug!", "Error while reading occurrence table from a file.");
//	        }
//        }
//    }
    
    private void doLoadRates()
    {
    	Session sesh = getActiveSession();
		Phylogeny main_tree = sesh.getModelBrowser().getSelectedTreeEntry().getTreeData().getContent();        				

		FancyFileDialog dialog = new FancyFileDialog(this,"Load rates",FileDialog.LOAD);
        dialog.setVisible(true);

        String file_name = dialog.getFile();
        if (file_name != null)
        {
            String directory = dialog.getDirectory();
            File rates_file = new File(directory,file_name);
            SwingWorker<DataFile<GammaInvariant>, Void> load_task = dialog.createTask(stream->
    			{                
    				BufferedReader ratesR = new BufferedReader(new InputStreamReader(stream));
    				GammaInvariant rate_model = RateVariationParser.readRates(ratesR, main_tree);
    				DataFile<GammaInvariant> rates_data = new DataFile<>(rate_model, rates_file);
    				return rates_data;
    			});
            load_task.addPropertyChangeListener(event->
    		{
    			if (! load_task.isCancelled()
    					&&
    					"state".equals(event.getPropertyName())
    	                 && SwingWorker.StateValue.DONE == event.getNewValue())   
    			{
    				DataFile<GammaInvariant> rates_data=null;
//    				System.out.println("#**AF.bIT task done");
    				try {rates_data = load_task.get();} 
    				catch (InterruptedException swallowed) {}
    				catch (ExecutionException swallowed) {}
    				if (rates_data!=null)
    				{
    					sesh.addRates(rates_data, true);    					
    				}
    			}
    		});
            load_task.execute();
        }
    }
    
//	private void foregroundLoadRates(BufferedReader R)
//    {
//        File rates_file = null;
//        if (R==null)
//        {
//            FancyFileDialog dialog = null;
//            dialog = new FancyFileDialog(this,"Load rates",FileDialog.LOAD);
//            dialog.setVisible(true);
//
//            String file_name = dialog.getFile();
//            String directory = dialog.getDirectory();
//            if (file_name != null)
//            {
//                try 
//                {
//                    rates_file = new File(directory,file_name);
////                    FileInputStream file_input = new FileInputStream(rates_file);
//                    R = dialog.getBufferedReader(); //  new InputStreamReader(file_input);
//                } catch (java.io.IOException E)
//                {
//                    exception_handler.handle(E, "I/O error", "Cannot open rates file.");
//                }
//            }
//        }
//        if (R != null)
//        {
//            GammaInvariant rate_model=null;
//        	Session sesh = getActiveSession();
//            try 
//            {
//// BUNDLE            	
////            	IndexedTree main_tree = sesh.getTreesBrowser().getSelectedTree().getTreeData().getContent();
//            	Phylogeny main_tree = sesh.getModelBrowser().getSelectedTreeEntry().getTreeData().getContent();
//// BUNDLE            	
//                rate_model = RateVariationParser.readRates(R, main_tree);
//                R.close();
//            } catch (java.io.InterruptedIOException E)
//            {
//                // canceled
//                rate_model = null;
//            } catch (RateVariationParser.FileFormatException E)
//            {
//                rate_model = null;
//                exception_handler.handle(E, "I/O error", "Badly formatted rates file.");
//            }
//            catch (java.io.IOException E)
//            {
//                rate_model = null;
//                exception_handler.handle(E, "I/O error", "File error while reading rates from a file.");
//            } catch (Exception E)
//            {
//                rate_model = null;
//                exception_handler.handle(E, "A bug maybe?", "Error while reading rates from a file.");
//            }
//            if (rate_model != null)
//            {
//                DataFile<GammaInvariant> rates_data = new DataFile<>(rate_model, rates_file);
//                sesh.addRates(rates_data, true);
//            } 
//        } // file_name 		
//    }
    
    private void checkTableEmpty(AnnotatedTable table)
    {
        String[] leaf_names = table.getTaxonNames();
        int has_zero = 0;
        StringBuilder zero_taxa = null;
        boolean has_spaces = false;
        for (int leaf_idx=0; leaf_idx<leaf_names.length; leaf_idx++)
        {
            int p = table.getNumFamiliesPresent(leaf_idx);
            if (p==0)
            {
                has_zero ++;
                has_spaces = has_spaces || (leaf_names[leaf_idx].indexOf(' ')!=-1);
                if (zero_taxa == null)
                {
                    zero_taxa = new StringBuilder();
                    zero_taxa.append("\"<tt>");
                    zero_taxa.append(leaf_names[leaf_idx]);
                    zero_taxa.append("</tt>\"");
                }
                else
                {
                    if (has_zero<4)
                    {
                        zero_taxa.append(", ");
                        zero_taxa.append("\"<tt>");
                        zero_taxa.append(leaf_names[leaf_idx]);
                        zero_taxa.append("</tt>\"");
                    }
                    else if (has_zero==4)
                        zero_taxa.append(",...");
                }
            }
        }
        if (has_zero!=0)
        {
            String msg = "<p><b>"+(has_zero==1?"Taxon ":"Taxa ")+zero_taxa.toString()+(has_zero==1?" has":" have");
            msg += " no members in any of the families.</b> <br />Maybe the names are misspelled in the tree.";
            if (has_spaces)
                msg += " <br />(Attention: in the Newick format, underscore is replaced by space unless the name is enclosed by quotation marks.)";
            msg += "</p>";
            JEditorPane EP = new JEditorPane("text/html",msg);
            EP.setEditable(false);
            JOptionPane.showMessageDialog(this,EP, "Missing data?", JOptionPane.WARNING_MESSAGE);
        }
    }
    
    
}
