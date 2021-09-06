/*
 * DealerCount.java
 *
 * Created on March 2, 2008, 4:10 PM
 *
 */

package ca.umontreal.iro.evolution.malin.ui.count;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Toolkit;
import java.awt.BorderLayout;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JDialog;
import javax.swing.JComponent;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.JSeparator;
import javax.swing.JScrollPane;
import javax.swing.KeyStroke;
import javax.swing.SwingWorker;
import javax.swing.UIManager;
import javax.swing.JPanel;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;

//import ca.umontreal.iro.banality.Verbose;

import ca.umontreal.iro.evolution.TreeNode;

import ca.umontreal.iro.evolution.malin.DataFile;

import ca.umontreal.iro.evolution.genecontent.NodeWithRates;
import ca.umontreal.iro.evolution.genecontent.OccurrenceTable;
import ca.umontreal.iro.evolution.genecontent.RateVariation;
import ca.umontreal.iro.evolution.genecontent.TreeWithRates;

import ca.umontreal.iro.evolution.malin.ui.Browser;
import ca.umontreal.iro.evolution.malin.ui.BrowserSelectionListener;
import ca.umontreal.iro.evolution.malin.ui.Dealer;
import ca.umontreal.iro.evolution.malin.ui.WorkSpace;
import ca.umontreal.iro.evolution.malin.ui.TreePanel;
//import org.xml.sax.SAXException;

/**
 *
 * Entry point for Count application.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class DealerCount extends Dealer
{
    private static final String EXECUTABLE_TITLE = "Count";
    private static final String EXECUTABLE_VERSION = "12.0127";

    @Override
    protected String getExecutableTitle()
    {
        return EXECUTABLE_TITLE;
    }
    
    @Override
    protected String getVersionNumber()
    {
        return EXECUTABLE_VERSION;
    }

    /**
     * Software package name + version.
     * 
     * @return software package name
     */
    public String getFullPackageName()
    {
        return getExecutableTitle()+" "+getVersionNumber();
    }
    
    public String getStandardHeader(Class C)
    {
        return "#| "+getFullPackageName()+"::"+C.getName();
    }

    public String getStandardHeader(String info)
    {
        return "#| "+info;
    
    }
    
    /** Creates a new instance of DealerCount 
     * @param frame th enclosing JFrame
     */
    public DealerCount(JFrame frame) 
    {
        super(frame);
        init();
    }

    private UsersGuide guide;
//    private static final String NO_SESSION_PANEL = "no session";

    private void init()
    {
        guide = new UsersGuide(top_frame);

        int w = top_frame.getWidth();
        int h = top_frame.getHeight();
        int x = top_frame.getX();
        int y = top_frame.getY();
        guide.setBounds(x+20,y+5,Math.min(w-40, 1000),600);

        doAbout();

//        JPanel empty_panel = new JPanel();
//        empty_panel.setLayout(new BoxLayout(empty_panel, BoxLayout.LINE_AXIS));
//        JLabel bats = new JLabel(this.getCountIcon());
//        empty_panel.setMaximumSize(new Dimension(w/10,h/10));
//        empty_panel.setPreferredSize(empty_panel.getMaximumSize());
//        bats.setMaximumSize(new Dimension(48,48));
//        bats.setPreferredSize(bats.getMaximumSize());
//        empty_panel.add(bats);
//        empty_panel.add(getAboutComponent());
//        JPanel enclosing = new JPanel(new BorderLayout());
//        enclosing.add(new JLabel("\u2b09  open a session or load a phylogeny."), BorderLayout.NORTH);
//        enclosing.add(empty_panel,BorderLayout.CENTER);
//        this.add(enclosing,NO_SESSION_PANEL);
//        show(NO_SESSION_PANEL);
    }

    
    /**
     * Allows access to this method within the same class.
     */
    @Override
    protected void initMenu()
    {
        super.initMenu();
    }

    private JMenuItem load_session;
    private JMenuItem save_session;
    
    @Override 
    protected void setSessionMenu()
    {
        if (load_session == null)
        {
            load_session = new JMenuItem("Open previously saved session(s) ...");
            load_session.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, 
                Toolkit.getDefaultToolkit().getMenuShortcutKeyMask())); 
            load_session.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {
                        doLoadSessions();
                    }
                });
            load_session.setToolTipText("This tooltip has no information");
            
            save_session = new JMenuItem("Save everything ...");
            save_session.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, 
                Toolkit.getDefaultToolkit().getMenuShortcutKeyMask())); 
            save_session.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {
                        doSaveSessions();
                    }
                }
            );
        }
        
        super.setSessionMenu();
    }
    
    @Override
    protected void addSessionMenuItems()
    {
        WorkSpaceCount wso = getActiveWorkSpaceCount();

        if (session_menu_selection_count != work_space_list.size())
        {

            session_menu.add(load_session);
            super.addSessionMenuItems();
            session_menu.addSeparator();
            session_menu.add(save_session);
        }

        save_session.setEnabled(wso!=null);
        load_session.setEnabled(work_space_list.isEmpty() || work_space_list.size()==0);
    }
    
    
    private JMenuItem data_load_table;
    private JMenuItem data_load_annotated;
    private JMenuItem data_load_annotation;
    private JMenuItem data_filter_lines;
    private JMenuItem data_binary;
    
    /**
     * Sets up the elements of the <q>Data</q> menu.
     */
    @Override
    protected void setDataMenu()
    {

        if (data_menu.getItemCount()==0) // initialization
        {


            data_load_table = new JMenuItem("Open table...");
            data_load_table.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {
                        doLoadDataFile(false);
                    }
                });
            data_load_table.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T,Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
            data_menu.add(data_load_table);       

            data_load_annotated = new JMenuItem("Open annotated table...");
            data_load_annotated.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {
                        doLoadDataFile(true);
                    }
                });
            data_menu.add(data_load_annotated);       

            data_load_annotation = new JMenuItem("Load family annotations...");
            data_load_annotation.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {
                        doLoadAnnotations();
                    }
                });
            data_menu.add(data_load_annotation);
            
            data_menu.add(new JSeparator());
         
            data_filter_lines = new JMenuItem("Extract selected families into a new table");
            data_filter_lines.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {
                        getActiveWorkSpaceCount().showFilteredFamilies();
                    }
                });

            data_menu.add(data_filter_lines);

            data_binary = new JMenuItem("Transform numerical profiles into binary (presence/absence) profiles ");
            data_binary.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {
                        getActiveWorkSpaceCount().showBinaryProfiles();
                    }
                });
            data_menu.add(data_binary);
        }
        
        // enabing/disabling menu items
        WorkSpaceCount wso = getActiveWorkSpaceCount();
        OccurrenceTable ot = null;
        File data_file=null;
        RateVariation rate_model = null;
        File model_file=null;
        if (wso!=null)
        {
            DataFile<OccurrenceTable> df = null;
            Object d0 = wso.getSelectedPrimaryDataDisplay();
            if (d0!=null && d0 instanceof FamilySizeTableDisplay)
                df = wso.getSelectedDataset();
            DataFile<RateVariation>          rf = null; // wso.getSelectedModel();
            if (df!=null)
            {
                ot = (OccurrenceTable) df.getData();
                data_file = df.getFile();
            }
            if (rf != null)
            {
            }
        }
        data_load_table.setEnabled(wso != null);
        data_load_annotated.setEnabled(wso != null);
        data_load_annotation.setEnabled(data_file!=null);
        if (wso==null)
            data_filter_lines.setEnabled(false);
        else if (WorkSpaceCount.DATA_PANEL.equals(wso.getSelectedWorkAreaElementName()))
        {
            Object d_selected = wso.getDatasetsBrowser().getSelectedItem();            
            data_filter_lines.setEnabled(d_selected instanceof FilterableFamilies);
        } else 
        {
            data_filter_lines.setEnabled(false);
            //System.out.println("#*DC.sDM filter '"+wso.getSelectedWorkAreaElementName()+"'");
        } 
        if (wso == null)
            data_binary.setEnabled(false);
        else
        {
            Object d0 = wso.getSelectedPrimaryDataDisplay();
            data_binary.setEnabled(d0!=null && (d0 instanceof FamilySizeTableDisplay));
        }
        
    }
    

    private JMenuItem rates_load;
    private JMenuItem rates_optimize;
    
    @Override
    protected void setRateMenu()
    {
        if (rate_menu.getItemCount()==0)
        {
            rates_load = new JMenuItem("Load rates...");
            rates_load.addActionListener(new ActionListener()
            {
                @Override
                public void actionPerformed(ActionEvent e)
                {
                    doLoadRates();
                }
            });
            rates_load.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_R,Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
            rate_menu.add(rates_load);
        
            rate_menu.addSeparator();
            rates_optimize = new JMenuItem("Optimize rates...");
            rates_optimize.addActionListener(new ActionListener()
            {
                @Override
                public void actionPerformed(ActionEvent e)
                {
                    getActiveWorkSpaceCount().optimizeRates(top_frame);
                }
            });
            rate_menu.add(rates_optimize);
        }

        WorkSpaceCount wso = getActiveWorkSpaceCount();
        OccurrenceTable ot = null;
        if (wso!=null)
        {
            DataFile<OccurrenceTable> df = null;
            Object d0 = wso.getSelectedPrimaryDataDisplay();
            if (d0!=null && d0 instanceof FamilySizeTableDisplay)
                df = wso.getSelectedDataset();
            if (df!=null)
                ot = (OccurrenceTable) df.getData();
        }

        rates_load.setEnabled(wso!=null);
        rates_optimize.setEnabled(ot!=null && !ot.hasMissingEntries());
    }
    
    private JMenuItem analysis_dollo;
    private JMenuItem analysis_wagner;
    private JMenuItem analysis_posterior;
    private JMenuItem analysis_PGL;
    
    @Override
    protected void setAnalysisMenu()
    {
        WorkSpaceCount wso = getActiveWorkSpaceCount();

        OccurrenceTable ot = null;
        DataFile<OccurrenceTable> df = null;
        DataFile<RateVariation> rf = null;
        RateVariation rate_model = null;
        if (wso!=null)
        {
            Object d0 = wso.getSelectedPrimaryDataDisplay();
            if (d0!=null && d0 instanceof FamilySizeTableDisplay)
                df = wso.getSelectedDataset();
            if (df!=null)
                ot = (OccurrenceTable) df.getData();
            rf = wso.getSelectedRateModel();
            if (rf != null)
                rate_model = (RateVariation) rf.getData();
        }

        if (analysis_menu.getItemCount()==0)
        {
            analysis_dollo = new JMenuItem("Family history by Dollo parsimony");
            analysis_menu.add(analysis_dollo);
            analysis_dollo.addActionListener(new java.awt.event.ActionListener() 
            {
                @Override
                public void actionPerformed(ActionEvent actionEvent) 
                {
                    getActiveWorkSpaceCount().showDollo();
                }
            });
        
            analysis_wagner = new JMenuItem("Family history by Wagner parsimony");
            analysis_menu.add(analysis_wagner);
            analysis_wagner.addActionListener(new java.awt.event.ActionListener() 
            {
                @Override
                public void actionPerformed(ActionEvent actionEvent) 
                {
                    getActiveWorkSpaceCount().showWagner();
                }
            });

            analysis_posterior = new JMenuItem("Family history by posterior probabilities");
            analysis_menu.add(analysis_posterior);
            analysis_posterior.addActionListener(new java.awt.event.ActionListener() 
            {
                @Override
                public void actionPerformed(ActionEvent actionEvent) 
                {
                    getActiveWorkSpaceCount().showPosteriors();
                }
            });    

            analysis_menu.addSeparator();
            analysis_PGL = new JMenuItem("PGL: propensity for gene loss (Krylov-Wolf-Rogozin-Koonin)");
            analysis_menu.add(analysis_PGL);
            analysis_PGL.addActionListener(new java.awt.event.ActionListener() 
            {
                @Override
                public void actionPerformed(ActionEvent actionEvent) 
                {
                    getActiveWorkSpaceCount().showPGL();
                }
            });
        }

        if (ot == null)
        {
            analysis_dollo.setEnabled(false);
            analysis_dollo.setToolTipText("A table must be selected in the Data browser to enable Dollo parsimony");
            analysis_wagner.setEnabled(false);
            analysis_wagner.setToolTipText("A table must be selected in the Data browser to enable Wagner parsimony");
            analysis_PGL.setEnabled(false);
            analysis_PGL.setToolTipText("A table must be selected in the Data browser to compute PGL");
        } else
        {
            analysis_dollo.setEnabled(true);
            analysis_dollo.setToolTipText("Dollo parsimony on table "+df.getFile().getName());
            analysis_wagner.setEnabled(!ot.hasMissingEntries());
            analysis_wagner.setToolTipText("Wagner parsimony on table "+df.getFile().getName());
            analysis_PGL.setEnabled(true);
            analysis_PGL.setToolTipText("PGL computation on table "+df.getFile().getName());
        }

        if (ot==null || rate_model == null)
        {
            analysis_posterior.setEnabled(false);
            if (ot==null)
            {
                if (rate_model==null)
                    analysis_posterior.setToolTipText("A table must be selected in the Data browser and a rate model in the Rates browser");
                else
                    analysis_posterior.setToolTipText("A rate model must be selected in the Rates browser");
            } else
                analysis_posterior.setToolTipText("A table must be selected in the Data browser");
        } else
        {
            analysis_posterior.setEnabled(!ot.hasMissingEntries());
            analysis_posterior.setToolTipText("Posterior reconstruction on table "+df.getFile().getName()
                    +" and rates "+rf.getFile().getName());
        }
        
        
        //System.out.println("#*DC.sAM pgl "+PGL.isEnabled());
        
        //
        //JMenuItem GLR = new JMenuItem("GLR: gene loss rate (Borenstein-Shlomi-Ruppin-Sharan)");
        //GLR.setEnabled(false);
        //analysis_menu.add(GLR);
        
    }
    
    @Override
    protected void setHelpMenu()
    {
        if (help_menu.getItemCount()==0)
        {
            super.setHelpMenu();

            JMenuItem guide = new JMenuItem("User's guide");
            help_menu.add(guide);

            guide.addActionListener(new ActionListener()
                {
                    @Override
                    public void actionPerformed(ActionEvent e)
                    {
                        doHelp();
                    }
            });
        }
    }
    
    private void doHelp()
    {
//        JEditorPane EP = new JEditorPane("text/html","<p><b>User's guide is not available in this version.</b></p>" +
//                "<p>Please consult the PDF version distributed with this package.</p>");
//        EP.setEditable(false);
//        EP.setBackground(top_frame.getBackground());
//        EP.setPreferredSize(new Dimension(200,100));
//        JScrollPane ep_scroll = new JScrollPane(EP);
//        JOptionPane.showMessageDialog(top_frame, ep_scroll, "User's guide", JOptionPane.INFORMATION_MESSAGE,
//                getCountIcon());
        Component activeC = null;
        WorkSpaceCount wsc = getActiveWorkSpaceCount();
        if (wsc != null)
            activeC = wsc.getDisplayedComponent();
        if (activeC != null)
        {
            if (activeC instanceof TreePanel)
            {
                guide.openPage(UsersGuide.SUBSECTION_PHYLOGENY_T);
            } else if (activeC instanceof Browser)
            {
                guide.openPage(UsersGuide.SUBSECTION_BROWSERS_T);
            } else if (activeC instanceof FamilySizeTableDisplay)
            {
                guide.openPage(UsersGuide.SUBSECTION_TABLE_T);
            } else if (activeC instanceof RateModelDisplay)
            {
                guide.openPage(UsersGuide.SUBSECTION_RATES_PANEL_T);
            } else if (activeC instanceof PGLDisplay)
            {
                guide.openPage(UsersGuide.SUBSECTION_ANALYSIS_PGL_T);
            } else if (activeC instanceof DolloDisplay)
            {
                guide.openPage(UsersGuide.SUBSECTION_ANALYSIS_DOLLO_T);
            } else if (activeC instanceof WagnerDisplay)
            {
                guide.openPage(UsersGuide.SUBSECTION_ANALYSIS_WAGNER_T);
            } else if (activeC instanceof PosteriorDisplay)
            {
                guide.openPage(UsersGuide.SUBSECTION_ANALYSIS_POSTERIORS_T);
//            } else if (activeC instanceof AlignmentSetDisplay)
//            {
//                guide.openPage(guide.SECTION_ALIGNMENTS_T);
//            } else if (activeC instanceof BootstrapRun)
//            {
//                guide.openPage(guide.SUBSECTION_ANALYSIS_BOOTSTRAP_T);
//            } else if (activeC instanceof HistoryDisplay)
//            {
//                guide.openPage(guide.SUBSECTION_ANALYSIS_SITES_T);
            } //else
              //  System.out.println("#**DI.dH called "+activeC.getClass().getName());

        }
        guide.setVisible(true);
    }
    
    @Override
    protected WorkSpace newWorkSpace(TreeNode root, File F)
    {
        
        WorkSpaceCount ws = new WorkSpaceCount(this, root, F);
        
        ws.getRatesBrowser().addBrowserSelectionListener(new BrowserSelectionListener() 
        {
            @Override
            public void valueChanged(Browser.PrimaryItemSelectionEvent E) {
                //System.out.println("#*DC.nWS.vC rates browser selection changed "+E);
                initMenu();
            }
        });
        ws.getDatasetsBrowser().addBrowserSelectionListener(new BrowserSelectionListener() 
        {
            @Override
            public void valueChanged(Browser.PrimaryItemSelectionEvent E) {
                initMenu();
            }
        });
        
        ws.addChangeListener(new ChangeListener() // when tab selection changes...
        {
            @Override
            public void stateChanged(ChangeEvent E)
            {
                initMenu();
                //System.out.println("#**DC.nWS.CL tab selection changed "+E);
            }
            
        });

        return ws;
    }
    
    private void doLoadDataFile(boolean is_annotated)
    {
        FileDialog dialog = null;
        String dialog_title = (is_annotated?"Open annotated family size table":"Open family size table");
        dialog = new FileDialog(top_frame,dialog_title,FileDialog.LOAD);
        dialog.setVisible(true);
       
        String file_name = dialog.getFile();
        String directory = dialog.getDirectory();
        OccurrenceTable ot = null;
        if (file_name != null)
        {
            
            try 
            {
                TreeWithRates TR = new TreeWithRates((NodeWithRates)getActiveWorkSpace().getPhylogeny());
                ot = new OccurrenceTable(TR.getLeaves());
                InputStreamReader R = new InputStreamReader(new FileInputStream(directory+file_name));
                ot.readTable(R,is_annotated);
                R.close();
                if (ot != null)
                {
                    // check table columns
                    {
                        String[] leaf_names = ot.getTerminalTaxonNames();
                        int[] num_families_present = new int[leaf_names.length];
                        int has_zero = 0;
                        StringBuffer zero_taxa = null;
                        boolean has_spaces = false;
                        for (int leaf_idx=0; leaf_idx<leaf_names.length; leaf_idx++)
                        {
                            int p = ot.getNumFamiliesPresent(leaf_idx);
                            num_families_present[leaf_idx] = p;
                            if (p==0)
                            {
                                has_zero ++;
                                has_spaces = has_spaces || (leaf_names[leaf_idx].indexOf(' ')!=-1);
                                if (zero_taxa == null)
                                {
                                    zero_taxa = new StringBuffer();
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
                            JOptionPane.showMessageDialog(top_frame,EP, "Missing data?", JOptionPane.WARNING_MESSAGE);
                        }
                    }
                    getActiveWorkSpaceCount().addDataSet(ot, new File(directory, file_name));
                }
            } catch (java.io.InterruptedIOException E)
            {
                // Cancel: nothing to do 
                ot = null;
            } catch (java.io.IOException E)
            {
                ot = null;
                exceptionCaught(E, "I/O error", "File error while reading occurrence table from a file.");
            } catch (IllegalArgumentException E)
            {
                ot = null;
                exceptionCaught(E, "Parsing error", "Cannot parse occurrence table in file.");
            } catch (Exception E)
            {
                ot = null;
                exceptionCaught(E, "A bug!", "Error while reading occurrence table from a file.");
            }
        }
    }
    
    private void doLoadAnnotations()
    {
        FamilySizeTableDisplay TD = getActiveWorkSpaceCount().getSelectedRootTableDisplay();//getSelectedTableDisplay();
        TD.loadAnnotations(top_frame);
    }
    
    
    private void doLoadRates()
    {
        FileDialog dialog = null;
        dialog = new FileDialog(top_frame,"Load rates",FileDialog.LOAD);
        dialog.setVisible(true);
       
        String file_name = dialog.getFile();
        String directory = dialog.getDirectory();
        
        RateVariation rate_model = null;
        if (file_name != null)
        {
            try 
            {
                File selected_file = new File(directory,file_name);
                FileInputStream file_input = new FileInputStream(selected_file);
                
                TreeWithRates TR = new TreeWithRates((NodeWithRates)getActiveWorkSpace().getPhylogeny());
                InputStreamReader R = new InputStreamReader(file_input);
                
                rate_model = RateVariation.read(R, TR);
                
                R.close();
            } catch (java.io.InterruptedIOException E)
            {
                // canceled
                rate_model = null;
            } catch (java.io.IOException E)
            {
                rate_model = null;
                exceptionCaught(E, "I/O error", "File error while reading rates from a file.");
            } catch (Exception E)
            {
                rate_model = null;
                exceptionCaught(E, "A bug maybe?", "Error while reading rates from a file.");
            }
            
            
            if (rate_model != null)
            {
                getActiveWorkSpaceCount().addRates(rate_model, new File(directory,file_name));
            } 
        } // file_name 
        
    }
    
    private void doSaveSessions()
    {
        FileDialog dialog = null;
        dialog = new FileDialog(top_frame,"Save all work",FileDialog.SAVE);
        dialog.setVisible(true);
       
        final String file_name = dialog.getFile();
        final String directory = dialog.getDirectory();
        
        if (file_name != null)
        {
            
            final JDialog wait_until_saved = new JDialog(this.top_frame,"Saving sessions",true);
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
                        out.println("<?xml version=\"1.0\" encoding=\"us-ascii\" ?>");
                        out.println("<"+WorkSpaceCount.EMT_TOP+" " +
                                "activesession=\""+"S"+Integer.toString(active_work_space)+"\" " +
                                ">");

                        for (int session_idx=0; session_idx<work_space_list.size(); session_idx++)
                        {
                            WorkSpaceCount ws = (WorkSpaceCount) work_space_list.get(session_idx);
                            String session_identifier = "S"+Integer.toString(session_idx);
                            ws.saveSession(out, session_identifier);
                        }

                        out.println("</"+WorkSpaceCount.EMT_TOP+">");
                        out.close();
                    } catch (java.io.IOException E)
                    {
                        exceptionCaught(E, "I/O error", "File error while saving sessions.");
                    } catch (Exception E)
                    {
                        exceptionCaught(E, "A bug maybe?", "Error while saving the sessions.");
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
            java.awt.Dimension frameD = top_frame.getSize();
            wait_until_saved.setBounds(frameD.width/3, frameD.height/3, frameD.width/3, frameD.height/6);
            wait_until_saved.setVisible(true);            
            
            
        } // file_name 
    }
    
    private void doLoadSessions()
    {
        FileDialog dialog = null;
        dialog = new FileDialog(top_frame,"Load sessions",FileDialog.LOAD);
        dialog.setVisible(true);
       
        final String file_name = dialog.getFile();
        final String directory = dialog.getDirectory();
        
        if (file_name != null)
        {
            final JDialog wait_until_loaded = new JDialog(this.top_frame,"Loading sessions",true);
            BoxLayout load_layout = new BoxLayout(wait_until_loaded.getContentPane(),BoxLayout.PAGE_AXIS);
            wait_until_loaded.setLayout(load_layout);
            wait_until_loaded.add(Box.createVerticalGlue());
            wait_until_loaded.add(new JLabel("Please be patient until session data load ..."));
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
                        WorkSpaceCount.loadSessions(selected_file, DealerCount.this);
                    } catch (java.io.IOException E)
                    {
                        exceptionCaught(E, "I/O error", "File error while loading sessions.");
                    } catch(org.xml.sax.SAXException E)
                    {
                        exceptionCaught(E, "File format problem", "The file is not formatted correctly");
                    } catch (Exception E)
                    {
                        exceptionCaught(E, "A bug maybe?", "Error while loading sessions.");
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
            java.awt.Dimension frameD = top_frame.getSize();
            wait_until_loaded.setBounds(frameD.width/3, frameD.height/3, frameD.width/3, frameD.height/6);
            wait_until_loaded.setVisible(true);

            //this.active_work_space = this.work_space_list.size()-1;
            initMenu();
        } // file_name         
    }
    
           ;
    
    private static String ABOUT_TEXT 
            = "<h1>The "+EXECUTABLE_TITLE+" v"+EXECUTABLE_VERSION+"</h1>" +
            "<p>The Count is a software package " +
            "for the evolutionary analysis of phylogenetic profiles, and numerical characters in general, " +
            "written entirely in Java. " +
            "</p>" +
            "<p>Author: Mikl&oacute;s Cs&#369;r&ouml;s http://www.iro.umontreal.ca/~csuros/</p> " +
            "<p>Some numerical optimization routines were adapted from " +
            "<em>Numerical Recipes in C: The Art of Scientific Computing</em> " +
            "[W. H. Press, S. A. Teukolsky, W. V. Vetterling and B. P. Flannery; " +
            "Second edition, Cambridge University Press, 1997].</p>" +
            "<p>Some code for MacOS X integration is based on Eirik Bjorsnos' Java package <tt>org.simplericity.macify.eawt</tt>, " +
            "distributed under the terms of the Apache License http://www.apache.org/licenses/LICENSE-2.0</p>" +
            "<p>On systems other than Mac OS and Windows, the JGoodies Look&amp;Feel package is used http://www.jgoodies.com/</p>" +
            "<p>The background image of the Count logo (Moon and Mars) is used with permission from the copyright owner, John Harms.</p>"+
            
            "</p>" +
            //"<p>Citing Malin in your work: "+UsersGuide.MALIN_REFERENCE+"</p>"+
            "<p>Algorithmic ideas uderlying the Count package were described in the following publications.</p>" +
            "<ul>" +
            UsersGuide.METHOD_REFERENCES+
            "</ul>" +
            "<p>Montr&eacute;al/Amsterdam/Szentendre, 2010 April 2</p>";


    private JComponent getAboutComponent()
    {
        JEditorPane EP = new JEditorPane("text/html",ABOUT_TEXT);
        EP.setEditable(false);
        EP.setBackground(top_frame.getBackground());
        EP.setPreferredSize(new Dimension(600,600));
        JScrollPane ep_scroll = new JScrollPane(EP);
        return ep_scroll;
    }

    @Override
    protected void doAbout()
    {
        JOptionPane.showMessageDialog(top_frame, getAboutComponent(), "About Count", JOptionPane.INFORMATION_MESSAGE,
                getCountIcon());
    }
    
    public WorkSpaceCount getActiveWorkSpaceCount()
    {
        return (WorkSpaceCount)super.getActiveWorkSpace();
    }
    
    public static ImageIcon getCountIcon()
    {
        java.net.URL icon_url = ClassLoader.getSystemResource("img/count-icon.jpg");
        try 
        {
            if (icon_url == null)
            {
                icon_url = new java.net.URL("http://www.iro.umontreal.ca/~csuros/gene_content/images/count-icon.jpg");
            }
        
            return new ImageIcon(icon_url);
        } catch (Exception E)
        {
            return null;
        }
    }
    
    
    public static void main(String[] args)
    {
        //Verbose.setVerbose(true);
        try
        {
            String system = System.getProperty("os.name", "[unknown OS]");

            //if (system.startsWith("Windows"))
            //{
            //    javax.swing.UIManager.setLookAndFeel("com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
            //} else 
            if (system.startsWith("Mac OS"))
            {
                // MacOS X 10.4+ --- Java 1.5 runtime properties; no effect on other systems 
                System.setProperty("apple.laf.useScreenMenuBar", "true");
                System.setProperty("com.apple.mrj.application.apple.menu.about.name", EXECUTABLE_TITLE);
                System.setProperty("apple.awt.showGrowBox","true");
                //System.setProperty("apple.awt.brushMetalLook","true");
            } 
            //else
            //{
            //    //com.jgoodies.looks.plastic.PlasticXPLookAndFeel.setPlasticTheme(new com.jgoodies.looks.plastic.theme.ExperienceBlue());
            //    javax.swing.UIManager.setLookAndFeel("com.jgoodies.looks.plastic.PlasticXPLookAndFeel");
            //}
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception E) // ClassNotFoundException, InstantiationException
        {
            System.out.println("Ooops - no L&F");
            E.printStackTrace();
        }

        
        javax.swing.SwingUtilities.invokeLater(new Runnable() 
            {
                @Override
                public void run() 
                {
                    createAndShowGUI();
                }
            });
    }
    
    
    private static void createAndShowGUI()
    {
        JFrame top = new JFrame();

        top.setIconImage(getCountIcon().getImage());        

        top.setBounds(25,25,
            (int)Toolkit.getDefaultToolkit().getScreenSize().getWidth()-50,
            (int)Toolkit.getDefaultToolkit().getScreenSize().getHeight()-50);
        
        top.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        DealerCount D=new DealerCount(top);
        D.makeThreadExceptionHandler(Thread.currentThread());
        D.setVisible(true);
    }    
    
}


