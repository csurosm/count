package count.gui;
/*
 * Copyright 2023 Mikl&oacute;s Cs&#369;r&ouml;s.
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
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;


import java.io.File;
import java.io.PrintStream;

import java.util.HashMap;
import java.util.Map;
import java.util.Properties;

import java.text.DecimalFormat;
import java.text.NumberFormat;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.InputVerifier;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;


import javax.swing.text.JTextComponent;

import count.Count;
import count.ds.AnnotatedTable;
import count.ds.Phylogeny;
//import count.gui.kit.FancyFileDialog;
import count.gui.kit.ShmancyFileDialog;
import count.io.CommandLine;
import count.io.DataFile;

import count.model.MixedRateModel;

/**
 * 
 * A modal dialog for constructing command-line execution.
 *
 */
public class CommandBuilderDialog extends JDialog
{
	public CommandBuilderDialog(AppFrame frame, Class<?> main_class, boolean save_tree, boolean want_rates)
	{
    	super(frame, "Command-line launch builder for "+Count.APP_TITLE+" ("+main_class.getCanonicalName()+")", true); // modal
		this.app = frame;
		this.command_main_class = main_class;
		initComponents(save_tree, want_rates);
		pack();
		this.setLocationRelativeTo(null);	
	}
	
	private final AppFrame app;
	private final Class<?> command_main_class;
    
    private static final int COLUMNS_FILE = 30;
    private static final int COLUMNS_SMALLINT = 4;
    private static final int COLUMNS_FLOAT = 16;
    private static final int COLUMNS_STRING = 12;
    private static final int COLUMNS_LONGSTRING = 50;
    private static final int COLUMNS_DOUBLE = COLUMNS_FLOAT;

    private static final Color COLOR_DECORATION = new Color(255,255,204);
    
    
    /*
     * Status when last dialog is disposed of.
     */
    public static final int CANCEL_OPTION = 0;
    public static final int SAVE_OPTION = 1;
    public static final int COPY_OPTION = 2;
    public static final int EXEC_OPTION = 3;
    public static final int PRINT_OPTION = 4;
    
    /**
     * Types for formatting and input verification.
     */
    public static enum FieldType
    {
        FILE, SMALL_INT, FLOAT, DOUBLE, STRING, LONGSTRING;
    }

    private static void initTextField(JFormattedTextField field, FieldType kind)
    {
        switch(kind)
        {
            case SMALL_INT:
                {
                    field.setColumns(COLUMNS_SMALLINT);
                    field.setToolTipText("An integer value");
                    field.setHorizontalAlignment(JTextField.LEFT);
                    field.setValue(0);

                    break;
                }
            case FLOAT:
            case DOUBLE:
                {

                    field.setColumns(kind==FieldType.FLOAT?COLUMNS_FLOAT:COLUMNS_DOUBLE);
                    field.setToolTipText("An integer or floating-point value"
                    		+(kind==FieldType.DOUBLE?", in scientific notation":""));
                    field.setHorizontalAlignment(JTextField.LEFT);
                    field.setValue(0.0);

                    break;
                }
            case STRING:
            case LONGSTRING:
            {
                field.setColumns(kind == FieldType.STRING?COLUMNS_STRING:COLUMNS_LONGSTRING);
                field.setToolTipText("");
                field.setHorizontalAlignment(JTextField.RIGHT);
                field.setValue("");
                break;
            }
            case FILE:
            {
                field.setColumns(COLUMNS_FILE);
                field.setToolTipText("A file");
                field.setHorizontalAlignment(JTextField.RIGHT);
                field.setValue(new File((File)null, ""));
                break;
            }
        }
        field.setEditable(true);
        field.setMaximumSize(field.getPreferredSize());

        InputVerifier verifier = new InputVerifier()
        {
            @Override
            public boolean verify(JComponent o)
            {
                JFormattedTextField f = (JFormattedTextField)o;
                JFormattedTextField.AbstractFormatter formatter = f.getFormatter();
                if (formatter != null) // ugly way of exploiting parse exceptions
                {
                    String text = f.getText();
                    try
                    {
                        formatter.stringToValue(text);
                        return true;
                    } catch (java.text.ParseException pe)
                    {
                        return false;
                    }
                }
                return true;
            }
            @Override
            public boolean shouldYieldFocus(JComponent input)
            {
                boolean yield = verify(input);
                return yield;
            }
        };
        field.setInputVerifier(verifier);
        field.setFocusLostBehavior(JFormattedTextField.COMMIT_OR_REVERT);

        ActionListener skippyskip = new ActionListener()
        {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                JFormattedTextField f = (JFormattedTextField) e.getSource();
                f.transferFocus();
            }
        };
        field.addActionListener(skippyskip); // when user presses Enter, go to next element
    }

    /**
     * Instantiates a text field with standard formats.
     * 
     * @param kind
     * @return
     */
    public static JFormattedTextField createTextField(FieldType kind)
    {
        JFormattedTextField field = null;

        switch(kind)
        {
            case SMALL_INT:
                {
                    NumberFormat format = NumberFormat.getIntegerInstance();
                    field = new JFormattedTextField(format);
                    break;
                }
            case FLOAT:
                {
                    NumberFormat format = new DecimalFormat(); //NumberFormat.getNumberInstance();
                    format.setMaximumFractionDigits(15);
                    field = new JFormattedTextField(format);
                    break;
                }
            case DOUBLE:
            	{
            		NumberFormat format = new DecimalFormat("0.############E0");
            		field = new JFormattedTextField(format);
            	}
            case STRING:
            case LONGSTRING:
            case FILE:
            {
                field = new JFormattedTextField();
                break;
            }
        }

        initTextField(field, kind);

        return field;
    }
    
    /**
     * Embeds the component in a scroll pane.
     * @param comp component to be scrolled 
     * @return scrollpane embedding comp
     */
    public static JScrollPane scrolled(JComponent comp)
    {
        JScrollPane scroll = new JScrollPane(comp);
        scroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        scroll.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        return scroll;
    }
    /**
     * Encloses a component in a box with padding on the right. 
     * @param component to be boxed
     * @return a box with the component and glue 
     */
    public static Box boxed(JComponent component)
    {
    	Box box = new Box(BoxLayout.LINE_AXIS);
    	box.add(component);
    	box.add(Box.createHorizontalGlue());
    	return box;
    }
    /**
     * Encloses 2 components in a box with padding on the right. 
     * @param comp1 to be boxed
     * @param comp2 second components
     * @return a box with the components and glue 
     */
    public static Box boxed(JComponent comp1, JComponent comp2)
    {
    	Box box = new Box(BoxLayout.LINE_AXIS);
    	box.add(comp1);
    	box.add(comp2);
    	box.add(Box.createHorizontalGlue());
    	return box;
    }

    private int closing_status = CANCEL_OPTION;
    
    private JFormattedTextField mainclass_TF;
    private JFormattedTextField javaexec_TF;
    private JFormattedTextField javathreads_TF;
    private JFormattedTextField options_TF; 
    
//    private FancyFileDialog jar_dialog;
    private JFormattedTextField jar_TF;
    
    private ShmancyFileDialog table_dialog;
    private ShmancyFileDialog tree_dialog;
    private ShmancyFileDialog rates_dialog;
    
    private ShmancyFileDialog load_bundle_dialog;
    
    private ShmancyFileDialog output_dialog;
    private ShmancyFileDialog output_tree_dialog;
    private ShmancyFileDialog save_bundle_dialog;
    
    private Properties extra_options;
    private Map<String, String> option_explanations;
    
    /**
     * Closing status: one of {@link #CANCEL_OPTION}, {@link #SAVE_OPTION}, {@link #COPY_OPTION}, {@link #EXEC_OPTION}.
     * @return
     */
    public int getClosedReason() { return closing_status;}
    private void initComponents(boolean save_tree, boolean want_rates)
	{
		extra_options= new Properties();
		option_explanations = new HashMap<>();

		Box front_box = new Box(BoxLayout.PAGE_AXIS);
		front_box.setAlignmentX(Box.LEFT_ALIGNMENT);
		
		front_box.add(initRuntimeComponent());
		front_box.add(initInputComponent(save_tree, want_rates));
		front_box.add(initOutputComponent(save_tree));
		
		Box button_box = new Box(BoxLayout.LINE_AXIS);
//		addTitledBorder(button_box, "Let's get to work");
        JButton seeB = new JButton("Show ▶︎");
        seeB.setToolTipText("The corresponding command-line arguments or shell script will be shown in a popup window");
		JButton cancelB = new JButton("Cancel");
		
		CommandDialog command_dialog = new CommandDialog();
		seeB.addActionListener(click->
		{
			dispose();
			command_dialog.command_panel.showCommands();
			command_dialog.pack();
			command_dialog.setVisible(true);
		});
		
		cancelB.addActionListener(click ->
		{
			dispose();
			closing_status = CANCEL_OPTION;
		});
		button_box.add(Box.createHorizontalGlue());
		button_box.add(cancelB);
		button_box.add(seeB);
		button_box.setBackground(COLOR_DECORATION);
		button_box.setOpaque(true);
		
		
		this.add(front_box, BorderLayout.CENTER);
		this.add(button_box, BorderLayout.PAGE_END);
	}
    
	/**
	 * Adds a titled border around the component with our preferred color {@link #COLOR_DECORATION}
	 * 
	 * @param comp component
	 * @param title border title
	 */
    private static void addTitledBorder(JComponent comp, String title)
    {
        javax.swing.border.TitledBorder titled = BorderFactory.createTitledBorder(BorderFactory.createLineBorder(COLOR_DECORATION), title);
        titled.setTitleFont(new Font("Serif",Font.BOLD,titled.getTitleFont().getSize()));
        comp.setBorder(titled);
    }
    
     
    
//    private JFormattedTextField jar_TF;
    
    private JComponent initRuntimeComponent()
    {
        Box runtime_panel = new Box(BoxLayout.PAGE_AXIS);
        runtime_panel.setAlignmentX(Box.LEFT_ALIGNMENT);
//        runtime_panel.setLayout(new BoxLayout(runtime_panel, BoxLayout.PAGE_AXIS));
        addTitledBorder(runtime_panel, "Runtime options");

        JLabel javaexecL = new JLabel("Java executable:");
        javaexec_TF = createTextField(FieldType.FILE);
        javaexec_TF.setValue(new File((File)null,"java"));
        javaexec_TF.setToolTipText("Path to java executable");
        javaexecL.setLabelFor(javaexec_TF);
        runtime_panel.add(boxed(javaexecL, javaexec_TF));
        
        
        JLabel jarL = new JLabel("Classpath:");
        jar_TF = createTextField(FieldType.LONGSTRING); // instead of FILE
        jarL.setLabelFor(jar_TF);
        jar_TF.setValue(new File((File)null,System.getProperty("java.class.path")));
        jar_TF.setToolTipText("Path to JAR bundle containing this executable.");
//        jar_dialog = new FancyFileDialog(app, "Where is the JAR file?","JAR file");
//        jar_dialog.setFile(System.getProperty("java.class.path"));
        runtime_panel.add(boxed(jarL,jar_TF));
        
//        javaB.add(jarL);
//        javaB.add(jar_TF);
//        
//        
//        javaB.add(Box.createHorizontalGlue());
        
//        
//        runtime_panel.add(jar_dialog.getEmbeddableChooser());
        

        Box mainB = new Box(BoxLayout.LINE_AXIS);
        JLabel mainclassL = new JLabel("Main class:");
        mainclass_TF = createTextField(FieldType.LONGSTRING);
        mainclass_TF.setEditable(false);
        mainclass_TF.setValue(command_main_class.getCanonicalName());
        mainclassL.setLabelFor(mainclass_TF);
        mainB.add(mainclassL);
        mainB.add(mainclass_TF);
        
        JLabel javathreadsL = new JLabel("Worker threads:");
        javathreads_TF = createTextField(FieldType.SMALL_INT);
        javathreads_TF.setValue(Count.THREAD_PARALLELISM);
        javathreads_TF.setToolTipText("Number of worker threads used in model optimization. 1 (and less) disables multithreading; you have maximum "+Runtime.getRuntime().availableProcessors()+" on this machine.");

        if (Count.UsesThreadpool.class.isAssignableFrom(command_main_class)) 
        {
	        mainB.add(javathreadsL);
	        mainB.add(javathreads_TF);
    	}
        mainB.add(Box.createHorizontalGlue());
        runtime_panel.add(mainB);
        
        JLabel optionsL = new JLabel("Options:");
        options_TF = createTextField(FieldType.LONGSTRING);
        options_TF.setToolTipText("Optional arguments passed to main class");
        options_TF.setAlignmentX(JFormattedTextField.LEFT_ALIGNMENT);
    	optionsL.setLabelFor(options_TF);
        runtime_panel.add(boxed(optionsL,options_TF));
        
        
        return runtime_panel;
    }
    
    private JComponent initInputComponent(boolean save_tree, boolean want_rates)
    {
        Box input_panel = new Box(BoxLayout.PAGE_AXIS);
        input_panel.setAlignmentX(Box.LEFT_ALIGNMENT);
//        input_panel.setLayout(new BoxLayout(input_panel, BoxLayout.PAGE_AXIS));
        addTitledBorder(input_panel, "Input");
        tree_dialog = new ShmancyFileDialog(app, "Input phylogeny", "Tree file:");
        input_panel.add(boxed(tree_dialog.getEmbeddableChooser()));
        table_dialog = new ShmancyFileDialog(app, "Input data table", "Table file:");
        input_panel.add(boxed(table_dialog.getEmbeddableChooser()));
        if (want_rates)
        {
	        rates_dialog = new ShmancyFileDialog(app, "Input rates model", "Rates file:");
	        input_panel.add(boxed(rates_dialog.getEmbeddableChooser()));
        } else
        {
        	rates_dialog = null;
        }
        if (save_tree)
        {
	        load_bundle_dialog = new ShmancyFileDialog(app, "Load session", "Session file (XML):");
	        input_panel.add(boxed(load_bundle_dialog.getEmbeddableChooser()));
        } else
        {
        	load_bundle_dialog = null;
        }
        return input_panel;
    }

   
    private JComponent initOutputComponent(boolean save_tree)
    {
    	Box output_panel = new Box(BoxLayout.PAGE_AXIS);
    	output_panel.setAlignmentX(Box.LEFT_ALIGNMENT);
    	addTitledBorder(output_panel, "Output");
//    	output_panel.setLayout(new BoxLayout(output_panel, BoxLayout.PAGE_AXIS));
    	output_dialog = new ShmancyFileDialog(app, "Save output", FileDialog.SAVE, "Output file:");
    	output_panel.add(boxed(output_dialog.getEmbeddableChooser()));
    	if (save_tree)
    	{
    		output_tree_dialog = new ShmancyFileDialog(app, "Save phylogeny", FileDialog.SAVE, "Tree output:");
    		output_panel.add(boxed(output_tree_dialog.getEmbeddableChooser()));
        	save_bundle_dialog = new ShmancyFileDialog(app, "Save all models", FileDialog.SAVE, "Save as session:");
        	output_panel.add(boxed(save_bundle_dialog.getEmbeddableChooser()));
    	} else
    	{
    		output_tree_dialog = null;
    		save_bundle_dialog = null;
    	}
    	return output_panel;
    }
    
    public void setTableData(DataFile<AnnotatedTable> table_data)
    {
    	table_dialog.setFile(table_data==null?".":table_data.getFile().getPath());
    }
    public void setRatesData(DataFile<MixedRateModel> rates_data)
    {
    	rates_dialog.setFile(rates_data==null?".":rates_data.getFile().getPath());
    }
    public void setTreeData(DataFile<Phylogeny> tree_data)
    {
    	tree_dialog.setFile(tree_data==null?".":tree_data.getFile().getPath());
    }
    
    public void addOption(String option_name, String option_value, String explanation)
    {
    	extra_options.setProperty(option_name,option_value);
    	option_explanations.put(option_name, explanation);
        StringBuilder options_sb = new StringBuilder();
    	for (Object opt: extra_options.keySet())
    	{
    		String key = (String) opt;
    		String value = extra_options.getProperty(key);
    		if (0<options_sb.length())
    			options_sb.append(" ");
    		options_sb.append("-").append(key).append(" ").append(value);
    	}
    	options_TF.setText(options_sb.toString());
    	options_TF.setCaretPosition(0);
    }
    
    private static final String TAB_CMD = "Command line";
    private static final String TAB_SCRIPT = "Script";

    /**
     * Class for showing the command-line arguments. 
     * 
     * @author csuros
     *
     */
    private class CommandPanel extends JTabbedPane //implements ActionListener, PropertyChangeListener
    {
        private JTextArea command_line_text;
        private JTextArea command_script_text;

        private CommandPanel()
        {
            super();
            setTabPlacement(JTabbedPane.BOTTOM);
            line_sb = new StringBuilder();
            preamble_sb = new StringBuilder();
            script_sb=new StringBuilder();
            initComponents();
            showCommands();
            //turnOn();
        }

        private StringBuilder line_sb;
        private StringBuilder preamble_sb;
        private StringBuilder script_sb;
        
        public JTextComponent getSelectedTextComponent()
        {
            showCommands();
            int selected_idx = this.getSelectedIndex();
            String tab_name = this.getTitleAt(selected_idx);
            if (TAB_CMD.equals(tab_name))
                return command_line_text;
            else
                return command_script_text;
        }

        private void initComponents()
        {
            JPanel command_line_panel = new JPanel();
            command_line_text = new JTextArea(10,200);
            command_line_text.setFont(new Font("Monospaced",Font.PLAIN,12));
            command_line_text.setEditable(false);
            command_line_text.setLineWrap(true);
            command_line_text.setWrapStyleWord(true);
            command_line_panel.add(command_line_text);

            JPanel command_script_panel = new JPanel();
            command_script_text = new JTextArea(20,200);
            command_script_text.setEditable(false);
            command_script_text.setLineWrap(true);
            command_script_text.setWrapStyleWord(true);
            command_script_text.setFont(command_line_text.getFont());
            command_script_panel.add(command_script_text);

            addTab(TAB_CMD, scrolled(command_line_panel));
            addTab(TAB_SCRIPT, scrolled(command_script_panel));
        }

        @Override
        public Dimension getPreferredSize()
        {
            Dimension pref = super.getPreferredSize();
            return
                    new Dimension(Math.min(pref.width, Toolkit.getDefaultToolkit().getScreenSize().width-40),
                    Math.min(pref.height,Toolkit.getDefaultToolkit().getScreenSize().height/4));
        }

        private String getScriptVariableName(String arg)
        {
            return arg.replace('-', '_').toUpperCase();
        }

        private boolean addSwitch(String switch_name, String value, String explanation)
        {
            boolean want = (value != null && value.length()>0);
            if (want)
            {
                value = quotify(value);
                String var = getScriptVariableName(switch_name);
                preamble_sb.append(var);
                preamble_sb.append("=");
                preamble_sb.append(value);
                if (explanation != null)
                {
                    preamble_sb.append("\t# ");
                    preamble_sb.append(explanation);
                }
                preamble_sb.append("\n");

                script_sb.append(" -");
                script_sb.append(switch_name);
                script_sb.append(" ${");
                script_sb.append(var);
                script_sb.append("}");

                line_sb.append(" -");
                line_sb.append(switch_name);
                line_sb.append(" ");
                line_sb.append(value);
            }
            return want;
        }

        private String quotify(String value)
        {
            boolean quote_needed =
                    (value.indexOf(' ')>=0)
                    || (value.indexOf('\"')>=0)
                    || (value.indexOf('\\')>=0)
                    || (value.indexOf('$')>=0)
                    || (value.indexOf('|')>=0)
                    || (value.indexOf('*')>=0)
                    || (value.indexOf('?')>=0)
                    || (value.indexOf('!')>=0)
                    || (value.indexOf('\'')>=0)
                    || (value.indexOf('&')>=0);
            if (quote_needed)
            {
                StringBuffer out = new StringBuffer();
                out.append("\"");
                for (char c: value.toCharArray())
                {
                    if ("\"\\$|*?!'&".indexOf(c)>=0)
                    {
                        out.append("\\");
                    }
                    out.append(c);
                }
                out.append("\"");
                return out.toString();
            } else
                return value;

        }

        private void addArgument(String var, String value, String explanation)
        {
            preamble_sb.append(var);
            preamble_sb.append("=");

            preamble_sb.append(quotify(value));
            if (explanation != null)
            {
                preamble_sb.append("\t# ");
                preamble_sb.append(explanation);
            }
            preamble_sb.append("\n");

            script_sb.append(" ${");
            script_sb.append(var);
            script_sb.append("}");

            String[] args = value.split("\\s+");
            boolean after_switch=false;
            for (int i=0; i<args.length; i++)
                if (args[i].trim().length()!=0)
                {
                    line_sb.append(' ');
                    if (args[i].startsWith("-"))
                    {
                        line_sb.append(args[i]);
                        after_switch = true;
                    } else
                    {
                        if (after_switch)
                        {
                            line_sb.append("\"");
                            line_sb.append(args[i]);
                            line_sb.append("\"");

                            after_switch = false;
                        } else
                        {
                            line_sb.append(args[i]);
                        }
                    }
                }
        }
        
        
        private void showCommands()
        {
            line_sb = new StringBuilder();
            script_sb = new StringBuilder();
            preamble_sb = new StringBuilder();

            preamble_sb.append("#!/bin/bash\n");

            addArgument("JAVAEXEC", javaexec_TF.getValue().toString(), "path to Java executable");
            addSwitch("cp", jar_TF.getValue().toString(), "path to JAR bundle for "+Count.APP_TITLE);
        	addArgument("JAVAMAIN", mainclass_TF.getValue().toString(), "Java main class");
        	if (Count.UsesThreadpool.class.isAssignableFrom(command_main_class))
        		addSwitch(CommandLine.OPT_THREADS, javathreads_TF.getValue().toString(), "worker threads used in model optimization");
        	for (Object key: extra_options.keySet())
        	{
        		String opt = (String) key;
        		addSwitch(opt, extra_options.getProperty(opt), option_explanations.get(key));
        	}
        	if (load_bundle_dialog != null && load_bundle_dialog.hasFile())
        	{
        		addSwitch(CommandLine.OPT_LOAD, load_bundle_dialog.getFile(), "Input session");
        	}
        	if (output_dialog.hasFile())
        	{
        		addSwitch(CommandLine.OPT_OUTPUT, output_dialog.getFile(), "Output file");
        	}
        	if (output_tree_dialog != null && output_tree_dialog.hasFile())
        	{
        		addSwitch(CommandLine.OPT_OUTPUT_TREE, output_tree_dialog.getFile(), "Output tree");
        	}
        	if (save_bundle_dialog != null && save_bundle_dialog.hasFile())
        	{
        		addSwitch(CommandLine.OPT_SAVE, save_bundle_dialog.getFile(), "Output session");
        	}
            String tree_file = tree_dialog.hasFile()?tree_dialog.getFile():".";
        	String table_file = table_dialog.hasFile()?table_dialog.getFile():".";
        	addArgument("TREE", tree_file, "path to input phylogeny (\".\" means none)");
        	addArgument("TABLE", table_file, "path to input table (\".\" means none)");
        	if (rates_dialog != null)
        	{
        		String rates_file = rates_dialog.hasFile()?rates_dialog.getFile():".";
            	addArgument("MODEL", rates_file, "path to input model (\".\" means none)");
        	}
            command_line_text.setText(line_sb.toString());
            preamble_sb.append("\n");
            preamble_sb.append(script_sb);
            command_script_text.setText(preamble_sb.toString());
        }
//        @Override
//        public void actionPerformed(ActionEvent ae)
//        {
//            showCommands();
//        }
//
//        @Override
//        public void propertyChange(PropertyChangeEvent pce)
//        {
//            showCommands();
//        }

    }    

    
    /**
     * Displayed JDialog after runtime options are selected 
     */
    private class CommandDialog extends JDialog
    {
    	CommandDialog()
    	{
        	super(app, "How to run "+Count.APP_TITLE+" ("+command_main_class.getCanonicalName()+") in the command line or in a script", true); // modal
        	initComponents();
        	pack();
        	this.setLocationRelativeTo(null); // center on screen
    	}
    	
        private CommandPanel command_panel;
        
        private void initComponents()
        {
			command_panel = new CommandPanel();
			
	        JPanel button_box = new JPanel();
	        button_box.setLayout(new BoxLayout(button_box, BoxLayout.LINE_AXIS));
	        button_box.setBorder(BorderFactory.createLineBorder(Color.BLACK));
	        JButton saveB = new JButton("Save into file...");
	        saveB.setActionCommand("save");
	        saveB.setToolTipText("Saves the selected text (command line or script) into a text file");
	        JButton printB = new JButton("Print onto console");
	        printB.setActionCommand("print");
	        printB.setToolTipText("Prints the selected text (command line or script) to the standard output.");
	        JButton copyB = new JButton("Copy text to clipboard");
	        copyB.setActionCommand("copy");
	
	        JButton execB = new JButton("◀︎◀︎ GUI launch");
	        execB.setActionCommand("exec");
	        execB.setToolTipText("Returns to execute the optimization in the graphical interface.");
	
	        JButton closeB = new JButton("Close window");
	        closeB.setActionCommand("close");
	        closeB.setToolTipText("Closes this window");
	
	        ActionListener button_listener = new ActionListener()
	        {
	            @Override
	            public void actionPerformed(ActionEvent ae)
	            {
	                String cmd = ae.getActionCommand();
	                JTextComponent selected_comp = command_panel.getSelectedTextComponent();
	                String selected_text = selected_comp.getText();
	                if ("save".equals(cmd))
	                {
	                    String selected_kind = command_panel.getTitleAt(command_panel.getSelectedIndex());
	                    FileDialog FD = new FileDialog((java.awt.Frame)null, "Save "+selected_kind+" into a text file", FileDialog.SAVE);
	                    FD.setVisible(true);
	
	                    String file_name = FD.getFile();
	                    String directory = FD.getDirectory();
	                    if (file_name == null)
	                    {
	                    	closing_status = CANCEL_OPTION;
	                    } else {	                    
	                        try
	                        {
	                            File file = new File(directory, file_name);
	                            PrintStream out = new PrintStream(file);
	                            out.println(selected_text);
	                            out.close();
	                        } catch (java.io.FileNotFoundException e)
	                        {
	                            throw new RuntimeException(e);
	                        }
	                        closing_status = SAVE_OPTION;
	                    }
	                } else if ("print".equals(cmd))
	                {
	                    System.out.println(selected_text);
	                    closing_status = PRINT_OPTION;
	                } else if ("copy".equals(cmd))
	                {
	                    selected_comp.selectAll();
	                    selected_comp.copy();
	                    closing_status = COPY_OPTION;
	                } else if ("close".equals(cmd))
	                {
	                	closing_status = CANCEL_OPTION;
	                } else if ("exec".equals(cmd))
	                {
	                	closing_status = EXEC_OPTION;
	//                    setVisible(false);
	//                    gui.getFrame().dispose();
	//                    String exec_cmd = command_panel.command_line_text.getText();
	//                    try
	//                    {
	//                        Process exec_proc = Runtime.getRuntime().exec(exec_cmd);
	//                        exec_proc.waitFor();
	//                    } catch (IOException e)
	//                    {
	//                        //
	//                    } catch (InterruptedException e)
	//                    {
	//                        // 
	//                    } finally
	//                    {
	//                        gui.doQuit();
	//                    }
	                }
	                setVisible(false);
	            }
	            
	        };
	        saveB.addActionListener(button_listener);
	        printB.addActionListener(button_listener);
	        closeB.addActionListener(button_listener);
	        copyB.addActionListener(button_listener);
	        execB.addActionListener(button_listener);
	
	        button_box.add(execB);
	        button_box.add(saveB);
//	        button_box.add(printB); // printing onto console is useless
	        button_box.add(copyB);
	        button_box.add(closeB);
	        button_box.add(Box.createHorizontalGlue());
	
	        button_box.setBackground(COLOR_DECORATION);
	        
	        this.add(command_panel, BorderLayout.CENTER);
	        this.add(button_box, BorderLayout.SOUTH);
        }
        
    	
    }
}
