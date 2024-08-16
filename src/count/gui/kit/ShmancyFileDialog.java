package count.gui.kit;
/*
 * Copyright 2024 Mikl&oacute;s Cs&#369;r&ouml;s.
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


import java.awt.FileDialog;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.FilterInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.JTextField;
import javax.swing.SwingWorker;
import javax.swing.Timer;

import static java.awt.FileDialog.LOAD;
import static java.awt.FileDialog.SAVE;

import java.awt.BorderLayout;
import java.awt.Container;

import count.gui.CommandBuilderDialog;
import count.io.GeneralizedFileReader;


/**
 * A front end to the OS-specific {@link FileDialog}. 
 * User can enter a path name or URL; or else 
 * can make a file chooser pop up. 
 */
public class ShmancyFileDialog extends JDialog implements ActionListener
{
	static
	{
		JDialog.setDefaultLookAndFeelDecorated(true);
	}
	/**
	 * option value when ok button was pressed 
	 */
	public final static String OK_OPTION = "ok";
	/**
	 * option value when choose-file button was pressed (keeps this dialog open)
	 */
	public final static String CHOOSE_OPTION = "choose";
	/**
	 * option value when cancel button was pressed 
	 */
	public final static String CANCEL_OPTION = "cancel";
	/**
	 * option value before user pressed a button
	 */
	public final static String UNINITIALIZED_VALUE = "uninitialized";
			

	private static final int TIME_POPUP_PROGRESS = 25; // 25 for development // popup in 0.5 sec
	
	
	
	/*
	 * Instance variables 
	 */
	/**
	 * SAVE or LOAD	
	 */
	private final int file_mode;
	/**
	 * For the embeddable component: label of text field
	 */
	private final String entry_label;
	/**
	 * Local copy of selected file, in sync with textfield's content
	 */
	private String filename="";
	
	/**
	 * Embeddable component: enter a path name or choose 
	 */
	private JComponent front_panel;
	
	/**
	 * Displayed file name. 
	 */
	private JTextField fileT;
	
	/**
	 * Choose button for embeddable component with a label
	 */
	private JButton fileB;
	
	/**
	 * dialog from last choose action
	 */
	private FileDialog file_dialog=null;
			
	/**
	 * Out ok button for enabling/disabling 
	 */
	private JButton okButton;
	
	/**
	 * optipn choosing value: ok, cancel, uninitialized, or null
	 */
	private Object value;
	
	private Icon icon;
	
	
	/**
	 * Full instantiation
	 * 
	 * @param parent_frame
	 * @param title dialog title 
	 * @param mode LOAD or SAVE
	 * @param label null if no embeddable chooser needed, or text for file to be chosen
	 */
	public ShmancyFileDialog(Frame parent_frame, String title, int mode, String label, Icon icon)
	{	
		super(parent_frame, title, true); // init as a modal dialog
		if (mode != LOAD && mode != SAVE)
		{
			throw new IllegalArgumentException("File mode must be FileDialog.LOAD or FiuleDialog.SAVE");
		}
		this.file_mode = mode;
		this.entry_label = label;
		this.icon = icon;
		initComponents();
	}
	
	public ShmancyFileDialog(Frame parent_frame, String title, int mode, String label)
	{
		this(parent_frame, title, mode, label, null);
	}
	
	/**
	 * Instantiation with no embeddable chooser needed
	 * 
	 * @param frame
	 * @param title
	 * @param mode
	 */
	public ShmancyFileDialog(Frame frame, String title, int mode)
	{
		this(frame, title, mode, (String) null);
	}

	/**
	 * Instantiation with no embeddable chooser needed
	 * 
	 * @param frame
	 * @param title
	 * @param mode
	 */
	public ShmancyFileDialog(Frame frame, String title, int mode, Icon icon)
	{
		this(frame, title, mode, null, icon);
	}
	
	
	/**
	 * Load dialog for a file
	 * 
	 * @param frame
	 * @param title
	 * @param label
	 */
	public ShmancyFileDialog(Frame frame, String title, String label)
	{
		this(frame,title, LOAD, label, null);
	}
	
	public int getMode() { return file_mode;}
	
	/**
	 * Sets the local file name, as well as the super's version.
	 */
	public void setFile(String file)
	{
		if (!fileT.getText().equals(file==null?"":file))
			fileT.setText(file==null?"":file);
		this.filename = file;
		okButton.setEnabled(file!=null && !"".equals(file));
	}
	
	public void setDirectory(String dir)
	{
		if (file_dialog!=null) file_dialog.setDirectory(dir);
	}
	
	/**
	 * Locally stored file name.
	 */
	public String getFile()
	{
		if (filename==null || "".equals(filename)) // maybe Enter was not pressed
			setFile(fileT.getText());
		return this.filename;
	}
	
	public String getDirectory()
	{
		// if file chooser was used
		if (file_dialog == null) return null;
		else return file_dialog.getDirectory();
	}
	
	public boolean hasFile()
	{
		String getfile = getFile();
		return getfile != null && !"".equals(getfile); 
	}
	
	/**
	 * Concatenated directory and file.
	 * @return null if file is null
	 */
	public String getPath()
	{
		String f = getFile();
		String d = getDirectory();
		return f==null?null:(d==null?f:d+f);
	}

	/*
	 * Reading the selected file
	 */
	public Reader getReader() throws IOException
	{
		return GeneralizedFileReader.guessReaderForInput(getPath());		
	}
	
	public InputStream getInputStream() throws IOException
	{
		return GeneralizedFileReader.guessInputStreamForInput(getPath());
	}
	
	public BufferedReader getBufferedReader() throws IOException 
	{
		return GeneralizedFileReader.guessBufferedReaderForInput(getPath());
	}
	/**
	 * Option value : uninitialized before user set a file; ok or cancel 
	 * after selection is done; retrieve the file via {@link #getFile()}, 
	 * {@link #getPath()}, {@link #getDirectory()}
	 * 
	 * @return
	 */
	public Object getValue()
	{
		return value;
	}
	
	private void initComponents()
	{
		final JButton chooseButton = new JButton("Choose ...");
		chooseButton.setActionCommand(CHOOSE_OPTION);
		chooseButton.addActionListener(this);
		if (entry_label==null)
			fileB = null;
		else
			fileB =chooseButton;

		Box button_box = new Box(BoxLayout.LINE_AXIS);
		final JButton cancelButton = new JButton("Cancel");
		cancelButton.setActionCommand(CANCEL_OPTION);
		cancelButton.addActionListener(this);
		okButton = new JButton((getMode()==LOAD?"Load":"Save")+" above named file");
		okButton.setActionCommand(OK_OPTION);
		okButton.addActionListener(this);
		button_box.add(cancelButton);
		button_box.add(okButton);
		if (fileB == null)
		{
			button_box.add(chooseButton);
		}
		
		final JLabel fileL = new JLabel(entry_label==null?"  File name or URL":"  "+entry_label);
		
		fileT = CommandBuilderDialog.createTextField(CommandBuilderDialog.FieldType.FILE); //   new JTextField(FILENAME_FIELD_LENGTH);
		fileL.setLabelFor(fileT);
		String file_tooltip =  "Enter a URL (http:// or ftp://) or a local file path, and press OK.";
		if (getMode()==LOAD)
			file_tooltip = file_tooltip+"Compressed files are recognized by .gz or .zip extension, and uncompressed on the fly.";
		fileT.setToolTipText(file_tooltip);
		
		setFile(filename);
		fileT.addActionListener(a->setFile(fileT.getText()));
		
		Box front_box = new Box(BoxLayout.LINE_AXIS);
		front_box.add(fileL);
		front_box.add(fileT);
		if (fileB != null)
		{
			front_box.add(fileB);
		}
		front_panel = front_box;	

		
		
		// using the content pane's BorderLayout
		Container contentPane = getContentPane();
		contentPane.add(front_box, BorderLayout.CENTER);
	    contentPane.add(button_box, BorderLayout.PAGE_END);	
	    
	    if (icon==null)
	    {
		    if (getMode()==LOAD)
		    	icon = CountActions.createOpenIcon(CountActions.SIZE_XL);
		    else
		    	icon = CountActions.createSaveIcon(CountActions.SIZE_XL);
	    }
	    
	    contentPane.add(new JLabel(icon), BorderLayout.LINE_START);
	    
		getRootPane().setDefaultButton(chooseButton);
		
		this.value = UNINITIALIZED_VALUE;
	    
//		Object[] selectionValues = new Object[3];
//		int i=0;
//		select_ichoose = i++;
//		select_iok = i++;
//		select_icancel = i++;
//		
//		selectionValues[select_ichoose] = "Choose file ...";
//		selectionValues[select_iok] = (getMode()==LOAD?"Load":"Save")+" above named file";
//		selectionValues[select_icancel] = "Cancel";
//		
//		JOptionPane optionPane = new JOptionPane(front_box, JOptionPane.QUESTION_MESSAGE, JOptionPane.YES_NO_CANCEL_OPTION,  null, selectionValues, selectionValues[select_ichoose]);
//		this.setContentPane(optionPane);
		this.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
		this.addWindowListener(new WindowAdapter() {
		    public void windowClosing(WindowEvent we) {
		        //System.out.println("#**SFF.iC/windowClose");
		        
		    }
		});
		this.pack();
		this.setLocationRelativeTo(getOwner());
		//this.setLocation(200, 200);
	}
	
	/**
	 * A small panel with a label+textfield, plus choose-file button; instead 
	 * of the JDialog.
	 * 
	 * @return
	 */
	public JComponent getEmbeddableChooser()
	{
		return front_panel;
	}
	
	
	/**
	 * When user clicks on choose-file button
	 * 
	 */ 
	private void chooseFile()
	{
//		System.out.println("#**SFF.cF ");
		
		file_dialog = new FileDialog(this,getTitle(),file_mode);
		file_dialog.setFile(this.getFile());
		//file_dialog.setUndecorated(false);
		
		file_dialog.setVisible(true);
		if (file_dialog.getFile()!=null)
		{
			this.setFile(file_dialog.getFile());
			getRootPane().setDefaultButton(okButton);
			okButton.requestFocusInWindow();
			if (fileB==null)
				okButton.doClick();
		}
//		System.out.println("#**SFF.cF done");

	}

	/**
	 * Handling our button events
	 */
	@Override
	public void actionPerformed(ActionEvent click)
	{
		String cmd = click.getActionCommand();
		if (CANCEL_OPTION.equals(cmd))
		{
			setFile(null);
			this.setVisible(false);
		} else if (OK_OPTION.equals(cmd))
		{
			String file = fileT.getText();
			if ("".equals(file))
			{
				JOptionPane.showMessageDialog(getParent(), "Empty file path: nothing to load", "No file name", JOptionPane.WARNING_MESSAGE);
				setFile(null);
			} else
				setFile(file);
			// OK
			this.setVisible(false);
		} else if (CHOOSE_OPTION.equals(cmd))
		{
			// want a FileDialog
			chooseFile();
		} else
		{
			throw new UnsupportedOperationException("Command "+cmd+" is not recognized:" + click);
		}
		value=cmd;
	}
	
	
	/**
	 * Acts on the file text field and choose button
	 */
	@Override
	public void setEnabled(boolean e)
	{
		this.fileT.setEnabled(e);
		if (fileB != null)
			fileB.setEnabled(e);
		super.setEnabled(e);
	}
	
	/**
	 * Creates a Swing task for reading a given type, with a progress monitor popup.
	 * 
	 * @param <T>
	 * @param want_to
	 * @return
	 */
	public <T> SwingWorker<T,Void> createTask(LoadTask<T> want_to)
	{
		SwingWorker<T,Void> createTask = new SwingWorker<>()
		{
			MonitoredStream in = null;
			
			@Override
			protected T doInBackground() throws Exception 
			{
				try
				{
					InputStream	stream = getInputStream();
					if (stream==null) return null;				
					in = new MonitoredStream(stream);
				} catch (IOException E)
				{
					if (in==null || in.is_cancelled)
					{
						this.cancel(true);
						return null;
					} else
						throw E;
				}
//				System.out.println("#**FFD.cT.dIB starting read "+Thread.currentThread());
				return want_to.read(in);
			}
			
			@Override 
			protected void done()
			{
				if (in==null) return;
				try
				{
//					System.out.println("#**FFD.cT.dIB done "+Thread.currentThread());
					in.close();
				} catch (IOException unlikely)
				{
					throw new RuntimeException(unlikely); 
				}
			}
			
		};
		return createTask;
	}
	
	/**
	 * Loads (reads) an object of a given type from an input stream
	 * 
	 * @param <T> return value 
	 */
	public static interface LoadTask<T>
	{
		public T read(InputStream in) throws IOException;
	}
	
	/**
	 * Redefines {@link MonitoredStream#close}, {@link MonitoredStream#read()} 
	 * and {@link MonitoredStream#read()} to follow progress. 
	 * 
	 * @author csuros
	 *
	 */
	private class MonitoredStream extends FilterInputStream
	{
		/**
		 * Instantiate on worker thread.
		 * 
		 * @param stream
		 */
		MonitoredStream(InputStream stream)
		{
			super(stream);
			
			try
			{
				int max = stream.available();
				progress = new JProgressBar(0, max);
				progress.setIndeterminate(false);
			} catch (IOException cannot_set_max)
			{
				progress = new JProgressBar(0,0);
				//progress.setIndeterminate(true);
				
			}
			progress.setStringPainted(true);
			String what_are_we_doing = "Reading and parsing";
			Object[] pane_elements = {
					"Loading "+getPath()+"...", 
					progress
			};
			cancel = new JButton("Cancel");
			Object[] pane_buttons = 
			{
					cancel
			};
			JOptionPane pane = new JOptionPane(pane_elements, JOptionPane.INFORMATION_MESSAGE);
			pane.setOptions(pane_buttons);
			dialog = pane.createDialog(what_are_we_doing);
			dialog.setModal(true);
			dialog.setResizable(true);
			dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
			is_cancelled = false;
			cancel.addActionListener(click->is_cancelled = true);

			num_bytes_read = 0;
			showProgress();
			dialog.pack();
//			System.out.println("#**FFD.MS() bf show "+Thread.currentThread());
			popup_timer = new Timer(TIME_POPUP_PROGRESS, null); 
			popup_timer.setRepeats(false);
			popup_timer.addActionListener(tick->{
				dialog.setVisible(true);
			});
			popup_timer.start();
//			System.out.println("#**FFD.MS() done ");
		}
		JProgressBar progress;
		JDialog dialog;
		JButton cancel;
		Timer popup_timer;
		boolean is_cancelled;
		int num_bytes_read;
		private void showProgress()
		{
			progress.setValue(num_bytes_read);
			
			progress.setString(Integer.toString(num_bytes_read/1000)+" kb");
		}
		private void interruptIfCanceled() throws java.io.InterruptedIOException
		{
			if (is_cancelled) throw new java.io.InterruptedIOException();
		}
		@Override
		public void close() throws IOException
		{
			popup_timer.stop();
			dialog.setVisible(false);
			super.close();
		}
		@Override
		public int read() throws IOException
		{
			int read = super.read();
			if (read != -1)
			{
				num_bytes_read++;
				showProgress();
				interruptIfCanceled();
			}
			return read;
		}
		@Override
		public int read(byte[] buf) throws IOException
		{
			return this.read(buf, 0, buf.length);
		}
		@Override
		public int read(byte[] buf, int offset, int length) throws IOException
		{
			int read = super.read(buf, offset, length);
			if (read != -1)
			{
				num_bytes_read+=read;
				showProgress();
				interruptIfCanceled();
			}
			return read;
		}
	}

	
	
}
