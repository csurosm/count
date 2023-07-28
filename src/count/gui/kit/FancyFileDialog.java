package count.gui.kit;

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

import java.awt.FileDialog;
import java.awt.Frame; // for FileDialog constructor

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.IOException;
import java.io.Reader;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JTextField;

import count.io.GeneralizedFileReader;

/**
 * A front end to the OS-specific {@link FileDialog}. 
 * User can enter a path name or URL; or else 
 * can make a file chooser pop up. 
 *
 */
public class FancyFileDialog extends FileDialog
{
	private static final int FILENAME_FIELD_LENGTH=40;

	public FancyFileDialog(Frame frame, String title, int mode)
	{
		this(frame, title, mode, null);
//		super(frame, title, mode);
//		this.title = title;
//		initComponents(null);
//		setFile(getFile()); // synchronize with super
	}
	
	public FancyFileDialog(Frame frame, String title, String label)
	{
		this(frame,title, LOAD, label);
	}
	
	private FancyFileDialog(Frame frame, String title, int mode, String entry_label)
	{
		super(frame, title, mode);
		this.title = title;
		initComponents(entry_label);
		setFile(getFile()); // synchronize with super
	}
	
	private String filename="";
	private String title;
	
	/**
	 * The first option pane's content: enter a path name.
	 */
	private JComponent front_panel;
	
	/**
	 * Displayed file name. 
	 */
	private JTextField fileT;
	
	private JButton fileB;
	
//	private JDialog front_dialog;
//	
	private void initComponents(String entry_label)
	{
		// JLabel[enter name] JTextField JButton[choose]
		final JLabel fileL = new JLabel(entry_label==null?"File name or URL":entry_label);
		fileT = new JTextField(FILENAME_FIELD_LENGTH);
		fileT.setText(filename);
		
		if (entry_label==null)
			fileB = null;
		else
			fileB = new JButton("Choose...");
		
		fileL.setLabelFor(fileT);
		String file_tooltip =  "Enter a URL or a local file path, and press OK.";
		if (getMode()==LOAD)
			file_tooltip = file_tooltip+"Compressed files are recognized by .gz or .zip extension, and uncompressed on the fly.";
		fileT.setToolTipText(file_tooltip);
		
		fileT.addActionListener(a->setFile(fileT.getText()));
		
		Box front_box = new Box(BoxLayout.LINE_AXIS);
		front_box.add(fileL);
		front_box.add(fileT);
		if (fileB != null)
		{
			fileB.addActionListener(click->chooseFile());
			front_box.add(fileB);
		}
		
		front_panel = front_box;
	}
	
	public JComponent getEmbeddableChooser()
	{
		return front_panel;
	}
	
	/**
	 * Sets the local file name, as well as the super's version.
	 */
	@Override
	public void setFile(String file)
	{
		super.setFile(file);
		if (!fileT.getText().equals(file==null?"":file))
			fileT.setText(file==null?"":file);
		this.filename = file;
//		System.out.println("#**FFD.setF "+file);
	}
	
	/**
	 * Locally stored file name.
	 */
	@Override
	public String getFile()
	{
		return this.filename;
	}
	
	/**
	 * Sets the locally stored title, and the super's title.
	 */
	@Override
	public void setTitle(String title)
	{
		super.setTitle(title);
		this.title = title;
	}
	
	/**
	 * Locally stored title.
	 */
	@Override
	public String getTitle()
	{
		return this.title;
	}
	
	/**
	 * Prior to popping up as a FileDialog, user can enter pathname/URL via {@link #front_panel}.
	 *  
	 */
	@Override
	public void setVisible(boolean visible)
	{
		if (visible)
		{
			// showInputDialog​(Component parentComponent, 
			// Object message, String title, int messageType, Icon icon, Object[] selectionValues, Object initialSelectionValue)			
			
			Object[] selectionValues = new Object[3];
			int i=0;
			int ichoose = i++;
			int iok = i++;
			int icancel = i++;
			
			selectionValues[ichoose] = "Choose file ...";
			selectionValues[iok] = (getMode()==LOAD?"Load":"Save")+" above named file";
			selectionValues[icancel] = "Cancel";
			
			int selected = JOptionPane.showOptionDialog(getParent(), front_panel, getTitle(), JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, selectionValues, selectionValues[ichoose]);

			if (selected==iok)
			{
				String file = fileT.getText();
				if ("".equals(file))
				{
					JOptionPane.showMessageDialog(getParent(), "Empty file path: nothing to load", "No file name", JOptionPane.WARNING_MESSAGE);
					setFile(null);
				} else
					setFile(file);
				// OK
			} else if (selected==ichoose)
			{
				// want a FileDialog
				chooseFile();
			} else // cancel or window closed 
			{
				setFile(null);
			}
		} else
		{
			super.setVisible(false);
		}
	}
	
	@Override
	public void setEnabled(boolean e)
	{
		this.fileT.setEnabled(e);
		if (fileB != null)
			fileB.setEnabled(e);
		super.setEnabled(e);
	}
	
	private void chooseFile()
	{
		super.setVisible(true);
		this.setFile(super.getFile());
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
	
	
}
