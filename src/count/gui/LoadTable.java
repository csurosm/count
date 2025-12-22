package count.gui;
/*
 * Copyright 2025 Mikl&oacute;s Cs&#369;r&ouml;s.
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
import java.io.File;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.Arrays;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JPanel;
import javax.swing.SwingWorker;

import count.ds.AnnotatedTable;
import count.gui.kit.ShmancyFileDialog;
import count.io.CogAnnotator;
import count.io.DataFile;
import count.io.TableParser;

/**
 * GUI for loading an {@link AnnotatedTable}
 */
public class LoadTable extends JPanel {
	
	public LoadTable(AppFrame app)
	{
		super(new BorderLayout());
		this.app = app;
		initComponents();
	}
	
	private final AppFrame app;
	
	
	private ShmancyFileDialog tableD;
	private JCheckBox annotationsCB;
	private JCheckBox cogCB; 
	
	private void initComponents() {
		this.setAlignmentX(LEFT_ALIGNMENT);
		Box main = new Box(BoxLayout.PAGE_AXIS);
		main.setAlignmentX(LEFT_ALIGNMENT);
		
		annotationsCB = new JCheckBox("Keep annotations");
		annotationsCB.setToolTipText("If unchecked, then only columns with matching leaf names are kept; if checked then nonmatching columns are added as annotations");
		cogCB = new JCheckBox("Add COG(2024)/arCOG(2022)/asCOG(2020) annotations");
		cogCB.setToolTipText("If checked, three extra annotation columns are added autmatically based on cog/arcog identifiers.");
	
		Box buttonB = new Box(BoxLayout.LINE_AXIS);
		buttonB.add(annotationsCB);
		buttonB.add(cogCB);
		
		tableD = new ShmancyFileDialog(app, "Open table", "Table file");
		main.add(buttonB);
		main.add(tableD.getEmbeddableChooser());
		
		this.add(main, BorderLayout.CENTER);
	}
	
	
	public SwingWorker<DataFile<AnnotatedTable>, Void> readTableTask(){
		SwingWorker<DataFile<AnnotatedTable>, Void> read_task; 

		
		File table_file = new File(tableD.getDirectory(), tableD.getFile());
		boolean with_annotation_columns = annotationsCB.isSelected();
		boolean with_cogs = cogCB.isSelected();
    	Session sesh = app.getActiveSession();
    	String[] terminal_names = sesh.getModelBrowser().getMainPhylogeny().getLeafNames();
		
		read_task  = tableD.createTask(stream-> {
			Reader tableR = new InputStreamReader(stream);
			AnnotatedTable table 
        	= TableParser.readTable(terminal_names,tableR, with_annotation_columns);
			if (with_cogs) addCogAnnotations(table);
			DataFile<AnnotatedTable> table_data = new DataFile<>(table, table_file);
			
			return table_data;
		});		
		
		return read_task;
	}
	
	private void addCogAnnotations(AnnotatedTable table) {
		String[] cogprops = CogAnnotator.getAvailableProperties();
		for (String p: cogprops)
			if (0<=table.getPropertyIndex(p)) {
				//System.out.println("#**LT.aCA already have poperty "+p+"\t// "+Arrays.toString(table.getKnownProperties()));
				return; // we already have that
			}
		
		String[] fprops = table.getKnownProperties();
		
		for (String p: cogprops)
			table.registerProperty(p);
		int nF = table.getFamilyCount();
		
		// we look for cog identifiers in properties and family names
		for (int f=0; f<nF; f++) {
			String name = table.getFamilyName(f);
			String[] annot=CogAnnotator.getAnnotation(name);
			if (annot==null) {
				for (String p: fprops) {
					annot = CogAnnotator.getAnnotation(table.getFamilyProperty(f, p));
					if (annot != null)
						break;
				}
			}
			if (annot == null) annot = CogAnnotator.getNullAnnotation();
			for (int i=0; i<cogprops.length; i++)
				table.setFamilyProperty(f, cogprops[i], annot[i]);
		}
	}

}
