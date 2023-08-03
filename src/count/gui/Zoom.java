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
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;

import javax.swing.Box;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import count.gui.TreePanel.LayoutStyle;
import count.gui.kit.MagnificationSpinner;

/**
 * Embeds the tree panel into a scrollable interface.
 * There is a control bar at the bottom, 
 * which (by default) contains a magnification spinner.
 */
public class Zoom<PANEL extends TreePanel> extends JPanel
{
    private final PANEL tree_panel;
    private final JScrollPane tree_scroll;

    private final Box control_bar;
    
    private double initial_magnification = Double.NaN;
    
    private final MagnificationSpinner zoom_spinner;
    
    public Zoom(PANEL tree_panel)
    {
        this.tree_panel = tree_panel;
        this.zoom_spinner = tree_panel.createMagnificationSpinner();
        this.tree_scroll = new JScrollPane(tree_panel);
        this.control_bar = Box.createHorizontalBox();
//        tree_panel.addPropertyChangeListener(TreePanel.MAGNIFICATION_PROPERTY, chg->
//        		{
//        			int label_font_size = tree_panel.getLabelFontSize();
//        			if (tree_panel.getTreeLayoutStyle()==TreePanel.LayoutStyle.NODE_TABLE)
//        			{
//            			tree_scroll.getVerticalScrollBar().setUnitIncrement(4*label_font_size);
//        				tree_scroll.getHorizontalScrollBar().setUnitIncrement(label_font_size);
//        			} else
//        			{
//            			tree_scroll.getHorizontalScrollBar().setUnitIncrement(4*label_font_size);
//        				tree_scroll.getVerticalScrollBar().setUnitIncrement(label_font_size);
//        			}
//        		});
//        tree_scroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
//        tree_scroll.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        initComponents();
    }
    
    /**
     * Retrieves the control bar at the bottom
     * (so that one can remove remove or add components as necessary).
     * 
     * @return control bar in SOUTH; initially containing only a magnification spinner.
     */
    public Box getControlBar(){ return control_bar;}

    public MagnificationSpinner getSpinner() {return zoom_spinner;}
    
    public PANEL getTreePanel(){ return tree_panel;}
    
    private void initComponents()
    {
        tree_scroll.getViewport().setBackground(tree_panel.getBackground());

        this.setLayout(new BorderLayout());
        this.add(tree_scroll, BorderLayout.CENTER);
        
        control_bar.add(Box.createHorizontalGlue());
        control_bar.add(zoom_spinner);
        this.add(control_bar, BorderLayout.SOUTH);
        
    }
    
    public void setInitialMagnification(Graphics g)
    {
    	tree_panel.calculateBoundingBoxes(g);

    	Dimension preferred_size = tree_panel.getPreferredSize();
        double mag=1.0;
        if (LayoutStyle.NODE_TABLE == tree_panel.getTreeLayoutStyle())
        {
        	mag = 1.0;
        } else
        {
//        	int w = tree_panel.getWidth();
//            double dw = w/preferred_size.getWidth();
//            mag = Math.sqrt(dw);
//            int h = tree_panel.getHeight();
//            double dh = h/preferred_size.getHeight();
//            
//            System.out.println("#**Z.sIM dw "+dw+" dh "+dh+"\tpref "+getPreferredSize()
//            	+" w "+w+" h " +h);
//            
//	        mag = Double.min(mag, Math.sqrt(dh));
        }
//        
//        
//        
//        int w = tree_panel.getWidth();
//        int h = tree_panel.getHeight();
//        double preferred_area = preferred_size.getWidth()*preferred_size.getHeight();
//        double actual_area = w*h;
//        double mag = actual_area/preferred_area;
//        
//        double dw = w/preferred_size.getWidth();
//        double dh = h/preferred_size.getHeight();
//        mag = Double.min(mag, Math.sqrt(dw));
//        mag = Double.min(mag, Math.sqrt(dh));
        mag = Double.max(mag, MagnificationSpinner.MAGNIFICATION_MIN);
        mag = Double.min(mag, MagnificationSpinner.MAGNIFICATION_MAX);
        
        this.initial_magnification = mag;
        tree_panel.setMagnification(mag);
    	
    }
    
    @Override
    protected void paintComponent(Graphics g)
    {
        if (Double.isNaN(this.initial_magnification))
        {
        	setInitialMagnification(g);
        }
        super.paintComponent(g);
    }
    
    /**
     * Adds a control widget to the control bar. 
     * 
     * @param idx
     * @param C
     */
    public void addControlAt(int idx, Component C)
    {
    	Component[] controls = control_bar.getComponents();
    	control_bar.removeAll();
    	int i=0;
    	while (i<idx && i<controls.length)
    		control_bar.add(controls[i++]);
    	control_bar.add(C);
    	while (i<controls.length)
    		control_bar.add(controls[i++]);
    	control_bar.validate();
    	control_bar.repaint();
    }
    
    public Component replaceControlAt(int idx, Component C)
    {
    	Component[] controls = control_bar.getComponents();
    	control_bar.removeAll();
    	int i=0;
    	while (i<idx && i<controls.length)
    		control_bar.add(controls[i++]);
    	Component replace = controls[i++];
    	control_bar.add(C);
    	while (i<controls.length)
    		control_bar.add(controls[i++]);
    	control_bar.validate();
    	control_bar.repaint();
    	return replace;
    }
    
    public Component removeControlAt(int idx)
    {
    	Component[] controls = control_bar.getComponents();
    	control_bar.removeAll();
    	int i=0;
    	while (i<idx && i<controls.length)
    		control_bar.add(controls[i++]);
    	Component remove = controls[i++];
    	while (i<controls.length)
    		control_bar.add(controls[i++]);
    	control_bar.validate();
    	control_bar.repaint();
    	return remove;
    }
    
    @Override 
    public String toString()
    {
    	return tree_panel.getTreeName();
    }
}
