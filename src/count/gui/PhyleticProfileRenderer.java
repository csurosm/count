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

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;

import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;

import count.ds.AnnotatedTable;
import count.gui.kit.ColoredValueRenderer;

/**
 * A class for rendering phylogenetic profiles which are  
 * {@link AnnotatedTable.PhyleticProfile} instances.
 */
public class PhyleticProfileRenderer extends DefaultTableCellRenderer
{
    public PhyleticProfileRenderer(int num_leaves)
    {
        super();
        //this.setBackground(Color.BLACK);
        this.setOpaque(true);
        this.num_leaves = num_leaves;
    }
    /*
     * Updated when getTableCellRendererComponent is called
     */
    private AnnotatedTable.PhyleticProfile profile;
    private boolean is_selected;
    private final int num_leaves;
    private Color[] leaf_colors = null;
    
    public void setLeafColors(Color[] leaf_colors)
    {
    	this.leaf_colors = leaf_colors;
    }
    
    /**
     * Returns the default table cell renderer.
     * 
     * @param table the JTable
     * @param profile_array a phyletic profile 
     * @param isSelected true if cell is selected
     * @param hasFocus true if cell has focus
     * @param row the row of the cell to render
     * @param column the column of the cell to render
     * @return the default table cell renderer
     */
    @Override
    public Component getTableCellRendererComponent(
                        JTable table, Object profile_array,
                        boolean isSelected, boolean hasFocus,
                        int row, int column) 
    {
        //System.out.println("#*OTM.PR.gTCRC "+profile_array+"\ttype "+profile_array.getClass());
        profile = (AnnotatedTable.PhyleticProfile) profile_array;
        is_selected = isSelected;
        Component C = super.getTableCellRendererComponent(table, profile_array, isSelected, hasFocus, row, column);
        setText("");
        setToolTipText(profile.getPatternString());
        return C;
    }
    
    @Override
    public void paintComponent(Graphics g)
    {
        int[] pattern = profile.getProfile();
        
        assert (pattern.length == num_leaves);
        
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g.create();
        
        Color Cs = Color.BLACK;
        Color Cm = Color.BLACK;//
        Color Cmiss = Color.DARK_GRAY;
                
        if (profile != null)
        {
            int w = getWidth();
            int h = getHeight();
            //g2.setColor(LookAndFeel.MULTI_PRESENCE_COLOR);
            //g2.drawRect(0,0,w,h);
            int max_membership = 16;
            double scale_y = h/Math.log(1+max_membership);
            int n = pattern.length;
            double n_1 = 1.0/n;
            for (int i = 0; i<n; i++)
            {
                double pos_left  = w*i*n_1;
                double pos_right = pos_left+w*n_1;
                int x_left = (int)pos_left;
                if (x_left != pos_left) x_left++; // ceiling
                int x_right = (int)pos_right;
                if (x_right != pos_right)
                    x_right--;
                
                //System.out.println("#*OTM.PR.pC ["+i+"] "+pattern[i]+"\txl "+x_left+" xr "+x_right+"\t// "+rendering_profile.getPatternString());
                double l0 = (pattern[i]==0?0.0:Math.log(1+Math.min(Math.abs(pattern[i]),max_membership)));
                if (pattern[i]!=0)
                {
                    //if (x_left<=x_right) // put a little block with height proportional to logarithm
                    {
                    	if (pattern[i]<0)
                            g2.setColor(Cmiss);
                    	else if (leaf_colors != null && !is_selected)
                    		g2.setColor(leaf_colors[i]);
                    	else if (is_selected)
                    	{
                    		g2.setColor(Color.WHITE);
                    	} else
                    	{
	                        if (pattern[i]==1)
	                            g2.setColor(Cs);//LookAndFeel.SINGLE_PRESENCE_COLOR);
	                        else
	                            g2.setColor(Cm);//LookAndFeel.MULTI_PRESENCE_COLOR);
                    	} 
                    		
                        double pos_top = scale_y*l0;
                        int y_top = Integer.min(1+(int)pos_top,h);
                        g2.fillRect(x_left, h-y_top, x_right-x_left+1, y_top);
                    } 
//                    else
//                    {
//                        // nothing plotted: too thin
//                    }
                }
                // and plot the remainder: just a vertical line for smooth transitions at the borders
//                if (pos_right>(int)pos_right) // fill in the gap between two little blocks
//                {
//                    int x_border = (int)pos_right;
//                    double wt0 = pos_right-(int)pos_right;
//                    double wt1 = 1.-wt0;
//                    int y_top = (int)(scale_y*l0);
//                    Color C0 = this.getBackground();
//                    if (pattern[i]==1)
//                        C0 = Cs;//LookAndFeel.SINGLE_PRESENCE_COLOR;
//                    else if (pattern[i]>1)
//                        C0 = Cm;//LookAndFeel.MULTI_PRESENCE_COLOR;
//                    else if (pattern[i]<0)
//                        C0 = Cmiss;
//                    Color C1 = this.getBackground();
//                    if (i!=n-1)
//                    {
//                        // weighing with the next profile entry
//                        double l1 = (pattern[i+1]==0?0.0:Math.log(2.*Math.min(Math.abs(pattern[i+1]),max_membership)));
//                        y_top = (int)(scale_y*Math.max(l0,l1));// (int)(scale_y*(wt0*l0+wt1*l1));
//                        if (pattern[i+1]==1)
//                            C1 = Cs;//LookAndFeel.SINGLE_PRESENCE_COLOR;
//                        else if (pattern[i+1]>1)
//                            C1 = Cm;//LookAndFeel.MULTI_PRESENCE_COLOR;
//                        else
//                            C1 = Cmiss;
//                    }
//                    g2.setColor(ColoredValueRenderer.intermediateColor(C0, C1, wt0));
//                    g2.drawLine(x_border, h-y_top, x_border, h);
//                }
                //g2.setColor(Color.LIGHT_GRAY);
                //g2.drawLine((int)pos_left, h-1, (int)(pos_left), h-2); // a little tick
            }
            
        }
        
    }
    
    public int getPreferredRendererWidth()
    {
        return num_leaves*4;
    }
    
    
    @Override
    public Dimension getMinimumSize()
    {
        return new Dimension(getPreferredRendererWidth(),20);
    }
}
