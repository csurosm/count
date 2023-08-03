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


import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.Graphics;
import java.util.Random;

import javax.swing.Icon;

import count.gui.RatesTreePanel;

/**
 * {@link Icon} for data table or ancestral reconstruction.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */
public class TableIcon implements Icon
{
	private static final int SIZE = 128; // should be a multiple of 8
//	private static final int CELL_WIDTH = SIZE/4;
//	private static final int CELL_HEIGHT = CELL_WIDTH/2;
	private static final Color BORDER_COLOR = Color.GRAY;
	private static final Color SELECTION_COLOR = Color.BLUE;
	private static final Color[] BACKGROUND_COLORS = null;

	/**
	 * 
	 * @param size
	 * @param cell_width
	 * @param cell_height
	 * @param background_colors
	 * @param border_color
	 * @param selection_color null if no selections illustrated
	 */
	public TableIcon(int size, int cell_width, int cell_height, Color[] background_colors, Color border_color, Color selection_color)
	{
		this.size = size;
		this.dw = cell_width;
		this.dh = cell_height;
		this.background_colors = background_colors;
		this.border_color = border_color;
		this.selection_color = selection_color;
		if (selection_color == null)
		{
			is_row_selected = null;
		} else
		{
			Random RND = new Random();
			int nrows = getIconHeight()/dh;
			is_row_selected = new boolean[nrows];
			int num_selected = 0;
			for (int i=0; i<nrows; i++)
			{
				boolean select = (num_selected<i-1) && ((num_selected< i/4) || (RND.nextDouble()<0.333));
				is_row_selected[i] = select;
				if (select)
					++num_selected;
			}
		}
	}
	
	public TableIcon(int icon_size, int cell_width, int cell_height, boolean want_selection)
	{
		this(icon_size, cell_width, cell_height, BACKGROUND_COLORS, BORDER_COLOR, want_selection?SELECTION_COLOR:null);
	}
	
	public TableIcon(int icon_size, boolean want_selection)
	{
		this(icon_size, preferredWidth(icon_size), preferredHeight(icon_size), want_selection);
	}
	
	public TableIcon(boolean want_selection)
	{
		this(SIZE, want_selection);
	}
	
	private static final int SIZE_M=32;
	
	private static final int preferredWidth(int icon_size)
	{
		return icon_size>SIZE_M?icon_size/8:icon_size/6;
	}
	
	private static final int preferredHeight(int icon_size)
	{
		return icon_size>SIZE_M?icon_size/16:icon_size/8;
	}
	
	public void setBinary(boolean isbinary)
	{
		this.isbinary = isbinary;
	}
	
	private final int size;
	private final int dw;
	private final int dh;
	private final Color border_color;
	private final Color selection_color;
	private final Color[] background_colors;
	private boolean isbinary=false;
	private final boolean[] is_row_selected;
	
	private boolean wantSelection() { return selection_color != null;}
	
    @Override
    public int getIconHeight()
    {
        return size;
    }

    @Override
    public int getIconWidth()
    {
        return size;
    }
    
    @Override
    public void paintIcon(Component C, Graphics g, int x, int y)
    {
        int w = getIconWidth();
        int h = getIconHeight();
        Graphics gg = g.create();
        
        gg.translate(x,y);
        if (background_colors == null)
        {
	        gg.setColor(Color.WHITE);
	        gg.fillRect(0,0,w,h);
        }
        gg.setColor(border_color);
        gg.drawRect(0,0,w,h);
        int num_selected = 0;
        for (int i=0; (i+1)*dh<=h; i++)
        {
            int row_y = dh*i;
            boolean is_selected = wantSelection() && is_row_selected[i];
//            		&& (num_selected<i-1) && ((num_selected< i/6) || (RND.nextDouble()<0.333));
            if (is_selected)
            {
                num_selected++;
                gg.setColor(selection_color);
                gg.fillRect(0, row_y, w, dh);
                
            }
            if (background_colors != null)
            {
                for (int j=0, cell_x=0; cell_x<w; j++, cell_x += dw)
                {
                	Color cell_color = ColoredValueRenderer.intermediateColor(Color.WHITE, background_colors[j%background_colors.length], 0.5);
                	gg.setColor(cell_color);
                	gg.fillRect(cell_x, row_y, dw, dh);
                }
            }
            gg.setColor(border_color);
            if (i!=0) 
                gg.drawLine(0, row_y, w, row_y);
            for (int j=1; j*dw<w; j++)
            {
            	int cell_x = j*dw;
        		gg.drawLine(cell_x, row_y, cell_x, row_y+dh);
            }
        }
        if (isbinary)
        {
        	gg.setFont(new Font("Serif",Font.PLAIN,3*size/4));
        	gg.setColor(Color.BLACK);
        	DrawString.drawCentered(gg, "01", size/2, size/2+gg.getFont().getSize()/2);
        }
        
    }	
    
    public static TableIcon historyIcon()
    {
    	return historyIcon(SIZE_M);
    }
    
    public static TableIcon historyIcon(int size)
    {
    	int num_shades = 2; //size>SIZE_M?3:2;
    	Color[] background = new Color[3*num_shades];
    	int ci = 0;
    	for (int i=0; i<num_shades; i++, ci++)
    	{
    		background[ci] = ColoredValueRenderer.intermediateColor(Color.WHITE, RatesTreePanel.RATES_GAIN_COLOR, (i+1.0)/(num_shades+1.0));
    	}
    	for (int i=0; i<num_shades; i++, ci++)
    	{
    		background[ci] = ColoredValueRenderer.intermediateColor(Color.WHITE, RatesTreePanel.RATES_LOSS_COLOR, (i+1.0)/(num_shades+1.0));
    	}
    	for (int i=0; i<num_shades; i++, ci++)
    	{
    		background[ci] = ColoredValueRenderer.intermediateColor(Color.WHITE, RatesTreePanel.RATES_DUPLICATION_COLOR,(i+1.0)/(num_shades+1.0));
    	}
    	TableIcon historyIcon = new TableIcon(size, preferredWidth(size), preferredHeight(size), background, BORDER_COLOR, null);
    	return historyIcon;
    }
    
    public static TableIcon columnFilterIcon(int size)
    {
//    	int ncol = (int)Math.ceil(size/((double)preferredWidth(size)));
    	Color[] background = new Color[4];
    	int c=0; 
    	while (c<3)
    	{
    		background[c++] = Color.WHITE;
    	}
    	while (c<background.length)
    		background[c++] = Color.DARK_GRAY;
    	TableIcon columnFilterIcon = new TableIcon(size, preferredWidth(size), preferredHeight(size), background, BORDER_COLOR, null);
    	return columnFilterIcon;
    }	
}
