package count.gui.kit;
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
import java.awt.Graphics;
import java.awt.Rectangle;
/**
 * Box-style icon for points in a scatterplot and such. 
 * 
 * @since November 14, 2007, 12:40 AM
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public final class BoxIcon extends Rectangle implements PointIcon
{
    private static final double SQRT2=Math.sqrt(2.);

    private int     side       = 0;
    private boolean filled     = false;
    private Color   fill_color = null;
    private Color   draw_color = null;
    
    private Color   cross_color = null;

    public BoxIcon(int point_size, boolean filled) 
    {
        setFilled(filled);
        setSize(point_size);
    }

    @Override
    public int getIconHeight(){return side;}
    @Override
    public int getIconWidth(){return side;}

    @Override
    public void  paint(Graphics g, int x, int y, Color draw_color, Color fill_color)
    {
        Color c=g.getColor();

        if (fill_color != null)
        {
            g.setColor(fill_color);
            g.fillRect(x-side/2,y-side/2,side,side);
        }

        if (draw_color != null)
        {
            g.setColor(draw_color);
            g.drawRect(x-side/2,y-side/2,side,side);
        }
        
        if (cross_color != null)
        {
            g.setColor(cross_color);
            g.drawLine(x-side/2, y+side/2, x+side/2, y-side/2);
            g.drawLine(x-side/2, y-side/2, x+side/2, y+side/2);
        }
        
        g.setColor(c);
    }

    /**
     * If unfilled, then only outline is drawn. If filled, then default
     * filling color is a bit darker than the drawing color.
     */
    @Override
    public void paint(Graphics g, int x, int y)
    {
        Color dc=(draw_color==null?g.getColor():draw_color);
        Color fc=null;
        if (filled)
            fc=(fill_color==null?dc.darker():fill_color);
        paint(g,x,y,dc,fc);
    }

    /**
     * Paints this BoxIcon.
     * Drawing and filling colors are determined by the color values set
     * - if none, then defaults are taken from C.
     * Default filling color for unfilled points is the background color of C.
     *
     * @param C a component, not null
     * @param g graphics context
     * @param x X coordinate for upper left corner of icon.
     * @param y Y coordinate for upper left corner of icon.
     */
    @Override
    public void paintIcon(Component C, Graphics g, int x, int y)
    {
        Color dc=(draw_color==null?C.getForeground():draw_color);
        Color fc=C.getBackground();
        if (filled)
            fc=(fill_color==null?dc:fill_color);
        paint(g,x+side/2,y+side/2,dc,fc);
    }

    @Override
    public void setFillColor(Color c)
    {
        fill_color=c;
    }
    
    @Override
    public Color getFillColor()   
    {
    	return fill_color;
    }

    @Override
    public void setDrawColor(Color c)
    {
        draw_color=c;
    }
    
    public void setCrossing(Color c)
    {
        this.cross_color = c;
    }

    @Override
    public void setFilled(boolean is_filled)
    {
        filled=is_filled;
    }

    /**
     * Sets the icon size for area 2*point_size*point_size
     * (matching the area of a diamond with same point size). 
     * 
     * @param point_size half the length of the diagonal. 
     */
    @Override
    public void setSize(int point_size)
    {
        side = (int)(point_size*SQRT2+.5);

        // center Rectangle at (0,0)
        x=-side/2; 
        y=-side/2; 
        width=side;
        height=side;
    }   

}
