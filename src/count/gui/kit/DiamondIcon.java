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
import java.awt.Polygon;

/**
 * Diamond-shaped icon. 
 * 
 * @since November 14, 2007, 12:45 AM
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public final class DiamondIcon extends Polygon implements PointIcon
{
    private int     point_size= 0;
    private boolean filled    = false;
    private Color draw_color;
    private Color fill_color;
    
    private final int[] xbuf = new int[4];
    private final int[] ybuf = new int[4];

    
    
    public DiamondIcon(int point_size, boolean filled) 
    {
        this.filled=filled;
        setSize(point_size);
    }
    
    
    
    @Override
    public int getIconHeight()
    {
        return 2*point_size;
    }

    
    
    @Override
    public int getIconWidth()
    {
        return 2*point_size;
    }

    
    
    @Override
    public void paint(Graphics g, int x, int y, Color draw_color, Color fill_color)
    {
        translate(x,y);
        Color c=g.getColor();
        if (fill_color != null)
        {
            g.setColor(fill_color);
            g.fillPolygon(this);
        }

        if (draw_color != null)
        {
            g.setColor(draw_color);
            g.drawPolygon(this);
        }
        g.setColor(c);

        // reset coordinates
        for (int i=0; i<4; i++)
        {
          xpoints[i]=xbuf[i];
          ypoints[i]=ybuf[i];
        }
    }
    
    

    /**
     * Paints this DiamondIcon.
     * Drawing and filling colors are determined by the values set for this PointDisplay
     * - if none, then defaults are taken from C.
     *
     * @param C a non-null component
     * @param g graphics context
     * @param x X coordinate for upper left corner of icon.
     * @param y Y coordinate for upper left corner of icon.
     *
     */
    @Override
    public void paintIcon(Component C, Graphics g, int x, int y)
    {
        Color dc=(draw_color==null?C.getForeground():draw_color);
        Color fc=C.getBackground();
        if (filled)
            fc=(fill_color==null?dc:fill_color);
        paint (g, x+point_size,y+point_size,dc,fc);
    }
    
    

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
     * Sets a 45-degree rotated square with area 2*point_size*point_size.
     */
    @Override
    public void setSize(int point_size)
    {
        this.point_size=point_size;
        xbuf[0]=-point_size;  ybuf[0]=0;
        xbuf[1]=0;            ybuf[1]=-point_size;
        xbuf[2]=point_size;   ybuf[2]=0;
        xbuf[3]=0;            ybuf[3]=point_size;

        if (npoints==4)
        {
            // Shape already initialized
            for (int i=0; i<4; i++)
            {
                xpoints[i]=xbuf[i];
                ypoints[i]=ybuf[i];
            }
        }
        else
        {
          for (int i=0; i<4; i++)
            addPoint(xbuf[i], ybuf[i]);
        }
    }

    
    
    @Override
    public void setFilled(boolean is_filled){this.filled=is_filled;}
    
    
    
    @Override
    public void setDrawColor(Color c){draw_color=c;}
    
    
    
    @Override
    public void setFillColor(Color c){fill_color=c;}  
    
    @Override
    public Color getFillColor()   
    {
    	return fill_color;
    }
}

