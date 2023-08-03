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
import java.awt.geom.Ellipse2D;

/**
 * Disc-shaped icon. 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */

public class DiscIcon extends Ellipse2D.Double implements PointIcon
{
	
    private boolean filled    = false;
    private Color draw_color;
    private Color fill_color;
    
    public DiscIcon(int point, boolean filled)
    {
    	super(-RADIUS_FACTOR*point, -RADIUS_FACTOR*point, 
    			2.0*RADIUS_FACTOR*point, 2.0*RADIUS_FACTOR*point);
    	this.filled = filled;
    }

	@Override
	public int getIconWidth() 
	{
		return (int)getWidth();
	}

	@Override
	public int getIconHeight() 
	{
		return (int)getHeight();
	}

	private static double RADIUS_FACTOR = Math.sqrt(2.0/Math.PI);

	@Override
	public void setSize(int point_size) 
	{
		// 2*s^2 = r^2 * pi
		// 2*s^2/pi = r^2
		// r = s* sqrt(2.0/pi)
		double radius = RADIUS_FACTOR*point_size;
		this.setFrame(-radius, -radius, 2.0*radius, 2.0*radius);
	}

	@Override
	public void setDrawColor(Color c) 
	{
		this.draw_color = c;
		
	}

	@Override
	public void setFillColor(Color c) 
	{
		this.fill_color = c;
	}

	@Override
	public Color getFillColor() 
	{
		return fill_color;
	}

	@Override
	public void setFilled(boolean is_filled) 
	{
		this.filled = is_filled;
		
	}

	@Override
	public void paint(Graphics g, int x, int y, Color draw_color, Color fill_color) 
	{
        Color c=g.getColor();
        int diam = this.getIconWidth();
        
        if (fill_color != null)
        {
            g.setColor(fill_color);
            g.fillOval(x-diam/2, y-diam/2, diam, diam);
        }

        if (draw_color != null)
        {
            g.setColor(draw_color);
            g.drawOval(x-diam/2, y-diam/2, diam, diam);
        }
        g.setColor(c);
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
	
	@Override
	public void paintIcon(Component C, Graphics g, int x, int y) 
	{
        Color dc=(draw_color==null?C.getForeground():draw_color);
        Color fc=C.getBackground();
        if (filled)
            fc=(fill_color==null?dc:fill_color);
        int radius = this.getIconWidth()/2;
        paint (g, x+radius,y+radius,dc,fc);
	}

	

}
