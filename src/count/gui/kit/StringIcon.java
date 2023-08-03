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
import java.awt.Font;
import java.awt.Graphics;
import javax.swing.Icon;

/**
 * An icon that shows a String. 
 * 
 */
public class StringIcon implements Icon
{
    public StringIcon(int width, int height, String text_shown)
    {
        this.width = width;
        this.height = height;
        this.icon_string = text_shown;
    }

    private final int width;
    private final int height;
    private final String icon_string; 
    
    private boolean in_center = true;
    
    public void setIndexing(boolean as_subscript)
    {
    	this.in_center = !as_subscript;
    }


    @Override
    public int getIconWidth(){return width;}
    @Override
    public int getIconHeight(){return height;}
    @Override
    public void paintIcon(Component c, Graphics g, int x, int y) // x y is top left 
    {
        Graphics myg = g.create();
        if (in_center)
        {
            myg.setColor(Color.BLACK);
        	myg.setFont(myg.getFont().deriveFont(1.2f*height).deriveFont(Font.BOLD));
            int w = myg.getFontMetrics().stringWidth(icon_string);
            myg.drawString(icon_string, x+width/2-w/2, y+height);
        } else
        {
        	myg.setColor(Color.WHITE);
        	myg.fillArc(x+width/2, y+height/2, width, height, 90, 90);
            myg.setColor(Color.BLACK);
        	myg.setFont(myg.getFont().deriveFont(0.6f*height).deriveFont(Font.BOLD));
        	DrawString.drawRight(myg, icon_string, x+width, y+height);
        }
    }
    
    private static final int DEF_WIDTH=30;
    private static final int DEF_HEIGHT = 20;
    
    public static StringIcon createRightPointingFinger()
    { 
        return new StringIcon(DEF_WIDTH,DEF_HEIGHT,"\u261e");
    }
    
    public static StringIcon superposedStringIcon(final Icon under, String letter)
    {
    	StringIcon superposed = 
    			new StringIcon(under.getIconWidth(), under.getIconHeight(), letter)
    			{
    				@Override
    			    public void paintIcon(Component c, Graphics g, int x, int y)  
    			    {
    					under.paintIcon(c, g, x, y);
    					super.paintIcon(c, g, x, y);
    			    }
    			};
    	return superposed;
    }
    
    public static StringIcon indexingStringIcon(final Icon under, String letter)
    {
    	StringIcon indexing = superposedStringIcon(under, letter);
    	indexing.setIndexing(true);
    	return indexing;
    }

}