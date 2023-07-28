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
import java.awt.Graphics;
import java.awt.Shape;

import javax.swing.Icon;

/**
 * Display style for a point.
 *
  *
 * A point (like in a scatterplot) can be displayed by using any shape. Must be able to
 * handle filled or unfilled points (diamonds, boxes, ...)
 * 
 * @since November 14, 2007, 12:21 AM
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 */
public interface PointIcon extends Shape, Icon
{
    /**
     * Sets the size of this point display.
     * 
     * @param point_size radius of the circle in which this guy has to fit.
     */
    public void setSize(int point_size);

    /**
     * Sets the preferred Color used for drawing the outline.
     * @param c if null, then default values are used.
     */
    public void setDrawColor(Color c);

    /**
     * Sets the preferred Color used for the inside of filled points.
     * @param c if null, then default values are used.
     */
    public void setFillColor(Color c);
    
    public Color getFillColor();

    /**
     * Sets the filling attribute. 
     * 
     * @param is_filled new value for whether this is filled 
     */
    public void setFilled(boolean is_filled);

    /**
     * Paints a point at the given location with the given colors.
     * (Does not change the drawing color of the graphics context.)
     * 
     * @param g graphics context
     * @param x X coordinate of the center.
     * @param y Y coordinate of the center.
     * @param draw_color Color used for the outline: if null, then no outline.
     * @param fill_color Color used for the inside: if null, then no filling.
     *
     */
    public void paint(Graphics g, int x, int y, Color draw_color, Color fill_color);

    /**
     * Paints a point at the given location with the given colors.
     * Drawing and filling colors are determined by the values set for this style
     * - if none, then defaults are taken from the graphics context.
     * 
     * @param g graphics context
     * @param x X coordinate of the center.
     * @param y Y coordinate of the center.
     */
    public void paint(Graphics g, int x, int y);
}