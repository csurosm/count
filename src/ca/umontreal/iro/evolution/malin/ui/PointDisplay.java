/*
 * PointDisplay.java
 *
 * Created on November 14, 2007, 12:21 AM
 */

package ca.umontreal.iro.evolution.malin.ui;

/**
 * Display style for a point.
 *
 * A point can be displayed by using any shape. The display must be able to
 * handle filled or unfilled points.
 *
 * @author  csuros
 */


import java.awt.Color;
import java.awt.Graphics;
import java.awt.Shape;

public interface PointDisplay extends Shape 
{

  /**
   * Sets the size of this point display.
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

  public void setFilled(boolean is_filled);

  /**
   * Paints a point at the given location with the given colors.
   * @param x X coordinate of the center.
   * @param y Y coordinate of the center.
   * @param draw_color Color used for the outline: if null, then no outline.
   * @param fill_color Color used for the inside: if null, then no filling.
   *
   * (Does not change the drawing color of g.)
   */
  public void  paint(Graphics g, int x, int y, Color draw_color, Color fill_color);

  /**
   * Paints a point at the given location with the given colors.
   * Drawing and filling colors are determined by the values set for this style
   * - if none, then defaults are taken from g.
   * @param x X coordinate of the center.
   * @param y Y coordinate of the center.
   */
  public void paint(Graphics g, int x, int y);

}