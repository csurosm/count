/*
 * BoxIcon.java
 *
 * Created on November 14, 2007, 12:40 AM
 */

package ca.umontreal.iro.evolution.malin.ui;

/**
 *
 * @author  csuros
 */

import javax.swing.Icon;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Rectangle;

public class BoxIcon extends Rectangle implements Icon, PointDisplay
{

  private int side = 0;
  private boolean filled;

  public BoxIcon(int point_size, boolean filled) {
    setFilled(filled);
    setSize(point_size);
  }

  /**
   * Unfilled BoxIcon with the given drawing color.
   */
  public BoxIcon(int point_size, Color draw_color){
    this(point_size, false);
    drawColor=draw_color;
  }

  public int getIconHeight(){return side;}
  public int getIconWidth(){return side;}

  public void  paint(Graphics g, int x, int y, Color draw_color, Color fill_color){
    Color c=g.getColor();

    if (fill_color != null){
      g.setColor(fill_color);
      g.fillRect(x-side/2,y-side/2,side,side);
    }

    if (draw_color != null){
      g.setColor(draw_color);
      g.drawRect(x-side/2,y-side/2,side,side);
    }
    g.setColor(c);
  }

  /**
   * If unfilled, then only outline is drawn. If filled, then default
   * filling color is a bit darker than the drawing color.
   */
  public void paint(Graphics g, int x, int y){
    Color dc=(drawColor==null?g.getColor():drawColor);
    Color fc=null;
    if (filled)
      fc=(fillColor==null?dc.darker():fillColor);
    paint(g,x,y,dc,fc);
  }

  /**
   * Paints this BoxIcon.
   * Drawing and filling colors are determined by the color values set
   * - if none, then defaults are taken from C.
   * Default filling color for unfilled points is the background color of C.
   *
   * @param x X coordinate for upper left corner of icon.
   * @param y Y coordinate for upper left corner of icon.
   *
   */
  public void paintIcon(Component C, Graphics g, int x, int y){
    Color dc=(drawColor==null?C.getForeground():drawColor);
    Color fc=C.getBackground();
    if (filled)
      fc=(fillColor==null?dc:fillColor);
    paint(g,x+side/2,y+side/2,dc,fc);
  }

  private Color fillColor;
  private Color drawColor;

  public void setFillColor(Color c){
    fillColor=c;
  }

  public void setDrawColor(Color c){
    drawColor=c;
  }

  public void setFilled(boolean is_filled){
    filled=is_filled;
  }

  private static final double SQRT2=Math.sqrt(2.);

  public void setSize(int point_size){
    side = (int)(point_size*SQRT2+.5);

    x=-side/2;
    y=-side/2;
    width=side;
    height=side;
  }


}