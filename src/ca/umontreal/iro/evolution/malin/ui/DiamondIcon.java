/*
 * DiamondIcon.java
 *
 * Created on November 14, 2007, 12:45 AM
 */

package ca.umontreal.iro.evolution.malin.ui;

/**
 * Diamond shaped point display.
 *
 * @author  csuros
 */


import javax.swing.Icon;
import java.awt.Polygon;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Component;

public class DiamondIcon extends Polygon implements Icon, PointDisplay 
{

  private int point_size=0;
  private boolean filled=false;

  private DiamondIcon(){}

  private int[] xbuf;
  private int[] ybuf;

  public DiamondIcon(int point_size, boolean filled) 
  {
    this.filled=filled;

    xbuf = new int[4];
    ybuf = new int[4];

    setSize(point_size);
  }

  public DiamondIcon(int point_size, Color c)
  {
    this(point_size,false);
  }

  public int getIconHeight()
  {
    return 2*point_size;
  }

  public int getIconWidth()
  {
    return 2*point_size;
  }

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
   * @param x X coordinate for upper left corner of icon.
   * @param y Y coordinate for upper left corner of icon.
   *
   */
   public void paintIcon(Component C, Graphics g, int x, int y)
   {
    Color dc=(drawColor==null?C.getForeground():drawColor);
    Color fc=C.getBackground();
    if (filled)
      fc=(fillColor==null?dc:fillColor);
    paint (g, x+point_size,y+point_size,dc,fc);
  }

  public void paint(Graphics g, int x, int y)
  {
    Color dc=(drawColor==null?g.getColor():drawColor);
    Color fc=null;
    if (filled)
      fc=(fillColor==null?dc.darker():fillColor);
    paint(g,x,y,dc,fc);
  }

  public void setSize(int point_size)
  {
    this.point_size=point_size;

    xbuf[0]=-point_size;  ybuf[0]=0;
    xbuf[1]=0;            ybuf[1]=-point_size;
    xbuf[2]=point_size;   ybuf[2]=0;
    xbuf[3]=0;            ybuf[3]=point_size;

    if (npoints==4)
      // already initialized
      for (int i=0; i<4; i++)
      {
        xpoints[i]=xbuf[i];
        ypoints[i]=ybuf[i];
      }
    else
      for (int i=0; i<4; i++)
        addPoint(xbuf[i],ybuf[i]);
  }

  public void setFilled(boolean is_filled){this.filled=is_filled;}

  Color drawColor;
  Color fillColor;

  public void setDrawColor(Color c){drawColor=c;}
  public void setFillColor(Color c){fillColor=c;}

}