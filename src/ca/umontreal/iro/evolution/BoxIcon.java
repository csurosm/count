package ca.umontreal.iro.evolution;


/**
 * Title:
 * Description:
 * Copyright:    Copyright (c)
 * Company:
 * @author
 * @version 1.0
 */

import javax.swing.Icon;
import java.awt.*;

public class BoxIcon extends Rectangle implements Icon {

  private int side = 0;
  private boolean filled;

  public BoxIcon(int point_size, boolean filled) {
    side = (int)(point_size*Math.sqrt(2)+.5);
    this.filled=filled;

    x=-side/2;
    y=-side/2;
    width=side;
    height=side;
  }

  public BoxIcon(int point_size, Color c){
    this(point_size, false);
    paint_color=c;
  }

  public int getIconHeight(){return side;}
  public int getIconWidth(){return side;}

  Color paint_color=null;
  public Color fillColor=Color.lightGray;
  public Color emptyColor=Color.white;

  public void  paint(Graphics g, int x, int y){
    Color c=g.getColor();
    g.setColor(filled?fillColor:emptyColor);
    g.fillRect(x-side/2,y-side/2,side,side);

    if (paint_color != null){
      g.setColor(paint_color);
      g.drawRect(x-side/2,y-side/2,side,side);
      g.setColor(c);
    } else {
      g.setColor(c);
      g.drawRect(x-side/2,y-side/2,side,side);
    }
  }


  public void paintIcon(Component unused_argument, Graphics g, int x, int y){
    paint(g,x+side/2,y+side/2);
  }

  static String RCS_ID="$Id: BoxIcon.java,v 1.1 2001/03/09 17:31:04 mcsuros Exp mcsuros $";
  //
  // $Log: BoxIcon.java,v $
  // Revision 1.1  2001/03/09 17:31:04  mcsuros
  // Initial revision
  //
  //


}