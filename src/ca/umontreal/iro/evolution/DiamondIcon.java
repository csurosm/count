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

public class DiamondIcon extends Polygon implements Icon {

  private int point_size=0;
  private boolean filled=false;

  private DiamondIcon(){}

  private int[] xbuf;
  private int[] ybuf;

  public DiamondIcon(int point_size, boolean filled) {
    this.point_size=point_size;
    this.filled=filled;

    this.addPoint(-point_size,0);
    this.addPoint(0,-point_size);
    this.addPoint(point_size,0);
    this.addPoint(0,point_size);

    xbuf = new int[4];
    ybuf = new int[4];
    for (int i=0; i<4; i++){ // values before translation
      xbuf[i]=this.xpoints[i];
      ybuf[i]=this.ypoints[i];
    }
  }

  public DiamondIcon(int point_size, Color c){
    this(point_size,false);
    paint_color=c;
  }

  Color paint_color=null;

  public int getIconHeight(){
    return 2*point_size;
  }

  public int getIconWidth(){
    return 2*point_size;
  }

  public Color fillColor=Color.lightGray;
  public Color emptyColor=Color.white;

  public void paint(Graphics g, int x, int y){
    this.translate(x,y);

    Color c=g.getColor();
    g.setColor(filled?fillColor:emptyColor);
    g.fillPolygon(this);

    if (paint_color != null){
      g.setColor(paint_color);
      g.drawPolygon(this);
      g.setColor(c);
    } else {
      g.setColor(c);
      g.drawPolygon(this);
    }
    for (int i=0; i<4; i++){
      xpoints[i]=xbuf[i];
      ypoints[i]=ybuf[i];
    }
  }

  public void paintIcon(Component unused_parameter, Graphics g, int x, int y){
    paint (g, x+point_size,y+point_size);
  }

  static String RCS_ID="$Id: DiamondIcon.java,v 1.1 2001/03/09 17:31:04 mcsuros Exp mcsuros $";
  //
  // $Log: DiamondIcon.java,v $
  // Revision 1.1  2001/03/09 17:31:04  mcsuros
  // Initial revision
  //
  //

}