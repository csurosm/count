/*
 * IndexedPoint.java
 *
 * Created on November 14, 2007, 12:02 AM
 */

package ca.umontreal.iro.evolution.malin;

/**
 * A Point with an integer index value.
 *
 * @author  csuros
 */


import java.awt.geom.Point2D;

public class IndexedPoint extends Point2D.Double 
{
    private int index;
    private static int next_id=1;

    /**
    * An IndexedPoint at (0,0).
    */
    public IndexedPoint(int index) 
    {
        init(index);
    }

    public IndexedPoint(int index, double x, double y)
    {
        super(x,y);
        init(index);
    }

    public IndexedPoint(int index, Point2D p)
    {
        super(p.getX(),p.getY());
        init(index);
    }

    private void init(int index)
    {
        this.index=index;
    }

    public void setIndex(int new_index){index=new_index;}
    public int getIndex(){return index;}

    protected String paramString()
    {
        StringBuffer sb=new StringBuffer();
        sb.append("idx ");
        sb.append(index);
        sb.append(" (");
        sb.append(x);
        sb.append(", ");
        sb.append(y);
        sb.append(")");
        return sb.toString();
    }

    public String toString()
    {
        StringBuffer sb=new StringBuffer("IPoint");//getClass().getName());
        sb.append('[');
        sb.append(paramString());
        sb.append(']');
        return sb.toString();
    }
}
