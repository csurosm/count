/*
 * AbstractPointSet.java
 *
 * Created on November 14, 2007, 12:05 AM
 */

package ca.umontreal.iro.evolution.malin;

import ca.umontreal.iro.evolution.malin.ui.*;

/**
 * Skeletal implementation of IndexedPointSet.
 *
 * @author  csuros
 */

public abstract class AbstractPointSet 
    extends java.util.AbstractSet
    implements IndexedPointSet 
{

  /**
   * Implementation based on add(IndexedPoint)
   */
  public boolean add(int point_index, double x, double y)
  {
    return add(new IndexedPoint(point_index,x,y));
  }

  /**
   * Implementation based on add(int, double, double);
   * @parm p mustn't be null.
   */
  public boolean add(IndexedPoint p)
  {
    return add(p.getIndex(),p.getX(),p.getY());
  }

  /**
   * @return false if not an IndexedPoint, otherwise normal behavior.
   */
  public boolean add(Object o)
  {
    if (o instanceof IndexedPoint)
      return add((IndexedPoint)o);
    else
      return false;
  }

  /**
   * Implementation based on closestPoint(int, int, double)
   */
  public IndexedPoint closestPoint(double x, double y)
  {
    return closestPoint(x,y,Double.POSITIVE_INFINITY);
  }
}