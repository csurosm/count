/*
 * IndexedPointSet.java
 *
 * Created on November 14, 2007, 12:03 AM
 */

package ca.umontreal.iro.evolution.malin;

import java.util.Set;


public interface IndexedPointSet extends Set 
{

  /**
   * Adds a point to the set with the given index.
   *
   * @return false iff there is an IndexedPoint in the set already with the same index.
   */
  public boolean add(int point_index, double x, double y);

  /**
   * @return false iff there is an IndexedPoint in the set already with the same index.
   */
  public boolean add(IndexedPoint P);


  /**
   * Finds a point by index.
   * @return null if there is no point with that index.
   */
  public IndexedPoint get(int i);

  /**
   * Nearest-neighbor query.
   * @return the closest point to (x,y) or null if the set is empty.
   */
  public IndexedPoint closestPoint(double x, double y);

  /**
   * Nearest-neighbor query.
   * @return the closest point to (x,y) or null if there is no point within
   * max_distance.
   */
  public IndexedPoint closestPoint(double x, double y, double max_distance);

  /**
   * @return List of IndexedPoints falling into the rectangle.
   */
  public java.util.List<IndexedPoint> withinRectangle(java.awt.Rectangle r);

}