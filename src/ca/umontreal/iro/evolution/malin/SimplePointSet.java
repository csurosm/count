/*
 * SimplePointSet.java
 *
 * Created on November 14, 2007, 12:06 AM
 */

package ca.umontreal.iro.evolution.malin;

import java.util.Vector;
import java.util.List;
import java.util.Iterator;
import java.util.ListIterator;
import java.awt.Rectangle;

public class SimplePointSet extends AbstractPointSet  
{

  // point w/ index i is stored at data[i]
  // some entries in data may be null
  private IndexedPoint[] data;
  private int number_of_elements; //every non-null entry in data[] is before this one
  private int capacity;

  private int number_of_points; // number of non-null entries in data[]

  private int number_of_iterators=0;
  private boolean modified_while_iterating;

  public SimplePointSet()
  {
    this(10);
  }

  public SimplePointSet(int initial_capacity) 
  {
    this.data=new IndexedPoint[initial_capacity];
    this.capacity=initial_capacity;
    this.number_of_elements=0;
    this.number_of_points=0;
  }

  /**
   * @param p a point with <em>non-negative</em> index for a new point.
   * If p is null or the index is
   * negative, then p is not added.
   *
   * @return true iff the point was added to the set.
   *
   * <strong>The space occupied by this class is determined by
   * the maximum index in the set.</strong>
   *
   */
  public boolean add(IndexedPoint p)
  {
    int idx;
    if (p==null || (idx=p.getIndex())<0) return false;
    ensureNewElementCanBeAdded(idx);
    if (data[idx] != null)
      return false;
    checkConcurrency();
    data[idx]=p;
    number_of_points++;
    return true;
  }

  private void checkConcurrency()
  {
    modified_while_iterating=(number_of_iterators>0);
  }

  public IndexedPoint get(int idx)
  {
    return data[idx];
  }

  /**
   * Removes all points from the set. (Does not release memory.)
   */
  public void clear()
  {
    if (number_of_points>0)
    {
      checkConcurrency();
      java.util.Arrays.fill(data,null);
      number_of_points=0;
    }
    number_of_elements=0;
  }

  public int size(){return number_of_points;}

  /**
   * Finds the closest point.
   *
   * Note: the search time of this method is proportional to
   * the maximum index in the set.
   */
  public IndexedPoint closestPoint(double x, double y, double within)
  {
    double d_min=within;
    int i_min=-1;
    for (int i=0; i<number_of_elements; i++)
      if (data[i] != null){
        double d=data[i].distanceSq(x,y);
        if (d<d_min){
          d=d_min;
          i_min=i;
        }
      }
    if (i_min>=0)
      return data[i_min];
    else
      return null;
  }

  public java.util.List<IndexedPoint> withinRectangle(Rectangle r)
  {
    Vector<IndexedPoint> v=new Vector<IndexedPoint>();
    for (int i=0; i<number_of_elements; i++)
      if (data[i] != null && r.contains(data[i]))
        v.add(data[i]);
    return v;
  }


  public Iterator iterator()
  {
    return new SimplePointSetIterator();
  }

  /**
   * Removes a point.
   * @return true iff p was in the set.
   */
  public boolean remove(IndexedPoint p)
  {
    if (!contains(p))
      return false;
    checkConcurrency();
    data[p.getIndex()]=null;
    number_of_points--;
    return true;
  }

  /**
   * Removes a point.
   * @return true iff o is an IndexedPoint in the set.
   */
  public boolean remove(Object o)
  {
    if (o instanceof IndexedPoint)
      return remove((IndexedPoint)o);
    else
      return false;
  }

  /**
   * Removes a point by index.
   * @return true iff there was a point with such index.
   */
  public boolean remove(int point_index)
  {
    if (!contains(point_index))
      return false;
    checkConcurrency();
    data[point_index]=null;
    number_of_points--;
    return true;
  }

  /**
   * Tests if a point is in the set.
   * @return true iff it's there.
   */
  public boolean contains(IndexedPoint p)
  {
    if (p==null) return false;
    int idx=p.getIndex();
    return (contains(idx) && p.equals(data[idx]));
  }

  /**
   * Tests if a point is in the set.
   * @return true iff o is an IndexedPoint and it's in the set.
   */
  public boolean contains(Object o)
  {
    if (o instanceof IndexedPoint)
      return contains((IndexedPoint)o);
    else
      return false;
  }

  /**
   * Tests if there is a point int the set with the given index.
   */
  public boolean contains(int point_index)
  {
    return (number_of_elements>point_index && point_index>=0 && data[point_index]!=null);
  }

  private void ensureNewElementCanBeAdded(int i)
  {
    if (i>=capacity){
      // expand array
      expandData(2*i);
    }
    number_of_elements = Math.max(i+1,number_of_elements);
  }

  private void expandData(int newsize)
  {
    IndexedPoint[] newarray=new IndexedPoint[newsize];
    System.arraycopy(data,0,newarray,0,capacity);
    data=newarray;
    capacity=newsize;
  }

  public class SimplePointSetIterator implements Iterator 
  {
    private int at_index=-1;
    private int num_points_left=number_of_points;

    public SimplePointSetIterator()
    {
      number_of_iterators++;
    }

    public boolean hasNext()
    {
      checkModified();
      return (num_points_left>0);
    }

    public Object next()
    {
      checkModified();
      if (num_points_left==0)
        throw new java.util.NoSuchElementException("There are no more points in this set.");
      // find next non-null
      do {at_index++;} while(at_index<number_of_elements && data[at_index]==null);
      num_points_left--;
      return data[at_index];
    }

    public void remove()
    {
      checkModified();
      if (at_index<0)
        throw new IllegalStateException("Cannot remove() before starting the iteration.");
      if (data[at_index]==null)
        throw new IllegalStateException("Point already removed.");

      data[at_index]=null;
      number_of_points--;
    }

    private void checkModified()
    {
      if (modified_while_iterating)
        throw new java.util.ConcurrentModificationException("Set was modified.");
    }

    public void finalize()
    {
      // this iterator is not used by anyone anymore
      number_of_iterators--;
    }

  }
}