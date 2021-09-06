/*
 * PointSetPanel.java
 *
 * Created on November 14, 2007, 12:10 AM
 */

package ca.umontreal.iro.evolution.malin.ui;

/**
 * A JPanel tracking a PointSet.
 *
 * @author  csuros
 */


import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;

import java.util.List;

import ca.umontreal.iro.evolution.malin.IndexedPoint;
import ca.umontreal.iro.evolution.malin.IndexedPointSet;

public abstract class PointSetPanel extends JPanel 
{

  private static int next_id;
  private int id;

  protected IndexedPointSet point_set;

  protected boolean points_moved;
  protected boolean window_resized;

  public PointSetPanel() 
  {
      this(true);
  }
  
  public PointSetPanel(boolean selectable_points)
  {
      this(selectable_points,selectable_points);
  }
  
  public PointSetPanel(boolean selectable_points, boolean selectable_area)
  {
    super(true);
    this.id=next_id++;
    initPointSet();

    if (selectable_points)
    {
        addMouseListener(new PointSelector());
    }
    if (selectable_area)
    {
        AreaSelectionTracker areaSelector = new AreaSelectionTracker();
        addMouseListener(areaSelector);
        addMouseMotionListener(areaSelector);
    }
    setToolTipText("");

    points_moved=true;
    window_resized=false;
      
  }


  public IndexedPoint closestPoint(int x, int y, double max_distance)
  {
    if (points_moved)
      rebuildPointSet();
    else if (window_resized)
      stretchPointSet();
    points_moved=window_resized=false;

    if (point_set==null)
      return null;
    return point_set.closestPoint(x,y,max_distance);
  }

  /**
   * Allocates <code>point_set</code>.
   * (Subclasses implement this.)
   */
  protected abstract void initPointSet();

  /**
   * Sets up <code>point_set</code> after points moved.
   */
  protected abstract void rebuildPointSet();

  /**
   * Sets up <code>point_set</code> after window is resized
   * but points haven't moved otherwise. Default implementation calls rebuildPointSet().
   */
  protected void stretchPointSet()
  {
    rebuildPointSet();
  }

  // ------------------------------------------
  // ---
  // ---  Area selection
  // ---
  // ------------------------------------------
  /**
   * Rectangle describing the currently selected area.
   */
  protected Rectangle currentRect = null;

  /**
   * Version of <code>currentRect</code> for drawing.
   * (width and height are always positive and within the enclosing component).
   */
  protected Rectangle rectToDraw = null;

  /**
   * Color used for displaying the selection.
   */
  private Color selectedRectangleColor=Color.gray;

  public void setSelectionAreaColor(Color c)
  {
    selectedRectangleColor=c;
  }

  /**
   * Draws the rectangle specified by <code>drawToRect</code>.
   */
  protected void drawSelectedRectangle(Graphics g)
  {
    Color c=g.getColor();
    g.setColor(selectedRectangleColor);
    g.fillRect(rectToDraw.x, rectToDraw.y,
                  rectToDraw.width - 1, rectToDraw.height - 1);
    g.setColor(c);
  }

  // window size
  private int prevWidth=-1;
  private int prevHeight=-1;

  /**
   * Draws the background (including the selection rectangle).
   * Sets window_resized.
   */
  public void paintComponent(Graphics g) 
  {
    super.paintComponent(g); //paints the background and image

    //If currentRect exists, paint a box on top.
    if (currentRect != null) 
      drawSelectedRectangle(g);

    // check if the window is resized
    if (getWidth()!=prevWidth || getHeight() != prevHeight)
    {
      window_resized=true;
      if (prevWidth==0 && prevHeight==0)
        points_moved=true; // setup if window wasn't init'd
      prevWidth=getWidth();
      prevHeight=getHeight();
    }
  }

  /**
   * Removes the area selection.
   */
  public void deSelectRectangle()
  {
    if (currentRect != null)
    {
      Rectangle totalRepaint = currentRect;
      currentRect=null;
      this.repaint(totalRepaint);
    }
  }

  /**
   * @return null if no area selection.
   */
  public Rectangle getSelectedRectangle()
  {
    return currentRect;
  }

  /**
   * @return <code>List</code> of </code>IndexedPoints</code> within the selected area.
   * (may return <code>null</code>)
   */
  public List<IndexedPoint> getSelectedPoints()
  {
    if (currentRect==null || point_set == null)
      return null;
    return point_set.withinRectangle(currentRect);
  }

  private double close_point_radius_sqr=9.;

  /**
   * Sets the radius for finding closest points (used for
   * mouse clicks and tool tip text).
   */
  public void setCloseRadius(double r)
  {
    close_point_radius_sqr=r*r;
  }
  
  public double getClosePointRadius()
  {
      return Math.sqrt(close_point_radius_sqr);
  }


  /**
   * Handles point selection by calling &quot;hooks&quot;
   * such as <code>selectPoint</code>, <code>selectPoints</code>
   * and <code>removeSelection()</code>
   */
  class PointSelector extends java.awt.event.MouseAdapter 
  {
    public void mouseClicked(MouseEvent e){
      if (!SwingUtilities.isLeftMouseButton(e)) return;
      int x=e.getX();
      int y=e.getY();

      IndexedPoint p=closestPoint(x,y,close_point_radius_sqr);
      {
          //System.out.println("#**PSP.PS.mC "+p+"\t"+e);
      }
      if (p!=null)
        selectPoint(p,e.getClickCount(),e.isControlDown());
      else if (!e.isControlDown())
        removeSelection();
    }

    public void mouseReleased(MouseEvent e)
    {
      if (!SwingUtilities.isLeftMouseButton(e)) return;
      Rectangle r=getSelectedRectangle();
      if (r==null)
          return;

      if (r.getWidth()<0)
      {
        r.setLocation((int)(r.getX()+r.getWidth()), (int)r.getY());
        r.setSize(-(int)r.getWidth(),(int)r.getHeight());
      }
      if (r.getHeight()<0)
      {
        r.setLocation((int)r.getX(), (int)(r.getY()+r.getHeight()));
        r.setSize((int)r.getWidth(),-(int)r.getHeight());
      }

      if (r.getWidth()!=0. && r.getHeight()!=0.)
      {
        selectPoints(getSelectedPoints(),e.isControlDown());
        deSelectRectangle();
      }
    }
  }

  /**
   * Hook for subclasses: called when a single point is selected.
   * @param P is not <code>null</code>.
   */
  protected abstract void selectPoint(IndexedPoint P, int num_mouse_clicks, boolean isCtrlDown);

  /**
   * Hook for subclasses: called when several points are selected in an area.
   * @param L <code>List</code> of <code>IndexedPoint</code>s [may be <code>null</code>.
   */
  protected abstract void selectPoints(List<IndexedPoint> L, boolean is_Ctrl_down);

  /**
   * Hook for subclasses: called when selection is removed.
   */
  protected abstract void removeSelection();

  private int previous_tool_tip_x=-1;
  private int previous_tool_tip_y=-1;
  private String previous_tool_tip_text="";

  public String getToolTipText(MouseEvent e)
  {
    int x=e.getX(), y=e.getY();
    if (x==previous_tool_tip_x && y==previous_tool_tip_y)
    {
      // called more than once even though mouse isn't moving
      return previous_tool_tip_text;
    } else 
    {
      previous_tool_tip_x=x;
      previous_tool_tip_y=y;
      String retval = "";

      //System.out.println("#**GP.gTT "+e.toString());

      IndexedPoint p=closestPoint(e.getX(),e.getY(),close_point_radius_sqr);

      retval = getToolTipTextForPoint(x,y,p);

      previous_tool_tip_text=retval;
      return retval;
    }
  }

  /**
   * Hook for subclasses: called for displaying tool tip.
   * @param x mouse location
   * @param y mouse location
   * @param p closest point [may be <code>null</code>]
   */
  protected abstract String getToolTipTextForPoint(int x, int y, IndexedPoint p);

  /**
   * This class tracks <code>MouseEvent</code>s
   * to update <code>currentRect</code> and <code>rectToDraw</code>,
   * and paint whenever necessary.
   */
  class AreaSelectionTracker extends javax.swing.event.MouseInputAdapter 
  {
    private Rectangle previousRectDrawn = new Rectangle();

    public void mousePressed(MouseEvent e) 
    {
      if (!SwingUtilities.isLeftMouseButton(e)) return;

      //System.out.println("#**GP.GPAS.mP "+e.toString()+" w "+getWidth()+" h "+getHeight());
      int x = e.getX();
      int y = e.getY();
      currentRect = new Rectangle(x, y, 0, 0);
      updateDrawableRect(getWidth(), getHeight());
      repaint();
    }

    public void mouseDragged(MouseEvent e) 
    {
      //System.out.println("#**GP.GPAS.mD "+e.toString());
      if (!SwingUtilities.isLeftMouseButton(e)) return;
      updateSize(e);
    }

    public void mouseReleased(MouseEvent e) 
    {
      if (!SwingUtilities.isLeftMouseButton(e)) return;
      updateSize(e);
    }

    /*
     * Update the size of the current rectangle
     * and call repaint.  Because currentRect
     * always has the same origin, translate it
     * if the width or height is negative.
     *
     * For efficiency (though
     * that isn't an issue for this program),
     * specify the painting region using arguments
     * to the repaint() call.
     *
     */

    /**
     * Called from AreaSelector.
     */
    private void updateSize(MouseEvent e) 
    {
      int x = e.getX();
      int y = e.getY();
      currentRect.setSize(x - currentRect.x,y - currentRect.y);
      updateDrawableRect(getWidth(), getHeight());
      Rectangle totalRepaint = rectToDraw.union(previousRectDrawn);


      // if just got focus beacuse of this selection, then repaint()
      // is not enough
      if ((e.getModifiers() & MouseEvent.MOUSE_RELEASED) != 0)
        paintImmediately(totalRepaint.x, totalRepaint.y,
                totalRepaint.width, totalRepaint.height);
      else
        repaint(totalRepaint.x, totalRepaint.y,
                totalRepaint.width, totalRepaint.height);
      //System.out.println("#**GP.GPAS.uS tR="+totalRepaint.toString()+"; e="+e.toString());
    }

    private void updateDrawableRect(int compWidth, int compHeight) 
    {
      int x = currentRect.x;
      int y = currentRect.y;
      int width = currentRect.width;
      int height = currentRect.height;

      //requestFocus(); // just in case

      /*System.out.println("#**GP.uDR "+currentRect.toString()
        +"; w "+width+", h "+height+", x "+x+", y "+y
        +"; cW "+compWidth+", cH "+compHeight);*/

      //Make the width and height positive, if necessary.
      if (width < 0) 
      {
          width = 0 - width;
          x = x - width + 1;
          if (x < 0) 
          {
              width += x;
              x = 0;
          }
      }
      if (height < 0) 
      {
          height = 0 - height;
          y = y - height + 1;
          if (y < 0) 
          {
              height += y;
              y = 0;
          }
      }

      //The rectangle shouldn't extend past the drawing area.
      if ((x + width) > compWidth) 
      {
          width = compWidth - x;
      }
      if ((y + height) > compHeight) 
      {
          height = compHeight - y;
      }

      //Update rectToDraw after saving old value.
      if (rectToDraw != null) 
      {
          previousRectDrawn.setBounds(
                      rectToDraw.x, rectToDraw.y,
                      rectToDraw.width, rectToDraw.height);
          rectToDraw.setBounds(x, y, width, height);
      } else 
      {
          rectToDraw = new Rectangle(x, y, width, height);
      }

      //System.out.println("#**GP.uDR rTD="+rectToDraw.toString());
    }

  }
}