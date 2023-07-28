package count.gui.kit;
/*
 * Copyright 2022 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.Set;
import java.util.stream.Collectors;

import javax.swing.DefaultListSelectionModel;
import javax.swing.JPanel;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.MouseInputAdapter;

/**
*
* A generic panel for a set of indexed points.
* 
* @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
* @since November 14, 2007, 12:10 AM
*/
public class PointSetPanel extends JPanel implements ListSelectionListener
{
    /**
     * A data structure for storing the points. 
     * 
     * (More sophisticated implementation would 
     * have a proper geometric data structure here.)
     */
    private final Set<IndexedPoint> point_set;

    /**
     * Rectangle describing the currently selected area.
     */
    private Rectangle current_selection;
    /**
     * Version of {@link #current_selection} for drawing.
     * Width and height are always positive and within the enclosing component.
     */
    private Rectangle drawn_selection;
    /**
     * Color used for displaying the selection.
     */
    private Color selection_area_color=Color.GRAY;

    /**
     * Default squared distance for finding points near a mouse click. 
     */
    private double close_point_radius_sqr=9.0;

    private ListSelectionModel point_selection;
    private MouseListener point_selection_handler=null;
    private AreaSelection area_selection_handler = null;
    
    private final TooltipSaver tooltip_info;
    
    
    /**
     * Instantiation with full set of arguments.  
     * 
     * @param point_set set of points stored here.
     * @param selection_mode one of the constants from {@link javax.swing.ListSelectionModel}
     * @param point_selection_enabled whether points can be selected by mouse clicks
     * @param area_selection_enabled whether mouse dragging is tracked by a rectangular selection area
     */
    public PointSetPanel(Set<IndexedPoint> point_set, int selection_mode, boolean point_selection_enabled, boolean area_selection_enabled)
    {
        super(true); // double buffered JPanel
        this.point_set = point_set;
        {
            ListSelectionModel selection_model = new DefaultListSelectionModel();
            selection_model.setSelectionMode(selection_mode);
            setSelectionModel(selection_model);
        }
        
        setPointSelectionEnabled(point_selection_enabled);
        setAreaSelectionEnabled(area_selection_enabled);
        setToolTipText(""); // so that tool tip text will be queried by mouse coordinates
        tooltip_info = new TooltipSaver();
    }
    
    
    /**
     * Instantiation with an empty point set.  
     * Add points later using {@link #point_set } and {@link java.util.Set#add(java.lang.Object) }, then repaint.
     * 
     * @param selection_mode one of the constants from {@link javax.swing.ListSelectionModel}
     * @param point_selection_enabled whether points can be selected by mouse clicks
     * @param area_selection_enabled whether mouse dragging is tracked by a rectangular selection area
     */
    public PointSetPanel(int selection_mode, boolean point_selection_enabled, boolean area_selection_enabled)
    {
        this(new java.util.HashSet<>(),selection_mode, point_selection_enabled, area_selection_enabled);
    }    

    
    /**
     * Paints the area within this panel. 
     * Points are not plotted (that is for the extending class), only the selection. 
     * 
     * @param g graphics context
     */
    @Override
    protected void paintComponent(Graphics g)
    {
        super.paintComponent(g);
        if (getSelectedArea()!=null)
            drawSelectedArea(g);
    }
    
    /**
     * Gets tooltip text via cached {@link #tooltip_info}.
     * @param e mouse position
     * @return tooltip text based on closest point 
     */
    @Override
    public String getToolTipText(MouseEvent e)
    {
        return tooltip_info.getText(e);
    }
    
    @Override
    public void valueChanged(ListSelectionEvent e)
    {
        repaint(); // show selected nodes
    }
    
    
    public Rectangle getSelectedArea()
    {
        return current_selection;
    }
    
    public final ListSelectionModel getSelectionModel()
    {
        return point_selection;
    }
    
    public final void setSelectionModel(ListSelectionModel model)
    {
        model.addListSelectionListener(this);
        this.point_selection = model;
    }

    
    /**
     * Finds the closest element to some coordinate.
     * Points outside the default radius {@link #close_point_radius_sqr} are ignored.
     * 
     * @param x X coordinate
     * @param y Y coordinate
     * @return the closest point in the set that is close enough, or null if everybody is too far
     */
    public IndexedPoint getClosestPoint(double x, double y)
    {
        double       min_dist  = this.close_point_radius_sqr;
        IndexedPoint min_point = null;
        for (IndexedPoint ipoint:point_set)
        {
            double d = ipoint.distanceSq(x, y);
            if (d<=min_dist)
            {
                min_dist = d;
                min_point = ipoint;
            }
        }
        return min_point;
    }
    
    public Set<IndexedPoint> coveredPoints(Shape shape)
    {
        return
            point_set.stream().filter((ipoint) -> (shape.contains(ipoint))).collect(Collectors.toSet()); // playing with map-reduce
    }
    
    /**
     * Points covered by the currently selected area.
     * 
     * @return list of points within the selected area, or null if no area selected.
     * 
     */
    public Set<IndexedPoint> coveredPoints()
    {
        if (this.current_selection==null || point_set == null)
            return null;
        else 
            return coveredPoints(current_selection);
    }
    
    /**
     * Adds a point to the underlying set.
     * 
     * @param point
     * @return value by {@link Set#add(Object)}.
     */
    public boolean addPoint(IndexedPoint point)
    {
    	return point_set.add(point);
    }
    
    /**
     * Whether the point with the given index is selected.
     * 
     * @param idx
     * @return 
     */
    public boolean isSelected(int idx)
    {
        return point_selection.isSelectedIndex(idx);
    }

    /**
     * Sets the color for the selection area.
     * 
     * @param color favorite color
     */
    public void setSelectionAreaColor(Color color)
    {
        this.selection_area_color=color;
    }
    
    /**
     * Sets the radius for finding closest points (used for
     * mouse clicks and tool tip text).
     * 
     * @param r (unsquared Euclidean) distance for being close
     */
    public void setCloseRadius(double r)
    {
        close_point_radius_sqr=r*r;
    }
    
    /**
     * Selects or adds a point to current selection.
     * 
     * @param ipoint
     * @param e if {@link java.awt.event.MouseEvent#isControlDown() } , then extends the current selection; double-click sets the selection
     */
    protected void selectPoint(IndexedPoint ipoint, MouseEvent e)
    {
        int idx = ipoint.getIndex();
        if (e.getClickCount()>1)
        {
            point_selection.setSelectionInterval(idx, idx);
        } else if (e.isControlDown())
        {
            point_selection.addSelectionInterval(idx, idx);
        } else
        {
            point_selection.setSelectionInterval(idx, idx);
        }
    }

    /**
     * Selects a set of points.
     * 
     * @param selected_point_set
     * @param e 
     */
    protected void selectPoints(Set<IndexedPoint> selected_point_set, MouseEvent e)
    {
        point_selection.setValueIsAdjusting(true);
        if (!e.isControlDown())
            point_selection.clearSelection();
        for (IndexedPoint ipoint: selected_point_set)
        {
            int idx = ipoint.getIndex();
            point_selection.addSelectionInterval(idx, idx);
        }
        point_selection.setValueIsAdjusting(false);
    }
    
    private void deselectArea()
    {
        if (current_selection!=null)
        {
            Rectangle repaint_area = current_selection;
            current_selection=null;
            repaint(repaint_area);
        }
    }
    
    
    public final void setPointSelectionEnabled(boolean enable)
    {
        if (enable)
        {
            if (point_selection_handler == null)
            {
            	point_selection_handler = new MouseAdapter()
        		{
                    /**
                     * Reacts to left-button clicks.
                     * 
                     * @param e 
                     */
                    @Override
                    public void mouseClicked(MouseEvent e)
            		{
	                    if (!SwingUtilities.isLeftMouseButton(e)) return;
	                    
	                    int x=e.getX();
	                    int y=e.getY();
	
	                    IndexedPoint p=getClosestPoint(x,y);
	                    if (p!=null)
	                    {
	                        selectPoint(p, e);
	                    }
	                    else if (!e.isControlDown())
	                    {
	                        point_selection.clearSelection();
	                    }
            		}
        		};
                this.addMouseListener(point_selection_handler);
            } else
            {
                // nothing to do if already enabled
            }
        } else
        {
            if (point_selection_handler!=null)
            {
                this.removeMouseListener(point_selection_handler);
                point_selection_handler=null;
            } else
            {
                // nothing to do if already disabled
            }
        }
    }
    
    public final void setAreaSelectionEnabled(boolean enable)
    {
        if (enable)
        {
            if (area_selection_handler==null)
            {
                area_selection_handler = new AreaSelection();
                this.addMouseListener(area_selection_handler);
                this.addMouseMotionListener(area_selection_handler);
            }
        } else
        {
            if (area_selection_handler != null)
            {
                this.removeMouseListener(area_selection_handler);
                this.removeMouseMotionListener(area_selection_handler);
                area_selection_handler = null;
            }
        }
    }
    
    private void drawSelectedArea(Graphics g)
    {
        Graphics2D g2 = (Graphics2D)g.create();
        g2.setColor(this.selection_area_color);
        g2.fill(drawn_selection);
    }
    
    /**
     * Hook for subclasses: called for displaying tool tip at a point. 
     * 
     * @param x mouse location
     * @param y mouse location
     * @param p closest point [may be <code>null</code>]
     * @return text to be displayed in tooltip
     */
    protected String getTooltipText(int x, int y, IndexedPoint p)
    {
        return Integer.toString(p.getIndex());
    }
     
    /**
     * Class for updating the selection area as the mouse is dragged.
     * 
     * @author csuros
     *
     */
    private class AreaSelection extends MouseInputAdapter
    {
        final Rectangle previous_drawn = new Rectangle();

        @Override
        public void mousePressed(MouseEvent e) 
        {
            if (!SwingUtilities.isLeftMouseButton(e)) return;
            int x = e.getX();
            int y = e.getY();
            current_selection = new Rectangle(x, y, 0, 0);
            updateDrawableRect(getWidth(), getHeight());
            repaint();
        }
        
        @Override
        public void mouseDragged(MouseEvent e) 
        {
            if (!SwingUtilities.isLeftMouseButton(e)) return;
            updateSize(e);
        }

        @Override
        public void mouseReleased(MouseEvent e) 
        {
            if (!SwingUtilities.isLeftMouseButton(e)) return;
            
            updateSize(e);

            Rectangle r=getSelectedArea();

            if (r==null)
                return;

            if (r.getWidth()<0) // flip width
            {
                r.setLocation((int)(r.getX()+r.getWidth()), (int)r.getY());
                r.setSize(-(int)r.getWidth(),(int)r.getHeight());
            }
            if (r.getHeight()<0) // flip height
            {
                r.setLocation((int)r.getX(), (int)(r.getY()+r.getHeight()));
                r.setSize((int)r.getWidth(),-(int)r.getHeight());
            }

            if (r.getWidth()!=0. && r.getHeight()!=0.) // for a regular mouse click, mousePressed sets a 0x0 rectangle mouseReleased brings us here
            {
                selectPoints(coveredPoints(),e);
                deselectArea();
            }

        }


        /**
         * Called from AreaSelector.
         */
        private void updateSize(MouseEvent e) 
        {
            /*
             * Update the size of the current rectangle
             * and call repaint.  Because {@link #current_selection}
             * always has the same origin, translate it
             * if the width or height is negative.
             *
             * For efficiency (though that isn't an issue here),
             * specify the painting region using arguments
             * to the repaint() call.
             *
             */
            int x = e.getX();
            int y = e.getY();
            
            
            current_selection.setSize(x - current_selection.x,y - current_selection.y);
            updateDrawableRect(getWidth(), getHeight());
            Rectangle area_repainted = drawn_selection.union(previous_drawn);


            // if just got focus because of this selection, then repaint()
            // is not enough
            if ((e.getModifiersEx() & MouseEvent.MOUSE_RELEASED) != 0)
                paintImmediately(area_repainted.x, area_repainted.y,
                    area_repainted.width, area_repainted.height);
            else
                repaint(area_repainted.x, area_repainted.y,
                    area_repainted.width, area_repainted.height);
      //System.out.println("#**GP.GPAS.uS tR="+totalRepaint.toString()+"; e="+e.toString());
        }

        private void updateDrawableRect(int compWidth, int compHeight) 
        {
            int x = current_selection.x;
            int y = current_selection.y;
            int width = current_selection.width;
            int height = current_selection.height;

            //requestFocus(); // just in case

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
            if (drawn_selection == null) 
            {
                drawn_selection = new Rectangle(x, y, width, height);
            } else
            {
                previous_drawn.setBounds(
                            drawn_selection.x, drawn_selection.y,
                            drawn_selection.width, drawn_selection.height);
                drawn_selection.setBounds(x, y, width, height);
            } 
        }
    }

    
    
    /**
     * Generates and caches toolip text based on closest 
     * point.
     * 
     * @author csuros
     *
     */
	private class TooltipSaver
	{
	    private int x=-1;
	    private int y=-1;
	    private String text="";
	    private IndexedPoint point = null;
	    
	    String getText(MouseEvent e)
	    {
	        int ex = e.getX();
	        int ey = e.getY();
	        
	        if (ex==x && ey==y)
	        {
	            return text;
	        } else
	        {
	            x=ex;
	            y=ey;
	            IndexedPoint pt = getClosestPoint(ex, ey);
	            if (pt==point) // same selection as before 
	                return text;
	            else // different point selected 
	            {
	                point = pt;
	                return text = getTooltipText(x, y, pt);
	            }
	        }
	    }
	}    
}
