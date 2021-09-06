
package ca.umontreal.iro.evolution.malin.ui;

/**
 * Static utility methods for displaying 
 * and aligning text in different manners. 
 *
 * @author csuros
 */

import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;

import java.awt.geom.Rectangle2D;
import java.awt.geom.AffineTransform;
import java.awt.Shape;
import java.awt.geom.NoninvertibleTransformException;

public class DrawString 
{

    /**
     * No instantiation.
     */
    private DrawString() {}
    

    /**
     * Draws a centered  String at a given coordinate.
     * 
     * @param g the Graphics context for drawing
     * @param text the String to be drawn
     * @param x X coordinate for the center of the drawn string
     * @param y Y coordinate for the baseline 
     */
    public static void drawCenteredString(Graphics g, String text, int x, int y)
    {
        FontMetrics fm = g.getFontMetrics();    
        int w = fm.stringWidth(text);
        g.drawString(text,x-w/2,y);
    }

    /**
     * Draws a rotated String at a given coordinate.
     * 
     * @param g the Graphics context for drawing
     * @param text the String to be drawn
     * @param theta angle for the drawing (in radians, clockwise, 0 is horizontal in g)
     * @param x X coordinate for the anchor point of the drawn string
     * @param y Y coordinate for the baseline at the anchor point
     * @param anchor_point relative position of the anchor point with respect to (x,y), given 
     *       as a fraction of the width (so 1.0f means right-aligned and 0.5f is centered)
     */
    public static void drawRotatedString(Graphics2D g, String text, int x, int y, double theta, float anchor_point)
    {
        //Rectangle2D R = getBoundingBoxForRotatedString(g,text,x,y,theta,anchor_point);
        //g.draw(R);
        AffineTransform old_transform = new AffineTransform(g.getTransform());
        g.translate(x,y);
        g.rotate(theta);
        int w = g.getFontMetrics().stringWidth(text);
        g.drawString(text,-anchor_point*w,0);
        g.setTransform(old_transform);
    }

    /**
     * Draws a rotated String at a coordinate specifying the lower left corner.
     * 
     * @param g the Graphics context for drawing
     * @param text the String to be drawn
     * @param theta angle for the drawing (in radians, clockwise, 0 is horizontal in g)
     * @param x X coordinate for the left-hand side of the drawn string
     * @param y Y coordinate for the baseline on the left-hand side
     */
    public static void drawRotatedStringLeft(Graphics2D g, String text, int x, int y, double theta)
    {
        drawRotatedString(g, text, x, y, theta, 0.0f);
    }

    /**
     * Draws a rotated String at a coordinate specifying the lower right corner.
     * 
     * @param g the Graphics context for drawing
     * @param text the String to be drawn
     * @param theta angle for the drawing (in radians, counterclockwise, 0 is horizontal in g)
     * @param x X coordinate for the right side of the drawn string
     * @param y Y coordinate for the baseline on the right-hand side
     */
    public static void drawRotatedStringRight(Graphics2D g, String text, int x, int y, double theta)
    {
        drawRotatedString(g, text, x, y, theta, 1.0f);
    }   
    
    /**
     * Draws a rotated String at a coordinate specifying the center bottom.
     * 
     * @param g the Graphics context for drawing
     * @param text the String to be drawn
     * @param theta angle for the drawing (in radians, clockwise, 0 is horizontal in g)
     * @param x X coordinate for the center of the drawn string
     * @param y Y coordinate for the baseline at the center
     */
    public static void drawRotatedStringCenter(Graphics2D g, String text, int x, int y, double theta)
    {
        drawRotatedString(g, text, x, y, theta, 0.5f);
    }   
    
    /**
     * Computes the bounding box a rotated String at a given coordinate.
     * 
     * @param g the Graphics context for drawing
     * @param text the String to be drawn
     * @param theta angle for the drawing (in radians, clockwise, 0 is horizontal in g)
     * @param x X coordinate for the anchor point of the drawn string
     * @param y Y coordinate for the baseline at the anchor point
     * @param anchor_point relative position of the anchor point with respect to (x,y), given 
     *       as a fraction of the width (so 1.0f means right-aligned and 0.5f is centered)
     * @return bounding box (in g's coordinates) for the rotated string
     */
    public static Rectangle2D getBoundingBoxForRotatedString(Graphics2D g, String text, int x, int y, double theta, float anchor_point)
    {
        FontMetrics fm = g.getFontMetrics();
        return getBoundingBoxForRotatedString(g, fm, text, x, y, theta, anchor_point);
    }
    /**
     * Computes the bounding box a rotated String at a given coordinate.
     * 
     * @param g the Graphics context for drawing
     * @param fm FontMetrics for the string
     * @param text the String to be drawn
     * @param theta angle for the drawing (in radians, clockwise, 0 is horizontal in g)
     * @param x X coordinate for the anchor point of the drawn string
     * @param y Y coordinate for the baseline at the anchor point
     * @param anchor_point relative position of the anchor point with respect to (x,y), given 
     *       as a fraction of the width (so 1.0f means right-aligned and 0.5f is centered)
     * @return bounding box (in g's coordinates) for the rotated string
     */
    public static Rectangle2D getBoundingBoxForRotatedString(Graphics2D g, FontMetrics fm, String text, int x, int y, double theta, float anchor_point)
    {
        Rectangle2D R = fm.getStringBounds(text, g);
        
        double w = R.getWidth();
        double h = R.getHeight();
        double Rx = R.getX();
        double Ry = R.getY();
        
        double sin_t = Math.sin(theta);
        double cos_t = Math.cos(theta); 
        
        double BB_x = 0.0; 
        double BB_y = 0.0; 
        double BB_w = 0.0;
        double BB_h = 0.0;
        
        if (sin_t<=0.)
        {
            if (cos_t>=0.)
            {
                // upper right quadrant
                BB_x = x - (anchor_point*w-Rx)*cos_t - Ry*sin_t; 
                BB_y = y + (w-anchor_point*w+Rx)*sin_t + Ry*cos_t;
                BB_w = w*cos_t - h*sin_t;
                BB_h = -w*sin_t + h*cos_t;
            } else
            {
                // upper left quadrant
                BB_x = x + (w-anchor_point*w+Rx)*cos_t - Ry*sin_t; 
                BB_y = y + (w-anchor_point*w+Rx)*sin_t + (h+Ry)*cos_t;
                BB_w = -w*cos_t - h*sin_t;
                BB_h = -w*sin_t - h*cos_t;
            }
        } else
        {
            if (cos_t<0)
            {
                // lower left quadrant
                BB_x = x + (w-anchor_point*w+Rx)*cos_t - (h+Ry)*sin_t; 
                BB_y = y - (anchor_point*w-Rx)*sin_t + (h+Ry)*cos_t;
                BB_w = -w*cos_t + h*sin_t;
                BB_h = w*sin_t - h*cos_t;
            } else
            {
                // lower right quadrant
                BB_x = x - (anchor_point*w-Rx)*cos_t - (h+Ry)*sin_t; 
                BB_y = y - (anchor_point*w-Rx)*sin_t + Ry*cos_t;
                BB_w = w*cos_t + h*sin_t;
                BB_h = w*sin_t + h*cos_t;
            }
        }
        
        //System.out.println("#**DS.gBB '"+text+"'\t@ ("+x+", "+y+")\tR "+R+"\tBB [x="+BB_x+" y="+BB_y+" w="+BB_w+" h="+BB_h+"]"+"\t//theta "+cos_t+", "+sin_t);
        
        Rectangle2D BB = new Rectangle2D.Double(BB_x, BB_y, BB_w, BB_h);
        //System.out.println("#**DS.gBB '"+text+"'\tBB "+BB);
        
        return BB;
    }
}
