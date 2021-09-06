/*
 * CircularSector.java
 *
 * Created on November 15, 2007, 9:49 AM
 */

package ca.umontreal.iro.evolution.malin;

/**
 *
 * @author  csuros
 */

import java.awt.geom.Arc2D;


public class CircularSector extends Arc2D.Double
{
    
    public CircularSector()
    {
        this(0.,0.,0.,0.,0.);
    }
    
    /**
     * Constructs a new circular sector, with the given center, radius, and angle range
     */
    public CircularSector(double x, double y, double radius, double start_angle, double end_angle) 
    {
        super(x-radius,y-radius,2.*radius,2.*radius,start_angle, end_angle-start_angle, Arc2D.PIE);
    }
    
    public void setCircularSector(double x, double y, double radius, double start_angle, double end_angle)
    
    {
        setArc(x-radius, y-radius, 2*radius, 2*radius, start_angle, end_angle-start_angle,Arc2D.PIE);
    }
    

    /**
     * Sets the circle radius
     */
    public void setRadius(double R)
    {
        double x0=x+0.5*width;
        double y0=y+0.5*height;
        width=height=2.*R;
        x=x0-R;
        y=y0-R;
    }
    

    /**
     * Sets the circle center
     */
    public void setCenterX(double x)
    {
        double diff = x-getCenterX();
        this.x += diff;
    }

    /**
     * Sets the circle center
     */
    public void setCenterY(double y)
    {
        double diff = y-getCenterY();
        this.y += diff;
    }

    
    /**
     * Returns the center of the circle this sector is inscribed to
     */
    public double getCenterX()
    {
        return this.x+getRadius();
    }
    
    /**
     * Returns the center of the circle this sector is inscribed to
     */
    public double getCenterY()
    {
        return this.y+getRadius();
    }

    /**
     * Returns the radius of the circle
     */
    public double getRadius()
    {
        return 0.5*this.width;
    }

    /**
     * Returns the end point of the angle range
     */
    public double getAngleEnd()
    {
        return getAngleStart()+getAngleExtent();
    }
    
    /**
     * Sets the end point of the angle range; start remains the same
     */
    public void setAngleEnd(double deg)
    {
        double xt = deg-getAngleStart();
        super.setAngleExtent(xt);
    }
    
    /**
     * Sets the starting point of the angle range; end remains the same
     */
    public void setAngleStart(double deg)
    {
        double e = getAngleEnd();
        super.setAngleStart(deg);
        setAngleEnd(e);
    }
    
    /**
     * Cosine of an angle given in degrees
     */
    public static final double cos(double deg)
    {
        return Math.cos(Math.toRadians(deg));
    }
    
    /**
     * Arc cosine in degrees.
     */
    public static final double acos(double x)
    {
        return Math.toDegrees(Math.acos(x));
    }
    
    /**
     * Sine of an angle given in degrees
     */
    public static final double sin(double deg)
    {
        return Math.sin(Math.toRadians(deg));
    }
    
    /**
     * Computes the center's distance and angle from a given point.
     * @return two-element array of (distance, angle); angle is in degrees, between -180 and 180 
     */
    public double[] centerDirectionFrom(double x, double y)
    {
        double dx = getCenterX()-x;
        double dy = getCenterY()-y;
        double l = Math.sqrt(dx*dx+dy*dy);
        if (l==0.0)
        {
            return new double[2];
        } else
        {
            double cos_a = dx/l;
            double sin_a = dy/l;
            double angle = acos(cos_a);
            if (sin_a<0.)
                angle= -angle;
            double[] retval = new double[2];
            retval[0]=l;
            retval[1]=angle;
            return retval;
        }
    }
        
    /**
     * Finds the intersection of the underlying circle (ie, within which this sector is inscribed),
     * and a half-line. Notice that only the circle's radius is affecting the computations (x,y,angles do not).
     *
     * @param line_start_distance how far the half-line's starting point is from the underlying circle's center
     * @param line_start_angle direction of the circle center from the line's starting point
     * @param angle direction of the half-line (angle measured in degrees, counter-clockwise, with y==0 line being 0 degrees)
     * @return a zero- one- or two-element array (aneg, apos) for the angles of the intersection points (0, 1 or two) viewed from the underlying circle's center
     *     for two intersection points; values are between -180 and +180 degrees
     */
    public double[] intersectionHalfLine(double line_start_distance, double line_start_angle, double angle)
    {
        double R = getRadius();
        double diff_angle = line_start_angle-angle;
        double b = line_start_distance*cos(diff_angle);
        double sin_sq = sin(diff_angle);
        double sq = R*R-line_start_distance*line_start_distance*sin_sq*sin_sq;
            
        // Now the distance t from the line start
        // is b +- sqrt(sq)
        if (sq<0.0)
        {
            // no hit here : guy fits entirely at the end angle 
            return new double[0];
        }
        else if (sq>0.0)
        {
            double t1 = b-Math.sqrt(sq); // the farther hit
            double t2 = b+Math.sqrt(sq); // the closer hit
            
            double cos_a1 = (t1*cos(angle)-line_start_distance*cos(line_start_angle))/R;
            double cos_a2 = (t2*cos(angle)-line_start_distance*cos(line_start_angle))/R;
            double sin_a1 = (t1*sin(angle)-line_start_distance*sin(line_start_angle))/R;
            double sin_a2 = (t2*sin(angle)-line_start_distance*sin(line_start_angle))/R;
            
            double a1 = acos(cos_a1);

            //System.out.println("#**CS.iHL lsd "+line_start_distance+" lsa "+line_start_angle+" a "+angle+";\tt1 "+t1+" t2 "+t2+"\tca1 "+cos_a1+" ca2 "+cos_a2);
            
            
            if (sin_a1<0.)
            {
                // negative sine...
                a1 = -a1;
            }
            double a2 = acos(cos_a2);
            if (sin_a2<0.)
            {
                // negative sine...
                a2 = -a2;
            }

            double[] retval = new double[2];
            retval[0]=a1;
            retval[1]=a2;
            //System.out.println("#**CS.iHL ["+cos_a1+", "+sin_a1+"] "+a1+"\t["+cos_a2+", "+sin_a2+"] "+a2+"\t// "+angle+"\t"+sin(angle)+", "+sin(line_start_angle));            
            
            return retval;
        } else
        {
            // sq==0
            double t = b;
            double cos_a = (t*cos(angle)-line_start_distance*cos(line_start_angle))/R;
            double a = acos(cos_a);
            if (t*sin(angle)<line_start_distance*sin(line_start_angle))
            {
                // negative sine...
                a = 360.0-a;
            }
            double[] retval = new double[1];
            retval[0]=a;
            return retval;
        }
    }
    
    public double interSectionPoint(double line_start_distance, double line_start_angle, double angle, double angle_from_center)
    {
        double cos_a = cos(angle);
        if (cos_a != 0.0)
        {
            double dx = line_start_distance*cos(line_start_angle)+getRadius()*cos(angle_from_center);
            double a = dx / cos(angle);
            return a;
        } else 
        {
            double dy = line_start_distance*sin(line_start_angle)+getRadius()*sin(angle_from_center);
            return dy/sin(angle);
        }
        
    }
    
    public double intersectionSector(double distance, double angle, CircularSector sector)
    {
        double R = sector.getRadius();
        if (R>distance)
        {
            if (angle<90.0)
            {
                double beta = sector.getAngleStart();
                double t = Math.sqrt(R*R+distance*distance-2.*R*distance*cos(beta-angle));
                double cos_phi = (R*cos(beta)-distance*cos(angle))/t;
                double sin_phi = (R*sin(beta)-distance*sin(angle))/t;
                double phi = acos(cos_phi);
                if (sin_phi<0.0)
                    phi=-phi;
                return phi;
            } else 
            {
                double beta = sector.getAngleEnd();
                double t = Math.sqrt(R*R+distance*distance-2.*R*distance*cos(beta-angle));
                double cos_phi = (R*cos(beta)-distance*cos(angle))/t;
                double sin_phi = (R*sin(beta)-distance*sin(angle))/t;
                double phi = acos(cos_phi);
                if (sin_phi<0.0)
                    phi=-phi;
                return phi;                
            }
        }
        double cos_delta = R/distance;
        double delta = acos(cos_delta);
        double alpha = delta-angle; // angle to the tangent
        if (angle>90.0)
            alpha = -alpha;
        if (alpha<sector.getAngleStart())
        {
            double beta = sector.getAngleStart();
            double t = Math.sqrt(R*R+distance*distance-2.*R*distance*cos(beta-angle));
            double cos_phi = (R*cos(beta)-distance*cos(angle))/t;
            double sin_phi = (R*sin(beta)-distance*sin(angle))/t;
            double phi = acos(cos_phi);
            if (sin_phi<0.0)
                phi=-phi;
            return phi;
        } else if (alpha>sector.getAngleEnd())
        {
            double beta = sector.getAngleEnd();
            double t = Math.sqrt(R*R+distance*distance-2.*R*distance*cos(beta-angle));
            double cos_phi = (R*cos(beta)-distance*cos(angle))/t;
            double sin_phi = (R*sin(beta)-distance*sin(angle))/t;
            double phi = acos(cos_phi);
            if (sin_phi<0.0)
                phi=-phi;
            return phi;
        } else
        {
            if (angle>90.0)
                return alpha - 90.0;
            else 
                return 180.0-alpha;
        }
    }
    
    protected String paramString()
    {
        StringBuffer sb = new StringBuffer();
        sb.append("(");
        sb.append(getCenterX());
        sb.append(",");
        sb.append(getCenterY());
        sb.append(" R=");
        sb.append(getRadius());
        sb.append(", ");
        sb.append((int)getAngleStart());
        sb.append(":");
        sb.append((int)getAngleEnd());
        sb.append(" (");
        sb.append(Math.toRadians(getAngleStart()));
        sb.append(":");
        sb.append(Math.toRadians(getAngleEnd()));
        sb.append(")");
        return sb.toString();
    }
    
    public String toString()
    {
        StringBuffer sb = new StringBuffer();
        sb.append("CSector[");
        sb.append(paramString());
        sb.append("]");
        return sb.toString();
    }
}
