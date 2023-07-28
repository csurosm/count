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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Stroke;
import javax.swing.Icon;

import count.matek.DiscreteGamma;
import count.matek.Functions;



/**
 * A plot (Icon) illustrating the density function for a gamma distribution and
 * its discrete approximation.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 */
public final class DiscreteGammaPlot implements Icon
{
    public static int GAMMA_PLOT_WIDTH = 320;
    public static int GAMMA_PLOT_HEIGHT = 60;
    public static int GAMMA_PLOT_PARTITIONS = 8;

    public DiscreteGammaPlot(DiscreteGamma distribution, int num_categories)
    {
        this(distribution, num_categories, GAMMA_PLOT_WIDTH, GAMMA_PLOT_HEIGHT);
    }
    
    
    public DiscreteGammaPlot(DiscreteGamma distribution, int num_categories, Color category_color)
    {
        this(distribution, num_categories);
        setCategoryColor(category_color);
        setPlottingColor(category_color.brighter().brighter());
    }

    public DiscreteGammaPlot(DiscreteGamma distribution, int num_categories, int width, int height)
    {
        setDistribution(distribution);
        setNumCategories(num_categories);
        setWidth(width);
        setHeight(height);
        background_color=legend_color=plotting_color= null;
    }

    private DiscreteGamma distribution;
    private int num_categories;

    private int width;
    private int height;

    /**
     * Background color
     */
    private Color background_color;
    /**
     * Color used for drawing the axes and parameter info
     */
    private Color legend_color;
    /**
     * Color used for plotting the distribution
     */
    private Color plotting_color;

    /**
     * Color used for plotting the categories
     */
    private Color category_color;

    public final void setWidth(int width){ this.width = width;}
    public final void setHeight(int height){this.height = height;}

    @Override
    public int getIconWidth(){return width;}
    @Override
    public int getIconHeight(){return height;}


    /**
     * Sets the bakground color for the plot. By default, this is null, meaning that 
     * no background is painted (transparent icon).
     * 
     * @param c a color (null means default behavior)
     */
    public void setBackground(Color c)
    {
        this.background_color = c;
    }

    /**
     * Sets the color for the axes, and parameter info. By 
     * default, this is the normal color for the Graphics 
     * context when paintComponent() is called.  
     * 
     * @param c a color (null means default behavior)
     */
    public void setLegendColor(Color c)
    {
        this.legend_color = c;
    }

    /**
     * Sets the plotting color for the distribution. By 
     * default, this is the same as the legend color   
     * 
     * @param c
     */
    public void setPlottingColor(Color c)
    {
        this.plotting_color = c;
    }

    /**
     * Sets the plotting color for the categories. By 
     * default, this is the same as the plotting color
     * 
     * @param c plotting color
     */
    public void setCategoryColor(Color c)
    {
        this.category_color = c;
    }

    /**
     * Sets the distribution that is to be plotted here.
     * 
     * @param D a discrete Gamma distribution 
     */
    public final void setDistribution(DiscreteGamma D){distribution = D;}    

    /**
     * Sets the number of partitions used in the discretization
     * 
     * @param num_categories number of partitions
     */
    public final void setNumCategories(int num_categories){this.num_categories=num_categories;}

    /**
     * Plots the distribution at a given location. 
     * 
     * @param ignored this parameter is ignored, may be even null (required parameter for Icon interface)
     * @param g graphics context (should be convertible to Graphics2D)
     * @param x coordinate for the left-hand side
     * @param y coordinate for the top
     */
    @Override
    public void paintIcon(java.awt.Component ignored,
               Graphics g,
               int x,
               int y)
    {
        Graphics2D g2 = (Graphics2D)g.create(); // we will work with a copy
        Color original_color = g2.getColor();
        Stroke original_stroke=g2.getStroke();
        // draw background
        if (background_color != null)
        {
            g2.setColor(background_color);
            g2.fillRect(x, y, width, height);
            g2.setColor(original_color);
        }

        // now comes the legend: X, Y axes and tics
        if (legend_color != null)
            g2.setColor(legend_color);

        // parameters for the plot
        int font_size = g2.getFont().getSize();
        int tic_length = 3;

        // how much pace is needed outside the axes
        int axis_offset_x = tic_length+g2.getFontMetrics().stringWidth("9e-9")+1;
        int axis_offset_y = tic_length+font_size*6/5;

        // position of the origin
        int origin_x = x+axis_offset_x+1;
        int origin_y = y+height-axis_offset_y+1;

        // size of the the actual plot area 
        int plot_width  = width - axis_offset_x - tic_length-1;
        int plot_height = height- axis_offset_y - tic_length-1;

        // draw the axes
        g2.drawLine(origin_x-1, origin_y+1, origin_x-1, origin_y-plot_height);
        g2.drawLine(origin_x-1, origin_y+1, origin_x+plot_width, origin_y+1);

        // get the discretization for the distribution 
        double[] partition_boundary = new double[num_categories-1];
        double[] partition_center = distribution.getPartitionMeans(num_categories, partition_boundary);

        // horizontal mapping for plot area 
        double max_value_x = partition_center[num_categories-1]*1.05; // leave a little bit after the last quantile's mean
        double scale_x = plot_width/max_value_x; // x is plotted at coordinate origin_x+scale_x*x

        // draw the tics
        double unit_x = Functions.roundToMostSignificantDigit(max_value_x/6.0,null, true);  // just a few tics 

        for (int i=1; i*unit_x<max_value_x; i++)
        {
            double r = i*unit_x; // this is the position of this tic
            int tic_x = (int)(r*scale_x); // transform to display coordinates
            g2.drawLine(origin_x+tic_x, origin_y, origin_x+tic_x, origin_y+tic_length);
            if (i==1)
                DrawString.drawCentered(g2, Float.toString((float)r), origin_x+tic_x, origin_y+tic_length+font_size);
        }


        // mode is at x=(alpha-1)/alpha, but here the actual plotted maximum is computed instead
        double max_value_y = 0.0; // will be positive

        // computed the points to be plotted
        double[] pdf = new double[plot_width];
        {
            double alpha = distribution.getAlpha();
            double d = alpha * Math.log(alpha)-Functions.gammln(alpha);
            for (int i=0; i<plot_width; i++)
            {
                // conversion from display coordinate to x
                double r = (i+1.0)/scale_x;
                double ln_f = d-alpha*r+(alpha-1.0)*Math.log(r);
                double f = Math.exp(ln_f); // this is the value of the density function at r
                pdf[i] = f;
                //System.out.println("#*RMD.DGP pI "+i+"\tr "+r+"\tf "+f);
                if (f>max_value_y) max_value_y=f;
            }
        }
        double scale_y = plot_height/max_value_y;


        // plot the distributions
        if (plotting_color != null)
            g2.setColor(plotting_color);
        g2.setStroke(new BasicStroke(1.f));
        for (int i=0; i<plot_width; i++)
        {
            int bar_height = (int)(scale_y*pdf[i]);
            //System.out.println("#*RMD.DGP pI "+i+"\tf "+pdf[i]+"\tbh "+bar_height+"\t// scale "+scale_y+"\tph "+plot_height+"\tmaxy "+max_value_y);
            g2.drawLine(origin_x+i, origin_y, origin_x+i, origin_y-bar_height);
        }

        // plot categories
        if (category_color != null)
            g2.setColor(category_color);

        g2.setStroke(new BasicStroke(2.f));
        for (int i=0; i<partition_center.length; i++)
        {
            int loc_x = (int)(partition_center[i]*scale_x);
            g2.drawLine(origin_x+loc_x, origin_y, origin_x+loc_x, origin_y-2*plot_height/4);
        }

        g2.setStroke(original_stroke);
        g2.setColor(original_color);

        if (legend_color != null)
            g2.setColor(legend_color);

        g2.drawString(num_categories+" categories, \u03b1="+distribution.getAlpha() //DoubleRoundedForDisplay.toString(distribution.getAlpha())+"\u22ef"
                , origin_x+50, origin_y-plot_height+font_size);
    }
}