package count.gui.kit;
/*
 * Copyright 2021 Mikl&oacute;s Cs&#369;r&ouml;s.
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
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Stroke;

import javax.swing.Icon;

import count.matek.DiscreteDistribution;
import count.matek.NegativeBinomial;
import count.matek.PointDistribution;
import count.matek.Poisson;
import count.matek.ShiftedGeometric;

/**
 * A plot (Icon) illustrating the point mass function for a discrete distribution.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 */
public class DiscreteDistributionPlot implements Icon
{
	/**
	 * Default width
	 */
    public static final int DISTRIBUTION_PLOT_WIDTH = 40;
    /**
     * Default height
     */
    public static final int DISTRIBUTION_PLOT_HEIGHT = 100;
    /**
     * Default range of values to plot
     */
    public static final int DISTRIBUTION_PLOT_RANGE = 6;
    /**
     * Default background color
     */
    public static final Color DISTRIBUTION_PLOT_BACKGROUND = new Color(204, 255, 102); // Honeydew // new Color(255,204,102); // Canteloupe
    /**
     * Fully parametrized instantiation. 
     * 
     * @param D a discrete distribution what we want to plot
     * @param width width of the plot (pixels)
     * @param height height of the plot (pixels)
     */
    public DiscreteDistributionPlot(DiscreteDistribution D, int width, int height)
    {
        setWidth(width);
        setHeight(height);
        setDistribution(D);
        background_color=legend_color=plotting_color= null;
    }
    
    /**
     * Instantiation width default size. 
     * 
     * @param D a discrete distribution what we want to plot
     */
    public DiscreteDistributionPlot(DiscreteDistribution D)
    {
        this(D, DISTRIBUTION_PLOT_WIDTH, DISTRIBUTION_PLOT_HEIGHT);
    }
    
    private DiscreteDistribution distribution;
    
    private int range_max=DISTRIBUTION_PLOT_RANGE;
    private int range_min=0;

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
     * Color used for plotting
     */
    private Color plotting_color;
    
    
    
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
     * default, this is null, so no legend is painted.
     * 
     * @param c a color (null means default behavior)
     */
    public void setLegendColor(Color c)
    {
        this.legend_color = c;
    }

    /**
     * Sets the plotting color. By 
     * default, this is the same as the legend color   
     * 
     * @param c
     */
    public void setPlottingColor(Color c)
    {
        this.plotting_color = c;
    }
    /**
     * Sets the distribution that is to be plotted here.
     * 
     * @param D a discrete distribution (Poisson, NegativeBinomial or PointDistribution)
     */
    public final void setDistribution(DiscreteDistribution D){distribution = D;}
    
    /**
     * Sets the range for the plot: probabilities for n=range_min,range_min+1,...,range_max are plotted.
     * 
     * @param range_min minimum value (inclusive) for which probability is plotted.
     * @param range_max maximum value (inclusive) for which probability is plotted.
     */
    public void setRange(int range_min, int range_max){this.range_min = range_min; this.range_max=range_max;}


    /**
     * Plots the distribution at a given location. 
     * 
     * @param ignored this parameter is ignored, may be even null (required parameter for Icon interface)
     * @param g graphics context (should be convertible to Graphics2D)
     * @param x coordinate for the left-hand side
     * @param y coordinate for the top
     */
    @Override
    public void paintIcon(Component ignored,
               Graphics g,
               int x,
               int y)
    {
        double[] prob = distribution.getPointMassFunction(range_max);
        double max_prob = 0.0; 
        for (int i=range_min; i<=range_max; i++) if (prob[i]>max_prob) max_prob = prob[i];
        
        max_prob = 1.0; // scale is always the same
        
        double scale_y = (max_prob!=0.?(height-1.0)/max_prob:height);
        double scale_x = (width-1.0) / (range_max-range_min+1.25);
        
        Graphics2D g2 = (Graphics2D)g.create(); 
        Color original_color = g2.getColor();
        Stroke original_stroke=g2.getStroke();
        // draw background
        if (background_color != null)
        {
            g2.setColor(background_color);
            g2.fillRect(x, y, width, height);
            g2.setColor(original_color);
        }
        // draw axes
        if (legend_color != null)
        {
            
            g2.setColor(legend_color);
            g2.setStroke(new BasicStroke(0.5f)); // thin
            g2.drawLine(x, y, x, y+height);
            g2.drawLine(x, y+height, x+width, y+height);
            for (int i=1; i<4; i++)
            {
                double p=0.25*i; 
                int py = (int)(scale_y*p);
                g2.drawLine(x, y+height-py, x+2, y+height-py);
            }
        }
        
            
        
        // plot the distribution
        if (plotting_color != null)
            g2.setColor(plotting_color);

        g2.setStroke(new BasicStroke(2.0f)); // thick
        for (int i=range_min; i<=range_max; i++)
        {
            int px = (int)(scale_x*(i-range_min));
            int py = (int)(scale_y*prob[i]);
            if (py>0)
                g2.drawLine(x+1+px, y+height, x+1+px, y+height-py);
        }
        
        if (legend_color != null)
        {
            g2.setColor(legend_color);
            g2.setStroke(new BasicStroke(0.5f));            
            int font_size = g.getFont().getSize();
            String distribution_name = distribution.getClass().getSimpleName();
            int num_parameters = distribution.getNumParameters();
            g2.drawString(distribution_name, x+6, y+height-(font_size+1+num_parameters*font_size*12/10));
            String[] parameter_name = new String[num_parameters];
            if (distribution instanceof Poisson)
            {
                parameter_name[0] = "r"; 
            } else if (distribution instanceof NegativeBinomial)
            {
                parameter_name[0] = "κ"; // kappa
                parameter_name[1] = "q";
            } else if (distribution instanceof PointDistribution)
            {
                parameter_name[0] = "p(0)";
            } else if (distribution instanceof ShiftedGeometric)
            {
                parameter_name[0]="p0";
                parameter_name[1]="q";

            } else
            {
                for (int i=0; i<num_parameters; i++) {parameter_name[i]="θ"+Integer.toString(i+1);} // theta
            }
        
            double[] parameter_value = distribution.getParameters();
            for (int par_idx=num_parameters-1; par_idx>=0; par_idx--)
            {
                int par_y = height-(font_size+1+par_idx*font_size*12/10);
                g2.drawString(parameter_name[par_idx]+"="+Double.toString(parameter_value[par_idx]), x+10, y+par_y);
            }
        }
    }
 }
