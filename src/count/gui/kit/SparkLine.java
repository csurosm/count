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
import java.awt.Component;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.List;

import javax.swing.Icon;

/**
 * A sparkline for displaying the tail of a sequence of values (history). 
 * 
 * @author csuros
 *
 */
public class SparkLine implements Icon
{
	/**
	 * Instantiation with dimensions and history tail length.
	 * 
	 * @param width history display width (icon width also includes the spark)
	 * @param height
	 * @param history_length nonnegative; use 0 to include all history
	 */
    public SparkLine(int width, int height, int history_length)
    {
        this.plot_width = width+SPARK_DIAMETER;
        this.plot_height = height;
        this.history = new ArrayList<>();
        setHistoryLength(history_length);         // 0 means all of history
    }

    private final int plot_width;
    private final int plot_height;
    private int max_history;
    private String info;

    private boolean stopped_history=false;

    private final List<Double> history;

    private Scaler scale; // set by paintIcon
    
    private boolean want_line = true;
    private boolean want_dots = true;
    private boolean want_bars = false;
    private boolean want_legend=true;

    private Color spark_color=Color.RED;
    private static final int SPARK_DIAMETER = 8;
    private static final Color DRAW_COLOR = Color.DARK_GRAY;
    private static final Color DOT_COLOR = new java.awt.Color(128, 0, 0); 
    

    /**
     * Color for displaying the accent (the "spark") at the end of the sparkline.
     * 
     * @param c null for no accent 
     */
    public void setSparkColor(Color c)
    {
        this.spark_color = c;
    }

    public void setHistoryLength(int max_history)
    {
        this.max_history = max_history;
    }

    /**
     * Displayed info about this sparkline. 
     * @param txt null if nothing to show 
     */
    public void setText(String txt)
    {
        this.info = txt;
    }
    
    public void setDrawLine(boolean want_line)
    {
    	this.want_line = want_line;
    }
    
    public void setDrawDots(boolean want_dots)
    {
    	this.want_dots = want_dots;
    }
    
    public void setDrawBars(boolean want_bars)
    {
    	this.want_bars = want_bars;
    }
    
    public void setDrawLegend(boolean want_legend)
    {
    	this.want_legend = want_legend;
    }

    
    public int historySize()
    {
        return history.size();
    }

    /**
     * Further calls to {@link #setValue(int, double)} have no effect
     * after this.
     */
    public void stopRecordingHistory()
    {
        stopped_history=true;
    }
    
    /**
     * The most recent scale used by {@link #paintIcon(Component, Graphics, int, int)}.
     * 
     * @return
     */
    public Scaler getScale()
    {
    	return scale;
    }
    
    

    @Override
    public void paintIcon(Component C, Graphics g, int x, int y)
    {
        Graphics2D g2 = (Graphics2D) g.create();


        int num_steps = historySize();

        double recent_max = Double.NEGATIVE_INFINITY;
        double recent_min = Double.POSITIVE_INFINITY;

        int idx = num_steps-1;
        int cnt = 0;

        while (idx>0 && cnt < max_history)
        {
            Double D = history.get(idx);
            if (D!=null)
            {
                double d = D.doubleValue();
                if (d<recent_min)
                    recent_min = d;
                if (d>recent_max)
                    recent_max = d;

                ++cnt;
            }
            --idx;
        }

        boolean has_data = !(Double.isInfinite(recent_min) || Double.isInfinite(recent_max)
                || recent_min==recent_max);

        int last_x = Integer.MAX_VALUE;
        int last_y = Integer.MAX_VALUE;

        int hgt = getIconHeight()-SPARK_DIAMETER;
        g2.setColor(DRAW_COLOR);

        int y0 = y+hgt+SPARK_DIAMETER/2;
        int legend_x=x+1;


        if (has_data)
        {
        	scale = new Scaler(recent_min, recent_max);
//            double recent_mid = 0.5*(recent_max+recent_min);
//            double recent_half_range = 0.5*(recent_max-recent_min);

            int legend_width = SPARK_DIAMETER;

            double scale_x =
                    max_history == 0
                    ?(getIconWidth()-legend_width-SPARK_DIAMETER)/((double)(num_steps))
                    :(getIconWidth()-legend_width-SPARK_DIAMETER)/((double)(max_history));

//            double half_height = getIconHeight()*0.5-ACCENT_DIAMETER;


            if (want_legend)
            {
                g2.drawLine(legend_x, y0-scale.getProjection(recent_min, hgt),
                            legend_x, y0-scale.getProjection(recent_max, hgt));

                double u = scale.getScaleUnit();
                double s0 = scale.getMin();

                int tic_idx = 1+(int)(s0/u);
                while (tic_idx*u<recent_min) tic_idx++;
                double ty = tic_idx*u;

                {
                    Font f = g2.getFont();
                    g2.setFont(C.getFont().deriveFont(9f));
                    g2.drawString(scale.toScaledString(tic_idx), legend_x+2, y0-scale.getProjection(ty, hgt)+2);

                    while (ty <= recent_max)
                    {
                        int tic_y = y0-scale.getProjection(ty, hgt);
                        g2.drawLine(legend_x, tic_y, legend_x+1, tic_y);
//                        if (ty % 3==0)
//                            g2.drawString(scale.toScaledString(tic_idx), legend_x+2, y0-scale.getProjection(ty, hgt)+2);

                        ty += u;
                        tic_idx++;
                    }
                    g2.setFont(f);
                }
            } else
            {
                g2.drawLine(legend_x, SPARK_DIAMETER/2, legend_x, y0);
            }

            int idx0=idx++;

            while (idx<num_steps)
            {
                Double D = history.get(idx);
                if (D!=null)
                {
                    int current_x = x+legend_width/2+(int)(scale_x*(idx-idx0));
                    int current_y = y+SPARK_DIAMETER/2+
                            hgt-scale.getProjection(D, hgt);


                    if (want_line)
                    {
                        if (last_x != Integer.MAX_VALUE)
                            g2.drawLine(last_x, last_y, current_x, current_y);
                    }
                    if (want_bars)
                        g2.drawLine(current_x, current_y, current_x, y+getIconHeight());
                    if (want_dots)
                    {
                    	g2.setColor(DOT_COLOR);
                        g2.fillOval(current_x-1, current_y-1, 1, 1);
                        g2.setColor(DRAW_COLOR);
                    }

                    last_x = current_x;
                    last_y = current_y;
                }
                idx++;
            }
        } else
        {
            last_x = 2;
            last_y =getIconHeight()/3;
        }


        if (info != null)
        {
                g2.setColor(new Color(0.9f,0.9f,0.9f,0.7f));
                java.awt.geom.Rectangle2D info_bounds = C.getFontMetrics(C.getFont()).getStringBounds(info, 0, info.length(), g2);

                int txt_x = (int)(getIconWidth()/2+last_x)/2;
                int txt_y = (getIconHeight()/2+last_y)/2;

                // centering
                txt_x -= (int)(0.5*info_bounds.getWidth());
                txt_y += (int)(0.5*info_bounds.getHeight());

                // whitened background
                g2.fillRoundRect((int)(txt_x+info_bounds.getX()), (int)(txt_y+info_bounds.getY()), (int)info_bounds.getWidth(), (int)info_bounds.getHeight(), 2, 2);

                g2.setColor(Color.BLACK);

                g2.drawString(info, txt_x, txt_y); // small offset for text height
        }

        if (has_data)
        { // draw the spark
            if (spark_color != null)
            {
                g2.setColor(spark_color);
                g2.fillOval(last_x-SPARK_DIAMETER/2, last_y-SPARK_DIAMETER/2, SPARK_DIAMETER, SPARK_DIAMETER);
            }
        }
    }


    public int setValue(int step, double value)
    {
        int cnt = 0;
        if (!stopped_history)
        {
            int s = historySize();
            while (s<=step)
            {
                history.add(null);
                ++s;
                ++cnt;
            }
            history.set(step, value);
        }
        return cnt;
    }

    @Override
    public int getIconWidth()
    {
        return plot_width;
    }

    @Override
    public int getIconHeight()
    {
        return plot_height;
    }
    

    
    

}
