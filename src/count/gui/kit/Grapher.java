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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.beans.IndexedPropertyChangeEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.text.NumberFormat;

import javax.swing.JLabel;
import javax.swing.Timer;

/**
 * A component (JLabel) for tracking a stream of double values
 * in a {@link SparkLine}, 
 * delivered via {@link IndexedPropertyChangeEvent}s.   
 * 
 * 
 * @author csuros
 *
 */
public class Grapher extends JLabel implements PropertyChangeListener
{
    public Grapher(int width, int height)
    {
        super();
        this.history_plot = new SparkLine(width, height, 100);
        setIcon(history_plot);
        setOpaque(true);
    }
    
    private SparkLine history_plot;


    public void setInfo(String txt)
    {
        history_plot.setText(txt);
    }

    public void setHistoryLength(int length)
    {
        history_plot.setHistoryLength(length);
    }

    public void setDrawLine(boolean want_line)
    {
    	history_plot.setDrawLine(want_line);
    }
    
    public void setDrawDots(boolean want_dots)
    {
    	history_plot.setDrawDots(want_dots);
    }
    
    public void setDrawBars(boolean want_bars)
    {
    	history_plot.setDrawBars(want_bars);
    }
    
    public void setDrawLegend(boolean want_legend)
    {
    	history_plot.setDrawLegend(want_legend);
    }


    @Override
    public String getToolTipText(MouseEvent e)
    {
        String s = super.getToolTipText(e);
        //System.out.println("#*T.G.gTT "+s+"\t// "+e);


        Scaler scale = history_plot.getScale();
        if (scale==null)
        	return s;
        else 
        {        
        	if (s==null) s="";
        	else s=s+"; ";

        	return s+"one tic on the axis is "+scale.getScaleUnit();
        }
    }

    /**
     *
     * @param pce must be an {@link java.beans.IndexedPropertyChangeEvent} or else ignored
     */
    @Override
    public void propertyChange(PropertyChangeEvent pce)
    {
        if (pce instanceof IndexedPropertyChangeEvent) // otherwise ignored
        {
            int step = ((IndexedPropertyChangeEvent)pce).getIndex();
            double d = ((Number)pce.getNewValue()).doubleValue();
            int num_steps_advanced = history_plot.setValue(step,d);
            if (num_steps_advanced>0)
            {
                if (d==0.0)
                    history_plot.setText(pce.getPropertyName()+"=0.0");
                else
                    {
                        int dec = (int)Math.log10(Math.abs(d));

                        int f = (dec>1?1:-(dec-2));
                        NumberFormat NF = NumberFormat.getNumberInstance();
                        NF.setMaximumFractionDigits(f);
                        if (f>0)
                            NF.setMinimumFractionDigits(f);
                        history_plot.setText(pce.getPropertyName()+"="+NF.format(d));
                    }
            }
            repaint();
        }
    }

    public void stopHistory()
    {
        history_plot.stopRecordingHistory();
        repaint();
    }


    private Color glow_color;
    private double glow_phase;

    private Timer glow_timer;

    /**
     * Pulsating spark for the current (last) value.
     *
     * @param color The full color for the glow
     * @param frequency in Hertz
     */
    public void glow(Color color, double frequency)
    {
        this.glow_color = color;
        this.glow_phase = 0.0;

        class GlowPhaser implements ActionListener
        {
            GlowPhaser(double freq){ this.freq = freq;}

            private double freq;

            @Override
            public void actionPerformed(ActionEvent ae)
            {
                Timer T = (Timer) ae.getSource();

                glow_phase += freq * T.getDelay()*0.001* 2.0 * Math.PI;

                repaint();
                
            }
        }

        if (glow_timer != null)
            glow_timer.stop();

        if (frequency == 0.0)
        {
            glow_timer = null;
            history_plot.setSparkColor(glow_color);
            repaint();
        } else
        {
            glow_timer = new Timer(20, new GlowPhaser(frequency));
            glow_timer.start();
        }
    }



    @Override
    protected void paintComponent(Graphics gr)
    {
        if (glow_color != null)
        {
            double phi = 0.5*(1.0+Math.cos(glow_phase)); // b/w 0 and 1 

            float r = (float) (glow_color.getRed()/256.0*phi);
            float g = (float) (glow_color.getGreen()/256.0*phi);
            float b = (float) (glow_color.getBlue()/256.0*phi);

            history_plot.setSparkColor(new Color(r,g,b,glow_color.getAlpha()/256f));
        }
        super.paintComponent(gr);
    }
}
