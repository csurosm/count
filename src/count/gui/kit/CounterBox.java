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
import java.awt.Font;
import java.awt.event.MouseListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.SwingConstants;

/**
 * A component with a JLabel for displaying a count, 
 * and another for maximum (if defined). As a PropertyChangeListener, 
 * the new vales update the counter value. 
 * 
 * @since Wed 29 Mar 2012
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 *
 */
public class CounterBox extends Box implements PropertyChangeListener
{
	/**
	 * 
	 * 
	 * @param name name of this counter (shown), or null if no name 
	 * @param max_value info on value we are counting up to (shown), or null if indeterminate
	 */
    public CounterBox(String name, String max_value)
    {
        super(BoxLayout.LINE_AXIS);
        this.counter_name = name;
        this.max_value = (max_value==null?null:"/"+max_value);
        initComponents();
    }

    private String counter_name;
    private String max_value;

    private JLabel counter_label;
    private JLabel counter_value;

    private void initComponents()
    {
        if (counter_name != null)
        {
            counter_label = new JLabel(counter_name);
            add(counter_label);
        }

        counter_value = new JLabel("      "); //JTextField("");
        counter_value.setHorizontalTextPosition(SwingConstants.RIGHT);

        add(counter_value);
    }

    public void setMaximumInfo(String max)
    {
        this.max_value = max;
    }

    @Override
    public void setFont(Font f)
    {
        super.setFont(f);
        if (counter_label != null)
            counter_label.setFont(f);
        counter_value.setFont(f);
    }

    public void setTrackersOpaque(boolean b)
    {
        if (counter_label != null)
            counter_label.setOpaque(b);
        counter_value.setOpaque(b);
    }

    @Override
    public void setBackground(Color c)
    {
        super.setBackground(c);
        if (counter_label != null)
            counter_label.setBackground(c);
        counter_value.setBackground(c);
    }
    
    void setValue(String value)
    {
        if (max_value != null)
            value = value +max_value;
        counter_value.setText(value);
        repaint();
    }
    
    @Override
    public void setEnabled(boolean b)
    {
        super.setEnabled(b);
        counter_value.setEnabled(b);
    }

    /**
     * As a PropertyChangeListeners, the new values are used to set the counter's value
     */
    @Override
    public void propertyChange(PropertyChangeEvent pce) 
    {
        Object v = pce.getNewValue();
        setValue(v.toString());
    }
    
    @Override
    public String toString()
    {
    	StringBuilder sb = new StringBuilder();
    	sb.append(getClass().getSimpleName())
    	.append("[")
    	.append(counter_name)
    	.append(",")
    	.append(max_value)
    	.append("; label ")
    	.append(counter_label)
    	.append("; value ")
    	.append(counter_value);
    	return sb.toString();
    }



    @Override
    public void addMouseListener(MouseListener listener)
    {
        counter_value.addMouseListener(listener);
        if (counter_label != null)
            counter_label.addMouseListener(listener);
    }
}
