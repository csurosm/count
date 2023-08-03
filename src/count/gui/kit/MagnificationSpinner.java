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


import java.awt.Dimension;
import java.awt.Font;

import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

/**
 * A spinner with variable increments:
 * possible values are 1.25,1.5,...,10,12.5,..,100,125,... 
 * and 0.9,0.8,...0.1,0.09,0.08,...,0.01,...
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 */
public class MagnificationSpinner extends JSpinner
{
	/** 
	 * Default minimum value
	 */
    public static double MAGNIFICATION_MIN = 0.1;
    /**
     * Default maximum value
     */
    public static double MAGNIFICATION_MAX = 10.;

    public MagnificationSpinner(double initial_value, double min_value, double max_value)
	{
		super(new Model(initial_value, min_value, max_value));

        JSpinner.NumberEditor magnification_editor = new JSpinner.NumberEditor(this, "#%");
        magnification_editor.getTextField().setFont(new Font("SansSerif", Font.PLAIN, 10));
        magnification_editor.getTextField().setBorder(null);
        magnification_editor.getTextField().setEditable(false);
        magnification_editor.getTextField().setBackground(this.getBackground()); // setEditable() may have changed the background color

        magnification_editor.getTextField().setToolTipText("Current magnification value");
        this.setToolTipText("Change magnification");
        this.setEditor(magnification_editor);
        this.setMaximumSize(new Dimension(80,30));
	}
	/**
	 * Instantiation with default range and initiali value 1.0. 
	 */
	public MagnificationSpinner()
	{
		this(1.0, MAGNIFICATION_MIN, MAGNIFICATION_MAX);
	}
	
    
	/**
	 * Our own {@link SpinnerNumberModel}. 
	 * 
	 * @author csuros
	 *
	 */
	private static class Model extends SpinnerNumberModel
	{
		public Model(double value, double minimum, double maximum)
	    {
	        super(value, minimum, maximum, 1.0);
	    }
		
	        
	    @Override
	    public Object getPreviousValue()
	    {
	        double v = ((Number)getValue()).doubleValue();
	        double m = ((Number)getMinimum()).doubleValue();
	        if (m>=v)
	            return null;
	        double d = increment();
	        double w= d*((int)(v/d+0.5))-d;
//	        System.out.println("#*MS.M.gPV "+v+"\td "+d+"\tw "+w);
	        return Double.valueOf(w);
	    }

	    @Override
	    public Object getNextValue()
	    {
	        double v = ((Number)getValue()).doubleValue();
	        double m = ((Number)getMaximum()).doubleValue();
	        if (m<=v)
	            return null;
	        double d = increment();

	        double w= d*((int)(v/d+0.5))+d;
//	        System.out.println("#*MS.M.gNV "+v+"\td "+d+"\tw "+w);
	        
	        return Double.valueOf(w);
	    }

		private double increment()
		{
	        Number V = (Number)getValue();
	        double v = V.doubleValue();
	        
	        double d = 0.1;
	        if (v<1.0)
	        {
	            double z = v; // scale
	            while (z<0.1)
	            {
	                z*=10.0;
	                d*=0.1;
	            }
	        } else
	        {
	            double z=v; // scale
	            while(z>=1.0)
	            {
	                z*=0.1;
	                d*=10.0;
	            }
	        }
	        
	        
	        if (d>=1.0)
	            d/=4.0;
	        
	        return d;
			
		}	    
	}
	


}
