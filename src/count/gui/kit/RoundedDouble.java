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

import java.awt.Component;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;


/**
 * A wrapper class four double values: toString() does rounding to about three digits.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 */

public class RoundedDouble extends Number implements Comparable<Number>
{
    public RoundedDouble(double value)
    {
        this.value = value;
    }
    private final double value;

    /**
     * Comparison is based on the full double values
     * 
     * @param z
     * @return
     */
    @Override
    public int compareTo(Number z)
    {
        return Double.compare(value, z.doubleValue());
    }
    

    @Override
    public double doubleValue(){ return value;}
    
    @Override
    public float floatValue(){return (float)value;}
    
    @Override
    public long longValue(){return (long)value;}
    
    @Override
    public int intValue(){return (int)value;}

    @Override
    public String toString()
    {
        return toString(value);
    }
            
    /**
     * Produces a nicely rounded value for the argument. 
     * 
     * @param value a double value
     * @return the String representation
     */
    public static String toString(double value)    
    {
        if (Double.isNaN(value))
            return "\u2639"; // sad face
        if (Double.isInfinite(value))
        {
            if (value == Double.POSITIVE_INFINITY)
                return "\u221e"; // infinity sign
            else 
                return "-\u221e";
        }
        String sign="";
        if (value<0)
        	sign = "-";
        
        value = Math.abs(value);
        
        if (value == 0) return " ";
        else if (value < 0.0005)
            return ".";
        else if (value == (int)value)
            return sign+Integer.toString((int)value);
        else if (value<0.0095)
        { // rounded as 0.00x
            int z = (int)(value*1000.0+0.5);
            return sign+"0.00"+Integer.toString(z);
        } else if (value<0.095)
        { // rounded at last valuable digit
            int z = (int)(value*100.0+0.5);
            return sign+"0.0"+Integer.toString(z);
        } else if (value<0.95)
        { // rounded at two valuable digits
            int z = (int)(value*100.0+0.5);
            return sign+"0."+Integer.toString(z);
        } else if (value<9.95)
        { // rounded at two valuable digits
            int z = (int)(value*10.0+0.5);
            return sign+Integer.toString(z/10)+"."+Integer.toString(z%10);
        } else
        { // rounded as an integer
            int z = (int)(value+0.5);
            return sign+Integer.toString(z);
        }
    }
    
    @Override
    public boolean equals(Object o)
    {
        boolean retval 
            =   (o instanceof Number)
                &&
                (((Number)o).doubleValue() == value);

        return retval;
    }
    
    @Override
    public int hashCode()
    {
        return Double.hashCode(value);
    }
    
    /**
     * An appropriate class for displaying these numbers in a Swing JTable: 
     * the rounded value appears in the cell but the tool tip text gives the true value. 
     */
    public static class Renderer extends DefaultTableCellRenderer
    {
        @Override
        public Component getTableCellRendererComponent(JTable table,
                                               Object value,
                                               boolean isSelected,
                                               boolean hasFocus,
                                               int row,
                                               int column)        
        {
            Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
            if (value == null)
            {
                setToolTipText("");
            } else
            {
                if (value instanceof RoundedDouble) // it'd better be...
                {
                    RoundedDouble D = (RoundedDouble) value;
                    setToolTipText(Double.toString(D.value));
                } else
                {
                    setToolTipText("[not a DoubleRoundedForDisplay] "+value);
                }
            }
            
            return comp;
        }
    }
}