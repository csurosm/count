
package ca.umontreal.iro.evolution.malin;

/**
 * A wrapper class four double values: toString() does rounding to about three digits.
 * 
 * @author csuros
 */

import java.awt.Component;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;

public class DoubleRoundedForDisplay extends Number implements Comparable<Number>
{
    public DoubleRoundedForDisplay(double value)
    {
        this.value = value;
    }
    private double value;

    /**
     * Comparison is based on the ful double values
     * 
     * @param z
     * @return
     */
    public int compareTo(Number z)
    {
        return Double.compare(value, z.doubleValue());
    }

    public double doubleValue(){ return value;}
    
    public float floatValue(){return (float)value;}
    
    public long longValue(){return (long)value;}
    
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
        if (value < 0.0005)
            return ".";
        else if (value == (int)value)
            return Integer.toString((int)value);
        else if (value<0.0095)
        { // rounded as 0.00x
            int z = (int)(value*1000.0+0.5);
            return "0.00"+Integer.toString(z);
        } else if (value<0.095)
        { // rounded at last valuable digit
            int z = (int)(value*100.0+0.5);
            return "0.0"+Integer.toString(z);
        } else if (value<0.95)
        { // rounded at two valuable digits
            int z = (int)(value*100.0+0.5);
            return "0."+Integer.toString(z);
        } else if (value<9.95)
        { // rounded at two valuable digits
            int z = (int)(value*10.0+0.5);
            return Integer.toString(z/10)+"."+Integer.toString(z%10);
        } else
        { // rounded as an integer
            int z = (int)(value+0.5);
            return Integer.toString(z);
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
                if (value instanceof DoubleRoundedForDisplay) // it'd better be...
                {
                    DoubleRoundedForDisplay D = (DoubleRoundedForDisplay) value;
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
    