package count.ds;
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
import java.util.HashMap;
import java.util.Map;




/**
 * Class for nonnegative integers + missing entries (any negative integer)
 * 
 * @author csuros
 */
public class IntegerOrMissing extends Number implements Comparable<IntegerOrMissing>
{
    private final int value;
    private static Map<Integer, IntegerOrMissing> cached_values;
    
    /**
     * Displayed string for missing entries.
     */
    public static String MISSING_ENTRY = "?";

    public IntegerOrMissing(int n)
    {
        this.value = n;
    }

    @Override
    public int compareTo(IntegerOrMissing z)
    {
        return Integer.compare(value, z.value);
    }

    @Override
    public double doubleValue()
    {
        return (double) value;
    }

    @Override
    public int intValue()
    {
        return (int) value;
    }

    @Override
    public float floatValue()
    {
        return (float) value;
    }

    @Override
    public long longValue()
    {
        return (long) value;
    }

    @Override
    public String toString()    
    {
        if (value<0)
            return MISSING_ENTRY; // sad face
        else
            return Integer.toString(value);
    }    
    
    public boolean isAmbiguous()
    {
        return value<0;
    }
    
    /**
     * Cached instances. 
     * 
     * @param value
     * @return
     */
    public static IntegerOrMissing getInteger(int value)
    {
    	if (cached_values==null)
    		cached_values=new HashMap<>();
    	if (cached_values.containsKey(value))
    		return cached_values.get(value);
    	else
    	{
    		IntegerOrMissing x = new IntegerOrMissing(value);
    		cached_values.put(value, x);
    		return x;
    	}
    }
    
    /**
     * An appropriate class for displaying these numbers in a Swing JTable: 
     * the rounded value appears in the cell but the tool tip text gives the true value. 
     */
    public static class Renderer extends javax.swing.table.DefaultTableCellRenderer
    {
        @Override
        public Component getTableCellRendererComponent(javax.swing.JTable table,
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
                if (value instanceof IntegerOrMissing) // it'd better be...
                {
                	IntegerOrMissing I = (IntegerOrMissing) value;
                	if (I.isAmbiguous())
                		setToolTipText("n/a");
                	else 
                		setToolTipText(I.toString());
                } else
                {
                    setToolTipText("[not and IntegerOrMissing] "+value);
                }
            }
            
            return comp;
        }
    }
    
}