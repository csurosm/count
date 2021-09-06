package ca.umontreal.iro.evolution.malin.ui;

/**
 * Table cell rendering for numerical values (table cell value of type Number) with variable background.
 * 
 * The cell content is expected to fall within a specified range: 
 * background color is set on a linear scale. 
 * A special color may be specified for 0. If this renderer is used with non-numerical 
 * cell content, the background is set to RED (would generate too many Exceptions otherwise).
 * 
 * @author csuros
 */
import java.awt.Color;
import java.awt.Component;
import java.awt.Font;

import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.SwingConstants;

import javax.swing.table.DefaultTableCellRenderer;

public class ColoredValueRenderer extends DefaultTableCellRenderer
{
    
    /**
     * Instantiates a new colored renderer with zero color equal to min color, and min value equal to 0.0. 
     * 
     * @param zero_color color corresponding to the minimum value of 0.0
     * @param max_color color corresponding to the maximum value
     * @param max_value numerical value of the range maximum
     */
    public ColoredValueRenderer(Color zero_color, Color max_color, double max_value)
    {
        this(zero_color, zero_color, max_color, 0.0, max_value);
    }

    /**
     * Instantiates a new colored renderer with zero color equal to min color. 
     * 
     * @param min_color color corresponding to the minimum value
     * @param max_color color corresponding to the maximum value
     * @param min_value numerical value of the range minimum
     * @param max_value numerical value of the range maximum
     */
    public ColoredValueRenderer(Color min_color, Color max_color, double min_value, double max_value)
    {
        this(min_color, min_color, max_color, min_value, max_value);
    }

    /**
     * Instantiates a new colored renderer.
     * 
     * @param zero_color color corresponding to a zero value
     * @param min_color color corresponding to the minimum value
     * @param max_color color corresponding to the maximum value
     * @param min_value numerical value of the range minimum
     * @param max_value numerical value of the range maximum
     */
    public ColoredValueRenderer(Color zero_color, Color min_color, Color max_color, double min_value, double max_value)
    {
        super();
        this.Cmin = min_color;
        this.Czero = zero_color;
        this.Cmax = max_color;
        this.max_value = max_value;
        this.min_value = min_value;
        this.setHorizontalTextPosition(SwingConstants.RIGHT);
        this.setHorizontalAlignment(SwingConstants.RIGHT);
    }
    
    
    private Color Cmin;
    private Color Cmax;
    private Color Czero;
    private double max_value;
    private double min_value;

    @Override
    public Component getTableCellRendererComponent(JTable table,
                                           Object value,
                                           boolean isSelected,
                                           boolean hasFocus,
                                           int row,
                                           int column)  
    {
        Component retval = super.getTableCellRendererComponent(table, value, false, hasFocus, row, column);
        if (row != -1)// && !isSelected)
        {
            if (value instanceof Number)
            {
                double x = ((Number)value).doubleValue();
                Color Cx = Czero;
                if (x!=0.0)
                {
                    double mid = Math.min((x-min_value)/(max_value-min_value),1.0);
                    Cx = intermediateColor(Cmin,Cmax,mid);
                }
                retval.setBackground(Cx);
                if (isSelected)
                {
                    JLabel label = (JLabel)retval;
                    label.setBorder(BorderFactory.createLineBorder(Cmax, 1));
                }
            } else
                retval.setBackground(Color.RED); // error!
        }
        return retval;
    }
    
    /**
     * Computes an intermediate color between two extreme values
     *
     * @param col0 color for 0
     * @param col1 color for 1
     * @param med a value between 0 and 1 (if outside that range, then the return value will be col0 or col1
     * @return a color that is the linear combination of col0 and col1, but it is enforced to fall between the two colors (fixing rounding problems).
     */
    public static Color intermediateColor(Color col0, Color col1, double med)
    {
        float[] rgb0 = col0.getRGBColorComponents(null);
        float[] rgb1 = col1.getRGBColorComponents(null);
        float[] rgb = new float[3];
        for (int j=0; j<3; j++)
        {
            rgb[j] = (float)(med*rgb1[j]+(1.-med)*rgb0[j]);
            float rmin = rgb0[j];
            float rmax = rgb1[j];
            if (rmin>rmax)
            {
                rmin = rmax;
                rmax = rgb0[j];
            }
            if (rgb[j]<rmin)
                rgb[j]=rmin;
            if (rgb[j]>rmax)
                rgb[j]=rmax;
        }
        return new Color(rgb[0], rgb[1], rgb[2]);
    }
    
}    
    
