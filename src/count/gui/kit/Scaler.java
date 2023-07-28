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


/**
 * A class for aligning a min..max range with 
 * multiples of a power of 10, and to 
 * project ingtermediate values onto an 0..w interval.  
 * 
 * @author csuros
 *
 */
public class Scaler 
{

    public Scaler(double min_value, double max_value)
    {
        assert (!Double.isInfinite(max_value) && !Double.isInfinite(min_value));
        assert (max_value >= min_value);

        this.min_value = min_value;
        this.max_value = max_value;
        double range = (max_value==min_value)?min_value:max_value-min_value;
        
        setScaling(range);
//
//        assert (range > 0.0);
//        
//        unit_magnitude = getDecimalMagnitude(range/3.0);
//        double unit = getPower10(unit_magnitude);
//
//        double steps = range/unit;
//        coin = 1;
//
//        // 1..7 stays by units of 1
//        if (steps>7.0)
//        {
//            if (steps>15)
//                coin = 5;
//            else
//                coin = 2;
//        }
//
//        // align min..max with scale units 
//        double scale_unit = getScaleUnit();
//        this.min_value = Math.floor(min_value/scale_unit)*scale_unit;
//        this.max_value = Math.ceil(max_value/scale_unit)*scale_unit;  
    }

    private double min_value;
    private double max_value;

    private static final double[] POSITIVE_POWERS =
        { 1, 10.0, 100.0, 1000.0, 10000.0, 100000.0,
                1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
                1e12, 1e13, 1e14, 1e15, 1e16, 1e17,
                1e18, 1e19, 1e20, 1e21, 1e22, 1e23,
                1e24, 1e25, 1e26, 1e27, 1e28, 1e29,
                1e30, 1e31
        };

    private static final double[] NEGATIVE_POWERS =
        {1., 0.1, 0.01, 0.001, 0.0001, 0.00001,
            1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11,
            1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17,
            1e-18, 1e-19, 1e-20, 1e-21, 1e-22, 1e-23,
            1e-24, 1e-25, 1e-26, 1e-27, 1e-28, 1e-29,
            1e-30, 1e-31

        };

    private static final String[] POSITIVE_MAGNITUDES =
        {
            "", "k", "M", "G", "T", "P", "E", "Z", "Y"
        };

    private static final String[] NEGATIVE_MAGNITUDES =
        {
            "", "m", "\u03bc", "n", "p", "f", "a", "z", "y"
        };

    private static final double getPower10(int exp)
    {
        return
                (exp<0)
                ?(-exp>=NEGATIVE_POWERS.length?NEGATIVE_POWERS[-exp]:Math.pow(10.0, exp))
                :(exp>=POSITIVE_POWERS.length?POSITIVE_POWERS[exp]:Math.pow(10.0, exp));
    }

    private int unit_magnitude;
    private int coin;

    public static int getDecimalMagnitude(double x)
    {
        return (int)Math.floor(Math.log10(x));
    }

    public void setScaling(double range)
    {

        if (range==0.0)
            return;
        
        unit_magnitude = getDecimalMagnitude(range/3.0);
        double unit = getPower10(unit_magnitude);

        double steps = range/unit;
        coin = 1;

        // 1..7 stays by units of 1
        if (steps>7.0)
        {
            if (steps>15)
                coin = 5;
            else
                coin = 2;
        }

        double scale_unit = getScaleUnit();
        this.min_value = Math.floor(min_value/scale_unit)*scale_unit;
        this.max_value = Math.ceil(max_value/scale_unit)*scale_unit;
    }

    public double getScaleUnit()
    {
        return coin * getPower10(unit_magnitude);
    }

    public double getMin()
    {
        return min_value;
    }

    public double getMax()
    {
        return max_value;
    }

    /**
     * Projects a value x within the min..max range to an
     * interger value in a 0-based interval. 
     * 
     * @param x 
     * @param length interval length; 0=min et length=max
     * @return 0 for {@link #getMin()}, length for {@link #getMax()}; linear interpolation for other x
     */
    public int getProjection(double x, int length)
    {

        int z =(int)(length*(x-min_value)/(max_value-min_value));
        return z; //Math.min(Math.max(0, z),width);
    }

    /**
     *
     * @param q this many units (multiplied by {@link #getScaleUnit() })
     * @return
     */
    public String toScaledString(int q)
    {
        if (q==0) return "0";

        int u3 = (unit_magnitude+(q<0?-1:1))/3;
        String prefix = "";
        if (u3<0)
        {
            if (-u3<NEGATIVE_MAGNITUDES.length)
            {
                prefix = NEGATIVE_MAGNITUDES[-u3];
                u3 *= 3;
            }
            else
            {
                u3 = unit_magnitude;
                prefix = Integer.toString(unit_magnitude);
            }
        } else
        {
            if (u3<POSITIVE_MAGNITUDES.length)
            {
                prefix = POSITIVE_MAGNITUDES[u3];
                u3 *= 3;
            }
            else
            {
                u3 = unit_magnitude;
                prefix = "+"+Integer.toString(unit_magnitude);
            }
        }

        String sign = "";
        if (q<0) {q=-q; sign="-";}

        int qint = q*coin;
        int qfrac = 0;
        if (u3 != unit_magnitude)
        {
            if (u3==unit_magnitude+1)
            {
                qfrac = qint % 10;
                qint = qint / 10;
            } else
            {
                qint *= 10;
            }
        }
        String dv = Integer.toString(qint);
        if (qfrac != 0)
            dv = dv + "." + Integer.toString(qfrac);
        return sign+dv+prefix;
    }

}
