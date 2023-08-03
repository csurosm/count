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


import javax.swing.InputVerifier;
import javax.swing.JComponent;
import javax.swing.JFormattedTextField;
import javax.swing.JOptionPane;

//import java.awt.Toolkit;

import java.text.ParseException;

/**
 * Some useful inpute verifiers for JFormattedTextField instances.
 * @author csuros
 */
public class InputVerifiers 
{

    private static PositiveInputVerifier positive_verifier;

    private static NonnegativeInputVerifier nonnegative_verifier;
    
    /**
     * Input must be a positive number.
     * 
     * @return
     */
    public static PositiveInputVerifier getPositiveInstance()
    {
        if (positive_verifier == null)
            positive_verifier = new PositiveInputVerifier();
        return positive_verifier;
    }

    /**
     * Input must be a non-negative number.
     * @return
     */
    public static NonnegativeInputVerifier getNonnegativeInstance()
    {
        if (nonnegative_verifier == null)
            nonnegative_verifier = new NonnegativeInputVerifier();
        return nonnegative_verifier;
    }


    /**
     * An input verifier for value within a given range.
     *
     */
    public static class RangeVerifier extends InputVerifier
    {
        private double min;
        private double max;

        /**
         * Instantiation with inclusive ranges
         *
         * @param min minimum value
         * @param max maximum value
         */
        public RangeVerifier(double min, double max)
        {
            this.min=min;
            this.max=max;
        }

        /**
         * The current value passes this test if it is not null,
         * it is a correctly formatted number, and falls within the range.
         *
         * @param C a
         * @return
         */
        @Override
        public boolean verify(JComponent C)
        {
            assert (C instanceof JFormattedTextField);
            //System.out.println("#**OUI.RV.v ["+min+","+max+"]\t"+C);

            try 
            {
                return verify((JFormattedTextField)C);
            } catch (ParseException PE)
            {
                return false;
            }
        }

        public boolean verify(JFormattedTextField F) throws ParseException
        {
            Object V = F.getFormatter().stringToValue(F.getText());
            if (V==null)
            {
                return false;
            } else
            {
                Number N = (Number)V;
                //System.out.println("#**OUI.RV.v ["+min+","+max+"]\t"+N+"\t"+C+"\t"+N.doubleValue());
                return valueWithinRange(N.doubleValue());
            }
        }
        
        @Override
        public boolean shouldYieldFocus(JComponent input)
        {
            return shouldYieldFocus((JFormattedTextField)input);
        }

        public boolean shouldYieldFocus(JFormattedTextField F)
        {
            String msg = null;
            boolean retval = false;
            try
            {
                Number N = (Number) F.getFormatter().stringToValue(F.getText());
                if (N==null)
                {
                    F.setFocusLostBehavior(JFormattedTextField.REVERT);
                    retval = true;
                } else
                {
                    msg = valueWithinRange(N);
                    retval  = (msg == null);
                }
            } catch (ParseException PE)
            {
                msg = "Badly formed value.";
                retval = false;
            }
            if (retval)
                F.setFocusLostBehavior(JFormattedTextField.COMMIT);
            else
            {
                msg = msg+"\nWould you like to keep the entered value and continue editing?";
                int edit_or_revert = JOptionPane.showConfirmDialog(F, msg , "Invalid entry", JOptionPane.YES_NO_OPTION, JOptionPane.ERROR_MESSAGE);
                if (edit_or_revert!=JOptionPane.YES_OPTION)
                {
                    F.setFocusLostBehavior(JFormattedTextField.REVERT);
                    retval=true;
                }
            }
            return retval;
        }

        protected String valueWithinRange(Number N)
        {
            boolean value_ok = valueWithinRange(N.doubleValue());
            if (!value_ok)
                return "Value outside of allowed range ("+min+".."+max+").";
            else
                return null;
        }
        protected boolean valueWithinRange(double x)
        {
            return (x>=min && x<=max);
        }

    }

    public static class NonnegativeInputVerifier extends RangeVerifier
    {
        public NonnegativeInputVerifier(){ super(0.0,Double.POSITIVE_INFINITY);}
        @Override
        protected String valueWithinRange(Number N)
        {
            double x = N.doubleValue();
            if (x<0.)
                return "Negative values are not allowed";
            else
                return null;
        }
    }
    
    public static class PositiveInputVerifier extends RangeVerifier
    {
        public PositiveInputVerifier(){super(0.0,Double.POSITIVE_INFINITY);}
        @Override
        protected boolean valueWithinRange(double x)
        {
            return (x>0.0);
        }
        @Override
        protected String valueWithinRange(Number N)
        {
            double x = N.doubleValue();
            if (x<0.)
                return "Negative values are not allowed";
            else if (x==0.0)
                return "Zero value is not allowed.";
            else
                return null;
        }
    }
}
