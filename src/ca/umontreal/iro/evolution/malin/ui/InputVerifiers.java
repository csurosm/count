
package ca.umontreal.iro.evolution.malin.ui;

import javax.swing.InputVerifier;
import javax.swing.JComponent;
import javax.swing.JFormattedTextField;
import javax.swing.JOptionPane;

//import java.awt.Toolkit;

import java.text.ParseException;

/**
 *
 * @author csuros
 */
public class InputVerifiers 
{
    private static final String MSG_SYNTAX = "Badly formed value";

    private static PositiveInputVerifier positive_verifier;

    private static NonnegativeInputVerifier nonnegative_verifier;
    
    public static PositiveInputVerifier getPositiveInstance()
    {
        if (positive_verifier == null)
            positive_verifier = new PositiveInputVerifier();
        return positive_verifier;
    }

    public static NonnegativeInputVerifier getNonnegativeInstance()
    {
        if (nonnegative_verifier == null)
            nonnegative_verifier = new NonnegativeInputVerifier();
        return nonnegative_verifier;
    }


        
//    public static class PositiveInputVerifier extends InputVerifier
//    {
//        @Override
//        public boolean verify(JComponent C)
//        {
//            try
//            {
//                JFormattedTextField F = (JFormattedTextField)C;
//                Object V = F.getFormatter().stringToValue(F.getText());
//                if (V==null)
//                {
//                    return false;
//                } else
//                {
//                    Number N = (Number)V;
//                    //System.out.println("#**OUI.PIV.v "+N+"\t"+C+"\t"+N.doubleValue());
//                    return N.doubleValue()>0;
//                }
//            } catch (ParseException PE)
//            {
//                return false;
//            }
//        }
//
//        @Override
//        public boolean shouldYieldFocus(JComponent input)
//        {
//            String msg = null;
//            boolean retval = false;
//            try
//            {
//                JFormattedTextField F = (JFormattedTextField)C;
//                Object V = F.getFormatter().stringToValue(F.getText());
//                if (V==null)
//                {
//                    F.setFocusLostBehavior(JFormattedTextField.REVERT);
//                    retval=true;
//                } else
//                {
//                    Number N = (Number)V;
//                    //System.out.println("#**OUI.PIV.v "+N+"\t"+C+"\t"+N.doubleValue());
//                    double v = N.doubleValue();
//                    if (v<0.)
//                    {
//                        msg = "Negative value";
//                    } else if (v==0.0)
//                    {
//                        msg = "Zero value";
//                    }
//                }
//            } catch (ParseException PE)
//            {
//                msg=MSG_SYNTAX;
//            }
//
//
////            boolean retval = verify(input);
////            if (!retval)
////                Toolkit.getDefaultToolkit().beep();
//            return retval;
//        }
//
//    }

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
                msg = msg+"\nWould you like to continue editing?";
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
