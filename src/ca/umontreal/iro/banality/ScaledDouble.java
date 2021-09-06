/*
 * ScaledDouble.java
 *
 * Created on January 23, 2005, 5:32 PM
 */

package ca.umontreal.iro.banality;

/**
 *
 * @author  csuros
 */
public class ScaledDouble implements Cloneable {
    
    /** Creates a new instance of ScaledDouble */
    public ScaledDouble() {
        this(0.0);
    }
    
    public ScaledDouble(double d){
        this(d, 0);
    }
    
    public ScaledDouble(double mantissa, int exponent){
        this.mantissa=mantissa;
        this.exponent=exponent;
        scale();
    }
    
    private double mantissa;
    private int exponent;
    
    private static double[] EXP10;
    private static final double LOG10=Math.log(10.0);
    
    private static final int EXP_SHIFT=300; // 307 is the maximum that can be handled with double
    static {
        int size=2*EXP_SHIFT+1;
        EXP10=new double[size];
        EXP10[EXP_SHIFT]=1.0;
        for (int i=EXP_SHIFT+1; i<size; i++){
            int e=i-EXP_SHIFT;
            if ((e & 1) == 0){
                int j=EXP_SHIFT+e/2;
                EXP10[i]=EXP10[j]*EXP10[j];
            } else 
                EXP10[i]=10.0*EXP10[i-1];
        }
        for (int i=EXP_SHIFT-1; i>=0; i--){
            int e=EXP_SHIFT-i;
            if ((e & 1) == 0){
                int j=EXP_SHIFT-e/2;
                EXP10[i]=EXP10[j]*EXP10[j];
            } else 
                EXP10[i]=EXP10[i+1]/10.0;
        }
        /*
        for (int j=0; j<size; j++)
            System.out.println("#**SD.static 10^"+(j-EXP_SHIFT)+"="+EXP10[j]);
         */
    }
    
    /**
     * Adds x to the current value.
     */
    public void addTo(ScaledDouble x){
        if (exponent==x.exponent){
            mantissa += x.mantissa;
        } else {
            int diff=exponent-x.exponent;
            if (diff<=EXP_SHIFT && diff >=-EXP_SHIFT){
                double sh=EXP10[EXP_SHIFT-diff];
                //System.out.println("#**SD.aT "+this+" + "+x+" diff "+diff+" shift "+sh);
                mantissa += x.mantissa * sh;
            } else {
                if (diff>EXP_SHIFT){
                    mantissa = x.mantissa;
                    exponent = x.exponent;
                }
            }
        }
    }
    
    public void addTo(double x){
        addTo(new ScaledDouble(x));
    }
    
    public Object clone(){
        ScaledDouble x=new ScaledDouble();
        x.mantissa=mantissa;
        x.exponent=exponent;
        return x;
    }
    
    public ScaledDouble add(ScaledDouble x){
        ScaledDouble b=(ScaledDouble)clone();
        b.addTo(x);
        return b;
    }

    public ScaledDouble add(double x){
        ScaledDouble b=(ScaledDouble)clone();
        b.addTo(x);
        return b;
    }
    /**
     * Multiplies the current value with x.
     */
    public void multiplyWith(ScaledDouble x){
        mantissa*=x.mantissa;
        exponent+=x.exponent;
        scale();
    }
    
    public void multiplyWith(double x){
        mantissa *= x;
        scale();
    }
    
    public ScaledDouble multiply(ScaledDouble x){
        ScaledDouble b=(ScaledDouble)clone();
        b.multiplyWith(x);
        return b;
    }
    
    public ScaledDouble multiply(double x){
        ScaledDouble b=(ScaledDouble)clone();
        b.multiplyWith(x);
        return b;
    }
    
    public void divideWith(double x){
        mantissa /= x;
        scale();
    }
    
    public void divideWith(ScaledDouble x){
        mantissa /= x.mantissa;
        exponent-=x.exponent;
        scale();
    }
    
    public ScaledDouble divide(double x){
        ScaledDouble b=(ScaledDouble)clone();
        b.divideWith(x);
        return b;
    }
    
    public ScaledDouble divide(ScaledDouble x){
        ScaledDouble b=(ScaledDouble)clone();
        b.divideWith(x);
        return b;
    }

    private void scale(){
        if (mantissa==0.0){
            exponent=0;
            return;
        }
        if (mantissa<1.0 || mantissa>10.0){
            //System.out.print("#**SD.s was "+mantissa+" E "+exponent);
            int e=(int)(Math.log(mantissa)/LOG10)-1;
            exponent+=e;
            mantissa *= EXP10[EXP_SHIFT-e];
            //System.out.println("; now "+mantissa+" E "+exponent);
        }
    }
        
    
    public double doubleValue(){
        if (exponent>EXP_SHIFT){
            if (mantissa<0)
                return Double.NEGATIVE_INFINITY;
            else 
                return Double.POSITIVE_INFINITY;
        } else if (exponent<-EXP_SHIFT)
            return 0.0;
        return   mantissa*EXP10[EXP_SHIFT+exponent];
    }
    
    public String toString(){
        scale();
        StringBuffer sb=new StringBuffer();
        if (exponent<4 && exponent >-3){
            sb.append(doubleValue());
        } else {
            sb.append(mantissa);
            sb.append("e");
            if (exponent>0)
                sb.append("+");
            sb.append(exponent);
        }
        return sb.toString();
    }
    
    public static void main(String[] args){
        ScaledDouble v=new ScaledDouble(1.0);
        System.out.println(v);
        v.multiplyWith(2e-100);
        System.out.println(v);
        v.addTo(3e-100);
        System.out.println(v);
        v.multiplyWith(1e-100);
        v.multiplyWith(1e-100);
        v.multiplyWith(1e-100);
        System.out.println(v);
        v.multiplyWith(2e200);
        System.out.println(v);
        v.divideWith(new ScaledDouble(2.0,-200));
        System.out.println(v);
    }
}
