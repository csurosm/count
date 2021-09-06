/*
 * ScaledDoubleArray.java
 *
 * Created on October 4, 2005, 9:45 PM
 */

package ca.umontreal.iro.banality;

/**
 * A class for handling vectors of (mantissa, exponent) 
 * pairs needed to avoid underflow in 
 * likelihood computations. Mantissa is stored as 
 * double (normalized between 1 and 10); 
 * exponent is stored as long.
 * It is expected that all numbers are non-negative: 
 * they are probabilities. 
 *
 * Memory usage for n elements is 48+n*16 bytes.
 * 
 * @author  csuros
 */
public class ScaledProbabilities implements Cloneable {
    
    private static final double MIN_MANTISSA=1.0; // inclusive
    private static final double MAX_MANTISSA=10.0; // exclusive
    
    
    /** Creates a new instance of ScaledProbabilities, 
     * for given number of elements. 
     * All elements are initialized to 0.
     *
     * @param size number of elements held in this array.
     */
    public ScaledProbabilities(int size) {
        mantissa=new double[size];
        exponent=new long[size];
    }
    
    /**
     * Creates a new ScaledProbabilities, with elements 
     * initialized to the given values.
     * 
     * @param values initial values of the elements.
     */
    public ScaledProbabilities(double[] values){
        this(values.length);
        System.arraycopy(values,0,mantissa,0,values.length);
        normalize();
    }
    
    private double[] mantissa;
    private long[] exponent;
    
    /**
     * Returns the number of elements in this array.
     */
    public int size(){return mantissa.length;}
    
    /**
     * Returns the mantissa of the element in the given position.
     */
    public double getMantissa(int pos){
        return mantissa[pos];
    }
    
    /**
     * Returns the exponent of the element in the given position.
     */
    public long getExponent(int pos){
        return exponent[pos];
    }
    
    /** 
     * Returns a double value for this guy
     */
    public double doubleValue(int pos){
        if (exponent[pos]<-OFFSET_POWER10)
            return 0.0;
        else if (exponent[pos]>OFFSET_POWER10)
            return Double.POSITIVE_INFINITY;
        else
            return mantissa[pos]*POWER10[(int)exponent[pos]+OFFSET_POWER10];
    }
    
    /** 
     * Returns the base-10 logarithm for the value 
     */
    public double lgValue(int pos){
        return exponent[pos]+log10(mantissa[pos]);
    }
    
        
    /**
     * Returns a clone of this scaled 
     * array (another ScaledProbabilities object),
     * with identical size and elements.
     */
    public Object clone(){
        ScaledProbabilities A=new ScaledProbabilities(mantissa.length);
        copyInto(A);
        return A;
    }
    
    /**
     * Copies the elements into another ScaledProbabilities object.
     *
     * @param A where the element are copied (must be of the same size as this guy)
     */
    public void copyInto(ScaledProbabilities A){
        int size=mantissa.length;
        System.arraycopy(mantissa,0,A.mantissa,0,size);
        System.arraycopy(exponent, 0, A.exponent, 0,size);
    }
    
    /**
     * Multiplies every element of this array with the elements of another one.
     *
     * @param factor another scaled prob array: elements will be set to e[i]*factor[i]
     */
    public void multiplyWith(ScaledProbabilities factor) {
        int size=mantissa.length;
        double[] fm=factor.mantissa;
        long[] fe=factor.exponent;
        for (int j=0; j<size; j++){
            mantissa[j]*=fm[j];
        }
        for (int j=0; j<size; j++){
            exponent[j]+=fe[j];
        }
        normalize();
    }
    
    /**
     * Element-wise multiplication of two arrays.
     *
     * @param factor1
     * @param factor2
     * @param result where the product is put: result[i]=factor1[i]*factor2[i]
     */
    public static void multiply(ScaledProbabilities factor1, ScaledProbabilities factor2, ScaledProbabilities result){
        {
            double[] a=factor1.mantissa;
            double[] b=factor2.mantissa;
            double[] c=result.mantissa;
            int n=factor1.size();
            for (int j=0; j<n; j++)
                c[j]=a[j]*b[j];
        }
        {
            long[] a=factor1.exponent;
            long[] b=factor2.exponent;
            long[] c=result.exponent;
            int n=factor1.size();
            for (int j=0; j<n; j++)
                c[j]=a[j]+b[j];
        }
        result.normalize();
    }
    
    private static final long DIGITS=15L;
    /**
     * Adds term[i] to every element emt[i].
     */
    public void add(ScaledProbabilities term){
        int size=mantissa.length;
        double[] tm=term.mantissa;
        long[] te=term.exponent;
        for (int j=0; j<size; j++){
            if (tm[j]==0.0) continue;
            if (mantissa[j]==0.0) {
                mantissa[j]=tm[j];
                exponent[j]=te[j];
                continue;
            }
            long diff = te[j]-exponent[j];
            if (diff>-DIGITS) {
                if (diff>DIGITS){
                    mantissa[j]=tm[j];
                    exponent[j]=te[j];
                } else { 
                    double x=tm[j]*POWER10[OFFSET_POWER10+(int)diff];
                    mantissa[j] += x;
                    normalize(j);
                }
            }
        }
    }

    
    /**
     * Subtracts term[i] from every element emt[i].
     */
    public void sub(ScaledProbabilities term){
        int size=mantissa.length;
        double[] tm=term.mantissa;
        long[] te=term.exponent;
        for (int j=0; j<size; j++){
            if (tm[j]==0.0) continue;
            if (mantissa[j]==0.0) {
                mantissa[j]=-tm[j];
                exponent[j]=te[j];
                continue;
            }
            long diff = te[j]-exponent[j];
            if (diff>-DIGITS) {
                if (diff>DIGITS){
                    mantissa[j]=tm[j];
                    exponent[j]=te[j];
                } else { 
                    double x=tm[j]*POWER10[OFFSET_POWER10+(int)diff];
                    mantissa[j] -= x;
                    normalize(j);
                }
            }
        }
    }    
    
    /**
     * Element-wise addition of two arrays.
     *
     * @param term1 first term
     * @param term2 second term
     * @param result where the sum is put: result[i]=term1[i]+term2[i]
     */
    public static void add(ScaledProbabilities term1, ScaledProbabilities term2, ScaledProbabilities result){
        double [] m1=term1.mantissa;
        double [] m2=term2.mantissa;
        double [] m=result.mantissa;
        
        long [] e1=term1.exponent;
        long [] e2=term2.exponent;
        long [] e=result.exponent;
        
        int n=term1.size();
        for (int j=0; j<n; j++){
            if (e1[j]==0.0){
                e[j] = e2[j];
                m[j] = m2[j];
                continue;
            }
            if (e2[j]==0.0){
                e[j] = e1[j];
                m[j] = m1[j];
                continue;
            }
            long diff = e2[j]-e1[j];
            if (diff>-DIGITS){
                m[j] = m1[j]+m2[j]*POWER10[OFFSET_POWER10+(int)diff];
                e[j] = e1[j];
                result.normalize(j);
            }
        }
    }
    

    /**
     * Sets the element at the given position.
     * @param pos position at which the element is set
     * @param M mantissa
     * @param E exponent (base 10)
     */
    public void set(int pos, double M, long E){
        mantissa[pos]=M;
        exponent[pos]=E;
        normalize(pos);
    }
    
    /**
     * Sets the element at the given position.
     * @param pos position at which the element is set
     * @param val new value
     */
    public void set(int pos, double val){
        set(pos,val,0L);
    }
    
    /**
     * Adds M*10^E to the element in the given position.
     * @param pos position at which the element is changed
     * @param M mantissa (if you are careful, this can even be negative, but if emt[pos] 
     *      becomes negative, weird things will happen)
     * @param E exponent (base 10)
     */
    public void add(int pos, double M, long E){
        if (M==0.)
            return;
        if (mantissa[pos]==0.0){
            mantissa[pos]=M;
            exponent[pos]=E;
        } else {

            long parE = E+(long)log10(M);

            //System.out.println("add-000 "+mantissa[pos]+"E"+exponent[pos]+"+= "+M+"E"+E+"; parE "+parE);

            if (parE<exponent[pos]-DIGITS)
                return;
            if (parE>exponent[pos]+DIGITS){
                mantissa[pos]=M;
                exponent[pos]=E;
            } else {
                long diff=E-exponent[pos];
                //System.out.println("#**SP.a "+pos+" "+mantissa[pos]+"E"+exponent[pos]+"+="+M+"E"+E+"; diff "+diff);
                double x = M*POWER10[OFFSET_POWER10+(int)diff];
                mantissa[pos]+=x;
            }
        }

        normalize(pos);
    }
    
    public void add(int pos, double value){
        add(pos, value, 0L);
    }
        
    /**
     * Multiplies with M*10^E to the element in the given position.
     * @param pos position at which the element is changed
     * @param M mantissa
     * @param E exponent (base 10)
     */
    public void multiply(int pos, double M, long E){
        mantissa[pos]*=M;
        exponent[pos]+=E;
        
        normalize(pos);
        
    }
    
    public void multiply(int pos, double val){
       multiply(pos, val, 0L); 
    }

    private static final double LOG10 = Math.log(10.0);
    private static final double LOG10_1 = 1.0/LOG10;

    /**
     * Computes the base-10 logarithm of a double value.
     */
    public static final double log10(double v){
        return Math.log(v)*LOG10_1;
    }

    /**
     * powers of 10 precomputed
     */
    private static double[] POWER10;
    
    /**
     * max and min of exponents for double
     */
    private static final int OFFSET_POWER10=307;
    
    static {
        // fill up POWER10[] array
        POWER10=new double[2*OFFSET_POWER10+1];
        POWER10[OFFSET_POWER10]=1.0;
        for (int j=-1; j+OFFSET_POWER10>=0; j--){
            if (j % 2 == 0){
                int idx_sqrt = OFFSET_POWER10+j/2;
                POWER10[OFFSET_POWER10+j]=POWER10[idx_sqrt]*POWER10[idx_sqrt];
            } else {
                POWER10[OFFSET_POWER10+j]=POWER10[OFFSET_POWER10+j+1]*0.1;
            }
        }
        for (int j=1; j<=OFFSET_POWER10; j++){
            if (j % 2 == 0){
                int idx_sqrt = OFFSET_POWER10+j/2;
                POWER10[OFFSET_POWER10+j]=POWER10[idx_sqrt]*POWER10[idx_sqrt];
            } else {
                POWER10[OFFSET_POWER10+j]=POWER10[OFFSET_POWER10+j-1]*10.0;
            }
        }
    }
    
//    /**
//     * Shifts all elements by one position down. Element 1 goes 
//     * to position 0, element 2 goes to position 1, etc. 
//     * Element 0 is lost, and the last element is set to 0.0.
//     */
//    public void shiftDown(){
//        int size=mantissa.length;
//        for (int j=0; j<size-1; j++){
//            mantissa[j]=mantissa[j+1];
//            exponent[j]=exponent[j+1];
//        }
//        mantissa[size-1]=0.;
//        exponent[size-1]=0L;
//    }
//    
//    /**
//     * Shifts all elements by one position up. Element 0 goes 
//     * to position 1, element 1 goes to position 2, etc. 
//     * Element 0 is set to 0.0, and the last element is lost.
//     */
//    public void shiftUp(){
//        int size=mantissa.length;
//        for (int j=size-1; j>0; j--){
//            mantissa[j]=mantissa[j-1];
//            exponent[j]=exponent[j-1];
//        }
//        mantissa[0]=0.;
//        exponent[0]=0L;
//    }

    /**
     * Computes 10^e.
     */
    public static final double power10(int e){
        return POWER10[OFFSET_POWER10+e];
    }
    
    /**
     * String representation of the floating-point value at a given index.
     */
    public String toString(int pos){
        StringBuffer sb = new StringBuffer();
        sb.append(mantissa[pos]);
        sb.append('E');
        sb.append(exponent[pos]);
        return sb.toString();
    }

    public String toString(){
        int max_idx=mantissa.length;
        if (max_idx>10) max_idx=10;
        StringBuffer sb = new StringBuffer();
        sb.append("ScaledProbs[");
        for (int j=0; j<max_idx; j++){
            if (j>0) sb.append(", ");
            sb.append(toString(j));
        }
        sb.append("]");
        return sb.toString();
    }
    
    /**
     * Scales the values so that all mantissas are between 1 an 10.
     */
    private void normalize(){
        int size=mantissa.length;
        for (int j=0; j<size; j++)
            normalize(j);
    }
    
    /**
     * Scales the element at this position so that its mantissa is between 1 and 10.
     */
    private void normalize(int pos){
        if (mantissa[pos]==0.0){
            exponent[pos]=0L;
            return;
        }
            
        if (mantissa[pos]<MIN_MANTISSA || mantissa[pos]>=MAX_MANTISSA){
            double e=log10(mantissa[pos]);
            int shift = (int)e; if (e<0.0) shift--;
            double scale = POWER10[OFFSET_POWER10-shift];
            if (scale==0.0){
                System.out.println("$$$ SP scale "+scale+" shift "+shift+" M "+mantissa[pos]+" E "+exponent[pos]+" e "+e);
            }
            mantissa[pos]*=scale;
            exponent[pos]+=shift;
            
        }
    }

    public static void main(String[] args){
        int s=10;
        double[] x=new double[s];
        for (int j=0; j<s; j++)
            x[j]=Math.pow(2.0,j-s/2);
        ScaledProbabilities A=new ScaledProbabilities(x);
        for (int j=0;j<s; j++){
            System.out.println("A["+j+"] "+x[j]+"; "+A.mantissa[j]+"E"+A.exponent[j]);
        }
        double[] y=new double[s];
        for (int j=0; j<s; j++)
            y[j]=Math.pow(400.0,j);
        ScaledProbabilities B=new ScaledProbabilities(y);
        ScaledProbabilities C=new ScaledProbabilities(s);
        multiply(A,B,C);
        for (int j=0;j<s; j++){
            System.out.println("A*B["+j+"] "+x[j]*y[j]+"; "+C.mantissa[j]+"E"+C.exponent[j]);
        }
        A.add(C);
        for (int j=0;j<s; j++){
            System.out.println("A+A*B["+j+"] "+(x[j]+x[j]*y[j])+"; "+A.mantissa[j]+"E"+A.exponent[j]);
        }
        A.add(s/2, 1.0, 100L);
        System.out.println("A["+(s/2)+"] "+A.mantissa[s/2]+"E"+A.exponent[s/2]);
        
    }
}

