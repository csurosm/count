/*
 * AuxAlgebra.java
 *
 * Created on March 28, 2002, 2:35 PM
 */

package ca.umontreal.iro.matek;

/**
 * Collection of static  methods for double vectors and arrays.
 *
 * @author  csuros
 */
public class AuxAlgebra {
   /**
    * Efficient cloning of a double[][] array.
    * (since clone() is shallow)
    */
   public static final double[][] arrayClone(double[][] A){
    if (A==null) return null;
    int nr=A.length;
    int nc=A[0].length;
    int i;
    double[][] M=new double[nr][nc];
    for (i=0;i<nr;i++)
      System.arraycopy(A[i],0,M[i],0,nc);
    return M;
  }

   /**
    * A wrapper for clone() allowing null.
    */
  public static final double[] arrayClone(double[] v){
    if (v==null) return null;
    return (double[])v.clone();
//    int len=v.length;
//    double[] w=new double[len];
//
//    System.arraycopy(v,0,w,0,len);
//    return w;
  }


  /**
   * A wrapper for clone() allowing null.
   */
  public static final int[] arrayClone(int[] v){
    if (v==null) return null;
    return (int[]) v.clone();
//    int len=v.length;
//    int[] w=new int[len];
//
//    System.arraycopy(v,0,w,0,len);
//    return w;
  }



  /**
   * Calculates squared distance between two vectors as
   * sum_{j=0}^d (u[j]-v[j])^2
   */
  public static final double Sqr_distance_between_vectors(double[] u, double[] v){
    //if (u.length != v.length) throw new java.lang.IllegalArgumentException();
    int d=u.length;
    double D=0.;
    for (int i=0; i<d; i++){
      double x=u[i]-v[i];
      D+= x*x;
    }
    return D;
  }

  // ===================================
  // ==
  // == Constant vectors
  // ==
  // ===================================
  /**
   * Fast initialization of constant vectors.
   * @return  double vector of length len, with all entries set to val
   */
  public static final double[] Constant_vector(int len, double val){
    if (len<=0) return null;

    double v[] = new double[len];

    if (len<8){
      int i;
      for (i=0;i<len;i++) v[i]=val;
    } else {
      int len2 = len/2;
      double[] v_half=Constant_vector(len2, val);
      System.arraycopy(v_half,0,v,0,len2);
      System.arraycopy(v_half,0,v,len2,len2);
      v_half=null;
      v[len-1]=val;
    }
    return v;
  }

  /**
   * Fast initialization of all-zero vector.
   */
  public static final double[] Zero_vector(int len){
    return new double[len];
    //return Constant_vector(len,0);
  }

  /**
   * Fast initalization of all-ones vector.
   */
  public static final double[] One_vector(int len){
    return Constant_vector(len,1);
  }
  
  public static final double average(double[] x)
  {
      double m = 0.0;
      int n = x.length;
      for (int i=0; i<n; i++)
          m += x[i];
      return m/n;
  }
  
  public static final double[] add(double a, double [] x)
  {
      int n = x.length;
      double[] y = new double[n];
      for (int i=0; i<n; i++)
          y[i] = x[i]+a;
      return y;
  }
  
  public static final double[] multiply(double a, double [] x)
  {
      int n = x.length;
      double[] y = new double[n];
      for (int i=0; i<n; i++)
          y[i] = x[i]*a;
      return y;
  }

  
  public static double[] MeanSubtract(double[][] A){ return _MeanSubtract(A);}
  private static double[] _MeanSubtract(double[][] A){
    // returns the vector of row-wise means
    double[] mean=new double[A.length];

    // set 0 mean in each row of A
    int i;
    for(i=0;i<A.length;i++){
      mean[i]=average(A[i]);
      A[i]=add(-mean[i], A[i]);
    }
    return mean;
  }
  
  public static final double innerProduct(double[] x, double[] y)
  {
      double v = 0.0;
      int n = x.length;
      for (int i=0; i<n; i++)
          v += x[i]*y[i];
      return v;
  }

  // ===================================
  // ==
  // == debug stuff
  // ==
  // ===================================

  public static String ArraytoString(double [][] A, String name, String prefix){
    StringBuffer s=new java.lang.StringBuffer(
      prefix+"{//matrix "+name+" size:"+A.length+"x"+A[0].length);

    int i,j;

    int jmax=Math.min(A[0].length,20);
    if (A[0].length>20) s=s.append(" -- truncated at 20 columns");
    s=s.append("\n");

    for (i=0; i<A.length; i++){
      s = s.append(prefix+" { ");
      for (j=0; j<jmax; j++){
        s = s.append(A[i][j]).append(",");
      }
      s=s.append("},\n");
    }
    s=s.append(prefix+"}");

    return s.toString();
  }

  public static String ArraytoString(double[][] A){
    return ArraytoString(A,"","");
  }

  public static String ArraytoString(double A[][], String name){
    return ArraytoString(A,name,"");
  }

  public static String ArraytoString(double v[], String name, String prefix){
    double[][] A = new double[1][];
    A[0]=v;
    return ArraytoString(A,name,prefix);
  }

  public static String ArraytoString(double v[], String name){
    return ArraytoString(v,name,"");
  }

  public static String ArraytoString(double v[]){
    return ArraytoString(v,"","");
  }

  public static String ArraytoString(int [][] A, String name, String prefix){
    StringBuffer s=new java.lang.StringBuffer(
      prefix+"{//matrix "+name+" size:"+A.length+"x"+A[0].length);

    int i,j;

    int jmax=Math.min(A[0].length,20);
    if (A[0].length>20) s=s.append(" -- truncated at 20 columns");
    s=s.append("\n");

    for (i=0; i<A.length; i++){
      s = s.append(prefix+" { ");
      for (j=0; j<jmax; j++){
        s = s.append(A[i][j]).append(",");
      }
      s=s.append("},\n");
    }
    s=s.append(prefix+"}");

    return s.toString();
  }

  public static String ArraytoString(int[][] A){
    return ArraytoString(A,"","");
  }

  public static String ArraytoString(int[][] A, String name){
    return ArraytoString(A,name,"");
  }
  public static String ArraytoString(int v[], String name, String prefix){
    int[][] A = new int[1][];
    A[0]=v;
    return ArraytoString(A,name,prefix);
  }

  public static String ArraytoString(int v[], String name){
    return ArraytoString(v,name,"");
  }

  public static String ArraytoString(int v[]){
    return ArraytoString(v,"","");
  }


  public static String getFormattedStr(double d){return getFormattedStr(d,3);}

  /**
   * @return IEEE format string rounded to the given number of decimals after the point
   */
  public static String getFormattedStr(double d, int decimals){
    double d0=d;
    String retval;

    double rounder=1;
    if (decimals == 1){
      rounder = 10.;
    } else if (decimals == 2){
      rounder = 100.;
    } else if (decimals == 3) {
      rounder = 1000.;
    } else {
      rounder = Math.exp(decimals*Math.log(10.));
    }

    if (d0==Double.NEGATIVE_INFINITY)
      retval = "\u2212\u221E";
    else if (d0==Double.POSITIVE_INFINITY)
      retval = "\u221E";
    else if (d0==Double.NaN)
      retval = "NaN";
    else {
      if (decimals == 0)
        retval = Integer.toString((int)(d0+0.5));
      else {
        d0=(((int)(d0*rounder+0.5))/rounder);
        retval = Double.toString(d0);
      }
    }

    return retval;
  }

}