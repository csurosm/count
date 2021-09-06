/*
 * SVD.java
 *
 * Created on March 28, 2002, 2:27 PM
 */

package ca.umontreal.iro.matek;

/**
 *
 * @author  csuros
 */

  // ===================================
  // ==
  // == SVD
  // ==
  // ===================================

   /** Singular Value Decomposition.
    *
   <p><em> (based on the JAMA package from <a href="http://www.nist.gov/">NIST</a>)</em>
   <P>
   For an m-by-n matrix A with m >= n, the singular value decomposition is
   an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
   an n-by-n orthogonal matrix V so that A = U*S*V'.
   <P>
   The singular values, sigma[k] = S[k][k], are ordered so that
   sigma[0] >= sigma[1] >= ... >= sigma[n-1].
   <P>
   The singular value decompostion always exists, so the constructor will
   never fail.  The matrix condition number and the effective numerical
   rank can be computed from this decomposition.
   */


import java.util.Vector;
import java.io.BufferedReader;
import java.io.BufferedInputStream;
import java.io.FileReader;
import java.io.FileInputStream;

public class SVD {

/* ------------------------
   Class variables
 * ------------------------ */

 public static double eps=2.2204460492503e-24;

   /** Arrays for internal storage of U and V.
   @serial internal storage of U.
   @serial internal storage of V.
   */
   public double[][] U, V;

   /** Array for internal storage of singular values.
   @serial internal storage of singular values.
   */
   public double[] s;

   /** Row and column dimensions.
   @serial row dimension.
   @serial column dimension.
   */
   private int m, n;
   private boolean wantu;
   private boolean wantv;


/* ------------------------
   Constructor
 * ------------------------ */

   /** Construct the singular value decomposition
   @param A    Rectangular matrix --- DESTROYED!
   @return     Structure to access U, S and V.
   */


   public SVD(double[][] A){
     this(A, true, true);
   }

   public SVD (double[][] A, boolean WantU, boolean WantV) {


      // Derived from LINPACK code.
      // Initialize.
      wantu=WantU;
      wantv=WantV;

      m = A.length;
      n = A[0].length;

      int nu = Math.min(m,n);
      s = new double [Math.min(m+1,n)];
      U = (wantu?new double [m][nu]:null);
      V = (wantv?new double [n][n]:null);
      double[] e = new double [n];
      double[] work = new double [m];

      // Reduce A to bidiagonal form, storing the diagonal elements
      // in s and the super-diagonal elements in e.

      int nct = Math.min(m-1,n);
      int nrt = Math.max(0,Math.min(n-2,m));
      for (int k = 0; k < Math.max(nct,nrt); k++) {
         if (k < nct) {

            // Compute the transformation for the k-th column and
            // place the k-th diagonal in s[k].
            // Compute 2-norm of k-th column without under/overflow.
            s[k] = 0;
            for (int i = k; i < m; i++) {
               s[k] = hypot(s[k],A[i][k]);
            }
            if (s[k] != 0.0) {
               if (A[k][k] < 0.0) {
                  s[k] = -s[k];
               }
               for (int i = k; i < m; i++) {
                  A[i][k] /= s[k];
               }
               A[k][k] += 1.0;
            }
            s[k] = -s[k];
         }
         for (int j = k+1; j < n; j++) {
            if ((k < nct) & (s[k] != 0.0))  {

            // Apply the transformation.

               double t = 0;
               for (int i = k; i < m; i++) {
                  t += A[i][k]*A[i][j];
               }
               t = -t/A[k][k];
               for (int i = k; i < m; i++) {
                  A[i][j] += t*A[i][k];
               }
            }

            // Place the k-th row of A into e for the
            // subsequent calculation of the row transformation.

            e[j] = A[k][j];
         }
         if (wantu & (k < nct)) {

            // Place the transformation in U for subsequent back
            // multiplication.

            for (int i = k; i < m; i++) {
               U[i][k] = A[i][k];
            }
         }
         if (k < nrt) {

            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e[k].
            // Compute 2-norm without under/overflow.
            e[k] = 0;
            for (int i = k+1; i < n; i++) {
               e[k] = hypot(e[k],e[i]);
            }
            if (e[k] != 0.0) {
               if (e[k+1] < 0.0) {
                  e[k] = -e[k];
               }
               for (int i = k+1; i < n; i++) {
                  e[i] /= e[k];
               }
               e[k+1] += 1.0;
            }
            e[k] = -e[k];
            if ((k+1 < m) & (e[k] != 0.0)) {

            // Apply the transformation.

               for (int i = k+1; i < m; i++) {
                  work[i] = 0.0;
               }
               for (int j = k+1; j < n; j++) {
                  for (int i = k+1; i < m; i++) {
                     work[i] += e[j]*A[i][j];
                  }
               }
               for (int j = k+1; j < n; j++) {
                  double t = -e[j]/e[k+1];
                  for (int i = k+1; i < m; i++) {
                     A[i][j] += t*work[i];
                  }
               }
            }
            if (wantv) {

            // Place the transformation in V for subsequent
            // back multiplication.

               for (int i = k+1; i < n; i++) {
                  V[i][k] = e[i];
               }
            }
         }
      }

      // Set up the final bidiagonal matrix or order p.


      int p = Math.min(n,m+1);
      if (nct < n) {
         s[nct] = A[nct][nct];
      }
      if (m < p) {
         s[p-1] = 0.0;
      }
      if (nrt+1 < p) {
         e[nrt] = A[nrt][p-1];
      }
      e[p-1] = 0.0;

      // If required, generate U.

      if (wantu) {
         for (int j = nct; j < nu; j++) {
            for (int i = 0; i < m; i++) {
               U[i][j] = 0.0;
            }
            U[j][j] = 1.0;
         }
         for (int k = nct-1; k >= 0; k--) {
            if (s[k] != 0.0) {
               for (int j = k+1; j < nu; j++) {
                  double t = 0;
                  for (int i = k; i < m; i++) {
                     t += U[i][k]*U[i][j];
                  }
                  t = -t/U[k][k];
                  for (int i = k; i < m; i++) {
                     U[i][j] += t*U[i][k];
                  }
               }
               for (int i = k; i < m; i++ ) {
                  U[i][k] = -U[i][k];
               }
               U[k][k] = 1.0 + U[k][k];
               for (int i = 0; i < k-1; i++) {
                  U[i][k] = 0.0;
               }
            } else {
               for (int i = 0; i < m; i++) {
                  U[i][k] = 0.0;
               }
               U[k][k] = 1.0;
            }
         }
      }

      // If required, generate V.

      if (wantv) {
         for (int k = n-1; k >= 0; k--) {
            if ((k < nrt) & (e[k] != 0.0)) {
               for (int j = k+1; j < nu; j++) {
                  double t = 0;
                  for (int i = k+1; i < n; i++) {
                     t += V[i][k]*V[i][j];
                  }
                  t = -t/V[k+1][k];
                  for (int i = k+1; i < n; i++) {
                     V[i][j] += t*V[i][k];
                  }
               }
            }
            for (int i = 0; i < n; i++) {
               V[i][k] = 0.0;
            }
            V[k][k] = 1.0;
         }
      }

      // Main iteration loop for the singular values.


      int pp = p-1;
      int iter = 0;


      while (p > 0) {
         int k,kase;


         // This section of the program inspects for
         // negligible elements in the s and e arrays.  On
         // completion the variables kase and k are set as follows.

         // kase = 1     if s(p) and e[k-1] are negligible and k<p
         // kase = 2     if s(k) is negligible and k<p
         // kase = 3     if e[k-1] is negligible, k<p, and
         //              s(k), ..., s(p) are not negligible (qr step).
         // kase = 4     if e(p-1) is negligible (convergence).


         // iter is increased by one when kase==3 and reset to 0 when kase==4
//         if (iter>1000){
//          System.out.println("#**SVD() too many iteration steps.");
//          break;
//        }

         for (k = p-2; k >= -1; k--) {
            if (k == -1) {
               break;
            }
            if (Math.abs(e[k]) <= eps*(Math.abs(s[k]) + Math.abs(s[k+1]))) {
               e[k] = 0.0;
               break;
            }
         }
         if (k == p-2) {
            kase = 4;
         } else {
            int ks;
            for (ks = p-1; ks >= k; ks--) {
               if (ks == k) {
                  break;
               }
               double t = (ks != p ? Math.abs(e[ks]) : 0.) +
                          (ks != k+1 ? Math.abs(e[ks-1]) : 0.);
               if (Math.abs(s[ks]) <= eps*t)  {
                  s[ks] = 0.0;
                  break;
               }
            }
            if (ks == k) {
               kase = 3;
            } else if (ks == p-1) {
               kase = 1;
            } else {
               kase = 2;
               k = ks;
            }
         }
         k++;

         // Perform the task indicated by kase.


         switch (kase) {

            // Deflate negligible s(p).

            case 1: {
               double f = e[p-2];
               e[p-2] = 0.0;
               for (int j = p-2; j >= k; j--) {
                  double t = hypot(s[j],f);
                  double cs = s[j]/t;
                  double sn = f/t;
                  s[j] = t;
                  if (j != k) {
                     f = -sn*e[j-1];
                     e[j-1] = cs*e[j-1];
                  }
                  if (wantv) {
                     for (int i = 0; i < n; i++) {
                        t = cs*V[i][j] + sn*V[i][p-1];
                        V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
                        V[i][j] = t;
                     }
                  }
               }
            }
            break;

            // Split at negligible s(k).

            case 2: {
               double f = e[k-1];
               e[k-1] = 0.0;
               for (int j = k; j < p; j++) {
                  double t = hypot(s[j],f);
                  double cs = s[j]/t;
                  double sn = f/t;
                  s[j] = t;
                  f = -sn*e[j];
                  e[j] = cs*e[j];
                  if (wantu) {
                     for (int i = 0; i < m; i++) {
                        t = cs*U[i][j] + sn*U[i][k-1];
                        U[i][k-1] = -sn*U[i][j] + cs*U[i][k-1];
                        U[i][j] = t;
                     }
                  }
               }
            }
            break;

            // Perform one qr step.

            case 3: {

               // Calculate the shift.

               double scale = Math.max(Math.max(Math.max(Math.max(
                       Math.abs(s[p-1]),Math.abs(s[p-2])),Math.abs(e[p-2])),
                       Math.abs(s[k])),Math.abs(e[k]));
               double sp = s[p-1]/scale;
               double spm1 = s[p-2]/scale;
               double epm1 = e[p-2]/scale;
               double sk = s[k]/scale;
               double ek = e[k]/scale;
               double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
               double c = (sp*epm1)*(sp*epm1);
               double shift = 0.0;
               if ((b != 0.0) | (c != 0.0)) {
                  shift = Math.sqrt(b*b + c);
                  if (b < 0.0) {
                     shift = -shift;
                  }
                  shift = c/(b + shift);
               }
               double f = (sk + sp)*(sk - sp) + shift;
               double g = sk*ek;

               // Chase zeros.

               for (int j = k; j < p-1; j++) {
                  double t = hypot(f,g);
                  double cs = f/t;
                  double sn = g/t;
                  if (j != k) {
                     e[j-1] = t;
                  }
                  f = cs*s[j] + sn*e[j];
                  e[j] = cs*e[j] - sn*s[j];
                  g = sn*s[j+1];
                  s[j+1] = cs*s[j+1];
                  if (wantv) {
                     for (int i = 0; i < n; i++) {
                        t = cs*V[i][j] + sn*V[i][j+1];
                        V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
                        V[i][j] = t;
                     }
                  }
                  t = hypot(f,g);
                  cs = f/t;
                  sn = g/t;
                  s[j] = t;
                  f = cs*e[j] + sn*s[j+1];
                  s[j+1] = -sn*e[j] + cs*s[j+1];
                  g = sn*e[j+1];
                  e[j+1] = cs*e[j+1];
                  if (wantu && (j < m-1)) {
                     for (int i = 0; i < m; i++) {
                        t = cs*U[i][j] + sn*U[i][j+1];
                        U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
                        U[i][j] = t;
                     }
                  }
               }
               e[p-2] = f;
               iter = iter + 1;

//                System.out.println("#**SVD() 3 p="+p+", k="+k+", iter="+iter);



            }
            break;

            // Convergence.

            case 4: {

               // Make the singular values positive.

               if (s[k] <= 0.0) {
                  s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
                  if (wantv) {
                     for (int i = 0; i <= pp; i++) {
                        V[i][k] = -V[i][k];
                     }
                  }
               }

               // Order the singular values.

               while (k < pp) {
                  if (s[k] >= s[k+1]) {
                     break;
                  }
                  double t = s[k];
                  s[k] = s[k+1];
                  s[k+1] = t;
                  if (wantv && (k < n-1)) {
                     for (int i = 0; i < n; i++) {
                        t = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = t;
                     }
                  }
                  if (wantu && (k < m-1)) {
                     for (int i = 0; i < m; i++) {
                        t = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = t;
                     }
                  }
                  k++;
               }
               iter = 0;
               p--;
            }
            break;
         }
      }


   }

/* ------------------------
   Public Methods
 * ------------------------ */


   public double[][] getU() {
    return (wantu?AuxAlgebra.arrayClone(this.U):null);
   }

   public double[][] getV() {
     return (wantv?AuxAlgebra.arrayClone(this.V):null);
   }

   /** Return the one-dimensional array of singular values
   * @return     diagonal of S.
   */
   public double[] getS () {
       return AuxAlgebra.arrayClone(s);
  }

   /** Two norm
   @return     max(S)
   */

   public double norm2 () {
      return s[0];
   }

   /** Two norm condition number
   @return     max(S)/min(S)
   */

   public double cond () {
      return s[0]/s[Math.min(m,n)-1];
   }

   /** Effective numerical matrix rank
   @return     Number of nonnegligible singular values.
   */

   public int rank () {
      double tol = Math.max(m,n)*s[0]*eps;
      int r = 0;
      for (int i = 0; i < s.length; i++) {
         if (s[i] > tol) {
            r++;
         }
      }
      return r;
   }

   /**
    * @return sqrt(a*a+b*b) calculated with numerically sensible method
   */
   public static double hypot(double a, double b) {
      double r;
      if (Math.abs(a) > Math.abs(b)) {
         r = b/a;
         r = Math.abs(a)*Math.sqrt(1+r*r);
      } else if (b != 0) {
         r = a/b;
         r = Math.abs(b)*Math.sqrt(1+r*r);
      } else {
         r = 0.0;
      }
      return r;
   }

   public static double[][] diag(double[] v){
    int l=v.length;
    double[][] X= new double[l][l];
    for (int i = 0; i < l; i++) {
    //  for (int j = 0; j < l; j++) {
    //    X[i][j] = 0.0;
    //  }
      X[i][i] = v[i];
    }
    return X;
   }
   
  // ===================================
  // ==
  // == Covariance
  // ==
  // ===================================
  public static double[][] Covariance(double[][] A)
  {
    return _Covariance(AuxAlgebra.arrayClone(A));
  }
  private static double[][] _Covariance(double[][] A)
  {
      
      // computes Covariance of A --- updates A by subtracting the mean

      AuxAlgebra.MeanSubtract(A);

      double[][] At=null;
      double[][] dot_matrix=null;
      boolean transpose_ok=true;

      int nr=A.length;
      int nc=A[0].length;
      try 
      {
          dot_matrix= new double[nr][nr];
      } catch (OutOfMemoryError x) 
      {
        int missing=nr*nr/(128*1024);
        throw new OutOfMemoryError(
        "Covariance: not enough memory for covariance matrix."
        +" --- need about "+missing+"M.");
        //return null;
      }
      for (int i=0; i<nr; i++)
      {
        for (int j=i+1; j<nr; j++)
          dot_matrix[j][i]=dot_matrix[i][j]=AuxAlgebra.innerProduct(A[i],A[j]);
        dot_matrix[i][i]=AuxAlgebra.innerProduct(A[i],A[i]);
      }
      double fact=1.0/A[0].length;
      for (int i=0;i<dot_matrix.length;i++)
        dot_matrix[i]=AuxAlgebra.multiply(fact,dot_matrix[i]);

      return dot_matrix;
  }
  
  public static SVD PCA(double[][] M)
  {
    return _PCA(AuxAlgebra.arrayClone(M));
  }
  private static SVD _PCA(double[][] M)
  {
    // updates eigenvalues[], returns eigenvectors
    // subtracts the means in M!!
    double[][] Cov  = _Covariance(M);

    //System.out.println(DoubleArraytoString(Cov,"Cov"));
    SVD SCov = new SVD(Cov,true,false);//does not want V
    
    return SCov;
  }   
  
  public static void main(String[] args) throws Exception
  {
      String file = args[0];
      Vector<Vector<Double>> input_matrix = new Vector<Vector<Double>>();
      String line = null;
      BufferedReader BR = new BufferedReader(new FileReader(file));
      do
      {
          line = BR.readLine();
          if (line == null)
              break;
          String[] e = line.split("\\t");
          Vector<Double> row = new Vector<Double>();
          for (int i=0; i<e.length; i++)
              row.add(Double.parseDouble(e[i]));
          input_matrix.add(row);
      } while (line != null);
      int nr = input_matrix.size();
      double[][] input_array = new double[nr][];
      for (int i=0; i<nr; i++)
      {
          Vector<Double> row = input_matrix.get(i);
          int nc = row.size();
          double[] x = new double[nc];
          for (int j=0; j<nc; j++)
          {
              Double D = row.get(j);
              x[j] =D.doubleValue();
          }
          input_array[i] = x;
      }
      SVD svd = PCA(input_array);
      
      double[][] U = svd.getU();
      double[] S = svd.getS();
      for (int i=0; i<U.length; i++)
      {
          for (int j=0; j<U[i].length; j++)
          {
              if (j>0)
                  System.out.print('\t');
              System.out.print(U[i][j]);
          }
          System.out.println();
      }
      System.out.println("Eigenvalues");
      for (int i=0; i<S.length; i++)
          System.out.println(i+"\t"+S[i]);
  }
}