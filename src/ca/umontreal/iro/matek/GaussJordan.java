/*
 * GaussJordan.java
 *
 * Created on August 8, 2006, 1:05 PM
 */

package ca.umontreal.iro.matek;

/**
 * Linear equation solution by Gauss-Jordan elimination.
 *
 * @author  csuros
 */
public class GaussJordan 
{
    
    /** Not used */
    private GaussJordan() {}

    /**
     *
     *
     * Linear equation solution by Gauss-Jordan elimination.
     *
     * Implementation based on Numerical Recipes 2.1
     *
     * @param a  is an n x n input matrix a[0..n-1][0..n-1]
     * @param b  is an n x m input containing the m right-hand side vectors b[0..n-1][0..m-1]
     * @return On output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution vectors. 
     */
    public static void gaussj(double[][] a,double[][] b)
    {
        int n = a.length;
        int m = b[0].length;
        
        // The integer arrays ipiv, indxr, and indxcare used for bookkeeping on the pivoting.
        int[] indxc=new int[n]; 
        int[] indxr=new int[n]; 
        int[] ipiv =new int[n]; 
        double[] tmp_n = new double[n];
        double[] tmp_m = new double[m];
        for(int j=0;j<n;j++)
            ipiv[j]=0; 
        for(int i=0;i<n;i++)
        { // This is the main loop over the columns to be reduced.
            double big=0.0;
            int irow=-1, icol=-1;
            for(int j=0;j<n;j++) // This is the outer loop of the search for a pivot element.
                if(ipiv[j]!=1)
                    for(int k=0;k<n;k++)
                    { 
                        if(ipiv[k]==0)
                        { 
                            double abs_ajk = Math.abs(a[j][k]);
                            if(abs_ajk>=big)
                            { 
                                big=abs_ajk; 
                                irow=j; 
                                icol=k; 
                            } 
                        } else if(ipiv[k]> 1)
                            throw new ArithmeticException("gaussj:SingularMatrix-1"); 
                    }
            ++(ipiv[icol]); 
            // We now have the pivot element, so we interchange rows, if needed, to put the pivot 
            // element on the diagonal. The columns are not physically interchanged, only relabeled: 
            // indxc[i], the column of the ith pivot element, is the ith column that is reduced, while 
            // indxr[i] is the row in which that pivot element was originally located. 
            // If indxr[i] <> indxc[i] there is an implied column interchange. 
            // With this form of bookkeeping, the 
            // solution b's will end up in the correct order, and the inverse matrix will be scrambled
            // by columns. 
            if(irow!=icol)
            { 
                System.arraycopy(a[irow],0,tmp_n, 0, n);
                System.arraycopy(a[icol],0,a[irow],0,n);
                System.arraycopy(tmp_n,0,a[icol],0,n);
                System.arraycopy(b[irow],0,tmp_m,0,m);
                System.arraycopy(b[icol],0,b[irow],0,m);
                System.arraycopy(tmp_m,0,b[icol],0,m);
            }
            // We are now ready to divide the pivot row by the 
            // pivot element, located at irow and icol.
            indxr[i]=irow; 
            indxc[i]=icol; 
            if(a[icol][icol]==0.0)
                throw new ArithmeticException("gaussj:SingularMatrix-2"); 
            double pivinv=1.0/a[icol][icol]; 
            a[icol][icol]=1.0; 
            for(int l=0;l<n;l++) a[icol][l]*=pivinv; 
            for(int l=0;l<m;l++) b[icol][l]*=pivinv;    
            for(int ll=0;ll<n;ll++) 
                // Next, we reduce the rows... 
                if(ll!=icol) // ...except for the pivot one, of course. 
                { 
                    double dum=a[ll][icol]; 
                    a[ll][icol]=0.0; 
                    for(int l=0;l<n;l++)a[ll][l]-=a[icol][l]*dum; 
                    for(int l=0;l<m;l++)b[ll][l]-=b[icol][l]*dum; 
                } 
        } 
        // This is the end of the main loop over columns of the reduction. 
        // It only remains to unscramble the solution in view of the column
        // interchanges. We do this by interchanging pairs of columns in the
        // reverse order that the permutation was built up. 
        for(int l=n-1;l>=0;l--)
        { 
            if(indxr[l]!=indxc[l]) 
                for(int k=0;k<n;k++) 
                {
                    double tmp = a[k][indxr[l]];
                    a[k][indxr[l]]=a[k][indxc[l]];
                    a[k][indxc[l]]=tmp;
                }
        } 
        // And we are done. 
    }
    
    
    /**
     * Test with random values
     */
    public static void main(String[] args)
    {
        java.util.Random RND = new java.util.Random();
        double[][] A = new double[3][3];
        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
                A[i][j] = RND.nextInt(20);
        double[][] b = new double[3][2];
        for (int i=0; i<3; i++)
            for (int j=0; j<2; j++)
                b[i][j] = RND.nextInt(10);
        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
                System.out.print("\t"+A[i][j]);
            System.out.print("\t|");
            for (int j=0; j<2; j++)
                System.out.print("\t"+b[i][j]);
            System.out.println();
        }
        System.out.println(" solution :");
        gaussj(A,b);
        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
                System.out.print("\t"+A[i][j]);
            System.out.print("\t|");
            for (int j=0; j<2; j++)
                System.out.print("\t"+b[i][j]);
            System.out.println();
        }
    }
    
}
