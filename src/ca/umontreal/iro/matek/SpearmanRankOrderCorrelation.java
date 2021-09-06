
package ca.umontreal.iro.matek;

import ca.umontreal.iro.banality.IndexSort;
import ca.umontreal.iro.banality.Functions;
/**
 * Spearman Rank-Order Correlation Coefficient 
 * 
 * Implementation from NR 14.6
 *  
 * @author csuros
 */
public class SpearmanRankOrderCorrelation 
{
    private double[] data1;
    private double[] data2;
    
    public SpearmanRankOrderCorrelation(double[] data1, double[] data2)
    {
        this.data1 = data1;
        this.data2 = data2;
        spear();
    }
    
    private void spear()
    {
        int n = data1.length;
        double[] wskap1 = new double[n];
        double[] wskap2 = new double[n];
        System.arraycopy(data1, 0, wskap1, 0, n);
        System.arraycopy(data2, 0, wskap2, 0, n);
        sort2(wskap1,wskap2);
        double sf = crank(wskap1);
        sort2(wskap2,wskap1);
        double sg = crank(wskap2);
        
        corr_d = 0.0;
        for (int i=0; i<n; i++)
        {
            double d = wskap1[i]-wskap2[i];
            //System.out.println("#**SROC.sp "+i+"\t"+wskap1[i]+"\t"+wskap2[i]+"\t//"+data1[i]+"\t"+data2[i]);
            corr_d += d*d;
        }
        double en = (double)n;
        double en3n = en*en*en-en;
        System.out.println("#**SROC.sp sf "+sf+"\tsg "+sg+"\tn3n "+en3n);
        corr_aved=en3n/6.0-(sf+sg)/12.0; // Expectation value of D, 
        double fac=(1.0-sf/en3n)*(1.0-sg/en3n); 
        corr_vard=((n-1.0)*n*n*(n+1.0)*(n+1.0)/36.0)*fac; // and variance of D give 
        corr_zd=(corr_d-corr_aved)/Math.sqrt(corr_vard); // number of standard deviations and significance.
        corr_probd=Functions.erfcc(Math.abs(corr_zd)/Math.sqrt(2.0)); 
        corr_rs=(1.0-(6.0/en3n)*(corr_d+(sf+sg)/12.0))/Math.sqrt(fac); // Rank correlation coefficient, 
        fac=(corr_rs+1.0)*(1.0-(corr_rs)); 
        if(fac>0.0)
        { 
            double t=(corr_rs)*Math.sqrt((n-2.0)/fac); // and its t value, 
            double df=n-2.0; 
            corr_probrs=Functions.betai(0.5*df,0.5,df/(df+t*t)); // give its significance. 
        } else 
            corr_probrs=0.0; 
    }
    
    public String toString()
    {
        StringBuffer sb = new StringBuffer("Spearman rank-order correlation [");
        sb.append("D=");
        sb.append(corr_d);
        sb.append(" (aveD "+corr_aved+", varD "+corr_vard+"), z="+corr_zd+" P(D)="+corr_probd);
        sb.append("; rs="+corr_rs);
        sb.append(", P(rs)="+corr_probrs);
        sb.append(']');
        return sb.toString();
    }
    
    private double corr_d;
    private double corr_aved;
    private double corr_vard;
    private double corr_zd;
    private double corr_probd;
    private double corr_rs;
    private double corr_probrs;
    
    private void sort2(double[] x, double[] y)
    {
        int[] xi = IndexSort.quickSort(x);
        double[] w = IndexSort.permute(x, xi);
        //for (int i=0; i<xi.length; i++)
        //{
        //    System.out.println("#**SROC.s2 "+i+"\t"+xi[i]+"\t"+w[i]);
        //}
        System.arraycopy(w,0,x,0,w.length);
        w = IndexSort.permute(y, xi);
        System.arraycopy(w,0,y,0,w.length);
    }
    
    /**
     * Given a sorted array w, replaces the elements by their rank, 
     * including mid-ranking of ties, 
     * and returns the sum of <var>f</var><sup>3</sup>-<var>f</var>,
     * where <var>f</var> is the number of elements in each tie. 
     * 
     * @param w input array of sorted values
     * @return 
     */
    private double crank(double w[]) 
    { 
        double total = 0.0;
        for (int j=0; j<w.length; )
        {
            int k = j+1;
            while (k<w.length && w[k]==w[j]) k++;
            double rank = (j+k-1.0)/2.0;
            for (int i=j; i<k; i++)
                w[i]=rank;
            int m = k-j;
            double em = (double)m; // m^2 could overflow with int calculations
            double m3m = em*em*em-em ; // 
            total += m3m;
            
            System.out.println("#**SROC.cr block "+j+".."+(k-1)+"\t"+m+"/"+m3m+"\tmid "+rank);
            j=k;
        }
        return total;
    }

}
