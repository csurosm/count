/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ca.umontreal.iro.banality;

/**
 * Experiments for comparing probabilistic counters
 * 
 * @author csuros
 */

import java.util.Random;

import java.util.Vector;

public class probCountExperiment 
{
    private enum DataType {ESTIMATE, VALUE, EXPONENT};
            
    private static final DataType DEFAULT_DATA_TYPE = DataType.ESTIMATE;
    
    // sampling the estimates at this many points between two magnitudes (base 10)
    // checkpoints are approximately equidistant in log-space
    private static final int RESOLUTION = 33;
    
    private static void go(int Nmag, int precision, int num_repeats, DataType data_type)
    {
        int N = 1;
        for (int i=0; i<Nmag; i++) N *= 10;
        
        Vector<Integer> check_pointsV = new Vector<Integer>();
        
        for (int i=0, d=1; i<Nmag-1; i++, d*=10)
        {
            for (int step=0; step<RESOLUTION; step++)
            {
                double x = Math.pow(10.0, ((double)step)/RESOLUTION);
                int sampled_value = (int)(x*10.0*d);
                if (check_pointsV.size()==0 || check_pointsV.get(check_pointsV.size()-1)!=sampled_value)
                    check_pointsV.add(sampled_value);
            }
            
        }

        check_pointsV.add(N);
        
        int num_checkpoints = check_pointsV.size();
        int[] checkpoints = new int[num_checkpoints];
        for (int i=0; i<num_checkpoints; i++)
            checkpoints[i] = check_pointsV.get(i);
        
        
        int[][] FPest = new int[num_checkpoints][num_repeats];
        int[][] Qest  = new int[num_checkpoints][num_repeats];
        Random RND = new Random(2009);
        
        // calculate 2^precision
        int M = 1<<precision;
                
        double LN2 = Math.log(2);
        for (int rep=1; rep<=num_repeats; rep++)
        {
            int Qcnt = 1;
            int FPcnt = 0;
            int cp_idx = 0;
            for (int n=1; n<=N; n++)
            {
                // update Q-counter
                double x = ((double)Qcnt)/((double)M);
                // random number
                double r = -Math.log(RND.nextDouble())/LN2;
                if (r>x)
                    Qcnt++;

                // update FPcounter
                int y = FPcnt/M;
                if (y==0 || RND.nextInt(1<<y)==0)
                    FPcnt++;


                // checkpoint?
                if (n==checkpoints[cp_idx])
                {
                    FPest[cp_idx][rep-1] = FPcnt;
                    Qest [cp_idx][rep-1] = Qcnt;
                    cp_idx++;
                }
            }

        }
        
        System.out.println("[ % -------------------------- Q-counter");
        for (int cp_idx=0; cp_idx<num_checkpoints; cp_idx++)
        {
            System.out.print("    "+checkpoints[cp_idx]);
            
            double[] estimates = new double[num_repeats];
            if (data_type == DataType.ESTIMATE)
                for (int i=0; i<num_repeats; i++)
                    estimates[i] = invertQCounter(Qest[cp_idx][i],M)/checkpoints[cp_idx];
            else if (data_type == DataType.VALUE)
                for (int i=0; i<num_repeats; i++)
                    estimates[i] = -Qest[cp_idx][i]+(M*Math.log(checkpoints[cp_idx])/LN2);
            else if (data_type == DataType.EXPONENT)
                for (int i=0; i<num_repeats; i++)
                    estimates[i] = (1 << ((Qest[cp_idx][i]/M)))/((double)checkpoints[cp_idx])*M;
            
            printEstimates(estimates);
            
            System.out.println();
        }
        System.out.println("] /data.Q exch def");
        
        System.out.println("[ % -------------------------- FP-counter");
        for (int cp_idx=0; cp_idx<num_checkpoints; cp_idx++)
        {
            System.out.print("    "+checkpoints[cp_idx]);

            double[] estimates = new double[num_repeats];
            if (data_type == DataType.ESTIMATE)
                for (int i=0; i<num_repeats; i++)
                    estimates[i] = invertFPCounter(FPest[cp_idx][i],precision)/checkpoints[cp_idx];
            else if (data_type == DataType.VALUE)
                for (int i=0; i<num_repeats; i++)
                    estimates[i] = -FPest[cp_idx][i]+(M*Math.log(checkpoints[cp_idx])/LN2);
            else if (data_type == DataType.EXPONENT)
                for (int i=0; i<num_repeats; i++)
                {
                    estimates[i] = (1 << ((FPest[cp_idx][i]/M)))/((double)checkpoints[cp_idx]+2.*M)*M;
                    //int t = FPest[cp_idx][i]/M;
                    //int t_wanted = (int)( Math.log(checkpoints[cp_idx]/M)/LN2);
                    //estimates[i] = Math.pow(2.0,t-t_wanted);
                }
                
            printEstimates(estimates);

            System.out.println();
        }
        System.out.println("] /data.FP exch def");
        
    }
    

    private static void printEstimates(double[] estimates)
    {
            double[] stats = computeStatistics(estimates);
            double mean = stats[0];
            double stdev = stats[1];
            
            System.out.print("\t"+mean+"\t"+stdev);
            
            double E2sd_neg = stats[0]-2.0*stats[1];
            double Esd_neg  = stats[0]-    stats[1];
            double E2sd_pos = stats[0]+2.0*stats[1];
            double Esd_pos  = stats[0]+    stats[1];
            
            
            System.out.print("\t[");
            for (int i=0; i<estimates.length; i++)
            {
                if (estimates[i]<E2sd_neg || estimates[i]>E2sd_pos )//|| (estimates[i]>Esd_neg && estimates[i]<Esd_pos))
                    System.out.print(" "+estimates[i]);
            }
            System.out.print("]");
    }
    
    /**
     * computes average and stdev on an array
     * 
     * @param A
     * @return {avg, stdev} pair
     */
    public static double[] computeStatistics(double[] A)
    {
        double s = 0.0;
        double s2 = 0.0;
        
        for (int j=0; j<A.length; j++)
        {
            s += A[j];
            s2 += A[j]*A[j];
        }
        
        double mean = s/A.length;
        
        double z = s2-A.length*mean*mean;
        if (z<0.0) // roundoff error
            z=0.0;
        
        double stdev = Math.sqrt(z/(A.length-1.0));
        
        double[] retval = new double[2];
        
        retval[0]=mean;
        retval[1]=stdev;
        return retval;
    }
    
    private static double invertQCounter(int x, int r)
    {
        double q = Math.pow(2.0, 1.0/r);
        double inv = (Math.pow(q, x)-1.0)/(q-1.0);
        inv = ((int)(10.0*inv+0.5))/10.0;
        return (inv-1.0);
    }
    

    private static double invertFPCounter(int x, int d)
    {
        int t = x >>> d;
        int u = x ^ (t<<d);
        int M = 1<<d;
        
        double w = Math.pow(2.0,t)*(u+M)-M;

        return w;
    }
    
    
    public static void main(String[] args)
    {
        if (args.length < 2)
        {
            System.err.println("Call as pCE <precision> <runs>" +
                    "\n\tCounts until 1 million  using Q-base and FP-counters of the same precision");
            System.exit(2009);
        }
        
        DataType data_type = DEFAULT_DATA_TYPE;
        int precision = Integer.parseInt(args[0]);
        int num_repeats = Integer.parseInt(args[1]);
        if (args.length > 2)
        {
            if ("estimate".equals(args[2]))
            {
                data_type = DataType.ESTIMATE;
            } else if ("counter".equals(args[2]))
            {
                data_type = DataType.VALUE;
            } else if ("exponent".equals(args[2]))
            {
                data_type = DataType.EXPONENT;
            } else
                throw new IllegalArgumentException("Unrecognized datatype (must be counter,exponent, or estimate");
        }
        
        
        go (4,precision, num_repeats, data_type);
    }
    
    

}
