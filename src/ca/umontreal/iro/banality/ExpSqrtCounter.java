
package ca.umontreal.iro.banality;
import java.util.Vector;
import java.util.Random;

/**
 *
 * @author csuros
 */
public class ExpSqrtCounter
{
    // sampling the estimates at this many points between two magnitudes (base 10)
    // checkpoints are approximately equidistant in log-space
    private static final int RESOLUTION = 15;

    public ExpSqrtCounter()
    {
        counter = 0;
        setProbability();
        RND = new Random();
    }

    private int counter;
    private double prob_increment;
    private Random RND;

    private void setProbability()
    {
        double e1= Math.exp(Math.sqrt(counter+1.0));
        double e0 = Math.exp(Math.sqrt((double)counter));
        double p = 1.0/(e1-e0);

        //double e = Math.exp(Math.sqrt(-(double)counter));
        //double d = Math.exp(Math.sqrt(counter+1.0)-Math.sqrt((double)counter))-1.0;

        //if (e/d != p)
        //    System.out.println("#*** INCREMENT "+p+"\t"+(e/d));
        prob_increment =p;
    }

    public void increment()
    {
        if (RND.nextDouble()<prob_increment)
        {
            counter++;
            setProbability();
        }
    }

    public int getCounterValue()
    {
        return counter;
    }

    public double getN()
    {
        return getN(counter);
    }

    public static double getN(int counter)
    {
        return Math.exp(Math.sqrt(counter))-1;
    }

    public void setCounterValue(int counter)
    {
        this.counter = counter;
        setProbability();
    }

    public void reset()
    {
        setCounterValue(0);
    }

    private enum DataType {ESTIMATE, VALUE};

    private static final DataType DEFAULT_DATA_TYPE = DataType.ESTIMATE;

    private static void go(int Nmag, int num_repeats, DataType data_type)
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


        int[][] SqEst = new int[num_checkpoints][num_repeats];

        for (int rep=1; rep<=num_repeats; rep++)
        {
            ExpSqrtCounter Ctr = new ExpSqrtCounter();
            int cp_idx = 0;
            for (int n=1; n<=N; n++)
            {
                Ctr.increment();
                // checkpoint?
                if (n==checkpoints[cp_idx])
                {
                    SqEst[cp_idx][rep-1] = Ctr.getCounterValue();
                    cp_idx++;
                }
            }
        } // for rep

        System.out.println("[ % -------------------------- ExpSqrt counter");
        for (int cp_idx=0; cp_idx<num_checkpoints; cp_idx++)
        {
            System.out.print("    "+checkpoints[cp_idx]);

            double[] estimates = new double[num_repeats];
            if (data_type == DataType.ESTIMATE)
                for (int i=0; i<num_repeats; i++)
                    estimates[i] = getN(SqEst[cp_idx][i])/checkpoints[cp_idx];
            else if (data_type == DataType.VALUE)
                for (int i=0; i<num_repeats; i++)
                    estimates[i] = SqEst[cp_idx][i];

            printEstimates(estimates);

            System.out.println();
        }
        System.out.println("] /data.ExpSqrt exch def");

    } // go

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


    public static void main(String[] args)
    {
        if (args.length < 1)
        {
            System.err.println("Call as ESC <runs> [<data type>]" +
                    "\n\tCounts until 1 million using ExpSqrt counter");
            System.exit(2009);
        }

        DataType data_type = DEFAULT_DATA_TYPE;
        int num_repeats = Integer.parseInt(args[0]);
        if (args.length > 1)
        {
            if ("estimate".equals(args[2]))
            {
                data_type = DataType.ESTIMATE;
            } else if ("counter".equals(args[1]))
            {
                data_type = DataType.VALUE;
            } else
                throw new IllegalArgumentException("Unrecognized datatype (must be counter,exponent, or estimate");
        }


        go (7,num_repeats, data_type);
    }


}
