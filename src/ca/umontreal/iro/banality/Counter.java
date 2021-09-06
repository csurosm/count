/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ca.umontreal.iro.banality;

/**
 *
 * @author csuros
 */

import java.util.Random;

public class Counter 
{
    private int d;
    public Counter(int precision)
    {
        this.d = precision;
        RND = new Random();
        value = significant_digits=nulls=0;
        power2 = 1<<d;
        max_n = 1;
        System.out.println("#** init pow2="+power2);
    }
    
    private int value;
    private Random RND;
    private int significant_digits;
    private int nulls;
    private int power2;
    private int max_n;
    
    public void increment()
    {
        if (nulls==0 && significant_digits<power2)
        {
            significant_digits++;
        } else
        {
            int x = RND.nextInt(max_n);
            if (x==0)
            {
                significant_digits++;
                //System.out.println("#**record "+significant_digits+"/"+power2+" E"+nulls);
            }
        }
        if (significant_digits==power2)
        {
            significant_digits/=2;
            max_n *= 2;
            nulls++;
        }
        
        { // Morris' counter
            double x = RND.nextDouble();
            double t = Math.pow(2.0, -((double)value)/((double)(power2/2)));
            if (x<t)
                value++;
        }
    }
    
    private void report()
    {
        double est_Flajolet = estimate(value, d-1);
        int est_mine = significant_digits<<nulls;
        //int x = nulls<<d+significant_digits;
        System.out.println("#** Estimates: "+est_Flajolet+" (Morris/Flajolet);\t"+est_mine+" (mine)"+"\t"+value);//+"; "+x);

    }
    
    private static final double estimate (int x, int delta)
    {
        double q = Math.pow(2.0, -1.0/(1<<delta));
        double q1=1.0/q;
        double est_Flajolet = (Math.pow(q, -x)-q1)/(q1-1.0);
        return est_Flajolet;
    }
    
    public static void main(String[] args)
    {
       int prec = Integer.parseInt(args[0]);
       int N = Integer.parseInt(args[1]);
       Counter C = new Counter(prec);
       
       for (int j=0; j<N; j++)
           C.increment();
       C.report();
       
    }

}
