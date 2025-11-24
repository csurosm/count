package count.matek;

/*
 * Copyright 2023 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Numerical tests used during code development: how to compute ln(1-exp(-x))? 
 * 
 * @author csuros
 *
 */
public class PrecisionTest 
{
	
	private void logTest(PrintStream out)
	{
		
		List<Double> X = new ArrayList<>();
		final double MIN_EPSI = 1e-24;
		final double DELTA_EPSI = 0.01;
		final double LONG_EPSI = 3.0;
		final double BIG_EPSI = 10.0;
		final double MAX_EPSI = 100.0;
		
		double epsi = MIN_EPSI;
		
		
		while (epsi<DELTA_EPSI)
		{
			X.add(epsi);
			epsi *= 10.0;
		}
		epsi = DELTA_EPSI;
		while (epsi<LONG_EPSI)
		{
			X.add(epsi);
			epsi += DELTA_EPSI;
		}
		epsi=BIG_EPSI;
		while (epsi<MAX_EPSI)
		{
			X.add(epsi);
			epsi += BIG_EPSI;
		}
		X.add(epsi);
		
		out.println("# x\tb=e^{-x}\tc=1-b\tc'=-expm1(-x)\t|\td=ln(c)\td'=ln(c')\td''=log1p(-b)");
		for (double x:X)
		{
			double exp_x = Math.exp(-x);
			double exp_x1 = 1.0-Math.exp(-x);
			double expm1_x = -Math.expm1(-x);
			double log_exp_x1 = Math.log(exp_x1);
			double log_expm1_x = Math.log(expm1_x);
			double log1p_exp_x = Math.log1p(-exp_x); 
			
			out.printf("%g\t%g\t%g\t%g\t|\t%g\t%g\t%g\n", x, 
					exp_x, // b
					exp_x1, // c
					expm1_x, // c'
					log_exp_x1, 
					log_expm1_x,
					log1p_exp_x);
			
			
		}
		
	}
	
	private void productTest(PrintStream out) {
		int n = 1000;
		int REP = 4096;
		double[] a = new double[n];
		double[] b = new double[n];
		
		long T1 = 0L;
		long T16 = 0L;
		
		double ss=0.0;
		
		for (int r=0; r<REP; r++) {
			for (int i=0; i<n; i++){
				a[i] = Math.sin(r+i);
				b[i] = Math.log((r+1)*(i+1));
			}
			
			long T0 = System.nanoTime();
			double s = 0.0;
			for (int i=0; i<n; i++) {
				s += a[i]*b[i];
			}
			long dT1 = System.nanoTime()-T0;
			T0 = System.nanoTime();
			double t = 0.0;
			double m = n % 16;
			int j=0;
			while (j<m) {
				t += a[j]*b[j];
				j++;
			}
			while (j < n) {
				double t0 = a[j]*b[j];
				final int j1 = j+1;
				double t1 = a[j1]*b[j1];
				final int j2 = j+2;
				double t2 = a[j2]*b[j2];
				final int j3 = j+3;
				double t3 = a[j3]*b[j3];
				final int j4 = j+4;
				double t4 = a[j4]*b[j4];
				final int j5 = j+5;
				double t5 = a[j5]*b[j5];
				final int j6 = j+6;
				double t6 = a[j6]*b[j6];
				final int j7 = j+7;
				double t7 = a[j7]*b[j7];
				final int j8 = j+8;
				double t8 = a[j8]*b[j8];
				final int j9 = j+9;
				double t9 = a[j9]*b[j9];
				final int j10 = j+10;
				double t10 = a[j10]*b[j10];
				final int j11 = j+11;
				double t11 = a[j11]*b[j11];
				final int j12 = j+12;
				double t12 = a[j12]*b[j12];
				final int j13 = j+13;
				double t13 = a[j13]*b[j13];
				final int j14 = j+14;
				double t14 = a[j14]*b[j14];
				final int j15 = j+15;
				double t15 = a[j15]*b[j15];
				t0 += t1;
				t2 += t3;
				t4 += t5;
				t6 += t7;
				t8 += t9;
				t10 += t11;
				t12 += t13;
				t14 += t15;
				t0 += t2;
				t4 += t6;
				t8 += t10;
				t12 += t14;
				t0 += t4;
				t8 += t12;
				t += t0+t8;
				
				j+= 16;
			}
			long dT16 = System.nanoTime() - T0;
			
			T1 += dT1;
			T16 += dT16;
			
			ss += (s-t);
		}
		
		double nano = 1e-9;
		double t1 = T1*nano/REP;
		double t16 = T16*nano/REP;
		
		out.printf("#DOTPRDUCT timing : t1 %.3g t16 %.3g (speedup %f)\t(ss %f)\n", t1, t16, t1/t16, ss);
		
		
		
		
		
		
		
		
	}
	
	public static void main(String[] args)
	{
		PrecisionTest T = new PrecisionTest();
		
		System.out.println("# Best way to calculate ln(1-b)=ln(1-exp(-x)): Math.log(-Math.expm1(-x)) when x is small (b near 1); Math.log1p(-Math.exp(-x)) when x is big (b near 0)");
		T.logTest(System.out);
		
		T.productTest(System.out);
		
		
	}

}
