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
	
	public static void main(String[] args)
	{
		PrecisionTest T = new PrecisionTest();
		
		System.out.println("# Best way to calculate ln(1-b)=ln(1-exp(-x)): Math.log(-Math.expm1(-x)) when x is small (b near 1); Math.log1p(-Math.exp(-x)) when x is big (b near 0)");
		T.logTest(System.out);
		
		
	}

}
