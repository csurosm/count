 /* Copyright 2024 Mikl&oacute;s Cs&#369;r&ouml;s.
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

package count.matek;

import java.util.Arrays;

/**
 * Utility methods for computing with double[] arrays (distributions or counts).
 * 
 * @author csuros
 *
 */
public class ArraySum 
{
	private ArraySum() {} 
	
	/**
	 * Adds b to a, cell by cell, reusing a and/or b if possible. 
	 * If a is null, then b is returned.   
	 * 
	 * @param a expanded if necessary to match length of b
	 * @param b same length; untouched if a is not null
	 * @return cell-by-cell sums in array a (or b if a is null), or the resized version  of a
	 */
	public static double[] addCells(double[] a, double[] b, double bmul)
	{
		if (b==null) return a;
		if (a==null)
		{
			a = b.clone();
			for (int i=0; i<b.length; i++)
				b[i]*=bmul;
			return b;
		}
		
		int i=0;
		while (i<a.length && i<b.length)
		{
			a[i] = Math.fma(bmul, b[i], a[i]);// == bmul*b[i]+a[i]
			++i;
		}
		if (a.length<b.length)
			a = Arrays.copyOf(a, b.length); // filled with 0.0
		while (i<b.length)
		{
			a[i] = bmul*b[i];
			++i;
		}
//		for (int i=0; i<b.length; i++) 
//			a[i] = Math.fma(bmul, b[i], a[i]); // == bmul*b[i]+a[i]
		return a;
	}

	/**
	 * Adds b to a, cell by cell, reusing a and/or b if possible. 
	 * If a is null, then b is returned.   
	 * 
	 * @param a expanded if necessary to match length of b
	 * @param b untouched if a is not null
	 * @return cell-by-cell sums in array a (or b if a is null), or the resized version  of a
	 */
	public static double[] addCells(double[] a, double[] b)
	{
		if (a==null)
			return b;
		if (a.length<b.length)
			a = Arrays.copyOf(a, b.length);
		for (int i=0; i<b.length; i++)
			a[i] += b[i];
		return a;
	}	
	
	
	/**
	 * Array of tail sums: t[i]=sum_{j&gt;i} p[j]   
	 * 
	 * @param p may be null
	 * @return array of tail sums 
	 */
	public static double[] tail(double[] p)
	{
		if (p==null) return new double[1];
		int j = p.length;
		double[] t = new double[j];
		double x = 0.0;
		while (j>0)
		{
			--j;
			t[j] = x;
			x+=p[j];
		}
		return t;
	}
	
	/**
	 * Sum across the array. 
	 * 
	 * @param x
	 * @return sum of x[i]
	 */
	public static double sum(double[] x)
	{
		double sum=0.0;
		if (x!=null)
		{
			int i=x.length;
			while (0<i)
			{
				--i;
				double y = x[i];
				sum += y;
			}
		}
		return sum;
	}
	

}
