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


import java.util.Random;

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

import java.util.Random;

/** 
* Pseudorandom number generation with various distributions.
* 
* @since November 12, 2006, 9:04 PM
*/
public class PseudoRandom 
{
    /** Creates a new instance of PseudoRandom */
    public PseudoRandom(Random RND) 
    {
        this.RND=RND;
    }
    
    public PseudoRandom()
    {
        this(new Random());
    }
    
    public PseudoRandom(long seed)
    {
    	this(new Random(seed));
    }
    
    private Random RND;
    

    /**
     * A random value by Exponential(lambda). Uses the classic transformation method. 
     *
     * @param lambda the rate parameter (expected value 1/lambda)
     *
     */
    public final double nextExponential(double lambda)
    {
    	return nextExponential(lambda, nextUniform());
    }
    
    
    /**
     * A random value by Exponential(lambda). Uses the classic transformation method. 
     *
     * @param lambda the rate parameter (expected value 1/lambda)
     * @param RND random uniform[0,1] number generator 
     *
     */
    public static final double nextExponential(double lambda, Random RND)
    {
        return nextExponential(lambda, RND.nextDouble());        
    }

    
    private static final double nextExponential(double lambda, double U)
    {
    	return -Math.log(U)/lambda;   
    }
    
    /**
     * Geometric random variable Pr{X=k} = p*(1-p)^{k-1}, k=1,2,...
     */
    public final int nextGeometric(double p)
    {
        return nextGeometric(p, nextUniform());
    }
    
    /**
     * Geometric r.v. transformation.
     * 
     * Pr{x=k} = p*(1-p)^(k-1) and Pr{x&le;k} = 1-(1-p)^k.
     * In order to generate a r.v. with this distribution,
     * we use ceiling(log U/log (1-p)), where U is uniform(0,1).
     * 
     * @param p distribution parameter
     * @param U uniform[0,1] r.v.
     * @return
     */
    private static final int nextGeometric(double p, double U)
    {

    	double z = Math.log(1.-p);
        double R = Math.log(U)/z;
        int retval = (int)(R+1.); // poor man's ceiling 
        return retval;
    }


    /**
     * Geometric random variable Pr{X=k} = p*(1-p)^{k-1}, k=1,2,...
     */
    public static final int nextGeometric(double p, Random RND)
    {
    	return nextGeometric(p, RND.nextDouble());
    }

    
    /**
     * Uniform[0,1] random variable (using Random.nextDouble())
     */
    public final double nextUniform()
    {
        return RND.nextDouble();
    }
    

    /**
     * Uniform[0..range] integer-valued random variable (using Random.nextInt(n))
     */
    public final int nextUniform(int range)
    {
        return RND.nextInt(range);
    }
    
    
    public final long nextLong()
    {
    	return RND.nextLong();
    }
    
    /**
     * A random permutation of the first <var>n</var> natural numbers
     * @param n
     * @return
     */
    public int[] nextPermutation(int n)
    {
    	int[] nextPermutation = new int[n];
    	for (int i=0; i<n; i++) nextPermutation[i]=i;
    	
    	// as if selection sort
    	int tail = n;
    	int i=0; 
    	while (1<tail)
    	{
    		int j = i+RND.nextInt(tail);
    		// exchange
    		int xi = nextPermutation[i];
    		nextPermutation[i] = nextPermutation[j];
    		nextPermutation[j] = xi;
    		++i;
    		--tail;
    	}
    	return nextPermutation;
    	
    }
    
}
