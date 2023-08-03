/*
 * Copyright 2021 Mikl&oacute;s Cs&#369;r&ouml;s.
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

/**
 * Common interface for working with probabilistic distributions.
 *
 * @author csuros
 * @since March 20, 2008, 10:41 AM
 */
public interface DiscreteDistribution 
{
    /**
     * Returns the first few elements of the probability mass function
     *
     * @param n limit on how far the probabilities are computed
     * @return array of <var>n</var>+1 elements: [<var>p</var><sub>0</sub> <var>p</var><sub>1</sub> ... <var>p</var><sub>n</sub>]
     */
    
    public double[] getPointMassFunction(int n);
    
    
    /**
     * Returns the distribution's parameters: depends on the implementation what they mean 
     */
    public abstract double[] getParameters();
    
    
    /**
     *  Sets one of the distribution parameters
     *
     * @param parameter_idx index of the parameter (parameters are indexed the same order as with getParameters())
     * @param value New value of the parameter
     */
    public abstract void setParameter(int parameter_idx, double value);
    

    /**
     * Returns the number of parameters for defining the distribution: same as the length for the array returned by getParameters()
     */
    public abstract int getNumParameters();
    
}
