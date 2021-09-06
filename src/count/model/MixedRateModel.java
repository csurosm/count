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
package count.model;

import count.ds.IndexedTree;

/**
 * Rate variation
 */
public interface MixedRateModel
{
    /**
     * Returns the number of discrete classes
     *
     * @return number of classes
     */
	public abstract int getNumClasses();
	
	
    /**
     * Prior probability for a rate class
     * @param class_idx 0..{@link #getNumClasses()}-1
     */
	public abstract double getClassProbability(int class_idx);
	
	/**
	 * Class-specific rate model. 
	 * 
	 * @param class_idx  0..{@link #getNumClasses()}-1
	 * @return
	 */
	public abstract RateModel.GLD getClassModel(int class_idx);
	
	
	

}
