package count.model;
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
 * Experimental code. 
 * 
 * @author csuros
 *
 */
public class ParametricBootstrap extends SimulatedEvolution 
{
	protected ParametricBootstrap(GammaInvariant rates_model, Random RND)
	{
		super(rates_model, RND);
	}
	
	private void oneShot(int nfam, int min_observed)
	{
		Table sim_tbl = table(nfam, min_observed);
		
		// optimize model on this data; compare with likelihood of true rates 
		// compare with original distribution parameters (for base)
		// compare posterior reconstructions: true history, estimated rates, true rates 
		
		TreeWithRates starting_rates = new TreeWithRates(getTree(), RND);
		
		for (int f=0; f<nfam; f++)
		{
			Table.ObservedProfile P = sim_tbl.getObservedProfile(f);
			
		}
	}
	
}
