package count.model;
/*
 * Copyright 2022 Mikl&oacute;s Cs&#369;r&ouml;s.
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

import count.matek.DiscreteDistribution;
import count.matek.NegativeBinomial;
import count.matek.PointDistribution;
import count.matek.Poisson;
import count.matek.ShiftedGeometric;

/**
 * Common interface for values associated with model parameters. 
 * 
 * @author Miklós Csűrös
 *
 */
public interface GLDParameters 
{
    public static final int PARAMETER_GAIN=0;
    public static final int PARAMETER_LOSS=1;
    public static final int PARAMETER_DUPLICATION=2;
	public static final int PARAMETER_LENGTH=3;
	public static final int PARAMETER_DUPLICATION_COMPLEMENT = 3;
	public static final int PARAMETER_LOSS_COMPLEMENT = 4;
	public static final int PARAMETER_RELATIVE_DUPLICATION = PARAMETER_LENGTH;
	
	
	public enum Type 
	{
		GAIN, LOSS, DUPLICATION, LENGTH
	};
	
	
	/**
	 * Gain parameter for transition probabilities. 
	 * 
	 * @param node edge or root prior
	 * @return a non-negative value
	 */
	public abstract double getGainParameter(int node);
	/**
	 * Loss parameter for transition probabilities. 
	 * 
	 * @param node node this edge is leading to (or root prior)
	 * @return a value between 0.0 and 1.0
	 */
	public abstract double getLossParameter(int node);
	/**
	 * Duplication parameter for transition probabilities. 
	 * 
	 * @param node node this edge is leading to (or root prior)
	 * @return a value between 0.0 and 1.0
	 */
	public abstract double getDuplicationParameter(int node);
	
	public default double getLossParameterComplement(int node)
	{
		return 1.0-getLossParameter(node);
	}
	
	public default double getDuplicationParameterComplement(int node)
	{
		return 1.0-getDuplicationParameter(node);
	}

	
	public default DiscreteDistribution getGainDistribution(int node)
	{
	    DiscreteDistribution D; // return value
	    
	    double q = getDuplicationParameter(node);
	    
	    if (q==0.0)
	    {
	    	// Poisson
	    	double r = getGainParameter(node);
	    	if(r==0.0)
	    		D=new PointDistribution(1.0);
	    	else
	    		D=new Poisson(r);
	    } else
	    {
	    	// Pólya
	    	double κ = getGainParameter(node);
	    	if (κ==0.0)
	    		D=new PointDistribution(1.0);
	    	else
	    		D=new NegativeBinomial(κ, q);
	    }
	    return D;
	}
	
	public default DiscreteDistribution getDuplicationDistribution(int node)
	{
		DiscreteDistribution D; // return value;
		double q = getDuplicationParameter(node);
		double p = getLossParameter(node);
		if (q==0.0)
			D = new PointDistribution(p);
		else 
			D = new ShiftedGeometric(p,q);
		return D;
    }
	
//	public static double getLogParameter(double param, double param_complement)
//	{
//		double log_param = (param<param_complement?Math.log(param):Math.log1p(-param_complement));
//		return log_param;
//		
//	}
	/**
	 * log(1-x) if we know both x and 1.0-x at double precision. Uses the smaller 
	 * value for calculation. 
	 * 
	 * @param param  variable  x at double precision
	 * @param param_complement 1.0-x at double precision 
	 * @return
	 */
	public static double log1m(double param, double param_complement)
	{
		double log_complement = (param<param_complement?Math.log1p(-param):Math.log(param_complement));
		return log_complement;
	}
	
	
	public static String paramName(int type)
	{
		if (type == PARAMETER_LENGTH) return "length";
		else if (type == PARAMETER_GAIN) return "gain";
		else if (type == PARAMETER_LOSS) return "loss";
		else if (type == PARAMETER_DUPLICATION) return "duplication";
		else return null;
	}
}
