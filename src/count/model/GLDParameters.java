package count.model;

import count.matek.DiscreteDistribution;
import count.matek.NegativeBinomial;
import count.matek.PointDistribution;
import count.matek.Poisson;
import count.matek.ShiftedGeometric;

public interface GLDParameters 
{
	
    public static final int PARAMETER_GAIN=0;
    public static final int PARAMETER_LOSS=1;
    public static final int PARAMETER_DUPLICATION=2;
	
	/**
	 * Gain parameter for transition probabilities. 
	 * 
	 * @param node_idx edge or root prior
	 * @return a non-negative value
	 */
	public abstract double getGainParameter(int node);
	/**
	 * Loss parameter for transition probabilities. 
	 * 
	 * @param node_idx node this edge is leading to (or root prior)
	 * @return a value between 0.0 and 1.0
	 */
	public abstract double getLossParameter(int node);
	/**
	 * Duplication parameter for transition probabilities. 
	 * 
	 * @param node_idx node this edge is leading to (or root prior)
	 * @return a value between 0.0 and 1.0
	 */
	public abstract double getDuplicationParameter(int node);

	
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
}
