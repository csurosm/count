package count.model;

import java.util.Arrays;

public class FreeMixedModel implements MixedRateModel
{
	public FreeMixedModel(MixedRateModel mixed_model)
	{
		int nc = 0;
		cat_rates = new TreeWithRates[mixed_model.getNumClasses()];
		cat_probs = new double[cat_rates.length];
		for (int c=0; c<cat_rates.length; c++)
		{
			double pc = mixed_model.getClassProbability(c);
			if (pc!=0.0)
			{
				cat_rates[nc] = mixed_model.getClassModel(c);
				cat_probs[nc] = pc;
				nc++;
			}
		}
		cat_rates = Arrays.copyOf(cat_rates, nc);
		cat_probs = Arrays.copyOf(cat_probs, nc);
	}
	
	public FreeMixedModel(TreeWithRates[] cat_rates, double[] cat_weights)
	{
		this.cat_rates = cat_rates;
		this.cat_probs = new double[cat_rates.length];
		if (cat_weights == null)
		{
			double pc = 1.0/cat_rates.length;
			Arrays.fill(cat_probs, pc);
		} else
		{
			assert (cat_weights.length == cat_rates.length);
			double sum = 0.0;
			for (double w: cat_weights) sum += w;
			for (int c=0; c<cat_rates.length; c++)
			{
				cat_probs[c] = cat_weights[c]/sum;
			}
		}
	}

	private TreeWithRates[] cat_rates;
	private double[] cat_probs;
	@Override
	public int getNumClasses() 
	{
		return cat_rates.length;
	}
	@Override
	public double getClassProbability(int class_idx) 
	{
		return cat_probs[class_idx];
	}
	
	public double[] getClassProbabilities()
	{
		return cat_probs.clone();
	}

	@Override
	public TreeWithRates getClassModel(int class_idx) 
	{
		return cat_rates[class_idx];
	}
	
	public void setClassProbabilities(double[] p)
	{
		assert (p.length == cat_probs.length);
		System.arraycopy(p, 0, cat_probs, 0, getNumClasses());
	}
}
