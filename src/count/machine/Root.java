package count.machine;

import java.util.Arrays;

public class Root extends Mutate
{
	Root(Insert Iv)
	{
		super(null, Iv);
	}
	
	@Override
	void computeWriteProbabilities()
	{
		final int m = Iv.Read.length-1;
		Write = new double[m+1][A+2];
		WriteR = new double[m+1][A+2];
		for (double[] row: Write)
			Arrays.fill(row, Double.NEGATIVE_INFINITY) ; // =log(0)
		for (double[] row: WriteR)
			Arrays.fill(row,  Double.NEGATIVE_INFINITY);
		
		Write[0][0]    =
		WriteR[0][0]   = 0.0; // turn on machine
		Write[1][A+1]  = 
		WriteR[1][A+1] = 0.0; // insert coin to play  

		if (PRINT_LIKELIHOODS)
			for (int s=1; s<=m; s++)
				System.out.println("#*M.cWL "+this+"\tW["+s+"] = "+Arrays.toString(Write[s]));
	
	}
	@Override
	public String toString()
	{
		return machine_id +" -> "+Iv.machine_id+"\tRoot";
	}
	
	@Override
	void computeSurvivalParameter()
	{
		this.r = 1.0; // nothing is known above the root
	}
	
	/**
	 * Nothing to mutate above the root. 
	 */
	@Override
	double defaultMutationProbability()
	{
		return 0.0; 
	}
}
