package count.machine;

import java.util.Arrays;
import java.util.List;

class Mutate extends Machine
{
	private double DEFAULT_MUT = 0.1;
	
	/**
	 * Instantiation called by upstream Difference machine.
	 * @param Du upstream difference machine
	 * @param Iv downstream insert machine
	 */
	Mutate(Differ Du, Insert Iv)
	{
		this.A = Iv.A;
		this.REPLACED=A+2;
		assert (Du==null || Du.A==Iv.A); // they use the same alphabet
		this.Du = Du;
		this.Iv = Iv;
		// connect the downstream Insert machine to us
		Iv.Mv = this;
	
		this.mut = new double[A+1][A+1];
		double pmut = defaultMutationProbability();
		for (int z=1; z<=A; z++)
		{
			mut[z][z]=1.0;
			int y = 1;
			while (y<z)
			{
				mut[z][y++] = pmut;
				mut[z][z] -= pmut;
			}
			assert (z==y);
			y++; 
			while (y<=A)
			{
				mut[z][y++] =pmut;
				mut[z][z]-=pmut;
			}
		}
	}
	
	/**
	 * Default mutation probability to a different residue, by the Neyman (Jukes-Cantor)
	 * model ={@link #DEFAULT_MUT}/((@link #A}-1). 
	 * Serves as a hook for the extending class {@link Root}
	 * 
	 * @return {@link Double#POSITIVE_INFINITY} if alphabet size is 1
	 */
	double defaultMutationProbability()
	{
		return DEFAULT_MUT/(A-1.0);
	}
	
	/**
	 * Alphabet size
	 */
	final int A;
	
	/**
	 * Meta-character for replaced residue = {@link #A}+2
	 */
	final int REPLACED;
	
	/**
	 * upstream Difference Machine
	 */
	final Differ Du;
	
	/**
	 * downstream Insert Machine
	 */
	final Insert Iv;
	
	/**
	 * Mutation matrix
	 */
	double[][] mut;
	
	/** 
	 * Survival replacement parameter
	 */
	double r;
	
	// mutator machine
	double[/*read s*/][/*sym*/] Read;
	//
	double[/*read s*/][/*sym*/] ReadW;
	
	double[/*write s*/][/*sym*/] Write;

	double[/*write s*/][/*sym*/] WriteR;

	@Override
	boolean isLeaf() { return false;}
	
	@Override 
	boolean isRoot() { return false;}

	@Override
	<T extends Machine> List<T> postOrder(List<T> machines, Class<T> machineClass) 
	{
		machines = Iv.postOrder(machines, machineClass);
		if (machineClass.isInstance(this))
			machines.add((T) this); // meaningless type cast (T), but we are sure of the instance
		return machines;
	}
	
	/**
	 * Survival non-homologous replacement probability; called by {@link #computeReadLikelihoods()}
	 */
	void computeSurvivalParameter()
	{
		double oq = Iv.orig_q;
		double op = Iv.Dv.orig_p;
		
		double mu_t;
		if (oq==op)
		{
			mu_t = op/(1.0-op);
		} else 
		{
			mu_t = Math.log((1.0-oq)/(1.0-op))/(1.0-oq/op);
		}
		double p = Iv.Dv.p;
		double e = Iv.Dv.e;
		
		this.r = (1.0-p-Math.exp(-mu_t)*(1.0-e))/(1.0-p);
	}
	
	@Override
	void computeReadLikelihoods()
	{
		this.computeSurvivalParameter();
		
		final int m = Iv.Read.length-1;
		Read = new double[m+1][A+2];
		ReadW = new double[m+2][A+2];
		
		for (double[] row: Read)
			Arrays.fill(row, Double.NEGATIVE_INFINITY) ; // =log(0)
		for (double[] row: ReadW)
			Arrays.fill(row,  Double.NEGATIVE_INFINITY); 

		final double log_r = Math.log(r);
		final double log1_r = Math.log(1.0-r);
		

		Read[0][0]=0.0; // alive
		ReadW[0][0]=0.0; // no replacement

		// precompute the  replacement and substitution probabilities 
		final double[][] rM = new double[A+1][A+1];
		final double[] rPi = new double[A+1];
		
		for (int y=1; y<=A; y++)
		{
			rPi[y] = log_r + Math.log(Iv.pi[y]);
			for (int z=1; z<=A; z++)
				rM[z][y] = log1_r + Math.log(mut[z][y]);
		}
		
		for (int s=1; s<Read.length; s++)
		{
			Read[s][0] = Iv.Read[s][0];
			ReadW[s][0] = Read[s][0] + log1_r;
			// same replacement prob for all residues z=1..A
			
			double replace = Double.NEGATIVE_INFINITY;
			for (int y=1; y<=A; y++)
			{
				replace = Logarithms.add(replace, rPi[y]+Iv.Read[s][y]);
			}
			
			for (int z=1; z<=A; z++)
			{
				double mutate = rM[z][1]+Iv.Read[s][1];
				for (int y=2; y<=A; y++)
					mutate = Logarithms.add(mutate,rM[z][y]+Iv.Read[s][y]);
				
				
				Read[s][z] = Logarithms.add(mutate, replace);
				ReadW[s][z] = mutate;
			}
			// end-of-sequence z = A+1 is copied
			Read[s][A+1] = ReadW[s][A+1]=Iv.Read[s][A+1]; // stop
			
			if (PRINT_LIKELIHOODS)
				System.out.println("#*M.cRL "+this+"\tR["+s+"] = "+Arrays.toString(Read[s]));
		}
	}
	
	@Override
	void computeWriteProbabilities() 
	{
		final int m = Iv.Read.length-1;
		Write = new double[m+1][A+2];
		WriteR = new double[m+1][A+2];
		
		// find out if we are left or right child 
		double[][] Du_Write;
		if (Du.Mv == this)
		{
			Du_Write = Du.Writev;
		} else
		{
			Du_Write = Du.Writew;
		}
		
		double log_r = Math.log(r);
		double log1_r = Math.log(1.0-r);

		// precompute the  substitution and replacement probabilities 
		double[][] wM = new double[A+1][A+1];
		double[] wPi = new double[A+1];		
		for (int y=1; y<=A; y++)
		{
			wPi[y] = Math.log(Iv.pi[y])+log_r;
			for (int z=1; z<=A; z++)
				wM[z][y] = Math.log(mut[z][y])+log1_r;
		}

		for (double[] row: Write)
			Arrays.fill(row, Double.NEGATIVE_INFINITY) ; // =log(0)
		for (double[] row: WriteR)
			Arrays.fill(row,  Double.NEGATIVE_INFINITY);
		
		Write[0][0]  =
		WriteR[0][0] = 0.0; // init
		for (int s=1; s<=m; s++)
		{
			Write[s][0] = Du_Write[s][0]; 
			WriteR[s][0]  = Write[s][0]+ log1_r;
			for (int y=1; y<=A; y++)
			{
				double mutate = Double.NEGATIVE_INFINITY;
				for (int z=1; z<=A; z++)
				{
					mutate = Logarithms.add(mutate, Du_Write[s][z]+wM[z][y]);
				}
				double replace = Du_Write[s][0]+wPi[y];
				Write[s][y] = Logarithms.add(mutate, replace);
				WriteR[s][y] = mutate;
			}
			Write[s][A+1]=WriteR[s][A+1] = Du_Write[s][A+1]; // stop
			
			if (PRINT_LIKELIHOODS)
				System.out.println("#*M.cWL "+this+"\tW["+s+"] = "+Arrays.toString(Write[s]));
		}
	}
	
	@Override
	public String toString()
	{
		return machine_id +" -> "+Iv.machine_id+"\t <- "+Du.machine_id;
	}
}
