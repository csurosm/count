package count.machine;

import java.util.Arrays;
import java.util.List;

import count.matek.Logarithms;

class Insert extends Machine
{
	private static final double DEFAULT_DUP = 0.1;

	Insert(int[] X)
	{
		this(new Node(X));
	}
	
	private Insert(Node Dv)
	{
		
		// copy alphabet size
		this.A = Dv.A; 
		this.INSERTED = A+2; // A+1=end=of-sequences
		// connect to child machine
		this.Dv = Dv;
		Dv.Iu = this;
		this.Mv = null; // no upstream connection yet

		
		// initialize residue frequencies
		this.pi = new double[A+1];
		for (int y=1; y<=A; y++)
			pi[y]=1.0/A;
		
	}
	
	Insert(Insert Iv, Insert Iw)
	{
		this(new Differ(Iv, Iw));
	}

	/**
	 * upstream Mutator; filled in by the instance of Mv when it is created.  
	 */
	Mutate Mv=null;
	
	/**
	 * Downstream difference or leaf
	 */
	final Node Dv;
	
	/**
	 * Alphabet size
	 */
	final int A;
	
	/**
	 * Meta-character for inserted slices = {@link #A}+2
	 */
	final int INSERTED;
	
	/**
	 * Non-conserved duplication probability. 
	 */
	double orig_q=DEFAULT_DUP;

	/**
	 * Symbol frequencies.
	 */
	double[] pi;
	
	/**
	 * Survival process duplication parameter
	 */
	double q;

	/**
	 * probability of empty profile
	 */
	double L0;
			
	/**
	 * Read likelihoods 
	 */
	double[/*read s*/][/*write ell*/][/*sym*/] ReadW; 
	double[/*read s*/][/*sym*/] Read;
	/**
	 * Write probabilities
	 */
	double[/*write ell*/][/*read s*/][/*sym*/] WriteR;
	double[/*write ell*/][/*sym*/] Write;
	

	@Override
	boolean isLeaf() { return false;}
	
	@Override 
	boolean isRoot() { return Mv==null;}

	
	/**
	 * Survival insert probability; called by {@link #computeReadLikelihoods()}
	 */
	void computeSurvivalParameters()
	{
		double e = Dv.e;
		this.q   = orig_q*(1.-e)/(1.-orig_q*e);
		
		this.L0 = 1.0-q;
		if (!Dv.isLeaf())
		{
			Differ D = (Differ) Dv;
			double L0v = D.Mv.Iv.L0;
			double L0w = D.Mw.Iv.L0;
			this.L0 *= L0v*L0w;
		}
	}
	
	
	@Override
	void computeReadLikelihoods()
	{
		this.computeSurvivalParameters();
		// allocate the likelihood arrays
		ReadW = new double[Dv.Read.length][Dv.Read.length][A+3]; // 1..A, A+1=end, A+2=INSERTED 
		Read = new double[Dv.Read.length][A+2];
		
		double logq = Math.log(q);
		double log1_q = Math.log(1.0-q);
		
		int m = ReadW.length-1; //maximum input length = maximum output length
		for (double[][] row: ReadW)
			for (double[] column: row)
				Arrays.fill(column, Double.NEGATIVE_INFINITY); // = log(0)
		for (double[] column: Read)
			Arrays.fill(column, Double.NEGATIVE_INFINITY);

		// precomputed insertprobs
		double[/*ell*/] insert = new double[m+1];
		for (int ell=1; ell<=m; ell++)
		{
			insert[ell] = Double.NEGATIVE_INFINITY; // log(0.0) 
			for (int x=1; x<=A; x++)
				insert[ell] = Logarithms.add(insert[ell],  Math.log(pi[x]) + Dv.Read[ell][x]-Dv.Read[ell-1][0]);
			insert[ell] += logq; 
		}
		
		// s=0: initialization
		ReadW[0][0][0]=0.0; // no insertions before position 0
		Read[0][0]=0.0;
		for (int s=1; s<=m; s++)
		{
			for (int ell=s; ell<=m; ell++)
			{
				// read likelihoods
				double slice_so_far = Logarithms.add(ReadW[s-1][ell-1][0], ReadW[s][ell-1][INSERTED]);
				for (int x=0; x<=A+1; x++)
				{
					double copy = log1_q + Dv.Read[ell][x]-Dv.Read[ell-1][0];
					double R =
					ReadW[s][ell][x] = slice_so_far + copy;
					Read[s][x] = Logarithms.add(Read[s][x], R);
				}
				ReadW[s][ell][INSERTED] = slice_so_far + insert[ell];
			} // ell
			if (PRINT_LIKELIHOODS)
				System.out.println("#*I.cRL "+this+"\tR["+s+"] = "+Arrays.toString(Read[s]));
		}
	}
	
	@Override
	void computeWriteProbabilities()
	{
		final int m = Dv.Read.length-1;
		// allocate the likelihood array
		WriteR = new double[m+1][m+1][2*A+3];
		// 0..A, A+1=end, A+2=inserted/0, ... 2A+2=inserted/A
		Write = new double[m+1][A+2];
		
		final double logq = Math.log(q);
		final double log1_q = Math.log(1.0-q);
		
		// precomputed insert probabilities
		final double[] logqPi = new double[A+1];
		logqPi[0] = logq;
		for (int x=1; x<=A; x++)
		{
			logqPi[x] = logq+Math.log(pi[x]);
		}
		
		// precomputed copy probabilities 
		final double[][] copy  = new double[m+1][A+2];
		for (int s=1; s<=m; s++)
		{
			for (int y=0; y<=A+1; y++)
			{
				if (Mv.Write[s][y]==Double.NEGATIVE_INFINITY) // log(0.0)
					copy[s][y] = Double.NEGATIVE_INFINITY; // before we do log(0.0)-log(0.0)=NaN
				else	
					copy[s][y] = log1_q + Mv.Write[s][y]-Mv.Write[s-1][0];
			}
		}
		
		for (double[][] row: WriteR)
			for (double[] column: row)
				Arrays.fill(column, Double.NEGATIVE_INFINITY); // = log(0)
		for (double[] column: Write)
			Arrays.fill(column, Double.NEGATIVE_INFINITY);

		WriteR[0][0][0] = 0.0; // start
		Write[0][0] = 0.0; // log(1)
		for (int ell=1; ell<=m; ell++)
		{
			for (int s=1; s<=ell; s++)
			{
				final double slice_so_far = Logarithms.add(WriteR[ell-1][s-1][0], WriteR[ell-1][s][INSERTED]);
				{ // x=0
					double W = WriteR[ell][s][INSERTED] = slice_so_far + logq; 
					Write[ell][0] = Logarithms.add(Write[ell][0],W) ;  
				}
				for (int x=1; x<=A; x++)
				{
					double W = 
					WriteR[ell][s][x+INSERTED] = slice_so_far + logqPi[x];
					Write[ell][x] = Logarithms.add(Write[ell][x],W);
				}
				for (int x=0; x<=A+1; x++)
				{
					double W = 
					WriteR[ell][s][x] = slice_so_far + copy[s][x];
					Write[ell][x] = Logarithms.add(Write[ell][x], W);
				}
			} // for s
			if (PRINT_LIKELIHOODS)
				System.out.println("#*I.cWL "+this+"\tW["+ell+"] = "+Arrays.toString(Write[ell]));
		}
	}
	
	@Override
	public String toString()
	{
		return machine_id + " -> "+Dv.machine_id;
	}

	@Override
	<T extends Machine> List<T> postOrder(List<T> machines, Class<T> machineClass) 
	{
		machines = Dv.postOrder(machines, machineClass);
		if (machineClass.isInstance(this))
			machines.add((T) this); // meaningless type cast (T), but we are sure of the instance
		return machines;
	}
	
}
