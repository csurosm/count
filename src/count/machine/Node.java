package count.machine;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

/**
 * A leaf or a difference machine.
 * 
 * @author csuros
 *
 */
class Node extends Machine
{
	private static final double DEFAULT_LOSS=0.17;
	
	private Node(int A, int[] X)
	{
		this.A = A;
		this.X = X;
	}

	/**
	 * Instantiator used by the Insert machine: creates a leaf for a sequence
	 * 
	 * @param X
	 */
	Node(int[] X)
	{
		this (X[X.length-1]-1,X); // last one is the end character A+1
	}
	
	/**
	 * Instantiator used by the Difference Machine
	 */
	Node(int A)
	{
		this(A,null);
	}
	
	/**
	 * alphabet size
	 */
	final int A;
	/**
	 * leaf sequence
	 */
	final int[] X;

	
	/**
	 * upstream Insert machine; filled in by the instance Iu when created 
	 */
	Insert Iu = null;

	
	/**
	 * extinction probability  (used in Difference machine)
	 */
	double e=0.0;
	
	/**
	 * Original loss probability for the edge leading here.
	 */
	double orig_p=DEFAULT_LOSS;
	
	/**
	 * Survival loss probability (used in Difference machine)
	 */
	double p;
	
	/**
	 * Read likelihoods.
	 */
	double [/*read ell*/][/*sym*/] Read;

	@Override
	boolean isLeaf() { return true;}
	
	@Override 
	boolean isRoot() { return false;}
	
	
	/**
	 * Computes {@link #p} from {@link #e}; called by {@link #computeReadLikelihoods()}
	 */
	void computeSurvivalParameters()
	{
		double oq = Iu.orig_q;
		this.p = (orig_p*(1.-e)+e*(1.0-oq))/(1.0-oq*e);
		assert (e>0.0 || p==orig_p);  // no extinction at a leaf 
	}
	
	@Override
	void computeReadLikelihoods()
	{
		this.computeSurvivalParameters();
		
		Read = new double[X.length][A+2];
		Read[0][0] = 0.0; // log(1)
		for (int i=1; i<X.length; i++)
		{
			Read[i][0] = Double.NEGATIVE_INFINITY; // log(0), unless we find a residue at this position
			for (int sym=1; sym<=A; sym++)
				if (X[i]==sym)
				{
					Read[i][sym]=0.0; // =log(1)
					Read[i][0]=0.0;   // =log(1)
				}
				else
					Read[i][sym]=Double.NEGATIVE_INFINITY; // =log(0)
			if (X[i]==A+1)
				Read[i][A+1]=0.0; // log(1)
			else
				Read[i][A+1]=Double.NEGATIVE_INFINITY; // log(0)
		}
		if (PRINT_LIKELIHOODS)
			for (int ell=0; ell<Read.length; ell++)
			{
				System.out.println("#*D.cRL "+this+"\tR["+ell+"] = "+Arrays.toString(Read[ell]));
			}
	}
	
	@Override
	void computeWriteProbabilities()
	{
		throw new NoSuchMethodError("Nowhere to write.");
	}
	
	@Override
	public String toString()
	{
		return machine_id+" -> " + Arrays.toString(X);
	}

}
