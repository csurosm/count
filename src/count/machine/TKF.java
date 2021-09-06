package count.machine;

import java.util.Map;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.List;
import java.io.PrintStream;

import count.ds.Phylogeny;
import count.ds.TreeTraversal;
import count.io.GeneralizedFileReader;
import count.util.Executable;
public class TKF 
{
	private TKF() 
	{
		this.sequence_headers = new ArrayList();
	}
	public TKF(String[] leaf_sequences)
	{
		this();
		initSequences(leaf_sequences);
	}
	/**
	 * Input sequences are stored at the leaves
	 */
	Insert[] leaves;
	
	/**
	 * Alphabet size. 
	 */
	private int A; 
	/**
	 * Inverted index from our encoding to characters in input sequences
	 */
	private String[] encoding=new String[6]; // init for DNA with integer codes 1..4; 5=end-of-sequece
	
	
	private void initSequences(final String[] leaf_sequences)
	{
		initSequences(leaf_sequences, true, 1);
	}

	/**
	 * Computes the integer encoding of 
	 * the input alphabet, and stores the 
	 * encoded sequences. Only  alphanumerical 
	 * characters are converted into 
	 * input symbols; spaces, dashes etc. are ignored.  
	 * 
	 * @param leaf_sequences 
	 * @param respect_case whether lower and upper case letters are the same
	 * @param alpha_block_size number of consecutive symbols constituting an input letter (eg, =3 for codons)
	 */
	private void initSequences(final String[] leaf_sequences, final boolean respect_case, final int alpha_block_size)
	{
		this.leaves = new Insert[leaf_sequences.length];
		int[][] X = new int[leaves.length][];
		// figure out input alphabet and sequence lengths
		this.A = 0;
		int leaf = 0;
		Map<String, Integer> input_alphabet = new HashMap<>();
		char[] input_symbol = new char[alpha_block_size];
		for (String seq: leaf_sequences)
		{
			ArrayList<Integer> Xleaf=new ArrayList<>(); // we don't know the input sequence length if gaps are allowed
			int n = seq.length(); 
			// for gapless imput and block length 1:
			// seq.charAt(0) -> position 1
			// seq.charAt(n-1) -> position n
			// end-of-sequence -> position n+1 = m
			int j = 0; // sequence position
			while (j<n)
			{
				int b = 0;
				while (b<alpha_block_size && j<n)
				{
					char c = seq.charAt(j);
					if (Character.isLetterOrDigit(c)) // increment b only if c is alphanumerical 
						input_symbol[b++] = respect_case?c:Character.toUpperCase(c);
					j++; 
				}
				if (b!=alpha_block_size)
					throw new IllegalArgumentException("Input sequence's ungapped length is not a multiple of alphabet block length "+alpha_block_size
								+" for sequence "+(1+leaf));
				String input = new String(input_symbol); // copy input symbols into a String
				int x; // our encoding
				if (input_alphabet.containsKey(input))
					x=input_alphabet.get(input);
				else
				{
					// new code
					A++;
					// store
					input_alphabet.put(input, A);
					if (encoding.length==A)
					{
						encoding = Arrays.copyOf(encoding, 2*A);
					}
					encoding[A] = input;
					System.out.println("#* Alphabet: "+input+"->"+A);
					x=A;
				}
				Xleaf.add(x);
			}
			// copy Xleaf into X[leaf] with position shift and Integer unboxing
			X[leaf] = new int[Xleaf.size()+2];
			int i=1; // position in X[leaf]
			for (int x: Xleaf)
				X[leaf][i++]=x;
			// next leaf
			++leaf;
		}
		System.out.println("#* Alphabet size: A="+A+"; start symbol=0, end symbol="+(A+1));

		// add end-of-sequence to each
		for (leaf=0; leaf<leaves.length; leaf++)
		{
			int seqlen = X[leaf].length-1;
			X[leaf][seqlen] = A+1; // end-of-sequence
			// cells 1..m used 
			leaves[leaf] = new Insert(X[leaf]);
		}
	}
	
	private ArrayList<String> sequence_headers;
	
	private void mainmain(String[] args) throws Exception
	{
		if (args.length<2)
			throw new IllegalArgumentException("Call as "+getClass().getName()+" tree sequences");
		
		int arg_idx=0;
		String tree_file = args[arg_idx++];
		String seq_file = args[arg_idx++];
		
        Phylogeny tree = count.io.NewickParser.readTree(new java.io.FileReader(tree_file));
	    java.io.BufferedReader reader    = new GeneralizedFileReader(seq_file);
	    
	    
        PrintStream out = System.out;
        
        out.println(Executable.getStandardHeader(this.getClass()));
        out.println(Executable.getStandardRuntimeInfo());
        out.println(Executable.getStandardHeader("Tree file:  "+tree_file));
        out.println(Executable.getStandardHeader("Table file: "+seq_file));
	    
		
	    StringBuilder current_sequence = null;
	    ArrayList<String> all_sequences = new ArrayList<>();
	    String line;
	    do
	    {
	    	line = reader.readLine();
	    	if (line != null)
	    	{
	    		if (line.startsWith(">"))
	    		{
	    			// new sequence
	    			if (current_sequence != null)
	    				all_sequences.add(current_sequence.toString());
	    			current_sequence=new StringBuilder();

	    			sequence_headers.add(line);
	    		} else
	    		{
	    			for (char c: line.toCharArray())
	    			{
	    				if (Character.isAlphabetic(c))
	    	    			current_sequence.append(c);
	    			}
	    		}
	    	}
	    } while (line != null);
	    if (sequence_headers.size()==all_sequences.size()+1)
	    	all_sequences.add(current_sequence.toString());
	    
	    
	    reader.close();
	    
	    initSequences(all_sequences.toArray(new String[0]));
	    
	    
	    
	    // build the tree
	    Insert[] machines = new Insert[tree.getNumNodes()];
	    
	    for (int node=0; node<tree.getNumNodes(); node++)
	    {
	    	if (tree.isLeaf(node))
	    		machines[node] = leaves[node];
	    	else
	    	{
	    		int v = tree.getChild(node, 0);
	    		int w = tree.getChild(node, 1);
	    		
	    		machines[node] = new Insert(machines[v], machines[w]);
	  	    }
	    }
	    Root root_machine = new Root(machines[machines.length-1]);
	    
	    List<Machine> all_machines = root_machine.postOrder(null, Machine.class);
	    for (Machine M: all_machines)
	    {
	    	out.println("#* Machine "+M);
	    	M.computeReadLikelihoods();
	    }
	    double LL = root_machine.Iv.Read[1][A+1]; // log-likelihood
	    
	    double L0 = root_machine.Iv.L0;
	    double corrLL = LL-Math.log(1.0-L0);
	    System.out.println("Log-likelihood "+LL+";\tL0 "+L0+";\tcorrected log-likelihood "+corrLL);
//	    Collections.reverse(all_machines);
//	    for (Machine M: all_machines)
//	    {
//	    	if (!M.isLeaf())
//	    	{
//		    	out.println("#* Machine "+M);
//		    	M.computeWriteProbabilities();
//	    	}
//	    }
	}
	
	
	public static void main(String[] args) throws Exception
	{
		TKF worker = new TKF();
		worker.mainmain(args);
	}
}
