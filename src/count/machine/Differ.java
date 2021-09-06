package count.machine;

import java.util.Arrays;
import java.util.List;

class Differ extends Node
{	
	Differ(Insert Iv, Insert Iw)
	{
		super(Iv.A);
		assert (Iw.A==Iv.A); // must have same alphabet at the children
		this.LOST = A+2;
		
		// make the mutator machines
    	this.Mv = new Mutate(this,Iv);
    	this.Mw = new Mutate(this,Iw);
	}

	/**
	 * Downstream left
	 */
	final Mutate Mv;
	/**
	 * Downstream right
	 */
	final Mutate Mw;
	
	/**
	 * Meta-character for not inheriting.
	 */
	final int LOST;
	
	

	@Override
	boolean isLeaf() { return false;}
	
	@Override 
	boolean isRoot() { return false;}

	@Override
	<T extends Machine> List<T> postOrder(List<T> machines, Class<T> machineClass) 
	{
		machines = Mv.postOrder(machines, machineClass);
		machines = Mw.postOrder(machines, machineClass);
		return super.postOrder(machines, machineClass);
	}
	
	
	/**
	 * Extinction probability {@link #e} and survival loss probability {@link #p}.
	 */
	@Override
	void computeSurvivalParameters()
	{
		Node v = Mv.Iv.Dv;
		Node w = Mw.Iv.Dv;

		double pv = v.p;
		double pw = w.p;
		this.e = pv*pw;
		double oq = Iu.orig_q;
		this.p = (orig_p*(1.-e)+e*(1.0-oq))/(1.0-oq*e);
		
	}
	
	
	
	
	double[/*read ell*/][/*s*/][/* sym */] FrameV; // FrameV is for use with v
	double[/*read ell*/][/*t*/][/* sym */] FrameW; // FrameW is for use with w
	double[/*read ell*/][/*write s in v*/][/* sym */] ReadWv;
	double[/*read ell*/][/*write t in w*/][/* sym */] ReadWw;  // we want both ReadWv and ReadWw for the posteriors 
	// inherited: double[/*read ell*/][/* sym */] Read;
	
	double[/*write s in v*/][/*read ell*/][/* sym */] WriteRv;
	double[/*write t in w*/][/*read ell*/][/* sym */] WriteRw;
	double[/*write s in v*/][/* sym */] Writev;
	double[/*write t in w*/][/* sym */] Writew;
			
	@Override
	void computeReadLikelihoods()
	{
		this.computeSurvivalParameters();
		
		int mv = Mv.Read.length-1;
		int mw = Mw.Read.length-1;
		int m = mv+mw-1; // maximum input length 
		FrameV = new double[m+1][Mv.Read.length][A+2];
		FrameW = new double[m+1][Mw.Read.length][A+2];
		// Frame[][]: 0=any, 1..A=residues, A+1=stop 
		
		ReadWv = new double[m+1][mv+1][LOST+A+1]; // LOST+0=any/0, LOST+A=A/0
		ReadWw = new double[m+1][mw+1][LOST+A+1];

		Read = new double[m+1][A+2];
		
		computeFrame(FrameV, Mw);
		computeReadChild(ReadWv, FrameV, Mv, Mw);
		
		computeFrame(FrameW, Mv);
		computeReadChild(ReadWw, FrameW, Mw, Mv);
		
		for (double[] row: Read)
			Arrays.fill(row, Double.NEGATIVE_INFINITY);
		Read[0][0]=0.0; // log(1)
		for (int ell=1; ell<=m; ell++)
		{
			for (int x=0; x<=A; x++)
			{
				Read[ell][x] = ReadWv[ell][0][x+LOST];
			}
			Read[ell][A+1] = ReadWv[ell][0][A+1]; 
			
			for (int s=1; s<=mv && s<=ell; s++)
			{
				for (int x=0; x<=A; x++)
				{
					double R = Logarithms.add(ReadWv[ell][s][x], ReadWv[ell][s][x+LOST]);
					Read[ell][x] = Logarithms.add(Read[ell][x], R);
				}
				Read[ell][A+1]=Logarithms.add(Read[ell][A+1], ReadWv[ell][s][A+1]);
			}
			if (PRINT_LIKELIHOODS)
				System.out.println("#*D.cRL "+this+"\tR["+ell+"] = "+Arrays.toString(Read[ell]));
		}
	}
	
	@Override
	void computeWriteProbabilities()
	{
		int mv = Mv.Read.length-1;
		int mw = Mw.Read.length-1;
		
		int m = mv+mw-1; // maximum input length
		WriteRv = new double[mv+1][m+1][LOST+A+1];
		Writev = new double[mv+1][A+2];
		WriteRw = new double[mw+1][m+1][LOST+A+1];
		Writew = new double[mw+1][A+2];
		
		computeWriteChild(WriteRv, Writev, FrameV, Mv, Mw);
		
		computeWriteChild(WriteRw, Writew, FrameW, Mw, Mv);
		
	}
	
	private void computeWriteChild(double[][][] WriteR, double[][] Write, double[][][] Frame, Mutate Mchild, Mutate Msibling)
	{
		double pv = Mchild.Iv.Dv.p;
		double pw = Msibling.Iv.Dv.p;
		final double log1_pv = Math.log((1.0-pv)/(1.0-pv*pw));
		final double log_pv = Math.log(pv*(1.0-pw)/(1.0-pv*pw));
		final double log1_pw = Math.log(1.0-pw); //Math.log1p(-pw);
		final double log_pw = Math.log(pw);

		final int mv = Write.length-1; // maximum output length
		assert (mv==WriteR.length-1);
		final int mu = Iu.Write.length-1; // maximum input length
		assert (mu==WriteR[0].length-1);
		
		// precomputed inheritance probabilities 
		double[][] inherit = new double[mu+1][A+2];
		for (int ell=1; ell<=mu; ell++)
		{
			for (int x=0; x<=A+1; x++)
				inherit[ell][x] = Iu.Write[ell][x]-Iu.Write[ell-1][0];
		}

		double[] Config = new double[mu+1]; // head configurations [s][ell] with current s
		Arrays.fill(Config, Double.NEGATIVE_INFINITY); // = log(0.0)

		for (double[][] row: WriteR)
			for (double[] col: row)
				Arrays.fill(col, Double.NEGATIVE_INFINITY); // log(0)
		for (double[] row: Write)
			Arrays.fill(row, Double.NEGATIVE_INFINITY); // log(0)

		
		Config[0] = 0.0; // start
		{ // s=0; 
			for (int ell=1; ell<=mu; ell++)
			{
				Config[ell] = log_pv+Config[ell-1]; // only lose 
			}
		}
		double prevConfig = Config[0]; // [s-1][ell-1]
				
		Write[0][0]=0.0;
		WriteR[0][0][0]=0.0; // init

		for (int s=1; s<=mv; s++)
		{
			for (int ell=s; ell<=mu; ell++)
			{
				double step_both = log1_pv + prevConfig; // [s-1][ell-1]
				double step_loss = log_pv + Config[ell-1]; // [s][ell-1]
				prevConfig = Config[ell]; // save [s-1][ell-1] for next iteration
				Config[ell] = Logarithms.add(step_both, step_loss);
				
				for (int x=0; x<=A; x++)
				{
					WriteR[s][ell][x+LOST]  = step_loss 
							+ inherit[ell][x] 
							+ Frame[ell][s][x]; // frame must match 
					double W = 
					WriteR[s][ell][x] = step_both + inherit[ell][x] 
							+ Logarithms.add(log1_pw+Frame[ell][s-1][x], log_pw);
					Write[s][x] = Logarithms.add(Write[s][x], W);
				}
				double stop =  
				WriteR[s][ell][A+1] = prevConfig + inherit[ell][A+1] 
						+ Frame[ell][s][A+1]; // frame must match
				Write[s][A+1] = Logarithms.add(Write[s][A+1], stop);
			}
			if (PRINT_LIKELIHOODS)
				System.out.println("#*D.cWC "+this+"->"+Mchild+"\tW["+s+"] = "+Arrays.toString(Write[s])
						+"\tlog1_pw "+log1_pw +"\tlogpw "+log_pw
						+"\tlog_pv "+log_pv+"\tlog1_pv "+log1_pv);
		}
	}
		
		
//		for (double[][] row: WriteR)
//			for (double[] col: row)
//				Arrays.fill(col, Double.NEGATIVE_INFINITY); // log(0)
//		for (double[] row: Write)
//			Arrays.fill(row, Double.NEGATIVE_INFINITY); // log(0)
//
//		Write[0][0]=0.0;
//		WriteR[0][0][0]=0.0; // start
//		for (int s=1; s<=mv; s++)
//		{
//			for (int ell=s; ell<=mu; ell++)
//			{
//				double slice_so_far = Logarithms.add(WriteR[s-1][ell-1][0], WriteR[s][ell-1][LOST]);
//				for (int x=0; x<=A; x++)
//				{
//					WriteR[s][ell][x+LOST]  = slice_so_far 
//							+ log_pv + inherit[ell][x] + Frame[ell][s][x]; // frame must match 
//					double W = 
//					WriteR[s][ell][x] = slice_so_far 
//							+ log1_pv + inherit[ell][x] + 
//							Logarithms.add(log1_pw+Frame[ell][s-1][x],log_pw);
//					Write[s][x] = Logarithms.add(Write[s][x], W);
//				}
//				double stop =  
//				WriteR[s][ell][A+1] = slice_so_far + inherit[ell][A+1] + Frame[ell][s][A+1]; // frame must match
//				Write[s][A+1] = Logarithms.add(Write[s][A+1], stop);
//			}
//			if (PRINT_LIKELIHOODS)
//				System.out.println("#*D.cWC "+this+"->"+Mchild+"\tW["+s+"] = "+Arrays.toString(Write[s])
//						+"\tlog1_pw "+log1_pw +"\tlogpw "+log_pw
//						+"\tlog_pv "+log_pv+"\tlog1_pv "+log1_pv);
//		}
//	}
	
	private void computeReadChild(double[][][] ReadW, double[][][] Frame, Mutate Mchild, Mutate Msibling)
	{
		double pv = Mchild.Iv.Dv.p;
		double pw = Msibling.Iv.Dv.p;
		
		final double log1_pv = Math.log((1.0-pv)/(1-pv*pw)); //Math.log1p(-pv)-Math.log1p(-pv*pw); // log ( (1-pv)/(1-pv*pw))
		final double log_pv =  Math.log(pv*(1.0-pw)/(1.0-pv*pw));  // Math.log(pv)+Math.log1p(-pw)+Math.log1p(-pv*pw); // log (pv*(1-pw)/(1-pvpw))
		
		final double log1_pw = Math.log(1.0-pw); //Math.log1p(-pw);
		final double log_pw = Math.log(pw);
		
		for (double[][] row: ReadW)
			for (double[] column: row)
				Arrays.fill(column, Double.NEGATIVE_INFINITY); // =log(0)
		
		int mv = Mchild.Read.length-1;
		int mw = Msibling.Read.length-1; // 0<= t=ell-s <= mw ; so s>= ell-mw
		int m = mv+mw-1;
		assert (ReadW.length == m+1);
		
		
		double[] Config = new double[mv+1]; // head configurations [ell][s] with current ell
		Arrays.fill(Config, Double.NEGATIVE_INFINITY); // = log(0.0)
		
		ReadW[0][0][0] = 0.0; // init
		Config[0] = 0.0; // start 
		for (int ell=1; ell<=m; ell++)
		{
			int s = 0; // only lose
			double lose = Config[s] + log_pv;
			double prevConfig = Config[s]; // Config[ell-1][s-1]
			Config[s]=lose;
			for (int x=0; x<=A; x++)
				ReadW[ell][s][x+LOST] = lose + Frame[ell][s][x];
			++s;
			do
			{
				assert (s<=mv && s<=ell);
				// stop 
				ReadW[ell][s][A+1] = prevConfig + Mchild.Read[s][A+1] + Frame[ell][s][A+1];

				lose = Config[s] + log_pv;
				double inherit = prevConfig + log1_pv;
				prevConfig = Config[s];
				Config[s] = Logarithms.add(inherit, lose);
	
				for (int x=0; x<=A; x++)
				{
					ReadW[ell][s][x+LOST] = lose + Frame[ell][s][x];
					ReadW[ell][s][x] = inherit + Mchild.Read[s][x]
							+ Logarithms.add(log1_pw+Frame[ell][s-1][x], log_pw);
				}
				++s;
			} while (s<=ell && s<=mv);
		}
		
	}
	
	/**
	 * Calculates the sibling's frame for connecting a child
	 * @param Frame will be filled with the frame likelihoods
	 * @param M sibling's mutator machine
	 */
	private void computeFrame(double[][][] Frame, Mutate Msib)
	{
		double pw = Msib.Iv.Dv.p;
		double log_pw = Math.log(pw);
		double log1_pw = Math.log(1.0-pw); //Math.log1p(-pw);
		
		int mw = Msib.Read.length-1; // max output pos at sibling
		int mu = Frame.length-1; // max for reading head pos ell
		int mv = mu-mw+1;
		assert (mv == Frame[0].length-1); // already allocated for 0<= s <= ell, mv
		
		for (double[][] row: Frame)
			for (double[] column: row)
				Arrays.fill(column, Double.NEGATIVE_INFINITY); // log(0)

		for (int ell=1; ell<=mu; ell++)
		{
			int s=0;
			if (ell<=mw)
			{
				for (int x=0; x<=A; x++)
					Frame[ell][s][x] = Msib.Read[ell][x];
				// Frame[ell][0][A+1] = log(0.0) already  
 			}
			++s;
			if (ell<=mw) 
			{
				Frame[ell][s][A+1] = Msib.Read[ell][A+1]; // stop with s=1 
			}
			do  // exit in the middle 
			{
				assert (s<=ell && s<=mv);

				for (int x=0; x<=A; x++)
				{
					double F_inherit = log1_pw + Frame[ell][s-1][x];
					double F_lose = log_pw + Frame[ell-1][s-1][x];
					
					double F=
					Frame[ell][s][x] = Logarithms.add(F_inherit, F_lose);
//					if (PRINT_LIKELIHOODS)
//						System.out.println("#*D.cF "+this+"\tF["+ell+"]["+s+"]["+x+"]=\t"+F+"\tnhrt "+F_inherit+"\tlose "+F_lose);
				}
					
				++s;
				// calculate stop frame
				if (s>ell || s>mv) break; // exit from loop here
				{
					double S_inherit = log1_pw + Frame[ell][s-1][A+1];
					double S_lose = log_pw + Frame[ell-1][s-1][A+1];
					double S = 
					Frame[ell][s][A+1] = Logarithms.add(S_inherit, S_lose);
//					if (PRINT_LIKELIHOODS)
//						System.out.println("#*D.cF "+this+"\tF["+ell+"]["+s+"]["+(A+1)+"]=\t"+S+"\tnhrt "+S_inherit+"\tlose "+S_lose);
				}
			} while (true);
		} // for ell
	}		
		
	@Override
	public String toString()
	{
		return machine_id+" -> "+Mv.machine_id+", "+Mw.machine_id;
	}


}
