package count.model;


import count.ds.AnnotatedTable;
import count.ds.ProfileTable;
import count.io.CommandLine;
import count.matek.Logarithms;
import count.matek.Functions;
import count.matek.Functions.RisingFactorial;

/**
 * 
 * @author csuros
 */
public class LikelihoodParametrized extends Likelihood 
{
	private static final int PARAMETER_EXTINCTION = PARAMETER_LENGTH;
	
	public LikelihoodParametrized(TreeWithLogisticParameters rates, ProfileTable table)
	{
		super(rates, table);
		this.num_nodes = rates.getTree().getNumNodes();
		
		this.logit_survival_parameters = new double[num_nodes][4];
		this.log_survival_parameters = new double[num_nodes][4];
		this.log_survival_complements = new double[num_nodes][4];
		
		for (int node=0; node<num_nodes; node++)
		{
			computeNodeParameters(node); // now we can 
		}
		this.parameter_factory = null;
	}
	
	public LikelihoodParametrized(TreeWithRates rates, ProfileTable table)
	{
		this(new TreeWithLogisticParameters(rates), table);
	}
	
	/**
	 * Instantiation with a linked factory using the same parameter set. 
	 * 
	 * @param that
	 * @param table
	 */
	public LikelihoodParametrized(LikelihoodParametrized that, ProfileTable table) 
	{
		super(that.rates,table);
		this.num_nodes = that.num_nodes;
		
		this.logit_survival_parameters = new double[num_nodes][];
		this.log_survival_parameters = new double[num_nodes][];
		this.log_survival_complements = new double[num_nodes][];
		this.parameter_factory = that;
		this.copyParametersFromFactory();
	}		
	
//	private TreeWithLogisticParameters getLogitModel() 
//	{
//		return (TreeWithLogisticParameters) rates;
//	}
	
	private final LikelihoodParametrized parameter_factory;
	private int num_nodes = 0; // set at instantiation; used as flag for super's instantiation
	
	private final double[][] logit_survival_parameters;
	private final double[][] log_survival_parameters;
	private final double[][] log_survival_complements;

	@Override
	public double getLogLossParameter(int node) { return log_survival_parameters[node][PARAMETER_LOSS];}
	@Override
	public double getLogDuplicationParameter(int node) { return log_survival_parameters[node][PARAMETER_DUPLICATION];}
	@Override
	public double getLogExtinction(int node) { return log_survival_parameters[node][PARAMETER_EXTINCTION];}
	@Override
	public double getLogLossComplement(int node) { return log_survival_complements[node][PARAMETER_LOSS];}
	@Override
	public double getLogDuplicationComplement(int node) { return log_survival_complements[node][PARAMETER_DUPLICATION];}
	@Override
	public double getLogExtinctionComplement(int node) { return log_survival_complements[node][PARAMETER_EXTINCTION];}
	
	@Override 
	public double getLossParameter(int node)
	{
		if (num_nodes==0) return super.getLossParameter(node);
		else return Math.exp(this.getLogLossParameter(node));
	}
	@Override 
	public double getDuplicationParameter(int node)
	{ 
		if (num_nodes==0) return super.getDuplicationParameter(node);
		else return Math.exp(this.getLogDuplicationParameter(node));
	}
	@Override 
	public double getLossParameterComplement(int node)
	{
		if (num_nodes==0) return super.getLossParameterComplement(node);
		else return Math.exp(this.getLogLossComplement(node));
	}
	@Override 
	public double getDuplicationParameterComplement(int node)
	{ 
		if (num_nodes==0) return super.getDuplicationParameterComplement(node);
		else return Math.exp(this.getLogDuplicationComplement(node));
	}
	@Override
	public double getExtinction(int node) 
	{
		if (num_nodes==0) return super.getExtinction(node);
		else return Math.exp(this.getLogExtinction(node));
	}
	@Override
	public double getExtinctionComplement(int node) 
	{ 
		if (num_nodes==0) return super.getExtinctionComplement(node);
		else return Math.exp(this.getLogExtinctionComplement(node));
	}
	
	private void setLogitSurvivalParameter(int node, int parameter_idx, double x)
	{
		double logp, log1_p;
		if (0.0<=x)
		{
			logp = Logarithms.logitToLogValue(x);
			log1_p = logp - x;
		} else
		{
			log1_p = Logarithms.logitToLogComplement(x);
			logp = log1_p + x;
		}
		logit_survival_parameters[node][parameter_idx] = x;
		log_survival_parameters[node][parameter_idx] = logp;
		log_survival_complements[node][parameter_idx] = log1_p;		
	}
	
	private void copyParametersFromFactory()
	{
		System.arraycopy(parameter_factory.logit_survival_parameters, 0, this.logit_survival_parameters, 0, num_nodes);
		System.arraycopy(parameter_factory.log_survival_parameters, 0, this.log_survival_parameters, 0, num_nodes);
		System.arraycopy(parameter_factory.log_survival_complements, 0, this.log_survival_complements, 0, num_nodes);
	}
	
	@Override
	public void computeNodeParameters(int node)
	{
		
		if (num_nodes == 0) // ==0 if super's instantiation calls here
		{
			super.computeNodeParameters(node);
		} else
		{
//			super.computeNodeParameters(node);
//			double old_logp = super.getLogLossParameter(node);
//			double old_log1_p = super.getLogLossComplement(node);
//			double oldp = super.getLossParameter(node);
//			double old_logq = super.getLogDuplicationParameter(node);
//			double old_log1_q = super.getLogDuplicationComplement(node);
//			double oldq = super.getDuplicationParameter(node);
			
			
//			System.out.println("#**LP.cNP "+node+"\tbefore "+toString(node)+"\t"+rates.toString(node));
			
			
			TreeWithLogisticParameters logit_model = (TreeWithLogisticParameters) rates;
			// calculate extinction
			double logit_e;
			
			int num_children = tree.getNumChildren(node);
			if (num_children==0) // leaf
			{
				logit_e = Double.NEGATIVE_INFINITY; // ==log(0/1)
			} else
			{
				logit_e = Double.POSITIVE_INFINITY; // == log(1/0)
				for (int c=0; c<num_children; c++)
				{
					int v = tree.getChild(node, c);
					double logit_p = logit_survival_parameters[v][PARAMETER_LOSS];
					logit_e = Logarithms.mulLogit(logit_e, logit_p);
				}
			}
			setLogitSurvivalParameter(node, PARAMETER_EXTINCTION, logit_e);
			// copy to super's 
			super.setExtinction(node, Math.exp(this.getLogExtinction(node)), Math.exp(this.getLogExtinctionComplement(node)));
			
			setNodeParametersLogit(node, logit_e);
			
//			// calculate duplication
//			double y = logit_model.logit_parameters[node][PARAMETER_DUPLICATION];
//			double log1_e = log_survival_complements[node][PARAMETER_EXTINCTION];
//			double logit_survival_q = y+log1_e;
//			setLogitSurvivalParameter(node, PARAMETER_DUPLICATION, logit_survival_q);
//			// copy to super's 
//			setDuplicationParameter(node, Math.exp(this.getLogDuplicationParameter(node)), Math.exp(this.getLogDuplicationComplement(node)));
//	
////			System.out.println("#**LP.cNP "+node
////					+"\tp "+getLossParameter(node)+"/"+oldp
////					+"\t("+(oldp-getLossParameter(node))+")"
////					+"\tq "+getDuplicationParameter(node)+"/"+oldq
////					+"\t("+(oldq-getDuplicationParameter(node))+")"
////					+"\tlq "+getLogDuplicationParameter(node)+"\tl1q "+getLogDuplicationComplement(node)+"\twas "+old_logq+"\t"+old_log1_q);
//			
//			
//			// calculate gain
//			if (y==Double.NEGATIVE_INFINITY) // duplication q==0.0
//			{
//				double r = rates.getGainParameter(node);
//				
//				
//				r *= Math.exp(log1_e);
//				this.gain_factorials[node] = null;
//	//			System.out.println("#*L.cP "+node+"\tq "+q+"\tpoisson "+r);
//	
//				setGainParameter(node, r);			
//			} else
//			{
//				double κ = rates.getGainParameter(node);
//				if (κ!=0.0)
//					this.gain_factorials[node] = new RisingFactorial(κ, 1+max_family_size[node]);
//				setGainParameter(node, κ);
//			}
//			
//			// calculate loss		
//			double x = logit_model.logit_parameters[node][PARAMETER_LOSS];
//			
//			double z = Logarithms.mulLogit(logit_e, -logit_survival_q); // epsilon*(1-q~)
//			
//			double logit_survival_p = -Logarithms.mulLogit(-x, -z);
//			setLogitSurvivalParameter(node, PARAMETER_LOSS, logit_survival_p);	
//			// copy to super's 
//			setLossParameter(node, Math.exp(this.getLogLossParameter(node)),Math.exp(this.getLogLossComplement(node)));
		}
	}	

	@Override 
	protected void setNodeParameters(int node, double epsi)
	{
		// ln (e/(1-e)) = log(e)-log1p(-e)
		if (num_nodes==0) 
		{
			// called from super's instantiation: ignore 
			// super.setNodeParameters(node, epsi);
		}
		else
		{
			double logit_e = Math.log(epsi)-Math.log1p(-epsi);
			setNodeParametersLogit(node, logit_e);
		}
	}
	
	protected void setNodeParametersLogit(int node, double logit_e)
	{
		TreeWithLogisticParameters logit_model = (TreeWithLogisticParameters) rates;
		// calculate duplication
		double y = logit_model.getLogitDuplicationParameter(node);//   logit_parameters[node][PARAMETER_DUPLICATION];
		double log_e, log1_e;
		if (0.0 <=logit_e)
		{
			log_e = Logarithms.logitToLogValue(logit_e);
			log1_e = log_e-logit_e;
		} else
		{
			log1_e = Logarithms.logitToLogComplement(logit_e);
			log_e = log1_e-logit_e;
		}
		double logit_survival_q = y+log1_e;
		setLogitSurvivalParameter(node, PARAMETER_DUPLICATION, logit_survival_q);
		// copy to super's 
		setDuplicationParameter(node, Math.exp(this.getLogDuplicationParameter(node)), Math.exp(this.getLogDuplicationComplement(node)));

//		System.out.println("#**LP.cNP "+node
//				+"\tp "+getLossParameter(node)+"/"+oldp
//				+"\t("+(oldp-getLossParameter(node))+")"
//				+"\tq "+getDuplicationParameter(node)+"/"+oldq
//				+"\t("+(oldq-getDuplicationParameter(node))+")"
//				+"\tlq "+getLogDuplicationParameter(node)+"\tl1q "+getLogDuplicationComplement(node)+"\twas "+old_logq+"\t"+old_log1_q);
		
		
		// calculate gain
		if (y==Double.NEGATIVE_INFINITY) // duplication q==0.0
		{
			double r = rates.getGainParameter(node);
			
			
			r *= Math.exp(log1_e);
			this.gain_factorials[node] = null;
//			System.out.println("#*L.cP "+node+"\tq "+q+"\tpoisson "+r);

			setGainParameter(node, r);			
		} else
		{
			
			double κ = rates.getGainParameter(node);
			double log_kappa = rates.getLogGainParameter(node);
			
			if (κ==0.0 )
			{
				double log_q = rates.getLogDuplicationParameter(node);
//				System.out.println("#**LP.sNP "+node+"\tkappa "+κ+"\tlogq "+log_q+"\t"+toString(node));

			}
//			if (κ!=0.0)
//				this.gain_factorials[node] = new RisingFactorial(κ, 1+max_family_size[node]);
			this.gain_factorials[node] = Functions.newRisingFactorialForLog(log_kappa, 1+max_family_size[node]); 
			setGainParameter(node, κ);
		}
		
		// calculate loss		
		double x = logit_model.getLogitLossParameter(node); // logit_parameters[node][PARAMETER_LOSS];
		
		double z = Logarithms.mulLogit(logit_e, -logit_survival_q); // epsilon*(1-q~)
		
		double logit_survival_p = -Logarithms.mulLogit(-x, -z);
		setLogitSurvivalParameter(node, PARAMETER_LOSS, logit_survival_p);	
		// copy to super's 
				
	}

	@Override
    protected Likelihood.Profile getProfileLikelihood(int family_idx)
    {
    	return this.new Profile(family_idx);
    }
	
	private class Profile extends Likelihood.Profile
	{
		Profile(int family_idx)
		{
			super(family_idx);
		}

		/**
		 * High-precision variant exploiting logistic scale
		 */
		@Override
		protected double[] computeNodeFromChildren(int node)
		{
			double[] C = null;
			int num_children = tree.getNumChildren(node);
			
			if (logit_survival_parameters[node][PARAMETER_EXTINCTION] == Double.POSITIVE_INFINITY)
			{ // sure extinction
				C = new double[1];
				C[0]=0.0; // log(1.0)
				for (int c=0; c<num_children; c++)
				{
					int junior = tree.getChild(node, c); // current sibling
					double[] K = getEdgeLikelihoods(junior); //  edge_likelihoods[junior];
					if (K.length>0)
						C[0] += K[0];
				}
			} else
			{
				double logit_eps = Double.POSITIVE_INFINITY; // extinction probability product across siblings
				for (int c=0; c<num_children; c++)
				{
					int junior = tree.getChild(node, c); // current child
					C = computeSiblingLogit(junior, C, logit_eps);
					
					
					double logit_pj = logit_survival_parameters[junior][PARAMETER_LOSS];
					logit_eps = Logarithms.mulLogit(logit_eps, logit_pj);
				} // loop across children 
				
			}
			
			
//			double[] oldC = super.computeNodeFromChildren(node);
//			System.out.println("#**LP.P.cNFC "+family_idx+"\t"+node+"\tC "+Arrays.toString(C)+"\toldC "+Arrays.toString(oldC));
			
			return C;
		}
		
		@Override
		protected double[] computeSiblingLogit(int junior, double[] C, double logit_eps)
		{
			double[] K =  getEdgeLikelihoods(junior);
			if (C==null || C.length==0)
			{
				// ambiguous on the right
				C = K.clone();
			} else if (K.length==0)
			{
				// ambiguous on the left 
				// nothing to do C is OK as is
			} else
			{
				int combined_family_size = (K.length-1)+(C.length-1);

				double loge, log1_e;
				if (0<=logit_eps)
				{
					loge = Logarithms.logitToLogValue(logit_eps);
					log1_e = loge-logit_eps;
				} else
				{
					log1_e = Logarithms.logitToLogComplement(logit_eps);
					loge = log1_e+logit_eps;
				}
				double logit_pj = logit_survival_parameters[junior][PARAMETER_LOSS];
				
				double logit_p2 = logit_pj + log1_e;
				double logp1,logp2;
				if (0<=logit_p2)
				{
					logp2 = Logarithms.logitToLogValue(logit_p2);
					logp1 = logp2-logit_p2;
				} else
				{
					logp1 = Logarithms.logitToLogComplement(logit_p2);
					logp2 = logp1 + logit_p2;
				}
				

//				{ // DEBUG
//					System.out.println("#**LP.P.cS "+family_idx+"/"+junior
//								+"\tlp1 "+logp1+"\tlp2 "+logp2
//								+"\tlogp "+getLogLossParameter(junior)+"\tlog1_p "+Math.log(one_minus_p)+"\tla "+loga+"\teps "+eps+"/loge "+loge+"\tlog1_e "+log1_e);
//				}
				
				
				
				double[] C2 = new double[combined_family_size+1]; // return value
				double[] terms = new double[K.length]; // auxiliary array to store summing terms 
				for (int ell=0; ell<=combined_family_size; ell++)
				{
					// calculate C(ell, t): t is at most as large as the max copy number in C[]  
					int t = Math.min(ell,  C.length);
					double x= // keeps the value of K2(ell,t+1) in the loop
						(t==C.length?Double.NEGATIVE_INFINITY:C[t]); // = C(ell,ell) while ell<C.length
					// when t starts with t==ell, C(t, t) = C[t] still, does not get touched
					// when t starts with t==C.length, all those values are log(0)
					int s= ell-t;
					while (t>0 && s<K.length-1)
					{
						--t;
						++s;
						x += log1_e;
						double y = C[t]+loge; // at this point C[t] = C(ell-1, t)
						x = C[t] = Logarithms.add(x, y); // recursion: C(ell,t) = (1-e)*C(ell,t+1)+e*C(ell-1,t)
					}
					
					// calculate C2[ell]
					double log_ellfact = factorials.factln(ell);

					assert (t==ell-s);
					int tmin = t; 
					while (s>=0 && t<C.length) // t<=ell by loop design 
					{
						double binom = log_ellfact - factorials.factln(s)-factorials.factln(t); // ell chose s 
						double slogp1 = s==0?0.0:s*logp1; // avoids 0*Double.NEGATIVE_INFINITY
						double tlogp2 = t==0?0.0:t*logp2; // avoids 0*Double.NEGATIVE_INFINITY
						terms[t-tmin] // fill lower cells
								= K[s] + C[t] + binom +  slogp1 +tlogp2;
//						if (Double.isNaN(terms[t-tmin]))
//						{
//							System.out.println("#**L.P.cS "+family_idx+"/"+junior+"\tell "+ell+"\tC2 "+Arrays.toString(C2)+"\tt "+t+"\ttmin "+tmin
//										+"\ts "+s+"\tK "+K[s]+"\tC "+C[t]
//										+"\tlp1 "+logp1+"\tlp2 "+logp2
//										+"\tp "+p+"\t1p "+one_minus_p+"\tla "+loga+"\teps "+eps+"/"+log1_e);
//						}
						--s;
						++t;
					} 
					C2[ell] = Logarithms.sum(terms, t-tmin);
//					if (Double.isNaN(C2[ell]))
//					{
//						System.out.println("#**L.P.cS "+family_idx+"/"+junior+"\tell "+ell+"\tC2 "+Arrays.toString(C2)+"\tt "+t+"\ttmin "+tmin+"// "+tree.toString(junior)+"\t// "+rates.toString(junior)
//								+"\tterms "+Arrays.toString(terms));
//					}
				} // for ell 
				C = C2;
			} // recurrence term for this child
			return C;
		}
	}
	

	public static void main(String[] args) throws Exception
	{
		CommandLine cli = new CommandLine(args, Likelihood.class);
        AnnotatedTable table = cli.getTable();
        TreeWithRates rates = cli.getRates();

        
		LikelihoodParametrized factory = new LikelihoodParametrized(rates, table);
		factory.mainmain(cli);
	}

}
