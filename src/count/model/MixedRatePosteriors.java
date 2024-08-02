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



import java.io.PrintStream;
import java.util.Arrays;
import java.util.function.IntFunction;

import count.ds.IndexedTree;
import count.ds.ProfileTable;
import count.io.CommandLine;
import count.matek.Logarithms;
import count.model.Posteriors.FamilyEvent;

import static count.io.CommandLine.OPT_ANCESTRAL;
import static count.io.CommandLine.OPT_HISTORY;


/**
 * Algorithms for computing posteriors in a rate-variation model.
 * 
 * @author csuros
 *
 */
public class MixedRatePosteriors 
{
	public MixedRatePosteriors(MixedRateModel mixed_model, ProfileTable table)
	{
		this.mixed_model = mixed_model;
		this.table = table;
		this.class_posteriors = new Posteriors[mixed_model.getNumClasses()];
		this.class_index = new int[class_posteriors.length];
		Arrays.fill(class_index,-1);
		this.log_class_prob = new double[class_posteriors.length]; // precomputed
		this.num_classes = 0;
		this.min_copies = Integer.min(table.minCopies(), 2);
	}
	private final MixedRateModel mixed_model;
	private final ProfileTable table;
	
	private int min_copies;
	
	private final Posteriors[] class_posteriors; 
	/**
	 * Mapping from our class indices to {@link #mixed_model} class indices.
	 */
	private final int[] class_index;
	private final double[] log_class_prob;
	
	public MixedRateModel getRateModel()
	{
		return mixed_model;
	}
	
	public ProfileTable getTable()
	{
		return table;
	}
	
	public Profile getProfile(int family_idx)
	{
		Profile P=null;
		P = new Profile(family_idx);
		
		return P;
	}
	
	private int threshold_width_absolute = Integer.MAX_VALUE;
	private double threshold_width_relative = Double.POSITIVE_INFINITY;
	/**
	 * Sets up truncated likelihood calculations: 
	 * at ancestors, the maximum copy number assumed
	 * is given by the minimum between the absolute
	 * threshold and the product of the relative threshold 
	 * times maximum copy number observed. 
	 * 
	 * @param absolute
	 * @param relative
	 */
	public void setCalculationWidthThresholds(int absolute, double relative)
	{
		this.threshold_width_absolute = absolute;
		this.threshold_width_relative = relative;
	}
	
	
	public int getCalculationWidthAbsolute()
	{
		return this.threshold_width_absolute; 
	}
	
	public double getCalculationWidthRelative()
	{
		return this.threshold_width_relative;
	}
	
	private double ancestor_deviation = 0.0; // set to 0 for don't care
	
	
	public void setAncestorWidthThreshold(double deviation)
	{
		this.ancestor_deviation = deviation;
	}
	
	/**
	 * Non-0 classes. 
	 */
	private int num_classes;
	
	/**
	 * Calculates class probabilities, indexing and instantiates 
	 * the {@link #class_posteriors} entries; 
	 * call if the underlying model changes.
	 */
	public void computeClasses()
	{
		Arrays.fill(class_index,-1);
		Arrays.fill(log_class_prob,Double.NEGATIVE_INFINITY);
		num_classes=0;
		for (int c = 0; c<class_posteriors.length; c++)
		{
			double pc = mixed_model.getClassProbability(c);
			if (pc != 0.0)
			{
				TreeWithRates rates =mixed_model.getClassModel(c);
				Posteriors post = new Posteriors(rates, table);
				post.setCalculationWidthThresholds(threshold_width_absolute, threshold_width_relative);
				post.setAncestorWidthThreshold(ancestor_deviation);
				class_posteriors[num_classes] = post;
				
				class_index[num_classes] = c;
				log_class_prob[num_classes] = Math.log(pc);
				num_classes++;
			}
		}
	}
	
	public int getClassCount()
	{
		if (num_classes==0) computeClasses();
		return num_classes;
	}
	
	public double getEmptyLL()
	{
		if (num_classes==0) computeClasses();
		double[] empty = new double[num_classes];	
		for (int c=0; c<num_classes; c++)
			empty[c] = class_posteriors[c].factory.getEmptyLL();
		return calculateLL(empty);
	}
	
	public double getSingletonLL()
	{
		if (num_classes==0) computeClasses();
		double[] single = new double[num_classes];
		for (int c=0; c<num_classes; c++)
			single[c] = class_posteriors[c].factory.getSingletonLL();
		return calculateLL(single);
	}
	
	/**
	 * Guides how to correct the likelihood in {@link Profile#getCorrectedLL()}
	 * 
	 * @param min_copies
	 */
	public void setMinimumObserved(int min_copies)
	{
		if (min_copies<0 || min_copies>2) throw new UnsupportedOperationException("Min copies 0,1,2 are implemented");
		assert (min_copies<=table.minCopies());
		this.min_copies = min_copies;
	}
	
	public int getMinimumObserved() { return  min_copies;}
	
	public double getUnobservedLL()
	{
		double unobservedLL = Double.NEGATIVE_INFINITY;
		if (min_copies>0)
		{
			unobservedLL = getEmptyLL();
			if (min_copies>1)
			{
				assert (min_copies==2);
				unobservedLL = Logarithms.add(unobservedLL, getSingletonLL());
			}
		}
		return unobservedLL;
	}
	
	/**
	 * Array for empty and possibly singleton profiles 
	 * 
	 * @return
	 */
	public Profile[] getUnobservedProfiles()
	{
		if (num_classes==0) computeClasses();
		
		Profile[] unobserved;		
		if (min_copies==0)
		{
			unobserved = new Profile[0];
		} else
		{			
			MixedRateModel rates_model = getRateModel();
			IndexedTree phylo = getTree();
			
			MixedRatePosteriors empty_posteriors = 
					new MixedRatePosteriors(rates_model, 
							ProfileTable.emptyProfile(phylo));
			empty_posteriors.setCalculationWidthThresholds(threshold_width_absolute, threshold_width_relative);
			empty_posteriors.setAncestorWidthThreshold(ancestor_deviation);
			empty_posteriors.computeClasses();
			if (min_copies == 1)
			{
				unobserved = new MixedRatePosteriors.Profile[1];
				unobserved[0] = empty_posteriors.getProfile(0);
				
			} else
			{
				MixedRatePosteriors singleton_posteriors = new MixedRatePosteriors(rates_model, 
						ProfileTable.singletonTable(phylo));
				singleton_posteriors.setCalculationWidthThresholds(threshold_width_absolute, threshold_width_relative);
				singleton_posteriors.setAncestorWidthThreshold(ancestor_deviation);
				singleton_posteriors.computeClasses();
				
				unobserved = new MixedRatePosteriors.Profile[1+phylo.getNumLeaves()];
				unobserved[0] = empty_posteriors.getProfile(0);
				for (int leaf=0; leaf<phylo.getNumLeaves(); leaf++)
				{
					unobserved[1+leaf] = singleton_posteriors.getProfile(leaf);
				}
			}
		}
		
		return unobserved;
	}
	
	private double calculateLL(double[] classLL)
	{
		double[] sum=new double[classLL.length];
		for (int c=0; c<num_classes; c++)
		{
			sum[c] = classLL[c]+log_class_prob[c];
		}
		double LL=Logarithms.sum(sum, num_classes);
		return LL;
	}
	
	public IndexedTree getTree()
	{
		return class_posteriors[0].factory.tree;
	}


	public class Profile
	{
		Profile(int family_idx)
		{
			if (num_classes==0) computeClasses();
//			System.out.println("#**MRP.P() "+family_idx);
			post = new Posteriors.Profile[num_classes];
			for (int c=0; c<num_classes; c++)
			{
				post[c] = class_posteriors[c].getPosteriors(family_idx);
				post[c].computeLikelihoods();
			}
			class_LL = new double[num_classes];
			pc = new double[num_classes];
//			System.out.println("#**MRP.P() "+family_idx+"\talloc");
			initDataStructures();
//			System.out.println("#**MRP.P() "+family_idx+"\tdone");
		}
		private final Posteriors.Profile[] post;
		private double LL;
		private double[] class_LL;
		/**
		 * Posterior class probabilities
		 */
		private double[] pc;
//		/**
//		 * Node survival means 
//		 */
//		private double[] expN;
//		/**
//		 * Edge survival means 
//		 */
//		private double[] expS;
//		
		private double[][] events;
		
		private void initDataStructures()
		{
			for (int c=0; c<num_classes; c++)
			{
				Likelihood.Profile cLik = post[c].inside;
				double ll = cLik.getLogLikelihood();
				class_LL[c]=ll;
			}
			LL = calculateLL(class_LL);

			for (int c=0; c<num_classes; c++)
			{
				pc[c]=Math.exp(class_LL[c]+log_class_prob[c]-LL);
			}
			int num_nodes = getTree().getNumNodes();
			events = new double[num_nodes][];
			
//			System.out.println("#**MRP.P.iDS classes ");
			for (int node=0; node<num_nodes; node++)
			{
				events[node]=new double[Posteriors.FamilyEvent.values().length];
				for (int c=0; c<num_classes; c++)
				{
					Posteriors.Profile P = post[c];
					double[] ev = P.getFamilyEventPosteriors(node);
					assert (ev.length==events[node].length);
					for (int j=0; j<ev.length; j++)
					{
						events[node][j] += pc[c]*ev[j];
					}
				}
			}
//			System.out.println("#**MRP.P.iDS events ");
		}
		
		
		public double getLL() { return LL;}
		
		
		
		public double getCorrectedLL()
		{
			double not_observed = getUnobservedLL();
			double pobs = -Math.expm1(not_observed); // 1-e^loglik
			return LL-Math.log(pobs);
		}
		
		/**
		 * Posterior class probabilities
		 * 
		 * @param class_idx
		 * @return
		 */
		public double getClassProbability(int class_idx)
		{
			return pc[class_idx];
		}
		
		private double sumClasses(IntFunction<Double> class_stats)
		{
			double sum = 0.0;
			for (int c=0; c<num_classes; c++)
			{
				
				double f = class_stats.apply(c);
				sum += pc[c]*f;
			}
			return sum;
		}
		
		public double getNodeMean(int node)
		{
			return sumClasses(c->post[c].getNodeMean(node));
		}
		
		public double getEdgeMean(int node)
		{
			return sumClasses(c->post[c].getEdgeMean(node));
		}
		
		/**
		 * Expected number of non-inherited ancestral copies at parent.
		 * 
		 * @param node
		 * @return
		 */
		public double getDeathMean(int node)
		{
			IndexedTree tree = mixed_model.getBaseModel().getTree();
			if (tree.isRoot(node))
				return 0.0;
			else
			{
				int parent = tree.getParent(node);
				double N = getNodeMean(parent);
				double S = getEdgeMean(node);
				return Double.max(0.0, N-S);
			}
		}		
		
		/**
		 * Expected number of new ancestral copies at node (by duplication and gain)
		 * 
		 * @param node
		 * @return
		 */
		public double getBirthMean(int node)
		{
			double S = getEdgeMean(node);
			double N = getNodeMean(node);
			
			return Double.max(0.0,N-S);
		}		
		
//		public double getConservedGain(int node)
//		{
//			return getNodeMean(node)-getEdgeMean(node);
//		}
//		
//		public double getConservedLoss(int node)
//		{
//			IndexedTree tree = getTree();
//			if (tree.isRoot(node)) return 0.0;
//			int parent = tree.getParent(node);
//			return getNodeMean(parent)-getEdgeMean(node);
//		}
		
		/**
		 * Posterior probabilities for copy numbers 0,1,...
		 * 
		 * @param node
		 * @return
		 */
		public double[] getNodeAncestor(int node)
		{
			
			double[] pN=new double[2];
			for (int c=0; c<num_classes; c++)
			{
				Posteriors.Profile P = post[c];
				double[] cN = P.getNodeAncestorPosteriors(node);
				if (pN.length<cN.length)
					pN = Arrays.copyOf(pN, cN.length);
				for (int n=0; n<cN.length; n++)
				{
					pN[n] += pc[c]*cN[n];
				}
			}
			return pN;
		}
		
		/** 
		 * Posterior probability for family-level events.
		 * 
		 * @param node
		 * @param event_type gain, loss, expansion or contraction
		 * @return
		 */
		public double getFamilyEvent(int node, Posteriors.FamilyEvent event_type)
		{
			final int event_idx = event_type.ordinal();
			return sumClasses(c->events[node][event_idx]);
		}
		
		
		public void reportPosteriors(PrintStream out)
		{
			int num_nodes = getTree().getNumNodes();
			IndexedTree phylo = getTree();
			int family_idx = post[0].inside.family_idx;
			for (int node=0; node<num_nodes; node++)
			{
				out.printf("%d\t%s", family_idx,phylo.getIdent(node));
				out.printf("\t%g", getNodeMean(node));
				out.printf("\t%g", getEdgeMean(node));
				double[] pN = getNodeAncestor(node);
				for (int n=0; n<pN.length; n++)
				{
					out.printf("\tp%d:%g", n, pN[n]);
				}
				for (FamilyEvent e: FamilyEvent.values())
				{
					out.printf("\t%s:%g", e, getFamilyEvent(node, e));
				}
				out.println();
			}
		}
		
	}
	
	
	private void printPosteriors(PrintStream out, boolean want_families)
	{
		// init unobserved[]
		Profile[] unobserved = getUnobservedProfiles();
		// compute correction factors
		double L0 = Double.NEGATIVE_INFINITY; // log(0)
		for (int i=0; i<unobserved.length; i++)
		{
			L0 = Logarithms.add(L0, unobserved[i].getLL());
		}
		// double p_unobs = Math.exp(L0);
		double p_obs = -Math.expm1(L0); // 1-exp(L0)
		double Lobs = Math.log(p_obs);
		
		// expected number of times unobserved profile comes up for each observed one
		double[] factor = new double[unobserved.length];
		for (int i=0; i<unobserved.length; i++)
		{
			double Ldiff = unobserved[i].getLL()-Lobs;
			factor[i] = Math.exp(Ldiff);
		}			
		
		IndexedTree phylo = getTree();
		int n = phylo.getNumNodes();
		// unobserved corrections
		double[] delta_family_present = new double[n];
		double[] delta_family_multi = new double[n];
		double[] delta_family_gain = new double[n];
		double[] delta_family_loss = new double[n];
		double[] delta_family_expand = new double[n];
		double[] delta_family_contract = new double[n];
		double[] delta_copies_surviving = new double[n];
		double[] delta_copies_true = new double[n];
		
		
		for (int f=0; f<unobserved.length; f++)
		{
			Profile P = unobserved[f];
			for (int node=0; node<n; node++)
			{
				double[] pN = P.getNodeAncestor(node);
				double fp = 1.0-pN[0];
				double fm = 1.0-(pN[0]+pN[1]);
				double fgain = P.getFamilyEvent(node, FamilyEvent.GAIN);
				double floss = P.getFamilyEvent(node, FamilyEvent.LOSS);
				double fexpand = P.getFamilyEvent(node, FamilyEvent.EXPAND);
				double fcontract = P.getFamilyEvent(node,  FamilyEvent.CONTRACT);
				double fsurviving = P.getNodeMean(node);
				double truecopies = mean(pN);
				

				delta_family_present[node] += factor[f]*fp;
				delta_family_multi[node] += factor[f]*fm;
				delta_family_gain[node] += factor[f]*fgain;
				delta_family_loss[node] += factor[f]*floss;
				delta_family_expand[node] += factor[f]*fexpand;
				delta_family_contract[node] += factor[f]*fcontract;
				delta_copies_surviving[node] += factor[f]*fsurviving;
				delta_copies_true[node] += factor[f]*truecopies;
			}
		}

		String prefix_ancestral = OPT_ANCESTRAL.toUpperCase();
		String prefix_families = OPT_HISTORY.toUpperCase();
		
		if (want_families)
		{
			out.println(prefix_families+"\tfamily\tnodeidx\tpresence:p\tmulti:m\tcorrected:p.\tgain:g\tloss:l\texpand:++\tcontract:--\tscopies:n\ttcopies:n.\tvarp:vp");
		}
		
		// lineage totals
		double[] lineage_family_present = new double[n];
		double[] lineage_family_present_corrected = new double[n];
		double[] lineage_family_multi = new double[n];
		double[] lineage_family_gain = new double[n];
		double[] lineage_family_loss = new double[n];
		double[] lineage_family_expand = new double[n];
		double[] lineage_family_contract = new double[n];
		double[] lineage_copies_surviving = new double[n];
		double[] lineage_copies_true = new double[n];
		double[] var_family_present = new double[n];
		
		for (int f=0; f<table.getFamilyCount(); f++)
		{
			Profile P = getProfile(f);
			for (int node=0; node<n; node++)
			{
				double[] pN = P.getNodeAncestor(node);
				
				double sum1;
				double sum2=0.0;
				for (int i=2; i<pN.length; i++)
					sum2+=pN[i];
				sum1 = sum2+pN[1];
				
				double fp = Double.max(1.0-pN[0],sum1);
				double fm = Double.max(1.0-(pN[0]+pN[1]),sum2);
				double fgain = P.getFamilyEvent(node, FamilyEvent.GAIN);
				double floss = P.getFamilyEvent(node, FamilyEvent.LOSS);
				double fexpand = P.getFamilyEvent(node, FamilyEvent.EXPAND);
				double fcontract = P.getFamilyEvent(node,  FamilyEvent.CONTRACT);
				double fsurviving = P.getNodeMean(node);
				double fcopies = mean(pN);
				double ftruecopies = fcopies +  delta_copies_true[node];
				
				double fpc = fp + delta_family_present[node];
				
				double varp = fp*pN[0];
				
				if (want_families)
				{
					out.printf("%s\t%d\t%d", prefix_families, f, node);
					out.printf("\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n" 
								, fp, fm, fpc
								, fgain, floss, fexpand, fcontract
								, fsurviving, ftruecopies, varp);
				}
				
				
				lineage_family_present[node] += fp;
				lineage_family_multi[node] += fm;
				lineage_family_present_corrected[node] += fpc;
				lineage_family_gain[node] += fgain;
				lineage_family_loss[node] += floss;
				lineage_family_expand[node] += fexpand;
				lineage_family_contract[node] += fcontract;
				
				lineage_copies_surviving[node] += fsurviving;
				lineage_copies_true[node] += ftruecopies;
				var_family_present[node] += varp;
			}
		}
		
		out.println(prefix_ancestral+"\tnode\tpresence:p\tmulti:m\tcorrected:p.\tgain:g\tloss:l\texpand:++\tcontract:--\tscopies:n\ttcopies:n.\tsdevp:sdevp");
		for (int node = 0; node<n; node++)
		{
			double fp = lineage_family_present[node];
			double fm = lineage_family_multi[node];
			double fpc = lineage_family_present_corrected[node];
			double fgain = lineage_family_gain[node];
			double floss = lineage_family_loss[node];
			double fexpand = lineage_family_expand[node];
			double fcontract = lineage_family_contract[node];
			
			double fsurviving = lineage_copies_surviving[node];
			double ftruecopies = lineage_copies_true[node];
			
			double varp = var_family_present[node];
			double sdevp = Math.sqrt(varp);
			
			out.printf("%s\t%s", prefix_ancestral, phylo.getIdent(node));
			out.printf("\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n" 
						, fp, fm, fpc
						, fgain, floss, fexpand, fcontract
						, fsurviving, ftruecopies, sdevp);
		}
	}
	
	private static double mean(double[] p)
	{
		double[] t = tail(p);
		double m = sum(t);
		return m;
	}
	
	/**
	 * Array of tail sums: t[i]=sum_{j&gt;i} p[j]   
	 * 
	 * @param p may be null
	 * @return array of tail sums 
	 */
	private static double[] tail(double[] p)
	{
		if (p==null) return new double[1];
		int j = p.length;
		double[] t = new double[j];
		double x = 0.0;
		while (j>0)
		{
			--j;
			t[j] = x;
			x+=p[j];
		}
		return t;
	}
	
	/**
	 * Sum across the array. 
	 * 
	 * @param x
	 * @return sum of x[i]
	 */
	private static double sum(double[] x)
	{
		double sum=0.0;
		if (x!=null)
		{
			int i=x.length;
			while (0<i)
			{
				--i;
				double y = x[i];
				sum += y;
			}
		}
		return sum;
	}
	
	public static void main(String[] args) throws Exception
	{
		Class<?> us = java.lang.invoke.MethodHandles.lookup().lookupClass();
		CommandLine cli = new CommandLine(args, us);
		PrintStream out = System.out;
		
		count.ds.Phylogeny tree = cli.getTree();
		count.ds.AnnotatedTable input_table = cli.getTable();
		MixedRateModel input_model = cli.getMixedrateModel();
		//GammaInvariant input_model = cli.getGammaModel();
		
		MixedRatePosteriors posteriors = new MixedRatePosteriors(input_model, input_table);
		
		// computation parameters
    	int absolute = cli.getOptionTruncateAbsolute();
    	double relative = cli.getOptionTruncateRelative();
        if (cli.getOptionValue(CommandLine.OPT_TRUNCATE)!=null)
        {
    		out.println(CommandLine.getStandardHeader("Truncated computation (absolute,relative)="
        				+absolute+","+relative));
        	posteriors.setCalculationWidthThresholds(absolute, relative);
		}
		double ancestor_width = 3.0;
		posteriors.setAncestorWidthThreshold(ancestor_width);
		
        int min_copies = Integer.min(2, input_table.minCopies());
        min_copies = cli.getOptionInt(CommandLine.OPT_MINCOPY, min_copies);
		posteriors.setMinimumObserved(min_copies);

		// now calculate the posterior reconstruction
		boolean want_families = cli.getOptionBoolean(OPT_HISTORY, false);
		
		posteriors.printPosteriors(out, want_families);
		
//		
//		
//		int nF = input_table.getFamilyCount();
//		for (int f=0; f<nF; f++)
//		{
//			Profile P = posteriors.getProfile(f);
//			P.reportPosteriors(System.out);
//		}
	}	
	
	
	
}
