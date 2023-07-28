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

/**
 * Algorithms for computing posteriors in a rate-variation model.
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
		
		public IndexedTree getTree()
		{
			return class_posteriors[0].factory.tree;
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
	
	public static void main(String[] args) throws Exception
	{
		CommandLine cli = new CommandLine(args, MixedRateGradient.class);
//        if (args.length!=3)
//        {
//            System.err.println("Call as $0 phylogeny table rates");
//            System.exit(2008);
//        }
//        String tree_file = args[0];
//        String table_file = args[1];
//        String rates_file = args[2];
//        count.ds.Phylogeny tree = NewickParser.readTree(
//        		GeneralizedFileReader.guessReaderForInput(tree_file));
//        count.ds.AnnotatedTable table = TableParser.readTable(tree.getLeafNames(),
//        		GeneralizedFileReader.guessReaderForInput(table_file),true);
//        GammaInvariant input_model = RateVariationParser.readRates(
//        		GeneralizedFileReader.guessReaderForInput(rates_file)
//        		, tree);
		count.ds.Phylogeny tree = cli.getTree();
		count.ds.AnnotatedTable input_table = cli.getTable();
		GammaInvariant input_model = cli.getModel();
		
		MixedRatePosteriors posteriors = new MixedRatePosteriors(input_model, input_table);
		
		int nF = input_table.getFamilyCount();
		for (int f=0; f<nF; f++)
		{
			Profile P = posteriors.getProfile(f);
			P.reportPosteriors(System.out);
		}
	}	
	
	
	
}
