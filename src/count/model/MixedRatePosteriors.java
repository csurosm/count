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
import static count.io.CommandLine.OPT_MINCOPY;


/**
 * Algorithms for computing posteriors in a rate-variation model.
 * 
 * @author csuros
 *
 */
public class MixedRatePosteriors 
{
	
	private static boolean FAMILY_PRESENCE_BY_SURVIVAL = true;

	public MixedRatePosteriors(MixedRateModel mixed_model, ProfileTable table)
	{
		this.mixed_model = mixed_model;
		this.table = table;
		this.class_posteriors = new Posteriors[mixed_model.getNumClasses()];
		this.class_index = new int[class_posteriors.length];
		Arrays.fill(class_index,-1);
		this.log_class_prob = new double[class_posteriors.length]; // precomputed
		this.num_classes = 0; // need to compute classes
		this.min_copies = Integer.min(table.minCopies(), 7);
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
		int nonzerop = 0;
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
				
//				System.out.println("#**MRP.cC class#"+num_classes
//						+"("+c+")\tpc "+pc);
				nonzerop++;
				num_classes++;
			}
		}
//		System.out.println("#**MRP.cC "+num_classes+" classes; nonzero prob "+nonzerop); 
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
//		if (min_copies<0 || min_copies>2) throw new UnsupportedOperationException("Min copies 0,1,2 are implemented");
		if (table.minCopies()<min_copies)
			throw new IllegalArgumentException("Minimum must be at least "+table.minCopies()+" in this table; requested "+min_copies);
		assert (min_copies<=table.minCopies());
		this.min_copies = min_copies;
	}
	
	public int getMinimumObserved() { return  min_copies;}
	
//	public double getUnobservedLL()
//	{
//		double unobservedLL = Double.NEGATIVE_INFINITY;
//		if (min_copies>0)
//		{
//			unobservedLL = getEmptyLL();
//			if (min_copies>1)
//			{
//				assert (min_copies==2);
//				unobservedLL = Logarithms.add(unobservedLL, getSingletonLL());
//			}
//		}
//		return unobservedLL;
//	}
	
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
		} else if (min_copies <= 2)
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
		} else { // 2<min_copies
			unobserved = new Profile[1];
			unobserved[0] = getProfile(-1);
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
			
			if (0<=family_idx){
				for (int c=0; c<num_classes; c++)
				{
					post[c] = class_posteriors[c].getPosteriors(family_idx);
					post[c].computeLikelihoods();
				}
			} else { // unbserved
				for (int c=0; c<num_classes; c++)
				{
					post[c] = class_posteriors[c].getUnobservedPosteriors(min_copies-1);
				}
			}
			cat_LL = new double[num_classes];
			cat_post = new double[num_classes];
			cat_log_post = new double[num_classes];
//			System.out.println("#**MRP.P() "+family_idx+"\talloc");
			initDataStructures();
//			System.out.println("#**MRP.P() "+family_idx+"\tdone");
		}
		private final Posteriors.Profile[] post;
		private double LL;
		private double[] cat_LL;
		/**
		 * Posterior class probabilities
		 */
		private final double[] cat_post;
		private final double[] cat_log_post;
		
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
		private double[][] log_survival_events;
		
		private void initDataStructures()
		{
			for (int c=0; c<num_classes; c++)
			{
				//Likelihood.Profile cLik = post[c].inside;
				double ll = post[c].getLogLikelihood();
				cat_LL[c]=ll;
			}
			LL = calculateLL(cat_LL);

			for (int c=0; c<num_classes; c++)
			{
				double log_p =
				cat_log_post[c] = cat_LL[c]+log_class_prob[c]-LL;
				cat_post[c]=Math.exp(log_p);
			}
			int num_nodes = getTree().getNumNodes();
			events = new double[num_nodes][];
			log_survival_events = new double[num_nodes][];
			
//			System.out.println("#**MRP.P.iDS classes ");
//			for (int node=0; node<num_nodes; node++)
//			{
//				events[node]=new double[Posteriors.FamilyEvent.values().length];
//				log_survival_events[node]=new double[Posteriors.FamilyEvent.values().length];
//				Arrays.fill(log_survival_events[node], Double.NEGATIVE_INFINITY);
//				for (int c=0; c<num_classes; c++)
//				{
//					Posteriors.Profile P = post[c];
//					double[] ev = P.getFamilyEventPosteriors(node);
//					double[] sv = P.getFamilyLogSurvivalEventPosteriors(node);
//					assert (ev.length==events[node].length);
//					assert (sv.length==ev.length);
//					for (int j=0; j<ev.length; j++)
//					{
//						events[node][j] += cat_post[c]*ev[j];
//						log_survival_events[node][j] = Logarithms.add(log_survival_events[node][j], cat_log_post[c]+sv[j]);
//					}
//				}
//			}
//			System.out.println("#**MRP.P.iDS events ");
		}
		
		
		public double getLL() { return LL;}
		
		
		
//		public double getCorrectedLL()
//		{
//			double not_observed = getUnobservedLL();
//			double pobs = -Math.expm1(not_observed); // 1-e^loglik
//			return LL-Math.log(pobs);
//		}
		
		/**
		 * Posterior class probabilities
		 * 
		 * @param class_idx
		 * @return
		 */
		public double getClassProbability(int class_idx)
		{
			return cat_post[class_idx];
		}
		
		private double sumClasses(IntFunction<Double> class_stats)
		{
			double sum = 0.0;
			for (int c=0; c<num_classes; c++)
			{
				
				double f = class_stats.apply(c);
				sum += cat_post[c]*f;
			}
			return sum;
		}
		private double sumPostLog(IntFunction<Double> cat_log_stats)
		{
			double sumLog = Double.NEGATIVE_INFINITY;
			int ncat = cat_log_post.length;
			for (int k=0; k<ncat; k++)
			{
				sumLog = Logarithms.add(sumLog, cat_log_post[k]+cat_log_stats.apply(k));
			}
			return sumLog;
		}
		
		
		
		public double getNodeMean(int node)
		{
//			return sumClasses(c->post[c].getNodeMean(node));
			double logN = sumPostLog(k->post[k].getLogNodeMean(node));
			
			return Math.exp(logN);	
		}
		
		public double getEdgeMean(int node)
		{
			double logS = sumPostLog(k->post[k].getLogEdgeMean(node));
			return Math.exp(logS);
//			return sumClasses(c->post[c].getEdgeMean(node));
		}
		
		/**
		 * Expected number of non-inherited ancestral copies at parent.
		 * 
		 * @param node
		 * @return
		 */
		public double getDeathMean(int node)
		{
			double logN_S = this.sumPostLog(k->post[k].getLogEdgeDecrease(node));
			return Math.exp(logN_S);
//
//			IndexedTree tree = mixed_model.getBaseModel().getTree();
//			if (tree.isRoot(node))
//				return 0.0;
//			else
//			{
//				int parent = tree.getParent(node);
//				double N = getNodeMean(parent);
//				double S = getEdgeMean(node);
//				return Double.max(0.0, N-S);
//			}
		}		
		
		/**
		 * Expected number of new ancestral copies at node (by duplication and gain)
		 * 
		 * @param node
		 * @return
		 */
		public double getBirthMean(int node)
		{
			double logN_S = this.sumPostLog(k->post[k].getLogNodeIncrease(node));
			return Math.exp(logN_S);
//			double S = getEdgeMean(node);
//			double N = getNodeMean(node);
//			
//			return Double.max(0.0,N-S);
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
		 * Logarithm of the posterior probabilities for ancestral copy numbers
		 * 
		 * @param node
		 * @return at least 3-element array
		 */
		public double[] getLogNodePosteriors(int node)
		{
			double[] log_pN = new double[3];
			Arrays.fill(log_pN, Double.NEGATIVE_INFINITY);
			
			
			// P{N=n} = sum_k p_k * P{N=n|k}
			
			
			for (int k=0; k<post.length; k++)
			{
				double[] log_cN = post[k].getLogNodePosteriors(node);
				int n=0;
				do
				{
					log_pN[n] = Logarithms.add(log_pN[n], log_cN[n]+cat_log_post[k]);
					n++;
				} while (n<log_pN.length && n<log_cN.length);
				if (n<log_cN.length)
				{
					log_pN = Arrays.copyOf(log_pN, log_cN.length);
					// no addition in the new cells: just copy 
					do
					{
						log_pN[n] = log_cN[n]+cat_log_post[k];
						n++;
					} while (n<log_cN.length);
				}
			}
			return log_pN;
		}
		
		
		/**
		 * Posterior distribution of ancestral copy numbers
		 * 
		 * @param node
		 * @return at least 3-element array
		 */
		public double[] getNodePosteriors(int node)
		{
			double[] log_pN = getLogNodePosteriors(node);
			double[] getNodePosteriors = new double[log_pN.length];
			for (int n=0; n<log_pN.length; n++)
				getNodePosteriors[n] = Math.exp(log_pN[n]);
			return getNodePosteriors;
		}
		
		public double[] getLogEdgePosteriors(int node)
		{
			double[] log_pS = new double[3];
			Arrays.fill(log_pS, Double.NEGATIVE_INFINITY);
			for (int k=0; k<post.length; k++)
			{
				double[] log_cS = post[k].getLogEdgePosteriors(node);
				int s=0;
				do
				{
					log_pS[s] = Logarithms.add(log_pS[s], log_cS[s]+cat_log_post[k]);
					s++;
				} while (s<log_pS.length && s<log_cS.length);
				if (s<log_cS.length)
				{
					int len = log_pS.length;
					log_pS = Arrays.copyOf(log_pS, log_cS.length);
					// no addition in the new cells: just copy 
					do
					{
						log_pS[s] = log_cS[s]+cat_log_post[k];
						s++;
					} while (s<log_cS.length);
				}
			} // for k 
			return log_pS;
		}
		
		
		private double[][] sumPostLogMatrix(IntFunction<double[][]> func)
		{
			double[][] logsum = new double[0][];
			for (int k=0; k<post.length; k++)
			{
				double[][] log_tNS =func.apply(k);
				int n=0; 
				while (n<logsum.length && n<log_tNS.length)
				{
					assert (logsum[n]!=null);
					int s = 0; 
					while (s<logsum[n].length && s<log_tNS[n].length)
					{
						logsum[n][s]=Logarithms.add(logsum[n][s],log_tNS[n][s]+cat_log_post[k]);
						++s;
					}
					if (s<log_tNS[n].length)
					{
						logsum[n] = Arrays.copyOf(logsum[n], log_tNS[n].length);
						do 
						{
							logsum[n][s]=log_tNS[n][s]+cat_log_post[k];
							++s;
						} while (s<logsum[n].length);
					}
					++n;
				}
				if (n<log_tNS.length)
				{
					logsum = Arrays.copyOf(logsum, log_tNS.length);
					do 
					{
						logsum[n]=log_tNS[n];
						for (int s=0; s<logsum[n].length; s++)
							logsum[n][s] += cat_log_post[k];
						++n;
					} while (n<logsum.length);
				}
			}
			return logsum;
		}
		
		/**
		 * Posterior probabilities for death (loss) changes
		 * on logarithmic scale : post[<var>n</var>][<var>s</var>] is the 
		 * logarithm of the posterior probability
		 * for <var>n</var>&rarr;<var>s</var> transition  
		 * to have <var>s</var> surviving copies at node
		 * with <var>n</var>  at parent.
		 * 
		 * @param node
		 * @return post[][]
		 */
		public double[][] getLogEdgeTransitionPosteriors(int node)
		{
			return sumPostLogMatrix(k->post[k].getLogEdgeTransitionPosteriors(node));
		}
		
		/**
		 * Posterior probabilities for birth (gain+duplication) changes
		 * on logarithmic scale : post[<var>n</var>][<var>s</var>] is the 
		 * logarithm of the posterior probability
		 * for <var>s</var>&rarr;<var>n</var> transition  
		 * to have <var>n</var> surviving copies at node
		 * with <var>s</var> inherited among them. 
		 * 
		 * @param node
		 * @return
		 */
		public double[][] getLogNodeTransitionPosteriors(int node)
		{
			return sumPostLogMatrix(k->post[k].getLogNodeTransitionPosteriors(node));
		}
		
		public double[] getEdgePosteriors(int node)
		{
			double[] log_pS = getLogEdgePosteriors(node);
			double[] getEdgePosteriors = new double[log_pS.length];
			for (int s=0; s<log_pS.length; s++)
			{
				getEdgePosteriors[s] = Math.exp(log_pS[s]);
			}
			return getEdgePosteriors;
		}
		
		/**
		 * Posterior probabilities for ancestor copy numbers 0,1,...
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
					pN[n] += cat_post[c]*cN[n];
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
			if (events[node]==null)
			{
				events[node]=new double[Posteriors.FamilyEvent.values().length];
				for (int c=0; c<num_classes; c++)
				{
					Posteriors.Profile P = post[c];
					double[] ev = P.getFamilyEventPosteriors(node);
					assert (ev.length==events[node].length);
					for (int j=0; j<ev.length; j++)
					{
						events[node][j] += cat_post[c]*ev[j];
					}
				}
			}
			final int event_idx = event_type.ordinal();
			return events[node][event_idx];
		}
		
		public double getFamilyLogSurvivalEvent(int node, Posteriors.FamilyEvent event_type)
		{
			if (log_survival_events[node]==null)
			{
				log_survival_events[node]=new double[Posteriors.FamilyEvent.values().length];
				Arrays.fill(log_survival_events[node], Double.NEGATIVE_INFINITY);
				for (int c=0; c<num_classes; c++)
				{
					Posteriors.Profile P = post[c];
					double[] sv = P.getFamilyLogSurvivalEventPosteriors(node);
					assert (sv.length==log_survival_events[node].length);
					for (int j=0; j<sv.length; j++)
					{
						log_survival_events[node][j] = Logarithms.add(log_survival_events[node][j], cat_log_post[c]+sv[j]);
					}
				}
			}
			final int event_idx = event_type.ordinal();
			return log_survival_events[node][event_idx];
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
	
	private void printPosteriors(PrintStream out, boolean want_families, int selected_node)
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
		double Lobs = Logarithms.logToLogComplement(L0);
//		double Lobs = Math.log(p_obs);
		
		// expected number of times unobserved profile comes up for each observed one
		double[] factor = new double[unobserved.length];
		double[] log_factor = new double[unobserved.length];
		for (int i=0; i<unobserved.length; i++)
		{
			double Ldiff = unobserved[i].getLL()-Lobs;
			log_factor[i] = Ldiff;
			factor[i] = Math.exp(Ldiff);
		}			
		
		IndexedTree phylo = getTree();
		int num_nodes = phylo.getNumNodes();
		// unobserved corrections
		double[] delta_family_present = new double[num_nodes];
		double[] delta_family_multi = new double[num_nodes];
		double[] delta_family_gain = new double[num_nodes];
		double[] delta_family_loss = new double[num_nodes];
		double[] delta_family_expand = new double[num_nodes];
		double[] delta_family_contract = new double[num_nodes];
		double[] delta_copies_surviving = new double[num_nodes];
		double[] delta_copies_true = new double[num_nodes];
		
		double[] delta_copies_birth = new double[num_nodes];
		double[] delta_copies_death = new double[num_nodes];
		
		
		int nodemin, nodemax;
		if (selected_node==-1) {
			nodemin=0; nodemax=num_nodes;
		} else {
			nodemin=selected_node; nodemax=selected_node+1;
		}
		
		for (int f=0; f<unobserved.length; f++)
		{
			Profile P = unobserved[f];
			for (int node=nodemin; node<nodemax; node++)
			{
				if (FAMILY_PRESENCE_BY_SURVIVAL)
				{
					double[] log_pn = P.getLogNodePosteriors(node);
					double log_present = Double.NEGATIVE_INFINITY;
					double log_mean = Double.NEGATIVE_INFINITY;
					
					int n=log_pn.length-1;
					while (1<n)
					{
						log_present = Logarithms.add(log_present, log_pn[n]);
						log_mean = Logarithms.add(log_mean, log_pn[n]+Math.log(n));
						--n;
					}
					double log_multi = log_present;
					if (0<n)
					{
						log_present = Logarithms.add(log_present, log_pn[n]);
						log_mean = Logarithms.add(log_mean, log_pn[n]);
						--n;
					}
					
					delta_family_present[node] +=  Math.exp(log_factor[f]+log_present);
					delta_family_multi[node] += Math.exp(log_factor[f]+log_multi);
					delta_family_gain[node] += Math.exp(log_factor[f]+P.getFamilyLogSurvivalEvent(node, FamilyEvent.GAIN));
					delta_family_loss[node] += Math.exp(log_factor[f]+P.getFamilyLogSurvivalEvent(node, FamilyEvent.LOSS));
					delta_family_expand[node] += Math.exp(log_factor[f]+P.getFamilyLogSurvivalEvent(node, FamilyEvent.EXPAND));
					delta_family_contract[node] += Math.exp(log_factor[f]+P.getFamilyLogSurvivalEvent(node, FamilyEvent.CONTRACT));
					delta_copies_surviving[node] += Math.exp(log_factor[f]+log_mean);
				} else
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
				double birth = P.getBirthMean(node);
				delta_copies_birth[node] += factor[f]*birth;
				if (!phylo.isRoot(node))
				{
					double death = P.getDeathMean(node);
					delta_copies_death[node] += factor[f]*death;
				}
			} // for node 
		
		}

		String prefix_ancestral = OPT_ANCESTRAL.toUpperCase();
		String prefix_families = OPT_HISTORY.toUpperCase();
		
		if (want_families)
		{
			out.println(prefix_families+"\tfamily\tnodeidx\tpresence:p\tmulti:m\tcorrected:p.\tgain:g\tloss:l\texpand:++\tcontract:--\tscopies:n\ttcopies:n.\tvarp:vp\tbirth:+\tdeath:-");
		}
		
		// lineage totals
		double[] lineage_family_present = new double[num_nodes];
		double[] lineage_family_present_corrected = new double[num_nodes];
		double[] lineage_family_multi = new double[num_nodes];
		double[] lineage_family_gain = new double[num_nodes];
		double[] lineage_family_loss = new double[num_nodes];
		double[] lineage_family_expand = new double[num_nodes];
		double[] lineage_family_contract = new double[num_nodes];
		double[] lineage_copies_surviving = new double[num_nodes];
		double[] lineage_copies_true = new double[num_nodes];
		double[] var_family_present = new double[num_nodes];
		double[] lineage_copies_birth = new double[num_nodes];
		double[] lineage_copies_death = new double[num_nodes];
		
		int[] max_family_sizes = table.getMaxFamilySizes(phylo);
		
		double[][] lineage_copies_distribution = new double[num_nodes][];
		for (int node=nodemin; node<nodemax; node++) {
			int max = max_family_sizes[node];
			lineage_copies_distribution[node]=new double[max+1];
			Arrays.fill(lineage_copies_distribution[node], Double.NEGATIVE_INFINITY);
		}
			
		
		for (int f=0; f<table.getFamilyCount(); f++)
		{
			Profile P = getProfile(f);
			double fp,fm,fgain,floss,fexpand,fcontract,fsurviving,ftruecopies,fpc,varp;
			
			for (int node=nodemin; node<nodemax; node++)
			{
				if (FAMILY_PRESENCE_BY_SURVIVAL)
				{
					double[] log_pn = P.getLogNodePosteriors(node);
					double log_present = Double.NEGATIVE_INFINITY;
					double log_mean = Double.NEGATIVE_INFINITY;
					
					int n=log_pn.length-1;
					if (lineage_copies_distribution[node].length<=n) {
						// got n=2 len=1 ? no worries 
						int len = lineage_copies_distribution[node].length;
						lineage_copies_distribution[node] = Arrays.copyOf(lineage_copies_distribution[node], n+1);
						for (int c=len; c<=n; c++)
							lineage_copies_distribution[node][c] = Double.NEGATIVE_INFINITY;
					}
					
					while (1<n)
					{
						lineage_copies_distribution[node][n] = Logarithms.add(lineage_copies_distribution[node][n], log_pn[n]);
						log_present = Logarithms.add(log_present, log_pn[n]);
						log_mean = Logarithms.add(log_mean, log_pn[n]+Math.log(n));
						--n;
					}
					double log_multi = log_present;
					if (0<n) // == 1
					{
						lineage_copies_distribution[node][n] = Logarithms.add(lineage_copies_distribution[node][n], log_pn[n]);
						log_present = Logarithms.add(log_present, log_pn[n]);
						log_mean = Logarithms.add(log_mean, log_pn[n]);
						--n;
					}
					
					fp = Math.exp(log_present);
					fm = Math.exp(log_multi);
					fgain = Math.exp(P.getFamilyLogSurvivalEvent(node, FamilyEvent.GAIN));
					floss = Math.exp(P.getFamilyLogSurvivalEvent(node, FamilyEvent.LOSS));
					fexpand = Math.exp(P.getFamilyLogSurvivalEvent(node, FamilyEvent.EXPAND));
					fcontract = Math.exp(P.getFamilyLogSurvivalEvent(node,  FamilyEvent.CONTRACT));
					fsurviving =  Math.exp(log_mean);
					ftruecopies = fsurviving+delta_copies_surviving[node];
					
					fpc = fp + delta_family_present[node];
					
					varp = Math.exp(log_present+log_pn[0]);
					
				} else
				{
					double[] pN = P.getNodeAncestor(node);
					
					double sum1;
					double sum2=0.0;
					for (int i=2; i<pN.length; i++)
						sum2+=pN[i];
					sum1 = sum2+pN[1];
					
					fp = Double.max(1.0-pN[0],sum1);
					fm = Double.max(1.0-(pN[0]+pN[1]),sum2);
					fgain = P.getFamilyEvent(node, FamilyEvent.GAIN);
					floss = P.getFamilyEvent(node, FamilyEvent.LOSS);
					fexpand = P.getFamilyEvent(node, FamilyEvent.EXPAND);
					fcontract = P.getFamilyEvent(node,  FamilyEvent.CONTRACT);
					fsurviving = P.getNodeMean(node);
					double fcopies = mean(pN);
					ftruecopies = fcopies +  delta_copies_true[node];
					
					fpc = fp + delta_family_present[node];
					
					varp = fp*pN[0];
				} // whether counting families by surviving members 
				
				double fbirth = P.getBirthMean(node)+delta_copies_birth[node];
				double fdeath;
				if (phylo.isRoot(node))
				{
					fdeath = 0.0;
				} else 
				{
					fdeath = P.getDeathMean(node)+delta_copies_death[node];
				}
				if (want_families)
				{
					out.printf("%s\t%d\t%d", prefix_families, f, node);
//					out.printf("%s\t%s\t%d", prefix_families, table.getFamilyName(f), node);
					out.printf("\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n" 
								, fp, fm, fpc
								, fgain, floss, fexpand, fcontract
								, fsurviving, ftruecopies, varp, fbirth, fdeath);
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
				lineage_copies_birth[node] += fbirth;
				lineage_copies_death[node] += fdeath;
			} // for node 
		}
		
		out.println(prefix_ancestral+"\tnode\tpresence:p\tmulti:m\tcorrected:p.\tgain:g\tloss:l\texpand:++\tcontract:--\tscopies:n\ttcopies:n.\tsdevp:sdevp\tbirth:+.\tdeath:-.\tsharedfam:p^\tsharedcopy:n^");
		for (int node = nodemin; node<nodemax; node++)
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
			double sdevp = Math.sqrt(varp/table.getFamilyCount());
			
			double fbirth = lineage_copies_birth[node];
			double fdeath = lineage_copies_death[node];
			
			
			
			double fshared = fp;
			double nshared = fsurviving;
			
			for (int c=0; c<phylo.getNumChildren(node); c++) {
				int child = phylo.getChild(node, c);
				nshared -= lineage_copies_death[child];
				fshared -= lineage_family_loss[child];
			}
			
			out.printf("%s\t%s", prefix_ancestral, phylo.getIdent(node));
			out.printf("\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n" 
						, fp, fm, fpc
						, fgain, floss, fexpand, fcontract
						, fsurviving, ftruecopies, sdevp, fbirth, fdeath
						, fshared, nshared);
		}
	}
	
	
	private void printTransitions(PrintStream out, int node)
	{
		
		computeClasses();
		
		
		double[][] log_death = null;
		double[][] log_birth = null;
		
		int nF = table.getFamilyCount();
		for (int f = 0; f<nF; f++)
		{
			Profile P = getProfile(f);
			double[][] fdeath = P.getLogEdgeTransitionPosteriors(node);
			log_death = Logarithms.add(log_death, fdeath);
			double[][] fbirth = P.getLogNodeTransitionPosteriors(node);
			log_birth = Logarithms.add(log_birth, fbirth);
		}
		
		String prefix_node = "#NODE\t";
		out.print(prefix_node+"node\tident");
		if (1<num_classes)
			out.print("\tclass\tpclass");
		out.println("\tp\tq\tr\tp~\tq~\tr~\t1-ext(parent)");
		String node_ident=null; 
		for (int ci=0; ci<num_classes; ci++)
		{
			assert (log_class_prob[ci]!=Double.NEGATIVE_INFINITY);
			{
				Posteriors P = class_posteriors[ci];
				node_ident = P.factory.tree.getIdent(node);
				out.print(prefix_node+node+"\t"+node_ident);
				if (1<num_classes)
					out.printf("\t%d\t%g", ci, Math.exp(log_class_prob[ci]));
				double logp = P.factory.rates.getLogLossParameter(node);
				double logq = P.factory.rates.getLogDuplicationParameter(node);
				double r = P.factory.rates.getGainParameter(node);
				if (logp!=Double.NEGATIVE_INFINITY) r*=Math.exp(logq);
				double logps = P.factory.getLogLossParameter(node);
				double logqs = P.factory.getLogDuplicationParameter(node);
				double rs = P.factory.getGainParameter(node);
				if (logqs != Double.NEGATIVE_INFINITY)
					rs *= Math.exp(logqs);
				double log_e1 ;
				if (P.factory.tree.isRoot(node))
				{
					log_e1 = Double.NEGATIVE_INFINITY;
				} else
				{
					int parent =  P.factory.tree.getParent(node);
					log_e1 = P.factory.getLogExtinctionComplement(parent);
				}
				
				out.printf("\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n"
						, Math.exp(logp), Math.exp(logq), r
						, Math.exp(logps), Math.exp(logqs), rs
						, Math.exp(log_e1)
						);
			}
		}
		

		
		int nmax = log_birth.length; // filter out < ulp probabilities 
		int smax = 0; // inclusive
		while (0<nmax)
		{
			double max=0.0;
			for (int s=0; s<log_birth[nmax-1].length; s++)
			{
				double lb = log_birth[nmax-1][s];
				if (lb != Double.NEGATIVE_INFINITY) smax= Integer.max(smax, s);
				max = Double.max(max, Math.exp(lb));
			}
			if (0.0<max) break;
			--nmax;
		}
		String prefix_birth = "#BIRTH\t";
		out.print(prefix_birth+"node\tident\ts /(n-s)=");
		for (int t=0; t<nmax; t++)
			out.print("\t"+t);
		out.println();
		
		for (int s=0; s<=smax; s++)
		{
			out.print(prefix_birth+node+"\t"+node_ident+"\t"+s);
			int t=0, n=s+t;
			while (n<nmax)
			{
				double pbirth = Math.exp(log_birth[n][s]);
				out.printf("\t%g", pbirth);
				t++;
				n++;
			}
			out.println();
		}
//		
		nmax = log_death.length; // exclusive
		smax = 0; // inclusive
		while (0<nmax)
		{
			double max=0.0;
			for (int s=0; s<log_death[nmax-1].length; s++)
			{
				double lb = log_death[nmax-1][s];
				if (lb != Double.NEGATIVE_INFINITY) smax= Integer.max(smax, s);
				max = Double.max(max, Math.exp(lb));
			}
			if (0.0<max) break;
			--nmax;
		}
		
		String prefix_death = "#DEATH\t";
		out.print(prefix_death+"node\tident\tn /(n-s)=");
		for (int t=0; t<=smax; t++)
			out.print("\t"+t);
		out.println();

		for (int n=0; n<nmax; n++)
		{
			out.print(prefix_death+node+"\t"+node_ident+"\t"+n);
			int t = 0, s=n;
			while (log_death[n].length<=s)
			{
				out.printf("\t%g", 0.0);
				++t; --s;
			}
			do 
			{
				out.printf("\t%g", Math.exp(log_death[n][s]));
				++t; --s;
			} while (0<=s);
			out.println();
		}
	}
	
	
	
	private void printNodePosteriors(PrintStream out, int node) {
		
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
		out.println(CommandLine.getStandardHeader("Minimum observed copies: -"+OPT_MINCOPY+" "+min_copies));
		
		
		
		int stats = cli.getOptionInt(CommandLine.OPT_STATISTICS, -1);
		if (0<=stats)
		{
			posteriors.printTransitions(out, stats);
		} else
		{
			boolean want_families = cli.getOptionBoolean(OPT_HISTORY, false);
			out.println(CommandLine.getStandardHeader("Detailed family posteriors: -"+OPT_HISTORY+" "+want_families));
			
			boolean family_survival = cli.getOptionBoolean(CommandLine.OPT_SURVIVAL, FAMILY_PRESENCE_BY_SURVIVAL);
			out.println(CommandLine.getStandardHeader("Family presence by survival: -"+CommandLine.OPT_SURVIVAL+" "+family_survival));
			FAMILY_PRESENCE_BY_SURVIVAL = family_survival;
			
			int selected_node = cli.getOptionInt(CommandLine.OPT_NODE, -1);
			if (selected_node != -1) {
				out.println(CommandLine.getStandardHeader("Selected node: -"+CommandLine.OPT_NODE+" "+selected_node+"\t ("+tree.toString(selected_node)+")"));
			}

			
			// now calculate the posterior reconstruction
			posteriors.printPosteriors(out, want_families, selected_node);
			
		}
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
