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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Set;
import java.util.function.ToIntFunction;

import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.Phylogeny;
import count.ds.TreeTraversal;
import count.io.CommandLine;
import count.io.TableParser;
import count.model.Posteriors.FamilyEvent;

import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_RND;

import static count.io.CommandLine.OPT_ANCESTRAL;
import static count.io.CommandLine.OPT_HISTORY;
import static count.io.CommandLine.OPT_STATISTICS;

/**
 * Precise, simulation of copy number evolution based 
 * on copy-level discrete events. 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */

public class SimulatedEvolution 
{
//	private static long RND_SEED=0L;
	
	private static boolean FAMILY_PRESENCE_BY_SURVIVAL = true;
	
	
	protected SimulatedEvolution(MixedRateModel rates_model, Random RND)
	{
		this (rates_model, 
				(RND==null?(new Random()).nextLong():RND.nextLong())
			);
	}
	
	protected SimulatedEvolution(MixedRateModel rates_model, long rnd_seed)
	{
		this.rnd_seed = rnd_seed;
		this.RND = new Random(rnd_seed);
		this.mixed_model = rates_model;
		phylo=rates_model.getClassModel(0).getTree();
		post_order = TreeTraversal.postOrder(phylo); // get it only once
		pre_order = TreeTraversal.preOrder(phylo);
		this.num_classes = rates_model.getNumClasses();
		this.class_rates = new TreeWithRates[num_classes];
		this.class_indexes = new int[num_classes];
		this.class_probs = new double[num_classes]; // precomputed
		this.class_root_cdfs = new double[num_classes][];
		initDataStructures();
	}
	
	public SimulatedEvolution(MixedRateModel rates_model)
	{
		this(rates_model, null);
	}
	
	private final MixedRateModel mixed_model;
	
	protected MixedRateModel getModel() { return mixed_model;}
	protected IndexedTree getTree() { return phylo;}
	
	private final long rnd_seed;
	
	protected final Random RND;
	
	
	private final IndexedTree phylo;
	
	/**
	 * Non-zero classes
	 */
	private int num_classes; 
	/**
	 * Original class index of active rate classes 
	 */
	private final int[] class_indexes;
	/**
	 * Cumulative distribution function for active rate classes 
	 */
	private final double[] class_probs;
	private final TreeWithRates[] class_rates;
	private final double[][] class_root_cdfs;
	private final int[] post_order;
	private final int[] pre_order;
	
	/**
	 * The random seed used to initialize the 
	 * {@link java.util.Random} generator. 
	 */
//	public long getRandomSeed() {return this.rnd_seed;}
	
	private void initDataStructures()
	{
		this.num_classes = 0;
		double cdf=0.0;
		for (int c = 0; c<class_rates.length; c++)
		{
			double pc = mixed_model.getClassProbability(c);
			if (pc != 0.0)
			{
				TreeWithRates rates =mixed_model.getClassModel(c);
				class_rates[num_classes] = rates;
				class_indexes[num_classes] = c;
				cdf = class_probs[num_classes] = cdf + pc;
				class_root_cdfs[num_classes]=getRootCDF(rates, 2); // initialized for 0,1 but will expand as necessary as random root copies are drawn 
				num_classes++;
			}
		}
		if (num_classes == 0)
			throw new IllegalArgumentException("All classes have 0.0 probability");
		class_probs[num_classes-1]=1.0; 
	}
	
	private synchronized void expandRootCDF(int class_idx)
	{
		TreeWithRates rates = class_rates[class_idx];
		double[] cdf = class_root_cdfs[class_idx];
		
		class_root_cdfs[class_idx] = getRootCDF(rates, 2*(cdf.length-1)); // doubling the maximum value		
	}
	
    private double[] getRootCDF(TreeWithRates rates, int max_value)
    {
        double[] pmf = rates.getRootDistribution().getPointMassFunction(max_value);
        double[] cdf = new double[max_value+1];
        cdf[0] = pmf[0];
        for (int i=1; i<=max_value; i++)
            cdf[i]=cdf[i-1]+pmf[i];
        return cdf;
    }
	
	
	private double nextRNDExponential(double lambda)
	{
		return -Math.log(RND.nextDouble())/lambda;        
	}
	
	/**
	 * Next random category 
	 * 
	 * @return
	 */
	private int nextRNDClass()
	{
		double p = RND.nextDouble();
        int x=Arrays.binarySearch(class_probs, 0, num_classes, p);
        if (x<0) // binarysearch does not find this value: very likely
            x = -(x+1);
//        System.out.println("#**SE.nRC "+x+"\tp "+p+"\t"+Arrays.toString(class_probs));
        return x;
	}
	
	
	/**
	 * Member of a family at a tree node, with associated random birth/death events.  
	 * 
	 */
	private class Copy implements Comparable<Copy>
	{
		Copy progenitor;
		private double timer_start;
		private double time_to_die;
		private double time_to_spawn;
		int node;
		
		Copy(int node)
		{
			this(node, 0.0);
		}
		Copy(int node, double time)
		{
			this.timer_start = time;
			this.time_to_die = Double.POSITIVE_INFINITY;
			this.time_to_spawn = Double.POSITIVE_INFINITY;
			this.progenitor = this;
			this.node = node;
		}
		
		/**
		 * Creates a copy with same timer start.
		 * 
		 * @param parent
		 */
		Copy(Copy parent)
		{
			this.timer_start = parent.timer_start;
			this.time_to_die = Double.POSITIVE_INFINITY;
			this.time_to_spawn = Double.POSITIVE_INFINITY;
			this.progenitor = parent.progenitor;
			this.node = parent.node;
		}
		
		Copy(Copy progenitor, int node)
		{
			this.timer_start = 0.0;
			this.progenitor = progenitor;
			this.node = node;
			this.time_to_die = Double.POSITIVE_INFINITY;
			this.time_to_spawn = Double.POSITIVE_INFINITY;
		}
		
		/**
		 * Spawn time after timer start.
		 * 
		 * @param x
		 */
        void setSpawnTime(double x)
        {
            this.time_to_spawn = x;
        }
        
        /**
         * Death time after timer start.
         * 
         * @param x
         */
        void setDeathTime(double x)
        {
            this.time_to_die = x;
        }
        
        /**
         * Next event time (die or spawn). 
         * @return
         */
        double nextEventTime()
        {
        	return timer_start + Double.min(time_to_spawn, time_to_die);
        }
		
		/**
		 * Whether this copy has an event before the other copy.
		 */
		@Override
		public int compareTo(Copy other)
		{
			return Double.compare(this.nextEventTime(), other.nextEventTime());
		}
		
		/**
		 * Whether this copy dies before it spawns a duplicate.
		 * @return
		 */
        boolean dies()
        {
            return time_to_die <= time_to_spawn;
        }

        boolean isProgenitor()
        {
        	return this == progenitor;
        }
        
        Copy spawn()
        {
        	timer_start += time_to_spawn;
        	assert (!dies());
        	time_to_die -= time_to_spawn;
        	time_to_spawn = 0.0;
            Copy child = new Copy(this);
            return child;
        }
        
        Copy inherit(int node)
        {
        	Copy child = new Copy(this, node);
        	return child;
        }
        
		 
	}
	
	public Table table(int num_rows, int min_observed)
	{
		return new Table(num_rows, min_observed);
	}

	public static Table table(MixedRateModel rates_model, long rnd_seed, int num_rows, int min_observed)
	{
		SimulatedEvolution sim = new SimulatedEvolution(rates_model, rnd_seed);
		return sim.table(num_rows, min_observed);
	}
	
	public static Table table(MixedRateModel rates_model, int num_rows, int min_observed)
	{
		SimulatedEvolution sim = new SimulatedEvolution(rates_model);
		return sim.table(num_rows, min_observed);
	}
	
	
	/**
	 * A table with initially unfilled profiles. 
	 * Calling {@link #getObservedProfile(int)} 
	 * generates the random profile for a family. (Subsequent calls return the same.) 
	 * 
	 * @author csuros
	 *
	 */
	public class Table extends AnnotatedTable
	{
		private Table(int num_rows, int min_observed)
		{
			super(phylo.getLeafNames());
			this.row_count=num_rows;
			this.min_observed = min_observed;
			this.observed_profiles = new ArrayList<>(num_rows);
			initDataStructures();
		}
		private final int row_count;
		
		private final List<ObservedProfile> observed_profiles;
		private final int min_observed;
		
		public int getMinimumObserved() { return min_observed;}
		public int getNumClasses()
		{
			return num_classes;
		}
		
		@Override
		public boolean isBinaryTable()
		{
			return false;
		}
		
		/**
		 * The random seed used to initialize the 
		 * {@link java.util.Random} generator at instantiation. 
		 */
		public long getRandomSeed() { return rnd_seed;}
		
		private void initDataStructures()
		{
			int[][] copy_numbers = new int[row_count][this.getTaxonCount()];
			for (int[] rows: copy_numbers)
			{
				Arrays.fill(rows, -1); // missing
				observed_profiles.add(null);
			}
			setTable(copy_numbers, null);
			assert (observed_profiles.size()==row_count);
		}
		
		/**
		 * Fills in and returns the random profile for a table row.  
		 * 
		 * @param family
		 * @return
		 */
		public ObservedProfile getObservedProfile(int family)
		{
			ObservedProfile obs = observed_profiles.get(family);
			if (obs == null)
			{
				obs = new ObservedProfile(family);
				observed_profiles.set(family, obs);
			}
			return obs;
		}
		
		public void fillTable()
		{
			int tot_profiles = 0;
			for (int f=0; f<row_count; f++)
			{
				ObservedProfile obs = getObservedProfile(f);
				int s = obs.profiles.size();
				tot_profiles += s;
			}
//			System.out.println("#*SE.T.fT "+row_count+"\ttot "+tot_profiles+"\thash "+this.tableHashCode());
		}
		
		
		public class ObservedProfile
		{
			private ObservedProfile(int family_idx)
			{
				profiles = new ArrayList<>();
				this.family_idx = family_idx;
				this.event_queue = new PriorityQueue<>();
				simulateTree();
				calculateStatistics();
				profiles.clear();
				event_queue.clear(); // should be empty though already
				
				// TODO
				// compute statistics from the profiles 
				// do not keep RandomProfile and Copy instances : too much memory
			}
			private final int family_idx;
			private final List<RandomProfile> profiles;

			private final PriorityQueue<Copy> event_queue;
			
			/**
			 * Precalculated statistics 
			 */
			private NodeStatistics[] unobs_node_statistics, obs_node_statistics;
			private int unobs_history_event_count, obs_history_event_count;
			private int[] unobs_family_class;
			int obs_family_class;
			private int num_unobserved_profiles;
			
			
			private class NodeStatistics
			{
				NodeStatistics(List<RandomProfile> profiles, int node, boolean want_unobserved)
				{
					node_survival_count = sumProfiles(profiles, P->P.getMemberCount(node), want_unobserved);
					node_true_count = sumProfiles(profiles, P->P.getTrueMemberCount(node), want_unobserved);
					edge_survival_count = sumProfiles(profiles, P->P.getEdgeSurvival(node), want_unobserved);
					boolean observed = !want_unobserved;
					int[] family_events = new int[FamilyEvent.values().length];
					int start = observed?profiles.size()-1:0;
					for (int i=start; i<profiles.size(); i++)
					{
						RandomProfile P = profiles.get(i);
						int[] events = P.getFamilyEvents(node);
						for (int j=0; j<events.length; j++)
							family_events[j]+=events[j];
					}
					family_events_gain = family_events[FamilyEvent.GAIN.ordinal()];
					family_events_loss = family_events[FamilyEvent.LOSS.ordinal()];
					family_events_expand = family_events[FamilyEvent.EXPAND.ordinal()];
					family_events_contract = family_events[FamilyEvent.CONTRACT.ordinal()];
					
					if (FAMILY_PRESENCE_BY_SURVIVAL) {
						family_present = sumProfiles(profiles, P->(P.getMemberCount(node)>0?1:0), want_unobserved);	
						family_multi = sumProfiles(profiles, P->(P.getMemberCount(node)>1?1:0), want_unobserved);
					} else {
						family_present = sumProfiles(profiles, P->(P.getTrueMemberCount(node)>0?1:0), want_unobserved);	
						family_multi = sumProfiles(profiles, P->(P.getTrueMemberCount(node)>1?1:0), want_unobserved);
					}
				}
				
				
				int node_survival_count;
				int node_true_count;
				int edge_survival_count;
				int family_events_gain;
				int family_events_loss;
				int family_events_expand;
				int family_events_contract;
				int family_present;
				int family_multi;
			}
			
			private void calculateStatistics()
			{
				num_unobserved_profiles = profiles.size()-1;
				int num_nodes = getTree().getNumNodes();
				unobs_node_statistics = new NodeStatistics[num_nodes];
				obs_node_statistics = new NodeStatistics[num_nodes];
				
				for (int node=0; node<num_nodes; node++)
				{
					unobs_node_statistics[node] = new NodeStatistics(profiles, node, true);
					obs_node_statistics[node] = new NodeStatistics(profiles, node, false);
				}
				
				unobs_family_class=new int[num_classes];
				obs_family_class=-1;
				
				for (int i=0;i<profiles.size(); i++)
				{
					RandomProfile P = profiles.get(i);
					int c = P.class_idx;
					unobs_family_class[c]++;
					obs_family_class=c;
				}

				unobs_history_event_count = sum(P->P.getHistoryEventCount(), true);
				obs_history_event_count = sum(P->P.getHistoryEventCount(), false);
			}
			
			
			
			private int simulateTree()
			{
				profiles.clear();
				RandomProfile last_profile;
				do
				{
					last_profile = new RandomProfile(nextRNDClass(), event_queue);
					profiles.add(last_profile);
					
				} while (last_profile.size()<Table.this.min_observed);
				
				int[] copy_numbers = getFamilyProfile(family_idx); // int[] directly in the underlying table
				for (int leaf=0; leaf<phylo.getNumLeaves(); leaf++)
				{
					copy_numbers[leaf] = last_profile.getMemberCount(leaf);
				}
//				System.out.println("#**SE.T.OP.sT "+family_idx+"\tnprof "+profiles.size()+"\tlast "+last_profile
//						+"\t// "+Table.this.getLineageCount(family_idx)
//						+"\t"+Arrays.toString(Table.this.getFamilyProfile(family_idx)));
				
				return profiles.size();
			}
			
			
			
			public int getNodeSurvivalCount(int node, boolean want_unobserved)
			{
				if (profiles.isEmpty())
				{
					NodeStatistics S = want_unobserved?unobs_node_statistics[node]:obs_node_statistics[node];
					return S.node_survival_count;
				} else
					return sum(P->P.getMemberCount(node), want_unobserved);
//				boolean observed = !want_unobserved;
//				int member_count = 0;
//				int start = observed?profiles.size()-1:0;
//				for (int i=start; i<profiles.size(); i++)
//				{
//					RandomProfile P = profiles.get(i);
//					member_count += P.getMemberCount(node);
//				}
//				return member_count;
			}
			
			public int getNodeTrueCount(int node, boolean want_unobserved)
			{
				if (profiles.isEmpty())
				{
					NodeStatistics S = want_unobserved?unobs_node_statistics[node]:obs_node_statistics[node];
					return S.node_true_count;
				} else
					return sum(P->P.getTrueMemberCount(node), want_unobserved);
			}
			
			public int getEdgeSurvivalCount(int node,  boolean want_unobserved)
			{
				if (profiles.isEmpty())
				{
					NodeStatistics S = want_unobserved?unobs_node_statistics[node]:obs_node_statistics[node];
					return S.edge_survival_count;
				} else
					return sum(P->P.getEdgeSurvival(node), want_unobserved);
//				boolean observed = !want_unobserved;
//				int edge_count = 0;
//				int start = observed?profiles.size()-1:0;
//				for (int i=start; i<profiles.size(); i++)
//				{
//					RandomProfile P = profiles.get(i);
//					edge_count+= P.getEdgeSurvival(node);
//				}
//				return edge_count;
			}
			
			public int[] getFamilyEvents(int node, boolean  want_unobserved)
			{
				int[] family_events = new int[FamilyEvent.values().length];
				if (profiles.isEmpty())
				{
					NodeStatistics S = want_unobserved?unobs_node_statistics[node]:obs_node_statistics[node];
					family_events[FamilyEvent.GAIN.ordinal()] = S.family_events_gain;
					family_events[FamilyEvent.LOSS.ordinal()] = S.family_events_loss;
					family_events[FamilyEvent.EXPAND.ordinal()] = S.family_events_expand;
					family_events[FamilyEvent.CONTRACT.ordinal()] = S.family_events_contract;
				} else
				{
					boolean observed = !want_unobserved;
					int start = observed?profiles.size()-1:0;
					for (int i=start; i<profiles.size(); i++)
					{
						RandomProfile P = profiles.get(i);
						int[] events = P.getFamilyEvents(node);
						for (int j=0; j<events.length; j++)
							family_events[j]+=events[j];
					}
				}
				return family_events;
			}
			
			public int getHistoryEventCount(boolean  want_unobserved)
			{
				if (profiles.isEmpty())
				{
					return want_unobserved?unobs_history_event_count:obs_history_event_count;
				} else
					return sum(P->P.getHistoryEventCount(), want_unobserved);
//				boolean observed = !want_unobserved;
//				int event_count = 0;
//				int start = observed?profiles.size()-1:0;
//				for (int i=start; i<profiles.size(); i++)
//				{
//					RandomProfile P = profiles.get(i);
//					event_count += P.getHistoryEventCount();
//				}
//				return event_count;
			}
			
			public int getFamilyPresent(int node, boolean want_unobserved)
			{
				if (profiles.isEmpty())
				{
					NodeStatistics S = want_unobserved?unobs_node_statistics[node]:obs_node_statistics[node];
					return S.family_present;
				} else
					return sum(P->(P.getTrueMemberCount(node)>0?1:0), want_unobserved);
			}
			
			public int getFamilyMulti(int node, boolean want_unobserved)
			{
				if (profiles.isEmpty())
				{
					NodeStatistics S = want_unobserved?unobs_node_statistics[node]:obs_node_statistics[node];
					return S.family_multi;
				} else
					return sum(P->(P.getTrueMemberCount(node)>1?1:0), want_unobserved);
			}
			
			public int getFamilySurvivalPresent(int node, boolean want_unobserved)
			{
				if (profiles.isEmpty())
				{
					NodeStatistics S = want_unobserved?unobs_node_statistics[node]:obs_node_statistics[node];
					return S.family_present;
				} else
					return sum(P->(P.getMemberCount(node)>0?1:0), want_unobserved);
			}
			
			public int getFamilySurvivalMulti(int node, boolean want_unobserved)
			{
				if (profiles.isEmpty())
				{
					NodeStatistics S = want_unobserved?unobs_node_statistics[node]:obs_node_statistics[node];
					return S.family_multi;
				} else
					return sum(P->(P.getMemberCount(node)>1?1:0), want_unobserved);
			}

			
			public int getFamilyClass(int class_idx, boolean want_unobserved)
			{
				if (profiles.isEmpty())
				{
					if (want_unobserved)
						return unobs_family_class[class_idx];
					else
						return class_idx==obs_family_class?1:0;
				} else
				return sum(P->(P.class_idx==class_idx?1:0), want_unobserved);
			}
			
			
			private int sum(ToIntFunction<RandomProfile> stats, boolean want_unobserved)
			{
				boolean observed = !want_unobserved;
				int sum = 0;
				int start = observed?profiles.size()-1:0;
				for (int i=start; i<profiles.size(); i++)
				{
					RandomProfile P = profiles.get(i);
					sum += stats.applyAsInt(P);
				}
				return sum;
			}

			private int sumProfiles(List<RandomProfile> profiles, ToIntFunction<RandomProfile> stats, boolean want_unobserved)
			{
				boolean observed = !want_unobserved;
				int sum = 0;
				int start = observed?profiles.size()-1:0;
				for (int i=start; i<profiles.size(); i++)
				{
					RandomProfile P = profiles.get(i);
					sum += stats.applyAsInt(P);
				}
				return sum;
			}
			
			
		}
		
	}

	
	private class RandomProfile 
	{
		RandomProfile(int class_idx)
		{
			this(class_idx, new PriorityQueue<>());
		}
		
		RandomProfile(int class_idx, PriorityQueue event_queue)
		{
			this.class_idx = class_idx;
			this.event_queue = event_queue;
//			this.member_copies = new int[post_order.length];
			this.node_copies = new ArrayList<>(pre_order.length);
			this.surviving_copies = new ArrayList<>(post_order.length);
			this.simulateTree();
		}
		void simulateTree()
		{
			node_copies.clear();
			for (int node=0; node<pre_order.length; node++)
				node_copies.add(new HashSet<>());
			history_event_count = simulateLineage(phylo.getRoot()); // fills up the Copy sets at the nodes
			
			surviving_copies.clear();
			for (int node=0; node<post_order.length; node++)
				surviving_copies.add(new HashSet<>());
			calculateSurvivors(phylo.getRoot());
		}
		
		private final int class_idx;
		
		private final PriorityQueue<Copy> event_queue;
		
//		private int[] member_copies;
//		private int[][] family_events;
		
		private final List<Set<Copy>> node_copies;
		private final List<Set<Copy>> surviving_copies;
		
		private int history_event_count;
		
		
		/**
		 * Sum of copy numbers across the leaves.
		 * @return
		 */
		int size()
		{
			int size = 0;
			for (int leaf=0;leaf<phylo.getNumLeaves(); leaf++)
				size += getMemberCount(leaf);
			return size;
		}
		
//		
		/**
		 * Number of copies at a node which survive at subtree leaves 
		 * 
		 * @param node
		 * @return
		 */
		int getMemberCount(int node)
		{
			return surviving_copies.get(node).size();
		}
		

		int getEdgeSurvival(int node)
		{
			int s;
			if (phylo.isRoot(node))
			{
				TreeWithRates rates = class_rates[class_idx];
				s = (rates.getEdgeLength(node)==Double.POSITIVE_INFINITY)?0:1;
			} else
			{
				int parent = phylo.getParent(node);
				Set<Copy> parent_copies = surviving_copies.get(parent);
				Set<Copy> our_copies = surviving_copies.get(node);
				Set<Copy> survivors = new HashSet<>();
				for (Copy copy: our_copies)
				{
					if (parent_copies.contains(copy.progenitor))
						survivors.add(copy.progenitor);
				}
				s = survivors.size();
			}
			return s;
		}
		
		/**
		 * All copies at a node (event those without surviving descendants)
		 * 
		 * @param node
		 * @return
		 */
		int getTrueMemberCount(int node)
		{
			return node_copies.get(node).size();
		}
		
		
		
		int[] getFamilyEvents(int node)
		{
			int[] events = new int[FamilyEvent.values().length];
			Set<Copy> members = node_copies.get(node);
			int n = members.size();
			if (phylo.isRoot(node))
			{
				TreeWithRates rates = class_rates[class_idx];
		        if (Double.isFinite(rates.getEdgeLength(node)))
		        {
			        events[FamilyEvent.GAIN.ordinal()] = 0;
		        	events[FamilyEvent.LOSS.ordinal()] = (n==0?1:0);
		        	events[FamilyEvent.EXPAND.ordinal()] = (n>1)?1:0;	 
		        } else
		        { // usual case
			        events[FamilyEvent.GAIN.ordinal()] = n>0?1:0;
		        }		        
			} else
			{
				int parent = phylo.getParent(node);
				Set<Copy> parent_copies = node_copies.get(parent);
				int pn = parent_copies.size();
				
				events[FamilyEvent.EXPAND.ordinal()] =  (n>1 && pn==1)?1:0;
				events[FamilyEvent.CONTRACT.ordinal()] = (pn>1 && n==1)?1:0;
				boolean parent_survived = false;
				for (Copy copy: members)
				{
					if (parent_copies.contains(copy.progenitor))
					{
						parent_survived = true;
						break;
					}
				}
				events[FamilyEvent.GAIN.ordinal()] = (n>0 && !parent_survived)?1:0;
				events[FamilyEvent.LOSS.ordinal()] = (pn>0 && !parent_survived)?1:0;
			}
			return events;
		}
		
		int getHistoryEventCount() { return history_event_count;}
		
		@Override
		public String toString()
		{
			return "RP#"+class_idx+"[size="+size()+",histry="+history_event_count+"]";
		}
		
		
		/**
		 * Generates random copy history in a subtree. Fills up the Copy sets at each node 
		 * in preorder ({@link #node_copies}).  
		 * 
		 * @param node
		 * @return number of events
		 */
		private int simulateLineage(int node)
		{
			event_queue.clear();
			TreeWithRates rates = class_rates[class_idx];
			IndexedTree tree = phylo; // rates.getTree(); // same as phylo
			
			int event_count = 0;
			if (tree.isRoot(node))
			{
				
				// simulate at root
				double p = RND.nextDouble();
				double[] cdf;
		        int n = -1;
		        do
		        {
		        	cdf = class_root_cdfs[class_idx];
		            n=Arrays.binarySearch(cdf, p);
		            if (n<0) // binarysearch does not find this value exactly: very likely
		                n = -(n+1);
		            if (n==cdf.length) expandRootCDF(class_idx);
		            	//class_root_cdfs[class_idx] = getRootCDF(rates, 2*(cdf.length-1));
		        } while (n==cdf.length);
//				System.out.println("#**SE.RP.sL root "+p+"\tn "+n+"\t// "+Arrays.toString(cdf));
		        
		        
		        Set<Copy> root_copies = node_copies.get(node);
		        for (int i=0; i<n; i++)
		        {
		        	root_copies.add(new Copy(node));
		        	event_count++;
		        }
			} // for root
			else
			{
				int parent = tree.getParent(node);
				event_queue.clear();
				Copy gain_source;
				double gain_rate = rates.getGainRate(node);
				double duplication_rate = rates.getDuplicationRate(node);
				double loss_rate = rates.getLossRate(node);
				if (duplication_rate == 0.0)
					gain_rate *= loss_rate;
				else
					gain_rate *= duplication_rate;
				if (gain_rate==0.0)
					gain_source = null;
				else
				{
					gain_source = new Copy(node);
					gain_source.setSpawnTime(nextRNDExponential(gain_rate));
					event_queue.add(gain_source); // immortal
				}
				Set<Copy> parent_copies = node_copies.get(parent);
				for (Copy p: parent_copies)
				{
					Copy ancestor = p.inherit(node);
					ancestor.setSpawnTime(nextRNDExponential(duplication_rate));
					ancestor.setDeathTime(nextRNDExponential(loss_rate));
					event_queue.add(ancestor);
				}
//				System.out.println("#**SE.RP.sL "+node+"\tprnt "+parent+"\tq "+event_queue.size());
				
				
				double edge_length = rates.getEdgeLength(node);
				while (!event_queue.isEmpty())
				{
					Copy copy = event_queue.peek();
					if (copy.nextEventTime()>=edge_length)
						break;
					else
						event_queue.poll(); // remove from the event queue 
					event_count++;
					
					if (copy.dies())
					{
						// ok, nothing to do; do not refile this copy 
					} else // a birth!
					{
						Copy newcomer = copy.spawn();
						if (copy == gain_source) // gain event
						{
							copy.setSpawnTime(nextRNDExponential(gain_rate));
						} else
						{
							copy.setSpawnTime(nextRNDExponential(duplication_rate));
							copy.setDeathTime(nextRNDExponential(loss_rate)); // not necessary since same distribution as after spawning
						}
						newcomer.setSpawnTime(nextRNDExponential(duplication_rate));
						newcomer.setDeathTime(nextRNDExponential(loss_rate));
						event_queue.add(newcomer);
						event_queue.add(copy); // refile
					}
				}
				Set<Copy> our_copies = node_copies.get(node);
				// whichever Copy is still in the event queue stays in the history
				for (Copy copy: event_queue)
				{
					if (copy != gain_source)
						our_copies.add(copy);
				}
//				// 
//				if (our_copies.isEmpty())
//				{
//					node_copies.set(node, Collections.EMPTY_SET);
//				} else if (our_copies.size()==1)
//				{
//					node_copies.set(node, Collections.singleton(our_copies.toArray(new Copy[0])[0]));
//				}
				
				
//				System.out.println("#**SE.RP.sL "+node+"\tgot "+our_copies.size());
				
			}
			for (int ci=0; ci<tree.getNumChildren(node); ci++)
			{
				int child = tree.getChild(node, ci);
				event_count += simulateLineage(child);
			}
			
			event_queue.clear(); // tidy up the events beyond the horizon @ edge_length
			
			return event_count;
		}
		
		
		/**
		 * Computes the {@link #surviving_copies} at each node in postorder traversal.
		 * 
		 * @param node
		 */
		private void calculateSurvivors(int node)
		{
			Set<Copy> survivors = surviving_copies.get(node); // to be set 
			Set<Copy> members = node_copies.get(node); // already known 
			if (phylo.isLeaf(node))
			{
				survivors.addAll(members);
			} else
			{
				for (int ci=0; ci<phylo.getNumChildren(node); ci++)
				{
					int child = phylo.getChild(node, ci);
					calculateSurvivors(child);
					for (Copy copy: surviving_copies.get(child))
					{
						if (members.contains(copy.progenitor))
							survivors.add(copy.progenitor);
					}
				}
			}
			if (survivors.isEmpty())
			{
//				System.out.println("#**SE.rP.cS "+node+"\tempty");
				surviving_copies.set(node, Collections.EMPTY_SET);
			}
			else if (survivors.size()==1)
			{
//				System.out.println("#**SE.rP.cS "+node+"\tsingleon");
				surviving_copies.set(node, Collections.singleton(survivors.toArray(new Copy[0])[0]));
			}
		}
		
//	       @Override 
//	       public int hashCode()
//	       {
//	    	   return Arrays.hashCode(profile);
//	       }
//	       
//	       @Override
//	       public boolean equals(Object o)
//	       {
//	    	   if (o!= null && o instanceof PhyleticProfile)
//	    	   {
//	    		   PhyleticProfile other = (PhyleticProfile) o;
//	    		   return Arrays.equals(this.profile, other.profile);
//	    	   } else
//	    		   return super.equals(o);
//	       }
//		
	} // RandomProfile class
	
	
	public static void main(String[] args) throws Exception
	{
		PrintStream out = System.out;
		
		Class<?> our_class = java.lang.invoke.MethodHandles.lookup().lookupClass();
		CommandLine cli = new CommandLine(args, our_class, 2);
		if (cli.getMixedrateModel() == null)
			throw new IllegalArgumentException("Specify the rates model");
		int num_rows;
		if (cli.getTable() == null)
		{
			num_rows = 10;
		} else
		{
			num_rows = cli.getTable().getFamilyCount();
		}
		num_rows = cli.getOptionInt(CommandLine.OPT_N, num_rows);
		int min_obs;
		if (cli.getTable() == null)
		{
			min_obs = 1;
		} else
		{
			min_obs = Integer.min(2,cli.getTable().minCopies());
		}
		min_obs = cli.getOptionInt(OPT_MINCOPY, min_obs);
		
		out.println(CommandLine.getStandardHeader("Families: -"+CommandLine.OPT_N+" "+num_rows));
		out.println(CommandLine.getStandardHeader("Minimum observed: -"+OPT_MINCOPY+" "+min_obs));
		MixedRateModel input_model = cli.getMixedrateModel();
//		if (cli.getFreeModel()!=null)
//		{
//			input_model= cli.getFreeModel();
//			out.println("#SE.main: FreeModel ("+input_model.getNumClasses()+" classes)");
//		}
		
		
		long rndSeed = cli.getOptionLong(OPT_RND, 0L);
		out.println(CommandLine.getStandardHeader("Random initialization: -"+OPT_RND+" "+rndSeed
				+(rndSeed==0L?" (chosen randomly)":"")));    			
		
//		Random RND = cli.getOptionRND(out);
		SimulatedEvolution S = new SimulatedEvolution(input_model, rndSeed);
		
		Table table = S.table(num_rows, min_obs); //  SimulatedEvolution.table(input_model, num_rows, min_obs);
		table.fillTable();

//		out.println(TableParser.getFormattedTable(table, false));
		TableParser.printFormattedTable(out, table, false);
		
		
		// other info on the table 
		boolean want_stats = cli.getOptionBoolean(OPT_STATISTICS, cli.getTable()!=null);
		boolean want_ancestral = cli.getOptionBoolean(OPT_ANCESTRAL, false);
		boolean want_families = cli.getOptionBoolean(OPT_HISTORY, want_ancestral);
		
		
		Phylogeny phylo = cli.getTree();
		
		if (want_stats)
		{
			// by lineages 
			AnnotatedTable input_table = cli.getTable();
			double relative_sample_size = input_table == null?1.0:(num_rows+0.0)/(input_table.getFamilyCount()+0.0);
			
			int[] profile_distribution = new int[phylo.getNumLeaves()+1];
			
			for (int f=0; f<table.getFamilyCount(); f++)
			{
	//			Table.ObservedProfile obs = table.getObservedProfile(f);
				int profile_size = table.getLineageCount(f);
				profile_distribution[profile_size]++;
			}
			int[] template_distribution;
			
			if (cli.getTable()==null)
				template_distribution=null;
			else
			{
				
				template_distribution = new int[profile_distribution.length];
				for (int f=0; f<input_table.getFamilyCount(); f++)
				{
					int size = input_table.getLineageCount(f);
					template_distribution[size]++;
				}
			}
			out.println("#DISTRIBUTION\tsize\tsim"+(template_distribution==null?"":"\tinput")
				+(relative_sample_size==1.0?"":"\textrapolation sim/"+relative_sample_size)
					);
			
			for (int s=0; s<profile_distribution.length; s++)
			{
				out.print("#DISTRIBUTION\t"+s+"\t"
						+(relative_sample_size==1.0?profile_distribution[s]:profile_distribution[s]/relative_sample_size));
				if (template_distribution!=null)
				{
					double diff = profile_distribution[s]/relative_sample_size-template_distribution[s];
					if (template_distribution[s]!=0)
						diff/=template_distribution[s];
					out.printf("\t%d\t%6.3f", template_distribution[s], diff);
				}
				out.println();
			}
			
			{
				// by ncopies
				int max_copies = 0;
				if (input_table != null)
					for (int f=0; f<input_table.getFamilyCount(); f++) {
						int m = input_table.getMemberCount(f);
						if (max_copies<m) max_copies= m;
					}
				for (int f=0; f<table.getFamilyCount(); f++) {
					int m = table.getMemberCount(f);
					if (max_copies<m) max_copies= m;
				}
				int[] profile_ncopies = new int[max_copies+1];
				int[] template_ncopies;;
				for (int f=0; f<table.getFamilyCount(); f++) {
					int m = table.getMemberCount(f);
					profile_ncopies[m]++;
				}
				if (input_table == null) {
					template_ncopies = null;
				} else {
					template_ncopies = new int[profile_ncopies.length];
					for (int f=0; f<input_table.getFamilyCount(); f++) {
						int m = input_table.getMemberCount(f);
						template_ncopies[m]++;
					}
				}
				
				out.println("#NCOPY\tncopy\tsim"+(template_ncopies==null?"":"\tinput")
						+(relative_sample_size==1.0?"":"\textrapolation sim/"+relative_sample_size)
							);
					
				for (int s=0; s<profile_ncopies.length; s++)
				{
					out.print("#NCOPY\t"+s+"\t"
							+(relative_sample_size==1.0?profile_ncopies[s]:profile_ncopies[s]/relative_sample_size));
					if (template_ncopies!=null)
					{
						double diff = profile_ncopies[s]/relative_sample_size-template_ncopies[s];
						if (template_ncopies[s]!=0)
							diff/=template_ncopies[s];
						out.printf("\t%d\t%6.3f", template_ncopies[s], diff);
					}
					out.println();
				}
			}
			
			
			
			
			
			Arrays.fill(profile_distribution, 0);
			for (int f=0; f<table.getFamilyCount(); f++)
			{
				Table.ObservedProfile obs = table.getObservedProfile(f);
	
				for (int leaf=0; leaf<phylo.getNumLeaves(); leaf++)
				{
					int c = obs.getNodeSurvivalCount(leaf, false);
					if (c>0)
						profile_distribution[leaf]++;
				}
				
			}
			if (template_distribution!=null)
			{
				Arrays.fill(template_distribution, 0);
				for (int f=0; f<input_table.getFamilyCount(); f++)
				{
					int[] copies = input_table.getFamilyProfile(f);
					for (int leaf=0; leaf<phylo.getNumLeaves(); leaf++)
					{
						int c = copies[leaf];
						if (c>0)
							template_distribution[leaf]++;
					}
				}
			}
			
			out.println("#FAMILIES\tnode\tsim"+(template_distribution==null?"":"\tinput"));
			
			for (int leaf=0; leaf<phylo.getNumLeaves(); leaf++)
			{
				out.print("#FAMILIES\t"+phylo.getIdent(leaf)+"\t"
						+(relative_sample_size==1.0?profile_distribution[leaf]:profile_distribution[leaf]/relative_sample_size));
				if (template_distribution!=null)
				{
					double diff = profile_distribution[leaf]/relative_sample_size-template_distribution[leaf];
					if (template_distribution[leaf]!=0)
						diff/=template_distribution[leaf];
					out.printf("\t%d\t%.3f", template_distribution[leaf], diff);
				}
				out.println();
			}
		}

		
		if (want_ancestral || want_families)
		{
			String prefix_ancestral = "#"+OPT_ANCESTRAL.toUpperCase();
			String prefix_families = "#"+OPT_HISTORY.toUpperCase();
			
			int n = phylo.getNumNodes();
			// integer counts 
			int[] lineage_family_present = new int[n];
			int[] lineage_family_present_corrected = new int[n];
			int[] lineage_family_multi = new int[n];
			int[] lineage_family_gain = new int[n];
			int[] lineage_family_loss = new int[n];
			int[] lineage_family_expand = new int[n];
			int[] lineage_family_contract = new int[n];
			
			int[] lineage_copies_surviving = new int[n];
			int[] lineage_copies_true = new int[n];
			
			if (want_families)
			{
				out.println(prefix_families+"\tfamily\tnodeidx\tpresence:p\tmulti:m\tcorrected:p.\tgain:g\tloss:l\texpand:++\tcontract:--\tscopies:n\ttcopies:n.");
			}
			
			for (int f=0; f<table.getFamilyCount(); f++)
			{
				Table.ObservedProfile obs = table.getObservedProfile(f);
				for (int node=0; node<n; node++)
				{
					int fp, fpc, fm;
					if (FAMILY_PRESENCE_BY_SURVIVAL) {
						fp = obs.getFamilySurvivalPresent(node, false);
						fpc = obs.getFamilySurvivalPresent(node, true);
						fm = obs.getFamilySurvivalMulti(node, false);
					} else {
						fp = obs.getFamilyPresent(node, false);
						fpc = obs.getFamilyPresent(node, true);
						fm = obs.getFamilyMulti(node, false);
					}
					
					int[] fevents = obs.getFamilyEvents(node, false);
					int fgain = fevents[FamilyEvent.GAIN.ordinal()];
					int floss = fevents[FamilyEvent.LOSS.ordinal()];
					int fexpand =fevents[FamilyEvent.EXPAND.ordinal()];
					int fcontract = fevents[FamilyEvent.CONTRACT.ordinal()];
					
					int fsurviving = obs.getNodeSurvivalCount(node, false);
					int ftruecopies;
					if (FAMILY_PRESENCE_BY_SURVIVAL) 
						ftruecopies = obs.getNodeSurvivalCount(node, true);
					else 
						ftruecopies = obs.getNodeTrueCount(node, true);
					
					if (want_families)
					{
						out.printf("%s\t%d\t%d", prefix_families, f, node);
						out.printf("\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" 
									, fp, fm, fpc
									, fgain, floss, fexpand, fcontract
									, fsurviving, ftruecopies);
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
				} // for node
			} // for family
			if (want_ancestral)
			{
				out.println(prefix_ancestral+"\tnode\tpresence:p\tmulti:m\tcorrected:p.\tgain:g\tloss:l\texpand:++\tcontract:--\tscopies:n\ttcopies:n.");
				for (int node = 0; node<n; node++)
				{
					int fp = lineage_family_present[node];
					int fm = lineage_family_multi[node];
					int fpc = lineage_family_present_corrected[node];
					int fgain = lineage_family_gain[node];
					int floss = lineage_family_loss[node];
					int fexpand = lineage_family_expand[node];
					int fcontract = lineage_family_contract[node];
					
					int fsurviving = lineage_copies_surviving[node];
					int ftruecopies = lineage_copies_true[node];
					
					out.printf("%s\t%s", prefix_ancestral, phylo.getIdent(node));
					out.printf("\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" 
								, fp, fm, fpc
								, fgain, floss, fexpand, fcontract
								, fsurviving, ftruecopies);
				}
			}
		} // want_ancestral or want_families
	}
	
}
