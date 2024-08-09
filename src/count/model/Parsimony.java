package count.model;
/*
 * Copyright 2023 Mikl&oacute;s Cs&#369;r&ouml;s.
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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;

import count.Count;
import count.ds.IndexedTree;
import count.ds.Phylogeny;
import count.ds.Phylogeny.Node;
import count.ds.ProfileTable;
import count.ds.TreeTraversal;
import count.ds.UniqueProfileTable;
import count.io.CommandLine;

import static count.io.CommandLine.OPT_DUPLICATION;
import static count.io.CommandLine.OPT_GAIN;
import static count.io.CommandLine.OPT_LOSS;
import static count.io.CommandLine.OPT_OUTPUT;

/**
 * Parsimony algorithms.
 * 
 * @author csuros
 *
 */
public class Parsimony implements Count.UsesThreadpool //, GLDParameters
{
	public static double DEFAULT_GAIN_PENALTY = 2.0;
	public static double DEFAULT_DUPLICATION_PENALTY = 1.0;
	public static double DEFAULT_LOSS_PENALTY = 1.0;
	
	public Parsimony(IndexedTree phylogeny, ProfileTable table)
	{
		this.tree = phylogeny;
		this.table = table;
		this.post_order = TreeTraversal.postOrder(tree);
		this.utable = (table instanceof UniqueProfileTable)?(UniqueProfileTable)table:null;
	}
	private final IndexedTree tree;
	private final ProfileTable table;
	private final int[] post_order;
	private final UniqueProfileTable utable;
	
//	/**
//	 * Thread pool shared across instances.
//	 */
//	private static final ForkJoinPool pool;
//	static 
//	{
//		if (THREAD_PARALLELISM>1)
//			pool = new ForkJoinPool(THREAD_PARALLELISM);			
//		else
//			pool = null;
//	}

	private static ForkJoinPool thread_pool=null;
	/**
	 * Initializad only once, if {@link Count#THREAD_PARALLELISM} is greater than 1.
	 * 
	 * @return
	 */
	protected synchronized static ForkJoinPool threadPool()
	{
		if (thread_pool == null && 1<Count.THREAD_PARALLELISM) // && Count.THREAD_UNIT_TASK<Integer.MAX_VALUE)
		{
//			System.out.println("#**P.threadPool init: "+Count.THREAD_PARALLELISM+" threads on "+Thread.currentThread());
			thread_pool =  Count.threadPool(); // new ForkJoinPool(Count.THREAD_PARALLELISM);	
		}
		return thread_pool;
	}
	
	
	public Profile getProfile(int family)
	{
		return new Profile(family);
	}
	
	private double gain_penalty=1.0;
	private double loss_penalty = 1.0;
	private double duplication_penalty = 1.0;
	
	public void setPenalties(double gain_penalty, double loss_penalty, double duplication_penalty)
	{
		this.gain_penalty = gain_penalty;
		this.loss_penalty = loss_penalty;
		this.duplication_penalty = duplication_penalty;

	}
	protected double scoreSankoff(int f)
	{
		Profile P = new Profile(f, false);
//		P.computeSankoff(gain_penalty, loss_penalty, duplication_penalty, false);
//		double s = P.getSankoffScore();
		double s = P.scoreSankoff();
		return s;
	}
	
	protected final int getMultiplicity(int f)
	{
		return utable==null?1:utable.getMultiplicity(f);
	}
	
	public double getSankoffScore()
	{
		int nF = table.getFamilyCount();
		ForkJoinPool pool = threadPool();
		final int unit_task = Count.unitTask(nF, Count.THREAD_UNIT_TASK*4);
		
		class PartialScore extends RecursiveTask<Double>
		{
			/**
			 * First family index in {@link Parsimony#utable}.
			 */
			private final int minF;
			/**
			 * Last family index in {@link Parsimony#utable}, exclusive.
			 */
			private final int maxF;
			PartialScore(int min, int max)
			{
				this.minF = min;
				this.maxF = max;
			}
			
			@Override 
			protected Double compute()
			{
				try
				{
					if (maxF-minF>unit_task)
					{
						int midF = (minF+maxF)/2;
						PartialScore left = new PartialScore(minF, midF);
						PartialScore right = new PartialScore(midF, maxF);
						right.fork();
						
						return left.compute()+right.join();
					} else
					{
						double score = 0;
						for (int f=minF; f<maxF; f++)
						{
	//						Profile P = getProfile(f);
	////						P.computeSankoff(gain_penalty, loss_penalty, duplication_penalty, false);
	////						double s = P.getSankoffScore();
	//						double s = P.scoreSankoff(gain_penalty, loss_penalty, duplication_penalty);
	//						if (utable!=null)
	//							s *= utable.getMultiplicity(f);
							double s = scoreSankoff(f);
							s *= getMultiplicity(f);
							score += s;
						}
						return score;
					}
				} catch (Throwable t)
				{
					throw new RuntimeException(t);
				}			
			} // compute()
		} // class
		
		double score;
		PartialScore bigjob = new PartialScore(0,nF);
		
		try
		{
			if (nF > unit_task)
			{
				score = pool.invoke(bigjob);
			} else
			{
				score = bigjob.compute();
			}
		} catch (Throwable t)
		{
			throw new RuntimeException(t);
		}
		return score;
	}
	
//	/**
//	 * Asymmetric Wagner parsimony: 0→1 transitions and n→n+1 transition for 0&lt;n 
//	 * are scored with same <em>birth</em> penalty.
//	 * 
//	 * @param birth_penalty
//	 * @return
//	 */
//	public double getSankoffScore(double birth_penalty)
//	{
//		return getSankoffScore(birth_penalty, 1.0, birth_penalty);
//	}
	
	
	public class Profile
	{
		private Profile(int family)
		{
			this(family, true);
		}
		private Profile(int family, boolean allocate_subtrees)
		{
			this.profile = table.getFamilyProfile(family);
			if (allocate_subtrees)
				setSubtreeCopyNumbers(tree.getRoot(), true);
		}
		private final int[] profile; 
 		private final int[][] subtree_copy_numbers = new int[post_order.length][];
 		private final float[][] sankoff_score = new float[post_order.length][];
 		
 		protected final int[] leafProfile() {return profile;}
 		
 		/**
 		 * Allocates {@link #sankoff_score} and {@link #subtree_copy_numbers},
 		 * and sets {@link #subtree_copy_numbers} at the node (no recursion), 
 		 * or within the subtree (postorder traversal by recursion).
 		 * 
 		 * @param node root of the subtree
 		 * @param recur whether work in subtree or at node 
 		 */
 		private void setSubtreeCopyNumbers(int node, boolean recur)
 		{
 			int[] possible_values;
 			if (tree.isLeaf(node))
 			{
 				if (profile[node]<0)
 				{
 					possible_values = new int[0];
 				} else
 				{
 					possible_values = new int[1];
 					possible_values[0] = profile[node];
 				}
 			} else
 			{
 				int num_children = tree.getNumChildren(node);
 				Set<Integer> union_children = new HashSet<>();
 				for (int ci=0; ci<num_children; ci++)
 				{
 					int child = tree.getChild(node, ci);
 					if (recur) 
 					{
 						setSubtreeCopyNumbers(child, recur);
 					}
 					for (int n: subtree_copy_numbers[child])
 						union_children.add(n);
 				}
 				if (!union_children.isEmpty())
 				{
 					union_children.add(0);
 					union_children.add(1);
 				}
 				possible_values = new int[union_children.size()];
 				int si=0;
 				for (Integer n: union_children)
 				{
 					possible_values[si++] = n;
 				}
 				assert (si == possible_values.length);
 				Arrays.sort(possible_values);
 			}
 			subtree_copy_numbers[node] = possible_values;
 			sankoff_score[node] = new float[possible_values.length];
 		}
 		

 		
 		
// 		public int[] computeSankoff(double birth_penalty, boolean max_copies)
// 		{
// 			return computeSankoff(birth_penalty, 1.0, birth_penalty, max_copies);
// 		}
 		
 		/**
 		 * Call after {@link #computeSankoff(boolean)}
 		 * to retrieve the cached min score at the root.
 		 * @return
 		 */
 		public double getSankoffScore()
 		{
 			return cached_sankoff_score;
 		}
 			
 		/**
 		 * Set by {@link #scoreSankoff()} an {@link #computeSankoff(boolean)}
 		 * methods.
 		 */
 	    private double cached_sankoff_score=Double.NaN; 
 	    
 	    private void addMinScores(int child)
		{
 	    	
 	    	assert (!tree.isRoot(child));
			final int[] child_values = subtree_copy_numbers[child];
			if (child_values.length==0) return; // nothing to do; all ambiguous hence no penalty
			
 	    	int node = tree.getParent(child);
			final float[] penalties = sankoff_score[node];
			final int[] values = subtree_copy_numbers[node];
			assert (penalties.length==values.length);

// 	    	System.out.println("#**P.P.aMS ch "+child+"\tn "+node+"\t"+child_values.length+"/"+values.length);
			
			if (tree.getLength(child)==Double.POSITIVE_INFINITY)
			{
				float[] child_penalties = sankoff_score[child];
				int pv = 0; // penalized for 0-> transitions, no matter what the actual parent copy number is
				double min_score=Double.POSITIVE_INFINITY;
				for (int cvidx=0; cvidx<child_values.length; cvidx++)
				{
					int cv = child_values[cvidx];
					double change_pty = getSankoffPenalty(pv,cv,gain_penalty,loss_penalty,duplication_penalty);
					double subtree_score = change_pty + child_penalties[cvidx];
					if (subtree_score<min_score)
					{
						min_score = subtree_score;
					}
				}
				for (int i=0; i<penalties.length; i++)
					penalties[i] += min_score;
			} else
			{
				float[] child_penalties = sankoff_score[child];
				
				// track left- and right minimal indices in child_XXX: 
				// jl is left-minimal if for all j<=jl, pty[jl] <= change(x_jl->x_j)+pty[j]
				// jr is right-minimal if for all jr<=j, pty[jr] <= change(x_jr->x_j)+pty(j)
				// I. store right-minimal indices in a right-to-left scan in child
				// II. in a left-to-right scan  (i) in node_XXX
				//  II.1 update jl and jr when matching a child_XXX[j]
				//  II.2 minimum is either change(y_i->x_jl)+pty[jl], or change(y_i->x_jr)+pty[jr]
				//    where jl is the last left minimum with x_jl <=y_i 
				//        and jr is the first right minimum with y_i<x_jr

				int j = child_values.length-1;
				int jr = j;
				int[] right_minimum = new int[child_values.length];
				right_minimum[j]=jr;
				double ptyr = child_penalties[jr];
				while (j>0)
				{
					--j;
					int cv = child_values[j];
					double ptyj = child_penalties[j];
					if (ptyj<ptyr 
						|| ptyj < ptyr+getSankoffPenalty(cv, child_values[jr], gain_penalty, loss_penalty, duplication_penalty))
					{
						jr = j;
						ptyr = ptyj;
					}
					right_minimum[j]=jr;
				}
				assert (j==0); // can start the left-to-right scan 
				
				int jl = -1;
				double ptyl = Double.POSITIVE_INFINITY;
				
				int cv = child_values[j];
				for (int i=0; i<values.length; i++)
				{
					int pv = values[i];
					if (pv==cv)
					{
						double ptyj = child_penalties[j];
						if (ptyj<ptyl // kicks in with ptyl==infty
							|| ptyj<ptyl
								+getSankoffPenalty(cv, child_values[jl],gain_penalty, loss_penalty, duplication_penalty))
						{ // new left minimum
							jl=j;
							ptyl = ptyj;
						}
						j++;
						if (j<child_values.length)
						{
							cv = child_values[j];
							// recover right minimum for this position
							jr = right_minimum[j];
							ptyr = child_penalties[jr];
						} else
						{
							cv = values[values.length-1]+1; // always pv < cv from now on
							jr = child_values.length;
							ptyr = Double.POSITIVE_INFINITY;
						}
					} 
					assert (pv<cv); 
					double min_score = Double.POSITIVE_INFINITY;
					assert Double.isFinite(ptyl) || Double.isFinite(ptyr);
					if (ptyl < min_score) // not when jl==-1 and ptyl==infty
					{
						// II.2a use left-minimum at or before pv
						double change_pty = getSankoffPenalty(pv, child_values[jl], gain_penalty, loss_penalty, duplication_penalty);
						min_score = change_pty+ptyl; // Double.min(change_pty + ptyl, min_score);
					}
					if (ptyr < min_score) // not when jr==m and ptyr==infty
					{
						// II.2b or use right-minimum after pv
						//    jl==jr only when it is a global minimum 
						double change_pty =  getSankoffPenalty(pv, child_values[jr], gain_penalty, loss_penalty, duplication_penalty);
						min_score = Double.min(change_pty + ptyr, min_score);
					}
					assert Double.isFinite(min_score);
					penalties[i] += min_score;
				} // for i 
			}
 	    }
 	     	    
 	    /**
 	     * Allocates and sets {@link #sankoff_score} and {@link #subtree_copy_numbers}
 	     * at a node; to be called in postorder. 
 	     * 
 	     * @param node
// 	     * @param gain_penalty 0→1 copy number change
// 	     * @param loss_penalty 1→0 copy number change  
// 	     * @param duplication_penalty n→n+1 copy number change at 0&lt;n
 	     */
 	    private void scoreSankoff(int node)
 	    {
 	    	this.setSubtreeCopyNumbers(node, false);
			final float[] penalties = sankoff_score[node];
			Arrays.fill(penalties, 0.0f);
			
			if (tree.isLeaf(node))
			{
				final int[] values = subtree_copy_numbers[node];
				
				assert (penalties.length==values.length);
				if (profile[node]<0)
				{
					// nothing to do
					assert (values.length==0);
				} else
				{
					penalties[0]=0.0f; // the single possible value there
				}
			} else
			{
				final int[] values = subtree_copy_numbers[node];
				
				assert (penalties.length==values.length);
				int num_children = tree.getNumChildren(node);
				for (int ci=0; ci<num_children; ci++)
				{
					int child = tree.getChild(node, ci);
					int[] child_values = subtree_copy_numbers[child];
					if (child_values.length==0)
					{
						// all ambiguous, nothing to do 
					} else if (tree.getLength(child)==Double.POSITIVE_INFINITY
							|| child_values.length>4)
					{
						addMinScores(child);
					} else
					{
						// quadratic time with embedded loops, only when not too many values at child
						for (int vidx=0; vidx<values.length; vidx++)
						{ 
							int pv = values[vidx];
							float[] child_penalties = sankoff_score[child];
	
							double min_score=Double.POSITIVE_INFINITY;
							for (int cvidx=0; cvidx<child_values.length; cvidx++)
							{
								int cv = child_values[cvidx];
								
								double change_pty = getSankoffPenalty(pv,cv,gain_penalty,loss_penalty,duplication_penalty);
								double subtree_score = change_pty + child_penalties[cvidx];
								if (subtree_score<min_score)
								{
									min_score = subtree_score;
								}
							}
							assert (Double.isFinite(min_score));
							penalties[vidx]+=min_score;
						} // for all values at node
					} // if 
				} // for all children
			} // if internal node
 	    }
 	    
 	    private double scoreSankoffRoot()
 	    {
			int root = tree.getRoot();

			double min_score;
			if (sankoff_score[root].length==0)
			{
				min_score = 0;
			} else
			{
				min_score  = Double.POSITIVE_INFINITY;
				for (int vidx=0; vidx<subtree_copy_numbers[root].length; vidx++)
				{
					int rv = subtree_copy_numbers[root][vidx];
					double change_pty = getSankoffPenalty(0,rv,gain_penalty,loss_penalty,duplication_penalty); // root copy number: 0→rv
					double score = change_pty+sankoff_score[root][vidx];
					min_score = Double.min(min_score,  score);
				}
			}
			cached_sankoff_score = min_score;
			return min_score; 	    	
 	    }
 	    
 	    /**
 	     * Calculates the Sankoff parsimony score without reconstructing the history.
 	     * @return
 	     */
 	    public double scoreSankoff()
 	    {
			for (int node: post_order)
				scoreSankoff(node);
			return cached_sankoff_score = scoreSankoffRoot();
//			
//			int root = tree.getRoot();
//
//			double min_score;
//			if (sankoff_score[root].length==0)
//			{
//				min_score = 0;
//			} else
//			{
//				min_score  = Double.POSITIVE_INFINITY;
//				for (int vidx=0; vidx<subtree_copy_numbers[root].length; vidx++)
//				{
//					int rv = subtree_copy_numbers[root][vidx];
//					double change_pty = getSankoffPenalty(0,rv,gain_penalty,loss_penalty,duplication_penalty); // root copy number: 0→rv
//					double score = change_pty+sankoff_score[root][vidx];
//					min_score = Double.min(min_score,  score);
//				}
//			}
//			cached_sankoff_score = min_score;
//			return min_score;
 	    }
 	    
 	    /**
 	     * Calculates the history with minimum parsimony 
 	     * 
  	     * @param max_copies favor max (or min if false) copy number among solutions with equal score 
 	     * @return
 	     */
 		public int[] computeSankoff(final boolean max_copies)
 		{
			for (int node: post_order)
	        {
				scoreSankoff(node);
//				final double[] penalties = sankoff_score[node];
//				Arrays.fill(penalties, 0.0);
//				final int[] values = subtree_copy_numbers[node];
//				
//				assert (penalties.length==values.length);
//				
//				if (tree.isLeaf(node))
//				{
//					if (profile[node]<0)
//					{
//						// nothing to do
//						assert (values.length==0);
//					} else
//					{
//						penalties[0]=0.0;
//					}
//				} else
//				{
//					int num_children = tree.getNumChildren(node);
//					for (int vidx=0; vidx<values.length; vidx++)
//					{
//						int pv = values[vidx];
//						for (int ci=0; ci<num_children; ci++)
//						{
//							int child = tree.getChild(node, ci);
//							int[] child_values = subtree_copy_numbers[child];
//							if (child_values.length>0)
//							{
//								double[] child_penalties = sankoff_score[child];
//								int min_cvidx=-1;
//								double min_score=Double.POSITIVE_INFINITY;
//								for (int cvidx=0; cvidx<child_values.length; cvidx++)
//								{
//									int cv = child_values[cvidx];
//									double change_pty = getSankoffPenalty(pv,cv,gain_penalty,loss_penalty,duplication_penalty);
//									double subtree_score = change_pty + child_penalties[cvidx];
//									if (subtree_score<min_score || (max_copies && subtree_score==min_score))
//									{
//										min_cvidx = cvidx;
//										min_score = subtree_score;
//									}
//								}
//								penalties[vidx]+=min_score;
//							} // child is not ambiguous
//						} // for all children
//					} // for all values at node
//				} // if internal node
	        } // all nodes in postorder
			
			
			
			int[] copies = new int[post_order.length];
			int root = tree.getRoot();

			if (sankoff_score[root].length==0)
			{
				copies[root]=0;
			} else
			{
				double min_score = Double.POSITIVE_INFINITY;
				for (int vidx=0; vidx<subtree_copy_numbers[root].length; vidx++)
				{
					int rv = subtree_copy_numbers[root][vidx];
					double change_pty = getSankoffPenalty(0,rv,gain_penalty,loss_penalty,duplication_penalty);
					double score = change_pty+sankoff_score[root][vidx];
					if (score<min_score || (max_copies && score==min_score))
					{
						copies[root]=rv;
						min_score = score;
					}
				}
				cached_sankoff_score = min_score;
			}
			int node_idx = copies.length-1;
			assert (root==post_order[node_idx]);

			while (node_idx>0)
			{
				--node_idx;
				int node = post_order[node_idx];
				int parent = tree.getParent(node);
				int pv = (tree.getLength(node)==Double.POSITIVE_INFINITY)?0:copies[parent];
				if (subtree_copy_numbers[node].length==0)
				{
					copies[node]=pv;
				} else
				{
					double min_score = Double.POSITIVE_INFINITY;
					for (int cvidx=0; cvidx<subtree_copy_numbers[node].length; cvidx++)
					{
						int cv = subtree_copy_numbers[node][cvidx];
						double change_pty = getSankoffPenalty(pv,cv,gain_penalty,loss_penalty,duplication_penalty);
						double score = change_pty + sankoff_score[node][cvidx];
						if (score<min_score || (max_copies && score ==min_score))
						{
							min_score=score;
							copies[node]=cv;
						}
					}
				}
			}
			return copies;
 		}

 		private int cached_dollo_score = -1;
 		
 		public int[] computeDollo(boolean root_surely_present)
 		{
 			int[] history = new int[post_order.length];
 			int soft1 = 2; // different from 0 and 1
 			for (int node: post_order)
 			{
 				if (tree.isLeaf(node))
 				{
 					history[node]=profile[node]>0?1:0;
 				} else
 				{
 					int nhist = 0;
 					int num_children = tree.getNumChildren(node);
 					for (int ci=0; ci<num_children && nhist != 1 ; ci++)
 					{
 						int ch = history[tree.getChild(node, ci)];
 						if (ch == 1 || ch==soft1)
 						{
 							if (nhist==soft1)
 								nhist=1;
 							else
 							{
 								assert (nhist==0);
 								nhist=soft1;
 							}
 						}
 					}
 					history[node]=nhist;
 				}
 			}
 			
 			int num_losses = 0;
 			int num_gains = 0;
 			int root = tree.getRoot();
 			int rhist = history[root];
 			if (rhist == soft1)
 				rhist = history[root] = root_surely_present?1:0;

 			num_gains += rhist;
 			
 			int node_idx = root;
 			assert (node_idx==post_order.length-1);
 			while (node_idx>0)
 			{
 				--node_idx;
 				int node = post_order[node_idx];
 				int nhist = history[node];
				int parent = tree.getParent(node);
				int phist = history[parent];
 				if (nhist==soft1)
 				{
 					assert (phist == 0 || phist==1); 
 					nhist = history[node] = (phist==1)?1:0;
 				}
 				if (phist==1 && nhist==0)
 					num_losses++;
 				if (phist==0 && nhist==1)
 					num_gains++;
 			}
 			assert (num_gains<2);
 			this.cached_dollo_score = num_losses+num_gains;
 			
 			return history;
 		}
 		
 		/**
 		 * Number of losses in the most recent call 
 		 * to {@link #computeDollo(boolean)}. 
 		 * 
 		 * @return
 		 */
 		public int getDolloScore()
 		{
 			return cached_dollo_score;
 		}
	
	} // Profile
	
	/**
	 * Calculates the penalty for pv→cv transition (parent to child)
	 * 
	 * @param pv
	 * @param cv
	 * @param gain_penalty
	 * @param loss_penalty
	 * @param duplication_penalty
	 * @return
	 */
	private static float getSankoffPenalty(int pv, int cv, double gain_penalty, double loss_penalty, double duplication_penalty)
	{
		double change_pty;
		if (cv==pv)
		{
			change_pty=0.0; 
		} else if (cv<pv)
		{
			if (cv==0)
			{
				change_pty = loss_penalty + (pv-cv-1.0);
			} else
			{
				change_pty = (pv-cv); // so many losses
			}
		} else 
		{
			// assert (pv<cv)
			if (pv==0)
			{
				change_pty = gain_penalty + duplication_penalty*(cv-pv-1.0);
			} else
			{
				change_pty = duplication_penalty*(cv-pv);
			}
		}
		return (float)change_pty;
	}
	
	
//	private class ChangePenalty 
//	{
//		private double gain_pty;
//		private double loss_pty;
//		private double dup_pty;
//		ChangePenalty(double gain_penalty, double loss_penalty, double duplication_penalty)
//		{
//			setPenalties(gain_penalty, loss_penalty, duplication_penalty);
//		}
//		
//		void setPenalties(double gain_penalty, double loss_penalty, double duplication_penalty)
//		{
//			this.gain_pty = gain_penalty;
//			this.loss_pty = loss_penalty;
//			this.dup_pty = duplication_penalty;
//			this.change_pty = new double[2][];
//		}
//		
//		private double[][] change_pty;
//		
//		double get(int pv, int cv)
//		{
//			if (change_pty.length<=pv)
//				change_pty = Arrays.copyOf(change_pty, pv+1);
//			
//			
//			
//		}
//	}
//	
	public class WagnerProfile extends Profile
	{
		WagnerProfile(int f)
		{
			super(f, true);
		}
		
		/* 
		 * Wagner parsimony work variables, preallocated and reused with 
		 * variable gain penalty calls of {@link #computeWagner}. 
		 */
        private final double[][] subtree_slope = new double[post_order.length][];
        private final int[][] subtree_breakpoint = new int[post_order.length][];
        private double[] subtree_shift = new double[post_order.length];
        private double[][] stem_slope = new double[post_order.length][];
        private int[][] stem_breakpoint = new int[post_order.length][]; 
        private double[] stem_shift = new double[post_order.length];
        private int[] stem_left_bp = new int[post_order.length];
        private int[] stem_right_bp = new int[post_order.length];

    	/**
    	 * Algorithm for 
    	 * summing two piecewise linears, returning the result in the first pair of arguments;
    	 * uses no instance-linked variables. (Used in Wagner parsimony.)
    	 * 
    	 * @param slopes repopulated with result value
    	 * @param breakpoints repopulated with result
    	 * @param slopes2
    	 * @param breakpoints2
    	 */
        private void sumPiecewiseLinear(List<Double> slopes, List<Integer> breakpoints, double[] slopes2, int[] breakpoints2)
        {
            if (slopes.size()==0)
            { // first call
                for (int i=0; i<slopes2.length; i++)
                {
                    slopes.add(slopes2[i]);
                    breakpoints.add(breakpoints2[i]);
                }
            } else
            {
                int n1 = slopes.size();
                double[] slopes1 = new double[n1];
                int[] breakpoints1 = new int[n1];
                for (int i=0; i<n1; i++)
                {
                    slopes1[i] = slopes.get(i);
                    breakpoints1[i] = breakpoints.get(i);
                }
                slopes.clear();
                breakpoints.clear();

                slopes.add(slopes1[0]+slopes2[0]);
                breakpoints.add(0); // dummy placeholder
                int n2 = slopes2.length;
                
                int i1=1; 
                int i2=1;
                while (i1<n1 || i2<n2)
                {
                    if (i1==n1 || (i2<n2 && breakpoints2[i2]<breakpoints1[i1]))
                    {
                        // add breakpoint2
                        int x = breakpoints2[i2];
                        double a = slopes1[i1-1]+slopes2[i2];
                        breakpoints.add(x);
                        slopes.add(a);
                        i2++;
                    } else if (i2==n2 || (i1<n1 && breakpoints1[i1]<breakpoints2[i2]))
                    {
                        // add breakpoint1
                        int x = breakpoints1[i1];
                        double a = slopes1[i1]+slopes2[i2-1];
                        breakpoints.add(x);
                        slopes.add(a);
                        i1++;
                    } else 
                    { // equality
                        int x = breakpoints1[i1];
                        double a = slopes1[i1]+slopes2[i2];
                        breakpoints.add(x);
                        slopes.add(a);
                        i1++;
                        i2++;
                    }
                } // while 
            } // if not first call
        } // piecewise linear
        
        
        /**
         * Returns the array of parsimony-optimal 
         * copy number labeling on the underlying phylogeny.
         * 
         * @param birth_penalty
         * @return
         */
 		public int[] computeWagner(final double birth_penalty)
		{
 			int[] profile = leafProfile();
	        Arrays.fill(subtree_shift, 0.0);
	        Arrays.fill(stem_shift, 0.0);
	        Arrays.fill(stem_left_bp, 0);
	        Arrays.fill(stem_right_bp, 0);
			for (int node: post_order)
	        {
	            if (tree.isLeaf(node))
	            {
	                // compute stem weight function
	                if (profile[node]<0) // ambiguous
	                {
	                    stem_slope[node]=new double[1];
	                    stem_breakpoint[node]=new int[1];
	                    stem_slope[node][0]=0.0;
	//                    stem_slope[node_idx][1]=0.0;
	//                    stem_breakpoint[node_idx][1]=0;
	                    stem_shift[node]=0.0;
	                } else
	                {
	                    stem_slope[node]=new double[2];
	                    stem_breakpoint[node]=new int[2];
	                    stem_slope[node][0]=-birth_penalty;
	                    stem_slope[node][1]=1.0;
	
	                    stem_breakpoint[node][1]=profile[node];
	                    stem_shift[node]=birth_penalty*profile[node];
	                }
	                //Verbose.message("PP.cWP leaf "+node_idx+"/"+N.getTaxonName()+"\t"+profile[node_idx]);
	            } else
	            {
	                // compute subtree weight functions
	                int num_children = tree.getNumChildren(node);
	                List<Double> slopeV = new ArrayList<>();
	                List<Integer> breakpointV = new ArrayList<>();
	                subtree_shift[node]=0.0;
	                
	                for (int ci=0; ci<num_children; ci++)
	                {
	                    int child = tree.getChild(node,ci);
	                    sumPiecewiseLinear(slopeV, breakpointV, stem_slope[child], stem_breakpoint[child]);
	                    subtree_shift[node]+=stem_shift[node];
	                }
	                int k = slopeV.size();
	                subtree_slope[node]=new double[k];
	                subtree_breakpoint[node]=new int[k];
	                for (int i=0; i<k; i++)
	                {
	                    subtree_slope[node][i] = slopeV.get(i).doubleValue();
	                    subtree_breakpoint[node][i] = breakpointV.get(i).intValue();
	                }
	                if (!tree.isRoot(node))
	                {
	                    // compute stem weight function
	                    int i_left = 0;
	                    for (int i=1; i<k; i++)
	                        if (subtree_slope[node][i]>=-birth_penalty)
	                        {
	                            i_left=i;
	                            break;
	                        }
	                    
	                    // what is the shift here?
	                    double phi = subtree_shift[node];
	                    for (int i=1; i<=i_left; i++)
	                    {
	                        if (i==1)
	                            phi += subtree_slope[node][0]*subtree_breakpoint[node][1];
	                        else
	                            phi += subtree_slope[node][i-1]*(subtree_breakpoint[node][i]-subtree_breakpoint[node][i-1]);
	                    }
	                    stem_shift[node] = phi + birth_penalty * subtree_breakpoint[node][i_left];
	                    
	                    int i_right = k-1;
	                    while (i_right>=i_left && subtree_slope[node][i_right]>=1.0)
	                        i_right--;
	                    i_right++;
	                    //Verbose.message("PP.cWP stem @ "+i_left+".."+i_right+" ["+k+"] shift "+stem_shift[node_idx]);
	                    stem_left_bp[node]=i_left;
	                    stem_right_bp[node]=i_right;
	                    
	                    stem_slope[node]=new double[i_right-i_left+2];
	                    stem_breakpoint[node]=new int[i_right-i_left+2];
	                    stem_slope[node][0]=-birth_penalty;
	                    for (int i=i_left; i<=i_right; i++)
	                    {
	                        stem_slope[node][i-i_left+1]=subtree_slope[node][i];
	                        stem_breakpoint[node][i-i_left+1]=subtree_breakpoint[node][i];
	                    }
	                    stem_slope[node][i_right-i_left+1]=1.0;
	                    //if (Verbose.isVerbose())
	                    //{
	                    //   for (int i=0; i<stem_slope[node_idx].length; i++)
	                    //   {
	                    //       Verbose.message("PP.cWP stem "+node_idx+"/"+N.getTaxonName()+"\t"+i+"\t"+stem_slope[node_idx][i]+"\t"+stem_breakpoint[node_idx][i]);
	                    //   }
	                    //}
	
	                }
	            }
	        } // for all nodes
	        int[] retval = new int[post_order.length];
	        
	        // find minimum at root
	        int imin=1;
	        int root = tree.getRoot();
	        assert (root == post_order.length-1);
	        
	        while (subtree_slope[root][imin]<0) imin++;
	        retval[root] = subtree_breakpoint[root][imin];

	        int node_idx=post_order.length-1;
	        assert (root==post_order[node_idx]);
	        while(node_idx>0)
	        {
	        	--node_idx;
	        	int node = post_order[node_idx];
	            int parent= tree.getParent(node);
	            int y = retval[parent];
	            
	            if (tree.isLeaf(node))
	            {
	                if (profile[node]<0)
	                {
	                    retval[node] = retval[parent];
	                } else
	                {
	                    retval[node]=profile[node];
	                }
	            }
	            else
	            {
	                int x0 = subtree_breakpoint[node][stem_left_bp[node]];
	                int x1 = subtree_breakpoint[node][stem_right_bp[node]];
	                if (y<x0)
	                    retval[node]=x0;
	                else if (y>x1)
	                    retval[node]=x1;
	                else
	                    retval[node]=y;
	                //Verbose.message("PP.cWP solution "+node_idx+"/"+N.getTaxonName()+"\t"+retval[node_idx]+"\tprn "+y+"\t["+x0+", "+x1+"]");
	            }
	            
	        } // all nodes 
	        
	        return retval;
	    }
 		
		
		
	} // WagnerProfile
	
	
//	/**
//	 * Parsimony penalties for {@link setPenalties(double, double, double)}
//	 * calculated by fitting to log-odds scores from the argument.
//	 * 
//	 * @return [gain, loss/death, dup]
//	 */
//	private double[] fitParsimonyPenalty(double[][] log_trans)
//	{
//		// calcukate log-odds scores
//		double[][] w=new double[log_trans.length][];
//		for (int n=0; n<w.length; n++)
//		{
//			w[n]=new double[log_trans[n].length];
//			for (int m=0; m<w[n].length; m++)
//			{
//				w[n][m] = log_trans[n][n]-log_trans[n][m];
////				System.out.println("#**DEM.fPP n "+n+"\tm "+m+"\tw "+w[n][m]);
//			}
//		}
//		
//		double avg_m=0.0;
//		double avg_w0m = 0.0;
//		//double tsum = 0.0;
//		
//		double log_tsum = Double.NEGATIVE_INFINITY;
//		
//		
//		{
//			int m=w[0].length;
//			while (m>1)
//			{
//				
//			}
//		}
//		for (int m=1; m<w[0].length; m++)
//		{
//			double t = log_trans[0][m];
//			if (t!=Double.NEGATIVE_INFINITY)
//			{
//				log_tsum = Logarithms.add(log_tsum, t);
//				double p = Math.exp(t);
//				
//				
//				avg_m += p*m;
//				avg_w0m += p*w[0][m];
//			}
//		}
//		
//		
////		System.out.println("#**DEM.fPP summ "+avg_m+"\tsumw "+avg_w0m+"\ttsum "+tsum);
//		
//		avg_m /= tsum;
//		avg_w0m /= tsum;
//		
//		double inc_num = 0.0;
//		double inc_denom = 0.0;
//		for (int m=1; m<w[0].length; m++)
//		{
//			double t = trans[0][m];
//			if (t!=0.0)
//			{
//				double dm = (m-avg_m);
//				inc_num += t*w[0][m]*dm;
//				inc_denom += t*dm*dm;
//			}
//		}		
//		for (int n=1; n<w.length; n++)
//		{
//			for (int d=1, m=n+d; m<w[n].length; m++, d++)
//			{
//				double t = trans[n][m];
//				if (t!=0.0)
//				{
//					inc_num += t*w[n][m]*d;
//					inc_denom += t*d*d;
//				}
//			}
//		}
//		double opt_inc = inc_num/inc_denom;
//		double opt_gain = avg_w0m-opt_inc*avg_m;
////		System.out.println("#**DEM.fPP inc_num "+inc_num+"\tinc_den "+inc_denom
////				+"\topt_inc "+opt_inc+"\topt_gain "+opt_gain);
//		
//		tsum=0.0;
//		double avg_n=0.0;
//		double avg_wn0 = 0.0;
//		for (int n=1; n<w.length; n++)
//		{
//			double t = trans[n][0];
//			if (t!=0.0)
//			{
//				tsum += t;
//				avg_n += t*n;
//				avg_wn0 += t*w[n][0];
//			}
//		}
////		System.out.println("#**DEM.fPP sumn "+avg_n+"\tsumw "+avg_wn0+"\ttsum "+tsum);
//
//		avg_n /= tsum;
//		avg_wn0 /= tsum;
//		
//		double dec_num=0.0;
//		double dec_denom=0.0;
//		for (int n=1; n<w.length; n++)
//		{
//			double t = trans[n][0];
//			if (t!=0.0)
//			{
//				double dn = n-avg_n;
//				dec_num += t*w[n][0]*dn;
//				dec_denom += t*dn*dn;
//			}
//		}
//		for (int n=2; n<w.length; n++)
//		{
//			for (int d=1, m=n-d; m>0; d++, m--)
//			{
//				double t = trans[n][m];
//				if (t!=0.0)
//				{
//					dec_num += t*w[n][m]*d;
//					dec_denom += t*d*d;
//				}
//			}
//		}
//
//		
//		double opt_dec = dec_num/dec_denom;
//		double opt_loss = avg_wn0-opt_dec*avg_n;
////		System.out.println("#**DEM.fPP dec_num "+dec_num+"\tdec_den "+dec_denom
////				+"\topt_dec "+opt_dec+"\topt_loss "+opt_loss);
//		
//		double pty_dup = opt_inc/opt_dec;
//		double pty_gain = (opt_gain+opt_inc)/opt_dec;
//		double pty_loss = (opt_loss+opt_dec)/opt_dec;
//		
//		
//		Count.out.println("#**DEM.fitP gain "+pty_gain+"\tloss "+pty_loss+"\tdup "+pty_dup);
//		double[] pty = new double[3];
//		pty[PARAMETER_GAIN] = pty_gain;
//		pty[PARAMETER_LOSS] = pty_loss;
//		pty[PARAMETER_DUPLICATION] = pty_dup;
//		
////		for (int n=0; n<w.length; n++)
////		{
////			for (int m=0; m<w[n].length && m<n; m++)
////			{
////				double c =   (n-m)*opt_dec;
////				double v = (n-m)*1.0;
////				
////				if (m==0)
////				{
////					c += opt_loss;
////					v += pty_loss-1.0;
////				}
////				double e = (w[n][m]-c)/w[n][m];
////				System.out.println("#**DEM.fitP\t"+n+"\t"+m+"\t"+trans[n][m]+"\t"+w[n][m]+"\t"+c+"\t"+e+"\t"+(v*opt_dec)+"\t"+v);
////			}
////			for (int m=n+1; m<w[n].length; m++)
////			{
////				double c = (m-n)*opt_inc;
////				double v = (m-n)*pty_dup;
////				if (n==0)
////				{
////					c += opt_gain;
////					v += pty_gain-pty_dup;
////				}
////
////				double e = (w[n][m]-c)/w[n][m];
////				System.out.println("#**DEM.fitP\t"+n+"\t"+m+"\t"+trans[n][m]+"\t"+w[n][m]
////						+"\t"+c+"\t"+e+"\t"+(v*opt_dec)+"\t"+v);
////			}
////		}
//		
//		
//		return pty;
//	}
		
	
	public static class SPRExplorer
	{
		public SPRExplorer(Phylogeny tree, UniqueProfileTable table, double gain_penalty, double loss_penalty, double duplication_penalty)
		{
			this.phylo = tree;
			this.og_profiles = new Profile[table.getFamilyCount()];
			this.factory = new Parsimony(phylo, table)
					{
						// saves the profiles in og_profiles
						@Override
						protected double scoreSankoff(int f)
						{
							Profile P = getProfile(f);
							og_profiles[f] = P;
							double s = P.scoreSankoff();
							return s;
						}
					};
			factory.setPenalties(gain_penalty, loss_penalty, duplication_penalty);
			this.factory_score = factory.getSankoffScore(); // fills the og_profiles[] table
		}
		private final Parsimony factory;
		private final Profile[] og_profiles;
		private final Phylogeny phylo;
		private final double factory_score; 
		
		public double getSankoffScore()
		{
			return factory_score;
		}
		
		public double getGainPenalty()
		{
			return factory.gain_penalty;
		}
		
		public double getLossPenalty()
		{
			return factory.loss_penalty;
		}
		
		public double getDuplicationPenalty()
		{
			return factory.duplication_penalty;
		}
		
		private Phylogeny last_spr;
		
		/**
		 * Parsimony score change for a subtree-prune-and-regraft 
		 * operation.  
		 * 
		 * @param graft_node graft node (new sibling for prune node), original parent becomes grandparet
		 * @param prune_node pruned node, detached from original parent, and joins graft node under a new common parent node 
		 * @return
		 */
		public double getDiffScore(int graft_node, int prune_node)
		{
			Phylogeny spr = new Phylogeny(phylo);
			Map<Node, Integer> og_index = new HashMap<>();
			for (Node N: spr.getRootNode().listNodes(null, null))
			{
				og_index.put(N, N.getIndex());
			}
			Node P = spr.getNode(prune_node);
			assert (!P.isRoot());
			assert (P.getParent().getNumChildren()>=2);
			Node Psib;
			Psib = P.getParent().getChild(1);
			if (Psib==P) Psib = P.getParent().getChild(0);
			
			Node G = spr.getNode(graft_node);
			
			spr.pruneAndRegraft(P, G, false);
			final Node U = G.getParent();
			final Node V = Psib.isRoot()?Psib:Psib.getParent();
			assert (U==P.getParent()); // since that was the spr move
			
//			System.out.println("#**P.SPR.gDS prune "+prune_node+"/"+P
//					+"\tregraft "+graft_node+"/"+G);
			
			/**
			 * Class for carrying out the Sankoff recursions
			 * only over two nodes (U, V) and their ancestors.
			 * 
			 * @author csuros
			 *
			 */
			class Update extends Parsimony
			{
				Update()
				{
					super(spr, factory.utable.mappedToTree(spr));
					this.setPenalties(factory.gain_penalty, factory.loss_penalty, factory.duplication_penalty);
					assert (!U.isLeaf());
					assert (!V.isLeaf());
					initDataStructures();
				}
				
				Node[] uv_ancestors;
				Node[] uv_children;
				
				private void initDataStructures()
				{
//					System.out.println("#**P.SPR.gDS/U.iDS U "+U+"\tV "+V);
					
					Set<Node> ancestors = new HashSet<>();
					ancestors.add(U);
					ancestors.add(V);
					ancestors.addAll(U.listAncestors(null));
					ancestors.addAll(V.listAncestors(null));
					Node[] empty = new Node[0];
					uv_ancestors = ancestors.toArray(empty);
					Arrays.sort(uv_ancestors, (a,b)->Integer.compare(a.getIndex(),b.getIndex())); // sompostorder is respected
					List<Node> children = new ArrayList<>();
					for (Node N: uv_ancestors)
					{
						for (int ci=0; ci<N.getNumChildren(); ci++)
						{
							Node C = N.getChild(ci);
							if (!ancestors.contains(C))
								children.add(C);
						}
					}
					uv_children = children.toArray(empty);
//					System.out.println("#**P.SPR.gDS/U.iDS anc "+uv_ancestors.length
//							+"\t// "+Arrays.toString(uv_ancestors));
//					System.out.println("#**P.SPR.gDS/U.iDS cds "+uv_children.length
//							+"\t// "+Arrays.toString(uv_children));
				}
				
				/**
				 * Reevaluates the score only at the 2 nodes an their ancestors. 
				 */
				@Override 
				protected double scoreSankoff(int f) 
				{
					Profile og_profile = og_profiles[f];
					Profile spr_profile = this.new Profile(f, false);
					// hotstart from children's already computed scores
					for (Node C: uv_children)
					{
						assert (og_index.containsKey(C)); // since U an V are the prune and graft position parents 
						int child = C.getIndex();
						int og_child = og_index.get(C);
						spr_profile.subtree_copy_numbers[child] = og_profile.subtree_copy_numbers[og_child];
						spr_profile.sankoff_score[child] = og_profile.sankoff_score[og_child];
					}
					// the scores change only at the two nodes and their ancestors up to the root
					for (Node A: uv_ancestors)
					{
						spr_profile.scoreSankoff(A.getIndex());
					}
					// now the root is correctly computed 
					return spr_profile.scoreSankoffRoot();
				}
			}
			Update spr_update = new Update();
			double spr_score = spr_update.getSankoffScore();
			double diff_score = spr_score-factory_score;
//			System.out.println("#**P.SPR.gDS/U.sS score "+spr_score+"\tdiff "+diff_score
//					+"\twas "+factory_score
//					+"\t// U "+U.getNodeIdentifier()+"\tV "+V.getNodeIdentifier()
//					+"\tP "+P.getNodeIdentifier()+"\tG "+G.getNodeIdentifier());
			this.last_spr = spr;
			return diff_score;
		}
		
		/**
		 * The phylogeny constructed in the most recent call of {@link #getDiffScore(int, int)}. 
		 * 
		 * @return
		 */
		public Phylogeny getLastSPR()
		{
			return last_spr;
		}
	}

	
	
	/**
	 * For testing.
	 * 
	 * @param out
	 */
	protected void reportSankoff(PrintStream out)
    {
    	int nF = table.getFamilyCount();
    	out.println("# Row\tFamily\tNode\tMin\tMax\tProfile");
    	int num_nodes = tree.getNumNodes();
    	double score = 0.0;
    	for (int f=0; f<nF; f++)
    	{
    		Profile P = new Profile(f);
    		int[] min_copies = P.computeSankoff(false);
    		int[] max_copies = P.computeSankoff(true);
    		score += P.getSankoffScore();
    		for (int node=0; node<num_nodes; node++)
    		{
    			out.print(f+"\t"+table.getFamilyName(f)+"\t"+tree.getIdent(node));
    			out.print("\t"+min_copies[node]);
    			out.print("\t"+max_copies[node]);
    			if (tree.isLeaf(node))
    				out.print("\t"+P.profile[node]);
    			out.println();
    		}
    	}
    	double s = getSankoffScore();
    	out.println("#SCORE "+s+"\t("+score+")");
    }
    
    
    
    

    public static void main(String[] args) throws Exception
    {
    	CommandLine cli = new CommandLine(args, Profile.class, 2);
    	ProfileTable table = cli.getTable();
    	Phylogeny tree = cli.getTree();
    	
    	
    	PrintStream out = cli.getOutput(OPT_OUTPUT, System.out);

    	double gain_pty = cli.getOptionDouble(OPT_GAIN, DEFAULT_GAIN_PENALTY);
    	double loss_pty = cli.getOptionDouble(OPT_LOSS, DEFAULT_LOSS_PENALTY);
    	double duplication_pty = cli.getOptionDouble(OPT_DUPLICATION, DEFAULT_DUPLICATION_PENALTY);
    	
    	out.println(CommandLine.getStandardHeader("Parsimony penalties: -"+OPT_GAIN+" "+gain_pty)
    				+"\t-"+OPT_LOSS+" "+loss_pty+"\t-"+OPT_DUPLICATION+" "+duplication_pty
    				+"\t(contraction 1)");
//    	
//    	
//    	String arg_gain = cli.getOptionValue(OPT_GAIN);
//    	String arg_dup = cli.getOptionValue(OPT_DUPLICATION);
//    	String arg_loss = cli.getOptionValue(OPT_LOSS);
//    	
//    	if (arg_gain!=null)
//    		gain_pty = Double.parseDouble(arg_gain);
//    	if (arg_dup!=null)
//    		duplication_pty = Double.parseDouble(arg_dup);
//    	if (arg_loss!=null)
//    		loss_pty = Double.parseDouble(arg_loss);
    	
//    	UniqueProfileTable utable = new UniqueProfileTable(table);
//    	SPRExplorer spr = new SPRExplorer(tree, utable, gain_pty, loss_pty, duplication_pty);
//    	double sc = spr.getDiffScore(0, 6);
    	
    	Parsimony factory = new Parsimony(tree, table);
    	factory.setPenalties(gain_pty, loss_pty, duplication_pty);
    	
    	factory.reportSankoff(out);
    	
    	
    }
}
