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


import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.function.Function;


import count.ds.AnnotatedTable;
import count.ds.Heap;
import count.ds.IndexedTree;
import count.ds.Phylogeny;
import count.ds.ProfileTable;
import count.ds.Phylogeny.Node;
import count.ds.TreeComparator;
import count.ds.UPGMA;
import count.ds.UniqueProfileTable;
import count.io.CommandLine;
import count.io.DataFile;
import count.io.GeneralizedFileReader;
import count.io.NewickParser;
import count.io.RateVariationParser;
import count.io.CountXML;

import static count.io.CommandLine.OPT_BOOTSTRAP;
import static count.io.CommandLine.OPT_BUILD;
import static count.io.CommandLine.OPT_CONTRACT;
import static count.io.CommandLine.OPT_EPS;
import static count.io.CommandLine.OPT_LOAD;
import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_PLACE;
import static count.io.CommandLine.OPT_POLLARD;
import static count.io.CommandLine.OPT_REROOT;
import static count.io.CommandLine.OPT_RND;
import static count.io.CommandLine.OPT_ROUNDS;
import static count.io.CommandLine.OPT_SAVE;
import static count.io.CommandLine.OPT_SEARCH;
import static count.io.CommandLine.OPT_SPR;
import static count.io.CommandLine.OPT_TOUR;
import static count.io.CommandLine.OPT_TRUNCATE;
import static count.io.CommandLine.OPT_OUTPUT;
//import static count.io.CommandLine.OPT_OUTPUT_SAVEALL;
import static count.io.CommandLine.OPT_OUTPUT_TREE;
import static count.io.CommandLine.OPT_PARSIMONY_FIT;
import static count.io.CommandLine.OPT_SNAPSHOT;
import static count.io.CommandLine.OPT_TOP;

import count.io.ModelBundle;

import static count.model.GLDParameters.PARAMETER_GAIN;
import static count.model.GLDParameters.PARAMETER_LOSS;
import static count.model.GLDParameters.PARAMETER_DUPLICATION;
import static count.model.GLDParameters.PARAMETER_LENGTH;
import static count.model.ML.EPS;
import static count.model.TreeWithRates.DEFAULT_EDGE_LENGTH;



/**
 * Experimental code for phylogeny reconstruction from copy-number data.
 * The optimization proceeds by considering subtree-prune-and-regraft (SPR) 
 * operations as candidate <em>moves</em>, 
 * defined by pairs of nodes for prune and regraft positions. 
 * Candidate moves are first score heuristically, 
 * and then with maximum-likelihood.  
 *  
 * <p>The heuristic scoring is done using minimum parsimony
 * by default. SPR moves are stored in a min-heap by heuristic score 
 * with each prune node ({@link MLPhylogeny#node_moves}). Top-scoring candidates from each prune position
 * are scored with ML, and are stored in a min-heap ({@link MLPhylogeny#model_moves}). In particular, the ML-score is 
 * the difference in negative log-likelihood if the SPR move is made. 
 * The best top move is accepted if the score is below a given threshold.   
 * After an SPR move, the data structures are updated to include new SPR moves, and delete 
 * those that are not possible anymore. The candidate moves that stay valid are rescored lazily 
 * on the heaps, trying to avoid recalculation of all scores. 
 * 
 *  <ul>
 * <li>{@link #buildModel(double)} builds the phylogeny from scratch. It initializes a star phylogeny
 * with infinitely long edges leading to the leaves from a single root. 
 * Infinite edges encode a forest of finite-length phylogenies, each rooted at 
 * an infinite edge. A completely resolved phylogeny contains (at most) one infinite edge 
 * leading to the root (to define the root prior distribution for copy numbers). The starting model 
 * may have one or more infinite-length edges. </li>
 * 
 * <li>{@link #walkModel(double)} makes a series of SPR moves, selecting the next move greedily until   
 * no more moves are available with a score change below the threshold. </li>
 * 
 * <li>{@link #placeRoot(double)} scores possible alternative rerootings descending from the root. </li>
 * 
 * <li>{@link #contractEdges(double)} contracts edges that do not affect the model score much. </li>
 * 
 * <li>{@link #pollardEdges(double)} prunes high-loss edges, grafts them at the root, and recomputes
 * model rates.</li>
 * 
 * <li>{@link #main(String[])} gives a command-line interface to launching topology exploration 
 * with a standard model class (3 parameters per edge).  </li>
 * 
 * </ul>
 * 
 * <p>Possible topology transformations are stored as {@link CandidateMove} instances. </p>
 * 
 * <p><strong>Time complexity.</strong> ML-scoring of a model takes quadratic time 
 * (times a constant depending on convergence settings) in the tree size. Initially, and 
 * for every current phylogeny over <var>n</var> leaves, there are O(<var>n</var><sup>2</sup>)
 * possible SPR moves (identified as pairs of prune and regraft nodes). 
 * At every iteration of accepting a move, 
 * O(<var>n</var>) SPR moves are created/destroyed; O(<var>t n</var>) are 
 * scored or rescored heuristically, and 
 * O(<var>t</var>) are scored with full ML, 
 * where <var>t</var>={@link #top_scored_moves} which can be set via 
 * {@link #setTopScoredMoves(int)}. 
 * Heuristic scoring by parsimony
 * ({@link CandidateMove#scoreParsimony(count.model.Parsimony.SPRExplorer)})
 * is O(<var>n</var>), 
 * but by FastML ({@link CandidateMove#scoreFastML(int)}), it is O(<var>n</var><sup>2</sup>).
 * (The heuristic scoring is set with the flag {@link #heuristic_byParsimony}.)   
 * {@link #buildModel(double)} and {@link #walkModel(double)}   
 * take O(<var>n</var>)
 * iterations.
 * Hence the total time is O(<var>t n</var><sup>3</sup>) 
 * for heuristic scoring with parsimony, or 
 * O(<var>t n</var><sup>4</sup>) for heuristic scoring with FastML.
 * (The SPR update per iteration calls {@link #recomputeWalkMoves()} 
 * which takes O(<var>n</var><sup>2</sup>) bc it enumerates all SPR moves, 
 * but it scores heuristically only the O(<var>n</var>) newly created ones.)
 * {@link #placeRoot(double)}
 * is O(<var>t n</var><sup>2</sup>); {@link #pollardEdges(double)} is O(<var>n</var><sup>2</sup>). 
 * {@link #contractEdges(double)} scores O(<var>n</var>) 
 * contractions initially, and then at worst O(<var>n</var>) but typically O(1) 
 * are rescored in each of O(<var>n</var>) iterations. Thus O(<var>n</var><sup>3</sup>) 
 * for a sequence of contractions.  </p>
 * 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */
public class MLPhylogeny 
{
	
	public static boolean PRINT_OPTIMIZATION_MESSAGES = true;
	private static final boolean DEBUG_SCORE = true;
	
	/**
	 * Absolute threshold for truncated computation.
	 */
	private static final int DEFAULT_CALCULATION_WIDTH_ABSOLUTE = 12;
	/**
	 * Relative threshold for truncated computation
	 */
	private static final double DEFAULT_CALCULATION_WIDTH_RELATIVE = 3.0;
	/**
	 * Max iterations for numerical optimization of distribution parameters 
	 */
	private static final int DEFAULT_ROUNDS = 200;
	/**
	 * Convergence (provides bound for gradient)
	 */
	private static final double DEFAULT_CONVERGENCE = 1e-6;
	/**
	 * Number of ML-scored models per node. 
	 */
	private static final int TOP_MOVES = 12;
	
	public MLPhylogeny(AnnotatedTable table)
	{
		this.table = table;
		this.min_copies = Integer.min(2,table.minCopies());
		this.utable = new UniqueProfileTable(table);
		setRandom(null);
	}
	
	
	private final AnnotatedTable table;
	private final UniqueProfileTable utable;
	
	
	private int calculation_width_absolute = DEFAULT_CALCULATION_WIDTH_ABSOLUTE;
	private double calculation_width_relative = DEFAULT_CALCULATION_WIDTH_RELATIVE;
	
	private int opt_rounds = DEFAULT_ROUNDS;
	private double opt_eps = DEFAULT_CONVERGENCE;
	
	private int min_copies;
	
	private final boolean eager_node_moves = false;
	private final boolean explore_nni = false; // too much memory
	private final boolean exclude_multiple_prune = true;
	
	private Random RND; 
	final 
	void setRandom(Random RND) { this.RND = RND==null?new Random(2023):RND;}
	
	private String session_id=null;
	
	public void SetSession(String session_id)
	{
		this.session_id = session_id;
	}
	
	public void setCalculationWidth(int absolute, double relative)
	{
		this.calculation_width_absolute = absolute;
		this.calculation_width_relative = relative;
	}
	
	public void setOptimizationConvergence(int max_rounds, double eps)
	{
		this.opt_eps= eps;
		this.opt_rounds = max_rounds;
	}
	
	public void setMinCopies(int min_copies)
	{
		if (min_copies > table.minCopies())
			throw new IllegalArgumentException();
 		this.min_copies = Integer.min(2,min_copies);
	}
	
	public void setOptimizationHistory(List<Double> history)
	{
		this.history = history;
	}
	
//	public void setRates(TreeWithRates rates)
//	{
//		best_model = new CandidateMove(rates);
//		best_model.score=best_model.scoreML(0);
//		System.out.println("#**MLP.sRM set "+best_model);
//	}
	
	public TreeWithRates getRates()
	{
		return best_model.rates;
	}
	
	/**
	 * Convenience method to wrap {@link #getRates()} into 
	 * 1 class. 
	 * 
	 * @return the rate model as GammaInvariant
	 */
	public GammaInvariant getRateVariation()
	{
		GammaInvariant model = new GammaInvariant(getRates(), 1, 1, 1, 1);
		return model;
	}
	
	public Phylogeny getPhylogeny()
	{
		return best_model.forest;
	}
	
	public double getScore()
	{
		return best_model.score;
	}
	
	public void randomModel()
	{
//		System.out.println("#**MLP.rM start");
		Phylogeny tree = Phylogeny.randomTree(table.getTaxonNames(), RND, true);
		TreeWithRates rates = new TreeWithRates(tree,RND);
		System.out.println("#**MLP.rM starting tree "+NewickParser.printTree(tree));
		
		this.best_model = new CandidateMove(rates);
		best_model.score = best_model.score();

		if (history != null)
			history.add(best_model.score);
		
		System.out.println("#**MLP.rM init "+best_model+"\t"+NewickParser.printTree(best_model.forest));
	}
	
	/**
	 * Builds a starting tree using UPGMA, and optimizes the model parameters.
	 */
	public void upgmaModel()
	{
		Phylogeny tree = UPGMA.buildTree(utable, UPGMA.Dissimilarity.BRAY_CURTIS, UPGMA.ClusteringPolicy.UPGMA);
		TreeWithRates rates = new TreeWithRates(tree,RND);
		System.out.println("#**MLP.upgma starting tree "+NewickParser.printTree(tree));
//		System.out.println(RateVariationParser.printRates(rates));
		
		this.best_model = new CandidateMove(rates);

		best_model.score = best_model.score();

		if (history != null)
			history.add(best_model.score);
		if (near_best_models == null)
		{
			near_best_models = new NearBestModels();
		}
		near_best_models.addMove(best_model);
		System.out.println("#**MLP.upgma init "+best_model+"\t"+NewickParser.printTree(best_model.forest));		
		
	}
	
	
	/**
	 * Building a model: 
	 * <ol>
	 * <li> {@link #initModel(TreeWithRates)} sets up a forest for {@link #best_model} </li>
	 * <li> {@link #initBuildMoves()} adds all starting SPR moves</li>
	 * <li> {@link #refineModel(double)} iteration step for accepting one SPR at a time</li>
	 * </ol>
	 * 
	 * 
	 * @param refine_threshold maximum score change at a refinement step
	 * 
	 */
	public void buildModel(double refine_threshold)
	{
		initModel(null);
		
		if (near_best_models == null)
		{
			near_best_models = new NearBestModels(refine_threshold);
			near_best_models.addMove(best_model);
		} else
			near_best_models.setThreshold(refine_threshold);

		
		if (history != null)
			history.add(best_model.score);
		
		initBuildMoves();
		
		double delta = 0.0;
		int num_leaves = table.getTaxonCount();
		
		for (int i=0; i<num_leaves-1 && !model_moves.isEmpty(); i++)
		{
			
			
			double d = this.refineModel(refine_threshold);
			if (d<refine_threshold)
			{
				delta += d;
				System.out.println("#CANDIDATE\t"+i+"\t"+best_model+"\td "+d+"\tdelta "+delta+"\t"+NewickParser.printTree(best_model.forest));
				if (history != null)
					history.add(best_model.score);
				
				near_best_models.addMove(best_model);
			} else
			{
				System.out.println("#**MLP.bM no more steps under the threshold: d "+d);
				if (explore_nni) exploreNNI(Math.sqrt(best_model.varML()));
				break;
			}
		}
		
//		return best_model.model;
	}
	
	
	
	
//	public void parsimonyBuild(double refine_threshold)
//	{
//		node_moves = new HashMap<>();
//		spr_parsimony = new Parsimony.SPRExplorer(best_model.forest, utable, DEFAULT_PARSIMONY_GAIN, DEFAULT_PARSIMONY_LOSS, DEFAULT_PARSIMONY_DUPLICATION);	
//		move_heuristics = new CandidateHeuristics(m->m.scoreParsimony(getSPRExplorer()));
//		model_moves = new Heap<>();
//
//		Node[] trees = best_model.getTreeRoots();
//		for (Node T: trees)
//		{
//			addForestMoves(T);
//		}
//		
//		// choose the best
//		{
//			Set<CandidateMove> node_top_moves = new HashSet<>();
//			for (Node N: node_moves.keySet())
//			{
//				Heap<CandidateMove> Nmoves = node_moves.get(N);
//				if (!Nmoves.isEmpty())
//				{
//					CandidateMove top = Nmoves.peek();
//					if (node_top_moves.contains(top))
//					{
//						// mutually best : score it 
//						
//						
//					} else
//						node_top_moves.add(top);
//				}
//			}
//		}
//		
//	}
	
	public void tourModel(int num_walks, double walk_threshold, double top_increase, double convergence_decrease)
	{
		for (int w=num_walks; 0<w;)
		{
			--w;
			walkModel(walk_threshold);
			System.out.println("#TOUR "+w+"\tscore "+best_model.score
					+"\ttop "+top_scored_moves
					+"\tconv "+opt_eps
					+"\t"+NewickParser.printTree(best_model.forest));
			if (w!=0)
			{
				int t = (int)(top_scored_moves*top_increase);
				setTopScoredMoves(t);
				double eps = opt_eps * convergence_decrease;
				setOptimizationConvergence(opt_rounds, eps);
			}
		}
	}
	
	/**
	 * Stochastic search in model space by perturbations (jumps)
	 * via random NNI sequences, followed by model walking.  
	 * 
	 * @param max_unsuccessful
	 * @param search_threshold
	 */
	public void searchModels(int max_unsuccessful, double search_threshold)
	{
		if (near_best_models == null)
		{
			near_best_models = new NearBestModels(search_threshold);
			near_best_models.addMove(best_model.snapshot());
		} else // already there
			near_best_models.setThreshold(search_threshold);
		
		int jump_distance = (int)((best_model.forest.getNumLeaves()-2.0)*0.5); // IQTree's logic: half of inner edges
				
		int nochange = 0;
		int iter = 0;
		
		while (nochange<max_unsuccessful)
		{
			CandidateMove starting_model = best_model;
			CandidateMove[] pick_one = near_best_models.allMoves();
			double current_best_score = near_best_models.bestest_score;
			
			int random_i = RND.nextInt(Integer.min(top_scored_moves,pick_one.length));
			
			CandidateMove C = pick_one[random_i];
			this.best_model = C.snapshot(); // need to create copy that we can mess with
			
			if (pick_one.length!=1)
			{
				stochasticNNI(jump_distance);
				CandidateMove J = this.best_model;
	//			boolean goodJ = near_best_models.addMove(J);
				System.out.println("#**MLP.sM "+iter+"\t"+nochange+"\tpick "+C+"\t("+random_i+"/"+pick_one.length+"\tjump "+J+"\tfrom "+starting_model+"\tstart "+current_best_score);
				
		    	if (DEBUG_SCORE) debugScore(best_model, System.out, "MLP.sM/jump");
			}
	    	
			walkModel(search_threshold);
			debugScore(best_model, System.out, "MLP.sM/walk");
			
			if (best_model.score < current_best_score)
			{
				// this is really better?
				double nFam = utable.getTotalFamilyCount();
				double old_avg = current_best_score / nFam;
				double avg_score = best_model.score / nFam;
				
				double rel_chg = (old_avg-avg_score)/Math.max(old_avg,1.0);

				System.out.println("#**MLP.sM "+iter+"\t"+nochange
						+"\tstart "+current_best_score
						+"\tjump "+C+"\twalk "+best_model
						+"\trchg "+rel_chg+"\tclose "+(rel_chg<opt_eps)
						+"\tnbm "+near_best_models.bestest_model+"\t/"+near_best_models.bestest_score);
				
				if (rel_chg<opt_eps)
					nochange++;
				else if (fit_parsimony)
				{
					fitParsimonyPenalties(best_model.rates, utable.mappedToTree(best_model.forest));
				}

				current_best_score = best_model.score;
				// keep the best model
//				if (near_best_models.bestest_score != best_model.score)
//				{
//					current_best_score = near_best_models.bestest_score;
//				}
//				this.best_model = near_best_models.bestest_model;
				//assert (near_best_models.bestest_score == best_model.score);
			} else
			{
				System.out.println("#**MLP.sM "+iter+"\t"+nochange
						+"\tstart "+current_best_score
						+"\tjump "+C+"\twalk "+best_model
						+"\twalk "+best_model
						
						+"\tback "+starting_model
						+"\tnear "+near_best_models.contains(best_model));
				nochange++;
				this.best_model = starting_model;
				// not better 
			}
			
			if (this.session_id != null && nochange<max_unsuccessful)
			{
				// save snapshot of state
				String bundle_file = session_id+"_s"+iter+"_"+nochange+".xml";
				try
				{
					saveNearBestModels(bundle_file, session_id);
					System.out.println("#**MLP.sM snapshot in "+bundle_file+" (session "+session_id+")");
				} catch (Exception E)
				{
					System.out.println("#**MLP.sM snapshot failed "+E+" for "+bundle_file+" (session "+session_id+")");
				}
			}
			iter++;
		}			// end of iterations

	}
	
	/**
	 * Improves the model by SPR moves. 
	 * Call when {@link #best_model} is set, after {@link #initModel(TreeWithRates)}
	 * or {@link #buildModel(double)}.  
	 * 
	 * Improving a model: 
	 * <ol>
	 * <li> {@link #initWalkMoves()} adds all starting SPR moves </li>
	 * <li> {@link #improveModel(double)} iteration step for accepting one SPR move at a time</li>
	 * </ol>
	 * 
	 * @param walk_threshold max accepted score change (exclusive bound)
	 */
	public void walkModel(double walk_threshold)
	{
		CandidateMove best_ever = best_model.snapshot();
		
		if (near_best_models == null)
		{
			near_best_models = new NearBestModels(walk_threshold);
			near_best_models.addMove(best_ever);
		} else // already there
			near_best_models.setThreshold(walk_threshold);
		
		if (history != null)
			history.add(best_model.score);

		
//		if (optimize_EM)
//		{
//			String prefix = "#walkModel.TRANSITIONS\tM"+best_model.id+"."+best_model.version;
//			System.out.println(prefix+"\tscore "+best_model.score);
//			DirectEM M = new DirectEM(best_model.rates, utable.mappedToTree(best_model.forest), true);
//			setMLParameters(M);			
//			M.reportTransitionCounts(System.out, prefix);
//		}
		
		initWalkMoves();
		
		
		

//		exploreNNI(Math.sqrt(best_model.varML()));

		double delta = 0.0;
		int num_nodes = best_model.forest.getNumNodes();
		
		for (int i=0; i<num_nodes-1 && !model_moves.isEmpty(); i++)
		{
			double d = this.improveModel(walk_threshold);
			
			if (d<walk_threshold)
			{
				delta += d;
				System.out.println("#CANDIDATE\t"+i+"\t"+best_model+"\td "+d+"\tdelta "+delta+"\t"+NewickParser.printTree(best_model.forest));
				if (history!=null)
					history.add(best_model.score);
				if (best_model.score<best_ever.score)
				{
					best_ever = best_model.snapshot();
					System.out.println("#**MPL.wM bestever "+best_ever);
					near_best_models.addMove(best_ever);
				} else
				{
					near_best_models.addMove(best_model.snapshot()); // store a copy bc best_model will change
				}
			} else
			{
				System.out.println("#**MLP.wM no more steps under the threshold: d "+d);
				if (explore_nni) exploreNNI(Math.sqrt(best_model.varML()));
				break;
			}
		}
		
		if (best_model.version != best_ever.version)
		{
			System.out.println("#**MLP.wM current "+best_model+"\treverting to "+best_ever);
			best_model = best_ever; // was better
		}
		
		
		// optimize one last time --- returns after 1 iteration in principle
		best_model.score = best_model.score();
		// update associated score 
		near_best_models.update(best_model);
		
		
//		near_best_models.addMove(best_model.snapshot()); // updates score 
		
//		return best_model.model;
		
		
		
		
//		if (optimize_EM)
//		{
//			String prefix = "#TRANSITIONS\tM"+best_model.id+"."+best_model.version;
//			System.out.println(prefix+"\tscore "+best_model.score);
//			DirectEM M = (DirectEM)parameterOptimizer(best_model.rates, utable.mappedToTree(best_model.forest));			
//			M.reportTransitionCounts(System.out, prefix);
//		}
		
		
	}
	
	
	public void placeNodes(double walk_threshold, Set<String> wanted_leaf_names)
	{
		CandidateMove best_ever = best_model.snapshot();
		
		if (near_best_models == null)
		{
			near_best_models = new NearBestModels(walk_threshold);
			near_best_models.addMove(best_ever);
		} else // already there
			near_best_models.setThreshold(walk_threshold);
		
		if (history != null)
			history.add(best_model.score);
		
		initWalkMoves(wanted_leaf_names);
		
		
		if (PRINT_OPTIMIZATION_MESSAGES)
		{
			List<CandidateMove> all_moves = new ArrayList<>();
			for (CandidateMove M: model_moves)
				all_moves.add(M);
			Collections.sort(all_moves);
			for (int mi=0; mi<all_moves.size(); mi++)
			{
				CandidateMove M = all_moves.get(mi);
				System.out.printf("#PLACE\t%d\t%d\t%d\t%g\t%g\t%s\n"
						, mi
						, M.og_left.getIndex()
						, M.og_right.getIndex()
						, M.score
						, move_heuristics.get(M)
						, M.toString()
						);
				
			}
		}
		
		double delta = 0.0;
		for (int i=0; i<wanted_leaf_names.size() && !model_moves.isEmpty(); i++)
		{
			double d = this.improveModel(walk_threshold);
			
			if (d<walk_threshold)
			{
				delta += d;
				System.out.println("#CANDIDATE\t"+i+"\t"+best_model+"\td "+d+"\tdelta "+delta+"\t"+NewickParser.printTree(best_model.forest));
				if (history!=null)
					history.add(best_model.score);
				if (best_model.score<best_ever.score)
				{
					best_ever = best_model.snapshot();
					System.out.println("#**MPL.wM bestever "+best_ever);
					near_best_models.addMove(best_ever);
				} else
				{
					near_best_models.addMove(best_model.snapshot()); // store a copy bc best_model will change
				}
			} else
			{
				System.out.println("#**MLP.wM no more steps under the threshold: d "+d);
				break;
			}
		}
		
		if (best_model.version != best_ever.version)
		{
			System.out.println("#**MLP.wM current "+best_model+"\treverting to "+best_ever);
			best_model = best_ever; // was better
		}
		
		
		// optimize one last time --- returns after 1 iteration in principle
		best_model.score = best_model.score();
		// update associated score 
		near_best_models.update(best_model);
	}
	
	/**
	 * Tries all rootings and picks the one with maximum likelihood.
	 * Root placements are considered in the order of distance from 
	 * original root, trying up to {@link #top_scored_moves} choices. 
	 * The best root placement is accepted if its score is below the 
	 * refinement threshold.  
	 * 
	 * @param refine_threshold maximum score change (exclusive bound)
	 */
	public void placeRoot(double refine_threshold)
	{
		assert (best_model!=null);
		assert (best_model.getNumTrees()==0);
		
		//best_model.score = best_model.scoreML();
		
		Phylogeny phylo = best_model.forest;
		int num_nodes = phylo.getNumNodes();
		for (int node=0; node < num_nodes; node++)
		{
			if (!phylo.isLeaf(node))
				phylo.getNode(node).setName(null); // since we will rotate subtrees around 
		}
		
		// gather the closest edges to the root
		Map<Node, Double> node_distances = new HashMap<>();
		Heap<Node> explored_nodes = new Heap<>((a,b)->Double.compare(node_distances.get(a), node_distances.get(b)));
		Heap<CandidateMove> root_moves = new Heap<>(); 
		
		Node root = phylo.getRootNode();
		node_distances.put(root, 0.0);
		explored_nodes.add(root);
		
		while (!explored_nodes.isEmpty()  && root_moves.size()<top_scored_moves)
		{
			Node N = explored_nodes.deleteLeast();
			int node = N.getIndex();

			Node P = N.getParent();
			if (P!=null && (!P.isRoot() || P.getNumChildren()>2))
			{
				// possible root
				CandidateMove alt = best_model.reroot(node);
				alt.score = alt.score()-best_model.score;
				if (alt.score>=refine_threshold)
				{
					// don't add
					System.out.println("#**MLP.pR noadd "+alt+"\t"+NewickParser.printTree(alt.forest));
				} else
				{
					System.out.println("#**MLP.pR add "+alt+"\t"+NewickParser.printTree(alt.forest));
					root_moves.add(alt);
				}
			}
			// determine children's distances
			double d = node_distances.get(N);
			for (int ci=0; ci<N.getNumChildren(); ci++)
			{
				Node C = N.getChild(ci);
				int child = C.getIndex();
				double len = best_model.rates.getEdgeLength(child);
				node_distances.put(C, d+len);
				explored_nodes.add(C);
			}
		}
		
//		// examine all possible root placements
//		for (int node=0; node<phylo.getNumNodes(); node++)
//		{
//			if (!phylo.isRoot(node))
//			{
//				int parent = phylo.getParent(node);
//				if (!phylo.isRoot(parent) || phylo.getNumChildren(parent)>2) // if binary root, exclude its chiold edges bc they yield the same phylo after removing unary nodes
//				{
//					CandidateMove alt = best_model.reroot(node);
//					alt.score = alt.score()-best_model.score;
//					if (alt.score>=refine_threshold)
//					{
//						// don't add
//						System.out.println("#**MLP.pR noadd "+alt+"\t"+NewickParser.printTree(alt.forest));
//					} else
//					{
//						System.out.println("#**MLP.pR add "+alt+"\t"+NewickParser.printTree(alt.forest));
//						root_moves.add(alt);
//					}
//				}
//			}
//		}
		
		
		
		
		if (!root_moves.isEmpty())
		{
			CandidateMove alt = root_moves.deleteLeast();
			if (alt.score<refine_threshold)
			{
				System.out.println("#**MLP.pR root "+alt+"\t"+NewickParser.printTree(alt.forest));
				best_model = alt;
				best_model.score = best_model.scoreML(0);
				best_model.og_edge = null;
				best_model.og_model = null;
			}
			else
			{
				System.out.println("#**MLP.pR stays "+best_model+"\t"+NewickParser.printTree(best_model.forest));
			}
		} else
		{
			System.out.println("#**MLP.pR noalt "+best_model+"\t"+NewickParser.printTree(best_model.forest));
		}
		
		while (!root_moves.isEmpty())
		{
			CandidateMove alt = root_moves.deleteLeast();
			if (alt.score>=refine_threshold)
			{
				break;
			}
			System.out.println("#**MLP.pR alt.root "+alt+"\t"+NewickParser.printTree(alt.forest));
		}
	}
	
	/**
	 * Contracts edges in a greedy order for score change; 
	 * each step must change by at most the given contract threshold. 
	 * For instance, choosing 3/2*log(F) for F families 
	 * optimizes scoring by Bayesian Information Criterion
	 * with 3 parameters per edge. Note that with positive threshold, the 
	 * model score (negative log-likelihood) may increase, so follow with 
	 * {@link #walkModel(double)}. 
	 *
	 * UNTESTED 
	 * 
	 * @param contract_threshold maximum score change allowed
	 * 
	 * 
	 */
	public void contractEdges(double contract_threshold)
	{
		assert (best_model!=null);
		assert (best_model.getNumTrees()==0);
		
		if (near_best_models == null)
		{
			near_best_models = new NearBestModels(contract_threshold);
			near_best_models.addMove(best_model.snapshot());
		} else // already there
			near_best_models.increaseThreshold(contract_threshold);
		
		
		Heap<CandidateMove> contraction_moves = new Heap<>(); 
		Phylogeny phylo = best_model.forest;
		for (int node=0; node<phylo.getNumNodes(); node++)
		{
			if (!phylo.isLeaf(node) && !phylo.isRoot(node))
			{
				CandidateMove move = best_model.contract(node);
				move.score = move.score() - best_model.score;
				contraction_moves.add(move);
				System.out.println("#**MLP.cE add/"+contraction_moves.size()+"\t"+move);
			}
		}
		
		int i = 0;
		double delta = 0.0;
		while (!contraction_moves.isEmpty())
		{
			
			CandidateMove best = contraction_moves.peek();
			while (best.isStale())
			{
				double bscore = best.score;
				best.score = best.scoreFastML() - best_model.score;
				double bdiff = best.score-bscore;
				System.out.println("#**MLP.cE refresh "+i+"\t"+best+"\twas "+bscore+"\tdiff "+bdiff);
				contraction_moves.updateOrder(best);
				best = contraction_moves.peek();
			}
			best = contraction_moves.deleteLeast();
			if (best.score > contract_threshold)
			{
				System.out.println("#**MLP.cE no more moves: "+best);
				break;
			} else
			{
				System.out.println("#**MLP.cE best "+i+"\t"+best);
			}
			double score = best_model.score;
			best.acceptContract();
			best_model.score = best_model.score();
			
			near_best_models.addMove(best_model.snapshot());
			
			double d = best_model.score - score;
			delta += d;
			System.out.println("#CANDIDATE\t"+i+"\t"+best_model+"\td "+d+"/"+best.score+"\tdelta "+delta+"\t"+NewickParser.printTree(best_model.forest));
			i++;
		}	
	}
	
	/**
	 * Prunes high-loss edges and regrafts them at the root.
	 * 
	 * @param loss_prob_threshold (exclusive) lower threshold on loss probability for pruning edges
	 */
	public void pollardEdges(double loss_prob_threshold)
	{
		Set<Node> long_edges = new HashSet<>();
		for (int node=best_model.forest.getNumNodes(); 0<node; )
		{
			--node;
			double p = best_model.rates.getLossParameter(node);
			Node N = best_model.forest.getNode(node);
			if (!N.isRoot() && !N.getParent().isRoot() 
					&& loss_prob_threshold < p)
			{
				long_edges.add(N);
			}
		}
		if (!long_edges.isEmpty())
		{
			Node R = best_model.forest.getRootNode();
			Map<Node, double[]> tree_rates = modelRates(best_model.rates);
			for (Node N: long_edges)
			{
				if (!N.isRoot() && !N.getParent().isRoot())
				{
					System.out.println("#**MLP.polE "+best_model.forest.getIdent(N.getIndex())+"\t// "+N);
					best_model.forest.pruneAndRegraft(N, R, true);
				}
			}
			TreeWithRates rates = new TreeWithRates(best_model.forest);
			for (int node=best_model.forest.getNumNodes(); 0<node; )
			{
				--node;
				Node N = best_model.forest.getNode(node);
				if (tree_rates.containsKey(N))
				{
					double[] node_params = tree_rates.get(N);
					rates.setGainRate(node, node_params[PARAMETER_GAIN]);
					rates.setLossRate(node, node_params[PARAMETER_LOSS]);
					rates.setDuplicationRate(node, node_params[PARAMETER_DUPLICATION]);
					double len = node_params[PARAMETER_LENGTH];
					rates.setEdgeLength(node, len);
					N.setLength(Double.max(EPS,len));
				}
			}
			best_model.rates = rates;
			double old_score = best_model.score;
			best_model.score = best_model.score();
			double diff = best_model.score - old_score;
			System.out.println("#**MLP.polE "+best_model+"\tscore "+diff);
		}
	}
	
	
	private CandidateMove best_model;
	private Heap<CandidateMove> model_moves;
	private NearBestModels near_best_models=null;
	
	
	/**
	 * Candidate moves for each prune position; min-heaps ordered by heuristic score. 
	 * 
	 */
	private Map<Node,Heap<CandidateMove>> node_moves;	
	/**
	 * Candidate moves with their heuristic score.   
	 */
	private CandidateHeuristics move_heuristics;
	private List<Double> history=null;

	/*
	 * model walking
	 */
	private Set<Node> accepted_prune_moves=null;
	
	private Parsimony.SPRExplorer spr_parsimony=null;
	
	private Parsimony.SPRExplorer getSPRExplorer(){return spr_parsimony;}
	
	
	/**
	 * Initializes the search or walk. 
	 * 
	 * @param rates if null, constructs the starting model with a star phylogeny (=forest of single-leaf trees)
	 */
	public void initModel(TreeWithRates rates)
	{
		if (rates == null)
		{
			Phylogeny base_phylo = Phylogeny.starTree(table.getTaxonNames());
			base_phylo.hasLength(true);
			rates = new TreeWithRates(base_phylo);
			for (int leaf=0; leaf<base_phylo.getNumLeaves(); leaf++)
			{
				rates.setEdgeLength(leaf, Double.POSITIVE_INFINITY);
				base_phylo.getNode(leaf).setLength(Double.POSITIVE_INFINITY);
			}
		} 
		
		this.best_model = new CandidateMove(rates);
		best_model.score = best_model.score();
		
		System.out.println("#**MLP.iM init "+best_model+"\t"+NewickParser.printTree(best_model.forest));
	}
	
	
	
	/**
	 * Initialization and population of the data structures
	 * for model building. SPR moves are between leaf pairs.
	 * 
	 * <ol>
	 * <li> {@link #addForestMoves(Node)} iteration step for adding the SPR moves with a leaf </li>
	 * <li> {@link #updateModelMoves()} called once at the end to advance the best node moves to {@link #model_moves}</li>
	 * </ol>
	 */
	private void initBuildMoves()
	{
//		model_moves=pairedTrees(best_model);
		node_moves = new HashMap<>();
//		double hscore = best_model.scoreParsimony();
		if (heuristic_byParsimony) // very first call; utable is mapped to best_model
			spr_parsimony = new Parsimony.SPRExplorer(best_model.forest, utable, DEFAULT_PARSIMONY_GAIN, DEFAULT_PARSIMONY_LOSS, DEFAULT_PARSIMONY_DUPLICATION);	
//		double hscore = spr_parsimony.getSankoffScore(); 

		move_heuristics = new CandidateHeuristics(
				heuristic_byParsimony?m->m.scoreParsimony(getSPRExplorer()):new FastMLHeuristic());

		//move_heuristics.put(best_model, hscore);
		model_moves = new Heap<>();

		Node[] trees = best_model.getTreeRoots();
		for (Node T: trees)
		{
			addForestMoves(T);
		}
		
//		
//		{ // DEBUG
//			for (Node T: trees)
//			{
//				Heap<CandidateMove> Tmoves = node_moves.get(T);
//				CandidateMove M = Tmoves.peek();
//				
//				System.out.println("#**MLP.iBM node "+T.getIndex()+"\tnmoves "+Tmoves.size()+"\ttop "+M+"\t// "+T);
//			}
//			
//			if (2==1+1) throw new RuntimeException("*** DEBUG stop ***");
//		}
		updateModelMoves();
	}

	/**
	 * Initialization and population of the data structures 
	 * for model walking.
	 * <ol>
	 * <li> {@link #addWalkMoves(Node)} iteration step for adding the SPR moves with a node</li>
	 * <li> {@link #updateModelMoves()} called once at the end to advance the best node moves to {@link #model_moves}</li>
	 * </ol>
	 */
	private void initWalkMoves()
	{
		node_moves = new HashMap<>();
		accepted_prune_moves = new HashSet<>();
		if (heuristic_byParsimony) // 
		{
//			double pty_gain = DEFAULT_PARSIMONY_GAIN;
//			double pty_loss = DEFAULT_PARSIMONY_LOSS;
//			double pty_dup = DEFAULT_PARSIMONY_DUPLICATION;
			UniqueProfileTable mapped_utable = utable.mappedToTree(best_model.forest);
//			{
//				DirectEM M = new DirectEM(best_model.rates, mapped_utable, true);
//				setMLParameters(M);
//				double[] pty = M.fitParsimonyPenalty();
//				pty_gain = pty[PARAMETER_GAIN];
//				pty_loss = pty[PARAMETER_LOSS];
//				pty_dup = pty[PARAMETER_DUPLICATION];
//			}
			
			spr_parsimony = new Parsimony.SPRExplorer(best_model.forest, mapped_utable, pty_gain, pty_loss, pty_dup);	
		}
		
//		double hscore = best_model.scoreParsimony();
//		double hscore = spr_parsimony.getSankoffScore();
//		move_heuristics.put(best_model, hscore);
		move_heuristics = new CandidateHeuristics(
				heuristic_byParsimony?m->m.scoreParsimony(getSPRExplorer()):new FastMLHeuristic());
		
		
		//best_model.score = best_model.score(); // already scored by initModel 
		model_moves = new Heap<>();

		for (Node N: best_model.forest.getRootNode().listNodes(null, null)) // Predicate.not(Node::isRoot)))
		{
			addWalkMoves(N);
//			System.out.println("#**MLP.iWM "+N);
			//updateModelMoves(0.0);
		}
		updateModelMoves(); 
	}
	
	private void initWalkMoves(Set<String> wanted_leaf_names)
	{
		node_moves = new HashMap<>();
		accepted_prune_moves = new HashSet<>();
		if (heuristic_byParsimony) // 
		{
			UniqueProfileTable mapped_utable = utable.mappedToTree(best_model.forest);
			spr_parsimony = new Parsimony.SPRExplorer(best_model.forest, mapped_utable, pty_gain, pty_loss, pty_dup);	
		}
		move_heuristics = new CandidateHeuristics(
				heuristic_byParsimony?m->m.scoreParsimony(getSPRExplorer()):new FastMLHeuristic());
		model_moves = new Heap<>();

		for (Node N: best_model.forest.getRootNode().listNodes(null, null)) 
		{
			if (!N.isLeaf() || !wanted_leaf_names.contains(N.getName()))
			{
				node_moves.put(N,  move_heuristics.newHeap()); // empty heaps
			}
		}
		for (Node N: best_model.forest.getRootNode().listNodes(null, null))
		{
			if (N.isLeaf() && wanted_leaf_names.contains(N.getName()))
			addWalkMoves(N);
		}
		updateModelMoves(); 
	}
	
	// 
	private int top_scored_moves = TOP_MOVES;
	/**
	 * Sets the number of top (by heuristic scoring)
	 * models per node moves that are considered for {@link #model_moves}.
	 * 
	 * @param ml_scored_models_per_node number of moves scored heuristically at a node
	 */
	public void setTopScoredMoves(int ml_scored_models_per_node)
	{
		this.top_scored_moves = ml_scored_models_per_node;
	}
	
	/**
	 * Accepts a move and updates {@link #best_model} and 
	 * {@link #spr_parsimony} accordingly. 
	 * 
	 * @param alt_model
	 * @return the newly created node: joint parent of the regraft an prune nodes
	 */
	private Node acceptJoinMove(CandidateMove alt_model)
	{
		assert (alt_model.og_model == best_model);
		alt_model.acceptJoin();
		best_model = alt_model.og_model;
		Phylogeny joined_phylo = best_model.forest;
		Node P = alt_model.og_left.getParent();
//		
//		Node L = alt_model.og_left;
//		Node R = alt_model.og_right;
//		
//		if (alt_model.isStale())
//		{
//			System.out.println("#**MLP.aJM version diff "+alt_model+"\there "+best_model+"\talt.og "+alt_model.og_model);
//		}
//
////		Map<Node,double[]> model_rates = modelRates(best_model.model.getBaseModel());
//
////		// join the two nodes 
//		best_model.forest.pruneAndRegraft(R, L, false);
//		Phylogeny joined_phylo = best_model.forest;
//		Node P = L.getParent();
//		assert (R.getParent() == P);
//
//		// adjust the model for the phylogeny change
//		TreeWithRates base = new TreeWithRates(joined_phylo);
//		TreeWithRates alt_base = alt_model.rates;
//		
//		for (int node=0; node<joined_phylo.getNumNodes(); node++)
//		{
//			base.setGainRate(node, alt_base.getGainRate(node));
//			base.setLossRate(node, alt_base.getLossRate(node));
//			base.setDuplicationRate(node, alt_base.getDuplicationRate(node));
//			base.setEdgeLength(node, alt_base.getEdgeLength(node));
//			joined_phylo.getNode(node).setLength(alt_model.forest.getLength(node));
//		}
//						
//		best_model.rates = base;
////		best_model.og_model = best_model;
////		best_model.og_left = L;
////		best_model.og_right = R;
//		best_model.version++;
		
		if (heuristic_byParsimony)
		{
			double pty_gain = spr_parsimony == null?DEFAULT_PARSIMONY_GAIN:spr_parsimony.getGainPenalty();
			double pty_loss = spr_parsimony == null?DEFAULT_PARSIMONY_LOSS:spr_parsimony.getLossPenalty();
			double pty_dup = spr_parsimony == null?DEFAULT_PARSIMONY_DUPLICATION:spr_parsimony.getDuplicationPenalty();
			
			spr_parsimony = new Parsimony.SPRExplorer(joined_phylo, utable.mappedToTree(joined_phylo), pty_gain, pty_loss, pty_dup);
		}
		
		return P;
	}
	
	/**
	 * One iteration for model refinement with forest:
	 *    make the best model move and update the 
	 *    set of possible node moves after the SPR. 
	 * 
	 * @param max_score max accepted move score 
	 * @return 0.0 if no good move, &ge; max_score if only bad models, &lt; max_score if successful
	 */
	private double refineModel(double max_score)
	{
		// I. accept the move suggested by the lowest score in model_moves:
		//     identifies left node L and right node R to be joined, with parent P; P is root when the forest is fully resolved 
		//     (a) L is either a tree root in the forest: length is +Infty, and parent is forest root
		//     (b) or L is a a nonroot node: finite length [only if not join roots]
		//		R is surely a tree root in the forest: length is +Infty, and parent is forest root
		// II. remove obviated moves from model moves (*) and node moves 
		//	II.1 update node moves involving L
		//		(a) found on node moves for L: if join roots only, then remove from pair's node moves+*;
		//										else put pair on the right (since it may be tree root), and * if pair is not tree root  
		//		(b) remove looping joins : L's subtree + a node in R; found on node moves for L's subtree root
		//	II.2 update node moves involving R as in II.1 (a)
		// III. add new moves for P if P is not root
		//		(a) P is either a tree root (when L is)
		//			generate node moves for P: pair with other tree roots (if join roots), 
		//				or with nodes in other trees (if join not only roots); put also on pair's node moves 
		//		(b) or P is not a tree root (when L is not) : only when join not only roots
		//			generate node move for pair = other tree roots 
		// IV. select model moves as top models from node moves if P is not root
		//		score the top models among model moves by likelihood
		
		if (model_moves.isEmpty())
			return 0.0;
		
//		CandidateModel alt_model = model_moves.deleteLeast();
		CandidateMove alt_model=null;
		// make sure scoring is up to date for the first guy 
		System.out.println("#**MLP.iM nmoves "+model_moves.size());
		int rescored_moves=0;
		while (!model_moves.isEmpty() && rescored_moves<top_scored_moves)
		{
			// make sure top guy has fresh score
			while (alt_model != model_moves.peek()) // peek() gives null if empty
			{
				alt_model=model_moves.peek();
				if (!alt_model.isValidJoin())
				{ // should never get here
					System.out.println("#**MLP.rM  invalid "+alt_model);
					model_moves.deleteLeast(); // we don't need it
					deleteNodeMove(alt_model);
				} else
				{
					double ascore = alt_model.score;
					alt_model.score = alt_model.score()-best_model.score;
					model_moves.updateOrder(alt_model);
					double diff = alt_model.score - ascore;
					System.out.println("#**MLP.rM refresh "+alt_model+"\tscore "+alt_model.score+"/"+ascore+"\tdiff "+diff);
				}
			}
			if (eager_node_moves || alt_model!=null) 
			{
				alt_model = model_moves.deleteLeast();
				System.out.println("#**MLP.rM got"+(alt_model.score>0.0?"\tRISE\t":"\tdrop\t")+alt_model);
				
				if (alt_model.score<max_score)
				{
					break; // perfect!
				} else
				{
					// iterate again, refresh score for 2nd best, and so on  
					rescored_moves++;
				}
			} else
			{
				assert model_moves.isEmpty();
			}
		}
		
		
		if (alt_model==null)
		{
			System.out.println("#**MLP.rM no more moves.");		
			return 0.0;
		} else if (alt_model.score >= max_score)
		{
			System.out.println("#**MLP.rM only bad moves "+alt_model);		
			return alt_model.score;
		}
		
		Node L = alt_model.og_left;
		Node R = alt_model.og_right;
		
		double lenL = best_model.rates.getEdgeLength(L.getIndex());
		double lenR = best_model.rates.getEdgeLength(R.getIndex());
		
		assert (lenR==Double.POSITIVE_INFINITY);
		
		assert Double.isFinite(lenL)==Double.isFinite(L.getLength());
		assert Double.isFinite(lenR)==Double.isFinite(R.getLength());
				
		assert (lenL==Double.POSITIVE_INFINITY)==(node_moves.containsKey(L));
		assert (lenR==Double.POSITIVE_INFINITY)==(node_moves.containsKey(R));
		
				
		// I. accept alt_model
		assert (alt_model.og_model == best_model);
		
		Node P = acceptJoinMove(alt_model);

		double best_old_score = best_model.score;
//		double hscore = best_model.scoreParsimony();
//		double hscore = spr_parsimony.getSankoffScore();
//		move_heuristics.put(best_model, hscore);
		best_model.score = best_model.score();
		double best_score_diff = best_model.score-best_old_score;
		
		System.out.println("#**MLP.rM best "+best_model+"\tscore upd "+best_score_diff+"\tfrom "+alt_model //+"\thscore "+hscore
					+"\t// "+NewickParser.printTree(best_model.forest));
		
		// II. remove moot node moves
		Node ourT = best_model.getTreeRoot(P);
		
		if (node_moves.containsKey(L))
		{
			// II.1 (a)
			node_moves.get(L).remove(alt_model);
			for (CandidateMove moot: node_moves.get(L))
			{
				if (moot.og_right==L) // no more pruning with L
				{
					if (node_moves.containsKey(moot.og_left))
						node_moves.get(moot.og_left).remove(moot);
					move_heuristics.remove(moot);
				}
				model_moves.remove(moot); // not necessarily there
			}
			node_moves.remove(L); // no more voting rights
		} else 
		{
			// II.1 (b)
			assert lenR==Double.POSITIVE_INFINITY;
			assert (Double.isFinite(lenL));
			
			// remove looping joins
			if (!ourT.isRoot())
			{
				Set<Node> subtreeR = new HashSet<>(R.listNodes(null, null));
				subtreeR.remove(R); // we deal with R separately 
				List<CandidateMove> looping_joins = new ArrayList<>();
				Heap<CandidateMove> our_moves = node_moves.get(ourT);
				
				for (CandidateMove moot: our_moves)
				{
					if (moot.og_right == ourT && subtreeR.contains(moot.og_left))
					{
						looping_joins.add(moot);
					}
				}
				for (CandidateMove moot: looping_joins)
				{
					System.out.println("#**MLP.rM remove loop "+moot);
					deleteNodeMove(moot);
//					our_moves.remove(moot);
//					if (node_moves.containsKey(moot.og_left))
//					{
//						node_moves.get(moot.og_left).remove(moot);
//					}
					model_moves.remove(moot);
				}
			}
		}
			
		assert (node_moves.containsKey(R)); 
		// II.2 
		node_moves.get(R).remove(alt_model);
		move_heuristics.remove(alt_model);
		for (CandidateMove moot: node_moves.get(R))
		{
			if (moot.og_right==R)
			{
				// no more pruning with R
				if (node_moves.containsKey(moot.og_left))
				{
					node_moves.get(moot.og_left).remove(moot);
				}
				move_heuristics.remove(moot);
			}
			model_moves.remove(moot);
		} 
		node_moves.remove(R);
			
		// generate new pairings 
		if (!ourT.isRoot()) // otherwise the phylogeny is fully resolved, and we are done  
		{
			addForestMoves(P); // III
			updateModelMoves();	// IV	
		}
		return best_score_diff;
	
	}
	
	/**
	 * Refreshes the SPR phylogeny 
	 * when reference model changes, and the
	 * node moves heaps where this move appears. 
	 * 
	 * @param move
	 */
	private void refreshNodeMove(CandidateMove move)
	{
		assert (move.isValidJoin());
		
		//move.makeJoin(); // don't need the phylogeny bc spr_parsimony builds it on the fly
		// double mscore = spr_parsimony.getDiffScore(move.og_left.getIndex(), move.og_right.getIndex()); call_MP++;
		move_heuristics.put(move); // recomputes heuristic score // , mscore);
		
		move.version = best_model.version;
		if (move.forest != null) // new phylogeny here
		{
			move.forest = null;
			move.rates = null;
		}
		
		if (node_moves.containsKey(move.og_left))
		{
			node_moves.get(move.og_left).updateOrder(move);
		}
		if (node_moves.containsKey(move.og_right))
		{
			node_moves.get(move.og_right).updateOrder(move);
		}
	}

	/**
	 * Removes a candidate move from {@link #node_moves} and 
	 * {@link #move_heuristics}. 
	 * 
	 * @param moot
	 */
	private void deleteNodeMove(CandidateMove moot)
	{
		if (node_moves.containsKey(moot.og_left))
			node_moves.get(moot.og_left).remove(moot);
		if (node_moves.containsKey(moot.og_right))
			node_moves.get(moot.og_right).remove(moot);
		move_heuristics.remove(moot);
	}
	
	/**
	 * One iteration for model walking:
	 *    make the best SPR move and update the 
	 *    set of possible node moves. 
	 * 
	 * @param improve_threshold max accepted move score 
	 * @return 0.0 if no good move, &ge; max_score if only bad models, &lt; max_score if successful
	 */
	private double improveModel(double improve_threshold)
	{
		long T0 = System.currentTimeMillis();
		// I. accept the move 
		// II. remove obviated moves from model moves and node moves 
		// 	1. all moves involving the pruneParent (that may disappear if was binary)
		//  2. moves with prune in graft ancestors, n graft in prune's subtree 
		//  3. if pruneParent was root and deleted, then pruneSibling becomes root:
		//		all moves involving S

		CandidateMove alt_model=null;
		// make sure scoring is up to date for the first guy 
		System.out.println("#**MLP.iM nmoves "+model_moves.size());
		int rescored_moves=0;
		while (!model_moves.isEmpty() && rescored_moves<top_scored_moves)
		{
			// make sure top guy has fresh score
			while (alt_model != model_moves.peek()) // peek() gives null if empty
			{
				alt_model=model_moves.peek();
				if (!alt_model.isValidJoin()) // should not happen
				{
					System.out.println("#**MLP.iM invalid "+alt_model);
					model_moves.deleteLeast(); // we dont need it
				} else // valid join: alt_model.score is set, but alt_model.rates and .forest may be null
				{
					double ascore = alt_model.score;
					alt_model.score = alt_model.score()-best_model.score;
					model_moves.updateOrder(alt_model);
					double diff = alt_model.score - ascore;
					System.out.println("#**MLP.iM refresh "+alt_model+"\tscore "+alt_model.score+"/"+ascore+"\tdiff "+diff);
				}
			}
			if (alt_model != null)
			{
				alt_model = model_moves.deleteLeast();
				System.out.println("#**MLP.iM got"+(alt_model.score>0.0?"\tRISE\t":"\tdrop\t")+alt_model);
				
				if (alt_model.score<improve_threshold)
				{
					break;
				} else
				{
					// iterate again, refresh score for 2nd best, and so on  
					rescored_moves++;
				}
			}
		}
		
		if (alt_model==null)
		{
			System.out.println("#**MLP.iM no more moves.");		
			return 0.0;
		} else if (alt_model.score >= improve_threshold)
		{
			System.out.println("#**MLP.iM only bad moves "+alt_model);		
			return alt_model.score;
		}
		
		Node L = alt_model.og_left;  // regraft node 
		Node R = alt_model.og_right; // prune node 
		
		double lenL = L.getLength();
		double lenR = R.getLength();
		
		assert Double.isFinite(lenL)==
			Double.isFinite(best_model.rates.getEdgeLength(L.getIndex()));
		assert Double.isFinite(lenR)==
			Double.isFinite(best_model.rates.getEdgeLength(R.getIndex()));
				
		// I. accept alt_model
		assert (alt_model.og_model == best_model);
		
		
		Node pruneParent = R.getParent();
		Node pruneSibling =null;
		if (pruneParent.getNumChildren()==2) // spr with sibling at a binary parent gives this very same tree 
		{
			if (pruneParent.getChild(0)==R)
			{
				pruneSibling = pruneParent.getChild(1);
			} else
			{
				pruneSibling = pruneParent.getChild(0);
			}
		}		
		Node P = acceptJoinMove(alt_model); // graft parent
		if (exclude_multiple_prune)
			this.accepted_prune_moves.add(R);

		double best_old_score = best_model.score;
//		double hscore = best_model.scoreParsimony();
//		double hscore = spr_parsimony.getSankoffScore();
		//move_heuristics.put(best_model, hscore);
		best_model.score = best_model.score();
		double best_score_diff = best_model.score-best_old_score;
		
		System.out.println("#**MLP.iM best "+best_model+"\tscore upd "+best_score_diff+"\tfrom "+alt_model);//+"\thscore "+move_heuristics.get(alt_model));
//		
//		// II. remove moot node moves
//		
		// II.0 remove the selected alt model
		deleteNodeMove(alt_model);

//		// II.3 if prune sibling becomes root  
//		if (pruneSibling != null && pruneSibling.isRoot())
//		{
//			List<CandidateMove> siblingMoves = new ArrayList<>(node_moves.get(pruneSibling));
//			for (CandidateMove moot: siblingMoves)
//			{
//				if (moot.og_right == pruneSibling)
//				{
////					node_moves.get(moot.og_left).remove(moot);
////					node_moves.get(pruneSibling).remove(moot);
////					move_heuristics.remove(moot);
//					deleteNodeMove(moot);
//				}
//				model_moves.remove(moot);
//			}
//		} 
//		
////		System.out.println("#**MLP.iM pP "+pruneParent+"\troot "+pruneParent.isRoot()+"\tpS "+pruneSibling);
		// II.1 moves involving pruneParent
		List<CandidateMove> pruneMoves = new ArrayList<>(node_moves.get(pruneParent));
		for (CandidateMove moot: pruneMoves)
		{
			if (pruneSibling == null)
			{
				// still there 
				refreshNodeMove(moot);
			} else
			{
				deleteNodeMove(moot);
					// and pruneParent moves will be removed
			}
			model_moves.remove(moot);
		}
		if (pruneSibling != null)
		{
			node_moves.remove(pruneParent); // since it is not in the tree anymore 
		}
//		
		Set<Node> Rsubtree = new HashSet<>(R.listNodes(null, null));
//
//		// II.1b moves involving pruneSibling
//		if (pruneSibling != null && !pruneSibling.isRoot() && pruneSibling.getParent().getNumChildren()==2)
//		{
//			Node pruneNewParent = pruneSibling.getParent();
//			Node pruneNewSibling = pruneNewParent.getChild(pruneNewParent.getChild(0)==pruneSibling?1:0);
//			List<CandidateMove> siblingMoves = new ArrayList<>(node_moves.get(pruneSibling));
//			for (CandidateMove moot: siblingMoves)
//			{
//				Node paired = moot.og_left==pruneSibling?moot.og_right:moot.og_left;
//				if (paired == pruneNewParent || paired == pruneNewSibling)
//				{
////					node_moves.get(paired).remove(moot);
////					node_moves.get(pruneSibling).remove(moot);
//					deleteNodeMove(moot);
//				}
//				model_moves.remove(moot);
//			}
//		}
//		
		// II.2 moves pruning at graft ancestors
		for (Node A: P.listAncestors(null))
		{
			if (!A.isRoot())
			{
				List<CandidateMove> ancestorMoves = new ArrayList<>(node_moves.get(A));
				for (CandidateMove moot: ancestorMoves)
				{
					if (moot.og_right==A)
					{
						Node paired = moot.og_left;
						if (Rsubtree.contains(paired))
						{
//							node_moves.get(A).remove(moot);
//							node_moves.get(paired).remove(moot);
							deleteNodeMove(moot);
						} else
						{
							// pruning at ancestor; keep unchanged
						}
					} else
					{
						// grafting at A ; ok 
						if (Rsubtree.contains(moot.og_right))
						{
							refreshNodeMove(moot);
						}
					}
					model_moves.remove(moot);
				}
			} 
		}
//
//		// II.4 no more pruning for R
		List<CandidateMove> Rmoves = new ArrayList<>(node_moves.get(R));
		for (CandidateMove moot: Rmoves)
		{
			if (moot.og_right == R)
			{
//				Node paired = moot.og_left;
//				node_moves.get(paired).remove(moot);
//				node_moves.get(R).remove(moot);
				deleteNodeMove(moot);
			} else
			{
				refreshNodeMove(moot);
			}
			model_moves.remove(moot);
		}
//		

//		// III.1 add new moves for pruning at ancestors of the old pruning position
//		Set<Node> prunables;
//		if (pruneSibling == null)
//		{
//			prunables = new HashSet<>(pruneParent.listAncestors(null));
//			prunables.add(pruneParent);
//		} else
//		{
//			prunables = new HashSet<>(pruneSibling.listAncestors(null));
//		}
//		prunables.removeAll(P.listAncestors(null));
//		prunables.remove(P);
//		prunables.removeAll(accepted_prune_moves);
//		
////		System.out.println("#**MLP.iM R "+R+"\tsubt "+Rsubtree);
////		System.out.println("#**MLP.iM psib "+pruneSibling+"\tanc "+pruneSibling.listAncestors(null));
////		System.out.println("#**MLP.iM P "+P+"\tanc "+P.listAncestors(null));
////		System.out.println("#**MLP.iM prunables "+prunables);
//		
//		
//		for (Node N: Rsubtree)
//		{
//			Node Nparent = N.getParent();
//			Node Nsib = null;
//			if (Nparent.getNumChildren()==2)
//			{
//				Nsib = Nparent.getChild(Nparent.getChild(0)==N?1:0);
//				
//			}
//			for (Node paired: prunables)
//			{
//				if (Nsib == null || (N!=Nsib && N!=Nparent))
//				{
//	//				System.out.println("#**MLP.iM add spr("+N+","+paired+")");
//					
//					CandidateMove alt = best_model.join(N.getIndex(), paired.getIndex());
//	//				double ascore = alt.scoreParsimony()-hscore;
//					double ascore = spr_parsimony.getDiffScore(N.getIndex(), paired.getIndex()); call_MP++;
//					move_heuristics.put(alt, ascore);
//					node_moves.get(N).remove(alt);
//					node_moves.get(N).add(alt);
//					
//					node_moves.get(paired).remove(alt);
//					node_moves.get(paired).add(alt);
//					System.out.println("#**MLP.iM add "+alt+"\thscore "+ascore);
//				}
//			}
//		}
//		
//		// III.2 add new moves for graft parent P
//		//if (!P.isRoot())
//		{
//			addWalkMoves(P);
//		}
		
		// Ideally, we would have deleted all the moot moves by now. 
		// But that code is subtly buggy, so we take a safe 
		// route by enumerating all legal SPR moves, and 
		// keeping the scores that are already computed. 
		recomputeWalkMoves(); 
		
		
		updateModelMoves();
		
		
		long dT = System.currentTimeMillis()-T0;
		System.out.println("#**MLP.iM done "+dT/1000.0+" sec");
		return best_score_diff;
		
	}
	
	
	
	/**
	 * Adds moves for initial model building. 
	 * 
	 * @param P new node that is the parent of a prune-regraft join pair. 
	 */
	private void addForestMoves(Node P)
	{
		// update scored by heuristics
		Heap<CandidateMove> update = move_heuristics.newHeap(); //  new Heap<>((a,b)->Double.compare(move_heuristics.get(a), move_heuristics.get(b)));
		
		assert (best_model.forest.getNode(P.getIndex())==P);
		
		double plen = P.getLength();
		assert 0.0<=plen; 
		
		Collection<Node> pairs;
		int left=-1, right=-1;
		
		pairs = new ArrayList<>();
		Node ourT = best_model.getTreeRoot(P);
		if (!ourT.isRoot())
		{
			if (Double.isFinite(plen))
			{
				left = P.getIndex(); // this guy is not to be pruned ; we had an inner join
			} else
			{
				right = P.getIndex();
				for (Node T: node_moves.keySet()) // only for non-finite plen
					pairs.addAll(T.listNodes(null, null)); // all nodes in T's subtree 
//					assert (!node_moves.containsKey(P));
			}
		}
		for (Node T: pairs)
		{
			assert (!T.isRoot() && T!=P);
			
			if (plen == Double.POSITIVE_INFINITY)
			{
				// T is any node outside of P's subtree
				left = T.getIndex();
			} else 
			{
				right = T.getIndex();
			}
			CandidateMove join =  best_model.join(left, right);
			//double hscore = spr_parsimony.getDiffScore(left, right); call_MP++;
			move_heuristics.put(join); //, hscore);
			System.out.println("#**MLP.aFM add "
					+((Double.isFinite(T.getLength())
							|| Double.isFinite(plen))?"inner":"roots")
					+"("+left+","+right+")"
					+"\t"+join+"\thscore "+move_heuristics.get(join)); //+"\tphylo "+NewickParser.printTree(join.forest));
			update.add(join);
			if (node_moves.containsKey(T))
			{
				Heap<CandidateMove> Tmoves = node_moves.get(T);
				Tmoves.remove(join);
				Tmoves.add(join);
			}
			
		}
		if (!P.isRoot())
		{
			if (plen == Double.POSITIVE_INFINITY)
			{
				node_moves.put(P, update);
//				System.out.println("#**MLP.aFM add ndmoves "+P+"\tmoves "+update.size());
			} else
			{
//				System.out.println("#**MLP.aFM not ndmoves "+P+"\tmoves "+update.size());
			}
		}
	}
	
	
	/**
	 * Adds the moves for a node: either as prune node, or as regraft node; sets up 
	 * its entry in {@link #node_moves}
	 * 
	 * @param P
	 */
	private void addWalkMoves(Node P)
	{
		Set<Node> Psubtree = new HashSet<>(P.listNodes(null, null));
		Set<Node> Pancestors = new HashSet<>(P.listAncestors(null));
		
		Heap<CandidateMove> update = move_heuristics.newHeap(); //new Heap<>((a,b)->Double.compare(move_heuristics.get(a), move_heuristics.get(b)));
//			{ // favors prune moves over graft moves
//				return (a.og_right==b.og_right)
//						?Double.compare(move_heuristics.get(a), move_heuristics.get(b))
//						:(a.og_right==P?-1:+1);
//			});
		
		Node parent = P.getParent();
		Node sibling=null;
		
//		assert (!P.isRoot());
//		assert (!Pancestors.isEmpty()); // since P is not root

		if (!P.isRoot())
		{
			if (parent.getNumChildren()==2) // spr with sibling at a binary parent gives this very same tree 
			{
				if (parent.getChild(0)==P)
				{
					sibling = parent.getChild(1);
				} else
				{
					sibling = parent.getChild(0);
				}
			}
		}
		
		for (Node N: node_moves.keySet())
		{
			if (sibling==null || (N!=parent && N!=sibling)) // otherwise SPR creates identical topology after removing unary nodes 
			{
				if (!Psubtree.contains(N) && !accepted_prune_moves.contains(P))
				{
					// N can be the target of a regraft
					int left = N.getIndex(); 
					int right = P.getIndex();
					CandidateMove spr = best_model.join(left, right);
					//double hscore = spr_parsimony.getDiffScore(left, right); call_MP++;
					move_heuristics.put(spr);//, hscore);
//					System.out.println("#**MLP.aWM add/pg spr("+left+","+right+")"
//							+"\t"+spr+"\thscore "+hscore); //  +"\tphylo "+NewickParser.printTree(spr.forest));
					update.add(spr);
					node_moves.get(N).remove(spr);
					node_moves.get(N).add(spr);			
				}
				if (!Pancestors.contains(N) && !accepted_prune_moves.contains(N))
				{
					// N can be grafted here
					if (P != N.getParent() || P.getNumChildren()>2) // otherwise SPR makes no diff 
					{
						int left = P.getIndex();
						int right = N.getIndex();
						CandidateMove spr = best_model.join(left, right);
//						double hscore = spr_parsimony.getDiffScore(left, right); call_MP++;
//						move_heuristics.put(spr, hscore);
						move_heuristics.put(spr);
						
//						System.out.println("#**MLP.aWM add/gp spr("+left+","+right+")"
//								+"\t"+spr+"\thscore "+hscore); //+"\tphylo "+NewickParser.printTree(spr.forest));
						update.add(spr);
						node_moves.get(N).remove(spr);
						node_moves.get(N).add(spr);				
					}
				}
			}
		}
		
		node_moves.put(P, update);
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#**MLP.aWM add ndmoves "+P+"\tmoves "+update.size());
	}
	
	/** 
	 * Enumerates all SPR moves and 
	 * resets candidate move structures; copying 
	 * the current scores for all moves that are still valid.  
	 * Sets {@link #node_moves}, {@link #model_moves} and {@link #move_heuristics}.  
	 * 
	 */
	void recomputeWalkMoves()
	{
		long T0 = System.currentTimeMillis();
		final Map<Node, Heap<CandidateMove>>
			updated_node_moves = new HashMap<>();
		final CandidateHeuristics //Map<CandidateMove, Double>
			updated_move_heuristics = new CandidateHeuristics(
					heuristic_byParsimony?m->m.scoreParsimony(getSPRExplorer()):new FastMLHeuristic()); //  HashMap<>();
		final Heap<CandidateMove>
			updated_model_moves = new Heap<>();
		final List<Node> postorder = best_model.forest.getRootNode().listNodes(null, null);
		
		// init the heap on heuristic score at every node 
		for (Node N: postorder)
		{
			Heap<CandidateMove> moves = updated_move_heuristics.newHeap(); // new Heap<>((a,b)->Double.compare(updated_move_heuristics.get(a), updated_move_heuristics.get(b)));
			updated_node_moves.put(N, moves);
		}
		
		// enumerate all SPR moves
		int num_scored = 0;
		int num_joins = 0;
		for (Node R: postorder) // prune node 
		{
			Heap<CandidateMove> Rmoves = updated_node_moves.get(R);
			Node P = R.getParent();
			if (P!=null && !accepted_prune_moves.contains(R))
			{ // thus can be pruned 
				int right = R.getIndex();
				List<Node> graft_nodes = R.graftNeighbors(3, postorder.size()); // all distances 
				if (P.getNumChildren()>2) // with a binary parent,  grafting there or at sibling gives identical topologies
				{
					graft_nodes.add(P);
					for (int ci=0; ci<P.getNumChildren(); ci++)
					{
						Node C = P.getChild(ci);
						if (C!=R)
							graft_nodes.add(C);
					}
				}
				
				for (Node L: graft_nodes)
				{
					int left = L.getIndex();
					CandidateMove join = best_model.join(left, right);
					// check if this was scored before 
					CandidateMove old;
					if (node_moves.containsKey(R))
						old = node_moves.get(R).get(join);
					else if (node_moves.containsKey(L))
						old = node_moves.get(L).get(join);
					else
						old = null;
					
					if (old == null)
					{
						//double jscore = spr_parsimony.getDiffScore(left, right);
						updated_move_heuristics.put(join); //, jscore);
						Rmoves.add(join);
						updated_node_moves.get(L).add(join);
						++ num_scored;
					} else
					{
						// keep the old version
						double hscore = move_heuristics.get(old);
						updated_move_heuristics.put(old, hscore); // do not recompute hscore
						Rmoves.add(old);
						updated_node_moves.get(L).add(old);
						
						if (model_moves.contains(old))
							updated_model_moves.add(old);
					}
					++ num_joins;
				} // for L

			} // if R can be pruned
		} // for R 
		long dT = System.currentTimeMillis()-T0;
			System.out.println("#**MLP.rWM "+best_model+"\tnjoins "+num_joins+"\thscored "+num_scored+"\t"+dT/1000.0+" sec");
//		for (Node R: postorder)
//		{
//			Heap<CandidateMove> Rheap = updated_node_moves.get(R);
//			if (node_moves.containsKey(R))
//			{
//				System.out.println("#**MLP.rWM seen "+R+"\tnmoves "+Rheap.size()
//						+"/"+node_moves.get(R).size());
//				for (CandidateMove old: node_moves.get(R))
//				{
//					if (!Rheap.contains(old))
//					{
//						System.out.println("#**MLP.rWM badmove "+R+"\t"+old+"\tvalid "+old.isValidJoin());
//					}
//				}
//				Heap<CandidateMove> old_moves = node_moves.get(R);
//				for (CandidateMove join: Rheap)
//				{
//					if (!old_moves.contains(join))
//					{
//						System.out.println("#**MLP.rWM miss "+R+"\t"+join);					
//					}
//				}
//			} else
//			{
//				System.out.println("#**MLP.rWM unseen "+R+"\tnmoves "+Rheap.size());
//			}
//		}
		this.node_moves = updated_node_moves;
		this.move_heuristics = updated_move_heuristics;
		this.model_moves = updated_model_moves;
	}
	
	
	/**
	 * Heap ordered by heuristic score from top moves at each node.
	 * 
	 * @param hscore_threshold
	 * @param top_per_node
	 * @return
	 */
	private Heap<CandidateMove> getTopNodeMoves(double hscore_threshold, int top_per_node)
	{
		Heap<CandidateMove> top_moves = move_heuristics.newHeap();  
		int num_refreshed=0;
		for (Node T: node_moves.keySet())
		{
			if (accepted_prune_moves==null || !accepted_prune_moves.contains(T))
			{
				Heap<CandidateMove> models = node_moves.get(T);
				for (int m=0; m<top_scored_moves && !models.isEmpty(); m++)
				{
					CandidateMove moveT = null;
					while (moveT != models.peek())
					{
						moveT = models.peek();
						if (!moveT.isValidJoin())
						{
//							System.out.println("#**MLP.gTNM purging "+moveT);
							models.deleteLeast(); // throw it away
							deleteNodeMove(moveT);
							model_moves.remove(moveT); 
							moveT = null;
						} else if (eager_node_moves && moveT.isStale())
						{
//							System.out.println("#**MLP.uMM refresh "+moveT);
							refreshNodeMove(moveT); // updates order of moveT in its two heaps
							num_refreshed++;
						} 
					}
					if (moveT != null)
					{
						if (move_heuristics.get(moveT)>=hscore_threshold) break; // won't look at them with high hscore
						
						moveT = models.deleteLeast();
						top_moves.add(moveT);
						
						
//						moveT.score = moveT.fastScoreML()-best_model.score;
//						fast_moves.add(moveT);
						
						//System.out.println("#**MLP.uMM fast "+moveT);
					} else
					{
						//assert models.isEmpty();
					}
				}
			}
		} // for all in node_moves
		
		for (CandidateMove moveT: top_moves)
		{
			// refile on node heaps; duplicates are detected in add() 
			if (node_moves.containsKey(moveT.og_left))
				node_moves.get(moveT.og_left).add(moveT);
			if (node_moves.containsKey(moveT.og_right))
				node_moves.get(moveT.og_right).add(moveT);
		}

		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#**MLP.gTNM topmoves "+top_moves.size()+"\ttop "+top_moves.peek()+"\thscore "+move_heuristics.get(top_moves.peek())+"\trefreshed "+num_refreshed);
		
		return top_moves;
		
	}
	
	
	
	/**
	 * Adds and scores the top candidates for each node move.
	 */
	private void updateModelMoves()
	{
		updateModelMoves(Double.POSITIVE_INFINITY);
	}
	
	/**
	 * Adds and scores the top candidates for each node. 
	 * Going through movable nodes with entry in {@link #node_moves}, 
	 * adds at most {@link #top_scored_moves} best (by heuristic score) 
	 * from each to grow a collection 
	 * of top moves to consider, while purging out the invalid joins encountered. 
	 * At most twice {@link #top_scored_moves} top moves are added to 
	 * {@link #model_moves}, and scored properly. 
	 *  
	 * @param hscore_threshold
	 */
	private void updateModelMoves(double hscore_threshold)
	{
		
		Heap<CandidateMove> top_moves; // = move_heuristics.newHeap();  // new Heap<>((a,b)->Double.compare(move_heuristics.get(a), move_heuristics.get(b)));
//		Heap<CandidateMove> fast_moves = new Heap<>();
		
		top_moves = getTopNodeMoves(hscore_threshold, top_scored_moves);
//		
//		int num_refreshed = 0;
//		for (Node T: node_moves.keySet())
//		{
//			if (accepted_prune_moves==null || !accepted_prune_moves.contains(T))
//			{
//				Heap<CandidateMove> models = node_moves.get(T);
//				for (int m=0; m<top_scored_moves && !models.isEmpty(); m++)
//				{
//					CandidateMove moveT = null;
//					while (moveT != models.peek())
//					{
//						moveT = models.peek();
//						if (!moveT.isValidJoin())
//						{
////							System.out.println("#**MLP.uMM purging "+moveT);
//							models.deleteLeast(); // throw it away
//							deleteNodeMove(moveT);
////							Node paired = (moveT.og_left==T?moveT.og_right:moveT.og_left);
////							if (node_moves.containsKey(paired))
////								node_moves.get(paired).remove(moveT);
//							model_moves.remove(moveT); 
//							moveT = null;
//						} else if (eager_node_moves && moveT.isStale())
//						{
////							System.out.println("#**MLP.uMM refresh "+moveT);
//							refreshNodeMove(moveT); // updates order of moveT in its two heaps
//							num_refreshed++;
//						} 
//					}
//					if (moveT != null)
//					{
//						if (move_heuristics.get(moveT)>=hscore_threshold) break; // won't look at them with high hscore
//						
//						moveT = models.deleteLeast();
//						top_moves.add(moveT);
//						
//						
////						moveT.score = moveT.fastScoreML()-best_model.score;
////						fast_moves.add(moveT);
//						
//						//System.out.println("#**MLP.uMM fast "+moveT);
//					} else
//					{
//						//assert models.isEmpty();
//					}
//				}
//			}
//		}

		int num_refreshed=0; // if eager_node_moves, then getTopNodeMoves did the refreshing
		
		int num_scored = top_scored_moves*2;//  *top_scored_moves;
		int num_top_moves= top_moves.size();

//		System.out.println("#**MLP.uMM fast moves "+fast_moves.size());
		
//		while (num_scored>0 && !fast_moves.isEmpty()) // 
		while (num_scored>0 && !top_moves.isEmpty())
		{
			CandidateMove moveT;
			if (eager_node_moves)
			{
				moveT = top_moves.deleteLeast();
			} else
			{
				moveT = null;
				while (moveT != top_moves.peek())
				{
					moveT = top_moves.peek();
					if (moveT.isStale())
					{
						double mscore = move_heuristics.get(moveT);
						refreshNodeMove(moveT);
						top_moves.updateOrder(moveT);
						double diff = move_heuristics.get(moveT)-mscore;
//						System.out.println("#**MLP.uMM refresh/T "+num_scored+"\t"+moveT+"\thscore "+move_heuristics.get(moveT)+"/"+mscore+"\tdiff "+diff);
						num_refreshed++;
					}
				}
				
				// make sure there is no stale top model for the same nodes 
				List<CandidateMove> top_list = new ArrayList<>(top_moves);
				CandidateMove moveN = null;
				while (moveN != moveT)
				{
					for (CandidateMove move2: top_list)
						if (moveT.sharedJoinNodes(move2)>0 && move2.isStale()) // moveT is not stale
						{
							double mscore = move_heuristics.get(move2);
							refreshNodeMove(move2);
							top_moves.updateOrder(move2);
							double diff = move_heuristics.get(move2)-mscore;
//							System.out.println("#**MLP.uMM refresh/m "+num_scored+"\t"+move2+"\thscore "+move_heuristics.get(move2)+"/"+mscore+"\tdiff "+diff);
							num_refreshed++;
						}
					moveN = top_moves.peek();
					if (moveN != moveT)
					{
						if (PRINT_OPTIMIZATION_MESSAGES)
							System.out.println("#**MLP.uMM refresh/N "+num_scored+"\t"+moveN+"\twas "+moveT);
						
						moveT = moveN;
						moveN = null;
						assert !moveT.isStale();
					}
				}
				moveT = top_moves.deleteLeast();
			}
			
			if (move_heuristics.get(moveT)>=hscore_threshold) break; // won't look at them with high hscore
//			CandidateMove moveT = fast_moves.deleteLeast();
			
			//System.out.println("#**MLP.uMM topmove "+moveT);
			
				
			if (accepted_prune_moves==null 
					|| !accepted_prune_moves.contains(moveT.og_right))
			{
				--num_scored;
				
				if (model_moves.contains(moveT))
				{
//					if (model_moves.get(moveT).score>moveT.score)
//					{
//						moveT.score = moveT.score()-best_model.score;
//						model_moves.remove(moveT);
//						model_moves.add(moveT);
//						System.out.println("#**MLP.uMM rescoring "+num_scored+"\t"+moveT+"\thscore "+move_heuristics.get(moveT));
//					} else		
//					{
						// should we recalculate the score? 
						// nah: it should typically only get worse; top model moves will be rescored anyway
						if (PRINT_OPTIMIZATION_MESSAGES)
							System.out.println("#**MLP.uMM notscoring "+num_scored+"\t"+moveT+"\thscore "+move_heuristics.get(moveT));
//					}
				} else
				{
					// score this guy and add to main pool
					moveT.score = moveT.score()-best_model.score;
					model_moves.add(moveT);
					if (PRINT_OPTIMIZATION_MESSAGES)
						System.out.println("#**MLP.uMM scoring "+num_scored+"\t"+moveT+"\thscore "+move_heuristics.get(moveT));
				}
			} else
			{
				// should not get here 
				System.out.println("#**MLP.uMM no more pruning of "+moveT.og_right+"\t"+moveT);
			}
		}
		if (PRINT_OPTIMIZATION_MESSAGES)
			System.out.println("#**MLP.uMM mmoves "+model_moves.size()+"\ttop "+model_moves.peek()+"\trefreshed "+num_refreshed+"/"+num_top_moves);
	}
	
	
	/**
	 * Makes a a number of random NNI moves
	 * 
	 * @param num_steps
	 */
	private void stochasticNNI(int num_steps)
	{
		for (int s=0; s<num_steps; s++)
		{
			List<CandidateMove> nni_moves = best_model.getNearestNeighborMoves();
			int i = RND.nextInt(nni_moves.size());
			CandidateMove alt_model = nni_moves.get(i);
			alt_model.score = alt_model.scoreML(0); // no optimization
			System.out.println("#**MLP.sNNI "+s+"/"+num_steps+"\t"+alt_model);
			acceptJoinMove(alt_model);
		}
		best_model.score = best_model.score();
	}
	
	
	/**
	 * 
	 * 
	 * @param increase_threshold
	 * @deprecated
	 */
	private void exploreNNI(double increase_threshold)
	{
		class Steps 
		{
			
			Steps(CandidateMove og, Parsimony.SPRExplorer spr)
			{ 
				this.origin = og;
				this.spr = spr;
				initDataStructures();
			}
			
			Steps(CandidateMove og)
			{
				this(og, new Parsimony.SPRExplorer(og.forest, utable.mappedToTree(og.forest), DEFAULT_PARSIMONY_GAIN, DEFAULT_PARSIMONY_LOSS, DEFAULT_PARSIMONY_DUPLICATION));
			}
			
			
			private final CandidateMove origin;
			private final Parsimony.SPRExplorer spr;
			private Map<Node, Map<Node, CandidateMove>> nni_moves; // prune -> {regraft -> move}
			private CandidateHeuristics heuristic_scores; // Map<CandidateMove, Double> heuristic_scores;
			private Heap<CandidateMove> sorted_moves;
			
			private void initDataStructures()
			{
				this.heuristic_scores =  new CandidateHeuristics(
						heuristic_byParsimony? m->m.scoreParsimony(spr):new FastMLHeuristic()); //new HashMap<>();
				this.nni_moves = new HashMap<>();
				for (CandidateMove nni: origin.getNearestNeighborMoves())
				{
					Node R = nni.og_right;
					Map<Node, CandidateMove> Rmoves;
					if (nni_moves.containsKey(R))
						Rmoves = nni_moves.get(R);
					else
					{
						Rmoves = new HashMap<>();
						nni_moves.put(R,  Rmoves);
					}
					Node L = nni.og_left;
					Rmoves.put(L, nni);
				}
			}
			
			private CandidateMove before(CandidateMove after)
			{
				CandidateMove match;
				Node R = origin.getOriginal(after.og_right);
				if (R==null) match = null;
				else 
				{
					Map<Node, CandidateMove> Rmoves = nni_moves.get(R);
					if (Rmoves==null) match = null;
					else
					{
						Node L = origin.getOriginal(after.og_left);
						if (L==null) match = null;
						else match = Rmoves.get(L);
					}
				}
				return match;
			}
			
			void computeHscores(Steps from)
			{
				sorted_moves = heuristic_scores.newHeap(); // new Heap<>((a,b)->Double.compare(heuristic_scores.get(a), heuristic_scores.get(b)));
				
				int num_scored = 0;

				for (Node R: nni_moves.keySet())
				{
					
					for (CandidateMove move: nni_moves.get(R).values())
					{
						double mscore;
						if (from == null)
						{
							heuristic_scores.put(move);
							// mscore = spr.getDiffScore(move.og_left.getIndex(), R.getIndex());
							++ num_scored;
						} else
						{
							CandidateMove match = before(move);
							if (match == null)
							{
								++ num_scored;
								heuristic_scores.put(move); // mscore = spr.getDiffScore(move.og_left.getIndex(), R.getIndex());
							}
							else 
							{
								mscore = from.heuristic_scores.get(match);
								if (match.rates!=null)
									move.score = match.score;
								heuristic_scores.put(move, mscore);							}
						}
						sorted_moves.add(move);
					}
				}
				
				System.out.println("#**MLP.eNNI/S.cHs "+this+"\tscored "+num_scored);
			}
			@Override
			public String toString()
			{
				return "S["+origin.id+":"+sorted_moves.size()+"]";
			}
		}
		
		
		long T0 = System.currentTimeMillis();
		Map<CandidateMove, Steps> current_steps = new HashMap<>();
		// init with starting model
		Steps first = new Steps(best_model, spr_parsimony);
		first.computeHscores(null);
		current_steps.put(best_model,first);
		
		int moves_per_current = 32-Integer.numberOfLeadingZeros(best_model.forest.getNumNodes()-1); 
		// ~ log_2(N) for N nodes
		Heap<CandidateMove> overall_best = new Heap<>((a,b)->b.compareTo(a));
		
		int max_steps = moves_per_current;
		
		System.out.println("#**MLP.eNNI threshold "+increase_threshold+"\tcurn "+moves_per_current+"\tmaxs "+max_steps+"\tnmoves "+first.sorted_moves.size());
		
		
		int step = 0;
		while(step<max_steps && !current_steps.isEmpty())
		{
			step++;
			// current_steps.size()<= top_moves;
			Heap<CandidateMove> fast_moves = new Heap<>();
			for (Steps S: current_steps.values() )
			{
				for (int m=0; m<moves_per_current && !S.sorted_moves.isEmpty(); m++)
				{
					CandidateMove move = S.sorted_moves.deleteLeast();
					double mscore = move.score;
					move.score = move.scoreFastML()-best_model.score;
					fast_moves.add(move);
					double mdiff = move.score - mscore;
					System.out.println("#**MLP.eNNI step "+step+"\tS "+S+"/"+m+"\tfast "+move+"\tmscore "+mscore+"\tdiff "+mdiff);
				}
			}
			// select the best of them ; with different topologies 
			Map<CandidateMove, TreeComparator> comparators = new HashMap<>();
			List<CandidateMove> top_moves = new ArrayList<>();
			for (int t=0; t<top_scored_moves && !fast_moves.isEmpty() ;)
			{
				CandidateMove top = fast_moves.deleteLeast();
				for (CandidateMove better: comparators.keySet())
				{
					TreeComparator cmp = comparators.get(better);
					if (cmp.sameTopology(top.forest))
					{
						System.out.println("#**MLP.eNNI sametopo "+top+"\tbetter "+better);
						top = null;
						break;
					}
				}
				if (top != null)
				{
					System.out.println("#**MLP.eNNI step "+step+"\ttop/"+t+" "+top);
					
					TreeComparator cmp = new TreeComparator(top.forest);
					comparators.put(top, cmp);
					top_moves.add(top);
					t++;
				}
			}
			
			Map<CandidateMove, Steps> next_steps = new HashMap<>();
			for (CandidateMove top: top_moves)
			{
				top.score = top.score()-best_model.score;
				
				if (overall_best.size()<top_scored_moves 
						|| top.score<overall_best.peek().score)
				{
					if (overall_best.size()==top_scored_moves)
						overall_best.deleteLeast();
					overall_best.add(top);
				}  else if (top.score > overall_best.peek().score+increase_threshold)
				{
					System.out.println("#**MLP.eNNI step "+step+"\ttop "+top+"\tnonext");
					break;
				} 
				Steps next = new Steps(top);
				next.computeHscores(current_steps.get(top.og_model));
				next_steps.put(top, next);
				
				System.out.println("#**MLP.eNNI step "+step+"\ttop "+top+"\tnext "+next);
			}
			current_steps = next_steps;
		}
		long dT = System.currentTimeMillis()-T0;
		
		System.out.println("#**MLP.eNNI done "+dT/1000.0+" sec");
		
		while (!overall_best.isEmpty())
		{
			CandidateMove best = overall_best.deleteLeast();
			List<CandidateMove> best_chain = new ArrayList<>();
			CandidateMove move = best;
			while (move != best_model)
			{
				best_chain.add(move);
				move = move.og_model;
			}
			for (int j=best_chain.size(); j>0; )
			{
				--j;
				move = best_chain.get(j);
				System.out.println("#**MLP.eNNI best/"+overall_best.size()+"."+j+"\t"+move);
			}
		}
		
		// pick the best, and accept the moves leading there
		// rescore best_model
		// reset spr_parsimony
		// recompute SPR moves:
		//   updated_node_moves[]; updated move_heuristics
		// 
		// for all prune-regraft pairs: <- O(n^2) 
		//     generate new join() move 
		//     copy hscore if it is in node_moves[]
		//		else compute hscore  
		//     add to updated_node_moves and updated_move_heuristics 
		// swap XX for updated_XX 
		// clear model_moves 
		// updateModelMoves
	}
	
	
	
	
	/**
	 * Arranges rate parameters per edge in a double array.
	 * 
	 * @param model
	 * @return
	 */
	private static Map<Node, double[]> modelRates(TreeWithRates model)
	{
		Phylogeny forest = (Phylogeny)model.getTree();
		Map<Node, double[]> model_rates = new HashMap<>();
		for (int node=0; node<forest.getNumNodes(); node++)
		{
			Node N = forest.getNode(node);
			double[] rate_params = new double[4];
			rate_params[PARAMETER_GAIN] = model.getGainRate(node);
			rate_params[PARAMETER_LOSS] = model.getLossRate(node);
			rate_params[PARAMETER_DUPLICATION] = model.getDuplicationRate(node);
			rate_params[PARAMETER_LENGTH] = model.getEdgeLength(node);
			model_rates.put(N,  rate_params);
		}
		return model_rates;
	}
	
	private void setMLParameters(ML M)
	{
		M.setCalculationWidth(calculation_width_absolute, calculation_width_relative);
		M.setMinimumObservedCopies(min_copies);
	}
	
	
	private ML parameterOptimizer(TreeWithRates rates, ProfileTable table)
	{
		ML M;
		if (optimize_EM)
			M = new DirectEM(rates, table);
		else
			M = optimize_rates?new MLRates(rates, table):new MLDistribution(rates, table);
		
		setMLParameters(M);
		return M;
	}
	
	private boolean optimize_EM = false;
	
	private boolean optimize_rates = false; // if true, use MLRates, else MLDistribution
	
	private boolean scoring_byML = true; // false (parsimony scoring) untested 
	
	private boolean heuristic_byParsimony = true ; //false; // false (fastML heuristic) untested
	
	private boolean fit_parsimony = true;
	
	
	private static int candidate_model_count = 0;
	
	private static final double DEFAULT_PARSIMONY_GAIN = 3.0;
	private static final double DEFAULT_PARSIMONY_LOSS = 2.0;
	private static final double DEFAULT_PARSIMONY_DUPLICATION = 1.5;
	
	private double pty_gain = DEFAULT_PARSIMONY_GAIN;
	private double pty_loss = DEFAULT_PARSIMONY_LOSS;
	private double pty_dup = DEFAULT_PARSIMONY_DUPLICATION;

	/**
	 * Uses {@link DirectEM#fitParsimonyPenalty()} to set penalties for parsimony.
	 * Nothing happens if {@link #optimize_EM} is false. 
	 * 
	 * @param model_rates the model used to calculate likelihoods 
	 * @param mapped_table profile table mapped to the model_rates' phylogeny leaf order
	 */
	private void fitParsimonyPenalties(TreeWithRates model_rates, UniqueProfileTable mapped_table)
	{
		if (optimize_EM)
		{
			DirectEM M = new DirectEM(model_rates, mapped_table, true);
			setMLParameters(M);
			double[] pty = M.fitParsimonyPenalty();
			if (Double.isFinite(pty[0]) && Double.isFinite(pty[1]) && Double.isFinite(pty[2]))
			{
				pty_gain = pty[PARAMETER_GAIN];
				pty_loss = pty[PARAMETER_LOSS];
				pty_dup = pty[PARAMETER_DUPLICATION];
			} else
			{
				System.out.println("#**MLP.fitP failed: "+Arrays.toString(pty)
					+"; staying with "+pty_gain+","+pty_loss+","+pty_dup);
			}
		}
	}
	
	
	private long time_LL=0L;
	private int call_LL=0;
	private int call_fLL = 0;
	private long time_fLL = 0L;
	private int call_MP=0;
	private long time_gradient=0L;
	private int call_gradient = 0;
	
	/**
	 * Heuristic scores stored for candidate moves. 
	 * 
	 * @author csuros
	 *
	 */
	private class CandidateHeuristics extends HashMap<CandidateMove, Double>
	{
		CandidateHeuristics(Function<CandidateMove,Double> heuristic)
		{
			super();
			this.heuristic = heuristic;
		}
		
		private final Function<CandidateMove,Double> heuristic;
		
		public Double put(CandidateMove move)
		{
			double hscore = heuristic.apply(move);
			return super.put(move, hscore);
		}
		
//		@Override
//		public Double put(CandidateMove move, Double ignore)
//		{
//			throw new UnsupportedOperationException("Do not calculate the score externally.");
//		}
		
		/**
		 *  
		 * New heap ordered by scores stored here. 
		 *  
		 * @return
		 */
		public Heap<CandidateMove> newHeap()
		{
			Heap<CandidateMove> newHeap = new Heap<>((a,b)->Double.compare(this.get(a), this.get(b)))	; 			
			return newHeap;
		}
	}
	
	private class FastMLHeuristic implements Function<CandidateMove,Double>
	{
		@Override
		public Double apply(CandidateMove M)
		{
			double fs = M.scoreFastML();
			double diff = fs-best_model.score;
			return diff;
		}
	}
	
	/**
	 * Candidate moves: SPR (join), edge contraction, and rerooting. 
	 * The {@link CandidateMove#score} variable serves to cache
	 * the ML score (set externally}, and the 
	 * model {@link CandidateMove#rates} is set 
	 * if it is optimized by either {@link CandidateMove#scoreFastML(int)}
	 * or {@link CandidateMove#scoreML(int)}. 
	 * {@link CandidateMove#version} serves to catch the updates
	 * to the reference model. Accepting a 
	 * candidate move updates the reference model, and increases its version. 
	 * New SPR moves are instanciated by the {@link CandidateMove#join(int, int)},
	 * {@link CandidateMove#contract(int)} and {@link CandidateMove#reroot(int)}
	 * methods.   
	 * 
	 * @author csuros
	 *
	 */
	private class CandidateMove implements Comparable<CandidateMove>
	{
		CandidateMove(TreeWithRates rates)
		{
			this.rates = rates;
			this.forest = (Phylogeny)rates.getTree();
			this.id = candidate_model_count++;
			this.version = 0;
		}
		private CandidateMove()
		{
			this.forest = null;
			this.rates = null;
//			model_rates = modelRates(base);
			this.id = candidate_model_count++;
			this.version = 0;
		}
		CandidateMove(CandidateMove copy)
		{
			this.id = copy.id;
			this.version = copy.version;
			
			this.forest = new Phylogeny(copy.forest);
			Map<Node, double[]> model_rates = modelRates(copy.rates);
			this.rates = new TreeWithRates(this.forest);
			for (int node = forest.getNumNodes(); 0<node ;)
			{
				--node;
				double[] rate_params = model_rates.get(copy.forest.getNode(node));
				this.rates.setGainRate(node, rate_params[PARAMETER_GAIN]);
				this.rates.setLossRate(node, rate_params[PARAMETER_LOSS]);
				this.rates.setDuplicationRate(node, rate_params[PARAMETER_DUPLICATION]);
				this.rates.setEdgeLength(node, rate_params[PARAMETER_LENGTH]);
				
				double oldlen = this.forest.getLength(node);
				double copylen = copy.forest.getLength(node);
				if (oldlen != copylen)
				{
					this.forest.getNode(node).setLength(copylen);
					//System.out.println("#**MLP.CM init badlen "+node+"\told "+oldlen+"\tcopy "+copylen+"\t"+forest.getIdent(node));
				}
			}
			this.og_model = copy;
			this.score = copy.score;
		}
		
		/**
		 * Creates a frozen snapshot of this 
		 * move, with its own phylogeny {@link #forest} 
		 * and rates {@link #rates}. This 
		 * instance an be safely used and updated further. 
		 * The snapshot {@link #equals(Object)} 
		 * the original for SPR moves, but possibly 
		 * with different rates and score.  
		 * 
		 * @return
		 */
		CandidateMove snapshot()
		{
			CandidateMove snapshot = new CandidateMove(this);
			snapshot.og_model = this.og_model;
			snapshot.og_left = this.og_left;
			snapshot.og_right = this.og_right;
			snapshot.og_edge = this.og_edge;
			
			return snapshot;
		}
		private Phylogeny forest;
		private int id;
		private int version; 
		private TreeWithRates rates;
		
		/**
		 * Set externally by enclosing class's methods. 
		 */
		private double score; 
		
		/** 
		 * The model that created this by {@link #join(int, int)}
.		 */
		private CandidateMove og_model=null;
		/*
		 * Nodes within {@link #og_model} 
		 */
		/**
		 * Regraft node in join move
		 */
		private Node og_left =null;
		/**
		 * Prune node in join move
		 */
		private Node og_right =null;
		/**
		 * Set by {@link #placeRoot(double)} or {@link #contractEdges(double)}
		 */
		private Node og_edge = null;
		
		/**
		 *  Mapping from nodes in our {@link #forest} nodes to those in {@link #og_model};
		 *  set by {@link #makeJoinPhylo()}. 
		 */
		private Map<Node, Node> og_nodes;		
		
		/**
		 * Node in {@link #og_model} phylogeny that 
		 * corresponds to a node in our {@link #forest}, 
		 * created by {@link #makeJoinPhylo()}.
		 * 		  
		 * @param node_after_join
		 * @return
		 */
		Node getOriginal(Node node_after_join)
		{
			return og_nodes.get(node_after_join);
		}

		/**
		 * Number of trees in the forest model: root children with infinite edge length.
		 * @return
		 */
		public int getNumTrees()
		{
			int nc = 0;
			Node R = forest.getRootNode();
			for (int ci=0; ci< R.getNumChildren(); ci++)
			{
				Node C = R.getChild(ci);
				if (C.getLength()==Double.POSITIVE_INFINITY)
					nc++;
			}
			assert (nc==0 || nc==R.getNumChildren());
			return nc;
		}
		
		/**
		 * Finds the tree in the forest where a node belongs. 
		 * The tree root is one of the forest root's children 
		 * (with infinite edge length), or, if 
		 * the forest is fully resolved, then  
		 * the forest root. 
		 * 
		 * @param N
		 * @return N's tree
		 */
		public Node getTreeRoot(Node N)
		{
			if (N.isRoot()) return N;
			while (!N.getParent().isRoot()) N=N.getParent();
			// now check if the forest is fully resolved
			if (Double.isFinite(N.getLength())) N=N.getParent();
			return N;
		}
		
		/**
		 * Tree roots in the forest model: root children with infinite edge lengths. 
		 * Assumes that the move was scored ({@link #rates} is not null).  
		 *  
		 * @return never null, but 0-length array if not a forest anymore 
		 */
		public Node[] getTreeRoots()
		{
			Node[] tree_roots = new Node[getNumTrees()];
			Node R = forest.getRootNode();
			int ri=0; 
			for (int ci=0; ci<R.getNumChildren(); ci++)
			{
				Node C = R.getChild(ci);
				double clen = rates.getEdgeLength(C.getIndex());
				if (clen==Double.POSITIVE_INFINITY)
					tree_roots[ri++]=C;
			}
			return tree_roots;
		}
		
		/**
		 * Whether the stored phylogeny is still OK. 
		 * {@link #scoreParsimony()} and {@link #scoreML()} 
		 * call to update the {@link #forest} before 
		 * computing scores. 
		 * 
		 * @return
		 */
		public boolean isStale()
		{
			return og_model != null && og_model.version>this.version;
		}
		
		
	
		
		private double scoreParsimony()
		{
			return scoreParsimony(null);
//			if (og_left != null)
//				makeJoinPhylo();
//			Parsimony P = new Parsimony(forest, utable.mappedToTree(forest));
//			double score = P.getSankoffScore();
//			return score;
		}
		double scoreParsimony(Parsimony.SPRExplorer spr)
		{
			double score;
			if (spr==null)
			{
				if (og_left != null)
					makeJoinPhylo();
				Parsimony P = new Parsimony(forest, utable.mappedToTree(forest));
				score = P.getSankoffScore();
			} else
			{
				int left = og_left.getIndex();
				int right = og_right.getIndex();
				score = spr.getDiffScore(left, right); call_MP++;				
			}
			return score;
		}
		double score()
		{
			return scoring_byML?scoreML():scoreParsimony();
		}

		double scoreFastML()
		{
			long T0 = System.currentTimeMillis();
			double s =scoreFastML(opt_rounds);
			time_fLL += System.currentTimeMillis()-T0;
			call_fLL++;
			return s;
		}
		
		/**
		 * Optimizes just a few edges around the prune and regraft positions, or
		 * around the contraction position. 
		 * (It may be that many other edges should change for ML.)
		 * 
		 * @param iter
		 * @return
		 */
		private double scoreFastML(int iter)
		{
			//MLDistribution 
			ML M;
			if (og_left == null)
			{
				assert (og_edge != null);
				// contraction
				Node P = contract();
				M = new MLDistribution(rates, utable.mappedToTree(forest)); //parameterOptimizer(rates, utable.mappedToTree(forest));  //
				
				
				Set<Node> optimized_edges = new HashSet<>();
				optimized_edges.add(P);
				for (int ci=0; ci<P.getNumChildren(); ci++)
				{
					optimized_edges.add(P.getChild(ci));
				}
				for (int node=0; node<forest.getNumNodes(); node++)
				{
					Node N = forest.getNode(node);
					if (!optimized_edges.contains(N))
					{
						M.fixNodeParameters(node, true);
//						M.fixGain(node, true);
//						M.fixLoss(node, true);
//						M.fixDuplication(node, true);
					}
				}
			} else // regular SPR move
			{
				assert (!og_right.isRoot());
				assert (!og_right.isDisownedChild());
				Node og_Rparent = og_right.getParent();
				int og_Rparent_arity = og_Rparent.getNumChildren();
				Node og_Rsibling = og_Rparent.getChild(
						og_right.getIndexAtParent()==0?1:0);
				
				Map<Node,Node> og_join_nodes = makeJoinPhylo();
				makeJoinModel(og_join_nodes);
				
				M = parameterOptimizer(rates, utable.mappedToTree(forest)); // new MLDistribution(rates, utable.mappedToTree(forest));
				
				Set<Node> optimized_edges = new HashSet<>();
				for (int node=0; node<forest.getNumNodes(); node++)
				{
					Node N = forest.getNode(node);
					if (og_join_nodes.containsKey(N))
					{
						Node ogN = og_join_nodes.get(N);
						if (ogN == og_left) // left edge at regraft
						{
							//assert (Double.isFinite(og_left.getLength())); // not really necessary 
							optimized_edges.add(N);
						}
						if (ogN == og_right) // right edge at regraft
						{
							//assert (Double.isFinite(og_right.getLength())); // not really necessary 
							optimized_edges.add(N);
						}
						if (ogN == og_Rsibling && og_Rparent_arity==2) // got replaced by Rsibling
						{
							optimized_edges.add(N);
							Node P = N.getParent();
							if (P!=null )
								optimized_edges.add(P);
						}
						if (ogN == og_Rparent)
						{
							// since prune position's parent had more than 2 children, it was kept 
							for (int ci=0; ci<N.getNumChildren(); ci++)
							{
								Node C = N.getChild(ci);
								if (Double.isFinite(rates.getEdgeLength(C.getIndex()))) // avoid optimizing all tree roots in the forest model
									optimized_edges.add(C);
							}
							optimized_edges.add(N); // even if infinite edge length
						}
					} else
					{
						// new node (regraft parent)
						 optimized_edges.add(N);
						
						// add all nodes in the subtree 
						//optimized_edges.addAll(N.listNodes(null, null));
					}
				}
				for (int node=0; node<forest.getNumNodes(); node++)
				{
					Node N = forest.getNode(node);
					if (!optimized_edges.contains(N))
					{
						M.fixNodeParameters(node, true);
//						M.fixGain(node, true);
//						M.fixLoss(node, true);
//						M.fixDuplication(node, true);
					}
				}
			}
			M.setCalculationWidth(calculation_width_absolute, calculation_width_relative);
			M.setMinimumObservedCopies(min_copies);
			double score = M.optimize(opt_eps, iter);
			for (int node=0; node<forest.getNumNodes(); node++)
			{
				double len = rates.getEdgeLength(node);
				forest.getNode(node).setLength(Double.max(EPS,len));
			}
			forest.hasLength(true);
			
			System.out.println("#**MLP.CM.sFML npar "+M.getModelParameterCount()+"\tscore "+score
						+"\t("+(score-best_model.score)+")"
						+"\t"+this);
			
			return score;
		}
		
		/**
		 * Parameter used for NNI exploration (unused)
		 * @return
		 */
		private double varML()
		{
			Gradient G  = new Gradient(rates, utable.mappedToTree(forest));
			double var = G.bootstrapLLVariance();
			return var;
		}
		
		private double scoreML(int iter)
		{
//			double fscore = 0.0;
//			TreeWithRates frates = null;
			Node P = null; // new node in SPR
			if (og_left != null && rates == null) 
			{
//				fscore = fastScoreML(iter);	
//				frates = this.rates;
				P = makeJoinModel(null);
			} 
			
			ML M = parameterOptimizer(rates, utable.mappedToTree(forest));
//			MLDistribution M = new MLDistribution(rates, utable.mappedToTree(forest));
//			M.setCalculationWidth(calculation_width_absolute, calculation_width_relative);
//			M.setMinimumObservedCopies(min_copies);
			
			double score = M.optimize(opt_eps, iter);
			for (int node=0; node<forest.getNumNodes(); node++)
			{
				double len = rates.getEdgeLength(node);
				forest.getNode(node).setLength(Double.max(EPS,len));
			}
			forest.hasLength(true);
			
//			if (fscore != 0.0)
//			{
//				double diff = fscore-score;
//				System.out.println("#**MLP.CM.sML "+this.id+"\tdiff "+diff+"\tscore "+score+"\tfast "+fscore+"\tP "+P);
//				if (Math.abs(diff)>30.0)
//				{
//					for (int node=0; node<forest.getNumNodes(); node++)
//					{
//						double r = rates.getGainParameter(node);
//						double p = rates.getLossParameter(node);
//						double q = rates.getDuplicationParameter(node);
//						
//						double fr = frates.getGainParameter(node);
//						double fp = frates.getLossParameter(node);
//						double fq = frates.getDuplicationParameter(node);
//						
//						double dr = (fr-r)/(r==0.0?1.0:r);
//						double dp = (fp-p)/p;
//						double dq = (fq-q)/(q==0.0?1.0:q);
//						
//						System.out.println("#**MLP.CM.sML "+node+"\t"+forest.getIdent(node)
//								+"\tdp "+dp+"\tdq "+dq+"\tdr "+dr
//								+"\tp "+fp+"/"+p
//								+"\tq "+fq+"/"+q
//								+"\tr "+fr+"/"+r);
//					}
//				}
//			}
			return score;
		}
		
		
		double scoreML()
		{
			long T0 = System.currentTimeMillis();
			double s = scoreML(opt_rounds);
			time_LL += System.currentTimeMillis()-T0;
			call_LL++;
			return s;
		}
		
		// Example: prune root U162 and graft as sibling of U105
		// #**MLP.sNNI 24/40       M11827.24@M0(N[U162 len Infinity prnt - chld {U82, U161}, nl 82, nn 163],N[U105 len 0.018661442940508362 prnt U161 chld {U97, U104}, nl 24, nn 47])[nodes=163, score=464427.85260794265]
		//         o U162                               * P
		//        / \                                  / \
		//       /   \                                /   \
		//      /     \                              /     \
		// U82 o       o U161               L=U162' o       o R=U105'
		//            / \                          / \
		//           /   \                        /   \
		//          /     \                      /     \
		//       s o       o U105          U82' o       o s
		
		// Example: prune root sibling U82 and graft as sibling of U105  
		//         o U162                               * U161'
		//        / \                                  / \
		//       /   \                                /   \
		//      /     \                              /     \
		// U82 o       o U161                     s o       o * P
 		//            / \                                  / \
		//           /   \                                /   \
		//          /     \                              /     \
		//       s o       o U105                L=U82' o       o R=U105'

		private Node makeJoinModel(Map<Node,Node> og_join_nodes)
		{
			if (og_join_nodes==null) og_join_nodes = makeJoinPhylo();
			Map<Node,double[]> model_rates = modelRates(og_model.rates);
			
			TreeWithRates base = new TreeWithRates(forest);
			Node L = null;
			Node R = null;
			for (int node=0; node<forest.getNumNodes(); node++)
			{
				Node N = forest.getNode(node);
				if (og_join_nodes.containsKey(N))
				{
					Node ogN = og_join_nodes.get(N);
					double[] rate_params = model_rates.get(ogN);
					base.setGainRate(node, rate_params[PARAMETER_GAIN]);
					base.setLossRate(node, rate_params[PARAMETER_LOSS]);
					base.setEdgeLength(node, rate_params[PARAMETER_LENGTH]);
					base.setDuplicationRate(node, rate_params[PARAMETER_DUPLICATION]);
					
					if (ogN == og_left) L = N;
					if (ogN == og_right) R = N;
					
					if (ogN.isRoot() && !N.isRoot()) // DEBUG
					{
						System.out.println("#**MLP.CM.mJM fromroot "+N+"\tfrom "+ogN+"\t"+base.toString(node)+"\t// "+this);
					}
					if (N.isRoot() && !ogN.isRoot()) 
					{
						base.setEdgeLength(node, Double.POSITIVE_INFINITY);
						base.setDuplicationRate(node, TreeWithRates.DEFAULT_DUPLICATION_RATE);
						
//						// DEBUG
//						System.out.println("#**MLP.CM.mJM toroot "+N+"\tfrom "+ogN+"\t"+base.toString(node)+"\t// "+this);
						
					}
				} else
				{
					// newly created graft parent node : keep params
//					if (N.isRoot()) // DEBUG
//					{
//						System.out.println("#**MLP.CM.mJM graftroot "+N+"\t"+base.toString(node)+"\t// "+this);
//					}
				}
			}
			
			// now the rates are copied, but they are problematic if the root changed:
			// (1) prune node was root (or forest root with infinite edge length)
			// (2) prune node was sibling 
			// (3) regraft target was root 
			
			
			assert (L!= null); // we found it
			assert (R != null); // we found it 
			Node P = L.getParent(); // graft node 
			assert (R.getParent() == P);
			
			int left = L.getIndex();
			int right = R.getIndex();
			int parent = P.getIndex();
			
			double lenL = base.getEdgeLength(left);
			if (lenL == Double.POSITIVE_INFINITY)
			{
				base.setEdgeLength(left, DEFAULT_EDGE_LENGTH);
				base.setEdgeLength(parent, Double.POSITIVE_INFINITY);
				base.setDuplicationRate(parent, TreeWithRates.DEFAULT_DUPLICATION_RATE);
//				// DEBUG
//				{
//					System.out.println("#**MLP.CM.mJM setleft+parent "+left+"\t"+parent
//							+"\t"+base.toString(left)+"\t"+base.toString(parent)+"\t// "+this);
//				}
			} else
			{
				// place P at half of L's edge; copy L'rates to P's rates 
				//lenL = lenL/2.0; // --- rather, dont create too short edges
				base.setEdgeLength(left, lenL);
				double[] rate_params = model_rates.get(og_left);

				base.setEdgeLength(parent, lenL);
				base.setGainRate(parent, rate_params[PARAMETER_GAIN]);
				base.setLossRate(parent, rate_params[PARAMETER_LOSS]);
				base.setDuplicationRate(parent, rate_params[PARAMETER_DUPLICATION]);
				
//				if (forest.isRoot(parent)) // DEBUG
//				{
//					System.out.println("#**MLP.CM.mJM setroot "+parent+"\t"+base.toString(parent)+"\t// "+this);
//				}
			} 
			double lenR = base.getEdgeLength(right);
			if (lenR == Double.POSITIVE_INFINITY)
			{
				lenR = DEFAULT_EDGE_LENGTH;
				base.setEdgeLength(right, lenR);
//				// DEBUG
//				{
//					System.out.println("#**MLP.CM.mJM setright "+right+"\t"+base.toString(right)+"\t// "+this);
//				}
			}

			this.rates = base;
			
//			// DEBUG
//			System.out.println("#**MLP.CM.mJM "+this+"\troot "+rates.toString(forest.getRoot())+"\tlenL "+lenL+"\tlenR "+lenR+"\tpar "+P+"\t"+rates.toString(parent));
			
			return P;
		}
		
		/**
		 * Sets {@link #forest} to be a copy of the {@link #og_model}'s 
		 * phylogeny, in which the SPR move is carried out.  
		 * 
		 * 
		 * @return Mapping from nodes in our {@link #forest} nodes to those in {@link #og_model}
		 */
		private Map<Node, Node> makeJoinPhylo()
		{
			Phylogeny joined_phylo = new Phylogeny(og_model.forest);
			
			Node L=null;
			Node R=null;
			this.og_nodes = new HashMap<>();
			for (int node=0; node<joined_phylo.getNumNodes(); node++)
			{
				Node ogN = og_model.forest.getNode(node);
				Node N = joined_phylo.getNode(node);
				og_nodes.put(N, ogN);
				if (ogN == og_left) L = N;
				if (ogN == og_right) R = N;
			}
			assert (L!=null);
			assert (R!=null);
			joined_phylo.pruneAndRegraft(R, L, false);
			this.forest = joined_phylo;
			this.version = og_model.version;
			this.rates = null;
			
			return og_nodes;
		}
		
		void acceptJoin()
		{
			if (isStale())
			{
				System.out.println("#**MLP.aJM version diff "+this+"\tog "+og_model);
			}

//			Map<Node,double[]> model_rates = modelRates(best_model.model.getBaseModel());

//			// join the two nodes 
			og_model.forest.pruneAndRegraft(og_right, og_left, false);
			Phylogeny joined_phylo = og_model.forest;
			Node P = og_left.getParent();
			assert (og_right.getParent() == P);

			// adjust the model for the phylogeny change
			TreeWithRates base = new TreeWithRates(joined_phylo);
			
			for (int node=0; node<joined_phylo.getNumNodes(); node++)
			{
				base.setGainRate(node, rates.getGainRate(node));
				base.setLossRate(node, rates.getLossRate(node));
				base.setDuplicationRate(node, rates.getDuplicationRate(node));
				base.setEdgeLength(node, rates.getEdgeLength(node));
				joined_phylo.getNode(node).setLength(this.forest.getLength(node));
			}
							
			og_model.rates = base;
			og_model.version++;
			//og_model.og_model = og_model;
		}
		
		/**
		 * Joins the right node's subtree as the sibling of the left node. 
		 * 
		 * @param left regraft node
		 * @param right prune node 
		 * @return
		 */
		CandidateMove join(int left, int right)
		{
			CandidateMove join = new CandidateMove();
			join.og_left = forest.getNode(left);
			join.og_right = forest.getNode(right);
			join.og_model = this;
			join.version = this.version;
			//makeJoin();
			
			return join;
		}
		
		
		List<CandidateMove> getNearestNeighborMoves()
		{
			List<CandidateMove> nni = new ArrayList<>();
			for (int rnode=0; rnode<forest.getNumNodes(); rnode++)
			{
				Node R = forest.getNode(rnode); // prune
				if (!R.isRoot())
				{
					List<Node> neighbors = R.graftNeighbors(3, 3);
					for (Node L: neighbors) // includes all NNI moves
					{
						int lnode = L.getIndex();
						nni.add(join(lnode,rnode));
					}
				}
			}
//			System.out.println("#**MLP.CM.gNNM "+nni.size());
			
			return nni;
			
		}
		
		/**
		 * Whether this is still a valid SPR move 
		 * in og_model: precisely if 
		 * regraft node (og_left) is <em>not</em> 
		 * in the prune node's (og_right) subtree.  
		 * 
		 * @return false if the SPR move would create a loop (prune becomes its own ancestor),
		 * 	
		 */
		boolean isValidJoin()
		{
			// regraft must not be in the prune's subtree 
			Node L = og_left;
			int num_nodes = og_model.forest.getNumNodes();  
			for (int __=0; __<num_nodes; __++) // used to avoid infinite loop caused by invalid phylogeny
			{
				if (L==og_right) return false; // ancestor 
				if (L.isRoot()) break;
				L = L.getParent();
			} 
//			boolean same_tree = (L==og_model.forest.getRootNode());
//			Node R = og_right; 
//			if (same_tree)
//			{
//				
//				for (int __=0; __<num_nodes && !R.isRoot(); __++) // used to avoid infinite loop caused by invalid phylogeny
//				{
//					R = R.getParent();
//				}
//				same_tree = (R == L);
//			}
//			if (!same_tree)
//			{
//				System.out.println("#**MLP.iVJ "+this+"\tdetached nodes "+og_left+"/"+L+"\t"+og_right+"/"+R+"\togroot "+og_model.forest.getRootNode());
//			}
//			
//			return same_tree;
			return true;
		}

		/**
		 * Places the root on a given edge. 
		 * 
		 * @param edge
		 * @return move with set {@link #rates} and {@link #forest}
		 */
		CandidateMove reroot(int edge)
		{
			Map<Node,double[]> model_rates = modelRates(rates);
			Phylogeny rotated = new Phylogeny(forest);
			Map<Node, double[]> og_rates = new HashMap<>();
			for (int node=0; node<rotated.getNumNodes(); node++)
			{
				Node ogN = forest.getNode(node);
				Node N = rotated.getNode(node);
				og_rates.put(N, model_rates.get(ogN));
			}
			Node ogR = rotated.getRootNode();
			Node R = rotated.getNode(edge);
			rotated.rerootEdge(R);
			
			TreeWithRates base = new TreeWithRates(rotated);
			for (int node=0; node<rotated.getNumNodes(); node++)
			{
				Node N = rotated.getNode(node);
				if(N==ogR)
				{
					base.setEdgeLength(ogR.getIndex(), DEFAULT_EDGE_LENGTH); // in case it was infinite 
				} else if (og_rates.containsKey(N))
				{
					double[] rate_params = og_rates.get(N);
					base.setGainRate(node, rate_params[PARAMETER_GAIN]);
					base.setLossRate(node, rate_params[PARAMETER_LOSS]);
					base.setDuplicationRate(node, rate_params[PARAMETER_DUPLICATION]);
					base.setEdgeLength(node, rate_params[PARAMETER_LENGTH]);
				} else 
				{
					// stay with default values
				}
			}
			//
			int rroot = rotated.getRoot();
			double[] rate_params = og_rates.get(ogR);
			base.setGainRate(rroot, rate_params[PARAMETER_GAIN]);
			base.setLossRate(rroot, rate_params[PARAMETER_LOSS]);
			base.setDuplicationRate(rroot, rate_params[PARAMETER_DUPLICATION]);
			base.setEdgeLength(rroot, rate_params[PARAMETER_LENGTH]);
			
			CandidateMove reroot = new CandidateMove(base);
			reroot.og_edge = forest.getNode(edge);
			reroot.og_model = this;
			return reroot;
		}
		
		/**
		 * Copies the contract move's parameters to
		 * {@link #og_model}.
		 */
		void acceptContract()
		{
			og_model.forest.fuseIntoParent(og_edge);
			og_model.rates = new TreeWithRates(og_model.forest);
			for (int node=0; node<forest.getNumNodes(); node++)
			{
				og_model.rates.setGainRate(node, rates.getGainRate(node));
				og_model.rates.setLossRate(node, rates.getLossRate(node));
				og_model.rates.setDuplicationRate(node, rates.getDuplicationRate(node));
				og_model.rates.setEdgeLength(node, rates.getEdgeLength(node));
				og_model.forest.getNode(node).setLength(forest.getLength(node));
			}
			og_model.version++;
		}
		
		/**
		 * Executes the contraction of the edge
		 * used at instantiation in {@link #contract(int)}.
		 * 
		 * @return the parent of the contracted edge
		 */
		private Node contract()
		{
			Map<Node,double[]> model_rates = modelRates(og_model.rates);
			forest = new Phylogeny(og_model.forest);

			og_nodes = new HashMap<>();
			Node E=null;
			for (int node=0; node<forest.getNumNodes(); node++)
			{
				Node N = forest.getNode(node);
				Node ogN = og_model.forest.getNode(node);
				if (ogN == og_edge)
				{
					E = N;
				} else
				{
					og_nodes.put(N, ogN);
				}
			}
			if (E==null)
				throw new java.lang.IllegalStateException("Edge does not exist anymore.");
			Node P = E.getParent();
			forest.fuseIntoParent(E);
			
			rates = new TreeWithRates(forest);
			for (int node=forest.getNumNodes(); 0<node; )
			{
				--node;
				Node N = forest.getNode(node);
				Node ogN = og_nodes.get(N);
				double[] og_rates = model_rates.get(ogN);
				rates.setGainRate(node, og_rates[PARAMETER_GAIN]);
				rates.setLossRate(node,  og_rates[PARAMETER_LOSS]);
				rates.setDuplicationRate(node, og_rates[PARAMETER_DUPLICATION]);
				rates.setEdgeLength(node, og_rates[PARAMETER_LENGTH]);
			}
			this.version = og_model.version;
			return P;
		}
		
		/**
		 * Contracts an edge in the phylogeny: child is fused int its parent
		 * so that  
		 * the grandchildren become the parent's children.
		 * 
		 * @param edge child node 
		 * @return corresponding move; call {@link #contract()} to initialize the model parameters before scoring
		 */
		CandidateMove contract(int edge)
		{
			CandidateMove contract = new CandidateMove();
			contract.og_model = this;
			contract.og_edge = forest.getNode(edge);
					
			contract.contract();
			return contract;
		}
		
		/**
		 * Compares two join moves : how many common nodes (in sets of {prune,regraft})
		 * @param other
		 * @return
		 */
		int sharedJoinNodes(CandidateMove other)
		{
			int shared  = 0;
			if (other.og_left == og_left || other.og_right == og_left)
				shared++;
			if (other.og_right == og_right || other.og_left == og_left)
				shared++;
			return shared;
		}
		
		
		@Override
		public int hashCode()
		{
			if (og_left==null) return super.hashCode();
			return 31*og_left.hashCode()+og_right.hashCode(); 
		}
		
		@Override
		public boolean equals(Object o)
		{
			if (o!=null && o instanceof CandidateMove)
			{
				if (og_left==null)
					return this==o;
				CandidateMove other = (CandidateMove) o;
				return og_model == other.og_model
						&& og_left == other.og_left
						&& og_right == other.og_right
						;
			} else
				return super.equals(o);
		}
		
		@Override
		public int compareTo(CandidateMove o) 
		{
			return Double.compare(this.score,o.score);
		}
				
		@Override
		public String toString()
		{
			StringBuilder sb = new StringBuilder();
			sb.append("M").append(id).append(".").append(version);
			if (og_model!=null)
			{
				sb.append("@M").append(og_model.id);
				if (og_edge != null)
				{
					sb.append(" * ").append(og_edge);
				} else if (og_left != null)
				{
					sb.append("(").append(og_left).append(",").append(og_right).append(")");
				} 
			}
			sb.append("[");
			if (forest != null)
				sb.append("nodes=").append(forest.getNumNodes()).append(", ");
			sb.append("score=").append(score).append("]");
			return sb.toString();
		}
		
	}
	
	
	private NearBestModels fromBundle(ModelBundle bundle)
	{
		ModelBundle.Entry R = bundle.getRoot();
		String threshold = R.getAttributeValue("max_worse_by");
		double max_worse_by = Double.POSITIVE_INFINITY;
		if (threshold == null)
		{
			System.out.println("#**MLP.fromBundle cannot find max_worse_by tag");
			// problem
		} else
		{
			max_worse_by = Double.parseDouble(threshold);
		}
		String type = R.getAttributeValue(CountXML.ATT_TYPE);
		if (!NearBestModels.class.getCanonicalName().equals(type))
		{
			System.out.println("#**MLP.fromBundle type mismatch "+type+"; expected "+NearBestModels.class.getCanonicalName());
		}
		NearBestModels NBM = new NearBestModels(max_worse_by);
		
		for (ModelBundle.Entry E: bundle.allEntries())
		{
			if (E.isRatesEntry())
			{
				GammaInvariant model = E.getModel();
				TreeWithRates rates = model.getBaseModel();
				CandidateMove M = new CandidateMove(rates);
				String score_str = E.getAttributeValue("score");
				if (score_str != null)
				{
					M.score = Double.parseDouble(score_str);
				} else
				{
					System.out.println("#*MLP.fromBundle computing score"+E.getId());
					M.score = M.score();
				}
				NBM.addMove(M);
			}
		}
		
//		
//		
//		for (GammaInvariant model: bundle.allModels())
//		{
//			NBM.addMove(M);
//		}
		double max_score_seen = NBM.top_models_seen.peek().score;
		double score_diff = max_score_seen - NBM.bestest_score;
		
		if (Double.isInfinite(max_worse_by))
		{
			max_worse_by = score_diff+0.0001;
			System.out.println("#*MLP.fromBundle guessing threshold "+max_worse_by);
			NBM.setThreshold(max_worse_by);
		}
		
		return NBM;
	}

	private class NearBestModels
	{
		private final Heap<CandidateMove> top_models_seen; // max-heap
		private double max_worse_by;
		private double bestest_score;
		private CandidateMove bestest_model;
		
		private final Map<CandidateMove, TreeComparator> trees;
		
		NearBestModels(){ this(Double.POSITIVE_INFINITY);}
		
		NearBestModels(double threshold)
		{
			this.top_models_seen = new Heap<>((a,b)->Double.compare(b.score, a.score)); // max-heap
			this.max_worse_by = threshold;
			this.bestest_score = Double.POSITIVE_INFINITY;
			this.bestest_model = null;
			trees = new HashMap<>();
		}
		
		void update(CandidateMove M)
		{
			System.out.println("#**MLP.NBM.uM "+M+"\tbestest "+bestest_model+"/"+bestest_score);

			top_models_seen.updateOrder(M);
			
			
			if (M.score<bestest_score)
			{
				newBestMove(M);
			}
		}
		
		int size()
		{
			return top_models_seen.size();
		}
		
		TreeComparator treeComparator(CandidateMove M)
		{
			return trees.get(M);
		}
		
		boolean contains(CandidateMove M)
		{
			return top_models_seen.contains(M);
		}
		
		private void newBestMove(CandidateMove M)
		{
			assert (M.score < bestest_score);
			// remove all by new threshold
			bestest_score = M.score;
			bestest_model = M;
			double smax = top_models_seen.peek().score;
			while (smax > bestest_score+max_worse_by)
			{
				
				CandidateMove bigM = top_models_seen.deleteLeast();
				trees.remove(bigM);
				System.out.println("#**MLP.NBM.nBM del "+bigM+"\tbest "+bestest_score);
				smax = top_models_seen.peek().score;
			}
			
		}
		
		/**
		 * Adds a model if it has a sufficiently good score.
		 * 
		 * @param M
		 * @return false if the move is not added 
		 */
		boolean addMove(CandidateMove M)
		{
			if (DEBUG_SCORE) MLPhylogeny.this.debugScore(M, System.out, ("NBM.aM\tbestsc "+bestest_score));

			if (M.score<=bestest_score+max_worse_by)
			{

				// should be added 
				IndexedTree Mtree = M.rates.getTree();
				for (CandidateMove prevM: trees.keySet())
				{
					TreeComparator cmp = trees.get(prevM);
					if (cmp.sameRootedTopology(Mtree))
					{
						System.out.println("#**MLP.NBM.aM "+M+"\trepeat "+prevM);

						if (prevM.score<=M.score)
						{
							// nah, we have seen better 
							M = null; // signal that its no good
						} else
						{
							trees.remove(prevM);
							top_models_seen.remove(prevM);
							System.out.println("#**MLP.NBM.nBM del "+prevM+"\tbetter "+M.score);
						}
						break;
					}
				}
				if (M!=null)
				{
					TreeComparator Mcmp = new TreeComparator(Mtree);
					trees.put( M, Mcmp);
					boolean new_on_heap = top_models_seen.add(M);


					if (M.score<bestest_score)
					{
						newBestMove(M);
					}
					System.out.println("#**MLP.NBM.aM add "+M+"\t(new "+new_on_heap+")"
							+"\tbest "+bestest_score+"\tn "+top_models_seen.size());
					return true;
				} // if M to be added 
			} // if score not too high
			System.out.println("#**MLP.NBM.aM "+M+"\tnoadd ");
			
			return false;
		}
		
		/**
		 * All stored moves in increasing order of score
		 * 
		 * @return
		 */
		CandidateMove[] allMoves()
		{
			CandidateMove[] all_moves = new CandidateMove[top_models_seen.size()];
			int i=0; 
			for (CandidateMove M: top_models_seen)
			{
				all_moves[i++] = M;
			}
			Arrays.sort(all_moves);
			return all_moves;
		}
		
		public ModelBundle toBundle(String bundle_id)
		{
			ModelBundle bundle = new ModelBundle(bundle_id);
			ModelBundle.Entry root = bundle.getRoot();
			root.setAttribute(CountXML.ATT_TYPE, this.getClass().getCanonicalName());
			root.setAttribute("max_worse_by", Double.toString(max_worse_by));
			
			CandidateMove[] all_moves = allMoves();
			for(int mi=0; mi<all_moves.length; mi++)
			{
				CandidateMove M = all_moves[mi];
				String tree_file = bundle_id+"."+Integer.toString(mi)+".tre";
				String rates_file = bundle_id+"."+Integer.toString(mi)+".rates.txt";
				
				DataFile<Phylogeny> tree_data = new DataFile<>(M.forest,new File((File)null, tree_file));
				DataFile<GammaInvariant> rates_data = new DataFile<>(new GammaInvariant(M.rates,1,1,1,1), new File((File)null, rates_file));
				ModelBundle.Entry T = root.addTree(tree_data);
				ModelBundle.Entry R =T.addRates(rates_data);
				R.setAttribute("score", Double.toString(M.score));
			}
			return bundle;
		}
		
//		TreeWithRates[] allModels()
//		{
//			CandidateMove[] all_moves = allMoves();
//			TreeWithRates[] all_models = new TreeWithRates[all_moves.length];
//			for (int m=0; m<all_models.length; m++)
//			{
//				all_models[m] = all_moves[m].rates;
//			}
//			return all_models;
//		}
		
		void setThreshold(double new_threshold)
		{
			if (new_threshold<max_worse_by)
			{
				double smax = top_models_seen.peek().score;
				while (smax>bestest_score + new_threshold)
				{
					CandidateMove bigM = top_models_seen.deleteLeast();
					trees.remove(bigM);
					smax =top_models_seen.peek().score;
				}
			} 
			max_worse_by = new_threshold;
		}
		
		/**
		 * sets new threshold but only if it is larger than the current one.
		 * @param new_threshold
		 */
		void increaseThreshold(double new_threshold)
		{
			if (new_threshold>max_worse_by)
				max_worse_by = new_threshold;
			// else untouched 
		}
		
		private Map<CandidateMove,double[]> annotateByBootstrap(int num_samples)
		{
			if (num_samples == 0) return null;
			
			CandidateMove[] moves = allMoves();
			RELL bootstrapper = new RELL(MLPhylogeny.this.table);
			
			Map<CandidateMove, double[]> annotate = new HashMap<>();
			
			boolean select_by_bootstrap;
//			do 
//			{
				bootstrapper.clear();
				for (CandidateMove C: moves)
				{
					bootstrapper.addRates(C.rates);
				}
			
				int[] model_wins = bootstrapper.getBootstrapWinners(MLPhylogeny.this.RND, num_samples);
			
				Integer[] LLorder = new Integer[model_wins.length];
				for (int i=0; i<LLorder.length; i++) LLorder[i]=i;
				Integer[] bootstrap_order = LLorder.clone();
				final double[] modelLL = bootstrapper.getCorrectedLL();
			
				Arrays.sort(LLorder,
						new java.util.Comparator<Integer>()
						{
							@Override
							public int compare(Integer a, Integer b)
							{
								return Double.compare(-modelLL[a],-modelLL[b]);
							}
						});
				Arrays.sort(bootstrap_order, (a,b)->Double.compare(model_wins[b],model_wins[a]));
			
				System.out.println("#**MPL.NBM.aBB boostrap winners: "+Arrays.toString(model_wins));
				System.out.println("#**MPL.NBM.aBB LL order : "+Arrays.toString(LLorder));
				System.out.println("#**MPL.NBM.aBB bootstrap order : "+Arrays.toString(bootstrap_order));
				
				for (int m=0; m<moves.length; m++)
				{
					int ci = bootstrap_order[m]; 
					System.out.println("#**MPL.NBM.aBB cand "+ci+"\tLL "+modelLL[ci]+"\twin "+model_wins[ci]);
				}
				
				select_by_bootstrap = false && (bootstrap_order[0]!=0);// do not reorder ever

//				if (select_by_bootstrap)
//				{
//					CandidateMove[] reordered = new CandidateMove[moves.length];
//					for (int m=0; m<reordered.length; m++)
//						reordered[m] = moves[bootstrap_order[m]];
//					moves = reordered;
//				} else
//				{
					double n = num_samples;
					for (int m=0; m<moves.length; m++)
					{
						CandidateMove M = moves[m];
	
						TreeComparator cmp = trees.get(M);
						int[] support = bootstrapper.getBootstrapSupport(model_wins, cmp);
						double[] rel_support = new double[support.length]; 
						for (int node=0; node<support.length; node++)
						{
							rel_support[node] = support[node]/n;
							if (m==0)
							{
								System.out.println("#BOOTSTRAP\t"+node+"\t"+support[node]+"\t"+rel_support[node]+"\t// "+cmp.getReferenceTree().toString(node));
							}
						}
						annotate.put(M, rel_support);
					}
//				}
//			} while (select_by_bootstrap);
			
			best_model = moves[0];
			return annotate;
		}
			
		
		void reportTrees(PrintStream out, int num_bootstrap_samples)
		{
			Map<CandidateMove, double[]> node_support = annotateByBootstrap(num_bootstrap_samples);
			for (CandidateMove M: allMoves())
			{
				String[] bootstrap_annotations=null;
				if (node_support != null)
				{
					double[] support = node_support.get(M);
					if (support != null)
					{
						assert (support.length ==  M.forest.getNumNodes());
						bootstrap_annotations = new String[support.length];
						for (int node=M.forest.getNumLeaves(); node<support.length; node++)
						{
							bootstrap_annotations[node] = Double.toString(support[node]);
						}
					}
				}
//				IndexedTree Mtree = M.forest;
				String tree_description = NewickParser.printTree(M.forest, bootstrap_annotations);
				double ascore = M.score/table.getFamilyCount();
				out.println(tree_description+"\t[SCORE "+M.score+"\tAVGSCORE "+ascore+"]");
			}
		}
		
		
		void reportAllModels(String tree_file, String rates_file, int num_bootstrap_samples) throws java.io.IOException
		{
			boolean bootstrap_in_name = false;
			CandidateMove[] moves = allMoves();
			int[] model_wins = null;
			String[][] node_support = new String[moves.length][];
			
			if (num_bootstrap_samples!=0)
			{
				RELL bootstrapper = new RELL(MLPhylogeny.this.table);
				for (CandidateMove C: moves)
				{
					bootstrapper.addRates(C.rates);
				}
				model_wins = bootstrapper.getBootstrapWinners(MLPhylogeny.this.RND, num_bootstrap_samples);
				System.out.println("#**MPL.NBM.rAM boostrap winners: "+Arrays.toString(model_wins));				
				
				for(int mi=0; mi<moves.length; mi++)
				{
					CandidateMove M = moves[mi];
					TreeComparator cmp = trees.get(M);					
					int[] support = bootstrapper.getBootstrapSupport(model_wins, cmp);
					Phylogeny Mtree = M.forest;
					
					int num_nodes = M.forest.getNumNodes();
					String[] bootstrap_info = (bootstrap_in_name?null:new String[num_nodes]);
					for (int node=M.forest.getNumLeaves(); node<num_nodes; node++)
					{
						if (!Mtree.isRoot(node))
						{
							String node_info = Double.toString(support[node]/(double)num_bootstrap_samples);
							if (bootstrap_in_name)
							{
								Mtree.getNode(node).setName(node_info);
							} else
							{
								bootstrap_info[node] = node_info;
							}
						}
					}
					node_support[mi] = bootstrap_info;
				}
			}
			for(int mi=0; mi<moves.length; mi++)
			{
				CandidateMove M = moves[mi];
				Phylogeny Mtree = M.forest;
				String model_suffix = "."+Integer.toString(mi);		
				String out_tree_file = tree_file + model_suffix;
				PrintStream tree_out = new PrintStream(out_tree_file);
				String tree_description = NewickParser.printTree(Mtree, node_support[mi]);
				tree_out.println(tree_description);
				
//				for (int node=M.forest.getNumLeaves(); node<num_nodes; node++)
//				{
//					TreeComparator cmp = trees.get(M);					
//					int[] support = bootstrapper.getBootstrapSupport(model_wins, cmp);
//					out.println("#BOOTSTRAP\t"+node+"\t"+support[node]+"\t"+rel_support+"\t// "+cmp.getReferenceTree().toString(node));
				
				
				tree_out.close();
				
				TreeWithRates Mrates = M.rates;
				String out_rates_file = rates_file + model_suffix;
	    		PrintStream rates_out = new PrintStream(out_rates_file);
	    	    rates_out.println(CommandLine.getStandardHeader(MLPhylogeny.class));	
	    	    rates_out.println("#TREE "+tree_description);
	    	    
	    	    TreeComparator cmp = trees.get(M);	
	    	    StringBuilder tree_info = null;
	    	    for (int mj=0; mj<moves.length; mj++)
	    	    {
	    	    	if (mi != mj)
	    	    	{
	    	    		Phylogeny tree2 = moves[mj].forest;
	    	    		TreeComparator.NodeMap map = cmp.map(tree2);
	    	    		boolean same_topo = true;
	    	    		int[] ridx = map.toReference();
	    	    		for (int node: ridx)
	    	    		{
	    	    			if (node==-1)
	    	    			{
	    	    				same_topo = false;
	    	    				break;
	    	    			}
	    	    		}
	    	    		if (same_topo)
	    	    		{
	    	    			if (tree_info==null)
	    	    			{
	    	    				tree_info = new StringBuilder("# Same unrooted topology as ");
	    	    			} else
	    	    			{
	    	    				tree_info.append(", ");
	    	    			}
	    	    			tree_info.append("tree#").append(mj);
	    	    			tree_info.append("(root ").append(Mtree.getIdent(ridx[tree2.getRoot()])).append(")");
	    	    		}
	    	    	}
	    	    }
	    	    if (tree_info != null)
	    	    {
	    	    	rates_out.println(tree_info.toString());
	    	    }
	    	    
	    	    
	    	    if (num_bootstrap_samples!=0)
	    	    {
	    	    	double model_support = model_wins[mi]/(double)num_bootstrap_samples;
	    	    	rates_out.println("#BOOTSTRAP best\t"+model_support+"\t"+model_wins[mi]+"\t/"+num_bootstrap_samples);
	    	    }
	    	    rates_out.println(RateVariationParser.printRates(Mrates));

    	    	if (DEBUG_SCORE) MLPhylogeny.this.debugScore(M, rates_out, "MLP.NBM.rAM");
	    	    
	    		double ascore = M.score/table.getFamilyCount();
	        	rates_out.println("#SCORE "+M.score);
	        	rates_out.println("#AVGSCORE "+ascore);
	        	
	        	rates_out.close();
			}
		}
	}
	
	
	private void saveNearBestModels(String bundle_file, String bundle_id) throws java.io.IOException
	{
		ModelBundle bundle = near_best_models.toBundle(bundle_id);
		PrintStream out = new PrintStream(bundle_file);
		
		ModelBundle.printBundle(out, bundle);
		
//		String bundle_xml = bundle.toXMLString();
//		out.println(bundle_xml);
		out.close();
	}
	
	private void loadNearBestModels(String bundle_file) throws java.io.IOException, org.xml.sax.SAXException, javax.xml.parsers.ParserConfigurationException 
	{
		List<ModelBundle> input_bundles = ModelBundle.readBundle(GeneralizedFileReader.guessBufferedReaderForInput(bundle_file));
		if (input_bundles.size()>1)
		{
			System.out.println("#**MLP.loadNBM multiple ("+input_bundles.size()+") bundles in file; using the first");
		}
		ModelBundle bundle = input_bundles.get(0);
		this.near_best_models = fromBundle(bundle);
		this.best_model = near_best_models.bestest_model; // with score set 
		if (fit_parsimony)
		{
			fitParsimonyPenalties(best_model.rates, utable.mappedToTree(best_model.forest));
		}
	}
	
	private void debugScore(CandidateMove M, PrintStream out, String reason)
	{
		long T0 = System.currentTimeMillis();
		Gradient G = new Gradient(M.rates, utable.mappedToTree(M.forest));
		G.setMinimumObservedCopies(min_copies);
		G.setCalculationWidthThresholds(12, 3.0);
		double gscore = -G.getCorrectedLL();
		double gLL = G.getLL();
		double g0 = G.getUnobservedLL();
		
		long dT = System.currentTimeMillis()-T0;
		this.time_gradient += dT;
		this.call_gradient ++;
		
		double dg = (M.score-gscore);
		out.println("#**MLP.debugS "+reason+"\t"+M+"\tG "+gscore+"\tuncorr "+gLL+"\tgunobs "+g0
					+"\t// diff "+dg+"\trdiff "+dg/gscore+"\ttiming "+dT+" ms");
//		if (Math.abs(dg)>1.0)
//		{
//			out.println("#CHECK phylo\t"+NewickParser.printTree(M.forest));
//			out.println(RateVariationParser.printRates(M.rates));
//		}
	}
	
	private void debugScore(PrintStream out)
	{
		this.debugScore(best_model, out, "");
//		Gradient G = new Gradient(best_model.rates, utable.mappedToTree(best_model.forest));
//		G.setMinimumObservedCopies(min_copies);
//		G.setCalculationWidthThresholds(12, 3.0);
//		double gscore = -G.getCorrectedLL();
//		double gLL = G.getLL();
//		double g0 = G.getUnobservedLL();
//		double dg = (best_model.score-gscore);
//		out.println("#**MLP.checkS "+"\t"+best_model+"\tG "+gscore+"\tuncorr "+gLL+"\tgunobs "+g0
//					+"\t// diff "+dg+"\trdiff "+dg/gscore);
	}
	
	public static void main(String[] args) throws Exception
	{
//		PRINT_OPTIMIZATION_MESSAGES = true;
		
		PrintStream out = System.out;
	    out.println(CommandLine.getStandardHeader(MLPhylogeny.class));
	    out.println(CommandLine.getStandardRuntimeInfo(MLPhylogeny.class, args));

	    CommandLine cli = new CommandLine(args,MLPhylogeny.class, 0) ;
		AnnotatedTable table = cli.getTable();
		if (table==null)
		{
			throw new IllegalArgumentException("Specify table either with "+CommandLine.OPT_TABLE+" or "
					+ "as 2nd argument with empty tree: "+CommandLine.NO_FILE+" table");
		}
		
		MLPhylogeny P = new MLPhylogeny(table);
        if (cli.getOptionValue(OPT_TRUNCATE)!=null)
        {
        	int absolute = cli.getOptionTruncateAbsolute();
        	double relative = cli.getOptionTruncateRelative();
        	P.setCalculationWidth(absolute, relative);
        } 
        out.println(CommandLine.getStandardHeader("Truncated computation: -"
        		+OPT_TRUNCATE+" "+P.calculation_width_absolute+","+P.calculation_width_relative));
        int opt_rounds = cli.getOptionInt(OPT_ROUNDS, DEFAULT_ROUNDS);
        double opt_eps = cli.getOptionDouble(OPT_EPS, DEFAULT_CONVERGENCE);
        P.setOptimizationConvergence(opt_rounds, opt_eps);
        out.println(CommandLine.getStandardHeader("Convergence: -"+OPT_ROUNDS+" "+P.opt_rounds+" -"+OPT_EPS+" "+P.opt_eps));
		if (cli.getOptionValue(OPT_MINCOPY)!=null)
		{
			int min_copies = Integer.parseInt(cli.getOptionValue(OPT_MINCOPY));
			P.setMinCopies(min_copies);
		}
		out.println(CommandLine.getStandardHeader("Minimum observed copies: -"+OPT_MINCOPY+" "+P.min_copies));
		int top_moves = cli.getOptionInt(OPT_TOP, TOP_MOVES);
		P.setTopScoredMoves(top_moves);
		out.println(CommandLine.getStandardHeader("Top scored moves: -"+OPT_TOP+" "+P.top_scored_moves));
		
		Random RND = cli.getOptionRND(out);
		P.setRandom(RND);

		if (cli.getOptionValue("opt.rates")!=null )
		{
			boolean want_rates = cli.getOptionBoolean("opt.rates", false);
			P.optimize_rates = want_rates;
			out.println(CommandLine.getStandardHeader("ML optimization: -opt.rates "+want_rates+"\t("+(want_rates?MLRates.class:MLDistribution.class)+")"));
			out.println(CommandLine.getStandardHeader("Heuristic scoring by "+(P.heuristic_byParsimony?"parsimony":"fastML")));
		}
		
//		int gamma_length = DEFAULT_GAMMA_LENGTH;
//		if (cli.getOptionValue(OPT_MODEL_GAMMA_CAT)!=null)
//		{
//			gamma_length = Integer.parseInt(cli.getOptionValue(OPT_MODEL_GAMMA_CAT));
//		}
//		P.setGammaCategories(gamma_length);
				
		if (cli.getOptionValue("opt.EM")!=null)
		{
			boolean wantEM = cli.getOptionBoolean("opt.EM", false);
			// wantEM = wantEM && P.best_model.getNumTrees()==0;
			
			P.optimize_EM = wantEM;
			{
				out.println(CommandLine.getStandardHeader("ML optimization: -opt.EM "+wantEM+"\t("+(DirectEM.class)+")"));
			}
//			int maxiterEM = 256;
//			if (wantEM && P.opt_rounds > maxiterEM)
//			{
//				P.setOptimizationConvergence(maxiterEM, opt_eps);
//				out.println(CommandLine.getStandardHeader("Convergence reset: -"+OPT_ROUNDS+" "+maxiterEM));
//			}
		} 
		
		if (cli.getOptionValue(OPT_PARSIMONY_FIT)!=null)
		{
			boolean parsimony_fit = cli.getOptionBoolean(OPT_PARSIMONY_FIT, P.fit_parsimony);
			P.fit_parsimony=parsimony_fit;
			out.println(CommandLine.getStandardHeader("Parsimony penalty fitting: -"+OPT_PARSIMONY_FIT+" "+parsimony_fit));
		}
		
		if (cli.getOptionValue(OPT_LOAD)!=null)
		{
			String bundle_file = cli.getOptionValue(OPT_LOAD);
			out.println("#**MLP.main loading bundle: -"+OPT_LOAD+" "+bundle_file);
			P.loadNearBestModels(bundle_file);
		} else
		{
			// 
			if (cli.getTree()==null)
			{
				if (OPT_RND.equals(cli.getOptionValue(OPT_BUILD)))
				{
					out.println("#**MLP.main random starting tree: -"+OPT_BUILD+" "+cli.getOptionValue(OPT_BUILD));
					P.randomModel();
				} else if ("upgma".equalsIgnoreCase(cli.getOptionValue(OPT_BUILD)))
				{
					out.println("#**MLP.main UPGMA starting tree: -"+OPT_BUILD+" "+cli.getOptionValue(OPT_BUILD));
					P.upgmaModel();
				} else
				{
					double build_threshold = cli.getOptionDouble(OPT_BUILD, Double.POSITIVE_INFINITY);
					out.println("#**MLP.main build: -"+OPT_BUILD+" "+build_threshold);
					P.buildModel(build_threshold);
				}
			} else
			{
				if (cli.getModel()==null)
				{
					TreeWithRates default_rates = new TreeWithRates(cli.getTree());
					P.initModel(default_rates);
				} else
				{
					GammaInvariant rates_model = cli.getModel();
					if (rates_model.getNumActiveClasses()>1)
					{
						out.println("#**MLP.main ignore rate variation with "+rates_model.getNumActiveClasses()+" classes; continuing with base rates.");
					}
					P.initModel(cli.getModel().getBaseModel());
				}
			}
		}
		
		
		if (P.best_model.getNumTrees()==0)
		{
			double root_threshold = cli.getOptionDouble(OPT_REROOT, Double.NEGATIVE_INFINITY); // don't want it by default
			if (root_threshold != Double.NEGATIVE_INFINITY)
			{
				out.println("#**MLP.main reroot: -"+OPT_REROOT+" "+root_threshold);
				if (top_moves > P.best_model.forest.getNumNodes())
				{
					top_moves = P.best_model.forest.getNumNodes();
					out.println("#**MLP.main top reset to max : -"+OPT_TOP+" "+top_moves);
					P.setTopScoredMoves(top_moves);
				}
				P.placeRoot(root_threshold);
			}
		}
		
		double pollard_threshold = cli.getOptionDouble(OPT_POLLARD, 1.0);
		out.println("#**MLP.main pollard: -"+OPT_POLLARD+" "+pollard_threshold);
		if (pollard_threshold != 1.0)
		{
			P.pollardEdges(pollard_threshold);
		}
		
		String save_bundle = null;
		if (cli.getOptionValue(OPT_SAVE)!=null)
		{
			save_bundle = cli.getOptionValue(OPT_SAVE);
			if (save_bundle.endsWith(".xml"))
			{
				save_bundle = save_bundle.substring(0, save_bundle.length()-4);
			}
			out.println("#**MLP.main bundle name: -"+OPT_SAVE+" "+save_bundle);
		}
		boolean want_snapshots= cli.getOptionBoolean(OPT_SNAPSHOT, false);
		if (save_bundle != null)
		{
			out.println("#**MLP.main snapshots: -"+OPT_SNAPSHOT+" "+want_snapshots);
			if (want_snapshots)
			{
				P.SetSession(save_bundle);
			}
		} else if (want_snapshots)
		{
			out.println("#**MLP.main snapshots: -"+OPT_SNAPSHOT+" "+want_snapshots+" ignored, because no bundle name set with -"+OPT_SAVE);
			want_snapshots = false;
		}

		double walk_threshold = cli.getOptionDouble(OPT_SPR   , Double.NEGATIVE_INFINITY); // don't want it by default
		out.println("#**MLP.main walk: -"+OPT_SPR+" "+walk_threshold);
		if (walk_threshold != Double.NEGATIVE_INFINITY)
		{
			if (cli.getOptionValue(OPT_TOUR)==null 
					&& cli.getOptionValue(OPT_SEARCH)==null
					&& cli.getOptionValue(OPT_PLACE)==null)
			{
				System.out.println("#**MLP.main walkModel: -"+OPT_SPR+" "+walk_threshold);
				P.walkModel(walk_threshold);
//				if (want_snapshots)
//				{
//					// save snapshot of state
//					String bundle_file = P.session_id+"w.xml";
//					try
//					{
//						P.saveNearBestModels(bundle_file, P.session_id);
//					} catch (Exception E)
//					{
//						System.out.println("#**MLP.main walk snapshot failed "+E);
//					}
//				}
			}  
			
			if (cli.getOptionValue(OPT_PLACE)!=null)
			{
				String comma_list =  cli.getOptionValue(OPT_PLACE);
				Set<String> wanted_leaf_names = new HashSet<>();
				for (String s: comma_list.split(","))
				{
					wanted_leaf_names.add(s);
				}
				out.println("#**MLP.main placing nodes: -"+OPT_PLACE+" "+comma_list);
				P.placeNodes(walk_threshold, wanted_leaf_names);
			}
			if (cli.getOptionValue(OPT_TOUR)!=null)
			{
				int num_walks = cli.getOptionInt(OPT_TOUR, 1);
				out.println("#**MLP.main tour: -"+OPT_TOUR+" "+num_walks);
				double top_increase = Math.sqrt(2.0);
				double convergence_decrease = 0.25;
				P.tourModel(num_walks, walk_threshold, top_increase, convergence_decrease);
			}
			if (cli.getOptionValue(OPT_SEARCH)!=null)
			{
				int search_steps = cli.getOptionInt(OPT_SEARCH, 0);
				out.println("#**MLP.main search: -"+OPT_SEARCH+" "+search_steps);
				P.searchModels(search_steps, walk_threshold);
			}
			
		}
		
		double contract_threshold = cli.getOptionDouble(OPT_CONTRACT, Double.NEGATIVE_INFINITY); 
		out.println("#**MLP.main contract: -"+OPT_CONTRACT+" "+contract_threshold);
		if (contract_threshold != Double.NEGATIVE_INFINITY)
		{
			P.contractEdges(contract_threshold);
		}

		double dsec = cli.getMillisSinceStart()/1000.0;
		
		// save results
		if (save_bundle!=null)
		{
			String bundle_file = save_bundle+".xml";
			P.saveNearBestModels(bundle_file, save_bundle);
			out.println("#**MLP.main model bundle written to "+bundle_file);			
		}
		
		
    	String out_tree_file = cli.getOptionValue(OPT_OUTPUT_TREE);
    	String out_rates_file = cli.getOptionValue(OPT_OUTPUT);
    	
    	String tree_description = NewickParser.printTree(P.getPhylogeny());

    	
    	
//		boolean saveall = (out_tree_file != null) && !"-".equals(out_tree_file)
//						&& (out_rates_file != null) && !"-".equals(out_rates_file);
//		
//		saveall = saveall && cli.getOptionBoolean(OPT_OUTPUT_SAVEALL, saveall);
//		System.out.println("#**MLP.main saveall: -"+OPT_OUTPUT_SAVEALL+" "+saveall);
//		
//		
//		
//		if (saveall)
//		{
//			
//			int num_smp = cli.getOptionInt(OPT_BOOTSTRAP, 0);
//			if (num_smp !=0 && P.near_best_models.size()<2)
//			{
//				System.out.println("#**MPL.main cannot boostrap with 1 tree");
//				num_smp = 0;
//			} 
//			P.near_best_models.reportAllModels(out_tree_file, out_rates_file, num_smp);
//		} else  // still do bootstrap, if requested
//		{
    		int num_bootstrap_samples = 0;
			if (cli.getOptionValue(OPT_BOOTSTRAP)!=null)
			{
				num_bootstrap_samples = cli.getOptionInt(OPT_BOOTSTRAP, 1024);
				out.println(CommandLine.getStandardHeader("Boostrap/RELL sampling: -"+OPT_BOOTSTRAP+" "+num_bootstrap_samples));
				if (P.near_best_models.size()<2)
				{
					out.println("#**MPL.main cannot bootstrap with 1 tree");
					num_bootstrap_samples = 0;
				} else
				{
					// we will do bootstrap when writing out the trees	
				}
			}
//		}
		
		double ascore = P.getScore()/table.getFamilyCount();
		
    	out.println("#SCORE "+P.getScore());
    	out.println("#AVGSCORE "+ascore);
    	
    	if (DEBUG_SCORE) P.debugScore(System.out);

    	GammaInvariant model = P.getRateVariation();
    	
    	PrintStream tree_out = null;
    	
    	if (out_tree_file == null && num_bootstrap_samples != 0)
    	{
    		if (save_bundle!=null)
    			out_tree_file = save_bundle+".xml.trees";
    		else
    			out_tree_file = "-";
    		out.println("#**MLP.main bootstrap requested but no -"+OPT_OUTPUT_TREE+"; so tree file name set to "+out_tree_file);
    	}
    	
    	if (out_tree_file != null)
    	{
    		if (out_tree_file.equals("-"))
    		{
    			tree_out = System.out;
    		} else
    		{
    			tree_out = new PrintStream(out_tree_file);
    		}
//    		tree_out.println(tree_description);
    		P.near_best_models.reportTrees(tree_out, num_bootstrap_samples); 
//    		if (tree_out != System.out)
//    			tree_out.close();
    		tree_out.close();
    		out.println("#**MLP.main best tree(s) written to "+out_tree_file);    		
    	} 
    	
    	PrintStream rates_out = out;
    	if (out_rates_file!=null)
    	{
    		rates_out = new PrintStream(out_rates_file);
    	    rates_out.println(CommandLine.getStandardHeader(MLPhylogeny.class));
    	    rates_out.println(CommandLine.getStandardRuntimeInfo(MLPhylogeny.class, args));
    	}
    	rates_out.println("#TREE "+tree_description);
    	rates_out.println(RateVariationParser.printRates(model));
		if (rates_out != System.out)
		{
    		out.println("#**MLP.main best rates written to "+out_rates_file);
			rates_out.close();
		}
		
		
		double t_eval = ((double)P.time_LL)/P.call_LL;
		
		double t_feval = ((double)P.time_fLL)/(1L+P.call_fLL);
		
		
		out.printf("#**MLP.main done in %.3f seconds; %d MLcalls %.3f msec/call total %.3f s\t%d MPcalls\t%d fastML %.3f msec/call total %.3f s"
					,dsec,P.call_LL,t_eval,P.time_LL/1000.0, P.call_MP,
					P.call_fLL, t_feval, P.time_fLL/1000.0);

		double t_gradient = ((double)P.time_gradient)/P.call_gradient;
		out.printf(";\t%d Gradient %.3f msec/call total %.3f s",P.call_gradient,t_gradient,P.time_gradient/1000.0);

		out.println();
	}

}
