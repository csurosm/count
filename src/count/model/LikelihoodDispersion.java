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

import static count.io.CommandLine.OPT_MINCOPY;
import static count.io.CommandLine.OPT_ROUNDS;
import static count.io.CommandLine.OPT_ROWCOUNT;

import java.io.PrintStream;
import java.util.Random;

import count.io.CommandLine;
import count.ds.Phylogeny;

/**
 * Experimental code used during development.
 * 
 * @author csuros
 *
 */
public class LikelihoodDispersion 
{
	public LikelihoodDispersion(TreeWithRates rates, Random RND)
	{
		this.rates = rates;
		
		this.sim = new SimulatedEvolution(new GammaInvariant(rates,1,1,1,1), RND);
	}
	
	private final SimulatedEvolution sim;
	private final TreeWithRates rates;

	void collectStatistics(int num_rows, int min_observed, int repeats)
	{
		this.collectStatistics(num_rows, min_observed, repeats, this.rates);
	}

	
	void collectStatistics(int num_rows, int min_observed, int repeats, TreeWithRates rates)
	{
		double[] ascore = new double[repeats];
		
		double sum_a = 0.0;
		double sum_a2 = 0.0;
		
		double sum_oria = 0.0;
		
		int ori_wins=0;
		for (int rep=0; rep<repeats; rep++)
		{
			SimulatedEvolution.Table table = sim.table(num_rows, min_observed);
			table.fillTable();
			Gradient G = new Gradient(rates,table.mappedToTree(rates.getTree()));
			G.setMinimumObservedCopies(min_observed);
			G.computeParameters();
			double LL = G.getCorrectedLL();
			double bvar = G.bootstrapLLVariance()/num_rows;
			double bsd = Math.sqrt(bvar);
			double smp_bsd = Math.sqrt(bvar/num_rows); 
			double avg_score = -LL/num_rows;
			Gradient oriG = new Gradient(this.rates, table);
			oriG.setMinimumObservedCopies(min_observed);
			double oriLL = oriG.getCorrectedLL();
			
			boolean ori_better = LL<oriLL;
			if (ori_better) ori_wins++; 

			if (!ori_better || rep<100 || rep % 1000==0)
			{ 
				System.out.println("#**LD.cS "+rep+"\t"+avg_score+"\t"+bvar+"\t"+bsd+"\t"+smp_bsd
						+"\t"+(-LL)+"\t"+(-oriLL)
						+"\t"+(LL<oriLL?"+1":"-1")+"\t"+ori_wins);
			}
			ascore[rep] = avg_score;
			sum_a += avg_score;
			sum_a2 += avg_score*avg_score;
			
			sum_oria += -oriLL/num_rows;
		}
		
		double avg_a = sum_a/repeats;
		double avg_a2 = sum_a2/repeats;
		
		double avg_oria = sum_oria/repeats;
		// 
		double var_a = avg_a2 - avg_a*avg_a;
		double sd_a = Math.sqrt(var_a);
		
		double p_ori_err = 1.0-ori_wins/(repeats+0.);
		System.out.println("#SCOREDIST\t"+num_rows+"\t"+repeats+"\t"+avg_a+"\t"+sd_a+"\t"+ori_wins+"\t"+p_ori_err+"\t"+avg_oria);//+"\taa2 "+avg_a2+"\tv "+var_a);
	}
	
	
	public static void main(String[] args) throws Exception
	{
		PrintStream out = System.out;
		CommandLine cli = new CommandLine(args, LikelihoodDispersion.class, 1);
		if (cli.getModel() == null)
			throw new IllegalArgumentException("Specify the rates model");
		int num_rows;
		if (cli.getTable() == null)
		{
			num_rows = 10;
		} else
		{
			num_rows = cli.getTable().getFamilyCount();
		}
		num_rows = cli.getOptionInt(OPT_ROWCOUNT, num_rows);
		int min_obs;
		if (cli.getTable() == null)
		{
			min_obs = 1;
		} else
		{
			min_obs = Integer.min(2,cli.getTable().minCopies());
		}
		min_obs = cli.getOptionInt(OPT_MINCOPY, min_obs);
		
		out.println(CommandLine.getStandardHeader("Families: -"+OPT_ROWCOUNT+" "+num_rows));
		out.println(CommandLine.getStandardHeader("Minimum observed: -"+OPT_MINCOPY+" "+min_obs));
		GammaInvariant input_model = cli.getModel();
//		if (cli.getFreeModel()!=null)
//		{
//			input_model= cli.getFreeModel();
//			out.println("#SE.main: FreeModel ("+input_model.getNumClasses()+" classes)");
//		}
		
		Random RND = cli.getOptionRND(out);
		LikelihoodDispersion LD = new LikelihoodDispersion(input_model.getBaseModel(), RND);

		
		
		int maxiter = cli.getOptionInt(OPT_ROUNDS, 100);
		if (2<=cli.getExtraArgumentCount())
		{
			String tree_file2 = cli.getExtraArgument(0);
			Phylogeny tree2 = cli.getTreeFromArgument(tree_file2).getContent();
			String model_file2 = cli.getExtraArgument(1);
			out.println(CommandLine.getStandardHeader("Alt. tree file: "+tree_file2));
			out.println(CommandLine.getStandardHeader("Alt. rates file: "+model_file2));
			GammaInvariant model2 = cli.getModelFromArgument(tree2, model_file2).getContent();
			LD.collectStatistics(num_rows, min_obs, maxiter, model2.getBaseModel());
		} else
		{
			LD.collectStatistics(num_rows, min_obs, maxiter);
		}
	}
}
