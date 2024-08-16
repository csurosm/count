package count.model;
/*
 * Copyright 2021 Mikl&oacute;s Cs&#369;r&ouml;s.
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

import count.ds.AnnotatedTable;
import count.ds.ProfileTable;
import count.io.CommandLine;
import count.matek.Logarithms;


/**
 * Standalone class for testing Ancestor calculations during early code development. 
 *
 * 
 * @author csuros
 * @deprecated
 */
public class AncestorPosteriors
{
	public AncestorPosteriors(TreeWithRates rates, ProfileTable table)
	{
		this(new Posteriors(rates, table));
	}
	
	public AncestorPosteriors(Posteriors post)
	{
		this.post = post;
		this.ancestor_posteriors = new Posteriors[post.factory.tree.getNumNodes()];
		
		initDataStructures();
	}
	
	private final Posteriors post;
	private final Posteriors[] ancestor_posteriors;
	
	private int ancestor_width = 2;
	
	private int threshold_width_absolute = 3; //Integer.MAX_VALUE;
	private double threshold_width_relative = 1.0; //Double.POSITIVE_INFINITY;
	/**
	 * Sets up truncated likelihood calculations: 
	 * absolute + relative * sqrt(max+1); calculated in {@link #getCalculationWidth(int)}
	 * 
	 * @param absolute
	 * @param relative
	 */
	public void setCalculationWidthThresholds(int absolute, double relative)
	{
		this.threshold_width_absolute = absolute;
		this.threshold_width_relative = relative;
		
//		int m = getCalculationWidth(max_family_size[tree.getRoot()]);
//		factorials = new RisingFactorial(1.0, 1+m); 		
	}	
	
	/**
	 * Maximum ancestor copy numbers 
	 * given maximum observed value at leaves. 
	 * 
	 * @param max_value
	 * @return
	 */
	private int getCalculationWidth(int max_value)
	{
		double rel = Math.ceil(threshold_width_relative*Math.sqrt(max_value+1.0));
		double width =  max_value+rel+threshold_width_absolute;
		int cwidth = (int) Double.min(width,Integer.MAX_VALUE);
		return cwidth;
	}
	
	private void initDataStructures()
	{
		int num_nodes = post.factory.tree.getNumNodes();
		int num_leaves = post.factory.tree.getNumLeaves();
		
		for (int a=num_leaves; a<num_nodes; a++)
		{
			Ancestor anc = new Ancestor(post.factory, a);
			anc.setCalculationWidth(ancestor_width);
			Posteriors P = new Posteriors(anc);
			ancestor_posteriors[a] = P;
		}
	}
	
	
	protected Profile getAncestorProfile(int family_idx)
	{
		int m = post.factory.table.maxCopies(family_idx);
		int cwidth = getCalculationWidth(m);
		System.out.println("#**AP.gAP "+family_idx+"\tm "+m+"\tcw "+cwidth);
		Profile AP = new Profile(family_idx, cwidth+1);
		return AP;
	}
	
	protected class Profile 
	{
		private Profile(int family_idx, int profile_width)
		{
//			this.family_idx = family_idx;
			this.profile_width = profile_width;
			int num_nodes = post.factory.tree.getNumNodes();
			ancestor = new Posteriors.Profile[num_nodes];
			this.og_profile = post.getPosteriors(family_idx);
			computeAncestors();
		}
		
//		final int family_idx;
		final int profile_width;
		private final Posteriors.Profile[] ancestor;
		private final Posteriors.Profile og_profile;
		private double LL;
		
		private void computeAncestors()
		{
			this.LL = og_profile.inside.getLogLikelihood();
			
			int num_nodes = ancestor.length;
			int root = num_nodes-1; 
			assert (root == post.factory.tree.getRoot()); // since that's how nodes are indexed

			int node = root;
			while (node>=0 && ancestor_posteriors[node]!=null)
			{
				Posteriors AP = ancestor_posteriors[node];
				
				Posteriors.Profile A = AP.getPosteriors(og_profile.inside.family_idx);
				Ancestor anc = (Ancestor) A.inside.getOwner();
				anc.setCalculationWidth(profile_width);
				
				if (node==root)
				{
					A.inside.computeLikelihoods(); // 
				} else
				{
					Likelihood.Profile rootLP =  og_profile.inside; // ancestor[root].inside;
					A.inside.copyLikelihoods(rootLP); // no need to recalculate 
				}
				A.computeOutsideFromRoot(node);
				ancestor[node] = A;
				node--;
			} // while loop over non-leaf nodes 
		}
		
		
		/**
		 * 
		 * @param outsideLL outside ancestor log-likelihoods
		 * @param insideLL inside log-likelihoods
		 * @param extinct extinction probability
		 * @return minimum 2-member array (padded with 0.0 if not computed)
		 */
		private double[] computePosteriors(double[] outsideLL, double[] insideLL, double extinct)
		{
			double[] posteriors;
			Likelihood factory = post.factory;
			
			if (extinct==0.0)
			{
				posteriors = Posteriors.computePosteriors(outsideLL, insideLL);
			} else
			{
				posteriors = outsideLL.clone();
				if (insideLL.length!=0)
				{
	
					double log_e = Math.log(extinct);
					double log_1e = Math.log1p(-extinct);
					
					for (int n=0; n<posteriors.length; n++)
					{
						int ell = 0;
						double ins = insideLL[ell]+n*log_e;
						++ ell;
						while ( ell<=n && ell<insideLL.length)
						{
							double binom = factory.factorials.factln(n)
									-factory.factorials.factln(ell)
									-factory.factorials.factln(n-ell);
							double log_p = binom + ell*log_1e+(n-ell)*log_e;
							ins = Logarithms.add(ins, insideLL[ell]+log_p);
							++ ell;
						}
						posteriors[n] += ins;
					}
				}
				for (int n=0; n<posteriors.length; n++)
				{
					posteriors[n]=Math.exp(posteriors[n]-LL);
				}
			}
			if (posteriors.length==1)
			{
				posteriors = Arrays.copyOf(posteriors, 2);
			}
			
			return posteriors;
		}		
		

		public double[] getNodePosteriors(int node) 
		{
			double[] pN;
			if (post.factory.tree.isLeaf(node))
			{
				pN = og_profile.getNodePosteriors(node);
			}
			else
			{
				double[] B = ancestor[node].getNodeOutside(node);
				double[] C = ancestor[node].inside.getNodeLikelihoods(node);
				double epsi = post.factory.extinction[node];
	
				pN = computePosteriors(B, C, epsi);
			}
			return pN;
		}
		
		public double[] getEdgePosteriors(int node)
		{
			double[] pS;
			if (post.factory.tree.isLeaf(node))
			{
				pS = og_profile.getEdgePosteriors(node);
			} else
			{
				double[] J = ancestor[node].getEdgeOutside(node);
				double[] K = ancestor[node].inside.getEdgeLikelihoods(node);
				double nu = post.factory.extinction[node]
						*post.factory.getDuplicationParameterComplement(node); 
				pS = computePosteriors(J, K, nu);
			}
			return pS;
		}
		
		
	}
	
	private static double mean(double[] p)
	{
		double m = 0.0;
		for (int i=1; i<p.length; i++)
			m += i*p[i];
		return m;
	}
	
	private static double sum(double[] p)
	{
		double sum=0.0;
		for (double x: p) sum+=x;
		return sum;
	}

	private void testAncestors(PrintStream out)
	{
		int nF = post.factory.table.getFamilyCount();
		int num_nodes = post.factory.tree.getNumNodes();
		int num_leaves = post.factory.tree.getNumLeaves();
		
		for (int f=0; f<nF; f++)
		{
			Profile P = getAncestorProfile(f);
			out.println("#FAMILY\t"+f+"\tLL "+P.LL);
			for (int node=num_leaves; node<num_nodes; node++)
			{
				double[] pN = P.getNodePosteriors(node);
				double N = mean(pN);
				out.print(f+"\t"+node+"\tN "+N);
				for (int n=0; n<pN.length; n++)
					out.print("\t"+pN[n]);
				double missN = 1.0-sum(pN);
				out.println("\t// Nmiss "+missN+"\t// "+post.factory.rates.toString(node));
				double[] pS = P.getEdgePosteriors(node);
				double S = mean(pS);
				out.print(f+"\t"+node+"\tS "+S);
				for (int s=0; s<pS.length; s++)
					out.print("\t"+pS[s]);
				double missS = 1.0-sum(pS);
				out.println("\t// Smiss "+missS);
			}
		}
	}
	
	public static void main(String[] args) throws Exception
	{
		CommandLine cli = new CommandLine(args, AncestorPosteriors.class);
        AnnotatedTable table = cli.getTable();
        TreeWithRates rates = cli.getRates();

        PrintStream out = System.out;       
        
        AncestorPosteriors AP = new AncestorPosteriors(rates, table);
        
    	int absolute = 12;
    	double relative = 3.0;
        if (cli.getOptionValue(CommandLine.OPT_TRUNCATE)!=null)
        {
        	absolute = cli.getOptionTruncateAbsolute();
        	relative = cli.getOptionTruncateRelative();
		}
		out.println(CommandLine.getStandardHeader("Truncated computation  for Posteriors (absolute,relative)="
				+absolute+","+relative));
		
		AP.post.setCalculationWidthThresholds(absolute, relative);		
		
		AP.testAncestors(out);
	}

}
