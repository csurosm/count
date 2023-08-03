package count.ds;

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





import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.Map;

/**
 * Direct implementation of the UPGMA algorithm 
 * (unweighted pair-group method with arithmetic mean), from dissimilarity/distance.
 * 
 */
public class UPGMA 
{
	/**
	 * Whether set edge lengths in clustering (true), or rather set all of them to 1.0 (false).
	 */
	private static final boolean SET_EDGE_LENGTHS = true;
	
	private UPGMA(String[] leaf_names, double[][] distances)
	{
		this.leaf_names = leaf_names;
		this.leaf_distances = distances;
	}
	
	private final String[] leaf_names;
	private final double[][] leaf_distances;
	
	/**
	 * Different methods for calculating the distances after joining two clusters. 
	 * 
	 */
	public enum ClusteringPolicy
	{
		/**
		 * Unweighted pair-group with arithmetic average 
		 */
		UPGMA("UPGMA"),
		/**
		 * Weighted pair-group with arithmetic average 
		 */
		WPGMA("WPGMA"),
		/**
		 * Single linkage (min-distance/closest neighbor)
		 */
		SingleLinkage("Single-linkage"),
		/**
		 * Complete linkage (max-distance/farthest neighbor)
		 */
		CompleteLinkage("Complete-linkage");
		
		private ClusteringPolicy(String name)
		{
			this.name = name;
		}
		private final String name;
		@Override
		public String toString() { return name;}
	}
	
	/**
	 * 
	 * Pairwise distances to store on a heap in the clustering algorithm.
	 *
	 */
	private static class Distance implements Comparable<Distance>
	{
		Distance(Phylogeny a, Phylogeny b, double ab_distance)
		{
			this.our_phylo = a;
			this.other_phylo = b;
			this.distance = ab_distance;
		}
		
		Distance(Phylogeny a, Phylogeny b)
		{
			this(a,b,Double.NaN);
		}
		
		private final Phylogeny our_phylo;
		private final Phylogeny other_phylo;
		private final double distance;
		
		
		@Override
		public int hashCode()
		{
			return our_phylo.hashCode() ^ other_phylo.hashCode();
		}
		
		@Override
		public boolean equals(Object o)
		{
			if (o!= null && o instanceof Distance)
			{
				Distance that = (Distance) o;
				return (this.our_phylo==that.our_phylo && this.other_phylo==that.other_phylo)
						|| (this.our_phylo==that.other_phylo && this.other_phylo==that.our_phylo);
			} else
				return super.equals(o);
		}
		
		@Override
		public int compareTo(Distance that)
		{
			return Double.compare(this.distance, that.distance);
		}
		
	}
	
	public static Phylogeny buildTree(ProfileTable table, ClusteringPolicy clustering)
	{
		double[][] distances = dissimilaritySorensenDice(table);
		return buildTree(table.getTaxonNames(), distances, clustering);
	}
	
	public static Phylogeny buildTree(ProfileTable table, Dissimilarity distance, ClusteringPolicy clustering)
	{
		double[][] distances = dissimilarityMatrix(table, distance);
		String[] leaf_names = table.getTaxonNames();
		UPGMA factory = new UPGMA(leaf_names, distances);
		return factory.buildTree(clustering);
	}	
	
	public static Phylogeny buildTree(String[] leaf_names, double[][] distances, ClusteringPolicy clustering)
	{
		UPGMA factory = new UPGMA(leaf_names, distances);
		return factory.buildTree(clustering);
	}
	

	private Phylogeny buildTree(ClusteringPolicy clustering)
	{
		Map<Phylogeny, Heap<Distance>> node_distances = new HashMap<>();
		List<Phylogeny> nodes = new ArrayList<>(); 
		
		for (int leaf=0; leaf<leaf_names.length; leaf++)
		{
			Phylogeny leaf_only = new Phylogeny();
			leaf_only.getRootNode().setName(leaf_names[leaf]);
			nodes.add(leaf_only);
			
		}
		for (int our=0; our<nodes.size(); our++)
		{
			Phylogeny our_phylo = nodes.get(our);
			Heap<Distance> pairs = new Heap<>();
			double[] our_dist = leaf_distances[our];
			for (int other=0; other<our_dist.length; other++)
				if (our!=other)
				{
					Phylogeny other_phylo = nodes.get(other);
					
					Distance D = new Distance(our_phylo, other_phylo, our_dist[other]);
					pairs.add(D);
				}
			node_distances.put(our_phylo, pairs);
		}
		
		Phylogeny join_phylo = null; // last join will be kept when the loop ends 
		
		while (node_distances.size()>1) // number of remaining free nodes; 2 remove and 1 put per iteration
		{
			// 1.1 find closest pair
			Distance mindist = null;
			for (Heap<Distance> dist: node_distances.values())
			{
				Distance shortest = dist.peek();
				if (mindist==null || shortest.compareTo(mindist)<0)
				{
					mindist = shortest; // will be found twice, since stored at each member
				}
			}
			// 2. join the pair members 
			Phylogeny our_phylo = mindist.our_phylo;
			Phylogeny.Node our_root = our_phylo.getRootNode();
			int our_size = our_phylo.getNumLeaves();
			Phylogeny other_phylo = mindist.other_phylo;
			Phylogeny.Node other_root = other_phylo.getRootNode();
			int other_size = other_phylo.getNumLeaves();
			int join_size = our_size+other_size;
			
			join_phylo = new Phylogeny();
			Phylogeny.Node join_root = join_phylo.getRootNode();
			join_root.addChild(our_phylo.getRootNode());
			join_root.addChild(other_phylo.getRootNode());
			assert (join_phylo.getNumLeaves()==join_size);
			nodes.add(join_phylo);
			
			// 3. calculate the average distances from other nodes 
			Heap<Distance> our_dist = node_distances.remove(our_phylo);
			Heap<Distance> other_dist = node_distances.remove(other_phylo);
			Heap<Distance> join_dist = new Heap<>(); // to be filled in 
			for (Phylogeny third: node_distances.keySet()) // without our and other
			{
				Heap<Distance> third_dist = node_distances.get(third);
				Distance our_third = our_dist.get(new Distance(our_phylo,third));
				Distance other_third = other_dist.get(new Distance(other_phylo,third));
				assert (our_third != null);
				assert (other_third != null);
				
				third_dist.remove(our_third);
				third_dist.remove(other_third);
				
				
				double combined_distance=Double.NaN;
				
				switch (clustering)
				{
					case UPGMA:
						double usum = our_size * our_third.distance + other_size*other_third.distance;
						combined_distance = usum/join_size;
						break;
					case WPGMA:
						double wsum = our_third.distance+other_third.distance;
						combined_distance = wsum/2.0;
						break;
					case SingleLinkage:
						combined_distance = Double.min(our_third.distance, other_third.distance);
						break;
					case CompleteLinkage:
						combined_distance = Double.max(our_third.distance, other_third.distance);
						break;
				}
				Distance join_third = new Distance(join_phylo, third, combined_distance);
				join_dist.add(join_third);
				third_dist.add(join_third); // member order is immaterial 
				
			}
			node_distances.put(join_phylo, join_dist);
			
			// 4. set branch lengths
			if (SET_EDGE_LENGTHS)
			{
				double half_dist = mindist.distance/2.0;
				double child_height = 0.0;
				Phylogeny.Node child = our_root;
				while (!child.isLeaf())
				{
					Phylogeny.Node grandchild = child.getChild(0); // ultrametric tree, any child is ok 
					child_height += grandchild.getLength();
					child = grandchild;
				}
				double len = half_dist-child_height;
				our_root.setLength(len);
				other_root.setLength(len);
			} else
			{
				our_root.setLength(1.0);
				other_root.setLength(1.0);
			}
//			System.out.println("#**UPGMA.bT join/"+clustering+" "
//					+(leaf_names.length-node_distances.size())
//					+"\tmindist "+mindist.distance
//					+"\tedgelen "+len
//					+"\ttree "+NewickParser.printTree(join_phylo));
			// 
		}
		join_phylo.computeIndexes(); // probably no need to be eager here, bc computeIndexes is called the first time a node property is queried by node index 
		return join_phylo;
	}

	/**
	 * Statistics to measure the dissimilarity between genomes
	 * based on copy numbers. 
	 * 	 
	 */
	public enum Dissimilarity
	{
		/**
		 * Bray-Curtis dissimilarity : dist(i,j)= 1.0-2.0*(copies in common bw i and j)/(sum of copies in i and j); common copies = minimum copy number  
		 */
		BRAY_CURTIS("Bray-Curtis"),
		/**
		 * Sorensen-Dice dissimilarity : dist(i,j) = 1.0-2.0*(families in common bw i and j)/(sum of families in i and j)
		 */
		SORENSEN_DICE("Sørensen-Dice"),
		/**
		 * Jaccard distance: dist(i,j) = 1.0-(families in common bw i and j)/(families in i or j)
		 */
		JACCARD("Jaccard"),
		/**
		 * Weighted Jaccard distance: dist(i,j) = 1.0-(copies in common bw i and j)/(copies in i or j); common = min, union = max
		 */
		WEIGHTED_JACCARD("Weighted Jaccard");
		private Dissimilarity(String name)
		{
			this.name = name;
		}
		private final String name;
		@Override
		public String toString()
		{
			return name;
		}
		
		

	}
	
	/**
	 * Dissimilarity values computed by arbitrary method. 
	 * 
	 * @param table profile table of copy counts 
	 * @param method one of {@link Dissimilarity} values; if null, default is assumed
	 * 
	 * 
	 * @return symmetric dissimilarity matrix 
	 */
	public static double[][] dissimilarityMatrix(ProfileTable table, Dissimilarity method)
	{
		if (method == null)
			method = Dissimilarity.BRAY_CURTIS;
			
		int n = table.getTaxonCount();

		double[][] fam_intersect = new double[n][];
		double[][] fam_union = new double[n][];
		double[][] copy_intersect = new double[n][];
		double[][] copy_union = new double[n][];

		for (int i=0; i<n; i++)
		{
			fam_intersect[i] = new double[i+1];
			fam_union[i] = new double[i];
			copy_intersect[i] = new double[i+1];
			copy_union[i] = new double[i];
		}
		
		
		UniqueProfileTable utable = null;
		if (table instanceof UniqueProfileTable)
			utable = (UniqueProfileTable)table;
		int nF = table.getFamilyCount();
		
		for (int f=0; f<nF; f++)
		{
			int[] copies = table.getFamilyProfile(f);
			int mul = utable==null?1:utable.getMultiplicity(f); // family multiplicity
			
			for (int i=0; i<n; i++)
			{				
				boolean ipresent = 0<copies[i];
				if (ipresent) fam_intersect[i][i] += mul; // fam_count[i] = fam_intersect[i][i] 
				copy_intersect[i][i] += mul*copies[i]; // copy_count[i] = copy_intersect[i][i]
				
				for (int j=0; j<i; j++)
				{
					boolean jpresent = 0<copies[j];
					if (ipresent)
					{
						fam_union[i][j] += mul;
						if (jpresent)
							fam_intersect[i][j] += mul;
					} else if (jpresent)
					{
						fam_union[i][j] += mul;
					}
					copy_intersect[i][j] += mul*Integer.min(copies[i], copies[j]);
					copy_union[i][j] += mul*Integer.max(copies[i], copies[j]);
				} // for j < i
			}  // for i
		} // for families
		double[][] D = new double[n][n]; // return value
		
		
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<i; j++)
			{
				double sim;
				switch(method)
				{
					case BRAY_CURTIS:
						sim = 2.0*copy_intersect[i][j]/(copy_intersect[i][i]+copy_intersect[j][j]);
						break;
					case SORENSEN_DICE:
						sim = 2.0*fam_intersect[i][j]/(fam_intersect[i][i]+fam_intersect[j][j]);
						break;
					case JACCARD:
						sim = fam_intersect[i][j]/fam_union[i][j];
						break;
					case WEIGHTED_JACCARD:
						sim = copy_intersect[i][j]/copy_union[i][j];
						break;
					default:
						sim=Double.NaN;
				}
				D[i][j] = 1.0-sim;
				D[j][i] = D[i][j];
			}
		}
		
		return D;
	}
		
	
	/**
	 * Dissimilarity values computed by Sorensen-Dice (Bray-Curtis) over presence-absence.
	 * 
	 * @param table
	 * @return
	 */
	public static double[][] dissimilaritySorensenDice(ProfileTable table)
	{
		int n = table.getTaxonCount();
		double[][] D = new double[n][n]; // return value
		
		int[] nfam = new int[n];
		int[] ncopy = new int[n];
		
		UniqueProfileTable utable = null;
		if (table instanceof UniqueProfileTable)
			utable = (UniqueProfileTable)table;
		
		int nF = table.getFamilyCount();
		int[] present=new int[n];
		for (int f=0; f<nF; f++)
		{
			int[] copies = table.getFamilyProfile(f);
			int mul = utable==null?1:utable.getMultiplicity(f);
			
			// extract the families in the table
			int nlin = 0; 
			for (int leaf=0; leaf<copies.length; leaf++)
			{
				if (copies[leaf]>0)
				{
					present[nlin++]=leaf;
				}
			}
			
			for (int i=0; i<nlin; i++)
			{
				int li = present[i];
				nfam[li]+=mul;
				ncopy[li] += mul*copies[li];
				
				for (int j=0; j<i; j++)
				{
					int lj = present[j];
					D[li][lj] += mul; // D counts intersection size 
				}
			}
		}
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<i; j++)
			{
				double common = D[i][j]/(nfam[i]+nfam[j]);
				D[i][j] = 1.0-2.0*common;
				D[j][i] = D[i][j];
			}
		}
		
		return D;
	}
	
}
