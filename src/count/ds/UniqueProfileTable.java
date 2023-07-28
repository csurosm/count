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
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
/**
 * A compressed {@link ProfileTable} that stores the multiplicities for 
 * repeated copy-number profiles. 
 * 
 * @author csuros
 *
 */
public class UniqueProfileTable implements ProfileTable
{
	public UniqueProfileTable(ProfileTable full_table)
	{
		this.full_table = full_table;
		this.uniq_index = new HashMap<>();
		this.unique_profiles = new ArrayList<>();
		this.last_occurrence = new ArrayList<>();
		this.previous_occurrence = new int[full_table.getFamilyCount()];
		this.uniq_family = new int[previous_occurrence.length];
		initDataStructures();
	}
	
	private UniqueProfileTable(UniqueProfileTable same_families)
	{
		this.full_table = same_families.full_table;
		this.unique_profiles = same_families.unique_profiles;
		this.uniq_index = same_families.uniq_index;
		this.last_occurrence = same_families.last_occurrence;
		this.previous_occurrence = same_families.previous_occurrence;
		this.uniq_family = same_families.uniq_family;
	}
	
	private final ProfileTable full_table;
	private final List<Profile> unique_profiles;
	private final Map<Profile, Integer> uniq_index;
	private final List<Integer> last_occurrence;
	private final int[] previous_occurrence;
	private final int[] uniq_family; 
	
	private void initDataStructures()
	{
		int nF = full_table.getFamilyCount();
		
		for (int f=0; f<nF; f++)
		{
			Profile P = new Profile(f);
			int uniq_idx;
			if (uniq_index.containsKey(P))
			{
				uniq_idx = uniq_index.get(P);
				Profile first = unique_profiles.get(uniq_idx);
				first.multiplicity++;
				previous_occurrence[f] = last_occurrence.get(uniq_idx);
				last_occurrence.set(uniq_idx, f);
			} else
			{
				uniq_idx = unique_profiles.size();
				uniq_index.put(P, uniq_idx);
				unique_profiles.add(P);
				last_occurrence.add(f);
				previous_occurrence[f]=-1;
			}
			uniq_family[f]=uniq_idx;
		}
		
//		System.out.println("#*UPT.iDS full "+full_table.getFamilyCount()+"\tuniq "+unique_profiles.size());
	}
	
	@Override 
	public int getFamilyCount()
	{
		return unique_profiles.size();
	}
	
	public int getTotalFamilyCount()
	{
		return full_table.getFamilyCount();
	}
	
	@Override
	public int[] getFamilyProfile(int family)
	{
		return unique_profiles.get(family).get();
	}
	
	@Override
	public int getTaxonCount()
	{
		return full_table.getTaxonCount();
	}
	
	@Override 
    public String[] getTaxonNames()
    {
		return full_table.getTaxonNames();
    }
	

	public int getMultiplicity(int family)
	{
		return unique_profiles.get(family).multiplicity;
	}
	
	
	
//	public Iterator<Integer> getOriginalFamilies(int uniq_idx)
//	{
//		class OccurrenceIterator implements Iterator<Integer>
//		{
//			OccurrenceIterator()
//			{
//				next_occurrence = last_occurrence.get(uniq_idx);
//			}
//			private int next_occurrence;
//			@Override 
//			public boolean hasNext()
//			{
//				return (next_occurrence >= 0);
//			}
//			@Override 
//			public Integer next()
//			{
//				int next = next_occurrence;
//				if (next!=-1)
//					next_occurrence=previous_occurrence[next_occurrence];
//				return next;		
//			}
//		}
//		return new OccurrenceIterator();
//	}
	
	private class Profile
	{
		Profile(int family_idx)
		{
			this.family_idx = family_idx;
			this.hashCode =  Arrays.hashCode(get());
			this.multiplicity = 1;
		}
		
		private final int family_idx;
		private final int hashCode;
		private int multiplicity;
		
		@Override
		public int hashCode()
		{
			return hashCode;
		}
		
		int[] get()
		{
			return full_table.getFamilyProfile(family_idx);
		}
		
		@Override
		public boolean equals(Object o)
		{
			if (o instanceof Profile)
			{
				Profile that = (Profile) o;
				return Arrays.equals(this.get(), that.get());
			} else
				return super.equals(o);
		}
		
		
	}
	
	public UniqueProfileTable mappedToTree(IndexedTree tree)
	{
	   String[] table_taxon_names = getTaxonNames();
	   Map<String, Integer> table_column_index = new HashMap<>();
	   for (int t=0; t<table_taxon_names.length; t++)
		   table_column_index.put(table_taxon_names[t], t);
	   String[] tree_taxon_names = tree.getLeafNames();
//	   System.out.println("#**AT.mTT table "+Arrays.toString(table_taxon_names)
//	   	+"\ttree "+Arrays.toString(tree_taxon_names));

	   
	   //	   if (tree_taxon_names.length != table_taxon_names.length)
//	   {
//		   throw new IllegalArgumentException("Cannot match all columnns: requested "
//				   	+Arrays.toString(tree_taxon_names)
//				   	+";\tknown "+Arrays.toString(table_taxon_names));
//	   }
	   int[] column_order = new int[tree_taxon_names.length];
	   Arrays.fill(column_order, -1);
	   boolean has_unmapped_leaves = false;
	   for (int c=0; c<tree_taxon_names.length; c++)
	   {
		   String name = tree_taxon_names[c];
		   if (table_column_index.containsKey(name))
		   {	
			   if (column_order[c]!=-1)
				   has_unmapped_leaves=true;
			   column_order[c]=table_column_index.get(name);
		   }
		   else
		   {
			   has_unmapped_leaves = true;
			   column_order[c]=-1;
		   }
	   }
	   if (has_unmapped_leaves)
	   {
		   System.err.println("Cannot match all leaf names to columnns: requested "
				   	+Arrays.toString(tree_taxon_names)
				   	+";\tknown "+Arrays.toString(table_taxon_names));
		   throw new IllegalArgumentException();
	   }
	   
	   /**
	    * Same families, but with the columns reordered. 
	    * 
	    * @author csuros
	    *
	    */
	   class MappedTable extends UniqueProfileTable
	   {
		   private final int[] column_order;
		   /**
		    * 
		    * 
		    * @param column_order every entry must be a column in the parent table
		    */
		   MappedTable(int[] column_order)
		   {
			   super(UniqueProfileTable.this);
			   boolean same_order = true;
			   for (int c=0; c<column_order.length && same_order; c++)
				   same_order = column_order[c]==c;
			   this.column_order = same_order?null:column_order;
		   }
		   @Override
		   public String[] getTaxonNames()
		   {
			   String[] orig_taxons = super.getTaxonNames();
			   String[] our_taxons;
			   if (column_order==null)
				   our_taxons = orig_taxons;
			   else
			   {
				   our_taxons = new String[column_order.length];
				   for (int c=0; c<column_order.length; c++)
					   our_taxons[c] = orig_taxons[column_order[c]];
			   }
			   return our_taxons;
		   }
		   
		   @Override
		   public int getTaxonCount()
		   {
			   return column_order==null?super.getTaxonCount():column_order.length;
			   
		   }
		   @Override
		   public int[] getFamilyProfile(int f)
		   {
			   int[] orig_profile = super.getFamilyProfile(f);
			   int[] our_profile = null;
			   if (column_order==null)
				   our_profile = orig_profile;
			   else
			   {
				   our_profile = new int[column_order.length];
				   for (int c=0; c<column_order.length; c++)
					   our_profile[c]=orig_profile[column_order[c]];
			   }
			   return our_profile;
		   }
	   }
	   
	   return new MappedTable(column_order);
	}

}
