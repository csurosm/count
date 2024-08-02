package count.io;

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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.StringReader;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;

import java.util.Map;

import java.util.function.Predicate;

import javax.xml.parsers.SAXParserFactory;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.DefaultHandler;


import count.ds.Phylogeny;
import count.Count;
import count.ds.AnnotatedTable;
import count.ds.IndexedTree;
import count.ds.TreeComparator;
import count.gui.PosteriorsView;
import count.model.GammaInvariant;
//import count.model.GammaInvariant;
import count.model.MixedRateModel;

/**
 * Hierarchy data structure for storing a session's members: phylogenies, 
 * rate models, and profile tables. 
 *
 */
public class ModelBundle implements CountXML
{	
	private static final boolean WANT_SECTION_TAGS = false;
	private static final boolean DEBUG_PARSING = false;
	
	public ModelBundle(String identifier)
	{
		this.entry_identifiers = new HashMap<>();
		this.root_identifier = identifier==null?"":identifier;
		this.root = new Entry(null,null,null,identifier);
	}
	
	
	private final String root_identifier;
	private final Map<String, Entry> entry_identifiers; 
	
	private final Entry root; // root has no data, only an identifier 
	
	
	public Entry getRoot() { return root;}
	
	public Entry add(DataFile<Phylogeny> tree_data, DataFile<MixedRateModel> rates_data)
	{
		Entry T = getRoot().addTree(tree_data);
		Entry R = T.addRates(rates_data);
		return R;
	}	
	
	public Entry[] allEntries()
	{
		List<Entry> allEntries = collectAllEntries(getRoot(), null);
		return allEntries.toArray(new Entry[0]);
	}
	
	private List<Entry> collectAllEntries(Entry E, List<Entry> entry_list)
	{
		if (entry_list == null) entry_list = new ArrayList<>();
		entry_list.add(E);
		int nc = E.getNumChildren();
		for (int ci=0; ci<nc; ci++)
		{
			Entry C = E.getChild(ci);
			entry_list = collectAllEntries(C, entry_list);
		}
		return entry_list;
	}
	
	
	public MixedRateModel[] allModels()
	{
		List<MixedRateModel> allModels = collectAllModels(getRoot(), null);
		return allModels.toArray(new MixedRateModel[0]);
	}
	
	private List<MixedRateModel> collectAllModels(Entry E, List<MixedRateModel> model_list)
	{
		if (model_list==null) model_list = new ArrayList<>();
		if (E.isRatesEntry())
			model_list.add(E.rates_data.getContent());
		int nc = E.getNumChildren();
		for (int ci=0; ci<nc; ci++)
		{
			Entry C = E.getChild(ci);
			model_list = collectAllModels(C, model_list);
		}
		return model_list;
	}
	
	public Entry getEntry(String id)
	{
		return entry_identifiers.get(id);
	}
	
	/**
	 * Recursive code to delete ids in the subtree, called before 
	 * the physical removal of an Entry via {@link Entry#removeFromParent()}
	 * 
	 * @param node
	 */
	private void deleteId(Entry node)
	{
		entry_identifiers.remove(node.getId());
		for (int ci=node.getNumChildren(); 0<ci; )
		{
			--ci;
			Entry child = node.getChild(ci);
			deleteId(child);
		}
	}
	
	/**
	 * The main phylogeny for this bundle (defining the default leaf ordering):
	 * the first tree child node of the root.   
	 * 
	 * 
	 * @return null if undefined
	 */
	public Entry getMainTree()
	{
		Entry R = getRoot();
		Entry getMainTree = R.getChild(0, c->c.isTreeEntry());
		
//		System.out.println("#**MB.gMT root "+R+"\ttree "+getMainTree);
		
		return getMainTree;
//		
//		Entry main_tree = null;
//		if (R.children != null)
//		{
//			for (Entry C: R.children)
//			{
//				if (C.isTreeEntry())
//				{
//					main_tree = C;
//					break;
//				}
//			}
//		}
//		return main_tree;
	}
	
	
	public class Entry
	{
		/**
		 * At most one of {@link #tree_data} and {@link #rates_data} is null.
		 */
		private DataFile<Phylogeny> tree_data;
		private DataFile<MixedRateModel> rates_data;
		private DataFile<AnnotatedTable> table_data = null;
		private final String entry_id;
		
		private Entry parent;
		private List<Entry> children=null;
		
		private Properties extra_attributes;
		
		public MixedRateModel getModel()
		{
			return rates_data==null?null:rates_data.getContent();
		}
		
		public String getId() { return entry_id;}
		
		public Entry getRoot()
		{
			if (this.isRoot()) return this;
			else return parent.getRoot();
		}
		
		
//		private Entry(String id, Entry parent)
//		{
//			this.entry_id = id;
//			tree_data = null;
//			rates_data = null;
//					
//		}
		public void setAttribute(String att_key, String att_value)
		{
			if (extra_attributes==null)
				extra_attributes = new Properties();
			extra_attributes.setProperty(att_key, att_value);
		}
		
		
		private Entry(DataFile<Phylogeny> tree_data, DataFile<MixedRateModel> rates_data, Entry parent, String id)
		{
			if (rates_data==null)
			{
				this.tree_data = tree_data;
				this.rates_data = null;
				if (parent==null)
				{
					this.entry_id = id==null?root_identifier:id;
				} else
				{
					int tidx = parent.countChildren(E->E.isTreeEntry());
					this.entry_id = id==null?parent.entry_id+"."+PREFIX_TREEID+tidx:id;
				}
			} else 
			{ // add tree and rates, or only rates
				if (tree_data != null)
				{
					Entry T = new Entry(tree_data, null, parent,null);
					parent = T;
				}
				assert (parent != null);
				this.tree_data = null;
				this.rates_data = rates_data;
				int ridx = parent.countChildren(E->E.isRatesEntry());
				this.entry_id = id==null?parent.entry_id+"."+PREFIX_RATESID+ridx:id;
			}

			entry_identifiers.put(entry_id, this);
			
			if (parent != null)
				parent.addChild(this);
		}
		
		private Entry(DataFile<AnnotatedTable> table_data, Entry parent, String id)
		{
			this.tree_data = null;
			this.rates_data = null;
			this.table_data = table_data;
			
			int didx = parent.countChildren(E->E.isTableEntry());
			this.entry_id = id==null?parent.entry_id+"."+PREFIX_DATAID+didx:id;
			entry_identifiers.put(entry_id, this);
			parent.addChild(this);
			
			
		}
		
		
		public boolean isRatesEntry() { return this.rates_data != null;}
		public boolean isTreeEntry() { return this.tree_data != null;}
		public boolean isTableEntry() { return this.table_data != null;}
		
		public int getNumChildren()
		{
			return children==null?0:children.size();
		}
		
		public Entry getParent()
		{
			return parent;
		}
		
		public Entry getChild(int child_idx)
		{
			return children.get(child_idx);
		}
		
		public boolean isRoot()
		{
			return parent == null;
		}
		
		
		public String getAttributeValue(String att_key)
		{
			String att_val = null;
			if (extra_attributes!=null)
			{
				att_val = extra_attributes.getProperty(att_key);
			}
			return att_val;
		}
		/**
		 * Value of type attribute 
		 * @return
		 */
		public String getType()
		{
			return getAttributeValue(ATT_TYPE);
		}
		
		
		public DataFile<Phylogeny> getTreeData()
		{
			return this.tree_data;
		}
		
		public DataFile<MixedRateModel> getRatesData()
		{
			return this.rates_data;
		}
		
		public DataFile<AnnotatedTable> getTableData()
		{
			return this.table_data;
		}
		
		
		public void setTableData(DataFile<AnnotatedTable> table_data)
		{
			assert (this.isTableEntry());
			this.table_data = table_data;
		}

		public Entry addTree(DataFile<Phylogeny> tree_data)
		{
			return this.addTree(tree_data, null);
		}
		
		private Entry addTree(DataFile<Phylogeny> data, String id)
		{
			if (!this.isRoot() && !this.isTreeEntry())
			{
				throw new IllegalArgumentException("Add tree "+data+" under the root or another phylogeny");
			}
			
			Entry find = data.getContent()==null?null:getRoot().find(data.getContent());
			Entry addTree;
			if (find == null)
			{
				 addTree = new Entry(data, null, this, id);
			}
			else
			{
				System.out.println("#**MB.E.addTree\t"+this+"\ttree already there for "+tree_data.toString());
				addTree = find;
			}
			return addTree;
		}
		
		public Entry addRates(DataFile<MixedRateModel> rates_data)
		{
			return addRates(rates_data, null);
		}
		
		private Entry addRates(DataFile<MixedRateModel> data, String id)
		{			
			Entry T = this;
			while (T!= null && !T.isTreeEntry()) T = T.getParent();
			boolean got_tree = false;
			if (T != null && T.isTreeEntry())
			{
				if (data.getContent()==null) got_tree = true;
				else
				{
					IndexedTree tree = data.getContent().getBaseModel().getTree();
					IndexedTree parent_tree = T.tree_data.getContent();
					got_tree = TreeComparator.sameTopology(true, parent_tree, tree);
				} 
			}
			if (!got_tree)
			{
				throw new IllegalArgumentException("Add rates "+rates_data+" under its proper phylogeny");
			}
			Entry addRates = new Entry(null, data, this, id);
			return addRates;
		}
				
		private Entry addTable(DataFile<AnnotatedTable> data, String id)
		{
			if (!this.isRoot() && !this.isTableEntry())
				throw new IllegalArgumentException("Add table "+data+" under the root or another data item/table");
			
			Entry addTable = new Entry(data, this, id);
			if (data.getContent() != null)
			{
				addTable.setAttribute(ATT_BINARY, Boolean.toString(data.getContent().isBinaryTable()));
			}
			return addTable;
		}
		
		public Entry addTable(DataFile<AnnotatedTable> data)
		{
			return addTable(data, null);
		}
		
		
		private Entry addDataItem(DataFile<AnnotatedTable> emptycontent, String type, String id)
		{
			if (!this.isTableEntry())
				throw new IllegalArgumentException("Add ata item under another data item or table");
			Entry item = addTable(emptycontent, id);
			item.setAttribute(ATT_TYPE, type);
			return item;
		}
		
		public Entry addDataItem(DataFile<AnnotatedTable> emptycontent, String type)
		{
			return addDataItem(emptycontent, type, null);
		}
		
		private void addChild(Entry C)
		{
			if (children==null)
				children = new ArrayList<>();
			children.add(C);
			C.parent = this;
		}
		
		/**
		 * Removes this entry (with all its subtree) from the parent.
		 * 
		 * @return
		 */
		public Entry removeFromParent()
		{
			ModelBundle.this.deleteId(this);
			Entry P = getParent();
			P.children.remove(this);
			return P;
		}
		
		/**
		 * Counts children with a given predicate true. 
		 * 
		 * @param to_be_counted
		 * @return
		 */
		public int countChildren(Predicate<Entry> to_be_counted)
		{
			int count = 0;
			if (children != null)
			{
				for (Entry child: children)
				{
					if (to_be_counted.test(child)) count++;
				}
			}
			return count;
		}
		
		/**
		 * Child enumeration by some selection type.
		 * 
		 * @param child_idx between 0 and {#link {@link #countChildren(Predicate)}}-1
		 * @param child_passes test for counting the children 
		 * @return null if no such child, or the child with the given property 
		 */
		public Entry getChild(int child_idx, Predicate<Entry> child_passes)
		{
			Entry child=null;
			if (children != null)
			{
				for (int ci=0; ci<children.size() && child==null; ci++)
				{
					Entry C = children.get(ci);
					if (child_passes.test(C))
					{
						if (child_idx==0) child = C;
						else child_idx--;
					}
				}
			}
			return child;
			
		}
		
		/**
		 * Finds the most recent ancestor (inclusive, with this)  
		 * that satisfies a predicate.
		 * 
		 * @param ancestor_selection
		 * @return null if no such node on the path to the root
		 */
		public Entry getClosestAncestor(Predicate<Entry> ancestor_selection)
		{
			Entry node = this;
			while (node != null && !ancestor_selection.test(node))
				node = node.getParent();
			return node;		
		}
			
		
//		private int countTreeChildren()
//		{
//			int count = 0;
//			if (children != null)
//			{
//				for (Entry child: children)
//				{
//					if (child.isTreeEntry()) count++;
//				}
//			}
//			return count;
//		}
		
//		private int countRateChildren()
//		{
//			int count = 0;
//			if (children != null)
//			{
//				for (Entry child: children)
//				{
//					if (child.isRatesEntry()) count++;
//				}
//			}
//			return count;
//		}

		/**
		 * @param tree_or_rates_or_table should not be null
		 * @return
		 */
		private Entry find(Object tree_or_rates_or_table)
		{
			Entry find = null;
			if ((tree_data != null && tree_data.getContent().equals(tree_or_rates_or_table))
				|| (rates_data != null && rates_data.getContent().equals(tree_or_rates_or_table))
				|| (table_data !=null && tree_or_rates_or_table.equals(table_data.getContent()))
				)
			{
				find = this;
			} else  
			{ // traverse subtree
				int nc  = getNumChildren();
				for (int ci=0; ci<nc && find==null; ci++)
				{
					Entry C = children.get(ci);
					find = C.find(tree_or_rates_or_table);
				}
			}
			return find;
		}
		
		private StringBuilder appendExtraAttributes(StringBuilder sb)
		{
			if (extra_attributes!=null)
			{
				for (String att_key: extra_attributes.stringPropertyNames())
				{
					sb = appendAttribute(sb, att_key, extra_attributes.getProperty(att_key));
//					sb.append(" ").append(att_key).append("=\"")
//					.append(extra_attributes.getProperty(att_key)).append("\"");
				}
			}
			return sb;
		}
		
		private StringBuilder appendAttribute(StringBuilder sb, String att_key, String att_value)
		{
			sb.append(" ").append(att_key).append("=\"").append(att_value).append("\"");
			return sb;
		}
		
		/**
		 * Recursive function for constructing the XML string 
		 * in a StringBuilder. 
		 * 
		 * @param sb
		 * @return
		 */
		private StringBuilder buildXMLString(StringBuilder sb)
		{
			if (sb==null) sb = new StringBuilder();
			if (parent==null) // root 
			{
				sb.append("<").append(EMT_SESSION);
				sb = appendAttribute(sb, ATT_ID, entry_id);
				sb = appendExtraAttributes(sb);
				sb.append(">\n");
			} else 
			{
				String parent_id = parent.entry_id;
				String content;
				if (this.isTreeEntry())
				{
					sb.append("<").append(EMT_TREE);
					sb = appendAttribute(sb, ATT_ID, entry_id);
					sb = appendAttribute(sb, ATT_NAME, tree_data.getFile().getName());
					sb = appendAttribute(sb, ATT_PARENT, parent_id);
//					.append(" id=\"").append(entry_id).append("\"");
//					sb.append(" name=\"").append(tree_data.getFile().getName()).append("\"");
//					sb.append(" dirty=\"").append(tree_data.isDirty()).append("\"");
//					sb.append(" parent=\"").append(parent_id).append("\"");
					sb = appendExtraAttributes(sb);
					sb.append(" >");
					content = NewickParser.printTree(tree_data.getContent());
				} else if (this.isRatesEntry())
				{
					if (WANT_SECTION_TAGS && parent.isTreeEntry())
					{
						int dot = entry_id.lastIndexOf('.');
						String id = entry_id.substring(0, dot+2); // includes ".R"
						sb.append("<").append(EMT_RATES);
						sb = appendAttribute(sb, ATT_ID, id);
						sb = appendAttribute(sb, ATT_PARENT, parent_id);
//						.append(" id=\"").append(id).append("\"");
//						sb.append(" parent=\"").append(parent_id).append("\"");
						sb.append(" >");
						parent_id = id;
						sb.append("\n");
					}
					sb.append("<").append(EMT_RATEMODEL);
					sb = appendAttribute(sb, ATT_ID, entry_id);
					sb = appendAttribute(sb, ATT_NAME, rates_data.getFile().getName());
					sb = appendAttribute(sb, ATT_PARENT, parent_id);
//					.append(" id=\"").append(entry_id).append("\"");
//					sb.append(" name=\"").append(rates_data.getFile().getName()).append("\"");
////					sb.append(" dirty=\"").append(rates_data.isDirty()).append("\"");
//					sb.append(" parent=\"").append(parent_id).append("\"");
					sb = appendExtraAttributes(sb);
					sb.append(" >");
					content = RateVariationParser.printRates(rates_data.getContent());
				} else if (this.isTableEntry())
				{
					if (WANT_SECTION_TAGS && parent.isRoot())
					{
						int dot = entry_id.lastIndexOf('.');
						String id = entry_id.substring(0, dot+2); // includes ".R"
						sb.append("<").append(EMT_DATA);
						sb = appendAttribute(sb, ATT_ID, id);
						sb = appendAttribute(sb, ATT_PARENT, parent_id);
						sb.append(" >");
						parent_id = id;
						sb.append("\n");						
					} 
					if (table_data.getContent()==null) // data item
					{
						sb.append("<").append(EMT_DATAITEM);
						content = null;
					} else
					{
						sb.append("<").append(EMT_TABLE);
						content = TableParser.getFormattedTable(table_data.getContent(), true);
					}
					sb = appendAttribute(sb, ATT_ID, entry_id);
					sb = appendAttribute(sb, ATT_NAME, table_data.getFile().getName());
					sb = appendAttribute(sb, ATT_PARENT, parent_id);
					sb = appendExtraAttributes(sb);
					sb.append(" >");
				} else
				{
					throw new IllegalStateException("XMLString: cannot deduce Entry type for "+this);
				}
				sb.append("\n");
				if (content != null)
				{
					sb.append("<![CDATA[\n").append(content).append("\n]]>\n");
				}
			}
			int nc = getNumChildren();
			for (int ci=0; ci<nc; ci++)
			{
				Entry C = getChild(ci);
				sb = C.buildXMLString(sb);
			}
			sb.append("</");
			if (this.isRoot()) // root 
			{
				sb.append(EMT_SESSION);
			} else 
			{
				if (this.isTreeEntry())
				{
					sb.append(EMT_TREE);
				} else if (this.isRatesEntry())
				{
					if (WANT_SECTION_TAGS && parent.isTreeEntry())
					{
						sb.append(EMT_RATEMODEL).append(">").append("\n</");
						sb.append(EMT_RATES);
					} else
					{
						sb.append(EMT_RATEMODEL);
					}
				} else 
				{
					assert (this.isTableEntry());
					String emt = table_data.getContent()==null?EMT_DATAITEM:EMT_TABLE;
					if (WANT_SECTION_TAGS && parent.isRoot())
					{
						sb.append(emt).append(">").append("\n</");
						sb.append(EMT_DATA);
					} else
					{
						sb.append(emt);
					}
				}
			}
			sb.append(">\n");
			return sb;
		}
		
		
		@Override
		public String toString()
		{
			return "E["+this.entry_id+", parent="+(parent==null?null:parent.entry_id)+"]";
		}
	}
	
	public String toXMLString()
	{
		StringBuilder sb = getRoot().buildXMLString(null);
		return sb.toString();
	}
	
	
	public static List<ModelBundle> readBundle(BufferedReader BR) throws IOException, SAXException, javax.xml.parsers.ParserConfigurationException
	{
		SAXParserFactory factory = SAXParserFactory.newInstance();
		XMLReader R = factory.newSAXParser().getXMLReader();
		Handler H = new Handler();
		R.setContentHandler(H);
		R.parse(new InputSource(BR));
		
		H.fillContent();
		return H.bundle_list;
	}
	
	/**
	 * Writes a complete XML file for multiple bundles, 
	 * with declaration and root element.
	 * 
	 * @param out
	 * @param bundle_list
	 */
	public static void printBundles(PrintStream out, List<ModelBundle> bundle_list )
	{
		out.println(XML_DECLARATION);
		java.util.Date now = java.util.Calendar.getInstance().getTime();
		out.printf("<%s %s=\"%.4f\" %s=\"%s\">\n"
				, EMT_ROOT
				, ATT_VERSION, Count.getAppVersionNumber()
				, ATT_DATE, now.toString()
				);
		for (ModelBundle bundle: bundle_list)
		{
			out.println(bundle.toXMLString());
		}
		out.printf("</%s>\n", EMT_ROOT);
	}
	
	/**
	 * Writes a complete XML file for a single bundle, with declaration and root element.
	 * 
	 * @param out where to print
	 * @param bundle what to write into the XML file
	 */ 
	public static void printBundle(PrintStream out, ModelBundle bundle)
	{
		printBundles(out, Collections.singletonList(bundle));
//		out.println(XML_DECLARATION);
//		java.util.Date now = java.util.Calendar.getInstance().getTime();
//		out.printf("<%s %s=\"%s\" %s=\"%d\" %s=\"%s\">\n"
//				, EMT_ROOT
//				, ATT_DATE, now.toString()
//				, ATT_THREADS, Count.THREAD_PARALLELISM
//				, ATT_VERSION, Count.APP_VERSION);
//		out.println(bundle.toXMLString());
//		out.printf("</%s>\n", EMT_ROOT);
	}
	
	/**
	 * Parser for XML-encoded bundle, with downward compatibility of 2012 session save format.  
	 * 
	 * @author csuros
	 *
	 */
	private static class Handler extends DefaultHandler
	{
		
		private final Map<String, StringBuilder> element_content;
		private final Map<String, Entry> entries;
		private final Map<String, Entry> enclosing_entry;
		private Entry current; 
		private ModelBundle current_bundle;
		private final List<ModelBundle> bundle_list;
		
		private Handler()
		{
			this.element_content = new HashMap<>();
			this.current = null;
			this.bundle_list = new ArrayList<>();
			this.entries = new HashMap<>();
			this.enclosing_entry = new HashMap<>();
			this.current_bundle = null;
		}
		
        @Override
        public void startElement(String uri,
                  String localName,
                  String qName,
                  Attributes attributes)
                  throws SAXException      
        {
            String elementName="".equals(uri)?qName:localName;
            String id = attributes.getValue(ATT_ID);
            if (id==null)
            {
            	assert (EMT_ROOT.equals(elementName));
            	id = EMT_ROOT;
            }
            element_content.put(id,  new StringBuilder());
            Entry previous = current;
            
            if (EMT_ROOT.equals(elementName))
            {
            	String sthreads = attributes.getValue(CountXML.ATT_THREADS);
            	if (sthreads != null)
            		Count.THREAD_PARALLELISM = Integer.parseInt(sthreads);
//            	String struncate = attributes.getValue(CountXML.ATT_TRUNCATE);
//            	if (struncate!=null)
//            	{
//            		int absolute = CommandLine.parseTruncateAbsolute(struncate);
//            		double relative = CommandLine.parseTruncateRelative(struncate);
//            	}
            } else if (EMT_SESSION.equals(elementName))
            {
            	current_bundle = new ModelBundle(id);
            	current = current_bundle.getRoot();
            	entries.put(current.entry_id,  current);
            	bundle_list.add(current_bundle);
                enclosing_entry.put(current.entry_id,  null);
            } else if (EMT_TREE.equals(elementName))
            {
            	String filename = attributes.getValue(ATT_NAME);
            	DataFile<Phylogeny> tree_data = new DataFile<>(null, new File(filename));
//            	String dirty = attributes.getValue("dirty");
//            	tree_data.setDirty("true".equals(dirty));
            	Entry T = current.addTree(tree_data, id);
            	entries.put(T.entry_id,  T);
            	current = T;
                enclosing_entry.put(T.entry_id,  previous);
            } else if (EMT_RATES.equals(elementName))
            {
            	
            	String parent_id = attributes.getValue(ATT_PARENT);
            	Entry T;
            	if (parent_id == null)
            	{
            		T = current;
            		Entry P = T.getParent();
            		while (P!=null && T.tree_data==null) {T=P;P=T.getParent();}
            		if (T.tree_data == null)
            		{
            			assert (T.isRoot());
            			int ci=0;
            			while(ci<T.getNumChildren() && !T.getChild(ci).isTreeEntry()) ci++;
            			if (ci!=T.getNumChildren()) T=T.getChild(ci); // first tree
            			if (T.tree_data == null)
            				throw new SAXException("Cannot determine parent for <"+EMT_RATES+" id="+id+">");
            		}
            		parent_id = T.entry_id;
            	} else
            	{
            		T = entries.get(parent_id);
            		if (T==null || !T.isTreeEntry())
            		{
            			throw new SAXException("Cannot determine parent for <"+EMT_RATES+" id="+id+" parent="+parent_id+">");
            		}
            	}
            	entries.put(id, T); // store the  tree
            	current = T;
                enclosing_entry.put(id,  previous);
            } else if (EMT_RATEMODEL.equals(elementName))
            {
            	String parent_id = attributes.getValue(ATT_PARENT);
            	Entry P = entries.get(parent_id);
            	String filename = attributes.getValue(ATT_NAME);
            	DataFile<MixedRateModel> rates_data = new DataFile<>(null, new File(filename));
//            	String dirty = attributes.getValue("dirty");
//            	rates_data.setDirty("true".equals(dirty));
            	Entry R = P.addRates(rates_data, id);
            	entries.put(R.entry_id, R);
            	current = R;
                enclosing_entry.put(R.entry_id,  previous);
            } else if (EMT_DATA.equals(elementName))
            {
            	String parent_id = attributes.getValue(ATT_PARENT);
            	Entry P;
            	if (parent_id == null) P = current.getRoot();
            	else 
            	{
            		P = entries.get(parent_id);
	            	if (!P.isRoot())
	            	{
	            		throw new SAXException("Parent for <"+EMT_DATA+" id="+id+"> should be root ("+EMT_SESSION+")");
	            	}
            	}
            	entries.put(id,  P);
            	current = P;
                enclosing_entry.put(id,  previous);
            } else if (EMT_TABLE.equals(elementName))
            {
            	String parent_id = attributes.getValue(ATT_PARENT);
            	Entry P = entries.get(parent_id);
            	String filename = attributes.getValue(ATT_NAME);
            	DataFile<AnnotatedTable> table_data = new DataFile<>(null, new File(filename));
            	Entry T = P.addTable(table_data, id);
            	entries.put(T.entry_id, T);
            	current = T;
                enclosing_entry.put(T.entry_id,  previous); 	
            } else if (EMT_DATAITEM.equals(elementName) || EMT_POSTERIORS.equals(elementName))
            {
            	String parent_id = attributes.getValue(ATT_PARENT);
            	Entry P = entries.get(parent_id);
            	String filename = attributes.getValue(ATT_NAME);
            	String type;
            	if (EMT_POSTERIORS.equals(elementName)) // legacy format
            		type = PosteriorsView.class.getCanonicalName();
            	else
            		type = attributes.getValue(ATT_TYPE);
            	if (filename == null)
            	{
            		int dot = type.lastIndexOf('.');
            		if (dot == -1)
            			filename = type + " @ "+P.getTableData().getFile().getName();
            		else
            			filename = type.substring(dot+1) + " @ "+P.getTableData().getFile().getName();
            	}
            	
            	DataFile<AnnotatedTable> empty_data = new DataFile<>(null, new File(filename));
            	
            	Entry DI = P.addDataItem(empty_data, type, id);
            	entries.put(DI.entry_id, DI);
            	current = DI;
                enclosing_entry.put(DI.entry_id,  previous);
            }
            else
            {
            	throw new SAXException("Unrecognized start tag "+elementName);
            }
            
            
            if (!EMT_ROOT.equals(elementName) &&  !EMT_RATES.equals(elementName) && !EMT_DATA.equals(elementName))
            {
	            // process all extra attributes
	            for (int ai=0; ai<attributes.getLength(); ai++)
	            {
	            	String att_key = "".equals(attributes.getURI(ai))?attributes.getQName(ai):attributes.getLocalName(ai);
	            	if (!ATT_ID.equals(att_key)
	            			&& !ATT_NAME.equals(att_key)
	            			&& !ATT_PARENT.equals(att_key))
	            	{
	            		if (!current.isTableEntry() || !ATT_TYPE.equals(att_key)) // table types are set already
	            			current.setAttribute(att_key, attributes.getValue(ai));
	            	}
	            }
            }
            if (DEBUG_PARSING)
            { 
            	System.out.println("#**MB.H.start "+elementName+"\t"+current+"\tencl "+
            			(current==null?null:enclosing_entry.get(current.getId())));
            }
        }		

        @Override
        public void endElement(String uri,
                String localName,
                String qName)
                throws SAXException
        {
            String elementName="".equals(uri)?qName:localName;
            // DEBUG
//            System.out.println("#**MB.H.end "+elementName+"\tcurrent "+current);
            
            if (EMT_ROOT.equals(elementName))
            {
                if (DEBUG_PARSING)
                { 
                	System.out.println("#**MB.H.end "+elementName+"\t"+current+"\troot");
                }            	
            } else if (EMT_SESSION.equals(elementName)
            	|| EMT_TREE.equals(elementName)
            	|| EMT_RATEMODEL.equals(elementName)
            	|| EMT_TABLE.equals(elementName)
            	|| EMT_DATAITEM.equals(elementName)
            	|| EMT_POSTERIORS.equals(elementName))
            {
            	Entry backtrack = enclosing_entry.get(current.entry_id);
            	if (backtrack !=current.getParent() || DEBUG_PARSING) // parent references should reflect XML structure 
            	{
            		System.out.println("#**MB.H.end "+elementName+"\t"+current
            				+"\tenclosing "+backtrack+"\tparent "+current.getParent()
            				+"\t(parsing output from old Count?)");
            	}

            	current = backtrack;
//            	current = current.getParent();
            } else if (EMT_RATES.equals(elementName) || EMT_DATA.equals(elementName)) // accommodating legacy format
            {
                if (DEBUG_PARSING)
                { 
                	System.out.println("#**MB.H.end "+elementName+"\t"+current+"\tstay");
                }
            	// current = enclosing_entry.get(current.entry_id);
            	// nothing to do
            } else
            {
            	throw new SAXException("Unrecognized end tag "+elementName);
            }
        }
        
        @Override
        public void characters(char[] ch,
                int start,
                int length)
                throws SAXException
        {     
            // DEBUG
//        	System.out.println("#**MB.H.char "+start+"\t"+length+"\tcurrent "+current);
        	StringBuilder content;
        	if (current == null)
        	{
        		 // root element containing the sessions
        		content = element_content.get(EMT_ROOT); // there is no info stored here 
        	} else 
        	{
        		content = element_content.get(current.entry_id);
        	}
        	content.append(ch, start, length);
        }
        
        void fillContent() throws NewickParser.ParseException,RateVariationParser.FileFormatException
        {
        	try
        	{
        		for (ModelBundle bundle: bundle_list)
        		{
        			fillContent(bundle.getRoot());
        		}
        	} catch (NewickParser.ParseException badparse)
        	{
        		throw badparse;
        	} catch (RateVariationParser.FileFormatException badparse)
        	{
        		throw badparse;
        	} catch (IOException notreally)
        	{
        		// BuferedReaders over StringReaders
        		// never happens 
        	}
        }
        
        /**
         * Adds the data in the DataFile fields. 
         * 
         * @param E
         * @throws IOException
         */
        void fillContent(Entry E) throws IOException
        {
        	StringBuilder sb = element_content.get(E.entry_id);
        	if (E.isTreeEntry())
        	{
        		StringReader R = new StringReader(sb.toString());
        		Phylogeny P = NewickParser.readTree(R);
        		P.fixZeroEdges();
        		
        		E.tree_data.setData(P);
//        		System.out.println("#**MB.H.fillC tree leaves "+Arrays.toString(P.getLeafNames()));
        		
        	} else if (E.isRatesEntry())
        	{
        		StringReader R = new StringReader(sb.toString());
        		Entry T = E.getParent();
        		while (T!=null && !T.isTreeEntry()) T=T.getParent();
        		Phylogeny P = T.tree_data.getContent();
        		MixedRateModel model = RateVariationParser.readModel(new BufferedReader(R), P);
        		E.rates_data.setData(model);
//        		E.rates_data.setData(model);
        	} else if (E.isTableEntry())
        	{
        		String table_type = E.getAttributeValue(ATT_TYPE);
        		if (table_type == null)
        		{
        			StringReader R = new StringReader(sb.toString());
        			Phylogeny main_phylo = current_bundle.getMainTree().tree_data.getContent();
//        			System.out.println("#**MB.H.fillC leaves "+Arrays.toString(main_phylo.getLeafNames()));
        			AnnotatedTable table = TableParser.readTable(main_phylo.getLeafNames(), R, true);
        			E.table_data.setData(table);
        		} else
        		{
        			// content ignored for other types  
        		}
        	}
        	int nc = E.getNumChildren();
        	for (int ci=0; ci<nc; ci++)
        	{
        		Entry C = E.getChild(ci);
        		fillContent(C);
        	}
        }
	} // Handler class
	
	public static void main(String[] args) throws Exception
	{
		Class<?> our_class = ModelBundle.class;
	    System.out.println(CommandLine.getStandardHeader(our_class));
	    System.out.println(CommandLine.getStandardRuntimeInfo(our_class, args));

	    CommandLine cli = new CommandLine(args,our_class, 0) ;
	    
	    String out_file = cli.getOptionValue(CommandLine.OPT_OUTPUT);
		PrintStream out = System.out;
	    if (out_file != null && !"-".equals(out_file))
	    {
	    	out = new PrintStream(out_file);
	    }
	    
	    DataFile<Phylogeny> tree_data = cli.getTreeData();
	    
	    if (cli.getOptionValue(CommandLine.OPT_LOAD)!=null)
	    {
	    	String session_file = cli.getOptionValue(CommandLine.OPT_LOAD);
	    	BufferedReader BR = GeneralizedFileReader.guessBufferedReaderForInput(session_file);
	    	List<ModelBundle> input_bundles = ModelBundle.readBundle(BR);
	    	BR.close();
	    	out.println(XML_DECLARATION);
	    	for (ModelBundle bundle:input_bundles)
	    	{
	    		out.println(bundle.toXMLString());
	    	}
	    } else if(tree_data == null)
	    {
	    	throw new IllegalArgumentException("Specify input tree or -"+CommandLine.OPT_LOAD);
	    } else
	    {
		    // output test: generate some random models

	    	Phylogeny tree = tree_data.getContent();
		    String session_name = DataFile.chopFileExtension(tree_data.getFile().getName());
		    {
		    	int slash = session_name.lastIndexOf('/');
		    	if (slash >= 0)
		    	{
		    		session_name = session_name.substring(slash+1, session_name.length());
		    	}
		    }

		    ModelBundle bundle = new ModelBundle(session_name);
		    String[] taxa = tree.getLeafNames();
		    java.util.Random RND = new java.util.Random(2023);
		    
		    int num_models = 3;
		    for (int mi=0; mi<num_models; mi++)
		    {
		    	Phylogeny mtree;
		    	if (mi==0)
		    		mtree = tree;
		    	else
		    	{
		    		mtree = Phylogeny.randomTree(taxa, RND, true);
		    	}
		    	GammaInvariant model = GammaInvariant.nullModel(mtree, RND);
		    	DataFile<Phylogeny> mdata = new DataFile<>(mtree, new File(tree_data.getFile().getParent(), session_name+"."+Integer.toString(mi)+".tre"));
		    	DataFile<MixedRateModel> rates_data = new DataFile<>(model, new File(tree_data.getFile().getParent(), session_name+"."+Integer.toString(mi)+".rates.txt"));
		    	bundle.add(mdata, rates_data);
		    }
		    String bundle_str = bundle.toXMLString();
		    out.println(XML_DECLARATION);
		    out.println(bundle_str);
	    }
	    out.close();
	}
	

}
