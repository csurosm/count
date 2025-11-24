package count.gui;
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


import java.awt.Color;

import javax.swing.JTable;
//import javax.swing.event.TableModelEvent;
import javax.swing.table.AbstractTableModel;

import count.ds.AnnotatedTable;
import count.ds.IntegerOrMissing;
import count.gui.kit.RoundedDouble;

/**
 * Model for family table content. Columns are the following
 * <ul> 
 * <li> 1 family name, 
 * <li> 2... <var>x</var> properties
 * <li> <var>x</var>+1 number of lineages in which it is present
 * <li> <var>x</var>+2 total number of members
 * <li> <var>x</var>+3 .. members in each lineage (if detailed profiles)
 * <li> <var>x</var>+3 profile summary
 * </ul>
 * 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 */
public class AnnotatedTableModel extends AbstractTableModel
{
    /**
     * Instantiation of the model
     * 
     * @param data_table the underlying family table
     * @param want_detailed_profiles whether there should be a separate column for each organism, or just a 
     *          single column for profile summary
     */
    public AnnotatedTableModel(AnnotatedTable data_table, boolean want_detailed_profiles)
    {
        super();
        this.data_table = data_table;
        this.want_detailed_profiles = want_detailed_profiles;
        this.leaves = data_table.getTaxonNames();
        initDataStructures();
    }
    
    /**
     * The underlying data.
     */
    private final AnnotatedTable data_table;
    
    public AnnotatedTable getTable() { return data_table;}
    
    private final boolean want_detailed_profiles;
    
    private AnnotatedTable.PhyleticProfile[] occurrence_data;

    /**
     * Names of terminal taxa
     */
    private final String[] leaves;

    /**
     *  Number of lineages in which a famiy is present
     */
    private int[] present_leaves;
    
    /**
     * Number of genes in the family
     */
    private int[] member_count;
    
    private void initDataStructures()
    {
        int num_families = data_table.getFamilyCount();

        present_leaves=new int[num_families];
        member_count = new int[num_families];
        for(int family_idx=0; family_idx<num_families; family_idx++)
        {
        	present_leaves[family_idx] = data_table.getLineageCount(family_idx);
        	member_count[family_idx] = data_table.getMemberCount(family_idx);
//            int[] pattern = data_table.getFamilyProfile(family_idx);
//            present_leaves[family_idx]=0;
//            member_count[family_idx]=0;
//            for (int leaf_idx=0; leaf_idx<leaves.length; leaf_idx++)
//            {
//                member_count[family_idx]+=Math.max(0,pattern[leaf_idx]);
//                if (pattern[leaf_idx]>0)
//                {
//                    present_leaves[family_idx]++;
//                }
//            }
        }
        occurrence_data = new AnnotatedTable.PhyleticProfile[num_families]; // will be filled in during scrolling as necessary        
        this.addTableModelListener(e->
					        {
					        	if (e.getType() == e.UPDATE)
					        	{
					        		for (int f=e.getFirstRow(); f<=e.getLastRow(); f++)
					        		{
					        			if (f != e.HEADER_ROW)
			        					{	
					        				updateRowData(f);
			        					}
					        		}
					        	}
					        });
    }
    
    private void updateRowData(int f)
    {
    	present_leaves[f] = data_table.getLineageCount(f);
    	member_count[f] = data_table.getMemberCount(f);
    	occurrence_data[f] = data_table.getPhyleticProfile(f);
    	
//    	System.out.println("#**ATM.uRD "+f);
    }
//    protected AnnotatedTable getDataTable()
//    {
//    	return data_table;
//    }

    @Override
    public int getColumnCount()
    {
        int ncol = 1 // row number
                +data_table.getKnownPropertiesCount()
                +2 // lineage and member totals
                +(want_detailed_profiles?leaves.length:1); // profile
        //System.out.println("#*OTM.gCC "+ncol);
        return ncol;
    }
    
    @Override
    public int getRowCount()
    {
        return data_table.getFamilyCount();
    }
    
    public int firstProfileColumn()
    {
    	return 1+data_table.getKnownPropertiesCount()+2;
    }
    
    public double getMaximumValue(int column)
    {
    	if (getColumnClass(column).isAssignableFrom(Number.class))
    	{
    		double max = Double.NEGATIVE_INFINITY;
    		for (int row=0; row<getRowCount(); row++)
    		{
    			Number val = (Number) getValueAt(row, column);
    			max = Double.max(max, val.doubleValue());
    		}
    		return max;
    	} else 
    		return Double.NaN;
    }
    
    /**
     * Cell content. The underlying data is queried to copy the values 
     * into a local data structure as it becomes necessary.
     * 
     * @param row_idx row index
     * @param column_idx column index
     * @return cell content value
     */
    @Override
    public Object getValueAt(int row_idx, int column_idx)
    {
        if (column_idx==0)
            return Integer.valueOf(1+row_idx);
        column_idx--;

        int num_props = data_table.getKnownPropertiesCount();
        if (column_idx<num_props)
        {
            return data_table.getFamilyProperty(row_idx, column_idx);
        } else
            column_idx -= num_props;

        if (column_idx==0)
            return Integer.valueOf(present_leaves[row_idx]);
        column_idx--;
        if (column_idx==0)
            return Integer.valueOf(member_count[row_idx]);
        column_idx--;

        if (occurrence_data[row_idx]==null)
            occurrence_data[row_idx] = data_table.getPhyleticProfile(row_idx);
        
        Object retval = null;
        if (want_detailed_profiles)
        {
//            if (occurrence_data[row_idx][column_idx]>=0)
                retval = occurrence_data[row_idx].getValue(column_idx);
//            else
//                retval = "?";
        }
        else
            retval = occurrence_data[row_idx];        
        //System.out.println("#*OTM.gVA row "+row_idx+" col "+column_idx+" ret "+retval+"\tclass "+retval.getClass());
        return retval;
    }    
    
    
    @Override
    public String getColumnName(int column_idx)
    {
        //System.out.println("#*FSTD.OTM.gCN col "+column_idx);
        if (column_idx==0)
            return "Order";
        column_idx--;

        int num_props = data_table.getKnownPropertiesCount();
        if (column_idx<num_props)
        {
            return data_table.getPropertyName(column_idx);
        } else
            column_idx -= num_props;

        if (column_idx==0)
            return "#lin";
        else
            column_idx--;
        if (column_idx==0)
            return "#mem";
        else
            column_idx--;
        
        if (want_detailed_profiles)
        {
            return leaves[column_idx]+":n";
//            return "\u03a6"+leaves[column_idx];
        }
        else
            return "Profile";
    }

    /**
     * Specify column classes to use the correct comparator in sorting
     * 
     * @param column_idx index of the column
     * @return what class the colum belongs to
     */
    @Override
    public Class<?> getColumnClass(int column_idx)
    {
        if (column_idx == 0) return Integer.class; else column_idx --;

        int num_props = data_table.getKnownPropertiesCount();
        if (column_idx<num_props)
            return String.class;
        else
            column_idx -= num_props;
        if (column_idx<2)
            return Integer.class;
        else
            column_idx -= 2;

        if (want_detailed_profiles)
            return IntegerOrMissing.class;
        else
            return AnnotatedTable.PhyleticProfile.class;
    }
    
    public boolean hasDetailedProfiles()
    {
    	return want_detailed_profiles;
    }
    
    /**
     * Description for the column
     * 
     * @param column_idx
     * @return 
     */
    public String getColumnDescription(int column_idx)
    {
        if (column_idx==0)
            return "Family index";
        else column_idx--;
        int num_props = data_table.getKnownPropertiesCount();
        if (column_idx<num_props)
        {
            return data_table.getPropertyName(column_idx);
        } else
            column_idx -= num_props;

        if (column_idx==0)
            return "Number of terminal lineages in which the family has at least one member";
        column_idx--;
        if (column_idx==0)
            return "Total number of homologs, i.e., sum of family size across the terminal lineages";
        column_idx--;
        return "Profile: number of homologs at terminal node"+(want_detailed_profiles?" "+leaves[column_idx]:"s");
    }
    
    /**
     * Can be called for producing a proper tool tip.
     * 
     * @param row_idx row index by model
     * @param column_idx column index by model
     * @return a tool tip for this cell
     */
    public String getCellToolTip(int row_idx, int column_idx)
    {
        Object val = getValueAt(row_idx, column_idx);
        if (column_idx==0)
            return "Family order is "+val;
        else column_idx--;
        
        int num_props = data_table.getKnownPropertiesCount();
        if (column_idx<num_props)
        {
            return data_table.getPropertyName(column_idx)+": "+val;
        } else
            column_idx -= num_props;
        
        if (column_idx==0)
            return "Family "+data_table.getFamilyName(row_idx)+" has at least one member in "+val+" terminal lineages.";
        else
            column_idx--;
        if (column_idx==0)
            return "Family "+data_table.getFamilyName(row_idx)+" has a total of "+val+" members across the terminal lineages";
        else
            column_idx--;
        
        if (want_detailed_profiles)
        {
            IntegerOrMissing ival = (IntegerOrMissing)val;
            
            if (ival.isAmbiguous())
            {
                return "Nobody knows how many members the family "+data_table.getFamilyName(row_idx)+" has at terminal node "+leaves[column_idx];
            } else
            {
                int c = ival.intValue();

                return "Family "+data_table.getFamilyName(row_idx)+" has "+(c==0?"no":Integer.toString(c))
                        +" member"+(c==1?"":"s")+" at terminal node "+leaves[column_idx];
            }
        } else
        {
            StringBuilder sb = new StringBuilder("Family ");
            sb.append(data_table.getFamilyName(row_idx));
            AnnotatedTable.PhyleticProfile P = (AnnotatedTable.PhyleticProfile)val;
            sb.append(":");
            sb.append(P.getPatternString());

            return sb.toString();
        }
    }
    
    
    
    /**
     * Sets the renderer for profile summaries. The 
     * argument is supposed to use this table model.
     * 
     * @param tbl a JTable object in which the renderers need to be set
     */
    public void setDefaultRenderers(JTable tbl)
    {
        if (tbl.getModel() != this)
            throw new IllegalArgumentException(getClass().getName()+".setDefaultRenderer should be used with a table that has the same model.");
        tbl.setDefaultRenderer(AnnotatedTable.PhyleticProfile.class, new PhyleticProfileRenderer(leaves.length));
        tbl.setDefaultRenderer(RoundedDouble.class, new RoundedDouble.Renderer())  ;      
        tbl.setDefaultRenderer(IntegerOrMissing.class, new IntegerOrMissing.Renderer());
    }
    
    protected void setPhyleticProfileColoring(JTable tbl, Color[] leaf_colors)
    {
        if (tbl.getModel() != this)
            throw new IllegalArgumentException(getClass().getName()+".setPhyleticProfileColoring should be used with a table that has the same model.");
    	
        if (tbl.getDefaultRenderer(AnnotatedTable.PhyleticProfile.class) instanceof PhyleticProfileRenderer)
        {
	        PhyleticProfileRenderer renderer  = (PhyleticProfileRenderer) tbl.getDefaultRenderer(AnnotatedTable.PhyleticProfile.class);
	        renderer.setLeafColors(leaf_colors);
        } else
        {
//        	System.out.println("#**ATM.sPPC no PhyleticProfileRenderer: table "+tbl+"\trenderer "+tbl.getDefaultRenderer(AnnotatedTable.PhyleticProfile.class));
        }
    }
    
}
