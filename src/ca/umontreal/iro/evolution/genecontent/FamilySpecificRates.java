
package ca.umontreal.iro.evolution.genecontent;

/**
 *
 * @author csuros
 */
public class FamilySpecificRates
{
    public FamilySpecificRates(TreeWithRates main_tree, OccurrenceTable table)
    {
        setTree(main_tree);
        setTable(table);
    }
    
    private OccurrenceTable[] table_columns;
    private TreeWithRates main_tree;
    private TreeWithRates scaled_tree;
    private RateVariation[] family_rates;

    private void setTable(OccurrenceTable table)
    {
        this.table_columns = table.allTablesForFamilies();
    }

    private void setTree(TreeWithRates tree)
    {
        this.main_tree = tree;
        scaled_tree = new TreeWithRates(NodeWithRates.copyTree(main_tree.getRoot()));
    }
    
    private void initRates()
    {
        int num_families = table_columns.length;
        family_rates = new RateVariation[num_families];
        for (int family_idx=0; family_idx<num_families; family_idx++)
        {
            RateVariation rates = new RateVariation(scaled_tree, null, 1, 1, 1, 1);
        }
        
    }


}
