package ca.umontreal.iro.evolution.malin.ui.count;

/**
 *
 * @author csuros
 */

import java.awt.Color;

import ca.umontreal.iro.evolution.genecontent.NodeWithRates;
import ca.umontreal.iro.evolution.genecontent.TreeWithRates;
import ca.umontreal.iro.evolution.malin.ui.ColoredValueRenderer;

public class LookAndFeel 
{
    public static Color MULTI_PRESENCE_COLOR = new Color(128,0,255);//Color.RED;
    public static Color ABSENCE_COLOR  = Color.WHITE;
    public static Color SINGLE_PRESENCE_COLOR = ColoredValueRenderer.intermediateColor(ABSENCE_COLOR, MULTI_PRESENCE_COLOR, 0.5); 
    
    public static Color LOSS_COLOR = new Color(255,128,0); // Tangerine
    public static Color GAIN_COLOR = new Color(64,128,0); // Fern
    public static Color DUPLICATION_COLOR = new Color(104,118,231); // Evening Blue
    
    public static Color EXPANSION_COLOR = new Color(128,128,0);
    public static Color REDUCTION_COLOR = new Color(255,255,0);
    
    public static int TABLE_FONT_SIZE = 14;
    public static int TREE_PANEL_FONT_SIZE = 12;

    public static final Color SMOKY_BACKGROUND = new Color(120,120,200,50);
    
    public static final Color DISTRIBUTION_PLOT_BACKGROUND = new Color(204, 255, 102); // Honeydew // new Color(255,204,102); // Canteloupe
    public static final Color RATE_VARIATION_DISTRIBUTION_COLOR = new Color(204, 255, 102); // Honeydew
    public static final Color RATE_VARIATION_CATEGORY_COLOR = new Color(64,128,0); // Fern
 
    public static final Color TREE_PANEL_SELECTED_NODE_INFO_BACKGROUND = new Color(255,220,220,240);

    public static int TREE_LEGEND_INSET = 4;
    public static int TREE_LEGEND_BAR_HEIGHT = 80;
    
    /**
     * Constructs a short name for the node.
     * 
     * @param tree the tree topology with a fixed node traversal order
     * @param idx index of tree node (into the traversal array)
     * @return leaf name if leaf, or a positive integer for non-leaf nodes
     */
    public static String getShortNodeName(TreeWithRates tree, int idx)
    {
        NodeWithRates N = tree.getNode(idx);
        if (N.isLeaf())
            return N.getTaxonName();
        else
            return Integer.toString(idx-tree.getNumLeaves()+1);
    }
    
    /**
     * Constructs a long name for the node.
     * 
     * @param tree the tree topology with a fixed node traversal order
     * @param idx index of tree node (into the traversal array)
     * @return leaf name if leaf, or String of style "int [name]" for non-leaf nodes
     */
    public static String getLongNodeName(TreeWithRates tree, int idx)
    {
        String short_name = getShortNodeName(tree, idx);
        NodeWithRates N = tree.getNode(idx);
        if (N.isLeaf())
            return short_name;
        else
        {
            return short_name+" ["+N.newickName()+"]";
        }
    }
    
    
}
