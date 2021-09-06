package ca.umontreal.iro.evolution;

/**
 * Title:        TreeLafNode
 * Description:  Has a reference to a distance matrix so that distances are easily accessible
 * Copyright:    Copyright (c) 2001
 * Company:
 * @author
 * @version 1.0
 */

public class TreeLeafNode extends TreeNode  {
  int index;
  DistanceMatrix dm;

  public TreeLeafNode(int index, DistanceMatrix dm) {
    this.index=index;
    this.dm=dm;
  }

  public double getDistance(TreeLeafNode node_from)
  {
    return dm.getDistance(index,node_from.index);
  }

  public double getSimilarity(TreeLeafNode node_to)
  {
    return dm.getSimilarity(index,node_to.index);
  }

  public int getIndex(){return index;}
  public void setIndex(int index){this.index=index;}

  boolean is_in_tree=false;

  public boolean isInTree(){return is_in_tree;}
  public void setIsInTree(boolean b){is_in_tree=b;}

  public boolean equals(TreeLeafNode o)
  {
    if (o==null)
      return false;
    return (index==o.index);
  }

  public static boolean isQuartetOK(
    TreeLeafNode i,
    TreeLeafNode j,
    TreeLeafNode k,
    TreeLeafNode l)
  {
    return i.dm.isQuartetOK(i.getIndex(),j.getIndex(),k.getIndex(),l.getIndex());
  }

  static String RCS_ID="$Id: TreeLeafNode.java,v 1.1 2001/03/09 17:31:04 mcsuros Exp mcsuros $";
  //
  // $Log: TreeLeafNode.java,v $
  // Revision 1.1  2001/03/09 17:31:04  mcsuros
  // Initial revision
  //
  //

}