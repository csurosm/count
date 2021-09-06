package ca.umontreal.iro.evolution;

/**
 * Title:        Triplet
 * Description:  Triplet of external nodes with orientation
 * Copyright:    Copyright (c) 2001
 * Company:
 * @author      Miklos Csuros csuros AT iro.umontreal.ca
 * @version 1.0
 */

public class Triplet {

  TreeLeafNode[] node;
  private double[] distance;
  private double average_similarity;

  public static final int TOP_DIRECTION=0;
  public static final int LEFT_DIRECTION=1;
  public static final int RIGHT_DIRECTION=2;

  public Triplet(TreeLeafNode top, TreeLeafNode left, TreeLeafNode right)
  {
    node=new TreeLeafNode[3];
    node[TOP_DIRECTION]=top;
    node[LEFT_DIRECTION]=left;
    node[RIGHT_DIRECTION]=right;

    distance=new double[3];

    double distance_top_left=node[TOP_DIRECTION].getDistance(node[LEFT_DIRECTION]);
    double distance_top_right=node[RIGHT_DIRECTION].getDistance(node[TOP_DIRECTION]);
    double distance_left_right=node[LEFT_DIRECTION].getDistance(node[RIGHT_DIRECTION]);

    if (distance_left_right != Double.POSITIVE_INFINITY
      && distance_top_left>=0. && distance_top_right >= 0.)
      distance[TOP_DIRECTION]=
        (distance_top_left+distance_top_right-distance_left_right)/2.0;
    else
      distance[TOP_DIRECTION]=Double.NaN;

    if (distance_top_right != Double.POSITIVE_INFINITY
      && distance_top_left>=0. && distance_left_right >= 0.)
      distance[LEFT_DIRECTION]=
        (distance_top_left-distance_top_right+distance_left_right)/2.0;
    else
      distance[LEFT_DIRECTION]=Double.NaN;

    if (distance_top_left != Double.POSITIVE_INFINITY
      && distance_left_right>=0. && distance_top_right >= 0.)
      distance[RIGHT_DIRECTION]=
        (-distance_top_left+distance_top_right+distance_left_right)/2.0;
    else
      distance[RIGHT_DIRECTION]=Double.NaN;


    if (distance_top_left<=0. || distance_top_right<=0. || distance_left_right<=0.)
      average_similarity = 0.;
    else
    {
      double c0=Math.exp(distance_top_left);
      double c1=Math.exp(distance_top_right);
      double c2=Math.exp(distance_left_right);
      average_similarity = 3.0/(c0+c1+c2);
    }

  }

  public TreeLeafNode getTop(){return getLeaf(TOP_DIRECTION);}
  public TreeLeafNode getLeft(){return getLeaf(LEFT_DIRECTION);}
  public TreeLeafNode getRight(){return getLeaf(RIGHT_DIRECTION);}
  public TreeLeafNode getLeaf(int direction){return node[direction];}

  public double getDistanceFromCenterTotop(){return getDistanceFromCenter(TOP_DIRECTION);}
  public double getDistanceFromCenterToLeft(){return getDistanceFromCenter(LEFT_DIRECTION);}
  public double getDistanceFromCenterToRight(){return getDistanceFromCenter(RIGHT_DIRECTION);}
  public double getDistanceFromCenter(int direction){return distance[direction];}

  public double getDistanceFromCenter(TreeLeafNode n)
  {
    for (int i=2; i>=0 ; i--)
      if (node[i].equals(n))
        return getDistanceFromCenter(i);
    return Double.NaN;
  }

  public double getAverageSimilarity(){return average_similarity;}

  /**
   * can be a wrapper for other scoring methods
   */
  public double getScore(){return getAverageSimilarity();}

  public String toString()
  {
    StringBuffer sb=new StringBuffer(getClass().getName());
    sb.append("[score: ");
    sb.append(getScore());
    sb.append(" top: ");
    sb.append(node[TOP_DIRECTION].toString());
    sb.append(", left: ");
    sb.append(node[LEFT_DIRECTION].toString());
    sb.append(", right: ");
    sb.append(node[RIGHT_DIRECTION].toString());
    sb.append(']');
    return sb.toString();
  }

  public boolean equals(Triplet t)
  {
    return t.node[TOP_DIRECTION].equals(node[TOP_DIRECTION])
        && t.node[LEFT_DIRECTION].equals(node[LEFT_DIRECTION])
        && t.node[RIGHT_DIRECTION].equals(node[RIGHT_DIRECTION]);
  }


  static String RCS_ID="$Id: Triplet.java,v 1.1 2001/03/09 17:31:04 mcsuros Exp mcsuros $";

  //
  // $Log: Triplet.java,v $
  // Revision 1.1  2001/03/09 17:31:04  mcsuros
  // Initial revision
  //
  //
}