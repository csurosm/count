package ca.umontreal.iro.evolution;

/**
 * Title:
 * Description:
 * Copyright:    Copyright (c) 2001
 * Company:
 * @author
 * @version 1.0
 */

public class RelevantTuple implements Comparable {
  Triplet triplet;
  TreeNode edge_lower_end;

  int direction_for_shared_with_lower;
  int direction_for_shared_with_upper;
  int direction_for_new;

  // the edge lengths if this guy is used for inserting new nodes
  double edge_length_upper;
  double edge_length_lower;
  double distance_to_new_leaf;
  private boolean distances_calculated=false;

  public RelevantTuple(
      TreeNode edge_lower_end,
      TreeLeafNode leaf_shared_with_lower,
      TreeLeafNode leaf_shared_with_upper,
      TreeLeafNode new_leaf) {
    this.edge_lower_end=edge_lower_end;
    initTriplet(leaf_shared_with_lower, leaf_shared_with_upper, new_leaf);
  }

  public Triplet getTriplet(){return triplet;}
  public TreeNode getEdgeLowerEnd(){return edge_lower_end;}
  public TreeLeafNode getNewLeaf(){return triplet.getLeaf(direction_for_new);}

  // orientations for top,left,right on the triplet within initTriplet()
  // defining the 4 possible permutations for
  // shared_with_upper, shared_with_lower, new_leaf in a relevant tuple
  public static final int ORIENTATION_ULN=0;
  public static final int ORIENTATION_UNL=1;
  public static final int ORIENTATION_LNU=2;
  public static final int ORIENTATION_LUN=3;

  public RelevantTuple(int orientation,
      TreeNode edge_lower_end,
      TreeLeafNode leaf_shared_with_lower,
      TreeLeafNode leaf_shared_with_upper,
      TreeLeafNode new_leaf) {
    this.edge_lower_end=edge_lower_end;
    initTriplet(orientation, leaf_shared_with_lower, leaf_shared_with_upper, new_leaf);
  }


  public double getEdgeLengthFromLowerEnd()
  {
    if (!distances_calculated)
      calculateDistances();
    return edge_length_lower;
  }

  public double getEdgeLengthFromUpperEnd()
  {
    if (!distances_calculated)
      calculateDistances();
    return edge_length_upper;
  }

  public double getDistanceFromNewLeaf()
  {
    if (!distances_calculated)
      calculateDistances();
    return distance_to_new_leaf;
  }

  private void initTriplet(
      TreeLeafNode leaf_shared_with_lower,
      TreeLeafNode leaf_shared_with_upper,
      TreeLeafNode new_leaf)
  {
    boolean flipped=false; // whether shared_with_lower_end should be on the top in the triplet


    if (! edge_lower_end.isLeaf()) // otherwise shared_with_lower == edge_lower_end
    {
      HGTNode edge_upper_end=(HGTNode)edge_lower_end.getParent();

      if (!edge_upper_end.isExternal()) // otherwise shared_with_upper == edge_upper_end
      {
        // both upper and lower endpoints are internal nodes
        flipped=(leaf_shared_with_lower.equals(
          ((HGTNode) edge_lower_end).getDefiningTriplet().getTop()));
      }
    }

    int orientation=-1;
    if (flipped)
    {
      if (edge_lower_end.isOnTheLeft())
        orientation=ORIENTATION_LUN;
      else
        orientation=ORIENTATION_LNU;
    } else
    {
      if (edge_lower_end.isOnTheLeft())
        orientation=ORIENTATION_ULN;
      else
        orientation=ORIENTATION_UNL;
    }

    initTriplet(orientation, leaf_shared_with_lower, leaf_shared_with_upper, new_leaf);
  }

  private void initTriplet(
    int orientation,
      TreeLeafNode leaf_shared_with_lower,
      TreeLeafNode leaf_shared_with_upper,
      TreeLeafNode new_leaf)
  {

    if (orientation == ORIENTATION_ULN)
    {
      triplet = new Triplet(
          leaf_shared_with_upper,
          leaf_shared_with_lower,
          new_leaf);
      direction_for_shared_with_lower=Triplet.LEFT_DIRECTION;
      direction_for_shared_with_upper=Triplet.TOP_DIRECTION;
      direction_for_new              =Triplet.RIGHT_DIRECTION;
    } else if (orientation == ORIENTATION_UNL)
    {
      triplet = new Triplet(
        leaf_shared_with_upper,
        new_leaf,
        leaf_shared_with_lower);
      direction_for_shared_with_lower=Triplet.RIGHT_DIRECTION;
      direction_for_shared_with_upper=Triplet.TOP_DIRECTION;
      direction_for_new              =Triplet.LEFT_DIRECTION;
    } else if (orientation == ORIENTATION_LNU)
    {
      triplet = new Triplet(
        leaf_shared_with_lower,
        new_leaf,
        leaf_shared_with_upper);
      direction_for_shared_with_lower=Triplet.TOP_DIRECTION;
      direction_for_shared_with_upper=Triplet.RIGHT_DIRECTION;
      direction_for_new              =Triplet.LEFT_DIRECTION;
    } else if (orientation == ORIENTATION_LUN)
    {
      triplet = new Triplet(
        leaf_shared_with_lower,
        leaf_shared_with_upper,
        new_leaf);
      direction_for_shared_with_lower=Triplet.TOP_DIRECTION;
      direction_for_shared_with_upper=Triplet.LEFT_DIRECTION;
      direction_for_new              =Triplet.RIGHT_DIRECTION;
    }
  }


  private void calculateDistances()
  {
    // this is for sure
    distance_to_new_leaf = triplet.getDistanceFromCenter(direction_for_new);

    // the triplet center's distance from the shared leaves
    double center_lower_distance=
      triplet.getDistanceFromCenter(direction_for_shared_with_lower);
    double center_upper_distance=
      triplet.getDistanceFromCenter(direction_for_shared_with_upper);

    // lower edge endpoint's distance from the shared leaf
    double edge_lower_distance=0.;
    if (!edge_lower_end.isLeaf())
    {
      edge_lower_distance =
        ((HGTNode)edge_lower_end).getDefiningTriplet().getDistanceFromCenter(
          triplet.getLeaf(direction_for_shared_with_lower));
    }
    // signed edge length for edge after insertion
    double d_lower=center_lower_distance-edge_lower_distance;
    if (direction_for_shared_with_lower == Triplet.TOP_DIRECTION)
      d_lower = -d_lower;

    // upper edge endpoint's distance from the shared leaf
    TreeNode edge_upper_end=edge_lower_end.getParent();
    double edge_upper_distance=0.;
    if (!edge_upper_end.isExternal())
    {
      edge_upper_distance =
        ((HGTNode)edge_upper_end).getDefiningTriplet().getDistanceFromCenter(
          triplet.getLeaf(direction_for_shared_with_upper));
    }
    double d_upper=center_upper_distance-edge_upper_distance;
    if (direction_for_shared_with_lower == Triplet.TOP_DIRECTION)
      d_upper=-d_upper;

    edge_length_lower=(d_lower+edge_lower_end.getLength()-d_upper)/2.0;
    edge_length_upper=(d_upper+edge_lower_end.getLength()-d_lower)/2.0;

    distances_calculated=true;
  }

  public double getScore()
  {
    return triplet.getScore();
  }

  public int compareTo(Object o)
  {
    double s0=getScore();

    if (o==null) // 0 score is just as bad as nothing at all
      return (s0<=0.?0:1);

    double s1=((RelevantTuple)o).getScore();
    if (s0>s1)
      return 1;
    else if (s0<s1)
      return -1;
    else
      return 0;
  }

  public boolean areEdgeLengthsOK()
  {
    return true;
//
// maybe want to include this to avoid negative edge lengths
//
//    return (getEdgeLengthFromLowerEnd()>0.
//            && getEdgeLengthFromUpperEnd()>0.);
  }

  public boolean isBetterThan(RelevantTuple r)
  {
    return (compareTo(r)>0);
  }

  public String toString()
  {
    StringBuffer sb=new StringBuffer(getClass().getName());
    sb.append("[");
    sb.append(" triplet: ");
    sb.append(triplet.toString());
    sb.append(" edge:");
    sb.append(edge_lower_end.toString());
    sb.append(']');
    return sb.toString();
  }

  static String RCS_ID="$Id: RelevantTuple.java,v 1.1 2001/03/09 17:31:04 mcsuros Exp mcsuros $";
  //
  // $Log: RelevantTuple.java,v $
  // Revision 1.1  2001/03/09 17:31:04  mcsuros
  // Initial revision
  //
  //

}