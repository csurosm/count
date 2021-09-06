package ca.umontreal.iro.evolution;

/**
 * Title:
 * Description:
 * Copyright:    Copyright (c) 2001
 * Company:
 * @author
 * @version 1.0
 */

import java.io.BufferedReader;
import java.io.IOException;

import java.util.Vector;

public class Algorithm {
  public Algorithm(){};

  DistanceMatrix distance;
  TreeLeafNode[] external_node; // array of external nodes
  int num_external_nodes;
  RelevantTuple[] relevant;
  int num_missing_leaves;

  private TreeNode root;

  public TreeNode getRoot(){return root;}

  /**
   * @return "Neighbor Joining-flavored" topology (w/ 3 children at root)
   *
   */
  public TreeNode getRootDegree3()
  {
    TreeNode new_root = root.getChild(0);

    double l=new_root.getLength();
    new_root.setParent(null);
    new_root.setLength(0.);

    root.setLength(l);
    root.clearChildren();

    new_root.addChild(new_root.getChild(1));
    new_root.setChild(new_root.getChild(0),1);
    new_root.setChild(root,0);

    return new_root;
  }

  public int getNumMissingNodes(){return num_missing_leaves;}
  public int getNumExternalNodes(){return num_external_nodes;}
  public TreeLeafNode getExternalNode(int index){return external_node[index];}
  public TreeLeafNode[] getExternalNodes(){return external_node;}

  public boolean done(){return num_missing_leaves==0;}

  /**
   * initializes the tree as a star
   */
  public HGTNode init(int first_node_index)
    throws HGTException
  {
    HGTNode first_internal_node
      =initTreeWithTriplet(getStartingTriplet(first_node_index));
    root=first_internal_node.getParent();

    // allocate space for relevant tuples
    relevant = new RelevantTuple[num_external_nodes];

    // add all initial relevant tuples
    addNewRelevantTuples(first_internal_node);
    addNewRelevantTuples(first_internal_node.getChild(TreeNode.LEFT_CHILD));
    addNewRelevantTuples(first_internal_node.getChild(TreeNode.RIGHT_CHILD));

    num_missing_leaves = num_external_nodes-3;

    return first_internal_node;
  }

  public void iterationStep() throws HGTException
  {
    if (num_missing_leaves==0)
      return;

    //System.out.println("to go: "+num_missing_leaves);

    // find best relevant tuple
    RelevantTuple r_best=getBestTuple();
    if (r_best == null)
      throw new HGTException(2, "Cannot find suitable triplet ["+num_missing_leaves+" leaves are missing");

    // add new node using the best relevant tuple
    HGTNode new_node = growTree(r_best);

    // delete relevant tuples containing the deleted edge
    deleteRelevantTuplesForEdge(r_best.getEdgeLowerEnd());

    // add relevant tuples for newly created edges
    addNewRelevantTuples(new_node);
    addNewRelevantTuples(new_node.getChild(TreeNode.LEFT_CHILD));
    addNewRelevantTuples(new_node.getChild(TreeNode.RIGHT_CHILD));

    num_missing_leaves --;
  }

  /**
   * Fast-HGT/FP algorithm
   */
  public void buildTree(
    int first_node_index,
    String[] leaf_name,
    DistanceMatrix D) throws HGTException
  {
    initExternalNodes(leaf_name, D);
    init(first_node_index);
    while(!done())
      iterationStep();
  }

//  public static TreeNode FastHGT_FP(BufferedReader Phylip_symmetric_distance_input)
//    throws HGTException
//  {
//    DistanceMatrix D=new DistanceMatrix();
//    String[] leaf_name=null;
//    try {
//      leaf_name=D.readSymmetric(Phylip_symmetric_distance_input);
//    } catch (Exception e)
//    {
//      return null;
//    }
//
//    return FastHGT_FP(FIRST_NODE_INDEX, leaf_name,D);
//  }

  /**
   * hook for extensions
   */
  protected HGTNode newHGTNode(RelevantTuple r){return new HGTNode(r.getTriplet());}
  protected HGTNode newHGTNode(Triplet t){return new HGTNode(t);}

  HGTNode growTree(RelevantTuple r)
  {
    double d_lower = r.getEdgeLengthFromLowerEnd();
    double d_upper = r.getEdgeLengthFromUpperEnd();
    double d_new   = r.getDistanceFromNewLeaf();

    HGTNode new_node = newHGTNode(r);
    TreeNode edge_lower_end = r.getEdgeLowerEnd();
    TreeNode edge_upper_end = edge_lower_end.getParent();
    TreeLeafNode new_leaf = r.getNewLeaf();

    if (edge_lower_end.isOnTheLeft())
    {
      // insert new_node preserving tree orientation
      new_node.addChild(edge_lower_end);
      edge_lower_end.setLength(d_lower);

      new_node.addChild(new_leaf);
      new_leaf.setLength(d_new);

      edge_upper_end.setChild(new_node,TreeNode.LEFT_CHILD);
      new_node.setLength(d_upper);
    } else
    {
      new_node.addChild(new_leaf);
      new_leaf.setLength(d_new);

       new_node.addChild(edge_lower_end);
      edge_lower_end.setLength(d_lower);

      edge_upper_end.setChild(new_node,TreeNode.RIGHT_CHILD);
      new_node.setLength(d_upper);
    }

    new_leaf.setIsInTree(true);
    return new_node;
  }

  public static int FIRST_NODE_INDEX=0;

  /**
   * allocates the array of external nodes and sets the names
   * @param leaf_name can be null
   * @param D distance matrix from which te tree will be built
   */
  public void initExternalNodes(String[] leaf_name, DistanceMatrix D)
  {
    distance=D;
    num_external_nodes=num_missing_leaves=D.getNumNodes();
    external_node = new TreeLeafNode[num_external_nodes];

    for (int i=0; i<num_external_nodes; i++)
    {
      external_node[i]=new TreeLeafNode(i,D);
      if (leaf_name==null)
        external_node[i].setName(new String(Integer.toString(i)));
      else
        external_node[i].setName(leaf_name[i]);
    }
  }

  /**
   * finds the best starting triplet
   */
  private Triplet getStartingTriplet(int first_node_index)
    throws HGTException
  {
    int i_best2, i_best3;
    i_best2=i_best3=0;

    int i1=first_node_index;
    double best_score = 0.;

    for (int i2=0; i2<num_external_nodes; i2++)
      if (i2 != i1)
      {
        double similarity_12 = external_node[i1].getSimilarity(external_node[i2]);
        if (similarity_12 > best_score/3.0) // all of them have to be to beat the highest score
        for (int i3=i2+1; i3<num_external_nodes; i3++)
          if (i3 != i1)
          {
            double similarity_23 = external_node[i2].getSimilarity(external_node[i3]);
            double similarity_31 = external_node[i3].getSimilarity(external_node[i1]);
            if (similarity_23 > best_score/3.0 && similarity_31 > best_score/3.0)
            {
              double score = 3.0/(1.0/similarity_12+1.0/similarity_23+1.0/similarity_31);
              if (score > best_score)
              {
                i_best2=i2;
                i_best3=i3;
                best_score=score;
                //System.out.println("#**A.gST "+i2+" "+i3+" "+score);
              }
            }
          }
      }

    if (i_best2==i_best3)
      throw new HGTException(1, "Cannot find suitable starting triplet for idnex "+first_node_index);


    Triplet first_triplet
      =new Triplet(external_node[first_node_index],external_node[i_best2], external_node[i_best3]);

    //System.out.println("#**A.gST "+first_triplet.toString());

    return first_triplet;
  }

  /**
   * initializes the tree as a star formed by a triplet
   * @return the internal node in the center
   */
  private HGTNode initTreeWithTriplet(Triplet t)
  {
    TreeLeafNode x=t.getTop();
    TreeLeafNode y=t.getLeft();
    TreeLeafNode z=t.getRight();
    HGTNode c=newHGTNode(t);
    x.addChild(c); c.setLength(t.getDistanceFromCenterTotop());
    c.addChild(y); y.setLength(t.getDistanceFromCenterToLeft());
    c.addChild(z); z.setLength(t.getDistanceFromCenterToRight());

    x.setIsInTree(true);
    y.setIsInTree(true);
    z.setIsInTree(true);

    return c;
  }

  void addNewRelevantTuples (
    TreeNode edge_lower_end)
  {
    Vector new_leaf_cache=new Vector(num_external_nodes-3); // its maximal size

    int regular_orientation =
      (edge_lower_end.isOnTheLeft()
        ?RelevantTuple.ORIENTATION_ULN
        :RelevantTuple.ORIENTATION_UNL);

    for (int i=0; i<num_external_nodes; i++)
      if (!external_node[i].isInTree()
          && canInsertOnEdge(edge_lower_end, external_node[i]))
        new_leaf_cache.add(external_node[i]);
    if (edge_lower_end.getParent().isRoot())
    {
      // parent is root [and is also the top member of the defining triplet]
      // edge_lower_end is a HGTNode
      Triplet t = ((HGTNode) edge_lower_end).getDefiningTriplet();
      addNewRelevantTuplesForEdge(edge_lower_end, new_leaf_cache,
        t.getLeft(),
        (TreeLeafNode)edge_lower_end.getParent(),
        regular_orientation);
      addNewRelevantTuplesForEdge(edge_lower_end,new_leaf_cache,
        t.getRight(),
        (TreeLeafNode)edge_lower_end.getParent(),
        regular_orientation);
    } else if (edge_lower_end.isLeaf())
    {
      // edge_lower_end is a leaf [and is in  the defining triplet at parent]
      // edge_upper_end is a HGTNode
      HGTNode edge_upper_end = (HGTNode) edge_lower_end.getParent();
      addNewRelevantTuplesForEdge(edge_lower_end, new_leaf_cache,
        (TreeLeafNode)edge_lower_end,
        edge_upper_end.getDefiningTriplet().getTop(),
        regular_orientation);
      addNewRelevantTuplesForEdge(edge_lower_end,new_leaf_cache,
        (TreeLeafNode)edge_lower_end,
        (edge_lower_end.isOnTheLeft()
          ?edge_upper_end.getDefiningTriplet().getRight()
          :edge_upper_end.getDefiningTriplet().getLeft()),
        regular_orientation);
    } else
    {
      // both lower and upper endpoints are internal nodes
      int flipped_orientation =
      (edge_lower_end.isOnTheLeft()
        ?RelevantTuple.ORIENTATION_LUN
        :RelevantTuple.ORIENTATION_LUN);
      Triplet t_upper=((HGTNode)edge_lower_end.getParent()).getDefiningTriplet();
      Triplet t_lower=((HGTNode)edge_lower_end).getDefiningTriplet();

      addNewRelevantTuplesForEdge(edge_lower_end,new_leaf_cache,
        t_lower.getLeft(),
        t_upper.getTop(),
        regular_orientation);
      addNewRelevantTuplesForEdge(edge_lower_end,new_leaf_cache,
        t_lower.getRight(),
        t_upper.getTop(),
        regular_orientation);
      if (edge_lower_end.isOnTheLeft())
      {
        addNewRelevantTuplesForEdge(edge_lower_end,new_leaf_cache,
          t_lower.getLeft(),
          t_upper.getRight(),
          regular_orientation);
        addNewRelevantTuplesForEdge(edge_lower_end,new_leaf_cache,
          t_lower.getRight(),
          t_upper.getRight(),
          regular_orientation);
        addNewRelevantTuplesForEdge(edge_lower_end,new_leaf_cache,
          t_lower.getTop(),
          t_upper.getLeft(),
          flipped_orientation);
      } else
      { // right edge
        addNewRelevantTuplesForEdge(edge_lower_end,new_leaf_cache,
          t_lower.getLeft(),
          t_upper.getLeft(),
          regular_orientation);
        addNewRelevantTuplesForEdge(edge_lower_end,new_leaf_cache,
          t_lower.getRight(),
          t_upper.getLeft(),
          regular_orientation);
        addNewRelevantTuplesForEdge(edge_lower_end,new_leaf_cache,
          t_lower.getTop(),
          t_upper.getRight(),
          flipped_orientation);
      }
    }
  }

  /**
   * @return null if none left
   */
  RelevantTuple getBestTuple()
  {
    RelevantTuple r_best=null;
    for (int i=0; i<num_external_nodes; i++)
      if (relevant[i] != null && relevant[i].compareTo(r_best)==1)
        r_best=relevant[i];
    return r_best;
 }

  /**
   * @param new_leaves list of new leaves for which the topologies are ok
   */
  private void addNewRelevantTuplesForEdge(
    TreeNode edge_lower_end,
    java.util.List new_leaves,
    TreeLeafNode shared_with_lower,
    TreeLeafNode shared_with_upper,
    int orientation)
  {
    for (int i=0; i<new_leaves.size(); i++)
    {
      TreeLeafNode n=(TreeLeafNode)new_leaves.get(i);
      RelevantTuple r=new RelevantTuple(
        orientation,
        edge_lower_end,
        shared_with_lower,
        shared_with_upper,
        n);

//      if (r.compareTo(relevant[n.getIndex()])<=0
//          && r.getScore()>relevant[n.getIndex()].getScore()-1e-5
//          && !r.getEdgeLowerEnd().equals(relevant[n.getIndex()].getEdgeLowerEnd())
//          )
//        System.out.println("#**A.aNRE -"+r.toString()+"\n#**A.aNRE +"+relevant[n.getIndex()].toString());
      if (r.areEdgeLengthsOK() &&
          r.isBetterThan(relevant[n.getIndex()]))
          relevant[n.getIndex()]=r;
    }
  }

  void deleteRelevantTuplesForEdge(TreeNode edge_lower_end)
  {
    for (int i=0; i<num_external_nodes; i++)
      if (relevant[i] != null &&
          relevant[i].getEdgeLowerEnd().equals(edge_lower_end))
        relevant[i]=null;
  }



  /**
   * checks the quartet topologies for defining triplets and a possible new leaf
   */
  public static boolean canInsertOnEdge(TreeNode edge_lower_end,
    TreeLeafNode n)
  {
    if (edge_lower_end.isRoot())
      return false;
    boolean retval=true;
    if (!edge_lower_end.isLeaf())
    {
      HGTNode x=(HGTNode) edge_lower_end;
      retval = TreeLeafNode.isQuartetOK(
          x.defining_triplet.getLeft(),
          x.defining_triplet.getRight(),
          x.defining_triplet.getTop(),
          n);
      if (!retval)
        return false;
    }
    TreeNode edge_upper_end=edge_lower_end.getParent();
    if (!edge_upper_end.isRoot())
    {
      Triplet t=((HGTNode)edge_upper_end).getDefiningTriplet();
      if (edge_lower_end.isOnTheLeft())
        retval = TreeLeafNode.isQuartetOK(
          t.getTop(),
          t.getRight(),
          t.getLeft(),
          n);
      else
        retval = TreeLeafNode.isQuartetOK(
          t.getTop(),
          t.getLeft(),
          t.getRight(),
          n);
    }
    return retval;
  }

  public static class HGTException extends Exception
  {
    public HGTException(int error_id, String error_message)
    {
      super("Error #"+error_id+":"+error_message);
    }
  }

  static String RCS_ID="$Id: Algorithm.java,v 1.1 2001/03/09 17:31:04 mcsuros Exp mcsuros $";
  //
  // $Log: Algorithm.java,v $
  // Revision 1.1  2001/03/09 17:31:04  mcsuros
  // Initial revision
  //
  //

}