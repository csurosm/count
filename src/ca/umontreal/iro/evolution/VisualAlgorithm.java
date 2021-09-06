package ca.umontreal.iro.evolution;

/**
 * Title:
 * Description:
 * Copyright:    Copyright (c) 2001
 * Company:
 * @author
 * @version 1.0
 */

import java.awt.Point;

public class VisualAlgorithm extends Algorithm {

  Point[] node_location;
  int[] subtree_size;

  // everything will be in a 1.0x1.0 rectangle
  // but this determines the resolution
  final double height=2048.0;

  final double edge_length_display_factor=2048.0;

  int num_nodes;

  public VisualAlgorithm() {
  }

  public void initExternalNodes(String[] leaf_name, DistanceMatrix dm)
  {
    super.initExternalNodes(leaf_name, dm);
    node_location=new Point[2*num_external_nodes-2];
    subtree_size = new int[2*num_external_nodes-2];
    // VisualHGTNodes will have indices n...2*n-3
    num_nodes=num_external_nodes;
    for (int i=0; i<num_nodes; i++)
      node_location[i]=new Point();
    setMissingNodeLocations();
  }

  public HGTNode init(int first_node_index) throws HGTException
  {
    VisualHGTNode c=(VisualHGTNode)super.init(first_node_index);

    // set locations for the first 3 nodes
    TreeLeafNode root=(TreeLeafNode)c.getParent();
    TreeLeafNode x=(TreeLeafNode)c.getChild(TreeNode.LEFT_CHILD);
    TreeLeafNode y=(TreeLeafNode)c.getChild(TreeNode.RIGHT_CHILD);

    node_location[root.getIndex()].setLocation(0.,height/2.0);

    double c_loc = c.getLength()*edge_length_display_factor;
    node_location[c.getIndex()].setLocation(
      c_loc,height/2.0);

    double x_loc = c_loc+x.getLength()*edge_length_display_factor;
    node_location[x.getIndex()].setLocation(
      x_loc, 0);

    double y_loc = c_loc+y.getLength()*edge_length_display_factor;
    node_location[y.getIndex()].setLocation(
      y_loc, height);

    tree_width=Math.max(x_loc,y_loc);
    tree_height=height;

    setMissingNodeLocations();

    subtree_size[root.getIndex()]=2;
    subtree_size[c.getIndex()]=2;
    subtree_size[x.getIndex()]=subtree_size[y.getIndex()]=1;

    return c;
  }

  HGTNode growTree(RelevantTuple r)
  {
    VisualHGTNode c=(VisualHGTNode)super.growTree(r);

    // set locations for the first 3 nodes

    TreeNode papa=c.getParent();
    int idx_papa
      =(papa.isExternal()
        ? ((TreeLeafNode)papa).getIndex()
        : ((VisualHGTNode)papa).getIndex());
    int idx=c.getIndex();

    TreeNode old_child;
    TreeLeafNode new_child;

    if (c.isOnTheLeft())
    {
      // insertion was on the left
      old_child=c.getChild(TreeNode.LEFT_CHILD);
      new_child=(TreeLeafNode)c.getChild(TreeNode.RIGHT_CHILD);
    } else
    {
      old_child=c.getChild(TreeNode.RIGHT_CHILD);
      new_child=(TreeLeafNode)c.getChild(TreeNode.LEFT_CHILD);
    }

    int idx_old_child
    =(old_child.isExternal()
      ? ((TreeLeafNode)old_child).getIndex()
      : ((VisualHGTNode)old_child).getIndex());

    int idx_new_child = new_child.getIndex();

//    System.out.println("#**VA.gT "
//        +" c "+idx+", p "+idx_papa+", o "+idx_old_child+", n "+idx_new_child+"; "
//        +c.toString()+", papa "+papa.toString()+", old "+old_child.toString()+", new "+new_child.toString());

//    System.out.println("#**VA.gT >> "
//        +" c "+node_location[idx].toString()
//        +" p "+node_location[idx_papa].toString()
//        +" o "+node_location[idx_old_child].toString()
//        +" n "+node_location[idx_new_child].toString());
//
    double x=node_location[idx_papa].getX();
    x+=c.getLength()*edge_length_display_factor;

    if (papa.isRoot())
    {
      // added on right edge for sure
      node_location[idx_papa].setLocation(0.,tree_height);

      node_location[idx].setLocation(
        x,
        tree_height);
      tree_width=Math.max(x,tree_width);


      x+=new_child.getLength()*edge_length_display_factor;
      double y=tree_height-node_location[idx_old_child].getY();
      node_location[idx_new_child].setLocation(
        x,
        tree_height+y-10.);
      tree_width=Math.max(x,tree_width);

      tree_height += y;
    } else
    {
      node_location[idx].setLocation(
        x,
        (node_location[idx_papa].getY()+2.*node_location[idx_old_child].getY())/3.0);

      tree_width=Math.max(tree_width,x);
      x+=new_child.getLength()*edge_length_display_factor;

      node_location[idx_new_child].setLocation(
        x,
        (2.*node_location[idx_papa].getY()+node_location[idx_old_child].getY())/3.0);
      tree_width=Math.max(x,tree_width);
    }

    subtree_size[idx]=subtree_size[idx_old_child];
    increaseSubtreeSize(new_child);
    tree_height=calcY(getRoot(),-1)*height;

//    System.out.println("#**VA.gT << "
//        +" c "+node_location[idx].toString()
//        +" p "+node_location[idx_papa].toString()
//        +" o "+node_location[idx_old_child].toString()
//        +" n "+node_location[idx_new_child].toString());

    setMissingNodeLocations();

    return c;
  }

  private void increaseSubtreeSize(TreeNode n)
  {
    if (n.isExternal())
    {
      subtree_size[((TreeLeafNode)n).getIndex()]++;
    } else
    {
      subtree_size[((VisualHGTNode)n).getIndex()]++;
    }
    if (!n.isRoot())
      increaseSubtreeSize(n.getParent());
  }

  private int calcY(TreeNode n, int lastY)
  {
    if (n.isLeaf())
    {
      int idx=((TreeLeafNode)n).getIndex();
      node_location[idx].y=
        (int)((lastY+1.)*height);
//      System.out.println("#**VA.cY L "+idx+" "+lastY+" "+subtree_size[idx]+"; "+node_location[idx].toString());
      return lastY+1;
    }
    else
    {
      double y=0.;
      int retval=lastY;
      for (int i=0; i<n.getNumChildren(); i++)
      {
        retval=calcY(n.getChild(i),retval);
        TreeNode c=n.getChild(i);
        if (c.isLeaf())
          y+=node_location[((TreeLeafNode)c).getIndex()].getY()*getSubtreeSize(c);
        else
          y+=node_location[((VisualHGTNode)c).getIndex()].getY()*getSubtreeSize(c);
      }
      if (!n.isRoot())
      {
        int idx=((VisualHGTNode)n).getIndex();
        node_location[idx].y=(int)(y/subtree_size[idx]);
//        System.out.println("#**VA.cY N "+idx+" "+y+" "+subtree_size[idx]+"; "+node_location[idx].toString());
      } else
        node_location[((TreeLeafNode)n).getIndex()].y=(int)(retval*height/2.);
      return retval;
    }
  }

  protected HGTNode newHGTNode(RelevantTuple r)
  {
    HGTNode n=new VisualHGTNode(r.getTriplet());
    //System.out.println("#**VA.nHN "+r.toString()+" "+n.toString());
    return n;
  }

  protected HGTNode newHGTNode(Triplet t){
    HGTNode n=new VisualHGTNode(t);
    //System.out.println("#**VA.nHN "+t.toString()+" "+n.toString());
    return n;
  }

  void placeTreeNodes(int w, int h)
  {
    System.out.println("#**VA.pTN "+w+"x"+h);
  }

  double getX(TreeLeafNode n)
  {
    if (n.isInTree())
      return node_location[n.getIndex()].getX()/tree_width;
    else
      return 0.;
  }

  double getY(TreeLeafNode n)
  {
    if (n.isInTree())
      return node_location[n.getIndex()].getY()/tree_height;
    else
      return node_location[n.getIndex()].getY()/height;
  }

  double getX(HGTNode n)
  {
    return node_location[((VisualHGTNode)n).getIndex()].getX()/tree_width;
  }

  double getY(HGTNode n)
  {
    return node_location[((VisualHGTNode)n).getIndex()].getY()/tree_height;
  }

  double getX(TreeNode n)
  {
    if (n.isExternal())
      return getX((TreeLeafNode)n);
    else
      return getX((HGTNode)n);
  }

  double getY(TreeNode n)
  {
    if (n.isExternal())
      return getY((TreeLeafNode)n);
    else
      return getY((HGTNode)n);
  }
  int getSubtreeSize(TreeNode n)
  {
    if (n.isExternal())
      return subtree_size[((TreeLeafNode)n).getIndex()];
    else
      return subtree_size[((VisualHGTNode)n).getIndex()];
  }


  double tree_width=0.;
  double tree_height=height;

  //double getTreeWidth(){return tree_width;}

//  Point getLocation(HGTNode n){return node_location[((VisualHGTNode)n).getIndex()];}

  private void setMissingNodeLocations()
  {
    int j=1;
    double dy=height/(num_external_nodes+1.0);
    for (int i=0; i<num_external_nodes; i++)
    {
      if (!external_node[i].isInTree())
      {
        node_location[i].setLocation(0.,j*dy);
        //System.out.println("#**VA.sMN "+i+"@ "+node_location[i].toString());
        j++;
      }
    }
  }

  class VisualHGTNode extends HGTNode
  {
    int index;
    VisualHGTNode(Triplet t)
    {
      super(t);
      node_location[num_nodes]=new Point();
      index=num_nodes;
      num_nodes++;
    }

    int getIndex(){return index;}
  }

}