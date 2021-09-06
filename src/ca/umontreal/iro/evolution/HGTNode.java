package ca.umontreal.iro.evolution;

/**
 * Title:        HGTNode
 * Description:  defining triplets
 * Copyright:    Copyright (c) 2001
 * Company:
 * @author
 * @version 1.0
 */

public class HGTNode extends TreeNode  {

  Triplet defining_triplet;

  public HGTNode(Triplet defining_triplet)
  {
    //super(-1,null);
    this.defining_triplet=defining_triplet;
  }

//  HGTNode(int index, DistanceMatrix dm)
//  {
//      super(index,dm);
//      this.defining_triplet=null;
//  }

  public Triplet getDefiningTriplet(){return defining_triplet;}

  static String RCS_ID="$Id: HGTNode.java,v 1.1 2001/03/09 17:31:04 mcsuros Exp mcsuros $";
  //
  // $Log: HGTNode.java,v $
  // Revision 1.1  2001/03/09 17:31:04  mcsuros
  // Initial revision
  //
  //


}
