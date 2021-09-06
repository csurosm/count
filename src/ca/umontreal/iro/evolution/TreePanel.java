package ca.umontreal.iro.evolution;

/**
 * Title:
 * Description:
 * Copyright:    Copyright (c) 2001
 * Company:
 * @author
 * @version 1.0
 */
import javax.swing.JPanel;

import java.awt.Graphics;
import java.awt.Point;

public class TreePanel extends JPanel {
  VisualAlgorithm algo;

  public TreePanel(VisualAlgorithm A) {
    super(true);
    algo=A;
    window_resized=true;
  }

  private static int prevWidth=0;
  private static int prevHeight=0;
  private static boolean window_resized=true;

  public void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    if (prevWidth != getWidth() || prevHeight != getHeight())
    {
      window_resized=true;
      prevWidth=getWidth();
      prevHeight=getHeight();
    }

    plotTree(g);
    plotExternalNodes(g);
    g.setColor(java.awt.Color.red);
    g.drawLine(prevWidth-missing_nodes_width,0,prevWidth-missing_nodes_width,prevHeight);

    window_resized=false;
  }

  int NODE_SIZE=4;
  DiamondIcon node_icon =  new DiamondIcon(NODE_SIZE, false);
  BoxIcon leaf_icon = new BoxIcon(NODE_SIZE,false);

  private void plotExternalNodes(Graphics g)
  {
    for (int i=0; i<algo.getNumExternalNodes(); i++)
      plotExternalNode(g,algo.getExternalNode(i));
  }

  int missing_nodes_width=100;

  private int transformX(double x)
  {
    return NODE_SIZE+(int) (x*(prevWidth-missing_nodes_width-2*NODE_SIZE));
  }

  private int transformY(double y)
  {
    return 2*NODE_SIZE+(int) (y*(prevHeight-4*NODE_SIZE));
  }

  private void plotExternalNode(Graphics g, TreeLeafNode n)
  {
    int x_offset=NODE_SIZE/2+4;
    int y_offset=g.getFont().getSize()/2;

    int x=prevWidth-missing_nodes_width+2*NODE_SIZE;
    if (n.isInTree())
    {
      x=transformX(algo.getX(n));
      if (n.isRoot())
      {
        x_offset=0;
        y_offset=-NODE_SIZE/2-3;
      }
    }
    int y=transformY(algo.getY(n));
    leaf_icon.paint(g,x,y);
    g.drawString(n.getName(),x+x_offset,y+y_offset);
  }

  private void plotInternalNode(Graphics g, HGTNode n)
  {

    int x=transformX(algo.getX(n));
    int y=transformY(algo.getY(n));
    node_icon.paint(g,x,y);
    //g.drawString(n.shortDesc(),x+6,y);

  }

  private void plotTree(Graphics g){
    plotSubTree(g, algo.getRoot());
  }

  private void plotSubTree(Graphics g, TreeNode n)
  {
    if (n==null) return;
    for (int i=0; i<n.getNumChildren(); i++)
      plotEdge(g,n,n.getChild(i));

    if (n.isExternal())
      plotExternalNode(g,(TreeLeafNode)n);
    else
      plotInternalNode(g,(HGTNode)n);

    for (int i=0; i<n.getNumChildren(); i++)
      plotSubTree(g,n.getChild(i));
  }

  private void plotEdge(Graphics g, TreeNode node_from, TreeNode node_to)
  {
    int x0=transformX(algo.getX(node_from));
    int y0=transformY(algo.getY(node_from));
    int x1=transformX(algo.getX(node_to));
    int y1=transformY(algo.getY(node_to));

    g.drawLine(x0,y0,x0,y1);
    g.drawLine(x0,y1,x1,y1);
  }



}