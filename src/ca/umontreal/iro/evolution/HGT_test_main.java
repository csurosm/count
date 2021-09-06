package ca.umontreal.iro.evolution;

/**
 * Title:
 * Description:
 * Copyright:    Copyright (c) 2001
 * Company:
 * @author
 * @version 1.0
 */

import java.io.*;

import javax.swing.*;
import java.awt.BorderLayout;
import java.awt.event.*;
import java.awt.Toolkit;

public class HGT_test_main {
  public static void main(String[] args) {
    launch();
  }

  public static void launch()
  {
    JFrame F = new JFrame("Fast-HGT/FP");
    F.setBounds(25,25,
      (int)Toolkit.getDefaultToolkit().getScreenSize().getWidth()-50,
      (int)Toolkit.getDefaultToolkit().getScreenSize().getHeight()-50);

//    TreeNode root=new TreeNode();
//    TreeNode a=new TreeNode();
//    root.addChild(a);
//    a.setLength(Double.NaN);
//    a.setName("Bt");
//
//    System.out.println(a.toString());
//    System.out.println(a.newickName(true));
//
//    System.out.println(root.newickSubtree(false,true,false,false));

    try {

//      String tree = "(((2:-Inf, ((3:0.1, 7:):-1e9, (10:0.1, 6:NaN):Inf):0.1):0.1, ((4:0.1,8:0.1):0.1, 5:0.1):0.1):0.1, 9:0.1, 1:1);";
//
//      root = Parser.readNewick(new StringReader(tree));
//
//      System.out.println(root.newickTree(false,true,false,false));
//
//
//      String matrix="3\negy 11 12 13\nketto 12 22 23\nharom 13 23 33";
//      DistanceMatrix M=new DistanceMatrix();


      Reader input=
        //new FileReader(ph_file);
        new StringReader(proba_ph);

//      String name[]=M.readSymmetric(new BufferedReader(input));
//      //System.out.println(M.getMaxFinitValue());
//
//      input.close();

//      System.out.println(M.toString());
//      for (int i=0; i<name.length; i++)System.out.println(name[i]);

//      VisualAlgorithm A = new VisualAlgorithm();
      //A.buildTree(0,name,M);
//      A.initExternalNodes(name,M);

      AlgorithmPanel ap=new AlgorithmPanel(input);
//    tp.setExternalNodes(A.getExternalNodes());



    F.getContentPane().add(ap);

    F.setVisible(true);

    F.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
    F.addWindowListener(new WindowAdapter(){
      public void windowClosed(WindowEvent e) {
        System.exit(0);
      }
    });

//    A.init(0);
//    A.iterationStep();
//    A.iterationStep();
//    A.iterationStep();
//    while(!A.done())
//      A.iterationStep();
//      TreeNode r=A.getRootDegree3();
//      //System.out.println(r.subtreeParamString());
//
//      System.out.println(r.newickTree(false,true,false,false));


    } catch (Throwable t){
      System.out.println(t.toString());
      t.printStackTrace();
    }
  }
  static String proba_ph
  = "10\n"
+"2         0 0.4081 0.4211 0.3946 0.4273 0.5075 0.4901 0.3897 0.3959 6.0600\n"
+"3         0.4081 0 0.1980 0.3943 0.4227 0.7544 0.6945 0.6068 0.6184 5.9700\n"
+"7         0.4211 0.1980 0 0.3784 0.3960 0.7526 0.7148 0.6060 0.5973 6.0600\n"
+"10        0.3946 0.3943 0.3784 0 0.1929 0.7465 0.7142 0.5793 0.5447 6.0600\n"
+"6         0.4273 0.4227 0.3960 0.1929 0 0.7937 0.7419 0.5913 0.5526 5.9700\n"
+"4         0.5075 0.7544 0.7526 0.7465 0.7937 0 0.2089 0.3247 0.5022 6.0600\n"
+"8         0.4901 0.6945 0.7148 0.7142 0.7419 0.2089 0 0.2886 0.4595 6.0600\n"
+"5         0.3897 0.6068 0.6060 0.5793 0.5913 0.3247 0.2886 0 0.3806 5.9700\n"
+"9         0.3959 0.6184 0.5973 0.5447 0.5526 0.5022 0.4595 0.3806 0 6.0600\n"
+"1         6.0600 5.9700 6.0600 6.0600 5.9700 6.0600 6.0600 5.9700 6.0600 0";
  static String ph_file = "E:\\mcsuros\\log\\SSU\\maki.ph";

  //eve135\\eve135-8.1000.10k.ph";
}