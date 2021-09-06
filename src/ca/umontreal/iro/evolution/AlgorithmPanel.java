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
import java.io.Reader;
import java.io.FileReader;

import java.awt.BorderLayout;
import javax.swing.*;
import java.awt.Dimension;
import java.awt.event.*;

public class AlgorithmPanel extends JPanel {
  public AlgorithmPanel(Reader input)
  {
    super(new BorderLayout(),true);
    try {
      setup(input);
    } catch (Exception e)
    {
    }
  }



  DistanceMatrix distance;
  VisualAlgorithm algo;
  TreePanel tree_display;

  JButton step_button;

  boolean initialized=false;

  private void setup(Reader input) throws java.io.IOException, Algorithm.HGTException, Parser.ParseException
  {
    distance=new DistanceMatrix();
    String[] leaf_name=distance.readSymmetric(new BufferedReader(input));
    algo = new VisualAlgorithm();
    algo.initExternalNodes(leaf_name,distance);
    tree_display = new TreePanel(algo);
    step_button=new JButton(" > ");
    add(step_button,BorderLayout.NORTH);
    add(tree_display,BorderLayout.CENTER);

    step_button.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent event){
        algoStep();
    }});
    //step_button.setMaximumSize(new Dimension(100,10));

    setVisible(true);
  }

  void algoStep()
  {
    try {
      if (!initialized)
      {
        algo.init(0);
        tree_display.repaint();
        initialized=true;
      }
      else
      {
        algo.iterationStep();
        tree_display.repaint();
      }
    } catch (Exception e)
    {

    }

  }



}