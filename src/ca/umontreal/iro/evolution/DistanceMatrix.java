package ca.umontreal.iro.evolution;

/**
 * Title:        DistanceMatrix
 * Description:  stores a distance matrix internally as a lower triangular matrix
 * Copyright:    Copyright (c) 2001
 * Company:
 * @author       Miklos Csuros csuros AT iro.umontreal.ca
 * @version 1.0
 */

import java.io.*;
import java.util.Vector;

public class DistanceMatrix {
  int size;
  private int capacity;
  double[][] distance; // triangular array

  public DistanceMatrix() {size=capacity=0;}

  /**
   * @return distance between nodes i and j, for i,j = 0...size [inclusive]
   */
  public double getDistance(int i, int j)
  {
    if (i>j) return distance[i-1][j];
    if (i<j) return distance[j-1][i];
    return 0.;
  }

  public double getSimilarity(int i, int j)
  {
    double d=getDistance(i,j);
    if (d>0.) return Math.exp(-d);
    return 0.;
  }

  public int getNumNodes(){return size+1;}

  private void newRow(){
    if (capacity==0)
    {
      capacity=10;
      distance=new double[capacity][0];
    } else if (capacity==size)
    {
      int new_capacity=capacity*2;
      double[][] new_distance=new double[new_capacity][0];
      System.arraycopy(distance,0,new_distance,0,capacity);
      distance=new_distance;
      capacity=new_capacity;
    }
    size++;
  }

//  public void readUpperTriangular(BufferedReader r) throws IOException, Parser.ParseException
//  {
//    int row=0;
//    int column=0;
//
//    while(true)
//    {
//      String line=r.readLine();
//      if (line==null) break;
//
//      PushbackReader br=new PushbackReader(new java.io.StringReader(line));
//      int c;
//      do
//      {
//        if (row==0){
//          newRow();
//          distance[column]=new double[column+1];
//        } else if (row+column>=size)
//        {
//          row++;
//          throw new Parser.ParseException(43, "Upper triangular matrix has too many columns in row "+row);
//        }
//
//        distance[row+column][row]=Parser.parseDouble(br);
//        column++;
//        c=Parser.nextNonWhitespace(br);
//        if (c!=-1) br.unread(c);
//      } while (c != -1);
//      row++;
//      column=0;
//    }
//  }
//
//  public void readLowerTriangular(BufferedReader r) throws IOException, Parser.ParseException
//  {
//    int row=0;
//    int column=0;
//
//    while(true)
//    {
//      String line=r.readLine();
//      if (line==null) break;
//
//      newRow();
//      distance[row]=new double[row+1];
//
//      PushbackReader br=new PushbackReader(new java.io.StringReader(line));
//      int c;
//      do
//      {
//        if (column>row)
//        {
//          row++;
//          throw new Parser.ParseException(44, "Lower triangular matrix has too many columns in row "+row);
//        }
//
//        distance[row][column]=Parser.parseDouble(br);
//        column++;
//        c=Parser.nextNonWhitespace(br);
//        if (c!=-1) br.unread(c);
//      } while (c != -1);
//      row++;
//      column=0;
//    }
//  }


  public String[] readSymmetric(BufferedReader r) throws IOException, Parser.ParseException
  {
    int row=0;
    int column=0;

    int num_leaves=Integer.parseInt(r.readLine());
    Vector node_name_vector=new Vector();

    while(true)
    {
      String line=r.readLine();
      if (line==null) break;

      PushbackReader br=new PushbackReader(new java.io.StringReader(line));

      if (column==0 && node_name_vector.size()==row)
        node_name_vector.add(Parser.parseString(br,true));
      //System.out.println("#**DM.rS @("+row+", "+column+") `"+line+"'");

      // skip whitespace separating name from the distances
      int c = Parser.nextNonWhitespace(br);
      if (c!=-1)
      {
        br.unread(c);
        do
        {
          double d=Parser.parseDouble(br);
          //System.out.println("#**DM.rS "+node_name_vector.get(row)+" @("+row+", "+column+") :"+d);

          if (row==0) {
            if (column>0)
            {
              newRow();
              distance[column-1]=new double[column];
              distance[column-1][0]=d;
            }
          } else
          {
            if (column>row)
              distance[column-1][row]=d;
            else if (column<row)
              distance[row-1][column]
                =(distance[row-1][column]+d)/2.0;

          }
          column++;
          c=Parser.nextNonWhitespace(br);
          if (c!=-1) br.unread(c);
        } while (c != -1);
      }
      if (column == num_leaves)
      {
        if (row % 100==0)
          System.out.println("Row "+row+" done.");
        row++;
        column=0;
      }
    }

    if (size != num_leaves-1)
    {
      size++;
      throw new Parser.ParseException(45,
        "Number of rows read ["+size+"] does not equal number of nodes expected ["+num_leaves+"]");
    }

    return (String[]) node_name_vector.toArray(new String[0]);


  }

  public String toString()
  {
    String retval="[failed toString]";
    try
    {
      StringWriter w=new StringWriter(size*size*6);
      writeMatrix(w);
      retval=w.toString();
    } catch (IOException ex)
    {
      //whatever
    }
    return retval;
  }

  public void writeMatrix(Writer w) throws IOException
  {
    for (int i=0; i<=size; i++)
    {
      if (i>0)
        w.write('\n');
      for (int j=0; j<=size; j++)
      {
        w.write(TreeNode.toIEEEFormat(getDistance(i,j)));
        w.write(' ');
      }
    }
  }

  /**
   * checks whether the quartet topology is ij|kl using the four-point condition
   */
  public boolean isQuartetOK(int i, int j, int k, int l)
  {
    double dist_small = getDistance(i,j)+getDistance(k,l);
    double dist_big =getDistance(i,k)+getDistance(j,l);
    if (dist_small >= dist_big)
      return false;
    dist_big = getDistance(i,l)+getDistance(k,j);
    return (dist_small < dist_big);
  }

  public double getMaxFinitValue()
  {
    double max=Double.NEGATIVE_INFINITY;
    for (int i=0; i<size; i++)
      for (int j=0; j<i; j++)
        if (!Double.isInfinite(distance[i][j])
            && !Double.isNaN(distance[i][j])
            && distance[i][j]>max)
          max=distance[i][j];
    return max;
  }

  static String RCS_ID="$Id: DistanceMatrix.java,v 1.1 2001/03/09 17:31:04 mcsuros Exp mcsuros $";
  //
  // $Log: DistanceMatrix.java,v $
  // Revision 1.1  2001/03/09 17:31:04  mcsuros
  // Initial revision
  //
  //

}