package count.gui.kit;
/*
 * Copyright 2022 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import java.awt.geom.Point2D;

/**
 *
 * A point with an associated integer index, set at instantiation.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s csurosm@gmail.com
 * @since November 14, 2007, 12:02 AM
 */
public class IndexedPoint extends Point2D.Double 
{
   private final int index;

   /**
   * Instantiation (0,0).
   * @param index the index associated with this point
   */
   public IndexedPoint(int index) 
   {
       this(index,0.0,0.0);
   }

   /**
    * instantiation with coordinates.  
    * @param index 
    * @param x
    * @param y 
    */
   public IndexedPoint(int index, double x, double y)
   {
       super(x,y);
       this.index = index;
   }

   public IndexedPoint(int index, Point2D p)
   {
       this(index, p.getX(), p.getY());
   }

//   private void init(int index)
//   {
//       this.index=index;
//   }
//
//   public void setIndex(int new_index){index=new_index;}

   public int getIndex(){return index;}

   protected String paramString()
   {
       StringBuilder sb=new StringBuilder();
       sb.append("idx ");
       sb.append(index);
       sb.append(" (");
       sb.append(x);
       sb.append(", ");
       sb.append(y);
       sb.append(")");
       return sb.toString();
   }

   @Override
   public String toString()
   {
       StringBuilder sb=new StringBuilder("iP");//getClass().getName());
       sb.append('[');
       sb.append(paramString());
       sb.append(']');
       return sb.toString();
   }
   
   @Override
   public int hashCode()
   {
       return super.hashCode()*3+getIndex();
   }
   
   @Override
   public boolean equals(Object O)
   {
       if (O instanceof IndexedPoint)
       {
          IndexedPoint P = (IndexedPoint) O;
          return P.getIndex() == getIndex() && super.equals(P);
       } else
       {
           return super.equals(O);
       }
   }
}