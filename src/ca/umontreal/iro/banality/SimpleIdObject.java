/*
 * SimpleIdObject.java
 *
 * Created on April 10, 2003, 10:03 AM
 */

package ca.umontreal.iro.banality;

/**
 *
 * @author  csuros
 */
public class SimpleIdObject implements IdObject {
  private int id;
  private static int next_id=1;

  public SimpleIdObject() {
    id=next_id++;
  }

  protected String paramString(){
    return "";
  }

  /**
   * buffers identString()
   */
  private String ident_save=null;
  public String identString(){
    if (ident_save !=null) return ident_save;

    StringBuffer sb=new StringBuffer(getClass().getName());
    sb.append('#');
    sb.append(id);
    return (ident_save=sb.toString());
  }

  public String toString(){
    StringBuffer sb=new StringBuffer(identString());
    sb.append('[');
    sb.append(paramString());
    sb.append(']');
    return sb.toString();
  }

  public int getObjectId(){
    return id;
  }

  /**
   * supported but buffered ident must be changed: clone gets new id also
   */
  protected Object clone() throws java.lang.CloneNotSupportedException {
    SimpleIdObject o=(SimpleIdObject)super.clone();
    o.id=next_id++;
    o.ident_save=null;
    return o;
  }
   
}
