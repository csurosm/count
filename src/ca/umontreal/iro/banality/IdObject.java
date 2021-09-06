/*
 * IdObject.java
 *
 * Created on April 10, 2003, 10:01 AM
 */

package ca.umontreal.iro.banality;

/**
 *
 * @author  csuros
 */
public interface IdObject {
  /**
   * The standard toString() format is identString()``[''paramString()``]''.
   * <code>paramString()</code> is <code>protected</code>, so
   * this interface cannot define it. identString() should return
   * &lt;class name&gt;``#''getObjectId()
   *
   * @return class name and unique id within class.
   */
  public String identString();

  /**
   * @return unique id within class.
   */
  public int getObjectId();


    
}
