/*
 * BrowserSelectionListener.java
 *
 * Created on December 30, 2007, 1:49 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package ca.umontreal.iro.evolution.malin.ui;

/**
 *
 * @author csuros
 */
public interface BrowserSelectionListener  extends java.util.EventListener
{
    public void valueChanged(Browser.PrimaryItemSelectionEvent E);
}
