/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ca.umontreal.iro.evolution.malin.ui.count;

/**
 * A common interface for displays of families, where 
 * one can create a new display for selected families only. 
 * 
 * @author csuros
 */

import javax.swing.JComponent;

public interface FilterableFamilies 
{   
    /**
     * This should return a display of the same class, which will be attached 
     * as a dependent item in the browser.  
     * 
     * @return a new display for the data or analysis browsers
     */
    public JComponent newDisplayWithSelectedFamilies();
}
