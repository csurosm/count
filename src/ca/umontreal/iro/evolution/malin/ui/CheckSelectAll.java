/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ca.umontreal.iro.evolution.malin.ui;

import java.awt.event.ItemListener;
import java.awt.event.ItemEvent;

import java.util.Set;
import java.util.HashSet;
import java.util.Iterator;

import javax.swing.JCheckBox;
import javax.swing.AbstractButton;

/**
 * Swing component for a "select all" checkbox
 * for a group of other check boxes.  When this checkbox is
 * selected (programmatically or by an action), all
 * attached boxes get selected. If an attached box is deselected
 * then this box gets deselected too (unless it was already deselected).
 *
 * @param <B> type for attached boxes (e.g., JCheckBox)
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class CheckSelectAll<B extends AbstractButton> extends JCheckBox //implements Iterable<B>
{
    /**
     * Instantiation with default JCheckBox text and icon.
     */
    public CheckSelectAll(){super(); init();}
    public CheckSelectAll(String text){super(text); init();}
    public CheckSelectAll(javax.swing.Action action){super(action); init();}
    public CheckSelectAll(javax.swing.Icon icon){super(icon); init();}
    public CheckSelectAll(String text, boolean select){super(text); init(); setSelected(select);}
    public CheckSelectAll(String text, javax.swing.Icon icon, boolean selected)
    {
        super(icon);
        init();
        setText(text);
        setSelected(selected);
    }


    private Set<B> attached_boxes;
    private ItemListener attached_box_listener;

    private void init()
    {
        attached_boxes = new HashSet<B>();
        attached_box_listener = new ItemListener()
            {
                @Override
                public void itemStateChanged(ItemEvent event_in_attached_box)
                {
                    if (isSelected())
                    {
                        AbstractButton src = (AbstractButton) event_in_attached_box.getSource();
                        if (ItemEvent.DESELECTED == event_in_attached_box.getStateChange())
                            setSelected(false);
                    }
                }

            };
        addItemListener(new ItemListener()
            {
                @Override
                public void itemStateChanged(ItemEvent event)
                {
                    if (ItemEvent.SELECTED == event.getStateChange())
                    {
                        for (B cb : attached_boxes)
                            cb.setSelected(true);
                    }
                }

        });
    }

    /**
     * Attaches a check box to this one.
     *
     * @param cb check box
     * @return true if this check box was not attached yet
     */

    public boolean addCheckBox(B cb)
    {
        boolean not_seen = attached_boxes.add(cb);
        if (not_seen) // if this box was not added yet
        {
            attached_boxes.add(cb);
            if (isSelected())
                cb.setSelected(true);
            cb.addItemListener(attached_box_listener);
        }
        return not_seen;
    }

    /**
     * Attaches a set of buttons to the master check box.
     *
     * @param check_boxes collection
     * @return whether at least one of them was not added yet
     */
    public boolean addCheckBoxes(java.util.Collection<B> check_boxes)
    {
        boolean not_seen = false;
        for (B cb: check_boxes)
            not_seen |= addCheckBox(cb);
        return not_seen;
    }

    /**
     * Attaches a set of buttons to the master check box.
     *
     * @param check_boxes non-null array
     * @return whether at least one of them was not added yet
     */
    public boolean addCheckBoxes(B[] check_boxes)
    {
        boolean not_seen = false;
        for (B cb: check_boxes)
            not_seen |= addCheckBox(cb) ;
        return not_seen;
    }


    /**
     * Removes a check box from dependency.
     *
     * @param cb check box
     * @return true if the check box was really attached
     */
    public boolean removeCheckBox(B cb)
    {
        boolean was_attached = attached_boxes.remove(cb);
        if (was_attached)
            cb.removeItemListener(attached_box_listener);
        return was_attached;
    }

    /**
     * Copies the current state (selected or deselected) of his button
     * to all attached buttons.
     */
    public void copySelectedToAll()
    {
        final boolean sel = isSelected();
        for (B cb: attached_boxes)
            cb.setSelected(sel);
    }

    /**
     * Enables/disables all attached check boxes.
     *
     * @param enabled
     */
    public void setEnabledAll(boolean enabled)
    {
        for (B cb: attached_boxes)
            cb.setEnabled(enabled);
    }
    
//    /**
//     * Iterator over the attached buttons.
//     *
//     * @return iterator
//     */
//    @Override
//    public Iterator<B> iterator()
//    {
//        return attached_boxes.iterator();
//    }





}
