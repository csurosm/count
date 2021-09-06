
package ca.umontreal.iro.evolution.malin.ui;


import java.awt.BorderLayout;
//import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;


import javax.swing.Box;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SpinnerModel;


/**
 * The basic implementation of a zoomable tree display.
 * The tree (EmbellishedTreePanel)
 * is within a JScrollPane and there is a zoom spinner below
 * it.
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class ZoomableTreePanel extends JPanel
{
    public ZoomableTreePanel(EmbellishedTreePanel TP)
    {
        this.tree_panel = TP;
        initComponents();
    }
    
    private EmbellishedTreePanel tree_panel;
    
    private JScrollPane tree_scroll;
    private JSpinner    magnification_spin;
    
    private void initComponents()
    {
        this.setLayout(new BorderLayout());

        tree_scroll = new JScrollPane(tree_panel);
        tree_scroll.getViewport().setBackground(tree_panel.getBackground());
        tree_scroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        tree_scroll.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        this.add(tree_scroll, BorderLayout.CENTER);
        
        Box bottom_bar = Box.createHorizontalBox();
        
        initBottomBarElements(bottom_bar);

        this.add(bottom_bar, BorderLayout.SOUTH);
        
        setMagnification(1.0);
    }
    
    /**
     * Called from initComponents() to add the elements in the bottom bar. 
     * The default implementation adds a magnification spinner. 
     * 
     * @param bottom_bar the Box containing the bottom bar
     */
    protected void initBottomBarElements(Box bottom_bar)
    {
        bottom_bar.add(Box.createGlue());        
        
        SpinnerModel spinner_model = 
                //new SpinnerListModel((Object[])mag_values);
                new SpinnerMagnifierModel(new Double(1.0), new Double(0.1), new Double(10.0));

        magnification_spin = new JSpinner(spinner_model);

        JSpinner.NumberEditor magnification_editor = new JSpinner.NumberEditor(magnification_spin, "#%");
        magnification_editor.getTextField().setFont(new Font("SansSerif", Font.PLAIN, 10));
        magnification_editor.getTextField().setBorder(null);
        magnification_editor.getTextField().setEditable(false);
        magnification_editor.getTextField().setBackground(this.getBackground()); // setEditable() may change the background color


        magnification_editor.getTextField().setToolTipText("Current magnification value");
        magnification_spin.setToolTipText("Change magnification");
        magnification_spin.setEditor(magnification_editor);
        magnification_spin.setMaximumSize(new Dimension(80,30));
        
                
        magnification_spin.addChangeListener(new javax.swing.event.ChangeListener() 
            {
                public void stateChanged(javax.swing.event.ChangeEvent E) 
                {
                    double value = ((Number)magnification_spin.getValue()).doubleValue();
                    if (value!=tree_panel.getMagnification())
                    {
                        setMagnification(value);
                    }
                }
            });
        
        bottom_bar.add(magnification_spin);        

    }
            
    /**
     * Sets the magnification value for the tree and flips a bit to indicate 
     * that it needs to be laid out again. 
     * 
     * @param m a magnification level (non-negative, normal size is 1.0) 
     */
    public void setMagnification(double m)
    {
        tree_panel.setMagnification(m);
        valid_bounding_boxes = false;
        repaint();
    }
    
    public void setValidBoundingBoxes(boolean b)
    {
        this.valid_bounding_boxes = b;
    }
    
    private boolean valid_bounding_boxes;
    
    @Override
    public void paintComponent(Graphics g)
    {
        if (!valid_bounding_boxes)
        {
            computeBoundingBoxes(g);
            tree_panel.revalidate();
        }
        
        
        super.paintComponent(g);
    }

    /**
     * Hook for recomputing the bounding boxes of the tree. 
     * Called from within paintComponent, before the actual 
     * painting that is preceded by revalidate(), 
     * so that the tree layout can be recomputed.
     *  
     * @param g graphics context for redrawing the tree
     */
    protected void computeBoundingBoxes(Graphics g)
    {
        tree_panel.computeNodeLabelBoundingBoxes(g);
        tree_panel.computeEdgeLabelBoundingBoxes(g);
        this.valid_bounding_boxes = true;
    }
    
    /**
     * A spinner number model with variable increments:
     * possible values are 1.25,1.5,...,10,12.5,..,100,125,... and 0.9,0.8,...0.1,0.09,0.08,...,0.01,...
     */
    public static class SpinnerMagnifierModel extends SpinnerNumberModel
    {
        public SpinnerMagnifierModel(double value, double minimum, double maximum)
        {
            super(value, minimum, maximum, 1.0);
        }
        
        @Override
        public Object getPreviousValue()
        {
            Number V = (Number)getValue();
            double v = V.doubleValue();
            Number M = (Number)getMinimum();
            double m = M.doubleValue();
            if (m>=v)
                return null;
            double d = 0.1;
            if (v<1.0)
            {
                double z=v;
                while (z<=0.1)
                {
                    z*=10.0;
                    d*=0.1;
                }
            } else
            {
                double z=v;
                while(z>1.0)
                {
                    z*=0.1;
                    d*=10.0;
                }
            }
            if (d>=1.0)
                d/=4.0;
            
            double w = v-d;
            return new Double(w);
        }
        @Override
        public Object getNextValue()
        {
            Number V = (Number)getValue();
            double v = V.doubleValue();
            Number M = (Number)getMaximum();
            double m = M.doubleValue();
            if (m<=v)
                return null;
            double d = 0.1;
            if (v<1.0)
            {
                double z = v;
                while (z<0.1)
                {
                    z*=10.0;
                    d*=0.1;
                }
            } else
            {
                double z=v;
                while(z>=1.0)
                {
                    z*=0.1;
                    d*=10.0;
                }
            }
            if (d>=1.0)
                d/=4.0;
            double w=v+d;
            return new Double(w);
        }
        
    }

}
