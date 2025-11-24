package count.gui.kit;

import java.awt.Dimension;

import javax.swing.JComboBox;

/**
 * A simple combo box for selecting within a range of 
 * integers.
 */
public class SmallIntegerSelector extends JComboBox<Integer> {
	public SmallIntegerSelector(int min, int max) {
		super(range(min,max));
		this.setSelectedItem(max);
		this.setEditable(false);
		this.setMaximumSize(new Dimension(80,24)); // do not stretch horizontally
		
	}
	
	public int getSelectedValue() {
		return getItemAt(getSelectedIndex()).intValue();
	}
	
	private static Integer[] range(int min, int max) {
		Integer[] range = new Integer[max-min+1];
		for (int i=0; i<range.length; i++)
			range[i]=Integer.valueOf(min+i);
		return range;
	}
	
	
}
