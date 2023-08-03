package count.gui.kit;

import java.awt.Dimension;

import javax.swing.DefaultBoundedRangeModel;
import javax.swing.JComponent;
import javax.swing.event.ChangeListener;

public class AsynchronousProgressBar extends JComponent
{
	private static final int DEFAULT_WIDTH = 500;
	private static final int DEFAULT_HEIGHT = 30;
	
	public AsynchronousProgressBar(int maximum, int pixel_height)
	{
		super();
		this.height = height;
		this.range_model = new DefaultBoundedRangeModel();
		setMaximum(maximum);
	}
	public AsynchronousProgressBar(int maximum)
	{
		this(maximum,Integer.max(DEFAULT_HEIGHT,1+(maximum-1)/DEFAULT_WIDTH));
	}
	
	private int height;
	private boolean[] done;

	private final DefaultBoundedRangeModel range_model;
	
	public void setMaximum(int maximum)
	{
		assert (maximum>=0);
		range_model.setValue(0);
		range_model.setMaximum(maximum);
		done = new boolean[maximum];
	}
	
	public int getMaximum()
	{
		return range_model.getMaximum();
	}
	
	public int getValue()
	{
		return range_model.getValue();
	}
	
	public void done(int n)
	{
		if (!done[n])
		{
			done[n]=true;
			range_model.setValue(range_model.getValue()+1);
		}
	}
	
	public void addChangeListener(ChangeListener listener)
	{
		range_model.addChangeListener(listener);
	}
	
	public void removeChangeListener(ChangeListener listener)
	{
		range_model.removeChangeListener(listener);
	}
	
	@Override
	public Dimension getPreferredSize()
	{
		return new Dimension(1+(range_model.getMaximum()-1)/height,height);
	}
	
	
	
	
}
