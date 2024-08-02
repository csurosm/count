package count.gui;
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

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.SwingWorker;
import javax.swing.Timer;

import count.gui.kit.CounterBox;
import count.gui.kit.Grapher;
import count.io.DataFile;
import count.model.GammaInvariant;
import count.model.MixedRateModel;
import count.model.RateVariationModel;
import count.model.old.MLGamma;
import count.model.ML;

/**
 * Panel for tracking model optimization. 
 * 
 * @author csuros
 *
 */
public class RateOptimizationPanel extends RateVariationPanel 
{
    private static final int HISTORY_LENGTH = 50;
    private static final int GRAPH_WIDTH = 200;
    private static final int GRAPH_HEIGHT = 72;
    private static final int STAGE_WIDTH = 300;
    private static final int COUNTER_WIDTH = 120;
    private static final int COUNTER_HEIGHT = 48;
	
	public static final String PROPERTY_OPTIMIZATION_VALUE = "lnL"; 
    public static final String PROPERTY_OPTIMIZATION_ROUND = "round";
    public static final String PROPERTY_OPTIMIZATION_DROP = "delta";
    public static final String PROPERTY_OPTIMIZATION_DONE = "optimization.done";
    private static final String ACTION_SNAPSHOT = "snapshot";
    private static final String ACTION_CANCEL = "cancel";
    private static final String ACTION_LIVE_UPDATE = "gyia";
    
    private static final Color GLOW_COLOR = new Color(1.0f,0.8f,0.1f,0.84f);
    private static final Color GLOW_FROZEN_COLOR = new Color(0.09f, 0.03f, 0.5f, GLOW_COLOR.getAlpha()/256f);
    

    public RateOptimizationPanel(DataFile<? extends MixedRateModel> model_data, ML optimization)
	{
		super(model_data);
		this.factory = optimization;
		this.is_descendant_model = false;
	}

	private ML factory;
	
	private OptimizationTask optimization_task;
	
	private boolean is_descendant_model;
	
	public boolean isDescendantModel() { return this.is_descendant_model;}
	public void setDescendantModel(boolean is_descendant_model) { this.is_descendant_model = is_descendant_model;}
	
	
	private static class OptimizationStep
	{
		OptimizationStep(int step, double value, double delta)
		{
			this.step = step;
			this.value = value;
			this.delta = delta;
		}
		final int step;
		final double value;
		final double delta;
		@Override
		public String toString()
		{
			return "["+step+","+value+","+delta+"]";
		}
	}
	
	
	@Override
    public int getModelParameterCount()
    {
    	return factory.getModelParameterCount();
    }
	
	
	class OptimizationTask extends SwingWorker<Double, OptimizationStep>
		implements ActionListener
	{
		OptimizationTask(AppFrame app, double delta, int max_iter)
		{
			this.delta = delta;
			this.max_iter = max_iter;
			this.app = app;
			this.optimization_property_support = new PropertyChangeSupport(this);
			initComponents();
		}
		private final double delta;
		private final int max_iter;
		private final AppFrame app;
		Thread worker_thread;

		// control & feedback
		JButton cancel_button ;
	    JButton snapshot_button;
	    Timer optimization_timer;
		JLabel track_stage;
		Grapher track_convergence;
		Grapher track_drop;
		CounterBox track_rounds;
		
		private final PropertyChangeSupport optimization_property_support; // DIY 
	    /**
	     * Handles actions coming from the snapshot and cancel buttons.
	     * 
	     * @param E event triggering this call
	     */
	    @Override
	    public void actionPerformed(ActionEvent E)
	    {
	        String cmd = E.getActionCommand();
	        if (ACTION_CANCEL.equals(cmd))
	        {
//	        	worker_thread.interrupt();
	        	this.cancel(true);
	        	updateDisplay(true);

	        } else if (ACTION_SNAPSHOT.equals(cmd))
	        {
	            snapshot();
	        } else if (ACTION_LIVE_UPDATE.equals(cmd))
	        {
	            updateDisplay(false);
	        }
	    }
		
		private void initComponents()
		{
			List<Double> history = new ArrayList<>()
			{
				@Override
				public boolean add(Double val)
				{
					double delta = isEmpty()?Double.NaN:(get(size()-1)-val);  // function is being minimized
					OptimizationStep step = new OptimizationStep(size(), val, delta);
					publish(step);
					return super.add(val);
				}
			};
			factory.setOptimizationHistory(history);
			

	        // create the buttons
	        snapshot_button = new JButton("Snapshot");
	        snapshot_button.setToolTipText("Takes a snapshot of the current model");
	        snapshot_button.setActionCommand(ACTION_SNAPSHOT); 
	        snapshot_button.addActionListener(this);

	        cancel_button = new JButton("Cancel");
	        cancel_button.setToolTipText("Stops the optimization");
	        cancel_button.setActionCommand(ACTION_CANCEL);
	        cancel_button.addActionListener(this);

			// timer for live updates
	        optimization_timer = new Timer(100, this); // start after 0.1 seconds
	        optimization_timer.setDelay(100); // update every 0.1 seconds
	        optimization_timer.setActionCommand(ACTION_LIVE_UPDATE);
	        
		}
		
		Box createTrackerBox()
		{
			Box tracker_box = new Box(BoxLayout.LINE_AXIS);
			tracker_box.setBorder(BorderFactory.createEtchedBorder());
			
			int font_size = TreePanel.getTreeNormalFontSize();
	        Font tp_font_rm = new Font("Serif", Font.PLAIN, font_size);
	        Font tp_font_it = new Font("Serif", Font.ITALIC, 7*font_size/6);

	        track_stage = new JLabel("Optimizing model parameters ("+factory.getClass().getSimpleName()+")");
	        track_stage.setFont(tp_font_it);
	        //track_stage.setBackground(Color.WHITE);   //new Color(255,255,192)); // pale yellow
	        track_stage.setOpaque(true);
	        track_stage.setMinimumSize(new Dimension(STAGE_WIDTH,GRAPH_HEIGHT));
	        track_stage.setPreferredSize(track_stage.getMinimumSize());
	        track_stage.setMaximumSize(track_stage.getMinimumSize());
	        //track_stage.setBorder(BorderFactory.createLoweredBevelBorder());
	        track_stage.setToolTipText("Current optimization status");
	        
	        tracker_box.add(track_stage);
	        
	        // track the optimization stage - not useful
	        factory.addPropertyChangeListener(ML.PROPERTY_OPTIMIZATION_PHASE, 
	        		e->{
	        			String phase;
	        			if (e.getNewValue()==null)
	        			{
	        				phase = "";
	        			} else
	        			{
	        				phase = e.getNewValue().toString();
	        			}
	        			track_stage.setText(phase);
	        		});
	        
	        track_rounds = new CounterBox("Round ", Integer.toString(max_iter));
	        track_rounds.setFont(tp_font_rm);
	        track_rounds.setMinimumSize(new Dimension(COUNTER_WIDTH,COUNTER_HEIGHT));
	        track_rounds.setPreferredSize(track_rounds.getMinimumSize());
	        track_rounds.setMaximumSize(track_rounds.getPreferredSize());
	        track_rounds.setTrackersOpaque(true);
	        track_rounds.setOpaque(true);
	        //track_rounds.setBackground(Color.WHITE); //new Color(255,192,160)); // pale orange
	        track_rounds.setToolTipText("Current iteration round"); 
	        tracker_box.add(track_rounds);
	        optimization_property_support.addPropertyChangeListener(PROPERTY_OPTIMIZATION_ROUND, track_rounds);
	        
	        track_convergence = new Grapher(GRAPH_WIDTH, GRAPH_HEIGHT);
	        track_convergence.setFont(tp_font_rm);
	        track_convergence.setHistoryLength(HISTORY_LENGTH);
	        track_convergence.setDrawLine(true);
	        track_convergence.setDrawDots(true);
	        track_convergence.setDrawBars(false);
	        track_convergence.setToolTipText("Log-likelihood value");
	        track_convergence.glow(GLOW_COLOR, 1.4); // 84 beats/min
	        tracker_box.add(track_convergence);
	        optimization_property_support.addPropertyChangeListener(PROPERTY_OPTIMIZATION_VALUE, track_convergence);
	        
			track_drop = new Grapher(GRAPH_WIDTH, GRAPH_HEIGHT);
	        track_drop.setFont(track_convergence.getFont());
	        track_drop.setHistoryLength(HISTORY_LENGTH);
	        track_drop.setDrawLine(false);
	        track_drop.setDrawDots(true);
	        track_drop.setDrawBars(true);
	        track_drop.setToolTipText("Increase of log-likelihood in the most recent optimization round");
	        track_drop.glow(GLOW_COLOR, 1.4);
			tracker_box.add(track_drop);
			optimization_property_support.addPropertyChangeListener(PROPERTY_OPTIMIZATION_DROP, track_drop);
			
			
			
			return tracker_box;
		}
	    
		@Override 
		public Double doInBackground()
		{
			this.worker_thread = Thread.currentThread();
			optimization_timer.start();	
//			System.out.println("#**GOP.OT.dIB start on "+worker_thread);
			try
			{
				Double X=factory.optimize(delta, max_iter);
				
				return X;
			} catch (Throwable t)
			{
				t.printStackTrace();
				throw new RuntimeException(t);
			}
		}
		
		@Override
		protected void process(List<OptimizationStep> steps)
		{
			for (OptimizationStep s:steps)
			{
				int round = s.step;
				double value = s.value;
				double delta = s.delta;
				
				if (!Double.isNaN(delta) && delta<0.)
				{
//					System.out.println("#*ROP.OT.process delta "+delta+"\tval "+value);
				} else
				{
					optimization_property_support.firePropertyChange(PROPERTY_OPTIMIZATION_ROUND, round-1, round);
					optimization_property_support.fireIndexedPropertyChange(PROPERTY_OPTIMIZATION_VALUE, round, null, -value);
					if (!Double.isNaN(delta))
						optimization_property_support.fireIndexedPropertyChange(PROPERTY_OPTIMIZATION_DROP, round, null, delta);
				}
			}
		}
		
		
		void publish(OptimizationStep x) // make method accessible from history
		{
			super.publish(x);
		}
		
		@Override
		protected void done()
		{
			optimization_timer.stop();
			snapshot_button.setEnabled(false);
			snapshot_button.setVisible(false);
			cancel_button.setEnabled(false);
			cancel_button.setVisible(false);
			track_stage.setText("Optimization "+(isCancelled()?"canceled prematurely.":"done."));
			track_convergence.glow(GLOW_FROZEN_COLOR, 0.0);
			track_drop.glow(GLOW_FROZEN_COLOR, 0.0);
			updateDisplay(true); // recompute layout
			RateOptimizationPanel.this.firePropertyChange(PROPERTY_OPTIMIZATION_DONE, false, true);
		}

		
	}		

	/**
	 * Main entrance point 
	 * 
	 * @param app
	 * @param delta
	 * @param max_iter
	 */
	public void launchOptimization(AppFrame app, double delta, int max_iter)
	{
		
		optimization_task = new OptimizationTask(app, delta, max_iter);

		// set control bar 
		Box tracker_box = optimization_task.createTrackerBox();
		Zoom<RatesTreePanel> ratesP = getRatesPanel();
		Box control_bar = ratesP.getControlBar();
		Component[] control_widgets = control_bar.getComponents(); // superclasses added them
		control_bar.removeAll();
		control_bar.add(tracker_box);
		control_bar.add(Box.createHorizontalGlue());
		control_bar.add(optimization_task.snapshot_button);
		control_bar.add(optimization_task.cancel_button);
		control_bar.add(Box.createHorizontalGlue());
		for (Component w: control_widgets)
			control_bar.add(w);
		
		optimization_task.execute();
		repaint();
	}
    
    private void updateDisplay(boolean recompute_layout)
    {
    	Zoom<RatesTreePanel> ratesZ = getRatesPanel();
    	RatesTreePanel ratesP = ratesZ.getTreePanel();
    	if (recompute_layout)
    		ratesP.calculateNodeLocations();
    	else
    		ratesP.setValidBoundingBoxes(false);
    	
    	RatesTable ratesT = getRatesTable();
    	int selected_node = ratesT.getSelectedModelRow(); // fireTableDataChanged clears row selection
    	((RatesTable.Model)ratesT.getModel()).fireTableDataChanged();
    	if (selected_node != -1) ratesT.setSelectedModelRow(selected_node);
    	repaint();
    		
    }
    
    void snapshot()
    {
	    int round = factory.getOptimizationHistory().size();
	    DataFile<MixedRateModel> current_file = this.getDataFile();
	    MixedRateModel current_rates =  current_file.getContent(); 
	    MixedRateModel copy_rates;
	    if (current_rates instanceof GammaInvariant)
	    {
	    	copy_rates=((GammaInvariant)current_rates).copy();
	    } else
	    {
	    	assert (current_rates instanceof RateVariationModel);
	    	copy_rates = ((RateVariationModel)current_rates).copy();
	    }
	
	    File snapshot_file = new File((File)null,DataFile.chopFileExtension(current_file.getFile().getName())+"#"+Integer.toString(round));
	    DataFile<MixedRateModel> copy_file = new DataFile<>(copy_rates, snapshot_file);
	    Session our_sesh = Session.getSession(this);
	    our_sesh.addRates(copy_file, false);
	    
	    updateDisplay(true);
    }
}
	
