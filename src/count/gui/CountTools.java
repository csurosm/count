package count.gui;
/*
 * Copyright 2023 Mikl&oacute;s Cs&#369;r&ouml;s.
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

import java.awt.Component;
import java.awt.event.ActionListener;
import java.util.function.Predicate;
import java.util.Map;
import java.util.HashMap;
import javax.swing.AbstractButton;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JToolBar;

import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;

import count.gui.BundleTree.Node;
import count.gui.kit.CountActions;
import count.io.SavableData;
import count.io.ExportableData;
import count.io.Removable;
import count.io.ModelBundle;

/**
 * GUI component (JToolBar) for Count's actions. 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */

public class CountTools extends JToolBar implements ChangeListener
{
	public CountTools(AppFrame app, String name)
	{
		super(name);
		this.app = app;
		this.action_enabler = new HashMap<>();
	}
	
	private final AppFrame app;
	private final Map<Action, Predicate<AppFrame>> action_enabler;
	
	public static final Predicate<AppFrame> HAS_ACTIVE_SESSION = app->(app.getActiveSession() != null);
	public static final Predicate<AppFrame> HAS_NO_SESSIONS = app->(app.getSessionCount()==0);
	public static final Predicate<AppFrame> HAS_SELECTED_TABLE = new Predicate<>() {
		@Override
		public boolean test(AppFrame app)
		{
			Session sesh = app.getActiveSession();
			ModelBundle.Entry table_entry = (sesh==null)?null:sesh.getDataBrowser().getSelectedTableEntry();				
			return table_entry != null;
		}
	};
	public static final Predicate<AppFrame> HAS_SELECTED_RATES = new Predicate<>() {
		@Override
		public boolean test(AppFrame app)
		{
			Session sesh = app.getActiveSession();
            ModelBundle.Entry rates_entry 
			= (sesh==null)
			? null
			: sesh.getModelBrowser().getSelectedRatesEntry();
            return rates_entry != null;
		}
	};
	public static final Predicate<AppFrame> HAS_SELECTED_BINARY = new Predicate<>() 
	{
		@Override
		public boolean test(AppFrame app)
		{
			Session sesh = app.getActiveSession();
			ModelBundle.Entry table_entry = (sesh==null)?null:sesh.getDataBrowser().getSelectedTableEntry();				
			return table_entry != null && table_entry.getTableData().getContent().isBinaryTable();
		}
	};
	public static final Predicate<AppFrame> HAS_FILTERED_FAMILIES = new Predicate<>()
	{
		@Override
		public boolean test(AppFrame app)
		{
			Session sesh = app.getActiveSession();
			Object selected = sesh==null?null:sesh.getDataBrowser().getSelectedItem();
			boolean test = selected != null && selected instanceof Session.FamilySelection;
			return test;
		}
	};
	public static final Predicate<AppFrame> IS_SELECTED_SAVABLE = app->(app.getActiveSession()!=null && app.getActiveSession().getSelectedComponent() instanceof SavableData);
	public static final Predicate<AppFrame> IS_SELECTED_EXPORTABLE = app->(app.getActiveSession()!=null && app.getActiveSession().getSelectedComponent() instanceof ExportableData);
	public static final Predicate<AppFrame> IS_SELECTED_REMOVABLE = app->(app.getActiveSession()!=null && app.getActiveSession().isSelectedRemovable());
	
	private JButton addAction(Action action, Predicate<AppFrame> enabler)
	{
		JButton button = this.add(action);
		action_enabler.put(action, enabler);
		button.setEnabled(enabler.test(app));
		button.setHideActionText(false);
//		button.setBorder(BorderFactory.createRaisedSoftBevelBorder());
		button.setBorderPainted(button.isEnabled());
		return button;
	}
	
	
	public JButton addLoadSession(String name, ActionListener do_it)
	{
		Action loadSession = CountActions.createLoadSession(name, do_it);
		return addAction(loadSession, HAS_NO_SESSIONS);
	}
	
	public JButton addNewSession(String name, ActionListener do_it)
	{
		Action newSession = CountActions.createNewSession(name, do_it);
		return addAction(newSession, always->true);
	}
	public JButton addInitSession(String name, ActionListener do_it)
	{
		Action initSession = CountActions.createInitSession(name, do_it);
		return addAction(initSession, always->true);
	}
	
	public JButton addSaveSession(String name, ActionListener do_it)
	{
		Action saveSession = CountActions.createSaveAll(name, do_it);
		return addAction(saveSession, HAS_ACTIVE_SESSION);
	}
	
	public JButton addCloseSession(String name, ActionListener do_it)
	{
		Action closeSession = CountActions.createDelete(name, do_it);
		return addAction(closeSession, HAS_ACTIVE_SESSION);
	}
	
	public JButton addLoadTree(String name, ActionListener do_it)
	{
		Action loadTree = CountActions.createLoadTree(name, do_it);
		return addAction(loadTree, HAS_ACTIVE_SESSION);
	}
	public JButton addBuildTree(String name, ActionListener do_it)
	{
		Action buildTree = CountActions.createBuildTree(name, do_it);
		return addAction(buildTree,HAS_SELECTED_TABLE);
	}
	
	public JButton addEditTree(String name, ActionListener do_it)
	{
		Action editTree = CountActions.createEditTree(name, do_it);
		return addAction(editTree, HAS_ACTIVE_SESSION);
	}
	
	
	public JButton addLoadTable(String name, ActionListener do_it)
	{
		Action loadTable = CountActions.createLoadTable(name, do_it);
		return addAction(loadTable, HAS_ACTIVE_SESSION);
	}
	
	public JButton addLoadAnnotations(String name, ActionListener do_it)
	{
		Action loadAnnotations = CountActions.createLoadAnnotations(name, do_it);
		return addAction(loadAnnotations, HAS_ACTIVE_SESSION);
	}

	public JButton addImportTable(String name, ActionListener do_it)
	{
		Action importTable = CountActions.createImportTable(name, do_it);
		return addAction(importTable, always->true);
	}
	
	
	public JButton addSimulation(String name, ActionListener do_it)
	{
		Action simulation = CountActions.createSimulation(name, do_it);
		return addAction(simulation, HAS_SELECTED_TABLE.and(HAS_SELECTED_RATES));
	}
	
	public JButton addFilterRows(String name, ActionListener do_it)
	{
		Action filterRows = CountActions.createFilterRows(name, do_it);
		return addAction(filterRows, HAS_FILTERED_FAMILIES);
	}
	
	public JButton addBinaryTable(String name, ActionListener do_it)
	{
		Action binaryTable = CountActions.createBinaryTable(name, do_it);
		return addAction(binaryTable, HAS_SELECTED_TABLE.and(Predicate.not(HAS_SELECTED_BINARY)));
	}
	
	public JButton addLoadRates(String name, ActionListener do_it)
	{
		Action loadRates = CountActions.createLoadRates(name, do_it);
		return addAction(loadRates, HAS_ACTIVE_SESSION);
	}
	
	public JButton addOptimizeRates(String name, ActionListener do_it)
	{
		Action optimizeRates = CountActions.createOptimizeRates(name, do_it);
		return addAction(optimizeRates, HAS_SELECTED_TABLE.and(Predicate.not(HAS_SELECTED_BINARY)));
	}
	public JButton addDollo(String name, ActionListener do_it)
	{
		Action dollo = CountActions.createDollo(name, do_it);
		return addAction(dollo, HAS_SELECTED_BINARY);
	}
	public JButton addParsimony(String name, ActionListener do_it)
	{
		Action parsimony = CountActions.createParsimony(name, do_it);
		return addAction(parsimony, HAS_SELECTED_TABLE.and(Predicate.not(HAS_SELECTED_BINARY)));
	}
	public JButton addPosteriors(String name, ActionListener do_it)
	{
		Action posteriors = CountActions.createPosteriors(name, do_it);
		return addAction(posteriors, HAS_SELECTED_TABLE.and(Predicate.not(HAS_SELECTED_BINARY)).and(HAS_SELECTED_RATES));
	}
	public JButton addRemove(String name, ActionListener do_it)
	{
		Action remove = CountActions.createRemove(name, do_it);
		return addAction(remove, IS_SELECTED_REMOVABLE);
	}
	public JButton addSave(String name, ActionListener do_it)
	{
		Action save = CountActions.createSave(name, do_it);
		return addAction(save, IS_SELECTED_SAVABLE.or(IS_SELECTED_EXPORTABLE));
	}
	
	public JButton addQuit(String name, ActionListener do_it)
	{
		Action quit = CountActions.createQuit(name, do_it);
		return addAction(quit, always->true);
	}
	
	/**
	 * Enables/disables all the added buttons using the associated rules. 
	 * 
	 */
	public void updateEnabledStates()
	{
		for (Component c: this.getComponents())
		{
			if (c instanceof AbstractButton)
			{
				AbstractButton button = (AbstractButton) c;
				Action a = button.getAction();
				button.setEnabled(action_enabler.get(a).test(app));
//				button.setBorder(BorderFactory.createRaisedSoftBevelBorder());
				button.setBorderPainted(button.isEnabled());
//				System.out.println("#**CT.uES "+button);
			}
		}
		
//		for (Map.Entry<Action, Predicate<AppFrame>> entry: action_enabler.entrySet())
//		{
//			entry.getKey().setEnabled();
//		}
	}

	/**
	 * Buttons are enabled/disabled when something changes in the session browsers.
	 */
	@Override
	public void stateChanged(ChangeEvent e) 
	{
		updateEnabledStates();
		this.repaint();
	}
}
