package count.gui.kit;

import java.awt.Component;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.net.URL; // for ImageIcons
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.Icon;
import javax.swing.ImageIcon;

/**
 * 
 * Count actions:
 * <ul>
 * <li>Session
 * 	<ol>
 * 	<li>load session (load XML)<li>
 *  <li>+load session by URL</li>
 *  <li>start new session (load tree)</li>
 *  <li>start new session
 *  <li>add tree<li>
 *  <li>edit tree</li>
 *  <li>+ close session</li>
 *  <li>save everything<li>
 *  </ol>
 * </li>
 * <li>Data
 * 	<ol>
 * 	<li>open table, keep annotations</li>
 * 	<li>open table, ignore annotations</li>
 *  <li>+ open table, add star tree</li>
 *  <li>simulate data</li>
 *  <li>extract filtered</li>
 *  <li>convert to binary</li>
 *  </ol>
 * </li>
 * <li>Rates
 *  <ol>
 *  <li>load rates</li>
 *  <li>optimize rates</li>
 *  </ol>
 * </li>
 * <li>Analysis
 *  <ol>
 *  <li>Dollo</li>
 *  <li>numerical parsimony</li>
 *  <li>posteriors</li>
 *  </ol>
 * </li>
 * </ul>
 * 
 */
public class CountActions 
{
	public static final int SIZE_S = 24;
	public static final int SIZE_M = 32;
	public static final int SIZE_L = 48;
	public static final int SIZE_XL = 64;
	public static final int SIZE_XXL = 128;
	
	private static final int PREFERRED_SIZE = SIZE_M;
	
	
	private CountActions() {/* no instantiations */}
	
	private static class Common extends AbstractAction
	{
		Common(String name, Icon icon, ActionListener do_it)
		{
			super(name,icon);
			this.do_it = do_it;
		}
		private final ActionListener do_it;
		
		@Override
		public void actionPerformed(ActionEvent e)
		{
			do_it.actionPerformed(e);
		}
	}
	
	private static class SuperposedIcon implements Icon
	{
		SuperposedIcon(Icon bottom, Icon top)
		{
			this.bottom = bottom;
			this.top = top;
		}
		
		private final Icon bottom;
		private final Icon top;
		
		@Override
		public int getIconWidth()
		{
			return Integer.max(bottom.getIconWidth(), top.getIconWidth());
		}
		
		@Override
		public int getIconHeight()
		{
			return Integer.max(bottom.getIconHeight(), top.getIconHeight());
		}
		
		@Override
		public void paintIcon(Component c, Graphics g, int x, int y)
		{
			bottom.paintIcon(c, g, x, y);
			top.paintIcon(c, g, x, y);
		}
		
	}
	private static ImageIcon createImageIcon(int size, String name)
	{
//		System.out.println("#**CA.cII loading image "+name);
        URL icon_url = ClassLoader.getSystemResource("img/"+name);
		
		ImageIcon createImageIcon = new ImageIcon(icon_url);
		if (size != createImageIcon.getIconHeight() || size != createImageIcon.getIconWidth())
		{
			Image image = createImageIcon.getImage();
			createImageIcon = new ImageIcon(image.getScaledInstance(size, size, Image.SCALE_SMOOTH));
		}
		return createImageIcon;
	}
	
	public static Icon createOpenIcon(int size)
	{
		return createImageIcon(size, "Open24.gif");
	}
	
	public static Icon createLoadSessionIcon(int size)
	{
		ImageIcon bottom_icon = createImageIcon(size, "Open24.gif");
//		TreeIcon top_icon = new TreeIcon(size); // square icon
//		return new SuperposedIcon(bottom_icon, top_icon);
		return bottom_icon;
	}
	public static Action createLoadSession(String name, ActionListener do_it)
	{
		return new Common(name,createLoadSessionIcon(PREFERRED_SIZE), do_it);
	}	
	public static Icon createNewSessionIcon(int size)
	{
		ImageIcon bottom_icon = createImageIcon(size, "New24.gif");
		TreeIcon top_icon = new TreeIcon(size); // square icon
		return new SuperposedIcon(bottom_icon, top_icon);
	}
	public static Action createNewSession(String name, ActionListener do_it)
	{
		return new Common(name, createNewSessionIcon(PREFERRED_SIZE), do_it);
	}
	
	public static Icon createInitSessionIcon(int size)
	{
		ImageIcon bottom_icon = createImageIcon(size, "New24.gif");
		return StringIcon.superposedStringIcon(bottom_icon, "?");
	}
	
	public static Action createInitSession(String name, ActionListener do_it)
	{
		return new Common(name, createInitSessionIcon(PREFERRED_SIZE), do_it);
	}
	
	public static Icon createLoadTreeIcon(int size)
	{
		ImageIcon bottom_icon = createImageIcon(size, "Add24.gif");
		TreeIcon top_icon = new TreeIcon(size); // square icon
		return new SuperposedIcon(bottom_icon, top_icon);
	}
	public static Action createLoadTree(String name, ActionListener do_it)
	{
		return new Common(name, createLoadTreeIcon(PREFERRED_SIZE), do_it);
	}
	
	public static Icon createBuildTreeIcon(int size)
	{
		ImageIcon bottom_icon = createImageIcon(size, "New24.gif");
		TreeIcon top_icon = new TreeIcon(size); // square icon
		return new SuperposedIcon(bottom_icon, top_icon);
	}
	
	public static Action createBuildTree(String name, ActionListener do_it)
	{
		return new Common(name, createBuildTreeIcon(PREFERRED_SIZE), do_it);
	}
	
	
	
	public static Icon createEditTreeIcon(int size)
	{
//		ImageIcon bottom_icon = createImageIcon(size, "Cut24.gif");
		TreeIcon top_icon = new TreeIcon(size);
//		return new SuperposedIcon(bottom_icon, top_icon);
		return StringIcon.indexingStringIcon(top_icon, "\u2703"); // scissors
	}
	public static Action createEditTree(String name, ActionListener do_it)
	{
		return new Common(name, createEditTreeIcon(PREFERRED_SIZE), do_it);
	}

	public static Icon createLoadTableIcon(int size)
	{
		TableIcon bottom_icon = new TableIcon(size, false);
		return StringIcon.indexingStringIcon(bottom_icon, "+");
	}
	public static Action createLoadTable(String name, ActionListener do_it)
	{
		return new Common(name, createLoadTableIcon(PREFERRED_SIZE), do_it);
	}
	public static Icon createLoadAnnotationsIcon(int size)
	{
		TableIcon bottom_icon =TableIcon.columnFilterIcon(size);
		return StringIcon.indexingStringIcon(bottom_icon, "+");
	}
	public static Action createLoadAnnotations(String name, ActionListener do_it)
	{
		return new Common(name, createLoadAnnotationsIcon(PREFERRED_SIZE), do_it);
	}
	
	public static Icon createImportTableIcon(int size)
	{
		ImageIcon bottom_icon = createImageIcon(size, "Import24.gif");
		return bottom_icon;
	}
	public static Action createImportTable(String name, ActionListener do_it)
	{
		return new Common(name, createImportTableIcon(PREFERRED_SIZE), do_it);
	}
	
	
	public static Icon createSimulationIcon(int size)
	{
		TableIcon bottom_icon = new TableIcon(size, false);
		TreeIcon top_icon = TreeIcon.ratesIcon(size);
		Icon superposed = new SuperposedIcon(bottom_icon, top_icon);
		return StringIcon.indexingStringIcon(superposed, "✴"); // \u2734e
	}
	
	
	
	
	
	public static Action createSimulation(String name, ActionListener do_it)
	{
		return new Common(name, createSimulationIcon(PREFERRED_SIZE), do_it);
	}
	public static Icon createFilterRowsIcon(int size)
	{
		TableIcon bottom_icon =  new TableIcon(size, true);
		return StringIcon.indexingStringIcon(bottom_icon, "+");
	}
	public static Action createFilterRows(String name, ActionListener do_it)
	{
		return new Common(name, createFilterRowsIcon(PREFERRED_SIZE), do_it);
	}
	public static Icon createBinaryTableIcon(int size)
	{
		TableIcon bottom_icon = new TableIcon(size, false);
		bottom_icon.setBinary(true);
		return StringIcon.indexingStringIcon(bottom_icon, "+");
	}
	public static Action createBinaryTable(String name, ActionListener do_it)
	{
		return new Common(name, createBinaryTableIcon(PREFERRED_SIZE), do_it);
	}

			
	
	public static Icon createLoadRatesIcon(int size)
	{
//		ImageIcon bottom_icon = createImageIcon(size, "Add24.gif");
		TreeIcon top_icon = TreeIcon.ratesIcon(size); // square icon
//		return new SuperposedIcon(bottom_icon, top_icon);
		return StringIcon.indexingStringIcon(top_icon, "+");
	}
	public static Action createLoadRates(String name, ActionListener do_it)
	{
		return new Common(name, createLoadRatesIcon(PREFERRED_SIZE), do_it);
	}
	public static Icon createOptimizeRatesIcon(int size)
	{
//		ImageIcon bottom_icon = createImageIcon(size, "New24.gif");
		TreeIcon top_icon = TreeIcon.ratesIcon(size); // square icon
//		return new SuperposedIcon(bottom_icon, top_icon);
		return StringIcon.indexingStringIcon(top_icon, "✴");
	}
	public static Action createOptimizeRates(String name, ActionListener do_it)
	{
		return new Common(name, createOptimizeRatesIcon(PREFERRED_SIZE), do_it);
	}
	
	public static Icon createSaveAllIcon(int size)
	{
		ImageIcon createSaveAllIcon = createImageIcon(size, "SaveAll24.gif");
		return createSaveAllIcon;
	}
	public static Action createSaveAll(String name, ActionListener do_it)
	{
		return new Common(name, createSaveAllIcon(PREFERRED_SIZE), do_it);
	}
	public static Icon createExportIcon(int size)
	{
		ImageIcon createExportIcon = createImageIcon(size, "Export24.gif");
		return createExportIcon;
	}
	public static Action createExport(String name, ActionListener do_it)
	{
		return new Common(name, createExportIcon(PREFERRED_SIZE), do_it);
	}
	public static Icon createSaveIcon(int size)
	{
		ImageIcon createSaveIcon = createImageIcon(size, "Save24.gif");
		return createSaveIcon;
	}
	public static Action createSave(String name, ActionListener do_it)
	{
		return new Common(name, createSaveIcon(PREFERRED_SIZE), do_it);
	}
	public static Icon createDeleteIcon(int size)
	{
		ImageIcon createDeleteIcon = createImageIcon(size, "Delete24.gif");
		return createDeleteIcon;
	}
	public static Action createDelete(String name, ActionListener do_it)
	{
		return  new Common(name, createDeleteIcon(PREFERRED_SIZE), do_it); 
	}
	public static Icon createRemoveIcon(int size)
	{
		ImageIcon createRemoveIcon = createImageIcon(size, "Remove24.gif");
		return createRemoveIcon;
	}
	public static Action createRemove(String name, ActionListener do_it)
	{
		return  new Common(name, createRemoveIcon(PREFERRED_SIZE), do_it);
	}
	public static Icon createPosteriorsIcon(int size)
	{
		TableIcon bottom_icon = TableIcon.historyIcon(size);
		TreeIcon top_icon = TreeIcon.ratesIcon(size);
		Icon superposed = new SuperposedIcon(bottom_icon, top_icon);
		
		return superposed;
		
//		return StringIcon.indexingStringIcon(superposed, "+");
	}
	public static Action createPosteriors(String name, ActionListener do_it)
	{
		return new Common(name, createPosteriorsIcon(PREFERRED_SIZE), do_it);
	}
	public static Icon createParsimonyIcon(int size)
	{
		TableIcon bottom_icon = TableIcon.historyIcon(size);
		TreeIcon top_icon = new TreeIcon(size);
		Icon superposed = new SuperposedIcon(bottom_icon, top_icon);
		
		return superposed;
//		
//		return StringIcon.indexingStringIcon(new SuperposedIcon(bottom_icon, top_icon), "+");
	}
	public static Action createParsimony(String name, ActionListener do_it)
	{
		return new Common(name, createParsimonyIcon(PREFERRED_SIZE), do_it);
	}
	public static Icon createDolloIcon(int size)
	{
		TableIcon bottom_icon = TableIcon.historyIcon(size);
		bottom_icon.setBinary(true);
		TreeIcon top_icon = new TreeIcon(size);
		Icon superposed = new SuperposedIcon(bottom_icon, top_icon);
		
		return superposed;
//		return StringIcon.indexingStringIcon(new SuperposedIcon(bottom_icon, top_icon), "+");
	}
	public static Action createDollo(String name, ActionListener do_it)
	{
		return new Common(name, createDolloIcon(PREFERRED_SIZE), do_it);
	}
	
	public static Icon createQuitIcon(int size)
	{
		ImageIcon stop_icon = createImageIcon(size, "Stop24.gif");
		return stop_icon;
	}
	
	public static Action createQuit(String name, ActionListener do_it)
	{
		return new Common(name, createQuitIcon(PREFERRED_SIZE), do_it);
	}
	
	public static Icon createAboutIcon(int size)
	{
		ImageIcon about_icon = createImageIcon(size, "About24.gif");
		return about_icon;
	}
	
	public static Action createAbout(String name, ActionListener do_it)
	{
		return new Common(name, createAboutIcon(PREFERRED_SIZE), do_it);
	}
}


