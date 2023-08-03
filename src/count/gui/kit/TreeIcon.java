package count.gui.kit;
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

import java.awt.BasicStroke;

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

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;

import java.awt.geom.Point2D;


import javax.swing.Icon;

import count.gui.HistoryView;
import count.gui.RatesTreePanel;
import count.gui.TreePanel;

/**
 * An Icon for phylogeny or rate parameters.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s 
 *
 */
public class TreeIcon implements Icon
{
	private static final int SIZE = 24;
	private static final Color BACKGROUND_COLOR = Color.WHITE;
	
	public TreeIcon()
	{
		this(SIZE);
	}
	
	public TreeIcon(int size)
	{
		this(size, TreePanel.TREE_EDGE_COLOR, BACKGROUND_COLOR);
	}
	public TreeIcon(int size, Color color, Color background)
	{
		this(size, color, null, background);
	}

	public TreeIcon(int size, Color[] colors, Color background)
	{
		this(size,colors[0],colors, background);
	}
	
	private TreeIcon(int size, Color color, Color[] colors, Color background)
	{
		this.size = size;
		this.single_color = color;
		this.edge_colors = colors;
		initPoints();
	}
	
	
	private final int size;
    @Override
    public int getIconHeight()
    {
        return size;
    }

    @Override
    public int getIconWidth()
    {
        return size;
    }
	
	private Point2D[] tree_nodes;
	private int[] parent;
	
	private final Color[] edge_colors;
	private final Color single_color;
	
	private void initPoints()
	{
		int num_leaves = 4;

		int num_nodes = 2*num_leaves-1;
		parent = new int[num_nodes];
		
		int dw = size/(1+num_leaves);
		int dh = size/4; // since height==3
		
		tree_nodes = new Point2D[num_nodes];
		int node = 0;
		while (node<num_leaves)
		{
			int px = (node+1)*dw+dw/3;
			int py = dh;
			tree_nodes[node]  = new Point2D.Double(px,py);
			parent[node] = num_leaves + node/2;
			node++;
		}
		for (int i=0; i<2; i++) // 2 parents
		{
			parent[node] = 6;
			double px = dw * (1.5+2.0*i); // 1.5, 3.5
			double py = 2*dh;
			tree_nodes[node] = new Point2D.Double(px,py);
			node++;
		}
		int root = node++;
		parent[root] = -1;
		double rx = 2.5*dw;
		double ry = 3*dh;
		tree_nodes[root] = new Point2D.Double(rx,ry);
	}
	
	private void drawTree(Graphics2D g2, boolean curved)
	{
		for (int node=0; node<tree_nodes.length; node++)
		{
			int p = parent[node];
			if (p!=-1)
			{
				Point2D Pchild = tree_nodes[node];
				Point2D Pparent = tree_nodes[p];
				if (curved)
					TreePanel.drawCurvedLine(g2, Pchild, Pparent, Pchild.getY());
				else
					TreePanel.drawBentLineVertical(g2, Pchild, Pparent);
			}
		}
	}
	
    @Override
    public void paintIcon(Component C, Graphics g, int x, int y)
    {
    	Graphics2D gg = (Graphics2D) g.create();
        //gg.setFont(new Font("Serif",Font.PLAIN,6));
        gg.translate(x,y);
    	if (edge_colors==null)
    	{
    		gg.setStroke(new BasicStroke(2));
    		gg.setColor(single_color);
    		drawTree(gg, true);
    	} else
    	{
    		double mid = 0.5*(edge_colors.length-1);
    		for (int c=edge_colors.length-1; 0<=c; --c)
    		{
    			Graphics2D g2 = (Graphics2D) gg.create();
    			double shift = mid+c-3;
    			g2.translate(2*shift, 2*shift);
    			g2.setColor(edge_colors[c]);
    			drawTree(g2, false);
    		}
    	}
    }	

    public static TreeIcon ratesIcon()
    {
    	return ratesIcon(SIZE);
    }
    public static TreeIcon ratesIcon(int size)
    {
    	Color[] edge_colors = new Color[3];
    	edge_colors[0] = RatesTreePanel.RATES_LOSS_COLOR;
    	edge_colors[1] = RatesTreePanel.RATES_DUPLICATION_COLOR;
    	edge_colors[2] = RatesTreePanel.RATES_GAIN_COLOR;
    	TreeIcon ratesIcon = new TreeIcon(size, edge_colors, BACKGROUND_COLOR);
    	return ratesIcon;
    }
    
//    public static TreeIcon zebraIcon(int size)
//    {
//    	Color[] edge_colors = new Color[3];
//    	edge_colors[0] = Color.WHITE;
//    	edge_colors[1] = TreePanel.TREE_EDGE_COLOR;
//    	edge_colors[2] = TreePanel.TREE_EDGE_COLOR;
//    	TreeIcon zebraIcon = new TreeIcon(size, edge_colors, BACKGROUND_COLOR);
//    	return zebraIcon;
//    }
    
    
}
