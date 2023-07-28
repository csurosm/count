package count.gui.kit;
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


import java.awt.geom.Point2D;

import java.util.Comparator;

import java.util.Arrays;

import java.util.Stack;

/**
 * Class for calculate Convex hull using algorithm Graham Scan
 * 
 * @author	Amir shayeghi
 */
public class ConvexHull 
{
	
	private Stack<Point2D> point = new Stack<Point2D>();
	
	private Point2D[] points;
	

	public ConvexHull(Point2D[] pts) {
		this.points = pts;
		Arrays.sort(points, new PointComparator());
		point.push(points[0]);	
		Arrays.sort(points, new PointComparator(points[0]));	
		if (points.length>1)
		{
			point.push(points[1]);
			if (points.length>2)
			{
				point.push(points[2]);	
				for(int i = 3; i < points.length; ++i) {
					Point2D top = point.peek();
					Point2D nextToTop = point.get(point.size() - 2);		
					while(polarAngle(points[i], nextToTop, top) >= 0 && point.size() >= 3) {
						point.pop();
						top = point.peek();
						nextToTop = point.get(point.size() - 2);
					}
					
					point.push(points[i]);
				}
			}
		}
	}
	
    /**
     * Points of the convex hull in counterclockwise order. 
     * 
     * @return 
     */
	public Point2D[] getHull() {
		Point2D[] pointArray = new Point2D[]{};
		return point.toArray(pointArray);
	}


	private double polarAngle(Point2D p0, Point2D p1, Point2D p2) {
		//(p1 - p0) * (p2 - p0) = (x1 - x0)(y2 - y0) - (x2 - x0)(y1 - y0)
		//If this cross product is positive, then p0p1 is clockwise from p0p2; if negative, it is counter-clockwise.
		return (p1.getX() - p0.getX()) * (p2.getY() - p0.getY()) - (p2.getX() - p0.getX()) * (p1.getY() - p0.getY());
	}
}


 class PointComparator implements Comparator<Point2D> {
	private Point2D origin;

	public PointComparator() {
	}

	public PointComparator(Point2D origin) {
		this.origin = origin;
	}


	@Override
	public int compare(Point2D p1, Point2D p2) {
		double angle;
		if(origin != null)
			angle = (p1.getX() - origin.getX()) * (p2.getY() - origin.getY()) - (p2.getX() - origin.getX()) * (p1.getY() - origin.getY());
		else
			angle = p1.getX()*p2.getY() - p1.getY()*p2.getX();
		return (int) angle;
	}
}
