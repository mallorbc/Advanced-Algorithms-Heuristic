// This is a library to be used to represent a Graph and various measurments for a Graph
//  and to perform optimization using Particle Swarm Optimization (PSO)
//    Copyright (C) 2008, 2015 
//       Patrick Olekas - polekas55@gmail.com
//       Ali Minai - minaiaa@gmail.com
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
package psograph.measurements;

import java.io.Serializable;


import java.util.*;

import psograph.graph.Graph;
import psograph.graph.Node;

public class AvgNodeDegree implements Serializable,IGraphMeasument {

	/**
	 * 
	 */
	private static final long serialVersionUID = -6387430381422847299L;

	public AvgNodeDegree(Graph graph)
	{
		m_graph = graph;
		isCalculated = false;
	}
	
	@Override
	public double Measure() throws Exception {
		double res = 0;

		Vector<Node> nodes = new Vector<Node>(m_graph.getHeaderNodesMap().values());
		
		double numNodes = nodes.size();
		
		for(int i =0; i < numNodes; i++)
		{
			double t = nodes.get(i).getDegree() / numNodes;
			res = res + t;
		}
		
		isCalculated = true;
		return res;
	}

	@Override
	public boolean isCalculated() {

		return isCalculated;
	}
	Graph m_graph;
	boolean isCalculated;
}
