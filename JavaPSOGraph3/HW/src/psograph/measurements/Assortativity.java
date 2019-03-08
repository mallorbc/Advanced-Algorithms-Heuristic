// This is a library to be used to represent a Graph and various measurments for a Graph
//  and to perform optimization using Particle Swarm Optimization (PSO)
//    Copyright (C) 2008, 2015 Patrick Olekas
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

import java.util.Vector;

import psograph.graph.Graph;
import psograph.graph.Node;


public class Assortativity implements IGraphMeasument {

	public Assortativity(Graph g) 
	{
		m_graph = g;
		isCalculated = false;
	}

	@Override
	public double Measure() throws Exception 
	{
		double value = 0.0;
		
		int l = m_graph.getTotalEdges();
		// number of links
		
		// r(g) = (first term - second term) / 
		//        (third term - second term)
		
		// first term = SUM(i,j)(  di * dj)/l
		
		// second term =  SUM(i,j) (.5 ( di + dj)/l)^2
		
		//third term = SUM(i,j) (.5 ( di^2 + dj^2)/l)
		
		double first_term = 0;

		
		double second_term = 0;

		
		double third_term = 0;

		Vector<Node> nodes = new Vector<Node>(m_graph.getHeaderNodesMap().values() );
		
		int[][] adjMatrix = m_graph.getAdjecencyGraph2();
		
		for(int i = 0; i < m_graph.getNumberOfNodes(); i++)
		{
			for(int j = i+1; j < m_graph.getNumberOfNodes(); j++)
			{
				if(adjMatrix[i][j]== 1)
				{
					int node_i_id = nodes.get(i).getID();
					int node_j_id = nodes.get(j).getID();
					
					Node Nodei = m_graph.getNode(node_i_id);
					Node Nodej = m_graph.getNode(node_j_id);
					
					// first term = SUM(i,j)(  di * dj)/l
					
					// second term =  SUM(i,j) (.5 ( di + dj)/l)^2
					
					//third term = SUM(i,j) (.5 ( di^2 + dj^2)/l)
					
					first_term += ((double)Nodei.getDegree() * (double)Nodej.getDegree())/l;
					
					
					second_term += Math.pow((.5 * ((double)Nodei.getDegree() + (double)Nodej.getDegree())/l), 2);				
					
					third_term += .5* ( Math.pow((double)Nodei.getDegree(), 2) + Math.pow((double)Nodej.getDegree(), 2))/l ;
					
				}
			}
		}
		
		value = (first_term - second_term) / (third_term - second_term);
		
		
		
		isCalculated = true;
		return value;
	}

	@Override
	public boolean isCalculated() {

		return isCalculated;
	}

	Graph m_graph;
	boolean isCalculated;
	
}
