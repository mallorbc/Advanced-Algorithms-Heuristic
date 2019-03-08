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

import java.util.TreeMap;
import java.util.Vector;

import psograph.graph.Graph;
import psograph.graph.GraphConstants;
import psograph.graph.Node;
import psograph.graph.NonLinearCostFunction;


public class ConnectedNodesOutsideRadius  {

	double numEdges = 0;
	
	public ConnectedNodesOutsideRadius(Graph g) 
	{
		m_graph = g;
		isCalculated = false;
	}
	
	public void setRadius(double r)
	{
		m_radius = r;
	}

	TreeMap<Integer, Integer> m_resultById;
	
	public TreeMap<Integer, Integer> Measure() throws Exception 
	{
		
		m_resultById = new TreeMap<Integer, Integer>();
		TreeMap<Integer, Integer> result = new TreeMap<Integer, Integer>();
		
		Vector<Node> nodes = new Vector<Node>(m_graph.getHeaderNodesMap().values() );
		
		int[][] adjMatrix = m_graph.getAdjecencyGraph2();
		
		for(int i = 0; i < m_graph.getNumberOfNodes(); i++)
		{
			int node_i_id = nodes.get(i).getID();
			Node Nodei = m_graph.getNode(node_i_id);
			
			int connectedNodesWithinRadius = 0;
			
			for(int j = 0; j < m_graph.getNumberOfNodes(); j++)
			{
				if(adjMatrix[i][j]== 1)
				{
					int node_j_id = nodes.get(j).getID();
					Node Nodej = m_graph.getNode(node_j_id);
					
					boolean isConnected = m_graph.isNodeAConnectedToNodeB(Nodei, Nodej);
					
					if(!isConnected)
						throw new Exception("Connnected Nodes withing R - node mismatch");
					
					double cost = Nodei.getConnectionInfo(Nodej).getWeight();
					
					double distance = NonLinearCostFunction.getDistance(cost);
					
					if(Double.compare(distance, m_radius) > 0.0)
					{
						connectedNodesWithinRadius++;
						numEdges=numEdges+1.0;
					}
				}
			}
			
			m_resultById.put(node_i_id, connectedNodesWithinRadius);
			
		}

		Vector<Integer> vNodeIds = new Vector<Integer>(m_resultById.keySet());
		for(int i=0; i < vNodeIds.size(); i++)
		{
			int connectNodeWithinR = m_resultById.get(vNodeIds.get(i));
			
			int value = 0;
			if(result.containsKey(connectNodeWithinR))
			{
				value = result.get(connectNodeWithinR);
			}
			
			result.put(connectNodeWithinR, value +1);
		}
		
		isCalculated = true;
		return result;
	}

	public double getNumEdges() {

		double val = numEdges/2.0;
		
		return val;
	}

	public boolean isCalculated() {

		return isCalculated;
	}

	double m_radius = GraphConstants.MeasurementR;
	
	Graph m_graph;
	boolean isCalculated;
	public TreeMap<Integer, Integer> ById() 
	{
		return m_resultById;
	}

}
