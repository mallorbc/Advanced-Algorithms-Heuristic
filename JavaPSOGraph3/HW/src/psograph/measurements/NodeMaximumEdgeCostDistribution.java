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

import psograph.graph.*;




public class NodeMaximumEdgeCostDistribution 
{

	public NodeMaximumEdgeCostDistribution(Graph g) 
	{
		m_graph = g;
		isCalculated = false;
	}
	
TreeMap<Integer, Double> m_resultById;
	

	public TreeMap<Double, Integer> Measure() throws Exception 
	{
		m_resultById = new TreeMap<Integer, Double>();
		TreeMap<Double, Integer> result = new TreeMap<Double, Integer>();
		
		Vector<Node> nodes = new Vector<Node>(m_graph.getHeaderNodesMap().values() );
		
		int[][] adjMatrix = m_graph.getAdjecencyGraph2();
		
		for(int i = 0; i < m_graph.getNumberOfNodes(); i++)
		{
			int node_i_id = nodes.get(i).getID();
			Node Nodei = m_graph.getNode(node_i_id);
			
			double maxEdgeCostForNode = -1;
			
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
					
					//double distance = NonLinearCostFunction.getDistance(cost);
					
					if(Double.compare(cost, maxEdgeCostForNode) > 0.0)
					{
						maxEdgeCostForNode = cost;
					}
				}
			}

			
			m_resultById.put(node_i_id, maxEdgeCostForNode);

		}

		
		Vector<Integer> vNodeIds = new Vector<Integer>(m_resultById.keySet());
		for(int i=0; i < vNodeIds.size(); i++)
		{
			double maxEdgeCostForNode = m_resultById.get(vNodeIds.get(i));
			
			int value = 0;
			if(result.containsKey(maxEdgeCostForNode))
			{
				value = result.get(maxEdgeCostForNode);
			}
			
			result.put(maxEdgeCostForNode, value +1);
		}
		
		
		isCalculated = true;
		return result;
	}

	public boolean isCalculated() {

		return isCalculated;
	}


	Graph m_graph;
	boolean isCalculated;
	public TreeMap<Integer, Double> ById() {
		return m_resultById;
	}

	
	
	
}
