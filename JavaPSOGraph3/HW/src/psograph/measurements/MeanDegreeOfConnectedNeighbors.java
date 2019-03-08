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
import psograph.graph.Node;


public class MeanDegreeOfConnectedNeighbors 
{
	
	public MeanDegreeOfConnectedNeighbors(Graph g) 
{
	m_graph = g;
	isCalculated = false;
	m_resultById =null;
}

TreeMap<Integer, Double> m_resultById;


public TreeMap<Double, Integer> Measure() throws Exception 
{
	
	
	
	m_resultById = new TreeMap<Integer, Double>();
	TreeMap<Double, Integer> result = new TreeMap<Double, Integer>();
	
	Vector<Node> vNodes = new Vector<Node>(m_graph.getHeaderNodesMap().values());
	
	for(int i = 0; i< vNodes.size(); i++)
	{
		Node n = vNodes.get(i);
		double bucket = 0;
		
		if(n.getDegree() == 0)
		{
			m_resultById.put(n.getID(), 0.0);
			continue;
		}
		
		Vector<Integer> neighbors = new Vector<Integer>(n.getNeighbors().keySet());
		
		for(int j = 0; j < neighbors.size(); j++)
		{
			Node neighborNode = m_graph.getNode(neighbors.get(j));
			bucket += neighborNode.getDegree();
		}
		
		bucket = bucket/neighbors.size();
		
		m_resultById.put(n.getID(), bucket);
	}
	
	

	Vector<Integer> vNodeIds = new Vector<Integer>(m_resultById.keySet());
	for(int i=0; i < vNodeIds.size(); i++)
	{
		double meanDegreeOfConnectedNeigh = m_resultById.get(vNodeIds.get(i));
		
		int value = 0;
		if(result.containsKey(meanDegreeOfConnectedNeigh))
		{
			value = result.get(meanDegreeOfConnectedNeigh);
		}
		
		result.put(meanDegreeOfConnectedNeigh, value +1);
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
