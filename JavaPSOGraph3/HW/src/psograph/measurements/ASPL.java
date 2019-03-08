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

import java.io.Serializable;
import java.util.LinkedList;
import java.util.TreeMap;
import java.util.Vector;

import psograph.graph.Edge;
import psograph.graph.Graph;
import psograph.graph.Node;




public class ASPL implements Serializable, IGraphMeasument
{
	/**
	 * 
	 */
	private static final long serialVersionUID = -3199607653710066896L;
	public ASPL(Graph g)
	{
		m_graph = g;
		m_isCalculated = false;
		m_solns = new Graph[m_graph.getNumberOfNodes()];
	}
	

	double determineSPL(boolean inverse) throws Exception
	{
		double value = 0;
		Vector<Integer> vhn = new Vector<Integer>(m_graph.getHeaderNodesMap().keySet());
		double number_of_edges = m_graph.getNumberOfNodes() * (m_graph.getNumberOfNodes() - 1);
		
		for(int i =0; i < vhn.size(); i++)
		{
			initSolutionGraph();
			Graph totallyTemp = new Graph(m_graph);
			evaulateNodeWithBFS(vhn.get(i), totallyTemp);
			
			double val;
			if(inverse)
			{
				val = calcValues.sumInverseWeights()/2.0;
				//System.out.println("Sum all weights " + val);
			}
			else
			{
				val = calcValues.SumAllWeights()/2.0;
			}
			
			
			val =  val/number_of_edges;
			if(inverse)
			{
				//System.out.println("perform calc " + val);
			}
			
			value += val;
			
			//TODO: turn this back on 
			m_solns[i] = new Graph(calcValues);
			
			totallyTemp = null; 
		}
		
		m_isCalculated = true;
		return value;
		
	}
	
	void initSolutionGraph() throws Exception
	{
		calcValues = new Graph(0);
		
		//First add in nodes
		Vector<Integer> vhn = new Vector<Integer>(m_graph.getHeaderNodesMap().keySet());
		
		for(int i =0; i < vhn.size(); i++)
		{
			calcValues.addNode(vhn.get(i), m_graph.getHeaderNodesMap().get(vhn.get(i)).getX(), 
					m_graph.getHeaderNodesMap().get(vhn.get(i)).getY() );
		}
	}
	
	public double Measure() throws Exception
	{
		double value = -1;
		if(m_graph.isFullyConnected() == false)
		{
			//System.out.println("ASPL - This measure is not accurate for a network with any disconnected nodes");
			return value;
		}
		
		value = determineSPL(false);
		//calcValues.printWithWeights();
		//value = peformCalculation();
		return value;
	}
	
	void evaulateNodeWithBFS(int NodeId , Graph g) throws Exception
	{
		
		//initialize the visit map
		m_visited = new TreeMap<Integer,Boolean>();
		LinkedList<Integer> m_linkedList = new LinkedList<Integer>();	
		
		Node n =g.getHeaderNodesMap().get(NodeId);
			
		m_visited.put(n.getID(), true);
		n.setDepth(0);
		m_linkedList.add(n.getID());
			

		while (!m_linkedList.isEmpty()) 
		{
			
			
			int u = m_linkedList.remove();
			//visit node
			Node node = g.getHeaderNodesMap().get(u);  
			
				
			TreeMap<Integer,Edge> tci = node.getNeighbors();
			if(tci == null)
			{		
				continue;
			}
			int i =0;
			
			Vector<Integer> vci = new Vector<Integer>(tci.keySet());
			
			for(i=0; i < vci.size(); i++)
			{				
				if (!IfValidUnvisitedNode(vci.get(i)))
					continue;
				
				g.getHeaderNodesMap().get(vci.get(i)).setDepth(node.getDepth() +1);
				m_linkedList.add(vci.get(i));
				m_visited.put(vci.get(i), true);
				
				//System.out.println("u " + u + " | next_node "
				//		+ vci.get(i));
			}
		}
		
		//Now that we have a graph with depth, we can make connections for this node.
		Vector<Integer> vhn = new Vector<Integer>(g.getHeaderNodesMap().keySet());	
		for(int i =0; i < vhn.size(); i++)
		{
			if(g.getHeaderNodesMap().get(vhn.get(i)).getID() == NodeId)
				continue;
			
			if(calcValues.getHeaderNodesMap().get(NodeId).isConnectedTo(g.getHeaderNodesMap().get(vhn.get(i))) == false)
			{
				calcValues.addConnection(NodeId, g.getHeaderNodesMap().get(vhn.get(i)).getID(), g.getHeaderNodesMap().get(vhn.get(i)).getDepth());
			}
		}
		
	}
	
	boolean IfValidUnvisitedNode(int id) throws Exception
	{
		if( m_visited.get(id) == null )
			return true;
		else
			return false;
	}
	
	public boolean isCalculated()
	{
		return m_isCalculated;
	}
	
	public TreeMap<Integer, Integer> getSPLDistribution() throws Exception
	{
		if(m_graph.isFullyConnected() == false)
		{
			//System.out.println("ASPL - This measure is not accurate for a network with any disconnected nodes");
			return null;
		}

		TreeMap<Integer, Integer> result = new TreeMap<Integer, Integer>();
		for(int i = 0; i < m_graph.getNumberOfNodes(); i++)
		{
			TreeMap<Double, Integer> iteration = m_solns[i].getWeightDistribution();
			
			Vector<Double> vWeight = new Vector<Double>(iteration.keySet());
			for(int j =0; j < vWeight.size(); j++)
			{
				double val =vWeight.get(j);
				int count = iteration.get(val);
				
				
				if(result.containsKey((int)Math.rint(val)))
				{
					count += result.get((int)Math.rint(val));
				}
				
				int intVal = (int)Math.rint(val);
				
				result.put(intVal, count);
			}
		
		}
		return result;
	}
	

	public Graph calcValues; 
	TreeMap<Integer,Boolean> m_visited;
	Graph m_graph; // the original graph
	
	Graph m_solns[];
	
	boolean m_isCalculated;
	
}
