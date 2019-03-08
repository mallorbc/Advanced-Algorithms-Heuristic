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
import java.util.TreeMap;
import java.util.Vector;
import java.util.LinkedList;

import psograph.graph.Edge;
import psograph.graph.Graph;
import psograph.graph.Node;




public class Diameter implements Serializable, IGraphMeasument
{
	
	/* For each node in the graph we do BF Traversal.  Recording of the length of the longest path in the traversal 
	 * for each node.  After we do this for all nodes, the node with the longest path is equal to the diameter of the graph.  
	 * We then compute the change in diameter from the original network.  */
	
	
	public boolean isCalculated()
	{
		return m_isCalculated;
	}
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 6814045651581260925L;
	public Diameter(Graph g, Vector<Node> lcc)
	{
		m_graph = g;

		diameterVector = new Vector<Integer>();
		LCC=lcc;
		
		m_isCalculated = false;
	}
	
	
	boolean IfValidUnvisitedNode(int id)
	{

		if( m_visited.get(id) == null && 
				LCC.contains(m_graph.getHeaderNodesMap().get(id)))
		{
			return true;
		}

		return false;
	}
	
	
	int evaluate1() throws Exception
	{
		int value =0;
		
	
		
		int v = -1;
		for(int i =0; i < LCC.size(); i++)
		{
			if(v < LCC.get(i).getDepth())
				v = LCC.get(i).getDepth();			
		}
		
		if(v < 0)
			throw new Exception("Diameter::evaulate1 - cannot determine diamter");
		else
			value = v;
		
		return value;
		
	}
	
	double evaluate() throws Exception
	{
		double value =0;

		int v = -1;
		for(int i =0; i < diameterVector.size(); i++)
		{
			if(v < diameterVector.get(i))
				v = diameterVector.get(i);			
		}
		
		if(v < 0)
		{
			//System.out.println("Diamter::evaulate - cannot determine diamter");
			value = v;
		}
		else
			value = v;
		
		m_isCalculated=true;
		
		return value;
		
	}
	

	int Next(TreeMap<Integer,Edge> neighbors)
	{
		int val = -1;
		
		int i=0;
		int max = neighbors.size();
		Vector<Integer> vkeys = new Vector<Integer>(neighbors.keySet());
		
		
		while(val == -1 && i < max) 
		{

			if(LCC.contains(m_graph.getHeaderNodesMap().get(vkeys.get(i))) &&
							m_visited.get(vkeys.get(i)) == null)
				val = vkeys.get(i);
			else
				i++;
		}

		return val;
	}
	
	int evaulateNodeWithBFS(Node n) throws Exception
	{
		
		//initialize the visit map
		m_visited = new TreeMap<Integer,Boolean>();
		LinkedList<Integer> m_linkedList = new LinkedList<Integer>();	
			
		m_visited.put(n.getID(),true);
		n.setDepth(0);
		m_linkedList.add(n.getID());
			

		while (!m_linkedList.isEmpty()) 
		{
			
			
			int u = m_linkedList.remove();
			//visit node
			Node node = m_graph.getHeaderNodesMap().get(u);  
			
				
			TreeMap<Integer,Edge> tci = node.getNeighbors();
			if(tci == null)
			{		
				//System.out.println("Isolated node "+node.getID());
				//System.out.println("THIS IS PROBALLY A SEVERE ERROR");
				//System.out.println("Diameter measurement - this measure is not accurate for a network with any disconnected nodes");
				return -2;
			}
			int i =0;
			
			Vector<Integer> vci = new Vector<Integer>(tci.keySet());
			
			for(i=0; i < vci.size(); i++)
			{
				if (!IfValidUnvisitedNode(vci.get(i)))
					continue;
				
				m_graph.getHeaderNodesMap().get(vci.get(i)).setDepth(node.getDepth() +1);
				m_linkedList.add(vci.get(i));
				m_visited.put(vci.get(i), true);
				
			//	System.out.println("u " + u + " | next_node "
			//			+ vci.get(i));
			}
		}
		
				
		return evaluate1();
	
	}
	
	
	public double Measure() throws Exception
	{
		if (m_graph.getHeaderNodesMap().size() == 0)
		{
			throw new Exception("Diamter::Measure - No nodes in network.");
		}
		
		int i=0;
		for(i=0; i < LCC.size();i++)
		{
			
			int d = evaulateNodeWithBFS(LCC.elementAt(i));
			//System.out.println("Diameter BFS for node :"+LCC.elementAt(i).m_id+ " is: "+d);
			diameterVector.add(d);
			
		}
		return evaluate();
	}
	
	Graph m_graph;
	Vector<Integer> diameterVector;
	TreeMap<Integer,Boolean> m_visited;
	Vector<Node> LCC;
	
	boolean m_isCalculated;

	
}
