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

import psograph.graph.Edge;
import psograph.graph.Graph;
import psograph.graph.Node;




public class PercentageInLargestCluster implements Serializable, IGraphMeasument
{
	
	public boolean isCalculated()
	{
		return m_isCalculated;
	}
	
		/**
	 * 
	 */
	private static final long serialVersionUID = -8064060957650207979L;

		public PercentageInLargestCluster(Graph g)
		{
			m_graph = g;
			
			Collection<Node> c = m_graph.getHeaderNodesMap().values();
			masterNodeVector = new Vector<Node>(c);	
			
			NodesByCluster = new TreeMap<Integer,Vector<Node>>();
			
			m_isCalculated = false;
			
		}
		
		public Vector<Node> LCC() throws Exception
		{
			if(!m_isCalculated)
				throw new Exception("PercentageInLargestCluster not computed yet.");
			else
				return NodesByCluster.get(IdLCC);
			
		}
		
		public double Measure() throws Exception
		{
	
			DFS();
			return evaluate();		
		}
		
		double evaluate() throws Exception
		{
			double value =0;
			
			int num_nodes = m_graph.getHeaderNodesMap().size();
			
			int v = -1;
			int u = -1;
			
			Vector<Integer> cluster = new Vector<Integer>(NodesByCluster.keySet());
			for(int i =0 ; i < cluster.size();i++)
			{
				if( v < NodesByCluster.get(cluster.elementAt(i)).size())
				{
					v = NodesByCluster.get(cluster.elementAt(i)).size();
					u = cluster.elementAt(i);
				}
				
			}
			
			if(v < 0)
			{
				throw new Exception ("PercentageInLargestCluster could not evaluate");
			}
			else
				value = (double)v/(double)num_nodes;
			
			IdLCC = u;
			
			
			m_isCalculated = true;
			
			return value;
			
		}
		
		int findUnvisitedNode()
		{
			for(int i =0; i < m_graph.getHeaderNodesMap().size(); i++)
			{
				if( m_visited.get(masterNodeVector.get(i).getID()) == null )
					return masterNodeVector.get(i).getID();
				
			}
			
			return -1;
		}
		
		void DFS()
		{
			//int cluster_count = 1;
			Vector<Node> NodesInCluster = new Vector<Node>();
			
			m_visited = new TreeMap<Integer,Boolean>();
			Stack<Integer> m_stack = new Stack<Integer>();
			
			
			while(m_visited.size()< m_graph.getHeaderNodesMap().size())
			{
				//System.out.println("Cluster :"+m_cluster_id);
				
				int prev_node = findUnvisitedNode();
				NodesInCluster.add(m_graph.getHeaderNodesMap().get(prev_node));
				if(prev_node < 0)
				{
					//error/execption
					return;
				}
				
				m_visited.put(prev_node,true);
				
				//visit node
				Node node = m_graph.getHeaderNodesMap().get(prev_node);  
				
				TreeMap<Integer,Edge> tci = node.getNeighbors();
				if(tci == null)
				{					
					NodesByCluster.put(m_cluster_id,NodesInCluster);
					NodesInCluster = new Vector<Node>();
					m_cluster_id++;
				//	System.out.println("Isolated node "+node.getID());
					continue;
				}
				
				int next_node = Next(tci);
					
				
				while (next_node != -1 || (!m_stack.empty())) 
				{
					//System.out.println("prev_node " + prev_node + " | next_node "
					//		+ next_node);
					if (next_node != -1) 
					{
						m_stack.push(prev_node);
						m_visited.put(next_node, true);
						prev_node = next_node;
						NodesInCluster.add(m_graph.getHeaderNodesMap().get(next_node));
					} 
					else 
					{
						prev_node = m_stack.pop();
					}
					node = m_graph.getHeaderNodesMap().get(prev_node);
	
					tci = node.getNeighbors();
					next_node = Next(tci);
	
				}
				NodesByCluster.put(m_cluster_id,NodesInCluster);
				NodesInCluster = new Vector<Node>();
				m_cluster_id++;

			}
			
		}
		
		int Next(TreeMap<Integer,Edge> neighbors)
		{
			int val = -1;
			
			int i=0;
			int max = neighbors.size();
			Vector<Integer> vkeys = new Vector<Integer>(neighbors.keySet());
			
			
			while(val == -1 && i < max) 
			{
				if(m_visited.get(vkeys.get(i)) == null)
					val = vkeys.get(i);
				else
					i++;
			}
	
			return val;
		}
		
		
		/*
		To determine the percentage of nodes in the largest component we can use simple 
		modified version of DF Traversal.   We add in storage for buckets to hold the size 
	    of component, and an index counter to indicate which component we are traversing.  
	    As we visit nodes in the current component we increment the count.  When we run out of 
	    connected nodes in the current component to traverse (i.e. we have popped the stack and 
	    but still have nodes to visit), we increment the index of the components buckets.  The 
	    complexity of this algorithm is O(m+n), where m is number of edges, and n is number of 
	    nodes.
		*/
		Graph m_graph;

		TreeMap<Integer,Boolean> m_visited;
		int m_cluster_id = 0;
		int IdLCC = -1;
	
		
		
		Vector<Node> masterNodeVector;
		
		TreeMap<Integer,Vector<Node>> NodesByCluster;
		
		boolean m_isCalculated;

		public TreeMap<Integer, Vector<Node>> getNodesByCluster() 
		{
			return NodesByCluster;
		}

	
}
