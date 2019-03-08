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
package psograph.graph;

import java.io.Serializable;
import java.util.Vector;
/**
 * A cost Basis is used to benchmark to future costs.  A Cost basis is exponentially connected 
 * graph.
 * @author Patrick
 *
 */
public class CostBasis extends Graph implements Serializable{

	static final long serialVersionUID = 2765L;
	
	/**
	 * Constructor for graph.
	 * @param g
	 */
	public CostBasis(Graph g) 
	{	//TODO this probably should only copy the node and node locations and 
		// not copy the edges
		super(g);
	}
	
	/**
	 * Goes through each node in the graph and adds in the longest link possible.
	 * @throws Exception
	 */
	void addInLongestLinks() throws Exception
	{
		int i=0;
		
		int numNodes = getHeaderNodesMap().values().size();
		Vector<Node> v_Nodes2 = new Vector<Node>(getHeaderNodesMap().values());
		for(; i < numNodes ; )
		{

			int t_idFurthestNode = r.nextInt(v_Nodes2.size());
			Node furthestNeighborA = m_nodeLoc.getHeaderNodesMap().get(t_idFurthestNode);
				
			if(furthestNeighborA == null)
			{
				System.out.println("null node returned in CostBasis::Generate");
				throw new Exception("null node returned in CostBasis::Generate");
			}
			else
			{
				//Remove the node we just picked so we don't get pick duplicates
				v_Nodes2.remove(t_idFurthestNode);
				
				i++;
			    
			    Integer furthestNeighborId = m_nodeLoc.pickFurthestNeighbor(furthestNeighborA);
			    Node furthestNeighborB = m_nodeLoc.getHeaderNodesMap().get(furthestNeighborId);
			     
			    Edge ci = furthestNeighborA.getConnectionInfo(furthestNeighborB);
			    addConnection(furthestNeighborB.getID(), furthestNeighborA.getID(), ci.getWeight());
			    m_nodeLoc.removeConnection(furthestNeighborB.getID(), furthestNeighborA.getID());
			}
		}
		
	
	}
	
	/**
	 * Generate a Cost base with the specified number of edges.  YOu need to at least
	 * have as many edges as nodes, and not more edges than is possible in a complete graph
	 * @param numUniqueEdges
	 * @throws Exception
	 */
	public void generate(int numUniqueEdges) throws Exception
	{

		m_nodeLoc = new NodeLocationCalculator(this, false);
		m_nodeLoc.calculateResults();
		
		m_nodeLoc.setClosestPercentage(0.0); // We want all random
		int i=0;
		
		Vector<Node> v_Nodes = new Vector<Node>(m_nodeLoc.getHeaderNodesMap().values());
		
		if(numUniqueEdges < v_Nodes.size())
			throw new Exception("Minimum number of edges equals number of nodes.");
		else if(numUniqueEdges > v_Nodes.size() * (v_Nodes.size() -1) )
			throw new Exception("Maximum number of edges equals N*N-1.");

		//addInLongestLinks();
		i = v_Nodes.size();
		
		for(; i < numUniqueEdges ; i++)
		{
			Vector<Node> v_Nodes2 = new Vector<Node>(getHeaderNodesMap().values());
			int t_idA = r.nextInt(v_Nodes2.size());
			Node NodeA = v_Nodes2.get(t_idA);
				
			//Remove the node we just picked so we don't get pick duplicates
			v_Nodes2.remove(t_idA);
				

			Node n = m_nodeLoc.chooseNode(NodeA);
			if(n == null)
			{
				System.out.println("null node returned for choose random/close node");
			}
			else
			{
				Edge ci = n.getConnectionInfo(NodeA);
				this.addConnection(NodeA.getID(), n.getID(), ci.getWeight());
				m_nodeLoc.removeConnection(NodeA.getID(), n.getID());
			}
		}

	}

	NodeLocationCalculator m_nodeLoc;

}
