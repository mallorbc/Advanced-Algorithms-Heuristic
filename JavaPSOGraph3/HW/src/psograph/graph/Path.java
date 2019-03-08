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
 * Represents a walk on a graph.
 * @author Patrick
 *
 */
public class Path  implements Serializable 
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 5406576084452654454L;
	Vector<Node> m_path;
	Node m_startNode;
	Node m_endNode;
	
	
	/**
	 * Start a Path staring at the specified node.
	 * We do not add the start node to keep path to make it zero for the case
	 * there are no more nodes in the path.
	 * @param start
	 * @throws Exception
	 */
	public Path(Node start)throws Exception
	{
		m_startNode = start;
		m_path = new Vector<Node>();
	}
	
	/**
	 * Create a new Path appending the new node as the end.
	 * This method checks to make sure the new end node does not create a 
	 * cycle.  And exception will be thrown in that case.
	 * @param p
	 * @param newEnd
	 * @throws Exception
	 */
	public Path(Path p, Node newEnd) throws Exception
	{
		m_startNode = p.getStart();
		m_path = new Vector<Node>();
		for(int i =0; i < p.m_path.size(); i++)
		{
			if(p.m_path.get(i) == newEnd)
				throw new Exception("Path just created a cycle");
			
			m_path.add(p.m_path.get(i));
		}
		m_path.add(newEnd);
		m_endNode = newEnd;
	}
	
	/**
	 * Returns the Vector of Nodes representing the Path
	 * @return
	 */
	public Vector<Node> getPath()
	{
		return m_path;
	}
	
	/**
	 * Returns Start Node.
	 * @return
	 */
	public Node getStart()
	{
		return m_startNode;
	}
	
	/**
	 * Gets End Node of Path.
	 * @return
	 */
	public Node getEnd()
	{
		return m_endNode;
	}
	
	/**
	 * Get number of Nodes on Path.  If only one Node on Path this will return 1.  
	 * @return
	 */
	public int getNumberOfNodesOnPath()
	{
		return 1 + m_path.size() ;
	}
	
	/**
	 * Returns the length of the Path.
	 * Note that start node is not counted, so if Path only contains
	 * 1 node is of length zero.
	 * @return
	 */
	public int getLength()
	{
		return m_path.size() ;
	}

	/**
	 * Returns true if Node n is on the Path.
	 * @param n
	 * @return
	 */
	public boolean isNodeOnPath(Node n)
	{
		boolean result = false;
		if(n.getID() == m_startNode.getID() 
				|| n.getID() == m_endNode.getID())
		{
			return result;
		}
		else
		{   
			boolean stop = false;
			for(int i =0; i < m_path.size() && !stop; i++)
			{
				if(n.getID() == m_path.get(i).getID())
				{
					result = true;
					stop = true;
				}
			}
		}
		
		return result;
		
	}
	/**
	 * Debug routine to print out Path.
	 */
	public void print() 
	{
		System.out.print(m_startNode+" ");
		
		for(int i =0; i < m_path.size(); i++)
		{
			System.out.print(m_path.get(i)+" ");
		}
		System.out.println("");
		
	}
		

}
