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

/**
 * This represents a Edge.
 * 
 * There is some commented out code I believe in this file to support DAG and the concept 
 * of multiple edges between two nodes.
 * @author Patrick
 *
 */
public class Edge implements Serializable
{

	static final long serialVersionUID = 45L;

	/**
	 * Copy Constructor
	 * @param ci
	 */
	public Edge(Edge ci)
	{
		m_weight = ci.m_weight;
	}
	
	/**
	 * Constructor
	 * @param weight
	 */
	public Edge(double weight)
	{
		m_weight = weight;

	}
	
	/**
	 * Comparison of two objects
	 */
	public boolean equals (Object obj)
	{
		boolean ret = true;
		
		Edge e = (Edge)obj;
		
		if(Double.compare(m_weight, e.getWeight()) != 0)
		{
			ret = false;
		}
		return ret;
		
	}
	
	/**
	 * Mutator for weight value.
	 * @param weight
	 */
	public void modifyWeight(double weight)
	{
		m_weight = weight;
	}
	/**
	 * Accessor for weight.
	 * @return
	 */
	public double getWeight()
	{
		return m_weight;
	
	}
	
	private double m_weight; 	
	
	
	
	/* Only allow on weight per node to node connection
	
	ConnectionInfo(ConnectionInfo ci)
	{
		m_weight = new Vector<Integer>(ci.m_weight);
	}
	
	ConnectionInfo(int weight)
	{
		m_weight = new Vector<Integer>();
		m_weight.add(weight);

	}
	
	ConnectionInfo(int weight[])
	{
		m_weight = new Vector<Integer>();
		for(int i=0; i < weight.length; i++)
			m_weight.add(weight[i]);
	}

	void addWeight(int weight)
	{
		m_weight.add(weight);
	}

	void addWeights(int weight[])
	{
		m_weight = new Vector<Integer>();
		for(int i=0; i < weight.length; i++)
			m_weight.add(weight[i]);
	}
	
	void removeWeight(int weight)
	{
		m_weight.remove(new Integer(weight));
	}
	
	void removeWeights(int weight[])
	{
		for(int i=0; i < weight.length; i++)
			m_weight.remove(new Integer(weight[i]));
	}
	
	Vector<Integer> m_weight; 
	*/
}
