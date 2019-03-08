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
import java.util.Vector;

import psograph.graph.Graph;
import psograph.graph.Node;
import psograph.util.Util;


public class ClusteringCoefficient implements Serializable, IGraphMeasument {

	/**
	 * 
	 */
	private static final long serialVersionUID = 4600457450108515368L;

	public ClusteringCoefficient(Graph graph)
	{
		m_graph = graph;
		isCalculated = false;
	}
	
	public double Measure() throws Exception {
		
		
		Vector<Node> v_Nodes = new Vector<Node>(m_graph.getHeaderNodesMap().values());
		double ci[] = new double[v_Nodes.size()];
		
		for(int i = 0; i < v_Nodes.size(); i++)
		{
			Node n = v_Nodes.get(i);
			
			if(n.getDegree() == 0 || n.getDegree() == 1)
			{
				ci[i] = 0;
				continue;
			}
			
			
			Vector<Integer> neighbors = new Vector<Integer>( n.getNeighbors().keySet() );
			
			double Eij[] = new double[neighbors.size()];
			
			
			for(int j = 0 ; j < Eij.length; j++)
			{
				Eij[j] =  0;
				
				Node n1 = m_graph.getHeaderNodesMap().get(neighbors.get(j));
				for(int k =0; k < Eij.length; k++)
				{
					if(j==k)
						continue;
					
					Node n2 = m_graph.getHeaderNodesMap().get(neighbors.get(k));
					if(n1.isConnectedTo(n2) == true)
					{
						Eij[j] += 1;
					}
					
				}
			}
			
			ci[i] =  Util.sumArray(Eij) / (double)(n.getDegree() * (n.getDegree() - 1));
			
		}
		
		double res = Util.sumArray(ci)/ci.length;

		isCalculated = true;
		return res;
	}

	public boolean isCalculated() {

		return isCalculated;
	}
	
	Graph m_graph;
	boolean isCalculated;
}
