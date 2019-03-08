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

import psograph.graph.CalculatedGraph;
import psograph.graph.Node;


public class Closeness {
	TreeMap<Integer, Double> m_resultById;

	public Closeness(CalculatedGraph g) 
	{
		m_graph = g;
		isCalculated = false;
		m_resultById =null;
	}




	public TreeMap<Double, Integer> Measure() throws Exception 
	{
		m_resultById = new TreeMap<Integer, Double>();
		TreeMap<Double, Integer> result = new TreeMap<Double, Integer>();
		
	
		//SPL SPL_measure = new SPL(this);
		//Graph graphs[] = SPL_measure.Measure(); 

		
		if(m_graph.getSPLGraphs() == null)
			return result;
		
		
		
		for(int iGraph = 0; iGraph < m_graph.getNumberOfNodes(); iGraph++)
		{
			double splSum =0;
			
		   // System.out.println("");
			for(int i =0; i < m_graph.getSPLGraphs()[iGraph].getNumberOfNodes() ; i++)
			{
				Node n = m_graph.getSPLGraphs()[iGraph].getNode(i);
				//Vector<Path> vPaths = n.getPaths();
				
				splSum += n.getSPLength();

//				System.out.println("Number of paths from "+ iGraph +" to "+n.getID()+" is " + vPaths.size());
//
//				for(int j = 0; j < vPaths.size(); j++)
//				{
//					Path p = vPaths.get(j);
//
//				//	System.out.print("start of path "+ p.getStart().getID()+" ");
//					for(int k =0; k < p.getLength(); k++)
//					{
//						Node n2 = p.getPath().get(k);
//						int id = n2.getID();
//						System.out.print(" " + id);
//					}
//				//	System.out.println(": Path Length " + p.getLength());					
//				}

			}
			
			
			
			double value = ((double)(m_graph.getNumberOfNodes() -1))/splSum;
			
			
			m_resultById.put(iGraph , value);
		
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
	
	
	CalculatedGraph m_graph;
	boolean isCalculated;
	
	public TreeMap<Integer, Double> ById() {
		return m_resultById;
	}
	
}
