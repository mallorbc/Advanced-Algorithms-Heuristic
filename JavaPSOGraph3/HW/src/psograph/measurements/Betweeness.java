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


public class Betweeness {
	
	TreeMap<Integer, Double> m_resultById;

	public Betweeness(CalculatedGraph g) 
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
		
		
		
		//i node being checked for betweeness
		//j node to start from
		//k node to end to
		
		for(int i = 0; i < m_graph.getNumberOfNodes(); i++)
		{
			//System.out.println("Betweeness - iter"+i);
			
			double resultNode=0;
			double value =0;

		//System.out.println("");
			for(int j =0; j <  m_graph.getNumberOfNodes() ; j++)
			{
				if(i==j)
					continue;
				
				for(int k=0; k < m_graph.getNumberOfNodes(); k++)
				{
					if(j==k || i==k)
						continue;
					
					Node n = m_graph.getSPLGraphs()[j].getNode(k);
					
			//		Vector<Path> vPaths = n.getPaths();
					
					double splNum =0;
					double splIncludeThisNodeInBetween = 0;
					
//					System.out.println("Number of paths from "+ j +" to "+k+" is " + vPaths.size());
//					System.out.println("Id in between "+i);
//					for(int jj = 0; jj < vPaths.size(); jj++)
//					{
//						Path p = vPaths.get(jj);
//
//					//	System.out.print("start of path "+ p.getStart().getID()+" ");
//						for(int kk =0; kk < p.getLength(); kk++)
//						{
//							Node n2 = p.getPath().get(kk);
//							int id = n2.getID();
//							System.out.print(" " + id);
//						}
//				//		System.out.println(": Path Length " + p.getLength());					
//					}
					
					

						splNum = n.numOfSPLength();
						splIncludeThisNodeInBetween = n.numberOfSPLNodeParticipatesIn(k,m_graph.getSPLGraphs()[i].getNode(i));
						
						value += splIncludeThisNodeInBetween/splNum;
					//	System.out.println("value update to "+ value);

				}				
			}
			
			double denom = ((double)(m_graph.getNumberOfNodes() -1)* (double)(m_graph.getNumberOfNodes() -2));
			
		//	System.out.println("denon is " +denom);
			resultNode = value/denom;



			m_resultById.put(i,resultNode);
			
			
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
