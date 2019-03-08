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
package psograph.test;


import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.TreeMap;
import java.util.Vector;

import org.junit.Before;
import org.junit.Test;

import psograph.graph.CalculatedGraph;
import psograph.graph.Graph;
import psograph.graph.NonLinearCostFunction;


public class MeanTest {

	Graph m_graph;
	CalculatedGraph m_calcGraph;
	
	@Before
	public void setUp() throws Exception 
	{
			
	}
		
	@Test
	public void testMeanLength() throws Exception
	{
		m_graph = TestUtils.createSquareGraph();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Integer, Double> nodeMeanLenghtDistribution = m_calcGraph.getMeanEdgeLengthById();
		
		Vector<Integer> vNodes = new Vector<Integer>(nodeMeanLenghtDistribution.keySet());		
	
		if(vNodes.size() <= 0)
		{
			fail("The node degree distribution is empty");
		}
	
		for(int i=0; i < vNodes.size(); i++)
		{
			int nodeId = vNodes.get(i);
			double meanLength = nodeMeanLenghtDistribution.get(nodeId);
			
			if(nodeId == 0 )
			{
				double v = (.1 + .1 + Math.sqrt(Math.pow(.1, 2) +  Math.pow(.1, 2) ))/3.0;
				
				
				assertTrue("node id "+nodeId+
						"has mean lenght of "+ meanLength +" which doesn't match " +v,
						Double.compare(meanLength, v)==0);	
			}
			else if(nodeId == 1 )
			{
				assertTrue("node id "+nodeId+
						"has mean lenght of "+ meanLength,
						Double.compare(meanLength, .1)==0);	
			}
			else if(nodeId == 3 )
			{
				assertTrue("node id "+nodeId+
						"has mean lenght of "+ meanLength,
						Double.compare(meanLength,
								Math.sqrt(Math.pow(.1, 2) +  Math.pow(.1, 2) ))==0);	
			}
			else if(nodeId == 2 )
			{
				assertTrue("node id "+nodeId+
						"has mean lenght of "+ meanLength,
						Double.compare(meanLength, .1)==0);	
			}
			else
			{
				fail("A node with a degree not anticipated for");
			}
		//	System.out.println("Mean Avg Degree "+ avgDegree + " Nodes with neighbor of this avg degree "+
		//			numAvgDegree);
		}
	}
	
	@Test
	public void testMeanCost() throws Exception
	{
		m_graph = TestUtils.createSquareGraph();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Integer, Double> nodeMeanLenghtDistribution = m_calcGraph.getMeanEdgeCostById();
		
		Vector<Integer> vNodes = new Vector<Integer>(nodeMeanLenghtDistribution.keySet());		
	
		if(vNodes.size() <= 0)
		{
			fail("The node degree distribution is empty");
		}
		
		NonLinearCostFunction cf = new NonLinearCostFunction();
	
		for(int i=0; i < vNodes.size(); i++)
		{
			int nodeId = vNodes.get(i);
			double meanCost = nodeMeanLenghtDistribution.get(nodeId);
			
			if(nodeId == 0 )
			{
				double v = cf.calculate(.1) + cf.calculate(.1)
						+ cf.calculate(Math.sqrt(Math.pow(.1, 2) +  Math.pow(.1, 2) ));
				double v_cost = v/3.0;
				
			//	java.lang.AssertionError: node id 0has mean lenght of 0.020000000000000004 instead of 0.026666666666666672



				//Math round off error
				assertTrue("node id "+nodeId+
						"has mean length of "+ meanCost + " instead of "+v_cost,
						Double.compare(meanCost, v_cost)==0);	
			}
			else if(nodeId == 1 )
			{
				assertTrue("node id "+nodeId+
						"has mean lenght of "+ meanCost+ " instead of "+cf.calculate(.1),
						Double.compare(meanCost, cf.calculate(.1))==0);	
			}
			else if(nodeId == 3 )
			{
				assertTrue("node id "+nodeId+
						"has mean lenght of "+ meanCost+ " instead of "+cf.calculate(Math.sqrt(Math.pow(.1, 2) +  Math.pow(.1, 2) )),
						Double.compare(meanCost,
								cf.calculate(Math.sqrt(Math.pow(.1, 2) +  Math.pow(.1, 2) )))==0);	
			}
			else if(nodeId == 2 )
			{
				assertTrue("node id "+nodeId+
						"has mean lenght of "+ meanCost+ " instead of "+cf.calculate(.1),
						Double.compare(meanCost, cf.calculate(.1))==0);	
			}
			else
			{
				fail("A node with a degree not anticipated for");
			}
//			System.out.println("Node Id  "+ nodeId + " has mean cost of "+
//					meanCost);
		}
	}
	

}
