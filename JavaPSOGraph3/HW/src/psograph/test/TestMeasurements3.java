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

import static org.junit.Assert.*;

import java.util.TreeMap;
import java.util.Vector;

import org.junit.Before;
import org.junit.Test;

import psograph.graph.*;


public class TestMeasurements3 {

	
	Graph m_graph;
	CalculatedGraph m_calcGraph;
	
	@Before
	public void setUp() throws Exception 
	{
		
	}
	
	

	@Test
	public void testGetBetweeness()throws Exception 
	{
		m_graph = TestUtils.createTestGraph7();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Double, Integer> betweenessDist = m_calcGraph.getBetweeness();
		
		Vector<Double> vBetweeness = new Vector<Double>(betweenessDist.keySet());

		if(vBetweeness.size() <= 0)
		{
			fail("The weight distribution is empty");
		}

		for(int i=0; i < vBetweeness.size(); i++)
		{
			double betweeness = vBetweeness.get(i);
			double betweenessNumber = betweenessDist.get(betweeness);
			
			if(Double.compare(0.0, betweeness) == 0)
			{
				assertTrue("cost of "+betweeness+" there are "+ betweenessNumber,
						betweenessNumber==1);	
			}
			else if(Double.compare(0.08333333333333333, betweeness) == 0)
			{
				assertTrue("cost of "+betweeness+" there are "+ betweenessNumber,
						betweenessNumber==1);	
			}
			else if(Double.compare(0.16666666666666666, betweeness) == 0)
			{
				assertTrue("cost of "+betweeness+" there are "+ betweenessNumber,
						betweenessNumber==2);	
			}
			else if(Double.compare(0.5833333333333334, betweeness) == 0)
			{
				assertTrue("cost of "+betweeness+" there are "+ betweenessNumber,
						betweenessNumber==1);	
			}
			else
			{
				fail("A betweeness not anticipated for "+betweeness);
			}
		//	System.out.println(" betweeness "+ betweeness + 
		//			" Number of with this betweeness "+	betweenessNumber);
		}
	}

	@Test
	public void testGetBetweenessInCompleteLayout()throws Exception 
	{
		m_graph = TestUtils.createTestGraph7();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Double, Integer> betweenessInCompleteLayout = m_calcGraph.getBetweenessInCompleteLayout();
		
		Vector<Double> vBetweeness = new Vector<Double>(betweenessInCompleteLayout.keySet());

		if(vBetweeness.size() <= 0)
		{
			fail("The weight distribution is empty");
		}

		for(int i=0; i < vBetweeness.size(); i++)
		{
			double betweeness = vBetweeness.get(i);
			double betweenessNumber = betweenessInCompleteLayout.get(betweeness);
			
			if(Double.compare(1.0, betweeness) == 0)
			{
				assertTrue("cost of "+betweeness+" there are "+ betweenessNumber,
						betweenessNumber==5);	
			}
			else
			{
				fail("A closeness not anticipated for "+betweeness);
			}
			//System.out.println(" betweeness "+ betweeness + 
			//		" Number of with this betweeness "+	betweenessNumber);
		}
	}

	@Test
	public void testGetCloseness()throws Exception 
	{
		m_graph = TestUtils.createTestGraph7();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Double, Integer> closenessDist = m_calcGraph.getCloseness();
		
		Vector<Double> vCloseness = new Vector<Double>(closenessDist.keySet());

		if(vCloseness.size() <= 0)
		{
			fail("The weight distribution is empty");
		}

		for(int i=0; i < vCloseness.size(); i++)
		{
			double closeness = vCloseness.get(i);
			double closenessNumber = closenessDist.get(closeness);
			
			if(Double.compare(0.5, closeness) == 0)
			{
				assertTrue("cost of "+closeness+" there are "+ closenessNumber,
						closenessNumber==1);	
			}
			else if(Double.compare(0.5714285714285714, closeness) == 0)
			{
				assertTrue("cost of "+closeness+" there are "+ closenessNumber,
						closenessNumber==1);	
			}
			else if(Double.compare(0.6666666666666666, closeness) == 0)
			{
				assertTrue("cost of "+closeness+" there are "+ closenessNumber,
						closenessNumber==2);	
			}
			else if(Double.compare(0.8, closeness) == 0)
			{
				assertTrue("cost of "+closeness+" there are "+ closenessNumber,
						closenessNumber==1);	
			}
			else
			{
				fail("A closeness not anticipated for "+closeness);
			}
			//System.out.println(" closeness "+ closeness + 
			//		" Number of with this closeness "+	closenessNumber);
		}
	}

	@Test
	public void testGetConnectedNodesWithinRadius() throws Exception
	{
		m_graph = TestUtils.createTestGraph6();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Integer, Integer> connectedNodesWithinRadius = m_calcGraph.getConnectedNodesWithinRadius();
		
		Vector<Integer> vConnectedNodesWithinRadius = new Vector<Integer>(connectedNodesWithinRadius.keySet());		
	
		if(vConnectedNodesWithinRadius.size() <= 0)
		{
			fail("The vConnectedNodesWithinRadius distribution is empty");
		}
	
		for(int i=0; i < vConnectedNodesWithinRadius.size(); i++)
		{
			int conNodes = vConnectedNodesWithinRadius.get(i);
			int numconNodes = connectedNodesWithinRadius.get(conNodes);
			
			if(conNodes == 1)
			{
				assertTrue("of node degree "+conNodes+
						"there are "+ numconNodes,
						numconNodes==4);	
			}
			else if(conNodes == 2)
			{
				assertTrue("of node degree "+conNodes+
						"there are "+ numconNodes + "not " +1,
						numconNodes==1);	
			}
			else
			{
				fail("A node with a degree not anticipated for");
			}
			//System.out.println("Number of Nodes within r "+ conNodes +
			//		" Number of Nodes with Number of Nodes within r "+
			//		numconNodes);
		}		

	}

	@Test
	public void testGetCostOfDirectConnectivity() throws Exception
	{
		m_graph = TestUtils.createTestGraph6();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Double, Integer> costDistribution = m_calcGraph.getCostOfDirectConnectivityDist();  
		

		Vector<Double> vCosts = new Vector<Double>(costDistribution.keySet());		

		if(vCosts.size() <= 0)
		{
			fail("The weight distribution is empty");
		}

		for(int i=0; i < vCosts.size(); i++)
		{
			double costOfDirectConnectivity = vCosts.get(i);
			int numCostOfDirectConnectivity = costDistribution.get(costOfDirectConnectivity);
			
			if(Double.compare(3.1800000000000006,costOfDirectConnectivity) == 0)
			{
				assertTrue("weight of "+costOfDirectConnectivity+
						"there are "+ numCostOfDirectConnectivity,
						numCostOfDirectConnectivity==1);	
			}
			else if(Double.compare(5.900000000000001,costOfDirectConnectivity) == 0)
			{
				assertTrue("weight of  "+costOfDirectConnectivity+
						"there are "+ numCostOfDirectConnectivity,
						numCostOfDirectConnectivity==1);	
			}
			else if(Double.compare(8.660000000000002 ,costOfDirectConnectivity) == 0)
			{
				assertTrue("weight of  "+costOfDirectConnectivity+
						"there are "+ numCostOfDirectConnectivity,
						numCostOfDirectConnectivity==1);	
			}	
			else if(Double.compare(12.560000000000002 ,costOfDirectConnectivity) == 0)
			{
				assertTrue("weight of  "+costOfDirectConnectivity+
						"there are "+ numCostOfDirectConnectivity,
						numCostOfDirectConnectivity==1);	
			}	
			else if(Double.compare(7.200000000000003,costOfDirectConnectivity) == 0)
			{
				assertTrue("weight of  "+costOfDirectConnectivity+
						"there are "+ numCostOfDirectConnectivity,
						numCostOfDirectConnectivity==1);	
			}	
			else
			{
			//	fail("A weight not anticipated for"+vWeight.get(i));
				//System.out.println(" weight "+ vWeight.get(i) + " Number of with this weight "+
					//	weightDistribution.get(vWeight.get(i)));
			}
			//System.out.println(" weight "+ costOfDirectConnectivity + " Number of with this weight "+
			//		numCostOfDirectConnectivity);
		}
	}

	@Test
	public void testGetLinksLongerThan_r() throws Exception
	{
		m_graph = TestUtils.createTestGraph6();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Double, Integer> linksLongerThanRDist = m_calcGraph.getLinksLongerThan_r();
		
		Vector<Double> vLength = new Vector<Double>(linksLongerThanRDist.keySet());

		if(vLength.size() <= 0)
		{
			fail("The weight distribution is empty");
		}

		for(int i=0; i < vLength.size(); i++)
		{
			double cost = vLength.get(i);
			double costCount = linksLongerThanRDist.get(cost);
			
			if(Double.compare(0.8485281374238571, cost) == 0)
			{
				assertTrue("cost of "+cost+
						"there are "+ costCount,
						costCount==1);	
			}
			else
			{
				fail("A cost not anticipated for"+cost);

			}
			//System.out.println(" cost "+ cost + 
			//		" Number of with this weight "+	costCount);
		}
		

		
	}

	@Test
	public void testGetMeanDegreeOfConnectedNeighbors() throws Exception
	{
		m_graph = TestUtils.createTestGraph4();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Double, Integer> nodeMeanDegreeDistribution = m_calcGraph.getMeanDegreeOfConnectedNeighbors();
		
		Vector<Double> vAvgDegree = new Vector<Double>(nodeMeanDegreeDistribution.keySet());		
	
		if(vAvgDegree.size() <= 0)
		{
			fail("The node degree distribution is empty");
		}
	
		for(int i=0; i < vAvgDegree.size(); i++)
		{
			double avgDegree = vAvgDegree.get(i);
			int numAvgDegree = nodeMeanDegreeDistribution.get(avgDegree);
			
			if(Double.compare(avgDegree, 10.0/3.0) == 0)
			{
				assertTrue("of node degree "+avgDegree+
						"there are "+ numAvgDegree,
						numAvgDegree==4);	
			}
			else if(Double.compare(vAvgDegree.get(i), 3) == 0)
			{
				assertTrue("of node degree "+avgDegree+
						"there are "+ numAvgDegree + "not " +1,
						numAvgDegree==1);	
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
	public void testGetNodeMaximumEdgeCost() throws Exception
	{
		m_graph = TestUtils.createTestGraph6();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Double, Integer> linksMaxEdgeCostPerNode = m_calcGraph.getNodeMaximumEdgeCost();
		
		Vector<Double> vMaxCost = new Vector<Double>(linksMaxEdgeCostPerNode.keySet());

		if(vMaxCost.size() <= 0)
		{
			fail("The weight distribution is empty");
		}

		for(int i=0; i < vMaxCost.size(); i++)
		{
			double totalCost = vMaxCost.get(i);
			double totalCostCount = linksMaxEdgeCostPerNode.get(totalCost);
			
			


			
			if(Double.compare(0.01999999999999999, totalCost) == 0)
			{
				assertTrue("cost of "+totalCost+" there are "+ totalCostCount,
						totalCostCount==3);	
			}
			else if(Double.compare(1.4400000000000004, totalCost) == 0)
			{
				assertTrue("cost of "+totalCost+" there are "+ totalCostCount,
						totalCostCount==2);	
			}
			else
			{
				fail("A cost not anticipated for"+totalCost);
			}
		//	System.out.println(" totalCost "+ totalCost + 
		//			" Number of with this totalCost "+	totalCostCount);
		}
	}

	@Test
	public void testGetNodeMaximumEdgeLength()throws Exception 
	{
		m_graph = TestUtils.createTestGraph6();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Double, Integer> linksMaxEdgeLengthPerNode = m_calcGraph.getNodeMaximumEdgeLength();
		
		Vector<Double> vMaxLength = new Vector<Double>(linksMaxEdgeLengthPerNode.keySet());

		if(vMaxLength.size() <= 0)
		{
			fail("The weight distribution is empty");
		}

		for(int i=0; i < vMaxLength.size(); i++)
		{
			double totalLength = vMaxLength.get(i);
			double totalCostLength = linksMaxEdgeLengthPerNode.get(totalLength);
			
			if(Double.compare(0.8485281374238571, totalLength) == 0)
			{
				assertTrue("cost of "+totalLength+" there are "+ totalCostLength,
						totalCostLength==2);	
			}
			else if(Double.compare(0.09999999999999998, totalLength) == 0)
			{
				assertTrue("cost of "+totalLength+" there are "+ totalCostLength,
						totalCostLength==3);	
			}
			else
			{
				fail("A cost not anticipated for "+totalLength);
			}
			//System.out.println(" totalCost "+ totalLength + 
			//		" Number of with this totalCost "+	totalCostLength);
		}
	}

	@Test
	public void testGetNodesWithinRadiusById()throws Exception 
	{
		m_graph = TestUtils.createTestGraph6();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Integer, Integer> nodesWithinRadius = m_calcGraph.getNodesWithinRadiusById();
		
		Vector<Integer> vNodesWithinRadius = new Vector<Integer>(nodesWithinRadius.keySet());		
	
		if(vNodesWithinRadius.size() <= 0)
		{
			fail("The vConnectedNodesWithinRadius distribution is empty");
		}
	
		for(int i=0; i < vNodesWithinRadius.size(); i++)
		{
			int conNodes = vNodesWithinRadius.get(i);
			int numconNodes = nodesWithinRadius.get(conNodes);
			
			if(conNodes == 0)
			{
				assertTrue("node id "+conNodes+
						"has this many nodes within R "+ numconNodes,
						numconNodes==2);	
			}
			else if(conNodes ==1)
			{
				assertTrue("node id "+conNodes+
						"has this many nodes within R "+ numconNodes,
						numconNodes==2);	
			}
			else if(conNodes ==2)
			{
				assertTrue("node id "+conNodes+
						"has this many nodes within R "+ numconNodes,
						numconNodes==2);	
			}
			else if(conNodes ==3)
			{
				assertTrue("node id "+conNodes+
						"has this many nodes within R "+ numconNodes,
						numconNodes==1);	
			}
			else if(conNodes ==4)
			{
				assertTrue("node id "+conNodes+
						"has this many nodes within R "+ numconNodes,
						numconNodes==1);	
			}
			else
			{
				fail("A node with a degree not anticipated for "+ conNodes);
			}
//			System.out.println("Node Id  "+ conNodes +
//					" has this many Nodes within r"+
//					numconNodes);
		}		
	}

	@Test
	public void testGetNodesWithinRadius()throws Exception 
	{
		m_graph = TestUtils.createTestGraph6();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		

		
		TreeMap<Integer, Integer> nodesWithinRadius = m_calcGraph.getNodesWithinRadius();
		
		Vector<Integer> vNodesWithinRadius = new Vector<Integer>(nodesWithinRadius.keySet());		
	
		if(vNodesWithinRadius.size() <= 0)
		{
			fail("The vConnectedNodesWithinRadius distribution is empty");
		}
	
		for(int i=0; i < vNodesWithinRadius.size(); i++)
		{
			int conNodes = vNodesWithinRadius.get(i);
			int numconNodes = nodesWithinRadius.get(conNodes);
			
			if(conNodes == 1)
			{
				assertTrue("of node degree "+conNodes+
						"there are "+ numconNodes,
						numconNodes==2);	
			}
			else if(conNodes ==2)
			{
				assertTrue("of node degree "+conNodes+
						"there are "+ numconNodes + "not " +1,
						numconNodes==3);	
			}
			else
			{
				fail("A node with a degree not anticipated for "+ conNodes);
			}
//			System.out.println("Number of Connected Nodes within r "+ conNodes +
//					"Number of Nodes with Number of Connected Nodes within r"+
//					numconNodes);
		}		
	}
	
	

	@Test
	public void testGetNodeTotalEdgeCost() throws Exception
	{
		m_graph = TestUtils.createTestGraph6();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Double, Integer> linksTotalEdgeCostPerNode = m_calcGraph.getNodeTotalEdgeCost();
		
		Vector<Double> vTotalCost = new Vector<Double>(linksTotalEdgeCostPerNode.keySet());

		if(vTotalCost.size() <= 0)
		{
			fail("The weight distribution is empty");
		}

		for(int i=0; i < vTotalCost.size(); i++)
		{
			double totalCost = vTotalCost.get(i);
			double totalCostCount = linksTotalEdgeCostPerNode.get(totalCost);
			
			
			if(Double.compare(0.01999999999999999, totalCost) == 0)
			{
				assertTrue("cost of "+totalCost+" there are "+ totalCostCount,
						totalCostCount==3);	
			}
			else if(Double.compare(1.4600000000000004, totalCost) == 0)
			{
				assertTrue("cost of "+totalCost+" there are "+ totalCostCount,
						totalCostCount==1);	
			}
			else if(Double.compare(1.4800000000000004, totalCost) == 0)
			{
				assertTrue("cost of "+totalCost+" there are "+ totalCostCount,
						totalCostCount==1);	
			}
			else
			{
				fail("A cost not anticipated for"+totalCost);
			}
			//System.out.println(" totalCost "+ totalCost + 
			//		" Number of with this totalCost "+	totalCostCount);
		}
	}

	@Test
	public void testGetNodeTotalEdgeLength()throws Exception 
	{
		m_graph = TestUtils.createTestGraph6();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		TreeMap<Double, Integer> linksTotalEdgeLengthPerNode = m_calcGraph.getNodeTotalEdgeLength();
		
		Vector<Double> vTotalLength = new Vector<Double>(linksTotalEdgeLengthPerNode.keySet());

		if(vTotalLength.size() <= 0)
		{
			fail("The weight distribution is empty");
		}

		for(int i=0; i < vTotalLength.size(); i++)
		{
			double totalLength = vTotalLength.get(i);
			double totalCostLength = linksTotalEdgeLengthPerNode.get(totalLength);
			
			if(Double.compare(0.09999999999999998, totalLength) == 0)
			{
				assertTrue("cost of "+totalLength+" there are "+ totalCostLength,
						totalCostLength==3);	
			}
			else if(Double.compare(0.9485281374238571, totalLength) == 0)
			{
				assertTrue("cost of "+totalLength+" there are "+ totalCostLength,
						totalCostLength==1);	
			}
			else if(Double.compare(1.048528137423857, totalLength) == 0)
			{
				assertTrue("cost of "+totalLength+" there are "+ totalCostLength,
						totalCostLength==1);	
			}
			else
			{
				fail("A cost not anticipated for "+totalLength);
			}
			//System.out.println(" totalCost "+ totalLength + 
			//		" Number of with this totalCost "+	totalCostLength);
		}
	}

}
