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

import psograph.graph.CalculatedGraph;
import psograph.graph.Graph;


/**
 * @author Owner
 *
 */
public class CalculatedGraphTest2 {

	Graph m_graph;
	CalculatedGraph m_calcGraph;
	
	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception 
	{
		m_graph = TestUtils.createTestGraph3();
		
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(56);
		m_calcGraph.UpdateCalcuations();
		
	}

	/**
	 * Test method for {@link psograph.graph.CalculatedGraph#getAssortativity()}.
	 */
	@Test
	public void testGetAssortativity() throws Exception
	{
		m_graph = TestUtils.createTestGraph4();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		double value = m_calcGraph.getAssortativity();
		



		assertTrue("Assortativity is  "+value, Double.compare(value,0.9734660033167496)==0); 
	}
	

	/**
	 * Test method for {@link psograph.graph.CalculatedGraph#getClumpinessOfLayout()}.
	 */
	@Test
	public void testGetClumpinessOfLayout() throws Exception
	{
		m_graph = TestUtils.createTestGraph13();
		
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		
		TreeMap<Integer, Integer> clumpinessDist = m_calcGraph.getClumpinessOfLayout();

		
		Vector<Integer> vClumpiness = new Vector<Integer>(clumpinessDist.keySet());		

		if(vClumpiness.size() <= 0)
		{
			fail("The ISPL distribution is empty");
		}

		for(int i=0; i < vClumpiness.size(); i++)
		{
			int clumpiness = vClumpiness.get(i);
			int numOfClumpiness = clumpinessDist.get(clumpiness);
			
			System.out.println("NUmber of nodes in cell :"+clumpiness+" numOfCells with  "+numOfClumpiness);
			
			if(clumpiness == 0)
			{
				assertTrue("clumpiness of "+clumpiness+
						"there are "+ numOfClumpiness,
						numOfClumpiness==81);	
			}
			else if(clumpiness==1)
			{
				assertTrue("clumpiness"+clumpiness+
						"there are "+ numOfClumpiness,
						numOfClumpiness==11);	
			}
			else if(clumpiness==2)
			{
				assertTrue("clumpiness of  "+clumpiness+
						"there are "+ numOfClumpiness,
						numOfClumpiness==5);	
			}
			else if(clumpiness== 3)
			{
				assertTrue("clumpiness of  "+clumpiness+
						"there are "+ numOfClumpiness,
						numOfClumpiness==3);	
			}
			else
			{
				fail("A clumpiness not anticipated for"+clumpiness);
			}
			//System.out.println("clumpiness "+ clumpiness + "Number of SPL with this weight "+
			//		numOfClumpiness);
		}
	}

	/**
	 * Test method for {@link psograph.graph.CalculatedGraph#getFitnessOfSeedWithFullConnectivity()}.
	 */
	@Test
	public void testGetFitnessOfSeedWithFullConnectivity()  throws Exception
	{
		double value = m_calcGraph.getFitnessOfSeedWithFullConnectivity();
		


		assertTrue("FitnessOfSeedWithFullConnectivity is  "+value,
				Double.compare(value, -1.6263151721194329)==0); 

	}

	/**
	 * Test method for {@link psograph.graph.CalculatedGraph#getHubbiness()}.
	 */
	@Test
	public void testGetHubbiness() {
		
	//	double stdDev = m_calcGraph.getNodeDegreeStdDeviation();
	//	double variance = m_calcGraph.getNodeDegreeVariance();
	//	double mean = m_calcGraph.getAvgNodeDegree();
		
	//	double k = 3.0;
		
	//	double clumpCuttoff = k * stdDev;
		
		//TreeMap<Integer, Integer> nodeDistribution = m_calcGraph.getNodeDistribution();
		
		//Vector<Integer> vDegree = new Vector<Integer>(nodeDistribution.keySet());		
	
		//for(int i=0; i < vDegree.size(); i++)
		//{

			//System.out.println("Degree "+ vDegree.get(i) + "Nodes with this degree "+
			//		nodeDistribution.get(vDegree.get(i)));
		//}	
		
		double value = m_calcGraph.getHubbiness();
		
		assertTrue("Hubbiness is  "+value, Double.compare(value,1.0)==0); 
		
	}

	/**
	 * Test method for {@link psograph.graph.CalculatedGraph#getHubbiness()}.
	 */
	@Test
	public void testGetPercentNodesConnectedToHub() throws Exception {
		
	//	double stdDev = m_calcGraph.getNodeDegreeStdDeviation();
	//	double variance = m_calcGraph.getNodeDegreeVariance();
	//	double mean = m_calcGraph.getAvgNodeDegree();
		
	//	double k = 3.0;
		
	//	double clumpCuttoff = k * stdDev;
		
		//TreeMap<Integer, Integer> nodeDistribution = m_calcGraph.getNodeDistribution();
		
		//Vector<Integer> vDegree = new Vector<Integer>(nodeDistribution.keySet());		
	
		//for(int i=0; i < vDegree.size(); i++)
		//{

			//System.out.println("Degree "+ vDegree.get(i) + "Nodes with this degree "+
			//		nodeDistribution.get(vDegree.get(i)));
		//}	
		
		double value = m_calcGraph.getPercentNodesConnectedToHub();
		
		//m_calcGraph.print();
		
		assertTrue("PercentNodesConnectedToHub is  "+value, Double.compare(value,.7)==0); 
		
	}	
	
}
