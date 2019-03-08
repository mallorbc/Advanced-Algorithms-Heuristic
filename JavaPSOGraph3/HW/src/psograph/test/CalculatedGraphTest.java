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

import org.junit.Before;
import org.junit.Test;

import psograph.graph.CalculatedGraph;
import psograph.graph.Graph;

import java.util.TreeMap;
import java.util.Vector;




/**
 * @author Owner
 *
 */
public class CalculatedGraphTest  {

	Graph m_graph;
	CalculatedGraph m_calcGraph;
	
	@Before
	public void setUp() throws Exception 
	{

		m_graph = TestUtils.createTestGraph();
		
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
	}
	
	@Test
	public void testNodeDistribution() throws Exception
	{
		
		
		TreeMap<Integer, Integer> nodeDistribution = m_calcGraph.getNodeDistribution();
		
		Vector<Integer> vDegree = new Vector<Integer>(nodeDistribution.keySet());		
	
		if(vDegree.size() <= 0)
		{
			fail("The node degree distribution is empty");
		}
	
		for(int i=0; i < vDegree.size(); i++)
		{
			
			if(vDegree.get(i) == 2)
			{
				assertTrue("of node degree "+vDegree.get(i)+
						"there are "+ nodeDistribution.get(vDegree.get(i)),
						nodeDistribution.get(vDegree.get(i))==3);	
			}
			else if(vDegree.get(i) == 3)
			{
				assertTrue("of node degree "+vDegree.get(i)+
						"there are "+ nodeDistribution.get(vDegree.get(i)),
						nodeDistribution.get(vDegree.get(i))==3);	
			}
			else if(vDegree.get(i) == 4)
			{
				assertTrue("of node degree "+vDegree.get(i)+
						"there are "+ nodeDistribution.get(vDegree.get(i)),
						nodeDistribution.get(vDegree.get(i))==3);	
			}
			else if(vDegree.get(i) == 5)
			{
				assertTrue("of node degree "+vDegree.get(i)+
						"there are "+ nodeDistribution.get(vDegree.get(i)),
						nodeDistribution.get(vDegree.get(i))==1);	
			}
			else
			{
				fail("A node with a degree not anticipated for");
			}
			//System.out.println("Degree "+ vDegree.get(i) + "Nodes with this degree "+
			//		nodeDistribution.get(vDegree.get(i)));
		}
		
		
		
		
		
	}

	@Test
	public void testNodeDegreeVariance() throws Exception
	{
		
		
		double nodeDegreeVariance = m_calcGraph.getNodeDegreeVariance() ;
		
		


		assertTrue("node degree variance of " +nodeDegreeVariance, Double.compare(nodeDegreeVariance, 1.0666666666666673)==0);			
		
		
	}
	@Test
	public void testNodeDegreeStdDeviation() throws Exception
	{
		
		
		double nodeDegreeStdDeviation = m_calcGraph.getNodeDegreeStdDeviation() ;
		
		
		
		assertTrue("node degree std dev of " +nodeDegreeStdDeviation, Double.compare(nodeDegreeStdDeviation, 1.0327955589886448)==0);			
		
		
		
	}
	@Test
	public void testAvgNodeDegree() throws Exception
	{
		double m_avgNodeDegree = m_calcGraph.getAvgNodeDegree();
		

		assertTrue("Avg node degree is "+ m_avgNodeDegree,Double.compare(m_avgNodeDegree, 3.2)==0);			
		
		
	}
	
	
}
