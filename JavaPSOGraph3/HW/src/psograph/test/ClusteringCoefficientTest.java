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



public class ClusteringCoefficientTest 
{

	Graph m_graph;
	
	@Before
	public void setUp() throws Exception 
	{

	}
	
	@Test
	public void testCompleteGraph() throws Exception
	{
		m_graph = new Graph(0);
		int m_id1 = m_graph.addNode();
		int m_id2 = m_graph.addNode();
		int m_id3 = m_graph.addNode();
		int m_id4 = m_graph.addNode();
		
		m_graph.addConnection(m_id1, m_id2);
		m_graph.addConnection(m_id1, m_id3);
		m_graph.addConnection(m_id1, m_id4);
		
		m_graph.addConnection(m_id2, m_id3);
		m_graph.addConnection(m_id2, m_id4);
		
		m_graph.addConnection(m_id3, m_id4);
		
		CalculatedGraph calcGraph= new CalculatedGraph(m_graph);
		calcGraph.setCostBasis(5);
		calcGraph.UpdateCalcuations();
		
		double factor = calcGraph.getClusteringCoefficient();
		
		assertTrue("factor is not equal to 1 but is equal to "+factor, Double.compare(1.0, factor) == 0);

		
	}
	@Test
	public void testGraph1() throws Exception
	{
		m_graph = new Graph(0);
		int m_id1 = m_graph.addNode();
		int m_id2 = m_graph.addNode();
		int m_id3 = m_graph.addNode();
		int m_id4 = m_graph.addNode();
		
		m_graph.addConnection(m_id1, m_id2);
		m_graph.addConnection(m_id1, m_id3);
		m_graph.addConnection(m_id1, m_id4);
		
		
		m_graph.addConnection(m_id3, m_id4);
		
		CalculatedGraph calcGraph= new CalculatedGraph(m_graph);
		calcGraph.setCostBasis(5);
		calcGraph.UpdateCalcuations();
		

		
		double factor = calcGraph.getClusteringCoefficient();
		@SuppressWarnings("unused")
		double ans = (1.0/3.0 + 2.0)/4.0; //0.5833333333333334 


		assertTrue("factor is not equal to " + 0.5833333333333333 +" but is equal to "+factor, Double.compare(0.5833333333333333, factor) == 0);
				
		
		
	}
	@Test
	public void testUnconnectedGraph() throws Exception
	{
		m_graph = new Graph(0);
		@SuppressWarnings("unused")
		int m_id1 = m_graph.addNode();
		@SuppressWarnings("unused")
		int m_id2 = m_graph.addNode();
		@SuppressWarnings("unused")
		int m_id3 = m_graph.addNode();
		@SuppressWarnings("unused")
		int m_id4 = m_graph.addNode();
		
		CalculatedGraph calcGraph= new CalculatedGraph(m_graph);
		calcGraph.setCostBasis(5);
		calcGraph.UpdateCalcuations();
		
		double factor = calcGraph.getClusteringCoefficient();
		
		assertTrue("factor is not equal to 0 ", Double.compare(0, factor) == 0);
				
		
		
	}
	

}
