// This is a library to be used to represent a Graph and various measurments for a Graph
//  and to perform optimization using Particle Swarm Optimization (PSO)
//    Copyright (C) 2008, 2015 Patrick Olekas
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
import psograph.measurements.ASPL;


public class ASPLTest  {

	Graph m_graph;
	CalculatedGraph m_calcGraph;
	
	@Before
	public void setUp() throws Exception 
	{
		m_graph = TestUtils.createTestGraph2();

		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
	}
	
	@Test
	public void testASPLDisconnected()  throws Exception
	{
		
		boolean isFullyConnected = m_calcGraph.isFullyConnected();
		
		assertTrue("Is graph fully connected should be false, not " + isFullyConnected,
				isFullyConnected == false);		
		
		ASPL ASPL_measure = new ASPL(m_graph);
		double m_ASPL = ASPL_measure.Measure();
		TreeMap<Integer, Integer> m_distributionSPL2 = ASPL_measure.getSPLDistribution();		

		double value = m_calcGraph.getASPL();
		
		assertTrue("ASPL should be -1 not " + value, Double.compare(value, -1) ==0);

		TreeMap<Integer, Integer> SPLDist = m_calcGraph.getDistributionSPL();
		
		assertTrue("SPL distribution should be NULL " , m_distributionSPL2 == null);
		
		System.out.println("Printing out original Graph");
		this.m_graph.printWithLocationAndWeights();
		
		

	}
	
	@Test
	public void testASPL()  throws Exception
	{
		m_graph = TestUtils.createSimpleRingGraph();
		
		boolean isFullyConnected = m_graph.isFullyConnected();
		
		assertTrue("Is graph fully connected should be true, not " + isFullyConnected,
				isFullyConnected == true);	
	
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		double value = m_calcGraph.getASPL();
				
		assertTrue("Unit test ASPL of network should be 1.5, not " + value,
				Double.compare(value,1.5) == 0.0);
		

		TreeMap<Integer, Integer> SPLDist = m_calcGraph.getDistributionSPL();

	
		Vector<Integer> vSPL = new Vector<Integer>(SPLDist.keySet());		

		if(vSPL.size() <= 0)
		{
			fail("The SPL distribution is empty");
		}

		for(int i=0; i < vSPL.size(); i++)
		{
			
			if(vSPL.get(i) == 1)
			{
				assertTrue("SPL of "+vSPL.get(i)+
						"there are "+ SPLDist.get(vSPL.get(i)),
						SPLDist.get(vSPL.get(i))==10);	
			}
			else if(vSPL.get(i) == 2)
			{
				assertTrue("SPL of  "+vSPL.get(i)+
						"there are "+ SPLDist.get(vSPL.get(i)),
						SPLDist.get(vSPL.get(i))==10);	
			}
			else
			{
				fail("A SPL not anticipated for"+vSPL.get(i));
			}
			//System.out.println("SPL weight "+ vSPL.get(i) + "Number of SPL with this weight "+
			//		SPLDist.get(vSPL.get(i)));
			
			
			
			
			ASPL ASPL_measure = new ASPL(m_graph);
			double m_ASPL = ASPL_measure.Measure();
			TreeMap<Integer, Integer> m_distributionSPL2 = ASPL_measure.getSPLDistribution();		

		
						
			System.out.println("Printing out original Graph");
			this.m_graph.printWithLocationAndWeights();
			
			
			System.out.println("Printing out aspl Graph");
			ASPL_measure.calcValues.printWithLocationAndWeights();
			
		}
	}

}
