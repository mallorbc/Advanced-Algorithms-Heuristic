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


import org.junit.Before;
import org.junit.Test;

import psograph.graph.*;

import static org.junit.Assert.*;

public class FitnessValues 
{
	
	Graph m_graph;
	
	@Before
	public void setUp() throws Exception 
	{
		m_graph = TestUtils.createTestGraph2();
	}
	
	@Test
	public void testOldFitnessValueMethod() throws Exception
	{
		
		CalculatedGraph calcGraph = new CalculatedGraph(m_graph);
		calcGraph.setCostBasis(5);
		calcGraph.UpdateCalcuations();
			
		//System.out.println("----------------Begin Measurements-------------");
		@SuppressWarnings("unused")
		
		double Rpercentage = calcGraph.getRandomLCC();
		//System.out.println("Avg Robustness Measure for Random - Percentage in LCC "+Rpercentage);
		
		@SuppressWarnings("unused")
		double Rdiameter = calcGraph.getRandomDiameter();
		//System.out.println("Avg Robustness Measure for Random - diameter "+Rdiameter);
		
		
		@SuppressWarnings("unused")
		double Dpercentage = calcGraph.getDirectLCC();
		//System.out.println("Avg Robustness Measure for Random - Percentage in LCC "+Dpercentage);
		
		@SuppressWarnings("unused")
		double Ddiameter = calcGraph.getDirectDiameter();
	//	System.out.println("Avg Robustness Measure for Random - diameter "+Ddiameter);
		
		double aispl=0;
		aispl = calcGraph.getAISPL();
		
		//System.out.println("Connectivity Measure - AISPL "+ aispl);


		assertTrue("APISL  is " + aispl,  
				Double.compare(aispl, 0.3371212121212121) ==0);

		
		double cost = m_graph.SumAllWeights()/2;
		
		//m_graph.printWithLocationAndWeights();
		
		//System.out.println("Cost Measure - summation weight costs "+ cost);	
		assertTrue("cost should be 52, but is " + cost, cost == 52);
		//System.out.println("------------------End Measurements-------------");
		
		//System.out.println("Fitness Value - Ran "+copy.calcRFitnessValue());
		//System.out.println("Fitness Value - Dir "+copy.calDRFitnessValue());
		
		
		//copy.printWithWeights();
	}
	
	
}
