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

import psograph.graph.CostBasis;
import psograph.graph.Graph;



public class CostBasisTests 
{



	Graph m_graph;
	
	@Before
	public void setUp() throws Exception 
	{

		m_graph = TestUtils.createEmptyGraph();
		
		
		@SuppressWarnings("unused")
		int id0 = m_graph.addNodeLocationDataOnly(0,circleEquationPositive(0));
		@SuppressWarnings("unused")
		int id1 = m_graph.addNodeLocationDataOnly(.1,circleEquationPositive(.1));
		@SuppressWarnings("unused")
		int id2 = m_graph.addNodeLocationDataOnly(.2,circleEquationPositive(.2));
		@SuppressWarnings("unused")
		int id3 = m_graph.addNodeLocationDataOnly(.3,circleEquationPositive(.3));
		@SuppressWarnings("unused")
		int id4 = m_graph.addNodeLocationDataOnly(.4,circleEquationPositive(.4));
		@SuppressWarnings("unused")
		int id5 = m_graph.addNodeLocationDataOnly(.5,circleEquationPositive(.5));			
		@SuppressWarnings("unused")
		int id6 = m_graph.addNodeLocationDataOnly(.6,circleEquationPositive(.6));
		@SuppressWarnings("unused")
		int id7 = m_graph.addNodeLocationDataOnly(.7,circleEquationPositive(.7));
		@SuppressWarnings("unused")
		int id8 = m_graph.addNodeLocationDataOnly(.8,circleEquationPositive(.8));
		@SuppressWarnings("unused")
		int id9 = m_graph.addNodeLocationDataOnly(.9,circleEquationPositive(.9));
		@SuppressWarnings("unused")
		int id10 = m_graph.addNodeLocationDataOnly(1,circleEquationPositive(1));
		@SuppressWarnings("unused")
		int id11 = m_graph.addNodeLocationDataOnly(.5, 0);
		@SuppressWarnings("unused")
		int id12 = m_graph.addNodeLocationDataOnly(0, 0);
		
	}

	public double circleEquationPositive(double x)
	{
		
		// eq of circle 
		// (x-.5)^2 + (y-.5)^2 =.5^2
		// (y-.5)^2 = .5^2 - (x -.5)^2
		// y = SQRT(.5^2- (x-.5)^2) + .5
		// y = SQRT(.5^2- (x-.5)^2) + .5

		double y =0;
		
		y = Math.sqrt(.25- (x-.5)*(x-.5)) + .5;
		
		return y;
	}
	@Test
	public void testLocationUniqueness()
	{
		
		try 
		{
			m_graph.addNodeLocationDataOnly(0, 0);
			fail("Exception should have been thrown");

		}
		catch(Exception e)
		{
			//Should have an exception for this case
			//e.printStackTrace();
			
		}
		
	}
	@Test
	public void testCostBasisCalculationsInputChecking() throws Exception
	{
		try 
		{
			CostBasis cb = new CostBasis(m_graph);
			cb.generate(12);
			fail("exception should have been thrown");
		}
		catch(Exception e)
		{
			
		}
		try 
		{
			CostBasis cb = new CostBasis(m_graph);
			cb.generate((13*14) +1);
			fail("exception should have been thrown");
		}
		catch(Exception e)
		{
			
		}
	}
	
//	public void testCostBasisCalculationsMinimal() throws Exception
//	{
//		//cand 1
//		//int id13 = m_graph.addNodeLocationDataOnly(1, 0);
//
//		
//		CostBasis cb = new CostBasis(m_graph);
//		cb.generate(13);
//
//		
//		assertTrue("N is not equal to 13",cb.getN() == 13);
//		// 29.892622072208884
//		assertTrue("cost is not equal to 29.892622072208884, but "+cb.getCost(),Double.compare(cb.getCost(),29.892622072208884) == 0.0);
//		assertTrue("edges is not equal to 13",cb.getTotalEdges() == 13);
//		/*		
//		cb.printWithLocationAndWeights();
//		
//	
//		Node 0 loc (0.0,0.5) neighbors 10:2.0
//		Node 1 loc (0.1,0.7999999999999999) neighbors 10:1.805
//		Node 2 loc (0.2,0.9) neighbors 11:1.805 12:1.6928
//		Node 3 loc (0.3,0.9582575694955839) neighbors 12:2.0
//		Node 4 loc (0.4,0.9898979485566356) neighbors 12:2.2898
//		Node 5 loc (0.5,1.0) neighbors 11:2.0 12:2.5088000000000004
//		Node 6 loc (0.6,0.9898979485566356) neighbors 12:2.6912
//		Node 7 loc (0.7,0.958257569495584) neighbors 12:2.8322
//		Node 8 loc (0.8,0.8999999999999999) neighbors 12:2.88
//		Node 9 loc (0.9,0.7999999999999999) neighbors 12:2.88
//		Node 10 loc (1.0,0.5) neighbors 0:2.0 1:1.805 12:2.5088000000000004
//		Node 11 loc (0.5,0.0) neighbors 2:1.805 5:2.0
//		Node 12 loc (0.0,0.0) neighbors 2:1.6928 3:2.0 4:2.2898 5:2.5088000000000004 6:2.6912 7:2.8322 8:2.88 9:2.88 10:2.5088000000000004
//*/
//		
//		//Util.printMFile(cb);		
//		//First ensure we are getting all the long connections
//		//And then connect those nodes together
//		//
//		
//		
//	}
	@Test
	public void testCostBasisCalculationsTypical() throws Exception
	{
		//cand 1
		//int id13 = m_graph.addNodeLocationDataOnly(1, 0);


		
		CostBasis cb = new CostBasis(m_graph);
		cb.generate(13*2);
		
		assertTrue("N is not equal to 13",cb.getNumberOfNodes() == 13);
		//assertTrue("cost is not equal to 29.8936",cb.getCost() == 29.8936);
		int totalEdges = cb.getTotalEdges();
		assertTrue("edges is not equal to 26, not "+totalEdges,totalEdges == 26);
		
		// Util.printMFile(cb);		
		
		//TODO we should at least check that it did connected to the longest links possible.
		
		//First ensure we are getting all the long connections
		//And then connect those nodes together
		//
		
	}
	
}
