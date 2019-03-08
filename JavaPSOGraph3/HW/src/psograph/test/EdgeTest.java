package psograph.test;

import static org.junit.Assert.*;

import org.junit.Test;

import psograph.graph.Edge;

public class EdgeTest {



	@Test
	public void testEdgeEdge() 
	{
		Edge e = new Edge(.55);
		Edge e2 = new Edge(e);
		
		assertTrue(Double.compare(e2.getWeight(),.55) == 0);
	}

	@Test
	public void testEdgeDouble() 
	{
		Edge e = new Edge(.55);
		
		assertTrue(Double.compare(e.getWeight(),.55) == 0);
	}

	@Test
	public void testModifyWeight() 
	{
		Edge e = new Edge(.55);
		
		assertTrue(Double.compare(e.getWeight(),.55) == 0);
		
		e.modifyWeight(.76);
		
		assertTrue(Double.compare(e.getWeight(),.76) == 0);
	}

//	@Test
//	public void testGetWeight() {
//		fail("Not yet implemented");
//	}

}
