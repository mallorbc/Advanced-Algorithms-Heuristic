/**
 * 
 */
package psograph.test;

import static org.junit.Assert.*;

import java.util.TreeMap;
import java.util.Vector;

import org.junit.Test;

import psograph.graph.Edge;
import psograph.graph.Node;

/**
 * @author Patrick
 *
 */
public class NodeTest  
{


	
	/** The following methods do not have Unit Tests since 
	 * we use them to validate the other methods
	 * Test method for {@link psograph.graph.Node#getDegree()}.
	 * Test method for {@link psograph.graph.Node#getDepth()}.
	 * Test method for {@link psograph.graph.Node#getID()}.
	 * Test method for {@link psograph.graph.Node#getMeanEdgeCost()}.
	 * Test method for {@link psograph.graph.Node#getTotalEdgeCost()}.
	 * Test method for {@link psograph.graph.Node#getX()}.
	 * Test method for {@link psograph.graph.Node#getY()}.
	 * Test method for {@link psograph.graph.Node#getNeighborDistribution()}.
	 * Test method for {@link psograph.graph.Node#getNeighbors()}.
	 */	
	
	/** These will be tested in the SPL testcase
	 * Test method for {@link psograph.graph.Node#getPaths()}.
	 * Test method for {@link psograph.graph.Node#getLongestPathLength()}.
	 * Test method for {@link psograph.graph.Node#getSPLength()}.
	 * Test method for {@link psograph.graph.Node#getVisited()}.
	 * Test method for {@link psograph.graph.Node#hasSPLTo(int)}.
	 * Test method for {@link psograph.graph.Node#numberOfSPLNodeParticipatesIn(int, psograph.graph.Node)}.
	 * Test method for {@link psograph.graph.Node#numOfSPLength()}.
	 * Test method for {@link psograph.graph.Node#setDepth(int)}.
	 * Test method for {@link psograph.graph.Node#setPath(java.util.Vector)}.
	 * Test method for {@link psograph.graph.Node#setVisited(boolean)}.
	 * Test method for {@link psograph.graph.Node#addPath(psograph.graph.Path)}.
	 * */
	


	/**
	 * Test method for {@link psograph.graph.Node#Node(int, double, double)}.
	 */
	@Test
	public void testNodeIntDoubleDouble() 
	{
		Node n = new Node(1, .5, .3);
		
		assertTrue(n.getID() == 1);
		
		assertTrue(Double.compare(.5, n.getX()) == 0);
		
		assertTrue(Double.compare(.3, n.getY()) == 0);
		
		assertTrue(n.getDegree() == 0);
		
		assertTrue(n.getDepth() == -1);
		
		assertTrue(n.getLongestPathLength()==-1);
		
		assertTrue(n.getMeanEdgeCost() == 0);
		
		assertTrue(n.getNeighborDistribution() == null);
		
		assertTrue(n.getNeighbors() == null);
		
		assertTrue(n.getPaths().size() == 0);
		
		assertTrue(n.getSPLength() == -1);
		
		assertTrue(n.getTotalEdgeCost() == 0);
		
		assertTrue(n.getVisited()==false);
		
	}

	/**
	 * Test method for {@link psograph.graph.Node#Node(int, double, double, int[])}.
	 */
	@Test
	public void testNodeIntDoubleDoubleIntArray() throws Exception
	{
		int connectedTo[] = null;
		
		Node n = new Node(1, .5, .3, connectedTo);
		
		assertTrue(n.getID() == 1);
		
		assertTrue(Double.compare(.5, n.getX()) == 0);
		
		assertTrue(Double.compare(.3, n.getY()) == 0);
		
		assertTrue(n.getDegree() == 0);
		
		assertTrue(n.getDepth() == -1);
		
		assertTrue(n.getLongestPathLength()==-1);
		
		assertTrue(n.getMeanEdgeCost() == 0);
		
		assertTrue(n.getNeighborDistribution() == null);
		
		assertTrue(n.getNeighbors() == null);
		
		assertTrue(n.getPaths().size() == 0);
		
		assertTrue(n.getSPLength() == -1);
		
		assertTrue(n.getTotalEdgeCost() == 0);
		
		assertTrue(n.getVisited()==false);
	}		

	/**
	 * Test method for {@link psograph.graph.Node#Node(int, double, double, int[])}.
	 */
	@Test
	public void testNodeIntDoubleDoubleIntArray2() throws Exception
	{
		int connectedTo[] = {1,2,3,4};
		try
		{
			@SuppressWarnings("unused")
			Node n = new Node(1, .5, .3, connectedTo);
			fail("Test shoudl have thrown exception");
		}
		catch(Exception e)
		{
			
			//Good as expected
		}

	}		
	
	/**
	 * Test method for {@link psograph.graph.Node#Node(int, double, double, int[])}.
	 */
	@Test
	public void testNodeIntDoubleDoubleIntArray3() throws Exception
	{
		int connectedTo[] = {2,3,4};
		
		Node n = new Node(1, .5, .3, connectedTo);
		
		assertTrue(n.getID() == 1);
		
		assertTrue(Double.compare(.5, n.getX()) == 0);
		
		assertTrue(Double.compare(.3, n.getY()) == 0);
		
		assertTrue(n.getDegree() == 3);
		
		assertTrue(n.getDepth() == -1);
		
		assertTrue(n.getLongestPathLength()==-1);
		
		assertTrue(n.getMeanEdgeCost() == 1);
		
		TreeMap<Double,Vector<Integer>> neighborDist = n.getNeighborDistribution();
		
		assertTrue(neighborDist != null);
		
		assertTrue(neighborDist.keySet().size() == 1);
		
		Vector<Integer> nodeIds = neighborDist.get(1.0);
		
		assertTrue(nodeIds.contains(2)== true);
		assertTrue(nodeIds.contains(3)== true);
		assertTrue(nodeIds.contains(4)== true);
	
		TreeMap<Integer, Edge> neighbors = n.getNeighbors();
		
		assertTrue(neighbors.size() == 3);
		
		assertTrue(neighbors.get(2).getWeight() == 1);
		assertTrue(neighbors.get(3).getWeight() == 1);
		assertTrue(neighbors.get(4).getWeight() == 1);
		
		assertTrue(n.getPaths().size() == 0);
		
		assertTrue(n.getSPLength() == -1);
		
		assertTrue(n.getTotalEdgeCost() == 3);
		
		assertTrue(n.getVisited()==false);
	}		
	


	/**
	 * Test method for {@link psograph.graph.Node#addConnection(int, double)}.
	 * @throws Exception 
	 */
	@Test
	public void testAddConnection() throws Exception 
	{
		Node n = new Node(1, .5, .4);
		
		n.addConnection(2, .28);
		try
		{
			n.addConnection(2, .28);
		}
		catch (Exception e)
		{
			//from adding in with duplicate ID
		}
		
		n.addConnection(3, .18);
		
		n.addConnection(5, .08);
		
		assertTrue(n.getID() == 1);
		
		assertTrue(Double.compare(.5, n.getX()) == 0);
		
		assertTrue(Double.compare(.4, n.getY()) == 0);
		
		assertTrue(n.getDegree() == 3);
		
		assertTrue(n.getDepth() == -1);
		
		assertTrue(n.getLongestPathLength()==-1);
		
		assertTrue(Double.compare(n.getMeanEdgeCost(),0.18000000000000002)==0);
		
		TreeMap<Double,Vector<Integer>> neighborDist = n.getNeighborDistribution();
		
		assertTrue(neighborDist != null);
		
		assertTrue(neighborDist.keySet().size() == 3);
		
		Vector<Integer> nodeIds = neighborDist.get(.08);
		
		assertTrue(nodeIds.contains(5)== true);

		 nodeIds = neighborDist.get(.18);
		
		assertTrue(nodeIds.contains(3)== true);

		nodeIds = neighborDist.get(.28);
		
		assertTrue(nodeIds.contains(2)== true);
	
		
		TreeMap<Integer, Edge> neighbors = n.getNeighbors();
		
		assertTrue(neighbors.size() == 3);
		
		assertTrue(Double.compare(neighbors.get(5).getWeight() , .08)==0);
		assertTrue(Double.compare(neighbors.get(3).getWeight() , .18)==0);
		assertTrue(Double.compare(neighbors.get(2).getWeight() , .28)==0);
		
		assertTrue(n.getPaths().size() == 0);
		
		assertTrue(n.getSPLength() == -1);
		
		assertTrue(n.getTotalEdgeCost() == .54);
		
		assertTrue(n.getVisited()==false);		
		
		
	}



	/**
	 * Test method for {@link psograph.graph.Node#equals(java.lang.Object)}.
	 * @throws Exception 
	 */
	@Test
	public void testEqualsObject() throws Exception 
	{
		Node n1 = new Node(1, .5, .4);
		n1.addConnection(3, .28);
		n1.addConnection(4, .18);
		n1.addConnection(6, .08);
	
		Node n2 = new Node(1, .5, .4);
		n2.addConnection(3, .28);
		n2.addConnection(4, .18);
		n2.addConnection(6, .08);
		
		Node n3 = new Node(1, .4, .4);
		n3.addConnection(3, .28);
		n3.addConnection(4, .18);
		n3.addConnection(6, .08);
		
		Node n3a = new Node(1, .5, .2);
		n3a.addConnection(3, .28);
		n3a.addConnection(4, .18);
		n3a.addConnection(6, .08);	
		
		Node n4 = new Node(1, .5, .4);
		n4.addConnection(3, .28);
		n4.addConnection(6, .08);		
	
		Node n5 = new Node(1, .5, .4);
		n5.addConnection(3, .18);
		n5.addConnection(4, .18);
		n5.addConnection(6, .08);			
		
		Node n6 = new Node(1, .5, .4);
		n6.addConnection(3, .28);
		n6.addConnection(7, .18);
		n6.addConnection(6, .08);	
		
		assertTrue(n1.equals(n2));
		
		assertFalse(n1.equals(n3));
		
		assertFalse(n1.equals(n3a));
		
		assertFalse(n1.equals(n4));
		
		assertFalse(n1.equals(n5));
		
		assertFalse(n1.equals(n6));
	}

	/**
	 * Test method for {@link psograph.graph.Node#getConnectionInfo(psograph.graph.Node)}.
	 * @throws Exception 
	 */
	@Test
	public void testGetConnectionInfoNode() throws Exception 
	{
		Node n1 = new Node(1, .5, .4);
		n1.addConnection(3, .28);
		n1.addConnection(4, .18);
		n1.addConnection(6, .08);
		
		Node n3 = new Node(3,.1, .1);
		Node n4 = new Node(4,.1, .1);
		Node n6 = new Node(6,.1, .1);
		Node n7 = new Node(7,.1, .1);
		
		assertTrue(Double.compare(n1.getConnectionInfo(n3).getWeight(),.28)==0);
		assertTrue(Double.compare(n1.getConnectionInfo(n4).getWeight(),.18)==0);
		assertTrue(Double.compare(n1.getConnectionInfo(n6).getWeight(),.08)==0);
		
		assertTrue(n1.getConnectionInfo(n7) == null);		

	}


	/**
	 * Test method for {@link psograph.graph.Node#isConnectedTo(int)}.
	 * @throws Exception 
	 */
	@SuppressWarnings("deprecation")
	@Test
	public void testIsConnectedToInt() throws Exception 
	{
		Node n1 = new Node(1, .5, .4);
		n1.addConnection(3, .28);
		n1.addConnection(4, .18);
		n1.addConnection(6, .08);
		
		assertTrue(n1.isConnectedTo(3));
		assertTrue(n1.isConnectedTo(4));
		assertTrue(n1.isConnectedTo(6));
		
		assertFalse(n1.isConnectedTo(7));
		
	}

	/**
	 * Test method for {@link psograph.graph.Node#isConnectedTo(psograph.graph.Node)}.
	 * @throws Exception 
	 */
	@Test
	public void testIsConnectedToNode() throws Exception {
		Node n1 = new Node(1, .5, .4);
		n1.addConnection(3, .28);
		n1.addConnection(4, .18);
		n1.addConnection(6, .08);
		
		Node n3 = new Node(3,.1, .1);
		Node n4 = new Node(4,.1, .1);
		Node n6 = new Node(6,.1, .1);
		Node n7 = new Node(7,.1, .1);
		
		assertTrue(n1.isConnectedTo(n3));
		assertTrue(n1.isConnectedTo(n4));
		assertTrue(n1.isConnectedTo(n6));
		
		assertFalse(n1.isConnectedTo(n7));		
	}	
	
	/**
	 * Test method for {@link psograph.graph.Node#removeConnection(int)}.
	 * @throws Exception 
	 */
	@Test
	public void testRemoveConnection() throws Exception 
	{
		Node n = new Node(1, .5, .4);
		
		n.addConnection(2, .28);
		n.addConnection(3, .18);
    	n.addConnection(5, .08);
    	
    	try
    	{
    		n.removeConnection(1);
    	}
    	catch(Exception e)
    	{
    		//remove self
    	}
    	
    	try
    	{
    		n.removeConnection(7);
    	}
    	catch(Exception e)
    	{
    		//remove non-existant
    	}
    	
    	n.removeConnection(2);
    	
    	assertTrue(n.getDegree()== 2);
    	
    	n.removeConnection(3);
    	
    	assertTrue(n.getDegree()== 1);
    	
    	n.removeConnection(5);
    	
    	assertTrue(n.getDegree()== 0);
	
	}	

	/**
	 * Test method for {@link psograph.graph.Node#Node(psograph.graph.Node)}.
	 * @throws Exception 
	 */
	@Test
	public void testNodeNode() throws Exception 
	{
		Node n = new Node(1, .5, .4);
		
		n.addConnection(2, .28);
		n.addConnection(3, .18);
    	n.addConnection(5, .08);
    	
    	Node n2 = new Node(n);
    	
    	
		assertTrue(n.equals(n2));
		
	}

}
