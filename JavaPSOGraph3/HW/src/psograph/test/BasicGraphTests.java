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

import java.util.*;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import psograph.graph.Graph;
import psograph.graph.Node;



public class BasicGraphTests  {

	Graph m_graph;
	
	@Before
	public void setUp() throws Exception 
	{

		m_graph = TestUtils.createTestGraph();
	}
	
	@Test
	public void testGraph() throws Exception
	{
		Graph copy = new Graph(m_graph);
		
		
		copy.addNode();
	
		//Test orig
		Collection<Node> c = m_graph.getHeaderNodesMap().values();
		Vector<Node> v = new Vector<Node>(c);
		Vector<Integer> vci;
		
		if(v.size() != 10)
		{
			fail("initial test network does not have 10 nodes");
		}
		// Node 0 neighbors 1 3 5 7
		assertTrue("Node 0 does not have a degree of 4",v.get(0).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(0).getNeighbors().keySet() );
		assertTrue("Node 0 does not have a connection to 1",vci.get(0) ==1);
		assertTrue("Node 0 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 0 does not have a connection to 5",vci.get(2) == 5);
		assertTrue("Node 0 does not have a connection to 7",vci.get(3) == 7);
		
		// Node 1 neighbors 0 2 3 4 6
		assertTrue("Node 1 does not have a degree of 5",v.get(1).getDegree() == 5);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(1).getNeighbors().keySet() );
		assertTrue("Node 1 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 1 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 1 does not have a connection to 3",vci.get(2) == 3);
		assertTrue("Node 1 does not have a connection to 4",vci.get(3) == 4);
		assertTrue("Node 1 does not have a connection to 6",vci.get(4) == 6);	
		
		//Node 2 neighbors 1 3 6 9
		assertTrue("Node 2 does not have a degree of 4",v.get(2).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(2).getNeighbors().keySet() );
		assertTrue("Node 2 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 2 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 2 does not have a connection to 6",vci.get(2) == 6);
		assertTrue("Node 2 does not have a connection to 9",vci.get(3) == 9);
		
		//Node 3 neighbors 0 1 2 8
		assertTrue("Node 3 does not have a degree of 4",v.get(3).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(3).getNeighbors().keySet() );
		assertTrue("Node 3 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 3 does not have a connection to 1",vci.get(1) == 1);
		assertTrue("Node 3 does not have a connection to 2",vci.get(2) == 2);
		assertTrue("Node 3 does not have a connection to 8",vci.get(3) == 8);		
		
		//Node 4 neighbors 1 5 9
		assertTrue("Node 4 does not have a degree of 3",v.get(4).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(4).getNeighbors().keySet() );
		assertTrue("Node 4 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 4 does not have a connection to 5",vci.get(1) == 5);
		assertTrue("Node 4 does not have a connection to 9",vci.get(2) == 9);

		//Node 5 neighbors 0 4 8
		assertTrue("Node 5 does not have a degree of 3",v.get(5).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(5).getNeighbors().keySet() );
		assertTrue("Node 5 does not have a connection to 1",vci.get(0) == 0);
		assertTrue("Node 5 does not have a connection to 5",vci.get(1) == 4);
		assertTrue("Node 5 does not have a connection to 9",vci.get(2) == 8);
		
		//Node 6 neighbors 1 2 7
		assertTrue("Node 6 does not have a degree of 3",v.get(6).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(6).getNeighbors().keySet() );
		assertTrue("Node 6 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 6 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 6 does not have a connection to 7",vci.get(2) == 7);
		
		//Node 7 neighbors 0 6
		assertTrue("Node 7 does not have a degree of 3",v.get(7).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(7).getNeighbors().keySet() );
		assertTrue("Node 7 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 7 does not have a connection to 6",vci.get(1) == 6);		
		
		//Node 8 neighbors 3 5
		assertTrue("Node 8 does not have a degree of 3",v.get(8).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(8).getNeighbors().keySet() );
		assertTrue("Node 8 does not have a connection to 3",vci.get(0) == 3);
		assertTrue("Node 8 does not have a connection to 5",vci.get(1) == 5);	
		
		//Node 9 neighbors 2 4
		assertTrue("Node 9 does not have a degree of 3",v.get(9).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(9).getNeighbors().keySet() );
		assertTrue("Node 9 does not have a connection to 2",vci.get(0) == 2);
		assertTrue("Node 9 does not have a connection to 4",vci.get(1) == 4);
		
		
		// now test the copy to make sure it is different
		c = copy.getHeaderNodesMap().values();
		v = new Vector<Node>(c);

		
		if(v.size() != 11)
		{
			fail("initial test network does not have 10 nodes");
		}
		// Node 0 neighbors 1 3 5 7
		assertTrue("Node 0 does not have a degree of 4",v.get(0).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(0).getNeighbors().keySet() );
		assertTrue("Node 0 does not have a connection to 1",vci.get(0) ==1);
		assertTrue("Node 0 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 0 does not have a connection to 5",vci.get(2) == 5);
		assertTrue("Node 0 does not have a connection to 7",vci.get(3) == 7);
		
		// Node 1 neighbors 0 2 3 4 6
		assertTrue("Node 1 does not have a degree of 5",v.get(1).getDegree() == 5);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(1).getNeighbors().keySet() );
		assertTrue("Node 1 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 1 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 1 does not have a connection to 3",vci.get(2) == 3);
		assertTrue("Node 1 does not have a connection to 4",vci.get(3) == 4);
		assertTrue("Node 1 does not have a connection to 6",vci.get(4) == 6);	
		
		//Node 2 neighbors 1 3 6 9
		assertTrue("Node 2 does not have a degree of 4",v.get(2).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(2).getNeighbors().keySet() );
		assertTrue("Node 2 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 2 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 2 does not have a connection to 6",vci.get(2) == 6);
		assertTrue("Node 2 does not have a connection to 9",vci.get(3) == 9);
		
		//Node 3 neighbors 0 1 2 8
		assertTrue("Node 3 does not have a degree of 4",v.get(3).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(3).getNeighbors().keySet() );
		assertTrue("Node 3 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 3 does not have a connection to 1",vci.get(1) == 1);
		assertTrue("Node 3 does not have a connection to 2",vci.get(2) == 2);
		assertTrue("Node 3 does not have a connection to 8",vci.get(3) == 8);		
		
		//Node 4 neighbors 1 5 9
		assertTrue("Node 4 does not have a degree of 3",v.get(4).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(4).getNeighbors().keySet() );
		assertTrue("Node 4 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 4 does not have a connection to 5",vci.get(1) == 5);
		assertTrue("Node 4 does not have a connection to 9",vci.get(2) == 9);

		//Node 5 neighbors 0 4 8
		assertTrue("Node 5 does not have a degree of 3",v.get(5).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(5).getNeighbors().keySet() );
		assertTrue("Node 5 does not have a connection to 1",vci.get(0) == 0);
		assertTrue("Node 5 does not have a connection to 5",vci.get(1) == 4);
		assertTrue("Node 5 does not have a connection to 9",vci.get(2) == 8);
		
		//Node 6 neighbors 1 2 7
		assertTrue("Node 6 does not have a degree of 3",v.get(6).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(6).getNeighbors().keySet() );
		assertTrue("Node 6 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 6 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 6 does not have a connection to 7",vci.get(2) == 7);
		
		//Node 7 neighbors 0 6
		assertTrue("Node 7 does not have a degree of 3",v.get(7).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(7).getNeighbors().keySet() );
		assertTrue("Node 7 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 7 does not have a connection to 6",vci.get(1) == 6);		
		
		//Node 8 neighbors 3 5
		assertTrue("Node 8 does not have a degree of 3",v.get(8).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(8).getNeighbors().keySet() );
		assertTrue("Node 8 does not have a connection to 3",vci.get(0) == 3);
		assertTrue("Node 8 does not have a connection to 5",vci.get(1) == 5);	
		
		//Node 9 neighbors 2 4
		assertTrue("Node 9 does not have a degree of 3",v.get(9).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(9).getNeighbors().keySet() );
		assertTrue("Node 9 does not have a connection to 2",vci.get(0) == 2);
		assertTrue("Node 9 does not have a connection to 4",vci.get(1) == 4);	
		
		//Node 10 neighbors 
		assertTrue("Node 10 does not have a degree of 3",v.get(10).getDegree() == 0);
		assertTrue( v.get(10).getNeighbors() == null );	
		
		
		
	}
	
	
	@Test
	public void testRemoveNode() 
	{
		m_graph.removeNode(0);
		m_graph.removeNode(5);
		
		Collection<Node> c = m_graph.getHeaderNodesMap().values();
		Vector<Node> v = new Vector<Node>(c);
		Vector<Integer> vci;
		
		if(v.size() != 8)
		{
			fail("initial test network does not have 10 nodes");
		}
	
		// Node 1 neighbors  2 3 4 6
		assertTrue("Node 1 does not have a degree of 5",v.get(0).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(0).getNeighbors().keySet() );
		assertTrue("Node 1 does not have a connection to 2",vci.get(0) == 2);
		assertTrue("Node 1 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 1 does not have a connection to 4",vci.get(2) == 4);
		assertTrue("Node 1 does not have a connection to 6",vci.get(3) == 6);	
		
		//Node 2 neighbors 1 3 6 9
		assertTrue("Node 2 does not have a degree of 4",v.get(1).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(1).getNeighbors().keySet() );
		assertTrue("Node 2 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 2 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 2 does not have a connection to 6",vci.get(2) == 6);
		assertTrue("Node 2 does not have a connection to 9",vci.get(3) == 9);
		
		//Node 3 neighbors 1 2 8
		assertTrue("Node 3 does not have a degree of 3",v.get(2).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(2).getNeighbors().keySet() );
		assertTrue("Node 3 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 3 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 3 does not have a connection to 8",vci.get(2) == 8);		
		
		//Node 4 neighbors 1 9
		assertTrue("Node 4 does not have a degree of 2",v.get(3).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(3).getNeighbors().keySet() );
		assertTrue("Node 4 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 4 does not have a connection to 9",vci.get(1) == 9);

	
		//Node 6 neighbors 1 2 7
		assertTrue("Node 6 does not have a degree of 3",v.get(4).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(4).getNeighbors().keySet() );
		assertTrue("Node 6 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 6 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 6 does not have a connection to 7",vci.get(2) == 7);
		
		//Node 7 neighbors 6
		assertTrue("Node 7 does not have a degree of 3",v.get(5).getDegree() == 1);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(5).getNeighbors().keySet() );
		assertTrue("Node 7 does not have a connection to 6",vci.get(0) == 6);		
		
		//Node 8 neighbors 3 
		assertTrue("Node 8 does not have a degree of 3",v.get(6).getDegree() == 1);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(6).getNeighbors().keySet() );
		assertTrue("Node 8 does not have a connection to 3",vci.get(0) == 3);

		
		//Node 9 neighbors 2 4
		assertTrue("Node 9 does not have a degree of 3",v.get(7).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(7).getNeighbors().keySet() );
		assertTrue("Node 9 does not have a connection to 2",vci.get(0) == 2);
		assertTrue("Node 9 does not have a connection to 4",vci.get(1) == 4);		
	}

	@Test
	public void test() 
	{
		Collection<Node> c = m_graph.getHeaderNodesMap().values();
		Vector<Node> v = new Vector<Node>(c);
		Vector<Integer> vci;
		
		if(v.size() != 10)
		{
			fail("initial test network does not have 10 nodes");
		}
		// Node 0 neighbors 1 3 5 7
		assertTrue("Node 0 does not have a degree of 4",v.get(0).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(0).getNeighbors().keySet() );
		assertTrue("Node 0 does not have a connection to 1",vci.get(0) ==1);
		assertTrue("Node 0 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 0 does not have a connection to 5",vci.get(2) == 5);
		assertTrue("Node 0 does not have a connection to 7",vci.get(3) == 7);
		
		// Node 1 neighbors 0 2 3 4 6
		assertTrue("Node 1 does not have a degree of 5",v.get(1).getDegree() == 5);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(1).getNeighbors().keySet() );
		assertTrue("Node 1 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 1 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 1 does not have a connection to 3",vci.get(2) == 3);
		assertTrue("Node 1 does not have a connection to 4",vci.get(3) == 4);
		assertTrue("Node 1 does not have a connection to 6",vci.get(4) == 6);	
		
		//Node 2 neighbors 1 3 6 9
		assertTrue("Node 2 does not have a degree of 4",v.get(2).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(2).getNeighbors().keySet() );
		assertTrue("Node 2 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 2 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 2 does not have a connection to 6",vci.get(2) == 6);
		assertTrue("Node 2 does not have a connection to 9",vci.get(3) == 9);
		
		//Node 3 neighbors 0 1 2 8
		assertTrue("Node 3 does not have a degree of 4",v.get(3).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(3).getNeighbors().keySet() );
		assertTrue("Node 3 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 3 does not have a connection to 1",vci.get(1) == 1);
		assertTrue("Node 3 does not have a connection to 2",vci.get(2) == 2);
		assertTrue("Node 3 does not have a connection to 8",vci.get(3) == 8);		
		
		//Node 4 neighbors 1 5 9
		assertTrue("Node 4 does not have a degree of 3",v.get(4).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(4).getNeighbors().keySet() );
		assertTrue("Node 4 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 4 does not have a connection to 5",vci.get(1) == 5);
		assertTrue("Node 4 does not have a connection to 9",vci.get(2) == 9);

		//Node 5 neighbors 0 4 8
		assertTrue("Node 5 does not have a degree of 3",v.get(5).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(5).getNeighbors().keySet() );
		assertTrue("Node 5 does not have a connection to 1",vci.get(0) == 0);
		assertTrue("Node 5 does not have a connection to 5",vci.get(1) == 4);
		assertTrue("Node 5 does not have a connection to 9",vci.get(2) == 8);
		
		//Node 6 neighbors 1 2 7
		assertTrue("Node 6 does not have a degree of 3",v.get(6).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(6).getNeighbors().keySet() );
		assertTrue("Node 6 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 6 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 6 does not have a connection to 7",vci.get(2) == 7);
		
		//Node 7 neighbors 0 6
		assertTrue("Node 7 does not have a degree of 3",v.get(7).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(7).getNeighbors().keySet() );
		assertTrue("Node 7 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 7 does not have a connection to 6",vci.get(1) == 6);		
		
		//Node 8 neighbors 3 5
		assertTrue("Node 8 does not have a degree of 3",v.get(8).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(8).getNeighbors().keySet() );
		assertTrue("Node 8 does not have a connection to 3",vci.get(0) == 3);
		assertTrue("Node 8 does not have a connection to 5",vci.get(1) == 5);	
		
		//Node 9 neighbors 2 4
		assertTrue("Node 9 does not have a degree of 3",v.get(9).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(9).getNeighbors().keySet() );
		assertTrue("Node 9 does not have a connection to 2",vci.get(0) == 2);
		assertTrue("Node 9 does not have a connection to 4",vci.get(1) == 4);	
		
	}	
	
	@Test
	public void testInputAddRemoveConnections()
	{
		try
		{
			Node n = m_graph.getHeaderNodesMap().get(0);
			n.removeConnection(0);
			fail("should have thrown an exception of removing an non-existant connection on a Node");
			
		}
		catch(Exception e)
		{
			
		}
		
		try
		{
			Node n = m_graph.getHeaderNodesMap().get(0);
			n.removeConnection(9);
			fail("should have thrown an exception of removing an non-existant connection on a Node");
			
		}
		catch(Exception e)
		{
			
		}
		
		try
		{
			Node n = m_graph.getHeaderNodesMap().get(0);
			n.addConnection(0,1);
			fail("should have thrown an exception of adding an already existant connection on a Node");
			
		}
		catch(Exception e)
		{
			
		}	
		
		try
		{
			Node n = m_graph.getHeaderNodesMap().get(0);
			n.addConnection(3,1);
			fail("should have thrown an exception of adding an already existant connection on a Node");
			
		}
		catch(Exception e)
		{
			
		}
		
	}
	
	@Test
	public void testAddConnectionIntInt() throws Exception
	{
		m_graph.addConnection(4,8);
		Collection<Node> c = m_graph.getHeaderNodesMap().values();
		Vector<Node> v = new Vector<Node>(c);
		Vector<Integer> vci;
		
		if(v.size() != 10)
		{
			fail("initial test network does not have 10 nodes");
		}
		// Node 0 neighbors 1 3 5 7
		assertTrue("Node 0 does not have a degree of 4",v.get(0).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(0).getNeighbors().keySet() );
		assertTrue("Node 0 does not have a connection to 1",vci.get(0) ==1);
		assertTrue("Node 0 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 0 does not have a connection to 5",vci.get(2) == 5);
		assertTrue("Node 0 does not have a connection to 7",vci.get(3) == 7);
		
		// Node 1 neighbors 0 2 3 4 6 8
		assertTrue("Node 1 does not have a degree of 5",v.get(1).getDegree() == 5);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(1).getNeighbors().keySet() );
		assertTrue("Node 1 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 1 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 1 does not have a connection to 3",vci.get(2) == 3);
		assertTrue("Node 1 does not have a connection to 4",vci.get(3) == 4);
		assertTrue("Node 1 does not have a connection to 6",vci.get(4) == 6);	
		
		//Node 2 neighbors 1 3 6 9
		assertTrue("Node 2 does not have a degree of 4",v.get(2).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(2).getNeighbors().keySet() );
		assertTrue("Node 2 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 2 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 2 does not have a connection to 6",vci.get(2) == 6);
		assertTrue("Node 2 does not have a connection to 9",vci.get(3) == 9);
		
		//Node 3 neighbors 0 1 2 8
		assertTrue("Node 3 does not have a degree of 4",v.get(3).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(3).getNeighbors().keySet() );
		assertTrue("Node 3 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 3 does not have a connection to 1",vci.get(1) == 1);
		assertTrue("Node 3 does not have a connection to 2",vci.get(2) == 2);
		assertTrue("Node 3 does not have a connection to 8",vci.get(3) == 8);		
		
		//Node 4 neighbors 1 5 8 9
		assertTrue("Node 4 does not have a degree of 3",v.get(4).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(4).getNeighbors().keySet() );
		assertTrue("Node 4 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 4 does not have a connection to 5",vci.get(1) == 5);
		assertTrue("Node 4 does not have a connection to 8",vci.get(2) == 8);
		assertTrue("Node 4 does not have a connection to 9",vci.get(3) == 9);

		//Node 5 neighbors 0 4 8
		assertTrue("Node 5 does not have a degree of 3",v.get(5).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(5).getNeighbors().keySet() );
		assertTrue("Node 5 does not have a connection to 1",vci.get(0) == 0);
		assertTrue("Node 5 does not have a connection to 5",vci.get(1) == 4);
		assertTrue("Node 5 does not have a connection to 9",vci.get(2) == 8);
		
		//Node 6 neighbors 1 2 7
		assertTrue("Node 6 does not have a degree of 3",v.get(6).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(6).getNeighbors().keySet() );
		assertTrue("Node 6 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 6 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 6 does not have a connection to 7",vci.get(2) == 7);
		
		//Node 7 neighbors 0 6
		assertTrue("Node 7 does not have a degree of 3",v.get(7).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(7).getNeighbors().keySet() );
		assertTrue("Node 7 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 7 does not have a connection to 6",vci.get(1) == 6);		
		
		//Node 8 neighbors 3 4 5
		assertTrue("Node 8 does not have a degree of 3",v.get(8).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(8).getNeighbors().keySet() );
		assertTrue("Node 8 does not have a connection to 3",vci.get(0) == 3);
		assertTrue("Node 8 does not have a connection to 4",vci.get(1) == 4);
		assertTrue("Node 8 does not have a connection to 5",vci.get(2) == 5);	
		
		//Node 9 neighbors 2 4
		assertTrue("Node 9 does not have a degree of 3",v.get(9).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(9).getNeighbors().keySet() );
		assertTrue("Node 9 does not have a connection to 2",vci.get(0) == 2);
		assertTrue("Node 9 does not have a connection to 4",vci.get(1) == 4);			
		
	}

	@Test
	public void testAddConnectionIntIntInt()  throws Exception
	{
		m_graph.addConnection(4, 8, 10);
		Collection<Node> c = m_graph.getHeaderNodesMap().values();
		Vector<Node> v = new Vector<Node>(c);
		Vector<Integer> vci;
		
		
		if(v.size() != 10)
		{
			fail("initial test network does not have 10 nodes");
		}
		// Node 0 neighbors 1 3 5 7
		assertTrue("Node 0 does not have a degree of 4",v.get(0).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(0).getNeighbors().keySet() );
		assertTrue("Node 0 does not have a connection to 1",vci.get(0) ==1);
		assertTrue("Node 0 does not have a connection to 1 with weight", v.get(0).getNeighbors().get(vci.get(0)).getWeight() == 1);
		assertTrue("Node 0 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 0 does not have a connection to 3 with weight", v.get(0).getNeighbors().get(vci.get(1)).getWeight() == 2.0);
		assertTrue("Node 0 does not have a connection to 5",vci.get(2) == 5);
		assertTrue("Node 0 does not have a connection to 5 with weight", v.get(0).getNeighbors().get(vci.get(2)).getWeight() == 2);
		assertTrue("Node 0 does not have a connection to 7",vci.get(3) == 7);
		assertTrue("Node 0 does not have a connection to 7 with weight", v.get(0).getNeighbors().get(vci.get(3)).getWeight() == 6);
		
		// Node 1 neighbors 0 2 3 4 6 8
		assertTrue("Node 1 does not have a degree of 5",v.get(1).getDegree() == 5);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(1).getNeighbors().keySet() );
		assertTrue("Node 1 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 1 does not have a connection to 0 with weight", v.get(1).getNeighbors().get(vci.get(0)).getWeight() == 1);
		assertTrue("Node 1 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 1 does not have a connection to 2 with weight", v.get(1).getNeighbors().get(vci.get(1)).getWeight() == 4);
		assertTrue("Node 1 does not have a connection to 3",vci.get(2) == 3);
		assertTrue("Node 1 does not have a connection to 3 with weight", v.get(1).getNeighbors().get(vci.get(2)).getWeight() == 4);
		assertTrue("Node 1 does not have a connection to 4",vci.get(3) == 4);
		assertTrue("Node 1 does not have a connection to 4 with weight", v.get(1).getNeighbors().get(vci.get(3)).getWeight() == 5);
		assertTrue("Node 1 does not have a connection to 6",vci.get(4) == 6);	
		assertTrue("Node 1 does not have a connection to 6 with weight", v.get(1).getNeighbors().get(vci.get(4)).getWeight() == 2);
		
		//Node 2 neighbors 1 3 6 9
		assertTrue("Node 2 does not have a degree of 4",v.get(2).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(2).getNeighbors().keySet() );
		assertTrue("Node 2 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 2 does not have a connection to 1 with weight", v.get(2).getNeighbors().get(vci.get(0)).getWeight() == 4);
		assertTrue("Node 2 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 2 does not have a connection to 3 with weight", v.get(2).getNeighbors().get(vci.get(1)).getWeight() == 3);
		assertTrue("Node 2 does not have a connection to 6",vci.get(2) == 6);
		assertTrue("Node 2 does not have a connection to 6 with weight", v.get(2).getNeighbors().get(vci.get(2)).getWeight() == 2);
		assertTrue("Node 2 does not have a connection to 9",vci.get(3) == 9);
		assertTrue("Node 2 does not have a connection to 9 with weight", v.get(2).getNeighbors().get(vci.get(3)).getWeight() == 4);
		
		//Node 3 neighbors 0 1 2 8
		assertTrue("Node 3 does not have a degree of 4",v.get(3).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(3).getNeighbors().keySet() );
		assertTrue("Node 3 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 3 does not have a connection to 0 with weight", v.get(3).getNeighbors().get(vci.get(0)).getWeight() == 2);
		assertTrue("Node 3 does not have a connection to 1",vci.get(1) == 1);
		assertTrue("Node 3 does not have a connection to 1 with weight", v.get(3).getNeighbors().get(vci.get(1)).getWeight() == 4);
		assertTrue("Node 3 does not have a connection to 2",vci.get(2) == 2);
		assertTrue("Node 3 does not have a connection to 2 with weight", v.get(3).getNeighbors().get(vci.get(2)).getWeight() == 3);
		assertTrue("Node 3 does not have a connection to 8",vci.get(3) == 8);
		assertTrue("Node 3 does not have a connection to 8 with weight", v.get(3).getNeighbors().get(vci.get(3)).getWeight() == 5);
		
		//Node 4 neighbors 1 5 8 9
		assertTrue("Node 4 does not have a degree of 3",v.get(4).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(4).getNeighbors().keySet() );
		assertTrue("Node 4 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 4 does not have a connection to 1 with weight", v.get(4).getNeighbors().get(vci.get(0)).getWeight() == 5);
		assertTrue("Node 4 does not have a connection to 5",vci.get(1) == 5);
		assertTrue("Node 4 does not have a connection to 5 with weight", v.get(4).getNeighbors().get(vci.get(1)).getWeight() == 4);
		assertTrue("Node 4 does not have a connection to 8",vci.get(2) == 8);
		assertTrue("Node 4 does not have a connection to 8 with weight", v.get(4).getNeighbors().get(vci.get(2)).getWeight() == 10);
		assertTrue("Node 4 does not have a connection to 9",vci.get(3) == 9);
		assertTrue("Node 4 does not have a connection to 9 with weight", v.get(4).getNeighbors().get(vci.get(3)).getWeight() == 22);

		//Node 5 neighbors 0 4 8
		assertTrue("Node 5 does not have a degree of 3",v.get(5).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(5).getNeighbors().keySet() );
		assertTrue("Node 5 does not have a connection to 1",vci.get(0) == 0);
		assertTrue("Node 5 does not have a connection to 1 with weight", v.get(5).getNeighbors().get(vci.get(0)).getWeight() == 2);
		assertTrue("Node 5 does not have a connection to 5",vci.get(1) == 4);
		assertTrue("Node 5 does not have a connection to 5 with weight", v.get(5).getNeighbors().get(vci.get(1)).getWeight() == 4);
		assertTrue("Node 5 does not have a connection to 9",vci.get(2) == 8);
		assertTrue("Node 5 does not have a connection to 9 with weight", v.get(5).getNeighbors().get(vci.get(2)).getWeight() == 4);
		
		//Node 6 neighbors 1 2 7
		assertTrue("Node 6 does not have a degree of 3",v.get(6).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(6).getNeighbors().keySet() );
		assertTrue("Node 6 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 6 does not have a connection to 1 with weight", v.get(6).getNeighbors().get(vci.get(0)).getWeight() == 2);
		assertTrue("Node 6 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 6 does not have a connection to 2 with weight", v.get(6).getNeighbors().get(vci.get(1)).getWeight() == 2);
		assertTrue("Node 6 does not have a connection to 7",vci.get(2) == 7);
		assertTrue("Node 6 does not have a connection to 7 with weight", v.get(6).getNeighbors().get(vci.get(2)).getWeight() == 2);
		
		//Node 7 neighbors 0 6
		assertTrue("Node 7 does not have a degree of 3",v.get(7).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(7).getNeighbors().keySet() );
		assertTrue("Node 7 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 7 does not have a connection to 0 with weight", v.get(7).getNeighbors().get(vci.get(0)).getWeight() == 6);
		assertTrue("Node 7 does not have a connection to 6",vci.get(1) == 6);		
		assertTrue("Node 7 does not have a connection to 6 with weight", v.get(7).getNeighbors().get(vci.get(1)).getWeight() == 2);
		
		//Node 8 neighbors 3 4 5
		assertTrue("Node 8 does not have a degree of 3",v.get(8).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(8).getNeighbors().keySet() );
		assertTrue("Node 8 does not have a connection to 3",vci.get(0) == 3);
		assertTrue("Node 8 does not have a connection to 3 with weight", v.get(8).getNeighbors().get(vci.get(0)).getWeight() == 5);
		assertTrue("Node 8 does not have a connection to 4",vci.get(1) == 4);
		assertTrue("Node 8 does not have a connection to 4 with weight", v.get(8).getNeighbors().get(vci.get(1)).getWeight() == 10);
		assertTrue("Node 8 does not have a connection to 5",vci.get(2) == 5);	
		assertTrue("Node 8 does not have a connection to 5 with weight", v.get(8).getNeighbors().get(vci.get(2)).getWeight() == 4);
		
		//Node 9 neighbors 2 4
		assertTrue("Node 9 does not have a degree of 3",v.get(9).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(9).getNeighbors().keySet() );
		assertTrue("Node 9 does not have a connection to 2",vci.get(0) == 2);
		assertTrue("Node 9 does not have a connection to 2 with weight", v.get(9).getNeighbors().get(vci.get(0)).getWeight() == 4);
		assertTrue("Node 9 does not have a connection to 4",vci.get(1) == 4);
		assertTrue("Node 9 does not have a connection to 4 with weight", v.get(9).getNeighbors().get(vci.get(1)).getWeight() == 22);
	}

	@Test
	public void testRemoveConnection() throws Exception
	{
		m_graph.removeConnection(2,9);
		
		Collection<Node> c = m_graph.getHeaderNodesMap().values();
		Vector<Node> v = new Vector<Node>(c);
		Vector<Integer> vci;
		
		if(v.size() != 10)
		{
			fail("initial test network does not have 10 nodes");
		}
		// Node 0 neighbors 1 3 5 7
		assertTrue("Node 0 does not have a degree of 4",v.get(0).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(0).getNeighbors().keySet() );
		assertTrue("Node 0 does not have a connection to 1",vci.get(0) ==1);
		assertTrue("Node 0 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 0 does not have a connection to 5",vci.get(2) == 5);
		assertTrue("Node 0 does not have a connection to 7",vci.get(3) == 7);
		
		// Node 1 neighbors 0 2 3 4 6
		assertTrue("Node 1 does not have a degree of 5",v.get(1).getDegree() == 5);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(1).getNeighbors().keySet() );
		assertTrue("Node 1 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 1 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 1 does not have a connection to 3",vci.get(2) == 3);
		assertTrue("Node 1 does not have a connection to 4",vci.get(3) == 4);
		assertTrue("Node 1 does not have a connection to 6",vci.get(4) == 6);	
		
		//Node 2 neighbors 1 3 6
		assertTrue("Node 2 does not have a degree of 3, but of "+v.get(2).getDegree(),v.get(2).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(2).getNeighbors().keySet() );
		assertTrue("Node 2 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 2 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 2 does not have a connection to 6",vci.get(2) == 6);

		
		//Node 3 neighbors 0 1 2 8
		assertTrue("Node 3 does not have a degree of 4",v.get(3).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(3).getNeighbors().keySet() );
		assertTrue("Node 3 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 3 does not have a connection to 1",vci.get(1) == 1);
		assertTrue("Node 3 does not have a connection to 2",vci.get(2) == 2);
		assertTrue("Node 3 does not have a connection to 8",vci.get(3) == 8);		
		
		//Node 4 neighbors 1 5 9
		assertTrue("Node 4 does not have a degree of 3",v.get(4).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(4).getNeighbors().keySet() );
		assertTrue("Node 4 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 4 does not have a connection to 5",vci.get(1) == 5);
		assertTrue("Node 4 does not have a connection to 9",vci.get(2) == 9);

		//Node 5 neighbors 0 4 8
		assertTrue("Node 5 does not have a degree of 3",v.get(5).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(5).getNeighbors().keySet() );
		assertTrue("Node 5 does not have a connection to 1",vci.get(0) == 0);
		assertTrue("Node 5 does not have a connection to 5",vci.get(1) == 4);
		assertTrue("Node 5 does not have a connection to 9",vci.get(2) == 8);
		
		//Node 6 neighbors 1 2 7
		assertTrue("Node 6 does not have a degree of 3",v.get(6).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(6).getNeighbors().keySet() );
		assertTrue("Node 6 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 6 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 6 does not have a connection to 7",vci.get(2) == 7);
		
		//Node 7 neighbors 0 6
		assertTrue("Node 7 does not have a degree of 3",v.get(7).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(7).getNeighbors().keySet() );
		assertTrue("Node 7 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 7 does not have a connection to 6",vci.get(1) == 6);		
		
		//Node 8 neighbors 3 5
		assertTrue("Node 8 does not have a degree of 3",v.get(8).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(8).getNeighbors().keySet() );
		assertTrue("Node 8 does not have a connection to 3",vci.get(0) == 3);
		assertTrue("Node 8 does not have a connection to 5",vci.get(1) == 5);	
		
		//Node 9 neighbors 4
		assertTrue("Node 9 does not have a degree of 3",v.get(9).getDegree() == 1);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(9).getNeighbors().keySet() );
		assertTrue("Node 9 does not have a connection to 4",vci.get(0) == 4);			
		
		
	}



	@Test
	public void testAddNode() throws Exception
	{	
		m_graph.addNode();
		
		Collection<Node> c = m_graph.getHeaderNodesMap().values();
		Vector<Node> v = new Vector<Node>(c);
		Vector<Integer> vci;
		
		if(v.size() != 11)
		{
			fail("initial test network does not have 10 nodes");
		}
		// Node 0 neighbors 1 3 5 7
		assertTrue("Node 0 does not have a degree of 4",v.get(0).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(0).getNeighbors().keySet() );
		assertTrue("Node 0 does not have a connection to 1",vci.get(0) ==1);
		assertTrue("Node 0 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 0 does not have a connection to 5",vci.get(2) == 5);
		assertTrue("Node 0 does not have a connection to 7",vci.get(3) == 7);
		
		// Node 1 neighbors 0 2 3 4 6
		assertTrue("Node 1 does not have a degree of 5",v.get(1).getDegree() == 5);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(1).getNeighbors().keySet() );
		assertTrue("Node 1 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 1 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 1 does not have a connection to 3",vci.get(2) == 3);
		assertTrue("Node 1 does not have a connection to 4",vci.get(3) == 4);
		assertTrue("Node 1 does not have a connection to 6",vci.get(4) == 6);	
		
		//Node 2 neighbors 1 3 6 9
		assertTrue("Node 2 does not have a degree of 4",v.get(2).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(2).getNeighbors().keySet() );
		assertTrue("Node 2 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 2 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 2 does not have a connection to 6",vci.get(2) == 6);
		assertTrue("Node 2 does not have a connection to 9",vci.get(3) == 9);
		
		//Node 3 neighbors 0 1 2 8
		assertTrue("Node 3 does not have a degree of 4",v.get(3).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(3).getNeighbors().keySet() );
		assertTrue("Node 3 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 3 does not have a connection to 1",vci.get(1) == 1);
		assertTrue("Node 3 does not have a connection to 2",vci.get(2) == 2);
		assertTrue("Node 3 does not have a connection to 8",vci.get(3) == 8);		
		
		//Node 4 neighbors 1 5 9
		assertTrue("Node 4 does not have a degree of 3",v.get(4).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(4).getNeighbors().keySet() );
		assertTrue("Node 4 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 4 does not have a connection to 5",vci.get(1) == 5);
		assertTrue("Node 4 does not have a connection to 9",vci.get(2) == 9);

		//Node 5 neighbors 0 4 8
		assertTrue("Node 5 does not have a degree of 3",v.get(5).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(5).getNeighbors().keySet() );
		assertTrue("Node 5 does not have a connection to 1",vci.get(0) == 0);
		assertTrue("Node 5 does not have a connection to 5",vci.get(1) == 4);
		assertTrue("Node 5 does not have a connection to 9",vci.get(2) == 8);
		
		//Node 6 neighbors 1 2 7
		assertTrue("Node 6 does not have a degree of 3",v.get(6).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(6).getNeighbors().keySet() );
		assertTrue("Node 6 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 6 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 6 does not have a connection to 7",vci.get(2) == 7);
		
		//Node 7 neighbors 0 6
		assertTrue("Node 7 does not have a degree of 3",v.get(7).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(7).getNeighbors().keySet() );
		assertTrue("Node 7 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 7 does not have a connection to 6",vci.get(1) == 6);		
		
		//Node 8 neighbors 3 5
		assertTrue("Node 8 does not have a degree of 3",v.get(8).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(8).getNeighbors().keySet() );
		assertTrue("Node 8 does not have a connection to 3",vci.get(0) == 3);
		assertTrue("Node 8 does not have a connection to 5",vci.get(1) == 5);	
		
		//Node 9 neighbors 2 4
		assertTrue("Node 9 does not have a degree of 3",v.get(9).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(9).getNeighbors().keySet() );
		assertTrue("Node 9 does not have a connection to 2",vci.get(0) == 2);
		assertTrue("Node 9 does not have a connection to 4",vci.get(1) == 4);	
		
		//Node 10 neighbors 
		assertTrue("Node 10 does not have a degree of 3",v.get(10).getDegree() == 0);
		assertTrue( v.get(10).getNeighbors() == null );

	}

	@Test
	public void testAddNodeIntArray() throws Exception
	{
		int neighbors[];
		neighbors = new int[3];
		neighbors[0] = 0;
		neighbors[1] = 2;
		neighbors[2] = 4;
		m_graph.addNode(neighbors);
		
		Collection<Node> c = m_graph.getHeaderNodesMap().values();
		Vector<Node> v = new Vector<Node>(c);
		Vector<Integer> vci;
		
		m_graph.print();
		
		if(v.size() != 11)
		{
			fail("initial test network does not have 11 nodes");
		}
		// Node 0 neighbors 1 3 5 7 10
		assertTrue("Node 0 does not have a degree of 5",v.get(0).getDegree() == 5);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(0).getNeighbors().keySet() );
		assertTrue("Node 0 does not have a connection to 1",vci.get(0) ==1);
		assertTrue("Node 0 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 0 does not have a connection to 5",vci.get(2) == 5);
		assertTrue("Node 0 does not have a connection to 7",vci.get(3) == 7);
		assertTrue("Node 0 does not have a connection to 10",vci.get(4) == 10);
		
		// Node 1 neighbors 0 2 3 4 6
		assertTrue("Node 1 does not have a degree of 5",v.get(1).getDegree() == 5);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(1).getNeighbors().keySet() );
		assertTrue("Node 1 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 1 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 1 does not have a connection to 3",vci.get(2) == 3);
		assertTrue("Node 1 does not have a connection to 4",vci.get(3) == 4);
		assertTrue("Node 1 does not have a connection to 6",vci.get(4) == 6);	
		
		//Node 2 neighbors 1 3 6 9 10
		assertTrue("Node 2 does not have a degree of 4",v.get(2).getDegree() == 5);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(2).getNeighbors().keySet() );
		assertTrue("Node 2 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 2 does not have a connection to 3",vci.get(1) == 3);
		assertTrue("Node 2 does not have a connection to 6",vci.get(2) == 6);
		assertTrue("Node 2 does not have a connection to 9",vci.get(3) == 9);
		assertTrue("Node 2 does not have a connection to 10",vci.get(4) == 10);
		
		//Node 3 neighbors 0 1 2 8
		assertTrue("Node 3 does not have a degree of 4",v.get(3).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(3).getNeighbors().keySet() );
		assertTrue("Node 3 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 3 does not have a connection to 1",vci.get(1) == 1);
		assertTrue("Node 3 does not have a connection to 2",vci.get(2) == 2);
		assertTrue("Node 3 does not have a connection to 8",vci.get(3) == 8);		
		
		//Node 4 neighbors 1 5 9 10
		assertTrue("Node 4 does not have a degree of 3",v.get(4).getDegree() == 4);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(4).getNeighbors().keySet());
		assertTrue("Node 4 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 4 does not have a connection to 5",vci.get(1) == 5);
		assertTrue("Node 4 does not have a connection to 9",vci.get(2) == 9);
		assertTrue("Node 4 does not have a connection to 10",vci.get(3) == 10);

		
		//Node 5 neighbors 0 4 8
		assertTrue("Node 5 does not have a degree of 3",v.get(5).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(5).getNeighbors().keySet() );
		assertTrue("Node 5 does not have a connection to 1",vci.get(0) == 0);
		assertTrue("Node 5 does not have a connection to 5",vci.get(1) == 4);
		assertTrue("Node 5 does not have a connection to 9",vci.get(2) == 8);
		
		//Node 6 neighbors 1 2 7
		assertTrue("Node 6 does not have a degree of 3",v.get(6).getDegree() == 3);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(6).getNeighbors().keySet() );
		assertTrue("Node 6 does not have a connection to 1",vci.get(0) == 1);
		assertTrue("Node 6 does not have a connection to 2",vci.get(1) == 2);
		assertTrue("Node 6 does not have a connection to 7",vci.get(2) == 7);
		
		//Node 7 neighbors 0 6
		assertTrue("Node 7 does not have a degree of 3",v.get(7).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(7).getNeighbors().keySet() );
		assertTrue("Node 7 does not have a connection to 0",vci.get(0) == 0);
		assertTrue("Node 7 does not have a connection to 6",vci.get(1) == 6);		
		
		//Node 8 neighbors 3 5
		assertTrue("Node 8 does not have a degree of 3",v.get(8).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(8).getNeighbors().keySet() );
		assertTrue("Node 8 does not have a connection to 3",vci.get(0) == 3);
		assertTrue("Node 8 does not have a connection to 5",vci.get(1) == 5);	
		
		//Node 9 neighbors 2 4
		assertTrue("Node 9 does not have a degree of 3",v.get(9).getDegree() == 2);
		vci = new Vector<Integer>( (Collection<Integer>) v.get(9).getNeighbors().keySet() );
		assertTrue("Node 9 does not have a connection to 2",vci.get(0) == 2);
		assertTrue("Node 9 does not have a connection to 4",vci.get(1) == 4);	
	}



}
