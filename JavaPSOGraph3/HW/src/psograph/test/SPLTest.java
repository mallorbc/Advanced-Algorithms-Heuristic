/**
 * 
 */
package psograph.test;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import psograph.graph.CalculatedGraph;
import psograph.graph.Graph;

import psograph.measurements.SPL;


	/**
	 * @author Patrick
	 *
	 */
	public class SPLTest extends CompareTest
	{
		Graph m_graph;
		CalculatedGraph m_calcGraph;
		/**
		 * @throws java.lang.Exception
		 */
		@BeforeClass
		public static void setUpBeforeClass() throws Exception 
		{
		}
	
		/**
		 * @throws java.lang.Exception
		 */
		@Before
		public void setUp() throws Exception 
		{
		}
	
		/**
		 * Test method for {@link psograph.measurements.SPL#SPL(psograph.graph.Graph)}.
		 */
		@Test
		public void testSPL1() throws Exception
		{
			m_graph = TestUtils.createTestGraph17();
			
			SPL spl = new SPL(m_graph);
			Graph graphs[] = spl.measure();
			
			assertTrue(graphs.length > 0);
			
			m_graph.print();
			
			StringBuilder sb = spl.print();
			
			compareData("psograph/test/SPL1_master.txt", sb);
			
		
		}
	
		/**
		 * Test method for {@link psograph.measurements.SPL#SPL(psograph.graph.Graph)}.
		 */
		@Test
		public void testSPL2() throws Exception
		{
			System.out.println("testSPL2==================");
			m_graph = TestUtils.createSquareGraph();
			
			SPL spl = new SPL(m_graph);
			Graph graphs[] = spl.measure();
			assertTrue(graphs.length > 0);
			
			m_graph.print();
			
			StringBuilder sb = spl.print();
			
			System.out.println("==================");
			System.out.println(sb.toString());
			System.out.println("==================");
			
			
			compareData("psograph/test/SPL2_master.txt", sb);
			System.out.println("End testSPL2==================");
		
		}	

		/**
		 * Test method for {@link psograph.measurements.SPL#SPL(psograph.graph.Graph)}.
		 */
		@Test
		public void testSPL3() throws Exception
		{
			System.out.println("testSPL3==================");
			m_graph = TestUtils.createSquareGraph2();
			
			SPL spl = new SPL(m_graph);
			Graph graphs[] = spl.measure();
			assertTrue(graphs.length > 0);
			
			m_graph.print();
			
			StringBuilder sb = spl.print();
			
			System.out.println("==================");
			System.out.println(sb.toString());
			System.out.println("==================");
			
			
			compareData("psograph/test/SPL3_master.txt", sb);
			System.out.println("End testSPL3==================");
		
		}			
		

		/**
		 * Test method for {@link psograph.measurements.SPL#determineSPLs()}.
		 */
		@Test
		public void testDetermineSPLs() {
			fail("Not yet implemented");
		}
	
		/**
		 * Test method for {@link psograph.measurements.SPL#Measure()}.
		 */
		@Test
		public void testMeasure() {
			fail("Not yet implemented");
		}
	
		/**
		 * Test method for {@link psograph.measurements.SPL#evaulateNodeWithBFS(int, psograph.graph.Graph)}.
		 */
		@Test
		public void testEvaulateNodeWithBFS() {
			fail("Not yet implemented");
		}
	
		/**
		 * Test method for {@link psograph.measurements.SPL#isCalculated()}.
		 */
		@Test
		public void testIsCalculated() {
			fail("Not yet implemented");
		}
	
	}

	
