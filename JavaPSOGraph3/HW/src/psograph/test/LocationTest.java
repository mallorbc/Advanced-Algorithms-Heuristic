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

import psograph.graph.Location;

/**
 * TestSuite for Location
 * @author Patrick
 *
 */
public class LocationTest  
{

	Location m_location1;
	Location m_location2;
	Location m_location3;
	Location m_location4;
	
	
	@Before
	public void setUp() throws Exception 
	{
		m_location1 = new Location (.132432, .643456354);
		m_location2 = new Location (.132432, .1325);
		m_location3 = new Location (.4345, .643456354);
		m_location4 = new Location (.132432, .643456354);
	}
	
	/**
	 * Test method for {@link psograph.graph.Location#Location(double, double)}.
	 */
	@Test
	public void testLocation() throws Exception
	{
		try 
		{
			@SuppressWarnings("unused")
			Location badLocation = new Location (-1, .2);
			fail("No exception thrown.");
		}
		catch (Exception e)
		{
			//Good we caught an exception
		}
		
		try 
		{
			@SuppressWarnings("unused")
			Location badLocation = new Location (.5, 1.2);
			fail("No exception thrown.");
		}
		catch (Exception e)
		{
			//Good we caught an exception
		}	
		
		@SuppressWarnings("unused")
		Location goodLocation1 = new Location (0.0, 0.0);

		@SuppressWarnings("unused")
		Location goodLocation2 = new Location (1.0, 1.0);

		
		
	}

	/**
	 * Test method for {@link psograph.graph.Location#equals(java.lang.Object)}.
	 */
	@Test
	public void testEqualsObject() 
	{
		assert(m_location1 != m_location2);
		assert(m_location1 != m_location3);
		assert(m_location1 == m_location4);		
	}

	/**
	 * Test method for {@link psograph.graph.Location#getX()}.
	 */
	@Test
	public void testGetX() 
	{
		if(!(Double.compare(m_location1.getX(), .132432) == 0))
		{
			fail("m_location1.getX() ~= .132432");
		}
		
		if(!(Double.compare(m_location2.getX(), .132432) == 0))
		{
			fail("m_location1.getX() ~= .132432");
		}
		
		if(!(Double.compare(m_location3.getX(), .4345) == 0))
		{
			fail("m_location1.getX() ~= .4345");
		}
		
		if(!(Double.compare(m_location4.getX(), .132432) == 0))
		{
			fail("m_location1.getX() ~= .132432");
		}		
	}

	/**
	 * Test method for {@link psograph.graph.Location#getY()}.
	 */
	@Test
	public void testGetY() 
	{
		if(!(Double.compare(m_location1.getY(), .643456354) == 0))
		{
			fail("m_location1.getX() ~= .643456354");
		}
		
		if(!(Double.compare(m_location2.getY(), .1325) == 0))
		{
			fail("m_location1.getX() ~= .1325");
		}
		
		if(!(Double.compare(m_location3.getY(), .643456354) == 0))
		{
			fail("m_location1.getX() ~= .643456354");
		}
		
		if(!(Double.compare(m_location4.getY(), .643456354) == 0))
		{
			fail("m_location1.getX() ~= .643456354");
		}		
	}

}
