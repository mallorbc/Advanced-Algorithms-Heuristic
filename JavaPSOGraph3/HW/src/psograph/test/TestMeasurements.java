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

import java.util.Collection;
import java.util.TreeMap;
import java.util.Vector;

import org.junit.Before;
import org.junit.Test;

import psograph.graph.*;



public class TestMeasurements  {

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
	public void testSumWeights() 
	{
		
		m_calcGraph.print();
		assertTrue("Sum of weights is not equaul to 104, but is equal to "+m_calcGraph.SumAllWeights(),m_calcGraph.SumAllWeights() == 104);
	}


	@Test
	public void getLengthDistribution() throws Exception
	{

			m_graph = TestUtils.createTestGraph3();
			m_calcGraph = new CalculatedGraph(m_graph);
			m_calcGraph.setCostBasis(5);
			m_calcGraph.UpdateCalcuations();
			
			TreeMap<Double, Integer> lengthDistribution = m_calcGraph.getDistributionLinkLengths();  
			
	
			Vector<Double> vCost = new Vector<Double>(lengthDistribution.keySet());		
	
			if(vCost.size() <= 0)
			{
				fail("The cost distribution is empty");
			}
	
			for(int i=0; i < vCost.size(); i++)
			{
				double value = vCost.get(i);
				
				if (Double.compare(0.030000000000000027,value) == 0)
				{
					assertTrue("weight of "+value+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.05385164807134499 ,value) == 0)
				{
					assertTrue("weight of "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.0670820393249937,vCost.get(i)) == 0)
				{
					assertTrue("weight of "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare( 0.068,vCost.get(i)) == 0)
				{
					assertTrue("weight of "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.22627416997969516,vCost.get(i)) == 0)
				{
					assertTrue("weight of "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}

				else if(Double.compare(0.24738633753705966,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.2493270943960965,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}	
				else if(Double.compare(0.2549509756796392,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}					

				


				
				else if(Double.compare(0.25553864678361277,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.28733952042836014,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}	
				else if(Double.compare(0.296951174437819,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.30265491900843117,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.3114482300479488,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				







				
				else if(Double.compare(0.32557641192199416,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==2);	
				}
				else if(Double.compare(0.33941125496954283,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(.4,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.45354161881794264,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.5,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.5755866572463264,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.579827560572969,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.6403124237432849,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}				
				else if(Double.compare(0.6826419266350405,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.6923149572268391,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.7022819946431774,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.7280109889280518,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(1.145512985522207,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(0.35 ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}				
				else if(Double.compare(0.3687817782917155 ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}				
				else if(Double.compare( 0.4440720662234903 ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}				
				else if(Double.compare( 0.4738143096192854  ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}		
				else if(Double.compare( 0.4904079934095691  ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}			
				else if(Double.compare( 0.4964876634922563   ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}	
				else if(Double.compare( 0.49979995998399207   ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}	
				else if(Double.compare( 0.5261178575186363   ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}		  
				else if(Double.compare(0.5586591089385369    ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}	 
				else if(Double.compare(0.5728001396647875     ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}			  
				else if(Double.compare( 0.5953990258641679     ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}	
				else if(Double.compare(  0.6128001305482891     ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}	
				else if(Double.compare(  0.6664082832618455    ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}		 
				else if(Double.compare(  0.6870953354520752    ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else if(Double.compare(  0.688839603971781    ,vCost.get(i)) == 0)
				{
					assertTrue("weight of  "+vCost.get(i)+
							"there are "+ lengthDistribution.get(vCost.get(i)),
							lengthDistribution.get(vCost.get(i))==1);	
				}
				else
				{
					fail("A weight not anticipated for"+vCost.get(i));
					//System.out.println(" cost "+ vCost.get(i) + " Number of with this weight "+
					//		lengthDistribution.get(vCost.get(i)));
					
				}
				//System.out.println(" cost "+ vCost.get(i) + " Number of with this weight "+
				//		lengthDistribution.get(vCost.get(i)));

			}
	}
	
	@Test
	public void testCostOfDirectConnectivity() throws Exception
	{
		m_graph = TestUtils.createSquareGraph();
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		double costOfDirectConnectivity  = m_calcGraph.getCostOfDirectConnectivity();

		assertTrue("CostOfDirectConnectivity "+ costOfDirectConnectivity,Double.compare(costOfDirectConnectivity, 0.16000000000000006)==0);			
		
		
		
	}
	
	@Test
	public void getCostDistribution() throws Exception
	{

			m_graph = TestUtils.createTestGraph3();
			m_calcGraph = new CalculatedGraph(m_graph);
			m_calcGraph.setCostBasis(5);
			m_calcGraph.UpdateCalcuations();
			
			TreeMap<Double, Integer> weightDistribution = m_calcGraph.getDistributionLinkCosts();  
			
	
			Vector<Double> vWeight = new Vector<Double>(weightDistribution.keySet());		
	
			if(vWeight.size() <= 0)
			{
				fail("The weight distribution is empty");
			}
	
			for(int i=0; i < vWeight.size(); i++)
			{
				
				if(Double.compare(0.10239999999999995,vWeight.get(i)) == 0)
				{
					assertTrue("weight of "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.12240000000000002,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.13060000000000002,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.21200000000000005,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==2);	
				}
				else if(Double.compare(0.23040000000000002,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.32000000000000006,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.41139999999999993,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.5,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.6626000000000001,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.6724000000000001,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.8200000000000001,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}				
				else if(Double.compare(0.9320000000000002,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.9586000000000001,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.9863999999999997,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(1.06,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(2.6244,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.0018000000000000032 ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.005799999999999988  ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.009000000000000001   ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.009248000000000001   ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}	
				else if(Double.compare(0.12432800000000002   ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}		
				else if(Double.compare(0.12999999999999995   ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}		 
				else if(Double.compare(0.165128    ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}		 
				else if(Double.compare(0.17636     ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}				  
				else if(Double.compare(0.18320000000000006   ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}				  
				else if(Double.compare(0.1940000000000001    ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}			 
				else if(Double.compare(0.24499999999999997     ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}	
				else if(Double.compare(0.272   ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}		 
				else if(Double.compare(0.3944    ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}	 
				else if(Double.compare(0.44900000000000007   ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}	 
				else if(Double.compare(0.4809999999999999   ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}	 			  
				else if(Double.compare(0.4929999999999999   ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}	 			 
				else if(Double.compare(0.49960000000000016    ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}	 
				else if(Double.compare(0.5536000000000001    ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}	 	 
				else if(Double.compare(0.6242    ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(0.6562000000000001    ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				} 
				else if(Double.compare(0.7090000000000002    ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				} 
				else if(Double.compare(0.7510480000000003     ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}  
				else if(Double.compare(0.8882000000000003    ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}  
				else if(Double.compare(0.9441999999999996     ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}  
				else if(Double.compare(0.9490000000000002      ,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}   
				else
				{
					fail("A weight not anticipated for"+vWeight.get(i));
					//System.out.println(" weight "+ vWeight.get(i) + " Number of with this weight "+
						//	weightDistribution.get(vWeight.get(i)));
				}
				//System.out.println(" weight "+ vWeight.get(i) + " Number of with this weight "+
				//		weightDistribution.get(vWeight.get(i)));
			}

	
	}
	
	@Test
	public void getWeightDistribution() 
	{
		try
		{
			TreeMap<Double, Integer> weightDistribution = m_calcGraph.getWeightDistribution();  
			
	
			Vector<Double> vWeight = new Vector<Double>(weightDistribution.keySet());		
	
			if(vWeight.size() <= 0)
			{
				fail("The weight distribution is empty");
			}
	
			for(int i=0; i < vWeight.size(); i++)
			{
				
				if(Double.compare(1.0,vWeight.get(i)) == 0)
				{
					assertTrue("weight of "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==3);	
				}
				else if(Double.compare(2.0,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==3);	
				}
				else if(Double.compare(3.0,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else if(Double.compare(4.0,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==2);	
				}
				else if(Double.compare(5.0,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==2);	
				}
				else if(Double.compare(22.0,vWeight.get(i)) == 0)
				{
					assertTrue("weight of  "+vWeight.get(i)+
							"there are "+ weightDistribution.get(vWeight.get(i)),
							weightDistribution.get(vWeight.get(i))==1);	
				}
				else
				{
					fail("A weight not anticipated for"+vWeight.get(i));
				}
				//System.out.println(" weight "+ vWeight.get(i) + "Number of with this weight "+
				//		weightDistribution.get(vWeight.get(i)));
			}
		}
		catch(Exception e)
		{
			fail("with an exception");
		}
	
	}
	
	
	@Test
	public void testGetDegreeDistributionAdv()
	{
		
		
		TreeMap<Integer,Vector<Integer>> degreeDist = m_graph.getDegreeDistributionAdvanced();
		
	/*	
	Bucket - 0
	     13
	Bucket - 1
	     7 8 9 11 12
	Bucket - 2
	     10
	Bucket - 3
	     3 4 6
	Bucket - 4
	     1 2
		*/
		
		Collection<Integer> c = degreeDist.keySet();
		Vector<Integer> v = new Vector<Integer>(c);
		
		assertTrue("Node 13 does not degree of 0" , degreeDist.get(v.get(0)).get(0) == 13);
		assertTrue("Bucket of degree 0 should only have 1 node", degreeDist.get(v.get(0)).size() ==1  );
		
		
		assertTrue("Node 7 does not degree of 1" , degreeDist.get(v.get(1)).get(0) == 7);
		assertTrue("Node 8 does not degree of 1" , degreeDist.get(v.get(1)).get(1) == 8);
		assertTrue("Node 9 does not degree of 1" , degreeDist.get(v.get(1)).get(2) == 9);
		assertTrue("Node 11 does not degree of 1" , degreeDist.get(v.get(1)).get(3) == 11);
		assertTrue("Node 12 does not degree of 1" , degreeDist.get(v.get(1)).get(4) == 12);
		
		assertTrue("Node 10 does not degree of 2" , degreeDist.get(v.get(2)).get(0) == 10);
		
		assertTrue("Node 3 does not degree of 3" , degreeDist.get(v.get(3)).get(0) == 3);
		assertTrue("Node 4 does not degree of 3" , degreeDist.get(v.get(3)).get(1) == 4);
		assertTrue("Node 6 does not degree of 3" , degreeDist.get(v.get(3)).get(2) == 6);
		
		assertTrue("Node 1 does not degree of 4" , degreeDist.get(v.get(4)).get(0) == 1);
		assertTrue("Node 2 does not degree of 4" , degreeDist.get(v.get(4)).get(1) == 2);
		

	}

	@Test
	public void testMeasurePercentageInLargestCluster() throws Exception
	{

		
		double measure = m_calcGraph.getLCC();
		assertTrue("Unit test LCC of network should be .75" , Double.compare(measure, 0.75) ==0);
	}

	@Test
	public void testDiameter() throws Exception
	{

	
		double Diameter  = m_calcGraph.getDiameter();
		assertTrue("Unit test diamter of LCC of network should be 4" , Double.compare(Diameter, 4.0) ==0);
	}
	
	@Test
	public void testAISPL()  throws Exception
	{
		
		m_graph = TestUtils.createTestGraph2();
		
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(5);
		m_calcGraph.UpdateCalcuations();
		
		double value = m_calcGraph.getAISPL();
		
		assertTrue("Unit test AISPL of network should be 0.3371212121212121, not " + value ,
				Double.compare(value,0.3371212121212121)==0);
		
		TreeMap<Double, Integer> ISPLDist = m_calcGraph.getDistributionISPL();

		
		Vector<Double> vISPL = new Vector<Double>(ISPLDist.keySet());		

		if(vISPL.size() <= 0)
		{
			fail("The ISPL distribution is empty");
		}

		for(int i=0; i < vISPL.size(); i++)
		{
			
			if(Double.compare(vISPL.get(i),0.0) == 0.0)
			{
				assertTrue("SPL of "+vISPL.get(i)+
						"there are "+ ISPLDist.get(vISPL.get(i)),
						ISPLDist.get(vISPL.get(i))==58);	
			}
			else if(Double.compare(vISPL.get(i),0.25) == 0.0)
			{
				assertTrue("SPL of  "+vISPL.get(i)+
						"there are "+ ISPLDist.get(vISPL.get(i)),
						ISPLDist.get(vISPL.get(i))==6);	
			}
			else if(Double.compare(vISPL.get(i),0.3333333333333333) == 0.0)
			{
				assertTrue("SPL of  "+vISPL.get(i)+
						"there are "+ ISPLDist.get(vISPL.get(i)),
						ISPLDist.get(vISPL.get(i))==18);	
			}
			else if(Double.compare(vISPL.get(i),0.5) == 0.0)
			{
				assertTrue("SPL of  "+vISPL.get(i)+
						"there are "+ ISPLDist.get(vISPL.get(i)),
						ISPLDist.get(vISPL.get(i))==26);	
			}
			else if(Double.compare(vISPL.get(i),1.0) == 0.0)
			{
				assertTrue("SPL of  "+vISPL.get(i)+
						"there are "+ ISPLDist.get(vISPL.get(i)),
						ISPLDist.get(vISPL.get(i))==24);	
			}
			else
			{
				fail("A SPL not anticipated for"+vISPL.get(i));
			}
			//System.out.println("ISPL weight "+ vISPL.get(i) + "Number of SPL with this weight "+
			//		ISPLDist.get(vISPL.get(i)));
		}
		
	}

	@Test
	public void testAISPL2()  throws Exception
	{
		
		m_graph = TestUtils.createSimpleRingGraph();
			
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(10);

		m_calcGraph.UpdateCalcuations();

		
		double value = m_calcGraph.getAISPL();
		assertTrue("Unit test AISPL of network should be 0.75, not " + value ,
				Double.compare(value,0.75)==0);
		
		TreeMap<Double, Integer> ISPLDist = m_calcGraph.getDistributionISPL();

		
		Vector<Double> vISPL = new Vector<Double>(ISPLDist.keySet());		

		if(vISPL.size() <= 0)
		{
			fail("The ISPL distribution is empty");
		}

		for(int i=0; i < vISPL.size(); i++)
		{
			
			if(Double.compare(vISPL.get(i),0.5) == 0.0)
			{
				assertTrue("SPL of  "+vISPL.get(i)+
						"there are "+ ISPLDist.get(vISPL.get(i)),
						ISPLDist.get(vISPL.get(i))==10);	
			}
			else if(Double.compare(vISPL.get(i),1.0) == 0.0)
			{
				assertTrue("SPL of  "+vISPL.get(i)+
						"there are "+ ISPLDist.get(vISPL.get(i)),
						ISPLDist.get(vISPL.get(i))==10);	
			}
			else
			{
				fail("A SPL not anticipated for"+vISPL.get(i));
			}
			//System.out.println("ISPL weight "+ vISPL.get(i) + "Number of SPL with this weight "+
			//		ISPLDist.get(vISPL.get(i)));
		}
		
	}
	
	
	@Test
	public void testAISPL3()  throws Exception
	{
		
		m_graph = TestUtils.createSimpleRingGraphButBroken();
			
		m_calcGraph = new CalculatedGraph(m_graph);
		m_calcGraph.setCostBasis(10);
		m_calcGraph.UpdateCalcuations();
		
		
		double value = m_calcGraph.getAISPL();
		assertTrue("Unit test AISPL of network should be 0.43333333333333335, not " + value ,
				Double.compare(value,0.43333333333333335)==0);
		
		TreeMap<Double, Integer> ISPLDist = m_calcGraph.getDistributionISPL();

		
		Vector<Double> vISPL = new Vector<Double>(ISPLDist.keySet());		

		if(vISPL.size() <= 0)
		{
			fail("The ISPL distribution is empty");
		}

		for(int i=0; i < vISPL.size(); i++)
		{
			
			if(Double.compare(vISPL.get(i),0.0) == 0.0)
			{
				assertTrue("SPL of "+vISPL.get(i)+
						"there are "+ ISPLDist.get(vISPL.get(i)),
						ISPLDist.get(vISPL.get(i))==8);	
			}
			else if(Double.compare(vISPL.get(i),0.25) == 0.0)
			{
				assertTrue("SPL of  "+vISPL.get(i)+
						"there are "+ ISPLDist.get(vISPL.get(i)),
						ISPLDist.get(vISPL.get(i))==3);	
			}
			else if(Double.compare(vISPL.get(i),0.3333333333333333) == 0.0)
			{
				assertTrue("SPL of  "+vISPL.get(i)+
						"there are "+ ISPLDist.get(vISPL.get(i)),
						ISPLDist.get(vISPL.get(i))==2);	
			}
			else if(Double.compare(vISPL.get(i),0.5) == 0.0)
			{
				assertTrue("SPL of  "+vISPL.get(i)+
						"there are "+ ISPLDist.get(vISPL.get(i)),
						ISPLDist.get(vISPL.get(i))==4);	
			}
			else if(Double.compare(vISPL.get(i),1.0) == 0.0)
			{
				assertTrue("SPL of  "+vISPL.get(i)+
						"there are "+ ISPLDist.get(vISPL.get(i)),
						ISPLDist.get(vISPL.get(i))==6);	
			}
			else
			{
				fail("A SPL not anticipated for"+vISPL.get(i));
			}
		//	System.out.println("ISPL weight "+ vISPL.get(i) + "Number of SPL with this weight "+
		//			ISPLDist.get(vISPL.get(i)));
		}
		
	}
	
	@Test
	public void testASPLDisconnected()  throws Exception
	{
		
		boolean isFullyConnected = m_calcGraph.isFullyConnected();
		
		assertTrue("Is graph fully connected should be false, not " + isFullyConnected,
				isFullyConnected == false);		
		

		double value = m_calcGraph.getASPL();
		
		assertTrue("ASPL should be -1 not " + value, Double.compare(value, -1) ==0);

		TreeMap<Integer, Integer> SPLDist = m_calcGraph.getDistributionSPL();
		
		assertTrue("SPL distribution should be NULL " + value, SPLDist == null);
		
		
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
		}
	}


	
	
	@Test
	public void testASPL2()  throws Exception
	{
		m_graph = TestUtils.createSimpleLineGraph();
		boolean isFullyConnected = m_graph.isFullyConnected();
		
		assertTrue("Is graph fully connected should be true, not " + isFullyConnected,
				isFullyConnected == true);	
		
		CalculatedGraph calcGraph = new CalculatedGraph(m_graph);
		calcGraph.setCostBasis(5);
		calcGraph.UpdateCalcuations();


		double value = calcGraph.getASPL();
				
		assertTrue("Unit test ASPL of network should be 2.0, not " + value,
				Double.compare(value, 2.0) == 0);
		
		TreeMap<Integer, Integer> SPLDist = calcGraph.getDistributionSPL();

		
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
						SPLDist.get(vSPL.get(i))==8);	
			}
			else if(vSPL.get(i) == 2)
			{
				assertTrue("SPL of  "+vSPL.get(i)+
						"there are "+ SPLDist.get(vSPL.get(i)),
						SPLDist.get(vSPL.get(i))==6);	
			}
			else if(vSPL.get(i) == 3)
			{
				assertTrue("SPL of  "+vSPL.get(i)+
						"there are "+ SPLDist.get(vSPL.get(i)),
						SPLDist.get(vSPL.get(i))==4);	
			}
			else if(vSPL.get(i) == 4)
			{
				assertTrue("SPL of  "+vSPL.get(i)+
						"there are "+ SPLDist.get(vSPL.get(i)),
						SPLDist.get(vSPL.get(i))==2);	
			}
			else
			{
				fail("A SPL not anticipated for"+vSPL.get(i));
			}
			//System.out.println("SPL weight "+ vSPL.get(i) + "Number of SPL with this weight "+
			//		SPLDist.get(vSPL.get(i)));
		}
	}

	
	
	
}
