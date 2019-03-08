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
package psograph.measurements;

import java.io.Serializable;
import java.util.Random;
import java.util.TreeMap;
import java.util.Vector;

import psograph.graph.Graph;
import psograph.graph.Node;


public class DirectedDamage implements IGraphMeasument, Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 3158315854166756480L;

	
	public DirectedDamage(Graph graph, double per)
	{
		m_isCalculated = false;
		m_graph = graph;
		m_directLCC=0;
		m_directDiameter =0;
		percentageToDamage = per;
	}
	public double getDirectDiamter()
	{
		return m_directDiameter;
	}
	
	public double Measure() throws Exception {

		
		double value = measureDamageDirect(percentageToDamage);
		m_isCalculated = true;
		
		return value;
	}

	public boolean isCalculated()
	{
		return m_isCalculated;
	}
	
	private void damageDirected(double percentageToDamage) throws Exception
	{
		
		//System.out.println("-----------------------------------------------");
		//print();
		
		if(m_graph.getHeaderNodesMap().size() == 0)
			return;
		
		//we need to remove a portion of the nodes in the netword		
		int num_nodes = m_graph.getHeaderNodesMap().size();
		
		
		double x= Math.rint( num_nodes * percentageToDamage); 
		//System.out.println("Rand num of nodes to remove "+x+" out of "+ num_nodes);
		
		//If there is no nodes to remove, at least remove one.
		if(x == 0)
			x = 1;
		
		
		int i;
		for(;x > 0;x=x-1)
		{

			TreeMap<Integer,Vector<Integer>> NodeDistribution = m_graph.getDegreeDistributionAdvanced();
			Vector<Integer> bucketsDegree = new Vector<Integer>(NodeDistribution.keySet());
			
			int numbucketsDegree = bucketsDegree.size();
			
			boolean done = false;
			
			for(i= numbucketsDegree -1; i > 0 && !done;i=i-1)
			{
				//Find the bucket with the largest degree nodes
				//assumption is targets are picked at once then all attacked

							
				if(NodeDistribution.get(bucketsDegree.get(i)).size() != 0)
				{
					int target;
					if(NodeDistribution.get(bucketsDegree.get(i)).size() == 1)
					{
						target = 0;
					}
					else
					{
						Random rr = new Random();
						target = rr.nextInt(NodeDistribution.get(bucketsDegree.get(i)).size() );
						/*for(int test=0; test < 10; test++)
						{
							System.out.println("using r - with "+ (NodeDistribution.get(bucketsDegree.get(i)).size() )+" get us "+r.nextInt(NodeDistribution.get(bucketsDegree.get(i)).size() ));
							Random rr = new Random();
							System.out.println("using a new Random - with "+ (NodeDistribution.get(bucketsDegree.get(i)).size() )+" get us "+rr.nextInt(NodeDistribution.get(bucketsDegree.get(i)).size() ));
							
							
						}*/
						
					}
						
					//System.out.println("Removing node "+m_headerNodesMap.get(NodeDistribution.getID()get(bucketsDegree.get(i)).get(target)).m_id     );
					m_graph.removeNode(m_graph.getHeaderNodesMap().get(NodeDistribution.get(bucketsDegree.get(i)).get(target)).getID() );					
					
					done =true;					
				}
			}
		}
	//	System.out.println("----------------after damage      -------------");
	//	print();
	//	System.out.println("-----------------------------------------------");
		
	}
	
	public double measureDamageDirect (double percentageToDamage) throws Exception
	{
		double Dpercentage2a =0;
		double DavgDiameter2a=0;
		
		Graph backupGraph = new Graph(m_graph);
		
		for(int i =0; i < 10; i++)
		{

			m_graph = new Graph(backupGraph);
			damageDirected(percentageToDamage);
			
			//System.out.println("num in graph now" +copy2a.m_headerNodesMap.size());
			
			
			PercentageInLargestCluster measure2a = new PercentageInLargestCluster(m_graph);
			Dpercentage2a += measure2a.Measure();
			Vector<Node> NodesInLCC2a;
			NodesInLCC2a = measure2a.LCC();
			//System.out.println("Robustness Measure for Directed - Percentage in LCC "+measure2a.Measure());
		
			Diameter measure22a = new Diameter(m_graph, NodesInLCC2a);
			DavgDiameter2a += measure22a.Measure();
		//	System.out.println("Unit test for Diameter for Directed"+measure22a.Measure());
			
		}
		Dpercentage2a= Dpercentage2a/10.0;
		DavgDiameter2a = DavgDiameter2a/10.0;
		//System.out.println("Avg Robustness Measure for Directed - Percentage in LCC "+Dpercentage2a);
		//System.out.println("Avg Diameter for Directed"+DavgDiameter2a);		
		
		m_directLCC = Dpercentage2a;
		m_directDiameter = DavgDiameter2a;
		return m_directLCC;
	}
	
	Graph m_graph; // the original graph
	
	boolean m_isCalculated;

	double m_directLCC;
	double m_directDiameter;
	double percentageToDamage;
}
