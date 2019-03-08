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
package psograph.measurements;

import java.io.Serializable;
import java.util.Random;
import java.util.Vector;

import psograph.graph.Graph;
import psograph.graph.Node;


public class RandomDamage implements IGraphMeasument, Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 2115567055322269541L;

	public RandomDamage(Graph graph, double per)
	{
		m_isCalculated = false;
		m_graph = graph;
		m_randomLCC=0;
		m_randomDiameter =0;
		percentageToDamage = per;
	}
	
	public double getRandomDiameter()
	{
		return m_randomDiameter;
	}

	public double Measure() throws Exception {

		
		double value = measureDamageRandom(percentageToDamage);
		m_isCalculated = true;
		
		return value;
	}
	


	public boolean isCalculated()
	{
		return m_isCalculated;
	}
	
	private double measureDamageRandom(double percentageToDamage) throws Exception
	{
		int i =0;
		
		double RavgDiameter2a=0;
		double Rpercentage2a =0;
		
		Graph backupGraph = new Graph(m_graph);
		
		for(i =0; i < 10; i++)
		{

			m_graph = new Graph(backupGraph);
			damageRandom(percentageToDamage);
			
			PercentageInLargestCluster measure2a = new PercentageInLargestCluster(m_graph);
			Rpercentage2a += measure2a.Measure();
			Vector<Node> NodesInLCC2a;
			NodesInLCC2a = measure2a.LCC();
			//System.out.println("Robustness Measure for Random - Percentage in LCC "+measure2a.Measure());
		
			Diameter measure22a = new Diameter(m_graph, NodesInLCC2a);
			RavgDiameter2a += measure22a.Measure();
			//System.out.println("Unit test for Diameter for Random "+measure22a.Measure());
			
		}
		Rpercentage2a= Rpercentage2a/10.0;
		RavgDiameter2a = RavgDiameter2a/10.0;
	//	System.out.println("Avg Robustness Measure for Random - Percentage in LCC "+Rpercentage2a);
	//	System.out.println("Avg Diameter for Random "+RavgDiameter2a);	
		
		m_randomLCC = Rpercentage2a;
		m_randomDiameter = RavgDiameter2a;
		return m_randomLCC;
	}
	
	private void damageRandom(double percentageToDamage) throws Exception
	{
		
		if(m_graph.getHeaderNodesMap().size() == 0)
			return;
		
		//we need to remove a portion of the nodes in the netword
		
		int num_nodes = m_graph.getHeaderNodesMap().size();
		
		
		double x= Math.rint( num_nodes * percentageToDamage); 
		//System.out.println("RAN num of nodes to remove "+x+" out of "+ num_nodes);
		
		
		for(;x > 0;x=x-1)
		{
			num_nodes = m_graph.getHeaderNodesMap().size();
			Vector<Integer> vNodeId = new Vector<Integer>(m_graph.getHeaderNodesMap().keySet());
			
			Random rr = new Random();
			int target = rr.nextInt(num_nodes);

			

		//	System.out.println("Removing node "+m_headerNodesMap.get(vNodeId.get(target)).m_id);
			m_graph.removeNode(m_graph.getHeaderNodesMap().get(vNodeId.get(target)).getID());

		}

	}
	
	Graph m_graph; // the original graph
	
	boolean m_isCalculated;

	double m_randomLCC;
	double percentageToDamage;
	double m_randomDiameter;
}
