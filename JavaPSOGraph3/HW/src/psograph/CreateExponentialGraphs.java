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
package psograph;

import java.io.File;
import java.util.Random;
import java.util.Vector;

import psograph.graph.*;
import psograph.util.*;


public class CreateExponentialGraphs {

	
	int MinimumEdges = 5279;
	int MaxEdges = 5591;
	String BaseDir = "C:\\TestHeuristic3\\Seed0";	
	String OutDir = BaseDir + "\\exponential";
	String SeedDir = BaseDir + "\\seed";
	
	int numGraphsToMake = 10;
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		try
		{
		CreateExponentialGraphs createExponentialGraphs = new CreateExponentialGraphs();
		createExponentialGraphs.run();
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}

	}

	public CreateExponentialGraphs()
	{
		r = new Random();
		
	}
	
	private void run() throws Exception
	{
		File OutDirFile = new File(OutDir);
		
		Graph seed = Util.streaminSeed(new File(SeedDir));
		
		//Graph costBasisGraph = Util.streaminGraph(new File(SeedDir+"\\ExponentialCostBasis1.Graph"));
		double costBasis;// = costBasisGraph.getCost();
		
		//costBasis = 327252.8216994972;
		
		System.out.println("num nodes in seed "+seed.getNumberOfNodes());
		
		for(int i =0; i < numGraphsToMake; i++)
		{
			Graph graphCandidate = new Graph(seed);
			
			System.out.println("num edges"+graphCandidate.getTotalEdges());
			
			NodeLocationCalculator nodeLocationCalculator = new NodeLocationCalculator(seed);
			nodeLocationCalculator.calculateResults();
			
			costBasis = nodeLocationCalculator.getCost();
			
			System.out.println("cost of basis" + costBasis);
		
			int edgesToAdd = r.nextInt(1+ MaxEdges -  MinimumEdges );
			edgesToAdd +=MinimumEdges;
			
			System.out.println("Adding "+edgesToAdd+" edges");
			
			for(; graphCandidate.getTotalEdges() < edgesToAdd; )
			{
				//Olekas 15-May this was hardecode to use the constant, but this is wrong
				//It should be using the number of nodes in the graph.
				int node_from = r.nextInt(graphCandidate.getNumberOfNodes());//(GraphConstants.NUM_NODES);
				int node_to ;			
				
				Vector<Integer> vNodes = new Vector<Integer>(nodeLocationCalculator.getNode(node_from).getNeighbors().keySet());
				
				if(vNodes.size() == 0)
				{
					continue;
				}
				else
				{
					int pick = r.nextInt(vNodes.size());
					
					node_to = vNodes.get(pick);
					
				}
				
				if(node_to != node_from &&
						!graphCandidate.getNode(node_from).isConnectedTo(node_to))
				{
					graphCandidate.addConnection(node_from, node_to);
					nodeLocationCalculator.removeConnection(node_from, node_to);
					System.out.println("added an edge  totEdges:"+graphCandidate.getTotalEdges());
				}
				else
				{
					System.out.println("a duplicate edge pick again totEdges:"+graphCandidate.getTotalEdges());
				}
				
		
			}
			System.out.println(i+" num edges"+graphCandidate.getTotalEdges());
			
			
			CalculatedGraph calcGraph = new CalculatedGraph(graphCandidate);
			calcGraph.setCostBasis(costBasis);
			calcGraph.UpdatePSOCalculations();
			
			//System.out.println("lCC "+calcGraph.getLCC());
			
			if(calcGraph.getLCC()==1.0)
			{
				System.out.println("perLCC = 1.0");
				Util.streamOutCalculatedGraph(OutDirFile, "exponentialGraph"+i, calcGraph);
			}
			else
			{
				System.out.println("perLCC = "+calcGraph.getLCC()+" redoing");
				i = i -1;
			}
		}
		
	}

	Random r;
	
}
