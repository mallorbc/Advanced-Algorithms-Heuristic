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
import java.util.TreeMap;
import java.util.Vector;

import psograph.graph.CalculatedGraph;
import psograph.graph.Edge;
import psograph.graph.Graph;
import psograph.graph.Node;
import psograph.graph.NodeLocationCalculator;
import psograph.graph.NonLinearCostFunction;
import psograph.measurements.*;
import psograph.util.Util;



public class CreateHeuristicGraphs
{

//	int MinimumEdges = 892;
//	int MaxEdges = 936;
	double Limit = .17;
	String BaseDir = "C:\\TestHeuristic3\\Seed0";	
	String OutDir = BaseDir + "\\optimizedheuristic";
	String SeedDir = BaseDir + "\\seed";
	
	
	int numGraphsToMake = 1;
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		try
		{
			CreateHeuristicGraphs createHeuristicGraphs = new CreateHeuristicGraphs();
			createHeuristicGraphs.run();
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}

	}

	public CreateHeuristicGraphs()
	{
		r = new Random();
		
	}
	
	private void run() throws Exception
	{
		File OutDirFile = new File(OutDir);
		
		Graph seed = Util.streaminSeed(new File(SeedDir));
		
		Graph costBasisGraph = Util.streaminGraph(new File(SeedDir+"\\ExponentialCostBasis1.Graph"));
		double costBasis = costBasisGraph.getCost();
			//;70785.81725041081;//
		
		System.out.println("num nodes in seed "+seed.getNumberOfNodes());
		System.out.println("cost of basis" + costBasis);
		for(int i =0; i < numGraphsToMake; i++)
		{
			
			//We should also think about the following heuristic: Each node makes connections with
			//all nodes within some small radius, then we identify the connected components in 
			//the resulting network and find the shortest connections between each pair 
			//of components to connect these. Is the resulting network better or worse than the
			//ones produced by PSO? 
			
			Graph graphCandidate = new Graph(seed);
			
			System.out.println("num edges"+graphCandidate.getTotalEdges());
			
			PercentageInLargestCluster lcc = null;

			
			double perLCC = 0;
			

			NodeLocationCalculator calculator = new NodeLocationCalculator(seed);
			calculator.calculateResults();

			Vector<Node> vNodes = new Vector<Node>( calculator.getHeaderNodesMap().values());

			for(int j =0; j < vNodes.size(); j++)
			{
				System.out.println(j+" iter");
				TreeMap<Integer,Edge> neighbors = vNodes.get(j).getNeighbors();
				Vector<Integer> vNodeIds = new Vector<Integer> (vNodes.get(j).getNeighbors().keySet());
				for(int k = 0; k < vNodeIds.size(); k++)
				{
					double length = NonLinearCostFunction.getDistance(neighbors.get(vNodeIds.get(k)).getWeight());
					if(Double.compare(length, Limit) < 0)
					{
						graphCandidate.addConnection(vNodes.get(j).getID(), vNodeIds.get(k) );
						calculator.removeConnection(vNodes.get(j).getID(), vNodeIds.get(k) );
					}
				}			

			}

			lcc =  new PercentageInLargestCluster(graphCandidate);
			perLCC = lcc.Measure();
			System.out.println(i+" num edges"+graphCandidate.getTotalEdges());
			
			if(Double.compare(perLCC,1.0) != 0)
			{
				TreeMap<Integer, Vector<Node>> clusterMap = lcc.getNodesByCluster();
				Vector<Integer> vClusterId = new Vector<Integer>(clusterMap.keySet());
				
				//Find largest Cluster
				int clusterSize = -1;
				int largestClusterId = -1;
				for(int pat = 0; pat < vClusterId.size(); pat++)
				{
					System.out.println("Cluster ID: "+vClusterId.get(pat)+
							" has this many nodes "+clusterMap.get(pat).size());
					
					if(clusterSize < clusterMap.get(pat).size())
					{
						clusterSize = clusterMap.get(pat).size();
						largestClusterId = pat;
					}
					
				}
				
				Vector<Node> LrgclusterNodes = new Vector<Node>(clusterMap.get(largestClusterId));
				//Now connect clusters to largest cluster
				for(int pat = 0; pat < vClusterId.size() ; pat++)
				{
					
					System.out.println("Cluster ID: "+vClusterId.get(pat)+
							" has this many nodes "+clusterMap.get(pat).size());
					
					if(pat == largestClusterId)
						continue;
					
					double shortestConnectionSoFar = 20;
					int shortestNodeFrom = -1; //small cluster 
					int shortestNodeTo = -1;  // largestcluster node id
					
					Vector<Node> clusterNodes = new Vector<Node>(clusterMap.get(pat));
					for(int c = 0; c < clusterNodes.size(); c++)
					{
						Node potFrom = calculator.getNode(clusterNodes.get(c).getID());
						
						for(int lrg = 0; lrg < LrgclusterNodes.size(); lrg++)
						{
							Node potTo = calculator.getNode(LrgclusterNodes.get(lrg).getID());
							
							if(! graphCandidate.getNode(potFrom.getID()).isConnectedTo(graphCandidate.getNode(potTo.getID()).getID()))
							{
								if(Double.compare(
										potFrom.getConnectionInfo(LrgclusterNodes.get(lrg)).getWeight(),
										shortestConnectionSoFar) < 0)
								{
									shortestConnectionSoFar = potFrom.getConnectionInfo(potTo).getWeight();
									shortestNodeFrom =  potFrom.getID();
									shortestNodeTo = potTo.getID();
								}
							}
							
						}

					}
					graphCandidate.addConnection(shortestNodeFrom, shortestNodeTo);
					calculator.removeConnection(shortestNodeFrom, shortestNodeTo );
				}
			}
			
			lcc =  new PercentageInLargestCluster(graphCandidate);
			perLCC = lcc.Measure();
			System.out.println(i+" num edges"+graphCandidate.getTotalEdges());
			System.out.println(perLCC+" perLCC");
		




			CalculatedGraph calcGraph = new CalculatedGraph(graphCandidate);
			calcGraph.setCostBasis(costBasis);
			calcGraph.UpdatePSOCalculations();
			System.out.println(calcGraph.getFitnessValue()+" fitness");
			
			//System.out.println("lCC "+calcGraph.getLCC());
			
			Util.streamOutCalculatedGraph(OutDirFile, "heuristicGraph"+i, calcGraph);
			
			
		}
		
	}

	Random r;
	
}
