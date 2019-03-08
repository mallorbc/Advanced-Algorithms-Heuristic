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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Random;
import java.util.TreeMap;
import java.util.Vector;

import psograph.graph.CalculatedGraph;
import psograph.graph.Edge;
import psograph.graph.Graph;
import psograph.graph.Node;
import psograph.graph.NodeLocationCalculator;
import psograph.graph.NonLinearCostFunction;
import psograph.measurements.MeanDegreeOfConnectedNeighbors;
import psograph.measurements.PercentageInLargestCluster;
import psograph.util.Util;


public class CreateOptimizedHeuristicGraphs2 {

//	int MinimumEdges = 892;
//	int MaxEdges = 936;
	double Limit = .16;
	Vector<String> BaseDir;	
	Vector<String> OutDir ;
	Vector<String> SeedDir ;
	int NodeLimit = 59;
	double perIntial = .5;
	
	String BaseDir0 = "C:\\TestHuer3\\Seed0";
	String BaseDir1 = "C:\\TestHuer3\\Seed1";
	String BaseDir2 = "C:\\TestHuer3\\Seed2";
	String BaseDir3 = "C:\\TestHuer3\\Seed3";
	String BaseDir4 = "C:\\TestHuer3\\Seed4";
	String BaseDir5 = "C:\\TestHuer3\\Seed5";
	String BaseDir6 = "C:\\TestHuer3\\BlankSeed0";
	String BaseDir7 = "C:\\TestHuer3\\BlankSeed1";
	String BaseDir8 = "C:\\TestHuer3\\BlankSeed2";
	String BaseDir9 = "C:\\TestHuer3\\BlankSeed3";
	String BaseDir10 = "C:\\TestHuer3\\BlankSeed4";
	String BaseDir11 = "C:\\TestHuer3\\BlankSeed5";
	
	String OutDir0 = BaseDir0 + "\\optimizedheuristic2";
	String SeedDir0 = BaseDir0 + "\\seed";
	
	String OutDir1 = BaseDir1 + "\\optimizedheuristic2";
	String SeedDir1 = BaseDir1 + "\\seed";
	
	String OutDir2 = BaseDir2 + "\\optimizedheuristic2";
	String SeedDir2 = BaseDir2 + "\\seed";
	
	String OutDir3 = BaseDir3 + "\\optimizedheuristic2";
	String SeedDir3 = BaseDir3 + "\\seed";
	
	String OutDir4 = BaseDir4 + "\\optimizedheuristic2";
	String SeedDir4 = BaseDir4 + "\\seed";
	
	String OutDir5 = BaseDir5 + "\\optimizedheuristic2";
	String SeedDir5 = BaseDir5 + "\\seed";
	
	String OutDir6 = BaseDir6 + "\\optimizedheuristic2";
	String SeedDir6 = BaseDir6 + "\\seed";
	
	String OutDir7 = BaseDir7 + "\\optimizedheuristic2";
	String SeedDir7 = BaseDir7 + "\\seed";
	
	String OutDir8 = BaseDir8 + "\\optimizedheuristic2";
	String SeedDir8 = BaseDir8 + "\\seed";
	
	String OutDir9 = BaseDir9 + "\\optimizedheuristic2";
	String SeedDir9 = BaseDir9 + "\\seed";
	
	String OutDir10 = BaseDir10 + "\\optimizedheuristic2";
	String SeedDir10 = BaseDir10 + "\\seed";
	
	String OutDir11 = BaseDir11 + "\\optimizedheuristic2";
	String SeedDir11 = BaseDir11 + "\\seed";

	
	int numGraphsToMake = 10;
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		try
		{
			CreateOptimizedHeuristicGraphs2 createOptimizedHeuristicGraphs = new CreateOptimizedHeuristicGraphs2();
			createOptimizedHeuristicGraphs.run();
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}

	}

	public CreateOptimizedHeuristicGraphs2()
	{
		r = new Random();
		
	}
	
	private void initData() throws Exception
	{
		BaseDir = new Vector<String>();	
		OutDir = new Vector<String>();	
		SeedDir = new Vector<String>();	
		
		
		BaseDir.add(BaseDir0);
		BaseDir.add(BaseDir1);
		BaseDir.add(BaseDir2);
		BaseDir.add(BaseDir3);
		BaseDir.add(BaseDir4);
		BaseDir.add(BaseDir5);
//		BaseDir.add(BaseDir6);
//		BaseDir.add(BaseDir7);
//		BaseDir.add(BaseDir8);
//		BaseDir.add(BaseDir9);
//		BaseDir.add(BaseDir10);
//		BaseDir.add(BaseDir11);
		
		OutDir.add(OutDir0);
		OutDir.add(OutDir1);
		OutDir.add(OutDir2);
		OutDir.add(OutDir3);
		OutDir.add(OutDir4);
		OutDir.add(OutDir5);
//		OutDir.add(OutDir6);
//		OutDir.add(OutDir7);
//		OutDir.add(OutDir8);
//		OutDir.add(OutDir9);
//		OutDir.add(OutDir10);
//		OutDir.add(OutDir11);
		
		SeedDir.add(SeedDir0);
		SeedDir.add(SeedDir1);
		SeedDir.add(SeedDir2);
		SeedDir.add(SeedDir3);
		SeedDir.add(SeedDir4);
		SeedDir.add(SeedDir5);
//		SeedDir.add(SeedDir6);
//		SeedDir.add(SeedDir7);
//		SeedDir.add(SeedDir8);
//		SeedDir.add(SeedDir9);
//		SeedDir.add(SeedDir10);
//		SeedDir.add(SeedDir11);
		
		
		
	}
	
	@SuppressWarnings("deprecation")
	private void run() throws Exception
	{
		initData();
		
		FileWriter fileLoc = new FileWriter("C:\\TestHeuristic2\\Heur2Report.txt");
		BufferedWriter mout = new BufferedWriter(fileLoc);
		
		for(int iDir=0; iDir < BaseDir.size(); iDir++)
		// BaseDir.size(); iDir++)
		//1; iDir++)
		
		{
			int minEdges = 10000000;
			int maxEdges = 0;
			
			System.out.println("Working on "+BaseDir.get(iDir));
			mout.write("Working on "+BaseDir.get(iDir));
			mout.write("\n");
			
			
			File OutDirFile = new File(OutDir.get(iDir));

			Graph seed = Util.streaminSeed(new File(SeedDir.get(iDir)));

			Graph costBasisGraph = Util.streaminGraph(new File(SeedDir.get(iDir)+"\\ExponentialCostBasis1.Graph"));
			double costBasis = costBasisGraph.getCost();

		//	System.out.println("num nodes in seed "+seed.getN());
			System.out.println("cost of basis" + costBasisGraph.getCost());
			for(int i =0; i < numGraphsToMake; i++)
			{

				//We should also think about the following heuristic: Each node makes connections with
				//all nodes within some small radius, then we identify the connected components in 
				//the resulting network and find the shortest connections between each pair 
				//of components to connect these. Is the resulting network better or worse than the
				//ones produced by PSO? 

				Graph graphCandidate = new Graph(seed);

			//	System.out.println("num edges"+graphCandidate.getTotalEdges());

				PercentageInLargestCluster lcc = new PercentageInLargestCluster(graphCandidate);

				double perLCC = 0;

				NodeLocationCalculator calculatorMain = new NodeLocationCalculator(seed);
				calculatorMain.calculateResults();

				Vector<Node> vNodes = new Vector<Node>( calculatorMain.getHeaderNodesMap().values());

				for(int j =0; j < vNodes.size(); j++)
				{
					System.out.println(j+" iter");
					TreeMap<Integer,Edge> neighbors = vNodes.get(j).getNeighbors();
					Vector<Integer> vNodeIds = new Vector<Integer> (vNodes.get(j).getNeighbors().keySet());

					Vector<Node> candNodes = new Vector<Node>();

					for(int k = 0; k < vNodeIds.size(); k++)
					{
						double length = NonLinearCostFunction.getDistance(neighbors.get(vNodeIds.get(k)).getWeight());
						if(Double.compare(length, Limit) < 0)
						{
							//	System.out.println("added this id to cand nodes "+vNodeIds.get(k));
							candNodes.add(calculatorMain.getNode(vNodeIds.get(k)));
						}
					}	
					
					
					int node_from_id = vNodes.get(j).getID();

					int numAlreadyConnected = graphCandidate.getNode(node_from_id).getDegree();
					
					
					
					double numToAccept = (candNodes.size() + (double)numAlreadyConnected ) * perIntial;
					int intNumToAccept = (int) Math.rint(numToAccept);
					int numToAdd = intNumToAccept - numAlreadyConnected;
					
					int numNodesDiscard = candNodes.size() - numToAdd;
					
					System.out.println("alreadyConnect:"+graphCandidate.getNode(node_from_id).getDegree()+
							" numToAdd:"+numToAdd+
							" numNodesDiscard:"+numNodesDiscard);
					
					if(numToAdd <= 0)
					{
						//we're done!
						
						
					}
					else
					{
						for(int ii =0; ii < numNodesDiscard; ii++)
						{
							int pick = r.nextInt(candNodes.size());

//							System.out.println("Removing edge from calculator from" +vNodes.get(j).getID() +" to "+candNodes.get(pick).getID() );

							//We just need to remove the possible connection form the candidate list
							//calculator.removeConnection(vNodes.get(j).getID(), candNodes.get(pick).getID() );


////							for(int pat=0; pat < candNodes.size();pat++)
//							{
//							System.out.print(" "+candNodes.get(pat).getID());
//							}
//							System.out.println();

							candNodes.remove(pick);

//							for(int pat=0; pat < candNodes.size();pat++)
//							{
//							System.out.print(" "+candNodes.get(pat).getID());
//							}
//							System.out.println();

						}


						for(int ii=0; ii < candNodes.size(); ii++ )
						{

							graphCandidate.addConnection(vNodes.get(j).getID(), candNodes.get(ii).getID() );
							//	System.out.println("Just added Removing edge from calculator from" +vNodes.get(j).getID() +" to "+candNodes.get(ii).getID() );
							calculatorMain.removeConnection(vNodes.get(j).getID(), candNodes.get(ii).getID() );

						}
						
					}
					


					

				}

				lcc =  new PercentageInLargestCluster(graphCandidate);
				perLCC = lcc.Measure();
			//	System.out.println(i+" num edges"+graphCandidate.getTotalEdges());

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

					//	System.out.println("Cluster ID: "+vClusterId.get(pat)+
					//			" has this many nodes "+clusterMap.get(pat).size());

						if(pat == largestClusterId)
							continue;

						double shortestConnectionSoFar = 20000;
						int shortestNodeFrom = -1; //small cluster 
						int shortestNodeTo = -1;  // largestcluster node id

						Vector<Node> clusterNodes = new Vector<Node>(clusterMap.get(pat));
						for(int c = 0; c < clusterNodes.size(); c++)
						{
							Node potFrom = calculatorMain.getNode(clusterNodes.get(c).getID());

							for(int lrg = 0; lrg < LrgclusterNodes.size(); lrg++)
							{
								Node potTo = calculatorMain.getNode(LrgclusterNodes.get(lrg).getID());

								if(! graphCandidate.getNode(potFrom.getID()).isConnectedTo(graphCandidate.getNode(potTo.getID()).getID())
										&&
										potFrom.isConnectedTo(potTo))
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
						calculatorMain.removeConnection(shortestNodeFrom, shortestNodeTo );
					}
				}

				lcc =  new PercentageInLargestCluster(graphCandidate);
				perLCC = lcc.Measure();
			//	System.out.println(i+" num edges"+graphCandidate.getTotalEdges());
				System.out.println(perLCC+" perLCC");



				// the nodes with a degree higher than 10 are connected to the 
				//nodes with a degree greater than 10 within the radius of .5.

				Vector<Node> nodes = new Vector<Node>(graphCandidate.getHeaderNodesMap().values());
				Vector<Node> highDegreenodes = new Vector<Node>();

				for(int j=0; j < nodes.size(); j++)
				{
					if(nodes.get(j).getDegree() >= NodeLimit)
					{

						highDegreenodes.add(nodes.get(j));
					}
				}
				
				System.out.println("High Degree Node total "+highDegreenodes.size());

				for(int j=0; j < highDegreenodes.size(); j++)
				{
					//Connect to another high degree node if within .5
					
					System.out.println("Hight degree iteration 2 ?? - "+ j);
					//First decide on target degree, if node already exceeds target degree,
					// add at least on connection and quit.

					int targetDegree = NodeLimit + 2 + r.nextInt(14);
					
					

					Vector<Node> candNodes = new Vector<Node>();
					for(int k=0; k <  highDegreenodes.size();k++)
					{
						if(k==j)
						{
							continue;
						}

						int from_id = highDegreenodes.get(j).getID();
						int to_id = highDegreenodes.get(k).getID();


						Node from = calculatorMain.getNode(highDegreenodes.get(j).getID());
						Node to = calculatorMain.getNode(highDegreenodes.get(k).getID());

						@SuppressWarnings("unused")
						Edge ci2 = from.getConnectionInfo(to);
						@SuppressWarnings("unused")
						Edge ci3 = calculatorMain.getNode(from.getID()).getEdgeInfo(to.getID());

						Edge ci ;
						ci = calculatorMain.getNode(from_id).getEdgeInfo(to_id);

						@SuppressWarnings("unused")
						Edge control = graphCandidate.getNode(from_id).getEdgeInfo(to_id);

						if(calculatorMain.isNodeAConnectedToNodeB(calculatorMain.getNode(from_id), calculatorMain.getNode(to_id))
								== true)
						{

							double length = NonLinearCostFunction.getDistance(ci.getWeight());


							if(Double.compare(length ,.5) <= 0 )
							{
								candNodes.add(highDegreenodes.get(k));
							}
						}

					}

					if(candNodes.size() == 0)
						continue;

					
					boolean done = false;
					
					@SuppressWarnings("unused")
					int numCands = candNodes.size();
					
					for(;  candNodes.size() > 0 && !done; )
					{
						int pick = r.nextInt(candNodes.size());

						Node from = calculatorMain.getNode(highDegreenodes.get(j).getID());
						Node to = calculatorMain.getNode(candNodes.get(pick).getID());
						
						if(graphCandidate.getNode(from.getID()).getDegree() > targetDegree)
						{
							done = true;
							continue;
						}
						


						@SuppressWarnings("unused")
						Edge ci = from.getConnectionInfo(to);
						@SuppressWarnings("unused")
						Edge ci2 = calculatorMain.getNode(from.getID()).getEdgeInfo(to.getID());


						graphCandidate.addConnection(from.getID(), to.getID());

						calculatorMain.removeConnection(from.getID(), to.getID());
						
						candNodes.remove(pick);

					}

				}

				//For any node that has a degree of five or less, we will evaluate its mean 
				//degree of it connected neighbors.  If it less than six, we will add a 
				//connection to nodes within a radius of .5 until we reach a mean 
				//average degree of at least six.
				
				int minDegree = 5;
				
				Vector<Node> candMinDegreeNodes = new Vector<Node>();
				
				vNodes = new Vector<Node>( graphCandidate.getHeaderNodesMap().values());
				for(int k=0; k < vNodes.size(); k++)
				{
					//System.out.println("Node degree:"+vNodes.get(k).getDegree());
					if(vNodes.get(k).getDegree() <= minDegree)
					{
						//System.out.println("Got a node with a degree less than or equal to 5");
						candMinDegreeNodes.add(vNodes.get(k));
					}
				}
				

				
				int numCandMinDegreeNodes =candMinDegreeNodes.size();
				
				System.out.println("Modify degree totals - "+ numCandMinDegreeNodes);
				
				for(int k=0; k < numCandMinDegreeNodes; k++)
				{
					System.out.println("Modify degree totals iter "+k);
					
					MeanDegreeOfConnectedNeighbors meanDegreeOfConnectedNeighbors =
						new MeanDegreeOfConnectedNeighbors(graphCandidate);
					meanDegreeOfConnectedNeighbors.Measure();
					TreeMap<Integer, Double> meanDegreeData = meanDegreeOfConnectedNeighbors.ById();
					
					
					int pick = r.nextInt(candMinDegreeNodes.size());
					
					
					if( graphCandidate.getNode(candMinDegreeNodes.get(pick).getID()).getDegree() < 6 &&						
							Double.compare(meanDegreeData.get(candMinDegreeNodes.get(pick).getID()), 6.0 ) < 0)
					{

						Node from = graphCandidate.getNode(candMinDegreeNodes.get(pick).getID());
						int numNeighbors = from.getNeighbors().size();
						double degreeSum = 0;
						//I am worried about roundoff error so I am going to sum it again myself
						Vector<Integer> vNodeIds = new Vector<Integer>(from.getNeighbors().keySet()); 
						for(int q=0; q< numNeighbors; q++)
						{
							degreeSum += graphCandidate.getNode(vNodeIds.get(q)).getDegree();
						}
						
												
					//	System.out.println("do Something - mean degree is " +meanDegreeData.get(candMinDegreeNodes.get(pick).getID()));
					//	System.out.println("   mean sum is " +degreeSum +
					//			"   numNeigh "+numNeighbors);						
						
						boolean done = false;
						while (!done)
						{
							//Get list of nodes within .5 that we are not connected to
							
							Vector<Integer> vdataIds = 
								new Vector<Integer>(calculatorMain.getNode(from.getID()).getNeighbors().keySet());
							
							Vector<Integer> cand = new Vector<Integer>();
							for(int q=0; q< vdataIds.size(); q++)
							{
								if(Double.compare(calculatorMain.getNode(from.getID()).getNeighbors().get(vdataIds.get(q)).getWeight(),.5) <= 0  
										&&
										graphCandidate.getNode(vdataIds.get(q)).getDegree() >= 10
								)
								{
									cand.add(vdataIds.get(q));
								}
							}
							
							if(cand.size() <= 0)
							{
						//		System.out.println("no good nodes to connect to");
								done = true;
								continue;
								
							}
							
							int pick2= r.nextInt(cand.size());
							
							numNeighbors++;
							degreeSum += graphCandidate.getNode(cand.get(pick2)).getDegree();
							
						//	System.out.println("     Added a connection of "+ graphCandidate.getNode(cand.get(pick2)).getDegree());
							graphCandidate.addConnection(from.getID(), cand.get(pick2));
							calculatorMain.removeConnection(from.getID(), cand.get(pick2));
							
							if(Double.compare(degreeSum/(double)numNeighbors, 6.0) >= 0)
							{
								done = true;
							}
							else
							{
							//	System.out.println("   mean sum is " +degreeSum +
							//			"   numNeigh "+numNeighbors);
								double v= degreeSum/(double)numNeighbors;
							//	System.out.println("   mean sum div numNeigh" +v);
							
							}
						}
						
						//System.out.println("did Something");
						//System.out.println("   mean sum is " +degreeSum +
						//		"   numNeigh "+numNeighbors);		
						
						
					}
					else
					{
					//	System.out.println("Just fine");
						//candMinDegreeNodes.remove(pick);
					}
					candMinDegreeNodes.remove(pick);
					
					
					
				}	
				

				CalculatedGraph calcGraph = new CalculatedGraph(graphCandidate);
				calcGraph.setCostBasis(costBasis);

				
				//System.out.println("Write out graph ...");
				Util.streamOutCalculatedGraph(OutDirFile, "heuristicGraph"+i, calcGraph);
				
				System.out.println("calc graph ...");
				//calcGraph.UpdateCalcuations();
				calcGraph.UpdatePSOCalculations();
				System.out.println("fitness:"+calcGraph.getFitnessValue()+" numEdges:"+calcGraph.getTotalEdges());
				mout.write("fitness:"+calcGraph.getFitnessValue()+" numEdges:"+calcGraph.getTotalEdges());
				mout.write("\n");
				if(calcGraph.getTotalEdges() > maxEdges)
				{
					maxEdges = calcGraph.getTotalEdges();
				}
				
				if( calcGraph.getTotalEdges() <  minEdges)
				{
					minEdges = calcGraph.getTotalEdges();
				}

				//System.out.println("lCC "+calcGraph.getLCC());
			//	System.out.println("write out calc graph ...");
				Util.streamOutCalculatedGraph(OutDirFile, "heuristicGraph"+i, calcGraph);


			}
			System.out.println("MinEdges: "+minEdges+" MaxEdges: "+maxEdges);
			mout.write("MinEdges: "+minEdges+" MaxEdges: "+maxEdges);
			mout.write("\n");
		}
		
		mout.flush();
		mout.close();
		
	}

	Random r;
	
}

