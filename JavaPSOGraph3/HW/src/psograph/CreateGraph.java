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

import java.io.*;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.TreeMap;
import java.util.Vector;

import psograph.graph.*;
import psograph.measurements.PercentageInLargestCluster;
import psograph.util.Util;



public class CreateGraph 
{
	
	public Vector<Node> Remove_Node_Index_by_ID(int Id_to_remove,Vector<Node> Node_Vector) {
		Vector<Node> modified_node_vector = Node_Vector;
		Boolean found_Id = false;
		int id_index = 0;
		
		for(int i = 0; i<modified_node_vector.size(); i++) {
			if(modified_node_vector.get(i).getID() == Id_to_remove) {
				id_index = i;
				modified_node_vector.remove(i);
				break;
			}
			
		}
		
		System.out.println("Test function");
		return modified_node_vector;
		
	}
	


	//This is the directory where the Seed is generated.  The Seed is the Node configuration
	//with no edges
	File m_SeedDirectory;  
	
	//This is where the graphs will reside TODO need more
	File m_GraphDirectory;
	
	Graph m_graphSeed;
	Graph canditate;
	
	double m_basisCost;
	NodeLocationCalculator m_nodeLoc;
	
	int seed = 0;
	int candidate = 0;
	
	public CreateGraph()
	{
		m_SeedDirectory = new File("C:\\TestHeuristic3\\StandardStartingPoint");
	}
	
	private void generateSeed() throws Exception
	{
		/* commenting out to use standardized seed
		//Generate a new graph for as the Seed
		m_graphSeed = new Graph(GraphConstants.NUM_NODES);

		//graphSeed.printWithLocationAndWeights();
		System.out.println("--------------------------------------");
		System.out.println("Saving Graph Seed "+seed);
		// stream out seed
		
		Util.streamOutSeed(m_SeedDirectory, m_graphSeed);   
			*/	
		File m_SeedDirectory2 = Util.CreateSeedDirectory();
		m_GraphDirectory = Util.CreateCalculatedGraphDirectory(m_SeedDirectory2);

		
		m_graphSeed = Util.streaminSeed(m_SeedDirectory);
		
		System.out.println("--------------------------------------");
		
		m_nodeLoc = new NodeLocationCalculator(m_graphSeed, false);
		m_nodeLoc.calculateResults();
		
		
		
	}
	
	private void calculateCostBasis() throws Exception
	{
		m_basisCost =0;
		//CostBasis costBasis = new CostBasis(m_graphSeed);
			
		//System.out.println("-----------exponentialCostBasis-----------------");
		//costBasis.generate(GraphConstants.MAX_CONNECTIONS+300);
		//System.out.println("Total edges "+exponentialCostBasis.getTotalEdges());			
		//Util.streamOutExponentialGraph(m_SeedDirectory, costBasis,1);
		
		//Used Standard Cost Basis
		File CostBasisFile = new File(m_SeedDirectory+"\\ExponentialCostBasis1.Graph");
		CostBasis costBasis = Util.streaminCostBasis(CostBasisFile);
		m_basisCost = costBasis.getCost();


		PercentageInLargestCluster expoLCC = new PercentageInLargestCluster(costBasis);
		double valueCostBasisLCC = expoLCC.Measure();
		if(Double.compare(valueCostBasisLCC,1.0) != 0)
			System.out.println("valueExpoLCC is not equal to 1.0 : "+valueCostBasisLCC + "differ of :" + (1.0 - valueCostBasisLCC));	
		
		//Do we want to perform some measurements on this guy
		
	}
	
	private void connnectCandidate() throws Exception
	{
		
		Random pickNode = new Random();

		NodeLocationCalculator workingNodeLoc = new NodeLocationCalculator(m_nodeLoc, false);
		
		//workingNodeLoc.printWithLocationAndWeights();
		
		canditate = new Graph(m_graphSeed);
		
		Vector<Node> v_Nodes = new Vector<Node>(m_graphSeed.getHeaderNodesMap().values());
		int jj;

		//
		// We should really be picking each node randomly.
		// 
		//  three different ideas on how to connect initially
		//  n = number of nodes
		//  1)Put in n initial connections.  Where each node
		//  will have a connection to another node
		//
		//  2) have n-1 connections.  This will gives us a chance to make a tree connecting
		// all nodes
		//
		// 3) Make a MST of the nodes
		//
		
		int num_of_nodes = v_Nodes.size();				
		
		for(jj=0; jj < num_of_nodes ; jj++)
		{
			int t_id = pickNode.nextInt(v_Nodes.size());
			
			Node n = workingNodeLoc.chooseCloseNode(v_Nodes.get(t_id));
			if(n == null)
			{
				throw new Exception("ERROR - null node returned for choose close node");
			}
			else
			{
				//System.out.println("a real node returned for choose close node");
				Edge ci = n.getConnectionInfo(v_Nodes.get(t_id));
				canditate.addConnection(v_Nodes.get(t_id).getID(), n.getID(), ci.getWeight());
				
				//Remove from working NodeLoc so we don't hit it in random phase
				workingNodeLoc.removeConnection(v_Nodes.get(t_id).getID(), n.getID());
				
				//System.out.println("v_Nodes.get(t_id).m_id is "+v_Nodes.get(t_id).m_id);
				//System.out.println("t_id is "+t_id);
				v_Nodes.remove(t_id);
				//System.out.println("n.m_id is "+n.m_id);
				//int t2 = v_Nodes.indexOf(n);
				//System.out.println("t2 is "+t2);							
			}
		}
		
	//	System.out.println("Total edges after first connect "+canditate.getTotalEdges());
						
		//System.out.println("Print out of candiate after initial connectiveness");
		//canditate.printWithLocationAndWeights();
		//printMFile(canditate);
		//System.out.println("End of Print out of candiate after initial connectiveness");
		
		
		v_Nodes = new Vector<Node>(m_graphSeed.getHeaderNodesMap().values());
		//int MAX_CONNECTIONS = (6 * v_Nodes.size());
		
		//System.out.println("MAX_CONNECTIONS "+MAX_CONNECTIONS);
		
		v_Nodes = new Vector<Node>(m_graphSeed.getHeaderNodesMap().values());
						
		for (jj=0; jj < GraphConstants.MAX_CONNECTIONS; )
		{
			
			int t_id = pickNode.nextInt(v_Nodes.size());
			Node NodeToConnect = v_Nodes.get(t_id);

			Node n = workingNodeLoc.chooseNode(NodeToConnect);
			if(n == null)
			{
				System.out.println("null node returned for choose random/close node");
			}
			else
			{
				Edge ci = n.getConnectionInfo(NodeToConnect);
				canditate.addConnection(NodeToConnect.getID(), n.getID(), ci.getWeight());
				workingNodeLoc.removeConnection(NodeToConnect.getID(), n.getID());
				jj=jj+2;
			}
		}
		
	//	System.out.println("Random :"+workingNodeLoc.m_random+" closest :"+workingNodeLoc.m_closest);
		if(canditate.getTotalEdges() != GraphConstants.MAX_CONNECTIONS)
			throw new Exception("Total edges "+canditate.getTotalEdges()+ " does not equal "+GraphConstants.MAX_CONNECTIONS);
		
	}
	
	private void connnectCandidate2() throws Exception
	{

		NodeLocationCalculator workingNodeLoc = new NodeLocationCalculator(m_nodeLoc, false);
		
		Vector<Node> v_Nodes = new Vector<Node>(workingNodeLoc.getHeaderNodesMap().values());
		Graph MST_candidate = new Graph(m_graphSeed);
		
		//Gets get number of nodes
		int num_of_nodes = v_Nodes.size();
		//used to store information about the connection, may not be needed, may be thought as Nearest[]
		class Edge_Info{
			public int node_id;
			public double Edge_Lengths;
		}
		
		//stores information regarding nodes, such as whether or not in tree, the id, and the parent
		class Node_Edge_Info{
		    public int node_id;
		    public int nearest;
		    public double nearest_weight;
		    public Boolean InTree;
		    public int node_parent;
		}
		//used solely to make it easier to add ids to the Node_Edge_Info
		Vector<Integer> All_Ids = new Vector<>();
		//Vector of type Node_Edge_Info, holds the values found
		Vector <Node_Edge_Info> All_Node_data = new Vector<>();
		//solely used to make it easier to initialize the values of the parents of the nodes
		Vector<Integer> Node_Parent = new Vector<>();

		
		//Should be GOOD; Just added the makes a vector to hold node ids
		for(int k = 0; k< num_of_nodes; k++) {
			All_Ids.add(v_Nodes.get(k).getID());

		}
		//Should be GOOD; Just initializes the values of parents to -1
		for(int p = 0; p< num_of_nodes; p++) {
			Node_Parent.add(-1);
		}
		//Initializes data structure
		//PROBABLY GOOD
		//BLAKE: we should have a parent node that starts as being in the tree
		for(int y = 0; y< num_of_nodes; y++) {
			Node_Edge_Info Entry_Node = new Node_Edge_Info();
			Entry_Node.InTree = false;
			Entry_Node.node_parent = -1;
			Entry_Node.node_id = All_Ids.get(y);
			Entry_Node.nearest_weight = 100000000;
			All_Node_data.add(Entry_Node);
			
		}
		//used to store the id of the current node (start at 0)
		int current_node_id = 0;
		Node current_node;
		//a node needs to be in the tree to start, zero will be that node
		All_Node_data.get(0).InTree = true;
		All_Node_data.get(0).node_parent = -1;
		//loops for all nodes and finds there minimum connection to start
		int id = 0;
		for(int r = 0; r < num_of_nodes; r++) {
			double minimum = 100000;
			int nearest_node = 0;
			//get the node object for the ID we have
			current_node = v_Nodes.get(current_node_id);
			//get all the neighbors of this node
			TreeMap<Integer, Edge> node_connections = current_node.getNeighbors();
			//And create an iterator to go through them all
			Iterator<Integer> neighbors = node_connections.keySet().iterator();
			//Update array going through every neighbor
			System.out.println(current_node_id);
			while(neighbors.hasNext()) {
				
				id = neighbors.next();
				//Any neighbors that aren't in the tree, their parent is now current node
				if(!All_Node_data.get(id).InTree) {
					System.out.print(id + ", ");
					//All_Node_data.get(id).node_parent = current_node.getID();
					//Find the minimum edge
					if(current_node.getEdgeInfo(id).getWeight() < minimum) {
						minimum = current_node.getEdgeInfo(id).getWeight();
						nearest_node = id;
					}
					//Update nearest
					if(All_Node_data.get(id).nearest_weight > current_node.getEdgeInfo(id).getWeight()) {
						All_Node_data.get(id).nearest_weight = current_node.getEdgeInfo(id).getWeight();
						All_Node_data.get(id).nearest = current_node.getID();
						All_Node_data.get(id).node_parent = current_node.getID();
					}
				}
				
			}
			System.out.println();
			All_Node_data.get(nearest_node).InTree = true;
			current_node_id = nearest_node;
		}
		System.out.println("Arrays fully updated");
		int counter = 0;
		for(int q = 0; q < num_of_nodes; q++) {
			if(All_Node_data.get(q).node_parent != -1) {
				
				MST_candidate.addConnection(q, All_Node_data.get(q).node_parent);
				counter++;
			}
		}
		System.out.println("edges added: "+counter);
		canditate = MST_candidate;
	}
	
	
	public void doWork() throws Exception
	{
		try
		{
			for( seed=0; seed < GraphConstants.MAX_SEEDS; seed++)
			{
				
				generateSeed();

				//Now to make some graphs to be used for normalizing the cost
				calculateCostBasis();

				//nodeLoc.printWithLocationAndWeights();

				for(candidate=0; candidate < 1; candidate++)
				{
					//canditate.printWithLocationAndWeights();
					//connnectCandidate();
					connnectCandidate2();
					measureAndOutputCandidate();					  
				}
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw e;
		}
	}
	
	private void measureAndOutputCandidate() throws Exception
	{

		CalculatedGraph calculatedCanditate = new CalculatedGraph(canditate); 
		calculatedCanditate.setCostBasis(m_basisCost);
	//	calculatedCanditate.UpdateCalcuations();
		calculatedCanditate.UpdatePSOCalculations();
		
		System.out.println("------------------Begin Measurements-----------");
		System.out.println("Avg Robustness Measure for Random - Percentage in LCC "+calculatedCanditate.getRandomLCC());
		System.out.println("Avg Robustness Measure for Random - Diameter in LCC "+calculatedCanditate.getRandomDiameter());
		System.out.println("Avg Robustness Measure for Directed - Percentage in LCC "+calculatedCanditate.getDirectLCC());
		System.out.println("Avg Robustness Measure for Directed - Diameter in LCC "+calculatedCanditate.getDirectDiameter());
		System.out.println("Connectivity Measure - AISPL "+ calculatedCanditate.getAISPL());
		System.out.println("Cost Measure - summation weight costs "+ calculatedCanditate.getCost());	
		System.out.println("Cost Basis -                      "+ calculatedCanditate.getCostBasis());
	    double t = calculatedCanditate.getCost() / calculatedCanditate.getCostBasis();
		System.out.println("Cost Basis ratio - "+t );
		System.out.println("Fitness Value -  "+calculatedCanditate.getFitnessValue());
		System.out.println("Diameter Value -  "+calculatedCanditate.getDiameter());
		System.out.println("ClusteringCoefficient - "+calculatedCanditate.getClusteringCoefficient());
		System.out.println("Per LCC -  "+calculatedCanditate.getLCC());

		System.out.println("------------------End Measurements-------------");
		calculatedCanditate.printWithLocationAndWeights();

		Util.streamOutCalculatedGraph(m_GraphDirectory, candidate, calculatedCanditate);	

	}
	

	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception
	{
		
		try
		{
			CreateGraph createGraph = new CreateGraph();
			createGraph.doWork();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			throw e;
		}
	}

}
