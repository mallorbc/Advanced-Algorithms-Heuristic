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
		
		//Random pickNode = new Random();

		NodeLocationCalculator workingNodeLoc = new NodeLocationCalculator(m_nodeLoc, false);
		
		//workingNodeLoc.printWithLocationAndWeights();
		
		canditate = new Graph(m_graphSeed);
		
		Vector<Node> v_Nodes = new Vector<Node>(workingNodeLoc.getHeaderNodesMap().values());
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
		    public Boolean IsParent;
		    
		}
		//used solely to make it easier to add ids to the Node_Edge_Info
		Vector<Integer> All_Ids = new Vector<>();
		//Vector of type Node_Edge_Info, holds the values found
		Vector <Node_Edge_Info> All_Node_data = new Vector<>();
		//solely used to make it easier to intialize the values of the parents of the nodes
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
			Entry_Node.IsParent = false;
			All_Node_data.add(Entry_Node);
			
		}
		//used to store the id of the current node (start at 0)
		int current_node_id = 0;
		//used to store the id of the node that is found with chooseCloseNode; I think it is the closest node
		int connection_node_id = 0;
		//Not used yet but would hold the value of whether or nor the current node is in the tree
		Boolean current_node_InTree_status;
		//Not used currently but would hold the value of the current nodes parent
		int current_node_parent;
		//Fairly certain this is not needed; keeping here just in case
		//TreeMap<Integer,Edge> node_connections;
		//stores the edge data structure of the two nodes;  The weight of which can be found by calling getWeight();
		Edge connection_weight;
		//node that will be connected (if this original node is added we have an error)
		Node connection_node = new Node(0, -1, -1);
		//current node
		Node current_node;
		//a node needs to be in the tree to start, zero will be that node
		All_Node_data.get(0).InTree = true;
		All_Node_data.get(0).node_parent = 0;
		//loops for all nodes and finds there minimum connection to start
		for(int r = 0; r<num_of_nodes; r++) {
			//get the node object for the ID we have
			current_node = v_Nodes.get(current_node_id);
			//gets the node id
			current_node_id = current_node.getID();
			//get all the neighbors of this node
			TreeMap<Integer, Edge> node_connections = current_node.getNeighbors();
			//Now find the smallest neighbor that is not in the tree
			//Edge for connected node (start with super high weight for our minimum weight finding)
			Edge lowest_connection_edge = new Edge(1000000);
			//used to store the current_edge that will be tested
			Edge current_edge = new Edge(100000);
			//You can't get an edge with an id of the current node, this value will be used to prevent that
			int edge_index;
			//Start searching for closest neighbor not in tree
			for(int h = 0; h < node_connections.size(); h++) {
				edge_index = h;
//				if(h==r) {
//					edge_index = edge_index + 1;
//					if(edge_index>node_connections.size()) {
//						edge_index = 0;
//					}
//					
//				}
				//Update parent vector and make this the parent of that node
//				if(!All_Node_data.get(edge_index).InTree) {
//					All_Node_data.get(edge_index).node_parent = All_Node_data.get(current_node_id).node_id;
//				}
				//you can't an edge to a node with the same id
				current_edge = node_connections.get(edge_index);
				if(current_edge==null) {
					while(current_edge == null) {
						edge_index++;
						if(edge_index>node_connections.size()) {
							edge_index =0;
						}
						if(All_Node_data.get(edge_index).InTree) {
							edge_index++;
						}
						current_edge = node_connections.get(edge_index);	
					}
					System.out.println("Current edge is null, expect an error");
				}
				//Make sure it is NOT in the tree and is the minimum neighbor
				//update nearest data
				if((!All_Node_data.get(edge_index).InTree) && (current_edge.getWeight() < lowest_connection_edge.getWeight())) {
					//update nearest with the minimum that isn't in the tree
					All_Node_data.get(current_node_id).nearest = All_Node_data.get(edge_index).node_id;
					All_Node_data.get(current_node_id).nearest_weight = current_edge.getWeight();
					All_Node_data.get(edge_index).node_parent = current_node_id;
					lowest_connection_edge = current_edge;
					
				}
				
			}
			int parent_id_to_add;
			int node_id_to_add;
			double connection_weight_canidate = lowest_connection_edge.getWeight();
			node_id_to_add = All_Node_data.get(current_node_id).nearest;
//			parent_id_to_add = All_Node_data.get(r).node_parent;
			parent_id_to_add = All_Node_data.get(node_id_to_add).node_parent;
			if(!All_Node_data.get(node_id_to_add).InTree && !All_Node_data.get(node_id_to_add).IsParent) {
			canditate.addConnection(parent_id_to_add,node_id_to_add,connection_weight_canidate);
			All_Node_data.get(node_id_to_add).InTree = true;
			//adding or removing this radically changes graph
			//All_Node_data.get(parent_id_to_add).IsParent = true;
			All_Node_data.get(node_id_to_add).node_parent=current_node_id;
			workingNodeLoc.removeConnection(parent_id_to_add, node_id_to_add);
			current_node_id = node_id_to_add;
			}
			
		}
		System.out.println("test");
			//search through nearest to find smallest edge
//			double min_weight = 1000000;
//			for(int z = 0; z < All_Node_data.size(); z++) {
//				//go through every node and look at the weight to its nearest edge
//				connection_weight = v_Nodes.get(All_Node_data.get(z).node_id).getConnectionInfo(v_Nodes.get(All_Node_data.get(z).nearest));
//				if(connection_weight.getWeight() < min_weight) {
//					min_weight = connection_weight.getWeight();
//					//Choose our connection node to be the minimum
//					connection_edge = node_connections.get(z);
//					connection_node_id = All_Node_data.get(z).nearest;
//					current_node_id = All_Node_data.get(z).node_id;
//				}
//				
//			}
//			//gets edge data structure from complete graph; can be used to find weight
//			//Do we really need to do this?
//			connection_weight = connection_node.getConnectionInfo(v_Nodes.get(current_node_id));
//			//nice print for debugging
//			System.out.println("Current node: "+ current_node.getID() + " Connection_node:"+connection_node_id +" weight:"+connection_weight.getWeight());
//			//adds the connection to the graph
//			canditate.addConnection(current_node_id, connection_node_id, connection_weight.getWeight());
//			//adds the relevant data to the custom structure vector, build if statements off this data
//			All_Node_data.get(current_node_id).InTree = true;
//			All_Node_data.get(connection_node_id).InTree = true;
//			All_Node_data.get(connection_node_id).node_parent = current_node_id;
//			//removes connection from complete graph.
//			workingNodeLoc.removeConnection(current_node_id, connection_node_id);
//			//for the next part of the loop, start at this node
//			current_node_id = connection_node_id;
//		}
		
		
		
		
		
		
		//COMMENTED OUT REST OF CODE FOR DEBUGGING
		//System.out.println("End of builing data structure");
		
		
//		for(int i=0; i< num_of_nodes ; i++) {
//			//Node_Parent.add(-1);
//			//Don't want to chooose randomly, want to start at 0 and work our way along the tree
//			Node temp_node = workingNodeLoc.chooseCloseNode(v_Nodes.get(i));
//			//Node temp_node2 = null;
//			int num_of_connections = temp_node.getDegree();
//			TreeMap<Integer,Edge> all_nodes_edges = temp_node.getNeighbors();
//			
//			
//			Vector<Double> edge_weights = null;
//			for(int j = 0; j<num_of_connections; j++) {
//				//int tree_key = temp_tree.keySet().toArray()[j];
//				edge_weights.add((Double) all_nodes_edges.values().toArray()[j]);
//				 //edge_weights.add(temp_tree.get(j).getWeight());
//				//Node_data.
//				//Node_data.add(i,edge_weight);
//				
//				//.
//			}
//			//Edge_Lengths.add(workingNodeLoc.chooseCloseNode(v_Nodes.get(i)).);
//			
//		}
//		
//		
//		for(jj=0; jj < num_of_nodes ; jj++)
//		{
//			//int t_id = pickNode.nextInt(v_Nodes.size());
//			//will start at index 0 node
//			int t_id = jj;
//			
//			Node n = workingNodeLoc.chooseCloseNode(v_Nodes.get(t_id));
//			if(n == null)
//			{
//				throw new Exception("ERROR - null node returned for choose close node");
//			}
//			else
//			{
//				//System.out.println("a real node returned for choose close node");
//				Edge ci = n.getConnectionInfo(v_Nodes.get(t_id));
//				canditate.addConnection(v_Nodes.get(t_id).getID(), n.getID(), ci.getWeight());
//				
//				//Remove from working NodeLoc so we don't hit it in random phase
//				workingNodeLoc.removeConnection(v_Nodes.get(t_id).getID(), n.getID());
//				
//				//System.out.println("v_Nodes.get(t_id).m_id is "+v_Nodes.get(t_id).m_id);
//				//System.out.println("t_id is "+t_id);
//				v_Nodes.remove(t_id);
//				//System.out.println("n.m_id is "+n.m_id);
//				//int t2 = v_Nodes.indexOf(n);
//				//System.out.println("t2 is "+t2);							
//			}
//		}
//		
//	//	System.out.println("Total edges after first connect "+canditate.getTotalEdges());
//						
//		//System.out.println("Print out of candiate after initial connectiveness");
//		//canditate.printWithLocationAndWeights();
//		//printMFile(canditate);
//		//System.out.println("End of Print out of candiate after initial connectiveness");
//		
//		
//		v_Nodes = new Vector<Node>(m_graphSeed.getHeaderNodesMap().values());
//		//int MAX_CONNECTIONS = (6 * v_Nodes.size());
//		
//		//System.out.println("MAX_CONNECTIONS "+MAX_CONNECTIONS);
//		
//		v_Nodes = new Vector<Node>(m_graphSeed.getHeaderNodesMap().values());
//						
//		for (jj=0; jj < GraphConstants.MAX_CONNECTIONS; )
//		{
//			
//			int t_id = pickNode.nextInt(v_Nodes.size());
//			Node NodeToConnect = v_Nodes.get(t_id);
//
//			Node n = workingNodeLoc.chooseNode(NodeToConnect);
//			if(n == null)
//			{
//				System.out.println("null node returned for choose random/close node");
//			}
//			else
//			{
//				Edge ci = n.getConnectionInfo(NodeToConnect);
//				canditate.addConnection(NodeToConnect.getID(), n.getID(), ci.getWeight());
//				workingNodeLoc.removeConnection(NodeToConnect.getID(), n.getID());
//				jj=jj+2;
//			}
//		}
//		
//	//	System.out.println("Random :"+workingNodeLoc.m_random+" closest :"+workingNodeLoc.m_closest);
//		if(canditate.getTotalEdges() != GraphConstants.MAX_CONNECTIONS)
//			throw new Exception("Total edges "+canditate.getTotalEdges()+ " does not equal "+GraphConstants.MAX_CONNECTIONS);
//		
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
