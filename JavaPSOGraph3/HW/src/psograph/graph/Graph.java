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
package psograph.graph;

import java.io.*; 
import java.util.*;

import psograph.measurements.PercentageInLargestCluster;





/**
 * This is the class representing a graph.  This Nodes in the graph are stored in 
 * a linked list implementation to reduce the memory foot print.
 * 
 * Please note this class was implemented to support Di-graphs.  So each node contains
 * information on what node it is connected.  This information is duplicated in the case of
 * undirected graphs.   The plan was to support to Di-Graphs but I ended up not needing to
 * use that portion, so I didn't fully implement it.  However, it shouldnt' be too much work to
 * make this work.  As architecture is there, just the implementation.
 * @author Patrick
 *
 */
public class Graph implements Serializable
{
	//There are two ways of implementing a graph
	// a) adjency matrix of size n by n.   This is only useful in graphs that are dense.
	//    easier to add edges
	// b) adjacency matrix - linked list implementation.  Reduces size in memory.
	
	
	// a) 200 nodes of degree r=3 -> 600 edges . matrix of size 40,000
	// b) 200 nodes of degree r=3 -> 600 edges .  header nodes 200, and each node having
	// 3 edges off of it, totals (3+1)*200 = 800 .   The 3+1 is 3 edges connected to the node
	// plus the header node.
	
	static final long serialVersionUID = 25L;

	//Used to make sure we don't have nodes with the identical coordinates
	private Vector<Location> m_NodeLocationMap;	
	private TreeMap<Integer,Node> m_headerNodesMap;
	private int m_nextNodeId;
	protected Random r;
	
	
	@SuppressWarnings("unused")
	private int m_version = 1;
	
	/**
	 * Copy constructor.  This will make a deep copy of the graph.
	 * @param g
	 */
	public Graph (Graph g)
	{
		r = new Random();
		
		m_headerNodesMap = new TreeMap<Integer,Node>();
		
		Vector<Integer> vi = new Vector<Integer>(g.m_headerNodesMap.keySet());
		for(int i = 0; i < vi.size(); i++)
			m_headerNodesMap.put(vi.get(i),new Node(g.m_headerNodesMap.get(vi.get(i))));

		m_nextNodeId = g.m_nextNodeId;		

		m_NodeLocationMap = g.m_NodeLocationMap;
	}
	
	/**
	 * Make a graph consisting of disconnected nodes.  The
	 * Graph will place the nodes in X-Y space in a unit square.
	 * That is the X,Y values will be in the range of 0.0 - 1.0  
	 * @param nodes the number of nodes 
	 */
	public Graph(int nodes) throws Exception
	{
		r = new Random();
		m_nextNodeId =0;
		m_headerNodesMap = new TreeMap<Integer,Node>();

		m_NodeLocationMap = new Vector<Location>();
		
		for (m_nextNodeId=0; m_nextNodeId < nodes; )
		{
			addNode();
		}
	}	
	
	/**
	 * Add a connection to between two nodes.   The weight will be determined
	 * by Cost function defined for the weight.
	 * If it cannot create the connection an exception will be thrown
	 * @param from_id First Node Id
	 * @param to_id Second Node Id
	 * @throws Exception
	 */
	public void addConnection(int from_id, int to_id) throws Exception
	{
		try
		{
		Node from = m_headerNodesMap.get(from_id);
		Node to = m_headerNodesMap.get(to_id);
		
		if(from == null || to == null)
		{
			throw new Exception("null pointer for node");
		}
		
		NodeLocationCalculator nodeLocationCalculator = new NodeLocationCalculator(this, false);
		//Now that we have a valid node pair, accept into A and return it
		double weight = nodeLocationCalculator.calculateWeight(from.getX(), from.getY(), to.getX(), to.getY());
		
		from.addConnection(to_id, weight);
		to.addConnection(from_id, weight);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}

	}
	
	/**
	 * This will make a connection between two Nodes but the weight will
	 * be specified instead of calculated.   This method should be used with care 
	 * as the weights will be calculate outside the Cost Function.
	 * 
	 *  THis is used in recording Hops in the SPL measurements
	 * 
	 * @param from_id
	 * @param to_id
	 * @param weight
	 * @throws Exception
	 */
	public void addConnection(int from_id, int to_id, double weight) throws Exception
	{
		Node from = m_headerNodesMap.get(from_id);
		Node to = m_headerNodesMap.get(to_id);
		
		if(from == null || to == null)
		{
			throw new Exception("null pointer for node");
		}
		
		from.addConnection(to_id, weight);
		to.addConnection(from_id, weight);
	}
	
	/**
	 * This adds a node to the graph.  THis node will be generated in a random location.
	 * 
	 * This is useful when make a random test node, but when making a new graph.
	 * 
	 * @return the new node id.
	 * 
	 * 
	 */
	public int addNode() throws Exception
	{	
		int retNodeId = m_nextNodeId;
		
		boolean isUnique = false;
		Location loc = null;
		while(!isUnique)
		{
			loc = generateLocation();
			
			if(!m_NodeLocationMap.contains(loc)  )
				isUnique = true;		
			else
				System.out.println("not unique");
		}
		//System.out.println("Unique location"+loc.getX()+", "+loc.getY()+")");
		m_NodeLocationMap.add(loc);

		m_headerNodesMap.put(m_nextNodeId, new Node (m_nextNodeId, loc.getX(), loc.getY()));
		
		m_nextNodeId++;
		
		return retNodeId;
	}
	
	/**
	 * Adds a node with the specified id a certain location.  This was used in the case
	 * of setting up SPL calculations, where we didn't want to copy the edges.  Probably should
	 * have made a copy constructor with an option to copy the edges or not.
	 * @param id
	 * @return
	 * @throws Exception 
	 */
	public int addNode(int id, double x, double y) throws Exception
	{	
		int retNodeId = id;
	
		Location loc = new Location(x,y);
			
		if(m_NodeLocationMap.contains(loc)  )
		{
			System.out.println("not unique");
			throw new Exception("Location is not unique");
		}
		//System.out.println("Unique location"+loc.getX()+", "+loc.getY()+")");
		m_NodeLocationMap.add(loc);
		
		m_headerNodesMap.put(id, new Node (id, loc.getX(), loc.getY()));
		
		m_nextNodeId++;
		
		return retNodeId;
	}
	

	
	/**This adds in a new node connected to the list of Nodes.
	 * The value of the weights will be determined by cost function.
	 * 
	 * @param neighbors array of node Ids
	 * @return node id of newly created Node
	 * @throws Exception
	 */
	public int addNode(int[] neighbors) throws Exception
	{	int retNodeId = m_nextNodeId;
	
	    boolean isUnique = false;
	    
		Location loc = null;
		while(!isUnique)
		{
			loc = generateLocation();
			
			if(!m_NodeLocationMap.contains(loc)  )
				isUnique = true;		
			else
				System.out.println("not unique");
		}
	
		m_headerNodesMap.put(m_nextNodeId, new Node (m_nextNodeId, loc.getX(), loc.getY(), neighbors));
		
		
		NodeLocationCalculator nodeLocationCalculator = new NodeLocationCalculator(this, false);

		
		// now do the back connections
		int numNeighbors = neighbors.length;
		for (int i=0; i < numNeighbors; i++)
		{
			Node t = m_headerNodesMap.get(neighbors[i]);
			
			//Now that we have a valid node pair, accept into A and return it
			double weight = nodeLocationCalculator.calculateWeight(t.getX(), t.getY(), loc.getX(), loc.getY());
			
			t.addConnection(m_nextNodeId, weight);
		}
		
		m_nextNodeId++;
		
		return retNodeId;
	}	
	
	/** I believe this method will add in connections to the nodes ID with the specified weights.
	 *  i.e. make a connection to neighbors[i] with weights[i]
	 * 
	 * I have some notes this was for the dag case.  But I never used this except in some testcode.
	 */
	public int addNode(int[] neighbors, double weights[]) throws Exception
	{	
		int retNodeId = m_nextNodeId;
		
	    Location loc = generateLocation();
	
		m_headerNodesMap.put(m_nextNodeId, new Node (m_nextNodeId, loc.getX(), loc.getY(), neighbors));
		
		// now do the back connections
		int numNeighbors = neighbors.length;
		for (int i=0; i < numNeighbors; i++)
		{
			Node t = m_headerNodesMap.get(neighbors[i]);
			

			t.addConnection(m_nextNodeId, weights[i]);
		}
		
		m_nextNodeId++;
		
		return retNodeId;
	}
	
	/**
	 * Adds a node but only Location data only.  This test looks to be only used
	 * in Tests.  I think I wanted this method to create a node at a specified spot
	 * but not to add connections.
	 * @param x
	 * @param y
	 * @return
	 * @throws Exception
	 */
	public int addNodeLocationDataOnly(double x, double y) throws Exception
		{
		
			int retNodeId = m_nextNodeId;
			
			Location loc = new Location(x, y);
			//We want to make sure we don't have a duplicate entry
			
	
			if(m_NodeLocationMap.contains(loc)  )
				throw new Exception("location already is occupied")	;	
	
			//System.out.println("Unique location"+loc.getX()+", "+loc.getY()+")");
			m_NodeLocationMap.add(loc);
	
			m_headerNodesMap.put(m_nextNodeId, new Node (m_nextNodeId, loc.getX(), loc.getY()));
			
			m_nextNodeId++;
			
			return retNodeId;
			
		}
	/**
	 * This generates a random location inside the unit square.
	 * This is not checked if it is unique location.
	 * 
	 * Due to how nextDouble works, the range of location
	 * will be 0.0 (inclusive) to 1.0 (exclusive).
	 * 
	 */
	protected Location generateLocation() throws Exception
	{
		double xi = r.nextDouble();
		double yi = r.nextDouble();
		
		double x = xi;
		double y = yi;
		
		Location loc = new Location(x, y);
			
		return loc;	
	}	

	/**
	 * This generates an adjecency matrix.  Useful when making
	 * Matlab files from graphs.  I believe this version 
	 * has to return an Integer explictily because 
	 * we stream this out.  But this can probably be refactored.
	 * @return
	 */
	public Integer [][] getAdjecencyGraph()
	{
		
		
		Vector<Node> nodes = new Vector<Node>(m_headerNodesMap.values() );
		int n = nodes.size();
		
		if (n == 0)
			return null;

		Integer [][]adj = new Integer[n][n];
		
		int i,j;
		
		for( i =0; i < n; i++)
		{
			for(j=0; j < n ; j++)
				adj[i][j] =0;
		}

		for( i =0; i < n; i++)
		{
			Node node = nodes.get(i);
			if(node.getNeighbors() == null)
				continue;
			
			Vector<Integer> neighbors = new Vector<Integer>(node.getNeighbors().keySet());

			int n_neighbors = neighbors.size();
			for(j=0; j < n_neighbors; j++)
			{
				int tt = neighbors.get(j);
				Node temp = m_headerNodesMap.get(tt);
				int place = nodes.indexOf(temp);
				adj[place][i] = 1;
				adj[i][place] = 1;
			}		
		}
		return adj;
		
	}
	/**
	 * This generates an adjecency matrix.  Useful when making
	 * Matlab files from graphs.  This can probably be refactored.
	 */
public int [][] getAdjecencyGraph2()
{
	
	
	Vector<Node> nodes = new Vector<Node>(m_headerNodesMap.values() );
	int n = nodes.size();
	
	if (n == 0)
		return null;

	int [][]adj = new int[n][n];
	
	int i,j;
	
	for( i =0; i < n; i++)
	{
		for(j=0; j < n ; j++)
			adj[i][j] =0;
	}

	for( i =0; i < n; i++)
	{
		Node node = nodes.get(i);
		if(node.getNeighbors() == null)
			continue;
		
		Vector<Integer> neighbors = new Vector<Integer>(node.getNeighbors().keySet());

		int n_neighbors = neighbors.size();
		for(j=0; j < n_neighbors; j++)
		{
			int tt = neighbors.get(j);
			Node temp = m_headerNodesMap.get(tt);
			int place = nodes.indexOf(temp);
			adj[place][i] = 1;
			adj[i][place] = 1;
		}		
	}
	return adj;
	
}

/**
 * Adjecency matrix with weights as opposed is connected or not.
 * @return
 */
	public double [][] getAdjecencyGraphWeights2()
	{
		
		
		Vector<Node> nodes = new Vector<Node>(m_headerNodesMap.values() );
		int n = nodes.size();
		
		if (n == 0)
			return null;

		double [][]adj = new double[n][n];
		
		int i,j;
		
		for( i =0; i < n; i++)
		{
			for(j=0; j < n ; j++)
				adj[i][j] =0;
		}

		for( i =0; i < n; i++)
		{
			Node node = nodes.get(i);
			if(node.getNeighbors() == null)
				continue;
			
			Vector<Integer> neighbors = new Vector<Integer>(node.getNeighbors().keySet());

			int n_neighbors = neighbors.size();
			for(j=0; j < n_neighbors; j++)
			{
				int tt = neighbors.get(j);
				Node temp = m_headerNodesMap.get(tt);
				
				double weight = temp.getConnectionInfo(node).getWeight();
				
				int place = nodes.indexOf(temp);
				adj[place][i] = weight;
				adj[i][place] = weight;
			}		
		}
		return adj;
		
	}
	
	/**
	 * Sums all the weights.   Please note this method assumes non-dag.
	 * 
	 * As we count all edges weigts then divide by two.  If we supported DAG
	 * we would need to stop divide by 2.
	 * @return
	 */
	public double getCost()
	{
		return SumAllWeights()/2;
	}
	
	/**
	 * NOt sure why this is advanced :)
	 * But this returns a Map with the key being the node degree.  And the value
	 * being a vector of the nodes IDs with that degree.
	 * @return
	 */
	public TreeMap<Integer,Vector<Integer>> getDegreeDistributionAdvanced()
	{
		TreeMap<Integer,Vector<Integer>> degreeDist = new TreeMap<Integer,Vector<Integer>>();
		
		if( m_headerNodesMap != null && m_headerNodesMap.size() != 0)
		{
			Collection<Node> c = m_headerNodesMap.values();
			Vector<Node> v = new Vector<Node>(c);
					
			for(int i = 0; i < v.size(); i++)
			{
				Node n = v.get(i);
				int deg = n.getDegree();
				
				if(degreeDist.get(deg) == null)
				{
					degreeDist.put(deg, new Vector<Integer>());
				}
				degreeDist.get(deg).add(n.getID());
			}
		}
	
		return degreeDist;
		
	}
	
	/**
	 * Return the Node Map.  This is Map of the Nodes that are keyed by ID number.
	 * @return
	 */
	public TreeMap<Integer,Node> getHeaderNodesMap()
	{
		return m_headerNodesMap;
	}
	/**
	 * Returns the Inverse of the Weigth Distribution.  THis returns a Map with the key being
	 * the inverse weight.  And the value being the number of edges with that inverse weight.
	 * @return
	 * @throws Exception
	 */
	public TreeMap<Double, Integer> getInverseWeightDistribution()  throws Exception
	{
		TreeMap<Double, Integer> result = getWeightDistribution(true);
		
		return result;
	}		
	
	/**Return the number of Nodes in the graph
	 * 
	 * @return
	 */
	public int getNumberOfNodes()
	{
		return m_headerNodesMap.size();	
	}	
	
	/**
	 * Get the Node object for a specified ID.
	 * @param id
	 * @return
	 * @throws Exception
	 */
	public Node getNode(int id) throws Exception
	{
		return  m_headerNodesMap.get(id);
		
	}	

	/**
	 * Returns the total number of edges.  Please note if we support
	 * DAG this will need modification as we add up all edges and divide by two.
	 * @return
	 */
	public int  getTotalEdges()
	{
		int return_val =0;
		
		if( m_headerNodesMap != null && m_headerNodesMap.size() != 0)
		{
			Collection<Node> c = m_headerNodesMap.values();
			Vector<Node> v = new Vector<Node>(c);
					
			for(int i = 0; i < v.size(); i++)
			{
				Node n = v.get(i);
				return_val =return_val + n.getDegree();
			}
			
		}
		else
			System.out.println("Empty Graph!");
		
		return return_val/2;
		
	}
	
	/**Returns a Map with the key being the weight value.  And the value being the number
	 * of edges with that weigth value.
	 * 
	 * @return
	 * @throws Exception
	 */
	public TreeMap<Double, Integer> getWeightDistribution()  throws Exception
	{
		TreeMap<Double, Integer> result = getWeightDistribution(false);
		
		return result;
	}
	
	/**
	 * Method to get the weight distribution.  With the notable option of sorting by the weight
	 * or the inverse.
	 * @param inverse
	 * @return
	 * @throws Exception
	 */
	private TreeMap<Double,Integer> getWeightDistribution(boolean inverse) throws Exception
	{
		TreeMap<Double, Integer> result = new TreeMap<Double, Integer>();
		
		
		if( m_headerNodesMap != null && m_headerNodesMap.size() != 0)
		{
			Collection<Node> c = m_headerNodesMap.values();
			Vector<Node> v = new Vector<Node>(c);
					
			for(int i = 0; i < v.size(); i++)
			{
				Node n = v.get(i);
				
				TreeMap<Integer,Edge> tci = n.getNeighbors();
				if(tci != null && tci.size() != 0)
				{
					Vector<Integer> vi = new Vector<Integer>(tci.keySet());
					for(int j = 0; j < vi.size();j++)
					{
						double treeKey = tci.get(vi.get(j)).getWeight();
						if(inverse)
						{
							if(treeKey == -1)
								treeKey = 0;
							else
								treeKey = 1/treeKey;
						}
						
						if(result.containsKey(treeKey))
						{
							if(inverse)
							{
								double key = 1/tci.get(vi.get(j)).getWeight();
								if(key == -1)
									key =0;
								
								int prev_value = result.get(key);
								result.put(key, prev_value+1);
							}
							else
							{
								int prev_value = result.get(tci.get(vi.get(j)).getWeight());
								result.put(tci.get(vi.get(j)).getWeight(), prev_value+1);
							}
						}
						else
						{
							if(inverse)
							{
								double key = 1/ tci.get(vi.get(j)).getWeight();
								if(key == -1)
									key =0;
								
								result.put(key, 1);
							}
							else 
							{
								result.put(tci.get(vi.get(j)).getWeight(), 1);
							}
							
						}

					}				
				}
			}
		}
		else
		{
			System.out.println("Empty Graph!");
			return null;
		}
		
		//Now we need to half all the values.	
		Vector<Double> vWeights = new Vector<Double>(result.keySet());		
	
		if(vWeights.size() <= 0)
		{
			throw new Exception ("The weight distribution is empty");
		}
	
		for(int i=0; i < vWeights.size(); i++)
		{
			
			int prev_result = result.get(vWeights.get(i));
			int new_result = prev_result/2;
			result.put(vWeights.get(i), new_result);
		
		}

		return result;
	}
	

	/**
	 * This returns true if every node is connected to the each other.  Otherwise false.
	 * @return
	 * @throws Exception
	 */
	public boolean isFullyConnected() throws Exception
	{
		double LCC =0;
		
		PercentageInLargestCluster percentageInLargestCluster = new PercentageInLargestCluster(this);
		LCC = percentageInLargestCluster.Measure();
		
		
		if(Double.compare(LCC, 1.0)==0)
			return true;
		else
			return false;
		
	}
	
	/**
	 * Query method to see if two nodes have a edge in the graph
	 * @param from From Node
	 * @param to To Node
	 * @return returns true if connected, otherwise false
	 * @throws Exception
	 */
	public boolean isNodeAConnectedToNodeB(Node from, Node to) throws Exception
	{			
		return from.isConnectedTo(to);
	}
	
	/**This is debug routine to dump the graph's node list 
	 * with connections information.  This includes only includes 
	 * which nodes it is connected to.
	 *  
	 */
	public void print()
	{
		if( m_headerNodesMap != null && m_headerNodesMap.size() != 0)
		{
			Collection<Node> c = m_headerNodesMap.values();
			Vector<Node> v = new Vector<Node>(c);
					
			for(int i = 0; i < v.size(); i++)
			{
				Node n = v.get(i);
				System.out.print("Node "+n.getID()+" neighbors");
				TreeMap<Integer,Edge> tci = n.getNeighbors();
				if(tci != null && tci.size() != 0)
				{
					Vector<Integer> vi = new Vector<Integer>(tci.keySet());
					for(int j = 0; j < vi.size();j++)
					{					
						System.out.print(" "+vi.get(j));
					}				
				}
				System.out.println("");
			}		
		}
		else
			System.out.println("Empty Graph!");
	}
	
	/**This is debug routine to dump the graph's node list (with node coordinates)
	 * with connections information.  This includes which nodes it is connected to
	 *  and what the weight of the connection is.
	 */
	public void printWithLocationAndWeights()
	{
		if( m_headerNodesMap != null && m_headerNodesMap.size() != 0)
		{
			Collection<Node> c = m_headerNodesMap.values();
			Vector<Node> v = new Vector<Node>(c);
					
			for(int i = 0; i < v.size(); i++)
			{
				Node n = v.get(i);
				System.out.print("Node "+n.getID()+" loc ("+n.getX()+","+n.getY()+") neighbors");
				TreeMap<Integer,Edge> tci = n.getNeighbors();
				if(tci != null && tci.size() != 0)
				{
					Vector<Integer> vi = new Vector<Integer>(tci.keySet());
					for(int j = 0; j < vi.size();j++)
					{
						System.out.print(" "+vi.get(j)+":"+tci.get(vi.get(j)).getWeight());
					}				
				}
				System.out.println("");
			}
			
		}
		else
		{
			System.out.println("Empty Graph!");
		}
	}
	

	/**This is debug routine to dump the graph's node list
	 * with connections information.  This includes which nodes it is connected to
	 *  and what the weight of the connection is.
	 */
	public void printWithWeights()
	{
		if( m_headerNodesMap != null && m_headerNodesMap.size() != 0)
		{
			Collection<Node> c = m_headerNodesMap.values();
			Vector<Node> v = new Vector<Node>(c);
					
			for(int i = 0; i < v.size(); i++)
			{
				Node n = v.get(i);
				System.out.print("Node "+n.getID()+" neighbors");
				TreeMap<Integer,Edge> tci = n.getNeighbors();
				if(tci != null && tci.size() != 0)
				{
					Vector<Integer> vi = new Vector<Integer>(tci.keySet());
					for(int j = 0; j < vi.size();j++)
					{
						System.out.print(" "+vi.get(j)+":"+tci.get(vi.get(j)).getWeight());
					}				
				}
				System.out.println("");
			}
			
		}
		else
			System.out.println("Empty Graph!");
	}
	

	/**
	 * This method will sever the connection between two nodes.
	 * This method will throw an exception if it is not able to 
	 * remove the connections.  
	 * @param from_id The first node
	 * @param to_id the second node
	 * @throws Exception
	 */
	public void removeConnection(int from_id, int to_id) throws Exception
	{
		Node from = m_headerNodesMap.get(from_id);
		Node to = m_headerNodesMap.get(to_id);
		from.removeConnection(to_id);
		to.removeConnection(from_id);
	}
	
	/**
	 * Removes a Node with the specified ID and removes all connections to and from 
	 * that node.
	 * @param id
	 */
	public void removeNode(int id)
	{
		Node t = m_headerNodesMap.get(id);	
		TreeMap<Integer,Edge> tci = t.getNeighbors();
		
		//for each neighbor, pull out our reference
		if(tci != null && tci.size() != 0)
		{   
			Vector<Integer> vkeys = new Vector<Integer>(tci.keySet());
			
			for(int j = vkeys.size() -1 ; j >= 0 ;j=j-1)
			{				
				Node tt = m_headerNodesMap.get(vkeys.get(j));
				TreeMap<Integer,Edge> tci2 = tt.getNeighbors();
				if(tci2 != null && tci2.size() != 0)
				{   
					Edge ci2 = tci2.get(id);
					if(ci2 != null)
						tci2.remove(id);
				}
			}
		}
		m_headerNodesMap.remove(id);
	}
	
	/**
	 * Sums all the weights.  Please note since I support DAG, in the non-DAG case we
	 * have to divide by 2 this value.
	 * @return
	 */
	public double SumAllWeights()
	{
		double value = 0;
		
		if( m_headerNodesMap != null && m_headerNodesMap.size() != 0)
		{
			Collection<Node> c = m_headerNodesMap.values();
			Vector<Node> v = new Vector<Node>(c);
					
			for(int i = 0; i < v.size(); i++)
			{
				Node n = v.get(i);
				TreeMap<Integer,Edge> tci = n.getNeighbors();
				if(tci != null && tci.size() != 0)
				{
					Vector<Integer> vi = new Vector<Integer>(tci.keySet());
					for(int j = 0; j < vi.size();j++)
					{
						
						if(tci.get(vi.get(j)).getWeight() == -1)
						{
							value += 0;
						}
						else
						{
							value += tci.get(vi.get(j)).getWeight();
						}

					}				
				}
			}
			
		}
		return value;
	}
	/**
	 * Sums all the inverse weights.  Please note since I support DAG, in the non-DAG case we
	 * have to divide by 2 this value.
	 * @return
	 */
	public double sumInverseWeights()
	{
		double value = 0;
		
		if( m_headerNodesMap != null && m_headerNodesMap.size() != 0)
		{
			Collection<Node> c = m_headerNodesMap.values();
			Vector<Node> v = new Vector<Node>(c);
					
			for(int i = 0; i < v.size(); i++)
			{
				Node n = v.get(i);
				TreeMap<Integer,Edge> tci = n.getNeighbors();
				if(tci != null && tci.size() != 0)
				{
					Vector<Integer> vi = new Vector<Integer>(tci.keySet());
					for(int j = 0; j < vi.size();j++)
					{

						if(tci.get(vi.get(j)).getWeight() == -1)
						{
							value += 0;
						}
						else
						{
							value += 1/(tci.get(vi.get(j)).getWeight());
						}
						
					}				
				}
			}
			
		}
		return value;
	}
	
}
