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
package psograph.graph;

import java.io.*; 
import java.util.*;

import psograph.graph.GraphConstants;
import psograph.measurements.*;


/**
 * Calculate graph extends the Graph object.
 * This extends it the Graph as it stores the calculated values.
 * For example AISPL, ASPL, and so no.
 * This class can then be stream out to a binary object.  And then read in later and querried.
 * THis way avoiding having to re-calculate these unless really needed.
 * @author Patrick
 *
 */
public class CalculatedGraph extends Graph implements Serializable
{

	static final long serialVersionUID = 2285L;
		
	
	//Be careful removing these, as this can prevent you from streaming in previous CalculatedGraph
	// I did explore it, but I believe you can have the version key off the serialVersionUID
	// and have it deal gracefully with this situation.
	


	
	//This guy takes up a lot of space, so we should not store him
	// in the object written to disk.
	private Graph SPLGraphs[];
	
	private double m_fitness_value;
	private double m_randomDiameter;
	private double m_directDiameter;
	private double m_random_lcc;
	private double m_direct_lcc;
	private double m_AISPL;
	private double m_ASPL;
	private double m_costBasis;
	private double m_diameter;
	private double m_LCC;
	private double m_ClusteringCoefficient;
	private double m_avgNodeDegree;
	private double m_avgLinkCost;
	private double m_avgNodeLength;
	private double m_assortativity;
	/**
	 * @deprecated I had a bug in the SPL calculation, but I wanted to make
	 * sure I didn't use the bugged calucation, so I put in new member
	 * variable to hold that data.  However, since we stream these in 
	 * I need to keep this around, but do not use this variable.
	 */
	
	@SuppressWarnings("unused")
	private TreeMap<Double, Integer> m_distributionSPL;
	/**
	 * @deprecated
	 */
	double m_costOfDirectConnectivity;
	/**
	@deprecated
	*/
	@SuppressWarnings("unused")
	private double m_clumpinessOfLayout;
	
	private TreeMap<Integer, Integer> m_distributionSPL2;
	private TreeMap<Double, Integer> m_distributionISPL;
	private TreeMap<Double, Integer> m_nodeTotalEdgeLength;
	private TreeMap<Integer, Double> m_nodeTotalEdgeLengthById;
	private TreeMap<Integer, Double> m_nodeTotalEdgeCostById;
	
	private TreeMap<Double,Integer> m_nodeTotalEdgeCost;
	private TreeMap<Double, Integer> m_nodeMaximumEdgeLength;
	private TreeMap<Integer, Double> m_nodeMaximumEdgeLengthById;
	private TreeMap<Double, Integer> m_nodeMaximumEdgeCost;
	private TreeMap<Integer, Double> m_nodeMaximumEdgeCostById;
	private TreeMap<Integer, Integer> m_connectedNodesWithinRadius;
	private TreeMap<Integer, Integer> m_connectedNodesWithinRadiusById;
	private TreeMap<Integer, Integer> m_nodesWithinRadius;
	private TreeMap<Integer, Integer> m_nodesWithinRadiusById;
	private TreeMap<Double, Integer> m_costOfDirectConnectivity2;
	
	private TreeMap<Double, Integer> m_meanDegreeOfConnectedNeighbors;
	private TreeMap<Integer, Double> m_meanDegreeOfConnectedNeighborsById;
	private TreeMap<Double, Integer> m_betweeness;
	private TreeMap<Integer, Double> m_betweenessById;
	private TreeMap<Double, Integer> m_closeness;
	private TreeMap<Integer, Double> m_closenessById;
	


	

/**
 * Copy constructor
 * @param g
 */
	public CalculatedGraph(CalculatedGraph g) {
		super(g);


		this.SPLGraphs = g.SPLGraphs;
		this.m_fitness_value = g.m_fitness_value;
		this.m_randomDiameter = g.m_randomDiameter;
		this.m_directDiameter = g.m_directDiameter;
		this.m_random_lcc = g.m_random_lcc;
		this.m_direct_lcc = g.m_direct_lcc;
		this.m_AISPL = g.m_AISPL;
		this.m_ASPL = g.m_ASPL;
		this.m_costBasis = g.m_costBasis;
		this.m_diameter = g.m_diameter;
		this.m_LCC = g.m_LCC;
		this.m_ClusteringCoefficient = g.m_ClusteringCoefficient;
		this.m_avgNodeDegree = g.m_avgNodeDegree;
		this.m_avgLinkCost = g.m_avgLinkCost;
		this.m_avgNodeLength = g.m_avgNodeLength;
		this.m_assortativity = g.m_assortativity;
		this.m_distributionSPL2 = g.m_distributionSPL2;
		this.m_distributionISPL = g.m_distributionISPL;
		this.m_nodeTotalEdgeLength = g.m_nodeTotalEdgeLength;
		this.m_nodeTotalEdgeLengthById = g.m_nodeTotalEdgeLengthById;
		this.m_nodeTotalEdgeCost = g.m_nodeTotalEdgeCost;
		this.m_nodeMaximumEdgeLength = g.m_nodeMaximumEdgeLength;
		this.m_nodeMaximumEdgeCost = g.m_nodeMaximumEdgeCost;
		this.m_connectedNodesWithinRadius = g.m_connectedNodesWithinRadius;
		this.m_nodesWithinRadius = g.m_nodesWithinRadius;
		this.m_costOfDirectConnectivity2 = g.m_costOfDirectConnectivity2;
	}
	
	/**
	 * generate a CalculatedGraph from a Graph object
	 * @param graph
	 */
	public CalculatedGraph(Graph graph)
	{
		super(graph);
	}

	/**
	 * Returns the AISPL value
	 * @return
	 */
	public double getAISPL()
	{
		return m_AISPL;
	}
	/**
	 * Returns the ASPL value
	 * @return
	 */
	public double getASPL()
	{
		return m_ASPL;
	}
	/**
	 * Returns the Assortativity value
	 * @return
	 */
	public double getAssortativity() {
		return m_assortativity;
	}
	/**
	 * Return average Link cost/weight
	 * @return
	 */
	public double getAvgLinkCost() {
		return m_avgLinkCost;
	}
	/**
	 * returns average link length
	 * @return
	 */
	public double getAvgLinkLength() {
		return m_avgNodeLength;
	}
	/**
	 * Returns average node degree
	 * @return
	 */
	public double getAvgNodeDegree() {
		return m_avgNodeDegree;
	}
	/**
	 * Betweeness is the property that measures how often a node is on the shortest path between two nodes.   
	 * Where the key is the betweeness value, and the value is the node ID.
	 * 
	 * @return
	 * @throws Exception
	 */
	public TreeMap<Double, Integer> getBetweeness() throws Exception
	{
		return m_betweeness;
	}
	/**
	 * Betweeness is the property that measures how often a node is on the shortest path between two nodes.   
	 * Where the key is the betweeness value, and the value is the node ID.
	 * 
	 * I forget what Complete Layout means.  I think I used this for testing, and this connects the graph's nodes
	 * to every other node, and caculate the between value, which should be zero
	 * @return
	 * @throws Exception
	 */
	public TreeMap<Double, Integer> getBetweenessInCompleteLayout() throws Exception
	{
		
		NodeLocationCalculator nodeLocationCalculator = 
			new NodeLocationCalculator(this);
		nodeLocationCalculator.calculateResults();
		
		SPL SPL_measure = new SPL(nodeLocationCalculator);
		Graph graphs[] = SPL_measure.measure(); 
		
		if(graphs == null)
			return null;
		
		TreeMap<Double, Integer> result = new TreeMap<Double,Integer>();
		
		for(int iGraph = 0; iGraph < this.getNumberOfNodes(); iGraph++)
		{
			double splSum =0;
			
		   // System.out.println("");
			for(int i =0; i < graphs[iGraph].getNumberOfNodes() ; i++)
			{
				Node n = graphs[iGraph].getNode(i);
	//			Vector<Path> vPaths = n.getPaths();
				
				splSum += n.getSPLength();

//				System.out.println("Number of paths from "+ iGraph +" to "+n.getID()+" is " + vPaths.size());
//
//				for(int j = 0; j < vPaths.size(); j++)
//				{
//					Path p = vPaths.get(j);
//
//					System.out.print("start of path "+ p.getStart().getID()+" ");
//					for(int k =0; k < p.getLength(); k++)
//					{
//						Node n2 = p.getPath().get(k);
//						int id = n2.getID();
//						System.out.print(" " + id);
//					}
//					System.out.println(": Path Length " + p.getLength());					
//				}

			}
			
			
			
			double value = ((double)(this.getNumberOfNodes() -1))/splSum;
			
			int prevValue = 0;
			
			if(result.containsKey(value))
			{
				prevValue = result.get(value);
			}
			
			result.put(value, prevValue +1);
		
		}

		return result;
	}
	/**
	 * Closeness is a measurement of the mean hop distance from all other nodes.
	 * Where the key is the closeness value, and the value is the node ID.
	 * @return
	 * @throws Exception
	 */
	public TreeMap<Double, Integer> getCloseness() throws Exception
	{
		return m_closeness;

	}	
	/**
	 * Ideally, the nodes in each layout should be uniformly distributed, and this is what clumpiness 
	 * tries to measure. The layout area is divided into identical square cells and the distribution of 
	 * the count of nodes in each cell gives a measure of clumpiness.
	 * Where the key is the clumpiness (i.e. number of nodes in cell) value, and the value is the number of cells 
	 * with that clumpiness number
	 * @return
	 * @throws Exception
	 */
	public TreeMap<Integer, Integer> getClumpinessOfLayout() 
	{
		 TreeMap<Integer, Integer> result = new TreeMap<Integer, Integer>();
		 
		 int buckets [][] =  
			 new int[GraphConstants.ClumpinessDivisor][GraphConstants.ClumpinessDivisor];
		 
		 Vector<Node> vNodes = new Vector<Node>(getHeaderNodesMap().values());
		 
		 for(int i =0; i < vNodes.size(); i++)
		 {
			 double x = vNodes.get(i).getX();
			 double y = vNodes.get(i).getY();
			 
			 int x_bucket = (int)Math.rint(x/.1);
			 int y_bucket =  (int)Math.rint(y/.1);
			 
			 if(x_bucket == 10)
				 x_bucket = 9;
			 
			 if(y_bucket == 10)
				 y_bucket = 9;
			 buckets[x_bucket][y_bucket]++;
			 
			// System.out.println("Node ("+x+", "+y+")  in bucket["+x_bucket+"]["+y_bucket+"]");
		 }
		 
		 for(int i = 0; i < GraphConstants.ClumpinessDivisor;i++)
		 {
			 for(int j =0; j < GraphConstants.ClumpinessDivisor; j++)
			 {
				 int value = 0;
				 if(result.containsKey(buckets[i][j]))
				 {
					 value = result.get(buckets[i][j]);
				 }
				 
				 result.put(buckets[i][j], value+1);
			 }
		 }
		
		
		
		return result;
	}
	/**
	 * Returns the clustering coefficient of the graph.
	 * @return
	 */
	public double getClusteringCoefficient() 
	{
		return m_ClusteringCoefficient;
	}

	/**
	 * This gives back a Map where the key is node ID and the value is number of 
	 * nodes that are connected to with a weight edge less than or equal to
	 * GraphConstants.MeasurementR.
	 */
	public TreeMap<Integer, Integer> getConnectedNodesWithinRadius() {
		return m_connectedNodesWithinRadius;
	}
	
	/**
	 * Returns the value of the Cost Basis.
	 * 
	 * To generate a Cost Basis  we use a heuristic method to generate a sufficiently expensive network for the given node
	 * layout, and then use the total cost of this network as the cost basis to scale the costs of the networks being 
	 * evaluated. 
	 * @return
	 */
	public double getCostBasis()
	{
		return m_costBasis;
	}
	
	/**
	 *  Returns a Map to represent the distribution of weights 
	 *  if the graph was fully connected (a complete graph).
	 * @return
	 */
	public TreeMap<Double, Integer> getCostOfDirectConnectivityDist() {
		return m_costOfDirectConnectivity2;
	}
	
	/**
	 * This returns the cost of the all the edges if the graph was a
	 * complete graph.
	 * @return
	 * @throws Exception
	 */
	public double getCostOfDirectConnectivity() throws Exception
	{
		NodeLocationCalculator calc = new  NodeLocationCalculator(this);
		calc.calculateResults();
		
		return calc.getCost();
	}
	
	/**
	 * Returns the Diamter of the graph.
	 * The diameter of a graph is the maximum eccentricity of any vertex in the graph. 
	 * That is, it is the greatest distance between any pair of vertices. To find the 
	 * diameter of a graph, first find the shortest path between each pair of vertices. 
	 * The greatest length of any of these paths is the diameter of the graph.
	 * @return
	 */
	public double getDiameter() {
		return m_diameter;
	}

	/**
	 * Average Diameter after graph has been directed Attacked (highest degree nodes have
	 * been removed).
	 * @return
	 */
	public double getDirectDiameter() {
		return m_directDiameter;
	}
	
	/**
	 * Average Largest Connected Component after graph has been directed Attacked (highest degree nodes have
	 * been removed).
	 * @return
	 */
	public double getDirectLCC()
	{return m_direct_lcc;}
	
	/**
	 * Returns a Map of the Invese SPL.  Where the key is the ISPL, and the value 
	 * is the number with that ISPL. 
	 */
	public TreeMap<Double, Integer> getDistributionISPL() 
	{
		return m_distributionISPL;
	}
	
	/**
	 * Returns a map representing the distribution of Link Costs.  The key is the weight cost
	 * and the value is the number with that cost.
	 * @return
	 * @throws Exception
	 */
	public TreeMap<Double, Integer> getDistributionLinkCosts() throws Exception
	{
		return getWeightDistribution();
	}
	
	/**
	 * Returns a map representing the distribution of Link length.  The key is the edge's length                                                                                                                                                                                                                                                                                                                                               
	 * and the value is the number with that cost.
	 * @return
	 * @throws Exception
	 */
	public TreeMap<Double, Integer> getDistributionLinkLengths() throws Exception
	{
		TreeMap<Double, Integer> result = new TreeMap<Double, Integer>();
		TreeMap<Double, Integer> temp = getWeightDistribution();
		
		//Need to calclate the data to get the length
		
		Vector<Double> vCosts = new Vector<Double>(temp.keySet());		
		
		if(vCosts.size() <= 0)
		{
			throw new Exception("The weight distribution is empty");
		}

		for(int i=0; i < vCosts.size(); i++)
		{
			double cost = vCosts.get(i);
			double distance = NonLinearCostFunction.getDistance(cost);
			int count = temp.get(cost);
			
			result.put(distance, count);

		}
		
		return result;
	}
	
	/**
	 * Retuns the SPL distribution.   The key is the number of hops, and the value is the 
	 * number of SPL with that number of hops.
	 * @return
	 */
	public TreeMap<Integer, Integer> getDistributionSPL() {
		return m_distributionSPL2;
	}
	
	/**
	 * Returns the fitness of the graph will full connectivity.
	 * @return
	 * @throws Exception
	 */
	public double getFitnessOfSeedWithFullConnectivity() throws Exception
	{
		double value  = -1;
		
		NodeLocationCalculator temp = new NodeLocationCalculator(this);
		temp.calculateResults();
		
		//System.out.println(" Big cost is "+ temp.getCost());
		
		CalculatedGraph tt = new CalculatedGraph(temp);
		tt.setCostBasis(getCostBasis());
		
		tt.UpdateCalcuations();

		value = tt.getFitnessValue();
		
		return value;
	}
	
	/**
	 * Fitness value of the graph.
	 * @return
	 */
	public double getFitnessValue()
	{
		return m_fitness_value;
	}
	
	/**Returns the Hubbiness Cuttoff.  This is the mean node degree plus three times the standard 
	 * deviation of the node degree distribution.
	 * 
	 * @return
	 */
	public double getHubinessCutoff()
	{
		
		double k = 3.0;
		double stdDev = this.getNodeDegreeStdDeviation();
		
		double mean = this.getAvgNodeDegree();
		
		
		double limit = mean + ( k * stdDev);
		
		return limit;
	}
	/**Returns the Hubbiness.  This measure is of our own devising that reports how many 
	 * nodes have a degree greater than the mean node degree plus three times the standard 
	 * deviation of the node degree distribution.
	 * 
	 * @return
	 */
	public double getHubbiness() 
	{
		int hubbiness = 0;
		double k = 3.0;
		double stdDev = this.getNodeDegreeStdDeviation();
		
		double mean = this.getAvgNodeDegree();
		
		
		double limit = mean + ( k * stdDev);
		
		TreeMap<Integer, Integer> nodeDist = getNodeDistribution();
		
	//	int iLimit = (int)Math.rint(limit);
		
	//	System.out.println("iLimit is "+iLimit);
		
	
		Vector<Integer> vDegree = new Vector<Integer>(nodeDist.keySet());
		
		for(int i=0; i < vDegree.size(); i++)
		{
			int cmpVal = vDegree.get(i);
			if(Double.compare( (double)cmpVal, limit) > 0)
			{
				hubbiness += nodeDist.get(vDegree.get(i));
				
			}
			
		}
		
		return hubbiness;
	}
	/**
	 * Returns the percentage of nodes in the Largest Connected Component.
	 * @return
	 */
	public double getLCC() {
		return m_LCC;
	}
	
	/**
	 * Returns a distribution of Edges Longer than GraphConstants.MeasurementR.
	 * The key is the distance, and the value is the number of edges with that value.
	 * @return
	 * @throws Exception
	 */
	public TreeMap<Double, Integer> getLinksLongerThan_r() throws Exception
	{
		TreeMap<Double, Integer> result = new TreeMap<Double, Integer>();
		TreeMap<Double, Integer> temp = getWeightDistribution();
		
		//Need to calclate the data to get the length
		
		Vector<Double> vCosts = new Vector<Double>(temp.keySet());		
		
		if(vCosts.size() <= 0)
		{
			throw new Exception("The weight distribution is empty");
		}

		for(int i=0; i < vCosts.size(); i++)
		{
			double cost = vCosts.get(i);
			double distance = NonLinearCostFunction.getDistance(cost);
			int count = temp.get(cost);
			
			if(Double.compare(distance , GraphConstants.MeasurementR) > 0.0)
			{
				result.put(distance, count);
			}

		}
		
		return result;
	}
	
	/** Returns the distribution of mean degree of connected nodes. 
	 * For each node, the mean degree of all nodes to which it is directly connected to is
	 *  averaged to give this measure. 
	 *  The key is the mean degree, and the value is the number of nodes with that mean degree.
	 * @return
	 * @throws Exception
	 */
	public TreeMap<Double, Integer> getMeanDegreeOfConnectedNeighbors() throws Exception
	{
		return m_meanDegreeOfConnectedNeighbors;
	}


	/**
	 * Returns the Standard Deviation of the node degrees.
	 * @return
	 */
	public double getNodeDegreeStdDeviation() 
	{
		double nodeDegreeStdDeviation =0;
		
		Vector<Node> nodes = new Vector<Node>(getHeaderNodesMap().values());
		
		// calculate standard deviation
		double sum = 0.0;
		for(int i =0; i < nodes.size(); i++)
		{
			sum += Math.pow( nodes.get(i).getDegree(),2) - 
			(2 * getAvgNodeDegree() * nodes.get(i).getDegree()) +
			Math.pow( getAvgNodeDegree(),2);
		}

		sum = sum / (nodes.size() - 1); // divide by n-1
		nodeDegreeStdDeviation = Math.sqrt(sum);
		
		return nodeDegreeStdDeviation;
	}
	
	/**
	 * Returns the Node degree variance of the graph.
	 * @return
	 */
	public double getNodeDegreeVariance() 
	{
		double nodeDegreeVariance = Math.pow(getNodeDegreeStdDeviation(), 2);
		
		return nodeDegreeVariance;
	}

	/**
	 * Returns a Map with the key being the Node ID, and the value
	 * being the node degree.
	 * @return
	 */
	public TreeMap<Integer, Integer> getNodeDistributionId() 
	{
		//degree  # of degree
		TreeMap<Integer,Integer> nodeDistributionById = new TreeMap<Integer, Integer>();
		
		Vector<Node> vec = new Vector<Node>(this.getHeaderNodesMap().values());
		
		for(int i = 0; i < vec.size(); i++)
		{
			Node node = vec.get(i);
						
			nodeDistributionById.put(node.getID(), node.getDegree());
		}
		
		return nodeDistributionById;
	}

	/**
	 * Returns a Map with the key being the node degree, and the value
	 * being the number of nodes with that degree.
	 * @return
	 */
	public TreeMap<Integer, Integer> getNodeDistribution() 
	{
		//degree  # of degree
		TreeMap<Integer,Integer> nodeDistribution = new TreeMap<Integer, Integer>();
		
		
		Vector<Node> vec = new Vector<Node>(this.getHeaderNodesMap().values());
		
		
		for(int i = 0; i < vec.size(); i++)
		{
			
			Node node = vec.get(i);
			
			
			if(nodeDistribution.containsKey(node.getDegree()))
			{
				//need to remove previous key value and update it
				int prevValue = nodeDistribution.get(node.getDegree());			
				nodeDistribution.put(node.getDegree(), prevValue + 1);
			}
			else
			{
				//First entry
				nodeDistribution.put(node.getDegree(), 1);
			}
			
		}
		
		return nodeDistribution;
	}
	

	/**
	 * Returns a map, with the key being a Nodes maximumly edge cost,
	 * and the value being the number of nodes with that maximal edge cost.
	 * @return
	 */
	public TreeMap<Double, Integer> getNodeMaximumEdgeCost() {
		return m_nodeMaximumEdgeCost;
	}
	/**
	 * Returns a map, with the key being a Nodes maximumly edge length,
	 * and the value being the number of nodes with that maximal edge length.
	 * @return
	 */	
	public TreeMap<Double, Integer> getNodeMaximumEdgeLength() {
		return m_nodeMaximumEdgeLength;
	}
	
	/**
	 * Returns a Map with the key being the number of Nodes within in
	 * GraphConstants.MeasurementR, and the value being the number of nodes
	 * with that many nodes within GraphConstants.MeasurementR.
	 * @return
	 */
	public TreeMap<Integer, Integer> getNodesWithinRadius() {
		return m_nodesWithinRadius;
	}
	/**
	 * Returns a map, with the key being a Nodes total edge costs,
	 * and the value being the number of nodes with that maximal edge costs.
	 * @return
	 */	
	public TreeMap<Double, Integer> getNodeTotalEdgeCost() 
	{
		return m_nodeTotalEdgeCost;
	}
	/**
	 * Returns a map, with the key being a Nodes total edge length,
	 * and the value being the number of nodes with that maximal edge length.
	 * @return
	 */	
	public TreeMap<Double, Integer> getNodeTotalEdgeLength() {
		return m_nodeTotalEdgeLength;
	}
	

	/** Returns the percentage of nodes directly connected to a Hub Node.
	 * A hub node is one with degree equal to or greater than 
	 * three times the Standard deviation of all the nodes.
	 * @return
	 * @throws Exception
	 */
	public double getPercentNodesConnectedToHub() throws Exception
	{
		double perConnectedHubs = 0;
		double k = 3.0;
		double stdDev = this.getNodeDegreeStdDeviation();
		
		double limit = k * stdDev;
	
		Vector<Node> vNodes = new Vector<Node>(getHeaderNodesMap().values());
		TreeMap<Integer,Boolean> Hubs = new TreeMap<Integer,Boolean>();
		
		for(int i=0; i < vNodes.size(); i++)
		{
			Node n = vNodes.get(i);
			if(Double.compare( (double)n.getDegree(), limit) > 0)
			{
				Hubs.put(n.getID(),true);
				
				Vector<Integer> t =
					new Vector<Integer>(n.getNeighbors().keySet());
				
				for(int j = 0; j < t.size(); j++)
				{
					Hubs.put(t.get(j),true);
				}
			}
		}
		
		double size = Hubs.size();
		
		perConnectedHubs = size/getNumberOfNodes();

		return perConnectedHubs;
	}
	
	/**
	 * average diameter after random attack.
	 * @return
	 */
	public double getRandomDiameter() {
		return m_randomDiameter;
	}
	
	/**
	 * average LCC after random attack.
	 * @return
	 */
	public double getRandomLCC()
	{return m_random_lcc;}
	
	/**
	 * Mutator for cost basis
	 * @param acb
	 */
	public void setCostBasis(double acb)
	{
		m_costBasis = acb;
	}
	
	/**
	 * Update all calculations at once.  This can take a long time.
	 * 
	 * @throws Exception
	 */
	public void UpdateCalcuations() throws Exception
	{
		//This measurement needs to be first as it depended upon by other measurements
		generateSPLGraphs();
		
		RandomDamage randomDamage = new RandomDamage(this,.15);
		m_random_lcc = randomDamage.Measure();
		m_randomDiameter = randomDamage.getRandomDiameter();
		
		DirectedDamage directedDamage = new DirectedDamage(this,.15);
		m_direct_lcc = directedDamage.Measure();
		m_directDiameter = directedDamage.getDirectDiamter();
	  
		AISPL AISPL_measure = new AISPL(this);
		m_AISPL = AISPL_measure.Measure();
		m_distributionISPL = AISPL_measure.getISPLDistribution();
		
		ASPL ASPL_measure = new ASPL(this);
		m_ASPL = ASPL_measure.Measure();
		m_distributionSPL2 = ASPL_measure.getSPLDistribution();
		
		PercentageInLargestCluster percentageInLargestCluster = 
			new PercentageInLargestCluster(this);
		m_LCC = percentageInLargestCluster.Measure();
		
		Diameter diameter = new Diameter(this, percentageInLargestCluster.LCC());
		m_diameter = diameter.Measure();
		
		ClusteringCoefficient clusteringCoefficient = new ClusteringCoefficient(this);
		m_ClusteringCoefficient = clusteringCoefficient.Measure();
		
		Assortativity assortativity = new Assortativity(this);
		m_assortativity = assortativity.Measure();
		
		AvgNodeDegree avgNodeDegree = new AvgNodeDegree(this);
		m_avgNodeDegree = avgNodeDegree.Measure();
		
		AvgLinkCost avgLinkCost = new AvgLinkCost(this);
		m_avgLinkCost = avgLinkCost.Measure();
		
		AvgLinkLength avgLinkLength = new AvgLinkLength(this);
		m_avgNodeLength = avgLinkLength.Measure();
		
		ConnectedNodesWithinRadius connectedNodesWithinRadius = new ConnectedNodesWithinRadius(this);
		m_connectedNodesWithinRadius = connectedNodesWithinRadius.Measure();
		m_connectedNodesWithinRadiusById = connectedNodesWithinRadius.ById();
		
		NodesWithinRadius nodesWithinRadius = new NodesWithinRadius(this);
		m_nodesWithinRadius = nodesWithinRadius.Measure();		
		m_nodesWithinRadiusById = nodesWithinRadius.ById();
		
		CostOfDirectConnectivity costOfDirectConnectivity = new CostOfDirectConnectivity(this);
		m_costOfDirectConnectivity2 = costOfDirectConnectivity.Measure();
		
		NodeMaximumEdgeCostDistribution nodeMaximumEdgeCostDistribution = 
			new NodeMaximumEdgeCostDistribution(this);
		m_nodeMaximumEdgeCost = nodeMaximumEdgeCostDistribution.Measure();
		m_nodeMaximumEdgeCostById = nodeMaximumEdgeCostDistribution.ById();

		
		
		
		NodeMaximumEdgeLengthDistribution nodeMaximumEdgeLengthDistribution = 
			new NodeMaximumEdgeLengthDistribution(this);
		m_nodeMaximumEdgeLength = nodeMaximumEdgeLengthDistribution.Measure();
		m_nodeMaximumEdgeLengthById = nodeMaximumEdgeLengthDistribution.ById();
		
		NodeTotalEdgeLengthDistribution nodeTotalEdgeLengthDistribution = 
			new NodeTotalEdgeLengthDistribution(this);
		m_nodeTotalEdgeLength = nodeTotalEdgeLengthDistribution.Measure();
		m_nodeTotalEdgeLengthById = nodeTotalEdgeLengthDistribution.ByIdData();
		
		NodeTotalEdgeCostDistribution nodeTotalEdgeCostDistribution = 
			new NodeTotalEdgeCostDistribution(this);
		m_nodeTotalEdgeCost = nodeTotalEdgeCostDistribution.Measure();
		m_nodeTotalEdgeCostById = nodeTotalEdgeCostDistribution.ByIdData();
		
		MeanDegreeOfConnectedNeighbors meanDegreeOfConnectedNeighbors =
			new MeanDegreeOfConnectedNeighbors(this);
		m_meanDegreeOfConnectedNeighbors = meanDegreeOfConnectedNeighbors.Measure();
		m_meanDegreeOfConnectedNeighborsById = meanDegreeOfConnectedNeighbors.ById();
		
		Betweeness betweeness = new Betweeness(this);
		m_betweeness = betweeness.Measure();
		m_betweenessById = betweeness.ById();
		
		Closeness closeness = new Closeness(this);
		m_closeness = closeness.Measure();
		m_closenessById = closeness.ById();
		
		
		
		
		//Not sure if these really should be here
		double cost = getCost();
		@SuppressWarnings("unused")
		double t =  1 - (cost/m_costBasis);		
		//End not sure comment.  I *think* this could be removed.
		
		//This one should always be last because it depends on other values
		m_fitness_value = FitnessFunction.calculateFitness(m_AISPL, (getCost()/m_costBasis), m_LCC) ;
		
		//To reduce space when streaming out, null out array
		SPLGraphs = null;
	}
	/**
	 * Generate SPL graphs needed for other measurements.  
	 * 
	 * 
	 * @throws Exception
	 */
	private void generateSPLGraphs() throws Exception
	{
		SPL SPL_measure = new SPL(this);
		SPLGraphs = SPL_measure.measure(); 
	}

	/**
	 * This will only update the calculations used for the fitness calculation 
	 * in the PSO.  This way we only spend time update the calculations that are
	 * needed instead of doing all the measurements.
	 * @throws Exception
	 */
	public void UpdatePSOCalculations() throws Exception
	{
		AISPL AISPL_measure = new AISPL(this);
		m_AISPL = AISPL_measure.Measure();
		
		PercentageInLargestCluster percentageInLargestCluster = 
			new PercentageInLargestCluster(this);
		m_LCC = percentageInLargestCluster.Measure();
		
		double cost = getCost();
		@SuppressWarnings("unused")
		double t =  1 - (cost/m_costBasis);		
		m_fitness_value = FitnessFunction.calculateFitness(m_AISPL, (getCost()/m_costBasis), m_LCC) ;
		
	}
	
	/**
	 * Returns distribution of Total edge length with the key being the node ID.
	 * @return
	 */
	public TreeMap<Integer, Double> getTotalEdgeLengthById()
	{
		return m_nodeTotalEdgeLengthById;
	}
	/**
	 * Returns distribution of Total edge cost with the key being the node ID.
	 * @return
	 */
	public TreeMap<Integer, Double> getTotalEdgeCostById() {
		return m_nodeTotalEdgeCostById;
	}
	/**
	 * Returns distribution of Max edge cost with the key being the node ID.
	 * @return
	 */
	public TreeMap<Integer, Double> getMaximumEdgeCostById() 
	{
		return m_nodeMaximumEdgeCostById;
	}
	/**
	 * Returns distribution of mean edge cost with the key being the node ID.
	 * @return
	 */
	public TreeMap<Integer, Double> getMeanEdgeCostById() 
	{
		TreeMap<Integer,Double> result = new TreeMap<Integer,Double>();
		
		TreeMap<Integer, Integer> nodeDistId = getNodeDistributionId();
		Vector<Integer> vNodeId = new Vector<Integer>(nodeDistId.keySet());
		for(int i =0 ; i < vNodeId.size(); i++)
		{
			double value =  
				m_nodeTotalEdgeCostById.get(vNodeId.get(i))/
				(double)nodeDistId.get(vNodeId.get(i));
			
			//System.out.println("Nod id " +vNodeId.get(i)+ " mean cost "+value);
			result.put(vNodeId.get(i), value);
			
		}
		return result;
	}
	/**
	 * Returns distribution of Max edge length with the key being the node ID.
	 * @return
	 */
	public TreeMap<Integer, Double> getMaximumEdgeLengthById() 
	{
		return m_nodeMaximumEdgeLengthById;
	}
	/**
	 * Returns distribution of MEan edge length with the key being the node ID.
	 * @return
	 */
	public TreeMap<Integer, Double> getMeanEdgeLengthById() 
	{
		TreeMap<Integer,Double> result = new TreeMap<Integer,Double>();
		
		TreeMap<Integer, Integer> nodeDistId = getNodeDistributionId();
		Vector<Integer> vNodeId = new Vector<Integer>(nodeDistId.keySet());
		for(int i =0 ; i < vNodeId.size(); i++)
		{
			double value =  
				m_nodeTotalEdgeLengthById.get(vNodeId.get(i))/
				(double)nodeDistId.get(vNodeId.get(i));
			
			result.put(vNodeId.get(i), value);
			
		}
		
		
		return result;
	}

	/**
	 * Returns a Map with the key being the Node ID, and value being
	 * the number of nodes connected within GraphConstants.MeasurementR
	 * @return
	 */
	public TreeMap<Integer, Integer> getConnectedNodesWithinRadiusRById() 
	{
		return m_connectedNodesWithinRadiusById;
	}

	/**
	 * Returns a Map with the key being the Node ID, and value mean degree of connected
	 * neighbors
	 * @return
	 */
	public TreeMap<Integer, Double> getMeanDegreeOfConnectedNeighborsById() 
	{
		return m_meanDegreeOfConnectedNeighborsById;
	}

	/**
	 * Returns a Map with the key being the Node ID, and betweeness of the node
	 * @return
	 */
	public TreeMap<Integer, Double> getBetweenessById() 
	{
		return m_betweenessById;
	}
	
	/**
	 * Returns a Map with the key being the Node ID, and Closeness of the node
	 * @return
	 */
	public TreeMap<Integer, Double> getClosenessById() 
	{
		return m_closenessById;
	}

	/**
	 * Returns the generated SPL Graphs, with the index of the array matching the Node ID.
	 * 
	 * This really should be HashMap for when the ids are not continous.
	 * @return
	 */
	public Graph[] getSPLGraphs() {
		return SPLGraphs;
	}
	
	/**
	 * Returns a Map with the key being the Node ID, and value being
	 * the number of nodes within GraphConstants.MeasurementR
	 * @return
	 */
	public TreeMap<Integer, Integer> getNodesWithinRadiusById() 
	{
		return m_nodesWithinRadiusById;
	}

	/**
	 * Returns a Map with the key being the Node ID, and value being
	 * the vector with the edge lengths
	 * @return
	 */
	public TreeMap<Integer,Vector< Double>> getEdgeLengthById() throws Exception
	{
		//id  edge length  - to make the getEdgeLengthById work, we need to make an entry
		// for each edge coming out of this node
		TreeMap<Integer,Vector< Double>> edgeLengthId = new TreeMap<Integer,Vector< Double>>();
		
		Vector<Node> vec = new Vector<Node>(this.getHeaderNodesMap().values());
		
		for(int i = 0; i < vec.size(); i++)
		{
			Node node = vec.get(i);
			
			Vector<Edge> nei = new Vector<Edge>(node.getNeighbors().values());
			
			
			Vector<Double> lengths = new Vector<Double>();
			
			for(int j=0; j < nei.size(); j++)
			{
				Edge ci = nei.get(j);
				double val = ci.getWeight();
				double real_value = NonLinearCostFunction.getDistance(val);
				
				lengths.add( real_value);
			}
			edgeLengthId.put(node.getID(), lengths);
		}
		
		return edgeLengthId;
	}

	/**
	 * Returns a Map with the key being the Node ID, and value being
	 * the vector with the edge lengths greater than .17
	 * @return
	 */
	public TreeMap<Integer, Vector<Double>> getEdgeLengthByIdG17()throws Exception {
		//id  edge length  - to make the getEdgeLengthById work, we need to make an entry
		// for each edge coming out of this node
		TreeMap<Integer,Vector< Double>> edgeLengthId = new TreeMap<Integer,Vector< Double>>();
		
		Vector<Node> vec = new Vector<Node>(this.getHeaderNodesMap().values());
		
		for(int i = 0; i < vec.size(); i++)
		{
			Node node = vec.get(i);
			
			Vector<Edge> nei = new Vector<Edge>(node.getNeighbors().values());
			
			
			Vector<Double> lengths = new Vector<Double>();
			
			for(int j=0; j < nei.size(); j++)
			{
				Edge ci = nei.get(j);
				double val = ci.getWeight();
				double real_value = NonLinearCostFunction.getDistance(val);
				
				if(Double.compare(0.17, real_value) > 0)
				{
					lengths.add( real_value);
				}
			}
			edgeLengthId.put(node.getID(), lengths);
		}
		
		return edgeLengthId;
	}

	/**
	 * Returns a Map with the key being the Node ID, and value being
	 * the vector with the edge lengths less than or equal to .17
	 * @return
	 */
	public TreeMap<Integer, Vector<Double>> getEdgeLengthByIdLT17() throws Exception{
		//id  edge length  - to make the getEdgeLengthById work, we need to make an entry
		// for each edge coming out of this node
		TreeMap<Integer,Vector< Double>> edgeLengthId = new TreeMap<Integer,Vector< Double>>();
		
		Vector<Node> vec = new Vector<Node>(this.getHeaderNodesMap().values());
		
		for(int i = 0; i < vec.size(); i++)
		{
			Node node = vec.get(i);
			
			Vector<Edge> nei = new Vector<Edge>(node.getNeighbors().values());
			
			
			Vector<Double> lengths = new Vector<Double>();
			
			for(int j=0; j < nei.size(); j++)
			{
				Edge ci = nei.get(j);
				double val = ci.getWeight();
				double real_value = NonLinearCostFunction.getDistance(val);
				
				if(Double.compare(0.17, real_value) <= 0)
				{
					lengths.add( real_value);
				}
			}
			edgeLengthId.put(node.getID(), lengths);
		}
		
		return edgeLengthId;
	}


	/**
	 * Updates WithinR calculations.  This was done to specifically calculate 
	 * certain measurements.
	 */	
	public void UpdateWithinRInfo () throws Exception
	{
	
		ConnectedNodesWithinRadius connectedNodesWithinRadius = new ConnectedNodesWithinRadius(this);
		m_connectedNodesWithinRadius = connectedNodesWithinRadius.Measure();
		m_connectedNodesWithinRadiusById = connectedNodesWithinRadius.ById();
		
	}
	
	/**
	 * Updates WithinR calculations.  This was done to specifically calculate 
	 * certain measurements using a different value than GraphConstants.MeasurementR
	 * @return
	 */	
	public void UpdateWithinRInfo(double val)  throws Exception
	{
		ConnectedNodesWithinRadius connectedNodesWithinRadius = new ConnectedNodesWithinRadius(this);
		connectedNodesWithinRadius.setRadius(val);
		m_connectedNodesWithinRadius = connectedNodesWithinRadius.Measure();
		m_connectedNodesWithinRadiusById = connectedNodesWithinRadius.ById();
		
	}

	/**
	 * Returns a Map with the key being the Node ID, and value being the 
	 * the number of connected nodes outside of radius GraphConstants.MeasurementR.
	 * 
	 * @return
	 */	
	public TreeMap<Integer, Integer> getConnectedNodesOutsideRadiusRById(double val) throws Exception
	{
		ConnectedNodesOutsideRadius connectedNodesOutsideRadius = new ConnectedNodesOutsideRadius(this);
		connectedNodesOutsideRadius.setRadius(val);
		
		TreeMap<Integer, Integer> connectedNodesOutsideRadiusById;
		
		connectedNodesOutsideRadius.Measure();
		connectedNodesOutsideRadiusById = connectedNodesOutsideRadius.ById();
		
		return connectedNodesOutsideRadiusById;
	}
	
	/**REturns the number of edges that are shorter than specified value.
	 * 	
	 * @param d
	 * @return
	 * @throws Exception
	 */
	public double getNumberShortEdges(double d) throws Exception 
	{

		ConnectedNodesWithinRadius connectedNodesWithinRadius = new ConnectedNodesWithinRadius(this);
		connectedNodesWithinRadius.setRadius(d);
		m_connectedNodesWithinRadius = connectedNodesWithinRadius.Measure();
		m_connectedNodesWithinRadiusById = connectedNodesWithinRadius.ById();
		
		return connectedNodesWithinRadius.getNumEdges();
		
	}
	/**REturns the number of edges that are longer than specified value.
	 * 	
	 * @param d
	 * @return
	 * @throws Exception
	 */
	public double getNumberLongEdges(double d) throws Exception 
	{

		ConnectedNodesOutsideRadius connectedNodesOutsideRadius = new ConnectedNodesOutsideRadius(this);
		connectedNodesOutsideRadius.setRadius(d);

		connectedNodesOutsideRadius.Measure();
		connectedNodesOutsideRadius.ById();

		return connectedNodesOutsideRadius.getNumEdges();
	}

	/**
	 * Returns a Map with the key being the Node ID, and value being the 
	 * the number of connected nodes within radius of .5 to .16.
	 * 
	 * @return
	 */	
	public TreeMap<Integer, Integer> getConnectedNodesWithinRadiusR12ById() throws Exception
	{
		ConnectedNodesWithinRadius1Radius2 connectedNodesWithinRadius1Radius2 = new ConnectedNodesWithinRadius1Radius2(this);
		connectedNodesWithinRadius1Radius2.setRadius1(.16);
		connectedNodesWithinRadius1Radius2.setRadius2(.5);
		
		TreeMap<Integer, Integer> m_connectedwithRadius12ById;
		
		connectedNodesWithinRadius1Radius2.Measure();
		m_connectedwithRadius12ById = connectedNodesWithinRadius1Radius2.ById();
		
		return m_connectedwithRadius12ById;
	}

	/**
	 * Returns a Map with the key being the Node ID, and value being the 
	 * the number of nodes within radius of .5 to .16.
	 * 
	 * @return
	 */	
	public TreeMap<Integer, Integer> getNodesWithinRadius12ById() throws Exception
	{

		NodesWithinRadius1Radius2 nodesWithinRadius1Radius2 = new NodesWithinRadius1Radius2(this);
		nodesWithinRadius1Radius2.setRadius1(.16);
		nodesWithinRadius1Radius2.setRadius2(.5);
		
		TreeMap<Integer, Integer> m_nodesWithRadius12ById;
		
		nodesWithinRadius1Radius2.Measure();
		m_nodesWithRadius12ById = nodesWithinRadius1Radius2.ById();
		
		return m_nodesWithRadius12ById;
	}

	/**
	 * Returns a Map with the key being the Node ID, and degree of the node
	 * minus the number of connected nodes outside of .16.
     *
	 * @return
	 */	
	public TreeMap<Integer, Integer> getSpecialDegreeById() throws Exception
	{

		TreeMap<Integer,Integer> specialDegreebyId = new TreeMap<Integer,Integer>();
		
		ConnectedNodesOutsideRadius connectedNodesOutsideRadius = new ConnectedNodesOutsideRadius(this);
		connectedNodesOutsideRadius.setRadius(.16);
		
		TreeMap<Integer, Integer> connectedNodesOutsideRadiusById;
		
		connectedNodesOutsideRadius.Measure();
		connectedNodesOutsideRadiusById = connectedNodesOutsideRadius.ById();		
		
		//degree  # of degree
		TreeMap<Integer,Integer> nodeDistributionById = getNodeDistributionId();
		
		Vector<Integer> vkeys = new Vector<Integer>(nodeDistributionById.keySet());
		
		
		for(int i = 0; i < vkeys.size(); i++)
		{
			int degree = nodeDistributionById.get(vkeys.get(i)) 
					- connectedNodesOutsideRadiusById.get(vkeys.get(i));
						
			specialDegreebyId.put(vkeys.get(i), degree);
		}
		
		return nodeDistributionById;

	}

	/**
	 * Returns a Map with the key being the Node ID, and value being the 
	 * the number of nodes outside of specified radius.
	 * 
	 * @return
	 */	
	public TreeMap<Integer, Integer> getNodesOutsideRadius(double d) throws Exception
	{
		NodesOutsideRadius nodesOutsideRadius = new NodesOutsideRadius(this);
		nodesOutsideRadius.setRadius(d);
		
		TreeMap<Integer, Integer> m_connectedNodesOutsideRadiusById;
		
		nodesOutsideRadius.Measure();
		m_connectedNodesOutsideRadiusById = nodesOutsideRadius.ById();
		
		return m_connectedNodesOutsideRadiusById;
	}

	/**
	 * Returns a Map with the key being the Node ID, and degree of the node
	 * minus the number of connected nodes outside of specified radius
     *
	 * @return
	 */	
	public TreeMap<Integer, Integer> getSpecialDegreeById(double d) throws Exception
	{
		TreeMap<Integer,Integer> specialDegreebyId = new TreeMap<Integer,Integer>();
//		
//		ConnectedNodesWithinRadius connectedNodesWithinRadius = new ConnectedNodesWithinRadius(this);
//		connectedNodesWithinRadius.setRadius(.16);
//		TreeMap<Integer,Integer> connectedNodesWithinRadiusData = connectedNodesWithinRadius.Measure();
//		TreeMap<Integer,Integer> connectedNodesWithinRadiusById = connectedNodesWithinRadius.ById();
//		
		
		ConnectedNodesOutsideRadius connectedNodesOutsideRadius = new ConnectedNodesOutsideRadius(this);
		connectedNodesOutsideRadius.setRadius(d);
		
		TreeMap<Integer, Integer> connectedNodesOutsideRadiusById;
		
		connectedNodesOutsideRadius.Measure();
		connectedNodesOutsideRadiusById = connectedNodesOutsideRadius.ById();		
		
		//degree  # of degree
		TreeMap<Integer,Integer> nodeDistributionById = getNodeDistributionId();
		
		Vector<Integer> vkeys = new Vector<Integer>(nodeDistributionById.keySet());
		
		
		for(int i = 0; i < vkeys.size(); i++)
		{

			
			int degree = nodeDistributionById.get(vkeys.get(i)) 
//					- connectedNodesWithinRadiusById.get(vkeys.get(i))
					- connectedNodesOutsideRadiusById.get(vkeys.get(i));
						
			specialDegreebyId.put(vkeys.get(i), degree);
		}
		
		return nodeDistributionById;
		
		
	}
	
	/**
	 * This method was implemented at the last minute.
	 * 
	 * 
	 * plotConnectedNodesR1R2VsNodesWithinR1R2LargeConnected.
     *
	 * @return
	 */	
	public TreeMap<Integer, Vector<Double>> pat2() throws Exception
	{
		TreeMap<Integer, Vector<Double>> result =null;
		 
		EdgeLengthDistribution edgeLengthDistribution = 
			new EdgeLengthDistribution(this);
		edgeLengthDistribution.setR(.16);
		edgeLengthDistribution.Measure();
		result = edgeLengthDistribution.ByIdData();
		
		return result;
	}
	
	/**
	 * Returns a Map with the key being the Node ID, and value being the 
	 * the number of connected nodes within of specified radius .
	 * 
	 * @return
	 */	
	public TreeMap<Integer, Integer> getConnectedNodesWithinRadius(double d) throws Exception
	{
		ConnectedNodesWithinRadius connectedNodesWithinRadius = new ConnectedNodesWithinRadius(this);
		connectedNodesWithinRadius.setRadius(d);
		TreeMap<Integer, Integer> connectedNodesWithinRadiusData = connectedNodesWithinRadius.Measure();
		//m_connectedNodesWithinRadiusById = connectedNodesWithinRadius.ById();
		
		return connectedNodesWithinRadiusData;
	}

	
}
