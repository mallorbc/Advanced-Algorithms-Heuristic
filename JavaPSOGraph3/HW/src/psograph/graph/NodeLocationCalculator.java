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

import java.io.Serializable;
import java.util.Random;
import java.util.TreeMap;
import java.util.Vector;






/**
 * A node location calculator takes a Graph in, and then connects it to every
 * other node.  This is useful for determining what node to pick as
 * you can make the decision based up cost/length.
 * @author Patrick
 *
 */
public class NodeLocationCalculator  extends Graph implements Serializable
{
	static final long serialVersionUID = 85L;
	
	/**
	 * Copy constructor
	 * @param g
	 * @param isLinearCostFunction determines linear or non-linear cost function
	 */
	public NodeLocationCalculator (Graph g, boolean isLinearCostFunction)
	{	
		super(g);
		
		if(isLinearCostFunction == true)
		{
			m_CostFunction = new LinearCostFunction();
		}
		else
		{
			m_CostFunction = new NonLinearCostFunction();
		}
		
		m_choseOrNot = new Random();
	}
	
	/**
	 * Copy constructor
	 * @param g
	 */
	public NodeLocationCalculator (Graph g)
	{	
		super(g);
		

		m_CostFunction = new NonLinearCostFunction();

		
		m_choseOrNot = new Random();
	}
	/**
	 * Returns the node which is the furthest from specified node.
	 * If there are mutliple Nodes that tied forfurthest neighbor this picks a random
	 * Node.
	 * @param from
	 * @return
	 * @throws Exception
	 */
	public Integer pickFurthestNeighbor(Node from) throws Exception
	{
		Integer target =-1;
		boolean found = false;
		
		TreeMap<Double, Vector<Integer>> dist = from.getNeighborDistribution();
		if(dist == null)
		{
			return -1;
		}
		
		Vector<Double> keys = new Vector<Double>(dist.keySet());
		int keysize = keys.size();
		
		for(int i=keysize; (i > 0) && (found != true); i=i-1)
		{
		//	double costOfInterest = keys.get(i-1);
			Vector<Integer> vNodes2 = new Vector<Integer>(dist.get(keys.get(i-1)));
		//	int vNodesSize = vNodes2.size();
						
			int longCandidates = vNodes2.size();
			
			if(longCandidates == 1)
			{
				target =  vNodes2.get(0);
				found = true;
			}
			else
			{
				//System.out.println("pickFurthestNeighbor - Choices :"+longCandidates);
				int nodePicked = m_choseOrNot.nextInt(longCandidates);
				target = vNodes2.get(nodePicked);
				found = true;
			}	
		}
		if(found == false)
			throw new Exception ("No node to connect to. NodeLocationCalculator::pickFurtherstNeighbor");
		return target;
	}
	
	public Node chooseNode(Node n)throws Exception
	{
		return chooseNode(n.getID());
	}
	
	
	/**
	 * Picks a node that calculator knows we are not connected to.
	 * On a percentage it will be a close node or random node.
	 */
	public Node chooseNode(int id)throws Exception
	{
		Node n = null;
		
		double per = m_choseOrNot.nextDouble();
		
		//System.out.print("m_choseOrNot - choice "+per);
		
		if(per > m_chooseClosetPercentage)
		{
			//System.out.println("random");
			n = chooseRandomNode(id);
			m_random++;
		}
		else
		{
			//System.out.println("closest");
			n = chooseCloseNode(id);
			m_closest++;
		}
		return n;
	}
	/**
	 * Choose random node to connect the inputted Node to.  This will specifically 
	 * exclude the 5 closest nodes.
	 * @param node
	 * @return
	 * @throws Exception
	 */
	public Node chooseRandomNode(Node node)throws Exception
	{
		return chooseRandomNode(node.getID());
	}
	
	/**
	 * Choose random node to connect the inputted Node to.  This will specifically 
	 * exclude the 5 closest nodes.
	 * @param id
	 * @return
	 * @throws Exception
	 */
	public Node chooseRandomNode(int id) throws Exception
	{
		Node n = null;
		
		Node from = getHeaderNodesMap().get(id);
		
		Vector<Integer> neighbors = new Vector<Integer>( from.getNeighbors().keySet() );
		
		if(neighbors.size()==0)
			return null;
		
		//Need to exclude the 5 or so close nodes
		int N = 5;
		
		//System.out.println("Initial pool size "+neighbors.size());
		
		Vector<Integer> candNodes = getCloseNNodes(id, N);	
		
		if(candNodes.size() != neighbors.size())
		{
			for(int i = 0; i < candNodes.size(); i++)
			{
				if(neighbors.contains(candNodes.get(i)))
					neighbors.remove(candNodes.get(i));
				else
				{
					throw new Exception("NodeLocationCalc::chooseRandomNode  What does this mean???");
				}
			}
		}
		//System.out.println("adjusted pool size "+neighbors.size());
		
		if(neighbors.size() == 0)
		{
			System.out.println("Only close nodes to choose from");
			return null;
		}
		
		int n_id = m_choseOrNot.nextInt(neighbors.size());
		
		n = getHeaderNodesMap().get(neighbors.get(n_id));
		
	//	System.out.println("picked randome " + n.m_id +" for node "+id);
			
		return n;		
	}
	
	/** Choose close node.
	 * 
	 * @param node
	 * @return
	 */
	public Node chooseCloseNode(Node node)
	{
		return chooseCloseNode(node.getID());
	}
	
	
	/**
	 * Get N close nodes to specified ID node
	 * @param id
	 * @param N
	 * @return
	 */
	public Vector<Integer> getCloseNNodes(int id, int N)
	{
		Vector<Integer> closeNNodes = new Vector<Integer>();
		Node from = getHeaderNodesMap().get(id);		
		TreeMap<Double, Vector<Integer>> pat = from.getNeighborDistribution();		
		Vector<Double> keys = new Vector<Double>(pat.keySet());
		int keysize = keys.size();
		
		for(int i=0; (i < keysize) && (closeNNodes.size() != N); i++)
		{
			Vector<Integer> vNodes = new Vector<Integer>(pat.get(keys.get(i)));
			int vNodesSize = vNodes.size();
			
			if(closeNNodes.size() + vNodesSize > N)
			{
				for(; closeNNodes.size() != N; )
				{
					int vSize = vNodes.size();
					
					int nodePicked = m_choseOrNot.nextInt(vSize);

					
					closeNNodes.add(vNodes.get(nodePicked));	
					vNodes.remove(nodePicked);
					//System.out.println("w "+keys.get(i)+" k "+vNodes.get(j));
				}
			}
			else
			{
				for(int j =0; j < vNodesSize; j++)
				{
					closeNNodes.add(vNodes.get(j));		
					//System.out.println("w "+keys.get(i)+" k "+vNodes.get(j));
				}
			}
		}

		return closeNNodes;
	}
	
	
	/** Choose close node to specified ID
	 * 
	 * @param id
	 * @return
	 */
	public Node chooseCloseNode(int id)
	{
		Node n = null;
				
		//this section needs to be re done to use something the buckets scheme as we 
				
		Vector<Integer> candNodes = getCloseNNodes(id, m_chooseClosestNumber);
		
		int size = candNodes.size();
		int t = m_chooseClosestNumber;
		if(size < m_chooseClosestNumber)
			t = size;
		
		if(t==0)
		{
			System.out.println("t is "+t);
			return null;
		}
		
		int n_id = m_choseOrNot.nextInt(t);
		
		int node_id = candNodes.get(n_id);
		
		//System.out.println("picked closest node of " + node_id +" for node "+id);
		
		n = getHeaderNodesMap().get(node_id);
			
		return n;			
	}
	
	/**
	 * Connects every node to every other node for the whole graph.
	 * @throws Exception
	 */
	public void calculateResults() throws Exception
	{
		Vector<Node> m_nodeVector = new Vector<Node>(getHeaderNodesMap().values());
		
		int i,j;
		for(i=0; i < m_nodeVector.size(); i++)
		{
			for(j=i; j < m_nodeVector.size(); j++)
			{
				if(i != j)
				{

					
					Node n1 =m_nodeVector.get(i);
					Node n2 =m_nodeVector.get(j);
					
					if(isNodeAConnectedToNodeB(n1, n2))
					{
						continue;
					}
					
					double weight = calculateWeight(n1.getX(), n1.getY(), n2.getX(), n2.getY());
					addConnection(n1.getID(), n2.getID(), weight);
				}
			}
		}
	}

	/** Percentage of the time to make a close connection
	 * 
	 * @param chooseClosest
	 */
	public void setClosestPercentage(double chooseClosest)
	{
		m_chooseClosetPercentage = chooseClosest;
	}
	/**
	 * Returns the close percentage
	 * @return
	 */
	public double getClosestPercentage()
	{
		return m_chooseClosetPercentage;
	}
	/**
	 * Sets the number of Close connection to exclude 
	 */
	public void setClosestNumber(int chooseClosestNumber)
	{
		m_chooseClosestNumber = chooseClosestNumber;
	}
	/**
	 * Returns the number of Close connection to exclude 
	 */
	public double getClosestNumber()
	{
		return m_chooseClosestNumber;
	}
	
	/**
	 * Calculates teh Node to Node distance based upon coordinates
	 * @param x_1
	 * @param y_1
	 * @param x_2
	 * @param y_2
	 * @return
	 * @throws Exception
	 */
	private double calculateNodeToNodeDistance(double x_1, double y_1, double x_2, double y_2) throws Exception
	{
		double x_comp = 0;
		double y_comp = 0;
		
		//Calculate distance
		//            |
		//     1      |      2
		//            |
		//------------X2------------
		//            |
		//     3      |      4
		//            |
		//
		//--------
		// type 2
		if(x_1 > x_2)
		{
			x_comp = x_1 - x_2;
		}
		else if ( x_1 < x_2)
		{
			x_comp = x_2 - x_1;
		}
		else if(x_1 == x_2)
		{
			x_comp = 0;
		}
		
		
		if(y_1 > y_2)
		{
			y_comp = y_1 - y_2;
		}
		else if ( y_1 < y_2)
		{
			y_comp = y_2 - y_1;
		}
		else if(y_1 == y_2)
		{
			y_comp = 0;
		}
		
		if(Double.compare(x_comp,0.0) == 0 && Double.compare(y_comp,0.0) == 0)
		{
			System.out.println("Node 1 ("+x_1+","+y_1+")");
			System.out.println("Node 2 ("+x_2+","+y_2+")");
			
			throw new Exception("Somehow points are coincident");
		}
		//System.out.println("--------------------------------------------");
		//System.out.println("Node 1 ("+x_1+","+y_1+")");
		//System.out.println("Node 2 ("+x_2+","+y_2+")");
		//System.out.println("X Comp"+x_comp);
		//System.out.println("Y Comp"+y_comp);
		//System.out.println("--------------------------------------------");
		
		
		double t = Math.sqrt( Math.pow(x_comp,2.0) + Math.pow(y_comp,2.0));

		return t;
		
	}
	
	/**
	 * Calculates weight using cost function based upon X,Y coordinates
	 * @param x_1
	 * @param y_1
	 * @param x_2
	 * @param y_2
	 * @return
	 * @throws Exception
	 */
	public double calculateWeight(double x_1, double y_1, double x_2, double y_2) throws Exception
	{
		double t_sig = calculateNodeToNodeDistance(x_1, y_1, x_2, y_2);
		double cost_sig = m_CostFunction.calculate(t_sig);		
		return cost_sig;
	}
	
	double m_chooseClosetPercentage = .8;
	int m_chooseClosestNumber = 5;
	Random m_choseOrNot;
	
	private ICostFunction m_CostFunction;
	
	public int m_closest = 0;
	public int m_random = 0;

}
