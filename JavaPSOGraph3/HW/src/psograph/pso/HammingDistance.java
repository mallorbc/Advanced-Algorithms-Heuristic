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
package psograph.pso;

import java.util.Random;
import java.util.Vector;

import psograph.graph.*;


public class HammingDistance 
{
	public class NodePairs
	{
		NodePairs(int A_id, int B_id)
		{
			m_A_id = A_id;
			m_B_id = B_id;
		}
		
		public boolean equals (Object obj)
		{
			boolean ret = false;
			
			NodePairs n = (NodePairs)obj;
			if(n.m_A_id == m_A_id  && n.m_B_id == m_B_id)
				ret = true;
			else if (n.m_A_id == m_B_id  && n.m_B_id == m_A_id)
				ret = true;
				
			return ret;
		}
		
		int m_A_id;
		int m_B_id;
	}
	
	HammingDistance()
	{
		r = new Random();
		alreadyPickedPairs = new Vector<NodePairs>();
		
	}
	
	public static int CalculateHammingDistance(CalculatedGraph A,CalculatedGraph B)throws Exception
	{
		//System.out.println("--------------Start A -------------");
		//Util.printMFile(A);
		if(A == null)
			System.out.println("hhelo");
		
		int [][] A_adj = A.getAdjecencyGraph2();
		//System.out.println("--------------End   A -------------");
		
		//System.out.println("--------------Start B -------------");
		//Util.printMFile(B);
		if(B == null)
			System.out.println("hhelo");
		
		int [][] B_adj = B.getAdjecencyGraph2();
		//System.out.println("--------------End   B--------------");		
		
		return CalculateHammingDistance(A_adj, B_adj);
		
	}
	
	public static int CalculateHammingDistance(int[][] A,int[][] B)throws Exception
	{
		int result =0;
		
		int size = A[0].length;
		if(size != B[0].length)
			throw new Exception("CalculateHammingDistance - A and B are different size");
		
		int i,j;
		
		for(i=0; i < size; i++)
		{
			for(j=i+1; j < size; j++ )
			{
				int t = (A[i][j] ^ B[i][j]);
				//System.out.println("A = "+ A[i][j] + ", B = " +B[i][j] +" A ^B ="+t );
				
				result = result + t;				
			}		
		}

		return result;		
	}
	
	public Vector<NodePairs> IdentifyCandidateLinks(CalculatedGraph A,CalculatedGraph B)throws Exception
	{
		//System.out.println("--------------Start A -------------");
	///	Util.printMFile(A);
		int [][] A_adj = A.getAdjecencyGraph2();
		//System.out.println("--------------End   A -------------");
		
		//System.out.println("--------------Start B -------------");
		//Util.printMFile(B);
		int [][] B_adj = B.getAdjecencyGraph2();
		//System.out.println("--------------End   B--------------");		
		
		return IdentifyCandidateLinks(A_adj, B_adj);
		
	}
	
	public Vector<NodePairs> IdentifyCandidateLinks(int[][] A,int[][] B)throws Exception
	{
		Vector<NodePairs> result = new Vector<NodePairs>();
		
		int size = A[0].length;
		if(size != B[0].length)
			throw new Exception("CalculateHammingDistance - A and B are different size");
		
		int i,j;
		
		for(i=0; i < size; i++)
		{
			for(j=i+1; j < size; j++ )
			{
				int t = (A[i][j] ^ B[i][j]);
				//System.out.println("A = "+ A[i][j] + ", B = " +B[i][j] +" A ^B ="+t );
				if(t == 1)
					result.add(new NodePairs(i,j));
		
			}		
		}

		return result;
			
	}
	
	
	public CalculatedGraph AcceptEdge(CalculatedGraph A,CalculatedGraph B)throws Exception
	{
		boolean picked = false;
		CalculatedGraph result = null;
		
		Vector<NodePairs> candidates = IdentifyCandidateLinks(A, B);
		int pick;
		NodePairs pickNodePairs = null;
		
//		Vector<Integer> nodes_id = new Vector<Integer>(A.getHeaderNodesMap().keySet());
//		for(int pat = 0; pat < nodes_id.size(); pat++)
//		{
//			System.out.println("Key number "+ nodes_id.get(pat) + " coresponding Node id "+A.getHeaderNodesMap().get(nodes_id.get(pat)).getID());
//			
//		}
		
		while(!picked && candidates.size() != 0)
		{
			pick = r.nextInt(candidates.size());
			
			//Need to verify that we didn't pick a node pair, that we have already accepted
			
			pickNodePairs = candidates.get(pick);
			
			if(alreadyPickedPairs.contains(pickNodePairs))
			{
				//Pick another pair
				candidates.remove(pickNodePairs);
			}
			else
			{
//				Node node1 = B.getNode(pickNodePairs.m_A_id);
//				Node node2 = B.getNode(pickNodePairs.m_B_id);
//				
//			    boolean removeConnection = !node1.isConnectedTo(node2);
//			    int numEdges = A.getTotalEdges();
//				if(removeConnection && numEdges == 20  )
//				{
//					picked = false;
//					candidates.remove(pickNodePairs);
//				}
//				else
					picked = true;
			}
		}
		
		if(!picked)
			return A;
		
		
		result = new CalculatedGraph(A);

		
		
		Node node1 = B.getNode(pickNodePairs.m_A_id);
		Node node2 = B.getNode(pickNodePairs.m_B_id);
		
	    boolean removeConnection = false;
		if(!node1.isConnectedTo(node2))
			removeConnection = true;
		
		if(removeConnection)
		{
			
			result.removeConnection(pickNodePairs.m_A_id, pickNodePairs.m_B_id);
		}
		else
		{
			//NodeLocationCalculator nodeLocationCalculator = new NodeLocationCalculator(result, false);
			//Now that we have a valid node pair, accept into A and return it
			//double weight = nodeLocationCalculator.calculateWeight(node1.getX(), node1.getY(), node2.getX(), node2.getY());
		//	result.addConnection(pickNodePairs.m_A_id, pickNodePairs.m_B_id, weight);
			result.addConnection(pickNodePairs.m_A_id, pickNodePairs.m_B_id);
			alreadyPickedPairs.add(pickNodePairs);
		}
		
		//System.out.println("--------------modified A -------------");
		//Util.printMFile(result);
		
		int hammingDistance = 
			HammingDistance.CalculateHammingDistance(result, A);
		if(hammingDistance !=1)
			throw new Exception ("AcceptEdge modifed graph too much, new Hamming distance is" + hammingDistance);
		//System.out.println("--------------modified   A--------------");	
		
		return result;
		
	}

	Random r;
	Vector<NodePairs> alreadyPickedPairs;
	
}
