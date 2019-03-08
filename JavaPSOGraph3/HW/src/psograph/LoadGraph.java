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
import java.io.FilenameFilter;
import java.util.Vector;

import psograph.graph.CalculatedGraph;
import psograph.graph.Graph;
import psograph.graph.NodeLocationCalculator;
import psograph.util.Util; 


public class LoadGraph {

	

	
	/**
	 * @param args
	 */
	public static void main(String[] args) 
	{

		try
		{
			System.out.println("getting CG1");
			CalculatedGraph cg1 = Util.streaminCalculatedGraph(new File("C:\\TestHeuristic2\\Seed0\\optimizedheuristic1\\heuristicGraph0.CalculatedGraph"));
			
			CalculatedGraph cg2 = Util.streaminCalculatedGraph(new File("C:\\TestHeuristic2\\Seed0\\optimizedheuristic2\\heuristicGraph0.CalculatedGraph"));
			
			CalculatedGraph cg3 = Util.streaminCalculatedGraph(new File("C:\\TestHeuristic2\\Seed0\\optimizedheuristic3\\heuristicGraph0.CalculatedGraph"));
			
			Graph costBasisGraph = Util.streaminGraph(new File("C:\\TestHeuristic2\\Seed0\\seed\\ExponentialCostBasis1.Graph"));
			
			NodeLocationCalculator calc = new NodeLocationCalculator(costBasisGraph);
			calc.calculateResults();
			
			double costBasisCost = calc.getCost();  
			
			//327252.8216994972
			cg1.setCostBasis(costBasisCost);
			
			cg1.UpdatePSOCalculations();
			
			
			cg2.setCostBasis(costBasisCost);
			
			cg2.UpdatePSOCalculations();
			
			cg3.setCostBasis(costBasisCost);
			
			cg3.UpdatePSOCalculations();
			
			System.out.println("got CG1 and CG2");
			System.out.println("getting expo");
			//CalculatedGraph expoGraph = Util.streaminCalculatedGraph(new File("C:\\TestHeuristic2\\Seed0\\exponential\\exponentialGraph0.CalculatedGraph"));
			
		//	expoGraph.setCostBasis(costBasisCost);
			
		//	expoGraph.UpdatePSOCalculations();
			System.out.println("got expo");
			
			
			


			int cg1_TotalEdges = cg1.getTotalEdges();
			int cg2_TotalEdges = cg2.getTotalEdges();
		//	int expo_TotalEdges = expoGraph.getTotalEdges();
			//
			System.out.println("----------------------------------------");
			System.out.println("CG1 : numNodes: "+cg1.getNumberOfNodes()+" numEdges: "+cg1.getTotalEdges());
			System.out.println("CG1 :Fitness: "+ cg1.getFitnessValue()+
					": Total Number of edges : "+cg1.getTotalEdges());
			System.out.println("----------------------------------------");
			System.out.println("CG2 : numNodes: "+cg2.getNumberOfNodes()+" numEdges: "+cg2.getTotalEdges());
			System.out.println("CG2 :Fitness: "+ cg2.getFitnessValue()+
					": Total Number of edges : "+cg2.getTotalEdges());
			System.out.println("----------------------------------------");
			System.out.println("CG3 : numNodes: "+cg3.getNumberOfNodes()+" numEdges: "+cg3.getTotalEdges());
			System.out.println("CG3 :Fitness: "+ cg3.getFitnessValue()+
					": Total Number of edges : "+cg3.getTotalEdges());
			System.out.println("----------------------------------------");
//			System.out.println("Expo : Cost: "+ expoGraph.getCost()+
//					": Total Number of edges : "+expoGraph.getTotalEdges());
//			System.out.println("Expo :Fitness: "+ expoGraph.getFitnessValue()+
//					": Total Number of edges : "+expoGraph.getTotalEdges());


			System.out.println("----------------------------------------");
			int pat =0;
			pat = 2;
			
			
//		    // It is also possible to filter the list of returned files.
//		    // This example does not return any files that start with `.'.
//		    FilenameFilter filter = new FilenameFilter() {
//		        public boolean accept(File dir, String name) {
//		            return name.endsWith(".CalculatedGraph");
//		        }
//		    };
			
//			File dir = new File ("C:\\finalOptimizationRuns\\seed0_3\\PSO_40");
//	    	File calculatedGraphs[] = dir.listFiles(filter);
//			if (calculatedGraphs == null)
//				throw new Exception ("SeedDirectory is empty");
//			
//
//			int maxCalculatedGraphs = calculatedGraphs.length;
//
//			System.out.println("------C:\\finalOptimizationRuns\\seed0_3\\PSO_40-------");
//			for(int j=0; j < maxCalculatedGraphs; j++)
//			{
//				CalculatedGraph graph =  null;
//
//				//System.out.println("CalculatedGraph was retrieved Graph " + calculatedGraphs[j].toString());
//				graph = thesis.util.Util.streaminCalculatedGraph(calculatedGraphs[j]);
//				
//			//	System.out.println("Fitness: "+ graph.getFitnessValue()+
//			//			": Total Number of edges : "+graph.getTotalEdges());
//				
//
//			}
//			System.out.println("----------------------------------------");
//			
//			dir = new File ("C:\\finalOptimizationRuns\\seed0_3\\graphs");
//	    	calculatedGraphs = dir.listFiles(filter);
//			if (calculatedGraphs == null)
//				throw new Exception ("SeedDirectory is empty");
//			
//
//			maxCalculatedGraphs = calculatedGraphs.length;
//
//			System.out.println("---------C:\\finalOptimizationRuns\\seed0_3\\graphs--------");
//			for(int j=0; j < maxCalculatedGraphs; j++)
//			{
//				CalculatedGraph graph =  null;
//
//				//System.out.println("CalculatedGraph was retrieved Graph " + calculatedGraphs[j].toString());
//				graph = thesis.util.Util.streaminCalculatedGraph(calculatedGraphs[j]);
//				
//			//	System.out.println("Fitness: "+ graph.getFitnessValue()+
//			//			": Total Number of edges : "+graph.getTotalEdges());
//				
//
//			}
//			System.out.println("----------------------------------------");
		}
		catch (Exception e)
		{
			e.printStackTrace();
			
		}
		
		
	}

	

	


    	



    	

		
 
	
}
