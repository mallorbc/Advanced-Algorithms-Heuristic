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
import java.io.FilenameFilter;
import java.util.Vector;

import psograph.analysis.*;
import psograph.graph.CalculatedGraph;
import psograph.graph.Graph;
import psograph.graph.NodeLocationCalculator;
import psograph.util.Util;


public class AnalyzeGraphs extends LoadGraph {

	
	//settings
	static boolean DO_LAYOUT = false;
	static boolean DO_VALUES = false;
	static boolean DO_DIST = false;
	static boolean DO_NODE_VALUES = false;
	static boolean DO_OPTIMIZED_COMPARE = false;
	static boolean DO_UPDATE = false;
	static boolean DO_WRITE_OUT_UPDATED_GRAPHS = false;
	static boolean DO_COMBINED_PLOTS = false;
	static boolean DO_COMBINED_PLOTS2 = true;


	public void OptimizedCompare()
	{
		
		
	}
	
	public static void main(String[] args) 
	{
	    // It is also possible to filter the list of returned files.
	    // This example does not return any files that start with `.'.
	    FilenameFilter filter = new FilenameFilter() {
	        public boolean accept(File dir, String name) {
	            return name.endsWith(".CalculatedGraph");
	        }
	    };
		
		
		
		try
		{
			Vector<File> baseDir = new Vector<File>();
			
			baseDir.add(new File("C:\\Comparisions\\Seed0"));
			baseDir.add(new File("C:\\Comparisions\\Seed1"));
			baseDir.add(new File("C:\\Comparisions\\Seed2"));
			baseDir.add(new File("C:\\Comparisions\\Seed3"));
			baseDir.add(new File("C:\\Comparisions\\Seed4"));
			baseDir.add(new File("C:\\Comparisions\\Seed5"));
			
			
//			if(DO_OPTIMIZED_COMPARE == true)
//			{
//				
//				
//				Vector<CalculatedGraph> seeds = new Vector<CalculatedGraph>();
//				Vector<Vector<CalculatedGraph>> optimalSolns = new Vector<Vector<CalculatedGraph>>();
//				for(int i = 0; i < baseDir.size(); i++)
//					//1 ; i++)
//					//baseDir.size(); i++)
//				{
//					
//					
//					optimalSolns.add(i,new Vector<CalculatedGraph>());
//					System.out.println("Processing "+baseDir.get(i).getAbsoluteFile());
//
//
//					File optimalDir = new File(baseDir.get(i).getAbsolutePath()+"\\optimized");
//					File optimalFile[] = optimalDir.listFiles();
//
//					System.out.println("Loading Optimized"+i);
//					for(int j =0; j < optimalFile.length; j++)
//					{
//
//						if(optimalFile[j].isDirectory())
//						{
//
//							CalculatedGraph solnToExamine = 
//								Util.streaminCalculatedGraph(new File(optimalFile[j].getAbsolutePath()+"\\globalBest_CalculatedGraph.CalculatedGraph"));
//							//Make sure all calculations are up to date
//							if(DO_UPDATE)
//							{
//								solnToExamine.UpdateCalcuations();
//							}
//							optimalSolns.get(i).add(solnToExamine);
//							if(j==0)
//								seeds.add(solnToExamine);
//
//							System.out.println(optimalFile[j].getName() +" has fitness of "+solnToExamine.getFitnessValue()+" with # edges "+solnToExamine.getTotalEdges());
//							if(DO_WRITE_OUT_UPDATED_GRAPHS)
//							{
//								Util.streamOutCalculatedGraphWithExtension(optimalFile[j],
//										"globalBest_CalculatedGraph.CalculatedGraph",
//										solnToExamine);
//
//							}
//
//						}
//					}
//					System.out.println("Done Loading Optimized");
//				}


//				System.out.println("Performing optimized value Plots");
//				OptimizedValuePlotWriter optimizedValuePlotWriter = 
//					new OptimizedValuePlotWriter(optimalSolns, new File(baseDir.get(0).getParent()));
//				optimizedValuePlotWriter.generatePlots();
//				System.out.println("Done Performing optimized value Plots");
				

				
				
//				System.out.println("Performing Layout Plots");
//				LayoutPropertiesWriter layoutPropertiesWriter = 
//					new LayoutPropertiesWriter(seeds, new File(baseDir.get(0).getParent()));
//				layoutPropertiesWriter.generateCombinedPlots();
//				System.out.println("Done with Layout Plots");


//			}	
			
			Vector<CalculatedGraph> optimalSolnsAll = new Vector<CalculatedGraph>();
			Vector<CalculatedGraph> exponentialSolnsAll = new Vector<CalculatedGraph>();
			Vector<CalculatedGraph> heuristicSolnsAll = new Vector<CalculatedGraph>();
			
			for(int i = 0; i < baseDir.size(); i++)
			//1 ; i++)
				//baseDir.size(); i++)
			{
				System.out.println("Processing "+baseDir.get(i).getAbsoluteFile());
				
				Graph costBasisGraph = Util.streaminGraph(new File(baseDir.get(i)+"\\seed\\ExponentialCostBasis1.Graph"));
				double costBasis = costBasisGraph.getCost();
				
				System.out.println("CostBasis - "+costBasis);
				
				Vector<CalculatedGraph> optimalSolns = new Vector<CalculatedGraph>();
				File optimalDir = new File(baseDir.get(i).getAbsolutePath()+"\\optimized");
				File optimalFile[] = optimalDir.listFiles();
				
				System.out.println("Loading Optimized");
				for(int j =0; j < optimalFile.length; j++)
				{
					if(optimalFile[j].isDirectory())
					{
						
						CalculatedGraph solnToExamine = 
							Util.streaminCalculatedGraph(new File(optimalFile[j].getAbsolutePath()+"\\globalBest_CalculatedGraph.CalculatedGraph"));
						//Make sure all calculations are up to date
						if(DO_UPDATE)
						{
							
							solnToExamine.setCostBasis(costBasis);
							
						    //solnToExamine.UpdateCalcuations();
							solnToExamine.UpdatePSOCalculations();
						}
						optimalSolns.add(solnToExamine);
						optimalSolnsAll.add(solnToExamine);

						System.out.println(optimalFile[j].getName() +" has fitness of "+solnToExamine.getFitnessValue()+" with # edges "+solnToExamine.getTotalEdges());
					//	System.out.println("    cost:"+solnToExamine.getCost()+" AISPL:"+solnToExamine.getAISPL()+" pLCC:"+solnToExamine.getLCC());
						if(DO_WRITE_OUT_UPDATED_GRAPHS)
						{
							Util.streamOutCalculatedGraphWithExtension(optimalFile[j],
									"globalBest_CalculatedGraph.CalculatedGraph",
									solnToExamine);

						}
						
					}
				}
				System.out.println("Done Loading Optimized");
				
				if(optimalSolns.size() == 0 )
					continue;
					
				Vector<CalculatedGraph> exponentialSolns = new Vector<CalculatedGraph>();
				File exponentialDir = new File(baseDir.get(i).getAbsolutePath()+"\\exponential");
				
				//File exponentialDir = new File(baseDir.get(i).getAbsolutePath()+"\\exponential");
				File exponentialFile[] = exponentialDir.listFiles(filter);
				
				double cost = -1;
				String costExp =null;
				
				System.out.println("Loading Exponential");
				for(int j =0;  j < exponentialFile.length; j++)
				//j < 5; j++)
				//j < exponentialFile.length; j++)
				{

					CalculatedGraph solnToExamine = 
						Util.streaminCalculatedGraph(exponentialFile[j]);
					//Make sure all calculations are up to date
					if(DO_UPDATE)
					{	
						solnToExamine.setCostBasis(costBasis);
						
					    //solnToExamine.UpdateCalcuations();
						solnToExamine.UpdatePSOCalculations();
					}
					exponentialSolns.add(solnToExamine);
					exponentialSolnsAll.add(solnToExamine);

					System.out.println(exponentialFile[j].getName() +" has fitness of "+solnToExamine.getFitnessValue()+" with # edges "+solnToExamine.getTotalEdges());
					System.out.println("    cost:"+solnToExamine.getCost()+" AISPL:"+solnToExamine.getAISPL()+" pLCC:"+solnToExamine.getLCC());
				//	System.out.println(exponentialFile[j].getName() +" has cost of "+solnToExamine.getCost());
					if(solnToExamine.getCost() > cost)
					{
						cost = solnToExamine.getCost();
						costExp = exponentialFile[j].getName();
					}
					
					
					if(DO_WRITE_OUT_UPDATED_GRAPHS)
					{
						Util.streamOutCalculatedGraphWithExtension(exponentialFile[j].getParentFile(),
								exponentialFile[j].getName(),
								solnToExamine);
					}
					
				}
				System.out.println("Done Loading Exponential");
				System.out.println("Most expensive network - "+costExp+" with cost of "+cost);
				
			
				Vector<CalculatedGraph> heuristicSolns = new Vector<CalculatedGraph>();
				File heuristicDir = new File(baseDir.get(i).getAbsolutePath()+"\\optimizedheuristic1a");
				File heuristicFile[] = heuristicDir.listFiles(filter);
				
				System.out.println("Loading Hueristic");
				for(int j =0;   j < heuristicFile.length; j++)
				{

					CalculatedGraph solnToExamine = 
						Util.streaminCalculatedGraph(heuristicFile[j]);
					//Make sure all calculations are up to date
					if(DO_UPDATE)
					{
						solnToExamine.setCostBasis(costBasis);
						
					    //solnToExamine.UpdateCalcuations();
						solnToExamine.UpdatePSOCalculations();
					}
				
					
					heuristicSolns.add(solnToExamine);
					heuristicSolnsAll.add(solnToExamine);
					
					System.out.println(heuristicFile[j].getName() +" has fitness of "+solnToExamine.getFitnessValue()+" with # edges "+solnToExamine.getTotalEdges());
					System.out.println("    cost:"+solnToExamine.getCost()+" AISPL:"+solnToExamine.getAISPL()+" pLCC:"+solnToExamine.getLCC());


					if(DO_WRITE_OUT_UPDATED_GRAPHS)
					{
						Util.streamOutCalculatedGraphWithExtension(heuristicFile[j].getParentFile(),
								heuristicFile[j].getName(),
								solnToExamine);
					}
					
				}
				System.out.println("Done Loading Hueristic");
			
				if(DO_LAYOUT == true)
				{
					System.out.println("Performing Layout Plots");
					LayoutPropertiesWriter layoutPropertiesWriter = 
						new LayoutPropertiesWriter(optimalSolns.get(0), baseDir.get(i), "Layout Seed "+i);
					layoutPropertiesWriter.generatePlots();
					System.out.println("Done with Layout Plots");
				}
				
				if(DO_VALUES == true)
				{
					System.out.println("Performing value Plots");
					ValuePlotWriter valuePlotWriter = 
						new ValuePlotWriter(optimalSolns, exponentialSolns, heuristicSolns, baseDir.get(i));
					valuePlotWriter.generatePlots();
					System.out.println("Done Performing value Plots");
				}
				
				if(DO_DIST == true)
				{
					System.out.println("Performing dist Plots");
					DistPlotWriter distPlotWriter = 
						new DistPlotWriter(optimalSolns, exponentialSolns, heuristicSolns, baseDir.get(i));
					distPlotWriter.generatePlots();
					System.out.println("Done Performing dist Plots");
				}
			
				if(DO_NODE_VALUES == true)
				{
					System.out.println("Performing node value Plots");
					NodeValuePlotWriter nodeValuePlotWriter = 
						new NodeValuePlotWriter(optimalSolns, exponentialSolns, heuristicSolns, baseDir.get(i));
					nodeValuePlotWriter.generatePlots();
					System.out.println("Done Performing node value Plots");
				}
				

				

			
			}
		
			
		if(DO_COMBINED_PLOTS2 == true)
		{
//			System.out.println("Performing node value Plots");
//			NodeValuePlotWriter nodeValuePlotWriter = 
//				new NodeValuePlotWriter(optimalSolnsAll, exponentialSolnsAll, heuristicSolnsAll, baseDir.get(0).getParentFile());
//			//nodeValuePlotWriter.generateConnectWithinRPlots();
//			nodeValuePlotWriter.generatePlots2();
//			System.out.println("Done Performing node value Plots");
			
			
//			System.out.println("Performing node value Plots");
//			NodeValuePlotWriter nodeValuePlotWriter = 
//				new NodeValuePlotWriter(optimalSolnsAll, exponentialSolnsAll, heuristicSolnsAll, baseDir.get(0).getParentFile());
//			nodeValuePlotWriter.plotNodeDegreeVsConnectedNodesOutsideRadiusR();
//			System.out.println("Done Performing node value Plots");
			
			
			System.out.println("Performing value Plots");
			ValuePlotWriter valuePlotWriter = 
				new ValuePlotWriter(optimalSolnsAll, exponentialSolnsAll, heuristicSolnsAll,  baseDir.get(0).getParentFile());
			valuePlotWriter.generateHubbinessPlot();
			System.out.println("Done Performing value Plots");
			
			
			
		}
			
			
		if(DO_COMBINED_PLOTS==true)
		{

//				System.out.println("Performing value Plots");
//				ValuePlotWriter valuePlotWriter = 
//					new ValuePlotWriter(optimalSolnsAll, exponentialSolnsAll, heuristicSolnsAll, baseDir.get(0).getParentFile());
//				valuePlotWriter.generatePlots();
//				System.out.println("Done Performing value Plots");


//				System.out.println("Performing dist Plots");
//				DistPlotWriter distPlotWriter = 
//						new DistPlotWriter(optimalSolnsAll, exponentialSolnsAll, heuristicSolnsAll, baseDir.get(0).getParentFile());
//				distPlotWriter.generatePlots();
//				System.out.println("Done Performing dist Plots");
//				
//				System.out.println("Performing optimized value Plots");
//				PLDistPlotWriter optimizedValuePlotWriter =  
//					new PLDistPlotWriter(optimalSolnsAll, exponentialSolnsAll, heuristicSolnsAll, baseDir.get(0).getParentFile());
//				optimizedValuePlotWriter.generatePlots();
//				System.out.println("Done Performing optimized value Plots");

//
//				System.out.println("Performing node value Plots");
//				NodeValuePlotWriter nodeValuePlotWriter = 
//					new NodeValuePlotWriter(optimalSolnsAll, exponentialSolnsAll, heuristicSolnsAll, baseDir.get(0).getParentFile());
//				nodeValuePlotWriter.generatePlots();
//				System.out.println("Done Performing node value Plots");


		}
			
			

		}
		catch(Exception e)
		{
			e.printStackTrace();
		}

	}

}
