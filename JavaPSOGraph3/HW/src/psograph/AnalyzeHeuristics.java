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


public class AnalyzeHeuristics {

	//settings
	static boolean DO_LAYOUT = false;
	static boolean DO_VALUES = false;
	static boolean DO_DIST = false;
	static boolean DO_NODE_VALUES = false;
	static boolean DO_OPTIMIZED_COMPARE = false;
	static boolean DO_UPDATE = false;
	static boolean DO_WRITE_OUT_UPDATED_GRAPHS = false;
	static boolean DO_COMBINED_PLOTS = true;
	

	public void AnalyzeHeuristics()
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
			//TestHuer3
			baseDir.add(new File("C:\\TestHeuristic3\\Seed0"));
//			baseDir.add(new File("C:\\Comparisions\\Seed1"));
//			baseDir.add(new File("C:\\Comparisions\\Seed2"));
//			baseDir.add(new File("C:\\Comparisions\\Seed3"));
//			baseDir.add(new File("C:\\Comparisions\\Seed4"));
//			baseDir.add(new File("C:\\Comparisions\\Seed5"));
//			
			
			Vector<CalculatedGraph> exponentialSolnsAll = new Vector<CalculatedGraph>();
			Vector<CalculatedGraph> optSolnsAll = new Vector<CalculatedGraph>();
			Vector<CalculatedGraph> heur1All = new Vector<CalculatedGraph>();
			Vector<CalculatedGraph> heur2All = new Vector<CalculatedGraph>();
			Vector<CalculatedGraph> heur3All = new Vector<CalculatedGraph>();
			Vector<CalculatedGraph> heur4All = new Vector<CalculatedGraph>();
			
			for(int i = 0; i < baseDir.size(); i++)
			//1 ; i++)
				//baseDir.size(); i++)
			{
				System.out.println("Processing "+baseDir.get(i).getAbsoluteFile());
				
				Graph costBasisGraph = Util.streaminGraph(new File(baseDir.get(i)+"\\seed\\ExponentialCostBasis1.Graph"));
				double costBasis = costBasisGraph.getCost();
				
				System.out.println("CostBasis - "+costBasis);
				
				System.out.println("SHould be 26856.48335871886");
				
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
							
							solnToExamine.UpdateWithinRInfo(.16);
//							solnToExamine.setCostBasis(costBasis);
//							
//						    //solnToExamine.UpdateCalcuations();
//							solnToExamine.UpdatePSOCalculations();
						}
						optimalSolns.add(solnToExamine);
						optSolnsAll.add(solnToExamine);

						System.out.println(optimalFile[j].getName() +" has fitness of "+solnToExamine.getFitnessValue()+" with # edges "+solnToExamine.getTotalEdges());
					//	System.out.println("    cost:"+solnToExamine.getCost()+" AISPL:"+solnToExamine.getAISPL()+" pLCC:"+solnToExamine.getLCC());
						System.out.println("Number of edges greater than  .16 -" + solnToExamine.getNumberLongEdges(.16) );
						System.out.println("Number of edges less or equal to .16 - " + solnToExamine.getNumberShortEdges(.16));
						
						
						if(DO_WRITE_OUT_UPDATED_GRAPHS)
						{
							Util.streamOutCalculatedGraphWithExtension(optimalFile[j],
									"globalBest_CalculatedGraph.CalculatedGraph",
									solnToExamine);

						}
						
					}
				}
				System.out.println("Done Loading Heur1");
				
				if(optimalSolns.size() == 0 )
					continue;
				
				
				Vector<CalculatedGraph> exponentialSolns = new Vector<CalculatedGraph>();
				File exponentialDir = new File(baseDir.get(i).getAbsolutePath()+"\\exponential");
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
						solnToExamine.UpdateWithinRInfo(.16);
//						solnToExamine.setCostBasis(costBasis);
//						
//					    //solnToExamine.UpdateCalcuations();
//						solnToExamine.UpdatePSOCalculations();
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
				
				
				Vector<CalculatedGraph> heur1Solns = new Vector<CalculatedGraph>();
				File heur1Dir = new File(baseDir.get(i).getAbsolutePath()+"\\optimizedheuristic");
				File heur1File[] = heur1Dir.listFiles(filter);
				
				cost = -1;
				costExp =null;
				
				System.out.println("Loading Heur1"+ heur1Dir.toString());
				for(int j =0;  j < heur1File.length; j++)
				//j < 5; j++)
				//j < exponentialFile.length; j++)
				{

					CalculatedGraph solnToExamine = 
						Util.streaminCalculatedGraph(heur1File[j]);
					//Make sure all calculations are up to date
					if(DO_UPDATE)
					{
						solnToExamine.UpdateWithinRInfo(.16);
//						solnToExamine.setCostBasis(costBasis);
//						
//					    //solnToExamine.UpdateCalcuations();
//						solnToExamine.UpdatePSOCalculations();
					}
					heur1Solns.add(solnToExamine);
					heur1All.add(solnToExamine);

					System.out.println(heur1File[j].getName() +" has fitness of "+solnToExamine.getFitnessValue()+" with # edges "+solnToExamine.getTotalEdges());
				//	System.out.println("    cost:"+solnToExamine.getCost()+" AISPL:"+solnToExamine.getAISPL()+" pLCC:"+solnToExamine.getLCC());
				//	System.out.println(exponentialFile[j].getName() +" has cost of "+solnToExamine.getCost());
					if(solnToExamine.getCost() > cost)
					{
						cost = solnToExamine.getCost();
						costExp = heur1File[j].getName();
					}
					
					
					if(DO_WRITE_OUT_UPDATED_GRAPHS)
					{
						Util.streamOutCalculatedGraphWithExtension(heur1File[j].getParentFile(),
								heur1File[j].getName(),
								solnToExamine);
					}
					
				}
				System.out.println("Done Loading Heur1");
				
					
				Vector<CalculatedGraph> heur2Solns = new Vector<CalculatedGraph>();
				File heur2Dir = new File(baseDir.get(i).getAbsolutePath()+"\\optimizedheuristic1b");
				File heur2File[] = heur2Dir.listFiles(filter);
				

				
				System.out.println("Loading Heur2"+ heur2Dir.toString());
				for(int j =0;  j < heur2File.length; j++)
				//j < 5; j++)
				//j < exponentialFile.length; j++)
				{

					CalculatedGraph solnToExamine = 
						Util.streaminCalculatedGraph(heur2File[j]);
					//Make sure all calculations are up to date
					if(DO_UPDATE)
					{
						solnToExamine.UpdateWithinRInfo(.16);
//						solnToExamine.setCostBasis(costBasis);
//						
//					    //solnToExamine.UpdateCalcuations();
//						solnToExamine.UpdatePSOCalculations();
					}
					heur2Solns.add(solnToExamine);
					heur2All.add(solnToExamine);

					System.out.println(heur2File[j].getName() +" has fitness of "+solnToExamine.getFitnessValue()+" with # edges "+solnToExamine.getTotalEdges());
				//	System.out.println("    cost:"+solnToExamine.getCost()+" AISPL:"+solnToExamine.getAISPL()+" pLCC:"+solnToExamine.getLCC());
				//	System.out.println(exponentialFile[j].getName() +" has cost of "+solnToExamine.getCost());
					if(solnToExamine.getCost() > cost)
					{
						cost = solnToExamine.getCost();
						costExp = heur2File[j].getName();
					}
					
					
					if(DO_WRITE_OUT_UPDATED_GRAPHS)
					{
						Util.streamOutCalculatedGraphWithExtension(heur2File[j].getParentFile(),
								heur2File[j].getName(),
								solnToExamine);
					}
					
				}
				System.out.println("Done Loading Heur2");
				System.out.println("Most expensive network - "+costExp+" with cost of "+cost);
				
			
				Vector<CalculatedGraph> heuristicSolns = new Vector<CalculatedGraph>();
				File heuristic3Dir = new File(baseDir.get(i).getAbsolutePath()+"\\optimizedheuristic1b");
				File heuristic3File[] = heuristic3Dir.listFiles(filter);
				
				System.out.println("Loading Heur3 "+ heuristic3Dir.toString());
				for(int j =0;   j < heuristic3File.length; j++)
				{

					CalculatedGraph solnToExamine = 
						Util.streaminCalculatedGraph(heuristic3File[j]);
					//Make sure all calculations are up to date
					if(DO_UPDATE)
					{
						solnToExamine.UpdateWithinRInfo(.16);
//						solnToExamine.setCostBasis(costBasis);
//						
//					    //solnToExamine.UpdateCalcuations();
//						solnToExamine.UpdatePSOCalculations();
					}
				
					
					heuristicSolns.add(solnToExamine);
					heur3All.add(solnToExamine);

					System.out.println(heuristic3File[j].getName() +" has fitness of "+solnToExamine.getFitnessValue()+" with # edges "+solnToExamine.getTotalEdges());

					if(DO_WRITE_OUT_UPDATED_GRAPHS)
					{
						Util.streamOutCalculatedGraphWithExtension(heuristic3File[j].getParentFile(),
								heuristic3File[j].getName(),
								solnToExamine);
					}
					
				}
				System.out.println("Done Loading Hueristic");
				
				
				Vector<CalculatedGraph> heuristicSolns2 = new Vector<CalculatedGraph>();
				File heuristic4Dir = new File(baseDir.get(i).getAbsolutePath()+"\\optimizedheuristic1b");
				File heuristic4File[] = heuristic4Dir.listFiles(filter);
				
				System.out.println("Loading Heur4 - "+ heuristic4Dir.toString());
				for(int j =0;   j < heuristic4File.length; j++)
				{

					CalculatedGraph solnToExamine = 
						Util.streaminCalculatedGraph(heuristic4File[j]);
					//Make sure all calculations are up to date
					if(DO_UPDATE)
					{
						solnToExamine.UpdateWithinRInfo(.16);
//						solnToExamine.setCostBasis(costBasis);
//						
//					    //solnToExamine.UpdateCalcuations();
//						solnToExamine.UpdatePSOCalculations();
					}
				
					
					heuristicSolns.add(solnToExamine);
					heur4All.add(solnToExamine);

					System.out.println(heuristic4File[j].getName() +" has fitness of "+solnToExamine.getFitnessValue()+" with # edges "+solnToExamine.getTotalEdges());

					if(DO_WRITE_OUT_UPDATED_GRAPHS)
					{
						Util.streamOutCalculatedGraphWithExtension(heuristic4File[j].getParentFile(),
								heuristic4File[j].getName(),
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
						new ValuePlotWriter(optimalSolns, heur2Solns, heuristicSolns, baseDir.get(i));
					valuePlotWriter.generatePlots();
					System.out.println("Done Performing value Plots");
				}
				
				if(DO_DIST == true)
				{
					System.out.println("Performing dist Plots");
					DistPlotWriter distPlotWriter = 
						new DistPlotWriter(optimalSolns, heur2Solns, heuristicSolns, baseDir.get(i));
					distPlotWriter.generatePlots();
					System.out.println("Done Performing dist Plots");
				}
			
				if(DO_NODE_VALUES == true)
				{
					System.out.println("Performing node value Plots");
					NodeValuePlotWriter nodeValuePlotWriter = 
						new NodeValuePlotWriter(optimalSolns, heur2Solns, heuristicSolns, baseDir.get(i));
					nodeValuePlotWriter.generatePlots();
					System.out.println("Done Performing node value Plots");
				}
				

				

			
			}
		
			
		if(DO_COMBINED_PLOTS==true)
		{

//				System.out.println("Performing value Plots");
//				ValuePlotWriter valuePlotWriter = 
//					new ValuePlotWriter(optSolnsAll, exponentialSolnsAll,  heur1All, heur2All, heur3All, heur4All, baseDir.get(0).getParentFile());
//				valuePlotWriter.generatePlotsH();
//				System.out.println("Done Performing value Plots");
//
//			System.out.println("Performing value Plots");
			ValuePlotWriter valuePlotWriter = 
				new ValuePlotWriter(optSolnsAll, exponentialSolnsAll, heur1All, heur2All, heur3All, heur4All, baseDir.get(0).getParentFile());
			valuePlotWriter.generatePlots2();
//			valuePlotWriter.generateAISPLPlotVsAssortativity();
//			System.out.println("Done Performing value Plots");
//			
//
//				System.out.println("Performing dist Plots");
//				DistPlotWriter distPlotWriter = 
//						new DistPlotWriter(optSolnsAll, exponentialSolnsAll, heur1All, heur2All, heur3All, heur4All, baseDir.get(0).getParentFile());
//				distPlotWriter.generatePlotsH();
//				System.out.println("Done Performing dist Plots");
				
				
				
//				System.out.println("Performing optimized value Plots");
//				PLDistPlotWriter optimizedValuePlotWriter =  
//					new PLDistPlotWriter(optSolnsAll, exponentialSolnsAll, heur1All, heur2All, heur3All, heur4All, baseDir.get(0).getParentFile());
//				optimizedValuePlotWriter.generatePlotsH();
//				System.out.println("Done Performing optimized value Plots");
//
//
//				System.out.println("Performing node value Plots");
//				NodeValuePlotWriter nodeValuePlotWriter = 
//					new NodeValuePlotWriter(optSolnsAll, exponentialSolnsAll, heur1All, baseDir.get(0).getParentFile());
//				nodeValuePlotWriter.generatePlotsH();
//				System.out.println("Done Performing node value Plots");

			
			System.out.println("Performing node value Plots");
			NodeValuePlotWriter nodeValuePlotWriter = 
				new NodeValuePlotWriter(optSolnsAll, exponentialSolnsAll, heur1All, baseDir.get(0).getParentFile());
			nodeValuePlotWriter.plotConnectedNodesR1R2VsNodesWithinR1R2LargeConnected();
			
//			nodeValuePlotWriter.plotConnectedNodesR1R2VsNodesWithinR1R2Large();
//			nodeValuePlotWriter.plotConnectedNodesR1R2VsNodesWithinR1R2Connected();
//			nodeValuePlotWriter.plotConnectedNodesR1R2VsNodesWithinR1R2LargeConnected();
//			nodeValuePlotWriter.plotConnectedNodesR1R2VsNodesWithinR1R2Largea();
//			nodeValuePlotWriter.plotConnectedNodesR1R2VsNodesWithinR1R2LargeConnecteda();			
			

			
			System.out.println("Done Performing node value Plots");
			
			

		}
			
			

		}
		catch(Exception e)
		{
			e.printStackTrace();
		}

	}
	
}
