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
package psograph.analysis;

import java.io.File;
import java.util.TreeMap;
import java.util.Vector;

import psograph.graph.*;
import psograph.util.MFilePlotData;
import psograph.util.PlotPairDoubleData;
import psograph.util.ValueWriter;



public class NodeValuePlotWriter {

	File outputDirectory;
	
	Vector<CalculatedGraph> optimized;
	Vector<CalculatedGraph> exponential;
	Vector<CalculatedGraph> heuristic;
	
	Vector<CalculatedGraph> heur1;
	Vector<CalculatedGraph> heur2;
	Vector<CalculatedGraph> heur3;
	
	String layoutTitle;
	
	
	public NodeValuePlotWriter(Vector<CalculatedGraph> optimized, Vector<CalculatedGraph> exponential,
			Vector<CalculatedGraph> heuristic, File OutDir)
	{
		this.optimized = optimized;
		this.exponential = exponential;
		this.heuristic = heuristic;
		outputDirectory = OutDir;
	}

	public NodeValuePlotWriter(Vector<CalculatedGraph> optSolnsAll,
			Vector<CalculatedGraph> exponentialSolnsAll,
			Vector<CalculatedGraph> heur1All, Vector<CalculatedGraph> heur2All,
			Vector<CalculatedGraph> heur3All, Vector<CalculatedGraph> heur4All,
			File parentFile) 
	{
		this.optimized =optSolnsAll;
		this.exponential = exponentialSolnsAll;
		this.heur1 =heur1All;
		this.heur2 =heur2All;
		this.heur3 =heur3All;
		this.heuristic = heur4All;
		
		outputDirectory = parentFile;
	}

	public void generatePlots2() throws Exception
	{
		plotNodeDegreeVsEdgeLength();
		plotNodeDegreeVsEdgeLength2();
		
		
	}
	
	private void plotNodeDegreeVsEdgeLength2() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "Node Degree","Edge Length", "Node Degree vs Edge Length",
				"NodeDegreeVsEdgeLength2.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Vector<Double>> edgeLenghtById = optimized.get(i).getEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				int NodeId = nodeDistById.get(j);
				Vector<Double> lengths = edgeLenghtById.get(j);
				
				if(NodeId != lengths.size())
				{
					throw new Exception (" degree doesn't match size of neighbors");
					
				}
				double val = 0;
				for(int k = 0; k < lengths.size() ; k++)
				{
					//System.out.println(NodeId +", "+ lengths.get(k) );
					val += 	lengths.get(k) ;
				
				}
				
				optimizedData.add(new PlotPairDoubleData(
						NodeId ,
						val/lengths.size() ));
				
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < exponential.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
//			TreeMap<Integer, Vector<Double>> edgeLenghtById = exponential.get(i).getEdgeLengthById();
//			for(int j=0; j < exponential.get(i).getN(); j++)
//			{
//				int NodeId = nodeDistById.get(j);
//				Vector<Double> lengths = edgeLenghtById.get(j);
//				
//				for(int k = 0; k < lengths.size() ; k++)
//				{
//					exponentialData.add(new PlotPairDoubleData(
//							NodeId ,
//							lengths.get(k) ));
//				
//				}
//			}
//		}
//		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < heuristic.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
//			TreeMap<Integer, Vector<Double>> edgeLenghtById = heuristic.get(i).getEdgeLengthById();
//			for(int j=0; j < heuristic.get(i).getN(); j++)
//			{
//				int NodeId = nodeDistById.get(j);
//				Vector<Double> lengths = edgeLenghtById.get(j);
//				
//				for(int k = 0; k < lengths.size() ; k++)
//				{
//					heuristicData.add(new PlotPairDoubleData(
//							NodeId ,
//							lengths.get(k) ));
//				
//				}
//			}
//		}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	

	private void plotNodeDegreeVsEdgeLength() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "Node Degree","Edge Length", "Node Degree vs Edge Length",
				"NodeDegreeVsEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Vector<Double>> edgeLenghtById = optimized.get(i).getEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				int NodeId = nodeDistById.get(j);
				Vector<Double> lengths = edgeLenghtById.get(j);
				
				if(NodeId != lengths.size())
				{
					throw new Exception (" degree doesn't match size of neighbors");
					
				}
				
				for(int k = 0; k < lengths.size() ; k++)
				{
					//System.out.println(NodeId +", "+ lengths.get(k) );
					optimizedData.add(new PlotPairDoubleData(
							NodeId ,
							lengths.get(k) ));
				
				}
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < exponential.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
//			TreeMap<Integer, Vector<Double>> edgeLenghtById = exponential.get(i).getEdgeLengthById();
//			for(int j=0; j < exponential.get(i).getN(); j++)
//			{
//				int NodeId = nodeDistById.get(j);
//				Vector<Double> lengths = edgeLenghtById.get(j);
//				
//				for(int k = 0; k < lengths.size() ; k++)
//				{
//					exponentialData.add(new PlotPairDoubleData(
//							NodeId ,
//							lengths.get(k) ));
//				
//				}
//			}
//		}
//		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < heuristic.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
//			TreeMap<Integer, Vector<Double>> edgeLenghtById = heuristic.get(i).getEdgeLengthById();
//			for(int j=0; j < heuristic.get(i).getN(); j++)
//			{
//				int NodeId = nodeDistById.get(j);
//				Vector<Double> lengths = edgeLenghtById.get(j);
//				
//				for(int k = 0; k < lengths.size() ; k++)
//				{
//					heuristicData.add(new PlotPairDoubleData(
//							NodeId ,
//							lengths.get(k) ));
//				
//				}
//			}
//		}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	public void generatePlots() throws Exception
	{
		
		plotNodeDegreeVsTotalEdgeLength();
		plotNodeDegreeVsTotalEdgeCost();
		plotNodeDegreeVsMaximumEdgeLength();
		plotNodeDegreeVsMaximumEdgeCost();
		plotNodeDegreeVsConnectedNodesWithinRadiusR();
		plotNodeDegreeVsMeanDegreeOfConnectedNeighbors();
		plotNodeDegreeVsBetweeness();
		plotNodeDegreeVsCloseness();
		plotNodeDegreeVsMeanEdgeLength();
		plotNodeDegreeVsMeanEdgeCost();
		
		
		plotTotalEdgeLengthVsTotalEdgeCost();
		plotTotalEdgeLengthVsMaximumEdgeLength();
		plotTotalEdgeLengthVsMaximumEdgeCost();
		plotTotalEdgeLengthVsConnectedNodesWithinRadiusR();
		plotTotalEdgeLengthVsMeanDegreeOfConnectedNeighbors();
		plotTotalEdgeLengthVsBetweeness();
		plotTotalEdgeLengthVsCloseness();
		plotTotalEdgeLengthVsMeanEdgeLength();
		plotTotalEdgeLengthVsMeanEdgeCost();
		
		plotTotalEdgeCostVsMaximumEdgeLength();
		plotTotalEdgeCostVsMaximumEdgeCost();
		plotTotalEdgeCostVsConnectedNodesWithinRadiusR();
		plotTotalEdgeCostVsMeanDegreeOfConnectedNeighbors();
		plotTotalEdgeCostVsBetweeness();
		plotTotalEdgeCostVsCloseness();
		plotTotalEdgeCostVsMeanEdgeLength();
		plotTotalEdgeCostVsMeanEdgeCost();
		
		plotMaximumEdgeLengthVsMaximumEdgeCost();
		plotMaximumEdgeLengthVsConnectedNodesWithinRadiusR();
		plotMaximumEdgeLengthVsMeanDegreeOfConnectedNeighbors();
		plotMaximumEdgeLengthVsBetweeness();
		plotMaximumEdgeLengthVsCloseness();
		plotMaximumEdgeLengthVsMeanEdgeLength();
		plotMaximumEdgeLengthVsMeanEdgeCost();
		
		plotMaximumEdgeCostVsConnectedNodesWithinRadiusR();
		plotMaximumEdgeCostVsMeanDegreeOfConnectedNeighbors();
		plotMaximumEdgeCostVsBetweeness();
		plotMaximumEdgeCostVsCloseness();
		plotMaximumEdgeCostVsMeanEdgeLength();
		plotMaximumEdgeCostVsMeanEdgeCost();
		
		plotConnectedNodesWithinRadiusRVsMeanDegreeOfConnectedNeighbors();
		plotConnectedNodesWithinRadiusRVsBetweeness();
		plotConnectedNodesWithinRadiusRVsCloseness();
		plotConnectedNodesVsMeanEdgeLength();
		plotConnectedNodesVsMeanEdgeCost();
		plotConnectedNodesVsNodesWithinR();
		
		
		plotMeanDegreeOfConnectedNeighborsRVsBetweeness();
		plotMeanDegreeOfConnectedNeighborsRVsCloseness();
		plotMeanDegreeOfConnectedNeighborsVsMeanEdgeLength();
		plotMeanDegreeOfConnectedNeighborsVsMeanEdgeCost();

		
		plotBetweenessRVsCloseness();
		plotBetweenessVsMeanEdgeLength();
		plotBetweenessVsMeanEdgeCost();
		

		plotClosenessVsMeanEdgeLength();
		plotClosenessVsMeanEdgeCost();


		plotMeanEdgeLengthVsMeanEdgeCost();
		
		
		
		
		
	}

	private void plotConnectedNodesVsNodesWithinR() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"Nodes Within Radius="+GraphConstants.MeasurementR,
				"Connected Nodes Within Radius="+GraphConstants.MeasurementR+ 
				" vs Nodes Within Radius="+GraphConstants.MeasurementR,
				"ConnectedNodesWithinRadiusRVsNodesWithinR.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getNodesWithinRadiusById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getNodesWithinRadiusById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getNodesWithinRadiusById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	
//	public void plotConnectedNodesR1R2VsNodesWithinR1R2() throws Exception
//	{
//		MFilePlotData plotData = new MFilePlotData( "Connected Nodes Greater than Radius1=.16 but less then Radius2=.5",
//				"Nodes Greater than Radius1=.16 but less then Radius2=.5",
//				"",
//				"ConnectedNodesWithinRadiusR12VsNodesWithinR12.m", outputDirectory, null, 20);
//		
//		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < optimized.size() ; i++)
//		{
//			
//			
//			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < optimized.get(i).getN(); j++)
//			{
//				optimizedData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}
//		}
//		
//		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
////		for(int i =0 ; i < exponential.size() ; i++)
////		{
////			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusR12ById();
////			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getNodesWithinRadius12ById();
////			for(int j=0; j < exponential.get(i).getN(); j++)
////			{
////				exponentialData.add(new PlotPairDoubleData(
////						nodeDistById.get(j) ,
////						edgeLenghtById.get(j) ));
////			}}
//		
//		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
////		for(int i =0 ; i < heuristic.size() ; i++)
////		{
////			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusR12ById();
////			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getNodesWithinRadius12ById();
////			for(int j=0; j < heuristic.get(i).getN(); j++)
////			{
////				heuristicData.add(new PlotPairDoubleData(
////						nodeDistById.get(j) ,
////						edgeLenghtById.get(j) ));
////			}	
////			}
//		
//		ValueWriter plotter = new ValueWriter(plotData);
//		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
//		
//	}
	

	public void plotConnectedNodesR1R2VsNodesWithinR1R2() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "Connected Nodes Greater than Radius1=.16 but less then Radius2=.5",
				"Nodes Greater than Radius1=.16 but less then Radius2=.5",
				"",
				"ConnectedNodesWRR12VsNodesWithinR12.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			
			
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesWithinRadiusR12ById();
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getNodesWithinRadius12ById();
			
			TreeMap<Integer,Integer> specialDegreeById = optimized.get(i).getSpecialDegreeById();
			
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				double percentage = (double)nodeDistById.get(j) / (double)edgeLenghtById.get(j);
				
				double specialDegree = specialDegreeById.get(j); 
				
				optimizedData.add(new PlotPairDoubleData(
						percentage ,
						specialDegree ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < exponential.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < exponential.get(i).getN(); j++)
//			{
//				exponentialData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < heuristic.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < heuristic.get(i).getN(); j++)
//			{
//				heuristicData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}	
//			}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	
	public void plotConnectedNodesR1R2VsNodesWithinR1R2Connected() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "Connected Nodes Greater than Radius1=.16 but less then Radius2=.5",
				"Nodes Greater than Radius1=.16 but less then Radius2=.5",
				"",
				"ConnectedNodesWRR12VsNodesWithinR12Connected.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			
			
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesWithinRadiusR12ById();
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getNodesWithinRadius12ById();
			
			TreeMap<Integer,Integer> specialDegreeById = optimized.get(i).getSpecialDegreeById();
			
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				double percentage = (double)nodeDistById.get(j) ;
				
				double specialDegree = specialDegreeById.get(j); 
				
				optimizedData.add(new PlotPairDoubleData(
						percentage ,
						specialDegree ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < exponential.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < exponential.get(i).getN(); j++)
//			{
//				exponentialData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < heuristic.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < heuristic.get(i).getN(); j++)
//			{
//				heuristicData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}	
//			}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	
	public void plotConnectedNodesR1R2VsNodesWithinR1R2Large() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "Connected Nodes Greater than Radius .5",
				"Nodes Greater than Radius .5",
				"",
				"ConnectedNodesWRR12VsNodesWithinR12Large.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			
			
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesOutsideRadiusRById(.5);
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getNodesOutsideRadius(.5);
			
			TreeMap<Integer,Integer> specialDegreeById = optimized.get(i).getSpecialDegreeById();
			
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				double percentage = (double)nodeDistById.get(j) / (double)edgeLenghtById.get(j);
				
				double specialDegree = specialDegreeById.get(j); 
				
				optimizedData.add(new PlotPairDoubleData(
						percentage ,
						specialDegree ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < exponential.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < exponential.get(i).getN(); j++)
//			{
//				exponentialData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < heuristic.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < heuristic.get(i).getN(); j++)
//			{
//				heuristicData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}	
//			}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	
	public void plotConnectedNodesR1R2VsNodesWithinR1R2LargeConnected() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "",
				"",
				"",
				"special2.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			
			
			TreeMap<Integer, Vector<Double>> nodeDistById = optimized.get(i).pat2();
			
			TreeMap<Integer,Integer> specialDegreeById = optimized.get(i).getSpecialDegreeById();
			
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{

				
				double specialDegree = specialDegreeById.get(j); 
				
				Vector<Double> vec = nodeDistById.get(j);
				
				for(int patty=0 ; patty < vec.size() && vec.size()!=0; patty++)
				{
					optimizedData.add(new PlotPairDoubleData(
							specialDegree, (double)vec.get(patty) ));
				}
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Vector<Double>> nodeDistById = exponential.get(i).pat2();
			
			TreeMap<Integer,Integer> specialDegreeById = exponential.get(i).getSpecialDegreeById();
			
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{

				
				double specialDegree = specialDegreeById.get(j); 
				
				Vector<Double> vec = nodeDistById.get(j);
				
				for(int patty=0 ; patty < vec.size() && vec.size()!=0; patty++)
				{
					exponentialData.add(new PlotPairDoubleData(
							specialDegree, (double)vec.get(patty) ));
				}
			}
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < heuristic.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < heuristic.get(i).getN(); j++)
//			{
//				heuristicData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}	
//			}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	
	public void plotConnectedNodesR1R2VsNodesWithinR1R2Largea() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "Connected Nodes Greater than Radius .5",
				"Nodes Greater than Radius .5",
				"",
				"ConnectedNodesWRR12VsNodesWithinR12Largea.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			
			
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesOutsideRadiusRById(.5);
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getNodesOutsideRadius(.5);
			
			TreeMap<Integer,Integer> specialDegreeById = optimized.get(i).getSpecialDegreeById(.5);
			
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				double percentage = (double)nodeDistById.get(j) / (double)edgeLenghtById.get(j);
				
				double specialDegree = specialDegreeById.get(j); 
				
				optimizedData.add(new PlotPairDoubleData(
						percentage ,
						specialDegree ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < exponential.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < exponential.get(i).getN(); j++)
//			{
//				exponentialData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < heuristic.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < heuristic.get(i).getN(); j++)
//			{
//				heuristicData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}	
//			}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	
	public void plotConnectedNodesR1R2VsNodesWithinR1R2LargeConnecteda() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "Node Degree",
				"Connected Nodes Greater than Radius1 .5",
				"",
				"ConnectedNodesWRR12VsNodesWithinR12LargeConnecteda.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			
			
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesOutsideRadiusRById(.5);
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getNodesOutsideRadius(.5);
			
			TreeMap<Integer,Integer> specialDegreeById = optimized.get(i).getSpecialDegreeById(.5);
			
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{

				
				double specialDegree = specialDegreeById.get(j); 
				
				optimizedData.add(new PlotPairDoubleData(
						(double)nodeDistById.get(j) ,
						specialDegree ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < exponential.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < exponential.get(i).getN(); j++)
//			{
//				exponentialData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < heuristic.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < heuristic.get(i).getN(); j++)
//			{
//				heuristicData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}	
//			}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	
	
	private void plotBetweenessRVsCloseness() throws Exception{
		
		MFilePlotData plotData = new MFilePlotData( "Betweeness",
				"Closeness",
				"Betweeness vs Closeness",
				"BetweenessVsCloseness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getBetweenessById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getClosenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getBetweenessById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getClosenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getBetweenessById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getClosenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotBetweenessVsMeanEdgeCost() throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Betweeness",
				"Mean Edge Cost",
				"Betweeness vs Mean Edge Cost",
				"BetweenessVsMeanEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getBetweenessById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getBetweenessById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getBetweenessById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotBetweenessVsMeanEdgeLength()throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "Betweeness",
				"Mean Edge Length",
				"Betweeness vs Mean Edge Length",
				"BetweenessVsMeanEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getBetweenessById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getBetweenessById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getBetweenessById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotClosenessVsMeanEdgeCost()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Closeness",
				"Mean Edge Cost",
				"Closeness vs Mean Edge Cost",
				"ClosenessVsMeanEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getClosenessById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getClosenessById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getClosenessById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotClosenessVsMeanEdgeLength() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Closeness",
				"Mean Edge Length",
				"Closeness vs Mean Edge Length",
				"ClosenessVsMeanEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getClosenessById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getClosenessById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getClosenessById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotConnectedNodesVsMeanEdgeCost()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"Mean Edge Cost",
				"Connected Nodes Within Radius="+GraphConstants.MeasurementR+ 
				" vs Mean Edge Cost",
				"ConnectedNodesWithinRadiusRVsMeanEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotConnectedNodesVsMeanEdgeLength() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"Mean Edge Length",
				"Connected Nodes Within Radius="+GraphConstants.MeasurementR+ 
				" vs Mean Edge Length",
				"ConnectedNodesWithinRadiusRVsMeanEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotConnectedNodesWithinRadiusRVsBetweeness() throws Exception{
		
		
		MFilePlotData plotData = new MFilePlotData( "Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"Betweeness",
				"Connected Nodes Within Radius="+GraphConstants.MeasurementR+ 
				" vs Betweeness",
				"ConnectedNodesWithinRadiusRVsBetweenessVsBetweeness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getBetweenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getBetweenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getBetweenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotConnectedNodesWithinRadiusRVsCloseness() throws Exception{
		
		MFilePlotData plotData = new MFilePlotData( "Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"Closeness",
				"Connected Nodes Within Radius="+GraphConstants.MeasurementR+ 
				" vs Closeness",
				"ConnectedNodesWithinRadiusRVsCloseness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getClosenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getClosenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getClosenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotConnectedNodesWithinRadiusRVsMeanDegreeOfConnectedNeighbors()throws Exception {

		MFilePlotData plotData = new MFilePlotData( "Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"Mean Degree Of Connected Neighbors",
				"Connected Nodes Within Radius="+GraphConstants.MeasurementR+ " vs Mean Degree Of Connected Neighbors",
				"ConnectedNodesWithinRadiusRVsMeanDegreeOfConnectedNeighbors.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusRById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeCostVsBetweeness() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Cost",
				"Betweeness",
				"Maximum Edge Cost vs Betweeness",
				"MaximumEdgeCostVsBetweeness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getBetweenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getBetweenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getBetweenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeCostVsCloseness()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Cost",
				"Closeness",
				"Maximum Edge Cost vs Closeness",
				"MaximumEdgeCostVsCloseness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getClosenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getClosenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getClosenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeCostVsConnectedNodesWithinRadiusR()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Cost",
				"Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"Maximum Edge Cost vs Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"MaximumEdgeCostVsConnectedNodesWithinRadiusR.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeCostVsMeanDegreeOfConnectedNeighbors() throws Exception{

		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Cost",
				"Mean Degree Of Connected Neighbors",
				"Maximum Edge Cost vs Mean Degree Of Connected Neighbors",
				"MaximumEdgeCostVsMeanDegreeOfConnectedNeighbors.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeCostVsMeanEdgeCost() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Cost",
				"Mean Edge Cost ",
				"Maximum Edge Cost vs Mean Edge Cost",
				"MaximumEdgeCostVsMeanEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeCostVsMeanEdgeLength()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Cost",
				"Mean Edge Length ",
				"Maximum Edge Cost vs Mean Edge Length",
				"MaximumEdgeCostVsMeanEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeLengthVsBetweeness()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Length",
				"Betweeness",
				"Maximum Edge Length vs Betweeness",
				"MaximumEdgeLengthVsBetweeness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getBetweenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getBetweenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getBetweenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeLengthVsCloseness() throws Exception{
		
		
		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Length",
				"Closeness",
				"Maximum Edge Length vs Closeness",
				"MaximumEdgeLengthVsCloseness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getClosenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getClosenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getClosenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeLengthVsConnectedNodesWithinRadiusR()throws Exception {
		
		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Length",
				"Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"Maximum Edge Length vs Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"MaximumEdgeLengthVsConnectedNodesWithinRadiusR.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeLengthVsMaximumEdgeCost()throws Exception 
	{
		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Length",
				"Maximum Edge Cost",
				"Maximum Edge Length vs Maximum Edge Cost",
				"MaximumEdgeLengthVsMaximumEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMaximumEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMaximumEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMaximumEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeLengthVsMeanDegreeOfConnectedNeighbors()throws Exception 
	{
	
		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Length",
				"Mean Degree Of Connected Neighbors",
				"Maximum Edge Length vs Mean Degree Of Connected Neighbors",
				"MaximumEdgeLengthVsMeanDegreeOfConnectedNeighbors.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeLengthVsMeanEdgeCost()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Length",
				"Mean Edge Cost",
				"Maximum Edge Length vs Mean Edge Cost",
				"MaximumEdgeLengthVsMeanEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMaximumEdgeLengthVsMeanEdgeLength() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Maximum Edge Length",
				"Mean Edge Length",
				"Maximum Edge Length vs Mean Edge Length",
				"MaximumEdgeLengthVsMeanEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMaximumEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMeanDegreeOfConnectedNeighborsRVsBetweeness() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Mean Degree Of Connected Neighbors",
				"Betweeness",
				"Mean Degree Of Connected Neighbors vs Betweeness",
				"MeanDegreeOfConnectedNeighborsVsBetweeness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMeanDegreeOfConnectedNeighborsById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getBetweenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMeanDegreeOfConnectedNeighborsById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getBetweenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMeanDegreeOfConnectedNeighborsById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getBetweenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMeanDegreeOfConnectedNeighborsRVsCloseness() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Mean Degree Of Connected Neighbors",
				"Closeness",
				"Mean Degree Of Connected Neighbors vs Closeness",
				"MeanDegreeOfConnectedNeighborsVsCloseness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMeanDegreeOfConnectedNeighborsById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getClosenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMeanDegreeOfConnectedNeighborsById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getClosenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMeanDegreeOfConnectedNeighborsById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getClosenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMeanDegreeOfConnectedNeighborsVsMeanEdgeCost()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Mean Degree Of Connected Neighbors",
				"Mean Edge Cost",
				"Mean Degree Of Connected Neighbors vs Mean Edge Cost",
				"MeanDegreeOfConnectedNeighborsVsMeanEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMeanDegreeOfConnectedNeighborsById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMeanDegreeOfConnectedNeighborsById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMeanDegreeOfConnectedNeighborsById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotMeanDegreeOfConnectedNeighborsVsMeanEdgeLength()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Mean Degree Of Connected Neighbors",
				"Mean Edge Length",
				"Mean Degree Of Connected Neighbors vs Mean Edge Length",
				"MeanDegreeOfConnectedNeighborsVsMeanEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMeanDegreeOfConnectedNeighborsById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMeanDegreeOfConnectedNeighborsById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMeanDegreeOfConnectedNeighborsById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}


	private void plotMeanEdgeLengthVsMeanEdgeCost()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Mean Edge Length",
				"Mean Edge Cost",
				"Mean Edge Length vs Mean Edge Cost",
				"MeanEdgeLengthVsCMeanEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getMeanEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getMeanEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getMeanEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotNodeDegreeVsBetweeness() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Node Degree","Betweeness", "Node Degree vs Betweeness",
				"NodeDegreeVsBetweeness.m", outputDirectory, null, 20);
		
		
		System.out.println("How many graphs - "+optimized.size());
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getBetweenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getBetweenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getBetweenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotNodeDegreeVsCloseness() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Node Degree","Closeness", "Node Degree vs Closeness",
				"NodeDegreeVsCloseness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getClosenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getClosenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getClosenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotNodeDegreeVsConnectedNodesWithinRadiusR() throws Exception 
	{
		MFilePlotData plotData = new MFilePlotData( "Node Degree",
				"Connected Nodes Within Radius="+GraphConstants.MeasurementR, 
				"Node Degree vs Connected Nodes Within Radius R="+GraphConstants.MeasurementR,
				"NodeDegreeVsConnectedNodesWithinRadiusR.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	
	public void plotNodeDegreeVsConnectedNodesOutsideRadiusR() throws Exception 
	{
		MFilePlotData plotData = new MFilePlotData( "Node Degree",
				"Connected Nodes Outside Radius="+GraphConstants.MeasurementR, 
				"Node Degree vs Connected Nodes Outside Radius R="+GraphConstants.MeasurementR,
				"NodeDegreeVsConnectedNodesOutsideRadiusR.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimized.get(i).UpdateWithinRInfo(.2);
			
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getConnectedNodesOutsideRadiusRById(.2);
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponential.get(i).UpdateWithinRInfo(.2);
			
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getConnectedNodesOutsideRadiusRById(.2);
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristic.get(i).UpdateWithinRInfo(.2);
			
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getConnectedNodesOutsideRadiusRById(.2);
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotNodeDegreeVsMaximumEdgeCost() throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Node Degree","Maximum Edge Cost", "Node Degree vs Maximum Edge Cost",
				"NodeDegreeVsMaximumEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMaximumEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMaximumEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMaximumEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotNodeDegreeVsMaximumEdgeLength() throws Exception {
		
		MFilePlotData plotData = new MFilePlotData( "Node Degree","Maximum Edge Length", "Node Degree vs Maximum Edge Length",
				"NodeDegreeVsMaximumEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMaximumEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMaximumEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMaximumEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotNodeDegreeVsMeanDegreeOfConnectedNeighbors() throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Node Degree"," Mean Degree Of Connected Neighbors", "Node Degree vs Mean Degree Of Connected Neighbors",
				"NodeDegreeVsMeanDegreeOfConnectedNeighbors.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}


	private void plotNodeDegreeVsMeanEdgeCost()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Node Degree",
				"Mean Edge Cost",
				"Node Degree vs Mean Edge Cost",
				"NodeDegreeVsMeanEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotNodeDegreeVsMeanEdgeLength() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Node Degree",
				"Mean Edge Length",
				"Node Degree vs Mean Edge Length",
				"NodeDegreeVsMeanEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
		
	}

	private void plotNodeDegreeVsTotalEdgeCost() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Node Degree","Total Edge Cost", "Node Degree vs Total Edge Cost",
				"NodeDegreeVsTotalEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getTotalEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getTotalEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getTotalEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotNodeDegreeVsTotalEdgeLength() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "Node Degree","Total Edge Length", "Node Degree vs Total Edge Length",
				"NodeDegreeVsTotalEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getTotalEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getTotalEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getNodeDistributionId();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getTotalEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeCostVsBetweeness()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Total Edge Cost",
				"Betweeness",
				"Total Edge Cost vs Betweeness",
				"TotalEdgeCostVsBetweeness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getBetweenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getBetweenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getBetweenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeCostVsCloseness() throws Exception{
		
		MFilePlotData plotData = new MFilePlotData( "Total Edge Cost",
				"Closeness",
				"Total Edge Cost vs Closeness",
				"TotalEdgeCostVsCloseness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getClosenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getClosenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getClosenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeCostVsConnectedNodesWithinRadiusR()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Total Edge Cost",
				"Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"Total Edge Cost vs Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"TotalEdgeCostVsConnectedNodesWithinRadiusR.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeCostVsMaximumEdgeCost()throws Exception {

		MFilePlotData plotData = new MFilePlotData( "Total Edge Cost",
				"Maximum Edge Cost",
				"Total Edge Cost vs Maximum Edge Cost",
				"TotalEdgeCostVsMaximumEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMaximumEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMaximumEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMaximumEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeCostVsMaximumEdgeLength()throws Exception {
		
		MFilePlotData plotData = new MFilePlotData( "Total Edge Cost",
				"Maximum Edge Length",
				"Total Edge Cost vs Maximum Edge Length",
				"TotalEdgeCostVsMaximumEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMaximumEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMaximumEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMaximumEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeCostVsMeanDegreeOfConnectedNeighbors() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "Total Edge Cost",
				"Mean Degree Of Connected Neighbors",
				"Total Edge Cost vs Mean Degree Of Connected Neighbors",
				"TotalEdgeCostVsMeanDegreeOfConnectedNeighbors.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeCostVsMeanEdgeCost() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Total Edge Cost",
				"Mean Edge Cost",
				"Total Edge Cost vs Mean Edge Cost",
				"TotalEdgeCostVsMeanEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}


	private void plotTotalEdgeCostVsMeanEdgeLength() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Total Edge Cost",
				"Mean Edge Length",
				"Total Edge Cost vs Mean Edge Length",
				"TotalEdgeCostVsMeanEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeCostById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeLengthVsBetweeness() throws Exception 
	{
		MFilePlotData plotData = new MFilePlotData( "Total Edge Length",
				"Betweeness",
				"Total Edge Length vs Betweeness",
				"TotalEdgeLengthVsBetweeness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getBetweenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getBetweenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getBetweenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeLengthVsCloseness()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Total Edge Length",
				"Closeness",
				"Total Edge Length vs Closeness",
				"TotalEdgeLengthVsCloseness.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getClosenessById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getClosenessById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getClosenessById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeLengthVsConnectedNodesWithinRadiusR() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Total Edge Length","Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"Total Edge Length vs Connected Nodes Within Radius="+GraphConstants.MeasurementR,
				"TotalEdgeLengthVsConnectedNodesWithinRadiusR.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getConnectedNodesWithinRadiusRById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}



	private void plotTotalEdgeLengthVsMaximumEdgeCost()throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Total Edge Length","Maximum Edge Cost",
				"Total Edge Length vs Maximum Edge Cost",
				"TotalEdgeLengthVsMaximumEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMaximumEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMaximumEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMaximumEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeLengthVsMaximumEdgeLength()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Total Edge Length","Maximum Edge Length",
				"Total Edge Length vs Maximum Edge Length",
				"TotalEdgeLengthVsMaximumEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMaximumEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMaximumEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMaximumEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeLengthVsMeanDegreeOfConnectedNeighbors() throws Exception{
		
		
		MFilePlotData plotData = new MFilePlotData( "Total Edge Length",
				"Mean Degree Of Connected Neighbors",
				"Total Edge Length vs Mean Degree Of Connected Neighbors",
				"TotalEdgeLengthVsMeanDegreeOfConnectedNeighbors.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanDegreeOfConnectedNeighborsById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}



	private void plotTotalEdgeLengthVsMeanEdgeCost() throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Total Edge Length",
				"Mean Edge Cost",
				"Total Edge Length vs Mean Edge Cost",
				"TotalEdgeLengthVsMeanEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void plotTotalEdgeLengthVsMeanEdgeLength()throws Exception {
		MFilePlotData plotData = new MFilePlotData( "Total Edge Length",
				"Mean Edge Length",
				"Total Edge Length vs Mean Edge Length",
				"TotalEdgeLengthVsMeanEdgeLength.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getMeanEdgeLengthById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getMeanEdgeLengthById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getMeanEdgeLengthById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}



	private void plotTotalEdgeLengthVsTotalEdgeCost() throws Exception{
		MFilePlotData plotData = new MFilePlotData( "Total Edge Length","Total Edge Cost",
				"Total Edge Length vs Total Edge Cost",
				"TotalEdgeLengthVsTotalEdgeCost.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = optimized.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = optimized.get(i).getTotalEdgeCostById();
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				optimizedData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = exponential.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = exponential.get(i).getTotalEdgeCostById();
			for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
			{
				exponentialData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer, Double> nodeDistById = heuristic.get(i).getTotalEdgeLengthById();
			TreeMap<Integer, Double> edgeLenghtById = heuristic.get(i).getTotalEdgeCostById();
			for(int j=0; j < heuristic.get(i).getNumberOfNodes(); j++)
			{
				heuristicData.add(new PlotPairDoubleData(
						nodeDistById.get(j) ,
						edgeLenghtById.get(j) ));
			}	}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	public void generatePlotsH() throws Exception
	{
//		plotNodeDegreeVsTotalEdgeLength();
//		plotNodeDegreeVsTotalEdgeCost();
//		plotNodeDegreeVsMaximumEdgeLength();
//		plotNodeDegreeVsMaximumEdgeCost();
//		plotNodeDegreeVsConnectedNodesWithinRadiusR();
//		plotNodeDegreeVsMeanDegreeOfConnectedNeighbors();
		plotNodeDegreeVsBetweeness();
//		plotNodeDegreeVsCloseness();
//		plotNodeDegreeVsMeanEdgeLength();
//		plotNodeDegreeVsMeanEdgeCost();
		
	}

	public void generateConnectWithinRPlots() throws Exception
	{

		MFilePlotData plotData = new MFilePlotData( "Radius",
				"Connected Nodes Within Radius=", 
				"Radius vs Connected Nodes Within Radius R",
				"RadiusVsConnectedNodesWithinRadiusR.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedDataAverages = new Vector<PlotPairDoubleData>();
		Vector<PlotPairDoubleData> exponentialDataAverages = new Vector<PlotPairDoubleData>();
		
		
		
		Vector<Double> values = new Vector<Double>();
		values.add(.05);
		values.add(.06);
		values.add(.07);
		values.add(.08);
		values.add(.09);
		values.add(.10);
		values.add(.11);
		values.add(.12);
		values.add(.13);
		values.add(.14);
		values.add(.15);
		values.add(.16);
		values.add(.17);
		values.add(.18);
		values.add(.19);
		values.add(.20);
		values.add(.21);
		values.add(.22);
		values.add(.23);
		values.add(.24);
		values.add(.25);
		values.add(.26);
		values.add(.27);
		values.add(.28);
		values.add(.29);
		values.add(.30);
		values.add(.31);
		values.add(.32);
		values.add(.33);
		values.add(.34);
		values.add(.35);
		values.add(.36);
		values.add(.37);
		values.add(.38);
		values.add(.39);
		values.add(.40);
		values.add(.41);
		values.add(.42);
		values.add(.43);
		values.add(.44);
		values.add(.45);		
		values.add(.46);
		values.add(.47);
		values.add(.48);
		values.add(.49);
		values.add(.50);
		values.add(.51);
		values.add(.52);
		values.add(.53);
		values.add(.54);
		values.add(.55);			
		values.add(.56);
		values.add(.57);
		values.add(.58);
		values.add(.59);
		values.add(.60);
		values.add(.61);
		values.add(.62);
		values.add(.63);
		values.add(.64);
		values.add(.65);		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		double startValue = .05;
		double incValue = .01;
		double endValue = .25;
		double val=0;
		
		for(int k =0; k < values.size(); k++)
		{
			val = values.get(k);
			
			System.out.println("Value  is "+val);
			
			Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
			for(int i =0 ; i < optimized.size() ; i++)
			{
				optimized.get(i).UpdateWithinRInfo(val);
				
				TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getNodeDistributionId();
				TreeMap<Integer, Integer> connectedNodeWR = optimized.get(i).getConnectedNodesWithinRadiusRById();
				for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
				{
					optimizedData.add(new PlotPairDoubleData(
							nodeDistById.get(j) ,
							connectedNodeWR.get(j) ));
				}
			}
			//Now Calc avg value
			double avg =0;
			avg = CalcIntAverage(optimizedData);
			optimizedDataAverages.add(new PlotPairDoubleData(val, avg ));
			
			
			Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
			for(int i =0 ; i < exponential.size() ; i++)
			{
				exponential.get(i).UpdateWithinRInfo(val);
				
				TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getNodeDistributionId();
				TreeMap<Integer, Integer> connectedNodeWR = exponential.get(i).getConnectedNodesWithinRadiusRById();
				for(int j=0; j < exponential.get(i).getNumberOfNodes(); j++)
				{
					exponentialData.add(new PlotPairDoubleData(
							nodeDistById.get(j) ,
							connectedNodeWR.get(j) ));
				}
			}
			//Now Calc avg value
			avg =0;
			avg = CalcIntAverage(exponentialData);
			exponentialDataAverages.add(new PlotPairDoubleData(val, avg ));
	
		}
		
		
		
		


		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedDataAverages, exponentialDataAverages, new Vector<PlotPairDoubleData>());
		
	}

	private double CalcIntAverage(Vector<PlotPairDoubleData> data) 
	{
		double res =0;
		
		double sum = 0;
		
		for(int i =0; i < data.size(); i++)
		{
			sum += (double)data.get(i).m_y;
		}

		res = sum/data.size();
		
		return res;
	}

	public void generateHubbinessPlots() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "df",
				"", 
				"",
				"HubbinessCutoffValues.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedDataAverages = new Vector<PlotPairDoubleData>();
		Vector<PlotPairDoubleData> exponentialDataAverages = new Vector<PlotPairDoubleData>();
		

			
			Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
			for(int i =0 ; i < optimized.size() ; i++)
			{

				
				double cuttoff_limit = optimized.get(i).getHubinessCutoff();
				
					optimizedData.add(new PlotPairDoubleData(
							i ,
							cuttoff_limit));

			}
			
			
			Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
			for(int i =0 ; i < exponential.size() ; i++)
			{
				double cuttoff_limit = exponential.get(i).getHubinessCutoff();
				
				exponentialData.add(new PlotPairDoubleData(
						i ,
						cuttoff_limit));
				

			}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, new Vector<PlotPairDoubleData>());
		
		
	}

	public void plotLinkLengthGreaterThanRVsLocalNodeDegree() throws Exception
	{
		MFilePlotData plotData = new MFilePlotData( "Connected Nodes Greater than Radius1=.16 but less then Radius2=.5",
				"Nodes Greater than Radius1=.16 but less then Radius2=.5",
				"",
				"ConnectedNodesWRR12VsNodesWithinR12.m", outputDirectory, null, 20);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			
			
			TreeMap<Integer, Integer> nodeDistById = optimized.get(i).getConnectedNodesWithinRadiusR12ById();
			TreeMap<Integer, Integer> edgeLenghtById = optimized.get(i).getNodesWithinRadius12ById();
			
			TreeMap<Integer,Integer> specialDegreeById = optimized.get(i).getSpecialDegreeById();
			
			for(int j=0; j < optimized.get(i).getNumberOfNodes(); j++)
			{
				double percentage = (double)nodeDistById.get(j) / (double)edgeLenghtById.get(j);
				
				double specialDegree = specialDegreeById.get(j); 
				
				optimizedData.add(new PlotPairDoubleData(
						percentage ,
						specialDegree ));
			}
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < exponential.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = exponential.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = exponential.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < exponential.get(i).getN(); j++)
//			{
//				exponentialData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
//		for(int i =0 ; i < heuristic.size() ; i++)
//		{
//			TreeMap<Integer, Integer> nodeDistById = heuristic.get(i).getConnectedNodesWithinRadiusR12ById();
//			TreeMap<Integer, Integer> edgeLenghtById = heuristic.get(i).getNodesWithinRadius12ById();
//			for(int j=0; j < heuristic.get(i).getN(); j++)
//			{
//				heuristicData.add(new PlotPairDoubleData(
//						nodeDistById.get(j) ,
//						edgeLenghtById.get(j) ));
//			}	
//			}
		
		ValueWriter plotter = new ValueWriter(plotData);
		plotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}



}
