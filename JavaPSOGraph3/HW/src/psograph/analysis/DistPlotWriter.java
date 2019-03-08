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
import psograph.util.DistWriter;
import psograph.util.MFilePlotData;



public class DistPlotWriter 
{
	File outputDirectory;
	
	Vector<CalculatedGraph> optimized;
	Vector<CalculatedGraph> exponential;
	Vector<CalculatedGraph> heuristic;
	
	Vector<CalculatedGraph> heur1;
	Vector<CalculatedGraph> heur2;
	Vector<CalculatedGraph> heur3;
	
	String layoutTitle;
	
	Vector<Vector<CalculatedGraph>> networks;
	
	
	public DistPlotWriter(Vector<CalculatedGraph> optimized, Vector<CalculatedGraph> exponential,
			Vector<CalculatedGraph> heuristic, File OutDir)
	{
		this.optimized = optimized;
		this.exponential = exponential;
		this.heuristic = heuristic;
		
		networks = new Vector<Vector<CalculatedGraph>>();
		networks.add(optimized);
		networks.add(exponential);
		networks.add(heuristic);

		
		
		outputDirectory = OutDir;
	}
	
	public DistPlotWriter(Vector<CalculatedGraph> optSolnsAll,
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
		
		networks = new Vector<Vector<CalculatedGraph>>();
		networks.add(optSolnsAll);
		networks.add(exponentialSolnsAll);
		networks.add(heur1All);
		networks.add(heur2All);
		networks.add(heur3All);
		networks.add(heur4All);
		
		outputDirectory = parentFile;
	}

	public void generatePlots() throws Exception
	{
	
		if(exponential.size() == 0 || heuristic.size()==0 || optimized.size()==0)
			return;
		
//		distributions
//		solnToExamine.getLinksLongerThan_r()
		//	System.out.println("generateLinksLongerThan_rPlot");
		generateLinksLongerThan_rPlot();
		
		
//		solnToExamine.getBetweeness()
		//	System.out.println("generateBetweeness");
		generateBetweeness();
		
//		solnToExamine.getCloseness()
		//	System.out.println("generateCloseness");
		generateCloseness();
		
//		solnToExamine.getConnectedNodesWithinRadius()
		//	System.out.println("generateConnectedNodesWithinRadius");
		generateConnectedNodesWithinRadius();
		
		
//		solnToExamine.getDistributionISPL()
		//	System.out.println("generateISPL");
		generateISPL();
		
//		solnToExamine.getDistributionLinkCosts()
		//	System.out.println("generateLinkCosts");
		generateLinkCosts();
		
//		solnToExamine.getDistributionLinkLengths()
		//	System.out.println("generateLinkLengths");
		generateLinkLengths();
		
//		solnToExamine.getDistributionSPL()
		//	System.out.println("generateSPL");
		generateSPL();
		
//		solnToExamine.getNodeDegreeStdDeviation()
		
		
//		solnToExamine.getNodeDegreeVariance()
		
//		solnToExamine.getNodeDistribution()
		//	System.out.println("generateNodeDistribution");
		generateNodeDistribution();
		
//		solnToExamine.getNodeMaximumEdgeCost()
		//	System.out.println("generateNodeMaximumEdgeCost");
		generateNodeMaximumEdgeCost();
		
//		solnToExamine.getNodeMaximumEdgeLength()
		//	System.out.println("generateNodeMaxiumumEdgeLength");
		generateNodeMaxiumumEdgeLength();
		
//		solnToExamine.getNodeTotalEdgeCost()
		//	System.out.println("generateNodeTotalEdgeCost");
		generateNodeTotalEdgeCost();
		
//		solnToExamine.getNodeTotalEdgeLength()
		//	System.out.println("generateNodeTotalEdgeLength");
		generateNodeTotalEdgeLength();
		
		//	System.out.println("generateMeanDegreeOfConnectedNeighbors");
		generateMeanDegreeOfConnectedNeighbors();
		
	
		generateFitness();

		
	}

	private void generateMeanDegreeOfConnectedNeighbors() throws Exception {
		MFilePlotData NodeTotalEdgeCostData = new MFilePlotData( "","Mean Degree Of Connected Neighbors", "Mean Degree Of Connected Neighbors",
				"MeanDegreeOfConnectedNeighbors.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Double,Integer>> optimizedData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Double,Integer> dist = optimized.get(i).getMeanDegreeOfConnectedNeighbors();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> exponentialData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Double,Integer> dist = exponential.get(i).getMeanDegreeOfConnectedNeighbors();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> heuristicData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Double,Integer> dist = heuristic.get(i).getMeanDegreeOfConnectedNeighbors();
			heuristicData.add(dist);
			
		}
		
		DistWriter NodeTotalEdgeLengthPlotter = new DistWriter(NodeTotalEdgeCostData);
		NodeTotalEdgeLengthPlotter.PlotDist(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateFitness() throws Exception
	{
		MFilePlotData mFilePlotData = new MFilePlotData( "Fitness","", "Distribution of Fitness",
				"FitnessDist.m", outputDirectory, null, 20);
		
		Vector<Double> optimizedData = new Vector<Double>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(optimized.get(i).getFitnessValue());
		}
		
		Vector<Double> exponentialData = new Vector<Double>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(exponential.get(i).getFitnessValue());
		}
		
		Vector<Double> heuristicData = new Vector<Double>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(heuristic.get(i).getFitnessValue());
			
		}
		
		DistWriter writer = new DistWriter(mFilePlotData);
		writer.plotDoubleVectors(optimizedData, exponentialData, heuristicData);
		
	}
	
	private void generateNodeTotalEdgeLength() throws Exception
	{
		MFilePlotData NodeTotalEdgeCostData = new MFilePlotData( "","Node Total Edge Length", "Node Total Edge Length",
				"NodeTotalEdgeLength.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Double,Integer>> optimizedData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Double,Integer> dist = optimized.get(i).getNodeTotalEdgeLength();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> exponentialData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Double,Integer> dist = exponential.get(i).getNodeTotalEdgeLength();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> heuristicData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Double,Integer> dist = heuristic.get(i).getNodeTotalEdgeLength();
			heuristicData.add(dist);
			
		}
		
		DistWriter NodeTotalEdgeLengthPlotter = new DistWriter(NodeTotalEdgeCostData);
		NodeTotalEdgeLengthPlotter.PlotDist(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateNodeTotalEdgeCost() throws Exception{
		MFilePlotData getNodeTotalEdgeCostData = new MFilePlotData( "","Node Total Edge Cost", "Node Total Edge Cost",
				"NodeTotalEdgeCost.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Double,Integer>> optimizedData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Double,Integer> dist = optimized.get(i).getNodeTotalEdgeCost();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> exponentialData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Double,Integer> dist = exponential.get(i).getNodeTotalEdgeCost();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> heuristicData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Double,Integer> dist = heuristic.get(i).getNodeTotalEdgeCost();
			heuristicData.add(dist);
			
		}
		
		DistWriter NodeTotalEdgeCostPlotter = new DistWriter(getNodeTotalEdgeCostData);
		NodeTotalEdgeCostPlotter.PlotDist(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateNodeMaxiumumEdgeLength()throws Exception {
		MFilePlotData NodeMaximumEdgeData = new MFilePlotData( "","Node Maximum Edge Length", "Node Maximum Edge Length",
				"NodeMaximumEdgeLength.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Double,Integer>> optimizedData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Double,Integer> dist = optimized.get(i).getNodeMaximumEdgeLength();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> exponentialData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Double,Integer> dist = exponential.get(i).getNodeMaximumEdgeLength();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> heuristicData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Double,Integer> dist = heuristic.get(i).getNodeMaximumEdgeLength();
			heuristicData.add(dist);
			
		}
		
		DistWriter NodeMaximumEdgeLengthPlotter = new DistWriter(NodeMaximumEdgeData);
		NodeMaximumEdgeLengthPlotter.PlotDist(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateNodeMaximumEdgeCost()throws Exception 
	{
		MFilePlotData NodeMaximumEdgeData = new MFilePlotData( "","Node Maximum Edge Cost", "Node Maximum Edge Cost",
				"NodeMaximumEdgeCost.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Double,Integer>> optimizedData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Double,Integer> dist = optimized.get(i).getNodeMaximumEdgeCost();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> exponentialData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Double,Integer> dist = exponential.get(i).getNodeMaximumEdgeCost();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> heuristicData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Double,Integer> dist = heuristic.get(i).getNodeMaximumEdgeCost();
			heuristicData.add(dist);
			
		}
		
		DistWriter NodeMaximumEdgePlotter = new DistWriter(NodeMaximumEdgeData);
		NodeMaximumEdgePlotter.PlotDist(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateNodeDistribution() throws Exception
	{
		MFilePlotData NodeDistributionData = new MFilePlotData( "","Node Degree Distribution", "Node Degree Distribution",
				"NodeDistribution.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Integer,Integer>> optimizedData = new Vector<TreeMap<Integer,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer,Integer> dist = optimized.get(i).getNodeDistribution();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Integer,Integer>> exponentialData = new Vector<TreeMap<Integer,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer,Integer> dist = exponential.get(i).getNodeDistribution();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Integer,Integer>> heuristicData = new Vector<TreeMap<Integer,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer,Integer> dist = heuristic.get(i).getNodeDistribution();
			heuristicData.add(dist);
			
		}
		
		DistWriter NodeDistributionPlotter = new DistWriter(NodeDistributionData);
		NodeDistributionPlotter.PlotDistInt(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateSPL() throws Exception
	{
		MFilePlotData SPLData = new MFilePlotData( "","SPL", "SPL",
				"SPL.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Integer,Integer>> optimizedData = new Vector<TreeMap<Integer,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer,Integer> dist = optimized.get(i).getDistributionSPL();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Integer,Integer>> exponentialData = new Vector<TreeMap<Integer,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer,Integer> dist = exponential.get(i).getDistributionSPL();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Integer,Integer>> heuristicData = new Vector<TreeMap<Integer,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer,Integer> dist = heuristic.get(i).getDistributionSPL();
			heuristicData.add(dist);
			
		}
		
		DistWriter SPLPlotter = new DistWriter(SPLData);
		SPLPlotter.PlotDistInt(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateLinkLengths()throws Exception {
		
		MFilePlotData LinkLengthsData = new MFilePlotData( "","Link Lengths", "Link Lengths",
				"LinkLengths.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Double,Integer>> optimizedData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Double,Integer> dist = optimized.get(i).getDistributionLinkLengths();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> exponentialData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Double,Integer> dist = exponential.get(i).getDistributionLinkLengths();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> heuristicData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Double,Integer> dist = heuristic.get(i).getDistributionLinkLengths();
			heuristicData.add(dist);
			
		}
		
		DistWriter LinkLengthsPlotter = new DistWriter(LinkLengthsData);
		LinkLengthsPlotter.PlotDist(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateLinkCosts()throws Exception {
		
		
		MFilePlotData LinkCostsData = new MFilePlotData( "","Link Costs", "Link Costs",
				"LinkCosts.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Double,Integer>> optimizedData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Double,Integer> dist = optimized.get(i).getDistributionLinkCosts();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> exponentialData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Double,Integer> dist = exponential.get(i).getDistributionLinkCosts();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> heuristicData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Double,Integer> dist = heuristic.get(i).getDistributionLinkCosts();
			heuristicData.add(dist);
			
		}
		
		DistWriter LinkCostsPlotter = new DistWriter(LinkCostsData);
		LinkCostsPlotter.PlotDist(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateISPL()throws Exception 
	{
		MFilePlotData ISPLData = new MFilePlotData( "","ISPL", "ISPL",
				"ISPL.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Double,Integer>> optimizedData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Double,Integer> dist = optimized.get(i).getDistributionISPL();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> exponentialData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Double,Integer> dist = exponential.get(i).getDistributionISPL();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> heuristicData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Double,Integer> dist = heuristic.get(i).getDistributionISPL();
			heuristicData.add(dist);
			
		}
		
		DistWriter ISPLPlotter = new DistWriter(ISPLData);
		ISPLPlotter.PlotDist(optimizedData, exponentialData, heuristicData);
		
	}

	

	private void generateConnectedNodesWithinRadius()throws Exception 
	{
		MFilePlotData ConnectedNodesWithinRadiusData = new MFilePlotData( "","Connected Nodes Within R="+GraphConstants.MeasurementR, "Connected Nodes Within Radius"+GraphConstants.MeasurementR,
				"ConnectedNodesWithinRadius_questionable.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Integer,Integer>> optimizedData = new Vector<TreeMap<Integer,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Integer,Integer> dist = optimized.get(i).getConnectedNodesWithinRadius();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Integer,Integer>> exponentialData = new Vector<TreeMap<Integer,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Integer,Integer> dist = exponential.get(i).getConnectedNodesWithinRadius();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Integer,Integer>> heuristicData = new Vector<TreeMap<Integer,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Integer,Integer> dist = heuristic.get(i).getConnectedNodesWithinRadius();
			heuristicData.add(dist);
			
		}
		
		DistWriter ConnectedNodesWithinRadiusPlotter = new DistWriter(ConnectedNodesWithinRadiusData);
		ConnectedNodesWithinRadiusPlotter.PlotDistInt(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateCloseness() throws Exception
	{
		MFilePlotData ClosenessData = new MFilePlotData( "","Closeness", "Closeness",
				"Closeness.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Double,Integer>> optimizedData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Double,Integer> dist = optimized.get(i).getCloseness();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> exponentialData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Double,Integer> dist = exponential.get(i).getCloseness();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> heuristicData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Double,Integer> dist = heuristic.get(i).getCloseness();
			heuristicData.add(dist);
			
		}
		
		DistWriter ClosenessPlotter = new DistWriter(ClosenessData);
		ClosenessPlotter.PlotDist(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateBetweeness() throws Exception
	{
		MFilePlotData BetweenessData = new MFilePlotData( "","Betweeness", "Betweeness",
				"Betweeness.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Double,Integer>> optimizedData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Double,Integer> dist = optimized.get(i).getBetweeness();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> exponentialData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			System.out.println("exp iter"+i);
			TreeMap<Double,Integer> dist = exponential.get(i).getBetweeness();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> heuristicData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Double,Integer> dist = heuristic.get(i).getBetweeness();
			heuristicData.add(dist);
			
		}
		
		DistWriter BetweenessPlotter = new DistWriter(BetweenessData);
		BetweenessPlotter.PlotDist(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateLinksLongerThan_rPlot() throws Exception
	{
		MFilePlotData getLinksLongerThanRData = new MFilePlotData( "","Links longer than R="+GraphConstants.MeasurementR, "Links longer than R="+GraphConstants.MeasurementR,
				"LinksLongerThanR.m", outputDirectory, null, 20);
		
		Vector<TreeMap<Double,Integer>> optimizedData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			TreeMap<Double,Integer> dist = optimized.get(i).getLinksLongerThan_r();
			optimizedData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> exponentialData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			TreeMap<Double,Integer> dist = exponential.get(i).getLinksLongerThan_r();
			exponentialData.add(dist);
		}
		
		Vector<TreeMap<Double,Integer>> heuristicData = new Vector<TreeMap<Double,Integer>>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			TreeMap<Double,Integer> dist = heuristic.get(i).getLinksLongerThan_r();
			heuristicData.add(dist);
			
		}
		
		DistWriter LinksLongerThanRPlotter = new DistWriter(getLinksLongerThanRData);
		LinksLongerThanRPlotter.PlotDist(optimizedData, exponentialData, heuristicData);
		
	}

	public void generatePlotsH() throws Exception
	{
		
		//generateNodeDistributionH();
		//generateConnectedNodesWithinRadiusH();
		generateLinkLengthsH();
		
		
	}

	private void generateLinkLengthsH() throws Exception {
		
		MFilePlotData NodeDistributionData = new MFilePlotData( "","Node Degree Distribution", "Node Degree Distribution",
				"DistLinkLengths.m", outputDirectory, null, 20);
		
		Vector<Vector<TreeMap<Double,Integer>>> data = new Vector<Vector<TreeMap<Double,Integer>>>();
		
		for(int i = 0; i < 1;
		//networks.size();
		i++)
		{
			Vector<TreeMap<Double,Integer>> curNetworkData = new Vector<TreeMap<Double,Integer>>();
			
			Vector<CalculatedGraph> networkSeries = networks.get(i);
			
			for(int j =0 ; j < networkSeries.size() ; j++)
			{
				TreeMap<Double,Integer> dist = networkSeries.get(j).getDistributionLinkLengths();
				curNetworkData.add(dist);
			}
			
			data.add(curNetworkData);
		}
		
		for(int i = 3; i < 4; i++)
		{
			Vector<TreeMap<Double,Integer>> curNetworkData = new Vector<TreeMap<Double,Integer>>();
			
			Vector<CalculatedGraph> networkSeries = networks.get(i);
			
			for(int j =0 ; j < networkSeries.size() ; j++)
			{
				TreeMap<Double,Integer> dist = networkSeries.get(j).getDistributionLinkLengths();
				curNetworkData.add(dist);
			}
			
			data.add(curNetworkData);
		}

		
		DistWriter NodeDistributionPlotter = new DistWriter(NodeDistributionData);
		NodeDistributionPlotter.PlotDist(data);
		
		

		
	}

	private void generateConnectedNodesWithinRadiusH()  throws Exception{
		MFilePlotData NodeDistributionData = new MFilePlotData( "","Node Degree Distribution", "Node Degree Distribution",
				"DistNodesWithinR.m", outputDirectory, null, 20);
		
		Vector<Vector<TreeMap<Integer,Integer>>> data = new Vector<Vector<TreeMap<Integer,Integer>>>();
		
		for(int i = 0; i < networks.size(); i++)
		{
			Vector<TreeMap<Integer,Integer>> curNetworkData = new Vector<TreeMap<Integer,Integer>>();
			
			Vector<CalculatedGraph> networkSeries = networks.get(i);
			
			for(int j =0 ; j < networkSeries.size() ; j++)
			{
				TreeMap<Integer,Integer> dist = networkSeries.get(j).getConnectedNodesWithinRadius(.16);
				curNetworkData.add(dist);
			}
			
			data.add(curNetworkData);
		}

		
		DistWriter NodeDistributionPlotter = new DistWriter(NodeDistributionData);
		NodeDistributionPlotter.PlotDistInt(data);
		
	}

	private void generateNodeDistributionH() throws Exception
	{
		MFilePlotData NodeDistributionData = new MFilePlotData( "","Node Degree Distribution", "Node Degree Distribution",
				"NodeDistribution.m", outputDirectory, null, 20);
		
		Vector<Vector<TreeMap<Integer,Integer>>> data = new Vector<Vector<TreeMap<Integer,Integer>>>();
		
		for(int i = 0; i < networks.size();i++)
		{
			Vector<TreeMap<Integer,Integer>> curNetworkData = new Vector<TreeMap<Integer,Integer>>();
			
			Vector<CalculatedGraph> networkSeries = networks.get(i);
			
			System.out.println("How many graphs - "+networkSeries.size());
			
			for(int j =0 ; j < networkSeries.size() ; j++)
			{
				TreeMap<Integer,Integer> dist = networkSeries.get(j).getNodeDistribution();
				
				curNetworkData.add(dist);
			}
			
			data.add(curNetworkData);
		}

		
		DistWriter NodeDistributionPlotter = new DistWriter(NodeDistributionData);
		NodeDistributionPlotter.PlotDistInt(data);
		
	}


	
}
