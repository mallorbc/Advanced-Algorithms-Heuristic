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

import psograph.graph.CalculatedGraph;
import psograph.graph.GraphConstants;
import psograph.util.HistogramWriter;
import psograph.util.MFilePlotData;



public class LayoutPropertiesWriter 
{
	File outputDirectory;
	CalculatedGraph layout;
	String layoutTitle;
	Vector<CalculatedGraph> seeds;
	
	
	public LayoutPropertiesWriter(CalculatedGraph Layout, File OutDir, String layoutTitle)
	{
		layout = Layout;
		outputDirectory = OutDir;
		this.layoutTitle = layoutTitle;
	}


	public LayoutPropertiesWriter(Vector<CalculatedGraph> seeds, File file) 
	{
		this.seeds = seeds;
		outputDirectory = file;
		this.layoutTitle = null;
	}


	public void generatePlots() throws Exception
	{
		
		//ClumpinessPlot
		MFilePlotData clumpinessData = new MFilePlotData("", "Clumpiness", "Clumpiness of Layout Seed",
				"clumpiness.m", outputDirectory,layoutTitle);
		TreeMap<Integer, Integer> clump = layout.getClumpinessOfLayout();
		HistogramWriter clumpinessPlotter = new HistogramWriter(clumpinessData);
		clumpinessPlotter.plotIntegerInteger(clump);
		
		
		//double fullFitness = layout.getFitnessOfSeedWithFullConnectivity();
		
		
		//betweenessInLayout
		System.out.println("Skipping BetweennessOfLayoutSeed");
//		MFilePlotData betweennessData = new MFilePlotData("", "BetweennessOfLayoutSeed", "Betweenness of Layout Seed",
//				"BetweennessOfLayoutSeed.m",outputDirectory,layoutTitle );
//		TreeMap<Double, Integer> betweenessInLayout = layout.getBetweenessInCompleteLayout();
//		HistogramWriter betweenessPlotter = new HistogramWriter(betweennessData);
//		betweenessPlotter.plotDoubleInteger(betweenessInLayout);
		
		
		
		
		//NodesWithinR
		MFilePlotData nodesWithinRData = new MFilePlotData("", "Nodes Within R", "Nodes Within radius R="+GraphConstants.MeasurementR,
		"nodesWithinR.m",outputDirectory,layoutTitle);
		TreeMap<Integer, Integer> nodesWithinRDist = layout.getNodesWithinRadius();
		HistogramWriter nodesWithinRPlotter = new HistogramWriter(nodesWithinRData);
		nodesWithinRPlotter.plotIntegerInteger(nodesWithinRDist);
		
		
		//NodesWithinR
		MFilePlotData CostOfDirectConnectivityData = 
			new MFilePlotData("","Cost Of Direct Connectivity", 
					"Cost Of Direct Connectivity",
					"CostOfDirectConnectivity.m",outputDirectory,layoutTitle);
		TreeMap<Double, Integer> costOfDirectConnectivityData = layout.getCostOfDirectConnectivityDist();
		HistogramWriter costOfDirectConnectivityDataPlotter = new HistogramWriter(CostOfDirectConnectivityData);
		costOfDirectConnectivityDataPlotter.plotDoubleInteger(costOfDirectConnectivityData);

		
	}
	
	
	public void generateCombinedPlots() throws Exception
	{
		
		//ClumpinessPlot
		MFilePlotData clumpinessData = new MFilePlotData("", "Clumpiness", "Clumpiness of Layout Seed",
				"clumpiness.m", outputDirectory,layoutTitle);
		TreeMap<Integer, Integer> clump = layout.getClumpinessOfLayout();
		HistogramWriter clumpinessPlotter = new HistogramWriter(clumpinessData);
		clumpinessPlotter.plotIntegerInteger(clump);
		
		
		//double fullFitness = layout.getFitnessOfSeedWithFullConnectivity();
		
		
		//betweenessInLayout
		System.out.println("Skipping BetweennessOfLayoutSeed");
//		MFilePlotData betweennessData = new MFilePlotData("", "BetweennessOfLayoutSeed", "Betweenness of Layout Seed",
//				"BetweennessOfLayoutSeed.m",outputDirectory,layoutTitle );
//		TreeMap<Double, Integer> betweenessInLayout = layout.getBetweenessInCompleteLayout();
//		HistogramWriter betweenessPlotter = new HistogramWriter(betweennessData);
//		betweenessPlotter.plotDoubleInteger(betweenessInLayout);
		
		
		
		
		//NodesWithinR
		MFilePlotData nodesWithinRData = new MFilePlotData("", "Nodes Within R", "Nodes Within radius R="+GraphConstants.MeasurementR,
		"nodesWithinR.m",outputDirectory,layoutTitle);
		TreeMap<Integer, Integer> nodesWithinRDist = layout.getNodesWithinRadius();
		HistogramWriter nodesWithinRPlotter = new HistogramWriter(nodesWithinRData);
		nodesWithinRPlotter.plotIntegerInteger(nodesWithinRDist);
		
		
		//NodesWithinR
		MFilePlotData CostOfDirectConnectivityData = 
			new MFilePlotData("","Cost Of Direct Connectivity", 
					"Cost Of Direct Connectivity",
					"CostOfDirectConnectivity.m",outputDirectory,layoutTitle);
		TreeMap<Double, Integer> costOfDirectConnectivityData = layout.getCostOfDirectConnectivityDist();
		HistogramWriter costOfDirectConnectivityDataPlotter = new HistogramWriter(CostOfDirectConnectivityData);
		costOfDirectConnectivityDataPlotter.plotDoubleInteger(costOfDirectConnectivityData);

		
	}

	
	
}
