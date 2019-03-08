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
import java.util.Vector;

import psograph.graph.CalculatedGraph;
import psograph.util.MFilePlotData;
import psograph.util.PlotPairDoubleData;
import psograph.util.ValueWriter;

public class ValuePlotWriter 
{
	File outputDirectory;
	
	Vector<CalculatedGraph> optimized;
	Vector<CalculatedGraph> exponential;
	Vector<CalculatedGraph> heuristic;
	

	Vector<CalculatedGraph> heur1;
	Vector<CalculatedGraph> heur2;
	Vector<CalculatedGraph> heur3;
	
	Vector<Vector<CalculatedGraph>> networks;
	
	String layoutTitle;
	
	
	public ValuePlotWriter(Vector<CalculatedGraph> optimized, Vector<CalculatedGraph> exponential,
			Vector<CalculatedGraph> heuristic, File OutDir)
	{
		this.optimized = optimized;
		this.exponential = exponential;
		this.heuristic = heuristic;
		outputDirectory = OutDir;
	}

	public ValuePlotWriter(Vector<CalculatedGraph> optSolnsAll,
			Vector<CalculatedGraph> heur1All, Vector<CalculatedGraph> heur2All,
			Vector<CalculatedGraph> heur3All, File parentFile) 
	{
		this.optimized =optSolnsAll;
		this.heur1 =heur1All;
		this.heur2 =heur2All;
		this.heur3 =heur3All;
		outputDirectory = parentFile;
	}

	public ValuePlotWriter(Vector<CalculatedGraph> optSolnsAll,
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
		
		networks = new Vector<Vector<CalculatedGraph>>();
		networks.add(optSolnsAll);
		networks.add(exponentialSolnsAll);
		networks.add(heur1All);
		networks.add(heur2All);
		networks.add(heur3All);
		networks.add(heur4All);

		
	}

	private void generateAISPLPlot() throws Exception
	{
	
		
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "AISPL","Fitness", "AISPL vs Fitness",
				"AISPLVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAISPL() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAISPL() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAISPL() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	
	public void generateAISPLPlotVsAssortativity()throws Exception {
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "AISPL","Assortativity", "AISPL vs Assortativity",
				"AISPLvsAssortativity.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAISPL() , optimized.get(i).getAssortativity() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAISPL() , exponential.get(i).getAssortativity() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAISPL() , heuristic.get(i).getAssortativity() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	
	private void generateAISPPLlotVsAvgLinkCost()throws Exception {
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "AISPL","Avgerage Link Cost", "AISPL vs Avgerage Link Cost",
				"AISPLvsAvgLinkCost.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAISPL() , optimized.get(i).getAvgLinkCost() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAISPL() , exponential.get(i).getAvgLinkCost() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAISPL() , heuristic.get(i).getAvgLinkCost() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	
	private void generateAISPPLlotVsAvgLinkLength()throws Exception {
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "AISPL","Avgerage Link Length", "AISPL vs Avgerage Link Length",
				"AISPLvsAvgLinkLength.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAISPL() , optimized.get(i).getAvgLinkLength() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAISPL() , exponential.get(i).getAvgLinkLength() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAISPL() , heuristic.get(i).getAvgLinkLength() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}



	private void generateAISPPLlotVsAvgNodeDegree() throws Exception{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "AISPL","Avgerage Node Degree", "AISPL vs Avgerage Node Degree",
				"AISPLvsAvgNodeDegree.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAISPL() , optimized.get(i).getAvgNodeDegree() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAISPL() , exponential.get(i).getAvgNodeDegree() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAISPL() , heuristic.get(i).getAvgNodeDegree() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAISPPLlotVsClusteringCoefficient()throws Exception {
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "AISPL","Cluster Coefficient", "AISPL vs Cluster Coefficient",
				"AISPLvsClusterCoeff.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAISPL() , optimized.get(i).getClusteringCoefficient() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAISPL() , exponential.get(i).getClusteringCoefficient() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAISPL() , heuristic.get(i).getClusteringCoefficient() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}


	private void generateAISPPLlotVsDiameter() throws Exception {
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "AISPL","Diameter", "AISPL vs Diamter",
				"AISPLvsDiameter.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAISPL() , optimized.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAISPL() , exponential.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAISPL() , heuristic.get(i).getDiameter() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAISPPLlotVsHubbiness() throws Exception
	{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "AISPL","Hubbiness", "AISPL vs Hubbiness",
				"AISPLvsHubbiness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAISPL() , optimized.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAISPL() , exponential.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAISPL() , heuristic.get(i).getHubbiness() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAISPPlotVsPercentNodesConnectedToHubPlot() throws Exception {
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "AISPL","Percent Nodes Connected To Hubs", "AISPL vs Percent Nodes Connected To Hubs",
				"AISPLvsPercentNodesConnectedToHub.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAISPL() , optimized.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAISPL() , exponential.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAISPL() , heuristic.get(i).getPercentNodesConnectedToHub() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateASPLPlot() throws Exception
	{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "ASPL","Fitness", "ASPL vs Fitness",
				"ASPLVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getASPL() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getASPL() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getASPL() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateASPLPlotVsAssortativity()throws Exception{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "ASPL","Assortativity ", "ASPL vs Assortativity",
				"ASPLVsAssortativity.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getASPL() , optimized.get(i).getAssortativity() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getASPL() , exponential.get(i).getAssortativity() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getASPL() , heuristic.get(i).getAssortativity() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateASPPLlotVsAvgLinkCost() throws Exception{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "ASPL"," Avgerage Link Cost ", "ASPL vs Avgerage Link Cost",
				"ASPLVsAvgLinkCost.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getASPL() , optimized.get(i).getAvgLinkCost() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getASPL() , exponential.get(i).getAvgLinkCost() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getASPL() , heuristic.get(i).getAvgLinkCost() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateASPPLlotVsAvgLinkLength()throws Exception {
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "ASPL"," Avgerage Link Length ", "ASPL vs Avgerage Link Length",
				"ASPLVsAvgLinkLength.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getASPL() , optimized.get(i).getAvgLinkLength() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getASPL() , exponential.get(i).getAvgLinkLength() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getASPL() , heuristic.get(i).getAvgLinkLength() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}




	private void generateASPPLlotVsAvgNodeDegree() throws Exception{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "ASPL"," Avgerage Node Degree ", "ASPL vs vgerage Node Degree",
				"ASPLVsAvgNodeDegree.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getASPL() , optimized.get(i).getAvgNodeDegree() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getASPL() , exponential.get(i).getAvgNodeDegree() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getASPL() , heuristic.get(i).getAvgNodeDegree() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateASPPLlotVsClusteringCoefficient()throws Exception {
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "ASPL"," Clustering Coefficient ", "ASPL vs Clustering Coefficient",
				"ASPLVsClusteringCoefficient.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getASPL() , optimized.get(i).getClusteringCoefficient() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getASPL() , exponential.get(i).getClusteringCoefficient() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getASPL() , heuristic.get(i).getClusteringCoefficient() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateASPPLlotVsDiameter()throws Exception {
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "ASPL"," Diameter ", "ASPL vs Diameter",
				"ASPLVsDiameter.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getASPL() , optimized.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getASPL() , exponential.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getASPL() , heuristic.get(i).getDiameter() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateASPPLlotVsHubbiness() throws Exception{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "ASPL"," Hubbiness ", "ASPL vs Hubbiness",
				"ASPLVsHubbiness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getASPL() , optimized.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getASPL() , exponential.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getASPL() , heuristic.get(i).getHubbiness() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateASPPlotVsPercentNodesConnectedToHubPlot() throws Exception{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "ASPL"," Percent Nodes Connected To Hub ", "ASPL vs Percent Nodes Connected To Hub",
				"ASPLVsPercentNodesConnectedToHub.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getASPL() , optimized.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getASPL() , exponential.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getASPL() , heuristic.get(i).getPercentNodesConnectedToHub() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAssortativityPlot() throws Exception
	{

		MFilePlotData AISPLData = new MFilePlotData( "Assortativity", "Fitness", "Assortativity vs. Fitness",
				"AssortativityVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAssortativity() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAssortativity() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAssortativity() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAssortativityVsAvgLinkCost() throws Exception
	{
		MFilePlotData AISPLData = new MFilePlotData( "Assortativity", "Average Link Cost", "Assortativity vs Average Link Cost",
				"AssortativityVsAvgLinkCost.m", outputDirectory, null);
		 
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAssortativity() , optimized.get(i).getAvgLinkCost()));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAssortativity() , exponential.get(i).getAvgLinkCost() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAssortativity() , heuristic.get(i).getAvgLinkCost() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAssortativityVsAvgLinkLength() throws Exception{
		MFilePlotData AISPLData = new MFilePlotData( "Assortativity", "Average Link Length", "Assortativity vs Average Link Length",
				"AssortativityVsAvgLinkLength.m", outputDirectory, null);
		 
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAssortativity() , optimized.get(i).getAvgLinkLength()));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAssortativity() , exponential.get(i).getAvgLinkLength() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAssortativity() , heuristic.get(i).getAvgLinkLength() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAssortativityVsAvgNodeDegree()throws Exception {
		MFilePlotData AISPLData = new MFilePlotData( "Assortativity", "Average Node Degree", "Assortativity vs Average Node Degree",
				"AssortativityVsAvgNodeDegree.m", outputDirectory, null);
		 
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAssortativity() , optimized.get(i).getAvgNodeDegree()));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAssortativity() , exponential.get(i).getAvgNodeDegree() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAssortativity() , heuristic.get(i).getAvgNodeDegree() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAssortativityVsClusteringCoefficient() throws Exception{
		MFilePlotData AISPLData = new MFilePlotData( "Assortativity", "ClusteringCoefficient", "Assortativity vs. ClusteringCoefficient",
				"AssortativityVsClusteringCoefficient.m", outputDirectory, null);
		 
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAssortativity() , optimized.get(i).getClusteringCoefficient()));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAssortativity() , exponential.get(i).getClusteringCoefficient() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAssortativity() , heuristic.get(i).getClusteringCoefficient() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAssortativityVsDiameter()throws Exception {
		MFilePlotData AISPLData = new MFilePlotData( "Assortativity", "Diameter", "Assortativity vs. Diameter",
				"AssortativityVsDiameter.m", outputDirectory, null);
		 
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAssortativity() , optimized.get(i).getDiameter()));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAssortativity() , exponential.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAssortativity() , heuristic.get(i).getDiameter() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAssortativityVsHubbiness() throws Exception{
		MFilePlotData AISPLData = new MFilePlotData( "Assortativity", "Hubbiness", "Assortativity vs. Hubbiness",
				"AssortativityVsHubbiness.m", outputDirectory, null);
		 
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAssortativity() , optimized.get(i).getHubbiness()));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAssortativity() , exponential.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAssortativity() , heuristic.get(i).getHubbiness() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAssortativityVsPercentNodesConnectedToHubPlot() throws Exception{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "Assortativity"," Percent Nodes Connected To Hub ", "Assortativity vs Percent Nodes Connected To Hub",
				"AssortativityVsPercentNodesConnectedToHub.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAssortativity() , optimized.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAssortativity() , exponential.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAssortativity() , heuristic.get(i).getPercentNodesConnectedToHub() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgLinkCostPlot() throws Exception {
		MFilePlotData getAvgLinkCostData = new MFilePlotData( "Average Link Cost","Fitness", "Average Link Cost vs Fitness",
				"AvgLinkCostVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgLinkCost() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgLinkCost() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgLinkCost() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter AvgLinkCostPlotter = new ValueWriter(getAvgLinkCostData);
		AvgLinkCostPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgLinkCostPlotVsAvgNodeDegree()throws Exception {
		MFilePlotData getAvgLinkCostData = new MFilePlotData( "Average Link Cost","Averge Node Degree", "Average Link Cost vs Avgerage Node Degree",
				"AvgLinkCostVsAvgNodeDegree.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgLinkCost() , optimized.get(i).getAvgNodeDegree() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgLinkCost() , exponential.get(i).getAvgNodeDegree() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgLinkCost() , heuristic.get(i).getAvgNodeDegree() ));
		}
		
		ValueWriter AvgLinkCostPlotter = new ValueWriter(getAvgLinkCostData);
		AvgLinkCostPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgLinkCostPlotVsClusteringCoefficient() throws Exception{
		MFilePlotData getAvgLinkCostData = new MFilePlotData( "Average Link Cost","ClusteringCoefficient", "Average Link Cost vs ClusteringCoefficient",
				"AvgLinkCostVsClusteringCoefficient.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgLinkCost() , optimized.get(i).getClusteringCoefficient() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgLinkCost() , exponential.get(i).getClusteringCoefficient() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgLinkCost() , heuristic.get(i).getClusteringCoefficient() ));
		}
		
		ValueWriter AvgLinkCostPlotter = new ValueWriter(getAvgLinkCostData);
		AvgLinkCostPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgLinkCostPlotVsDiameter() throws Exception{
		MFilePlotData getAvgLinkCostData = new MFilePlotData( "Average Link Cost","Diameter", "Average Link Cost vs Diameter",
				"AvgLinkCostVsDiameter.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgLinkCost() , optimized.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgLinkCost() , exponential.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgLinkCost() , heuristic.get(i).getDiameter() ));
		}
		
		ValueWriter AvgLinkCostPlotter = new ValueWriter(getAvgLinkCostData);
		AvgLinkCostPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgLinkCostPlotVsHubbiness()throws Exception {
		MFilePlotData getAvgLinkCostData = new MFilePlotData( "Average Link Cost","Hubbiness", "Average Link Cost vs Hubbiness",
				"AvgLinkCostVsHubbiness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgLinkCost() , optimized.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgLinkCost() , exponential.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgLinkCost() , heuristic.get(i).getHubbiness() ));
		}
		
		ValueWriter AvgLinkCostPlotter = new ValueWriter(getAvgLinkCostData);
		AvgLinkCostPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgLinkCostVsPercentNodesConnectedToHubPlot() throws Exception{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "Avgerage Link Cost"," Percent Nodes Connected To Hub ", "Avgerage Link Cost vs Percent Nodes Connected To Hub",
				"AvgLinkCostVsPercentNodesConnectedToHub.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgLinkCost() , optimized.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgLinkCost() , exponential.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgLinkCost() , heuristic.get(i).getPercentNodesConnectedToHub() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgLinkLengthPlot() throws Exception 
	{
		MFilePlotData AvgLinkLengthData = new MFilePlotData( "Average Link Length","Fitness", "Average Link Length vs Fitness",
				"AvgLinkLengthVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgLinkLength() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgLinkLength() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgLinkLength() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter AvgLinkLengthPlotter = new ValueWriter(AvgLinkLengthData);
		AvgLinkLengthPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgLinkLengthPlotVsAvgNodeDegree()throws Exception {
		MFilePlotData AvgLinkLengthData = new MFilePlotData( "Average Link Length","Average Node Degree", "Average Link Length vs Average Node Degree",
				"AvgLinkLengthVsAvgNodeDegree.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgLinkLength() , optimized.get(i).getAvgNodeDegree() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgLinkLength() , exponential.get(i).getAvgNodeDegree() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgLinkLength() , heuristic.get(i).getAvgNodeDegree() ));
		}
		
		ValueWriter AvgLinkLengthPlotter = new ValueWriter(AvgLinkLengthData);
		AvgLinkLengthPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgLinkLengthPlotVsClusteringCoefficient() throws Exception{
		MFilePlotData AvgLinkLengthData = new MFilePlotData( "Average Link Length","ClusteringCoefficient", "Average Link Length vs ClusteringCoefficient",
				"AvgLinkLengthVsClusteringCoefficient.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgLinkLength() , optimized.get(i).getClusteringCoefficient() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgLinkLength() , exponential.get(i).getClusteringCoefficient() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgLinkLength() , heuristic.get(i).getClusteringCoefficient() ));
		}
		
		ValueWriter AvgLinkLengthPlotter = new ValueWriter(AvgLinkLengthData);
		AvgLinkLengthPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgLinkLengthPlotVsDiameter()throws Exception {
		MFilePlotData AvgLinkLengthData = new MFilePlotData( "Average Link Length","Diameter", "Average Link Length vs Diameter",
				"AvgLinkLengthVsDiameter.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgLinkLength() , optimized.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgLinkLength() , exponential.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgLinkLength() , heuristic.get(i).getDiameter() ));
		}
		
		ValueWriter AvgLinkLengthPlotter = new ValueWriter(AvgLinkLengthData);
		AvgLinkLengthPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgLinkLengthPlotVsHubbiness() throws Exception{
		MFilePlotData AvgLinkLengthData = new MFilePlotData( "Average Link Length","Hubbiness", "Average Link Length vs Hubbiness",
				"AvgLinkLengthVsHubbiness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgLinkLength() , optimized.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgLinkLength() , exponential.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgLinkLength() , heuristic.get(i).getHubbiness() ));
		}
		
		ValueWriter AvgLinkLengthPlotter = new ValueWriter(AvgLinkLengthData);
		AvgLinkLengthPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgLinkLengthVsPercentNodesConnectedToHubPlot() throws Exception{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "Avgerage Link Length"," Percent Nodes Connected To Hub ", "Avgerage Link Length vs Percent Nodes Connected To Hub",
				"AvgLinkLengthVsPercentNodesConnectedToHub.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgLinkLength() , optimized.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgLinkLength() , exponential.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgLinkLength() , heuristic.get(i).getPercentNodesConnectedToHub() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgNodeDegree() throws Exception 
	{
		MFilePlotData AvgNodeDegreeData = new MFilePlotData( "Average Node Degree","Fitness", "Average Node Degree vs Fitness",
				"AvgNodeDegreeVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgNodeDegree() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgNodeDegree() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgNodeDegree() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter AvgNodeDegreePlotter = new ValueWriter(AvgNodeDegreeData);
		AvgNodeDegreePlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgNodeDegreeVsClusteringCoefficient() throws Exception{
		MFilePlotData AvgNodeDegreeData = new MFilePlotData( "Average Node Degree","ClusteringCoefficient", "Average Node Degree vs ClusteringCoefficient",
				"AvgNodeDegreeVsClusteringCoefficient.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgNodeDegree() , optimized.get(i).getClusteringCoefficient() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgNodeDegree() , exponential.get(i).getClusteringCoefficient() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgNodeDegree() , heuristic.get(i).getClusteringCoefficient() ));
		}
		
		ValueWriter AvgNodeDegreePlotter = new ValueWriter(AvgNodeDegreeData);
		AvgNodeDegreePlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgNodeDegreeVsDiameter()throws Exception {
		MFilePlotData AvgNodeDegreeData = new MFilePlotData( "Average Node Degree","Diameter", "Average Node Degree vs Diameter",
				"AvgNodeDegreeVsDiameter.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgNodeDegree() , optimized.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgNodeDegree() , exponential.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgNodeDegree() , heuristic.get(i).getDiameter() ));
		}
		
		ValueWriter AvgNodeDegreePlotter = new ValueWriter(AvgNodeDegreeData);
		AvgNodeDegreePlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgNodeDegreeVsHubbiness() throws Exception{
		MFilePlotData AvgNodeDegreeData = new MFilePlotData( "Average Node Degree","Hubbiness", "Average Node Degree vs Hubbiness",
				"AvgNodeDegreeVsHubbiness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgNodeDegree() , optimized.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgNodeDegree() , exponential.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgNodeDegree() , heuristic.get(i).getHubbiness() ));
		}
		
		ValueWriter AvgNodeDegreePlotter = new ValueWriter(AvgNodeDegreeData);
		AvgNodeDegreePlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateAvgNodeDegreeVsPercentNodesConnectedToHubPlot() throws Exception{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "Avgerage Node Degree"," Percent Nodes Connected To Hub ", "Avgerage Node Degree vs Percent Nodes Connected To Hub",
				"AvgNodeDegreeVsPercentNodesConnectedToHub.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getAvgNodeDegree() , optimized.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getAvgNodeDegree() , exponential.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getAvgNodeDegree() , heuristic.get(i).getPercentNodesConnectedToHub() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateClusteringCoefficient() throws Exception 
	{
		MFilePlotData ClusteringCoefficientData = new MFilePlotData( "Average Clustering Coefficient","Fitness", "Average Clustering Coefficient vs Fitness",
				"ClusteringCoefficientVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getClusteringCoefficient() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getClusteringCoefficient() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getClusteringCoefficient() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter ClusteringCoefficientPlotter = new ValueWriter(ClusteringCoefficientData);
		ClusteringCoefficientPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateClusteringCoefficientVsDiameter()throws Exception {
		MFilePlotData ClusteringCoefficientData = new MFilePlotData( "Average Clustering Coefficient","Diameter", "Average Clustering Coefficient vs Diameter",
				"ClusteringCoefficientVsDiameter.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getClusteringCoefficient() , optimized.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getClusteringCoefficient() , exponential.get(i).getDiameter() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getClusteringCoefficient() , heuristic.get(i).getDiameter() ));
		}
		
		ValueWriter ClusteringCoefficientPlotter = new ValueWriter(ClusteringCoefficientData);
		ClusteringCoefficientPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateClusteringCoefficientVsHubbiness() throws Exception{
		MFilePlotData ClusteringCoefficientData = new MFilePlotData( "Average Clustering Coefficient","Hubbiness", "Average Clustering Coefficient vs Hubbiness",
				"ClusteringCoefficientVsHubbiness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getClusteringCoefficient() , optimized.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getClusteringCoefficient() , exponential.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getClusteringCoefficient() , heuristic.get(i).getHubbiness() ));
		}
		
		ValueWriter ClusteringCoefficientPlotter = new ValueWriter(ClusteringCoefficientData);
		ClusteringCoefficientPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateClusteringCoefficientVsPercentNodesConnectedToHubPlot() throws Exception{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "Clustering Coefficient"," Percent Nodes Connected To Hub ", "Clustering Coefficient vs Percent Nodes Connected To Hub",
				"ClusteringCoefficientVsPercentNodesConnectedToHub.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getClusteringCoefficient() , optimized.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getClusteringCoefficient() , exponential.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getClusteringCoefficient() , heuristic.get(i).getPercentNodesConnectedToHub() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateCostPlot() throws Exception 
	{
		MFilePlotData CostData = new MFilePlotData( "Cost","Fitness", "Cost vs Fitness",
				"CostVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getCost() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getCost() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getCost() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter CostPlotter = new ValueWriter(CostData);
		CostPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	private void generateCostBasisPlot() throws Exception 
	{
		MFilePlotData CostData = new MFilePlotData( "Cost Basis","Fitness", "Cost Basis vs Fitness",
				"CostBasisVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getCostBasis() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();

		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();

		
		ValueWriter CostPlotter = new ValueWriter(CostData);
		CostPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}
	
	private void generateDiameterPlot()throws Exception 
	{
		MFilePlotData DiameterData = new MFilePlotData( "Diameter","Fitness", "Diameter vs Fitness",
				"DiameterVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getDiameter() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getDiameter() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getDiameter() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter DiameterPlotter = new ValueWriter(DiameterData);
		DiameterPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateDiameterVsHubbiness() throws Exception{
		MFilePlotData DiameterData = new MFilePlotData( "Diameter","Hubbiness", "Diameter vs Hubbiness",
				"DiameterVsHubbiness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getDiameter() , optimized.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getDiameter() , exponential.get(i).getHubbiness() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getDiameter() , heuristic.get(i).getHubbiness() ));
		}
		
		ValueWriter DiameterPlotter = new ValueWriter(DiameterData);
		DiameterPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateDiameterVsPercentNodesConnectedToHubPlot() throws Exception{
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "Diameter"," Percent Nodes Connected To Hub ", "Diameter vs Percent Nodes Connected To Hub",
				"DiameterVsPercentNodesConnectedToHub.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getDiameter() , optimized.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getDiameter() , exponential.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getDiameter() , heuristic.get(i).getPercentNodesConnectedToHub() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateDirectDiameterPlot()throws Exception 
	{
		MFilePlotData DirectDiameterData = new MFilePlotData( "Direct Diameter","Fitness", "Direct Diameter vs Fitness",
				"DirectDiameterVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getDirectDiameter() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getDirectDiameter() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getDirectDiameter() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter DirectDiameterPlotter = new ValueWriter(DirectDiameterData);
		DirectDiameterPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateDirectLCCPlot()throws Exception 
	{
		MFilePlotData DirectLCCData = new MFilePlotData( "Direct LCC","Fitness", "Direct LCC vs Fitness",
				"DirectLCCVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getDirectLCC() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getDirectLCC() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getDirectLCC() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter DirectLCCPlotter = new ValueWriter(DirectLCCData);
		DirectLCCPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	public void generateHubbinessPlot() throws Exception 
	{
		MFilePlotData HubbinessData = new MFilePlotData( "Hubbiness","Fitness", "Hubbiness vs Fitness",
				"HubbinessVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getHubbiness() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getHubbiness() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getHubbiness() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter HubbinessPlotter = new ValueWriter(HubbinessData);
		HubbinessPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateHubbinessVsPercentNodesConnectedToHubPlot() throws Exception {
		//solnToExamine.getAISPL()
		MFilePlotData AISPLData = new MFilePlotData( "Hubbiness"," Percent Nodes Connected To Hub ", "Hubbiness vs Percent Nodes Connected To Hub",
				"HubbinessVsPercentNodesConnectedToHub.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getHubbiness() , optimized.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getHubbiness() , exponential.get(i).getPercentNodesConnectedToHub() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getHubbiness() , heuristic.get(i).getPercentNodesConnectedToHub() ));
		}
		
		ValueWriter AISPLPlotter = new ValueWriter(AISPLData);
		AISPLPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateLCCPlot() throws Exception 
	{
		MFilePlotData HubbinessData = new MFilePlotData( "LCC","Fitness", "LCC vs Fitness",
				"LCCVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getLCC() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getLCC() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getLCC() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter LCCPlotter = new ValueWriter(HubbinessData);
		LCCPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generatePercentNodesConnectedToHubPlot() throws Exception 
	{
		MFilePlotData PercentNodesConnectedToHubData = new MFilePlotData( "Percent Nodes Connected To Hubs","Fitness", "Percent Nodes Connected To Hubs vs Fitness",
				"PercentNodesConnectedToHubPlotVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getPercentNodesConnectedToHub() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getPercentNodesConnectedToHub() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getPercentNodesConnectedToHub() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter PercentNodesConnectedToHubPlotter = new ValueWriter(PercentNodesConnectedToHubData);
		PercentNodesConnectedToHubPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	public void generatePlots() throws Exception
	{
	
		if(exponential.size() == 0 || heuristic.size()==0 || optimized.size()==0)
			return;
		
		//Only value vs. fitness function
//		solnToExamine.getRandomDiameter()
		generateRandomDiameterPlot();
		
//		solnToExamine.getRandomLCC()
		generateRandomLCCPlot();
		
//		solnToExamine.getLCC()
		generateLCCPlot();
		
//		solnToExamine.getCost()
		generateCostPlot();
		
//		solnToExamine.getDirectDiameter()
		generateDirectDiameterPlot();
		
//		solnToExamine.getDirectLCC()
		generateDirectLCCPlot();
		
		
		//solnToExamine.getAISPL()
		generateAISPLPlot();
		generateAISPLPlotVsAssortativity();
		generateAISPPLlotVsAvgLinkCost();
		generateAISPPLlotVsAvgLinkLength();
		generateAISPPLlotVsAvgNodeDegree();
		generateAISPPLlotVsClusteringCoefficient();
		generateAISPPLlotVsDiameter();
		generateAISPPLlotVsHubbiness();
		generateAISPPlotVsPercentNodesConnectedToHubPlot();
		
		
		
		//solnToExamine.getASPL()
		generateASPLPlot();
		generateASPLPlotVsAssortativity();
		generateASPPLlotVsAvgLinkCost();
		generateASPPLlotVsAvgLinkLength();
		generateASPPLlotVsAvgNodeDegree();
		generateASPPLlotVsClusteringCoefficient();
		generateASPPLlotVsDiameter();
		generateASPPLlotVsHubbiness();
		generateASPPlotVsPercentNodesConnectedToHubPlot();
		
//		solnToExamine.getAssortativity()
		generateAssortativityPlot();
		generateAssortativityVsAvgLinkCost();
		generateAssortativityVsAvgLinkLength();
		generateAssortativityVsAvgNodeDegree();
		generateAssortativityVsClusteringCoefficient();
		generateAssortativityVsDiameter();
		generateAssortativityVsHubbiness();
		generateAssortativityVsPercentNodesConnectedToHubPlot();
		
//		solnToExamine.getAvgLinkCost()
		generateAvgLinkCostPlot();
		generateAvgLinkCostPlotVsAvgNodeDegree();
		generateAvgLinkCostPlotVsClusteringCoefficient();
		generateAvgLinkCostPlotVsDiameter();
		generateAvgLinkCostPlotVsHubbiness();
		generateAvgLinkCostVsPercentNodesConnectedToHubPlot();
		
//		solnToExamine.getAvgLinkLength()
		generateAvgLinkLengthPlot();
		generateAvgLinkLengthPlotVsAvgNodeDegree();
		generateAvgLinkLengthPlotVsClusteringCoefficient();
		generateAvgLinkLengthPlotVsDiameter();
		generateAvgLinkLengthPlotVsHubbiness();
		generateAvgLinkLengthVsPercentNodesConnectedToHubPlot();
		
//		solnToExamine.getAvgNodeDegree()
		generateAvgNodeDegree();
		generateAvgNodeDegreeVsClusteringCoefficient();
		generateAvgNodeDegreeVsDiameter();
		generateAvgNodeDegreeVsHubbiness();
		generateAvgNodeDegreeVsPercentNodesConnectedToHubPlot();
		
//		solnToExamine.getClusteringCoefficient()
		generateClusteringCoefficient();
		generateClusteringCoefficientVsDiameter();
		generateClusteringCoefficientVsHubbiness();
		generateClusteringCoefficientVsPercentNodesConnectedToHubPlot();
		

		
//		solnToExamine.getDiameter()
		generateDiameterPlot();
		generateDiameterVsHubbiness();
		generateDiameterVsPercentNodesConnectedToHubPlot();

		
//		solnToExamine.getFitnessValue()
		
//		solnToExamine.getHubbiness()	
		generateHubbinessPlot();
		generateHubbinessVsPercentNodesConnectedToHubPlot();
			
//		solnToExamine.getPercentNodesConnectedToHub()
		generatePercentNodesConnectedToHubPlot();
		

		
		//double fullFitness = layout.getFitnessOfSeedWithFullConnectivity();
		generateFitnessOfSeedWithFullConnectivityVsFitness();
		
		generateCostOfFullConnectivityVsFitness();
		generateCostBasisPlot();

	

		
	}

	

	private void generateCostOfFullConnectivityVsFitness() throws Exception{
		MFilePlotData HubbinessData = new MFilePlotData( "Cost Of Seed With Full Connectivity",
				"Fitness",
				"Cost Of Seed With Full Connectivity vs Fitness",
				"CostOfSeedWithFullConnectivityVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size(); i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getCostOfDirectConnectivity() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();

		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();

		
		ValueWriter LCCPlotter = new ValueWriter(HubbinessData);
		LCCPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
		
	}

	private void generateFitnessOfSeedWithFullConnectivityVsFitness() throws Exception
	{
		MFilePlotData HubbinessData = new MFilePlotData( "Fitness Of Seed With Full Connectivity",
				"Fitness",
				"Fitness Of Seed With Full Connectivity vs Fitness",
				"FitnessOfSeedWithFullConnectivityVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size(); i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getFitnessOfSeedWithFullConnectivity() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();

		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();

		
		ValueWriter LCCPlotter = new ValueWriter(HubbinessData);
		LCCPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateRandomDiameterPlot2() throws Exception
	{
		MFilePlotData RandomDiameterData = new MFilePlotData("", "Fitness"," Fitness",
				"MFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getRandomDiameter() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getRandomDiameter() , exponential.get(i).getFitnessValue() ));
		}
				
		Vector<PlotPairDoubleData> heuristic1Data = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heur1.size() ; i++)
		{
			heuristic1Data.add(new PlotPairDoubleData(heur1.get(i).getRandomDiameter() , heur1.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristic2Data = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heur2.size() ; i++)
		{
			heuristic2Data.add(new PlotPairDoubleData(heur2.get(i).getRandomDiameter() , heur2.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristic3Data = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heur3.size() ; i++)
		{
			heuristic3Data.add(new PlotPairDoubleData(heur3.get(i).getRandomDiameter() , heur3.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristic4Data = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heur3.size() ; i++)
		{
			heuristic4Data.add(new PlotPairDoubleData(heur3.get(i).getRandomDiameter() , heur3.get(i).getFitnessValue() ));
		}
		
		ValueWriter RandomDiameterPlotter = new ValueWriter(RandomDiameterData);
		RandomDiameterPlotter.plotDouble2(optimizedData, exponentialData, heuristic1Data, heuristic2Data,heuristic3Data, heuristic4Data);
		
	}
	
	private void generateRandomDiameterPlot() throws Exception
	{
		MFilePlotData RandomDiameterData = new MFilePlotData("Random Diameter", "Fitness","Random Diameter vs Fitness",
				"RandomDiameterVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getRandomDiameter() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getRandomDiameter() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getRandomDiameter() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter RandomDiameterPlotter = new ValueWriter(RandomDiameterData);
		RandomDiameterPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	private void generateRandomLCCPlot() throws Exception 
	{
		MFilePlotData RandomLCCData = new MFilePlotData( "Random LCC","Fitness", "Random LCC vs Fitness",
				"RandomLCCVsFitness.m", outputDirectory, null);
		
		Vector<PlotPairDoubleData> optimizedData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < optimized.size() ; i++)
		{
			optimizedData.add(new PlotPairDoubleData(optimized.get(i).getRandomLCC() , optimized.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> exponentialData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < exponential.size() ; i++)
		{
			exponentialData.add(new PlotPairDoubleData(exponential.get(i).getRandomLCC() , exponential.get(i).getFitnessValue() ));
		}
		
		Vector<PlotPairDoubleData> heuristicData = new Vector<PlotPairDoubleData>();
		for(int i =0 ; i < heuristic.size() ; i++)
		{
			heuristicData.add(new PlotPairDoubleData(heuristic.get(i).getRandomLCC() , heuristic.get(i).getFitnessValue() ));
		}
		
		ValueWriter RandomLCCPlotter = new ValueWriter(RandomLCCData);
		RandomLCCPlotter.plotDouble(optimizedData, exponentialData, heuristicData);
		
	}

	public void generatePlots2() throws Exception
	{

		
		//Only value vs. fitness function
//		solnToExamine.getRandomDiameter()
		generateRandomDiameterPlot2();

		
	}

//	public void generatePlotsH() 
//	{
//		//Only value vs. fitness function
////		solnToExamine.getRandomDiameter()
//		generateRandomDiameterPlotH();
//		
////		solnToExamine.getRandomLCC()
//		generateRandomLCCPlotH();
//		
////		solnToExamine.getLCC()
//		generateLCCPlotH();
//		
////		solnToExamine.getCost()
//		generateCostPlotH();
//		
////		solnToExamine.getDirectDiameter()
//		generateDirectDiameterPlotH();
//		
////		solnToExamine.getDirectLCC()
//		generateDirectLCCPlotH();
//		
//		
//		//solnToExamine.getAISPL()
//		generateAISPLPlotH();
//		generateAISPLPlotVsAssortativityH();
//		generateAISPPLlotVsAvgLinkCostH();
//		generateAISPPLlotVsAvgLinkLengthH();
//		generateAISPPLlotVsAvgNodeDegreeH();
//		generateAISPPLlotVsClusteringCoefficientH();
//		generateAISPPLlotVsDiameterH();
//		generateAISPPLlotVsHubbinessH();
//		generateAISPPlotVsPercentNodesConnectedToHubPlotH();
//		
//		
//		
//		//solnToExamine.getASPL()
//		generateASPLPlotH();
//		generateASPLPlotVsAssortativityH();
//		generateASPPLlotVsAvgLinkCostH();
//		generateASPPLlotVsAvgLinkLengthH();
//		generateASPPLlotVsAvgNodeDegreeH();
//		generateASPPLlotVsClusteringCoefficientH();
//		generateASPPLlotVsDiameterH();
//		generateASPPLlotVsHubbinessH();
//		generateASPPlotVsPercentNodesConnectedToHubPlotH();
//		
////		solnToExamine.getAssortativity()
//		generateAssortativityPlotH();
//		generateAssortativityVsAvgLinkCostH();
//		generateAssortativityVsAvgLinkLengthH();
//		generateAssortativityVsAvgNodeDegreeH();
//		generateAssortativityVsClusteringCoefficientH();
//		generateAssortativityVsDiameterH();
//		generateAssortativityVsHubbinessH();
//		generateAssortativityVsPercentNodesConnectedToHubPlotH();
//		
////		solnToExamine.getAvgLinkCost()
//		generateAvgLinkCostPlotH();
//		generateAvgLinkCostPlotVsAvgNodeDegreeH();
//		generateAvgLinkCostPlotVsClusteringCoefficientH();
//		generateAvgLinkCostPlotVsDiameterH();
//		generateAvgLinkCostPlotVsHubbinessH();
//		generateAvgLinkCostVsPercentNodesConnectedToHubPlotH();
//		
////		solnToExamine.getAvgLinkLength()
//		generateAvgLinkLengthPlotH();
//		generateAvgLinkLengthPlotVsAvgNodeDegreeH();
//		generateAvgLinkLengthPlotVsClusteringCoefficientH();
//		generateAvgLinkLengthPlotVsDiameterH();
//		generateAvgLinkLengthPlotVsHubbinessH();
//		generateAvgLinkLengthVsPercentNodesConnectedToHubPlotH();
//		
////		solnToExamine.getAvgNodeDegree()
//		generateAvgNodeDegreeH();
//		generateAvgNodeDegreeVsClusteringCoefficientH();
//		generateAvgNodeDegreeVsDiameterH();
//		generateAvgNodeDegreeVsHubbinessH();
//		generateAvgNodeDegreeVsPercentNodesConnectedToHubPlotH();
//		
////		solnToExamine.getClusteringCoefficient()
//		generateClusteringCoefficientH();
//		generateClusteringCoefficientVsDiameterH();
//		generateClusteringCoefficientVsHubbinessH();
//		generateClusteringCoefficientVsPercentNodesConnectedToHubPlotH();
//		
//
//		
////		solnToExamine.getDiameter()
//		generateDiameterPlotH();
//		generateDiameterVsHubbinessH();
//		generateDiameterVsPercentNodesConnectedToHubPlotH();
//
//		
////		solnToExamine.getFitnessValue()
//		
////		solnToExamine.getHubbiness()	
//		generateHubbinessPlotH();
//		generateHubbinessVsPercentNodesConnectedToHubPlotH();
//			
////		solnToExamine.getPercentNodesConnectedToHub()
//		generatePercentNodesConnectedToHubPlotH();
//		
//
//		
//		//double fullFitness = layout.getFitnessOfSeedWithFullConnectivity();
//		generateFitnessOfSeedWithFullConnectivityVsFitnessH();
//		
//		generateCostOfFullConnectivityVsFitnessH();
//		generateCostBasisPlotH();
//		
//	}


}
