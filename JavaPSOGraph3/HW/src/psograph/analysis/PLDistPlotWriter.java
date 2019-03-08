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
import psograph.util.DistWriter2;
import psograph.util.MFilePlotData;


public class PLDistPlotWriter
{
File outputDirectory;
	
	Vector<CalculatedGraph> optimized;
	Vector<CalculatedGraph> exponential;
	Vector<CalculatedGraph> heuristic;
	
	Vector<CalculatedGraph> heur1;
	Vector<CalculatedGraph> heur2;
	Vector<CalculatedGraph> heur3;
	
	String layoutTitle;
	
	
	public PLDistPlotWriter(Vector<CalculatedGraph> optimized, Vector<CalculatedGraph> exponential,
			Vector<CalculatedGraph> heuristic, File OutDir)
	{
		this.optimized = optimized;
		this.exponential = exponential;
		this.heuristic = heuristic;
		outputDirectory = OutDir;
	}
	
	public PLDistPlotWriter(Vector<CalculatedGraph> optSolnsAll,
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

	public void generatePlots() throws Exception
	{
	
		if(exponential.size() == 0 || heuristic.size()==0 || optimized.size()==0)
			return;
		

	generateISPL();
		

//		generateSPL();
		

		
	}

	private void generateSPL() throws Exception
	{
		MFilePlotData SPLData = new MFilePlotData( "","SPL", "SPL",
				"SPL.m", outputDirectory, null, 0, 20);
		
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
		
		DistWriter2 SPLPlotter = new DistWriter2(SPLData);
		SPLPlotter.PlotDistInt(optimizedData, exponentialData, heuristicData);
		
	}


	private void generateISPL()throws Exception 
	{
		MFilePlotData ISPLData = new MFilePlotData( "","ISPL", "ISPL",
				"ISPL.m", outputDirectory, null, 0.0, 1.0);
		
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
		
		DistWriter2 ISPLPlotter = new DistWriter2(ISPLData);
		ISPLPlotter.PlotDist(optimizedData, exponentialData, heuristicData);
		
	}

	public void generatePlotsH() {
		// TODO Auto-generated method stub
		
	}

	

}
