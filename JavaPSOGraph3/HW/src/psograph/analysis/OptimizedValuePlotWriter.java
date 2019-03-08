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
import psograph.util.DistWriter;
import psograph.util.MFilePlotData;


public class OptimizedValuePlotWriter  {

	File outputDirectory;
	
	Vector<Vector<CalculatedGraph>> optimalSolns;
	
	String layoutTitle;
	
	
	
	public OptimizedValuePlotWriter(
			Vector<Vector<CalculatedGraph>> optimalSolns, File parent) 
	{
		this.optimalSolns= optimalSolns;
		outputDirectory = parent;
	}

	public void generatePlots() throws Exception{
		generateLinkCosts();
		generateLinkLengths();
		generateSPLS();
		generateISPLS();


		
	}

	private void generateISPLS() throws Exception{
		MFilePlotData LinkCostsData = new MFilePlotData( "Inverse Shortest Path Length",
				"Number of Inverse Shortest Path Length", 
				"Histogram Plots of Inverse Shortest Path Lengths Distribution of Optimized Networks",
				"OptimizedISPLs.m", outputDirectory, null, 8);
		
		Vector<Vector<TreeMap<Double,Integer>>> data = new Vector<Vector<TreeMap<Double,Integer>>>();
		//System.out.println("optimal Solns size "+optimalSolns.size());
		for(int i =0; i < optimalSolns.size(); i++)
		{
			data.add(i, new Vector<TreeMap<Double,Integer>>());
			for(int j =0 ; j < optimalSolns.get(i).size() ; j++)
			{
				TreeMap<Double,Integer> dist = optimalSolns.get(i).get(j).getDistributionISPL();
				data.get(i).add(j,dist);
			}
		}
		

		
		
		DistWriter LinkCostsPlotter = new DistWriter(LinkCostsData);
		LinkCostsPlotter.PlotDist(data);
		
	}

	private void generateSPLS()throws Exception 
	{
		MFilePlotData LinkCostsData = new MFilePlotData( "","Shortest Path Lengths", "Shortest Path  Lengths",
				"OptimizedSPLs.m", outputDirectory, null, 20);
		
		Vector<Vector<TreeMap<Integer,Integer>>> data = new Vector<Vector<TreeMap<Integer,Integer>>>();
		//System.out.println("optimal Solns size "+optimalSolns.size());
		for(int i =0; i < optimalSolns.size(); i++)
		{
			data.add(i, new Vector<TreeMap<Integer,Integer>>());
			for(int j =0 ; j < optimalSolns.get(i).size() ; j++)
			{
				TreeMap<Integer,Integer> dist = optimalSolns.get(i).get(j).getDistributionSPL();
				data.get(i).add(j,dist);
			}
		}
		

		
		
		DistWriter LinkCostsPlotter = new DistWriter(LinkCostsData);
		LinkCostsPlotter.PlotDistInt(data);
		
	}

	private void generateLinkLengths()throws Exception {
		MFilePlotData LinkCostsData = new MFilePlotData( "","Link Lengths", "Link Lengths",
				"OptimizedLinkLengths.m", outputDirectory, null, 20);
		
		Vector<Vector<TreeMap<Double,Integer>>> data = new Vector<Vector<TreeMap<Double,Integer>>>();
		//System.out.println("optimal Solns size "+optimalSolns.size());
		for(int i =0; i < optimalSolns.size(); i++)
		{
			data.add(i, new Vector<TreeMap<Double,Integer>>());
			for(int j =0 ; j < optimalSolns.get(i).size() ; j++)
			{
				TreeMap<Double,Integer> dist = optimalSolns.get(i).get(j).getDistributionLinkLengths();
				data.get(i).add(j,dist);
			}
		}
		

		
		
		DistWriter LinkCostsPlotter = new DistWriter(LinkCostsData);
		LinkCostsPlotter.PlotDist(data);
		
	}

	private void generateLinkCosts()throws Exception {
		
		
		MFilePlotData LinkCostsData = new MFilePlotData( "","Link Costs", "Link Costs",
				"OptimizedLinkCosts.m", outputDirectory, null, 20);
		
		Vector<Vector<TreeMap<Double,Integer>>> data = new Vector<Vector<TreeMap<Double,Integer>>>();
		System.out.println("optimal Solns size "+optimalSolns.size());
		for(int i =0; i < optimalSolns.size(); i++)
		{
			data.add(i, new Vector<TreeMap<Double,Integer>>());
			for(int j =0 ; j < optimalSolns.get(i).size() ; j++)
			{
				TreeMap<Double,Integer> dist = optimalSolns.get(i).get(j).getDistributionLinkCosts();
				data.get(i).add(j,dist);
			}
		}
		

		
		
		DistWriter LinkCostsPlotter = new DistWriter(LinkCostsData);
		LinkCostsPlotter.PlotDist(data);
		
	}
	
	
}
