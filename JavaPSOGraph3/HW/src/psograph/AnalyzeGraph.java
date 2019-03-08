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

import psograph.analysis.CreateFitnessCompHistograms;
import psograph.analysis.CreateHistogram;
import psograph.analysis.FitnessFunctionGraphs;
import psograph.analysis.GenerateFitnessFunctionDistribution;


public class AnalyzeGraph 
{


    public static void main(String[] args) 
    {
    	try
    	{
    		
//    		System.out.println("Performing BasicPlots");
//    		BasicDataPlots basicDataPlots = new BasicDataPlots();
//    		basicDataPlots.generatePlots();
//    		System.out.println("Done Performing BasicPlots");
    		
//    		System.out.println("Performing CreateFitnessCompHistograms");
//      	CreateFitnessCompHistograms createFitnessCompHistograms = new CreateFitnessCompHistograms();
//			createFitnessCompHistograms.generateHistorgrams();
//			System.out.println("Done Performing CreateFitnessCompHistograms");
//			
//			System.out.println("Performing CreateHistogram");
//			CreateHistogram createHistogram = new CreateHistogram();
//			createHistogram.generateHistorgrams();
//			System.out.println("Done Performing CreateHistogram");
			
//			System.out.println("Performing FitnessFunctionGraphs");
//    		FitnessFunctionGraphs fitnessFunctionGraphs = new FitnessFunctionGraphs();
//    		fitnessFunctionGraphs.generatePlots();
//    		System.out.println("Done Performing FitnessFunctionGraphs");
//
    		System.out.println("Performing GenerateFitnessFunctionDistribution");
    		GenerateFitnessFunctionDistribution generateFitnessFunctionDistribution = new GenerateFitnessFunctionDistribution();
    		generateFitnessFunctionDistribution.generateDist();
    		System.out.println("Done Performing GenerateFitnessFunctionDistribution");
    		

    		
    		
    		
    		
    	}
    	catch (Exception e)
    	{
    		e.printStackTrace();
    	}
    }

	


}

