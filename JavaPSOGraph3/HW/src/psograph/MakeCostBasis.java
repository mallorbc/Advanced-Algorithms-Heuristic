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
import java.util.Vector;

import psograph.graph.CalculatedGraph;
import psograph.graph.Graph;
import psograph.graph.NodeLocationCalculator;
import psograph.util.Util;



public class MakeCostBasis {


	
	/**
	 * @param args
	 */
	public static void main(String[] args) 
	{

		String BaseDir0 = "C:\\Comparisions\\Seed0";
		String BaseDir1 = "C:\\Comparisions\\Seed1";
		String BaseDir2 = "C:\\Comparisions\\Seed2";
		String BaseDir3 = "C:\\Comparisions\\Seed3";
		String BaseDir4 = "C:\\Comparisions\\Seed4";
		String BaseDir5 = "C:\\Comparisions\\Seed5";
		
		String SeedDir0 = BaseDir0 + "\\seed";
		
		String SeedDir1 = BaseDir1 + "\\seed";
		
		String SeedDir2 = BaseDir2 + "\\seed";
		
		String SeedDir3 = BaseDir3 + "\\seed";
		
		String SeedDir4 = BaseDir4 + "\\seed";
		
		String SeedDir5 = BaseDir5 + "\\seed";
		
		Vector<String> dirs = new Vector<String>();
		dirs.add(SeedDir0);
		dirs.add(SeedDir1);
		dirs.add(SeedDir2);
		dirs.add(SeedDir3);
		dirs.add(SeedDir4);
		dirs.add(SeedDir5);
		
		
		
		try
		{
			System.out.println("getting CG1");

			for(int i=0; i < dirs.size(); i++)
			{
				System.out.println("Working on "+dirs.get(i));
				
				Graph costBasisGraph = Util.streaminGraph(new File(dirs.get(i)+"\\ExponentialCostBasis1.Graph"));

				NodeLocationCalculator calc = new NodeLocationCalculator(costBasisGraph);
				calc.calculateResults();
				
				System.out.println("Cost of cost Basis is "+calc.getCost());
				System.out.println("num of edges is "+calc.getTotalEdges());
				
				
				Util.streamOutExponentialGraph(new File(dirs.get(i)), calc, 0);

			}


		}
		catch (Exception e)
		{
			e.printStackTrace();

		}

	}
		
}

	

