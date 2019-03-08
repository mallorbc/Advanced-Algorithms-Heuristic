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
package psograph.util;


import java.io.File;
import java.io.FilenameFilter;

import psograph.graph.CalculatedGraph;




public class GraphUpdater {

    // It is also possible to filter the list of returned files.
    // This example does not return any files that start with `.'.
    FilenameFilter filter = new FilenameFilter() {
        public boolean accept(File dir, String name) {
            return name.endsWith(".CalculatedGraph");
        }
    };
    public void processSeedDirectory(File seedDirectory) throws Exception
    {
		File graphDirectory = new File(seedDirectory.getAbsolutePath()+"\\graphs");
		
		if(!graphDirectory.exists())
		{
			return;			
		}
    	
		File calculatedGraphs[] = graphDirectory.listFiles(filter);
		if (calculatedGraphs == null)
			throw new Exception ("SeedDirectory is empty");
		
		System.out.println("Working on SeedDirectory - "+ graphDirectory.getAbsolutePath());
    	
    	int j;

	

		int maxCalculatedGraphs = calculatedGraphs.length;
		//for loop to go over each .calculatedGraph file
		
		for(j=0; j < maxCalculatedGraphs; j++)
		{
			CalculatedGraph graph =  null;

			System.out.println("CalculatedGraph was retrieved Graph " + calculatedGraphs[j].toString());
			graph = psograph.util.Util.streaminCalculatedGraph(calculatedGraphs[j]);
			
			//graph.UpdateCalcuations();
			graph.UpdatePSOCalculations();
		
			Util.streamOutCalculatedGraphWithExtension( graphDirectory,  calculatedGraphs[j].getName(), graph);

		}
		
    }
}
