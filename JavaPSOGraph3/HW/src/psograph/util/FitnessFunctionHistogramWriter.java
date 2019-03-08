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

import java.io.BufferedWriter;
import java.util.Iterator;
import java.util.Vector;

import psograph.graph.CalculatedGraph;



public class FitnessFunctionHistogramWriter extends MfileWriter
{

	int numBins = 0; 
	

	
	public FitnessFunctionHistogramWriter(int nBins)
	{
		super();
		m_baseFileNames = new String[1];
		m_baseFileNames[0] = "fitnessHistogram";
		
		
    	m_labels = new Vector<GraphLabel>();
    	m_labels.add(new GraphLabel("", "Distribution of Fitness Values"));
		
    	numBins = nBins;
		
	}
	

	
	void printDataEntryCalculatedMfile(Vector<BufferedWriter> writers, CalculatedGraph calcGraph, boolean is_last) throws Exception
	{		
		
		
		//Note Order is very important do not move seq around without considering how the 
		// writers are added to the Vector
		BufferedWriter writer;
		
		writer = writers.get(0);
		if(!is_last)
			writer.write(  calcGraph.getFitnessValue() +",");
		else
			writer.write(" "+  calcGraph.getFitnessValue() );
		
		
		
	}
	void printEndOfCalculatedMFile(Vector<BufferedWriter> writers, int numGraphs) throws Exception
	{
		Iterator<BufferedWriter> iter = writers.iterator();
		
		int i =0;
	    for (; iter.hasNext(); ) {
	    	BufferedWriter writer = iter.next();
			writer.write("]; \n");

			writer.write("hist(xy," + numBins +  "); \n");
			
			writer.write("xlabel('fitness value'); \n");
			writer.write("ylabel('"+ m_labels.get(i).getYLabel() +"'); \n");
			writer.write("title('"+ m_labels.get(i).getGraphLabel() +"'); \n");
			
			//m_randomDataMFile.write("gplot(A,xyz,'r+');  \n");		
			
			writer.close();
			i++;
	    }

	}
	

	

	
	
	
}
