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


/** 
 * @deprecated
 * @author Owner
 *
 */
public class BasicMeasurementsPlotWriter extends MfileWriter 
{

	public BasicMeasurementsPlotWriter()
	{
		super();
		
		m_baseFileNames = new String[12];
		m_baseFileNames[0] ="m_randomDiameter";	
		m_baseFileNames[1] ="m_directDiameter";
		m_baseFileNames[2] ="m_directLCC";
		m_baseFileNames[3] ="m_randomLCC" ;
		m_baseFileNames[4] ="m_AISPLFile" ;    	
		m_baseFileNames[5] ="m_costBasisFile" ;
		m_baseFileNames[6] ="m_diameterFile";
		m_baseFileNames[7] ="m_LCCFile";
		m_baseFileNames[8] ="m_clusteringCoefficientFile";
		m_baseFileNames[9] ="m_avgNodeDegree";
		m_baseFileNames[10] ="m_avgLinkCost";
		m_baseFileNames[11] ="m_avgEdgeLength";
    	
		
		
		m_labels = new Vector<GraphLabel>();
		m_labels.add(new GraphLabel("Diameter", "Diameter after random attack 15%"));
		m_labels.add(new GraphLabel("Diameter", "Diameter after direct attack 15%"));
		m_labels.add(new GraphLabel("LCC", "Largest Connected Component (LCC) after direct attack 20%"));
		m_labels.add(new GraphLabel("LCC", "Largest Connected Component (LCC) after random attack 20%"));
		m_labels.add(new GraphLabel("AISPL", "Average inverse shortest path length"));
		m_labels.add(new GraphLabel("Cost ratio (graph cost/cost basis)","Cost Ratio"));
		m_labels.add(new GraphLabel("Diameter", "Diameter"));
		m_labels.add(new GraphLabel("LCC", "Largest Connected Component (LCC)"));	    	
		m_labels.add(new GraphLabel("Clustering Coefficient", "Clustering Coefficient Largest"));	  
		m_labels.add(new GraphLabel("Avg. Node Degree", "Avg. Node Degree"));
		m_labels.add(new GraphLabel("Avg. Link Cost", "Avg. Link Cost"));	    	
		m_labels.add(new GraphLabel("Avg. Link Length", "Avg. Link Length"));	
		
	}

	void printDataEntryCalculatedMfile(Vector<BufferedWriter> writers, CalculatedGraph calcGraph, boolean is_last) throws Exception
	{		
		double fitness = calcGraph.getFitnessValue();
		
		//Note Order is very important do not move seq around without considering how the 
		// writers are added to the Vector
		BufferedWriter writer;
		
		writer = writers.get(0);
		writer.write(fitness +", " + calcGraph.getRandomDiameter() +";\n");
		
		writer = writers.get(1);
		writer.write(fitness +", " + calcGraph.getDirectDiameter() +";\n");		
		
		writer = writers.get(2);
		writer.write(fitness +", " + calcGraph.getRandomLCC() +";\n");	    	

		writer = writers.get(3);
		writer.write(fitness +", " + calcGraph.getDirectLCC() +";\n");	    	

    	
		writer = writers.get(4);
		writer.write(fitness +", " + calcGraph.getAISPL() +";\n");	    	

		writer = writers.get(5);
		writer.write(fitness +", " + Util.calculateCostBasisRation(calcGraph) +";\n");	    	

		writer = writers.get(6);
		writer.write(fitness +", " + calcGraph.getDiameter() +";\n");	    	

		writer = writers.get(7);
		writer.write(fitness +", " + calcGraph.getLCC() +";\n");	 	

		writer = writers.get(8);
		writer.write(fitness +", " + calcGraph.getClusteringCoefficient() +";\n");
		
		writer = writers.get(9);
		writer.write(fitness +", " + calcGraph.getAvgNodeDegree() +";\n");	    	

		writer = writers.get(10);
		writer.write(fitness +", " + calcGraph.getAvgLinkCost() +";\n");	 	

		writer = writers.get(11);
		writer.write(fitness +", " + calcGraph.getAvgLinkLength() +";\n");
		
	}
	
	void printEndOfCalculatedMFile(Vector<BufferedWriter> writers, int numGraphs) throws Exception
	{
		Iterator<BufferedWriter> iter = writers.iterator();
		
		int i =0;
	    for (; iter.hasNext(); ) {
	    	BufferedWriter writer = iter.next();
			writer.write("]; \n");
			
			writer.write("x = xy(:,1) \n");
			writer.write("y = xy(:,2) \n");
			
			writer.write("plot(x,y,'r+'); \n");
			
			writer.write("xlabel('fitness value'); \n");
			writer.write("ylabel('"+ m_labels.get(i).getYLabel() +"'); \n");
			writer.write("title('"+ m_labels.get(i).getGraphLabel() +"'); \n");
			
			//m_randomDataMFile.write("gplot(A,xyz,'r+');  \n");		
			
			writer.close();
			i++;
	    }

	}
	
	
}
