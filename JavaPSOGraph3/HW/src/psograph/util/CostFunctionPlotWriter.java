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
public class CostFunctionPlotWriter extends MfileWriter {

		public CostFunctionPlotWriter()
		{
			super();
			
			m_baseFileNames = new String[3];
			m_baseFileNames[0] ="m_AISPLVsFitness";	
			m_baseFileNames[1] ="m_CostFactorVsFitness";
			m_baseFileNames[2] ="m_AISPLVsCostFactor";

	    	
			
			
			m_labels = new Vector<GraphLabel>();
			m_labels.add(new GraphLabel("Fitness", "AISPL vs. Fitness"));
			m_labels.add(new GraphLabel("Fitness", "(1- Cost Factor)^3 vs. Fitness"));
			m_labels.add(new GraphLabel("(1- Cost Factor)^3", "Scatterplot of AISPL vs. (1- Cost Factor)^3  "));
  
			
		}

		void printDataEntryCalculatedMfile(Vector<BufferedWriter> writers, CalculatedGraph calcGraph, boolean is_last) throws Exception
		{		
			double fitness = calcGraph.getFitnessValue();
			
			//Note Order is very important do not move seq around without considering how the 
			// writers are added to the Vector
			BufferedWriter writer;
			
			writer = writers.get(0);
			writer.write(calcGraph.getAISPL()+", " + fitness +";\n");
			
			
			double cost = calcGraph.getCost();
			double cost_factor =  1 - (cost/calcGraph.getCostBasis());		
			
			//To make the values spread out more, we are taking the (1-CF)^3
			double tt = Math.pow(cost_factor,3);
		
			writer = writers.get(1);
			writer.write(tt +", " + fitness +";\n");		
			
			writer = writers.get(2);
			writer.write(calcGraph.getAISPL()+", " + tt +";\n");
			
    	

		
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
				
				writer.write("plot(x,y,'d'); \n");
				
				
				
				if(i==0)
					writer.write("xlabel('AISPL'); \n");
				else if(i==1)
					writer.write("xlabel('1 - cost factor'); \n");
				else if(i==2)
					writer.write("xlabel('AISPL'); \n");
				
				writer.write("ylabel('"+ m_labels.get(i).getYLabel() +"'); \n");
				writer.write("title('"+ m_labels.get(i).getGraphLabel() +"'); \n");
				
				//m_randomDataMFile.write("gplot(A,xyz,'r+');  \n");		
				
				writer.close();
				i++;
		    }

		}
		
		
}
