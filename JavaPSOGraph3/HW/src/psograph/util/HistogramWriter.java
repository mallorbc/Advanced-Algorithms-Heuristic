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
import java.io.FileWriter;
import java.util.TreeMap;
import java.util.Vector;


public class HistogramWriter 
{
	MFilePlotData plotData;
	BufferedWriter m_bufWriter;
	String m_baseFileNames[];
	Vector<GraphLabel> m_labels;
	int numBins = 0;
	
	public HistogramWriter(MFilePlotData PlotData)
	{
		plotData = PlotData;
	}


	
	void printEndOfMFile() throws Exception
	{
		m_bufWriter.write("]; \n");
		
		//m_bufWriter.write("x = xy(:,1) \n");
		//m_bufWriter.write("y = xy(:,2) \n");
		
		
		
		m_bufWriter.write("hist(xy, "+numBins+"); \n");
		
		m_bufWriter.write("xlabel('"+plotData.xAxisLabel+"'); \n");
		m_bufWriter.write("ylabel('"+plotData.yAxisLabel+"'); \n");
		if(plotData.ofWhat ==null)
		{
			m_bufWriter.write("title('"+plotData.plotLabel+"'); \n");
		}
		else
		{
			m_bufWriter.write("title('"+plotData.plotLabel+plotData.ofWhat+"'); \n");
		}
		
	
	}
	
	void printStartOfMFile() throws Exception
	{
    	m_bufWriter.write("xy = [ \n");
	}
	
	private void BaseShutDown() throws Exception
	{
		
		m_bufWriter.close();
	}
	
	private void BaseSetupForPlot() throws Exception
	{
		m_bufWriter = 
			new BufferedWriter(new FileWriter(plotData.outputDirectory.getAbsolutePath()+"\\"+plotData.fileName));
		
		
	}

	public void plotIntegerInteger(TreeMap<Integer, Integer> data) throws Exception
	{
		BaseSetupForPlot();
		printStartOfMFile();
		printDataIntegerInteger(data);
		printEndOfMFile();
		BaseShutDown();
		
	}

	private void printDataIntegerInteger(TreeMap<Integer, Integer> data) throws Exception
	{
		Vector<Integer> vKeys = new Vector<Integer>(data.keySet());
		for(int i =0; i < vKeys.size(); i++)
		{
			for(int j =0 ; j < data.get(vKeys.get(i));j++)
			{
				m_bufWriter.write(vKeys.get(i) +";");
			}
		}
		
		numBins = vKeys.lastElement() + 1;
		
	}




	public void plotDoubleInteger(TreeMap<Double, Integer> data) throws Exception
	{
		
		BaseSetupForPlot();
		printStartOfMFile();
		printDataDoubleInteger(data);
		printEndOfMFile();
		BaseShutDown();
		
	}






	private void printDataDoubleInteger(TreeMap<Double, Integer> data) throws Exception{
		
		Vector<Double> vKeys = new Vector<Double>(data.keySet());
		for(int i =0; i < vKeys.size(); i++)
		{   
			for(int j =0 ; j < data.get(vKeys.get(i));j++)
			{
				m_bufWriter.write(vKeys.get(i) +";");
			}
		}
		numBins = vKeys.size() + 1;
		
	}



}
