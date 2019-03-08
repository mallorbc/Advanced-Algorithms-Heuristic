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

import java.util.Vector;


public class ValueWriter {
	MFilePlotData plotData;
	BufferedWriter m_bufWriter;
	String m_baseFileNames[];
	Vector<GraphLabel> m_labels;
	int numBins = 0;
	
	public ValueWriter(MFilePlotData PlotData)
	{
		plotData = PlotData;
	}

	private void printStartOfMFile2() {
		// TODO Auto-generated method stub
		
	}

	void printStartOfMFile() throws Exception
	{
    	
	}
	
	void printStartOfDataArray(String name) throws Exception
	{
		m_bufWriter.write(name+" = [ \n");
	}
	
	void printEndOfDataArray() throws Exception
	{
		m_bufWriter.write("]; \n");
	}
	
	public void plotDouble(Vector<PlotPairDoubleData> optimizedData,
			Vector<PlotPairDoubleData> exponentialData,
			Vector<PlotPairDoubleData> heuristicData) throws Exception
	{
		BaseSetupForPlot();
		printStartOfMFile();
		
		printStartOfDataArray("A");
		printDataDoublePair(optimizedData);
		printEndOfDataArray();
		
		printStartOfDataArray("B");
		printDataDoublePair(exponentialData);
		printEndOfDataArray();
		
		printStartOfDataArray("C");
		printDataDoublePair(heuristicData);
		printEndOfDataArray();
		
		printEndOfMFile("A , B, C");
		BaseShutDown();
		
	}

	private void printDataDoublePair(Vector<PlotPairDoubleData> data) throws Exception
	{

		for(int i =0; i < data.size(); i++)
		{   
			m_bufWriter.write(data.get(i).m_x +" , " +data.get(i).m_y+ " ; \n");
		}
		
		
	}

	private void BaseShutDown() throws Exception
	{
		m_bufWriter.close();
		
	}

	private void BaseShutDown2() throws Exception{
		m_bufWriter.close();
		
	}

	void printEndOfMFile(String arrayNames) throws Exception
	{
		
		
		m_bufWriter.write("Ax = A(:,1) \n");
		m_bufWriter.write("Ay = A(:,2) \n");
		m_bufWriter.write("Bx = B(:,1) \n");
		m_bufWriter.write("By = B(:,2) \n");		
		m_bufWriter.write("Cx = C(:,1) \n");
		m_bufWriter.write("Cy = C(:,2) \n");		
		
		m_bufWriter.write("plot(Ax,Ay,'rd',Bx, By,'b+',Cx,Cy,'gx' ); \n");
		m_bufWriter.write("legend('Optimized','Exponential','Heuristic'); \n");
		
		
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

	private void printEndOfMFile2() throws Exception{

		
		m_bufWriter.write("Ax = A(:,1) \n");
		m_bufWriter.write("Ay = A(:,2) \n");
		m_bufWriter.write("Bx = B(:,1) \n");
		m_bufWriter.write("By = B(:,2) \n");		
		m_bufWriter.write("Cx = C(:,1) \n");
		m_bufWriter.write("Cy = C(:,2) \n");	
		m_bufWriter.write("Dx = D(:,1) \n");
		m_bufWriter.write("Dy = D(:,2) \n");	
		
		m_bufWriter.write("plot(Ax,Ay,'rd',Bx, By,'b+',Cx,Cy,'gx', Dx, Dy, 'ro' ); \n");
		m_bufWriter.write("legend('Optimized','Heuristic1','Heuristic2','Heuristic3'); \n");
		
		
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

	private void BaseSetupForPlot2() throws Exception{
		
		m_bufWriter = 
			new BufferedWriter(new FileWriter(plotData.outputDirectory.getAbsolutePath()+"\\"+plotData.fileName));
		
	}

	private void BaseSetupForPlot() throws Exception
	{
		m_bufWriter = 
			new BufferedWriter(new FileWriter(plotData.outputDirectory.getAbsolutePath()+"\\"+plotData.fileName));
		
		System.out.println("file is at"+ plotData.outputDirectory.getAbsolutePath()+"\\"+plotData.fileName);
		
	}

	public void plotDouble2(Vector<PlotPairDoubleData> optimizedData,
			Vector<PlotPairDoubleData> exponentialData,
			Vector<PlotPairDoubleData> heuristic1Data,
			Vector<PlotPairDoubleData> heuristic2Data,
			Vector<PlotPairDoubleData> heuristic3Data,
			Vector<PlotPairDoubleData> heuristic4Data)  throws Exception
	{
		BaseSetupForPlot2();
		printStartOfMFile2();
		
		printStartOfDataArray("A");
		printDataDoublePair(optimizedData);
		printEndOfDataArray();
		
		printStartOfDataArray("B");
		printDataDoublePair(exponentialData);
		printEndOfDataArray();
		
		printStartOfDataArray("C");
		printDataDoublePair(heuristic1Data);
		printEndOfDataArray();
		
		printStartOfDataArray("D");
		printDataDoublePair(heuristic2Data);
		printEndOfDataArray();
		
		
		printStartOfDataArray("E");
		printDataDoublePair(heuristic3Data);
		printEndOfDataArray();
		
		printStartOfDataArray("F");
		printDataDoublePair(heuristic4Data);
		printEndOfDataArray();
		
		printEndOfMFile2a();
		BaseShutDown2();
		
	}

	private void printEndOfMFile2a() throws Exception 
	{
		m_bufWriter.write("Ax = A(:,1) \n");
		m_bufWriter.write("Ay = A(:,2) \n");
		m_bufWriter.write("Bx = B(:,1) \n");
		m_bufWriter.write("By = B(:,2) \n");		
		m_bufWriter.write("Cx = C(:,1) \n");
		m_bufWriter.write("Cy = C(:,2) \n");	
		m_bufWriter.write("Dx = D(:,1) \n");
		m_bufWriter.write("Dy = D(:,2) \n");	
		m_bufWriter.write("Ex = E(:,1) \n");
		m_bufWriter.write("Ey = E(:,2) \n");	
		m_bufWriter.write("Fx = F(:,1) \n");
		m_bufWriter.write("Fy = F(:,2) \n");	
		
		m_bufWriter.write("plot(Ax,Ay,'rd',Bx, By,'b+',Cx,Cy,'gx', Dx, Dy, 'ro' , Ex, Ey, 'ro', Fx, Fy, 'ro'); \n");
		m_bufWriter.write("legend('Optimized','Heuristic1','Heuristic2','Heuristic3'); \n");
		
		
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

	public void plotDouble2(Vector<PlotPairDoubleData> optimizedData,
			Vector<PlotPairDoubleData> heuristic1Data,
			Vector<PlotPairDoubleData> heuristic2Data,
			Vector<PlotPairDoubleData> heuristic3Data) throws Exception
	{
		BaseSetupForPlot2();
		printStartOfMFile2();
		
		printStartOfDataArray("A");
		printDataDoublePair(optimizedData);
		printEndOfDataArray();
		
		printStartOfDataArray("B");
		printDataDoublePair(heuristic1Data);
		printEndOfDataArray();
		
		printStartOfDataArray("C");
		printDataDoublePair(heuristic2Data);
		printEndOfDataArray();
		
		
		printStartOfDataArray("D");
		printDataDoublePair(heuristic3Data);
		printEndOfDataArray();
		
		printEndOfMFile2();
		BaseShutDown2();
		
	}
}
