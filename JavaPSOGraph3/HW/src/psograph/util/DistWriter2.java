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

public class DistWriter2 {
	MFilePlotData plotData;
	BufferedWriter m_bufWriter;
	String m_baseFileNames[];
	Vector<GraphLabel> m_labels;
	
	TreeMap<Double, Integer> optEdges;
	TreeMap<Double, Integer> expEdges;
	TreeMap<Double, Integer> heuEdges;
	
	double minValue;
	boolean maxMinValueSet = false;
	double maxValue;
	
	boolean isIntPlot = false;

	
	int numFiles =6;
	
	public DistWriter2(MFilePlotData PlotData)
	{
		plotData = PlotData;
		minValue =0;
	    maxValue =0;
	    
	    optEdges = new TreeMap<Double, Integer>();
		expEdges = new TreeMap<Double, Integer>();
		heuEdges = new TreeMap<Double, Integer>();
	}

	
	
	void printStartOfMFile() throws Exception
	{
    	
	}
	
	void printStartOfDataArray(String name) throws Exception
	{
		m_bufWriter.write(name+" = [ \n");
		m_bufWriter.flush();
	}
	
	void printEndOfDataArray() throws Exception
	{
		m_bufWriter.write("]; \n");
		
		m_bufWriter.flush();
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

//	public void plotIntegerInteger(TreeMap<Integer, Integer> data) throws Exception
//	{
//		BaseSetupForPlot();
//		printStartOfMFile();
//		printDataIntegerInteger(data);
//		printEndOfMFile();
//		BaseShutDown();
//		
//	}

//	private void printDataIntegerInteger(TreeMap<Integer, Integer> data) throws Exception
//	{
//		Vector<Integer> vKeys = new Vector<Integer>(data.keySet());
//		for(int i =0; i < vKeys.size(); i++)
//		{
//			for(int j =0 ; j < data.get(vKeys.get(i));j++)
//			{
//				m_bufWriter.write(vKeys.get(i) +";");
//			}
//		}
//		
//		numBins = vKeys.lastElement() + 1;
//		
//	}




//	public void plotDoubleInteger(TreeMap<Double, Integer> data) throws Exception
//	{
//		
//		BaseSetupForPlot();
//		printStartOfMFile();
//		printDataDoubleInteger(data);
//		printEndOfMFile();
//		BaseShutDown();
//		
//	}






//	private void printDataDoubleInteger(TreeMap<Double, Integer> data) throws Exception{
//		
//		Vector<Double> vKeys = new Vector<Double>(data.keySet());
//		for(int i =0; i < vKeys.size(); i++)
//		{   
//			for(int j =0 ; j < data.get(vKeys.get(i));j++)
//			{
//				m_bufWriter.write(vKeys.get(i) +";");
//			}
//		}
//		numBins = vKeys.size() + 1;
//		
//	}

	
	private void dealOptEdges(TreeMap<Double, Integer> data)throws Exception
	{

		Vector<Double> vKeys = new Vector<Double>(data.keySet());

		for(int i =0; i < vKeys.size(); i++)
		{   
			optEdges.put(vKeys.get(i), 1);
		}
	
	}
	

	private void printDataDoubleInteger(TreeMap<Double, Integer> data, TreeMap<Double, Integer> Edges) throws Exception
	{

		Vector<Double> vKeys = new Vector<Double>(data.keySet());
		
		Vector<Double> vEdges = new Vector<Double>(Edges.keySet());
		
		for(int i =0; i < vEdges.size(); i++)
		{
			if(data.containsKey(vEdges.get(i)))
			{
				m_bufWriter.write(vEdges.get(i)  + " "+  data.get(vEdges.get(i))+" ;");
			}
			else
			{
				m_bufWriter.write(vEdges.get(i)  + " "+  0 +" ;");
			}
			
		}
		
		
		
		
		//numBins = vKeys.size() + 1;
		
	}
	
	private void printDataDoubleIntegerOld(TreeMap<Double, Integer> data) throws Exception
	{
		
		
		
		Vector<Double> vKeys = new Vector<Double>(data.keySet());
		
		
		
		
		for(int i =0; i < vKeys.size(); i++)
		{   
			if(maxMinValueSet == false)
			{
				maxValue = vKeys.get(i);
				minValue = vKeys.get(i);
				maxMinValueSet =true;
			}
			else if(Double.compare(vKeys.get(i), maxValue) > 0)
			{
				maxValue = vKeys.get(i);
			}
			else if(Double.compare(vKeys.get(i), minValue) < 0)
			{
				minValue = vKeys.get(i);
			}
						
//			for(int j =0 ; j < data.get(vKeys.get(i));j++)
//			{
				m_bufWriter.write(vKeys.get(i)  + " "+  data.get(vKeys.get(i))+" ;");
//			}
		}
		
		

		
		
		
		
		//numBins = vKeys.size() + 1;
		
	}
	
//	private void printDataIntegerInteger(TreeMap<Integer, Integer> data) throws Exception
//	{
//		Vector<Integer> vKeys = new Vector<Integer>(data.keySet());
//		for(int i =0; i < vKeys.size(); i++)
//		{
//			for(int j =0 ; j < data.get(vKeys.get(i));j++)
//			{
//				m_bufWriter.write(vKeys.get(i) +";");
//			}
//		}
//		
//		numBins = vKeys.lastElement() + 1;
//		
//	}

	public void PlotDist(Vector<TreeMap<Double, Integer>> optimizedData,
			Vector<TreeMap<Double, Integer>> exponentialData,
			Vector<TreeMap<Double, Integer>> heuristicData) throws Exception
	{
		
		
		

		
		BaseSetupForPlot();
		printStartOfMFile();
		
		
		for(int i =0; i < optimizedData.size() ; i++)
		{
			dealOptEdges(optimizedData.get(i));
		
		}
		printOptEdges();
		
		
		for(int i =0; i < optimizedData.size() ; i++)
		{
			printStartOfDataArray("opt"+i);
			printDataDoubleInteger(optimizedData.get(i), optEdges);
			printEndOfDataArray();
			
		}
		
		for(int i =0; i < exponentialData.size() ; i++)
		{
			dealExpEdges(exponentialData.get(i));
		
		}
		printExpEdges();
		
		for(int i =0; i < exponentialData.size() ; i++)
		{
			
			printStartOfDataArray("exp"+i);
			printDataDoubleInteger(exponentialData.get(i), expEdges);
			printEndOfDataArray();
			
		}
		
		
		for(int i =0; i < heuristicData.size() ; i++)
		{
			dealHeuEdges(heuristicData.get(i));
		
		}
		printHeuEdges();
		
		for(int i =0; i < heuristicData.size() ; i++)
		{
			
			printStartOfDataArray("heu"+i);
			printDataDoubleInteger(heuristicData.get(i), heuEdges);
			printEndOfDataArray();
			
		}
		
		
		
		
		
		printEndOfMFile(optimizedData.size(), exponentialData.size(), heuristicData.size() );
		BaseShutDown();
		
	}

	private void dealExpEdges(TreeMap<Double, Integer> treeMap) {
		Vector<Double> vKeys = new Vector<Double>(treeMap.keySet());

		for(int i =0; i < vKeys.size(); i++)
		{   
			expEdges.put(vKeys.get(i), 1);
		}
		
	}



	private void printExpEdges() throws Exception {
		Vector<Double> vEdges = new Vector<Double>(expEdges.keySet());
		
		m_bufWriter.write("exp_edges = [ \n");
		
		
		for(int i =0; i < vEdges.size(); i++)
		{

			m_bufWriter.write(vEdges.get(i)  + "  ;");

			
		}
		m_bufWriter.write("];\n");
		
	}



	private void printHeuEdges() throws Exception {
		Vector<Double> vEdges = new Vector<Double>(heuEdges.keySet());
		
		m_bufWriter.write("heu_edges = [ \n");
		
		
		for(int i =0; i < vEdges.size(); i++)
		{

			m_bufWriter.write(vEdges.get(i)  + "  ;");

			
		}
		m_bufWriter.write("];\n");
		
	}



	private void dealHeuEdges(TreeMap<Double, Integer> treeMap) {
		Vector<Double> vKeys = new Vector<Double>(treeMap.keySet());

		for(int i =0; i < vKeys.size(); i++)
		{   
			heuEdges.put(vKeys.get(i), 1);
		}
		
	}



	private void printOptEdges() throws Exception
	{
		
		Vector<Double> vEdges = new Vector<Double>(optEdges.keySet());
		
		m_bufWriter.write("opt_edges = [ \n");
		
		
		for(int i =0; i < vEdges.size(); i++)
		{

			m_bufWriter.write(vEdges.get(i)  + "  ;");

			
		}
		m_bufWriter.write("];\n");
}



	private void printEndOfMFile(int optimizedSize, int exponentialSize, int heuristicSize) throws Exception
	{
		double step = (maxValue - minValue)/plotData.nBins;
		int i=0;
		
		m_bufWriter.write("edges = [ \n");
		if(isIntPlot)
		{
			for(i = 0; i < plotData.maxBin; i++)
			{
				double v = plotData.minBin + i;
				m_bufWriter.write(v+ "; \n");
				
			}
		}
		else
		{
			for(i = 0; i < plotData.nBins; i++)
			{
				double v = minValue + (step*i);
				m_bufWriter.write(v+ "; \n");
				
			}
		}
		m_bufWriter.write(plotData.maxBin+" ]; \n");
		
		m_bufWriter.write("edges_combined_opt = [ \n");
		for(i =0; i < optimizedSize ;i++)
		{
			m_bufWriter.write(" edges ");
		}
		m_bufWriter.write(" ]; \n");
		
				
		
		for(i =0; i < optimizedSize; i++)
		{
			m_bufWriter.write("opt"+i+"_x = opt"+i+"(:,2); \n");

		}
		
		m_bufWriter.write("opt_avg = [");
		for(i =0; i < optimizedSize; i++)
		{
			m_bufWriter.write("opt"+i+"_x ");

		}
		m_bufWriter.write("]; \n");
		m_bufWriter.write("opt_mean = mean(opt_avg,2);");
		
		for(i =0; i < exponentialSize; i++)
		{
			m_bufWriter.write("exp"+i+"_x = exp"+i+"(:,2); \n");
			
		}
		
		m_bufWriter.write("exp_avg = [");
		for(i =0; i < exponentialSize; i++)
		{
			m_bufWriter.write("exp"+i+"_x ");

		}
		m_bufWriter.write("]; \n");
		m_bufWriter.write("exp_mean = mean(exp_avg,2);");
		
		for(i =0; i < heuristicSize; i++)
		{
			m_bufWriter.write("heu"+i+"_x = heu"+i+"(:,2); \n");
			
		}
		
		m_bufWriter.write("heu_avg = [");
		for(i =0; i < heuristicSize; i++)
		{
			m_bufWriter.write("heu"+i+"_x ");

		}
		m_bufWriter.write("]; \n");
		m_bufWriter.write("heu_mean = mean(heu_avg,2);");
		
				
		
		
		m_bufWriter.write("plot(edges,opt_mean,'r', edges,exp_mean,'b',edges,heu_mean,'g' ); \n");
		
		
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

	

	public void PlotDistInt(Vector<TreeMap<Integer, Integer>> optimizedData,
			Vector<TreeMap<Integer, Integer>> exponentialData,
			Vector<TreeMap<Integer, Integer>> heuristicData) throws Exception
	{

		
		BaseSetupForPlot();
		printStartOfMFile();
		
		for(int i =0; i < optimizedData.size() ; i++)
		{
			printStartOfDataArray("opt"+i);
			printDataIntegerInteger(optimizedData.get(i));
			printEndOfDataArray();
			
		}
		
		for(int i =0; i < exponentialData.size() ; i++)
		{
			
			printStartOfDataArray("exp"+i);
			printDataIntegerInteger(exponentialData.get(i));
			printEndOfDataArray();
			
		}
		
		for(int i =0; i < heuristicData.size() ; i++)
		{
			
			printStartOfDataArray("heu"+i);
			printDataIntegerInteger(heuristicData.get(i));
			printEndOfDataArray();
			
		}
		
		printEndOfMFile(optimizedData.size(), exponentialData.size(), heuristicData.size() );
		BaseShutDown();
		
	}

	

	



	private void printDataIntegerInteger(TreeMap<Integer, Integer> data) throws Exception
	{
		isIntPlot = true;
		
		Vector<Integer> vKeys = new Vector<Integer>(data.keySet());
		for(int i =0; i <= plotData.maxBin; i++)
		{   
			if(data.containsKey(i))
			{
				m_bufWriter.write(i +" "+data.get(i) +";");
			}
			else
			{
				m_bufWriter.write(i+" 0;");
			}
		}
		
	}



	



	
}
