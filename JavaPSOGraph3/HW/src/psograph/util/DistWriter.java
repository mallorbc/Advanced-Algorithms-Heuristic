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

public class DistWriter 
{
	MFilePlotData plotData;
	BufferedWriter m_bufWriter;
	String m_baseFileNames[];
	Vector<GraphLabel> m_labels;
	
	double minValue;
	boolean maxMinValueSet = false;
	double maxValue;
	
	boolean isIntPlot = false;

	
	public DistWriter(MFilePlotData PlotData)
	{
		plotData = PlotData;
		minValue =0;
	    maxValue =0;
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

	
	private void printDataDoubleVector(Vector<Double> data) throws Exception{
		
		for(int i =0; i < data.size(); i++)
		{   
			if(maxMinValueSet == false)
			{
				maxValue = data.get(i);
				minValue = data.get(i);
				maxMinValueSet =true;
			}
			else if(Double.compare(data.get(i), maxValue) > 0)
			{
				maxValue = data.get(i);
			}
			else if(Double.compare(data.get(i), minValue) < 0)
			{
				minValue = data.get(i);
			}
			
		//	for(int j =0 ; j < data.size();j++)
		//	{
			m_bufWriter.write(data.get(i) +";");
		//	}
		}
		//numBins = vKeys.size() + 1;
		
	}
	

	private void printDataDoubleInteger(TreeMap<Double, Integer> data) throws Exception{
		
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
			
			for(int j =0 ; j < data.get(vKeys.get(i));j++)
			{
				m_bufWriter.write(vKeys.get(i) +";");
			}
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
			printStartOfDataArray("opt"+i);
			printDataDoubleInteger(optimizedData.get(i));
			printEndOfDataArray();
			
		}
		
		for(int i =0; i < exponentialData.size() ; i++)
		{
			
			printStartOfDataArray("exp"+i);
			printDataDoubleInteger(exponentialData.get(i));
			printEndOfDataArray();
			
		}
		
		for(int i =0; i < heuristicData.size() ; i++)
		{
			
			printStartOfDataArray("heu"+i);
			printDataDoubleInteger(heuristicData.get(i));
			printEndOfDataArray();
			
		}
		
		
		
		
		
		printEndOfMFile(optimizedData.size(), exponentialData.size(), heuristicData.size() );
		BaseShutDown();
		
	}

	private void printEndOfMFile(int optimizedSize, int exponentialSize, int heuristicSize) throws Exception
	{
		double step = (maxValue - minValue)/plotData.nBins;
		int i=0;
		
		m_bufWriter.write("edges = [ \n");
		if(isIntPlot)
		{
			for(i = 0; i < plotData.nBins; i++)
			{
				double v = minValue + i;
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
		m_bufWriter.write(maxValue+" ]; \n");
		
		m_bufWriter.write("edges_combined_opt = [ \n");
		for(i =0; i < optimizedSize ;i++)
		{
			m_bufWriter.write(" edges ");
		}
		m_bufWriter.write(" ]; \n");
		
		m_bufWriter.write("edges_combined_exp = [ \n");
		for(i =0; i < exponentialSize;i++)
		{
			m_bufWriter.write(" edges ");
		}
		m_bufWriter.write(" ]; \n");
		
		m_bufWriter.write("edges_combined_heu = [ \n");
		for(i =0; i <  heuristicSize;i++)
		{
			m_bufWriter.write(" edges ");
		}
		m_bufWriter.write(" ]; \n");
		
		m_bufWriter.write("edges_combined = [ \n");
		for(i =0; i < optimizedSize +exponentialSize+ heuristicSize;i++)
		{
			m_bufWriter.write(" edges ");
		}
		m_bufWriter.write(" ]; \n");
		
		
		for(i =0; i < optimizedSize; i++)
		{
			m_bufWriter.write("opt"+i+"_x = histc(opt"+i+",edges); \n");

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
			m_bufWriter.write("exp"+i+"_x = histc(exp"+i+",edges); \n");
			
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
			m_bufWriter.write("heu"+i+"_x = histc(heu"+i+",edges); \n");
			
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

	private void printEndOfMFile2(Vector<String> varName, Vector<Vector<TreeMap<Double, Integer>>> data ) throws Exception
	{
		Vector<String> lineSpec = new Vector<String>();
		lineSpec.add("ro");
		lineSpec.add("bs");
		lineSpec.add("rd");
		lineSpec.add("b^");
		lineSpec.add("rv");
		lineSpec.add("b>");

		
		double step = (maxValue - minValue)/plotData.nBins;
		int i=0;
		
		m_bufWriter.write("edges = [ \n");
		if(isIntPlot)
		{
			for(i = 0; i < plotData.nBins; i++)
			{
				double v = minValue + i;
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
		m_bufWriter.write(maxValue+" ]; \n");
		
		System.out.println("data size is "+ data.size());
		for( i =0; i < data.size() ;i++)
		{
			for(int j = 0; j < data.get(i).size(); j++)
			{
//				printStartOfDataArray(varName.get(i)+i);
//				printDataDoubleInteger(data.get(i).get(j));
//				printEndOfDataArray();
				
				m_bufWriter.write(varName.get(i)+j+"_x = histc("+varName.get(i)+j+",edges); \n");
				
			}
			
			m_bufWriter.write("opt"+varName.get(i)+"_avg = [");
			for(int j = 0; j < data.get(i).size(); j++)
			{
				m_bufWriter.write(varName.get(i)+j+"_x ");

			}
			m_bufWriter.write("]; \n");
			m_bufWriter.write("opt"+varName.get(i)+"_mean = mean(opt"+varName.get(i)+"_avg,2);");
			
		}
		
		m_bufWriter.write("plot(");
		for( i =0; i < data.size() ; i++)
		{
			if(i == data.size() -1)
			{
				m_bufWriter.write("edges, opt"+varName.get(i)+"_mean,'"+lineSpec.get(i)+"'");
			}
			else
			{
				m_bufWriter.write("edges, opt"+varName.get(i)+"_mean,'"+lineSpec.get(i)+"',");
			}
		}
		m_bufWriter.write("	); \n");
	
		
		
		
		m_bufWriter.write("legend('Seed0','Seed1','Seed2','Seed3','Seed4','Seed5'); \n");
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

	public void PlotDistInt(Vector<Vector<TreeMap<Integer, Integer>>> data) throws Exception
	{
		
		Vector<String> varName = new Vector<String>();
		
		varName.add("A");
		varName.add("B");
		varName.add("C");
		varName.add("D");
		varName.add("E");
		varName.add("F");
		

		
		BaseSetupForPlot();
		printStartOfMFile();
		
		for(int i =0; i < data.size() ; i++)
		{
			for(int j = 0; j < data.get(i).size(); j++)
			{
				printStartOfDataArray(varName.get(i)+j);
				printDataIntegerInteger(data.get(i).get(j));
				printEndOfDataArray();
			}
			
		}

		printEndOfMFile2Int(varName ,data);
		BaseShutDown();
		
	}

	private void printEndOfMFile2Int(Vector<String> varName,
			Vector<Vector<TreeMap<Integer, Integer>>> data) throws Exception{
		
		Vector<String> lineSpec = new Vector<String>();
		lineSpec.add("ro");
		lineSpec.add("bs");
		lineSpec.add("rd");
		lineSpec.add("b^");
		lineSpec.add("rv");
		lineSpec.add("b>");

		
		double step = (maxValue - minValue)/plotData.nBins;
		int i=0;
		
		m_bufWriter.write("edges = [ \n");
		if(isIntPlot)
		{
			for(i = 0; i < plotData.nBins; i++)
			{
				double v = minValue + i;
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
		m_bufWriter.write(maxValue+" ]; \n");
		
		System.out.println("data size is "+ data.size());
		for( i =0; i < data.size() ; i++)
		{
			for(int j = 0; j < data.get(i).size(); j++)
			{
//				printStartOfDataArray(varName.get(i)+i);
//				printDataDoubleInteger(data.get(i).get(j));
//				printEndOfDataArray();
				
				m_bufWriter.write(varName.get(i)+j+"_x = histc("+varName.get(i)+j+",edges); \n");
				
			}
			
			m_bufWriter.write("opt"+varName.get(i)+"_avg = [");
			for(int j = 0; j < data.get(i).size(); j++)
			{
				m_bufWriter.write(varName.get(i)+j+"_x ");

			}
			m_bufWriter.write("]; \n");
			m_bufWriter.write("opt"+varName.get(i)+"_mean = mean(opt"+varName.get(i)+"_avg,2);");
			
		}
		
		m_bufWriter.write("plot(");
		for( i =0; i < data.size() ; i++)
		{
			if(i == data.size() -1)
			{
				m_bufWriter.write("edges, opt"+varName.get(i)+"_mean,'"+lineSpec.get(i)+"'");
			}
			else
			{
				m_bufWriter.write("edges, opt"+varName.get(i)+"_mean,'"+lineSpec.get(i)+"',");
			}
		}
		m_bufWriter.write("	); \n");
	
		
		
		
		m_bufWriter.write("legend('Seed0','Seed1','Seed2','Seed3','Seed4','Seed5'); \n");
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



	private void printDataIntegerInteger(TreeMap<Integer, Integer> data) throws Exception
	{
		
		if(data == null)
		{
			return;
		}
		
		try
		{
		isIntPlot = true;
		
		Vector<Integer> vKeys = new Vector<Integer>(data.keySet());
		for(int i =0; i < vKeys.size(); i++)
		{   
			if(maxMinValueSet == false)
			{
				maxValue = vKeys.get(i);
				minValue = vKeys.get(i);
				maxMinValueSet =true;
			}
			else if(vKeys.get(i) > maxValue)
			{
				maxValue = vKeys.get(i);
			}
			else if(vKeys.get(i) < minValue) 
			{
				minValue = vKeys.get(i);
			}
			
			for(int j =0 ; j < data.get(vKeys.get(i));j++)
			{
				m_bufWriter.write(vKeys.get(i) +";");
			}
		}
		//numBins = vKeys.size() + 1;
		}
		catch (Exception e)
		{
			throw e;
			
		}
	}



	public void PlotDist(Vector<Vector<TreeMap<Double, Integer>>> data) throws Exception
	{
		
		Vector<String> varName = new Vector<String>();
		
		varName.add("A");
		varName.add("B");
//		varName.add("C");
//		varName.add("D");
//		varName.add("E");
//		varName.add("F");
		

		
		BaseSetupForPlot();
		printStartOfMFile();
		
		for(int i =0; i < data.size() ; i++)
		{
			for(int j = 0; j < data.get(i).size(); j++)
			{
				printStartOfDataArray(varName.get(i)+j);
				printDataDoubleInteger(data.get(i).get(j));
				printEndOfDataArray();
			}
			
		}
		
	
		
		printEndOfMFile2(varName ,data);
		BaseShutDown();
		
	}



	public void plotDoubleVectors(Vector<Double> optimizedData,
			Vector<Double> exponentialData, Vector<Double> heuristicData) throws Exception
	{
		BaseSetupForPlot();
		printStartOfMFile();
		
	//	for(int i =0; i < optimizedData.size() ; i++)
	//	{
			printStartOfDataArray("opt"+0);
			printDataDoubleVector(optimizedData);
			printEndOfDataArray();
			
	//	}
		
	//	for(int i =0; i < exponentialData.size() ; i++)
	//	{
			
			printStartOfDataArray("exp"+0);
			printDataDoubleVector(exponentialData);
			printEndOfDataArray();
			
	//	}
		
	//	for(int i =0; i < heuristicData.size() ; i++)
	//	{
			
			printStartOfDataArray("heu"+0);
			printDataDoubleVector(heuristicData);
			printEndOfDataArray();
			
	//	}

		
		printEndOfMFile(1,
						1,
						1 );
		
		BaseShutDown();
		
	}



	




	
}
