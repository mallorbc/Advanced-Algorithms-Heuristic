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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Random;
import java.util.TreeMap;
import java.util.Vector;

import psograph.graph.CalculatedGraph;
import psograph.graph.CostBasis;
import psograph.graph.Graph;
import psograph.graph.GraphConstants;
import psograph.graph.Node;
import psograph.graph.NonLinearCostFunction;
import psograph.pso.PSOConstants;


public class Util {

	public static String baseDirectory = new String("C:\\TestHeuristic3\\");

	//public static String baseDirectory = new String("C:\\ThesisBaseDir\\");
	public static int i =0;
	
	public static String readTextFile(String  file)throws Exception
	{
		String everything;
		
		try(BufferedReader br = new BufferedReader(new FileReader(file))) 
		{
		    StringBuilder sb = new StringBuilder();
		    String line = br.readLine();

		    while (line != null) {
		        sb.append(line);
		        sb.append(System.lineSeparator());
		        line = br.readLine();
		    }
		    everything = sb.toString();
		}
		return everything;
	}
	
	public static File CreateCalculatedGraphDirectory(File seedDirectory) throws Exception
	{
	    File GraphDirectory = new File(seedDirectory.getAbsolutePath()+"\\graphs");
	    boolean success = GraphDirectory.mkdirs();
	    if (!success) {
	        throw new Exception ("Cannot make sub directory to store seeds.");
	    }
	    
	    return GraphDirectory;
	
	}
	
	public static File CreateSeedDirectory() throws Exception
	{
		File SeedDirectory;
		
		
	    // Get today's date
	    Date date = new Date();
	    Format formatter = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss");
	    String directory = formatter.format(date);
	    // 2002.01.29.08.36.33
	
	    // Create a directory; all non-existent ancestor directories are
	    // automatically created
	    SeedDirectory = new File(baseDirectory+directory+"_"+i);
	    i++;
	    boolean success = SeedDirectory.mkdirs();
	    if (!success) {
	        throw new Exception ("Cannot make directory to store seed.");
	    }
	    

	    
	    return SeedDirectory;
		
	}

	public static void streamOutExponentialGraph(File SeedDirectory,Graph ExponentialCostBasis,int i) throws Exception
	{
		// Serialize to a file
		ObjectOutput out = new ObjectOutputStream(new FileOutputStream(SeedDirectory+"\\ExponentialCostBasis"+i+".Graph"));
		out.writeObject(ExponentialCostBasis);
		out.close();  
		
		//Util.streamOutExponentialGraphAsMFile(new FileWriter(SeedDirectory+"\\ExponentialCostBasis"+i+".m") ,ExponentialCostBasis);
	}
	
	public static void streamOutSeed(File SeedDirectory,Graph graphSeed) throws Exception
	{
		// Serialize to a file
		ObjectOutput out = new ObjectOutputStream(new FileOutputStream(SeedDirectory+"\\Seed.Graph"));
		out.writeObject(graphSeed);
		out.close();  
		
		//Util.streamOutSeedGraphAsMFile(new FileWriter(SeedDirectory+"\\Seed.m") ,graphSeed);
	}

	public static void streamOutGraph(File SeedDirectory,Graph graph, int id) throws Exception
	{
		// Serialize to a file
		Integer t_id = new Integer(id);
		ObjectOutput out = new ObjectOutputStream(new FileOutputStream(SeedDirectory+"\\graph"+t_id.toString()+".Graph"));
		out.writeObject(graph);
		out.close();      
	}
	
	public static void streamOutFitnessDistribution(File SeedDirectory, TreeMap<Double, Vector<String>> m_graphFitnessFunctionMap) throws Exception
	{
		// Serialize to a file
		ObjectOutput out = new ObjectOutputStream(new FileOutputStream(SeedDirectory+"\\FitnessDistribution.FitnessDistribution"));
		out.writeObject(m_graphFitnessFunctionMap);
		out.close();      
	}
	
	@SuppressWarnings("unchecked")
	public static TreeMap<Double, Vector<String>> streaminFitnessDistribution(File FitnessDistribution) throws Exception
	{
		TreeMap<Double, Vector<String>> fitDist = null;
		// Serialize to a file
		
		ObjectInput in = new ObjectInputStream(new FileInputStream(FitnessDistribution));
		fitDist = (TreeMap<Double, Vector<String>>)in.readObject();
		in.close();   
		 
		return fitDist;
		
		
	}
	

	public static void printMFile(Graph graph) throws Exception
	{
		
		//Now lets write out a M-Lab file
		Integer adjecency[][] = graph.getAdjecencyGraph();
		if(adjecency == null)
			throw new Exception("Adjecency graph is NULL");
		
		int n =   adjecency.length;
		int ii, jj;
		System.out.println("A = [");
		for( ii =0; ii< n; ii++)
		{
			for( jj = 0; jj <n; jj++)
			{
				System.out.print(adjecency[ii][jj] + " ");
			}
			System.out.println(";");
		}
		System.out.println("];");
		
		Vector<Node> nodes = new Vector<Node>(graph.getHeaderNodesMap().values());
		
		System.out.println("xy = [");
		for( ii =0; ii< nodes.size(); ii++)
		{
			System.out.println(nodes.get(ii).getX() +" "+ nodes.get(ii).getY()+" 0;");
		}
		System.out.println("];");
		
		System.out.println("gplot(A,xy);");
	}
	
	public static void printMFileSpecial(CalculatedGraph cg1, double limit, File file) throws Exception
	{
		//Now lets right out a M-Lab file
		double adjecency[][] = cg1.getAdjecencyGraphWeights2();
		if(adjecency == null)
			throw new Exception("Adjecency graph is NULL");
		
		FileWriter fileLoc = new FileWriter(file);
		BufferedWriter mout = new BufferedWriter(fileLoc);
		
		int n =   adjecency.length;
		int ii, jj;
		
		mout.write("A1 = [ \n");
		
		//System.out.println("A = [");
		for( ii =0; ii< n; ii++)
		{
			for( jj = 0; jj <n; jj++)
			{
				double value = NonLinearCostFunction.getDistance(adjecency[ii][jj]);
				
				if(Double.compare(value, limit ) >= 0 )
				{
					mout.write(1 + " ");
				}
				else
				{
					mout.write(0 + " ");
				}
				//System.out.print(value + " ");
			}
			mout.write("; \n");
			//System.out.println(";");
		}
		mout.write("]; \n");
		//System.out.println("];\n");
		
		Vector<Node> nodes = new Vector<Node>(cg1.getHeaderNodesMap().values());
		
		mout.write("xy1 = [ \n");
		//System.out.println("xy = [");
		for( ii =0; ii< nodes.size(); ii++)
		{
			mout.write(nodes.get(ii).getX() +" "+ nodes.get(ii).getY()+"; \n");
			//System.out.println(nodes.get(ii).m_x +" "+ nodes.get(ii).m_y+";");
		}
		mout.write("]; \n");
		//System.out.println("];");
		
		mout.write("gplot(A1,xy1,'rd'); \n");
		
		mout.write("hold on \n;");
		
		mout.write("A2 = [ \n");
		
		//System.out.println("A = [");
		for( ii =0; ii< n; ii++)
		{
			for( jj = 0; jj <n; jj++)
			{
				double value = NonLinearCostFunction.getDistance(adjecency[ii][jj]);
				
				if(adjecency[ii][jj] == 0)
				{
					mout.write(0 + " ");
				}
				else if(Double.compare(value, limit ) < 0 )
				{
					mout.write(1 + " ");
				}
				else
				{
					mout.write(0 + " ");
				}
				//System.out.print(adjecency[ii][jj] + " ");
			}
			mout.write("; \n");
			//System.out.println(";");
		}
		mout.write("]; \n");
		//System.out.println("];");
		
		Vector<Node> nodes2 = new Vector<Node>(cg1.getHeaderNodesMap().values());
		
		mout.write("xy2 = [ \n");
		//System.out.println("xy = [");
		for( ii =0; ii< nodes2.size(); ii++)
		{
			mout.write(nodes2.get(ii).getX() +" "+ nodes2.get(ii).getY()+"; \n");
			//System.out.println(nodes.get(ii).m_x +" "+ nodes.get(ii).m_y+";");
		}
		mout.write("]; \n");
		
		mout.write("gplot(A2,xy2,':bd'); \n");

		mout.close();
		
	}

	public static void streamOutSeedGraphAsMFile(FileWriter fileLoc, Graph graph) throws Exception
	{
		
		BufferedWriter mout = new BufferedWriter(fileLoc);
	    
		//Now lets right out a M-Lab file
		
		int ii=0;
		Vector<Node> nodes = new Vector<Node>(graph.getHeaderNodesMap().values());
		
		mout.write("x = [ \n");
		//System.out.println("xy = [");
		for( ii =0; ii< nodes.size(); ii++)
		{
			mout.write(nodes.get(ii).getX() +"; \n");
			//System.out.println(nodes.get(ii).m_x +" "+ nodes.get(ii).m_y+";");
		}
		mout.write("]; \n");
		
		nodes = new Vector<Node>(graph.getHeaderNodesMap().values());
		
		mout.write("y = [ \n");
		//System.out.println("xy = [");
		for( ii =0; ii< nodes.size(); ii++)
		{
			mout.write(nodes.get(ii).getY()+"; \n");
			//System.out.println(nodes.get(ii).m_x +" "+ nodes.get(ii).m_y+";");
		}
		
		mout.write("]; \n");
		//System.out.println("];");
		mout.write("subplot(2,2,1); \n");
		mout.write("plot(x, y,'d'); \n");
		//System.out.println("gplot(A,xy);");
		mout.close();
		
	}
	
	public static void streamOutExponentialGraphAsMFile(FileWriter fileLoc, Graph graph) throws Exception
{
		
		BufferedWriter mout = new BufferedWriter(fileLoc);
	    
		//Now lets right out a M-Lab file
		Integer adjecency[][] = graph.getAdjecencyGraph();
		if(adjecency == null)
			throw new Exception("Adjecency graph is NULL");
		
		int n =   adjecency.length;
		int ii, jj;
		
		mout.write("A = [ \n");
		
		//System.out.println("A = [");
		for( ii =0; ii< n; ii++)
		{
			for( jj = 0; jj <n; jj++)
			{
				mout.write(adjecency[ii][jj] + " ");
				//System.out.print(adjecency[ii][jj] + " ");
			}
			mout.write("; \n");
			//System.out.println(";");
		}
		mout.write("]; \n");
		//System.out.println("];");
		
		Vector<Node> nodes = new Vector<Node>(graph.getHeaderNodesMap().values());
		
		mout.write("xy = [ \n");
		//System.out.println("xy = [");
		for( ii =0; ii< nodes.size(); ii++)
		{
			mout.write(nodes.get(ii).getX() +" "+ nodes.get(ii).getY()+"; \n");
			//System.out.println(nodes.get(ii).m_x +" "+ nodes.get(ii).m_y+";");
		}
		mout.write("]; \n");
		//System.out.println("];");
		mout.write("subplot(2,2,2); \n");
		mout.write("gplot(A,xy); \n");
		//System.out.println("gplot(A,xy);");
		mout.close();
		
	}

	public static void streamOutGraphAsMFile(FileWriter fileLoc, Graph graph, int iterID) throws Exception
	{
		
		BufferedWriter mout = new BufferedWriter(fileLoc);
	    
		//Now lets right out a M-Lab file
		Integer adjecency[][] = graph.getAdjecencyGraph();
		if(adjecency == null)
			throw new Exception("Adjecency graph is NULL");
		
		int n =   adjecency.length;
		int ii, jj;
		
		mout.write("A = [ \n");
		
		//System.out.println("A = [");
		for( ii =0; ii< n; ii++)
		{
			for( jj = 0; jj <n; jj++)
			{
				mout.write(adjecency[ii][jj] + " ");
				//System.out.print(adjecency[ii][jj] + " ");
			}
			mout.write("; \n");
			//System.out.println(";");
		}
		mout.write("]; \n");
		//System.out.println("];");
		
		Vector<Node> nodes = new Vector<Node>(graph.getHeaderNodesMap().values());
		
		mout.write("xy = [ \n");
		//System.out.println("xy = [");
		for( ii =0; ii< nodes.size(); ii++)
		{
			mout.write(nodes.get(ii).getX() +" "+ nodes.get(ii).getY()+"; \n");
			//System.out.println(nodes.get(ii).m_x +" "+ nodes.get(ii).m_y+";");
		}
		mout.write("]; \n");
		//System.out.println("];");
		
		if(iterID==0 || iterID==1)
		{
			int temp = iterID+3;
			mout.write("subplot(2,2,"+temp+"); \n");
		}
		
		mout.write("gplot(A,xy); \n");
		//System.out.println("gplot(A,xy);");
		mout.close();
		
	}

	public static double calculateCostBasisRation(CalculatedGraph graph)
	{
		double costBasisRatio = graph.getCost()/graph.getCostBasis();
		return costBasisRatio;
	}
	
	public static void printCalculatedGraphVitals(CalculatedGraph graph)
	{
		System.out.println("graph.getN " + graph.getNumberOfNodes());
		System.out.println("getTotalEdges "+graph.getTotalEdges());
		System.out.println("getAvgCostBasis " + graph.getCostBasis());
		System.out.println("getCost " + graph.getCost());
		
		
		System.out.println("getCost/getAvgCostBasis "+ Util.calculateCostBasisRation(graph) );
		System.out.println("1 - getCost/getAvgCostBasis "+ (1 - Util.calculateCostBasisRation(graph))); 
		System.out.println("directLCC "+ graph.getDirectLCC());
		System.out.println("randomLCC "+ graph.getRandomLCC());
		System.out.println("AISPL " + graph.getAISPL());
		
		System.out.println("Graph1(cost score, AISPL, directLCC) "+"("+(1 - Util.calculateCostBasisRation(graph)) +", "
				+graph.getAISPL() + ", "+graph.getRandomLCC()+")");
		System.out.println("Graph2(cost score, AISPL, directLCC) "+"("+(1 - Util.calculateCostBasisRation(graph)) +", "
				+graph.getAISPL() + ", "+graph.getDirectLCC()+")");
	}
	

	
	public static void streamOutCalculatedGraph(File SeedDirectory, int i, CalculatedGraph graph) throws Exception
	{
		// Serialize to a file
		ObjectOutput out = new ObjectOutputStream(new FileOutputStream(SeedDirectory+"\\Graph"+i+".CalculatedGraph"));
		out.writeObject(graph);
		out.close();      
		
		streamOutGraphAsMFile(new FileWriter(SeedDirectory+"\\Graph"+i+".m") ,graph, i);		
	
	}
	
	public static void streamOutCalculatedGraphWithExtension( File SeedDirectory, String fileName, CalculatedGraph graph) throws Exception
	{
		// Serialize to a file
		ObjectOutput out = new ObjectOutputStream(new FileOutputStream(SeedDirectory+"\\"+fileName));
		out.writeObject(graph);
		out.close();      
		
		//streamOutGraphAsMFile(new FileWriter(SeedDirectory+"\\Graph"+i+".m") ,graph, i);		
	
	}
	public static void streamOutCalculatedGraph(File SeedDirectory, String fileName, CalculatedGraph graph) throws Exception
	{
		// Serialize to a file
		ObjectOutput out = new ObjectOutputStream(new FileOutputStream(SeedDirectory+"\\"+fileName+".CalculatedGraph"));
		out.writeObject(graph);
		out.close();      
		
		streamOutGraphAsMFile(new FileWriter(SeedDirectory+"\\"+fileName+".m") ,graph, i);		
	
	}

	public static CalculatedGraph streaminCalculatedGraph(String SeedDirectoryStr, int id) throws Exception
	{

		return streaminCalculatedGraph(new File(SeedDirectoryStr) , id);
			
	}		
	

	public static CalculatedGraph streaminCalculatedGraph(File calcGraphFile) throws Exception
	{
		CalculatedGraph graph = null;
		// Serialize to a file
		
		ObjectInput in = new ObjectInputStream(new FileInputStream(calcGraphFile));
		graph = (CalculatedGraph)in.readObject();
		in.close();   
		
		return graph;
		
		
	}
	
	public static CalculatedGraph streaminCalculatedGraph(File SeedDirectory, int id) throws Exception
	{
		CalculatedGraph graph = null;
		// Serialize to a file
		Integer t_id = new Integer(id);
		String tt = SeedDirectory+"\\graph"+t_id.toString()+".CalculatedGraph";
		File f = new File(tt);
		@SuppressWarnings("unused")
		boolean b = f.isFile();
		
		ObjectInput in = new ObjectInputStream(new FileInputStream(SeedDirectory+"\\graph"+t_id.toString()+".CalculatedGraph"));
		graph = (CalculatedGraph)in.readObject();
		in.close();   
		
		return graph;
		
		
	}
	
	public static Graph streaminGraph(String SeedDirectoryStr, int id) throws Exception
	{

		return streaminGraph(new File(SeedDirectoryStr) , id);
			
	}		
	
	public static Graph streaminGraph(File SeedDirectory,int id) throws Exception
	{
		
		Graph graph  =null;
		// Serialize to a file
		Integer t_id = new Integer(id);
		ObjectInput in = new ObjectInputStream(new FileInputStream(SeedDirectory+"\\graph"+t_id.toString()+".graph"));
		graph = (Graph)in.readObject();
		in.close();     
		
		
		return graph;
		/*System.out.println("Proof of steaming input of graph");
		graph.printWithLocationAndWeights();
		System.out.println("End Proof of steaming input of graph");*/
	}

	public static Graph streaminSeed(File SeedDirectory) throws Exception
	{
		Graph graphSeed = null;
		
		// Serialize to a file
		ObjectInput in = new ObjectInputStream(new FileInputStream(SeedDirectory+"\\Seed.graph"));
		graphSeed = (Graph)in.readObject();
		in.close(); 
		
		/*System.out.println("Proof of steaming input of seed");
		graphSeed.printWithLocationAndWeights();
		System.out.println("End Proof of steaming input of seed");*/
		return graphSeed;
		
	}
	
	public static double sumArray(double inputArray[])
	{
		double res = 0;
		
		for(int i = 0; i < inputArray.length ; i++)
		{
			res += inputArray[i];
		}
		
		return res;
	}

	public static CostBasis streaminCostBasis(File file) throws Exception
	{
		CostBasis graph  =null;
		// Serialize to a file
		ObjectInput in = new ObjectInputStream(new FileInputStream(file));
		graph = (CostBasis)in.readObject();
		in.close();     
		
		
		return graph;
	}
	
	public static Graph streaminGraph(File file) throws Exception
	{
		Graph graph  =null;
		// Serialize to a file
		ObjectInput in = new ObjectInputStream(new FileInputStream(file));
		graph = (Graph)in.readObject();
		in.close();     
		
		
		return graph;
	}
	
	public static void perturbParticle(Graph particle) throws Exception
	{
		Random r = new Random();
		
		int base = r.nextInt(1+ PSOConstants.PERTURB_MAXIMUM - PSOConstants.PERTURB_MINIMUM);
		
		int perturbAmount = PSOConstants.PERTURB_MINIMUM + base;
		
		for(int i =0; i < perturbAmount; i++)
		{
			int pickA, pickB;
			
			//Olekas May-19-2012 had this using the constant which is wrong
			//should be indicated by the number of Nodes in the graph
			pickA = r.nextInt(particle.getNumberOfNodes());
			pickB = r.nextInt(particle.getNumberOfNodes());
			while(pickA == pickB)
			{
				pickB = r.nextInt(GraphConstants.NUM_NODES);				
			}
			
			if(!particle.isNodeAConnectedToNodeB(particle.getNode(pickA), particle.getNode(pickB)))
			{
				particle.addConnection(pickA, pickB);
			}
			else
			{
				particle.removeConnection(pickA, pickB);
			}
		}
	}
	
	

}
