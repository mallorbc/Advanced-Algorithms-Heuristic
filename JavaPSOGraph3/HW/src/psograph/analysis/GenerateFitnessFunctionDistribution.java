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
package psograph.analysis;

import java.io.File;
import java.io.FilenameFilter;
import java.util.TreeMap;
import java.util.Vector;

import psograph.graph.CalculatedGraph;




public class GenerateFitnessFunctionDistribution 
{

	
    public GenerateFitnessFunctionDistribution() throws Exception
    {
    	m_mainDirectory= new File(psograph.util.Util.baseDirectory);
		File m_seedDirectories[] = m_mainDirectory.listFiles();
		
		Vector<File> vFiles = new Vector<File>();
		
		for(int i =0 ; i < m_seedDirectories.length; i++)
		{
			File t[] = m_seedDirectories[i].listFiles();
			for(int j =0; j < t.length; j++)
			{
				vFiles.add(t[j]);
			}
			
		}
		
		m_graphDirectories = new File[vFiles.size()];
		for(int i = 0; i < vFiles.size(); i++)
		{
			m_graphDirectories[i]=vFiles.get(i);
		}

		  
    }
    
    public GenerateFitnessFunctionDistribution(File canditateDir) throws Exception
    {
    	m_graphDirectories = new File[1];
    	m_graphDirectories[0] = canditateDir;
    }
    
    public void generateDist() throws Exception
    {
		
		if(m_graphDirectories == null)
			throw new Exception ("No seed directories !!!!");
		
		int i;
		int maxSeedDirectories = m_graphDirectories.length;
		

		
		//For loop to go over each directory
		for( i=0; i <  maxSeedDirectories; i++)
		{
			if(m_graphDirectories[i].isDirectory())
			{
				//System.out.println("Output fitness distribution of "+ m_graphDirectories[i].getAbsolutePath());
				processSeedDirectory(m_graphDirectories[i]);
			}
		}
		

    }
	/**
	 * @param args
	 */
    // It is also possible to filter the list of returned files.
    // This example does not return any files that start with `.'.
    FilenameFilter filter = new FilenameFilter() {
        public boolean accept(File dir, String name) {
            return name.endsWith(".CalculatedGraph");
        }
    };
	
    private TreeMap<Double, Vector<String>> m_graphFitnessFunctionMap;
    
    public void processSeedDirectory(File graphDirectory) throws Exception
    {  	
		
		if(!graphDirectory.exists())
		{
			throw new Exception("Directory does not exist "+graphDirectory.getAbsolutePath());
		}
    	
    	
		File calculatedGraphs[] = graphDirectory.listFiles(filter);
		if (calculatedGraphs == null)
			throw new Exception ("SeedDirectory is empty");
		
		System.out.println("Generating Fitness distribution for - "+ graphDirectory.getAbsolutePath());

    	int j;

    	m_graphFitnessFunctionMap = new TreeMap<Double, Vector<String>>();
		
		int maxCalculatedGraphs = calculatedGraphs.length;
		//for loop to go over each .calculatedGraph file

		
		for(j=0; j < maxCalculatedGraphs; j++)
		{
			CalculatedGraph graph =  null;

			//System.out.println("CalculatedGraph was retrieved Graph " + calculatedGraphs[j].toString());
			graph = psograph.util.Util.streaminCalculatedGraph(calculatedGraphs[j]);
			
			addToMap(graph.getFitnessValue(), calculatedGraphs[j].getName());
		}
	
		
//		Vector<Double> keyset = new Vector<Double>( m_graphFitnessFunctionMap.keySet());	
//		for(int k=0; k < keyset.size(); k++)
//		{
//			System.out.println("fitness value - "+ keyset.get(k));
//
//		}
		
		psograph.util.Util.streamOutFitnessDistribution(graphDirectory, m_graphFitnessFunctionMap);
		

    }
    
    private void addToMap(double fitnessValue, String graphFileName) throws Exception
    {
    	//m_graphFitnessFunctionMap
    	//System.out.println("fitness value - "+ fitnessValue + ", graphName - " + graphFileName);
    	Vector<String> v_strings = m_graphFitnessFunctionMap.get(fitnessValue);
    	if(v_strings == null)
    	{
    		v_strings = new Vector<String>();	
    	}
    	v_strings.add(graphFileName);
    	m_graphFitnessFunctionMap.put(fitnessValue, v_strings);
    }
    
	File m_mainDirectory;
	File m_graphDirectories[];
    


}
