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
package psograph;

import java.io.File;
import java.io.FilenameFilter;
import java.util.TreeMap;
import java.util.Vector;

import psograph.util.BasicMeasurementsPlotWriter;



public class NbVsNgPlots {

	TreeMap<Double, Vector<String>> m_graphFitnessFunctionMap ;
		
		@SuppressWarnings("unused")
		private BasicMeasurementsPlotWriter m_basicMeasurementsPlotWriter;
		
	    public NbVsNgPlots() throws Exception
	    {
	    	
	    	m_basicMeasurementsPlotWriter = new BasicMeasurementsPlotWriter();
	    	m_numBad = 100;
	 	    m_numGood = 20;


	    }
	    // It is also possible to filter the list of returned files.
	    // This example does not return any files that start with `.'.
	    FilenameFilter filter = new FilenameFilter() {
	        public boolean accept(File dir, String name) {
	            return name.endsWith(".FitnessDistribution");
	        }
	    };
	    public void processSeedDirectory(File seedDirectory) throws Exception
	    {
			File fitnessDist[] = seedDirectory.listFiles(filter);
			if (fitnessDist == null || fitnessDist.length > 1)
				throw new Exception ("SeedDirectory is empty or has too many distrib files");
			
			System.out.println("Working on SeedDirectory - "+ seedDirectory.getAbsolutePath());

	    	m_graphFitnessFunctionMap = psograph.util.Util.streaminFitnessDistribution(fitnessDist[0] );
			
	    	Vector<Double> keySet = new Vector<Double>(m_graphFitnessFunctionMap.keySet());
	    	
	    	
	    	if(keySet.size() < m_numBad || keySet.size() < m_numGood )
	    	{
	    		throw new Exception("number of samples is too small");
	    	}
	    	
	    	

	    }
	    
	    
	    public void generatePlots() throws Exception
	    {
	    	
	    	File mainDirectory= new File(psograph.util.Util.baseDirectory);
			File seedDirectories[] = mainDirectory.listFiles();
			
			if(seedDirectories == null)
				throw new Exception ("No seed directories !!!!");
			
			int i;
			int maxSeedDirectories = seedDirectories.length;
			
			
			
			//For loop to go over each directory
			for( i=0; i <  maxSeedDirectories; i++)
			{
				
				processSeedDirectory(seedDirectories[i]);		
			}
			

	    }
		/**
		 * @param args
		 */


	    public static void main(String[] args) 
	    {
	    	try
	    	{
	    		
	    		NbVsNgPlots nbVsNgPlots = new NbVsNgPlots();
	    		nbVsNgPlots.generatePlots();

	    	}
	    	catch (Exception e)
	    	{
	    		e.printStackTrace();
	    	}
	    }

	    int m_numBad;
	    int m_numGood;


	}


