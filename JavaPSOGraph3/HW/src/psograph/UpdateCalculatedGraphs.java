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

import psograph.util.GraphUpdater;


public class UpdateCalculatedGraphs 
{

	private GraphUpdater m_GraphUpdater;
	
    public UpdateCalculatedGraphs() throws Exception
    {
    	
    	m_GraphUpdater = new GraphUpdater();


    }
    
    public void update() throws Exception
    {
    	
    	File mainDirectory= new File(psograph.util.Util.baseDirectory);
		File seedDirectories[] = mainDirectory.listFiles();
		
		if(seedDirectories == null)
			throw new Exception ("No seed directories !!!!");
		
		int i;
		int maxSeedDirectories = seedDirectories.length;
		

		
		//For loop to go over each directory
		for( i=0; i <  maxSeedDirectories;
		i++)
		{
			m_GraphUpdater.processSeedDirectory(seedDirectories[i]);		
		}
		

    }
	/**
	 * @param args
	 */


    public static void main(String[] args) 
    {
    	try
    	{
    		
    		UpdateCalculatedGraphs graphUpdater = new UpdateCalculatedGraphs();
    		graphUpdater.update();

    	}
    	catch (Exception e)
    	{
    		e.printStackTrace();
    	}
    }

}
