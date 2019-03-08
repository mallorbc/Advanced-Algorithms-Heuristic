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

package psograph.test;

import static org.junit.Assert.fail;

import java.io.File;
import java.util.List;

import difflib.Delta;
import difflib.DiffUtils;
import difflib.Patch;

/**
 * This is Base TestCase utility, to 
 * allow the testcase to compare known output(master)
 * against new data.  The idea is 
 * we can output detailed information and compare it
 * to master, as opposed to writing a lot of code to manually
 * go through the solution (in the case of SPL solution analysis)
 *
 * This solution depends on https://code.google.com/p/java-diff-utils/wiki/SampleUsage
 * Which means the tests now have a dependency on this as well now.   Shout out to 
 * this project for allowing me to focus on adding more functionality instead of writing
 * a system to do these differences myself.
 */
public class CompareTest 
{

	public File getMasterFile(String masterFile)
	{
		ClassLoader classLoader = getClass().getClassLoader();
		File file = new File(classLoader.getResource(masterFile).getFile());	
		
		return file;
	}
	
	public void compareData(String masterFileLocation, StringBuilder sb) throws Exception
	{
		File masterFile = getMasterFile(masterFileLocation);

        List<String> originalData = TestUtils.fileToLines(masterFile);
        List<String> newData  = TestUtils.stringToLines(sb.toString());
                
        // Compute diff. Get the Patch object. Patch is the container for computed deltas.
        Patch patch = DiffUtils.diff(originalData, newData);
        
        List<Delta> listDeltas = patch.getDeltas();

        for (Delta delta: listDeltas) 
        {
        	System.out.println(delta);
        }
		if(listDeltas.size()  != 0)
		{
			fail("Master is not the same as new.");
		}	
		
	}
	
}
