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


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;


import psograph.graph.Graph;


public class TestUtils 
{
	
    // Helper method for get the file content
    static List<String> fileToLines(File file) 
    {
            List<String> lines = new LinkedList<String>();
            String line = "";
            try {
                    BufferedReader in = new BufferedReader(new FileReader(file));
                    while ((line = in.readLine()) != null) {
                            lines.add(line);
                    }
                    in.close();
            } catch (IOException e) {
                    e.printStackTrace();
            }
            
            return lines;
    }	
    
    static List<String> stringToLines(String str) throws IOException
    {
        BufferedReader rdr = new BufferedReader(new StringReader(str));
        List<String> lines = new ArrayList<String>();
        for (String line = rdr.readLine(); line != null; line = rdr.readLine()) {
            lines.add(line);
        }
        rdr.close();
        
        return lines;
    }
    

	
    // Helper method for get the file content
    static List<String> fileToLines(String filename) 
    {
            List<String> lines = new LinkedList<String>();
            String line = "";
            try {
                    BufferedReader in = new BufferedReader(new FileReader(filename));
                    while ((line = in.readLine()) != null) {
                            lines.add(line);
                    }
                    in.close();
            } catch (IOException e) {
                    e.printStackTrace();
            }
            return lines;
    }

	static public Graph createTestGraph3()  throws Exception
	{
		Graph tt = new Graph(0);
		
		tt.addNodeLocationDataOnly(.1, .1);
		tt.addNodeLocationDataOnly(.34, .34);
		tt.addNodeLocationDataOnly(.1, 1);
		tt.addNodeLocationDataOnly(.23, .78);
		tt.addNodeLocationDataOnly(.4, .1);
		tt.addNodeLocationDataOnly(.1, .5);
		tt.addNodeLocationDataOnly(.5, .5);
		tt.addNodeLocationDataOnly(.91, .91);
		tt.addNodeLocationDataOnly(.8, .7);
		tt.addNodeLocationDataOnly(.34, .78);
		tt.addNodeLocationDataOnly(.12, .1);
		tt.addNodeLocationDataOnly(.34, .342);
		tt.addNodeLocationDataOnly(.32, .83);
		tt.addNodeLocationDataOnly(.4, .34);
		tt.addNodeLocationDataOnly(.4, .31);
		tt.addNodeLocationDataOnly(.1, .35);
		tt.addNodeLocationDataOnly(.15, .95);
		tt.addNodeLocationDataOnly(.11, .91);
		tt.addNodeLocationDataOnly(.81, .71);
		tt.addNodeLocationDataOnly(.134, .178);
		
		
		
		tt.addConnection(0,1);
		tt.addConnection(0,3);
	
		tt.addConnection(1,2);
		tt.addConnection(1,3);
		tt.addConnection(1,4);
	

		tt.addConnection(2,3);
			
		tt.addConnection(5,0);
		tt.addConnection(5,4);
		
		tt.addConnection(6,1);
		tt.addConnection(6,2);
		
		tt.addConnection(7,6);
		tt.addConnection(7,0);
		
	
		tt.addConnection(8,3);
		tt.addConnection(8,5);
		
		tt.addConnection(9,2);
		tt.addConnection(9,4);		
		tt.addConnection(9,12);
		tt.addConnection(9,14);
		tt.addConnection(9,13);
		tt.addConnection(9,16);
		
		
		tt.addConnection(10,14);
		tt.addConnection(10,1);
		tt.addConnection(10,3);
		
		tt.addConnection(11,14);
		tt.addConnection(11,4);
		tt.addConnection(11,5);
		tt.addConnection(11,17);
		
		tt.addConnection(12,14);
		tt.addConnection(12,1);
		tt.addConnection(12,13);
		tt.addConnection(12,7);


		tt.addConnection(13,14);
		tt.addConnection(13,10);
		tt.addConnection(13,19);
		
		tt.addConnection(14,1);
		tt.addConnection(14,8);
		tt.addConnection(14,3);
		
		tt.addConnection(15,14);
		
		tt.addConnection(16,14);
		
		tt.addConnection(17,14);
		
		tt.addConnection(18,14);
		
		tt.addConnection(19,14);
		

		
		return tt;
	}
	
	
	static public Graph createTestGraph()  throws Exception
	{
		Graph tt = new Graph(10);
		
		tt.addConnection(0,1,1);
		tt.addConnection(0,3, 2);
	
		//tt.addConnection(1,0,3);
		tt.addConnection(1,2,4);
		tt.addConnection(1,3,4);
		tt.addConnection(1,4,5);
	
	//	tt.addConnection(2,1,6);
		tt.addConnection(2,3,3);
	
		//tt.addConnection(3,0,3);
	//	tt.addConnection(3,1,3);
		//tt.addConnection(3,2,3);
	
		//tt.addConnection(4,1,3);
		
		tt.addConnection(5,0,2);
		tt.addConnection(5,4,4);
		
		tt.addConnection(6,1,2);
		tt.addConnection(6,2,2);
		
		tt.addConnection(7,6,2);
		tt.addConnection(7,0,6);
		
	
		tt.addConnection(8,3,5);
		tt.addConnection(8,5,4);
		
		tt.addConnection(9,2,4);
		tt.addConnection(9,4,22);
		
		return tt;
	}
	static public Graph createEmptyGraph() throws Exception
	{
		Graph tt = new Graph(0);
		
		return tt;
	}
	static public Graph createTestGraph2() throws Exception
	{
		Graph tt = createTestGraph();
		
		double weights[];
		
		int neighbors[];
		neighbors = new int[3];
		neighbors[0] = 0;
		neighbors[1] = 2;
		neighbors[2] = 4;
		weights = new double[3];
		weights[0] = weights[1] =weights[2] =1;
		tt.addNode(neighbors, weights);
	
	
		tt.addNode();
	
	
		tt.removeNode(0);
		tt.removeNode(5); 
		
		neighbors = new int[1];
		neighbors[0] = 11;		
		weights = new double[1];
		weights[0] =  1;
		tt.addNode(neighbors, weights);
		tt.addNode();
	
		tt.removeConnection(2,9);
	
		return tt;
	}

	static public Graph createSimpleLineGraph() throws Exception
	{
		Graph tt = new Graph(5);
		
		tt.addConnection(0,1,1);
		tt.addConnection(1, 2,1);
	
		tt.addConnection(2,3,1);
		tt.addConnection(3,4,1);

		
		return tt;
	}
	
	static public Graph createSimpleRingGraph() throws Exception
	{
		Graph tt = new Graph(5);
		
		tt.addConnection(0,1,1);
		tt.addConnection(1, 2,1);
	
		tt.addConnection(2,3,1);
		tt.addConnection(3,4,1);
		tt.addConnection(4,0,1);
		
		return tt;
	}


	public static Graph createTestGraph4() throws Exception
	{
		
		Graph tt = new Graph(0);
		
		tt.addNodeLocationDataOnly(.2, .2);
		tt.addNodeLocationDataOnly(.2, .8);
		tt.addNodeLocationDataOnly(.8, .8);
		tt.addNodeLocationDataOnly(.8, .2);
		tt.addNodeLocationDataOnly(.5, .5);
		
		
		tt.addConnection(0, 1);
		tt.addConnection(0, 3);
		tt.addConnection(0, 4);
		
		//tt.addConnection(1, 0);
		tt.addConnection(1, 2);
		tt.addConnection(1, 4);
		
		//tt.addConnection(2, 1);
		tt.addConnection(2, 3);
		tt.addConnection(2, 4);
		
		//tt.addConnection(3, 0);
		//tt.addConnection(3, 2);
		tt.addConnection(3, 4);
		
		return tt;
		
	}


	public static Graph createTestGraph6() throws Exception
	{
		Graph tt = new Graph(0);
		
		tt.addNodeLocationDataOnly(.2, .2);
		tt.addNodeLocationDataOnly(.2, .3);
		tt.addNodeLocationDataOnly(.3, .2);
		tt.addNodeLocationDataOnly(.8, .8);
		tt.addNodeLocationDataOnly(.8, .9);
		
		
		tt.addConnection(0, 1);
		tt.addConnection(0, 2);
		tt.addConnection(0, 3);
		
		//tt.addConnection(3, 0);
		//tt.addConnection(3, 2);
		tt.addConnection(3, 4);
		
		return tt;
	}


	public static Graph createTestGraph5() throws Exception{
		Graph tt = new Graph(0);
		
		tt.addNodeLocationDataOnly(.2, .2);
		tt.addNodeLocationDataOnly(.2, .8);
		tt.addNodeLocationDataOnly(.8, .2);
		tt.addNodeLocationDataOnly(.8, .8);

		
		
		tt.addConnection(0, 1);
		tt.addConnection(0, 2);
		tt.addConnection(1, 3);
		
		//tt.addConnection(3, 0);
		//tt.addConnection(3, 2);
		tt.addConnection(2, 3);
		
		return tt;
	}


	public static Graph createSimpleRingGraphButBroken() throws Exception {
		Graph tt = new Graph(5);
		
		tt.addConnection(0,1,1);
		tt.addConnection(1, 2,1);
	
		tt.addConnection(2,3,1);


		
		return tt;
	}


	public static Graph createTestGraph7() throws Exception {

			Graph tt = new Graph(0);
			
			tt.addNodeLocationDataOnly(.2, .2);
			tt.addNodeLocationDataOnly(.2, .8);
			tt.addNodeLocationDataOnly(.8, .2);
			tt.addNodeLocationDataOnly(.8, .8);
			tt.addNodeLocationDataOnly(.5, .8);
			
			
			tt.addConnection(0, 1);
			tt.addConnection(0, 2);
			tt.addConnection(1, 3);
			
			//tt.addConnection(3, 0);
			//tt.addConnection(3, 2);
			tt.addConnection(2, 3);
			
			tt.addConnection(3,4);
			
			return tt;
		}


	public static Graph createTestGraph13() throws Exception {
		Graph tt = new Graph(0);
		
		tt.addNodeLocationDataOnly(.01, .01);
		tt.addNodeLocationDataOnly(.11, .11);
		tt.addNodeLocationDataOnly(.21, .21);
		tt.addNodeLocationDataOnly(.31, .31);
		tt.addNodeLocationDataOnly(.41, .41);
		tt.addNodeLocationDataOnly(.51, .51);
		tt.addNodeLocationDataOnly(.61, .61);
		tt.addNodeLocationDataOnly(.71, .71);
		tt.addNodeLocationDataOnly(.81, .81);
		tt.addNodeLocationDataOnly(.92, .92);
		tt.addNodeLocationDataOnly(.1, .1);
		tt.addNodeLocationDataOnly(.34, .34);
		tt.addNodeLocationDataOnly(.1, 1);
		tt.addNodeLocationDataOnly(.23, .78);
		tt.addNodeLocationDataOnly(.4, .1);
		tt.addNodeLocationDataOnly(.1, .5);
		tt.addNodeLocationDataOnly(.5, .5);
		tt.addNodeLocationDataOnly(.91, .91);
		tt.addNodeLocationDataOnly(.8, .7);
		tt.addNodeLocationDataOnly(.34, .78);
		tt.addNodeLocationDataOnly(.12, .1);
		tt.addNodeLocationDataOnly(.34, .342);
		tt.addNodeLocationDataOnly(.32, .83);
		tt.addNodeLocationDataOnly(.4, .34);
		tt.addNodeLocationDataOnly(.4, .31);
		tt.addNodeLocationDataOnly(.1, .35);
		tt.addNodeLocationDataOnly(.15, .95);
		tt.addNodeLocationDataOnly(.11, .91);
		tt.addNodeLocationDataOnly(.81, .71);
		tt.addNodeLocationDataOnly(.134, .178);
		
		
		
		tt.addConnection(0,1);
		tt.addConnection(0,3);
	
		tt.addConnection(1,2);
		tt.addConnection(1,3);
		tt.addConnection(1,4);
	

		tt.addConnection(2,3);
			
		tt.addConnection(5,0);
		tt.addConnection(5,4);
		
		tt.addConnection(6,1);
		tt.addConnection(6,2);
		
		tt.addConnection(7,6);
		tt.addConnection(7,0);
		
	
		tt.addConnection(8,3);
		tt.addConnection(8,5);
		
		tt.addConnection(9,2);
		tt.addConnection(9,4);		
		tt.addConnection(9,12);
		tt.addConnection(9,14);
		tt.addConnection(9,13);
		tt.addConnection(9,16);
		
		
		tt.addConnection(10,14);
		tt.addConnection(10,1);
		tt.addConnection(10,3);
		
		tt.addConnection(11,14);
		tt.addConnection(11,4);
		tt.addConnection(11,5);
		tt.addConnection(11,17);
		
		tt.addConnection(12,14);
		tt.addConnection(12,1);
		tt.addConnection(12,13);
		tt.addConnection(12,7);


		tt.addConnection(13,14);
		tt.addConnection(13,10);
		tt.addConnection(13,19);
		
		tt.addConnection(14,1);
		tt.addConnection(14,8);
		tt.addConnection(14,3);
		
		tt.addConnection(15,14);
		
		tt.addConnection(16,14);
		
		tt.addConnection(17,14);
		
		tt.addConnection(18,14);
		
		tt.addConnection(19,14);
		

		
		return tt;
	}


	public static Graph createTestGraph17() throws Exception{
		Graph tt = new Graph(0);
		
		tt.addNode();
		tt.addNode();
		tt.addNode();
		tt.addNode();
		tt.addNode();
		tt.addNode();
		tt.addNode();
		tt.addNode();

		
		tt.addConnection(0, 1);
		tt.addConnection(0, 2);
		tt.addConnection(0, 3);
		
		tt.addConnection(1, 7);
		
		//tt.addConnection(3, 0);
		//tt.addConnection(3, 2);
		tt.addConnection(2, 7);
		
		tt.addConnection(3,5);
		
		tt.addConnection(5,6);
		
		tt.addConnection(6,7);
		
		tt.addConnection(7,4);
		
		return tt;
	}

	/**
	 * This will create a square graph where each node 0 is connected to
	 * the other 3 nodes
	 * @return
	 * @throws Exception
	 */
	public static Graph createSquareGraph() throws Exception 
	{
		Graph tt = new Graph(0);
		
		tt.addNodeLocationDataOnly(.1,.1);
		tt.addNodeLocationDataOnly(.1,.2);
		tt.addNodeLocationDataOnly(.2, .1);
		tt.addNodeLocationDataOnly(.2, .2);


		
		tt.addConnection(0, 1);
		tt.addConnection(0, 2);
		tt.addConnection(0, 3);
		

		
		return tt;

	}

	/**
	 * This will create a square graph where each node is connected to
	 * 2 other nodes, thus making a classic square.
	 * @return
	 * @throws Exception
	 */
	public static Graph createSquareGraph2() throws Exception 
	{
		Graph tt = new Graph(0);
		
		tt.addNodeLocationDataOnly(.1,.1);
		tt.addNodeLocationDataOnly(.1,.2);
		tt.addNodeLocationDataOnly(.2, .1);
		tt.addNodeLocationDataOnly(.2, .2);

		tt.addConnection(0, 1);
		tt.addConnection(0, 2);
		
		tt.addConnection(1, 3);
		
		tt.addConnection(2, 3);
		
		return tt;

	}

}