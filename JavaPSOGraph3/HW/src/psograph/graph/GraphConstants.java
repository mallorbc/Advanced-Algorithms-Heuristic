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
package psograph.graph;

/**
 * This is to hold Graph Constanst that are used in multiple places.  But I wanted to be be tweak the values easily.
 * @author Patrick
 *
 */
public class GraphConstants 
{

	public static final int ClumpinessDivisor = 10;
	//static numbers used to run simulation
	/**
	 * Max number of Seeds
	 */
	public static int MAX_SEEDS = 1;
	/**
	 * Number of iterations in PSO
	 */
	public static int MAX_ITERATIONS = 1;
	/**This the number of nodes to generated when making a new graph*/
	public static int NUM_NODES = 200;//500;
	/**Max connections for starting graphs */
	public static int MAX_CONNECTIONS = (2 * (NUM_NODES));
	
	/** This used in generated graphs in Matlab that show the edges in different colors
	 * based upon the weight of the edge.
	 */
	public static double MeasurementR = .17;
	
}
