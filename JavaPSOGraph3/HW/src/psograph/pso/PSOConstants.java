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
package psograph.pso;

public class PSOConstants {

	public static int NUM_PARTICLES = 40;
	
	public static int PARTICLES_BAD = 10;
	public static int PARTICLES_GOOD = 15;
	public static int PARTICLES_RANDOM = 15;
	
	public static int LOCAL_SIZE = 4;  
	
	public static int MAX_THREADS = 2;
	
	public static boolean PERTURB_RUN = true;
	
	public static int PERTURB_MINIMUM = 5;
	public static int PERTURB_MAXIMUM = 15;
	
	public static int PERTURB_ITERATIONS = 100;

	public static int PERTURB_STALENESS = 10;
	
	public static int STALENESS = 30;
	
	public static double V_MAX = 2;
	//Vmaxx set to 10,000 (no vmax) 
	// 2
	// 4 seems to produce good results
	//  8 
	//  6
	
	public static double WEIGHT_GLOBAL = .2;
	public static double WEIGHT_LOCAL = .3;
	public static double WEIGHT_PERSONAL = .4;
	
	
	public static int ITERATIONS = 1500;
	
	public enum LinkType { PERSONAL, LOCAL, GLOBAL }


}
