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

import java.io.Serializable;


/**
 * Nonlinear cost function.  y(x) = 2x^2
 * @author Patrick
 *
 */
public class NonLinearCostFunction implements ICostFunction, Serializable
{

	private static final long serialVersionUID = -2786664952929279473L;

	/**
	 * Constructor
	 */
	public NonLinearCostFunction() {}
	
	/**
	 * Calculate the cost.  If the cost comes to be zero, we make the cost
	 * be .01
	 */
	public double calculate (double t_sig) throws Exception
	{
		double cost = 2.0 * Math.pow(t_sig,2);
		double cost_sig = cost;
		//System.out.println("Cost is "+cost);

		if( Double.compare(cost_sig , 0.0) == 0)
		{
			System.out.println(" cost_sig is zero!!!! Setting cost to .01");
			cost_sig = .01;
		}
		
		return cost_sig;
		
	}
	/**
	 * Calculates the distance based upon the cost.
	 * @param cost
	 * @return
	 * @throws Exception
	 */
	public static double getDistance(double cost)throws Exception
	{
		double first_term = cost/2.0;
		double final_term = Math.sqrt(first_term);
		
		return final_term;
	}
	
}
