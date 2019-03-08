// This is a library to be used to represent a Graph and various measurments for a Graph
//  and to perform optimization using Particle Swarm Optimization (PSO)
//    Copyright (C) 2008, 2015 Patrick Olekas
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
package psograph.measurements;

import java.io.Serializable;
import java.util.TreeMap;
import java.util.Vector;

import psograph.graph.Graph;





public class AISPL extends ASPL implements Serializable, IGraphMeasument {

	
	/**
	 * 
	 */
	private static final long serialVersionUID = -5536904660766078822L;

	public AISPL(Graph g)
	{
		super(g);
	}
	
	
	public double Measure() throws Exception
	{
		double value =0;
		value = determineSPL(true);
		//value = peformCalculation();
		//calcValues.printWithWeights();
		
		return value;
	}

	public TreeMap<Double, Integer> getISPLDistribution() throws Exception
	{
		TreeMap<Double, Integer> result = new TreeMap<Double, Integer>();
		int i =0;
		try
		{

			for( i = 0; i < m_graph.getNumberOfNodes(); i++) 
			{
				TreeMap<Double, Integer> iteration = m_solns[i].getInverseWeightDistribution();

				Vector<Double> vWeight = new Vector<Double>(iteration.keySet());
				for(int j =0; j < vWeight.size(); j++)
				{
					double val = vWeight.get(j);
					int count = iteration.get(val);


					if(result.containsKey(val))
					{
						count += result.get(val);
					}
					result.put(val, count);
				}

			}

		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		return result;
	}
	

	

}
