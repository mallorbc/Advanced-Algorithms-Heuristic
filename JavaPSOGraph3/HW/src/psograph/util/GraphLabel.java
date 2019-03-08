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

public class GraphLabel
{

	public GraphLabel(String YLabel, String GraphLabel)
	{
		setYLabel(YLabel);
		setGraphLabel(GraphLabel);
	}
	public void setGraphLabel(String m_GraphLabel) {
		this.m_GraphLabel = m_GraphLabel;
	}
	public String getGraphLabel() {
		return m_GraphLabel;
	}
	public void setYLabel(String m_YLabel) {
		this.m_YLabel = m_YLabel;
	}
	public String getYLabel() {
		return m_YLabel;
	}
	private String m_YLabel;
	private String m_GraphLabel;
}
