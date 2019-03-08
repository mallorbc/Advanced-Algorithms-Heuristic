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

import java.io.File;

public class MFilePlotData 
{
	String xAxisLabel;
	String yAxisLabel;
	String plotLabel;
	String fileName;
	File outputDirectory;
	String ofWhat;
	int nBins;
	
	double minBin ;
	double maxBin;
	
	//Only to be used for ValuePlots
	boolean reverseValues = true;
	
	public String getXAxisLabel() {
		return xAxisLabel;
	}
	public void setXAxisLabel(String axisLabel) {
		xAxisLabel = axisLabel;
	}
	public String getYAxisLabel() {
		return yAxisLabel;
	}
	public void setYAxisLabel(String axisLabel) {
		yAxisLabel = axisLabel;
	}
	public String getPlotLabel() {
		return plotLabel;
	}
	public void setPlotLabel(String plotLabel) {
		this.plotLabel = plotLabel;
	}
	public String getFileName() {
		return fileName;
	}
	public void setFileName(String fileName) {
		this.fileName = fileName;
	}
	
	public boolean getReverseValues() {
		return this.reverseValues;
	}
	public void setReverseValues(boolean reverseValues) {
		this.reverseValues = reverseValues;
	}
	
	public File getOutputDirectory() {
		return outputDirectory;
	}
	public void setOutputDirectory(File outputDirectory) {
		this.outputDirectory = outputDirectory;
	}
	public MFilePlotData(String axisLabel, String axisLabel2, String plotLabel,
			String fileName, File outputDirectory, String OfWhat) {
		super();
		xAxisLabel = axisLabel;
		yAxisLabel = axisLabel2;
		this.plotLabel = plotLabel;
		this.fileName = fileName;
		this.outputDirectory = outputDirectory;
		if(OfWhat == null)
			this.ofWhat = OfWhat;
		else
			this.ofWhat = "of "+OfWhat;
		
		nBins =0;
	}
	public MFilePlotData(String axisLabel, String axisLabel2, String plotLabel,
			String fileName, File outputDirectory, String OfWhat, int bins) {
		super();
		xAxisLabel = axisLabel;
		yAxisLabel = axisLabel2;
		this.plotLabel = plotLabel;
		this.fileName = fileName;
		this.outputDirectory = outputDirectory;
		if(OfWhat == null)
			this.ofWhat = OfWhat;
		else
			this.ofWhat = "of "+OfWhat;
		
		nBins = bins;
	}
	public MFilePlotData(String axisLabel, String axisLabel2, String plotLabel,
			String fileName, File outputDirectory, String OfWhat, double d,
			double e) 
	{
		super();
		xAxisLabel = axisLabel;
		yAxisLabel = axisLabel2;
		this.plotLabel = plotLabel;
		this.fileName = fileName;
		this.outputDirectory = outputDirectory;
		if(OfWhat == null)
			this.ofWhat = OfWhat;
		else
			this.ofWhat = "of "+OfWhat;
		
		nBins =0;
		
		minBin = d;
		maxBin = e;
	}
	

}
