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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Random;
import java.util.TreeMap;
import java.util.Vector;

import psograph.analysis.GenerateFitnessFunctionDistribution;
import psograph.graph.CalculatedGraph;
import psograph.graph.GraphConstants;
import psograph.util.Util;




public class PSO extends Thread
{
	public PSO(File seedDirectory, int id) throws Exception
	{
		m_seedDirectory = seedDirectory;
		
		m_graphsDirectory = new File(m_seedDirectory.getCanonicalPath()+"\\graphs");
		
		m_particles = new CalculatedGraph[PSOConstants.NUM_PARTICLES];
		
		m_localBest = new CalculatedGraph[PSOConstants.NUM_PARTICLES/PSOConstants.LOCAL_SIZE];
		m_personalBest = new CalculatedGraph[PSOConstants.NUM_PARTICLES];	
		
		m_stalenessCount =0;
		m_peturb_iteration = 0;
		
		m_id = id;
		
		r = new Random();
		
		bestFitnessFromPreviousPerturbRun =0;
		
		
	}

	public void run() 
	{
		try
		{
			
			//File graphDir= new File(m_seedDirectory.getCanonicalPath()+"\\graphs");
			
			int m_peturb_iteration;
			if(PSOConstants.PERTURB_RUN)
				m_peturb_iteration=0;
			else
				m_peturb_iteration = PSOConstants.PERTURB_ITERATIONS -1;
			
			for(; m_peturb_iteration < PSOConstants.PERTURB_ITERATIONS; m_peturb_iteration++)
			{
				m_canditateDir = new File (m_seedDirectory.getAbsolutePath() + "\\PSO_"+m_peturb_iteration);
				m_canditateDir.mkdir();
				
				doWork();
				finalReport();
				
				if(PSOConstants.PERTURB_RUN)
				{
					m_graphsDirectory = m_canditateDir;
					
					if(m_globalBest.getFitnessValue() > bestFitnessFromPreviousPerturbRun)
					{
						System.out.println("Thd_id - " + m_id + "Perturb Run found a better solution");
						System.out.println("Thd_id - " + m_id + "old - "+bestFitnessFromPreviousPerturbRun);
						System.out.println("Thd_id - " + m_id + "new - "+m_globalBest.getFitnessValue());
						bestFitnessFromPreviousPerturbRun = m_globalBest.getFitnessValue();
					}
					else
					{
						System.out.println("Thd_id - " + m_id + "Older solution better or equal to previous solution, stopping run");
					
						m_peturb_iteration = PSOConstants.PERTURB_ITERATIONS;

					}

					m_stalenessCount=0;
					m_i = 0;
				}
				
			}

		}
		catch (Exception e)
		{
			e.printStackTrace();
			
		}
		
	}
	

	/**
	 * This method will perturb the particles.  Specifically it will
	 * Pick two random Nodes.  And if they are connected it will sever the connection.
	 * If they are not connected it will connect them.
	 * 
	 * It add/sever edges for a random number of times.  The number of times will be in the 
	 * range between PSOConstants.PERTURB_MINIMUM and PSOConstants.PERTURB_MAXIMUM
	 * @param particle Graph to perturbed
	 * @throws Exception
	 */
	void perturbParticle(CalculatedGraph particle) throws Exception
	{
		
		int base = r.nextInt(1+ PSOConstants.PERTURB_MAXIMUM - PSOConstants.PERTURB_MINIMUM);
		
		int perturbAmount = PSOConstants.PERTURB_MINIMUM + base;
		
		for(int i =0; i < perturbAmount; i++)
		{
			int pickA, pickB;
			
			//Olekas May-19-2012 had this using the constant which is wrong
			//should be indicated by the number of Nodes in the graph
			pickA = r.nextInt(particle.getNumberOfNodes());
			pickB = r.nextInt(particle.getNumberOfNodes());
			while(pickA == pickB)
			{
				pickB = r.nextInt(GraphConstants.NUM_NODES);				
			}
			
			if(!particle.isNodeAConnectedToNodeB(particle.getNode(pickA), particle.getNode(pickB)))
			{
				particle.addConnection(pickA, pickB);
			}
			else
			{
				particle.removeConnection(pickA, pickB);
			}
		}
	}
	
	void perturbParticles() throws Exception
	{
		if(!PSOConstants.PERTURB_RUN)
			return;
		
		InitCalculateBests();
		
		for(int i =0; i < PSOConstants.NUM_PARTICLES; i++)
		{   
			perturbParticle(m_particles[i]);
		}
		
		updatePSOParticlesInfo();
		
	}
	
	void alternativeInititialize(Vector<String> valueSetStrings ) throws Exception
	{
		
		int i =0;
		TreeMap<Integer,Boolean> picked_so_far = new TreeMap<Integer,Boolean>();
		
		int pick = 0;
		for(;picked_so_far.size() < PSOConstants.NUM_PARTICLES;)
		{
			pick = r.nextInt(40);
			if(picked_so_far.get(pick)!=null)
			{
				continue;
			}
			else
			{
				File t = new File(m_graphsDirectory.getCanonicalPath() +"\\"+valueSetStrings.get(i));
				
				if(!t.exists())
					throw new Exception("Non-existant calculated graph - " +t.getCanonicalPath());
				
				m_particles[pick] = Util.streaminCalculatedGraph(t);
				
				i++;
				picked_so_far.put(pick, true);
			}
		}

		//InitCalculateBests();

	}
	
	void initialize() throws Exception
	{
		//Need to stream in graphs and fill out initial particles.
		
		TreeMap<Double, Vector<String>> fitnessDistribution;
		
		fitnessDistribution = Util.streaminFitnessDistribution(new File(m_graphsDirectory.getCanonicalPath() + "\\FitnessDistribution.FitnessDistribution"));
		
		Vector<Double> keySet = new Vector<Double>(fitnessDistribution.keySet());
		

		
		int sizeKeySet = keySet.size();
		if(sizeKeySet < PSOConstants.NUM_PARTICLES)
		{
			Vector<Vector<String>> valueSet = new Vector<Vector<String>>(fitnessDistribution.values() );
			int valueSetSize = valueSet.size();

			Vector<String> valueSetStrings = new Vector<String>();
			
			for(int i = 0; i < valueSetSize ; i++)
			{
				Vector<String> val = valueSet.get(i);
				for(int j=0; j < val.size(); j++)
				{
					valueSetStrings.add(val.get(j));
				}				
			}
			
			int sizeStrings = valueSetStrings.size();
			if(sizeStrings < PSOConstants.NUM_PARTICLES)
				throw new Exception("Too few particles to do a PSO - skipping");
			else
			{
				 alternativeInititialize(valueSetStrings);
				 return;
				
			}
		}
		
		System.out.println("Thd_id - " + m_id + " number of particle candidates for "+m_graphsDirectory.getName()+ " is "+sizeKeySet);
		
	
		TreeMap<Integer,Boolean> picked_so_far = new TreeMap<Integer,Boolean>();
		TreeMap<Integer,Boolean> picked_so_far_random = new TreeMap<Integer,Boolean>();
		
		int i =0;
		int pick = 0;
		for(;picked_so_far.size() < PSOConstants.PARTICLES_BAD;)
		{
			pick = r.nextInt(40);
			if(picked_so_far.get(pick)!=null)
			{
				continue;
			}
			else
			{
				Vector<String> m = fitnessDistribution.get(keySet.get(i));
				File t = new File(m_graphsDirectory.getCanonicalPath() +"\\"+m.get(0));
				
				if(!t.exists())
					throw new Exception("Non-existant calculated graph - " +t.getCanonicalPath());
				
				m_particles[pick] = 
					Util.streaminCalculatedGraph(t);
				i++;
				picked_so_far.put(pick, true);
			}
		}
		
		i =1;
		pick = 0;
		for(;picked_so_far.size() < PSOConstants.PARTICLES_GOOD + PSOConstants.PARTICLES_BAD;)
		{
			pick = r.nextInt(40);
			if(picked_so_far.get(pick)!=null)
			{
				continue;
			}
			else
			{
				Vector<String> m = fitnessDistribution.get(keySet.get(keySet.size() - i));
				File t = new File(m_graphsDirectory.getCanonicalPath() +"\\"+m.get(0));
				
				if(!t.exists())
					throw new Exception("Non-existant calculated graph - " +t.getCanonicalPath());
				
				this.m_particles[pick] = 
					Util.streaminCalculatedGraph(t);
				i++;
				picked_so_far.put(pick, true);
			}
		}
		
		pick = 0;
		for(;picked_so_far.size() < PSOConstants.NUM_PARTICLES;)
		{
			pick = r.nextInt(40);
			if(picked_so_far.get(pick)!=null)
			{
				continue;
			}
			else
			{
				int pick2 = r.nextInt(keySet.size() - PSOConstants.PARTICLES_GOOD - PSOConstants.PARTICLES_BAD);
				if(picked_so_far_random.get(pick2)!= null)
					continue;
				
				Vector<String> m = fitnessDistribution.get(keySet.get(PSOConstants.PARTICLES_BAD +pick2));
				File t = new File(m_graphsDirectory.getCanonicalPath() +"\\"+m.get(0));
				
				if(!t.exists())
					throw new Exception("Non-existant calculated graph - " +t.getCanonicalPath());
				
				this.m_particles[pick] = 
					Util.streaminCalculatedGraph(t);
				i++;
				picked_so_far.put(pick, true);
				picked_so_far_random.put(pick2,true);
			}
		}
		
		//System.out.println("Print out fitness values to see distributions");
		for(i=0;i < PSOConstants.NUM_PARTICLES;i++)
		{
			//System.out.println(this.m_particles[i].fitnessValue());
			m_personalBest[i] = new CalculatedGraph(m_particles[i]);
		}
		
		
		//InitCalculateBests();
	}
	
	void updatePSOParticlesInfo() throws Exception
	{
		//Need to update the fields that matter 
		for(int i=0;i < PSOConstants.NUM_PARTICLES;i++)
		{
			m_particles[i].UpdatePSOCalculations();
				
		}
	}
	
	void InitCalculateBests() throws Exception
	{
		int i;
		
		//System.out.println("Print out fitness values to see distributions");
		for(i=0;i < PSOConstants.NUM_PARTICLES;i++)
		{
			//System.out.println(this.m_particles[i].fitnessValue());
			m_personalBest[i] = new CalculatedGraph(m_particles[i]);
		}

		
		TreeMap<Double,Integer> fitnessDistOfParticles = new TreeMap<Double,Integer>();
		for(i=0;i < PSOConstants.NUM_PARTICLES;i++)
		{
			fitnessDistOfParticles.put(m_particles[i].getFitnessValue(), i);
		}
		
		Vector<Double> keySet = new Vector<Double>(fitnessDistOfParticles.keySet());
		
		int maxSize = PSOConstants.NUM_PARTICLES;
		//If we have graphs with the same fitness we need to lower the number
		if(keySet.size() < PSOConstants.NUM_PARTICLES)
			maxSize = keySet.size();
		

		int Global_Best = fitnessDistOfParticles.get(keySet.get(maxSize -1));
		m_globalBest = new CalculatedGraph(m_particles[Global_Best]);
		
		i=0;
		double best_fitness_seen;
		int best_fitness_index;
		for(i=0; i < PSOConstants.NUM_PARTICLES/PSOConstants.LOCAL_SIZE; i++)
		{
			best_fitness_seen =-1;
			best_fitness_index =-1;

			int j;
			for(j=0; j < PSOConstants.LOCAL_SIZE; j++)
			{
				if(m_particles[(i*PSOConstants.LOCAL_SIZE)  +j].getFitnessValue() > best_fitness_seen)
				{
					best_fitness_seen = m_particles[(i*PSOConstants.LOCAL_SIZE)  +j].getFitnessValue();
					best_fitness_index = (i*PSOConstants.LOCAL_SIZE)  +j;
				}
				m_localBest[i] = new CalculatedGraph(m_particles[best_fitness_index]);
				
			}
		}
		
	//	System.out.println();
		//System.out.println("---------------- Done updating ----now reporting-----");
		System.out.println("Thd_id - " + m_id + " Global best fitness value "+ m_globalBest.getFitnessValue());
		//System.out.println();
		
//		System.out.println("Print out fitness values of local bests");
//		for(i=0;i < PSOConstants.NUM_PARTICLES/PSOConstants.LOCAL_SIZE;i++)
//		{
//			System.out.println("Personal best of local particle "+i+" is fitness of "+m_localBest[i].fitnessValue());
//		}
		
//		System.out.println();
//		System.out.println("Print out fitness values of personal bests");
//		for(i=0;i < PSOConstants.NUM_PARTICLES;i++)
//		{
//			System.out.println("Personal best of particle "+i+" is fitness of "+m_particles[i].fitnessValue());
//		}
		
		//System.out.println("---------------- Done reporting new values ---------------------");
	//	System.out.println();
		
		System.out.println("Thd_id - " + m_id + " Total number of edges "+ m_globalBest.getTotalEdges());
		
	}
	
	//This will update the Global and Local bests
	//The personal bests will be handled during the particle velocity update
	void CalculateBests() throws Exception
	{
		int i;
		

		updatePSOParticlesInfo();

	
		
		TreeMap<Double,Integer> fitnessDistOfParticles = new TreeMap<Double,Integer>();
		for(i=0;i < PSOConstants.NUM_PARTICLES;i++)
		{
			fitnessDistOfParticles.put(m_particles[i].getFitnessValue(), i);
		}
		
		Vector<Double> keySet = new Vector<Double>(fitnessDistOfParticles.keySet());
		
		int Global_Best = fitnessDistOfParticles.get(keySet.lastElement());
		
		if(m_particles[Global_Best].getFitnessValue() > m_globalBest.getFitnessValue())
		{
			m_globalBest = new CalculatedGraph(m_particles[Global_Best]);
			m_stalenessCount =0;
		}
		else
		{
			m_stalenessCount++;
		}
		
		i=0;
		double best_fitness_seen;
		int best_fitness_index;
		for(i=0; i < PSOConstants.NUM_PARTICLES/PSOConstants.LOCAL_SIZE; i++)
		{
			best_fitness_seen =-1;
			best_fitness_index =-1;

			int j;
			for(j=0; j < PSOConstants.LOCAL_SIZE; j++)
			{
				if(m_particles[(i*PSOConstants.LOCAL_SIZE)  +j].getFitnessValue() > best_fitness_seen)
				{
					best_fitness_seen = m_particles[(i*PSOConstants.LOCAL_SIZE)  +j].getFitnessValue();
					best_fitness_index = (i*PSOConstants.LOCAL_SIZE)  +j;
				}
				
				if(m_particles[best_fitness_index].getFitnessValue() > m_localBest[i].getFitnessValue())
					m_localBest[i] = new CalculatedGraph(m_particles[best_fitness_index]);
				
			}
		}
		
		
		for(i=0;i < PSOConstants.NUM_PARTICLES;i++)
		{
			
			if(m_particles[i].getFitnessValue() > m_personalBest[i].getFitnessValue())
				m_personalBest[i] = new CalculatedGraph(m_particles[i]);
				

		}
		
		
	    //System.out.println();
		//System.out.println("---------------- Done updating ----now reporting-----");
		
		double cost = m_globalBest.getCost();	
		double CF = (cost/m_globalBest.getCostBasis());
		
		System.out.println("Thd_id - " + m_id + "	"+m_i+" "+ m_globalBest.getFitnessValue()+
				"  ("+m_globalBest.getAISPL()+
				","+CF+
				","+ m_globalBest.getLCC() +") ");
		
		m_progressText.write(m_i+" "+ m_globalBest.getFitnessValue()+
				"  ("+m_globalBest.getAISPL()+
				","+CF+
				","+ m_globalBest.getLCC() +") \n");
		
		
		
		
		//System.out.println();
		
//		System.out.println("Print out fitness values of local bests");
//		for(i=0;i < PSOConstants.NUM_PARTICLES/PSOConstants.LOCAL_SIZE;i++)
//		{
//			System.out.println("Personal best of local particle "+i+" is fitness of "+m_localBest[i].fitnessValue());
//		}
		
//		System.out.println();
//		System.out.println("Print out fitness values of personal bests");
//		for(i=0;i < PSOConstants.NUM_PARTICLES;i++)
//		{
//			System.out.println("Personal best of particle "+i+" is fitness of "+m_particles[i].fitnessValue());
//		}
		
		//System.out.println("---------------- Done reporting new values ---------------------");
	//	System.out.println();
	}
	
	void UpdateParticles() throws Exception
	{
		int i;
		for (i=0; i< PSOConstants.NUM_PARTICLES;i++)
		{
			UpdateParticle(i);			
		}
		
		
	}
	
	void resetPickLinkType()
	{
		allowPickingPersonalLinks =true;
		allowPickingLocalLinks = true;
		allowPickingGlobalLinks = true;
		
	}
	
	void disablePickingPersonalLinks()
	{
		allowPickingPersonalLinks=false;
	}
	
	void disablePickingLocalLinks()
	{
		allowPickingLocalLinks=false;
	}
	
	void disablePickingGlobalLinks()
	{
		allowPickingGlobalLinks=false;
	}
	
	void UpdateParticle(int i)throws Exception
	{
		
//		CalculatedGraph copy_Particle = new CalculatedGraph(m_particles[i]);
		
		int hammingDistanceLocal = 
			HammingDistance.CalculateHammingDistance(m_particles[i],
					m_localBest[i / PSOConstants.LOCAL_SIZE]);
		
		int hammingDistancePersonal = 
			HammingDistance.CalculateHammingDistance(m_particles[i], m_personalBest[i]);	
	

		
		int hammingDistanceGlobal = 
			HammingDistance.CalculateHammingDistance(m_particles[i], m_globalBest);
		
		double K_p = PSOConstants.WEIGHT_PERSONAL * hammingDistancePersonal;
		double K_l = PSOConstants.WEIGHT_LOCAL * hammingDistanceLocal;
		double K_g = PSOConstants.WEIGHT_GLOBAL * hammingDistanceGlobal;
		
		if(K_p > 0 && K_p < 1.0)
			K_p = 1.0;
		
		if(K_l > 0 && K_l < 1.0)
			K_l = 1.0;
		
		if(K_g > 0 && K_g < 1.0)
			K_g = 1.0;
		
		
		
		
		int NumLinksAccepted = 0;

		double NumPersonalLinks = Math.rint(K_p);
		double NumGlobalLinks = Math.rint(K_g);
		double NumLocalLinks = Math.rint(K_l);
		double maxLinks = PSOConstants.V_MAX;
		
		if(NumPersonalLinks > PSOConstants.V_MAX)
		{
			NumPersonalLinks = PSOConstants.V_MAX;
//			System.out.println("Wow we are accepting less than VMAX " +maxLinks);
		}
		
		if(NumGlobalLinks > PSOConstants.V_MAX)
		{
			NumGlobalLinks = PSOConstants.V_MAX;
//			System.out.println("Wow we are accepting less than VMAX " +maxLinks);
		}
		
		if(NumLocalLinks > PSOConstants.V_MAX)
		{
			NumLocalLinks = PSOConstants.V_MAX;
//			System.out.println("Wow we are accepting less than VMAX " +maxLinks);
		}
			
		
		HammingDistance hammingDistance = new HammingDistance();
		
		maxLinks = NumPersonalLinks + NumGlobalLinks + NumLocalLinks;
		
		resetPickLinkType();
		
		while(NumLinksAccepted < maxLinks)
		{
			PSOConstants.LinkType pickLinkType = PickLinkType();
			
			// now we pick randomly from the three different types of link
			switch(pickLinkType)
			{
				case PERSONAL:
					if(NumPersonalLinks == 0)
					{
						//pick another link type
						break;
					}
					
					//Otherwise accept a personal link;
					m_particles[i] = hammingDistance.AcceptEdge(m_particles[i], m_personalBest[i]);
					
					NumPersonalLinks--;
					NumLinksAccepted++;
					
					if(NumPersonalLinks ==0)
						disablePickingPersonalLinks();
					break;
				case LOCAL:
			
					if(NumLocalLinks == 0)
					{
						//pick another link type
						break;
					}
				
					//Otherwise accept a local link;
					m_particles[i] = hammingDistance.AcceptEdge(m_particles[i],
							m_localBest[i / PSOConstants.LOCAL_SIZE]);
				
					NumLocalLinks--;
					NumLinksAccepted++;
					
					if(NumLocalLinks ==0)
						disablePickingLocalLinks();
					break;

				case GLOBAL:
					if(NumGlobalLinks == 0)
					{
						//pick another link type
						break;
					}
					
					//Otherwise accept a global link;
					m_particles[i] = hammingDistance.AcceptEdge(m_particles[i], m_globalBest);
					
					NumGlobalLinks--;
					NumLinksAccepted++;
					
					if(NumGlobalLinks ==0)
						disablePickingGlobalLinks();
					break;

					
				
			}
		}
		
		//int hammingDistanceTotal = HammingDistance.CalculateHammingDistance(m_particles[i], copy_Particle);
		
		//System.out.println("Hamming distance after update links -" + hammingDistanceTotal);
		
	//	int pdakfj = 0;
		
		
		
	}
	
	PSOConstants.LinkType  PickLinkType() throws Exception
	{
		PSOConstants.LinkType result =PSOConstants.LinkType.PERSONAL;
		
		if(!allowPickingPersonalLinks && !allowPickingLocalLinks && !allowPickingGlobalLinks)
			throw new Exception("trying to pick a link type and not types are valid");
		
		if(allowPickingPersonalLinks && !allowPickingLocalLinks && !allowPickingGlobalLinks)
			result = PSOConstants.LinkType.PERSONAL;
		else if(!allowPickingPersonalLinks && allowPickingLocalLinks && !allowPickingGlobalLinks)
			result = PSOConstants.LinkType.LOCAL;
		else if(!allowPickingPersonalLinks && !allowPickingLocalLinks && allowPickingGlobalLinks)
			result = PSOConstants.LinkType.GLOBAL;
		else
		{
		
			boolean legal = false;
			
			while(!legal)
			{
				int pick = r.nextInt(3);
				
				if(pick == 0)// PSOConstants.WEIGHT_PERSONAL)
				{
					result = PSOConstants.LinkType.PERSONAL;
					if(allowPickingPersonalLinks)
						legal = true;
					
				}
				else if(pick == 1) // PSOConstants.WEIGHT_PERSONAL + PSOConstants.WEIGHT_LOCAL)
				{
					result = PSOConstants.LinkType.LOCAL;
					if(allowPickingLocalLinks)
						legal = true;
				}
				else //by default this means pick > PSOConstants.WEIGHT_PERSONAL + PSOConstants.WEIGHT_LOCAL
				{
					result = PSOConstants.LinkType.GLOBAL;
					if(allowPickingGlobalLinks)
						legal = true;
				}
				
			}
			
		}
		return result;
		
	}
	
	public void doWork() throws Exception
	{
		try
		{
			initialize();
			

		
			FileWriter m_progressTextFW = new FileWriter(m_canditateDir.getAbsolutePath()+"\\progress.txt");
			m_progressText = new BufferedWriter(m_progressTextFW);
		    
			m_progressText.write("iteration Global_best_fitness (AISPL, CF, LCC) \n");
			
			perturbParticles();
			InitCalculateBests();
			
			

			for(m_i=0; m_i < PSOConstants.ITERATIONS && (m_stalenessCount < PSOConstants.STALENESS ); m_i++)
			{
				
				UpdateParticles();
				CalculateBests();
			}
			
			m_progressText.close();
			
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.out.println("Thd_id - " + m_id + "Error on PSO, going to next Seed");
		}
		
		
	}

	public void finalReport() throws Exception
	{
		System.out.println("Thd_id - " + m_id + "Print out global best candidate");
		
		
		Util.i = 6;
		//m_globalBest.UpdateCalcuations();
		Util.streamOutCalculatedGraph(m_canditateDir, "globalBest_CalculatedGraph", m_globalBest);
		
		
		
		
		System.out.println("Thd_id - " + m_id + " Total number of edges "+ m_globalBest.getTotalEdges());
		
		//Now stream out the personal best 
		
		System.out.println("Thd_id - " + m_id + " location of final graphs -"+m_canditateDir.getCanonicalPath());
		for(int i =0; i <  PSOConstants.NUM_PARTICLES;i++)
		{
			//m_personalBest[i].UpdateCalcuations();
			Util.streamOutCalculatedGraph(m_canditateDir, "personalBest"+i, m_personalBest[i]);
		}
		//We need to generate the FitnessDistribution
		
		GenerateFitnessFunctionDistribution generateFitnessFunctionDistribution = 
			new GenerateFitnessFunctionDistribution(m_canditateDir);
		generateFitnessFunctionDistribution.generateDist();
		
	}
	
    public static void main(String[] args)
    {
    	try
    	{

    		Thread currentThreads[] = new Thread[PSOConstants.MAX_THREADS];
    		
    		
    		File mainDirectory= new File(psograph.util.Util.baseDirectory);
    		File seedDirectories[] = mainDirectory.listFiles();
    		
    		if(seedDirectories == null)
    			throw new Exception ("No seed directories !!!!");
    		
    		int i;
    		int maxSeedDirectories = seedDirectories.length;
    		
    		PSO pso = null;
    		//For loop to go over each directory 
    		for( i=0; i < maxSeedDirectories  && seedDirectories[i].isDirectory(); )
    		{
    			if(currentThreads[0] == null)
    			{
    			    System.out.println("Starting new thread");
	    			pso = new PSO(seedDirectories[i],i);
	    			pso.start();
	    			currentThreads[0] = pso;
	    			i++;
    			}
    			else if (currentThreads[1] == null)
    			{
    				System.out.println("Starting new thread");
	    			pso = new PSO(seedDirectories[i],i);
	    			pso.start();
	    			currentThreads[1] = pso;
	    			i++;
    			}
    			else if (!currentThreads[0].isAlive())
    			{
    				System.out.println("Starting new thread");
	    			pso = new PSO(seedDirectories[i],i);
	    			pso.start();
	    			currentThreads[0] = pso;
    				i++;
    			}
    			else if (!currentThreads[1].isAlive())
    			{
    				System.out.println("Starting new thread");
	    			pso = new PSO(seedDirectories[i],i);
	    			pso.start();
	    			currentThreads[1] = pso;
    				i++;
    			}
    			else
    				Thread.sleep(5000);
    				

    		}
    		
    		boolean working = true;
    		while(working)
    		{
    			boolean workA = true;
    			boolean workB = true;
    			
    			if(currentThreads[0] == null  && currentThreads[1] == null)
    			{
    				working = false;
    			}
    			
    			if(currentThreads[0] != null && !currentThreads[0].isAlive() )
    			{
    				workA = false;
    			}
    			else if(currentThreads[0] == null  )
    				workA = false;
    				
    				
    			if(currentThreads[1] != null && !currentThreads[1].isAlive() )
    			{
    				workB = false;	
    			}	
    			else if(currentThreads[1] == null  )
    				workB = false;	
    				
    			
    			working = workA || workB;
    				
    			Thread.sleep(5000);
    		}
    		
    		
    		System.out.println("done with all threads");

    	}
    	catch (Exception e)
    	{
    		e.printStackTrace();
    	}
    }
	
    int m_i;
    
	CalculatedGraph[] m_particles;
	
	CalculatedGraph m_globalBest;
	CalculatedGraph[] m_localBest;
	CalculatedGraph[] m_personalBest;	
	File m_seedDirectory;
	File m_graphsDirectory;
	File m_canditateDir;
	BufferedWriter m_progressText;
	
	int m_stalenessCount ;
	
	double bestFitnessFromPreviousPerturbRun;
	
	boolean allowPickingPersonalLinks, allowPickingLocalLinks, allowPickingGlobalLinks;
	
	int m_id;
	int m_peturb_iteration;
	int m_peturb_iteration_staleness;
	
	Random r;
}
