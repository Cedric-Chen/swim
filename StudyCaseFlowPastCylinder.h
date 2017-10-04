#pragma once

// Std
#include<math.h>
#include<iostream>
#include<fstream>
#include <stdio.h>
#include <sys/stat.h>

// Fluids2D
#include "MRAGProfiler.h"
#include "ArgumentParser.h"
#include "Environment.h"
#include "LayerToVTK.h"
#include "Layer.h"
#include "Saxpby_TBB.h"
#include "OperatorsTBB.h"
#include "Operators_DFT.h"
#include "ClearBoundary.h"
#include "ComputeVorticity_TBB.h"
#include "ComputeBaroclinicTerm_TBB.h"
#include "ComputeDiffusion_TBB.h"
#ifdef _m_FLOAT_AGGRESSIVE
#include "MoveGridPointsRK2_TBB_Aggressive.h"
#include "P2M_MP4_TBB_Aggressive.h"
#else
#include "MoveGridPointsRK2_TBB.h"
#include "P2M_MP4_TBB.h"
#endif
#include "ComputeRotationVelocity_TBB.h"
#include "ComputeImplicitlyUPenalized_TBB.h"
#include "DumpVorticity_TBB.h"
#include "DumpVorticityHalo_TBB.h"
#include "ComputeDensity_TBB.h"
#include "ComputeMomenta.h"
#include "KillVorticityAtBoundarySin.h"
#include "UpdateVorticityPenalization_TBB.h"
#include "UpdateVorticityDiffusionRK2_TBB.h"
#include "ComputeForces.h"
#include "DumpRaw.h"
#include "ComputeInsideOutside.h"
#include "ComputeForcesPK.h"
#include "ComputeForcesPKvort.h"
#include "CheckExistenceFile.h"
#include "CopyFileOnTheFly.h"

// Shapes
#include "Ellipse_TBB.h"
#include "Disk_TBB.h"

namespace StudyCaseFlowPastCylinder
{	
	template<int sizeX, int sizeY>
	void StudyCaseFlowPastCylinder()
	{
		// Instantiate profiler
		MRAG::Profiler profiler;
		profiler.reset();
		FILE * profilerFile;
		profilerFile = fopen("profiler.txt","w");

		// Set up writing flags
		const bool append = true;
		const bool write = false;

		// Open report file
		FILE * pFile;
		pFile = fopen("report.txt","w");
		
		// Non-dimensional quantities
		double lambda = 1e4;
		double LCFL = 0.1;
		double CFL = LCFL;
		double Re = 1000.0;
		
		// Quantities in real space
		double Ut_infinity = -0.1;
		double D = 0.2;
		double nu = fabs(Ut_infinity)*D/Re;
		
		// Additional settings (time)
		unsigned int framePerSecond = 30;
		int istep = 0;
		double time = 0.0;
		double dt = 0.0;
		double terminalTime = 6.0;
		double fotoTimer = 1.0/(double)framePerSecond;
		
		// Additional quantities (shape)
		double CM[2] = {0.5,(double)sizeY/(double)sizeX-1.0*D};
		//double CM[2] = {0.5,0.5};
		double V[2]={0.0,Ut_infinity};
		double F[2]={0.0,0.0};
		double volume = 0.0;
		
		// Quantitities that do not need to be converted
		const double dx = 1./sizeX;
		const double mollificationFactor = 2.0;
		const double epsilon = mollificationFactor*sqrt(2.0)*dx;
		
		// Instantiate layers
		Layer<sizeX, sizeY,1> & omegastar = *(new Layer<sizeX, sizeY,1>);
		Layer<sizeX, sizeY,1> & vorticity = *(new Layer<sizeX, sizeY,1>);
		Layer<sizeX, sizeY,2> & ustar = *(new Layer<sizeX, sizeY,2>);
		Layer<sizeX, sizeY,2> & ustar2 = *(new Layer<sizeX, sizeY,2>);
		Layer<sizeX, sizeY,2> & old_ustar  = *(new Layer<sizeX, sizeY,2>);
		Layer<sizeX, sizeY,1> & charFunction = *(new Layer<sizeX, sizeY,1>);
		Layer<sizeX, sizeY,1> & density = *(new Layer<sizeX, sizeY,1>);
		Layer<sizeX, sizeY,1> & new_vorticity = *(new Layer<sizeX, sizeY,1>);
		Layer<sizeX, sizeY,1> & new_vorticity2 = *(new Layer<sizeX, sizeY,1>);
		Layer<sizeX, sizeY,2> & new_velocity = *(new Layer<sizeX, sizeY,2>);
		Layer<sizeX, sizeY,2> & deformationVel = *(new Layer<sizeX, sizeY,2>);
		Layer<sizeX, sizeY,2> & rotationVel = *(new Layer<sizeX, sizeY,2>);
		Layer<sizeX, sizeY,3> & particles = *(new Layer<sizeX, sizeY,3>);
				
		// Initialize fields to zero
		omegastar = 0;
		vorticity = 0;
		ustar = 0;
		ustar2 = 0;
		old_ustar  = 0;
		charFunction = 0;
		density = 1.0;
		new_vorticity = 0;
		new_vorticity2 = 0;
		new_velocity = 0;
		deformationVel = 0;
		rotationVel = 0;
		particles = 0;
		
		// Initial statements
		printf("\n\n\n#----------------------------------------------\n");
		printf("#             SIMULATION STARTED!              \n");
		printf("#----------------------------------------------\n");
		printf("Size: %dx%d\n", sizeX, sizeY);
		printf("Lambda=%f\n", lambda);
		printf("LCFL=%f\n", LCFL);
		printf("SCALE=%f\n", 1.0);
		printf("nu_Real=%f\n", nu);
		printf("D_Real=%f\n", D);
		printf("epsilon=%f\n", epsilon);
		printf("mollificationFactor=%f\n", mollificationFactor);
#ifdef _m_DOUBLE
		printf("Precision=DOUBLE\n");
#endif
#ifdef _m_FLOAT
		printf("Precision=FLOAT\n");
#endif
#ifdef _m_FLOAT_AGGRESSIVE
		printf("Precision=FLOAT AGGRESSIVE\n");
#endif
		printf("\n");		
		fprintf(pFile,"\nSIMULATION STARTED!\n\n");
		fprintf(pFile,"Size: %dx%d\n", sizeX, sizeY);
		fprintf(pFile,"Lambda=%f\n", lambda);
		fprintf(pFile,"LCFL=%f\n", LCFL);
		fprintf(pFile,"SCALE=%f\n", 1.0);
		fprintf(pFile,"nu_Real=%f\n", nu);
		fprintf(pFile,"D_Real=%f\n", D);
		fprintf(pFile,"epsilon=%f\n", epsilon);
		fprintf(pFile,"mollificationFactor=%f\n", mollificationFactor);
#ifdef _m_DOUBLE
		fprintf(pFile,"Precision=DOUBLE\n");
#endif
#ifdef _m_FLOAT
		fprintf(pFile,"Precision=FLOAT\n");
#endif
#ifdef _m_FLOAT_AGGRESSIVE
		fprintf(pFile,"Precision=FLOAT AGGRESSIVE\n");
#endif
		fprintf(pFile,"\n");
		fflush(pFile);
		
		VelocitySolver_Unbounded_FFTW velsolver(sizeX, sizeY);
		ComputeImplicitlyUPenalized_TBB<sizeX,sizeY> & implicit_u = *(new ComputeImplicitlyUPenalized_TBB<sizeX,sizeY>);
		UpdateVorticityPenalization_TBB<sizeX,sizeY> & updateomega_1 = *(new UpdateVorticityPenalization_TBB<sizeX,sizeY>);
		UpdateVorticityDiffusionRK2_TBB<sizeX,sizeY> & updateomega_2=  *(new UpdateVorticityDiffusionRK2_TBB<sizeX,sizeY>);
#ifdef _m_FLOAT_AGGRESSIVE
		MoveGridPointsRK2_TBB_Aggressive<sizeX,sizeY> & movegridpoints =  *(new MoveGridPointsRK2_TBB_Aggressive<sizeX,sizeY>);
		P2M_MP4_TBB_Aggressive<sizeX,sizeY> & p2m = *(new P2M_MP4_TBB_Aggressive<sizeX,sizeY>);
#else
		MoveGridPointsRK2_TBB<sizeX,sizeY> & movegridpoints =  *(new MoveGridPointsRK2_TBB<sizeX,sizeY>);
		P2M_MP4_TBB<sizeX,sizeY> & p2m = *(new P2M_MP4_TBB<sizeX,sizeY>);
#endif
		Disk_TBB<sizeX,sizeY> & shape = *(new Disk_TBB<sizeX,sizeY>);
		KillVorticityAtBoundarySin antiperiodic;
		ComputeForces computeForces;

		// Set up data dumper
		DumpRaw dumpRaw("dumpRaw");

		// Print header simulation
		DumpRaw dumpHeader("header.txt");
		std::vector<double> headerData;
		headerData.push_back(sizeX);
		headerData.push_back(sizeY);
		headerData.push_back(V[1]);
		headerData.push_back(lambda);
		headerData.push_back(LCFL);
		headerData.push_back(D);
		headerData.push_back(nu);
		dumpHeader.lightDump(headerData,write);
		std::cout << std::endl << std::endl;
		
		// Main loop
		while(true)
		{
			// Create characteristic function
			profiler.getAgent("mollification [TBB] ").start();
			charFunction = 0.0;
			double R = D/2.0;
			shape.setup(R,0.0,epsilon,CM,charFunction);
			processTBB<sizeX, sizeY>(shape);
			profiler.getAgent("mollification [TBB] ").stop();
			
			// Kill vorticity at the boundary
			profiler.getAgent("antiperiodic [no] ").start();
			antiperiodic(vorticity, 10);
			profiler.getAgent("antiperiodic [no] ").stop();
			
			// Solve Poisson equation for stream function
			profiler.getAgent("streamFunction [TBB] ").start();
			omegastar = vorticity;
			velsolver(omegastar, ustar);
			profiler.getAgent("streamFunction [TBB] ").stop();
			
			// Set dt
			vector<double> maxVort = vorticity.computeNorms();
			vector<double> maxVel = ustar.computeNorms();
			double limitDiffusion = 0.9*(dx)*(dx)/(4*nu);
			double limitLCFL = LCFL/maxVort[0];
			double limitCFL = dx/max(maxVel[0],maxVel[1])*CFL;
			dt = min(limitLCFL,limitDiffusion);
			if(istep<200){ dt = 1e-6; }
			
			// Penalize velocity
			profiler.getAgent("penalization [TBB] ").start();
			implicit_u.setup(lambda, dt, charFunction, ustar, new_velocity, deformationVel, V, rotationVel);
			processTBB<sizeX, sizeY>(implicit_u);
			profiler.getAgent("penalization [TBB] ").stop();
			
			// Update vorticity: given penalized velocity
			profiler.getAgent("vorticityPenal [TBB] ").start();
			updateomega_1(vorticity, new_velocity, ustar, new_vorticity);
			profiler.getAgent("vorticityPenal [TBB] ").stop();
			
			// Update vorticity: diffusion
			profiler.getAgent("diffusion [TBB] ").start();
			updateomega_2(dt, nu, new_vorticity, new_vorticity, new_vorticity2);
			profiler.getAgent("diffusion [TBB] ").stop();
			
			// Advect particles
			profiler.getAgent("advect [TBB] ").start();
			movegridpoints.setup(dt, new_vorticity2, new_velocity, particles);
			processTBB<sizeX, sizeY>(movegridpoints);
			profiler.getAgent("advect [TBB] ").stop();
			
			// Remesh
			profiler.getAgent("remesh [TBB] ").start();
			p2m.setup(particles, vorticity);
			p2m();
			profiler.getAgent("remesh [TBB] ").stop();
			
			// Additional calculations
			double modV = sqrt(V[0]*V[0]+V[1]*V[1]);
			double Re = modV*D/nu;
			double Re_h_vort = maxVort[0]*dx*dx/nu;
			double Re_h_vel = max(maxVel[0],maxVel[1])*dx/nu;
			double dimensionlessTime = 2.0*fabs(V[1])*time/D;

			// Print out stuff
			if( (istep%25) == 0 )
			{
				// Compute forces
				profiler.getAgent("forces [no] ").start();
				computeForces(F,V,volume,lambda,charFunction,density,deformationVel,rotationVel,new_velocity);
				const double cD = 2.0*F[0]/(V[1]*V[1]*D);
				const double cL = 2.0*F[1]/(V[1]*V[1]*D);
				profiler.getAgent("forces [no] ").stop();

				profiler.getAgent("print stuff [no] ").start();
				std::vector<double> stuff;
				stuff.push_back(time);
				stuff.push_back(CM[0]);
				stuff.push_back(CM[1]);
				stuff.push_back(V[0]);
				stuff.push_back(V[1]);
				stuff.push_back(0.0);
				stuff.push_back(F[0]);
				stuff.push_back(F[1]);
				stuff.push_back(Re);
				stuff.push_back(Re_h_vort);
				stuff.push_back(Re_h_vel);
				stuff.push_back(cD);
				stuff.push_back(cL);
				stuff.push_back(dimensionlessTime);
								
				dumpRaw.lightDump(stuff,append);
				
				printf("step=%d --> time=%f, Re=%f, Re_h=%f, dt=%e\n",istep,time,Re,Re_h_vort,dt);
				printf("\tLCFL=%e,CFL=%e,diff=%e,maxVort=%e,maxVel=%e,mindx=%e,dx=%e\n\n",limitLCFL,limitCFL,limitDiffusion,maxVort[0],max(maxVel[0],maxVel[1]),pow(Re,-0.5)/6.0,dx);
				fprintf(pFile,"step=%d --> time=%f, Re=%f, Re_h=%f, dt=%e\n",istep,time,Re,Re_h_vort,dt);
				fprintf(pFile,"\tLCFL=%e,CFL=%e,diff=%e,maxVort=%e,maxVel=%e,mindx=%e,dx=%e\n\n",limitLCFL,limitCFL,limitDiffusion,maxVort[0],max(maxVel[0],maxVel[1]),pow(Re,-0.5)/6.0,dx);
				fflush(pFile);
				profiler.getAgent("print stuff [no] ").stop();
			}

			if(istep > 0 && istep%500 == 0)
			{
				profiler.printSummary();
				profiler.printSummaryToFile(profilerFile);
			}

			// Visualization
			profiler.getAgent("visualization [no] ").start();
			fotoTimer += dt;
			if(fotoTimer > 1.0/(double)framePerSecond)
			{
				fotoTimer = 0.0;
				dumpLayer2VTK(istep, "vorticity", vorticity);
			}
			profiler.getAgent("visualization [no] ").stop();
			
			// Terminal time
			if(dimensionlessTime>=terminalTime){ exit(0); }

			// Update time
			time += dt;
			//CM[0] += V[0]*dt;
			CM[1] += V[1]*dt;
			
			// Update counter
			istep++;
		}
	}
	
}


