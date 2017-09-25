#pragma once

#include <iostream>
using namespace std;

#include "Layer.h"

template <int sizeX, int sizeY, bool bPeriodic=true>
class ComputeViscousForces_TBB
{
	const Layer<sizeX, sizeY, 2>* m_velocity;
	Layer<sizeX, sizeY, 2>* m_viscousforces;
    double m_nu;

public:
	void setup(const Layer<sizeX, sizeY, 2>& velocity, 
            Layer<sizeX, sizeY, 2>& viscousforces, double nu)
	{
		m_velocity = &velocity;
		m_viscousforces = &viscousforces;
        m_nu = nu;
	}

	void operator() (int start[2]=NULL, int end[2]=NULL)
	{
		const Layer<sizeX, sizeY, 2>& velocity = *m_velocity;
		Layer<sizeX, sizeY, 2>& viscousforces = *m_viscousforces;

		const int startX = start==NULL? (bPeriodic? 0 : 1) : start[0];
		const int startY = start==NULL? (bPeriodic? 0 : 1) : start[1];
		const int endX = end==NULL? (bPeriodic? sizeX: (sizeX-1)) : end[0];
		const int endY = end==NULL? (bPeriodic? sizeY: (sizeY-1)) : end[1];

		const Real invH2 = 1./(velocity.getH0()*velocity.getH0());
		
		if (bPeriodic)
		{
            for (int idir=0; idir<2; idir++)
			for (int dy=startY; dy<endY; dy++)
			for (int dx=startX; dx<endX; dx++)
			{
				const Real un[4] = {
					velocity.read((dx+1+sizeX) % sizeX,dy,idir),
					velocity.read((dx-1+sizeX) % sizeX,dy,idir),
					velocity.read(dx,(dy+1+sizeY) % sizeY,idir),
					velocity.read(dx,(dy-1+sizeY) % sizeY,idir)
				};
				
				const Real u = velocity.read(dx, dy, idir);

				viscousforces(dx, dy, idir) = 
                    m_nu * (un[0]+un[1]+un[2]+un[3] - 4*u)*invH2;

//            //debug
//            cout <<"location: " <<dx <<"," <<dy <<"\t\t";
//            cout <<viscousforces(dx,dy,0) <<"\t";
//            cout <<viscousforces(dx,dy,1) <<"\t\n";

			}
		}
		else
		{
			Real un[4], u;

            for(int idir=0; idir<2; idir++)
			for (int dy=startY; dy<endY; dy++)
			{
				un[1] = velocity.read(startX-1 ,dy, idir);				
				u = velocity.read(startX, dy, idir);

				for (int dx=startX; dx<endX; dx++)
				{
					assert(un[1] == velocity.read(dx-1, dy, idir));
					assert(u == velocity.read(dx, dy, idir));

					un[0] = velocity.read(dx+1, dy, idir);
					un[2] = velocity.read(dx, dy+1, idir);
					un[3] = velocity.read(dx, dy-1, idir);

					viscousforces(dx, dy, idir) = 
                        m_nu * (un[0]+un[1]+un[2]+un[3] - 4*u)*invH2;
					
					un[1] = u; 
					u = un[0];
				}
			}
		}
	}
};
