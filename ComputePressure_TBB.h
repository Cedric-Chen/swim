#pragma once

#include <iostream>

#include "Layer.h"
#include "Operators_DFT.h"
using namespace std;

template <int sizeX, int sizeY, int idir=0, bool bPeriodic=true>
class ComputePressure_TBB
{
private:
    double m_dt;
    double m_nu;
    double m_g[2];
//    Layer<sizeX, sizeY, 1>* m_density;
//    Layer<sizeX, sizeY, 1>* m_vorticity;
    Layer<sizeX, sizeY, 2>* m_velocityOld;
    Layer<sizeX, sizeY, 2>* m_velocityNew;
    Layer<sizeX, sizeY, 1>* m_pressure;

public:
//void setup( double dt, Layer<sizeX, sizeY, 1>& density, 
//        Layer<sizeX, sizeY, 1>& vorticity, 
//        Layer<sizeX, sizeY, 2>& velocityOld, 
//        Layer<sizeX, sizeY, 2>& velocityNew, const double * g, double nu)
void setup( double dt, const double * g, double nu, 
        Layer<sizeX, sizeY, 2>& velocityOld, 
        Layer<sizeX, sizeY, 2>& velocityNew,
        Layer<sizeX, sizeY, 1>& pressure)
{
	m_dt = dt;
	m_nu = nu;
	m_g[0] = g[0];
	m_g[1] = g[1];
//	m_density = &density;
//	m_vorticity = &vorticity;
	m_velocityOld = &velocityOld;
	m_velocityNew = &velocityNew;
    m_pressure = &pressure;
}

void operator()(int start[2]=NULL, int end[2]=NULL)
{
    cout <<"hahaha\n";
	const double dt = m_dt;
	const double nu = m_nu;
	const double g[2] = {m_g[0],m_g[1]};
//	Layer<sizeX, sizeY, 1>& density = *m_density;
//	Layer<sizeX, sizeY, 1>& vorticity = *m_vorticity;
	Layer<sizeX, sizeY, 2>& velocityOld = *m_velocityOld;
	Layer<sizeX, sizeY, 2>& velocityNew = *m_velocityNew;
    Layer<sizeX, sizeY, 1>& pressure = *m_pressure;

    Layer<sizeX, sizeY, 2> gradPressure{};
    Layer<sizeX, sizeY, 1> minusLaplacePressure{};
    gradPressure = 0;
    minusLaplacePressure = 0;

	const int startX = start==NULL?(bPeriodic?0:1):start[0];
	const int startY = start==NULL?(bPeriodic?0:1):start[1];
	const int endX = end==NULL?(bPeriodic?sizeX:(sizeX-1)):end[0];
	const int endY = end==NULL?(bPeriodic?sizeY:(sizeY-1)):end[1];
	
	
	const double factor = 0.5/(double)velocityNew.getH0();
	const double invH2 = 1.0/((double)velocityNew.getH0()*(double)velocityNew.getH0());
	
//	double gradRhoX = 0.0;
//	double gradRhoY = 0.0;
	double du_dt_X = 0.0;
	double du_dt_Y = 0.0;
	double uX = 0.0;
	double uY = 0.0;
	double vX = 0.0;
	double vY = 0.0;
	double laplaceUX = 0.0;
	double laplaceUY = 0.0;
	
	for (int dy=startY; dy<endY; dy++)
		for (int dx=startX; dx<endX; dx++)
		{
//			const double rho_left = (double)density.read((dx-1+sizeX) % sizeX,dy,0);
//			const double rho_right = (double)density.read((dx+1) % sizeX, dy, 0);
//			const double rho_down = (double)density.read(dx,(dy-1+sizeY) % sizeY, 0);
//			const double rho_up = (double)density.read(dx, (dy+1) % sizeY, 0);
//			const double rho = (double)density.read(dx,dy,0);
			
//			gradRhoX = (rho_right-rho_left)*factor;
//			gradRhoY = (rho_up-rho_down)*factor;
			
			du_dt_X = ((double)velocityNew.read(dx,dy,0) - (double)velocityOld.read(dx,dy,0)) / dt;
			du_dt_Y = ((double)velocityNew.read(dx,dy,1) - (double)velocityOld.read(dx,dy,1)) / dt;
			
			const double u_left = (double)velocityNew.read((dx-1+sizeX) % sizeX,dy,0);
			const double u_right = (double)velocityNew.read((dx+1) % sizeX, dy, 0);
			const double u_down = (double)velocityNew.read(dx,(dy-1+sizeY) % sizeY, 0);
			const double u_up = (double)velocityNew.read(dx, (dy+1) % sizeY, 0);
			const double u = (double)velocityNew.read(dx,dy,0);
			
			uX = (u_right-u_left)*factor;
			uY = (u_up-u_down)*factor; 
			
			const double v_left = (double)velocityNew.read((dx-1+sizeX) % sizeX,dy,1);
			const double v_right = (double)velocityNew.read((dx+1) % sizeX, dy, 1);
			const double v_down = (double)velocityNew.read(dx,(dy-1+sizeY) % sizeY, 1);
			const double v_up = (double)velocityNew.read(dx, (dy+1) % sizeY, 1);
			const double v = (double)velocityNew.read(dx,dy,1);
			
			vX = (v_right-v_left)*factor;
			vY = (v_up-v_down)*factor;
			
			const double unX[4] = {
				(double)velocityNew.read((dx+1+sizeX) % sizeX,dy,0),
				(double)velocityNew.read((dx-1+sizeX) % sizeX,dy,0),
				(double)velocityNew.read(dx,(dy+1+sizeY) % sizeY,0),
				(double)velocityNew.read(dx,(dy-1+sizeY) % sizeY,0)
			};
			const double uxx = (double)velocityNew.read(dx,dy,0);
			laplaceUX = (unX[0]+unX[1]+unX[2]+unX[3] - 4*uxx)*invH2;
			
			const double unY[4] = {
				(double)velocityNew.read((dx+1+sizeX) % sizeX,dy,1),
				(double)velocityNew.read((dx-1+sizeX) % sizeX,dy,1),
				(double)velocityNew.read(dx,(dy+1+sizeY) % sizeY,1),
				(double)velocityNew.read(dx,(dy-1+sizeY) % sizeY,1)
			};
			const double uyy = (double)velocityNew.read(dx,dy,1);
			laplaceUY = (unY[0]+unY[1]+unY[2]+unY[3] - 4*uyy)*invH2;

			gradPressure(dx,dy,0) = du_dt_X + u*uX + v*uY - nu*laplaceUX - g[0];
			gradPressure(dx,dy,1) = du_dt_Y + u*vX + v*vY - nu*laplaceUY - g[1];

            //debug
//            cout <<"location: " <<dx <<"," <<dy <<"\t\t";
//            cout <<gradPressure(dx,dy,0) <<"\t";
//            cout <<gradPressure(dx,dy,1) <<"\t\n";
		
//			const double a = gradRhoX / rho;
//			const double b = gradRhoY / rho;
//			const double c = du_dt_X + u*uX + v*uY - nu*laplaceUX - g[0];
//			const double d = du_dt_Y + u*vX + v*vY - nu*laplaceUY - g[1];
		
// 			vorticity(dx,dy,0) += (Real)(-(a*d-b*c)*dt);

		}

	for (int dy=startY; dy<endY; dy++)
		for (int dx=startX; dx<endX; dx++){
            // compute lapace(p)
			const double u_left = 
                (double)gradPressure.read((dx-1+sizeX) % sizeX,dy,0);
			const double u_right = 
                (double)gradPressure.read((dx+1) % sizeX, dy, 0);
			const double v_down = 
                (double)gradPressure.read(dx,(dy-1+sizeY) % sizeY, 1);
			const double v_up = 
                (double)gradPressure.read(dx, (dy+1) % sizeY, 1);
			
            minusLaplacePressure(dx,dy,0) = 
                (u_right - u_left) * factor + (v_up - v_down) * factor;

            //debug
            cout <<"location: " <<dx <<"," <<dy <<"\t\t";
            cout <<minusLaplacePressure(dx,dy,0) <<"\t\n";

        }

    VelocitySolver_Unbounded_FFTW poisson_solver(sizeX, sizeY);
    poisson_solver(minusLaplacePressure, pressure);
}
};
