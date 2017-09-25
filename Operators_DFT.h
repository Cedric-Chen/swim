#pragma once

#define _USE_MATH_DEFINES

#include <math.h>
#include <typeinfo>
#include "Layer.h"
#include "Matrix2D.h"

class VelocitySolver_Unbounded
{
protected:
	int m_nSize[2], m_nBigSize[2];
	bool m_bReady;
	
	Matrix2D<Real[2]> * m_matFourierGF;
	Matrix2D<Real[2]> * m_matFourierPsi;
	Matrix2D<Real> * m_matPsi;
	
	static inline Real _physGF(int x_, int y_, Real h)
	{
		if (!(x_==0 && y_==0))
		{
			const Real x = (x_+0.0)*h;
			const Real y = (y_+0.0)*h;
			const Real r = sqrt(x*x+y*y);
			return log(r)/(M_PI*2);
		}
		else
		{
			const Real R = h/sqrt(M_PI);
			return (2*log(R)-1)/(4*M_PI);
		}
	}
	
	virtual void _init() =0;
	virtual void _solve()=0;
	
	template <int sizeX, int sizeY>
	void _copyRHS(const Layer<sizeX, sizeY, 1>& vorticity, Matrix2D<Real>& RHS )
	{		
		RHS = 0.0;
		for(int dy=0; dy<sizeY; dy++)
			for(int dx=0; dx<sizeX; dx++)
				RHS(dx,dy) = -vorticity.read(dx,dy);
	}
	
	template <int sizeX, int sizeY>
	void _copyResults(Layer<sizeX,sizeY,1>& results )
	{	
		Matrix2D<Real>& PSI = *m_matPsi;
		for(int dy=0; dy<sizeY; dy++)
			for(int dx=0; dx<sizeX; dx++)
				results(dx,dy,0) = PSI(dx,dy);
	}
	
	template <int sizeX, int sizeY>
	void _computeVelocityFromStreamFunction(Layer<sizeX, sizeY, 2>& velocity)
	{
		const double h= velocity.getH0();
		const double factor = 0.5/h;
		Matrix2D<Real>& PSI = *m_matPsi;
		
		for(int dy=1; dy<sizeY-1; dy++)
			for(int dx=1; dx<sizeX-1; dx++)
			{
				const int dxp = (dx + 1);
				const int dxm = (dx - 1);
				const int dyp = (dy + 1);
				const int dym = (dy - 1);
				
				velocity(dx,dy,0) = (PSI(dx,dyp) - PSI(dx,dym))*factor;
				velocity(dx,dy,1) = -(PSI(dxp,dy) - PSI(dxm,dy))*factor;
			}
	}
	
	template <int sizeX, int sizeY>
	void _computeVelocityFromVelocityPotential(Layer<sizeX, sizeY, 2>& velocity)
	{
		const double h= velocity.getH0();
		const double factor = 0.5/h;
		Matrix2D<Real>& PSI = *m_matPsi;
		
		for(int dy=0; dy<sizeY; dy++)
			for(int dx=0; dx<sizeX; dx++)
			{
				const int dxp = (dx + 1);
				const int dxm = (dx - 1);
				const int dyp = (dy + 1);
				const int dym = (dy - 1);
				
				velocity(dx,dy,0) += (PSI(dxp,dy) - PSI(dxm,dy))*factor;
				velocity(dx,dy,1) += (PSI(dx,dyp) - PSI(dx,dym))*factor;
			}
	}
	
	
public:
	
	VelocitySolver_Unbounded(const int nX, const int nY): 
	m_matFourierGF(NULL), m_matPsi(NULL), m_matFourierPsi(NULL), m_bReady(false)
	{		
		if(!(nX%2))
		{
			m_nSize[0] = nX;
			m_nSize[1] = nY;
			m_nBigSize[0] = nX*2;
			m_nBigSize[1] = nY*2;
		}
		else
		{
			m_nSize[0] = nX;
			m_nSize[1] = nY;
			m_nBigSize[0] = nX*2-2;
			m_nBigSize[1] = nY*2-2;
		}
	}
	
	virtual ~VelocitySolver_Unbounded(){}
	
    template <int sizeX, int sizeY>
    void operator() (const Layer<sizeX, sizeY, 1>& input, 
            Layer<sizeX, sizeY, 1>& output)
    // lapace(output) = -input
    {
		if (!m_bReady)
		{
			_init();
			m_bReady = true;
		}
		
		_copyRHS(input, *m_matPsi);
		_solve();
		_copyResults(output);

    }

	template <int sizeX, int sizeY>
	void operator() (const Layer<sizeX, sizeY, 1>& vorticity, Layer<sizeX, sizeY, 2>& velocity)
	{
		if (!m_bReady)
		{
			_init();
			m_bReady = true;
		}
		
		_copyRHS(vorticity, *m_matPsi);
		_solve();
		_computeVelocityFromStreamFunction(velocity);
	}
	
	template <int sizeX, int sizeY>
	void operator() (const Layer<sizeX,sizeY,1>& vorticity, const Layer<sizeX,sizeY,1>& minusDivergence, Layer<sizeX,sizeY,1>& streamFunction,Layer<sizeX,sizeY,1>& velocityPotential)
	{
		if (!m_bReady)
		{
			_init();
			m_bReady = true;
		}
		
		_copyRHS(vorticity, *m_matPsi);
		_solve();
		_copyResults(streamFunction);
		//_computeVelocityFromStreamFunction(velocity);
		
		_copyRHS(minusDivergence, *m_matPsi);
		_solve();
		_copyResults(velocityPotential);
		//_computeVelocityFromVelocityPotential(velocity);
	}
};



class VelocitySolver_Unbounded_FFTW: public VelocitySolver_Unbounded
{
protected:
	
	void _init();
	void _solve();

public:
	
	VelocitySolver_Unbounded_FFTW(const int nX, const int nY): 
		VelocitySolver_Unbounded(nX, nY) {}
	
	~VelocitySolver_Unbounded_FFTW();
};

#include "Operators_DFT_FFTW.h"
