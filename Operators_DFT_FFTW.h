#include "fftw3.h"
#include <math.h>

#include "Operators_DFT.h"
#include "InterfaceFFTW.h"

struct MetadataFFWT_float
{
	fftwf_plan plan_forward;
	fftwf_plan plan_backward;
	fftwf_plan velocityplan_backward[2];

};
struct MetadataFFWT_double
{
	fftw_plan plan_forward;
	fftw_plan plan_backward;
	fftw_plan velocityplan_backward[2];
};


VelocitySolver_Unbounded_FFTW::
~VelocitySolver_Unbounded_FFTW()
{
	InterfaceFFTW::erasePlan("RHS: Real->Fourier");
	InterfaceFFTW::erasePlan("Psi: Fourier->Real");
	InterfaceFFTW::deallocateMatrix(m_matFourierGF);
	InterfaceFFTW::deallocateMatrix(m_matFourierPsi);
	InterfaceFFTW::deallocateMatrix(m_matPsi);
}

void
VelocitySolver_Unbounded_FFTW::
_init()
{	
	InterfaceFFTW::allocateMatrix(m_matFourierGF, m_nBigSize[0]/2+1, m_nBigSize[1]);
	InterfaceFFTW::allocateMatrix(m_matFourierPsi, m_nBigSize[0]/2+1, m_nBigSize[1]);
	InterfaceFFTW::allocateMatrix(m_matPsi, m_nBigSize[0], m_nBigSize[1]);
	
	InterfaceFFTW::createPlanForward("RHS: Real->Fourier", *m_matPsi, *m_matFourierPsi);
	InterfaceFFTW::createPlanBackward("Psi: Fourier->Real", *m_matFourierPsi, *m_matPsi);
	
	{
		Matrix2D<Real> * matGF  = NULL;
		InterfaceFFTW::allocateMatrix(matGF, m_nBigSize[0], m_nBigSize[1]);
		InterfaceFFTW::createPlanForward("GF: Real->Fourier", *matGF, *m_matFourierGF);
		
		const Real H = 1./m_nSize[0];
		for(int dy=0; dy<m_nBigSize[1]; dy++)
			for(int dx=0; dx<m_nBigSize[0]; dx++)
				matGF->Access(dx, dy) = _physGF(min(dx, m_nBigSize[0]-dx), min(dy, m_nBigSize[1]-dy), H);
		
		InterfaceFFTW::executePlan("GF: Real->Fourier");
		InterfaceFFTW::erasePlan("GF: Real->Fourier");	
		InterfaceFFTW::deallocateMatrix(matGF);
	}	
}

void
VelocitySolver_Unbounded_FFTW::
_solve()
{
	const Real h = 1./m_nSize[0];
	const int nX = m_nBigSize[0]/2+1;
	const int nY = m_nBigSize[1];
	
	Matrix2D<Real[2]>& PSI = *m_matFourierPsi;
	Matrix2D<Real[2]>& GF = *m_matFourierGF;
	
	InterfaceFFTW::executePlan("RHS: Real->Fourier");
	
	for(int dy=0; dy<nY; dy++)
	for(int dx=0; dx<nX; dx++)
	{
		const Real a[2] = {PSI(dx, dy)[0], PSI(dx, dy)[1] };
		const Real b[2] = {GF(dx, dy)[0], GF(dx, dy)[1] };
		
		PSI(dx, dy)[0] = (a[0]*b[0]-a[1]*b[1]);
		PSI(dx, dy)[1] = (a[1]*b[0]+a[0]*b[1]);
	}
	
      	InterfaceFFTW::executePlan("Psi: Fourier->Real");
	(*m_matPsi) *= h*h/(m_nBigSize[0]*m_nBigSize[1]);
}

