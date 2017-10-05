#pragma once 
#include <math.h>
#include <string>
#include <vector>
using namespace std;
#include <assert.h>
#include <cstdlib>
#include "Typedefs.h"

#ifdef _WIN32
#define ALIGN_ATTRIBUTE __declspec( align( 16 ) )
//
#else
#define ALIGN_ATTRIBUTE __attribute__((aligned(16))) 
#endif

template<int _sizeX, int _sizeY, int _nDim=2>
struct ALIGN_ATTRIBUTE Layer
{
	static const int sizeX = _sizeX;
	static const int sizeY = _sizeY;
	static const int nDim = _nDim;
	
	Real data[nDim][sizeY][sizeX];

	inline Real& operator ()(int ix = 0, int iy=0, int dim = 0)
	{
#ifndef NDEBUG
		assert(ix>=0 && ix<sizeX);
		assert(iy>=0 && iy<sizeY);
#endif
		
		return data[dim][iy][ix];
	}

	inline Real operator ()(int ix = 0, int iy=0, int dim = 0) const
	{
#ifndef NDEBUG
		assert(ix>=0 && ix<sizeX);
		assert(iy>=0 && iy<sizeY);
#endif
		
		return data[dim][iy][ix];
	}


	inline Real read(int ix = 0, int iy=0, int dim = 0) const
	{
#ifndef NDEBUG
		assert(ix>=0 && ix<sizeX);
		assert(iy>=0 && iy<sizeY);
#endif
		return data[dim][iy][ix];
	}

	const Layer& operator=(const Real val)
	{
		for(int idim = 0; idim<nDim; idim++)
		for(int iy = 0; iy<sizeY; iy++)
		for(int ix = 0; ix<sizeX; ix++)
			data[idim][iy][ix] = val;
		
		return *this;
	}

	const Layer& operator=(const Layer& l)
	{
		for(int idim = 0; idim<nDim; idim++)
			for(int iy = 0; iy<sizeY; iy++)
				for(int ix = 0; ix<sizeX; ix++)
					data[idim][iy][ix] = l.data[idim][iy][ix];
		
		return *this;
	}

	template<int dim>
	void clear(Real val)
	{
		for(int iy = 0; iy<sizeY; iy++)
		for(int ix = 0; ix<sizeX; ix++)
			data[dim][iy][ix] = val;
	}

	const vector<double> operator -(const Layer& layer)
	{
		vector<double> result;
		
		//compute linf distance
		{
			double LInf_diff = 0;
			for(int idim = 0; idim<nDim; idim++)
				for(int iy = 0; iy<sizeY; iy++)
					for(int ix = 0; ix<sizeX; ix++)
						LInf_diff = max(LInf_diff, (double)fabs(data[idim][iy][ix] - layer.data[idim][iy][ix]));
			
			result.push_back(LInf_diff);
		}
		
		//compute linf distance
		{
			double L2error = 0;
			for(int idim = 0; idim<nDim; idim++)
				for(int iy = 0; iy<sizeY; iy++)
					for(int ix = 0; ix<sizeX; ix++)
						L2error += pow(data[idim][iy][ix] - layer.data[idim][iy][ix], 2);
			
			result.push_back(sqrt((double)L2error/(sizeY*sizeX)));
		}
		
		return result;
	}

	const vector<double> computeNorms(int start_d=0, int end_d=nDim) const
	{
		vector<double> result;
		
		//compute linf distance
		{
			for(int idim = start_d; idim<end_d; idim++)
			{
				double LInf_diff = 0;
				for(int iy = 0; iy<sizeY; iy++)
					for(int ix = 0; ix<sizeX; ix++)
						LInf_diff = max(LInf_diff, (double)fabs(data[idim][iy][ix]));
				
				result.push_back(LInf_diff);
			}
			
			
		}
		
		//compute linf distance
		{
			for(int idim = start_d; idim<end_d; idim++)
			{
				double L2error = 0;
				
				for(int iy = 0; iy<sizeY; iy++)
					for(int ix = 0; ix<sizeX; ix++)
						L2error += pow(data[idim][iy][ix], 2);
				
				result.push_back(sqrtf(L2error/(sizeY*sizeX)));
			}
		}
		
		return result;
	}

	template<int iDim>
	Real * getPlane() 
	{
		return (Real*)&data[iDim][0][0];
	}

	static Real getH0() {return 1.0/sizeX;}
	static Real getH1() {return 1.0/sizeX;}

	void MatlabDelCazzo(string sFileName) const
	{
		FILE * f = fopen(sFileName.data(), "w");
		assert(f!=NULL);
		
		fprintf(f, "%d\n", sizeX);
		fprintf(f, "%d\n", sizeY);
		for (int iy=0; iy<sizeY; iy++)
		{
			for (int ix=0; ix<sizeX; ix++)
			{
				fprintf(f, "%e\t", data[0][iy][ix]);
			}
			fprintf(f, "\n");
		}
		
		fclose(f);
	}

	void serialize(string sFileName) const
	{
		vector<double> norms = computeNorms();
		
		FILE * f = fopen(sFileName.data(),"wb");
		fwrite(&norms[0], sizeof(double), 1, f);
		fwrite(&norms[1], sizeof(double), 1, f);
		fwrite(data, sizeof(Real), nDim*sizeX*sizeY, f);
		
		printf("serialized,  norm are : Linf norm  = %e, L2 norm = %e\n", norms[0], norms[1]);
		
		fclose(f);
	}

	void deserialize(string sFileName) 
	{
		FILE * f = fopen(sFileName.data(),"rb");
		if (f==NULL) abort();
		
		double n1=0, n2=0;
		fread(&n1, sizeof(double), 1, f);
		fread(&n2, sizeof(double), 1, f);
		
		fread(data, sizeof(Real), _nDim*_sizeX*_sizeY, f);
		fclose(f);
		
		vector<double>norms = computeNorms();
		
		printf("deserialized, checking norm errors (che palle): Linf norm error = %e, L2 norm error= %e\n", fabs(n1-norms[0]), fabs(n2-norms[1]));
			   
		if(n1 != norms[0]) abort();
		if(n2 != norms[1]) abort();
	}

    // added by Cedric
    bool operator ==(const Layer& layer)
    {
        bool isequal{true};

        for(int dim = 0; dim<nDim; dim++)
		for(int iy = 0; iy<sizeY; iy++)
		for(int ix = 0; ix<sizeX; ix++)
            if(data[dim][iy][ix] != layer.data[dim][iy][ix])
                isequal = false;

        return isequal;
    }
};
