#pragma once

#include <vector>
#include <cstdlib>
#include "Layer.h"
using namespace std;

template <int sizeX,int sizeY>
class Beam{
public:
	Beam(double width,double length,double epsilon,int numofElem)
		:mWidth{width},
        mLength{length},
        mEpsilon{epsilon},
		mNumofElem{numofElem},
		pCharFunc{NULL},
        pDeltaFunc{NULL},
        pVelBeam{NULL} 
		{
            mLengthofElem = mLength/mNumofElem;
        }
	~Beam() {}

    /// Update the position of beam.
    /// 
    /// This function update charFunc, deltaFunc, and velbeam 
    /// according to new position, velocity and shape, 
    /// which include location and velocity of head, 
    /// angle and angle velocity of the beam.
	void update(
        //input
        vector<double> locHead,
        vector<double> velHead,
		vector<double> angle,
        vector<double> velAngle,
        //output
		Layer<sizeX,sizeY,1>& charFunc, 
        Layer<sizeX,sizeY,1>& deltaFunc,
		Layer<sizeX,sizeY,2>& velBeam
        );

    void update_pressure(
        //input
        Layer<sizeX,sizeY,1>& pressure,
        //output
        vector<double>& pressureOnBody
        );

private:
	const double mWidth; //width of the beam
	const double mLength; //length of the beam
	const double mEpsilon; //thickness of buffer area near interface

	int mNumofElem;  //number of elements
	double mLengthofElem; // = mLength/mNumofElement;

	vector<double> mLocHead = vector<double>(2); //location of head
	vector<double> mVelHead = vector<double>(2); //veloctiy of head

	typedef vector< vector<double> > matrix;
	//do not include location and of velocity of head.
	matrix mLoc = matrix(mNumofElem+1, vector<double>(2)); 
	matrix mVel = matrix(mNumofElem+1, vector<double>(2));
	vector<double> mAngle = vector<double>(mNumofElem);
	vector<double> mVelAngle = vector<double>(mNumofElem); //couterclockwise is positive

	Layer<sizeX, sizeY, 1> * pCharFunc;
	Layer<sizeX, sizeY, 1> * pDeltaFunc;
	Layer<sizeX, sizeY, 2> * pVelBeam;

	//update helper function
	void update_loc();
	void update_vel();
	void update_charFunc();
	void update_deltaFunc();
	void update_velBeam();

	bool isBetweenVerticalLine(double x,double y,
			vector<double> p1,vector<double> p2);
	double getDistance(double x, double y, 
			vector<double> p1, vector<double> p2);

};

#include "Beam.cpp"
