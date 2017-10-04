#pragma once

#include <vector>
#include <cstdlib>
#include "Layer.h"
using namespace std;

template<int sizeX, int sizeY>
class AdvanceBeam{
public:
    AdvanceBeam(){}

    void operator()(
        const Layer<sizeX,sizeY,1> & currPressure,
        const vector<double> & oldAngle, 
        const vector<double> & oldVelAngle,
        const vector<double> & currAngle, 
        const vector<double> & currVelAngle,
        const vector<double> & currT,
        // output
        vector<double> & newAngle, vector<double> & newVelAngle,
        vector<double> & newT
    );

private:
    const Layer<sizeX,sizeY,1> * pCurrPressure;
    const vector<double> * pOldAngle;
    const vector<double> * pOldVelAngle;
    const vector<double> * pCurrAngle;
    const vector<double> * pCurrVelAngle;
    const vector<double> * pCurrT;

    vector<double> * pNewAngle;
    vector<double> * pNewVelAngle;
    vector<double> * pNewT;
   
    double ComputeF(const vector<double> & estimatedVelAngle);
    double ComputeG(const vector<double> & estimatedT);
    void Broyden(
        double (*Func)(vector<double> &)
        const vector<double> &estimated,
        vector<double> &result
    );
};

#include "AdavanceBeam.cpp"
