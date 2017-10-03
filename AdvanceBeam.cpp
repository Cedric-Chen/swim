#pragma once

template<int sizeX, int sizeY>
void AdvanceBeam::operator()(
    //input
    const Layer<int,int,1> & currPressure,
    const vector<double> & oldAngle,
    const vector<double> & oldVelAngle,
    const vector<double> & currAngle,
    const vector<double> & currVelAngle,
    const vector<double> & currT,
    
    //output
    const vector<double> & newAngle,
    const vector<double> & newVelAngle,
    const vector<double> & newT
){
    //input
    pCurrPressure = &currPressure;
    pOldAngle = &oldAngle;
    pOldVelAngle = &oldVelAngle;
    pCurrAngle = &currAngle;
    pCurrVelAngle = &currVelAngle;
    pCurrT = &currT;
    //output
    pNewAngle = &newAngle;
    pNewVelAngle = &newVelAngle;
    pNewT = &newT;

    Broyden(&ComputeF, *pCurrVelAngle, *pNewVelAngle);
    *pNewAngle = ----
}

template<int sizeX, int sizeY>
double AdvanceBeam::ComputeF(const vector<double> & estimatedVelAngle){
    const Layer<sizeX,sizeY,1> &currPressure = *pCurrPressure;
    const vector<double> &currAngle = *pCurrAngle;
    const vector<double> &currVelAngle = *pCurrVelAngle;
}

template<int sizeX, int sizeY>
double AdvanceBeam::ComputeG(const vector<double> & estimatedT){

}

template<int sizeX, int sizeY>
void AdvanceBeam::Broyden(
    double (*Func)(vector<double> &),
    const vector<double> &estimated,
    vector<double> &result
){

}
