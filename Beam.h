#pragma once

#include <vector>
#include <armadillo>
#include <cstdlib>
#include "Layer.h"
using namespace std;

template <int sizeX,int sizeY>
class Beam{
public:
    Beam(double width,double length,
        vector<double> locHead, vector<double> velHead,
        double epsilon,int numofElem);
    ~Beam() {}

    void advance(
        double oldDt,
        double currDt,
        const Layer<sizeX,sizeY,1> & pressure,
        const vector<double> & oldAngle,
        const vector<double> & oldVelAngle,
        const vector<double> & currAngle,
        const vector<double> & currVelAngle,
        vector<double> & newAngle,
        vector<double> & newVelAngle,
        const vector<double> & currT,
        vector<double> & newT,
        Layer<sizeX,sizeY,2> & velBeam,
        Layer<sizeX,sizeY,1> & charFunc
    );

private:
    // geometry
    const double mWidth; //width of the beam
    const double mLength; //length of the beam
    const vector<double> mLocHead; //location of head
    const vector<double> mVelHead; //veloctiy of head
    // physical properties of beam
    double mu;
    double eta;
    // computational setting
    const double mEpsilon; //thickness of buffer area near interface
    const int mNumofElem;  //number of elements
    const double mLengthofElem; // = mLength/mNumofElement;

    typedef vector< vector<double> > matrix;
    //couterclockwise is positive
    const vector<double> * pOldAngle;
    const vector<double> * pOldVelAngle; 
    const vector<double> * pCurrAngle;
    const vector<double> * pCurrVelAngle; 
    const vector<double> * pCurrT; 
    const Layer<sizeX, sizeY, 1> * pPressure;
    double oldDt;
    double currDt;

    // internal variable
    //do not include location and of velocity of head.
    matrix mLoc = matrix(mNumofElem+1, vector<double>(2)); 
    matrix mVel = matrix(mNumofElem+1, vector<double>(2));
    vector<double> mAngle = vector<double>(mNumofElem);
    vector<double> mVelAngle = vector<double>(mNumofElem);
    Layer<sizeX, sizeY, 1> mDeltaFunc;
    vector<double> mPressureOnBody;
    // internal shared calculation variable
    arma::vec dw1ds;
    arma::vec dw1ds_2;
    arma::vec dw1ds_3;
    arma::vec w2;
    arma::mat D;
    arma::mat I_adddown;

    //output
    vector<double> * pNewAngle;
    vector<double> * pNewVelAngle; 
    vector<double> * pNewT; 
    Layer<sizeX, sizeY, 1> * pCharFunc;
    Layer<sizeX, sizeY, 2> * pVelBeam;

    // deform helper function
    arma::vec compute_F(const arma::vec & estimatedVelAngle);
    arma::vec compute_G(const arma::vec & estimatedT);
    arma::vec compute_test(const arma::vec & estimated);
    arma::vec compute(char choose_function, const arma::vec &estimated);
    void broyden(
        char choose_function,
        const vector<double> &estimated,
        vector<double>::iterator result_begin,
        vector<double>::iterator result_end
//        vector<double> &result
    );

    //update helper function
    void update(
        //input
	vector<double> angle,
        vector<double> velAngle
        );
    void update_pressure();
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
