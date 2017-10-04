#pragma once

//Environment
//#include "Environment.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <armadillo>

using namespace std;

template<int sizeX,int sizeY>
Beam<sizeX, sizeY>::Beam(double width,double length,
    vector<double> locHead, vector<double> velHead,
    double epsilon,int numofElem)
    :mWidth{width},
    mLength{length},
    mLocHead{locHead},
    mVelHead{velHead},
    mEpsilon{epsilon},
    mNumofElem{numofElem},
    mLengthofElem{mLength/mNumofElem},

    pOldAngle{NULL},
    pOldVelAngle{NULL},
    pCurrAngle{NULL},
    pCurrVelAngle{NULL},
    pNewAngle{NULL},
    pNewVelAngle{NULL},

    pCurrT{NULL},
    pNewT{NULL},
    pPressure{NULL},

    pCharFunc{NULL},
    pVelBeam{NULL} 
{
    mLoc = matrix(mNumofElem+1, vector<double>(2)); 
    mVel = matrix(mNumofElem+1, vector<double>(2));
    mAngle = vector<double>(mNumofElem);
    mVelAngle = vector<double>(mNumofElem);
    mDeltaFunc = Layer<sizeX, sizeY, 1>{};
    mPressureOnBody = vector<double>(mNumofElem);
}

template<int sizeX,int sizeY>
void Beam<sizeX,sizeY>::advance(
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
    )
{
    pOldAngle = &oldAngle;
    pOldVelAngle = &oldVelAngle;
    pCurrAngle = &currAngle;
    pCurrVelAngle = &currVelAngle;
    pNewAngle = &newAngle;
    pNewVelAngle = &newVelAngle;

    pCurrT = &currT;
    pNewT = &newT;
    pPressure = &pressure;

    pCharFunc = &charFunc;
    pVelBeam = &velBeam;
    
    vector<double> est{0,0};
    vector<double> res{0,0};
    broyden('T',est,res);
    cout <<"result: " <<res[0] <<", " <<res[1] <<endl; 
//    broyden('F', *pCurrVelAngle, *pNewVelAngle);
//    *pNewAngle = -----;
}

template<int sizeX,int sizeY>
arma::vec Beam<sizeX,sizeY>::compute_F(const arma::vec & estimatedVelAngle){
    return arma::vec{};
}

template<int sizeX,int sizeY>
arma::vec Beam<sizeX,sizeY>::compute_G(const arma::vec & estimatedT){
    return arma::vec{};
}

template<int sizeX,int sizeY>
arma::vec Beam<sizeX,sizeY>::compute_test(const arma::vec & estimated){
    arma::mat res{arma::size(estimated),arma::fill::zeros};
    res[0] = estimated[0] + estimated[1] - 3;
    res[1] = estimated[0]^2 + estimated[1]^2 - 9;
    return res;
}

template<int sizeX,int sizeY>
arma::vec Beam<sizeX,sizeY>::compute(
    char choose_function, const arma::vec & estimated)
{
    switch(choose_function){
        case 'F':
            return compute_F(estimated);
        case 'G':
            return compute_G(estimated);
        case 'T':
            return compute_test(estimated);
        default:
            cout<<"choose_function is wrong";
            exit(0);
    }
}


template<int sizeX,int sizeY>
void Beam<sizeX,sizeY>::broyden(
    char choose_function,
    const vector<double> &estimated,
    vector<double> &result
){
    const double eps1 = 1e-8;
    const double eps2 = 1e-8;
    int dim = estimated.size();
    arma::vec x0{estimated};
    arma::vec y0{compute(choose_function, x0)};
    
    arma::mat S(dim,dim,arma::fill::eye);
    arma::mat err(arma::size(x0),arma::fill::ones);

    double ferr{};
    arma::vec x{};
    arma::vec y{};
    while(err >0){
        arma::vec d = -solve(S, y0);
        x = x0 + d;
        y = compute(choose_function,x);
        S = S + ((y-y0) - S*d) *d.t() / (d.t()*d);
        double temp = sqrt( (x-x0).t() * (x-x0));
        err = temp - eps1 * (1 + abs(x));
        ferr = sqrt(y.t() * y);
        x0 = x;
        y0 = y;
    }

    if(ferr >= eps2){
        cout<<"Broyden fails; ferr>= eps2";
        exit(0);
    }

    for(int i=0;i<dim;i++){
        result[i] = x[i];
    }
}

template<int sizeX,int sizeY>
void Beam<sizeX,sizeY>::update(
        vector<double> angle,
        vector<double> velAngle
        ){

	vector<double> mAngle = angle;
	vector<double> mVelAngle = velAngle;
//    //debug
//    cout <<"angle: ";
//    for(auto i = mAngle.begin(); i != mAngle.end(); i++)
//        cout <<" " <<*i;
//    cout <<"\n";
//    cout <<"velAngale: ";
//    for(auto i = mVelAngle.begin(); i != mVelAngle.end(); i++)
//        cout <<" " <<*i;
//    cout <<"\n";

	update_loc();
	update_vel();
	update_charFunc();
	update_deltaFunc();
	update_velBeam();
        update_pressure();
}

template<int sizeX,int sizeY>
void Beam<sizeX,sizeY>::update_pressure(){
    const Layer<sizeX,sizeY,1> & deltaFunc = mDeltaFunc;
    const Layer<sizeX,sizeY,1> & pressure = *pPressure;
    vector<double> & pressureOnBody = mPressureOnBody;

    const double H = (double)deltaFunc.getH0(); 
    const int sX{0}; const int eX{sizeX};
    const int sY{0}; const int eY{sizeY};
    double xx{0};
    double yy{0};

    for(int i=0;i<mNumofElem;i++){
        vector<double> p1{mLoc[i]};
        vector<double> p2{mLoc[i+1]};

        double pressure_up{0};
        double pressure_down{0};
        double ratio_up{0};
        double ratio_down{0};

        for(int iy=sY; iy<eY; iy++)
        for(int ix=sX; ix<eX; ix++){
                xx = ((double)ix+0.5)*H;
                yy = ((double)iy+0.5)*H;

        if(isBetweenVerticalLine(xx,yy,p1,p2)){
            if(deltaFunc(ix,iy,0) != 0){
                //whether point(xx,yy) is above or below line p1-p2
                //I assume that p2[0] > p1[0] all the time. 
                //line p1-p2 is characterized as Ax + y + C = 0
                double A = - ( (p2[1]-p1[1]) / (p2[0]-p1[0]) );
                double C = -A*p1[0] - p1[1];
                bool above_or_not = false;
                if(A*xx + yy + C >= 0)
                    above_or_not = true;

                if(above_or_not){
                    pressure_up += pressure(ix,iy,0)*deltaFunc(ix,iy,0);
                    ratio_up += deltaFunc(ix,iy,0);
                }
                else{
                    pressure_down += pressure(ix,iy,0)*deltaFunc(ix,iy,0);
                    ratio_down += deltaFunc(ix,iy,0);
                }
            }
        }
        }
        pressure_up /= ratio_up;
        pressure_down /= ratio_down;

        pressureOnBody[i] = pressure_down - pressure_up;
    }	
}

template<int sizeX,int sizeY>
void Beam<sizeX,sizeY>::update_loc(){
    mLoc[0] = mLocHead;
    for(int i=0;i<mNumofElem;i++){
        mLoc[i+1][0] = mLoc[i][0] + mLengthofElem*cos(mAngle[i]);
        mLoc[i+1][1] = mLoc[i][1] + mLengthofElem*sin(mAngle[i]);
    }

//    //debug
//    cout <<"location: ";
//    for(auto i = mLoc.begin(); i != mLoc.end(); i++){
//        cout <<" (";
//        for(auto j = (*i).begin(); j != (*i).end(); j++)
//            cout <<" " <<*j;
//        cout <<")";
//    }
//    cout <<"\n";
}

template<int sizeX,int sizeY>
void Beam<sizeX,sizeY>::update_vel(){
//    cout <<"length of elem = " <<mLengthofElem <<endl;
    mVel[0] = mVelHead;
    for(int i=0;i<mNumofElem;i++){
            mVel[i+1][0] = mVel[i][0] - 
                    mLengthofElem*mVelAngle[i]*sin(mAngle[i]);
            mVel[i+1][1] = mVel[i][1] + 
                    mLengthofElem*mVelAngle[i]*cos(mAngle[i]);
    }

    //debug
//    cout <<"velocity: ";
//    for(auto i = mVel.begin(); i != mVel.end(); i++){
//        cout <<" (";
//        for(auto j = (*i).begin(); j != (*i).end(); j++)
//            cout <<" " <<*j;
//        cout <<")";
//    }
//    cout <<"\n";

}

template<int sizeX,int sizeY>
void Beam<sizeX,sizeY>::update_charFunc(){
    Layer<sizeX,sizeY,1> & charFunc = *pCharFunc;
    charFunc = 0;

    const double width = mWidth;
    const double epsilon = mEpsilon;
    const double H = (double)charFunc.getH0(); 
    const int sX{0}; const int eX{sizeX};
    const int sY{0}; const int eY{sizeY};
    double xx{0};
    double yy{0};
    double dist{0};

    for(int i=0;i<mNumofElem;i++){
        vector<double> p1{mLoc[i]};
        vector<double> p2{mLoc[i+1]};
//        cout <<p1[0] <<" ," <<p1[1] <<"\n";

        for(int iy=sY; iy<eY; iy++)
        for(int ix=sX; ix<eX; ix++){
            xx = ((double)ix+0.5)*H;
            yy = ((double)iy+0.5)*H;

            if(isBetweenVerticalLine(xx,yy,p1,p2)){
                //distance with the line pass p1 and p2
                dist=getDistance(xx,yy,p1,p2);		
                if(dist<width/2-epsilon){
                    charFunc(ix,iy,0)=1;
                }
                else if(dist>width/2+epsilon){
                    ;
                    //charFunc(ix,iy,0)=
                        //std::max(charFunc(ix,iy,0),0);
                }
                else{
                    charFunc(ix,iy,0)=
                        std::max(charFunc(ix,iy,0),
                        0.5*(1.0 + (width/2-dist)/epsilon + 
                        1.0/M_PI*sin(M_PI*(width/2-dist)/epsilon)));
                }
            }
        }
    }
}

template<int sizeX,int sizeY>
void Beam<sizeX,sizeY>::update_deltaFunc(){
    Layer<sizeX,sizeY,1> & deltaFunc = mDeltaFunc;
    deltaFunc = 0;

    const double width = mWidth;
    const double epsilon = mEpsilon;
    const double H = (double)deltaFunc.getH0(); 
    const int sX{0}; const int eX{sizeX};
    const int sY{0}; const int eY{sizeY};
    double xx{0};
    double yy{0};
    double dist{0};

    for(int i=0;i<mNumofElem;i++){
        vector<double> p1{mLoc[i]};
        vector<double> p2{mLoc[i+1]};

        for(int iy=sY; iy<eY; iy++)
        for(int ix=sX; ix<eX; ix++){
            xx = ((double)ix+0.5)*H;
            yy = ((double)iy+0.5)*H;

            if(isBetweenVerticalLine(xx,yy,p1,p2)){
                //distance with the line pass p1 and p2		
                dist=getDistance(xx,yy,p1,p2);
                //distance with interface
                double r = abs(dist-width/2);				
                double r_p = r/epsilon;

                if(r_p<0.5){
                    deltaFunc(ix,iy,0)=
                        std::max(deltaFunc(ix,iy,0),
                        3.0/8.0 + M_PI/32.0 - r_p*r_p/4.0 );
                }
                else if(r_p<1.5){
                    deltaFunc(ix,iy,0)=
                        std::max(deltaFunc(ix,iy,0),
                        (0.25 + (1-r_p)/8.0*sqrt(-2.0+8.0*r_p-4*r_p*r_p) 
                        - 0.125*asin(sqrt(2)*(r_p-1))));
                }
                else if(r_p<2.5){
                    deltaFunc(ix,iy,0)=
                        std::max(deltaFunc(ix,iy,0),
                        (17.0/16.0 - M_PI/64.0 
                        - 3.0*r_p/4.0 + r_p*r_p/8.0 
                        + (r_p-2)/16.0*sqrt(16*r_p-4*r_p*r_p-14) 
                        + 1.0/16.0*asin(sqrt(2)*(r_p-2))));
                }
                else{
                    ;
                }
            }
        }
    }
}	


		
template<int sizeX,int sizeY>
void Beam<sizeX,sizeY>::update_velBeam(){
    Layer<sizeX,sizeY,2> & velBeam = *pVelBeam;
    velBeam = 0;
    Layer<sizeX,sizeY,1> & charFunc = *pCharFunc;

    const double width = mWidth;
    const double epsilon = mEpsilon;
    // the length of the domain is 1 by default.
    const double H = (double)charFunc.getH0(); 
    const int sX{0}; const int eX{sizeX};
    const int sY{0}; const int eY{sizeY};
    double xx{0};
    double yy{0};
    double dist{0};

    for(int i=0;i<mNumofElem;i++){
        vector<double> p1{mLoc[i]};
        vector<double> p2{mLoc[i+1]};

        for(int iy=sY; iy<eY; iy++)
        for(int ix=sX; ix<eX; ix++){
            xx = ((double)ix+0.5)*H;
            yy = ((double)iy+0.5)*H;

            if(isBetweenVerticalLine(xx,yy,p1,p2)){
                //distance with the line pass p1 and p2
                dist=getDistance(xx,yy,p1,p2);		
                if(dist<width/2+epsilon){
                    // square of distance between (xx,yy) and p1
                    double dist_hyp_square = pow(xx-p1[0],2) + pow(yy-p1[1],2);
                    // distance between foot point of (xx,yy) on p1-p2 and p1
                    double dist_leg = sqrt(dist_hyp_square - pow(dist,2));
                    // distance between p1 and p2
                    double dist_p1_p2 = 
                        sqrt( pow(p2[0]-p1[0],2) + pow(p2[1]-p1[1],2) );
                    double ratio = dist_leg/dist_p1_p2;

                    velBeam(ix,iy,0) = charFunc(ix,iy) * (mVel[i][0] 
                        + (mVel[i+1][0]-mVel[i][0]) * ratio);
                    velBeam(ix,iy,1) = charFunc(ix,iy) * (mVel[i][1] 
                        + (mVel[i+1][1]-mVel[i][1]) * ratio);
                }
            }
        }
    }
}

template<int sizeX,int sizeY>
bool Beam<sizeX,sizeY>::isBetweenVerticalLine(double x, double y, 
    vector<double> p1, vector<double> p2){
    double condition=0;

    if(p2[1]-p1[1] != 0 && p2[0]-p1[0] != 0){
        double slope=(p2[1] - p1[1])/(p2[0] - p1[0]);
        double slope_v = -1.0/slope; //slope of the vertical line
		
        condition = (slope_v*(x-p1[0])-(y-p1[1])) 
            * (slope_v*(x-p2[0])-(y-p2[1]));
    }
    else if(p2[1] - p1[1] == 0)
        condition = (x-p1[0]) * (x-p2[0]);
    else
        condition = (y-p1[1])*(y-p2[1]);

    if(condition>0)
        return 0;
    else
        return 1;
}

/// Returns the absolute distance between a point and a line
template<int sizeX,int sizeY>
double Beam<sizeX, sizeY>::getDistance(double x, double y, 
    vector<double> p1, vector<double> p2){
    if(p2[0]-p1[0] != 0){
        double slope=(p2[1] - p1[1])/(p2[0] - p1[0]);
        return abs((slope*(x-p1[0])-(y-p1[1]))/sqrt(slope*slope+1));
    }
    else
        return abs(x-p1[0]);
}
