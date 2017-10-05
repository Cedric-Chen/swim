#pragma once

//Environment
//#include "Environment.h"
#include "LayerToVTK.h"

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
    
    mu{1.0},
    eta{1.0},

    pOldAngle{NULL},
    pOldVelAngle{NULL},
    pCurrAngle{NULL},
    pCurrVelAngle{NULL},
    pNewAngle{NULL},
    pNewVelAngle{NULL},
    oldDt{},
    currDt{},

    pCurrT{NULL},
    pNewT{NULL},
    pPressure{NULL},

    pCharFunc{NULL},
    pVelBeam{NULL},

    dw1ds{arma::vec{}},
    dw1ds_2{arma::vec{}},
    dw1ds_3{arma::vec{}},
    w2{arma::vec{}},
    D{arma::mat{}},
    I_adddown{arma::mat{}}
{
    mLoc = matrix(mNumofElem+1, vector<double>(2)); 
    mVel = matrix(mNumofElem+1, vector<double>(2));
    mAngle = vector<double>(mNumofElem);
    mVelAngle = vector<double>(mNumofElem);
    mDeltaFunc = Layer<sizeX, sizeY, 1>{};
    mPressureOnBody = vector<double>(mNumofElem);
    cout <<"constructor: " <<mLocHead[0]<<", " <<mLocHead[1]<<endl;
}

template<int sizeX,int sizeY>
void Beam<sizeX,sizeY>::advance(
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
    )
{

    this->oldDt = oldDt;
    this->currDt = currDt;
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

    vector<double> estimatedVelAngle
        {(*pCurrVelAngle).begin()+1,(*pCurrVelAngle).end()};
    broyden(
        'F', 
        estimatedVelAngle,
        (*pNewVelAngle).begin()+1,
        (*pNewVelAngle).end()
    );

    // from newCurrVelAngle to newCurrAngle
    const arma::vec w2{
         vector<double>{
            (*pNewAngle).begin()+1,
            (*pNewAngle).end()
            }
        };
    const arma::vec w1_0{
        vector<double>{
            (*pCurrAngle).begin()+1,
            (*pCurrAngle).end()
            }
        };
    arma::vec w1 = w1_0+ currDt*w2;
    w1.insert_rows(0,1,0);
    (*pNewAngle) = arma::conv_to<vector<double>>::from(w1);
}

template<int sizeX,int sizeY>
arma::vec Beam<sizeX,sizeY>::compute_F(const arma::vec & estimatedVelAngle){
    auto dim = estimatedVelAngle.size()+1;
    
    // w1 - estiamted angle, (dim-1)*1
    // w2 - estimated velAngle, (dim-1)*1
    w2 = estimatedVelAngle;
    const arma::vec w2_0{
        vector<double>{
            (*pCurrVelAngle).begin()+1,
            (*pCurrVelAngle).end()
            }
        };
    const arma::vec w2_minus1{
        vector<double>{
            (*pOldVelAngle).begin()+1,
            (*pOldVelAngle).end()
            }
        };
    const arma::vec w1_0{
        vector<double>{
            (*pCurrAngle).begin()+1,
            (*pCurrAngle).end()
            }
        };
    const arma::vec w1 = w1_0+ currDt*w2;

    // solve pressure at estimated geometry
    update(
        arma::conv_to<vector<double>>::from(w1),
        arma::conv_to<vector<double>>::from(w2)
    );
    const arma::vec pressureOnBody{mPressureOnBody};

    // build transforamtion and derivative matrix
    arma::mat I{dim, dim, arma::fill::eye};
    arma::mat I_addtop = I.cols(1,dim-1);
    I_adddown = I.cols(0,dim-2);
    const double H = (double)(*pPressure).getH0(); 
    D = -arma::diagmat(arma::vec{dim,1,arma::fill::ones}) +
        arma::diagmat(arma::vec{dim-1,1,arma::fill::ones}, 1);
    D = D.rows(0,dim-2) / H;

    dw1ds = (I_adddown * D) * (I_addtop * w1);
    dw1ds_2 = (I_adddown * D) * dw1ds;
    dw1ds_3 = D * dw1ds_2;
    dw1ds_3.insert_rows(0,1,0);

    // T - longitual force, dim * 1
    broyden(
        'G',
        vector<double>{(*pCurrT).begin(),(*pCurrT).end()-1},
        (*pNewT).begin(),
        (*pNewT).end()
    );
    const arma::vec T{(*pNewT)};

    dw1ds_3[0] = T[0] * dw1ds[0] - pressureOnBody[0];
    arma::vec dw1ds_4 = D * dw1ds_3; // (dim-1)*1
    arma::vec dpds = D * pressureOnBody; // (dim-1)*1
    arma::vec dTds = D * T; // (dim-1)*1

    arma::vec dw2dt = - dpds - eta*dw1ds_4 -
        (T.rows(1,dim-1) + eta * dw1ds.rows(1,dim-1) % dw1ds.rows(1,dim-1)) % 
        dw1ds_2.rows(1,dim-1) +
        2 * dTds % dw1ds.rows(1,dim-1);
    dw2dt /= mu;

    double xi = currDt / oldDt;
    arma::vec F = -w2 + pow(1+xi,2)/(1+2*xi) * w2_0 - 
        pow(xi,2)/(1+2*xi) * w2_minus1 + 
        (1+xi)/(1+2*xi) * currDt * dw2dt;
    F.print();
    return F;
}

template<int sizeX,int sizeY>
arma::vec Beam<sizeX,sizeY>::compute_G(const arma::vec & estimatedT){
    auto dim = estimatedT.size()+1;

    const arma::vec pressureOnBody{mPressureOnBody};
    pressureOnBody.print();

    const arma::vec & Tx = estimatedT;
    dw1ds_3[0] = Tx[0] * dw1ds[0] - pressureOnBody[0];

    arma::vec dTds{dim,1,arma::fill::ones};
    dTds.rows(1,dim-1) = D * I_adddown * Tx;
    dTds[0] = - dw1ds[0] * dw1ds_2[0];
    arma::vec dTds_2{D * dTds};
    
    arma::vec G = dTds_2 - 
        Tx % dw1ds.rows(0,dim-2) % dw1ds.rows(0,dim-2) +
        pressureOnBody.rows(0,dim-2) % dw1ds.rows(0,dim-2) +
        2*eta*dw1ds.rows(0,dim-2) % dw1ds_3.rows(0,dim-2) +
        eta * dw1ds_2.rows(0,dim-2) % dw1ds_2.rows(0,dim-2) +
        mu * w2 % w2;
    G.print();
    exit(0);
    return G;
}

template<int sizeX,int sizeY>
arma::vec Beam<sizeX,sizeY>::compute_test(const arma::vec & estimated){
    arma::mat res{arma::size(estimated),arma::fill::zeros};
    res[0] = 2*estimated[0] + 2*estimated[1] - estimated[2] -1;
    res[1] = 2*estimated[0] - 2*estimated[1] + 4*estimated[2] +2;
    res[2] = -estimated[0] + 0.5*estimated[1] - estimated[2];
    cout <<res[2]<<endl;
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
    vector<double>::iterator result_start,
    vector<double>::iterator result_end
//    vector<double> &result
){
    const double eps1 = 1e-8;
    const double eps2 = 1e-8;
    int dim = estimated.size();
    arma::vec x0{estimated};
    cout <<"broyden: "<<choose_function<<endl;
    arma::vec y0{compute(choose_function, x0)};
    
    arma::mat S(dim,dim,arma::fill::eye);
    arma::vec err(arma::size(x0),arma::fill::ones);

    double ferr{};
    arma::vec x{};
    arma::vec y{};

    x0.print();
    y0.print();

    bool state = true;
    while(state){
        arma::vec d = -solve(S, y0);
        x = x0 + d;
        y = compute(choose_function,x);
        S = S + ((y-y0) - S*d) *d.t() / arma::as_scalar((d.t()*d));
        double temp = arma::as_scalar(sqrt( (x-x0).t() * (x-x0)));
        err = temp - eps1 * (1 + abs(x));
        ferr = sqrt(arma::as_scalar(y.t() * y));

//        d.print();
//        S.print();
//        err.print();
//        cout <<ferr<<endl;
        x.print();
        y.print();

        x0 = x;
        y0 = y;

        state = true;
        for(arma::vec::const_row_iterator i = err.begin();
            i != err.end(); i++){
            if((*i)<=0){
                state = false;
                break;
            }
        }
    }

    if(ferr >= eps2){
        cout<<"Broyden fails; ferr>= eps2";
        exit(0);
    }

    int j{0};
    for(auto i=result_start; i != result_end; i++){
        *i = x[j];
        j++;
    }
}

template<int sizeX,int sizeY>
void Beam<sizeX,sizeY>::update(
        vector<double> angle,
        vector<double> velAngle
        ){

	vector<double> mAngle = angle;
	vector<double> mVelAngle = velAngle;

	update_loc();
	update_vel();
	update_charFunc();
	update_deltaFunc();
	update_velBeam();
        update_pressure();
        dumpLayer2VTK(0,"charFunc",*pCharFunc);
        dumpLayer2VTK(0,"deltaFunc",mDeltaFunc);
        dumpLayer2VTK(0,"velBeam",*pVelBeam);
        dumpLayer2VTK(0,"pressure",*pPressure);
}

template<int sizeX,int sizeY>
void Beam<sizeX,sizeY>::update_pressure(){
    Layer<sizeX,sizeY,1> & deltaFunc = mDeltaFunc;
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
    cout <<"locHead: "<<mLocHead[0]<<", " <<mLocHead[1]<<endl;
    mLoc[0] = mLocHead;
    for(int i=0;i<mNumofElem;i++){
        mLoc[i+1][0] = mLoc[i][0] + mLengthofElem*cos(mAngle[i]);
        mLoc[i+1][1] = mLoc[i][1] + mLengthofElem*sin(mAngle[i]);
    }

    //debug
    cout <<"location: ";
    for(auto i = mLoc.begin(); i != mLoc.end(); i++){
        cout <<" (";
        for(auto j = (*i).begin(); j != (*i).end(); j++)
            cout <<" " <<*j;
        cout <<")";
    }
    cout <<"\n";
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
