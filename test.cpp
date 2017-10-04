#include "Layer.h"
#include "LayerToVTK.h"
#include "Beam.h"
#include "ComputePressure_TBB.h"
#include "ComputeViscousForces_TBB.h"

#include <vector>
#include <iostream>
using namespace std;

void testBeam();
//void testCP();
//void testCVF();

const int grid_x = 1000;
const int grid_y = 100;
double width = 0.02;
double length = 0.5;
double epsilon = 0.002;
int numofElem = 51;
vector<double> locHead{0,0.05};
vector<double> velHead{0,0};

Layer<grid_x,grid_y,1> pressure = Layer<grid_x,grid_y,1>{};
vector<double> oldAngle = vector<double>(numofElem);
vector<double> oldVelAngle = vector<double>(numofElem);
vector<double> currAngle = vector<double>(numofElem);
vector<double> currVelAngle = vector<double>(numofElem);
vector<double> newAngle = vector<double>(numofElem);
vector<double> newVelAngle = vector<double>(numofElem);
vector<double> currT = vector<double>(numofElem);
vector<double> newT = vector<double>(numofElem);
Layer<grid_x,grid_y,2> velBeam = Layer<grid_x,grid_y,2>{};
Layer<grid_x,grid_y,1> charFunc = Layer<grid_x,grid_y,1>{};

int main()
{
    testBeam();
//    testCP();
//    testCVF();
}

void testBeam(){
    Beam<grid_x, grid_y>beam = Beam<grid_x, grid_y>
        {width,length,locHead,velHead,epsilon,numofElem};
    beam.advance(pressure,oldAngle,oldVelAngle,currAngle,
        currVelAngle,newAngle,newVelAngle,currT,newT,
        velBeam,charFunc);
}

//void testBeam(){
//    double theta1 = 0.3;
//    double theta2 = -0.4;
//    vector<double> angle{0};
//    for(int i = 1; i<numofElem; i++)
//        if(i%2 == 1)
//            angle.push_back(theta1);
//        else
//            angle.push_back(theta2);
//
//    double speed = -0.1;
//    vector<double> velAngle{0};
//    for(int i = 1; i<numofElem; i++)
//        velAngle.push_back( (i%2*2-1) * speed);
//
//    Layer<grid_x,grid_y,1> charFunc_old;
//    Layer<grid_x,grid_y,1> deltaFunc_old;
//    Layer<grid_x,grid_y,2> velBeam_old;
//    charFunc_old = 0;
//    deltaFunc_old = 0;
//    velBeam_old = 0;
//
//    Layer<grid_x,grid_y,1> charFunc;
//    Layer<grid_x,grid_y,1> deltaFunc;
//    Layer<grid_x,grid_y,2> velBeam;
//    charFunc = 0;
//    deltaFunc = 0;
//    velBeam = 0;
//
//    Beam<grid_x,grid_y> beam{width,length,epsilon,numofElem};
//    beam.update(locHead,velHead,angle,velAngle,
//            charFunc,deltaFunc,velBeam);
//
//    cout << "if charFunc changes: " << !(charFunc_old == charFunc) <<"\n";
//    cout << "if deltaFunc changes: " << !(deltaFunc_old == deltaFunc) <<"\n";;
//    cout << "if velBeam changes: " << !(velBeam_old == velBeam) <<"\n";;
//    dumpLayer2VTK(0,"charFunc",charFunc);
//    dumpLayer2VTK(0,"deltaFunc",deltaFunc);
//    dumpLayer2VTK(0,"vel",velBeam);
//
//    //test update_pressure function
//    Layer<grid_x,grid_y,1> pressure;
//    vector<double> pressureOnBody(numofElem);
//    for(int i=0;i<grid_x;i++)
//    for(int j=0;j<grid_y;j++)
//        pressure(i,j,0) = i;
//    beam.update_pressure(pressure,pressureOnBody);
//
//    dumpLayer2VTK(0,"pressure_beam_input",pressure);
//    cout <<"pressure_on_body: ";
//    for(auto i = pressureOnBody.begin(); i != pressureOnBody.end(); i++){
//        cout << " " <<*i <<",";
//    }
//    cout <<"\n";
//}
//
//void testCP()
//{
//    const int grid = 50;
//    ComputePressure_TBB<grid,grid> computepressure;
//    Layer<grid,grid,2> velold;
//    Layer<grid,grid,2> velnew;
//    Layer<grid,grid,1> pressure_old;
//    Layer<grid,grid,1> pressure;
//
//    velold = 0;
//    for(int xx=0; xx<grid; xx++)
//        for(int yy=0; yy<grid; yy++){
//            velold(xx,yy,0) = -xx;
//            velold(xx,yy,1) = -yy;
//        }
//    velnew = 0;
////    for(int xx=0; xx<grid; xx++)
////        for(int yy=0; yy<grid; yy++){
////            velnew(xx,yy,0) = -xx;
////            velnew(xx,yy,1) = -yy;
////        }
//
//    pressure_old = 0;
//    pressure = 0;
//    double nu = 1;
//    double dt=1;
//    double g[2] ={0,0};
//
//    computepressure.setup(dt,g,nu,velold,velnew,pressure);
//    computepressure();
//
//    cout << "if pressure changes: " << !(pressure_old == pressure) <<"\n";
//    dumpLayer2VTK(0,"pressure",pressure);
//}
//
//void testCVF()
//{
//    const int grid = 50;
//    ComputeViscousForces_TBB<grid,grid> computeviscousforces;
//    Layer<grid,grid,2> vel;
//    Layer<grid,grid,2> viscousForces;
//    Layer<grid,grid,2> viscousForces_old;
//
//    vel = 0;
//    for(int xx=0; xx<grid; xx++)
//        for(int yy=0; yy<grid; yy++){
//            vel(xx,yy,0) = pow(xx,3);
//            vel(xx,yy,1) = pow(yy,3);
//        }
//
//    viscousForces_old = 0;
//    viscousForces = 0;
//    double nu = 1;
//
//    computeviscousforces.setup(vel,viscousForces,nu);
//    computeviscousforces();
//
//    cout << "if viscous forces changes: " 
//        << !(viscousForces_old == viscousForces) <<"\n";
//    dumpLayer2VTK(0,"viscous_forces",viscousForces);
//}
