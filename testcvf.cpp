#include "ComputeViscousForces_TBB.h"
#include "Layer.h"

int main(){
    ComputeViscousForces_TBB<20,20> computeviscousforces;
    Layer<20,20,2> vel;
    Layer<20,20,2> forces;
    vel = 0;
    forces = 0;
    double nu = 1;

    computeviscousforces.setup(vel,forces,nu);
    return 1;
}
