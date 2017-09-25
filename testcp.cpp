#include "ComputeGradPressure_TBB.h"
#include "Layer.h"

int main(){
    ComputeGradPressure_TBB<20,20> computegradpressure;
    Layer<20,20,2> velold;
    Layer<20,20,2> velnew;
    Layer<20,20,2> gradPressure;
    velold = 0;
    velnew = 0;
    gradPressure = 0;
    double nu = 1;
    double dt=1;
    double g[2] ={0,1};

    computegradpressure.setup(dt,g,nu,velold,velnew,gradPressure);
    return 1;
}
