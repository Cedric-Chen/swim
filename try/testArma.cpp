#include "armadillo"
#include <iostream>
#include <vector>

using namespace arma;

int main(){
    std::vector<double> a{1,2,3};
    vec b{a};
    mat c(3,3,fill::eye);
    std::cout <<"haha"<<endl;    
}
