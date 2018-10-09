// C++ Script to simulate 2D incompressible flow
#include <iostream>

void updateCn(double* Cn,double dt, int nxc,int nyc){
for(int i = 1; i<(nyc-1); ++i){
 for(int j =1; j<(nxc-1); ++j){
    Cn[i*nxc+j] = Cn[i*nxc+j]/dt;
  }
 }
}
