
#include <iostream>

void corrector(double* ux, double *uy,
               int nxc, int nyc,double dt,
               double* gradxP, double *gradyP){
    
for(int i = 1; i<(nyc-1); ++i){
    for(int j =1; j<(nxc-1); ++j){
        ux[i*nxc+j] = ux[i*nxc+j]-dt*gradxP[i*nxc+j];
        uy[i*nxc+j] = uy[i*nxc+j]-dt*gradyP[i*nxc+j];
    }
}
    
}
