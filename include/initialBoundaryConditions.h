// C++ Script to simulate 2D incompressible flow
#include <iostream>

void boundaryCondition(double* ux, double *uy, double *p,
                  int nxc, int nyc,int nc){
// Boundary conditions
// North Boundary
for (int i = 1; i<nxc; i++)
{
    ux[i]= 1;
    uy[i]= 0;
    p[i]  = p[i+nxc];
}

// South Boundary
for (int i = nc-nxc; i<nc; i++)
{
    ux[i]= 0;
    uy[i]= 0;
    p[i]= p[i-nxc];
    
}

// West Boundary - Left end
for (int i = 0; i<nc; i=(i+nxc))
{
    ux[i]= 0;
    uy[i]= 0;
    p[i] = p[i+1];
}


// East Boundary - Right end
for (int i = nxc-1; i<nc; i=(i+nxc))
{
    ux[i]=0;
    uy[i]=0;
    p[i]  = p[i-1];
}

}
