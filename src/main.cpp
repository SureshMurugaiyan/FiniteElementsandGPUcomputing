// C++ Script to simulate 2D incompressible flow
#include <iostream>
#include <math.h>
#include <fstream>
#include "functions.h"
#include "mesh.h"
#include "initialBoundaryConditions.h"
#include "updateBoundaryConditions.h"
#include "predictor.h"
#include "corrector.h"
#include "printData.h"
#include "updateCn.h"

using namespace std;

int main(int argc, char **argv)
{
//USER INPUT
// Opposite boundaries should have equal number of points
const int nxc = 32;        // number of cells in north and south boundary
const int nyc = 32;        // number of cells in east and west boundary
double Re = 100.0;          // User Input instead of material property

const int nx  = nxc+1;      // number of grid points in x-direction
const int ny  = nyc+1;      // number of grid points in y-direction
const int nc  = nxc*nyc;    // total number of cells
const int nv  = nx*ny;      // total number of vertices


const int nyHfc = ny;         // centers of Horizontal coordinate lines // row
const int nxHfc = nx-1;       // centers of Horizontal coordinate lines // column
const int nHfc  = nxHfc*nyHfc;// total number of horizontal face centers

const int nyVfc = ny-1;       // center of vertical coordinate lines  // row
const int nxVfc = nx;         // center of vertical coordinate lines  // column
const int nVfc  = nxVfc*nyVfc;// total number of vertical face centers

double X[nx*ny];
double Y[nx*ny];
double Z[nx*ny];

mesh2Dsquare(X,Y,Z,nxc,nyc,nx,ny,nc);
double dx = 1.0/nxc;
double dy = 1.0/nyc;
double dz = 1.0; // For 2D unit depth

double dt = 0.01;//Co*dx;
double dV = dx*dy*dz;

double ux[nc]; // x component of velocity
double uy[nc]; // y component of velocity
double uz[nc]; // z component of velocity
double p[nc];  // Pressure


double uxOld[nc]; // x component of velocity
double uyOld[nc]; // y component of velocity
double uzOld[nc]; // z component of velocity
double pOld[nc];  // Pressure

double Dnx[nc]; // allocating space for Diffusion term
double Dny[nc]; // allocating space for Diffusion term
double Cnx[nc]; // allocating space for Convection term
double Cny[nc]; // allocating space for Convection term

double DnxOld[nc]; // allocating space for Diffusion term
double DnyOld[nc]; // allocating space for Diffusion term
double CnxOld[nc]; // allocating space for Convection term
double CnyOld[nc]; // allocating space for Convection term

double Cn[nc]; //  allocating space for Convection source term in pressure poisson
double gradxP[nc]; //  allocating space for gradient
double gradyP[nc]; //  allocating space for gradient

double vertP[nv];    //  Interpolating pressure at vertices
double hfaceP[nHfc]; //  Interpolating Pressure at horizontal face center
double vfaceP[nVfc]; //  Iterpolating Pressure at vertical face center



//L2-norm- L2 norm at each time step
double L2ux = 1.0;
double L2uy = 1.0;
double L2uz = 1.0;


//L2-norm - Initial Correction
double L2oux = 1.0;
double L2ouy = 1.0;
double L2ouz = 1.0;


InitializeField(ux,nyc,nxc);
InitializeField(uy,nyc,nxc);
InitializeField(uz,nyc,nxc);
InitializeField(p,nyc,nxc);

boundaryCondition(ux,uy,p,nxc,nyc,nc);

InitializeArray(vertP,nv);
InitializeArray(hfaceP,nHfc);
InitializeArray(vfaceP,nVfc);


double normux = 1.0;
double normuy = 1.0;
double normuz = 1.0;


Laplacian(Dnx,ux,nyc,nxc,dx,dy);      //  Dnx Diffusion term
Laplacian(Dny,uy,nyc,nxc,dx,dy);      //  Dny Diffusion term
Div(Cnx,ux,ux,uy,nyc,nxc,dx,dy);      //  Cnx Convection term
Div(Cny,uy,ux,uy,nyc,nxc,dx,dy);      //  Cny Convection term

int k = 1;

while (normuy>(1e-6)){
//--------------Store old values-------------------//
storeOldValue(ux,uxOld,nc);
storeOldValue(uy,uyOld,nc);
storeOldValue(uz,uzOld,nc);
storeOldValue(p,pOld,nc);

storeOldValue(Dnx,DnxOld,nc);
storeOldValue(Dny,DnyOld,nc);
storeOldValue(Cnx,CnxOld,nc);
storeOldValue(Cny,CnyOld,nc);



//--------------- Predictor step--------------------------------//
 Laplacian(Dnx,ux,nyc,nxc,dx,dy);      //  Dnx Diffusion term
 Laplacian(Dny,uy,nyc,nxc,dx,dy);      //  Dny Diffusion term
 Div(Cnx,ux,ux,uy,nyc,nxc,dx,dy);      //  Cnx Convection term
 Div(Cny,uy,ux,uy,nyc,nxc,dx,dy);      //  Cny Convection term


//eulerPredictor(ux,uy,nxc,nyc,dt,dV,Cnx,Cny,Dnx,Dny,Re);
adamPredictor(ux,uy,nxc,nyc,dt,dV,Cnx,Cny,Dnx,Dny,CnxOld,CnyOld,DnxOld,DnyOld,Re);



//--------------Solve poisson eqn for Pressure ---------------//
Divergence(Cn,ux,uy,nyc,nxc,dx,dy); //  source term in poisson equation

updateCn(Cn,dt,nxc,nyc);

Poisson(p,nyc,nxc,dx,dy,Cn,nc);

//--------------- Corrector Step ----------------------

vertexInterpolate(vertP,p,ny,nx,nxc);
vfaceInterpolate(vfaceP,vertP,p,nyVfc,nxVfc,nx,nxc);
hfaceInterpolate(hfaceP,vertP,p,nyHfc,nxHfc,nx,nxc);
gradient(gradxP,gradyP,vfaceP,hfaceP,nyc,nxc,dx,dy,nxVfc,nxHfc);
//printmatrix(gradxP,nyc, nxc);

//gradient(gradxP,gradyP,p,nyc,nxc,dx,dy);

corrector(ux,uy,nxc,nyc,dt,gradxP,gradyP);

//--------------- Update Boundary condition ----------------------/

updateBoundaryCondition(ux,uy,p,nxc,nyc,nc);

//--------------- Update reference Pressure ----------------------/

refPressure(p,nyc, nxc);
if(k==1){
L2norm(ux,uxOld, &L2oux,nc);
L2norm(uy,uyOld, &L2ouy,nc);
L2norm(uz,uzOld, &L2ouz,nc);

L2norm(ux,uxOld, &L2ux,nc);
L2norm(uy,uyOld, &L2uy,nc);
L2norm(uz,uzOld, &L2uz,nc);

} else {
L2norm(ux,uxOld, &L2ux,nc);
L2norm(uy,uyOld, &L2uy,nc);
L2norm(uz,uzOld, &L2uz,nc);

}

normux = L2ux/L2oux;
normuy = L2uy/L2ouy;
normuz = L2uz/L2ouz;


cout<<endl<<"Itr:"<<k<<"\t"<<"uRes:"<<normux<<"\t"<<"vRes:"<<normuy<<"\t"<<endl;
k = k+1;
}
printData(ux,uy,p,nxc,nyc);

return 0;
//END OF PROGRAM

}

