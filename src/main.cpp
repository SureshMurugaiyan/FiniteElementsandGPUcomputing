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

int main()
{
//USER INPUT
// Opposite boundaries should have equal number of points
const int nxc = 32;      // number of cells in north and south boundary
const int nyc = 32;      // number of cells in east and west boundary
double Re = 100.0;         // User Input instead of material property
double Co = 0.4; 	// Courant number

int noIteration=20000;


const int nx  = nxc+1; // number of grid points in x-direction
const int ny  = nyc+1; // number of grid points in y-direction
const int nc  = nxc*nyc; // total number of cells

const int nxHfc = nx;  // column // center of vertical lines
const int nyHfc = ny-1; //row
const int nHfc  = nxHfc*nyHfc;
const int nxVfc = nx-1; // center of horizontal lines // column
const int nyVfc = ny; // row
const int nVfc  = nxVfc*nyVfc;

double X[nx*ny];
double Y[nx*ny];
double Z[nx*ny];

mesh2Dsquare(X,Y,Z,nxc,nyc,nx,ny,nc);
double dx = 1.0/nxc;
double dy = 1.0/nyc;

double dt = 0.00001; //Co*dx;
double dV = dx*dy;

double ux[nc]; // x component of velocity
double uy[nc]; // y component of velocity
double uz[nc]; // z component of velocity
double p[nc];  // Pressure



double Dnx[nc]; // allocating space for Diffusion term
double Dny[nc]; // allocating space for Diffusion term
double Cnx[nc]; // allocating space for Convection term
double Cny[nc]; // allocating space for Convection term

double Cn[nc]; //  allocating space for Convection source term in pressure poisson
double gradxP[nc]; //  allocating space for gradient
double gradyP[nc]; //  allocating space for gradient

double vertP[ny*nx];
double hfaceP[nHfc]; //  allocating space for gradient
double vfaceP[nVfc]; //  allocating space for gradient
boundaryCondition(ux,uy,p,nxc,nyc,nc);

int iteration = 1;
while (iteration < noIteration){

//--------------- Predictor step--------------------------------//
 Laplacian(Dnx,ux,nyc,nxc,dx,dy);      //  Dnx Diffusion term
 Laplacian(Dny,uy,nyc,nxc,dx,dy);      //  Dny Diffusion term
 Div(Cnx,ux,ux,uy,nyc,nxc,dx,dy);      //  Cnx Convection term
 Div(Cny,uy,ux,uy,nyc,nxc,dx,dy);      //  Cny Convection term


 eulerPredictor(ux,uy,nxc,nyc,dt,dV,Cnx,Cny,Dnx,Dny,Re);

//--------------Solve poisson eqn for Pressure ---------------//
Divergence(Cn,ux,uy,nyc,nxc,dx,dy); //  source term in poisson equation

updateCn(Cn,dt,nxc,nyc);

Poisson(p,nyc,nxc,dx,dy,Cn);

//--------------- Corrector Step ----------------------

gradient(gradxP,gradyP,p,nyc,nxc,dx,dy);

corrector(ux,uy,nxc,nyc,dt,gradxP,gradyP);

//--------------- Update Boundary condition ----------------------/

updateBoundaryCondition(ux,uy,p,nxc,nyc,nc);

//--------------- Update reference Pressure ----------------------/

refPressure(p,nyc, nxc);
cout<<iteration<<endl;
iteration = iteration+1;

}

printData(ux,uy,p,nxc,nyc);


return 0;
//END OF PROGRAM

}

