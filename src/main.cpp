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
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

using namespace std;

int main(int argc, char **argv)
{
//USER INPUT
// Opposite boundaries should have equal number of points
const int nxc = 64;      // number of cells in north and south boundary
const int nyc = 64;      // number of cells in east and west boundary
double Re = 100.0;         // User Input instead of material property
double Co = 0.4; 	// Courant number

int noIteration=13000;


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
double dz = dx;

double dt = 0.003; //Co*dx;
double dV = dx*dy;

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

double vertP[ny*nx];
double hfaceP[nHfc]; //  allocating space for gradient
double vfaceP[nVfc]; //  allocating space for gradient
boundaryCondition(ux,uy,p,nxc,nyc,nc);
//L2-norm
double L2ux = 1.0;
double L2uy = 1.0;
double L2uz = 1.0;
double L2p  = 1.0;

//L2-norm
double L2oux = 1.0;
double L2ouy = 1.0;
double L2ouz = 1.0;
double L2op  = 1.0;

//--------------- Predictor step--------------------------------//
 Laplacian(Dnx,ux,nyc,nxc,dx,dy);      //  Dnx Diffusion term
 Laplacian(Dny,uy,nyc,nxc,dx,dy);      //  Dny Diffusion term
 Div(Cnx,ux,ux,uy,nyc,nxc,dx,dy);      //  Cnx Convection term
 Div(Cny,uy,ux,uy,nyc,nxc,dx,dy);      //  Cny Convection term


// Parallelizing
// Initialize MPI: always do this "early"
MPI_Init(&argc, &argv);
int rank;
int size;

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
printf("This is process %d of %d total processes\n",rank,size);

// collective barrier
MPI_Barrier(MPI_COMM_WORLD);  // all processor should enter this barrier function // synchrony

// send a "message" from rank 0 to rank 1
  if(rank==0){
    int dimension = 3;      // size of message out
     double *dxdydzout = (double*) 
      calloc(dimension, sizeof(double));

    
      dxdydzout[0] = dx;
      dxdydzout[1] = dy;
      dxdydzout[2] = dz;

    int dest = 1;
    int tag = 666;
    
    MPI_Request sendRequest;
    
    MPI_Isend(dxdydzout, 
	      dimension,
	      MPI_DOUBLE,
	      dest,
	      tag,
	      MPI_COMM_WORLD,
	      &sendRequest);

    // random task
    int k = 87;
    for(int j=0;j<10000;++j){
      k = (k^j)*3;
    }
    
    MPI_Status status;
    MPI_Wait(&sendRequest, &status);

    free(dxdydzout);

  }

  // recv a "message" from rank 1
  if(rank==1){
    int dimension = 3; 
double *dxdydzin = (double*) 
      calloc(dimension, sizeof(double));


    int source = 0;
    int tag = 666;
    
    MPI_Request recvRequest;
    
    MPI_Irecv(dxdydzin, 
	      dimension,
	      MPI_DOUBLE,
	      source,
	      tag,
	      MPI_COMM_WORLD,
	      &recvRequest);
    
    // random task
    int k = 99;
    for(int j=0;j<1000;++j){
      k = (k^j)*3;
    }
    
    MPI_Status status;
    MPI_Wait(&recvRequest, &status);

double dx = dxdydzin[0];
double dy = dxdydzin[1];
double dz = dxdydzin[2];
cout<<dx<<endl;

    for(int i=0;i<dimension;++i){
      printf("dxdydzin[%d] = %4.8f\n",
	     i, dxdydzin[i]);
   }
    
   free(dxdydzin);

  }

  // tear down MPI (potentially block for all proceses)
  MPI_Finalize();




/*

int iteration = 1;
while (iteration < noIteration){
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


 eulerPredictor(ux,uy,nxc,nyc,dt,dV,Cnx,Cny,Dnx,Dny,Re);
//adamPredictor(ux,uy,nxc,nyc,dt,dV,Cnx,Cny,Dnx,Dny,CnxOld,CnyOld,DnxOld,DnyOld,Re);

//--------------Solve poisson eqn for Pressure ---------------//
Divergence(Cn,ux,uy,nyc,nxc,dx,dy); //  source term in poisson equation

updateCn(Cn,dt,nxc,nyc);

Poisson(p,nyc,nxc,dx,dy,Cn,nc);

//--------------- Corrector Step ----------------------

gradient(gradxP,gradyP,p,nyc,nxc,dx,dy);

corrector(ux,uy,nxc,nyc,dt,gradxP,gradyP);

//--------------- Update Boundary condition ----------------------/

updateBoundaryCondition(ux,uy,p,nxc,nyc,nc);

//--------------- Update reference Pressure ----------------------/

refPressure(p,nyc, nxc);
if(iteration==1){
L2norm(ux,uxOld, &L2oux,nc);
L2norm(uy,uyOld, &L2ouy,nc);
L2norm(uz,uzOld, &L2ouz,nc);
L2norm(p,pOld, &L2op,nc);
}

L2norm(ux,uxOld, &L2ux,nc);
L2norm(uy,uyOld, &L2uy,nc);
L2norm(uz,uzOld, &L2uz,nc);
L2norm(p,pOld, &L2p,nc);

cout<<endl<<"iteration no:"<<iteration<<"\t"<<"u residual:"<<L2ux/L2oux<<"\t"
     <<"v residual:"<<L2uy/L2ouy<<"\t"<<"P residual:"<<L2p/L2op<<"\t"<<endl;
iteration = iteration+1;
}

printData(ux,uy,p,nxc,nyc);
*/

return 0;
//END OF PROGRAM

}

