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
const int nxc = 6;      // number of cells in north and south boundary
const int nyc = 6;      // number of cells in east and west boundary
double Re = 100.0;         // User Input instead of material property
double Co = 0.4; 	// Courant number

int noIteration=2;


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
double XC[(ny-1)*(nx-1)];
double YC[(ny-1)*(nx-1)];
double ZC[(ny-1)*(nx-1)];


double X_HFC[(ny)*(nx-1)];
double Y_HFC[(ny)*(nx-1)];
double Z_HFC[(ny)*(nx-1)];


double X_VFC[(ny-1)*(nx)];
double Y_VFC[(ny-1)*(nx)];
double Z_VFC[(ny-1)*(nx)];


mesh2Dsquare(X,Y,Z,XC,YC,ZC,X_HFC,Y_HFC,Z_HFC,X_VFC,Y_VFC,Z_VFC,nxc,nyc,nx,ny,nc);
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
//cout<<iteration<<endl;
iteration = iteration+1;

}

//printData(ux,uy,p,nxc,nyc);

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
//if(rank==0) printf("hello from bob\n");

// send a "message" from rank 0 to rank1

int ncL = (nc)/2;

if(rank==0){
  int messageN=ncL;
  int *messageOutXL = (int*)
    calloc(messageN, sizeof(int));
  for(int i=ncL+1; i<messageN;++i){
    XL[i] = X[i];
    YL[i] = Y[i];
    ZL[i] = Z[i];
    XCL[i] = XC[i];
    YCL[i] = YC[i];
    ZCL[i] = ZC[i];
    X_HFCL[i] = X_HFC[i];
    Y_HFCL[i] = Y_HFC[i];
    Z_HFCL[i] = Z_HFC[i];
    X_VFCL[i] = X_VFC[i];
    Y_VFCL[i] = Y_VFC[i];
    Z_VFCL[i] = Z_VFC[i];
  }
  int dest = 1;
  int tagX = 661;  // for vector X
  int tagY = 662;  // for vector Y
  int tagZ = 663;  // for vector Y
  int tagXC = 664; // for vector Y
  int tagYC = 665;  // for vector Y
  int tagZC = 666;  // for vector Y
  int tagX_HFC = 667; // for vector Y
  int tagY_HFC = 668;  // for vector Y
  int tagZ_HFC = 669;  // for vector Y
  int tagX_VFC = 670; // for vector Y
  int tagY_VFC = 671;  // for vector Y
  int tagZ_VFC = 672;  // for vector Y
  MPI_Request sendRequest;
  MPI_Isend(messageOutXL, messageN,MPI_INT,dest,tagX,MPI_COMM_WORLD,&sendRequest);
  MPI_Isend(messageOutYL, messageN,MPI_INT,dest,tagY,MPI_COMM_WORLD,&sendRequest);
  MPI_Isend(messageOutZL, messageN,MPI_INT,dest,tagY,MPI_COMM_WORLD,&sendRequest);
  MPI_Isend(messageOutXCL, messageN,MPI_INT,dest,tagXC,MPI_COMM_WORLD,&sendRequest);
  MPI_Isend(messageOutYCL, messageN,MPI_INT,dest,tagYC,MPI_COMM_WORLD,&sendRequest);
  MPI_Isend(messageOutZCL, messageN,MPI_INT,dest,tagZC,MPI_COMM_WORLD,&sendRequest);
  MPI_Isend(messageOutX_HFCL, messageN,MPI_INT,dest,tagX_HFC,MPI_COMM_WORLD,&sendRequest);
  MPI_Isend(messageOutY_HFCL, messageN,MPI_INT,dest,tagY_HFC,MPI_COMM_WORLD,&sendRequest);
  MPI_Isend(messageOutZ_HFCL, messageN,MPI_INT,dest,tagZ_HFC,MPI_COMM_WORLD,&sendRequest);
  MPI_Isend(messageOutX_VFCL, messageN,MPI_INT,dest,tagX_VFC,MPI_COMM_WORLD,&sendRequest);
  MPI_Isend(messageOutY_VFCL, messageN,MPI_INT,dest,tagY_VFC,MPI_COMM_WORLD,&sendRequest);
  MPI_Isend(messageOutZ_VFCL, messageN,MPI_INT,dest,tagZ_VFC,MPI_COMM_WORLD,&sendRequest);

// random task
int k = 87;
for(int j=0;j<10000;++j){
k = (k^j)*3;
}

MPI_Status status;
MPI_Wait(&sendRequest, &status);

  free(messageOutXL);  // blocking communication routine
}


// rank 1 receives a "message" from rank 0

if(rank==1){
  int tag =666;
  int messageN=ncL;
  int *messageIn = (int*)
    calloc(messageN, sizeof(int));
  int source=0;
  MPI_Request recvRequest;

  MPI_Irecv(messageIn,
            messageN,
            MPI_INT,
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
  for(int i=0; i<messageN;++i){
    printf("messageIn[%d]=%d\n",
           i,messageIn[i]);
  }
  free(messageIn);  // blocking communication routine
}



// tear down MPI( potentially block for all processes)
MPI_Finalize(); // Book keeping. All of the instances to wait
return 0;
//END OF PROGRAM

}

