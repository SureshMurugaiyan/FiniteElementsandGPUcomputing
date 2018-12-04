#include <iostream>
#include <math.h>
//#include "mpi.h"
#include "input.h"
#include "inputParallel.h"
#include <stdio.h>
#include <cuda.h>
using namespace std;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Function Declarations                                                    !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
__global__ void InitializeField(double *phi, int row, int col);
__global__ void Div(double* Dn, double* Phi, double* U, double* V, int row, int col,double delX,double delY);
__global__ void Laplacian(double* Ln, double *Phi, int row, int col, double delX, double delY);
__global__ void timeStep(double* delt,double* ux,double* uy);
__global__ void storeOldValue(double *phinew, double *phiOld,int totCell);
void eulerPredictor(double* ux, double *uy,double delT,
                    double* Cnx, double *Cny,
                    double* Dnx,double* Dny,
                    int col, int row,double vol);
__global__ void adamPredictor(double* ux, double *uy,double delT,
                      double* Cnx, double *Cny,
                      double* Dnx,double* Dny,
                      double* CnxOld, double *CnyOld,
                      double* DnxOld,double* DnyOld,
                      int col, int row,double vol);
__global__ void Divergence(double* Dn, double* U, double* V,int row, int col, double delX, double delY);
__global__ void updateCn(double* Cn,double dt, int col,int row);
__global__ void PoissonPressure(double* Phi, int row, int col,
             double delX,double delY,double* source,
             int totCell);
__global__ void corrector(double* ux, double* uy,
               double* gradxP,double* gradyP,
               double dt, int col, int row);
__global__ void gradient(double* gradxPhi,double* gradyPhi,double* Phi,
                        int row, int col, double delX, double delY);
void L2norm(double *Phinew, double *Phiold,double *L2Phi,int  totCell);
void vertexInterpolate(double* vertPhi, double* cellPhi, int vertRow, int vertCol, int cellCol);
void vfaceInterpolate(double* vfacePhi, double* vertPhi,
                      double* cellPhi, int vfaceRow, int vfaceCol,
                      int vertCol, int cellCol);
void hfaceInterpolate(double* hfacePhi, double* vertPhi,
                      double* cellPhi, int hfaceRow, int hfaceCol,
                      int vertCol, int cellCol);
void gradient(double* gradxPhi,double* gradyPhi,
              double* vfacePhi,double* hfacePhi,
			        int cellRow, int cellCol, double delX, double delY,
              int vfaceCol, int hfaceCol);
void sendLocal(double* sendPhi,double* PhiL,int col,int row);
void printData(double* ux, double *uy, double *p,
               int col, int row);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Main Function                                                            !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void doNavierParallel(double* h_ux,double* h_uy,double* h_p, 
                      double* h_uxOld,double* h_uyOld,double* h_pOld,
                      int itr,int rank){

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Variable Declarations                                                    !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Host vectors declaration and allocation is already done is previous func.

// Device vectors declaration and allocation
//  Interpolating pressure at vertices,horizontaland vertical face center 
// Device vectors declaration and allocation
// Size, in bytes, of each vector
//size_t bytes1 = nvGL*sizeof(double);
//double *vertP; cudaMalloc(&vertP, bytes1);  
//size_t bytes2 = nHfcGL*sizeof(double);
//double *hfaceP;cudaMalloc(&hfaceP, bytes2); 
//size_t bytes3 = nVfcGL*sizeof(double);
//double *vfaceP;cudaMalloc(&vfaceP, bytes3); 

// Declaration of arrays for storing Derived Variables
// allocating space for Diffusion term andConvection term
size_t bytes = ncGL*sizeof(double);

double *Dnx;cudaMalloc(&Dnx, bytes);
double *Dny;cudaMalloc(&Dny, bytes);
double *Cnx;cudaMalloc(&Cnx, bytes);
double *Cny;cudaMalloc(&Cny, bytes);
//allocating space for Convection term in poisson eqn & gradient
double *Cn;    cudaMalloc(&Cn, bytes);
double *gradxP;cudaMalloc(&gradxP, bytes);
double *gradyP;cudaMalloc(&gradyP, bytes);
// Storing previous timestep diffusion term and convection term
double *DnxOld;cudaMalloc(&DnxOld, bytes);
double *DnyOld;cudaMalloc(&DnyOld, bytes);
double *CnxOld;cudaMalloc(&CnxOld, bytes);
double *CnyOld;cudaMalloc(&CnyOld, bytes);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
double *ux;cudaMalloc(&ux, bytes);
double *uy;cudaMalloc(&uy, bytes);
double *p;cudaMalloc(&p, bytes);
double *uxOld;cudaMalloc(&uxOld, bytes);
double *uyOld;cudaMalloc(&uyOld, bytes);
double *pOld;cudaMalloc(&pOld, bytes);


// Copy host vectors to device
cudaMemcpy(ux,h_ux,bytes,cudaMemcpyHostToDevice);
cudaMemcpy(uy,h_uy,bytes,cudaMemcpyHostToDevice);
cudaMemcpy(p,h_p,bytes,cudaMemcpyHostToDevice);
cudaMemcpy(uxOld,h_uxOld,bytes,cudaMemcpyHostToDevice);
cudaMemcpy(uyOld,h_uyOld,bytes,cudaMemcpyHostToDevice);
cudaMemcpy(pOld, h_pOld,bytes,cudaMemcpyHostToDevice);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// GPU device
int blockSize, gridSize;
// Number of threads in each thread block
blockSize = 1024;
// Number of thread blocks in grid
gridSize = (int)ceil((float)ncGL/blockSize);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Initialize all the matrices                                              !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Execute the kernels
InitializeField<<<gridSize, blockSize>>>(uxOld,nycGL,nxcGL);
InitializeField<<<gridSize, blockSize>>>(uyOld,nycGL,nxcGL);
InitializeField<<<gridSize, blockSize>>>(pOld,nycGL,nxcGL);

InitializeField<<<gridSize, blockSize>>>(Dnx,nycGL,nxcGL);
InitializeField<<<gridSize, blockSize>>>(Dny,nycGL,nxcGL);
InitializeField<<<gridSize, blockSize>>>(Cnx,nycGL,nxcGL);
InitializeField<<<gridSize, blockSize>>>(Cny,nycGL,nxcGL);
InitializeField<<<gridSize, blockSize>>>(Cn,nycGL,nxcGL);
InitializeField<<<gridSize, blockSize>>>(gradxP,nycGL,nxcGL);
InitializeField<<<gridSize, blockSize>>>(gradyP,nycGL,nxcGL);
InitializeField<<<gridSize, blockSize>>>(DnxOld,nycGL,nxcGL);
InitializeField<<<gridSize, blockSize>>>(DnyOld,nycGL,nxcGL);
InitializeField<<<gridSize, blockSize>>>(CnxOld,nycGL,nxcGL);
InitializeField<<<gridSize, blockSize>>>(CnyOld,nycGL,nxcGL);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculate TimeStep at each iteration based on max velocity at each step  !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
double dt = 0.000000001; // Initializing time step
//timeStep<<<gridSize, blockSize>>>(&dt,ux,uy); // Calculate the actual timeStep    //check

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Things to be done for first iteration alone                              !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
if(itr==1){
// Calculate Laplacian and Divergence term from Initial conditions
Laplacian<<<gridSize, blockSize>>>(Dnx,ux,nycGL,nxcGL,dx,dy);      //  Dnx Diffusion term
Laplacian<<<gridSize, blockSize>>>(Dny,uy,nycGL,nxcGL,dx,dy);      //  Dny Diffusion term
Div<<<gridSize, blockSize>>>(Cnx,ux,ux,uy,nycGL,nxcGL,dx,dy);      //  Cnx Convection term
Div<<<gridSize, blockSize>>>(Cny,uy,ux,uy,nycGL,nxcGL,dx,dy);      //  Cny Convection term 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Store old values                                                         !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
storeOldValue<<<gridSize, blockSize>>>(ux,uxOld,ncGL);
storeOldValue<<<gridSize, blockSize>>>(uy,uyOld,ncGL);
storeOldValue<<<gridSize, blockSize>>>(p,pOld,ncGL);

storeOldValue<<<gridSize, blockSize>>>(Dnx,DnxOld,ncGL);
storeOldValue<<<gridSize, blockSize>>>(Dny,DnyOld,ncGL);
storeOldValue<<<gridSize, blockSize>>>(Cnx,CnxOld,ncGL);
storeOldValue<<<gridSize, blockSize>>>(Cny,CnyOld,ncGL);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculation of Laplacian and Divergence for Predictor step               !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
 Laplacian<<<gridSize, blockSize>>>(Dnx,ux,nycGL,nxcGL,dx,dy);      //  Dnx Diffusion term
 Laplacian<<<gridSize, blockSize>>>(Dny,uy,nycGL,nxcGL,dx,dy);      //  Dny Diffusion term
 Div<<<gridSize, blockSize>>>(Cnx,ux,ux,uy,nycGL,nxcGL,dx,dy);      //  Cnx Convection term
 Div<<<gridSize, blockSize>>>(Cny,uy,ux,uy,nycGL,nxcGL,dx,dy);      //  Cny Convection term
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Predictor                                                                !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
 //eulerPredictor(ux,uy,dt,Cnx,Cny,Dnx,Dny,nxcGL,nycGL,dV);
adamPredictor<<<gridSize, blockSize>>>(ux,uy,dt,Cnx,Cny,Dnx,Dny,CnxOld,CnyOld,DnxOld,DnyOld,nxcGL,nycGL,dV);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculation of Source term in Poisson Equation                           !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
Divergence<<<gridSize, blockSize>>>(Cn,ux,uy,nycGL,nxcGL,dx,dy); //  source term in poisson equation
updateCn<<<gridSize, blockSize>>>(Cn,dt,nxcGL,nycGL);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Solver Poisson Equation For Pressure                                     !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
PoissonPressure<<<gridSize, blockSize>>>(p,nycGL,nxcGL,dx,dy,Cn,ncGL);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculation of pressure gradient                                         !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
gradient<<<gridSize, blockSize>>>(gradxP,gradyP,p,nycGL,nxcGL,dx,dy);  // Gradient computation using Finite difference-testing
//vertexInterpolate<<<gridSize, blockSize>>>(vertP,p,nyGL,nxGL,nxcGL);
//vfaceInterpolate<<<gridSize, blockSize>>>(vfaceP,vertP,p,nyVfcGL,nxVfcGL,nxGL,nxcGL);
//hfaceInterpolate<<<gridSize, blockSize>>>(hfaceP,vertP,p,nyHfcGL,nxHfcGL,nxGL,nxcGL);
//gradient<<<gridSize, blockSize>>>(gradxP,gradyP,vfaceP,hfaceP,nycGL,nxcGL,dx,dy,nxVfcGL,nxHfcGL);   // Gradient computation using finite volume

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Corrector Step                                                           !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
corrector<<<gridSize, blockSize>>>(ux,uy,gradxP,gradyP,dt,nxcGL,nycGL);

// Copy host vectors to device
cudaMemcpy(h_ux,ux,bytes,cudaMemcpyDeviceToHost);
cudaMemcpy(h_uy,uy,bytes,cudaMemcpyDeviceToHost);
cudaMemcpy(h_p,p,bytes,cudaMemcpyDeviceToHost);
cudaMemcpy(h_uxOld,uxOld,bytes,cudaMemcpyDeviceToHost);
cudaMemcpy(h_uyOld,uyOld,bytes,cudaMemcpyDeviceToHost);
cudaMemcpy(h_pOld,pOld,bytes,cudaMemcpyDeviceToHost);

// Release device memory
//free(vertP);
//free(hfaceP);
//free(vfaceP);
cudaFree(Dnx);
cudaFree(Dny);
cudaFree(Cnx);
cudaFree(Cny);
cudaFree(Cn);
cudaFree(gradxP);
cudaFree(gradyP);
cudaFree(DnxOld);
cudaFree(DnyOld);
cudaFree(CnxOld);
cudaFree(CnyOld);

cudaFree(ux);
cudaFree(uy);
cudaFree(p);
cudaFree(uxOld);
cudaFree(uyOld);
cudaFree(pOld);

// Release host memory is done in mother program of this routine

}
