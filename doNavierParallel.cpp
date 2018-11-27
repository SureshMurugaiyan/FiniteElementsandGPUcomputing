#include <iostream>
#include <math.h>
#include "mpi.h"
#include "input.h"
#include "inputParallel.h"
#include <stdio.h>
using namespace std;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Function Declarations                                                    !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

void Div(double* Dn, double* Phi, double* U, double* V, int row, int col,double delX,double delY);
void Laplacian(double* Ln, double *Phi, int row, int col, double delX, double delY);
void timeStep(double* delt,double* ux,double* uy);
void storeOldValue(double *phinew, double *phiOld,int totCell);
void eulerPredictor(double* ux, double *uy,double delT,
                    double* Cnx, double *Cny,
                    double* Dnx,double* Dny,
                    int col, int row,double vol);
void adamPredictor(double* ux, double *uy,double delT,
                      double* Cnx, double *Cny,
                      double* Dnx,double* Dny,
                      double* CnxOld, double *CnyOld,
                      double* DnxOld,double* DnyOld,
                      int col, int row,double vol);
void Divergence(double* Dn, double* U, double* V,int row, int col, double delX, double delY);
void updateCn(double* Cn,double dt, int col,int row);
void PoissonPressure(double* Phi, int row, int col,
             double delX,double delY,double* source,
             int totCell);
void corrector(double* ux, double* uy,
               double* gradxP,double* gradyP,
               double dt, int col, int row);
void gradient(double* gradxPhi,double* gradyPhi,double* Phi,
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
void InitializeField(double *phi, int row, int col);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Main Function                                                            !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

void doNavierParallel(double* ux,double* uy,double* p, int* it, int* stop,FILE * FILE1,int rank){

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Variable Declarations                                                    !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Declaration and Initialization of variables for L2 norm calculation
double L2ux = 1.0;   //L2-norm- L2 norm at each time step
double L2uy = 1.0;
double L2p  = 1.0;

double L2oux = 1.0;  //L2-norm - Initial Correction
double L2ouy = 1.0;
double L2op  = 1.0;

double normux = L2ux/L2oux; //Normalized L2 norm
double normuy = L2uy/L2ouy;
double normp  = L2p/L2op;

double totNormUx = 1.0;
double totNormUy = 1.0;
double totNormP  = 1.0;


double vertP[nvGL];    //  Interpolating pressure at vertices
double hfaceP[nHfcGL]; //  Interpolating Pressure at horizontal face center
double vfaceP[nVfcGL]; //  Iterpolating Pressure at vertical face center

// Declaration of arrays for storing Derived Variables
double Dnx[ncGL]; // allocating space for Diffusion term
double Dny[ncGL]; // allocating space for Diffusion term
double Cnx[ncGL]; // allocating space for Convection term
double Cny[ncGL]; // allocating space for Convection term

double Cn[ncGL];     //  allocating space for Convection source term in pressure poisson
double gradxP[ncGL]; //  allocating space for gradient
double gradyP[ncGL]; //  allocating space for gradient

double uxOld[ncGL]; // x component of velocity
double uyOld[ncGL]; // y component of velocity
double pOld[ncGL];  // Pressure

double DnxOld[ncGL]; // allocating space for Diffusion term
double DnyOld[ncGL]; // allocating space for Diffusion term
double CnxOld[ncGL]; // allocating space for Convection term
double CnyOld[ncGL]; // allocating space for Convection term
double dt = 0.0001; // Initializing time step

double sendUx[nycL][nxcL]; // x component of velocity
double sendUy[nycL][nxcL]; // y component of velocity
double  sendP[nycL][nxcL];  // Pressure

double GlobalUx[nc]; // x component of velocity
double GlobalUy[nc]; // y component of velocity
double  GlobalP[nc];  // Pressure
double resultUx[nc]; // x component of velocity
double resultUy[nc]; // x component of velocity
double resultP[nc]; // x component of velocity
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Initialize all the matrices                                              !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
InitializeField(uxOld,nycGL,nxcGL);
InitializeField(uyOld,nycGL,nxcGL);
InitializeField(pOld,nycGL,nxcGL);


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculate TimeStep at each iteration based on max velocity at each step  !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
timeStep(&dt,ux,uy); // Calculate the actual timeStep


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Things to be done for first iteration alone                              !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
if(*it==1){
// Calculate Laplacian and Divergence term from Initial conditions
Laplacian(Dnx,ux,nycGL,nxcGL,dx,dy);      //  Dnx Diffusion term
Laplacian(Dny,uy,nycGL,nxcGL,dx,dy);      //  Dny Diffusion term
Div(Cnx,ux,ux,uy,nycGL,nxcGL,dx,dy);      //  Cnx Convection term
Div(Cny,uy,ux,uy,nycGL,nxcGL,dx,dy);      //  Cny Convection term
// L2 norm for initial correction
L2norm(ux,uxOld, &L2oux,ncGL);
L2norm(uy,uyOld, &L2ouy,ncGL);
L2norm(p,pOld,   &L2op,ncGL);

// Write Data= logfile into file
fprintf(FILE1,"2D Navier Stokes Equation Using Finite Volume Method\n");
fprintf(FILE1,"GRID Size:\t %d \t X %d \n",nxc,nyc);
fprintf(FILE1,"Time Step for the simulation: %f\n",dt);
fprintf(FILE1,"Reynolds number of the simulation:%f\n",Re);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Store old values                                                         !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
storeOldValue(ux,uxOld,ncGL);
storeOldValue(uy,uyOld,ncGL);
storeOldValue(p,pOld,ncGL);

storeOldValue(Dnx,DnxOld,ncGL);
storeOldValue(Dny,DnyOld,ncGL);
storeOldValue(Cnx,CnxOld,ncGL);
storeOldValue(Cny,CnyOld,ncGL);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculation of Laplacian and Divergence for Predictor step               !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
 Laplacian(Dnx,ux,nycGL,nxcGL,dx,dy);      //  Dnx Diffusion term
 Laplacian(Dny,uy,nycGL,nxcGL,dx,dy);      //  Dny Diffusion term
 Div(Cnx,ux,ux,uy,nycGL,nxcGL,dx,dy);      //  Cnx Convection term
 Div(Cny,uy,ux,uy,nycGL,nxcGL,dx,dy);      //  Cny Convection term
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Predictor                                                                !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
eulerPredictor(ux,uy,dt,Cnx,Cny,Dnx,Dny,nxcGL,nycGL,dV);
//adamPredictor(ux,uy,dt,Cnx,Cny,Dnx,Dny,CnxOld,CnyOld,DnxOld,DnyOld,nxcGL,nycGL,dV);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculation of Source term in Poisson Equation                           !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
Divergence(Cn,ux,uy,nycGL,nxcGL,dx,dy); //  source term in poisson equation
updateCn(Cn,dt,nxcGL,nycGL);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Solver Poisson Equation For Pressure                                     !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
PoissonPressure(p,nycGL,nxcGL,dx,dy,Cn,ncGL);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculation of pressure gradient                                         !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
gradient(gradxP,gradyP,p,nycGL,nxcGL,dx,dy);  // Gradient computation using Finite difference-testing
vertexInterpolate(vertP,p,nyGL,nxGL,nxcGL);
vfaceInterpolate(vfaceP,vertP,p,nyVfcL,nxVfcL,nxL,nxcL);
hfaceInterpolate(hfaceP,vertP,p,nyHfcL,nxHfcL,nxL,nxcL);
//gradient(gradxP,gradyP,vfaceP,hfaceP,nycGL,nxcGL,dx,dy,nxVfcGL,nxHfcGL);   // Gradient computation using finite volume

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Corrector Step                                                           !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
corrector(ux,uy,gradxP,gradyP,dt,nxcGL,nycGL);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// L2-Norm Calculation                                                      !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
L2norm(ux,uxOld, &L2ux,ncGL);
L2norm(uy,uyOld, &L2uy,ncGL);
L2norm(p,pOld, &L2p,ncGL);

int selectNorm = 1;    // choose 0 for normalized and 1 for direct norm

if(selectNorm ==0){
normux = L2ux/L2oux; // normalized norm wrt to initial correction
normuy = L2uy/L2ouy;
normp = L2p/L2op;
}else{
normux = L2ux;  // Actual norm without normalization wrt to initial correction
normuy = L2uy;
normp = L2p;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// L2-Norm Calculation-Sum of Norm over all processes                       !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Sum of Norm over all processes

MPI_Allreduce(&normux,&totNormUx,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(&normuy,&totNormUy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(&normp,&totNormP,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Stopping Criteria                                                        !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
if(500<*it){
*stop = 1;
}
if((totNormUx<(nproc*(1e-6))) &&(totNormUy<(nproc*(1e-6)))){
*stop = 1;
}

if(*stop==1){

// Sending Results to task 0
for(int i=0; i<nxcL;i++){
  for(int j=0;j<nycL;j++){
    sendUx[i][j]= ux[(i+1)*nxcGL+(j+1)];
    sendUy[i][j]= uy[(i+1)*nxcGL+(j+1)];
    sendP[i][j] =  p[(i+1)*nxcGL+(j+1)];
  }
}

 MPI_Gather(&sendUx,ncL,MPI_DOUBLE,&GlobalUx,ncL,MPI_DOUBLE,0,MPI_COMM_WORLD);
 MPI_Gather(&sendUy,ncL,MPI_DOUBLE,&GlobalUy,ncL,MPI_DOUBLE,0,MPI_COMM_WORLD);
 MPI_Gather(&sendP ,ncL,MPI_DOUBLE,&GlobalP ,ncL,MPI_DOUBLE,0,MPI_COMM_WORLD);

if(rank==0){
// we have gathered block by block into Global and is not in correct order
// Reordering of data before printing

int index =0;
for(int k=0;k<nprocx;k++){
  for( int l=0;l<nprocy;l++){
    for(int i =k*nxcL;i<(k+1)*nxcL;i++){
      for(int j=l*nycL;j<(l+1)*nycL;j++){
        resultUx[i*nxc+j]=GlobalUx[index];
        resultUy[i*nxc+j]=GlobalUy[index];
         resultP[i*nxc+j]=GlobalP[index];
        index++;
      }
    }
  }
}


}
if(rank==0){
printData(resultUx,resultUy,resultP,nxc,nyc);
//printData(ux,uy,p,nxcGL,nycGL);
}

}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Writing LogFile                                                          !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
if(rank==0){
fprintf(FILE1,"Iteration no:\t%d\t Ures: \t%.10e\t Vres: \t%.10e\t \n",*it,totNormUx,totNormUy);
}

}
