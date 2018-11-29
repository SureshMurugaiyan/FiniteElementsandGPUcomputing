#include <iostream>
#include <math.h>
#include "inputSerial.h"
#include <stdio.h>
using namespace std;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Function Declarations                                                    !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void InitializeField(double *phi, int row, int col);
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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Main Function                                                            !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

void doNavier(double* ux,double* uy,double* p, int* it, int* stop,FILE * FILE1){

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

//double vertP[nvG];    //  Interpolating pressure at vertices
//double hfaceP[nHfcG]; //  Interpolating Pressure at horizontal face center
//double vfaceP[nVfcG]; //  Iterpolating Pressure at vertical face center

double *vertP;
vertP = (double*) malloc(nvG * sizeof(double));
double *hfaceP;
hfaceP = (double*) malloc(nHfcG * sizeof(double));
double *vfaceP;
vfaceP = (double*) malloc(nVfcG * sizeof(double));

// Declaration of arrays for storing Derived Variables
//double Dnx[ncG]; // allocating space for Diffusion term
//double Dny[ncG]; // allocating space for Diffusion term
//double Cnx[ncG]; // allocating space for Convection term
//double Cny[ncG]; // allocating space for Convection term
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
double *Dnx;
Dnx = (double*) malloc(ncG * sizeof(double));
double *Dny;
Dny = (double*) malloc(ncG * sizeof(double));
double *Cnx;
Cnx = (double*) malloc(ncG * sizeof(double));
double *Cny;
Cny = (double*) malloc(ncG * sizeof(double));
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
//double Cn[ncG];     //  allocating space for Convection source term in pressure poisson
//double gradxP[ncG]; //  allocating space for gradient
//double gradyP[ncG]; //  allocating space for gradient
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
double *Cn;
Cn = (double*) malloc(ncG * sizeof(double));
double *gradxP;
gradxP = (double*) malloc(ncG * sizeof(double));
double *gradyP;
gradyP = (double*) malloc(ncG * sizeof(double));
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
//double uxOld[ncG]; // x component of velocity
//double uyOld[ncG]; // y component of velocity
//double pOld[ncG];  // Pressure
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
double *uxOld;
uxOld = (double*) malloc(ncG * sizeof(double));
double *uyOld;
uyOld = (double*) malloc(ncG * sizeof(double));
double *pOld;
pOld = (double*) malloc(ncG * sizeof(double));
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
//double DnxOld[ncG]; // allocating space for Diffusion term
//double DnyOld[ncG]; // allocating space for Diffusion term
//double CnxOld[ncG]; // allocating space for Convection term
//double CnyOld[ncG]; // allocating space for Convection term
double dt = 0.0001; // Initializing time step
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
double *DnxOld;
DnxOld = (double*) malloc(ncG * sizeof(double));
double *DnyOld;
DnyOld = (double*) malloc(ncG * sizeof(double));
double *CnxOld;
CnxOld = (double*) malloc(ncG * sizeof(double));
double *CnyOld;
CnyOld = (double*) malloc(ncG * sizeof(double));
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Initialize all the matrices                                              !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
InitializeField(uxOld,nycG,nxcG);
InitializeField(uyOld,nycG,nxcG);
InitializeField(pOld,nycG,nxcG);


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculate TimeStep at each iteration based on max velocity at each step  !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
timeStep(&dt,ux,uy); // Calculate the actual timeStep


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Things to be done for first iteration alone                              !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
if(*it==1){
// Calculate Laplacian and Divergence term from Initial conditions
Laplacian(Dnx,ux,nycG,nxcG,dx,dy);      //  Dnx Diffusion term
Laplacian(Dny,uy,nycG,nxcG,dx,dy);      //  Dny Diffusion term
Div(Cnx,ux,ux,uy,nycG,nxcG,dx,dy);      //  Cnx Convection term
Div(Cny,uy,ux,uy,nycG,nxcG,dx,dy);      //  Cny Convection term
// L2 norm for initial correction
L2norm(ux,uxOld, &L2oux,ncG);
L2norm(uy,uyOld, &L2ouy,ncG);
L2norm(p,pOld,   &L2op,ncG);

// Write Data= logfile into file
fprintf(FILE1,"2D Navier Stokes Equation Using Finite Volume Method\n");
fprintf(FILE1,"GRID Size:\t %d \t X %d \n",nxc,nyc);
fprintf(FILE1,"Time Step for the simulation: %f\n",dt);
fprintf(FILE1,"Reynolds number of the simulation:%f\n",Re);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Store old values                                                         !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
storeOldValue(ux,uxOld,ncG);
storeOldValue(uy,uyOld,ncG);
storeOldValue(p,pOld,ncG);

storeOldValue(Dnx,DnxOld,ncG);
storeOldValue(Dny,DnyOld,ncG);
storeOldValue(Cnx,CnxOld,ncG);
storeOldValue(Cny,CnyOld,ncG);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculation of Laplacian and Divergence for Predictor step               !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
 Laplacian(Dnx,ux,nycG,nxcG,dx,dy);      //  Dnx Diffusion term
 Laplacian(Dny,uy,nycG,nxcG,dx,dy);      //  Dny Diffusion term
 Div(Cnx,ux,ux,uy,nycG,nxcG,dx,dy);      //  Cnx Convection term
 Div(Cny,uy,ux,uy,nycG,nxcG,dx,dy);      //  Cny Convection term
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Predictor                                                                !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
//eulerPredictor(ux,uy,dt,Cnx,Cny,Dnx,Dny,nxcG,nycG,dV);
adamPredictor(ux,uy,dt,Cnx,Cny,Dnx,Dny,CnxOld,CnyOld,DnxOld,DnyOld,nxcG,nycG,dV);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculation of Source term in Poisson Equation                           !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
Divergence(Cn,ux,uy,nycG,nxcG,dx,dy); //  source term in poisson equation
updateCn(Cn,dt,nxcG,nycG);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Solver Poisson Equation For Pressure                                     !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
PoissonPressure(p,nycG,nxcG,dx,dy,Cn,ncG);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Calculation of pressure gradient                                         !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
//gradient(gradxP,gradyP,p,nycG,nxcG,dx,dy);  // Gradient computation using Finite difference-testing
vertexInterpolate(vertP,p,nyG,nxG,nxcG);
vfaceInterpolate(vfaceP,vertP,p,nyVfcG,nxVfcG,nxG,nxcG);
hfaceInterpolate(hfaceP,vertP,p,nyHfcG,nxHfcG,nxG,nxcG);
gradient(gradxP,gradyP,vfaceP,hfaceP,nycG,nxcG,dx,dy,nxVfcG,nxHfcG);   // Gradient computation using finite volume
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Corrector Step                                                           !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
corrector(ux,uy,gradxP,gradyP,dt,nxcG,nycG);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// L2-Norm Calculation                                                      !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
L2norm(ux,uxOld, &L2ux,ncG);
L2norm(uy,uyOld, &L2uy,ncG);
L2norm(p,pOld, &L2p,ncG);

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
// Stopping Criteria                                                        !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
if(500<*it){
*stop = 1;
}
if((normux<(1e-6)) &&(normuy<(1e-6))){
*stop = 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Writing LogFile                                                          !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
fprintf(FILE1,"Iteration no:\t%d\t Ures: \t%.10e\t Vres: \t%.10e\t \n",*it,normux,normuy);

}
