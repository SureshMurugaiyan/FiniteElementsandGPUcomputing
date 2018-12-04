#include <iostream>
#include <math.h>
using namespace std;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Function Declarations                                                    !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void storeOldValue(double *phinew, double *phiOld,int totCell);
void L2norm(double *Phinew, double *Phiold,double *L2Phi,int  totCell);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Main Function-->Poisson Solver for Pressure Finite Volume Solver         !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void PoissonPressureCPU(double* Phi, int row, int col,
             double delX,double delY,double* source,
             int totCell){
double lam = 1;
int itr = 0;
int stop = 0;
while (stop==0){
itr++;
for(int i=1; i<(row-1); i++){
 for(int j=1; j<(col-1); j++){
   int k = i*col+j;
   double PhiP = Phi[k];
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];
   double AP   = (-2*delY/delX)-(2*delX/delY);
   double AS   = (delX/delY);
   double AW   = (delY/delX);
   double AE   = (delY/delX);
   double AN   = (delX/delY);

   double  R     = source[k]- AP*PhiP-AE*PhiE-AW*PhiW-AN*PhiN-AS*PhiS;
   double delPhi = R/AP;
   Phi[k] = Phi[k]+lam*delPhi;
  }
 }
//L2norm2( R,L2R,col,row);// L2 norm 
if(itr>1000){stop=1;}
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
__global__ void PoissonPressure(double* Phi, int row, int col,
             double delX,double delY,double* source,
             int totCell){
double lam = 1;
int itr = 0;
int stop = 0;
while (stop==0){
itr++;

   // Get global thread ID
  int k = blockIdx.x*blockDim.x+threadIdx.x;
  int n  = (row-1)*(col-1);
    // Do for only inner points
    if (k >0 && k<n) {

   //int k = i*col+j;
   double PhiP = Phi[k];
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];
   double AP   = (-2*delY/delX)-(2*delX/delY);
   double AS   = (delX/delY);
   double AW   = (delY/delX);
   double AE   = (delY/delX);
   double AN   = (delX/delY);

   double  R     = source[k]- AP*PhiP-AE*PhiE-AW*PhiW-AN*PhiN-AS*PhiS;
   double delPhi = R/AP;
   Phi[k] = Phi[k]+lam*delPhi;
  }
 
//L2norm2( R,L2R,col,row);// L2 norm 
if(itr>1000){stop=1;}
  }
}
