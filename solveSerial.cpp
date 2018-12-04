/*---------------------*- C++ 2D Incompressible FLow -*-----------------------*
|  Solves the  2D incompressible Fluid Flow in 2D geometry                    |
|  User Input is located in include/input.h                                   |
|  Subroutine name: solveSerial-> Solves in Serial Mode                       |
*-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "inputSerial.h"
/*----------------------------------------------------------------------------*
|                    Function Declarations                                    |
*----------------------------------------------------------------------------*/
void InitializeField(double *phi, int row, int col);
void updateBoundaryCondition(double* ux, double *uy, double *p,int col,
                             int totCell);
void refValueUpdate(double* Phi, int row, int col, int refCell);
void doNavier(double* ux,double* uy,double* p,double* uxOld,double* uyOld,
                                        double* pOld, int itr);
void printData(double* ux, double *uy, double *p,int col, int row);
void normL2(double *Phi1new, double *Phi1old,
	    double *Phi2new, double *Phi2old,
            double *Phi3new, double *Phi3old,
            double *L2Phi,   int  totCell);
/*----------------------------------------------------------------------------*
|                      Main Subroutine                                        |
*----------------------------------------------------------------------------*/
void solveSerial(){

// Declaration of arrays for storing Primitive Variables
double *ux;ux = (double*) malloc(ncG * sizeof(double));
double *uy;uy = (double*) malloc(ncG * sizeof(double));
double *p;  p = (double*) malloc(ncG * sizeof(double));

// storing previous time step velocity,pressure
double *uxOld;uxOld = (double*) malloc(ncG * sizeof(double));
double *uyOld;uyOld = (double*) malloc(ncG * sizeof(double));
double *pOld; pOld  = (double*) malloc(ncG * sizeof(double));

//Storing L2 norm 
double *L2; L2  = (double*) malloc(nvar * sizeof(double));
double *L2o;L2o = (double*) malloc(nvar * sizeof(double));

// Declaration and Initialization of variables for L2 norm calculation
double normux=0.0;
double normuy=0.0;
double normp =0.0;

// Initialization of arrays for storing Primitive Variables
InitializeField(ux,nycG,nxcG);
InitializeField(uy,nycG,nxcG);
InitializeField(p ,nycG,nxcG);

InitializeField(uxOld,nycG,nxcG);
InitializeField(uyOld,nycG,nxcG);
InitializeField(pOld,nycG,nxcG);

InitializeField(L2,nvar,1);
InitializeField(L2o,nvar,1);

//------------------Update Boundary condition---------------------------------!
updateBoundaryCondition(ux,uy,p,nxcG,ncG);
//----------------------------------------------------------------------------!
FILE * FILE1;
FILE * FILE2;
FILE1 = fopen ("logSerial.txt" , "w");
FILE2 = fopen ("ResidualPlotting.txt" , "w");
// Write Data= logfile into file
fprintf(FILE1,"2D Navier Stokes Equation Using Finite Volume Method\n");
fprintf(FILE1,"GRID Size:\t %d \t X %d \n",nxc,nyc);
//fprintf(FILE1,"Time Step for the simulation: %f\n",dt);
fprintf(FILE1,"Reynolds number of the simulation:%f\n",Re);
//----------------------------------------------------------------------------!
// Looping Variables
int itr = 0;
int stop =0;
// Time Marching Loop
//while (itr<500){
while (stop==0){
itr++;
//--------------- Solve Navier Stokes equation -------------------/
doNavier(ux,uy,p,uxOld,uyOld,pOld,itr);
//--------------- Update reference Pressure ----------------------/
//refValueUpdate(p,nycG,nxcG,nxcG+2);
//--------------- Update Boundary condition ----------------------/
updateBoundaryCondition(ux,uy,p,nxcG,ncG);
//--------------------------------------------------------------------------!
// L2-Norm Calculation                                                      !
//--------------------------------------------------------------------------!
//L2norm1( ux,uxOld,uy,uyOld,p,pOld,L2,nxcG,nycG);// L2 norm 

normL2( ux,uxOld,uy,uyOld,p,pOld,L2,ncG);                        // L2 norm 
if(itr==1){L2o[0]=L2[0];L2o[1]=L2[1];L2o[2]=L2[2];}//initial correction

 normux=sqrt(L2[0]);
 normuy=sqrt(L2[1]);
 normp =sqrt(L2[2]);

printf("Iteration no:\t%d\t Ures: \t%.10e\t Vres: \t%.10e\t Pres: \t%.10e\t \n",itr,sqrt(L2[0]),sqrt(L2[1]),sqrt(L2[2]));
//--------------------------------------------------------------------------!
// Stopping Criteria                                                        !
//--------------------------------------------------------------------------!
if(MAXitr<=itr){
stop = 1;
}
if((normux<MAXnormux) &&(normuy<MAXnormuy)&&(normp<MAXnormp)){
stop = 1;
}
//--------------------------------------------------------------------------!
// Writing LogFile                                                          !
//--------------------------------------------------------------------------!
fprintf(FILE1,"Iteration no:\t%d\t Ures: \t%.10e\t Vres: \t%.10e\t Pres: \t%.10e\t \n",itr,sqrt(L2[0]/L2o[0]),sqrt(L2[1]/L2o[1]),sqrt(L2[2]/L2o[2]));
if(itr>0){fprintf(FILE2,"%d\t %6.10f\t  %6.10f\t  %6.10f\t \n",itr,sqrt(L2[0]/L2o[0]),sqrt(L2[1]/L2o[1]),sqrt(L2[2]/L2o[2]));}
}
//-----------END OF TIME STEP ^-----------------------------!
fprintf(FILE1,"Solution converged");
fclose(FILE1);fclose(FILE2);
//--------------------------------------------------------------------------!
// Printing cell center data
// Remove the ghost cell data 
double *uxInner;uxInner = (double*) malloc(nxc*nyc * sizeof(double));
double *uyInner;uyInner = (double*) malloc(nxc*nyc * sizeof(double));
double *pInner;pInner   = (double*) malloc(nxc*nyc * sizeof(double));
for(int i=0; i<nxc;i++){
  for(int j=0;j<nyc;j++){
    uxInner[i*nxc+j]= ux[(i+1)*nxcG+(j+1)];
    uyInner[i*nxc+j]= uy[(i+1)*nxcG+(j+1)];
    pInner[i*nxc+j] = p[(i+1)*nxcG+(j+1)];
  }
}

printData(uxInner,uyInner,pInner,nxc,nyc);

free(ux);
free(uy);
free(p);
free(uxOld);
free(uyOld);
free(pOld);
free(uxInner);
free(uyInner);
free(pInner);
free(L2o);
free(L2);
}
// * * * * * * * * * * END  OF SUB ROUTINE * * * * * * * * * * * * * * * * * //
