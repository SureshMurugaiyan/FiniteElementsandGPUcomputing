/*---------------------*- C++ 2D Incompressible FLow -*-----------------------*
|  Solves the  2D incompressible Fluid Flow in 2D geometry                    |
|  User Input is located in include/input.h                                   |
|  Subroutine name: solveSerial-> Solves in Serial Mode                       |
*-----------------------------------------------------------------------------*/
#include <stdio.h>
#include "inputSerial.h"
/*----------------------------------------------------------------------------*
|                    Function Declarations                                    |
*----------------------------------------------------------------------------*/
void InitializeField(double *phi, int row, int col);
void updateBoundaryCondition(double* ux, double *uy, double *p,int col,
                             int totCell);
void refValueUpdate(double* Phi, int row, int col, int refCell);
void doNavier(double* ux,double* uy,double* p, int* it, int* stop,
              FILE * pFile);
void printData(double* ux, double *uy, double *p,int col, int row);
/*----------------------------------------------------------------------------*
|                      Main Subroutine                                        |
*----------------------------------------------------------------------------*/
void solveSerial(){

// Declaration of arrays for storing Primitive Variables
double ux[ncG]; // x component of velocity
double uy[ncG]; // y component of velocity
double p[ncG];  // Pressure

// Initialization of arrays for storing Primitive Variables
InitializeField(ux,nycG,nxcG);
InitializeField(uy,nycG,nxcG);
InitializeField(p,nycG,nxcG);

// Looping Variables
int itr = 0;
int stop =0;
FILE * FILE1;
FILE1 = fopen ("logSerial.txt" , "w");
// Time Marching Loop
while (stop==0){
itr++;
//--------------- Update Boundary condition ----------------------/
updateBoundaryCondition(ux,uy,p,nxcG,ncG);
//--------------- Solve Navier Stokes equation -------------------/
doNavier(ux,uy,p,&itr,&stop,FILE1);
//--------------- Update reference Pressure ----------------------/
refValueUpdate(p,nycG,nxcG,pRefCell);
//--------------- Update Boundary condition ----------------------/
updateBoundaryCondition(ux,uy,p,nxcG,ncG);

}

fprintf(FILE1,"Solution converged");
fclose(FILE1);
printData(ux,uy,p,nxcG,nycG);

}
// * * * * * * * * * * END  OF SUB ROUTINE * * * * * * * * * * * * * * * * * //
