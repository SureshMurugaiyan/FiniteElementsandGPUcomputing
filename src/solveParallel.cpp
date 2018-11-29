#include <stdlib.h>
#include<stdio.h>
#include<math.h>
#include "mpi.h"

#include "input.h"
#include "inputParallel.h"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Function Declaration                                                     !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void InitializeField(double *phi, int row, int col);
void haloExchange(double* Phi, int* commMatrix, int procRank);
void updatePhysicalNorthBC(double* uxL,double* uyL,double* pL,int ProcRank);
void updatePhysicalSouthBC(double* uxL,double* uyL,double* pL,int ProcRank);
void  updatePhysicalEastBC(double* uxL,double* uyL,double* pL,int ProcRank);
void  updatePhysicalWestBC(double* uxL,double* uyL,double* pL,int ProcRank);
void  doNavierParallel(double* uxL,double* uyL,double* pL, int* itrNumber,
                       int* stopVal,FILE * pFileRank, int rank);
void printData(double* uxL, double *uyL, double *pL,
               int col, int row);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Main Function starts here                                                !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void solveParallel(int rank ,int* comm){

// Declaration of arrays for storing Primitive Variables
/**/
double uxL[ncGL]; // x component of velocity
double uyL[ncGL]; // y component of velocity
double pL[ncGL];  // Pressure

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
/*
double *uxL;
uxL = (double*) malloc(ncGL * sizeof(double));
double *uyL;
uyL = (double*) malloc(ncGL * sizeof(double));
double *pL;
pL = (double*) malloc(ncGL * sizeof(double));
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

// Initialization of arrays for storing Primitive Variables
InitializeField(uxL,nycGL,nxcGL);
InitializeField(uyL,nycGL,nxcGL);
InitializeField(pL,nycGL,nxcGL);

// Looping Variables
int itr = 0;
int stop =0;

FILE * FILE1;
if(rank==0){
FILE1 = fopen ("logParallel.txt" , "w");
}

// Time Marching Loop
while (stop ==0){
//while (itr<2001){
itr++;
//--------------- Update Ghost layers----------------------------/
haloExchange(uxL,comm,rank);
haloExchange(uyL,comm,rank);
haloExchange(pL,comm,rank);
//--------------- Update Physical boundary ----------------------/
updatePhysicalNorthBC(uxL,uyL,pL,rank);
updatePhysicalSouthBC(uxL,uyL,pL,rank);
updatePhysicalEastBC(uxL,uyL,pL,rank);
updatePhysicalWestBC(uxL,uyL,pL,rank);

//--------------- Solve Navier Stokes equation -------------------/
doNavierParallel(uxL,uyL,pL,&itr,&stop,FILE1,rank);

//--------------- Update Physical boundary ----------------------/
updatePhysicalNorthBC(uxL,uyL,pL,rank);
updatePhysicalSouthBC(uxL,uyL,pL,rank);
updatePhysicalEastBC(uxL,uyL,pL,rank);
updatePhysicalWestBC(uxL,uyL,pL,rank);

}
if(rank==0){
fprintf(FILE1,"Solution converged");
fclose(FILE1);
}

//printf("value of cell 1 is %6.3f for rank %d\n",uxL[19],rank);

}

