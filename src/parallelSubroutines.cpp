#include <iostream>
#include <mpi.h>
#include "input.h"
#include "inputParallel.h"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Function Declaration                                                     !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

void setOut(double* out,double* Phi ,int which);
void setIn(double* Phi,int which ,double* In);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Main Function starts here                                                !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void haloExchange(double* Phi, int* commMatrix, int procRank){

double outN[nxcL];
double outE[nycL];
double outS[nxcL];
double outW[nycL];

double inN[nxcL];
double inE[nycL];
double inS[nxcL];
double inW[nycL];

int partner;
int tag;
int nreq = 8;   // No of requests and No of status

MPI_Request requests[nreq];
MPI_Status status[nreq];

// Initialize requests
for( int i = 0; i<nreq; i++){
  requests[i] = MPI_REQUEST_NULL;
}
/* RECEIVING DATA TO GHOST CELLS*/
// Receive From North  Neighbor

if(commMatrix[0]==1){
  partner = procRank-nprocx;
  tag=0;
  MPI_Irecv(&inN,nxcL,MPI_DOUBLE,partner,tag,MPI_COMM_WORLD,&requests[0]);
}

// Receive From East  Neighbor

if(commMatrix[1]==1){
   partner = procRank+1;
   tag=1;
   MPI_Irecv(&inE,nycL,MPI_DOUBLE,partner,tag,MPI_COMM_WORLD,&requests[1]);
}

// Receive From South  Neighbor

if(commMatrix[2]==1){
   partner = procRank+nprocx;
   tag=2;
   MPI_Irecv(&inS,nxcL,MPI_DOUBLE,partner,tag,MPI_COMM_WORLD,&requests[2]);
}

// Receive From West  Neighbor

if(commMatrix[3]==1){
   partner = procRank-1;
   tag=3;
   MPI_Irecv(&inW,nycL,MPI_DOUBLE,partner,tag,MPI_COMM_WORLD,&requests[3]);
}

/* SENDING DATA TO GHOST CELLS */


// Send North Data to South boundary of neighbor

if(commMatrix[0]==1){
  partner = procRank-nprocx;
  tag=2;
  setOut(outN,Phi,0);
  MPI_Isend(&outN,nxcL,MPI_DOUBLE,partner,tag,MPI_COMM_WORLD,&requests[4]);
}

// Send East data to West boundary of neighbor

if(commMatrix[1]==1){
  partner = procRank+1;
  tag=3;
  setOut(outE,Phi,1);
  MPI_Isend(&outE,nycL,MPI_DOUBLE,partner,tag,MPI_COMM_WORLD,&requests[5]);
}

// Send South Data to North boundary of neighbor

if(commMatrix[2]==1){
  partner = procRank+nprocx;
  tag=0;
  setOut(outS,Phi,2);
  MPI_Isend(&outS,nxcL,MPI_DOUBLE,partner,tag,MPI_COMM_WORLD,&requests[6]);
}

// Send West Data to East boundary of neighbor

if(commMatrix[3]==1){
  partner = procRank-1;
  tag=1;
  setOut(outW,Phi,3);
  MPI_Isend(&outW,nycL,MPI_DOUBLE,partner,tag,MPI_COMM_WORLD,&requests[7]);
}
/* Wait for all communicationst to complete*/

MPI_Waitall(8,requests,status);

/* Copy Boundary values into Ghost cells, sent by neighbors into Matrix Phi */

if(commMatrix[0]==1){setIn(Phi,0,inN);}
if(commMatrix[1]==1){setIn(Phi,1,inE);}
if(commMatrix[2]==1){setIn(Phi,2,inS);}
if(commMatrix[3]==1){setIn(Phi,3,inW);}

return;
}

/* Pulls off the edge values of Phi to send to another task*/
void setOut(double* out,double* Phi ,int which){

  int i;
int j;
switch(which)
{
  case 0:{  // Pull of North boundary data
  i=1;
  for(j=1;j<((nxcL+2)-1);j++){out[j-1]=Phi[i*(nxcL+2)+j];}break;}

  case 1:{  // Pull of East boundary data
  j=(nxcL+2)-2;
  for(i=1;i<((nycL+2)-1);i++){out[i-1]=Phi[i*(nxcL+2)+j];}break;}

  case 2:{  // Pull of South boundary data
  i=(nycL+2)-2;
  for(j=1;j<((nxcL+2)-1);j++){out[j-1]=Phi[i*(nxcL+2)+j];}break;}

  case 3:{  // Pull of West boundary data
  j=1;
  for(i=1;i<((nycL+2)-1);i++){out[i-1]=Phi[i*(nxcL+2)+j];}break;}

}

return;
}
void setIn(double* Phi,int which ,double* In){

int i;
int j;
switch(which)
{
  case 0:{ // Copy data to North Ghost cell boundary
  i=0;
  for(j=1;j<((nxcL+2)-1);j++){Phi[i*(nxcL+2)+j]=In[j-1];}break;}

  case 1:{ // Copy data to East  Ghost cell boundary
  j=(nxcL+2)-1;
  for(i=1;i<((nycL+2)-1);i++){Phi[i*(nxcL+2)+j]=In[i-1];}break;}

  case 2:{ // Copy data to South Ghost cell boundary
  i=(nycL+2)-1;
  for(j=1;j<((nxcL+2)-1);j++){Phi[i*(nxcL+2)+j]=In[j-1];}break;}

  case 3:{ // Copy data to West  Ghost cell boundary
  j=0;
  for(i=1;i<((nycL+2)-1);i++){Phi[i*(nxcL+2)+j]=In[i-1];}break;}
}

return;
}

void sendLocal(double* sendPhi,double* PhiL,int col,int row){
for(int i=0;i<row;i++){
  for(int j=0;j<col;j++){
    sendPhi[i*col+j]=PhiL[(i+1)*col+(j+1)];
  }
}
}

