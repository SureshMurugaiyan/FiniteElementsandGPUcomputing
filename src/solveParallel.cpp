#include <stdlib.h>
#include<stdio.h>
#include<math.h>
#include "mpi.h"
#include "input.h"
#include "inputParallel.h"
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Function Declaration
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void mpiCheck( int procRank, int totalTasks, int totProc);
void setcomm( int procRank, int* comm, int totProc,int nchunkx,int nchunky);
void iterateParallel(int processRank ,int* commMatrix);


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Main Function starts here                                                !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void solveParallel(int argc, char **argv){

int ntasks;
int rank;
double wtime;
int nfriend =4;  // for 2D we have 4 neighbors
int comm[nfriend];
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// MPI Initializations                                                      !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
wtime = MPI_Wtime();
mpiCheck( rank, ntasks, nproc);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Figure out neighbors of each chunk                                       !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
/* Figure out the neighbors of each block */
if(rank==0){printf(" Setting up the list of neighbors.\n");}
setcomm(rank,comm,nproc,nprocx,nprocy);

//if(rank==3){printf(" Neighbors of rank 0 are: \nTOP\t%i\tRIGHT\t%i\tBOTTOM\t%i\tLEFT\t%i\t\n",comm[0],comm[1],comm[2],comm[3]);}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Solving in Parallel Mode                                                 !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
if(rank==0)
{ printf(" Begin Solving in Parallel Mode.\n");}
//iterateParallel(rank,comm);

// Report Wall time
wtime=MPI_Wtime()-wtime;

printf("Task %i took %6.3f seconds\n",rank, wtime);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Finalize MPI                                                             !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
/* Terminate MPI */
MPI_Finalize();
if(rank==0){
  printf("\n");
  printf("Normal End of Parallel exectution.\n");
}

}

