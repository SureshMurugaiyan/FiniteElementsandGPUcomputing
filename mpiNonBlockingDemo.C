#include <stdlib.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char **argv){
 // Initialize MPI: always do this early
 MPI_Init(&argc, &argv);

  int rank;
  int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("This is process %d of %d total processes\n",rank,size);

   MPI_Barrier(MPI_COMM_WORLD);
 if(rank==0) printf("hello from bob\n");

 // send a "message from rank 0 to rank 1
// for example data contained in a array

 if(rank ==0){
  int messageN = 10;
  int *messageOut = (int*)
   calloc(messageN, sizeof(int));
  for(int i = 0; i<messageN;++i){
    messageOut[i]=i;
  }
  int dest = 1;
  int tag = 666;
  MPI_Request sendRequest;

  MPI_Isend(messageOut,
           messageN,
           MPI_INT,
           dest,
           tag,
           MPI_COMM_WORLD,
           &sendRequest);
  //random task
   int k = 87;
   for(int j =0; j<10000;++j){
     k = (k^j)*3;
   }

   MPI_Status status;
   MPI_Wait(&sendRequest, &status);

  free(messageOut);  // blocking communicator- it copies data to some other buffer
 }

// You cannot send a large message to a receiver, he may not be able to receive. So have a receive function to receive large amount of data
// rank 1 receives a "message from rank 0


 if(rank ==1){
  int tag = 666;
  int messageN = 10;   // array size has to be same
  int *messageIn = (int*)
   calloc(messageN, sizeof(int));
                       // you dont have to initialize the array
                       // He just have to catch the message
  int source = 0;
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
    k=(k^j)*3;
  }

  MPI_Status status;
  MPI_Wait(&recvRequest, &status);
  for(int i = 0; i<messageN;++i){
    printf("messageIn[%d]=%d\n",
            i, messageIn[i]);
  }
  free(messageIn);  // blocking communicator- it copies data to some other buffer


 }
 // tear down MPI( potentially block for all processes)
 MPI_Finalize();

  return 0;
}
