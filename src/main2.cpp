// C++ Script to simulate 2D incompressible flow in parallel Mode with OpenMPI

// HeaderFiles for C++ 
#include <iostream>
#include <math.h>
#include <fstream>
//#include "functions.h"
#include "mesh.h"
#include "initialBoundaryConditions.h"
#include "updateBoundaryConditions.h"
#include "predictor.h"
#include "corrector.h"
#include "printData.h"
#include "updateCn.h"
#include "iterate2DNavier.h"
using namespace std;

// HeaderFiles for OpenMPI and C
# include <stdlib.h>
# include <stdio.h>
# include "mpi.h"

# define n 32            /* Total number of cells in each direction X and Y : Matrix= 32 X 32 */
# define nodeedge 16     /* Total number of cells in each direction X and Y for each chunks : Matrix 16 X 16 */
# define nblock n/nodeedge   /* number of tasks per row of matrix : Matrix 2 X 2 = 4 blocks */
# define nproc nblock*nblock /* total number of tasks (processors) = 4 procs/tasks     */

int main ( int argc, char **argv );
void doblack ( double w, double M[][nodeedge+2] );
void dored ( double w, double M[][nodeedge+2] );
void exchange ( double M[][nodeedge+2], int comm[], int rank );
void iterate ( double w, double M[][nodeedge+2], double result[][n], int rank, int comm[] );
void setcomm ( int rank, int comm[] );
void setex ( double ex[], double M[][nodeedge+2], int which );
void initialize_matrix ( double M[][nodeedge+2] );
void unpack ( double M[][nodeedge+2], int where, double in[] );

/******************************************************************************/

int main ( int argc, char **argv )

/******************************************************************************/
/*
  Purpose:

    SOLVER FOR 2D INCOMPRESSIBLE MHD EQUATIONS USING MPI.


  Reference: To do Reference


  Local parameters:

    Local, int COMM[4], contains a 0 (no) or 1 (yes) if
    communication is needed for the UP(0), RIGHT(1), DOWN(2)
    and LEFT(3) neighbors.

    Local, FILE *fp, a pointer to the output file.

    Local, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process.

    Local, double RESULT[n][n], the results for the complete problem,
    kept by process 0.

    Local, double W, the SOR factor, which must be strictly between 0 and 2.
*/ 
{
  int comm[4];
  FILE *fp;
  int i;
  int j;
  double M[nodeedge+2][nodeedge+2];
  int ntasks;
  int rank;
  double result[n][n];
  double w;
  double wtime;

const int nxc = 32;          // number of cells in north and south boundary
const int nyc = 32;          // number of cells in east and west boundary
const int nx  = nxc+1;       // number of grid points in x-direction
const int ny  = nyc+1;       // number of grid points in y-direction
const int nc  = nxc*nyc;     // total number of cells
const int nxcL = nodeedge+2;        // number of cells in north and south boundary Local
const int nycL = nodeedge+2;        // number of cells in east and west boundary Local
const int ncL  = (nodeedge+2)*(nodeedge+2);    // total number of cells
const int nv  = nx*ny;      // total number of vertices



double X[nx*ny];
double Y[nx*ny];
double Z[nx*ny];

double ux[nc]; // x component of velocity
double uy[nc]; // y component of velocity
double uz[nc]; // z component of velocity
double p[nc];  // Pressure

double uxL[ncL]; // x component of velocity
double uyL[ncL]; // y component of velocity
double uzL[ncL]; // z component of velocity
double pL[ncL];  // Pressure


double uxOld[nc]; // x component of velocity
double uyOld[nc]; // y component of velocity
double uzOld[nc]; // z component of velocity
double pOld[nc];  // Pressure





  MPI_Init ( &argc, &argv );      //done

  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );  //done

  MPI_Comm_size ( MPI_COMM_WORLD, &ntasks );  //done 

  wtime = MPI_Wtime ( );   //done

  if ( rank == 0 )     //done 
  {
    printf ( "\n" );
    printf ( "SOLVER FOR 2D INCOMPRESSIBLE MHD EQUATIONS:\n" );
    printf ( "  C/MPI version\n" );
    printf ( "  SOLVER USING MPI.\n" );
  }

  if ( ntasks != nproc )    //done 
  {
    if ( rank == 0 ) 
    {
      printf ( "\n" );
      printf ( "Fatal error!\n" );
      printf ( "  MP_PROCS should be set to %i!\n", nproc );
    }
    MPI_Finalize ( );
    exit ( 1 );
  }

  if ( rank == 0 )    //done 
  {
    printf ( "\n" );
    printf ( "  MPI has been set up.\n" );
  }
/*  Initialize the matrix M. */
  if ( rank == 0 )     //done 
  {
    printf ( "  Initialize the matrix M.\n" );
  }

  mesh2Dsquare(X,Y,Z,nxc,nyc,nx,ny,nc);
  double dx = 1.0/nxc;
  double dy = 1.0/nyc;
  double dz = 1.0; // For 2D unit depth
  double dt = 0.01;//Co*dx;
  double dV = dx*dy*dz;
  
  InitializeField(ux,nyc,nxc);
  InitializeField(uy,nyc,nxc);
  InitializeField(uz,nyc,nxc);
  InitializeField(p,nyc,nxc);

   boundaryCondition(ux,uy,p,nxc,nyc,nc);


  initialize_matrix ( M ); //test

/*  Figure out who I communicate with.*/
  if ( rank == 0 ) 
  {
    printf ( "  Set the list of neighbors.\n" );  //done 
  }
  setcomm ( rank, comm ); //done 
/* Update M, using SOR value W, until convergence.*/
  if ( rank == 0 ) 
  {
    printf ( "  Begin the iteration.\n" );
  }



  w = 1.2;  //test
  iterate ( w, M, result, rank, comm ); //test
  
  iterate2DNavier(uxL,uyL,pL,ux,uy,p,ncL,nxcL,nycL,dx,dy,dz,dV,dt,rank,comm);



/* Report timing */ 
  wtime = MPI_Wtime ( ) - wtime;

  printf ( "  Task %i took %6.3f seconds\n", rank, wtime );

/*Write the solution to a file.*/
  if ( rank == 0 )
  {
    fp = fopen ( "laplace_solution.txt", "w" );

    for ( i = 0; i < n; i++ ) 
    {
      for ( j = 0; j < n; j++ )
      {
        fprintf ( fp, "%f \n", result[i][j] );
      }  
    }
    fclose ( fp );
    printf ( "  Solution written to \"laplace_solution.txt\".\n" );
  }
/*
  Terminate MPI.
*/
  MPI_Finalize ( );
/*
  Terminate.
*/
  if ( rank == 0 )
  {
    printf ( "\n" );
    printf ( "LAPLACE_MPI:\n" );
    printf ( "  Normal end of execution.\n" );
  }
  return 0;
}
/******************************************************************************/

void doblack ( double w, double M[][nodeedge+2] )

/******************************************************************************/
/*
  Purpose:

    DOBLACK iterates on the upper right and lower left corners of my matrix.

  Parameters:

    Input, double W, the SOR factor, which must be strictly between 0 and 2.

    Input/output, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process.
*/
{
  int i;
  int j;
/*
  Upper right corner.
*/
  for ( i = 1; i <= nodeedge / 2; i++ )
  {
    for ( j = nodeedge / 2 + 1; j <= nodeedge; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
/*
  Lower left corner.
*/
  for ( i = nodeedge / 2 + 1; i <= nodeedge; i++ )
  {
    for ( j = 1; j <= nodeedge / 2; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
  return;
}
/******************************************************************************/

void dored ( double w, double M[][nodeedge+2] )

/******************************************************************************/   
/*
  Purpose:

    DORED iterates on the upper left and lower right corners of my matrix.

  Parameters:

    Input, double W, the SOR factor, which must be strictly between 0 and 2.

    Input/output, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process.
*/  
{
  int i;
  int j;
/*
  Upper left corner.
*/
  for ( i = 1; i <= nodeedge / 2; i++ )
  {
    for ( j = 1; j <= nodeedge / 2; j++ ) 
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
/*
  Lower right corner.
*/
  for ( i = nodeedge / 2 + 1; i <= nodeedge; i++ )
  {
    for ( j = nodeedge / 2 + 1; j <= nodeedge; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
  return;
}
/******************************************************************************/

void exchange ( double M[][nodeedge+2], int comm[], int rank )

/******************************************************************************/
/*
  Purpose:

   EXCHANGE trades edge values with up to four neighbors.

  Discussion:

    Up to 4 MPI sends are carried out, and up to 4 MPI receives.


  Parameters:

    Input/output, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process.

    Input, int COMM[4], contains a 0 (no) or 1 (yes) if
    communication is needed for the UP(0), RIGHT(1), DOWN(2)
    and LEFT(3) neighbors.

    Input, int RANK, the rank of this process.
*/
{
  double ex0[nodeedge];
  double ex1[nodeedge];
  double ex2[nodeedge];
  double ex3[nodeedge];
  int i;
  double in0[nodeedge];
  double in1[nodeedge];
  double in2[nodeedge];
  double in3[nodeedge];
  int partner;
  MPI_Request requests[8];
  MPI_Status status[8];
  int tag;
/* 
  Initialize requests.
*/
  for ( i = 0; i < 8; i++ ) 
  {
    requests[i] = MPI_REQUEST_NULL; 
  }
/* 
  Receive from UP neighbor (0).
*/
  if ( comm[0] == 1 )
  {
    partner = rank - nblock;
    tag = 0;
    MPI_Irecv ( &in0, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD, 
      &requests[0] );
  }
/*
  Receive from RIGHT neighbor (1).
*/
  if ( comm[1] == 1 )
  {
    partner = rank + 1;
    tag = 1;
    MPI_Irecv ( &in1, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[1] );
  }
/*
  Receive from DOWN neighbor (2).
*/
  if ( comm[2] == 1 )
  {
    partner = rank + nblock;
    tag = 2;
    MPI_Irecv ( &in2, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[2] );
  }
/*
  Receive from LEFT neighbor (3).
*/
  if ( comm[3] == 1 )
  {
    partner = rank - 1;
    tag = 3;
    MPI_Irecv ( &in3, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[3] );
  }
/*
  Send up from DOWN (2) neighbor.
*/
  if ( comm[0] == 1 )
  {
    partner = rank - nblock;
    tag = 2;
    setex ( ex0, M, 0 );
    MPI_Isend ( &ex0, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[4] );
  }
/*
  Send right form LEFT (3) neighbor.
*/
  if (comm[1] == 1 )
  {
    partner = rank + 1;
    tag = 3;
    setex ( ex1, M, 1 );
    MPI_Isend ( &ex1, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[5] );
  }
/*
  Send down from UP (0) neighbor.
*/
  if ( comm[2] == 1 )
  {
    partner = rank + nblock;
    tag = 0;
    setex ( ex2, M, 2 );
    MPI_Isend ( &ex2, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[6] );
  }
/*
  Send left from RIGHT (1) neighbor.
*/
  if ( comm[3] == 1 )
  {
    partner = rank - 1;
    tag = 1;
    setex ( ex3, M, 3 );
    MPI_Isend ( &ex3, nodeedge, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[7] );
  }
/* 
  Wait for all communication to complete.
*/ 
  MPI_Waitall ( 8, requests, status );
/*
  Copy boundary values, sent by neighbors, into M.
*/
  if ( comm[0] == 1 ) 
  {
    unpack ( M, 0, in0 );
  }
  if ( comm[1] == 1 ) 
  {
    unpack ( M, 1, in1 );
  }
  if ( comm[2] == 1 ) 
  {
    unpack ( M, 2, in2 );
  }
  if ( comm[3] == 1 ) 
  {
    unpack ( M, 3, in3 );
  }

  return;
}
/******************************************************************************/

void initialize_matrix ( double M[][nodeedge+2] )

/******************************************************************************/
/*
  Purpose:

    INITIALIZE_MATRIX initializes the partial results array M.

  Parameters:

    Output, double M[nodeedge+2][nodeedge+2], the initialized partial 
    results array.
*/
{
  double avg;
  double bv[4];
  int i;
  int j;

  bv[0] = 100.0;
  bv[1] = 0.0;
  bv[2] = 0.0;
  bv[3] = 0.0;
/* 
  Put the boundary values into M.
*/ 
  for ( i = 1; i <= nodeedge; i++ )
  { 
    M[0][i] =          bv[0];
    M[i][nodeedge+1] = bv[1];
    M[nodeedge+1][i] = bv[2];
    M[i][0] =          bv[3];
  }
/* 
  Set all interior values to be the average of the boundary values.
*/ 
  avg = ( bv[0] + bv[1] + bv[2] + bv[3] ) / 4.0;

  for ( i = 1; i <= nodeedge; i++ )
  {
    for ( j = 1; j <= nodeedge; j++ )
    {
      M[i][j] = avg;
    }
  }

  return;
}
/******************************************************************************/

void iterate ( double w, double M[][nodeedge+2], double result[][n], int rank, 
  int comm[] )

/******************************************************************************/
/*
  Purpose:

    ITERATE controls the iteration, including convergence checking.

  Parameters:

    Input, double W, the SOR factor, which must be strictly between 0 and 2.

    Input/output, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process.

    Output, double RESULT[n][n], the results for the complete problem,
    kept by process 0.

    Input, int RANK, the rank of the process.

    Input, int COMM[4], contains a 0 (no) or 1 (yes) if
    communication is needed for the UP(0), RIGHT(1), DOWN(2)
    and LEFT(3) neighbors.

  Local parameters:

    Local, int COUNT, the length, in elements, of messages.

    Local, double DIFF, the average absolute difference in elements
    of M since the last comparison.

    Local, int IT, the iteration counter.

    Local, double MM[n*n], a vector, used to gather the data from
    all processes.  This data is then rearranged into a 2D array.
*/
{
  int count;
  double diff;
  int done;
  double ediff;
  int i;
  double in;
  int index;
  int it;
  int j;
  int k;
  int l;
  double MM[n*n];
  double mold[nodeedge+2][nodeedge+2];
  double send[nodeedge][nodeedge];

  it = 0;
  done = 0;
  for ( i = 1; i <= nodeedge; i++ )
  {
    for ( j = 1; j <= nodeedge; j++ )
    {
      mold[i][j] = M[i][j];
    }
  }

  while ( done == 0 )
  {
    it++;
/*
  Exchange values with neighbors, update red squares, exchange values
  with neighbors, update black squares.
*/
    exchange ( M, comm, rank );
    dored ( w, M );
    exchange ( M, comm, rank );
    doblack ( w, M );
/*
  Check for convergence every 20 iterations.
  Find the average absolute change in elements of M.
  Maximum iterations is 5000.
*/
    if ( 5000 < it )
    {
      done = 1;
    }

    if ( ( it % 20 == 0.0 ) && ( done != 1 ) )
    { 
      diff = 0.0;
      for ( i = 1; i <= nodeedge; i++ )
      {
        for ( j = 1; j <= nodeedge; j++ )
        {
          ediff = M[i][j] - mold[i][j];
          if ( ediff < 0.0 ) 
          {
            ediff = - ediff;
          }
          diff = diff + ediff;
          mold[i][j] = M[i][j];
        }
      }
      diff = diff / ( ( double ) ( nodeedge * nodeedge ) );
/*
  IN = sum of DIFF over all processes.
*/
      MPI_Allreduce ( &diff, &in, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

      if ( in < ( double ) nproc * 0.0001 ) 
      {
        done = 1;
      }
    }

  }
/* 
  Send results to task 0.
*/ 
  for ( i = 0; i < nodeedge; i++ )
  {
    for ( j = 0; j < nodeedge; j++ )
    {
      send[i][j] = M[i+1][j+1];
    }
  }

  count = nodeedge * nodeedge;

  MPI_Gather ( &send, count, MPI_DOUBLE, &MM, count, MPI_DOUBLE, 0, 
    MPI_COMM_WORLD );

  printf ( "  ITERATE gathered updated results to process 0.\n" );
/* 
  Storage on task 0 has to be consistent with a NBLOCK x NBLOCK decomposition.
*/
  if ( rank == 0 ) 
  {
    printf ( "\n did %i iterations\n", it );

    index = 0;
    for ( k = 0; k < nblock; k++ )
    {
      for ( l = 0; l < nblock; l++ )
      {
        for ( i = k * nodeedge; i < ( k + 1 ) * nodeedge; i++ )
        {
          for ( j = l * nodeedge; j < ( l + 1 ) * nodeedge; j++ )
          {
            result[i][j] = MM[index];
            index++;
          }
        }
      }
    }
  }
  return;
}
/******************************************************************************/

void setcomm ( int rank, int comm[] )

/******************************************************************************/
/*
  Purpose:

    SETCOMM determines the active communication directions.

  Discussion:

    In this picture, we're assuming the RESULTS array is split among 
    four processes, numbered 0 through 3 and arranged as suggested by the 
    following:

        0  |  1
     ------+-------
        2  |  3

    Then process 0 must communicate with processes 1 and 2 only,
    so its COMM array would be { 0, 1, 1, 0 }.

  Parameters:

    Input, int RANK, the rank of the process.

    Output, int COMM[4], contains a 0 (no) or 1 (yes) if
    communication is needed for the UP(0), RIGHT(1), DOWN(2)
    and LEFT(3) neighbors.
*/
{
  int i;
/*
  Start out by assuming all four neighbors exist.
*/
  for ( i = 0; i < 4; i++ ) 
  {
    comm[i] = 1;
  }
/*
  Up neighbor?
*/
  if ( rank < nblock )
  {
    comm[0] = 0;    
  }
/*
  Right neighbor?
*/
  if ( ( rank + 1 ) % nblock == 0 )
  {
    comm[1] = 0;
  }
/*
  Down neighbor?
*/
  if ( rank > (nblock*(nblock-1)-1) )
  {
    comm[2] = 0;
  }
/*
  Left neighbor?
*/
  if ( ( rank % nblock ) == 0 )
  {
    comm[3] = 0;
  }

  return;
}
/******************************************************************************/

void setex ( double ex[], double M[][nodeedge+2], int which )

/******************************************************************************/
/*
  Purpose:

    SETEX pulls off the edge values of M to send to another task.

  Parameters:

    Output, double EX[NODEEDGE], the values to be exchanged.

    Input, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process. 

    Input, int WHICH, 0, 1, 2, or 3, indicates the edge from which
    the data is to be copied.
*/                  
{
  int i;

  switch ( which ) 
  {
    case 0:
    {
      for ( i = 1; i <= nodeedge; i++) 
      {
        ex[i-1] = M[1][i];
      }
      break;
    }
    case 1:
    {
      for ( i = 1; i <= nodeedge; i++)
      {
        ex[i-1] = M[i][nodeedge];
      }
      break;
    }
    case 2:
    {
      for ( i = 1; i <= nodeedge; i++)
      {
        ex[i-1] = M[nodeedge][i];
      }
      break;
    }
    case 3:
    {
      for ( i = 1; i <= nodeedge; i++)
      {
        ex[i-1] = M[i][1];
      }
      break;
    }
  }
  return;
}
/******************************************************************************/

void unpack ( double M[][nodeedge+2], int where, double in[] )

/******************************************************************************/
/*
  Purpose:

    UNPACK puts the vector of new edge values into the edges of M.

  Parameters:

    Output, double M[nodeedge+2][nodeedge+2], the part of the results 
    kept by this process.

    Input, int WHERE, 0, 1, 2, or 3, indicates the edge to which the 
    data is to be applied.

    Input, int IN[nodeedge], the boundary data.
*/
{
  int i;

  if ( where == 0 )
  {
    for ( i = 0; i < nodeedge; i++ )
    {
      M[0][i+1] = in[i]; 
    }
  }
  else if ( where == 1 )
  {
    for ( i = 0; i < nodeedge; i++ )
    {
      M[i+1][nodeedge+1] = in[i];
    }
  }
  else if ( where == 2 )
  {
    for ( i = 0; i < nodeedge; i++ )
    {
      M[nodeedge+1][i+1] = in[i];
    }
  }
  else if ( where == 3 )
  {
    for ( i = 0; i < nodeedge; i++ )
    {
      M[i+1][0] = in[i];
    }
  }

  return;
}

