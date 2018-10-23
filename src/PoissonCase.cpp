# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include "mpi.h"

// Total number of cells in X-Direction and same as Y-Direction = Matrix = 32 X 32
# define n 32
// Total number of Local cells in each Block in X-Direction and same as Y-Direction = Matrix = 16 X 16
# define chunkCell 16
// Total number of blocks in each direction Matrix = 2 X 2 = 4 blocks
# define nblock n/chunkCell
// Total number of processor = 4 processes
# define nproc nblock*nblock

// Function Declarations for the Poisson Equation

int main ( int argc, char **argv );
void doblock1 ( double w, double M[][chunkCell+2] );
void doblock2 ( double w, double M[][chunkCell+2] );
void exchange ( double M[][chunkCell+2], int comm[], int rank );
void iterate ( double w, double M[][chunkCell+2], double result[][n], int rank, int comm[] );
void setcomm ( int rank, int comm[] );
void setex ( double ex[], double M[][chunkCell+2], int which );
void initialize_matrix ( double M[][chunkCell+2] );
void unpack ( double M[][chunkCell+2], int where, double in[] );

// Main function
/******************************************************************************/
int main ( int argc, char **argv )
/******************************************************************************/

{
  int comm[4];
  FILE *fp;
  int i;
  int j;
  double M[chunkCell+2][chunkCell+2];
  int ntasks;
  int rank;
  double result[n][n];
  double w;
  double wtime;

  MPI_Init ( &argc, &argv );

  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

  MPI_Comm_size ( MPI_COMM_WORLD, &ntasks );

  wtime = MPI_Wtime ( );

  if ( rank == 0 )
  {
    printf ( "\n" );
    printf ( "2D_LAPLACE_SOLVER_MPI:\n" );
    printf ( "  C/MPI version\n" );
  }

  if ( ntasks != nproc )
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

  if ( rank == 0 )
  {
    printf ( "\n" );
    printf ( "  MPI has been set up.\n" );
  }
/*
  Initialize the matrix M.
*/
  if ( rank == 0 )
  {
    printf ( "  Initialize the matrix M.\n" );
  }
  initialize_matrix ( M );
/*
  Figure out who I communicate with.
*/
  if ( rank == 0 )
  {
    printf ( "  Set the list of neighbors.\n" );
  }
  setcomm ( rank, comm );
/*
  Update M, using SOR value W, until convergence.
*/
  if ( rank == 0 )
  {
    printf ( "  Begin the iteration.\n" );
  }
  w = 1.2;
  iterate ( w, M, result, rank, comm );
/*
  Report timing
*/
  wtime = MPI_Wtime ( ) - wtime;

  printf ( "  Task %i took %6.3f seconds\n", rank, wtime );
/*
  Write the solution to a file.
*/
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

void doblock1 ( double w, double M[][chunkCell+2] )

/******************************************************************************/
/* doblock1 iterates on the upper right and lower left corners of my matrix. */
{
  int i;
  int j;
/*
  Upper right corner.
*/
  for ( i = 1; i <= chunkCell / 2; i++ )
  {
    for ( j = chunkCell / 2 + 1; j <= chunkCell; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
/*
  Lower left corner.
*/
  for ( i = chunkCell / 2 + 1; i <= chunkCell; i++ )
  {
    for ( j = 1; j <= chunkCell / 2; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
  return;
}
/******************************************************************************/

void doblock2 ( double w, double M[][chunkCell+2] )

/******************************************************************************/
/* doblock2 iterates on the upper left and lower right corners of my matrix*/
{
  int i;
  int j;
/*
  Upper left corner.
*/
  for ( i = 1; i <= chunkCell / 2; i++ )
  {
    for ( j = 1; j <= chunkCell / 2; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
/*
  Lower right corner.
*/
  for ( i = chunkCell / 2 + 1; i <= chunkCell; i++ )
  {
    for ( j = chunkCell / 2 + 1; j <= chunkCell; j++ )
    {
      M[i][j] = w / 4.0 * ( M[i-1][j] + M[i][j-1] + M[i+1][j] + M[i][j+1] )
        + ( 1.0 - w ) * M[i][j];
    }
  }
  return;
}
/******************************************************************************/

void exchange ( double M[][chunkCell+2], int comm[], int rank )

/******************************************************************************/
/* EXCHANGE trades edge values with up to four neighbors.*/
{
  double ex0[chunkCell];
  double ex1[chunkCell];
  double ex2[chunkCell];
  double ex3[chunkCell];
  int i;
  double in0[chunkCell];
  double in1[chunkCell];
  double in2[chunkCell];
  double in3[chunkCell];
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
    MPI_Irecv ( &in0, chunkCell, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[0] );
  }
/*
  Receive from RIGHT neighbor (1).
*/
  if ( comm[1] == 1 )
  {
    partner = rank + 1;
    tag = 1;
    MPI_Irecv ( &in1, chunkCell, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[1] );
  }
/*
  Receive from DOWN neighbor (2).
*/
  if ( comm[2] == 1 )
  {
    partner = rank + nblock;
    tag = 2;
    MPI_Irecv ( &in2, chunkCell, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
      &requests[2] );
  }
/*
  Receive from LEFT neighbor (3).
*/
  if ( comm[3] == 1 )
  {
    partner = rank - 1;
    tag = 3;
    MPI_Irecv ( &in3, chunkCell, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
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
    MPI_Isend ( &ex0, chunkCell, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
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
    MPI_Isend ( &ex1, chunkCell, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
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
    MPI_Isend ( &ex2, chunkCell, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
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
    MPI_Isend ( &ex3, chunkCell, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD,
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

void initialize_matrix ( double M[][chunkCell+2] )

/******************************************************************************/

{
  double avg;
  double bv[4];
  int i;
  int j;

  bv[0] = 50.0;
  bv[1] = 0.0;
  bv[2] = 0.0;
  bv[3] = 0.0;
/*
  Put the boundary values into M.
*/
  for ( i = 1; i <= chunkCell; i++ )
  {
    M[0][i] =          bv[0];
    M[i][chunkCell+1] = bv[1];
    M[chunkCell+1][i] = bv[2];
    M[i][0] =          bv[3];
  }
/*
  Set all interior values to be the average of the boundary values.
*/
  avg = ( bv[0] + bv[1] + bv[2] + bv[3] ) / 4.0;

  for ( i = 1; i <= chunkCell; i++ )
  {
    for ( j = 1; j <= chunkCell; j++ )
    {
      M[i][j] = avg;
    }
  }

  return;
}
/******************************************************************************/

void iterate ( double w, double M[][chunkCell+2], double result[][n], int rank,
  int comm[] )

/******************************************************************************/
/* Main solving loop*/
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
  double mold[chunkCell+2][chunkCell+2];
  double send[chunkCell][chunkCell];

  it = 0;
  done = 0;
  for ( i = 1; i <= chunkCell; i++ )
  {
    for ( j = 1; j <= chunkCell; j++ )
    {
      mold[i][j] = M[i][j];
    }
  }

  while ( done == 0 )
  {
    it++;
/*
  Exchange values with neighbors, update block2, exchange values
  with neighbors, update block1 */
    exchange ( M, comm, rank );
    doblock2 ( w, M );
    exchange ( M, comm, rank );
    doblock1 ( w, M );
/*
  Check for convergence every 20 iterations.
*/
    if ( 2000 < it )
    {
      done = 1;
    }

    if ( ( it % 20 == 0.0 ) && ( done != 1 ) )
    {
      diff = 0.0;
      for ( i = 1; i <= chunkCell; i++ )
      {
        for ( j = 1; j <= chunkCell; j++ )
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
      diff = diff / ( ( double ) ( chunkCell * chunkCell ) );
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
  for ( i = 0; i < chunkCell; i++ )
  {
    for ( j = 0; j < chunkCell; j++ )
    {
      send[i][j] = M[i+1][j+1];
    }
  }

  count = chunkCell * chunkCell;

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
        for ( i = k * chunkCell; i < ( k + 1 ) * chunkCell; i++ )
        {
          for ( j = l * chunkCell; j < ( l + 1 ) * chunkCell; j++ )
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
/*   determines the active communication directions.

        0  |  1
     ------+-------
        2  |  3

   Process 0 must communicate with processes 1 and 2 and its , COMM array would be { 0, 1, 1, 0 }.

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

void setex ( double ex[], double M[][chunkCell+2], int which )

/******************************************************************************/
/* SETEX pulls off the edge values of M to send to another task. */
{
  int i;

  switch ( which )
  {
    case 0:
    {
      for ( i = 1; i <= chunkCell; i++)
      {
        ex[i-1] = M[1][i];
      }
      break;
    }
    case 1:
    {
      for ( i = 1; i <= chunkCell; i++)
      {
        ex[i-1] = M[i][chunkCell];
      }
      break;
    }
    case 2:
    {
      for ( i = 1; i <= chunkCell; i++)
      {
        ex[i-1] = M[chunkCell][i];
      }
      break;
    }
    case 3:
    {
      for ( i = 1; i <= chunkCell; i++)
      {
        ex[i-1] = M[i][1];
      }
      break;
    }
  }
  return;
}
/******************************************************************************/

void unpack ( double M[][chunkCell+2], int where, double in[] )

/******************************************************************************/
/* UNPACK puts the vector of new edge values into the edges of M.
WHERE, 0, 1, 2, or 3, indicates the edge to which the
    data is to be applied.*/
{
  int i;

  if ( where == 0 )
  {
    for ( i = 0; i < chunkCell; i++ )
    {
      M[0][i+1] = in[i];
    }
  }
  else if ( where == 1 )
  {
    for ( i = 0; i < chunkCell; i++ )
    {
      M[i+1][chunkCell+1] = in[i];
    }
  }
  else if ( where == 2 )
  {
    for ( i = 0; i < chunkCell; i++ )
    {
      M[chunkCell+1][i+1] = in[i];
    }
  }
  else if ( where == 3 )
  {
    for ( i = 0; i < chunkCell; i++ )
    {
      M[i+1][0] = in[i];
    }
  }

  return;
}

