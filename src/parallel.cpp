#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

double L=1.0;        /*Linear size of square region */
int N = 32;         /* number of interior points per dimension*/

double *u, *un_new      /* linear arrays to hold solution*/

/*macro to index into a 2-D (N+2) x(N+2) array */

#define INDEX(i,j) ((N+2)*(i)+(j))

int my_rank;    /*rank of this process */

int *proc            /* process indexed by vertex */
int *i_min, *i_max;  /* min, max vertex indices of processes */
int *left_proc, *right_proc;  /* processes to left and right */


/*
Functions:
*/

int main (int argc, char *argv[]);
void allocate_arrays ();
void jacobi ( int num_procs, double f[]);
void make_domains( int num_procs );
double *make_source();
void timestamp ();

/*****************************************************/
 int main(int argc, char *argv[])
 /****************************************************/

/*
 Purpose:
  MAIN is the main program for POISSON_MPI

  Discussion:

  This program solves Poisson's equation in a 2D region.

  The Jacobi iterative method is used to solve the linear system.

  MPI is used for parallel execution, with the domain divided into strips.

  Local Parameters:

  Local, double F[(N+2)x(N+2)], the source term.

  Local, int N, the number of interior vertices in one dimension.

  Local, int NUM_PROCS, the number of MPI processes

  Local, double U[(N+2)*(N+2)], a solution estimate

  Local, double U_NEW[(N+2)*(N+2)], a solution estimate.

  */
{

  double change;
  double epsilon = 1.0E-03;
  double *f;
  char file_name[100];
  int i;
  int j;

  double my_change;
  int my_n;
  int num_procs;
  int step;
  double *swap;
  double wall_time;

  /*  MPI Initialization */

MPI_Ini( &argc, &argv );

MPI_COMM_size( MPI_COMM_WORLD, &num_procs );
MPI_COMM_rank( MPI_COMM_WORLD, &my_rank);

/*
 Read commandline arguments, if present.
 */

if (1<argc)
{
  sscanf(argv[1], "%d", &N);
}
else
{
  N=32;
}
if(2<argc)
{
  sscanf( argv[2],"%lf", &epsilon );
}
else
{
  epsilon=1.0E-03;
}
if(3 <argc )
{
  strcpy (file_name, argv[3]);
}
else
{
  strcpy(file_name, "poisson_mpi.out");
}
/*
  print out initial information.
 */

if(my_rank==0)
{
  timestamp ();
  printf ( "\n" );
  printf ( "POISSON_MPI:\n" );
  printf ( " 2D Poisson equation using Jacobi algorithm\n" );
  printf("================================================\n");
  printf("MPI :1D domains, non-blocking send/receive\n");
  printf("Number of processes  = %d\n",num_procs );
  printf("Number of interior vertices = %d\n", N);
  printf("Desired accuracy=%f\n",epsilon );
  printf ("\n" );
}

allocate_arrays ();
f=make_source ();
make_domains ( num_procs );

step = 0;
/* Begin timing. */
wall_time = MPI_Wtime();
/* Begin iteration */
do
{
  jacobi (num_procs,f);
  ++step;
  /* Estimate error */

  change = 0.0;
  n =0 ;
  my_change = 0.0;
  my_n = 0;
for (i = i_min[my_rank]; i<=i_max[my_rank]; i++)
{
  for (j = 1; j<=N; j++)
  {
    if (u_new[INDEX(i,j)] !=0.0 )
    {
      my_change = my_change+fabs(1.0-u[INDEX(i,j)]/u_new[INDEX(i,j)]);

      my_n = my_n+1;

    }
  }
}

MPI_Allreduce( &my_change, &change, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce( &my_n, &n, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

if (n!=0)
{
  change = change / n;
}
if ( my_rank ==0 && (step % 10 ) == 0 )
{
  printf("N = %d, n = %d, my_n = %d, step %4d Error = %g\n", N,n,my_n,stp,change);
}

/* Interchange U and U_NEW. */

swap = u;
u = u_new;
u_new = swap;
}while (epsilon < change );

/* Copy the solution to process 0 and print to a file */

/* Report on wallclock time */

wall_time = MPI_Wtime()-wall_time;

if(my_rank ==0)
{
  printf( myrank ==0)
  {
    printf("\n");
    printf(" Wall clock time = %f secs\n", wall_time );
  }

  /* Terminate MPI here */

  MPI_Finalize ();

  free (f);

  /* Terminate */

  if (my_rank ==0)

  {
    printf ("\n");
    printf( "POISSON_MPI:\n");
    printf("Normal end of execution.\n");
    printf("\n");
    timestamp ();
  }
  return 0;
}

void allocate_arrays()
  /* ALLOCATE_ARRAYS creates and zeros out the arays U and U_NEW */

{
  int i;
  int ndof;
  ndof = (N+2)*(N+2);

  u = ( double * ) malloc ( ndof * sizeof (double ) );

  for ( i = 0 ; i < ndof; i++)
  {
    u[i] = 0.0;
  }

  u_new = (double * ) malloc ( ndof * sizeof ( double ));
  for ( i = 0; i< ndof; i++)
  {
    u_new[i] = 0.0;
  }
  return ;
}


























