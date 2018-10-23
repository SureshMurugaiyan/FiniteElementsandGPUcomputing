// C++ Script to simulate 2D incompressible flow
#include <iostream>
#include "functions.h"
 void iterate2DNavier(double* uxL,double* uyL,double* pL, 
		      double* ux,double* uy, double* p, 
                      int ncL, int nxcL, int nycL,
                      double dx,double dy,double dz,double dV,double dt,
		      int rank,int* comm){
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
  double MM[n*n];    // ux   and M = nodeedge+2 X nodeedge+2
  double mold[nodeedge+2][nodeedge+2];
  double send[nodeedge][nodeedge];
  it = 0;
  done = 0;




double Dnx[ncL]; // allocating space for Diffusion term
double Dny[ncL]; // allocating space for Diffusion term
double Cnx[ncL]; // allocating space for Convection term
double Cny[ncL]; // allocating space for Convection term

double DnxOld[ncL]; // allocating space for Diffusion term
double DnyOld[ncL]; // allocating space for Diffusion term
double CnxOld[ncL]; // allocating space for Convection term
double CnyOld[ncL]; // allocating space for Convection term

double Cn[ncL]; //  allocating space for Convection source term in pressure poisson
double gradxP[ncL]; //  allocating space for gradient
double gradyP[ncL]; //  allocating space for gradient


//L2-norm- L2 norm at each time step
double L2ux = 1.0;
double L2uy = 1.0;
double L2uz = 1.0;


//L2-norm - Initial Correction
double L2oux = 1.0;
double L2ouy = 1.0;
double L2ouz = 1.0;


double normux = 1.0;
double normuy = 1.0;
double normuz = 1.0;


Laplacian(Dnx,uxL,nycL,nxcL,dx,dy);      //  Dnx Diffusion term
Laplacian(Dny,uyL,nycL,nxcL,dx,dy);      //  Dny Diffusion term
Div(Cnx,uxL,uxL,uyL,nycL,nxcL,dx,dy);    //  Cnx Convection term
Div(Cny,uyL,uxL,uyL,nycL,nxcL,dx,dy);    //  Cny Convection term 



  for ( i = 1; i <= nodeedge; i++ ){
    for ( j = 1; j <= nodeedge; j++ ){
      mold[i][j] = M[i][j];
    }
  }

  for ( i = 1; i <= nodeedge; i++ ){
    for ( j = 1; j <= nodeedge; j++ ){
      mold[i][j] = M[i][j];
    }
  }


  while ( done == 0 ){
    it++;


/* Exchange values with neighbors, then update each block */
    exchange ( M, comm, rank );
    dored ( w, M );

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








int k = 1;

while (normuy>(1e-6)){
//--------------Store old values-------------------//
storeOldValue(ux,uxOld,nc);
storeOldValue(uy,uyOld,nc);
storeOldValue(uz,uzOld,nc);
storeOldValue(p,pOld,nc);

storeOldValue(Dnx,DnxOld,nc);
storeOldValue(Dny,DnyOld,nc);
storeOldValue(Cnx,CnxOld,nc);
storeOldValue(Cny,CnyOld,nc);



//--------------- Predictor step--------------------------------//
 Laplacian(Dnx,ux,nyc,nxc,dx,dy);      //  Dnx Diffusion term
 Laplacian(Dny,uy,nyc,nxc,dx,dy);      //  Dny Diffusion term
 Div(Cnx,ux,ux,uy,nyc,nxc,dx,dy);      //  Cnx Convection term
 Div(Cny,uy,ux,uy,nyc,nxc,dx,dy);      //  Cny Convection term


//eulerPredictor(ux,uy,nxc,nyc,dt,dV,Cnx,Cny,Dnx,Dny,Re);
adamPredictor(ux,uy,nxc,nyc,dt,dV,Cnx,Cny,Dnx,Dny,CnxOld,CnyOld,DnxOld,DnyOld,Re);



//--------------Solve poisson eqn for Pressure ---------------//
Divergence(Cn,ux,uy,nyc,nxc,dx,dy); //  source term in poisson equation

updateCn(Cn,dt,nxc,nyc);

Poisson(p,nyc,nxc,dx,dy,Cn,nc);

//--------------- Corrector Step ----------------------

//gradient(gradxP,gradyP,vfaceP,hfaceP,nyc,nxc,dx,dy,nxVfc,nxHfc);
//printmatrix(gradxP,nyc, nxc);

//gradient(gradxP,gradyP,p,nyc,nxc,dx,dy);

corrector(ux,uy,nxc,nyc,dt,gradxP,gradyP);

//--------------- Update Boundary condition ----------------------/

updateBoundaryCondition(ux,uy,p,nxc,nyc,nc);

//--------------- Update reference Pressure ----------------------/

refPressure(p,nyc, nxc);
if(k==1){
L2norm(ux,uxOld, &L2oux,nc);
L2norm(uy,uyOld, &L2ouy,nc);
L2norm(uz,uzOld, &L2ouz,nc);

L2norm(ux,uxOld, &L2ux,nc);
L2norm(uy,uyOld, &L2uy,nc);
L2norm(uz,uzOld, &L2uz,nc);

} else {
L2norm(ux,uxOld, &L2ux,nc);
L2norm(uy,uyOld, &L2uy,nc);
L2norm(uz,uzOld, &L2uz,nc);

}

normux = L2ux/L2oux;
normuy = L2uy/L2ouy;
normuz = L2uz/L2ouz;


cout<<endl<<"Itr:"<<k<<"\t"<<"uRes:"<<normux<<"\t"<<"vRes:"<<normuy<<"\t"<<endl;
k = k+1;
}
printData(ux,uy,p,nxc,nyc);








}
