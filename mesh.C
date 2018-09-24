// C++ Script to simulate 2D incompressible flow 
#include <iostream>
#include <math.h>
#include <fstream>
#include "functions.h"
#include <cstring>



using namespace std;

int main()
{
  // North boundary
  const double xmin_n = -0.5;  // left  co-ordinate
  const double ymin_n =  0.5;  // left  co-ordinate
  const double xmax_n =  0.5;  // right co-ordinate
  const double ymax_n =  0.5;  // right co-ordinate

  //South boundary
  const double xmin_s = -0.5;  // left  co-ordinate
  const double ymin_s = -0.5;  // left  co-ordinate
  const double xmax_s =  0.5;  // right co-ordinate
  const double ymax_s = -0.5;  // right co-ordinate

  //East boundary
  const double xmin_e = xmax_s; // left  co-ordinate
  const double ymin_e = ymax_s; // left  co-ordinate
  const double xmax_e = xmax_n; // right co-ordinate
  const double ymax_e = ymax_n; // right co-ordinate

  //West boundary
  const double xmin_w = xmin_s; // left  co-ordinate
  const double ymin_w = ymin_s; // left  co-ordinate
  const double xmax_w = xmin_n; // right co-ordinate
  const double ymax_w = ymin_n; // right co-ordinate

  //Length of the domain
  const double length_n = sqrt(pow((xmin_n -xmax_n),2)+pow((ymin_n-ymax_n),2));
  const double length_s = sqrt(pow((xmin_s -xmax_s),2)+pow((ymin_s-ymax_s),2));
  const double length_e = sqrt(pow((xmin_e -xmax_e),2)+pow((ymin_e-ymax_e),2));
  const double length_w = sqrt(pow((xmin_w -xmax_w),2)+pow((ymin_w-ymax_w),2));

  //Evaluation of boundary lines      
  // Opposite boundaries should have equal number of points
  const int nxc = 5;       // number of cells in north and south boundary
  const int nyc = 5;       // number of cells in east and west boundary
  const int nx  = nxc+1;   // number of grid points in x-direction
  const int ny  = nyc+1;   // number of grid points in y-direction
  const int nc  = nxc*nyc; // total number of cells


  // Initializing matrix with grid points
  Matrix x(nx, ny);
  Matrix y(nx, ny);

  for (int i = 0; i<nx; i++)
  {
    x(0, i) = xmin_n+((xmax_n-xmin_n)/(nx-1))*i;   //north boundary grid points
    y(0,i) = ymin_n+((ymax_n-ymin_n)/(nx-1))*i;   // arranging from left to right
    x(ny-1,i) = xmin_s+((xmax_s-xmin_s)/(nx-1))*i;   //south boundary grid points
    y(ny-1,i) = ymin_s+((ymax_s-ymin_s)/(nx-1))*i;
  }

  for (int i = 0; i<ny; i++)
  {
    x(i,nx-1) = xmax_e-((xmax_e-xmin_e)/(ny-1))*i;   //east boundary grid points
    y(i,nx-1) = ymax_e-((ymax_e-ymin_e)/(ny-1))*i;   //arranging from top to bottom
    x(i,0) = xmax_w-((xmax_w-xmin_w)/(ny-1))*i;   //west boundary grid points
    y(i,0) = ymax_w-((ymax_w-ymin_w)/(ny-1))*i;
  }



  // Transfinite Interpolation to find inner points
  Matrix xnew = x; // x co-ordinates of grid
  Matrix ynew = y; // y co-ordinates of grid

  /* Computational domain is unit square
  coordinate lines are Psi and Eta
  Setting Psi = s ; and Eta = t
  Setting dPsi = dX ; dEta = dY */
  double s,t,dX,dY;

  dX = 1/(nx-1);   // Numerator is unity => unit square
  dY = 1/(ny-1);   // ,,

  // Interpolation of inner points by TFI

  for (int j = 1; j<nx-1; j++)
  {
    s = j*dX;
    for (int i = 1; i<ny-1; i++)
    {
      t = i*dY;
      xnew(j,i) = (1-t) * x(ny-1,i) + t * x(0,i) + 
                (1-s) * x(j,0) + s * x(j,nx-1) -
                (s*t*x(0,nx-1)+ s*(1-t)*x(ny-1,nx-1) +
                t*(1-s)*x(0,0) + (1-s)*(1-t)*x(ny-1,0));

      ynew(j, i)= (1-t)*y(ny-1,i) + t*y(0,i) +
                  (1-s)*y(j,0) + s*y(j,nx-1) -
                  (s*t*y(0,nx-1) + s*(1-t)*y(ny-1,nx-1) +
                  t*(1-s)*y(0,0) + (1-s)*(1-t)*y(ny-1,0));
    }
  }

  // Calculation of cell face centers
  Matrix x_hfc(ny,nx-1);   //x-co-ordinates of horizontal lines
  Matrix y_hfc(ny,nx-1);   //y-co-ordinates of horizontal lines

  Matrix x_vfc(ny-1,nx);   //x-co-ordinates of vertical lines
  Matrix y_vfc(ny-1,nx);   //y-co-ordinates of vertical lines

  Matrix x_c(ny-1,nx-1);   //x-co-ordinates of cell centers
  Matrix y_c(ny-1,nx-1);   //y-co-ordinates of cell centers

  // Calculation of horizontal face centers
  // center of horizontal grid lines

  for (int j = 0; j<ny; j++)
  {
    for (int i = 0; i<nx-1; i++)
    {
      x_hfc(j,i)=0.5*(xnew(j,i)+xnew(j,i+1));
      y_hfc(j,i)=0.5*(ynew(j,i)+ynew(j,i+1));
    }
  }

  // Calculation of vertical face centers
  // center of vertical grid lines

  for (int j = 0; j<ny-1; j++)
  {
    for (int i = 0; i<nx; i++)
    {
      x_vfc(j,i)=0.5*(xnew(j,i)+xnew(j+1,i));
      y_vfc(j,i)=0.5*(ynew(j,i)+ynew(j+1,i));
    }
  }


  // Calculation of cell centers

  for (int j = 0; j<ny-1; j++)
  {
    for (int i = 0; i<nx-1; i++)
    {
      x_c(j,i)=0.25*(x_hfc(j,i)+x_hfc(j+1,i)+x_vfc(j,i)+x_vfc(j,i+1));
      y_c(j,i)=0.25*(y_hfc(j,i)+y_hfc(j+1,i)+y_vfc(j,i)+y_vfc(j,i+1));
    }
  }

  /* Variables*/
  double dx = length_n/nxc ;
  double dy = length_e/nyc;



  ofstream myfile1;
  myfile1.open ("xVertices.txt");
  //myfile1 << "Writing this to a file.\n";
  for(int i=0; i<ny; i++)
  {
    for(int j=0; j<nx; ++j)
    {
      myfile1 << xnew(i,j) << "\t";
    }
    myfile1<<endl;
  }
  myfile1.close();

  ofstream myfile2;
  myfile2.open ("yVertices.txt");
  //myfile1 << "Writing this to a file.\n";
  for(int i=0; i<ny; i++)
  {
    for(int j=0; j<nx; ++j)
    {
      myfile2 << ynew(i,j) << "\t";
    }
    myfile2<<endl;
  }
  myfile2.close();

  ofstream myfile3;
  myfile3.open ("xCellcenter.txt");
  //myfile1 << "Writing this to a file.\n";
  for(int i=0; i<nxc; i++)
  {
    for(int j=0; j<nxc; j++)
    {
      myfile3<<x_c(i,j) << "\t";
    }
    myfile3<<endl;
  }
  myfile3.close();

  ofstream myfile4;
  myfile4.open ("yCellcenter.txt");
  //myfile1 << "Writing this to a file.\n";
  for(int i=0; i<nyc; i++)
  {
    for(int j=0; j<nxc; ++j)
    {
      myfile4<<y_c(i,j) << "\t";
    }
    myfile4<<endl;
  }
  myfile4.close();
 cout<<x_c(1,1)<<endl; 


Vector X(nx*ny);
Vector Y(nx*ny);
Vector Z(nx*ny);

 ofstream myfile5;
 myfile5.open ("Mesh.txt");
 myfile5 << "#MeshFile" << endl;  
   for(int i=0; i<ny; i++){
      for(int j=0; j<nx; j++){
        X(i*nx+j)= xnew(i,j);
        Y(i*nx+j)= ynew(i,j);
    myfile5<< (i*nx+j) <<"\t"<< X(i*nx+j)<<"\t"<<Y(i*nx+j)<<endl;
    
   } 
  }










return 0;
  //END OF PROGRAM

} 

