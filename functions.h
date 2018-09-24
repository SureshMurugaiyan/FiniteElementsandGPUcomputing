#include <iostream>
#include "Matrix.h"
#include "Vector.h"
using namespace std;

Vector Laplacian( Vector &Phi, int row, int col, double dx, double dy)
{
  Vector Dn(row*col);
  for(int i = 1; i<(row-1); ++i)
  {
    for(int j =1; j<(col-1); ++j)
    {
      int k = i*col+j;
      Dn(k) = (Phi(k+1) - 2*Phi(k) + Phi(k-1))*dy/dx
            +(Phi(k-col) - 2*Phi(k) + Phi(k+col))*dx/dy;
    }
  }
  return Dn;
} 
void gradient(Vector& gradxPhi,Vector& gradyPhi,Vector& Phi,int row, int col, double dx, double dy)
{
  for(int i = 1; i<(row-1); ++i)
  {
    for(int j =1; j<(col-1); ++j)
    {
      gradxPhi(i*col+j) = (Phi(i*col+j+1)-Phi(i*col+j-1))/dx;
      gradyPhi(i*col+j) = (Phi(i*col+j-col)-Phi(i*col+j+col))/dy;
    }
  } 
}
void Div(Vector& Cnx, Vector& Phi, Vector& U, Vector& V, int row, int col,double dx,double dy){


for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){
      int k = i*col+j; 
         Cnx(k) = (0.5*(U(k+1)+U(k))*0.5*(Phi(k+1)+Phi(k))-0.5*(U(k-1)+U(k))*0.5*(Phi(k-1)+Phi(k)))*dy+\
                  (0.5*(V(k-col)+U(k))*0.5*(Phi(k-col)+Phi(k))-0.5*(V(k+col)+V(k))*0.5*(Phi(k+col)+Phi(k)))*dx;

      }
   }
}
 
void Divergence(Vector& Cn, Vector& U, Vector& V,int row, int col, double dx, double dy){

for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){
    int k = i*col+j; 
    Cn(k) = ((U(k+1)+U(k))/2-(U(k)+ U(k-1))/2)*dy
            + ((V(k-col)+V(k))/2-(V(k)+ V(k+col))/2)*dx;
   } 
 }
}

void Poisson(Vector& Phi, int row, int col, double dx,double dy,Vector& W){
double lam = 1;
int k = 1;
double A_P   = (-2*dy/dx)-(2*dx/dy);
double A_S   = -(dx/dy);
double A_W   = -(dy/dx);
double A_E   = -(dy/dx);
double A_N   = -(dx/dy);
while(k < 5000){
for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){
   double R     = W(i*col+j)- A_P*Phi(i*col+j)-A_E*Phi(i*col+j+1)-A_W*Phi(i*col+j-1)-A_N*Phi(i*col+j-col)-A_S*Phi(i*col+j+col);
   double delPhi = R/A_P;
   Phi(i*col+j) = Phi(i*col+j)+lam*delPhi;

      }
    }k = k+1;
  }
}

void refPressure(Vector& p, int row, int col){
for(int i=0;i<(row*col);i++){
 p(i) = p(i)-p(0);
  }
 }


