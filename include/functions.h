#include <iostream>
using namespace std;

void printmatrix(double *matrix, int row, int col){
     for(int i=0; i<row; i++){
	for(int j=0; j<col; ++j){
		cout << matrix[i*col+j] << "\t";}
      cout<<endl;
     }
}

 void storeOldValue(double *phi, double *phiOld,int nc){
   for(int i =0; i<nc; i++){
     phiOld[i]=phi[i];
   }
 }

void L2norm(double *Phinew, double *Phiold,double *L2Phi,int  nc){
*L2Phi = 0;
  for(int i = 0; i<nc;i++){
    *L2Phi= *L2Phi+pow((Phiold[i]-Phinew[i]),2);
   }
   *L2Phi=sqrt(*L2Phi/nc);
}

void Laplacian(double* Ln, double *Phi, int row, int col, double dx, double dy){
for(int i = 1; i<(row-1); i++){
 for(int j =1; j<(col-1); j++){
   int k = i*col+j;
   double PhiP = Phi[k];
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];

   double Ee  = (PhiE-PhiP)/dx;
   double Ew  = (PhiP-PhiW)/dx;
   double Fn  = (PhiN-PhiP)/dy;
   double Fs  = (PhiP-PhiS)/dy;
   Ln[k]      = dx*(Fn-Fs)+dy*(Ee-Ew);
     }
  }
}

void vertexInterpolate(double* vertPhi, double* cellPhi, int row, int col){
for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){

int k = i*col+j;		// K = 7->21 AND  COL =6
double Ta = cellPhi[k-col-1];   // 7-6-1 = 0
double Tb = cellPhi[k-col];    // 7-6 = 1
double Tc = cellPhi[k-2];     // 7-2  = 5
double Td = cellPhi[k-1];     // 7-1  = 6

vertPhi[k]=0.25*(Ta+Tb+Tc+Td);
   }
 }
}

void hfaceInterpolate(double* hfacePhi, double* vertPhi, double* cellPhi, int row, int col){

for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){

double Ta = cellPhi[i*col+j-2];
double Tb = cellPhi[i*col+j-1];
double Tc = vertPhi[i*col+j];
double Td = vertPhi[i*col+j+col+1];

hfacePhi[i*col+j]=0.25*(Ta+Tb+Tc+Td);
  }
}
}

void vfaceInterpolate(double* vfacePhi, double* vertPhi, double* cellPhi, int row, int col){
for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){

double Ta = vertPhi[i*col+j+1];
double Tb = vertPhi[i*col+j+2];
double Tc = cellPhi[i*col+j-col];
double Td = cellPhi[i*col+j];

vfacePhi[i*col+j]=0.25*(Ta+Tb+Tc+Td);
 }

}

}
void gradient(double* gradxPhi,double* gradyPhi,double* hfacePhi,double* vfacePhi,
			int row, int col, double dx, double dy){
for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){

   int k = i*col+j;
   double Phie = hfacePhi[k+2];
   double Phiw = hfacePhi[k+1];
   double Phin = vfacePhi[k];
   double Phis = vfacePhi[k+col];

   gradxPhi[k] = (Phie-Phiw)/dx;
   gradyPhi[k] = (Phin-Phis)/dy;
    }
  }
}

void gradient(double* gradxPhi,double* gradyPhi,double* Phi,
                        int row, int col, double dx, double dy){
for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){

   int k = i*col+j;
   double Phie = Phi[k+1];
   double Phiw = Phi[k-1];
   double Phin = Phi[k-col];
   double Phis = Phi[k+col];

   gradxPhi[k] = (Phie-Phiw)/dx;
   gradyPhi[k] = (Phin-Phis)/dy;
    }
  }
}

void Div(double* Dn, double* Phi, double* U, double* V, int row, int col,double dx,double dy){


for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){
    int k = i*col+j;

   double PhiP = Phi[k];
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];

   double UP = U[k];
   double UE = U[k+1];
   double UW = U[k-1];
   double UN = U[k-col];
   double US = U[k+col];

   double VP = V[k];
   double VE = V[k+1];
   double VW = V[k-1];
   double VN = V[k-col];
   double VS = V[k+col];

   double Ee  = 0.5*(UE*PhiE+UP*PhiP);
   double Ew  = 0.5*(UW*PhiW+UP*PhiP);
   double Fn  = 0.5*(VN*PhiN+VP*PhiP);
   double Fs  = 0.5*(VS*PhiS+VP*PhiP);
   Dn[k]      = dx*(Fn-Fs)+dy*(Ee-Ew);

      }
   }
}

void Divergence(double* Dn, double* U, double* V,int row, int col, double dx, double dy){

for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){
    int k = i*col+j;
   double UP = U[k];
   double UE = U[k+1];
   double UW = U[k-1];
   double UN = U[k-col];
   double US = U[k+col];

   double VP = V[k];
   double VE = V[k+1];
   double VW = V[k-1];
   double VN = V[k-col];
   double VS = V[k+col];

   double Ue = 0.5*(UE+UP);
   double Uw = 0.5*(UW+UP);
   double Un = 0.5*(UN+UP);
   double Us = 0.5*(US+UP);

   double Ve = 0.5*(VE+VP);
   double Vw = 0.5*(VW+VP);
   double Vn = 0.5*(VN+VP);
   double Vs = 0.5*(VS+VP);

  Dn[k] = (Ue-Uw)*dy+(Vn-Vs)*dx;
   }
 }
}

void Poisson(double* Phi, int row, int col, double dx,double dy,double* W, int nc){
double lam = 1.0;
double PhiOld[nc];
double L2oPhi = 1.0;
double L2Phi = 1.0;
int k = 1;
while((L2Phi/L2oPhi)>(1e-3)){
storeOldValue(Phi,PhiOld,nc);
for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){
   int k = i*col+j;
   double PhiP = Phi[k];
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];
   double AP   = (-2*dy/dx)-(2*dx/dy);
   double AS   = (dx/dy);
   double AW   = (dy/dx);
   double AE   = (dy/dx);
   double AN   = (dx/dy);
   double vol  = dx*dy;
   double R     = vol*W[k]- AP*PhiP-AE*PhiE-AW*PhiW-AN*PhiN-AS*PhiS;
   double delPhi = R/AP;
   Phi[k] = Phi[k]+lam*delPhi;

      }
    }
if(k==1){
  L2norm(Phi,PhiOld, &L2oPhi,nc);
}

  L2norm(Phi,PhiOld, &L2Phi,nc);

//cout<<"iteration no:"<<k<<"\t"<<"P inner Residual:"<<L2Phi/L2oPhi<<"\t"<<endl;

k = k+1;
  }


}


void refPressure(double* p, int row, int col){
for(int i=0;i<(row*col);i++){
 p[i] = p[i]-p[col];
  }
 }
