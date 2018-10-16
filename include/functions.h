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

void InitializeField(double *phi, int row, int col){
for(int i = 0; i<row; ++i){
 for(int j =0; j<col; ++j){
  int k = i*col+j;
   phi[k]=0.0;
    }
   }
}

 void InitializeArray(double *phi,int n){
   for(int i =0; i<n; i++){
     phi[i]=0.0;
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

void vertexInterpolate(double* vertPhi, double* cellPhi, int vertRow, int vertCol, int cellCol){
// Inner Points
for(int i = 1; i<(vertRow-1); ++i){
 for(int j =1; j<(vertCol-1); ++j){

double Ta = cellPhi[(i-1)*cellCol+(j-1)]; // NW
double Tb = cellPhi[(i-1)*cellCol+j];    // NE
double Tc = cellPhi[(i*cellCol)+(j-1)];      // SW
double Td = cellPhi[i*cellCol+j];      //SE

vertPhi[i*vertCol+j]=0.25*(Ta+Tb+Tc+Td);
   }
 }
// North boundary
// South boundary
// East boundary
// West boundary

}

void vfaceInterpolate(double* vfacePhi, double* vertPhi,
                      double* cellPhi, int vfaceRow, int vfaceCol,
                      int vertCol, int cellCol){

for(int i = 1; i<(vfaceRow-1); ++i){
 for(int j =1; j<(vfaceCol-1); ++j){

double Ta = cellPhi[i*cellCol+(j-1)]; // West cell
double Tb = cellPhi[i*cellCol+j];     // East cell
double Tc = vertPhi[i*vertCol+j];   // North Vertex
double Td = vertPhi[(i+1)*vertCol+j];         //South Vertex


vfacePhi[i*vfaceCol+j]=0.25*(Ta+Tb+Tc+Td);
  }

}
}

void hfaceInterpolate(double* hfacePhi, double* vertPhi,
                      double* cellPhi, int hfaceRow, int hfaceCol,
                      int vertCol, int cellCol){
for(int i = 1; i<(hfaceRow-1); ++i){
 for(int j =1; j<(hfaceCol-1); ++j){

double Ta = vertPhi[i*vertCol+j]; // West vertex
double Tb = vertPhi[i*vertCol+(j+1)];     // East vertex
double Tc = cellPhi[(i-1)*cellCol+j];     // North cell
double Td = cellPhi[i*cellCol+j]; //South cell

hfacePhi[i*hfaceCol+j]=0.25*(Ta+Tb+Tc+Td);
 }

}

}
void gradient(double* gradxPhi,double* gradyPhi,
              double* vfacePhi,double* hfacePhi,
			        int cellRow, int cellCol, double dx, double dy,
              int vfaceCol, int hfaceCol){
for(int i = 1; i<(cellRow-1); ++i){
 for(int j =1; j<(cellCol-1); ++j){


   double Phie = vfacePhi[i*vfaceCol+(j+1)];
   double Phiw = vfacePhi[i*vfaceCol+j];
   double Phin = hfacePhi[i*hfaceCol+j];
   double Phis = hfacePhi[(i+1)*hfaceCol+j];

   gradxPhi[i*cellCol+j] = (Phie-Phiw)/dx;
   gradyPhi[i*cellCol+j] = (Phin-Phis)/dy;
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

void Poisson(double* Phi, int row, int col,
             double dx,double dy,double* W,
             int nc){
double lam = 1.0;
double PhiOld[nc];
double L2oPhi = 1.0;
double L2Phi = 1.0;
double norm  = 1.0;
double normtarget  = 1.0e-8;;
int itr = 1;
while(itr<5000 && norm>normtarget){
storeOldValue(Phi,PhiOld,nc);
for(int i=1; i<(row-1); i++){
 for(int j=1; j<(col-1); j++){
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

   double R     = W[k]- AP*PhiP-AE*PhiE-AW*PhiW-AN*PhiN-AS*PhiS;
   double delPhi = R/AP;
   Phi[k] = Phi[k]+lam*delPhi; }
 }
if(itr==1) {
  L2norm(Phi,PhiOld, &L2oPhi,nc);
  L2norm(Phi,PhiOld, &L2Phi,nc);
 } else {
 L2norm(Phi,PhiOld, &L2Phi,nc);
 }

norm = L2Phi/L2oPhi;

itr = itr+1;
//cout<<"iteration no:"<<itr<<"\t"<<"P inner Residual:"<<norm<<"\t"<<endl;
  }
cout<<"no of inner itr for p:"<<itr<<"pnorm"<<norm<<endl;
}


void refPressure(double* p, int row, int col){
for(int i=0;i<(row*col);i++){
 p[i] = p[i]-p[col];
  }
 }
