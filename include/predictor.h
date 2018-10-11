#include <iostream>

void eulerPredictor(double* ux, double *uy,
                       int nxc, int nyc,double dt, double dV,
                        double* Cnx, double *Cny, double* Dnx,
                        double* Dny, double Re){
for(int i = 1; i<(nyc-1); ++i){
    for(int j =1; j<(nxc-1); ++j){
        ux[i*nxc+j] = ux[i*nxc+j]+(dt/dV)*(-Cnx[i*nxc+j]+(Dnx[i*nxc+j])/Re);
        uy[i*nxc+j] = uy[i*nxc+j]+(dt/dV)*(-Cny[i*nxc+j]+(Dny[i*nxc+j])/Re);
    }
}

}


void adamPredictor(double* ux, double *uy,
                      int nxc, int nyc,
                      double dt, double dV,
                      double* Cnx, double *Cny,
                      double* Dnx,double* Dny,
                      double* CnxOld, double *CnyOld,
                      double* DnxOld,double* DnyOld,
                      double Re){
for(int i = 1; i<(nyc-1); ++i){
    for(int j =1; j<(nxc-1); ++j){
        ux[i*nxc+j] = ux[i*nxc+j]+(dt/dV)*(
                      1.5*(-Cnx[i*nxc+j]+(Dnx[i*nxc+j])/Re)
                     -0.5*(-CnxOld[i*nxc+j]+(DnxOld[i*nxc+j])/Re)
                      );
        uy[i*nxc+j] = uy[i*nxc+j]+(dt/dV)*(
                      1.5*(-Cny[i*nxc+j]+(Dny[i*nxc+j])/Re)
                     -0.5*(-CnyOld[i*nxc+j]+(DnyOld[i*nxc+j])/Re)
                      );
    }
}
}

