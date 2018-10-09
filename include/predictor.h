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

void adamPredictor(double* ux, double *uy,double* ux_star, double *uy_star,
                       int nxc, int nyc,double dt, double dV,
                        double* Cnx, double *Cny,double* Cnx_star, double *Cny_star, 
			double* Dnx,double* Dny, double* Dnx_star,double* Dny_star,double Re){
for(int i = 1; i<(nyc-1); ++i){
    for(int j =1; j<(nxc-1); ++j){
        ux_star[i*nxc+j] = ux[i*nxc+j];
        uy_star[i*nxc+j] = uy[i*nxc+j];

        ux[i*nxc+j] = ux[i*nxc+j]+(dt/dV)*(1.5*(-Cnx_star[i*nxc+j]+(Dnx_star[i*nxc+j])/Re)-0.5*(-Cnx[i*nxc+j]+(Dnx[i*nxc+j])/Re));
        uy[i*nxc+j] = uy[i*nxc+j]+(dt/dV)*(1.5*(-Cny_star[i*nxc+j]+(Dny_star[i*nxc+j])/Re)-0.5*(-Cny[i*nxc+j]+(Dny[i*nxc+j])/Re));
    }
}

}
