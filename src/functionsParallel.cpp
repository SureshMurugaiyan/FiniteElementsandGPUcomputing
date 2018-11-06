#include "input.h"
#include "inputParallel.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Update Physical Boundary conditions for Parallel solver-North Boundary   !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void updatePhysicalNorthBC(double* uxL,double* uyL,double* pL,int ProcRank){
if(ProcRank<nprocx){
    for (int i = 0; i<nxcGL; i++){
        uxL[i]= 1;
        uyL[i]= 0;
        pL[i]  = pL[i+nxcGL];
    }
}
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Update Physical Boundary conditions for Parallel solver-South Boundary   !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void updatePhysicalSouthBC(double* uxL,double* uyL,double* pL,int ProcRank){
if(ProcRank>(nprocx*(nprocy-1)-1)){
    for (int i = ncGL-nxcGL; i<ncGL; i++){
        uxL[i]= 0;
        uyL[i]= 0;
        pL[i]= pL[i-nxcGL];
    }
}
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Update Physical Boundary conditions for Parallel solver-East  Boundary   !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void  updatePhysicalEastBC(double* uxL,double* uyL,double* pL,int ProcRank){
if((ProcRank+1)%nprocx==0){
    for (int i = nxcGL-1; i<ncGL; i=(i+nxcGL)){
        uxL[i]=0;
        uyL[i]=0;
        pL[i]  = pL[i-1];
    }
}
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Update Physical Boundary conditions for Parallel solver-West  Boundary   !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void  updatePhysicalWestBC(double* uxL,double* uyL,double* pL,int ProcRank){
if((ProcRank%nprocx)==0){
    for (int i = 0; i<ncGL; i=(i+nxcGL)){
        uxL[i]= 0;
        uyL[i]= 0;
        pL[i] = pL[i+1];
    }

}
}  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
