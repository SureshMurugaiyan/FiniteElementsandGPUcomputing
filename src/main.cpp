/*---------------------*- C++ 2D Incompressible FLow -*-----------------------*
|  Solves the  2D incompressible Fluid Flow in 2D geometry                    |
|  User Input is input.h File                                                 |
|  This is the main file of the solver                                        |
*-----------------------------------------------------------------------------*/
#include <iostream>
#include "input.h"
using namespace std;
/*----------------------------------------------------------------------------*
|                    Function Declarations                                    |
*----------------------------------------------------------------------------*/
void solveSerial();
void solveParallel(int argc,char **argv);
/*----------------------------------------------------------------------------*
|                      Main Function                                          |
*----------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
switch(Mode){
  case 'S':
    cout<<"Running in Serial Mode"<<endl;
    solveSerial();                                    // Calling Serial Solver
    break;
  case 'P':
    cout<<"Running in Parrallel Mode"<<endl;
    solveParallel(argc,argv);                       // calling parallel Solver
    break;
  default :
    cout<<"Running in Serial Mode"<<endl;
    solveSerial();
}
return 0;
}
// * * * * * * * * * * END  OF PROGRAM * * * * * * * * * * * * * * * * * * * //
