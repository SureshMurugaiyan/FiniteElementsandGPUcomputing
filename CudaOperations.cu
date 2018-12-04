/*---------------------*- C++ 2D Incompressible FLow -*-----------------------*
|  Solves the  2D incompressible Fluid Flow in 2D geometry                    |
|  User Input is input.h File                                                 |
|  This is the main file of the solver                                        |
*-----------------------------------------------------------------------------*/
#include <stdlib.h>
#include<stdio.h>
#include <cuda.h>
#include <iostream>
//#include <cuda_runtime.h>
using namespace std;


//*--------------------------------------------------------------------------*/
//check the return value of a cuda function for success
void checkError(cudaError_t err)
{
    if (err != cudaSuccess)
    {
        //print a human readable error message
        std::cout << cudaGetErrorString(err) << std::endl;
        exit(-1);
    }
    
        if (err == cudaSuccess)
    {
        //print a human readable error message
       std::cout << "Successfully loaded cuda" << std::endl;
        exit(-1);
    }
}











// * * * * * * * * * * END  OF PROGRAM * * * * * * * * * * * * * * * * * * * //
