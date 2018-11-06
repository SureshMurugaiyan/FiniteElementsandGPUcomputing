#include "input.h"
#define nx   (nxc+1)      // Number of vertices in North/South Boundary
#define ny   (nyc+1)      // Number of vertices in East/West   Boundary
#define nv   (nx*ny)      // Number of vertices in total


#define nxcG (nxc+2)       // Number of cells in North/South Boundary including 1 layer of Ghost cells
#define nycG (nyc+2)       // Number of cells in East/West   Boundary including 1 layer of Ghost cells
#define ncG  (nxcG)*(nycG) // Number of cells in Total including 1 layer of Ghost cells


#define nxG   (nx+2)       // Number of vertices in North/South Boundary including 1 layer of Ghost cells
#define nyG   (ny+2)       // Number of vertices in East/West   Boundary including 1 layer of Ghost cells
#define nvG   (nxG*nyG)    // Number of vertices in total including 1 layer of Ghost cells

#define nxHfc  (nx-1)       // Number of face center in North/South Boundary along Horizontal co-ordinate lines
#define nyHfc  ny           // Number of face center in East/West   Boundary along Horizontal co-ordinate lines
#define nHfc   (nxHfc*nyHfc)// Number of face center in Total along Horizontal co-ordinate lines 

#define nxVfc  nx           // Number of face center in North/South Boundary along Vertical co-ordinate lines
#define nyVfc  (ny-1)       // Number of face center in East/West   Boundary along Vertical co-ordinate lines
#define nVfc   (nxVfc*nyVfc)// Number of face center in Total along Vertical co-ordinate lines


#define nxHfcG  (nxHfc+2)       // No of fcenter in North/South  alg Horz codnate lines + 1 layr of Ghost cells
#define nyHfcG  (nyHfc+2)       // No of fcenter in East/West   alg Horz codnate lines + 1 layr of Ghost cells
#define nHfcG   (nxHfcG*nyHfcG) // No of fcenter in Total alg Horz codnate lines + 1 layr of Ghost cells 

#define nxVfcG  (nxVfc+2)        // No of fcenter in North/South alg Vert codnate lines + 1 layr of Ghost cells
#define nyVfcG  (nyVfc+2)        // No of fcenter in East/West alg Vert codnate lines + 1 layr of Ghost cells
#define nVfcG   (nxVfcG*nyVfcG)  // No of fcenter in Total alg Vert codnate lines + 1 layr of Ghost cells

#define dx (1.0/nxc)
#define dy (1.0/nyc)
#define dV (dx*dy)

/*

//USER INPUT
#define Re 100.0

// Opposite boundaries should have equal number of points
#define nxc 32      // number of cells in north and south boundary
#define nyc 32      // number of cells in east and west boundary
#define nc  nxc*nyc  // total number of cells
#define chunknxc 16  // number of cells in north and south of each chunk
#define chunknyc 16  // number of cells in each and west boundary of each chunk
#define nchunkx nxc/chunknxc  // number of chunks divided along row
#define nchunky nyc/chunknyc  // number of chunks divided along column
#define nproc nchunkx*nchunky // total number of processors or taks = no of subdomains
#define ncL chunknxc*chunknyc
#define ncGL (chunknxc+2)*(chunknyc+2)
#define nxcL chunknxc
#define nycL chunknyc
#define nxcGL (chunknxc+2)
#define nycGL (chunknyc+2)

#define dz 1
#define dV dx*dy*dz
#define nx (nxc+1)
#define ny (nyc+1)
#define nv nx*ny

#define nxGL (nxcGL+1)
#define nyGL (nycGL+1)
#define nvGL (nxGL*nyGL)
#define nyHfcGL nyGL
#define nxHfcGL (nxGL-1)
#define nHfcGL  (nxHfcGL*nyHfcGL)
#define nyVfcGL (nyGL-1)
#define nxVfcGL nxGL
#define nVfcGL  (nxVfcGL*nyVfcGL)
*/

