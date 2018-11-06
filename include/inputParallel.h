#define nc  nxc*nyc  // total number of cells
#define ncL nxcL*nycL
#define nxL   (nxcL+1)      // Number of vertices in North/South Boundary
#define nyL   (nycL+1)      // Number of vertices in East/West   Boundary
#define nvL   (nxL*nyL)      // Number of vertices in total

#define nprocx nxc/nxcL  // number of chunks divided along row ---- nchunkx
#define nprocy nyc/nycL  // number of chunks divided along column---- nchunky
#define nproc nprocx*nprocy // total number of processors or taks = no of subdomains

#define nxcGL (nxcL+2)
#define nycGL (nycL+2)
#define ncGL (nxcGL)*(nycGL)

#define nxGL   (nxL+2)       // Number of vertices in North/South Boundary including 1 layer of Ghost cells
#define nyGL   (nyL+2)       // Number of vertices in East/West   Boundary including 1 layer of Ghost cells
#define nvGL   (nxGL*nyGL)    // Number of vertices in total including 1 layer of Ghost cells

#define nxHfcL  (nxL-1)       // Number of face center in North/South Boundary along Horizontal co-ordinate lines
#define nyHfcL  nyL           // Number of face center in East/West   Boundary along Horizontal co-ordinate lines
#define nHfcL   (nxHfcL*nyHfcL)// Number of face center in Total along Horizontal co-ordinate lines 

#define nxVfcL  nxL           // Number of face center in North/South Boundary along Vertical co-ordinate lines
#define nyVfcL  (nyL-1)       // Number of face center in East/West   Boundary along Vertical co-ordinate lines
#define nVfcL   (nxVfcL*nyVfcL)// Number of face center in Total along Vertical co-ordinate lines

#define nxHfcGL  (nxHfcL+2)       // No of facecenter in North/South  alg Horz codnate lines + 1 layer of Ghost cells
#define nyHfcGL  (nyHfcL+2)       // No of facecenter in East/West   alg Horz codnate lines + 1 layer of Ghost cells
#define nHfcGL   (nxHfcGL*nyHfcGL) // No of facecenter in Total alg Horz codnate lines + 1 layer of Ghost cells 

#define nxVfcGL  (nxVfcL+2)        // No of fcenter in North/South alg Vert codnate lines + 1 layr of Ghost cells
#define nyVfcGL  (nyVfcL+2)        // No of fcenter in East/West alg Vert codnate lines + 1 layr of Ghost cells
#define nVfcGL   (nxVfcGL*nyVfcGL)  // No of fcenter in Total alg Vert codnate lines + 1 layr of Ghost cells

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

