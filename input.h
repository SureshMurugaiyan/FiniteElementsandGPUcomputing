//USER INPUT
#define Mode 'P'       // Choice S == Serial Mode and Choice P == Parallel mode
#define nxc 32         // Number of cells in North and South Boundary
#define nyc 32         // Number of cells in East  and West  Boundary
#define Re  100.0      // Number of cells in East  and West  Boundary
#define pRefCell 1     // Pressure Reference cell // Assume center of cell
#define MAXitr 10     // Maximum number of Main Loop iterations// Time marching iterations
#define MAXnormux (1e-8)     // Desired normalized L2 norm for U-Velocity
#define MAXnormuy (1e-8)     // Desired normalized L2 norm for V-Velocity
#define MAXnormp (1e-8)      // Desired normalized L2 norm for Pressure
#define MAXpresitr 3000      // Maximum number of pressure iterations

//USER INPUT FOR PARALLEL
#define nxcL 16  // number of cells in north and south of each chunk --- chunknxc
#define nycL 16  // number of cells in each and west boundary of each chunk------chunknyc


#define nvar 3         // u,v,p // dont change it // this in only 2D solver
