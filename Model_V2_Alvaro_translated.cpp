using namespace std;
#include <vector>
#include <math.h>       /* atan2 */
#include <gsl/gsl_rng.h>  // Uncomment when gsl works
#include <gsl/gsl_randist.h>
#include <iostream>     // std::cout
#include <fstream> // create files
#include <time.h>
#define PI 3.1415926535897932384626433832795


// create a rng in order to create random numbers
gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);


// SET SYSTEM PARAMETERS HERE
const int MAXSTEP = 2*365*24/10;  // Number of time steps to be taken
const int N0 = pow(10,1);         // Number of initial cells
const int Alt = 2;                // Number of alterations
const int L = pow(2,Alt);         // Number of different possible genotypes
const int Xmax = 21;              // Number of voxels per x-side
const int Ymax = 21;              // Number of voxels per y-side
const int Zmax = 1;               // Number of voxels per z-side
const int K = pow(10,6);          // Maximum number of cells allowed to be in a voxel


// DECLARE FUNCTIONS HERE

// Function to compute population increase due to reproduction
int rep(int Pgen, int Ptot);
// Function to compute population decrease due to cell death
int kill(int Pgen, int Ptot);
// Function to compute new genotypes due to mutation
int mut(int Popgen, int Poptot, int chosen);
// Function to compute how many cells will migrate from a given voxel,
// and where to migrate
int mig(int Popgen, int Poptot, int chosen);



// CREATE CLASSES
class submod {
  public:
    int N;            // Total population inside a voxel
    int clones[4];    // Vector of population for each genotype
    int existing[4];  // Vector indicating whether a given genotype exists or not
    int D;            // Number of necrotic cells
};


// MAIN CODE
int main(){

// Set a random seed each time the program is executed
gsl_rng_set (r, time(NULL));

// Create matrix of submod objects
submod submodel[Xmax][Ymax];  // This one stores each timestep information
submod subnew[Xmax][Ymax];    // This one is used for correct updating

// Initialize system with zeros
for (int i = 0; i < Xmax; i++){
  for (int j = 0; j < Ymax; j++){
    submodel[i][j].N = 0;
    submodel[i][j].D = 0;
    for (int e = 0; e < L; e++){
      submodel[i][j].clones[e] = 0;
      submodel[i][j].existing[e] = 0;
    }
  }
}

// PLACE INITIAL POPULATION
// Set coordinates for initial population to be placed
// We want them to be at the center (voxel [11,11]), so we
// use the previous index, as indexing starts at 0 in C++
int x0 = 10;
int y0 = 10;
submodel[x0][y0].clones[0] = N0;
submodel[x0][y0].existing[0] = 1;
for (int i = 0; i < Xmax; i++){
  for (int j = 0; j < Ymax; j++){
    for (int e = 0; e < L; e++){
      submodel[i][j].N = submodel[i][j].N+submodel[i][j].clones[e];
    }
  }
}

// INITIALIZE UPDATING SUBMODEL ARRAY
for (int i = 0; i < Xmax; i++){
  for (int j = 0; j < Ymax; j++){
    subnew[i][j].N = submodel[i][j].N;
    subnew[i][j].D = submodel[i][j].D;
    for (int e = 0; e < L; e++){
      subnew[i][j].clones[e] = submodel[i][j].clones[e];
      subnew[i][j].existing[e] = submodel[i][j].existing[e];
    }
  }
}


// Check that population at initial point is higher than 0
//cout << submodel[x0][y0].clones[0] << "\t" << submodel[x0][y0].N << endl;
//cout << subnew[x0][y0].clones[0] << "\t" << subnew[x0][y0].N << endl;


// CREATE VARIABLES TO KEEP TRACK OF SYSTEM EVOLUTION
int evo[MAXSTEP][2];  // evo is for storing system population
int area[MAXSTEP][2]; // area is for storing occuppied voxels
// Initialize them
for (int t = 0; t < MAXSTEP; t++){
  evo[t][0] = 0;
  evo[t][1] = 0;
  area[t][0] = 0;
  area[t][1] = 0;
}
for (int i = 0; i < Xmax; i++){
  for (int j = 0; j < Ymax; j++){
    evo[0][1] = evo[0][1]+submodel[i][j].N;
    area[0][1] = area[0][1]+submodel[i][j].existing[0];
  }
}

// Check if code above worked properly (evo should be 10, and area 1)
// cout << evo[0][1] << "\t" << area[0][1] << endl;

// CREATE VECTOR STATE VARIABLES FOR EACH VOXEL
int Pop[Xmax][Ymax];        // Total population at each voxel
int clones[Xmax][Ymax][L];   // Population for each genotype at each voxel
int existing[Xmax][Ymax][L]; // Array indicating if genotypes exists at each voxel
int dead[Xmax][Ymax];       // Necrotic cells at each voxels
// Initialize them
for (int i = 0; i < Xmax; i++){
  for (int j = 0; j < Ymax; j++){
    Pop[i][j] = submodel[i][j].N;
    dead[i][j] = 0;
    for (int e = 0; e < L; e++){
      clones[i][j][e] = submodel[i][j].clones[e];
      existing[i][j][e] = submodel[i][j].existing[e];
    }
  }
}


// TIME EVOLUTION
// Iterate along each time step
int Pborn;
int Pkill;
int Pnew;
for (int t = 0; t < MAXSTEP; t++){
  // Iterate along each voxel
  for (int i = 0; i < Xmax; i++){
    for (int j = 0; j < Ymax; j++){
      // Update submodel only if the voxel is populated
      if (submodel[i][j].N > 0){
        for (int e = 0; e < L; e++){
          Pnew = 0;
          // REPRODUCTION EVENT (CHECK PROPER FUNCTION)
          //------------------------------------------------------------------------------
          Pborn = rep(submodel[i][j].clones[e],submodel[i][j].N);
          Pnew = Pnew + Pborn;

          // MUTATION EVENT (CHECK PROPER FUNCTION)
          //------------------------------------------------------------------------------

          // MIGRATION EVENT (CHECK PROPER FUNCTION)
          //------------------------------------------------------------------------------
          // Laplacian variable (lap) will store, for each possible movement (4),
          // the indexes of the voxels available for migration, and their laplacian (3)
          double LAP[4][3];
          double MIG[3];
          double lap;
          double cumlap = 0;
          double Pmig = 0.0025;
          int index = 0;
          for (int movi = -1; movi < 2; movi++){
            for (int movj = -1; movj < 2; movj++){
              // Calculate neighbour spatial indexes
              int in = i+movi;
              int jn = j+movj;

              // Check that we do not get out of boundaries
              if (in < Xmax+1 && in > 0 && jn < Ymax+1 && jn > 0 && abs(movi)+abs(movj) != 0 && abs(movi)+abs(movj) != 2){
                lap = (submodel[i][j].N-submodel[in][jn].N)/submodel[i][j].N;
                if (lap < 0){
                  lap = 0;
                }
                cumlap = cumlap + lap;
                LAP[index][0] = in;
                LAP[index][1] = jn;
                LAP[index][2] = lap;
                //cout << movi << "\t" << movj << "\t" << "\t" << lap << endl;
                index = index + 1;
              }
            }
          }

          // Calculate cumulative probability of migratin to any given voxel available
          for (int l = 0; l < 4; l++){
            if (l == 0){
              LAP[l][2] = LAP[l][2]/cumlap;
            } else {
              LAP[l][2] = LAP[l-1][2] + LAP[l][2]/cumlap;
            }
            //cout << LAP[l][0]-i << "\t" << LAP[l][1]-j << "\t" << LAP[l][2] << endl;
          }

          // Select where to migrate
          double Rand =  gsl_rng_uniform (r);
          for (int l = 0; l < 4; l++){
            if (Rand < LAP[l][2]){
              MIG[0] = LAP[l][0];
              MIG[1] = LAP[l][1];
              MIG[2] = LAP[l][2];
              break;
            }
          }

          //cout << MIG[0]-i << "\t" << MIG[1]-j << "\t" << MIG[2] << endl;
          int xmig = int(MIG[0]);
          int ymig = int(MIG[1]);

          Rand =  gsl_rng_uniform (r);
          if (Rand < Pmig){
            subnew[i][j].N = subnew[i][j].N - 1;
            subnew[xmig][ymig].N = subnew[xmig][ymig].N + 1;
            subnew[i][j].clones[e] = subnew[i][j].clones[e] - 1;
            subnew[xmig][ymig].clones[e] = subnew[xmig][ymig].clones[e] + 1;
            subnew[xmig][ymig].existing[e] = 1;
          }

          // DEATH EVENT (CHECK PROPER FUNCTION)
          //------------------------------------------------------------------------------
          Pkill = kill(submodel[i][j].clones[e],submodel[i][j].N);
          Pnew = Pnew - Pkill;

          // UPDATE Population
          //------------------------------------------------------------------------------
          subnew[i][j].clones[e] = subnew[i][j].clones[e] + Pnew;
          subnew[i][j].N = subnew[i][j].N + Pnew;
          subnew[i][j].D = subnew[i][j].D + Pkill;
          //cout << "\t" << Pnew << endl;
        }
      }

      Pop[i][j] = subnew[i][j].N;

    }
  }


  // UPDATE VOXELS
  for (int i = 0; i < Xmax; i++){
    for (int j = 0; j < Ymax; j++){
      if (subnew[i][j].N > 0){
        submodel[i][j].N = subnew[i][j].N;
        submodel[i][j].D = subnew[i][j].D;
        for (int e = 0; e < L; e++){
          submodel[i][j].clones[e] = subnew[i][j].clones[e];
          submodel[i][j].existing[e] = subnew[i][j].existing[e];
        }
      }
    }
  }


  // Calculate total population and necrotic cells

  int poptot = 0;
  int deadtot = 0;
  int areatot = 0;
  for (int i = 0; i < Xmax; i++){
    for (int j = 0; j < Ymax; j++){
      poptot = poptot + submodel[i][j].N;
      deadtot = deadtot + submodel[i][j].D;
      areatot = areatot + submodel[i][j].existing[0];
    }
  }

  evo[t+1][0] = t;
  evo[t+1][1] = poptot;

  area[t+1][0] = t;
  area[t+1][1] = areatot;

  cout << evo[t][0] << "\t" << evo[t][1] << "\t" << area[t][1] << endl;

}

}

//##############################################################################
//                             MAIN FUNCTIONS
//##############################################################################

//------------------------------------------------------------------------------
//                              REPRODUCTION
//------------------------------------------------------------------------------
int rep(int Popgen, int Poptot){

  // Popgen stands for current genotype population, while
  // Poptot stands for total voxel population

  // Compute reproduction probability
  // Make it proportional to voxel occupancy
  double Prep = (1-Poptot/K);
  int Popnew = Popgen;

  for (int i = 0; i < 0.02*Popgen; i++){
    double Rand =  gsl_rng_uniform (r);
    if (Rand < Prep){
      Popgen = Popgen + 1;
    }
  }

  Popnew = Popgen-Popnew;
  return Popnew;

}

//------------------------------------------------------------------------------
//                                MUTATION
//------------------------------------------------------------------------------

int mut(int Popgen, int Poptot, int chosen){

  double Pmut = 0.1;
  int Popnew = Popgen;

}

//------------------------------------------------------------------------------
//                                 DEATH
//------------------------------------------------------------------------------

int kill(int Popgen, int Poptot){

  double Pkill = 0.01;
  int Popnew = Popgen;

  // Use only if current genotype exist
  if (Popgen > 0){
    for (int i = 0; i < Popgen; i++){
      double Rand =  gsl_rng_uniform (r);
      if (Rand < Pkill){
        Popgen = Popgen - 1;
      }
    }
  }

  Popnew = Popnew-Popgen;
  return Popnew;
  // Although Popnew is the population killed, it is positive
  // Beware of substracting it instead of adding to current population
}

//------------------------------------------------------------------------------
//                               MIGRATION
//------------------------------------------------------------------------------

int mig(int Popgen, int Poptot, int chosen){

  double Pmig = 0.1;
  int Popnew = Popgen;

}


//##############################################################################
//                          AUXILIARY FUNCTIONS
//##############################################################################

//------------------------------------------------------------------------------
//                     CONVERT FROM DECIMAL TO BINARY
//------------------------------------------------------------------------------
