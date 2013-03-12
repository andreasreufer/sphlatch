#include <iostream>
#include <sstream>
#include <vector>

//#define SPHLATCH_SINGLEPREC

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#include "typedefs.h"
typedef sphlatch::fType      fType;
typedef sphlatch::cType      cType;
typedef sphlatch::iType      iType;
typedef sphlatch::vect3dT    vect3dT;
typedef sphlatch::box3dT     box3dT;

typedef sphlatch::fvectT     fvectT;
typedef sphlatch::fmatrT     fmatrT;

const fType finf = sphlatch::fTypeInf;

#include "hdf5_io.cpp"
typedef sphlatch::HDF5File   HDF5File;

///
/// define the particle we are using
///
#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "io_particle.h"
#include "clump_particle.h"
#include "friend_particle.h"
class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::energyPart,
   public sphlatch::clumpPart,
   public sphlatch::ANEOSPart,
   public sphlatch::IOPart,
   public sphlatch::friendPart
{
public:

   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));

      vars.push_back(storeVar(clumpid, "clumpid"));
      vars.push_back(storeVar(orbit, "orbit"));
      vars.push_back(storeVar(ecc, "ecc"));
      vars.push_back(storeVar(pot, "pot"));

      vars.push_back(storeVar(rho, "rho"));
      vars.push_back(storeVar(mat, "mat"));

      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(friendid, "friendid"));
      return(vars);
   }

   iType clumpidsub;
};

typedef particle                               partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;

#include "clump_finder.cpp"
typedef sphlatch::Clumps<partT>                clumpsT;

#include "bhtree.cpp"
#include "bhtree_worker_grav.cpp"
typedef sphlatch::BHTree                       treeT;
typedef sphlatch::fixThetaMAC                  macT;
typedef sphlatch::GravityWorker<macT, partT>   gravT;

#include "disk_binner.cpp"
typedef sphlatch::DiskBinner<partT>            dbinT;

// particles are global
partSetT parts;
clumpsT  clumps;

int main(int argc, char* argv[])
{
   if (not ((argc == 7) || (argc == 10)))
   {
      std::cerr <<
      "usage: moon_postproc <inputdump> <diskfile> <clumpid> <rmin> <rmax> <nobins> (<rhoMin> <hMult> <mMin>)\n";
      return(1);
   }

   std::string inFilename = argv[1];
   std::string diskFilename = argv[2];

   std::istringstream ccidStr(argv[3]);
   iType ccid;
   ccidStr >> ccid;

   std::istringstream rminStr(argv[4]);
   fType rmin;
   rminStr >> rmin;

   std::istringstream rmaxStr(argv[5]);
   fType rmax;
   rmaxStr >> rmax;

   std::istringstream nobinsStr(argv[6]);
   size_t             nobins;
   nobinsStr >> nobins;

   fType rhomin, hmult, mmin;
   if (argc == 10)
   {
      std::istringstream rhominStr(argv[7]);
      rhominStr >> rhomin;

      std::istringstream hmultStr(argv[8]);
      hmultStr >> hmult;

      std::istringstream mminStr(argv[9]);
      mminStr >> mmin;
   }


   const size_t maxMat = 32;

   // load the particles
   parts.loadHDF5(inFilename);
   std::cerr << "particles loaded\n";
   const fType G = parts.attributes["gravconst"];

   dbinT diskBinner(parts);
   std::cerr << "bin disk no." << ccid << "\n";
   diskBinner.saveBins(ccid, rmin, rmax, nobins, diskFilename);

   if (argc == 10)
   {
      std::cerr << "looking for friends of friends ...\n";
      diskBinner.findFOF(ccid, rhomin, hmult, mmin);
      parts.saveHDF5(inFilename);
   }

   return(0);
}
