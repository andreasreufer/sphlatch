#include <iostream>
#include <sstream>
#include <vector>

//#define SPHLATCH_SINGLEPREC

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#include "typedefs.h"
typedef sphlatch::fType     fType;
typedef sphlatch::cType     cType;
typedef sphlatch::vect3dT   vect3dT;
typedef sphlatch::box3dT    box3dT;
typedef sphlatch::partsIndexListT   plistT;


const fType finf = sphlatch::fTypeInf;

#include "bhtree.cpp"
typedef sphlatch::BHTree    treeT;

///
/// define the particle we are using
///
#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "io_particle.h"

class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::energyPart,
   public sphlatch::IOPart,
   public sphlatch::ANEOSPart
{
public:
   ioVarLT getLoadVars()
   {
      ioVarLT vars;
      
      vars.push_back(storeVar(id,  "id"));
      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(rho, "rho"));
      vars.push_back(storeVar(mat, "mat"));
      vars.push_back(storeVar(u, "u"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(p, "p"));
      vars.push_back(storeVar(S, "S"));
      vars.push_back(storeVar(noneigh, "noneigh"));
      vars.push_back(storeVar(cost, "cost"));

      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(id,  "id"));
      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(rho, "rho"));
      vars.push_back(storeVar(mat, "mat"));
      vars.push_back(storeVar(u, "u"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(p, "p"));
      vars.push_back(storeVar(S, "S"));
      vars.push_back(storeVar(noneigh, "noneigh"));
      vars.push_back(storeVar(cost, "cost"));
      return(vars);
   }
};

typedef particle                       partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>   partSetT;

// particles are global
partSetT parts;

#include "misc_physics.cpp"
using sphlatch::addToCOM;


int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif

   if (not (argc == 3))
   {
      std::cerr <<
      "usage: h5part2ascii <input dump> <output dump>\n";
      return(1);
   }


   std::string iname = argv[1];
   std::string oname = argv[2];

   // load the particles
   parts.loadHDF5(iname);
   const size_t nop = parts.getNop();

   std::cerr << "loaded " << nop << " particles\n";

   std::fstream fout;
   fout.open(oname.c_str(), std::ios::out);


   for (size_t i = 0; i < nop; i++)
   {
     const partT part = parts[i];
     std::string bodstr = part.id > 2000000 ? "1" : "0";
     std::string matstr = "silc";
     if (part.mat == 5)
       matstr = "iron";

     fout << std::scientific \
       << part.pos[0] << "," \
       << part.pos[1] << "," \
       << part.pos[2] << "," \
       << part.vel[0] << "," \
       << part.vel[1] << "," \
       << part.vel[2] << "," \
       << matstr << bodstr << "," \
       << part.h << "\n";
   }

   std::cerr << "stored dump\n";

   fout.close();

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}
