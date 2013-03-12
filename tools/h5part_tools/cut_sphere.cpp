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

   if (not (argc == 5))
   {
      std::cerr <<
      "usage: cut_sphere <input dump> <output dump> <rmax> <vmax>\n";
      return(1);
   }


   std::string iname = argv[1];
   std::string oname = argv[2];

   std::istringstream rmaxStr(argv[3]);
   fType rmax;
   rmaxStr >> rmax;

   std::istringstream vmaxStr(argv[4]);
   fType vmax;
   vmaxStr >> vmax;

   // load the particles
   parts.loadHDF5(iname);
   const size_t nop = parts.getNop();

   std::cerr << "loaded " << nop << " particles\n";

   vect3dT comvel, compos;
   fType   mtot;
   addToCOM(parts, compos, comvel, mtot);

   plistT rmvlst;

   for (size_t i = 0; i < nop; i++)
   {
      const vect3dT rvec = compos - parts[i].pos;
      const fType   r    = sqrt(dot(rvec, rvec));
      
      const vect3dT vvec = comvel - parts[i].vel;
      const fType   v    = sqrt(fabs(dot(vvec, vvec)));

      if ( (r > rmax) or (v > vmax) )
         rmvlst.push_back(i);

   }

   plistT::const_reverse_iterator rmvItr;
   for (rmvItr = rmvlst.rbegin(); rmvItr != rmvlst.rend(); rmvItr++)
      partT popP = parts.pop(*rmvItr);

   const size_t rop = parts.getNop();
   
   for (size_t i = 0; i < rop; i++)
      parts[i].id = i;

   std::cerr << "keep   " << rop << " particles\n";

   parts.doublePrecOut();
   parts.saveHDF5(oname);
   std::cerr << "stored dump\n";

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}
