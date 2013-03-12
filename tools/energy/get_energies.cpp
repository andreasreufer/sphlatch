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

typedef sphlatch::fvectT    fvectT;

const fType finf = sphlatch::fTypeInf;

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
#ifdef SPHLATCH_ANEOS
   public sphlatch::ANEOSPart,
#endif
   public sphlatch::IOPart
{
public:

   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(u, "u"));
      vars.push_back(storeVar(id, "id"));

#ifndef SPHLATCH_CALC_POTENTIAL
      vars.push_back(storeVar(pot, "pot"));
#endif
      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;
#ifdef SPHLATCH_CALC_POTENTIAL
      vars.push_back(storeVar(pot, "pot"));
#endif
      return(vars);
   }
};

typedef particle                               partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;

#include "clump_finder.cpp"
typedef sphlatch::Clumps<partT>                clumpsT;

#ifdef SPHLATCH_CALC_POTENTIAL
 #include "bhtree.cpp"
 #include "bhtree_worker_grav.cpp"
typedef sphlatch::BHTree                       treeT;
typedef sphlatch::fixThetaMAC                  macT;
typedef sphlatch::GravityWorker<macT, partT>   gravT;
#endif

#include "hdf5_io.cpp"
typedef sphlatch::HDF5File H5FT;

// particles are global
partSetT parts;
clumpsT  clumps;

int main(int argc, char* argv[])
{
   if (not ((argc == 3) || (argc == 4)))
   {
      std::cerr << 
      "usage: get_energies <inputdump> <clumpsfile> (<numthreads>)\n";
      return(1);
   }

   std::string inFilename = argv[1];
   std::string clFilename = argv[2];

   if (argc == 4)
   {
      std::istringstream threadStr(argv[3]);
      int numThreads;
      threadStr >> numThreads;
      omp_set_num_threads(numThreads);
   }

   // load the particles
   parts.loadHDF5(inFilename);
   std::cerr << "particles loaded\n";
   const fType G    = parts.attributes["gravconst"];
   const size_t nop = parts.getNop();

#ifdef SPHLATCH_CALC_POTENTIAL
   treeT& Tree(treeT::instance());
   gravT  gravWorker(&Tree, G);

   Tree.setExtent(parts.getBox() * 1.1);

   const fType  costppart = 1. / nop;
   for (size_t i = 0; i < nop; i++)
   {
      parts[i].cost = costppart;
      Tree.insertPart(parts[i]);
   }
   Tree.update(0.8, 1.2);

   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

 #pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcPot(CZbottomLoc[i]);
   std::cerr << "potential calculated\n";
#endif

   fType Ethm = 0., Ekin = 0., Epot = 0.;
   for (size_t i = 0; i < nop; i++)
   {
      const fType   mi = parts[i].m;
      const fType   ui = parts[i].u;
      const vect3dT vi = parts[i].vel;
      const fType poti = parts[i].pot;

      Ekin += 0.5 * dot(vi, vi) * mi;
      Epot += 0.5 * poti * mi;
      Ethm += ui * mi;
   }

   const fType Etot = Ekin + Epot + Ethm;

   std::cerr << "Ekin =   " << Ekin << "\n"
             << "Epot =   " << Epot << "\n"
             << "U    =   " << Ethm << "\n"
             << "Etot =   " << Etot << "\n";

   parts.attributes["ekin"] = Ekin;
   parts.attributes["ethm"] = Ethm;
   parts.attributes["epot"] = Epot;
   parts.attributes["etot"] = Etot;

   std::cerr << "storing energies  in " << inFilename << "\n";
   parts.doublePrecOut();
   parts.saveHDF5(inFilename);

   H5FT clumpf(clFilename);
   clumpf.setNewRoot(parts.getStepName());
   std::cerr << "storing energies  in " << clFilename << "\n";
   clumpf.saveAttribute("ekin", Ekin);
   clumpf.saveAttribute("ethm", Ethm);
   clumpf.saveAttribute("epot", Epot);
   clumpf.saveAttribute("etot", Etot);

   return(0);
}

