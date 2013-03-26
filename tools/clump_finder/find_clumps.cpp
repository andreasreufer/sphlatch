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
typedef sphlatch::vect3dT    vect3dT;
typedef sphlatch::box3dT     box3dT;

typedef sphlatch::fvectT     fvectT;

#include "hdf5_io.cpp"
typedef sphlatch::HDF5File   H5FT;

const fType finf = sphlatch::fTypeInf;

#include "misc_helpers.cpp"

///
/// define the particle we are using
///
#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "io_particle.h"
#include "clump_particle.h"
class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::energyPart,
   public sphlatch::clumpPart,
#ifdef SPHLATCH_ANEOS
   public sphlatch::ANEOSPart,
#endif
   public sphlatch::IOPart
{
public:
   fType morig;

   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(u, "u"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));

      vars.push_back(storeVar(rho, "rho"));
#ifdef SPHLATCH_ANEOS
      vars.push_back(storeVar(mat, "mat"));
#endif
#ifndef SPHLATCH_CALC_POTENTIAL
      vars.push_back(storeVar(pot, "pot"));
#endif
      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(clumpid, "clumpid"));
      vars.push_back(storeVar(orbit, "orbit"));
      vars.push_back(storeVar(ecc, "ecc"));
      vars.push_back(storeVar(a, "a"));
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

#include "bhtree.cpp"
#include "bhtree_worker_grav.cpp"
typedef sphlatch::BHTree                       treeT;
typedef sphlatch::fixThetaMAC                  macT;
typedef sphlatch::GravityWorker<macT, partT>   gravT;

// particles are global
partSetT parts;
clumpsT  clumps;

int main(int argc, char* argv[])
{
   if (not ((argc == 5) || (argc == 6)))
   {
      std::cerr <<
      "usage: clump_finder <inputdump> <clumpsfile> <minmass> <minmassorbits> (<numthreads>)\n";
      return(1);
   }

   std::string inFilename = argv[1];
   std::string clFilename = argv[2];

   std::istringstream mMassStr(argv[3]);
   fType minMass;
   mMassStr >> minMass;

   std::istringstream mMassOrbitsStr(argv[4]);
   fType minMassOrbits;
   mMassOrbitsStr >> minMassOrbits;

   if (argc == 6)
   {
      std::istringstream threadStr(argv[5]);
      int numThreads;
      threadStr >> numThreads;
      omp_set_num_threads(numThreads);
   }

   // load the particles
   parts.loadHDF5(inFilename);
   const size_t nop = parts.getNop();

   std::cerr << nop << " particles loaded\n";

   if (minMass < 0.)
      if (parts.attributes.count("mminclump") == 0)
      {
         std::cerr << " attribute \"mminclump\" not found in "
                   << inFilename << "!\n";
         return(1);
      }
      else
       minMass = parts.attributes["mminclump"];

   if (minMassOrbits < 0.)
      if (parts.attributes.count("mminorbit") == 0)
      {
         std::cerr << " attribute \"mminorbit\" not found in "
                   << inFilename << "!\n";
         return(1);
      }
      else
         minMassOrbits = parts.attributes["mminorbit"];

   treeT& Tree(treeT::instance());

   Tree.setExtent(parts.getBox() * 1.1);

   const fType costppart = 1. / nop;
   for (size_t i = 0; i < nop; i++)
   {
      parts[i].cost = costppart;
      Tree.insertPart(parts[i]);
   }
   Tree.update(0.8, 1.2);
   
   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

#ifdef SPHLATCH_GRAVITY
   const fType G = parts.attributes["gravconst"];
   gravT       gravWorker(&Tree, G);
 #pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcPot(CZbottomLoc[i]);
   std::cout << "calculated potential\n";


   fType Ethm = 0., Ekin = 0.;
 #ifdef SPHLATCH_GRAVITY
   fType Epot = 0.;
 #endif

   //
   for (size_t i = 0; i < nop; i++)
   {
      const fType   mi = parts[i].m;
      const fType   ui = parts[i].u;
      const vect3dT vi = parts[i].vel;

      Ekin += 0.5 * dot(vi, vi) * mi;
      Ethm += ui * mi;

 #ifdef SPHLATCH_GRAVITY
      const fType poti = parts[i].pot;
      Epot += 0.5 * poti * mi;
 #endif
   }


 #ifdef SPHLATCH_ESCAPEES
   const size_t noep = escapees.getNop();
   for (size_t i = 0; i < noep; i++)
   {
      const fType   mi = escapees[i].m;
      const fType   ui = escapees[i].u;
      const vect3dT vi = escapees[i].vel;

      Ekin += 0.5 * dot(vi, vi) * mi;
      Ethm += ui * mi;

  #ifdef SPHLATCH_GRAVITY
      const fType poti = escapees[i].pot;
      Epot += 0.5 * poti * mi;
  #endif
   }
 #endif

   parts.attributes["ekin"] = Ekin;
   parts.attributes["ethm"] = Ethm;
 #ifdef SPHLATCH_GRAVITY
   parts.attributes["epot"] = Epot;
 #endif


 #ifdef SPHLATCH_FIND_CLUMPS
   // set the minimal clump mass to 10 times that of the lightest particles
   fType pMinMass = 0., totMass = 0.;
   for (size_t i = 0; i < nop; i++)
   {
      pMinMass = parts[i].m > pMinMass ? parts[i].m : pMinMass;
      totMass += parts[i].m;
   }
   const fType cMinMass       = 10. * pMinMass;
   const fType cMinMassOrbits = 0.1 * totMass;
   const fType cMinRho        = parts.attributes["rhominclump"];

   parts.attributes["mminclump"] = cMinMass;
   parts.attributes["mminorbit"] = cMinMassOrbits;

   clumps.getClumps(parts, cMinMass);

   const size_t noc = clumps.getNop();
   std::cout << "found " << noc - 1 << " clump(s) with m > "
                 << cMinMass << "\n";

   std::fstream cfile;
   cfile.open("clumps.txt", std::ios::app | std::ios::out);
   cfile << std::setw(18) << std::setprecision(6) << std::scientific;
   cfile << time << "   ";
   for (size_t i = 0; i < noc; i++)
      cfile << clumps[i].m << " ";
   for (size_t i = noc; i < 10; i++)
      cfile << 0. << " ";
   cfile << "\n";
   cfile.close();

  #ifdef SPHLATCH_GRAVITY
   // store original mass
   for (size_t i = 0; i < nop; i++)
      parts[i].morig = parts[i].m;
  #endif

   for (size_t i = 1; i < noc; i++)
      if (clumps[i].m > cMinMassOrbits)
      {
         clumps[i].getCentralBodyOrbits(cMinRho, G);
         std::cout << "got orbits for clump " << i << "\n";

  #ifdef SPHLATCH_GRAVITY
         const int cid = i;
         for (size_t j = 0; j < nop; j++)
         {
            if (parts[j].clumpid == cid)
               parts[j].m = parts[j].morig;
            else
               parts[j].m = 0.;
            parts[j].treeNode->update();
         }

         Tree.redoMultipoles();
         std::cout << "    Tree.redoMultipoles()\n";

   #pragma omp parallel for firstprivate(gravWorker)
         for (int j = 0; j < noCZbottomLoc; j++)
            gravWorker.calcPot(CZbottomLoc[j]);
         std::cout << "    Tree.calcPot()\n";

         fType EpotCC = 0.;
         for (size_t j = 0; j < nop; j++)
            if (parts[j].clumpid == cid)
               EpotCC += 0.5 * parts[j].pot * parts[j].m;

         clumps[i].Epot = EpotCC;

         std::cout << "    Epot = " << EpotCC
                       << ", Erot = " << clumps[i].Erot
                       << ", Ekin = " << clumps[i].Ekin << "\n";
  #endif
      }

  #ifdef SPHLATCH_GRAVITY
   // restore original mass
   for (size_t i = 0; i < nop; i++)
      parts[i].m = parts[i].morig;
  #endif

   clumps.doublePrecOut();
   clumps.saveHDF5("clumps.h5part");
 #endif
#endif


#ifdef SPHLATCH_FIND_CLUMPS
   H5FT clumpf("clumps.h5part");
   clumpf.setNewRoot(parts.getStepName());

   clumpf.saveAttribute("ekin", Ekin);
   clumpf.saveAttribute("ethm", Ethm);
 #ifdef SPHLATCH_GRAVITY
   clumpf.saveAttribute("epot", Epot);
 #endif
#endif

   Tree.clear();
   std::cout << "Tree.clear() \n";

#ifdef SPHLATCH_LRDISK
   const fType rmin = parts.attributes["diskrmin"];
   const fType rmax = parts.attributes["diskrmax"];

   dbinT diskBinner(parts);
   diskBinner.saveBins(1., rmin, rmax, 100, "disk.hdf5");
   std::cout << "binned disk and stored to disk.hdf5\n";
#endif

   std::cerr << "storing clump IDs in " << inFilename << "\n";
   parts.saveHDF5(inFilename);
   
   if (sphlatch::fileExists(clFilename))
   {
      H5FT        clFile(clFilename);
      std::string grpName = clumps.getStepName();
      if (clFile.groupExists(grpName))
      {
         std::cerr << "del old clumps    in " << clFilename << "\n";
         clFile.deleteGroup(grpName);
      }
   }
   clumps.doublePrecOut();
   std::cerr << "storing clumps    in " << clFilename << "\n";
   clumps.saveHDF5(clFilename);

   return(0);
}



