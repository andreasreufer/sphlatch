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
typedef sphlatch::HDF5File   h5FT;

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
   const fType  G   = parts.attributes["gravconst"];
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
   gravT  gravWorker(&Tree, G);

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

#ifdef SPHLATCH_CALC_POTENTIAL
 #pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcPot(CZbottomLoc[i]);
   std::cerr << "potential calculated\n";
#endif

   fType Ethm = 0., Ekin = 0., Epot = 0.;
   for (size_t i = 0; i < nop; i++)
   {
      const fType   mi   = parts[i].m;
      const fType   ui   = parts[i].u;
      const vect3dT vi   = parts[i].vel;
      const fType   poti = parts[i].pot;

      Ekin += 0.5 * dot(vi, vi) * mi;
      Epot += 0.5 * poti * mi;
      Ethm += ui * mi;
   }

   parts.attributes["ekin"] = Ekin;
   parts.attributes["ethm"] = Ethm;
   parts.attributes["epot"] = Epot;

   clumps.getClumps(parts, minMass);
   const size_t noc = clumps.getNop();
   std::cerr << noc - 1 << " clump(s) found!\n";

   for (size_t i = 0; i < nop; i++)
      parts[i].morig = parts[i].m;

   for (size_t i = 1; i < noc; i++)
      if (clumps[i].m > minMassOrbits)
      {
         fvectT mfrac = clumps[i].getMassFractions(32);

         fType rhoMin = 0.75;
         // if there is no rhominclump in attributes, use standard values below
         if (parts.attributes.count("rhominclump") == 0)
         {
#ifdef SPHLATCH_ANEOS
            // take rho0 from the lightes material with a mass
            // fraction of at least 10%
            const fType minfrac = 0.1;
            fType       rho0    = 1.;

            if (mfrac(2) > minfrac)
               rho0 = 1.11;
            else if (mfrac(1) > minfrac)
               rho0 = 2.65;
            else if (mfrac(4) > minfrac)
               rho0 = 3.32;
            else if (mfrac(5) > minfrac)
               rho0 = 7.85;

            rhoMin = 0.75 * rho0;
#endif
         }
         else
            rhoMin = parts.attributes["rhominclump"];

         /*rhoMin = 2.00;
         clumps[i].getCentralBodyOrbits(rhoMin, G);
         std::cerr << "got orbits for clump " << i
                   << " (rhoMin: " << rhoMin
                   << ", rho: "   << clumps[i].rho
                   << ", rclmp: " << clumps[i].rclmp
                   << ", rmean: " << clumps[i].rc << ")\n";*/

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
         std::cerr << "    Tree.redoMultipoles()\n";
         
         clumps[i].getCentralBodyOrbits(rhoMin, G);
         std::cerr << "got orbits for clump " << i
                   << " (rhoMin: " << rhoMin
                   << ", rho: "   << clumps[i].rho
                   << ", rclmp: " << clumps[i].rclmp
                   << ", rmean: " << clumps[i].rc << ")\n";


#pragma omp parallel for firstprivate(gravWorker)
         for (int j = 0; j < noCZbottomLoc; j++)
            gravWorker.calcPot(CZbottomLoc[j]);
         std::cerr << "    Tree.calcPot()\n";

         fType EpotCC = 0.;
         for (size_t j = 0; j < nop; j++)
            if (parts[j].clumpid == cid)
               EpotCC += 0.5 * parts[j].pot * parts[j].m;

         clumps[i].Epot = EpotCC;

         fType EpotCCCB = 0.;
         fType ErotCCCB = 0.;
         for (size_t j = 0; j < nop; j++)
         {
            if (parts[j].clumpid == cid and parts[j].orbit ==
                sphlatch::ORBITCLUMP)
            {
               const fType mj   = parts[j].m;
               const fType potj = parts[j].pot;

               const vect3dT rj  = parts[j].pos - clumps[i].posclmp;
               const vect3dT vj  = parts[j].vel - clumps[i].velclmp;
               const vect3dT rXv = cross(rj, vj);

               const fType rXvrXv = dot(rXv, rXv);
               const fType rr     = dot(rj, rj);
               ErotCCCB += 0.5 * mj * rXvrXv / rr;
               EpotCCCB += 0.5 * potj * mj;
            }
         }

         clumps[i].Epotclmp = EpotCCCB;
         clumps[i].Erotclmp = ErotCCCB;

         std::cerr << "    Epot = " << EpotCC
                   << ", Erot  = " << clumps[i].Erot
                   << ", EpotCB = " << clumps[i].Epotclmp
                   << ", ErotCB = " << clumps[i].Erotclmp
                   << ", Ekin   = " << clumps[i].Ekin << "\n";
      }
   
   std::cerr << "storing clump IDs in " << inFilename << "\n";
   parts.saveHDF5(inFilename);
   
   if (sphlatch::fileExists(clFilename))
   {
      h5FT        clFile(clFilename);
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
