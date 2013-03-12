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
#include "friend_particle.h"
class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::energyPart,
   public sphlatch::clumpPart,
   public sphlatch::friendPart,
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

      vars.push_back(storeVar(orbit, "orbit"));
      vars.push_back(storeVar(clumpid, "clumpid"));

      vars.push_back(storeVar(rho, "rho"));
#ifdef SPHLATCH_ANEOS
      vars.push_back(storeVar(mat, "mat"));
#endif
      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(friendid, "friendid"));
      //vars.push_back(storeVar(orbit, "orbit"));
      //vars.push_back(storeVar(ecc, "ecc"));
      return(vars);
   }
};

typedef particle                               partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;

#include "fof.cpp"
typedef sphlatch::FOF<partT>                   fofT;

#include "clump.cpp"
typedef sphlatch::Clump<partT>                 clumpT;

#include "bhtree.cpp"
#include "bhtree_worker_grav.cpp"
typedef sphlatch::BHTree                       treeT;
typedef sphlatch::fixThetaMAC                  macT;
typedef sphlatch::GravityWorker<macT, partT>   gravT;

// particles are global
partSetT parts;

int main(int argc, char* argv[])
{
   if (not ((argc == 6) || (argc == 7)))
   {
      std::cerr <<
      "usage: clump_finder <inputdump> <clumpsfile> <minmass> <rhomin> <hmult> (<numthreads>)\n";
      return(1);
   }

   std::string inFilename = argv[1];
   std::string clFilename = argv[2];

   std::istringstream mMassStr(argv[3]);
   fType minMass;
   mMassStr >> minMass;

   std::istringstream rhominStr(argv[4]);
   fType rhomin;
   rhominStr >> rhomin;

   std::istringstream hmultStr(argv[5]);
   fType hmult;
   hmultStr >> hmult;

   if (argc == 7)
   {
      std::istringstream threadStr(argv[6]);
      int numThreads;
      threadStr >> numThreads;
      omp_set_num_threads(numThreads);
   }

   // load the particles
   parts.loadHDF5(inFilename);
   const fType  G   = parts.attributes["gravconst"];
   const size_t nop = parts.getNop();
   const fType  ML  = 7.34770e25;
   const fType day = 86400;

   std::cerr << nop << " particles loaded\n";

   if (minMass < 0.)
      if (parts.attributes.count("mminclump") == 0)
      {
         std::cerr << " attribute \"mminclump\" not found in "
                   << inFilename << "!\n";
         return(1);
      }
   minMass = parts.attributes["mminclump"];

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
   std::cerr << "Tree() constructed\n";

   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   // get central clump
   clumpT cb;
   for (size_t i = 0; i < nop; i++)
   {
      if (parts[i].clumpid == 1)
         cb.addParticle(&(parts[i]));
   }
   cb.calcProperties();

   // get central clump mass fractions and determine surface density
   fvectT      mfrac   = cb.getMassFractions(32);
   const fType minfrac = 0.1;

   fType rho0 = 1.;
   if (mfrac(2) > minfrac)
      rho0 = 1.11;
   else if (mfrac(1) > minfrac)
      rho0 = 2.65;
   else if (mfrac(4) > minfrac)
      rho0 = 3.32;
   else if (mfrac(5) > minfrac)
      rho0 = 7.85;

   fType rhoClmpMin = 0.75 * rho0;
   std::cerr << "central body surface density: " << rhoClmpMin << "\n";

   // get central clump surface
   cb.getCentralBodyOrbits(rhoClmpMin, G);
   const fType rcb(cb.rclmp);
   const fType mcb(cb.m);
   const fType rhocb(cb.rho);

   std::cerr << "central body surface radius:  " << rcb << "\n";

   // now search friends of friends
   fofT fof(parts, Tree);
   fof.search(rhomin, hmult, minMass);
   const size_t noc = fof.getNop();
   std::cerr << noc - 1 << " clump(s) found!\n";
   
   // now get the clumps orbits
   for (size_t i = 2; i < noc; i++)
   {
      //clumps
      clumpT& cc(fof[i]);

      cc.calcProperties();
      cc.getCentralBodyOrbits(rhoClmpMin, G);

      const vect3dT r0v = cc.pos - cb.posclmp;
      const vect3dT v0v = cc.vel - cb.vel;

      const fType mi   = cc.m;
      const fType rhoi = cc.rho;

      const fType r0 = sqrt(dot(r0v, r0v));
      const fType v0 = sqrt(dot(v0v, v0v));

      const fType K  = (mcb + mi) * G;
      const fType k1 = r0 * v0 * v0 / K;

      const fType beta0  = asin(dot(r0v, v0v) / (r0 * v0));
      const fType theta0 = atan2(k1 * sin(beta0) * cos(beta0),
                                 (k1 * cos(beta0) * cos(beta0) - 1.));

      const vect3dT h = cross(r0v, v0v);

      const vect3dT ev = (cross(v0v, h) / K) - (r0v / r0);
      const fType   e  = sqrt(dot(ev, ev));
      const fType   a  = dot(h, h) / (K * (1. - e * e));
      const fType   T  = 2.*M_PI*pow( a, 3./2.)/sqrt(K);
      
      const fType   rp = r0 * ((1. + e * cos(theta0)) / (1. + e));
      const fType   ra = rp * (1. + e) / (1. - e);


      const fType rrche = 2.455* rcb* sqrt(rhocb / rhoi);
      std::cout << "FOF clump " << i << "\n"
                << "  m: " << mi << "   rho: " << cc.rho << "   rRoche: "
                << rrche << "\n"
                << "  e: " << e << "    rp: " << rp << "       ra: "
                << ra << "\n"
                << "  a: " << a << "   T: " << T /day << " days\n"
                << "  rp/rcb: " << rp / rcb << "    rp/rRoche: "
                << rp / rrche << "  m/mL: " << mi / ML << "\n"
                << "\n"
                << "  SiO2: " << cc.mmat[1]
                << "   (" << 100. * cc.mtarg[1] / cc.mmat[1] << " % targ)\n"
                << "  Fe:   " << cc.mmat[5]
                << "   (" << 100. * cc.mtarg[5] / cc.mmat[5] << " % targ)\n"
                << "  H2O:  " << cc.mmat[2]
                << "   (" << 100. * cc.mtarg[2] / cc.mmat[2] << " % targ)\n"
                << "\n"
                << "\n";
   }
   
   if (noc > 1)
     fof[1].getCentralBodyOrbits(rhoClmpMin, G);

   std::cerr << "storing clump IDs in " << inFilename << "\n";
   parts.saveHDF5(inFilename);

   if (sphlatch::fileExists(clFilename))
   {
      h5FT        clFile(clFilename);
      std::string grpName = fof.getStepName();
      if (clFile.groupExists(grpName))
      {
         std::cerr << "del old clumps    in " << clFilename << "\n";
         clFile.deleteGroup(grpName);
      }
   }
   fof.doublePrecOut();
   std::cerr << "storing clumps    in " << clFilename << "\n";
   fof.saveHDF5(clFilename);

   return(0);
}
