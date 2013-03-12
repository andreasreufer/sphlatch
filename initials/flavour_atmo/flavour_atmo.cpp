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

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(u, "u"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));
      vars.push_back(storeVar(mat, "mat"));

      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(rho, "rho"));
#ifdef SPHLATCH_MISCIBLE
      vars.push_back(storeVar(delta, "delta"));
#endif
      vars.push_back(storeVar(mat, "mat"));
      vars.push_back(storeVar(u, "u"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(S, "S"));
      vars.push_back(storeVar(p, "p"));
      vars.push_back(storeVar(noneigh, "noneigh"));
      vars.push_back(storeVar(cost, "cost"));
      return(vars);
   }
};

typedef particle                       partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>   partSetT;
typedef std::list<partT>               partLT;

#include "sph_algorithms.cpp"
#include "sph_kernels.cpp"
#include "bhtree_worker_sphsum.cpp"

typedef sphlatch::CubicSpline3D                krnlT;

typedef sphlatch::densSum<partT, krnlT>        densT;
typedef sphlatch::SPHsumWorker<densT, partT>   densSumT;

#include "bhtree_worker_cost.cpp"
typedef sphlatch::CostWorker<partT>            costT;


#include "lookup_table1D.cpp"
typedef sphlatch::InterpolateLinear            intLinT;
typedef sphlatch::InterpolateStepwise          intStepT;

typedef sphlatch::LookupTable1D<intLinT>       LUT1DlinT;
typedef sphlatch::LookupTable1D<intStepT>      LUT1DstepT;

#include "hdf5_io.cpp"
typedef sphlatch::HDF5File                     H5FT;

#include "lattice_hcp.h"
typedef sphlatch::LatticeHCP                   lattT;

#include "eos_super.cpp"
typedef sphlatch::SuperEOS<partT>              eosT;

// particles are global
partSetT parts;


int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif
   
   if (not (argc == 6))
      {
      std::cerr <<
      "usage: flavour_atmo <indump> <outdump> <profile> <relgaspartmass> <rhocutoff>\n";
      return(1);
      }

   std::string idname = argv[1];
   std::string odname = argv[2];
   std::string pname = argv[3];

   std::istringstream rgmsstr(argv[4]);
   fType mrelgas;
   rgmsstr >> mrelgas;
   
   std::istringstream rhocosstr(argv[5]);
   fType rhoco;
   rhocosstr >> rhoco;

   const fType shellWidth = 0.2;
   const fType shelleps = 1.e-3;


   parts.loadHDF5(idname);
   const size_t nocp = parts.getNop();

   // offset center of mass for core particles
   vect3dT com;
   com = 0., 0., 0.;
   fType mcore = 0., hmean = 0.;
   for (size_t i = 0; i < nocp; i++)
   {
      com   += parts[i].pos * parts[i].m;
      mcore += parts[i].m;
      hmean += parts[i].h;
   }
   com /= mcore;
   hmean /= static_cast<fType>( nocp );

   std::cout << "center of mass: ["
             << com[0] << ","
             << com[1] << ","
             << com[2] << "]\n";
   
   const fType micore = mcore / static_cast<fType>( nocp );
   const fType migas = micore*mrelgas;
   
   for (size_t i = 0; i < nocp; i++)
      parts[i].pos -= com;

   fType rcore = -sphlatch::fTypeInf;
   for (size_t i = 0; i < nocp; i++)
   {
     const fType ri = sqrt( dot( parts[i].pos, parts[i].pos ) );
     rcore = rcore > ri ? rcore : ri;
   }
   
   std::cout << "mcore:  " << mcore << " rcore: " << rcore << " hmean: " << hmean << "\n";

   fType       lmed = 0.;

   LUT1DlinT rhoLUT(pname, "r", "rho");
   LUT1DlinT uLUT(pname, "r", "u");
   
   H5FT profFile(pname);
   const fType rminprof = rhoLUT.getXmin();
   const fType rmincore = rcore + 2.*hmean;

   fType rmin = std::max(rminprof, rmincore);
   fType rmed = rmin, rmax;

   std::cout << "rmin: " << rmin << "  rminProf: " << rminprof << "  rcore+2*hmean: " << rmincore << "\n";

   // particle list for the atmosphere
   partLT partsL;
   size_t notp = nocp;
   while (rhoLUT(rmed) > rhoco)
   {
      rmed = rmin;
      fType rold = -rmin;
      while (fabs((rold - rmed) / rmed) > shelleps)
      {
         lmed = pow(2. * migas / rhoLUT(rmed), 1. / 3.);
         rold = rmed;
         rmed = rmin + 0.5 * shellWidth * lmed;
      }
      rmax = rmin + shellWidth * lmed; 

      const fType Vshell = (4. / 3.) * M_PI * (pow(rmax, 3.) - pow(rmin, 3.));

      lattT Lattice(lmed, 1.01 * rmax, 1.01 * rmax, 1.01 * rmax);

      // count the particles
      size_t pc = 0;
      Lattice.first();
      while (!Lattice.isLast)
      {
         if ((Lattice.rCur < rmax) && (Lattice.rCur > rmin))
            pc++;
         Lattice.next();
      }

      // assign the particles
      const fType Vpart = Vshell / pc;
      const fType h     = 0.85 * lmed;

      Lattice.first();
      while (!Lattice.isLast)
      {
         if ((Lattice.rCur < rmax) && (Lattice.rCur > rmin))
         {
            partT np;
            np.pos[0] = Lattice.xCur;
            np.pos[1] = Lattice.yCur;
            np.pos[2] = Lattice.zCur;

            np.h = h;
            np.m = Vpart * rhoLUT(Lattice.rCur);

            np.u   = uLUT(Lattice.rCur);
            np.mat = 0;
            np.id = notp;

            partsL.push_back(np);
            pc++;
            notp++;
         }
         Lattice.next();
      }

      std::cout << "new shell r = " << rmin << " ... " << rmax << "  with " << pc << " parts \n";
      rmin = rmax;
   }


   // inserting atmosphere particles
   partLT::const_iterator plItr;
   for (partLT::const_iterator plItr = partsL.begin();
        plItr != partsL.end();
        plItr++)
      parts.insert(*plItr);

   const size_t nop = parts.getNop();

   treeT&      Tree(treeT::instance());
   const fType costppart = 1. / nop;

   Tree.setExtent(parts.getBox() * 1.1);
   for (size_t i = 0; i < nop; i++)
   {
      parts[i].cost = costppart;
      Tree.insertPart(parts[i]);
   }
   std::cout << "created tree\n";

   Tree.update(0.8, 1.2);
   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   std::cout << "Tree.update() -> " << noCZbottomLoc << " CZ cells\n";

   densSumT densWorker(&Tree);
#pragma omp parallel for firstprivate(densWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      densWorker(CZbottomLoc[i]);
   std::cout << "Tree.densWorker()\n";

   costT costWorker(&Tree);
#pragma omp parallel for firstprivate(costWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      costWorker(CZbottomLoc[i]);
   std::cout << "Tree.costWorker()\n";

   std::cout << "init EOS\n";
   eosT& EOS(eosT::instance());
   std::cout << "load table\n";
   EOS.aneos.loadTableU("aneos_tables.hdf5", 1);
   EOS.aneos.loadTableU("aneos_tables.hdf5", 2);
   EOS.aneos.loadTableU("aneos_tables.hdf5", 4);
   EOS.aneos.loadTableU("aneos_tables.hdf5", 5);

   const fType gasGamma = profFile.loadAttribute("gamma");
   std::cout << "set gamma:      " << gasGamma << "\n";
   parts.attributes["gamma"] = gasGamma;
   EOS.idealgas.setGamma(gasGamma);
   
   const fType molMass  = profFile.loadAttribute("molarmass");
   std::cout << "set molar mass: " << molMass  << "\n";
   parts.attributes["molarmass"] = molMass;
   EOS.idealgas.setMolarmass(molMass);

   std::cout << "EOS()\n";
   for (size_t i = 0; i < nop; i++)
      EOS(parts[i]);


   Tree.clear();
   std::cout << "Tree.clear()\n";

   parts.saveHDF5(odname);


#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}
