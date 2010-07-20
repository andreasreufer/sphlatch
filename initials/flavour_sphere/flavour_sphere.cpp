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
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));

      vars.push_back(storeVar(u, "u"));
      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(mat, "mat"));
      vars.push_back(storeVar(u, "u"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(noneigh, "noneigh"));
      vars.push_back(storeVar(cost, "cost"));
      return(vars);
   }
};

typedef particle                       partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>   partSetT;

#include "sph_algorithms.cpp"
#include "sph_kernels.cpp"
#include "bhtree_worker_sphsum.cpp"

typedef sphlatch::CubicSpline3D                krnlT;

typedef sphlatch::densSum<partT, krnlT>        densT;
typedef sphlatch::SPHsumWorker<densT, partT>   densSumT;

#include "bhtree_worker_cost.cpp"
typedef sphlatch::CostWorker<partT>            costT;


#include "lookup_table1D.cpp"
typedef sphlatch::InterpolateLinear            intplT;
typedef sphlatch::LookupTable1D<intplT>        LUT1DlinT;

#include "hdf5_io.cpp"
typedef sphlatch::HDF5File                     H5FT;


// particles are global
partSetT parts;


int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif

   if (not (argc == 3))
   {
      std::cerr <<
      "usage: flavour_sphere <dump> <profile>\n";
      return(1);
   }


   std::string dname = argv[1];
   std::string pname = argv[2];

   // load the particles
   parts.loadHDF5(dname);
   const size_t nop = parts.getNop();

   std::cout << "loaded " << nop << "\n";

   LUT1DlinT rhoLUT(pname, "r", "rho");
   LUT1DlinT engLUT(pname, "r", "u");
   LUT1DlinT matLUT(pname, "r", "mat");

   std::cout << "1D profile " << pname << " loaded" << "\n";

   const fType rScale = rhoLUT.getXmax();
   const fType vScale = rScale * rScale * rScale;
   std::cout << " rScale: " << rScale << "\n";

   // find the center of mass
   vect3dT com;
   com = 0., 0., 0.;
   fType totM = 0.;
   for (size_t i = 0; i < nop; i++)
   {
      parts[i].pos *= rScale;
      parts[i].h   *= rScale;
      parts[i].m   *= vScale;

      com  += parts[i].pos * parts[i].m;
      totM += parts[i].m;
   }
   com /= totM;

   std::cout << "center of mass: ["
             << com[0] << ","
             << com[1] << ","
             << com[2] << "]\n";

   // scale positions, specific volume and smoothing length
   totM = 0.;
   for (size_t i = 0; i < nop; i++)
   {
      const vect3dT rveci = parts[i].pos - com;
      const fType   ri    = sqrt(dot(rveci, rveci));

      const fType rhoi = rhoLUT(ri);
      parts[i].m  *= rhoi;
      parts[i].rho = rhoi;

      parts[i].u   = engLUT(ri);
      parts[i].mat = lrint(matLUT(ri));

      totM += parts[i].m;
   }
   std::cout << " total mass: " << totM << "\n";

   treeT& Tree(treeT::instance());

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

   Tree.clear();
   std::cout << "Tree.clear()\n";

   // store some important attributes
   H5FT profile(pname);
   parts.attributes["G"]    = profile.loadAttribute("gravconst");
   parts.attributes["umin"] = profile.loadAttribute("umin");
   
   parts.attributes["courant"] = 0.3;
   parts.attributes["time"] = 0.;

   parts.doublePrecOut();
   parts.saveHDF5(dname);
   std::cout << "stored dump\n";

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}
