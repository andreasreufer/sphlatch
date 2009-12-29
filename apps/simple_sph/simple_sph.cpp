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
typedef sphlatch::vect3dT   vect3dT;
typedef sphlatch::box3dT    box3dT;

const fType finf = sphlatch::fTypeInf;

#define SPHLATCH_LOGGER
#include "logger.cpp"
typedef sphlatch::Logger   logT;

#include "bhtree.cpp"
typedef sphlatch::BHTree   treeT;

///
/// define the particle we are using
///
#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "io_particle.h"
#include "integrator_predcorr.cpp"

class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::energyPart,
   public sphlatch::IOPart
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   , public sphlatch::varHPart
#endif
#ifdef SPHLATCH_ANEOS
   , public sphlatch::ANEOSPart
#endif
{
public:

   sphlatch::PredictorCorrectorO2<vect3dT> posInt;
   sphlatch::PredictorCorrectorO1<fType>   energyInt;

   void bootstrap()
   {
      posInt.bootstrap(pos, vel, acc);
      energyInt.bootstrap(u, dudt);
   }

   void predict(const fType _dt)
   {
      posInt.predict(pos, vel, acc, _dt);
      energyInt.predict(u, dudt, _dt);
   }

   void correct(const fType _dt)
   {
      posInt.correct(pos, vel, acc, _dt);
      energyInt.correct(u, dudt, _dt);
   }

   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));

      vars.push_back(storeVar(u, "u"));
#ifdef SPHLATCH_ANEOS
      vars.push_back(storeVar(mat, "mat"));
#endif
#ifdef SPHLATCH_GRAVITY_EPSSMOOTHING
      vars.push_back(storeVar(eps, "eps"));
#endif
      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(acc, "acc"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));
      vars.push_back(storeVar(rho, "rho"));
      vars.push_back(storeVar(noneigh, "noneigh"));

      vars.push_back(storeVar(p, "p"));
      vars.push_back(storeVar(u, "u"));
      vars.push_back(storeVar(cs, "cs"));
#ifdef SPHLATCH_ANEOS
      vars.push_back(storeVar(mat, "mat"));
      vars.push_back(storeVar(T, "T"));
      vars.push_back(storeVar(phase, "phase"));
#endif
#ifdef SPHLATCH_GRAVITY_EPSSMOOTHING
      vars.push_back(storeVar(eps, "eps"));
#endif
      return(vars);
   }
};

typedef particle   partT;

class ghost :
   public sphlatch::treeGhost,
   public sphlatch::movingGhost
{ };

typedef ghost                                  ghstT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;

#ifdef SPHLATCH_GRAVITY
 #include "bhtree_worker_grav.cpp"
typedef sphlatch::fixThetaMAC                  macT;
typedef sphlatch::GravityWorker<macT, partT>   gravT;
#endif

#include "sph_algorithms.cpp"
#include "sph_kernels.cpp"
#include "bhtree_worker_sphsum.cpp"

typedef sphlatch::CubicSpline3D                  krnlT;

typedef sphlatch::densSum<partT, krnlT>          densT;
typedef sphlatch::SPHsumWorker<densT, partT>     densSumT;

typedef sphlatch::accPowSum<partT, krnlT>        accPowT;
typedef sphlatch::SPHsumWorker<accPowT, partT>   accPowSumT;

#ifdef SPHLATCH_ANEOS
 #include "eos_aneos.cpp"
typedef sphlatch::ANEOS<partT>                   eosT;
#else
 #include "eos_idealgas.cpp"
typedef sphlatch::IdealGas<partT>                eosT;
#endif

// particles are global
partSetT parts;

void derive()
{
   treeT& Tree(treeT::instance());
   logT&  Logger(logT::instance());
   
   const size_t nop       = parts.getNop();
   const fType  costppart = 1. / nop;

   Tree.setExtent(parts.getBox() * 1.1);

   for (size_t i = 0; i < nop; i++)
   {
      parts[i].cost = costppart;
      Tree.insertPart(parts[i]);
      parts[i].bootstrap();
   }
   Logger << "created tree";

   Tree.update(0.8, 1.2);
   Logger << "Tree.update()";

   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   //const size_t nop = parts.getNop();
   for (size_t i = 0; i < nop; i++)
   {
      parts[i].acc  = 0., 0., 0.;
      parts[i].dudt = 0.;
   }

#ifdef SPHLATCH_GRAVITY
   const fType G = parts.attributes["gravconst"];
   gravT       gravWorker(&Tree, G);
 #pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcGravity(CZbottomLoc[i]);
   Logger << "Tree.calcGravity()";
#endif

   densSumT densWorker(&Tree);
#pragma omp parallel for firstprivate(densWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      densWorker(CZbottomLoc[i]);
   Logger << "Tree.densWorker()";

   eosT& EOS(eosT::instance());
   for (size_t i = 0; i < nop; i++)
   {
      EOS(parts[i]);
   }
   Logger << "pressure";

   accPowSumT accPowWorker(&Tree);
#pragma omp parallel for firstprivate(accPowWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      accPowWorker(CZbottomLoc[i]);
   Logger << "Tree.accPowWorker()";

   Tree.clear();
   Logger << "Tree.clear()";
}

fType timestep()
{
   logT& Logger(logT::instance());

   const fType courantNumber = parts.attributes["courant"];
   const fType time          = parts.attributes["time"];

   fType dtA   = finf;
   fType dtCFL = finf;

#ifdef SPHLATCH_TIMEDEP_ENERGY
   fType dtU = finf;
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   fType dtH = finf;
#endif
   
   const size_t nop = parts.getNop();
   for (size_t i = 0; i < nop; i++)
   {
      const fType ai = sqrt(dot(parts[i].acc, parts[i].acc));
      if (ai > 0.)
      {
         const fType dtAi = sqrt(parts[i].h / ai);
         dtA = dtAi < dtA ? dtAi : dtA;
      }

      const fType dtCFLi = parts[i].h / parts[i].cs;
      dtCFL = dtCFLi < dtCFL ? dtCFLi : dtCFL;

#ifdef SPHLATCH_TIMEDEP_ENERGY
      const fType dudti = parts[i].dudt;
      if (dudti < 0.)
      {
         const fType dtUi = -parts[i].u / dudti;
         dtU = dtUi < dtU ? dtUi : dtU;
      }
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
      const fType dhdti = parts[i].dhdt;
      if (dhdti < 0.)
      {
         const fType dtHi = -parts[i].h / dhdti;
         dtH = dtHi < dtH ? dtHi : dtH;
      }
#endif
   }

   dtA   *= 0.5;
   dtCFL *= courantNumber;
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   dtH *= 0.15;
#endif

   //FIXME: globally minimize all dts
   fType dt = finf;

   dt = dtA < dt ? dtA : dt;
   dt = dtCFL < dt ? dtCFL : dt;

#ifdef SPHLATCH_TIMEDEP_ENERGY
   dt = dtU < dt ? dtU : dt;
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   dt = dtH < dt ? dtH : dt;
#endif
   Logger.stream << "dt: " << dt    << " "
                 << "dtA: " << dtA   << " "
                 << "dtCFL: " << dtCFL << " "
#ifdef SPHLATCH_TIMEDEP_ENERGY
                 << "dtU: " << dtU   << " "
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
                 << "dtH: " << dtH   << " "
#endif
                 << " ";
   Logger.flushStream();
   return(dt);
}

int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif

   omp_set_num_threads(4);

   if (argc != 4)
   {
      std::cerr <<
      "usage: simple_sph_XXXXXX <inputdump> <saveStepTime> <stopTime>\n";
      return(1);
   }

   std::string inFilename = argv[1];

   std::istringstream stepStr(argv[2]);
   fType stepTime;
   stepStr >> stepTime;

   std::istringstream stopStr(argv[3]);
   fType stopTime;
   stopStr >> stopTime;

   logT& Logger(logT::instance());

   ///
   /// log program compilation time
   ///
   Logger.stream << "executable compiled from " << __FILE__
                 << " on " << __DATE__
                 << " at " << __TIME__ << "\n\n"
                 << "    features: \n"
#ifdef SPHLATCH_GRAVITY
                 << "     gravity  \n"
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
                 << "     time dependent smoothing length\n"
#endif
#ifdef SPHLATCH_TIMEDEP_ENERGY
                 << "     time dependent specific energy\n"
#endif
#ifdef SPHLATCH_TILLOTSON
                 << "     Tillotson EOS\n"
#endif
#ifdef SPHLATCH_ANEOS
                 << "     ANEOS\n"
#endif
#ifdef SPHLATCH_INTEGRATERHO
                 << "     integrated density\n"
#endif
#ifdef SPHLATCH_FRICTION
                 << "     friction\n"
#endif
                 << "     basic SPH\n";
   Logger.flushStream();

   // load the particles
   parts.loadHDF5(inFilename);
   Logger << "loaded particles";

   treeT& Tree(treeT::instance());

   /*const size_t nop       = parts.getNop();
   const fType  costppart = 1. / nop;

   Tree.setExtent(parts.getBox() * 1.5);

   for (size_t i = 0; i < nop; i++)
   {
      parts[i].cost = costppart;
      Tree.insertPart(parts[i]);
      parts[i].bootstrap();
   }
   Logger << "created tree";*/

   fType& time( parts.attributes["time"] );

   // start the loop
   while (time < stopTime)
   {
      derive();
   
      const size_t nop = parts.getNop();
      const fType dt = timestep();

      for (size_t i = 0; i < nop; i++)
         parts[i].predict(dt);
      Logger.finishStep("predicted");

      derive();
      
      for (size_t i = 0; i < nop; i++)
         parts[i].correct(dt);
      Logger.finishStep("corrected");
      
      time += dt;
   }

   parts.doublePrecOut();
   parts.saveHDF5("out_both.h5part");
   Logger << "stored dump ";

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}
