#include <iostream>
#include <sstream>
#include <vector>

//#define SPHLATCH_SINGLEPREC

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#ifdef SPHLATCH_TIMEDEP_ENERGY
 #ifndef SPHLATCH_VELDIV
  #define SPHLATCH_VELDIV
 #endif
#endif

#ifdef SPHLATCH_TIMEDEP_SMOOTHING
 #ifndef SPHLATCH_VELDIV
  #define SPHLATCH_VELDIV
 #endif
 #ifndef SPHLATCH_NONEIGH
  #define SPHLATCH_NONEIGH
 #endif
#endif

#ifdef SPHLATCH_INTEGRATERHO
 #ifndef SPHLATCH_VELDIV
  #define SPHLATCH_VELDIV
 #endif
#endif

#ifdef SPHLATCH_ANEOS
 #ifndef SPHLATCH_NONEGPRESS
  #define SPHLATCH_NONEGPRESS
 #endif
#endif


#include "typedefs.h"
typedef sphlatch::fType     fType;
typedef sphlatch::cType     cType;
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
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   sphlatch::PredictorCorrectorO1<fType>   smolenInt;
#endif

   void bootstrap()
   {
      posInt.bootstrap(pos, vel, acc);
#ifdef SPHLATCH_TIMEDEP_ENERGY
      energyInt.bootstrap(u, dudt);
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
      smolenInt.bootstrap(h, dhdt);
#endif
   }

   void predict(const fType _dt)
   {
      posInt.predict(pos, vel, acc, _dt);
#ifdef SPHLATCH_TIMEDEP_ENERGY
      energyInt.predict(u, dudt, _dt);
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
      smolenInt.predict(h, dhdt, _dt);
#endif
   }

   void correct(const fType _dt)
   {
      posInt.correct(pos, vel, acc, _dt);
#ifdef SPHLATCH_TIMEDEP_ENERGY
      energyInt.correct(u, dudt, _dt);
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
      smolenInt.correct(h, dhdt, _dt);
#endif
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
      vars.push_back(storeVar(S, "S"));
      vars.push_back(storeVar(phase, "phase"));
#endif
#ifdef SPHLATCH_GRAVITY_EPSSMOOTHING
      vars.push_back(storeVar(eps, "eps"));
#endif
#ifdef SPHLATCH_VELDIV
      vars.push_back(storeVar(divv, "divv"));
#endif
#ifdef SPHLATCH_TIMEDEP_ENERGY
      vars.push_back(storeVar(dudt, "dudt"));
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
      vars.push_back(storeVar(dhdt, "dhdt"));
#endif

      vars.push_back(storeVar(cost, "cost"));
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

#include "bhtree_worker_cost.cpp"
typedef sphlatch::CostWorker<partT>              costT;

#ifdef SPHLATCH_ANEOS
 #include "eos_aneos.cpp"
typedef sphlatch::ANEOS<partT>                   eosT;
#else
 #include "eos_idealgas.cpp"
typedef sphlatch::IdealGas<partT>                eosT;
#endif

#ifdef SPHLATCH_KEEPENERGYPROFILE
#include "lookup_table1D.cpp"
typedef sphlatch::InterpolateLinear              intplT;
typedef sphlatch::LookupTable1D<intplT>          engLUT;

engLUT energyLUT("profile1D.hdf5", "r", "u");
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
   }
   Logger << "created tree";

   Tree.update(0.8, 1.2);

   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   Logger.stream << "Tree.update() -> " << noCZbottomLoc << " CZ cells";
   Logger.flushStream();

   //const size_t nop = parts.getNop();
   for (size_t i = 0; i < nop; i++)
   {
      parts[i].acc = 0., 0., 0.;
#ifdef SPHLATCH_TIMEDEP_ENERGY
      parts[i].dudt = 0.;
#endif
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

#ifdef SPHLATCH_TIMEDEP_ENERGY
   const fType uMin = parts.attributes["umin"];
   Logger.stream << "assure minimal spec. energy umin = " << uMin;
   Logger.flushStream();
#endif

   eosT& EOS(eosT::instance());
   for (size_t i = 0; i < nop; i++)
   {
#ifdef SPHLATCH_TIMEDEP_ENERGY
      if (parts[i].u < uMin)
         parts[i].u = uMin;
#endif
      EOS(parts[i]);
   }
   Logger << "pressure";

   accPowSumT accPowWorker(&Tree);
#pragma omp parallel for firstprivate(accPowWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      accPowWorker(CZbottomLoc[i]);
   Logger << "Tree.accPowWorker()";

#ifdef SPHLATCH_VELDIV
   for (size_t i = 0; i < nop; i++)
      parts[i].setDivvMax();
   Logger << "setDivvMax()";
#endif

#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   for (size_t i = 0; i < nop; i++)
      setDhDt(parts[i]);
#endif

#ifdef SPHLATCH_FRICTION
   const fType fricCoeff = 1. / parts.attributes["frictime"];
   for (size_t i = 0; i < nop; i++)
      parts[i].acc -= parts[i].vel * fricCoeff;
   Logger.stream << "friction (t_fric = " << 1. / fricCoeff << ")";
   Logger.flushStream();
#endif

#ifdef SPHLATCH_KEEPENERGYPROFILE
   vect3dT com;
   com = 0., 0., 0.;
   fType   totM = 0.;
   for (size_t i = 0; i < nop; i++)
   {
      com  += parts[i].pos * parts[i].m;
      totM += parts[i].m;
   }
   com /= totM;

   Logger.stream << "center of mass: ["
                 << com[0] << ","
                 << com[1] << ","
                 << com[2] << "]";
   Logger.flushStream();

   const fType thermFricCoeff = 1. / parts.attributes["frictime"];
   for (size_t i = 0; i < nop; i++)
   {
      const vect3dT rveci  = parts[i].pos - com;
      const fType   ri     = sqrt(dot(rveci, rveci));
      const fType   utheoi = energyLUT(ri);

      parts[i].dudt -= (parts[i].u - utheoi) * thermFricCoeff;
   }
#endif

   costT costWorker(&Tree);
#pragma omp parallel for firstprivate(costWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      costWorker(CZbottomLoc[i]);
   Logger << "Tree.costWorker()";

   Tree.clear();
   Logger << "Tree.clear()";
}

fType timestep(const fType _stepTime)
{
   logT& Logger(logT::instance());

   const fType courant = parts.attributes["courant"];
   const fType time    = parts.attributes["time"];

   fType dtSave = (floor((time / _stepTime) + 1.e-6) + 1.) * _stepTime - time;
   fType dtA    = finf;
   fType dtCFL  = finf;

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
   dtCFL *= courant;
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   dtH *= 0.15;
#endif

   //FIXME: globally minimize all dts
   fType dt = finf;

   dt = dtSave < dt ? dtSave : dt;
   dt = dtA < dt ? dtA : dt;
   dt = dtCFL < dt ? dtCFL : dt;

#ifdef SPHLATCH_TIMEDEP_ENERGY
   dt = dtU < dt ? dtU : dt;
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   dt = dtH < dt ? dtH : dt;
#endif
   Logger.stream << "dt: " << dt << " "
                 << "dtSave: " << dtSave << " "
                 << "dtA: " << dtA << " "
                 << "dtCFL: " << dtCFL << " "
#ifdef SPHLATCH_TIMEDEP_ENERGY
                 << "dtU: " << dtU << " "
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
                 << "dtH: " << dtH << " "
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

   if (not ((argc == 4) || (argc == 5)))
   {
      std::cerr <<
      "usage: simple_sph_XXXXXX <inputdump> <saveStepTime> <stopTime> (<numthreads>)\n";
      return(1);
   }


   std::string inFilename = argv[1];

   std::istringstream stepStr(argv[2]);
   fType stepTime;
   stepStr >> stepTime;

   std::istringstream stopStr(argv[3]);
   fType stopTime;
   stopStr >> stopTime;

   if (argc == 5)
   {
      std::istringstream threadStr(argv[4]);
      int numThreads;
      threadStr >> numThreads;
      omp_set_num_threads(numThreads);
   }

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

   Logger.stream << "working on " << omp_get_num_threads() << " threads";
   Logger.flushStream();

   // load the particles
   parts.loadHDF5(inFilename);
   Logger << "loaded particles";

   fType& time(parts.attributes["time"]);
   cType& step(parts.step);
   treeT& Tree(treeT::instance());

   const size_t nop = parts.getNop();

   parts[0].noneighOpt = 50;

   // first bootstrapping step
   derive();
   for (size_t i = 0; i < nop; i++)
      parts[i].bootstrap();

   Logger.finishStep("bootstrapped integrator");

   parts.doublePrecOut();
   parts.saveHDF5("bootstrap.h5part");

   // start the loop
   while (time < stopTime)
   {
      derive();

      const fType dt = timestep(stepTime);

      for (size_t i = 0; i < nop; i++)
         parts[i].predict(dt);
      Logger.finishStep("predicted");

      time += dt;
      derive();

      for (size_t i = 0; i < nop; i++)
         parts[i].correct(dt);
      step++;

      std::stringstream sstr;
      sstr << "corrected (t = " << time << ")";
      Logger.finishStep(sstr.str());

      const fType dumpTime = (floor((time / stepTime) + 1.e-9)) * stepTime;
      if (fabs(dumpTime - time) < 1.e-9)
      {
         std::stringstream dumpStr, stepStr, timeStr;

         dumpStr << "dump";

         stepStr << step;
         // pad step number to 7 digits
         for (size_t i = stepStr.str().size(); i < 7; i++)
            dumpStr << "0";
         dumpStr << stepStr.str();

         dumpStr << "_T";
         timeStr << std::setprecision(4) << std::scientific << time;
         // pad time string to 11 digits
         for (size_t i = timeStr.str().size(); i < 11; i++)
            dumpStr << "0";
         dumpStr << timeStr.str();
         dumpStr << ".h5part";

         parts.saveHDF5(dumpStr.str());

         Logger.stream << "write " << dumpStr.str();
         Logger.flushStream();
      }
   }

   Logger << "stored dump ";

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}
