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

#ifndef SPHLATCH_TIMEDEP_ENERGY
 #undef SPHLATCH_TRACK_UAV
#endif

#ifndef SPHLATCH_ANEOS
 #undef SPHLATCH_TRACK_TMAX
#endif

#ifdef SPHLATCH_ESCAPEES
 #ifndef SPHLATCH_FIND_CLUMPS
  #define SPHLATCH_FIND_CLUMPS
 #endif
#endif

#ifdef SPHLATCH_LRDISK
 #ifndef SPHLATCH_FIND_CLUMPS
  #define SPHLATCH_FIND_CLUMPS
 #endif
#endif


#include "typedefs.h"
typedef sphlatch::fType             fType;
typedef sphlatch::cType             cType;
typedef sphlatch::vect3dT           vect3dT;
typedef sphlatch::box3dT            box3dT;
typedef sphlatch::partsIndexListT   plistT;

const fType finf = sphlatch::fTypeInf;

#define SPHLATCH_LOGGER
#include "logger.cpp"
typedef sphlatch::Logger     logT;

#include "hdf5_io.cpp"
typedef sphlatch::HDF5File   H5FT;

#include "bhtree.cpp"
typedef sphlatch::BHTree     treeT;

///
/// define the particle we are using
///
#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "io_particle.h"
#include "integrator_predcorr.cpp"

#ifdef SPHLATCH_FIND_CLUMPS
 #include "clump_particle.h"
#endif

#ifdef SPHLATCH_LRDISK
 #include "friend_particle.h"
#endif

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
#ifdef SPHLATCH_FIND_CLUMPS
   , public sphlatch::clumpPart
#endif
#ifdef SPHLATCH_LRDISK
   , public sphlatch::friendPart
#endif
{
public:
   sphlatch::PredictorCorrectorO2<vect3dT> posInt;
   sphlatch::PredictorCorrectorO1<fType>   energyInt;
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   sphlatch::PredictorCorrectorO1<fType>   smolenInt;
#endif
#ifdef SPHLATCH_TRACK_UAV
   sphlatch::PredictorCorrectorO1<fType>   avengergyInt;
#endif
#ifdef SPHLATCH_INTEGRATERHO
   sphlatch::PredictorCorrectorO1<fType>   densInt;
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
#ifdef SPHLATCH_TRACK_UAV
      avengergyInt.bootstrap(uav, dudtav);
#endif
#ifdef SPHLATCH_INTEGRATERHO
      densInt.bootstrap(rho, drhodt);
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
#ifdef SPHLATCH_TRACK_UAV
      avengergyInt.predict(uav, dudtav, _dt);
#endif
#ifdef SPHLATCH_INTEGRATERHO
      densInt.predict(rho, drhodt, _dt);
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
#ifdef SPHLATCH_TRACK_UAV
      avengergyInt.correct(uav, dudtav, _dt);
#endif
#ifdef SPHLATCH_INTEGRATERHO
      densInt.correct(rho, drhodt, _dt);
#endif
   }

#ifdef SPHLATCH_TRACK_TMAX
   fType Tmax;
#endif
#ifdef SPHLATCH_TRACK_PMAX
   fType pmax;
#endif
#ifdef SPHLATCH_TRACK_UAV
   fType uav, dudtav;
#endif
#ifdef SPHLATCH_ESCAPEES
   fType rmvtime;
#endif
#ifdef SPHLATCH_TRACK_ACCP
   vect3dT accp;
#endif
#if defined SPHLATCH_FIND_CLUMPS && defined SPHLATCH_GRAVITY
   fType morig;
#endif

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
#ifdef SPHLATCH_TRACK_TMAX
      vars.push_back(storeVar(Tmax, "Tmax"));
#endif
#ifdef SPHLATCH_TRACK_PMAX
      vars.push_back(storeVar(pmax, "pmax"));
#endif
#ifdef SPHLATCH_TRACK_UAV
      vars.push_back(storeVar(uav, "uav"));
#endif
#ifdef SPHLATCH_TIMEDEP_ENTROPY
      vars.push_back(storeVar(S, "S"));
#endif
#ifdef SPHLATCH_INTEGRATERHO
      vars.push_back(storeVar(rho, "rho"));
#endif
      vars.push_back(storeVar(cost, "cost"));
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
      vars.push_back(storeVar(rhoL, "rhoL"));
      vars.push_back(storeVar(rhoH, "rhoH"));
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
#ifdef SPHLATCH_INTEGRATERHO
      vars.push_back(storeVar(drhodt, "drhodt"));
#endif
#ifdef SPHLATCH_GRAVITY
      vars.push_back(storeVar(pot, "pot"));
#endif
#ifdef SPHLATCH_FIND_CLUMPS
      vars.push_back(storeVar(clumpid, "clumpid"));
      vars.push_back(storeVar(orbit, "orbit"));
      vars.push_back(storeVar(ecc, "ecc"));
#endif
#ifdef SPHLATCH_TRACK_TMAX
      vars.push_back(storeVar(Tmax, "Tmax"));
#endif
#ifdef SPHLATCH_TRACK_PMAX
      vars.push_back(storeVar(pmax, "pmax"));
#endif
#ifdef SPHLATCH_TRACK_ACCP
      vars.push_back(storeVar(accp, "accp"));
#endif
#ifdef SPHLATCH_TRACK_UAV
      vars.push_back(storeVar(uav, "uav"));
      vars.push_back(storeVar(dudtav, "dudtav"));
#endif
#ifdef SPHLATCH_ESCAPEES
      vars.push_back(storeVar(rmvtime, "rmvtime"));
#endif
#ifdef SPHLATCH_LRDISK
      vars.push_back(storeVar(friendid, "friendid"));
#endif
#ifdef SPHLATCH_MISCIBLE
      vars.push_back(storeVar(delta, "delta"));
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

#ifndef SPHLATCH_INTEGRATERHO
typedef sphlatch::densSum<partT, krnlT>          densT;
typedef sphlatch::SPHsumWorker<densT, partT>     densSumT;
#endif
typedef sphlatch::accPowSum<partT, krnlT>        accPowT;
typedef sphlatch::SPHsumWorker<accPowT, partT>   accPowSumT;

#include "bhtree_worker_cost.cpp"
typedef sphlatch::CostWorker<partT>              costT;

#include "eos_super.cpp"
typedef sphlatch::SuperEOS<partT>                eosT;

#ifdef SPHLATCH_KEEPENERGYPROFILE
 #include "lookup_table1D.cpp"
typedef sphlatch::InterpolateLinear              intplT;
typedef sphlatch::LookupTable1D<intplT>          engLUT;

engLUT energyLUT("profile1D.hdf5", "r", "u");
#endif

#ifdef SPHLATCH_FIND_CLUMPS
 #include "clump_finder.cpp"
typedef sphlatch::Clumps<partT>   clumpsT;

 #ifdef SPHLATCH_ESCAPEES
partSetT escapees;
 #endif

 #ifdef SPHLATCH_LRDISK
  #include "disk_binner.cpp"
typedef sphlatch::DiskBinner<partT>   dbinT;
 #endif

#endif

#include "misc_physics.cpp"
using sphlatch::addToCOM;

// particles are global
partSetT parts;

#ifdef SPHLATCH_FIND_CLUMPS
clumpsT clumps;
#endif

void derive()
{
   treeT& Tree(treeT::instance());
   logT&  Logger(logT::instance());

   Logger << "start deriving step";

   const size_t nop = parts.getNop();

   Tree.setExtent(parts.getBox() * 1.1);

   // renormalize
   fType totCost = 0.;
   for (size_t i = 0; i < nop; i++)
      totCost += parts[i].cost;
   for (size_t i = 0; i < nop; i++)
      parts[i].cost /= totCost;

   // avoid extremes of relative cost
   const fType meanCost = 1. / nop;
   const fType minCost  = 0.01 * meanCost;
   const fType maxCost  = 20. * meanCost;

   for (size_t i = 0; i < nop; i++)
   {
      fType& costi(parts[i].cost);
      if (costi > maxCost)
         costi = maxCost;
      if (costi < minCost)
         costi = minCost;
   }

   for (size_t i = 0; i < nop; i++)
      Tree.insertPart(parts[i]);
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
      gravWorker.calcAcc(CZbottomLoc[i]);
   Logger << "Tree.calcAcc()";
#endif

#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   const fType hmin = parts.attributes["hmin"];
   Logger.stream << "assure minimal smth. length hmin = " << hmin;
   Logger.flushStream();
   for (size_t i = 0; i < nop; i++)
      if (parts[i].h < hmin)
         parts[i].h = hmin;
#endif

#ifndef SPHLATCH_INTEGRATERHO
   densSumT densWorker(&Tree);
 #pragma omp parallel for firstprivate(densWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      densWorker(CZbottomLoc[i]);
   Logger << "Tree.densWorker()";
#endif

   const fType pmin = parts.attributes["pmin"];
   Logger.stream << "assure minimal pressure     pmin = " << pmin;
   Logger.flushStream();

#ifdef SPHLATCH_TIMEDEP_ENERGY
   const fType uMin = parts.attributes["umin"];
   fType&      EthermAdded(parts.attributes["ethermadded"]);
   Logger.stream << "assure minimal spec. energy umin = " << uMin;
   Logger.flushStream();

   for (size_t i = 0; i < nop; i++)
   {
      fType& uCur(parts[i].u);
      if (uCur < uMin)
      {
         EthermAdded += (uMin - uCur) * parts[i].m;
         uCur         = uMin;
      }
   }
   Logger.stream << "               E_therm_added = " << EthermAdded;
   Logger.flushStream();
#endif

   eosT& EOS(eosT::instance());
   for (size_t i = 0; i < nop; i++)
   {
      EOS(parts[i]);

// FIXME: include ideal gas EOS for certain materials
#ifdef SPHLATCH_TRACK_TMAX
      fType& Tcur(parts[i].T);
      fType& Tmax(parts[i].Tmax);
      if (Tcur > Tmax)
         Tmax = Tcur;
#endif

      fType& pcur(parts[i].p);
      if (pcur < pmin)
         pcur = pmin;

#ifdef SPHLATCH_TRACK_PMAX
      fType& pmax(parts[i].pmax);
      if (pcur > pmax)
         pmax = pcur;
#endif
   }
   Logger << "pressure";

#ifdef SPHLATCH_TIMEDEP_ENERGY
#endif

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


#if defined SPHLATCH_KEEPENERGYPROFILE || defined SPHLATCH_SPINUP
   vect3dT com(0., 0., 0.);
   vect3dT vom(0., 0., 0.);
   fType   totM = 0.;
   for (size_t i = 0; i < nop; i++)
   {
      com  += parts[i].pos * parts[i].m;
      vom  += parts[i].vel * parts[i].m;
      totM += parts[i].m;
   }
   com /= totM;
   vom /= totM;

   Logger.stream << " center of mass: " << com;
   Logger.flushStream();
   Logger.stream << " vel of com.   : " << vom;
   Logger.flushStream();



 #ifdef SPHLATCH_KEEPENERGYPROFILE
   const fType thermFricCoeff = 1. / parts.attributes["frictime"];
 #endif
 #ifdef SPHLATCH_SPINUP
   const fType   spinupCoeff = 1. / parts.attributes["spinuptime"];
   const fType   omega       = (2 * M_PI) / parts.attributes["targetrotperiod"];
   const vect3dT omegavec(0., 0., omega);
   vect3dT       L(0., 0., 0.);
   fType         I = 0.;
 #endif


   for (size_t i = 0; i < nop; i++)
   {
      const vect3dT rveci = parts[i].pos - com;
      const vect3dT vveci = parts[i].vel - vom;
      const fType   ri    = sqrt(dot(rveci, rveci));

 #ifdef SPHLATCH_KEEPENERGYPROFILE
      const fType utheoi = energyLUT(ri);
      parts[i].dudt -= (parts[i].u - utheoi) * thermFricCoeff;
 #endif

 #ifdef SPHLATCH_SPINUP
      vect3dT vtarg = cross(omegavec, rveci);
      parts[i].acc += spinupCoeff * (vtarg - parts[i].vel);

      L += parts[i].m* cross(rveci, vveci);
      I += parts[i].m * ri * ri;
 #endif
   }

   vect3dT rotperact = (2. * M_PI * I) / (3600. * L);
 #ifdef SPHLATCH_KEEPENERGYPROFILE
   Logger << " enforced radial energy profile";
 #endif
 #ifdef SPHLATCH_SPINUP
   Logger.stream << " spin up to:    "
                 << ((2 * M_PI) / (omega / 3600.)) << " h";
   Logger.flushStream();
   Logger.stream << " current omega: " << rotperact << " h";
   Logger.flushStream();
 #endif
#endif

#ifdef SPHLATCH_ZONLY
   for (size_t i = 0; i < nop; i++)
      parts[i].acc[0] = parts[i].acc[1] = 0.;
#endif

   Tree.normalizeCost();
   costT costWorker(&Tree);
#pragma omp parallel for firstprivate(costWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      costWorker(CZbottomLoc[i]);
   Logger << "Tree.costWorker()";

   Tree.clear();
   Logger << "Tree.clear()";
}

fType timestep(const fType _stepTime, const fType _nextTime)
{
   logT& Logger(logT::instance());

   const fType courant = parts.attributes["courant"];
   const fType time    = parts.attributes["time"];

   fType dtSave = _nextTime - time;
   fType dtA    = finf;
   fType dtCFL  = finf;

#ifdef SPHLATCH_TIMEDEP_ENERGY
   fType dtU = finf;
   //const fType umin    = parts.attributes["umin"];
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   fType dtH = finf;
   //const fType hmin    = parts.attributes["hmin"];
#endif
#ifdef SPHLATCH_INTEGRATERHO
   fType dtRho = finf;
#endif

   const size_t nop = parts.getNop();
   for (size_t i = 0; i < nop; i++)
   {
      const fType ai = sqrt(dot(parts[i].acc, parts[i].acc));
      if (ai > 0.)
      {
         const fType dtAi = 0.25 * sqrt(parts[i].h / ai);
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
         //const fType dtHi = -( parts[i].h - hmin ) / dhdti;
         const fType dtHi = -(parts[i].h) / dhdti;
         dtH = dtHi < dtH ? dtHi : dtH;
      }
#endif
#ifdef SPHLATCH_INTEGRATERHO
      const fType drhodti = parts[i].drhodt;
      if (drhodti < 0.)
      {
         const fType dtRhoi = -(parts[i].rho) / drhodti;
         dtRho = dtRhoi < dtRho ? dtRhoi : dtRho;
      }
#endif
   }

   dtA   *= 0.5;
   dtCFL *= courant;
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   dtH *= 0.15;
#endif
#ifdef SPHLATCH_INTEGRATERHO
   dtRho *= 0.15;
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
#ifdef SPHLATCH_INTEGRATERHO
   dt = dtRho < dt ? dtRho : dt;
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
#ifdef SPHLATCH_INTEGRATERHO
                 << "dtRho: " << dtRho << " "
#endif
                 << " ";
   Logger.flushStream();
   return(dt);
}

void save(std::string _dumpPrefix)
{
   std::stringstream dumpStr, stepStr, timeStr;

   logT& Logger(logT::instance());

   const fType time = parts.attributes["time"];
   const cType step = parts.step;

   const size_t nop = parts.getNop();

   // create the dump filename
   dumpStr << _dumpPrefix;

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

   treeT& Tree(treeT::instance());
   Tree.setExtent(parts.getBox() * 1.1);

   for (size_t i = 0; i < nop; i++)
      Tree.insertPart(parts[i]);

   Logger << "created tree";

   Tree.update(0.8, 1.2);

   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

#ifdef SPHLATCH_GRAVITY
   const fType G = parts.attributes["gravconst"];
   gravT       gravWorker(&Tree, G);
 #pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcPot(CZbottomLoc[i]);
   Logger << "Tree.calcPot()";


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

  #ifdef SPHLATCH_GRAVITY
   const fType virialfact = parts.attributes["noclumpsvirialfactor"];
   parts.attributes["virialfactor"] = Ekin / (-0.5 * Epot);
   if (Ekin > -virialfact * 0.5 * Epot)
   {
      Logger << "too much kinetic energy, clumps search inhibited";
      clumps.noClumps(parts);
   }
   else
  #endif
   clumps.getClumps(parts, cMinMass);

   const size_t noc = clumps.getNop();
   Logger.stream << "found " << noc - 1 << " clump(s) with m > "
                 << cMinMass;
   Logger.flushStream();

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
         Logger.stream << "got orbits for clump " << i;
         Logger.flushStream();

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
         Logger << "    Tree.redoMultipoles()";

   #pragma omp parallel for firstprivate(gravWorker)
         for (int j = 0; j < noCZbottomLoc; j++)
            gravWorker.calcPot(CZbottomLoc[j]);
         Logger << "    Tree.calcPot()";

         fType EpotCC = 0.;
         for (size_t j = 0; j < nop; j++)
            if (parts[j].clumpid == cid)
               EpotCC += 0.5 * parts[j].pot * parts[j].m;

         clumps[i].Epot = EpotCC;

         Logger.stream << "    Epot = " << EpotCC
                       << ", Erot = " << clumps[i].Erot
                       << ", Ekin = " << clumps[i].Ekin;
         Logger.flushStream();
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
   Logger << "Tree.clear()";

#ifdef SPHLATCH_LRDISK
   const fType rmin = parts.attributes["diskrmin"];
   const fType rmax = parts.attributes["diskrmax"];

   dbinT diskBinner(parts);
   diskBinner.saveBins(1., rmin, rmax, 100, "disk.hdf5");
   Logger << "binned disk and stored to disk.hdf5";

 #ifdef SPHLATCH_LRDISKFOF
   const fType  hmultfof  = parts.attributes["hmultfof"];
   const fType  rhominfof = parts.attributes["rhomindiskfof"];
   const size_t nobins    = 100;

   diskBinner.findFOF(1., rhominfof, hmultfof, cMinMass);
   Logger.stream << "FOF clumps in disk with rho > " << rhominfof;
   Logger.flushStream();
 #endif
#endif

   parts.doublePrecOut();
   std::string pdumpFilename = dumpStr.str() + ".h5part";
   parts.saveHDF5(pdumpFilename);

   Logger.stream << "wrote " << pdumpFilename;
   Logger.flushStream();

#ifdef SPHLATCH_ESCAPEES
   if (clumps.getNop() > 1)
   {
      const fType rmaxubd = parts.attributes["rmaxunbound"];
      const fType rmaxbd  = parts.attributes["rmaxbound"];

      const vect3dT ccom = clumps[1].pos;
      plistT        esclst;
      for (size_t i = 0; i < nop; i++)
      {
         const vect3dT rvec = ccom - parts[i].pos;
         const fType   r    = sqrt(dot(rvec, rvec));

         if (((parts[i].clumpid == sphlatch::CLUMPNONE) && (r > rmaxubd)) or
             ((parts[i].clumpid > sphlatch::CLUMPNONE) && (r > rmaxbd)))
         {
            parts[i].rmvtime = time;
            esclst.push_back(i);
         }
      }

      plistT::const_reverse_iterator escItr;
      for (escItr = esclst.rbegin(); escItr != esclst.rend(); escItr++)
         escapees.insert(parts.pop(*escItr));

      escapees.step = parts.step;
      escapees.saveHDF5(dumpStr.str() + "_esc.h5part");

      Logger.stream << "removed " << esclst.size() << " escapees (" <<
      escapees.getNop() << " total, " << parts.getNop() << " particles left)";
      Logger.flushStream();
   }
#endif
}

int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif

   if (not ((argc == 5) || (argc == 6)))
   {
      std::cerr <<
      "usage: simple_sph_XXXXXX <inputdump> <saveStepTime> <stopTime> <dumpPrefix> (<numthreads>)\n";
      return(1);
   }


   std::string inFilename = argv[1];

   std::istringstream stepStr(argv[2]);
   fType stepTime;
   stepStr >> stepTime;

   std::istringstream stopStr(argv[3]);
   fType stopTime;
   stopStr >> stopTime;

   std::string dumpPrefix = argv[4];

   if (argc == 6)
   {
      std::istringstream threadStr(argv[5]);
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
#ifdef SPHLATCH_MISCIBLE
                 << "     miscible SPH\n"
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
                 << "     ANEOS EOS\n"
#endif
#ifdef SPHLATCH_MANEOS
                 << "     MANEOS EOS\n"
#endif
#ifdef SPHLATCH_NONEGPRESS
                 << "     no negative pressure\n"
#endif
#ifdef SPHLATCH_INTEGRATERHO
                 << "     integrated density\n"
#endif
#ifdef SPHLATCH_FRICTION
                 << "     friction\n"
#endif
#ifdef SPHLATCH_TRACK_TMAX
                 << "     track Tmax\n"
#endif
#ifdef SPHLATCH_TRACK_PMAX
                 << "     track Pmax\n"
#endif
                 << "     ideal gas EOS\n"
                 << "     basic SPH\n";
   Logger.flushStream();

   Logger.stream << "working on " << omp_get_num_threads() << " thread(s)";
   Logger.flushStream();

   // load the particles
   parts.loadHDF5(inFilename);
   const size_t nop = parts.getNop();


   Logger.stream << "loaded " << nop << " particles";
   Logger.flushStream();

   if (nop == 0)
   {
      Logger << "something's fishy, exiting ...";
#ifdef SPHLATCH_MPI
      MPI::Finalize();
#endif
      return(1);
   }

#ifdef SPHLATCH_ESCAPEES
   // load escapees, if a escapee file exists
   std::string inEscFilename = inFilename.substr(0,
                                                 inFilename.find(".h5part")) +
                               "_esc.h5part";
   struct stat statBuff;
   if (stat(inEscFilename.c_str(), &statBuff) == 0)
   {
      escapees.loadHDF5(inEscFilename);
      Logger.stream << "loaded " << escapees.getNop() << " escapees";
      Logger.flushStream();
   }
   else
   {
      Logger.stream << inEscFilename << " not found, continuing";
      Logger.flushStream();
   }
#endif

   // normalize cost
   fType totCost = 0.;
   for (size_t i = 0; i < nop; i++)
      totCost += parts[i].cost;

   if (totCost > 0.)
      for (size_t i = 0; i < nop; i++)
         parts[i].cost /= totCost;
   else
   {
      Logger << "re-normalize cost";
      const fType nopinv = 1. / static_cast<fType>(nop);
      for (size_t i = 0; i < nop; i++)
         parts[i].cost = nopinv;
   }

   if (parts.attributes.count("courant") == 0)
      parts.attributes["courant"] = 0.3;
#ifdef SPHLATCH_ESCAPEES
   if (parts.attributes.count("rmaxunbound") == 0)
      parts.attributes["rmaxunbound"] = 1.e10;

   if (parts.attributes.count("rmaxbound") == 0)
      parts.attributes["rmaxbound"] = 2.e10;

   if (parts.attributes.count("rhominclump") == 0)
      parts.attributes["rhominclump"] = 0.75 * 2.65;
#endif

#ifdef SPHLATCH_LRDISK
   if (parts.attributes.count("diskrmin") == 0)
      parts.attributes["diskrmin"] = 0.;
   if (parts.attributes.count("diskrmax") == 0)
      parts.attributes["diskrmax"] = 2.e9;
   if (parts.attributes.count("rhomindiskfof") == 0)
      parts.attributes["rhomindiskfof"] = 0.1;
   if (parts.attributes.count("hmultfof") == 0)
      parts.attributes["hmultfof"] = 1.5;
#endif
   if (parts.attributes.count("gamma") == 0)
      parts.attributes["gamma"] = 1.4;
   if (parts.attributes.count("molarmass") == 0)
      parts.attributes["molarmass"] = 1.007977; // atomic hydrogen

   if (parts.attributes.count("pmin") == 0)
      parts.attributes["pmin"] = 0.;

#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   if (parts.attributes.count("hmin") == 0)
      parts.attributes["hmin"] = 0.;
#endif

#ifdef SPHLATCH_FIND_CLUMPS
 #ifdef SPHLATCH_GRAVITY
   if (parts.attributes.count("noclumpsvirialfactor") == 0)
      parts.attributes["noclumpsvirialfactor"] = 3.;
 #endif
#endif

#ifdef SPHLATCH_XSPH
   if (parts.attributes.count("xsphfactor") == 0)
      parts.attributes["xsphfactor"] = 0.5;

   parts[0].xsphfactor = parts.attributes["xsphfactor"];
#endif

   fType& time(parts.attributes["time"]);
   cType& step(parts.step);
   //treeT& Tree(treeT::instance());

   parts[0].noneighOpt = 50;

   eosT& EOS(eosT::instance());
#ifdef SPHLATCH_ANEOS_TABLE
 #ifdef SPHLATCH_TIMEDEP_ENERGY
   EOS.aneos.loadTableU("aneos_tables.hdf5", 1);
   EOS.aneos.loadTableU("aneos_tables.hdf5", 2);
   EOS.aneos.loadTableU("aneos_tables.hdf5", 4);
   EOS.aneos.loadTableU("aneos_tables.hdf5", 5);
 #endif
 #ifdef SPHLATCH_TIMEDEP_ENTROPY
   EOS.aneos.loadTableS("aneos_tables.hdf5", 1);
   EOS.aneos.loadTableS("aneos_tables.hdf5", 2);
   EOS.aneos.loadTableS("aneos_tables.hdf5", 4);
   EOS.aneos.loadTableS("aneos_tables.hdf5", 5);
 #endif
#endif

   EOS.idealgas.setGamma(parts.attributes["gamma"]);
   EOS.idealgas.setMolarmass(parts.attributes["molarmass"]);

   // first bootstrapping step
   derive();
   for (size_t i = 0; i < nop; i++)
      parts[i].bootstrap();

   Logger.finishStep("bootstrapped integrator");

   fType nextTime = (floor(time / stepTime) + 1.) * stepTime;
   // start the loop

   while (time < stopTime)
   {
      derive();

      const fType dt = timestep(stepTime, nextTime);

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

      if (fabs(nextTime - time) < 1.e-9)
      {
         save(dumpPrefix);
         nextTime += stepTime;
      }
   }

   Logger << "simulation stopped";

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}
