// some defs

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

// enable parallel version
#define SPHLATCH_PARALLEL

// enable intensive logging for toptree global summation
//#define SPHLATCH_TREE_LOGSUMUPMP

//#define SPHLATCH_RANKSPACESERIALIZE

// enable checking of bounds for the neighbour lists
//#define SPHLATCH_CHECKNONEIGHBOURS

// enable selfgravity?
//#define SPHLATCH_GRAVITY

// time dependent energy?
#define SPHLATCH_TIMEDEP_ENERGY

// time dependent smoothing length?
#define SPHLATCH_TIMEDEP_SMOOTHING

// do we need the veolicty divergence?
#ifdef SPHLATCH_TIMEDEP_ENERGY
#ifndef SPHLATCH_VELDIV
#define SPHLATCH_VELDIV
#endif
#endif

#ifdef SPHLATCH_TIMEDEP_SMOOTHING
#ifndef SPHLATCH_VELDIV
#define SPHLATCH_VELDIV
#endif
#endif

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>

#include <boost/assign/std/vector.hpp>

#include <boost/mpl/vector_c.hpp>

namespace po = boost::program_options;

#include "typedefs.h"
typedef sphlatch::fType fType;
typedef sphlatch::valueRefType valueRefType;
typedef sphlatch::valvectType valvectType;
typedef sphlatch::valvectRefType valvectRefType;

typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::matrixRefType matrixRefType;
typedef sphlatch::partsIndexVectType partsIndexVectType;
typedef sphlatch::quantsType quantsType;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager comm_type;

#include "costzone.h"
typedef sphlatch::CostZone costzone_type;

#include "log_manager.h"
typedef sphlatch::LogManager log_type;

#include "integrator_verlet.h"
#include "integrator_predcorr.h"

#include <boost/progress.hpp>
#include <vector>

#include "kernel_cubicspline.h"

/// tree stuff
#include "bhtree.h"

/// neighbour search
#include "rankspace.h"

using namespace sphlatch::vectindices;
using namespace boost::assign;

fType timeStep()
{
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  log_type& Logger(log_type::instance());
  io_type& IOManager(io_type::instance());

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  valvectRefType p(PartManager.p);
  valvectRefType u(PartManager.u);
  valvectRefType rho(PartManager.rho);
  valvectRefType dudt(PartManager.dudt);
  valvectRefType dhdt(PartManager.dhdt);
  valvectRefType divv(PartManager.divv);

#ifdef SPHLATCH_GRAVITY
  valvectRefType eps(PartManager.eps);
#endif

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  valueRefType time(PartManager.attributes["time"]);
  int& step(PartManager.step);
  //int& substep(PartManager.substep);
  const size_t noParts = PartManager.getNoLocalParts();
  //const size_t myDomain = CommManager.getMyDomain();

  ///
  /// timestep criterion for acceleration
  /// ( see Wetzstein et. al 2008 )
  ///
  fType dtA = std::numeric_limits<fType>::max();
  for (size_t i = 0; i < noParts; i++)
    {
      const fType ai = sqrt(acc(i, X) * acc(i, X) +
                                acc(i, Y) * acc(i, Y) +
                                acc(i, Z) * acc(i, Z));

      if (ai > 0.)
        {
          const fType dtAi = sqrt(h(i) / ai);
          dtA = dtAi < dtA ? dtAi : dtA;
        }
    }
  dtA *= 0.5;
  CommManager.min(dtA);

#ifdef SPHLATCH_TIMEDEP_ENERGY
  ///
  /// limit oooling speed in integration time
  /// energy integration
  ///
  fType dtU = std::numeric_limits<fType>::max();
  for (size_t i = 0; i < noParts; i++)
    {
      if (dudt(i) < 0.)
        {
          const fType dtUi = -u(i) / dudt(i);
          dtU = dtUi < dtU ? dtUi : dtU;
        }
    }
  CommManager.min(dtU);
#endif

#ifdef SPHLATCH_TIMEDEP_SMOOTHING
  ///
  /// timestep criterion for smoothing length integration
  /// ( see Wetzstein et. al 2008 )
  ///
  fType dtH = std::numeric_limits<fType>::max();
  for (size_t i = 0; i < noParts; i++)
    {
      const fType absdtHi = fabs(h(i) / dhdt(i));

      if (absdtHi > 0.)
        {
          dtH = absdtHi < dtH ? absdtHi : dtH;
        }
    }
  dtH *= 0.15;
  CommManager.min(dtH);
#endif

  ///
  /// CFL condition
  ///
  fType dtCFL = std::numeric_limits<fType>::max();
  const fType gamma = PartManager.attributes["gamma"];
  const fType courantNumber = PartManager.attributes["courant"];
  for (size_t i = 0; i < noParts; i++)
    {
      const fType ci = sqrt(p(i) * gamma / rho(i));
      const fType dtCFLi = h(i) / ci;

      if (ci > 0.)
        {
          dtCFL = dtCFLi < dtCFL ? dtCFLi : dtCFL;
        }
    }
  dtCFL *= courantNumber;
  CommManager.min(dtCFL);

  ///
  /// distance to next saveItrvl
  ///
  const fType saveItrvl = 10.;
  std::string fileName = "saveDump000000.h5part";
  const fType nextSaveTime = (floor((time / saveItrvl) + 1.e-6)
                                  + 1.) * saveItrvl;
  fType dtSave = nextSaveTime - time;

  ///
  /// determine global minimum.
  /// by parallelly minimizing the timesteps, we
  /// can estimate which ones are dominant
  ///
  fType dtGlob = std::numeric_limits<fType>::max();

  dtGlob = dtA < dtGlob ? dtA : dtGlob;
  dtGlob = dtCFL < dtGlob ? dtCFL : dtGlob;
#ifdef SPHLATCH_TIMEDEP_ENERGY
  dtGlob = dtU < dtGlob ? dtU : dtGlob;
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
  dtGlob = dtH < dtGlob ? dtH : dtGlob;
#endif
  dtGlob = dtSave < dtGlob ? dtSave : dtGlob;

  Logger.stream << "dtA: " << dtA
#ifdef SPHLATCH_TIMEDEP_ENERGY
                << " dtU: " << dtU
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
                << " dtH: " << dtH
#endif
                << " dtCFL: " << dtCFL
                << " dtSave: " << dtSave
                << "   dtGlob: " << dtGlob;
  Logger.flushStream();

  ///
  /// define the quantities to save in a dump
  ///
  quantsType saveQuants;
  saveQuants.vects += &pos, &vel, &acc;
  saveQuants.scalars += &m, &rho, &u, &p, &h;
  saveQuants.ints += &id, &noneigh;

#ifdef SPHLATCH_TIMEDEP_ENERGY
  saveQuants.scalars += &dudt;
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
  saveQuants.scalars += &dhdt;
  saveQuants.scalars += &divv;
#endif
#ifdef SPHLATCH_GRAVITY
  saveQuants.scalars += &eps;
#endif

  const fType curSaveTime = (floor((time / saveItrvl) + 1.e-9))
                                * saveItrvl;
  if (fabs(curSaveTime - time) < 1.e-9)
    {
      Logger << "save dump";

      std::string stepString = boost::lexical_cast<std::string>(step);

      fileName.replace(fileName.size() - 7 - stepString.size(),
                       stepString.size(), stepString);

      IOManager.saveDump(fileName, saveQuants);
    }

  return dtGlob;
}

void derivate()
{
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  //io_type& IOManager(io_type::instance());
  costzone_type& CostZone(costzone_type::instance());
  log_type& Logger(log_type::instance());

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  valvectRefType p(PartManager.p);
  valvectRefType u(PartManager.u);
  valvectRefType rho(PartManager.rho);
  valvectRefType dudt(PartManager.dudt);
  valvectRefType dhdt(PartManager.dhdt);
  valvectRefType divv(PartManager.divv);

#ifdef SPHLATCH_GRAVITY
  valvectRefType eps(PartManager.eps);
#endif

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  ///
  /// exchange particle data
  ///
  CostZone.createDomainPartsIndex();
  Logger << "new domain decomposition";
  CommManager.exchange(CostZone.domainPartsIndex,
                       CostZone.getNoGhosts());

  ///
  /// prepare ghost sends
  ///
  CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());
  Logger.stream << "distributed particles: "
                << PartManager.getNoLocalParts() << " parts. & "
                << PartManager.getNoGhostParts() << " ghosts";
  Logger.flushStream();

  const size_t noParts = PartManager.getNoLocalParts();
  const size_t noTotParts = noParts + PartManager.getNoGhostParts();
  int& step(PartManager.step);
  int& substep(PartManager.substep);
  //const size_t myDomain = CommManager.getMyDomain();

  /// send ghosts to other domains
  CommManager.sendGhosts(pos);
  CommManager.sendGhosts(vel);
  CommManager.sendGhosts(id);
  CommManager.sendGhosts(m);
  CommManager.sendGhosts(h);
#ifdef SPHLATCH_GRAVITY
  CommManager.sendGhosts(eps); // << eps is not yet used for interacting partners!
#endif
  Logger << " sent to ghosts: pos, vel, id, m, h, eps";

  ///
  /// zero the derivatives
  ///
  for (size_t i = 0; i < noParts; i++)
    {
      acc(i, X) = 0.;
      acc(i, Y) = 0.;
      acc(i, Z) = 0.;
#ifdef SPHLATCH_TIMEDEP_ENERGY
      dudt(i) = 0.;
#endif

      ///
      /// shocktube boundary conditions
      ///
      vel(i, X) = 0.;
      vel(i, Y) = 0.;

      if (pos(i, Z) < 5. || pos(i, Z) > 295.)
        vel(i, Z) = 0.;
    }


#ifdef SPHLATCH_GRAVITY
  const fType gravTheta = PartManager.attributes["gravtheta"];
  const fType gravConst = PartManager.attributes["gravconst"];

  sphlatch::BHtree<sphlatch::Quadrupoles> Tree(gravTheta,
                                               gravConst,
                                               CostZone.getDepth(),
                                               CostZone.getCenter(),
                                               CostZone.getSidelength()
                                               );

  ///
  /// fill up tree, determine ordering, calculate multipoles
  ///
  for (size_t k = 0; k < noTotParts; k++)
    {
      Tree.insertParticle(k);
    }

  Tree.detParticleOrder();
  Tree.calcMultipoles();
  Logger << "Tree ready";

  for (size_t k = 0; k < noParts; k++)
    {
      const size_t i = Tree.particleOrder[k];
      Tree.calcGravity(i);
    }
  Logger << "gravity calculated";
#endif

  ///
  /// define kernel and neighbour search algorithm
  ///
  sphlatch::CubicSpline3D Kernel;
  sphlatch::Rankspace RSSearch;

  ///
  /// the size of the neighbour list doesn't really
  /// affect the performance of the neighbour search,
  /// so it can be choosen quite large
  ///
  RSSearch.prepare();
  RSSearch.neighbourList.resize(1024);
  RSSearch.neighDistList.resize(1024);
  Logger << "Rankspace prepared";

  ///
  /// 1st SPH sum: density
  /// I need    : pos, h, m
  /// I provide : rho
  ///
  for (size_t k = 0; k < noParts; k++)
    {
      ///
      /// find neighbours
      ///
#ifdef SPHLATCH_GRAVITY
      const size_t i = Tree.particleOrder[k];
#else
      const size_t i = k;
#endif
      const fType hi = h(i);
      const fType srchRad = 2. * hi;
      RSSearch.findNeighbours(i, srchRad);

      const size_t noNeighs = RSSearch.neighbourList[0];

      ///
      /// store the number of neighbours
      ///
      noneigh(i) = noNeighs;

      static fType rhoi;
      rhoi = 0.;

      ///
      /// sum over neighbours
      ///
      for (size_t curNeigh = 1; curNeigh <= noNeighs; curNeigh++)
        {
          const fType r = RSSearch.neighDistList[curNeigh];
          const size_t j = RSSearch.neighbourList[curNeigh];

          const fType hij = 0.5 * (hi + h(j));

          rhoi += m(j) * Kernel.value(r, hij);
        }
      rho(i) = rhoi;
    }
  Logger << "SPH sum: rho";
  CommManager.sendGhosts(rho);
  Logger << " sent to ghosts: rho";

#ifdef SPHLATCH_TIMEDEP_ENERGY
  ///
  /// lower temperature bound
  ///

  const fType uMin = PartManager.attributes["umin"];
  for (size_t i = 0; i < noParts; i++)
    {
      if (u(i) < uMin)
        {
          u(i) = uMin;
        }
    }
  Logger.stream << "assure minimal temperature (umin = " << uMin << ")";
  Logger.flushStream();
#endif
  CommManager.sendGhosts(u);
  Logger << " sent to ghosts: u";

  ///
  /// pressure
  /// I need    : u, rho
  /// I provide : p
  ///
  const fType gamma = PartManager.attributes["gamma"];
  p = (gamma - 1) * (boost::numeric::ublas::element_prod(u, rho));
  Logger << "pressure";

  ///
  /// 2st SPH sum: acceleration, specific power & velocity divergence
  /// I need    : pos, vel, h, m, rho, u
  /// I provide : acc, dudt, divv
  ///
  const fType alpha = 1;
  const fType beta = 2;

  fType curAccX = 0., curAccY = 0., curAccZ = 0.;
#ifdef SPHLATCH_TIMEDEP_ENERGY
  fType curPow = 0.;
#endif
#ifdef SPHLATCH_VELDIV
  fType curVelDiv = 0.;
  fType divvMax = std::numeric_limits<fType>::min();
#endif
  for (size_t k = 0; k < noParts; k++)
    {
#ifdef SPHLATCH_GRAVITY
      const size_t i = Tree.particleOrder[k];
#else
      const size_t i = k;
#endif
      const fType hi = h(i);
      const fType rhoi = rho(i);
      const fType pi = p(i);

      const fType piOrhoirhoi = pi / (rhoi * rhoi);

      /// find the neighbours
      const fType srchRad = 2. * hi;
      RSSearch.findNeighbours(i, srchRad);

      const size_t noNeighs = RSSearch.neighbourList[0];

      curAccX = 0.;
      curAccY = 0.;
      curAccZ = 0.;
#ifdef SPHLATCH_TIMEDEP_ENERGY
      curPow = 0.;
#endif
#ifdef SPHLATCH_VELDIV
      curVelDiv = 0.;
#endif
      const fType viX = vel(i, X);
      const fType viY = vel(i, Y);
      const fType viZ = vel(i, Z);

      const fType riX = pos(i, X);
      const fType riY = pos(i, Y);
      const fType riZ = pos(i, Z);

      const fType ci = sqrt(gamma * pi / rhoi);

      ///
      /// sum over the neighbours
      ///
      for (size_t curNeigh = 1; curNeigh <= noNeighs; curNeigh++)
        {
          const fType rij = RSSearch.neighDistList[curNeigh];
          const size_t j = RSSearch.neighbourList[curNeigh];

          const fType rhoj = rho(j);
          const fType pj = p(j);

          const fType hij = 0.5 * (hi + h(j));

          const fType rijX = riX - pos(j, X);
          const fType rijY = riY - pos(j, Y);
          const fType rijZ = riZ - pos(j, Z);

          const fType vijX = viX - vel(j, X);
          const fType vijY = viY - vel(j, Y);
          const fType vijZ = viZ - vel(j, Z);

          /// replace by scalar expressions?
          const fType vijrij = rijX * vijX + rijY * vijY + rijZ * vijZ;

          fType av = 0;

          /// AV
          if (vijrij < 0.)
            {
              const fType rijrij = rijX * rijX + rijY * rijY + rijZ * rijZ;
              const fType rhoij = 0.5 * (rhoi + rhoj);
              const fType cij = 0.5 * (ci + sqrt(gamma * pj / rhoj));
              const fType muij = hij * vijrij / (rijrij + 0.01 * hij * hij);

              av = (-alpha * cij * muij + beta * muij * muij) / rhoij;
            }

          const fType accTerm = piOrhoirhoi + (pj / (rhoj * rhoj)) + av;
          const fType mj = m(j);

          /// acceleration
          Kernel.derive(rij, hij, rijX, rijY, rijZ);

          curAccX -= mj * accTerm * Kernel.derivX;
          curAccY -= mj * accTerm * Kernel.derivY;
          curAccZ -= mj * accTerm * Kernel.derivZ;

#ifdef SPHLATCH_VELDIV
          ///
          /// m_j * v_ij * divW_ij
          ///
          const fType mjvijdivWij = mj * (vijX * Kernel.derivX +
                                              vijY * Kernel.derivY +
                                              vijZ * Kernel.derivZ);
#endif

#ifdef SPHLATCH_TIMEDEP_ENERGY
          ///
          /// pdV + AV heating
          ///
          curPow += (0.5 * accTerm * mjvijdivWij);
#endif

#ifdef SPHLATCH_VELDIV
          ///
          /// velocity divergence
          ///
          curVelDiv += mjvijdivWij;
#endif
        }
      acc(i, X) += curAccX;
      acc(i, Y) += curAccY;
      acc(i, Z) += curAccZ;

      ///
      /// shocktube boundary condition
      ///
      if (pos(i, Z) < 5. || pos(i, Z) > 295.)
        acc(i, Z) = 0.;

#ifdef SPHLATCH_TIMEDEP_ENERGY
      dudt(i) = curPow;
#endif
#ifdef SPHLATCH_VELDIV
      divv(i) = curVelDiv / rho(i);
      divvMax = divv(i) > divvMax ? divv(i) : divvMax;
#endif
    }
#ifdef SPHLATCH_VELDIV
  CommManager.max(divvMax);
#endif

  Logger << "SPH sum: acc, pow, divv";

#ifdef SPHLATCH_TIMEDEP_SMOOTHING
  ///
  /// time derivative of smoothing length, limit smoothing
  /// length if necessary
  /// I need    : divv, noneigh
  /// I provide : h, dhdt
  ///
  const fType noNeighOpt = PartManager.attributes["noneigh"];
  const fType noNeighMin = (2. / 3.) * noNeighOpt;
  const fType noNeighMax = (5. / 3.) * noNeighOpt;
  const fType cDivvMax = divvMax;
#endif

  const fType czAtomicLength = CostZone.getAtomicLength();
  for (size_t i = 0; i < noParts; i++)
    {
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
      const fType noNeighCur = static_cast<fType>(noneigh(i));

      const fType A = exp((noNeighCur - noNeighMin) / 5.);
      const fType B = exp((noNeighCur - noNeighMax) / 5.);

      const fType k1 = 1. / (A * (A + (1. / A)));
      const fType k2 = (A / (A + (1. / A)))
                           + 1. / (B * (B + (1. / B))) - 1.;
      const fType k3 = B / (B + (1. / B));

      dhdt(i) = (k1 * cDivvMax - k3 * cDivvMax
                 - k2 * static_cast<fType>(1. / 3.) * divv(i)) * h(i);
#endif
      ///
      /// hard upper limit
      ///
      if (2.5 * h(i) > czAtomicLength)
        {
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
          if (dhdt(i) > 0.)
            dhdt(i) = 0.;
#endif
          h(i) = czAtomicLength / 2.5;
        }
    }
  Logger.stream << "adapted smoothing length (2.5h < " << czAtomicLength
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
                << ", divvmax " << cDivvMax
#endif
                << ")";
  Logger.flushStream();
};


int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  ///
  /// instantate managers
  ///
  io_type& IOManager(io_type::instance());
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  costzone_type& CostZone(costzone_type::instance());
  log_type& Logger(log_type::instance());

  ///
  /// some simulation parameters
  /// attributes will be overwritten, when defined in file
  ///
  PartManager.attributes["gamma"] = 1.4;
  PartManager.attributes["courant"] = 0.3;

  PartManager.attributes["noneigh"] = 50.;
  PartManager.attributes["umin"] = 0.1;

  std::string loadDumpFile = "initial.h5part";
  const fType maxTime = 100.;

  ///
  /// define what we're doing
  ///
  PartManager.useGravity();
  PartManager.useBasicSPH();
  PartManager.useAVMonaghan();
  PartManager.useEnergy();

#ifdef SPHLATCH_TIMEDEP_ENERGY
  PartManager.useTimedepEnergy();
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
  PartManager.useTimedepH();
#endif

  ///
  /// some useful references
  ///
  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);
  valvectRefType u(PartManager.u);
#ifdef SPHLATCH_TIMEDEP_ENERGY
  valvectRefType dudt(PartManager.dudt);
#endif
  valvectRefType h(PartManager.h);
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
  valvectRefType dhdt(PartManager.dhdt);
#endif
#ifdef SPHLATCH_GRAVITY
  valvectRefType eps(PartManager.eps);
#endif

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  int& step(PartManager.step);
  valueRefType time(PartManager.attributes["time"]);

  const size_t myDomain = CommManager.getMyDomain();

  ///
  /// register the quantites to be exchanged
  ///
  CommManager.exchangeQuants.vects += &pos, &vel;
  CommManager.exchangeQuants.scalars += &m, &u, &h;
  CommManager.exchangeQuants.ints += &id;

#ifdef SPHLATCH_GRAVITY
  CommManager.exchangeQuants.scalars += &eps;
#endif

  ///
  /// instantate the MetaIntegrator
  ///
  //sphlatch::VerletMetaIntegrator Integrator(derivate, timestep);
  sphlatch::PredCorrMetaIntegrator Integrator(derivate, timeStep);

  ///
  /// register spatial, energy and smoothing length integration
  ///
  Integrator.regIntegration(pos, vel, acc);
#ifdef SPHLATCH_TIMEDEP_ENERGY
  Integrator.regIntegration(u, dudt);
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
  Integrator.regIntegration(h, dhdt);
#endif

  IOManager.loadDump(loadDumpFile);
  Logger.stream << "loaded " << loadDumpFile;
  Logger.flushStream();

  ///
  /// bootstrap the integrator
  ///
  Integrator.bootstrap();
  Logger << "integrator bootstrapped";

  ///
  /// the integration loop
  ///
  //fType nextSaveTime = time;
  while (time <= maxTime)
    {
      ///
      /// integrate
      ///
      Integrator.integrate();

      Logger.stream << "finished step " << step << ", now at t = " << time;
      Logger.flushStream();
      Logger.zeroRelTime();

      if (myDomain == 0)
        {
          std::cout << "t = " << std::fixed << std::right
                    << std::setw(12) << std::setprecision(6)
                    << time << " (" << step << ")\n";
        }
    }

  MPI::Finalize();
  return EXIT_SUCCESS;
}


