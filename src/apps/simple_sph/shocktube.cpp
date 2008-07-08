// some defs

#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
//#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

// enable parallel version
#define SPHLATCH_PARALLEL

// enable intensive logging for toptree global summation
//#define SPHLATCH_TREE_LOGSUMUPMP

//#define SPHLATCH_RANKSPACESERIALIZE


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
typedef sphlatch::valueType valueType;
typedef sphlatch::valueRefType valueRefType;
typedef sphlatch::valvectType valvectType;
typedef sphlatch::zerovalvectType zerovalvectType;
typedef sphlatch::valvectRefType valvectRefType;

typedef sphlatch::particleRowType particleRowType;

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

#include "kernel_cubicspline3d.h"

/// tree stuff
#include "bhtree.h"

/// neighbour search
#include "rankspace.h"

using namespace sphlatch::vectindices;
using namespace boost::assign;

valueType timeStep()
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

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  valueRefType time(PartManager.attributes["time"]);
  size_t& step(PartManager.step);


  const size_t noParts = PartManager.getNoLocalParts();

  ///
  /// don't go further than one smoothing length
  ///
  valueType dtA = std::numeric_limits<valueType>::max();

  for (size_t i = 0; i < noParts; i++)
    {
      const valueType ai = sqrt(acc(i, X) * acc(i, X) +
                                acc(i, Y) * acc(i, Y) +
                                acc(i, Z) * acc(i, Z));

      if (ai > 0.)
        {
          const valueType dtAi = sqrt(h(i) / ai);
          dtA = dtAi < dtA ? dtAi : dtA;
        }
    }

  ///
  /// energy integration
  ///
  valueType dtU = std::numeric_limits<valueType>::max();
  for (size_t i = 0; i < noParts; i++)
    {
      const valueType absdudti = dudt(i);
      const valueType dtUi = (u(i)) / absdudti;

      if (absdudti > 0.)
        {
          dtU = dtUi < dtU ? dtUi : dtU;
        }
    }

  ///
  /// smoothing length integration
  ///
  valueType dtH = std::numeric_limits<valueType>::max();
  for (size_t i = 0; i < noParts; i++)
    {
      const valueType absdhdti = dhdt(i);
      const valueType dtHi = (h(i)) / absdhdti;

      if (absdhdti > 0.)
        {
          dtH = dtHi < dtH ? dtHi : dtH;
        }
    }

  ///
  /// CFL condition
  ///
  valueType dtCFL = std::numeric_limits<valueType>::max();
  const valueType gamma = 1.4;
  for (size_t i = 0; i < noParts; i++)
    {
      const valueType ci = sqrt(p(i) * gamma / rho(i));
      const valueType dtCFLi = h(i) / ci;

      if (ci > 0.)
        {
          dtCFL = dtCFLi < dtCFL ? dtCFLi : dtCFL;
        }
    }

  ///
  /// distance to next saveItrvl
  ///
  const valueType saveItrvl = 10.;
  std::string fileName = "saveDump000000.h5part";

  const valueType nextSaveTime = (floor((time / saveItrvl) + 1.e-6)
                                  + 1.) * saveItrvl;
  valueType dtSave = nextSaveTime - time;

  ///
  /// determine global minimum.
  /// by parallelly minimizing the timesteps, we
  /// can estimate which ones are dominant
  ///
  valueType dtGlob = std::numeric_limits<valueType>::max();

  CommManager.min(dtA);
  CommManager.min(dtU);
  CommManager.min(dtH);
  CommManager.min(dtCFL);

  const valueType courantNumber = 0.3;
  dtGlob = dtA < dtGlob ? dtA : dtGlob;
  dtGlob = dtU < dtGlob ? dtU : dtGlob;
  dtGlob = dtH < dtGlob ? dtH : dtGlob;
  dtGlob = dtCFL < dtGlob ? dtCFL : dtGlob;
  dtGlob *= courantNumber;

  dtGlob = dtSave < dtGlob ? dtSave : dtGlob;

  Logger.stream << "dtA: " << dtA
                << " dtU: " << dtU
                << " dtH: " << dtH
                << " dtCFL: " << dtCFL
                << " dtSave: " << dtSave
                << "   dtGlob: " << dtGlob;
  Logger.flushStream();

  ///
  /// define the quantities to save in a dump
  ///
  quantsType saveQuants;
  saveQuants.vects += &pos, &vel, &acc;
  saveQuants.scalars += &m, &rho, &u, &p, &h, &divv, &dudt, &dhdt;
  saveQuants.ints += &id, &noneigh;


  const valueType curSaveTime = (floor((time / saveItrvl) + 1.e-9) )
                                * saveItrvl;
  if ( fabs(curSaveTime - time) < 1.e-9)
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
  //costzone_type& CostZone(costzone_type::instance());
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

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  const size_t noParts = PartManager.getNoLocalParts();
  const size_t noTotParts = noParts + PartManager.getNoGhostParts();
  
  //size_t& step(PartManager.step);
  const size_t myDomain = CommManager.getMyDomain();

  /// little helper vector to zero a 3D quantity
  const zerovalvectType zero(3);

  /// send ghosts to other domains
  CommManager.sendGhosts(pos);
  CommManager.sendGhosts(vel);
  CommManager.sendGhosts(id);
  CommManager.sendGhosts(m);
  CommManager.sendGhosts(h);

  ///
  /// boundary conditions
  ///
  for (size_t i = 0; i < noTotParts; i++)
    {
      vel(i, 0) = 0.;
      vel(i, 1) = 0.;
      
      if (pos(i, 2) < 5. || pos(i, 2) > 295.)
        {
          vel(i, 2) = 0.;
        }
    }

/*
   const valueType gravTheta = 1.0;
   const valueType gravConst = 0.0;

   sphlatch::BHtree<sphlatch::NoMultipoles> Tree(gravTheta,
                                               gravConst,
                                               CostZone.getDepth(),
                                               CostZone.getCenter(),
                                               CostZone.getSidelength()
                                               );


   ///
   /// fill up tree, determine ordering, calculate multipoles
   ///
   for (size_t i = 0; i < noTotParts; i++)
    {
      Tree.insertParticle(i);
    }

   Tree.detParticleOrder();
   Tree.calcMultipoles();
   Logger << "Tree ready";
 */
///
/// define kernel and neighbour search algorithm
///
  sphlatch::CubicSpline3D Kernel;
  sphlatch::Rankspace RSSearch;

  RSSearch.prepare();
  RSSearch.neighbourList.resize(512);
  RSSearch.neighDistList.resize(512);

  Logger << "Rankspace prepared";

  for (size_t i = 0; i < noParts; i++)
    {
      ///
      /// find neighbours
      ///
      const valueType hi = h(i);
      const valueType srchRad = 2. * hi;
      RSSearch.findNeighbours(i, srchRad);

      const size_t noNeighs = RSSearch.neighbourList[0];

      ///
      /// store the number of neighbours
      ///
      noneigh(i) = noNeighs;

      static valueType rhoi;
      rhoi = 0.;

      ///
      /// SPH density sum
      /// I need    : pos, h, m
      /// I provide : rho
      ///
      for (size_t curNeigh = 1; curNeigh <= noNeighs; curNeigh++)
        {
          const valueType r = RSSearch.neighDistList[curNeigh];
          const size_t j = RSSearch.neighbourList[curNeigh];

          const valueType hij = 0.5 * (hi + h(j));

          rhoi += m(j) * Kernel.value(r, hij);
        }
      rho(i) = rhoi;
    }
  Logger << "SPH sum: rho";
  CommManager.sendGhosts(rho);
  Logger << " density sent to ghosts";

  ///
  /// lower temperature bound
  ///
  const valueType uMin = 0.1;
  for (size_t i = 0; i < noParts; i++)
    {
      if (u(i) < uMin)
        {
          u(i) = uMin;
        }
    }
  Logger << "assure minimal temperature";
  CommManager.sendGhosts(u);
  Logger << " temperature sent to ghosts";

  ///
  /// pressure
  ///
  /// I need    : u, rho
  /// I provide : p
  ///
  const valueType gamma = 1.4;
  p = (gamma - 1) * (boost::numeric::ublas::element_prod(u, rho));
  Logger << "pressure";

  ///
  /// acceleration, power and velocity divergence
  ///
  const valueType alpha = 1;
  const valueType beta = 1;

  valvectType curAcc(3);

  valueType curPow = 0., curVelDiv = 0.;
  valueType divvMax = std::numeric_limits<valueType>::min();

  for (size_t i = 0; i < noParts; i++)
    {
      const valueType hi = h(i);
      const valueType rhoi = rho(i);
      const valueType pi = p(i);

      const valueType piOrhoirhoi = pi / (rhoi * rhoi);

      /// find the neighbours
      const valueType srchRad = 2. * hi;
      RSSearch.findNeighbours(i, srchRad);

      const size_t noNeighs = RSSearch.neighbourList[0];

      curAcc = zero;
      curPow = 0.;
      curVelDiv = 0.;

      const particleRowType veli(vel, i);
      const particleRowType Ri(pos, i);

      const valueType ci = sqrt(gamma * p(i) / rho(i));

      ///
      /// SPH acceleration and specific power sum
      ///
      /// I need    : pos, vel, h, m, rho, u
      /// I provide : acc, dhdu, divv
      ///
      for (size_t curNeigh = 1; curNeigh <= noNeighs; curNeigh++)
        {
          const valueType rij = RSSearch.neighDistList[curNeigh];
          const size_t j = RSSearch.neighbourList[curNeigh];

          const valueType rhoj = rho(j);
          const valueType pj = p(j);

          const valueType hij = 0.5 * (hi + h(j));

          const particleRowType velj(vel, j);
          const particleRowType Rj(pos, j);

          /// replace by scalar expressions?
          const valueType vijrij =
            boost::numeric::ublas::inner_prod(velj - veli, Rj - Ri);

          valueType av = 0;

          /// AV
          if (vijrij < 0.)
            {
              const valueType rijrij =
                boost::numeric::ublas::inner_prod(Ri - Rj, Ri - Rj);
              const valueType rhoij = 0.5 * (rhoi + rhoj);
              const valueType cij = 0.5 * (ci + sqrt(gamma * p(j) / rhoj));
              const valueType muij = hij * vijrij / (rijrij + 0.01 * hij * hij);

              av = (-alpha * cij * muij + beta * muij * muij) / rhoij;
            }

          const valueType accTerm = piOrhoirhoi + (pj / (rhoj * rhoj)) + av;

          /// acceleration
          curAcc -= (m(j) * accTerm * Kernel.derive(rij, hij, Ri - Rj));

          /// m_j * v_ij * divW_ij
          const valueType mjvijdivWij = m(j) * (
            boost::numeric::ublas::inner_prod(veli - velj,
                                              Kernel.derivative));

          /// pdV + AV heating
          curPow += (0.5 * accTerm * mjvijdivWij);

          /// velocity divergence
          curVelDiv += mjvijdivWij;
        }

      ///
      /// a_x = a_y = 0.
      /// hold both shocktube ends:
      /// a_z = 0. if z < z_min || z > z_max
      /// z_min = 5.
      /// z_max = 295.
      ///
      curAcc(0) = 0.;
      curAcc(1) = 0.;
      if (pos(i, 2) < 5. || pos(i, 2) > 295.)
        {
          curAcc(2) = 0.;
        }

      particleRowType(acc, i) += curAcc;
      dudt(i) = curPow;
      divv(i) = curVelDiv / rho(i);

      divvMax = divv(i) > divvMax ? divv(i) : divvMax;
    }
  CommManager.max(divvMax);

  Logger << "SPH sum: acc, pow, divv";

  /// define desired number of neighbours
  const size_t noNeighOpt = 50;

  const valueType noNeighMin = (2. / 3.) * static_cast<valueType>(noNeighOpt);
  const valueType noNeighMax = (5. / 3.) * static_cast<valueType>(noNeighOpt);
  const valueType cDivvMax = divvMax;

  for (size_t i = 0; i < noParts; i++)
    {
      const valueType noNeighCur = static_cast<valueType>(noneigh(i));

      const valueType A = exp((noNeighCur - noNeighMin) / 5.);
      const valueType B = exp((noNeighCur - noNeighMax) / 5.);

      const valueType k1 = 1. / (A * (A + (1. / A)));
      const valueType k2 = (A / (A + (1. / A)))
                           + 1. / (B * (B + (1. / B))) - 1.;
      const valueType k3 = B / (B + (1. / B));

      dhdt(i) = (k1 * cDivvMax - k3 * cDivvMax
                 - k2 * static_cast<valueType>(1. / 3.) * divv(i)) * h(i);
    }
  Logger << "adapted smoothing length";
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
  ///

  std::string loadDumpFile = "shocktube270k.h5part";
  //std::string loadDumpFile = "shocktube_t30.h5part";

  const valueType maxTime = 50.;
  const valueType saveItrvl = 10.;

  ///
  /// define what we're doing
  ///
  PartManager.useGravity();
  PartManager.useBasicSPH();
  PartManager.useAVMonaghan();
  PartManager.useEnergy();
  PartManager.useTimedepEnergy();
  PartManager.useTimedepH();

  ///
  /// some useful references
  ///
  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);
  valvectRefType rho(PartManager.rho);
  valvectRefType u(PartManager.u);
  valvectRefType dudt(PartManager.dudt);
  valvectRefType h(PartManager.h);
  valvectRefType dhdt(PartManager.dhdt);
  valvectRefType p(PartManager.p);
  valvectRefType eps(PartManager.eps);
  valvectRefType divv(PartManager.divv);

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  size_t& step(PartManager.step);
  valueRefType time(PartManager.attributes["time"]);

  const size_t myDomain = CommManager.getMyDomain();

  ///
  /// register the quantites to be exchanged
  ///
  CommManager.exchangeQuants.vects += &pos, &vel;
  CommManager.exchangeQuants.scalars += &m, &u, &h;
  CommManager.exchangeQuants.ints += &id, &noneigh;

  ///
  /// instantate the MetaIntegrator
  ///

  //sphlatch::VerletMetaIntegrator Integrator(derivate, timestep);
  sphlatch::PredCorrMetaIntegrator Integrator(derivate, timeStep);

  ///
  /// register spatial, energy and smoothing length integration
  ///
  Integrator.regIntegration(pos, vel, acc);
  Integrator.regIntegration(u, dudt);
  Integrator.regIntegration(h, dhdt);

  IOManager.loadDump(loadDumpFile);
  Logger.stream << "loaded " << loadDumpFile;
  Logger.flushStream();

  ///
  /// exchange particles
  ///
  CostZone.createDomainPartsIndex();
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

  ///
  /// bootstrap the integrator
  ///
  Integrator.bootstrap();
  Logger << "integrator bootstrapped";

  ///
  /// the integration loop
  ///
  //valueType nextSaveTime = time;
  while (time <= maxTime)
    {
      ///
      /// exchange particles and ghosts
      /// \todo: only do this, when necessary
      ///
      CostZone.createDomainPartsIndex();
      CommManager.exchange(CostZone.domainPartsIndex,
                           CostZone.getNoGhosts());

      CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());
      Logger.stream << "distributed particles: "
                    << PartManager.getNoLocalParts() << " parts. & "
                    << PartManager.getNoGhostParts() << " ghosts";
      Logger.flushStream();

      ///
      /// integrate
      ///
      Integrator.integrate();

      Logger.stream << "finished step " << step << ", now at t = " << time;
      Logger.flushStream();
      Logger.zeroRelTime();
      step++;

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


