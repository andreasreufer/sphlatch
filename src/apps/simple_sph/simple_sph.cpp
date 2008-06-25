// some defs

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
#define SPHLATCH_SINGLEPREC

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

// tree stuff
#include "bhtree.h"

#include "rankspace.h"

using namespace sphlatch::vectindices;

valueType timeStep()
{
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  log_type& Logger(log_type::instance());

  matrixRefType pos(PartManager.pos);
  //matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  //valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  valvectRefType p(PartManager.p);
  valvectRefType u(PartManager.u);
  valvectRefType rho(PartManager.rho);
  valvectRefType dudt(PartManager.dudt);
  valvectRefType dhdt(PartManager.dhdt);

  idvectRefType id(PartManager.id);

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

  Logger.stream << "dtA: " << dtA
                << " dtU: " << dtU
                << " dtH: " << dtH
                << " dtCFL: " << dtCFL
                << "   dtGlob: " << dtGlob;
  Logger.flushStream();

  return dtGlob;
}

void derivate()
{
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
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

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  const size_t noParts = PartManager.getNoLocalParts();
  size_t& step(PartManager.step);

  const size_t noTotParts = noParts + PartManager.getNoGhostParts();
  const size_t myDomain = CommManager.getMyDomain();

  //const valueType gravTheta = 0.7;
  const valueType gravTheta = 0.7;
  //const valueType gravConst = 6.674e-11;
  const valueType gravConst = 1.e-3;

  //sphlatch::BHtree<sphlatch::Monopoles>   Tree(gravTheta,
  sphlatch::BHtree<sphlatch::Quadrupoles> Tree(gravTheta,
  //sphlatch::BHtree<sphlatch::Octupoles>   Tree(gravTheta,
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
  ///
  /// calculate the accelerations
  ///
  const zerovalvectType zero(3);
  for (size_t i = 0; i < noParts; i++)
    {
      const size_t curIndex = Tree.particleOrder[i];
      particleRowType(acc, curIndex) = zero;
      Tree.calcGravity(curIndex);
      if (id(curIndex) == 4920)
        {
          std::cout << "a_z: " << acc(curIndex, 2) << "\n";
        }
    }
  Logger << "gravity calculated";


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
  Logger << "SPH density sum";

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

  ///
  /// pressure
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
      /// hold both shocktube ends
      ///

/*
      curAcc(0) = 0.;
      curAcc(1) = 0.;

      if (pos(i, 2) < 5. || pos(i, 2) > 95.)
        {
          curAcc(2) = 0.;
          vel(i, 2) = 0.;
        }
*/
      const valueType dist = sqrt( pos(i, 0)*pos(i, 0)
                                 + pos(i, 1)*pos(i, 1)
                                 + pos(i, 2)*pos(i, 2) );

      if ( dist > 400. )
      {
        vel(i, 0) = 0.;
        vel(i, 1) = 0.;
        vel(i, 2) = 0.;
      }

      particleRowType(acc, i) += curAcc;
      dudt(i) = curPow;
      divv(i) = curVelDiv / rho(i);

      divvMax = divv(i) > divvMax ? divv(i) : divvMax;
    }
  CommManager.max(divvMax);

  Logger << "acceleration, power & velocity divergence sum";

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

      //dhdt(i) = ( - 0.333333333333333 * h(i) * divv(i) );

      if (id(i) == 4920)
        {
          std::cout << "a_z: " << acc(i, 2) << "\n";
          std::cout << "z: " << pos(i, 2) << "   rho: " << rho(i) << "   divv: " << divv(i) << "   dhdt: " << dhdt(i) << "    h: " << h(i) << "    noneigh: " << noneigh(i) << "\n";
          //std::cout << "k1 " << k1 << "    k2 " << k2 << "   k3 " << k3 << "\n";
        }
    }
  Logger << "adapted smoothing length";
};

using namespace boost::assign;

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

  //IOManager.loadDump("shocktube.h5part");
  //Logger << "loaded shocktube.h5part";
  
  IOManager.loadDump("random010k.h5part");
  Logger << "loaded random010k.h5part";
  
  const size_t noParts = PartManager.getNoLocalParts();
  valueType xCom = 0., yCom = 0., zCom = 0., mTot = 0.;

  for (size_t i = 0; i < noParts; i++)
  {
    xCom += pos(i, X)*m(i);
    yCom += pos(i, Y)*m(i);
    zCom += pos(i, Z)*m(i);
    mTot += m(i);
    std::cout << pos(i, X) << " " << pos(i, Y) << "\n";
  }
  xCom /= mTot;
  yCom /= mTot;
  zCom /= mTot;

  Logger.stream << "center of mass: ["
                << xCom << ","
                << yCom << ","
                << zCom << "]\n";
  Logger.flushStream();

  for (size_t i = 0; i < noParts; i++)
  {
    pos(i, X) = 20.*( pos(i, X) - xCom );
    pos(i, Y) = 20.*( pos(i, Y) - yCom );
    pos(i, Z) = 20.*( pos(i, Z) - zCom );
    h(i)   *= 20.;
    eps(i) = .5;
    m(i)   = 1.;
    u(i)   = 1.;
  }
  Logger << "data ready";
  std::cout  << "center of mass: " << xCom << " " << yCom << " " << zCom << "\n";

  ///
  /// define the quantities to save in a dump
  ///
  quantsType saveQuants;
  saveQuants.vects += &pos, &vel, &acc;
  saveQuants.scalars += &m, &rho, &u, &p, &h;
  saveQuants.ints += &id, &noneigh;

  ///
  /// exchange particles
  ///
  CostZone.createDomainPartsIndex();
  CommManager.exchange(CostZone.domainPartsIndex,
                       CostZone.getNoGhosts());

  ///
  /// prepare ghost sends and send ghosts
  ///
  CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());

  CommManager.sendGhosts(pos);
  CommManager.sendGhosts(id);
  CommManager.sendGhosts(m);
  CommManager.sendGhosts(rho);
  CommManager.sendGhosts(u);


  Logger.stream << "distributed particles: "
                << PartManager.getNoLocalParts() << " parts. & "
                << PartManager.getNoGhostParts() << " ghosts";
  Logger.flushStream();

  ///
  /// define some simulation parameters
  ///
  //const valueType dt = 138.e6 / 50.;

  //const size_t maxSteps = 500000;
  //const size_t saveSteps =  5000;
  const size_t maxStep = 3000;
  const size_t saveSteps = 20;

  ///
  /// bootstrap the integrator
  ///
  Integrator.bootstrap();
  Logger << "integrator bootstrapped";

  //MPI::Finalize();
  //return EXIT_SUCCESS;
  ///
  /// the integration loop
  ///
  while (step < maxStep)
    {
      ///
      /// exchange particles and ghosts
      ///
      CostZone.createDomainPartsIndex();
      CommManager.exchange(CostZone.domainPartsIndex,
                           CostZone.getNoGhosts());

      CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());

      CommManager.sendGhosts(pos);
      CommManager.sendGhosts(id);
      CommManager.sendGhosts(m);

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


      std::cout << "t = " << std::fixed << std::right
                << std::setw(12) << std::setprecision(6)
                << time << " (" << step << ")\n";

      ///
      /// save a dump
      ///
      if ((step % saveSteps) == 0)
        {
          if (myDomain == 0)
            {
              std::cout << "save dump!\n";
            }
          IOManager.saveDump("saveDump.h5part", saveQuants);
        }
    }

  MPI::Finalize();
  return EXIT_SUCCESS;
}


