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

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  valueRefType time(PartManager.attributes["time"]);
  size_t& step(PartManager.step);
  size_t& substep(PartManager.substep);
  const size_t noParts = PartManager.getNoLocalParts();
  const size_t myDomain = CommManager.getMyDomain();

  ///
  /// don't go further than one smoothing length
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

  ///
  /// energy integration
  ///
  fType dtU = std::numeric_limits<fType>::max();
  for (size_t i = 0; i < noParts; i++)
    {
      const fType absdudti = dudt(i);
      const fType dtUi = (u(i)) / absdudti;

      if (absdudti > 0.)
        {
          dtU = dtUi < dtU ? dtUi : dtU;
        }
    }

  ///
  /// smoothing length integration
  ///
  fType dtH = std::numeric_limits<fType>::max();
  for (size_t i = 0; i < noParts; i++)
    {
      const fType absdhdti = dhdt(i);
      const fType dtHi = (h(i)) / absdhdti;

      if (absdhdti > 0.)
        {
          dtH = dtHi < dtH ? dtHi : dtH;
        }
    }

  ///
  /// CFL condition
  ///
  fType dtCFL = std::numeric_limits<fType>::max();
  const fType gamma = 1.4;
  for (size_t i = 0; i < noParts; i++)
    {
      const fType ci = sqrt(p(i) * gamma / rho(i));
      const fType dtCFLi = h(i) / ci;

      if (ci > 0.)
        {
          dtCFL = dtCFLi < dtCFL ? dtCFLi : dtCFL;
        }
    }

  ///
  /// distance to next saveItrvl
  ///
  //const fType saveItrvl = 0.1;
  const fType saveItrvl = 0.2;
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

  CommManager.min(dtA);
  CommManager.min(dtU);
  CommManager.min(dtH);
  CommManager.min(dtCFL);

  const fType courantNumber = 0.3;
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
  io_type& IOManager(io_type::instance());
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
  valvectRefType eps(PartManager.eps);

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  const size_t noParts = PartManager.getNoLocalParts();
  const size_t noTotParts = noParts + PartManager.getNoGhostParts();
  size_t& step(PartManager.step);
  size_t& substep(PartManager.substep);

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
  CommManager.sendGhosts(eps);
  Logger << " sent to ghosts: pos, vel, id, m, h, eps";

  const fType gravTheta = 0.7;
  //const fType gravConst = 1.0;
  const fType gravConst = 5.e-3;
  PartManager.attributes["gravTheta"] = gravTheta;
  PartManager.attributes["gravConst"] = gravConst;

  sphlatch::BHtree<sphlatch::Quadrupoles> Tree(gravTheta,
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

  for (size_t i = 0; i < noParts; i++)
    {
      const size_t curIdx = Tree.particleOrder[i];
      Tree.calcGravity(curIdx);
    }
  Logger << "gravity calculated";

///
/// define kernel and neighbour search algorithm
///
  sphlatch::CubicSpline3D Kernel;
  sphlatch::Rankspace RSSearch;

  RSSearch.prepare();
  RSSearch.neighbourList.resize(1024);
  RSSearch.neighDistList.resize(1024);

  Logger << "Rankspace prepared";

  for (size_t i = 0; i < noParts; i++)
    {
      ///
      /// find neighbours
      ///
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

      //if (id(i) == 125664)
      /*if (id(i) == 2292)
         {
          std::cerr << "part" << id(i) << ":  "
                    << h(i) << "   " << dhdt(i) << "   "
                    << noneigh(i) << "   "
                    << pos(i, X) << " " << pos(i, Y) << " " << pos(i, Z) << "   " << myDomain << "\n";
         }*/

      ///
      /// SPH density sum
      /// I need    : pos, h, m
      /// I provide : rho
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

  ///
  /// lower temperature bound
  ///
  //const fType uMin = 1000.;
  const fType uMin = .1;
  for (size_t i = 0; i < noParts; i++)
    {
      if (u(i) < uMin)
        {
          u(i) = uMin;
        }
    }
  Logger << "assure minimal temperature";
  CommManager.sendGhosts(u);
  Logger << " sent to ghosts: u";

  ///
  /// pressure
  ///
  /// I need    : u, rho
  /// I provide : p
  ///
  const fType gamma = 1.4;
  //const fType gamma = 1.66666666666666667;
  p = (gamma - 1) * (boost::numeric::ublas::element_prod(u, rho));
  PartManager.attributes["gamma"] = gamma;
  Logger << "pressure";

  ///
  /// acceleration, power and velocity divergence
  ///
  const fType alpha = 1;
  const fType beta = 1;

  valvectType curAcc(3);

  fType curPow = 0., curVelDiv = 0.;
  fType divvMax = std::numeric_limits<fType>::min();

  for (size_t i = 0; i < noParts; i++)
    {
      const fType hi = h(i);
      const fType rhoi = rho(i);
      const fType pi = p(i);

      const fType piOrhoirhoi = pi / (rhoi * rhoi);

      /// find the neighbours
      const fType srchRad = 2. * hi;
      RSSearch.findNeighbours(i, srchRad);

      const size_t noNeighs = RSSearch.neighbourList[0];

      curAcc = zero;
      curPow = 0.;
      curVelDiv = 0.;

      const particleRowType veli(vel, i);
      const particleRowType Ri(pos, i);

      const fType ci = sqrt(gamma * p(i) / rho(i));

      ///
      /// SPH acceleration and specific power sum
      ///
      /// I need    : pos, vel, h, m, rho, u
      /// I provide : acc, dhdu, divv
      ///
      for (size_t curNeigh = 1; curNeigh <= noNeighs; curNeigh++)
        {
          const fType rij = RSSearch.neighDistList[curNeigh];
          const size_t j = RSSearch.neighbourList[curNeigh];

          const fType rhoj = rho(j);
          const fType pj = p(j);

          const fType hij = 0.5 * (hi + h(j));

          const particleRowType velj(vel, j);
          const particleRowType Rj(pos, j);

          /// replace by scalar expressions?
          const fType vijrij =
            boost::numeric::ublas::inner_prod(velj - veli, Rj - Ri);

          fType av = 0;

          /// AV
          if (vijrij < 0.)
            {
              const fType rijrij =
                boost::numeric::ublas::inner_prod(Ri - Rj, Ri - Rj);
              const fType rhoij = 0.5 * (rhoi + rhoj);
              const fType cij = 0.5 * (ci + sqrt(gamma * p(j) / rhoj));
              const fType muij = hij * vijrij / (rijrij + 0.01 * hij * hij);

              av = (-alpha * cij * muij + beta * muij * muij) / rhoij;
            }

          const fType accTerm = piOrhoirhoi + (pj / (rhoj * rhoj)) + av;

          /// acceleration
          curAcc -= (m(j) * accTerm * Kernel.derive(rij, hij, Ri - Rj));

          /// m_j * v_ij * divW_ij
          const fType mjvijdivWij = m(j) * (
            boost::numeric::ublas::inner_prod(veli - velj,
                                              Kernel.derivative));

          /// pdV + AV heating
          curPow += (0.5 * accTerm * mjvijdivWij);

          /// velocity divergence
          curVelDiv += mjvijdivWij;
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

  const fType noNeighMin = (2. / 3.) * static_cast<fType>(noNeighOpt);
  const fType noNeighMax = (5. / 3.) * static_cast<fType>(noNeighOpt);
  const fType cDivvMax = divvMax;

  const fType czAtomicLength = CostZone.getAtomicLength();
  for (size_t i = 0; i < noParts; i++)
    {
      const fType noNeighCur = static_cast<fType>(noneigh(i));

      const fType A = exp((noNeighCur - noNeighMin) / 5.);
      const fType B = exp((noNeighCur - noNeighMax) / 5.);

      const fType k1 = 1. / (A * (A + (1. / A)));
      const fType k2 = (A / (A + (1. / A)))
                           + 1. / (B * (B + (1. / B))) - 1.;
      const fType k3 = B / (B + (1. / B));

      dhdt(i) = (k1 * cDivvMax - k3 * cDivvMax
                 - k2 * static_cast<fType>(1. / 3.) * divv(i)) * h(i);

      ///
      /// hard upper limit
      ///
      /*if (2.5 * h(i) > czAtomicLength)
         {
          dhdt(i) = 0.;
          h(i) = czAtomicLength / 2.5;
         }*/
    }
  Logger.stream << "adapted smoothing length (2.5h < " << czAtomicLength << ")";
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
  ///

  //std::string loadDumpFile = "shocktube270k.h5part";
  //std::string loadDumpFile = "shocktube_t30.h5part";
  //std::string loadDumpFile = "test1_small.h5part";
  std::string loadDumpFile = "initial.h5part";
  //const fType maxTime = 5.;
  const fType maxTime = 200.;

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
  /// "colour" particles
  ///
  size_t noParts = PartManager.getNoLocalParts();
  for (size_t i = 0; i < noParts; i++)
    {
      pos(i, X) -= 0.5;
      pos(i, X) *= 20.;
      pos(i, Y) -= 0.5;
      pos(i, Y) *= 20.;
      pos(i, Z) -= 0.5;
      pos(i, Z) *= 20.;

      h(i) *= 20.;
      eps(i) = 0.3;
      u(i) = 1.;

      if (pos(i, X) > 0)
        id(i) += 1e6;
    }

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
  //fType nextSaveTime = time;
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

      ///
      /// set center of mass to 0 again
      ///
      noParts = PartManager.getNoLocalParts();
      fType comX = 0., comY = 0., comZ = 0., totM = 0.;
      for (size_t i = 0; i < noParts; i++)
        {
          comX += pos(i, X) * m(i);
          comY += pos(i, Y) * m(i);
          comZ += pos(i, Z) * m(i);
          totM += m(i);
        }
      comX /= totM;
      comY /= totM;
      comZ /= totM;

      for (size_t i = 0; i < noParts; i++)
        {
          pos(i, X) -= comX;
          pos(i, Y) -= comY;
          pos(i, Z) -= comZ;
        }

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

