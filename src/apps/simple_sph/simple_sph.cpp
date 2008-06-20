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
#include <math.h>

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

#include "integrator_verlet.h"
#include "integrator_predcorr.h"

#include <boost/progress.hpp>
#include <vector>

#include "kernel_cubicspline3d.h"

// tree stuff
#include "bhtree.h"

#include "rankspace.h"

valueType timeStep()
{
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());

  //matrixRefType pos(PartManager.pos);
  //matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  //valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  valvectRefType p(PartManager.p);
  valvectRefType u(PartManager.u);
  valvectRefType rho(PartManager.rho);
  valvectRefType dudt(PartManager.dudt);

  //idvectRefType id(PartManager.id);

  const size_t noParts = PartManager.getNoLocalParts();

  ///
  /// don't go further than one smoothing length
  ///
  valueType dtA = VAL_MAX;

  for (size_t i = 0; i < noParts; i++)
    {
      const valueType ai = sqrt(acc(i, 0) * acc(i, 0) +
                                acc(i, 1) * acc(i, 1) +
                                acc(i, 2) * acc(i, 2));

      if (ai > 0.)
        {
          const valueType dtAi = sqrt(h(i) / ai);
          dtA = dtAi < dtA ? dtAi : dtA;
        }
    }

  ///
  /// energy integration
  ///
  const valueType uMin = 0.1;
  valueType dtU = VAL_MAX;
  for (size_t i = 0; i < noParts; i++)
    {
      const valueType absdudti = dudt(i);

      const valueType dtUi = (u(i) + uMin) / absdudti;

      if ( absdudti > 0. )
      {
        dtU = dtUi < dtU ? dtUi : dtU;
      }
    }

  ///
  /// CFL condition
  ///
  valueType dtCFL = VAL_MAX;
  const valueType gamma = 1.4;
  for (size_t i = 0; i < noParts; i++)
    {
      const valueType ci = sqrt(p(i) * gamma / rho(i));
      const valueType dtCFLi = h(i) / ci;

      if ( ci > 0. )
      {
        dtCFL = dtCFLi < dtCFL ? dtCFLi : dtCFL;
      }
    }

  ///
  /// determine global minimum
  /// by parallel minimizing the timesteps, we
  /// can estimate which ones are dominant
  ///
  valueType dtGlob = VAL_MAX;

  CommManager.min(dtA);
  CommManager.min(dtU);
  CommManager.min(dtCFL);
  
  dtGlob = dtA < dtGlob ? dtA : dtGlob;
  dtGlob = dtU < dtGlob ? dtU : dtGlob;
  dtGlob = dtCFL < dtGlob ? dtCFL : dtGlob;

  return dtGlob;
}

void derivate()
{
  part_type& PartManager(part_type::instance());
  //comm_type& CommManager(comm_type::instance());
  //costzone_type& CostZone(costzone_type::instance());

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  valvectRefType p(PartManager.p);
  valvectRefType u(PartManager.u);
  valvectRefType rho(PartManager.rho);
  valvectRefType dudt(PartManager.dudt);

  //idvectRefType id(PartManager.id);

  const size_t noParts = PartManager.getNoLocalParts();
  //const size_t noTotParts = noParts + PartManager.getNoGhostParts();
  //const size_t myDomain = CommManager.getMyDomain();
/*
   //const valueType gravTheta = 0.7;
   const valueType gravTheta = 0.7;
   const valueType gravConst = 6.674e-11;

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

   ///
   /// calculate the accelerations
   ///
   for (size_t i = 0; i < noParts; i++)
   {
    const size_t curIndex = Tree.particleOrder[i];
    ///
    /// don't calculate the acceleration for particle with ID = 1 (star)
    ///
    if ( lrint( id(curIndex) ) != 1 )
    {
      Tree.calcGravity( curIndex );
    }
   }
 */

  ///
  /// define kernel and neighbour search algorithm
  ///
  sphlatch::CubicSpline3D Kernel;
  sphlatch::Rankspace RSSearch;
  RSSearch.prepare();

  for (size_t i = 0; i < noParts; i++)
    {
      ///
      /// find neighbours
      ///
      const valueType hi = h(i);
      const valueType srchRad = 2. * hi;
      RSSearch.findNeighbours(i, srchRad);

      const size_t noNeighs = RSSearch.neighbourList[0];

      if (noNeighs < 10)
        std::cerr << noNeighs << " ";

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

          /*if ( fpclassify( rhoi ) == FP_NAN )
             {
             std::cerr << "nan in density sum  i:" << i << " j: " << j << " (Kernel: " <<  Kernel.value(r, hij) << " hij: " << hij << " r: " << r << ")\n";
             }*/
        }
      rho(i) = rhoi;
    }
  std::cerr << "density queue finished!\n";

  const valueType gamma = 1.4;

  ///
  /// lower temperature bound
  ///
  for (size_t i = 0; i < noParts; i++)
    {
      if (u(i) < 0.1)
        {
          u(i) = 0.1;
        }
    }
  p = (gamma - 1) * (boost::numeric::ublas::element_prod(u, rho));


  const valueType alpha = 1;
  const valueType beta = 1;

  valvectType curAcc(3);
  zerovalvectType zero(3);

  valueType curPow = 0.;

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
          const valueType rhoij = 0.5 * (rhoi + rhoj);
          const valueType cij = 0.5 * (ci + sqrt(gamma * p(j) / rhoj));

          const particleRowType velj(vel, j);
          const particleRowType Rj(pos, j);

          /// replace by scalar expressions?
          const valueType vijrij =
            boost::numeric::ublas::inner_prod(velj - veli, Rj - Ri);

          const valueType rijrij =
            boost::numeric::ublas::inner_prod(Ri - Rj, Ri - Rj);

          /// AV
          const valueType muij = hij * vijrij / (rijrij + 0.01 * hij * hij);
          const valueType av = (-alpha * cij * muij + beta * muij * muij) / rhoij;

          const valueType accTerm = piOrhoirhoi + (pj / (rhoj * rhoj)) + av;

          /// acceleration
          curAcc += (m(j) * accTerm * Kernel.derive(rij, hij, Ri - Rj));


          /// pdV + AV heating
          curPow += m(j) * accTerm *
                    boost::numeric::ublas::inner_prod(veli - velj,
                                                      Kernel.derivative); /// check sign!
          /*if ( fpclassify( curPow ) == FP_NAN )
             {
             std::cerr << "nan in power sum  i:" << i << " j: " << j << " accTerm: " << accTerm << " (Kernel deriv: " << boost::numeric::ublas::inner_prod(veli - velj,Kernel.derivative) << " hij: " << hij << " R: " << Ri - Rj << " Ri: " << Ri << " Rj: " << Rj << ")\n";
             std::cerr << "  piOrhoirhoi: " << piOrhoirhoi << " pj: " << pj << " rhoj: " << rhoj << " av: " << av << "\n";
             std::cerr << "  av stuff  cij: " << cij << " muij: " << muij << " rhoi: " << rhoi << "\n";
             std::cerr << "  cij stuff  ci: " << ci << " rhoj: " << rhoj << " gamma: " << gamma << " p(j): " << p(j) << " sqrt(gamma * p(j) / rhoj): " << sqrt(gamma * p(j) / rhoj) << " u(j): " << u(j) << "\n";

             for (size_t k = 0; k < noParts; k++)
             {
              //if (pos(k, 0) < 8 && pos(k, 0) > 2 &&
              //    pos(k, 1) < 8 && pos(k, 1) > 2)
                {
                  std::cout << pos(k, 2) << "\t" << rho(k) << "\t" << p(k) << "\t" << dudt(k) << "\t" << vel(k, 2) << "\t" << u(k) << "\n";
                }

             }

             MPI::Finalize();
             exit(EXIT_FAILURE);
             }*/
        }

      curAcc(0) = 0.;
      curAcc(1) = 0.;

      if (pos(i, 2) < 5. || pos(i, 2) > 95.)
        {
          curAcc(2) = 0.;
        }

      particleRowType(acc, i) = curAcc;
      dudt(i) = curPow;
    }
  std::cerr << "acc&pow queue finished!\n";
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

  ///
  /// define what we're doing
  ///
  PartManager.useGravity();
  PartManager.useBasicSPH();
  PartManager.useEnergy();
  PartManager.useTimedepEnergy();

  ///
  /// some useful references
  ///
  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);
  valvectRefType rho(PartManager.rho);
  valvectRefType u(PartManager.u);
  //valvectRefType dudt(PartManager.dudt);
  valvectRefType p(PartManager.p);

  idvectRefType id(PartManager.id);

  size_t& step(PartManager.step);

  const size_t myDomain = CommManager.getMyDomain();

  ///
  /// register the quantites to be exchanged
  ///
  CommManager.exchangeQuants.vects += &pos, &vel;
  CommManager.exchangeQuants.scalars += &m, &u;
  CommManager.exchangeQuants.ints += &id;

  ///
  /// instantate the MetaIntegrator and
  /// register spatial integration
  ///

  //sphlatch::VerletMetaIntegrator Integrator(derivate, timestep);
  sphlatch::PredCorrMetaIntegrator Integrator(derivate, timeStep);
  Integrator.regIntegration(pos, vel, acc);
  //Integrator.regIntegration(u, dudt);

  IOManager.loadDump("shocktube.hdf5");

  ///
  /// define the quantities to save in a dump
  ///
  quantsType saveQuants;
  saveQuants.vects += &pos, &vel, &acc;
  saveQuants.scalars += &m, &rho, &u, &p;
  saveQuants.ints += &id;

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

  ///
  /// define some simulation parameters
  ///
  //const valueType dt = 138.e6 / 50.;

  //const size_t maxSteps = 500000;
  //const size_t saveSteps =  5000;
  const size_t maxStep = 300;
  const size_t saveSteps = 30;

  ///
  /// bootstrap the integrator
  ///
  Integrator.bootstrap();

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

      ///
      /// integrate
      ///
      Integrator.integrate();
      step++;

      std::cerr << "now at step " << step << "\n";

      ///
      /// save a dump
      ///
      if ((step % saveSteps) == 0)
        {
          if (myDomain == 0)
            {
              std::cerr << "save dump!\n";
            }

          /*const size_t noParts = PartManager.getNoLocalParts();
             std::cout << "#z\trho\tp\tdudt\tvz\n";
             for (size_t i = 0; i < noParts; i++)
             {
              if (pos(i, 0) < 8 && pos(i, 0) > 2 &&
                  pos(i, 1) < 8 && pos(i, 1) > 2)
                {
                  std::cout << pos(i, 2) << "\t" << rho(i) << "\t" << p(i) << "\t" << dudt(i) << "\t" << vel(i, sphlatch::Z) << "\n";
                }
             }*/
          IOManager.saveDump("saveDump.h5part", saveQuants);
        }
    }

  MPI::Finalize();
  return EXIT_SUCCESS;
}


