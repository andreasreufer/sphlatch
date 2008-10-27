// some defs

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
//#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
#define SPHLATCH_SINGLEPREC

// enable parallel version
//#define SPHLATCH_PARALLEL

// enable intensive logging for toptree global summation
//#define SPHLATCH_TREE_LOGSUMUPMP

//#define SPHLATCH_RANKSPACESERIALIZE

// enable checking of bounds for the neighbour lists
#define SPHLATCH_CHECKNONEIGHBOURS

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

//#include "communication_manager.h"
//typedef sphlatch::CommunicationManager comm_type;

//#include "costzone.h"
//typedef sphlatch::CostZone costzone_type;

//#include "log_manager.h"
//typedef sphlatch::LogManager log_type;

#include "integrator_verlet.h"
#include "integrator_predcorr.h"

#include <boost/progress.hpp>
#include <vector>

#include "kernel_cubicspline.h"

/// tree stuff
//#include "bhtree.h"

/// neighbour search
#include "rankspace.h"

using namespace sphlatch::vectindices;
using namespace boost::assign;


void derivate()
{
  part_type& PartManager(part_type::instance());
  //comm_type& CommManager(comm_type::instance());
  io_type& IOManager(io_type::instance());
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
//  const size_t myDomain = CommManager.getMyDomain();

/// little helper vector to zero a 3D quantity
  const zerovalvectType zero(3);

/// send ghosts to other domains
/*  CommManager.sendGhosts(pos);
  CommManager.sendGhosts(vel);
  CommManager.sendGhosts(id);
  CommManager.sendGhosts(m);
  CommManager.sendGhosts(h);
  CommManager.sendGhosts(eps);
  Logger << " sent to ghosts: pos, vel, id, m, h, eps";*/


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
          const valueType r = RSSearch.neighDistList[curNeigh];
          const size_t j = RSSearch.neighbourList[curNeigh];

          const valueType hij = 0.5 * (hi + h(j));

          rhoi += m(j) * Kernel.value(r, hij);
        }
      rho(i) = rhoi;
    }
  //Logger << "SPH sum: rho";
  //CommManager.sendGhosts(rho);
  //Logger << " sent to ghosts: rho";

  ///
  /// lower temperature bound
  ///
  const valueType uMin = 1000.;
  for (size_t i = 0; i < noParts; i++)
    {
      if (u(i) < uMin)
        {
          u(i) = uMin;
        }
    }
  //Logger << "assure minimal temperature";
  //CommManager.sendGhosts(u);
  //Logger << " sent to ghosts: u";

  ///
  /// pressure
  ///
  /// I need    : u, rho
  /// I provide : p
  ///
  //const valueType gamma = 1.4;
  const valueType gamma = 1.66666666666666667;
  p = (gamma - 1) * (boost::numeric::ublas::element_prod(u, rho));
  //Logger << "pressure";

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

      particleRowType(acc, i) += curAcc;
      dudt(i) = curPow;
      divv(i) = curVelDiv / rho(i);

      divvMax = divv(i) > divvMax ? divv(i) : divvMax;
    }
  //CommManager.max(divvMax);

  //Logger << "SPH sum: acc, pow, divv";

};


int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  ///
  /// instantate managers
  ///
  io_type& IOManager(io_type::instance());
  part_type& PartManager(part_type::instance());
  //comm_type& CommManager(comm_type::instance());
  //costzone_type& CostZone(costzone_type::instance());
  //log_type& Logger(log_type::instance());

  ///
  /// some simulation parameters
  ///

  //std::string loadDumpFile = "shocktube270k.h5part";
  //std::string loadDumpFile = "shocktube_t30.h5part";
  //std::string loadDumpFile = "test1_small.h5part";
  std::string loadDumpFile = "initial.h5part";
  //const valueType maxTime = 5.;

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
  /*CommManager.exchangeQuants.vects += &pos, &vel;
  CommManager.exchangeQuants.scalars += &m, &u, &h;
  CommManager.exchangeQuants.ints += &id, &noneigh;*/

  IOManager.loadDump(loadDumpFile);
  //Logger.stream << "loaded " << loadDumpFile;
  //Logger.flushStream();

  ///
  /// exchange particles
  ///
  /*CostZone.createDomainPartsIndex();
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
*/


  MPI::Finalize();
  return EXIT_SUCCESS;
}


