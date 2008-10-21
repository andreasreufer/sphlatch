// some defs

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

// enable parallel version
#define SPHLATCH_PARALLEL

// enable Logger
#define SPHLATCH_LOGGER

// enable intensive logging for toptree global summation
//#define SPHLATCH_TREE_LOGSUMUPMP

//#define SPHLATCH_RANKSPACESERIALIZE

// enable checking of bounds for the neighbour lists
//#define SPHLATCH_CHECKNONEIGHBOURS

// enable selfgravity?
//#define SPHLATCH_GRAVITY

// time dependent energy?
//#define SPHLATCH_TIMEDEP_ENERGY

// time dependent smoothing length?
//#define SPHLATCH_TIMEDEP_SMOOTHING

// shocktube boundary conditions?
//#define SPHLATCH_SHOCKTUBE

// integrate rho instead of using the SPH sum?
//#define SPHLATCH_INTEGRATERHO

// linear velocity damping term?
//#define SPHLATCH_FRICTION

// remove particles on escaping orbits?
//#define SPHLATCH_REMOVEESCAPING

// do we need the velocity divergence?
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

#ifdef SPHLATCH_INTEGRATERHO
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
typedef sphlatch::valueType fType;
typedef sphlatch::valueRefType valueRefType;

typedef sphlatch::valvectType valvectType;
typedef sphlatch::valvectRefType valvectRefType;

typedef sphlatch::countsType countsType;

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

//#include "integrator_verlet.h"
//typedef sphlatch::VerletMetaIntegrator integrator_type;

#include "integrator_predcorr.h"
typedef sphlatch::PredCorrMetaIntegrator integrator_type;

#ifdef SPHLATCH_TILLOTSON
#include "eos_tillotson.h"
#define SPHLATCH_EOS_DEFINED
typedef sphlatch::Tillotson eos_type;
#endif

#ifdef SPHLATCH_ANEOS
#include "eos_aneos.h"
#define SPHLATCH_EOS_DEFINED
typedef sphlatch::ANEOS eos_type;
#endif

#ifndef SPHLATCH_EOS_DEFINED
#include "eos_idealgas.h"
#define SPHLATCH_EOS_DEFINED
typedef sphlatch::IdealGas eos_type;
#endif

#include <boost/progress.hpp>
#include <vector>

#include "kernel_cubicspline3d.h"
typedef sphlatch::CubicSpline3D kernel_type;

#include "bhtree.h"
typedef sphlatch::BHtree<sphlatch::Quadrupoles> tree_type;

#include "rankspace.h"
typedef sphlatch::Rankspace neighsearch_type;

#include "lookup_table1D.h"
typedef sphlatch::LookupTable1D<sphlatch::InterpolateLinear> lut_type;

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
  valvectRefType cs(PartManager.cs);
  valvectRefType rho(PartManager.rho);
  valvectRefType dudt(PartManager.dudt);
  valvectRefType dhdt(PartManager.dhdt);
  valvectRefType divv(PartManager.divv);

#ifdef SPHLATCH_GRAVITY
  valvectRefType eps(PartManager.eps);
#endif
#ifdef SPHLATCH_INTEGRATERHO
  valvectRefType drhodt(PartManager.drhodt);
#endif
  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);
#ifdef SPHLATCH_TILLOTSON
  idvectRefType mat(PartManager.mat);
#endif
#ifdef SPHLATCH_ANEOS
  idvectRefType mat(PartManager.mat);
  idvectRefType phase(PartManager.phase);
  valvectRefType T(PartManager.T);
#endif
#ifdef SPHLATCH_REMOVEESCAPING
  valvectRefType ecc(PartManager.ecc);
#endif

  valueRefType time(PartManager.attributes["time"]);
  valueRefType saveItrvl(IOManager.saveItrvl);

  int& step(PartManager.step);
  //int& substep(PartManager.substep);

  const size_t noParts = PartManager.getNoLocalParts();
  const size_t myDomain = CommManager.getMyDomain();

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
  const fType courantNumber = PartManager.attributes["courant"];
  for (size_t i = 0; i < noParts; i++)
    {
      const fType dtCFLi = h(i) / cs(i);
      dtCFL = dtCFLi < dtCFL ? dtCFLi : dtCFL;
      if (cs(i) > 0.)
        {
          dtCFL = dtCFLi < dtCFL ? dtCFLi : dtCFL;
        }
    }
  dtCFL *= courantNumber;
  CommManager.min(dtCFL);

  ///
  /// distance to next save time
  ///
  const fType dtSave = (floor((time / saveItrvl) + 1.e-6)
                            + 1.) * saveItrvl - time;

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

#ifdef SPHLATCH_REMOVEESCAPING
  ///
  /// estimate center of protoearth by calculating
  /// the center of mass of all particles with higher than
  /// minimal core density
  ///
  const fType minCoreDensity = PartManager.attributes["mincoredens"];
  fType comX = 0., comY = 0., comZ = 0., comM = 0.;
  for (size_t i = 0; i < noParts; i++)
    {
      if (rho(i) > minCoreDensity)
        {
          comX += pos(i, X) * m(i);
          comY += pos(i, Y) * m(i);
          comZ += pos(i, Z) * m(i);
          comM += m(i);
        }
    }
  CommManager.sum(comX);
  CommManager.sum(comY);
  CommManager.sum(comZ);
  CommManager.sum(comM);

  comX /= comM;
  comY /= comM;
  comZ /= comM;

  if (myDomain == 0)
    {
      std::fstream comFile;
      comFile.open("centerOfMass", std::ios::out | std::ios::app);
      comFile << time << "\t"
              << comX << "\t"
              << comY << "\t"
              << comZ << "\t"
              << comM << "\n";
      comFile.close();
    }

  ///
  /// determine maximal radius relative to center of mass
  ///
  fType maxRad = 0.;
  const size_t noMassBins = 10000;
  for (size_t i = 0; i < noParts; i++)
    {
      const fType curRad = sqrt((pos(i, X) - comX) * (pos(i, X) - comX) +
                                    (pos(i, Y) - comY) * (pos(i, Y) - comY) +
                                    (pos(i, Z) - comZ) * (pos(i, Z) - comZ));
      maxRad = curRad > maxRad ? curRad : maxRad;
    }
  CommManager.max(maxRad);

  ///
  /// bin particle masses according to radius, cumulate mass bins
  /// and sum them up globally
  ///
  valvectType massRad(noMassBins), massBins(noMassBins);
  for (size_t i = 0; i < noMassBins; i++)
    {
      massRad(i) = maxRad * (i / static_cast<fType>(noMassBins));
      massBins(i) = 0.;
    }

  for (size_t i = 0; i < noParts; i++)
    {
      const fType curRad = sqrt((pos(i, X) - comX) * (pos(i, X) - comX) +
                                    (pos(i, Y) - comY) * (pos(i, Y) - comY) +
                                    (pos(i, Z) - comZ) * (pos(i, Z) - comZ));

      const size_t curBin = std::max(std::min(
                                       lrint((curRad / maxRad) * noMassBins
                                             - 0.5),
                                       static_cast<long int>(noMassBins)),
                                     static_cast<long int>(0));
      massBins(curBin) += m(i);
    }

  for (size_t i = 1; i < noMassBins; i++)
    {
      massBins(i) += massBins(i - 1);
    }
  CommManager.sum(massBins);

  ///
  /// create lookup table
  ///
  lut_type massProf(massRad, massBins);

  ///
  /// finally blacklist particles going to escape
  ///
  /// particles with
  ///  e > 1 AND
  ///  ( rho < maxFreeDensity OR noneigh <= 1 OR outside of > 90% tot mass) AND
  ///  t > 0
  /// are considered escapees and blacklisted
  ///
  /// particles with
  ///  e > 1 AND
  ///  ( 0.5 t_orbit > 2. t_togo ) AND
  ///  outside of > 90% tot mass AND
  ///  t > 0
  /// are considered disk particles which will not reimpact until
  /// the end of the simulation
  ///
  const fType gravConst = PartManager.attributes["gravconst"];
  const fType maxFreeDensity = PartManager.attributes["maxfreedens"];
  const fType maxOrbitTime = PartManager.attributes["maxorbittime"];
  const fType rMax = PartManager.attributes["maxradius"];

  const fType totMass = massBins(noMassBins - 1);
  fType curEscapeesMass = 0.;
  countsType noEscapees = 0;
  fType curDiskRemovedMass = 0.;
  countsType noDiskRemoved = 0;

  std::string domString = boost::lexical_cast<std::string>(myDomain);
  std::fstream escFile, rmvFile;
  
  std::string escFilename = "escDom000";
  escFilename.replace(escFilename.size() - 0 - domString.size(),
                      domString.size(), domString);
  escFile.open(escFilename.c_str(), std::ios::out | std::ios::app);
  
  std::string rmvFilename = "rmvDom000";
  rmvFilename.replace(rmvFilename.size() - 0 - domString.size(),
                      domString.size(), domString);
  rmvFile.open(rmvFilename.c_str(), std::ios::out | std::ios::app);

  for (size_t i = 0; i < noParts; i++)
    {
      ///
      /// calculate eccentricity vector and its magnitude
      ///
      const fType rx = (pos(i, X) - comX);
      const fType ry = (pos(i, Y) - comY);
      const fType rz = (pos(i, Z) - comZ);

      const fType r = sqrt(rx * rx + ry * ry + rz * rz);

      const fType vx = vel(i, X);
      const fType vy = vel(i, Y);
      const fType vz = vel(i, Z);

      const fType vivi = vx * vx + vy * vy + vz * vz;
      const fType rivi = rx * vx + ry * vy + rz * vz;

      const fType M = massProf(r);
      const fType mu = M * gravConst;

      const fType ex = (vivi * rx / mu) - (rivi * vx / mu) - (rx / r);
      const fType ey = (vivi * ry / mu) - (rivi * vy / mu) - (ry / r);
      const fType ez = (vivi * rz / mu) - (rivi * vz / mu) - (rz / r);

      const fType e = sqrt(ex * ex + ey * ey + ez * ez);
      ecc(i) = e;

      ///
      /// specific energy, semi-major axis and orbital period
      ///
      const fType specE = (vivi / 2.) - (mu / r);
      const fType a = mu / (2. * specE);
      const fType Torbit = 2.*M_PI*sqrt(a * a * a / mu);

      if (r > rMax && time > 0.)
        {
          if (e > 1. &&
              (rho(i) < maxFreeDensity || noneigh(i) <= 1 || M > 0.9 * totMass)
              )
            {
              PartManager.blacklisted[i] = true;
              curEscapeesMass += m(i);
              noEscapees++;

              escFile << time << "\t"
                      << id(i) << "\t"
                      << pos(i, X) << "\t"
                      << pos(i, Y) << "\t"
                      << pos(i, Z) << "\t"
                      << vel(i, X) << "\t"
                      << vel(i, Y) << "\t"
                      << vel(i, Z) << "\t"
                      << m(i) << "\t"
                      << h(i) << "\t"
                      << u(i) << "\t"
                      << p(i) << "\t"
                      << noneigh(i) << "\t"
                      << e << "\t"
                      << a << "\t"
                      << Torbit << "\t"
#ifdef SPHLATCH_TILLOTSON
                      << mat(i) << "\t"
#endif
#ifdef SPHLATCH_ANEOS
                      << mat(i) << "\t"
                      << phase(i) << "\t"
                      << T(i) << "\t"
#endif
                      << "\n";
            }

          if (e < 1. && e > 0.7 &&
              M > 0.9 * totMass &&
              Torbit > maxOrbitTime)
            {
              PartManager.blacklisted[i] = true;
              curDiskRemovedMass += m(i);
              noDiskRemoved++;

              rmvFile << time << "\t"
                      << id(i) << "\t"
                      << pos(i, X) << "\t"
                      << pos(i, Y) << "\t"
                      << pos(i, Z) << "\t"
                      << vel(i, X) << "\t"
                      << vel(i, Y) << "\t"
                      << vel(i, Z) << "\t"
                      << m(i) << "\t"
                      << h(i) << "\t"
                      << u(i) << "\t"
                      << p(i) << "\t"
                      << noneigh(i) << "\t"
                      << e << "\t"
                      << a << "\t"
                      << Torbit << "\t"
#ifdef SPHLATCH_TILLOTSON
                      << mat(i) << "\t"
#endif
#ifdef SPHLATCH_ANEOS
                      << mat(i) << "\t"
                      << phase(i) << "\t"
                      << T(i) << "\t"
#endif
                      << "\n";
            }
        }
    }
  escFile.close();
  rmvFile.close();

  CommManager.sum(curEscapeesMass);
  CommManager.sum(noEscapees);
  PartManager.attributes["escapedmass"] += curEscapeesMass;

  CommManager.sum(curDiskRemovedMass);
  CommManager.sum(noDiskRemoved);
  PartManager.attributes["diskremovemass"] += curDiskRemovedMass;

  Logger.stream << noEscapees << " escapees with mass of "
                << curEscapeesMass << "   and "
                << noDiskRemoved << " disk parts. with a mass of "
                << curDiskRemovedMass << " removed";
  Logger.flushStream();
#endif

  ///
  /// define the quantities to save in a dump
  ///
  quantsType saveQuants;
  saveQuants.vects += &pos, &vel, &acc;
  saveQuants.scalars += &m, &rho, &u, &p, &h, &cs;
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
#ifdef SPHLATCH_TILLOTSON
  saveQuants.ints += &mat;
#endif
#ifdef SPHLATCH_ANEOS
  saveQuants.scalars += &T;
  saveQuants.ints += &mat, &phase;
#endif
#ifdef SPHLATCH_INTEGRATERHO
  saveQuants.scalars += &drhodt;
#endif
#ifdef SPHLATCH_REMOVEESCAPING
  saveQuants.scalars += &ecc;
#endif
  const fType curSaveTime = (floor((time / saveItrvl) + 1.e-9))
                                * saveItrvl;
  if (fabs(curSaveTime - time) < 1.e-9)
    {
      std::string fileName = "dump";

      std::ostringstream stepSS;
      stepSS << step;

      ///
      /// pad step number to 7 numbers
      ///
      for (size_t i = stepSS.str().size(); i < 7; i++)
        {
          fileName.append("0");
        }
      fileName.append(stepSS.str());

      fileName.append("_T");
      std::ostringstream timeSS;
      timeSS << std::setprecision(4) << std::scientific << time;

      ///
      /// pad to a size of 11 characters
      /// (sign 1, mantissa 1, '.' 1, prec 3, 'e+' 2, exponent 3)
      ///
      for (size_t i = timeSS.str().size(); i < 11; i++)
        {
          fileName.append("0");
        }
      fileName.append(timeSS.str());
      fileName.append(".h5part");

      IOManager.saveDump(fileName, saveQuants);
      Logger.stream << "dump saved: " << fileName;
      Logger.flushStream();
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
  eos_type& EOS(eos_type::instance());

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  valvectRefType p(PartManager.p);
  valvectRefType u(PartManager.u);
  valvectRefType cs(PartManager.cs);
  valvectRefType rho(PartManager.rho);
  valvectRefType dudt(PartManager.dudt);
  valvectRefType dhdt(PartManager.dhdt);
  valvectRefType divv(PartManager.divv);

#ifdef SPHLATCH_GRAVITY
  valvectRefType eps(PartManager.eps);
#endif
#ifdef SPHLATCH_INTEGRATERHO
  valvectRefType drhodt(PartManager.drhodt);
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
#ifdef SPHLATCH_SHOCKTUBE
      ///
      /// shocktube boundary condition
      ///
      vel(i, X) = 0.;
      vel(i, Y) = 0.;

      if (pos(i, Z) < 5. || pos(i, Z) > 295.)
        vel(i, Z) = 0;
#endif
#ifdef SPHLATCH_GRAVITY
      eps(i) = 0.7 * h(i);
#endif
    }

#ifdef SPHLATCH_GRAVITY
  const fType gravTheta = PartManager.attributes["gravtheta"];
  const fType gravConst = PartManager.attributes["gravconst"];

  tree_type Tree(gravTheta, gravConst,
                 CostZone.getDepth(), CostZone.getCenter(),
                 CostZone.getSidelength());

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
  kernel_type Kernel;
  neighsearch_type Nsearch;

  ///
  /// the size of the neighbour list doesn't really
  /// affect the performance of the neighbour search,
  /// so it can be choosen quite large
  ///
  Nsearch.prepare();
  Nsearch.neighbourList.resize(16384);
  Nsearch.neighDistList.resize(16384);
  Logger << "Rankspace prepared";

#ifndef SPHLATCH_INTEGRATERHO
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
      Nsearch.findNeighbours(i, srchRad);

      const size_t noNeighs = Nsearch.neighbourList[0];

      static fType rhoi;
      rhoi = 0.;

      ///
      /// sum over neighbours
      ///
      for (size_t curNeigh = 1; curNeigh <= noNeighs; curNeigh++)
        {
          const fType r = Nsearch.neighDistList[curNeigh];
          const size_t j = Nsearch.neighbourList[curNeigh];

          const fType hij = 0.5 * (hi + h(j));

          rhoi += m(j) * Kernel.value(r, hij);
        }
      rho(i) = rhoi;
    }
  Logger << "SPH sum: rho";
#endif
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

  ///
  /// pressure
  /// I need    : u, rho
  /// I provide : p, cs
  ///
  for (size_t k = 0; k < noParts; k++)
    {
      EOS(k, p(k), cs(k));

      if (p(k) < 0.)
        p(k) = 0.;
    }
  Logger << "calculated pressure";
  CommManager.sendGhosts(p);
  CommManager.sendGhosts(cs);
  Logger << " sent to ghosts: p, cs";

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
  fType curDrhoDt = 0.;
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
      Nsearch.findNeighbours(i, srchRad);

      ///
      /// store the number of neighbours
      ///
      const size_t noNeighs = Nsearch.neighbourList[0];
      noneigh(i) = noNeighs;

      curAccX = 0.;
      curAccY = 0.;
      curAccZ = 0.;
#ifdef SPHLATCH_TIMEDEP_ENERGY
      curPow = 0.;
#endif
#ifdef SPHLATCH_VELDIV
      curDrhoDt = 0.;
#endif
      const fType viX = vel(i, X);
      const fType viY = vel(i, Y);
      const fType viZ = vel(i, Z);

      const fType riX = pos(i, X);
      const fType riY = pos(i, Y);
      const fType riZ = pos(i, Z);

      const fType ci = cs(i);

      ///
      /// sum over the neighbours
      ///
      for (size_t curNeigh = 1; curNeigh <= noNeighs; curNeigh++)
        {
          const fType rij = Nsearch.neighDistList[curNeigh];
          const size_t j = Nsearch.neighbourList[curNeigh];

          const fType rhoj = rho(j);
          const fType pj = p(j);

          const fType hij = 0.5 * (hi + h(j));

          const fType rijX = riX - pos(j, X);
          const fType rijY = riY - pos(j, Y);
          const fType rijZ = riZ - pos(j, Z);

          const fType vijX = viX - vel(j, X);
          const fType vijY = viY - vel(j, Y);
          const fType vijZ = viZ - vel(j, Z);

          const fType vijrij = rijX * vijX + rijY * vijY + rijZ * vijZ;

          /// make that a static or define outside of loop?
          fType av = 0;

          /// AV
          if (vijrij < 0.)
            {
              const fType rijrij = rijX * rijX + rijY * rijY + rijZ * rijZ;
              const fType rhoij = 0.5 * (rhoi + rhoj);
              const fType cij = 0.5 * (ci + cs(j));
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
          curDrhoDt += mjvijdivWij;
#endif
        }
      acc(i, X) += curAccX;
      acc(i, Y) += curAccY;
      acc(i, Z) += curAccZ;
#ifdef SPHLATCH_TIMEDEP_ENERGY
      dudt(i) = curPow;
#endif
#ifdef SPHLATCH_VELDIV
      divv(i) = curDrhoDt / rho(i);
      divvMax = divv(i) > divvMax ? divv(i) : divvMax;
#endif
#ifdef SPHLATCH_INTEGRATERHO
      drhodt(i) = -curDrhoDt;
#endif
#ifdef SPHLATCH_FRICTION
      const fType fricCoeff = 1. / PartManager.attributes["frictime"];
      acc(i, X) -= vel(i, X) * fricCoeff;
      acc(i, Y) -= vel(i, Y) * fricCoeff;
      acc(i, Z) -= vel(i, Z) * fricCoeff;
#endif
#ifdef SPHLATCH_SHOCKTUBE
      ///
      /// shocktube boundary condition
      ///
      if (pos(i, Z) < 5. || pos(i, Z) > 295.)
        acc(i, Z) = 0.;
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

      ///
      /// k1: too few neighbours   k2: number ok   k3: too many neighbours
      ///
      //const fType k1 = 0.5 * (1 + tanh((noNeighCur - noNeighMin) / -5.));
      const fType k1 = 0.; // inhibits h oscillation above hard surfaces
      const fType k3 = 0.5 * (1 + tanh((noNeighCur - noNeighMax) / 5.));
      const fType k2 = 1. - k1 - k3;

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
  ///
  /// parse program options
  ///
  po::options_description Options(
    "<input-file> <save-time> <stop-time>\n ... or use options");

  Options.add_options()
  ("input-file,i", po::value<std::string>(),
   "input file")
  ("save-time,s", po::value<fType>(),
   "save dumps when (time) modulo (save time) = 0.")
  ("stop-time,S", po::value<fType>(),
   "stop simulaton at this time");

  po::positional_options_description posDesc;
  posDesc.add("input-file", 1);
  posDesc.add("save-time", 1);
  posDesc.add("stop-time", 1);

  po::variables_map poMap;
  po::store(po::command_line_parser(argc, argv).options(Options).positional(posDesc).run(), poMap);
  po::notify(poMap);

  if (!poMap.count("input-file") ||
      !poMap.count("save-time") ||
      !poMap.count("stop-time"))
    {
      std::cerr << Options << "\n";
      return EXIT_FAILURE;
    }

  ///
  /// everythings set, now start the parallel envirnoment
  ///
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
  std::string loadDumpFile = poMap["input-file"].as<std::string>();
  const fType maxTime = poMap["stop-time"].as<fType>();
  IOManager.saveItrvl = poMap["save-time"].as<fType>();

  PartManager.attributes["gamma"] = 5. / 3.;
  PartManager.attributes["gravconst"] = 1.0;
  PartManager.attributes["gravtheta"] = 0.7;
  PartManager.attributes["courant"] = 0.3;

  PartManager.attributes["noneigh"] = 50.;
  PartManager.attributes["umin"] = 1000.;
#ifdef SPHLATCH_FRICTION
  PartManager.attributes["frictime"] = 200.;
#endif
#ifdef SPHLATCH_REMOVEESCAPING
  PartManager.attributes["mincoredens"] = 5.0;
  PartManager.attributes["maxfreedens"] = 0.1;
  PartManager.attributes["maxorbittime"] = 432000.;
  PartManager.attributes["maxradius"] = 2.e9;
  PartManager.attributes["escapedmass"] = 0.0;
#endif

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
#ifdef SPHLATCH_TILLOTSON
  PartManager.useMaterials();
#endif
#ifdef SPHLATCH_ANEOS
  PartManager.useMaterials();
  PartManager.usePhase();
  PartManager.useTemperature();
#endif
#ifdef SPHLATCH_INTEGRATERHO
  PartManager.useIntegratedRho();
#endif
#ifdef SPHLATCH_REMOVEESCAPING
  PartManager.useEccentricity();
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
#ifdef SPHLATCH_TILLOTSON
  idvectRefType mat(PartManager.mat);
#endif
#ifdef SPHLATCH_ANEOS
  idvectRefType mat(PartManager.mat);
  valvectRefType T(PartManager.T);
#endif
#ifdef SPHLATCH_INTEGRATERHO
  valvectRefType rho(PartManager.rho);
  valvectRefType drhodt(PartManager.drhodt);
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
#ifdef SPHLATCH_TILLOTSON
  CommManager.exchangeQuants.ints += &mat;
#endif
#ifdef SPHLATCH_ANEOS
  CommManager.exchangeQuants.ints += &mat;
  CommManager.exchangeQuants.scalars += &T;
#endif
#ifdef SPHLATCH_INTEGRATERHO
  CommManager.exchangeQuants.scalars += &rho;
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
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
  Integrator.regIntegration(h, dhdt);
#endif
#ifdef SPHLATCH_TIMEDEP_ENERGY
  Integrator.regIntegration(u, dudt);
#endif
#ifdef SPHLATCH_INTEGRATERHO
  Integrator.regIntegration(rho, drhodt);
#endif

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
#ifdef SPHLATCH_REMOVEESCAPING
                << "     removal of particles on escaping orbits\n"
#endif
                << "     basic SPH\n";
  Logger.flushStream();

  ///
  /// load particles
  ///
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


