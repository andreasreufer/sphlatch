//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

// enable parallel version
#define SPHLATCH_PARALLEL

//#define SPHLATCH_RANKSPACESERIALIZE

// enable checking of bounds for the neighbour lists
//#define SPHLATCH_CHECKNONEIGHBOURS

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assign/std/vector.hpp>

namespace po = boost::program_options;

#include "typedefs.h"
typedef sphlatch::valueType valueType;
typedef sphlatch::identType identType;
typedef sphlatch::valvectType valvectType;

typedef sphlatch::valvectRefType valvectRefType;
typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::matrixRefType matrixRefType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager comm_type;

#include "costzone.h"
typedef sphlatch::CostZone costzone_type;

#include "kernel_cubicspline3d.h"
typedef sphlatch::CubicSpline3D kernel_type;

#include "rankspace.h"
typedef sphlatch::Rankspace neighsearch_type;

#include "eos_tillotson.h"
typedef sphlatch::Tillotson eos_type;

#include "lagrange_sphere1D_solver.h"
typedef sphlatch::LagrangeSphere1DSolver lg1D_solver_type;

#include "lookup_table.h"
typedef sphlatch::LookupTable<sphlatch::InterpolateStepwise> lut_type;

using namespace boost::assign;
using namespace sphlatch::vectindices;

void densPress()
{
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  costzone_type& CostZone(costzone_type::instance());

  matrixRefType pos(PartManager.pos);

  valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  valvectRefType u(PartManager.u);
  valvectRefType p(PartManager.p);
  valvectRefType cs(PartManager.cs);
  valvectRefType rho(PartManager.rho);

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  const size_t noParts = PartManager.getNoLocalParts();
  const size_t noTotParts = noParts + PartManager.getNoGhostParts();

  /// send ghosts to other domains
  CommManager.sendGhosts(pos);
  CommManager.sendGhosts(id);
  CommManager.sendGhosts(m);
  CommManager.sendGhosts(h);

  ///
  /// define kernel and neighbour search algorithm
  ///
  kernel_type Kernel;
  neighsearch_type Nsearch;
  Nsearch.prepare();
  Nsearch.neighbourList.resize(1024);
  Nsearch.neighDistList.resize(1024);

  ///
  /// 1st SPH sum: density
  /// I need    : pos, h, m
  /// I provide : rho
  ///
/*  for (size_t k = 0; k < noParts; k++)
    {
      ///
      /// find neighbours
      ///
      const size_t i = k;
      const valueType hi = h(i);
      const valueType srchRad = 2. * hi;
      Nsearch.findNeighbours(i, srchRad);

      const size_t noNeighs = Nsearch.neighbourList[0];

      ///
      /// store the number of neighbours
      ///
      noneigh(i) = noNeighs;

      static valueType rhoi;
      rhoi = 0.;

      ///
      /// sum over neighbours
      ///
      for (size_t curNeigh = 1; curNeigh <= noNeighs; curNeigh++)
        {
          const valueType r = Nsearch.neighDistList[curNeigh];
          const size_t j = Nsearch.neighbourList[curNeigh];

          const valueType hij = 0.5 * (hi + h(j));

          rhoi += m(j) * Kernel.value(r, hij);
        }
      rho(i) = rhoi;
    }*/


  ///
  /// calculate pressure
  ///
  eos_type& EOS(eos_type::instance());
  for (size_t k = 0; k < noParts; k++)
    {
      EOS(k, p(k), cs(k));
    }
}


int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  po::options_description Options("Global Options");
  Options.add_options()
  ("help,h", "Produces this Help")
  ("output-file,o", po::value<std::string>(), "output file")
  ("input-file,i", po::value<std::string>(), "input  file");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  if (!VMap.count("output-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  io_type&        IOManager(io_type::instance());
  part_type&      PartManager(part_type::instance());
  comm_type&      CommManager(comm_type::instance());
  costzone_type&  CostZone(costzone_type::instance());

  const valueType mTot = 5.3014e27;
  const valueType mCore = 0.30 * mTot;
  const valueType mMantle = 0.70 * mTot;

  const identType matCore = 5;
  const identType matMantle = 4;

  ///
  /// guess for the density
  ///
  const valueType rhoCore = 6.;
  const valueType rhoMantle = 3.;

  const valueType rCoreGuess = pow(0.75 * mCore / (M_PI * rhoCore), 1. / 3.);
  const valueType rMantleGuess = pow(0.75 * mMantle / (M_PI * rhoMantle)
                                + pow(rCoreGuess, 3.), 1. / 3.);
  const size_t noCells = 200;
  const size_t noEdges = noCells + 1;

  const size_t noCoreCells = lrint(noCells * (rCoreGuess / rMantleGuess));
  const size_t noMantleCells = noCells - noCoreCells;

  lg1D_solver_type Solver(noCells);

  const valueType drCore = rCoreGuess / noCoreCells;
  const valueType drMantle = (rMantleGuess - rCoreGuess) / noMantleCells;

  ///
  /// set the shell edges
  ///
  Solver.r(0) = 0.;
  Solver.v(0) = 0.;
  for (size_t i = 1; i < noEdges; i++)
    {
      if (i <= noCoreCells)
        Solver.r(i) = Solver.r(i - 1) + drCore;
      else
        Solver.r(i) = Solver.r(i - 1) + drMantle;

      Solver.v(i) = 0.;
    }

  ///
  /// set the shell centered values
  ///
  for (size_t i = 0; i < noCells; i++)
    {
      const valueType rI = Solver.r(i);
      const valueType rO = Solver.r(i + 1);
      const valueType rC = 0.5 * (rI + rO);

      const valueType vol = (4. * M_PI / 3.) * (rO * rO * rO - rI * rI * rI);

      if (rC < rCoreGuess)
        {
          Solver.m(i) = vol * rhoCore;
          Solver.mat(i) = 5;
        }
      else
        {
          Solver.m(i) = vol * rhoMantle;
          Solver.mat(i) = 4;
        }

      Solver.u(i) = 1.e10;
    }

  ///
  /// set some constants
  ///
  Solver.uMin = 1.e4;
  Solver.friction = 1.e-2;

  ///
  /// integrate for a certain physical time
  ///
  std::cerr << " start 1D Lagrange solver\n";
  Solver.integrateTo(5.e3); /// replace by option

  valvectType rCen;
  rCen.resize(noCells);

  for (size_t i = 0; i < noCells; i++)
    {
      rCen(i) = 0.5 * (Solver.r(i) + Solver.r(i + 1));
    }

  std::string dumpfilename = "dump.hdf5";
  IOManager.savePrimitive(Solver.r, "r", dumpfilename);
  IOManager.savePrimitive(rCen, "rCen", dumpfilename);
  IOManager.savePrimitive(Solver.v, "v", dumpfilename);
  IOManager.savePrimitive(Solver.m, "m", dumpfilename);
  IOManager.savePrimitive(Solver.u, "u", dumpfilename);
  IOManager.savePrimitive(Solver.p, "p", dumpfilename);
  IOManager.savePrimitive(Solver.mat, "mat", dumpfilename);
  IOManager.savePrimitive(Solver.rho, "rho", dumpfilename);

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  valvectRefType p(PartManager.p);
  valvectRefType u(PartManager.u);
  valvectRefType rho(PartManager.rho);

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);
  idvectRefType mat(PartManager.mat);

  PartManager.useBasicSPH();
  PartManager.useEnergy();
  PartManager.useTimedepEnergy();
  PartManager.useMaterials();

  ///
  /// load particles
  ///
  std::string inputFilename = VMap["input-file"].as<std::string>();
  IOManager.loadDump(inputFilename);
  std::cerr << " read " << inputFilename << "\n";

  ///
  /// exchange particle data
  ///
  CostZone.createDomainPartsIndex();
  CommManager.exchange(CostZone.domainPartsIndex,
                       CostZone.getNoGhosts());
  CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());

  ///
  /// scale position, specific volume and smoothing length
  ///
  const valueType rScale = Solver.r(noCells);
  valueType rCore = 0.;
  for (size_t i = 1; i < noCells; i++)
  {
    if ( Solver.mat(i-1) != Solver.mat(i) )
      rCore = Solver.r(i);
  }
  const valueType vScale = pow(rScale, 3.);

  const valueType rhoSet = 3.;

  const size_t noParts = PartManager.getNoLocalParts();
  for (size_t k = 0; k < noParts; k++)
    {
      pos(k, X) *= rScale;
      pos(k, Y) *= rScale;
      pos(k, Z) *= rScale;

      h(k) *= rScale;
      m(k) *= vScale;
    }


  ///
  /// set mass by looking up the density and multiplying
  /// it with the particle volume stored in the mass variable
  ///
  lut_type rhoLUT(rCen, Solver.rho);
  lut_type uLUT(rCen, Solver.u);

  for (size_t k = 0; k < noParts; k++)
    {
      const valueType curRad = sqrt(pos(k, X) * pos(k, X) +
                                    pos(k, Y) * pos(k, Y) +
                                    pos(k, Z) * pos(k, Z));

      const valueType curRho = rhoLUT(curRad);
      m(k) *= curRho;
      rho(k) = curRho;
      u(k) = uLUT(curRad);

      if (curRad < rCore)
        mat(k) = matCore;
      else
        mat(k) = matMantle;
    }

  std::cerr << " initial guess for particle masses, now calc. SPH density ...\n";
  densPress();

  ///
  /// correct, check, do again
  ///

  sphlatch::quantsType saveQuants;
  saveQuants.vects += &pos;
  saveQuants.scalars += &m, &rho, &h, &u, &p;
  saveQuants.ints += &id, &noneigh, &mat;

  PartManager.step = 0;
  PartManager.attributes["time"] = 0.;
  PartManager.attributes["gravconst"] = Solver.gravConst;
  PartManager.attributes["umin"] = Solver.uMin;

  std::string outputFilename = VMap["output-file"].as<std::string>();
  std::cerr << " -> " << outputFilename << "\n";
  IOManager.saveDump(outputFilename, saveQuants);
  std::cerr << "particles saved ... \n";

  MPI::Finalize();
  return EXIT_SUCCESS;
}

