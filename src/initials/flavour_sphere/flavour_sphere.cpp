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

void density()
{
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  //costzone_type& CostZone(costzone_type::instance());

  matrixRefType pos(PartManager.pos);

  valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  valvectRefType rho(PartManager.rho);

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  const size_t noParts = PartManager.getNoLocalParts();
  //const size_t noTotParts = noParts + PartManager.getNoGhostParts();

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
  for (size_t k = 0; k < noParts; k++)
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
    }


}

void pressure()
{
  part_type& PartManager(part_type::instance());
  eos_type& EOS(eos_type::instance());

  valvectRefType p(PartManager.p);
  valvectRefType cs(PartManager.cs);

  const size_t noParts = PartManager.getNoLocalParts();
  ///
  /// calculate pressure
  ///
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
  ("output-file,o", po::value<std::string>(), "output  file")
  ("profile-file,p", po::value<std::string>(), "profile file")
  ("input-file,i", po::value<std::string>(), "input   file");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  if (!VMap.count("output-file") ||
      !VMap.count("profile-file") ||
      !VMap.count("input-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  io_type&        IOManager(io_type::instance());
  part_type&      PartManager(part_type::instance());
  comm_type&      CommManager(comm_type::instance());
  costzone_type&  CostZone(costzone_type::instance());

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
  std::cerr << " read particle positions from: " << inputFilename << "\n";

  ///
  /// exchange particle data
  ///
  CostZone.createDomainPartsIndex();
  CommManager.exchange(CostZone.domainPartsIndex,
                       CostZone.getNoGhosts());
  CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());

  ///
  /// load profile
  ///
  std::string profileFilename = VMap["profile-file"].as<std::string>();
  valvectType rProf;
  IOManager.loadPrimitive(rProf, "r", profileFilename);
  valvectType uProf;
  IOManager.loadPrimitive(uProf, "u", profileFilename);
  valvectType rhoProf;
  IOManager.loadPrimitive(rhoProf, "rho", profileFilename);
  valvectType matProf;
  IOManager.loadPrimitive(matProf, "mat", profileFilename);
  std::cerr << " read profile from:            " << profileFilename << "\n";

  ///
  /// scale position, specific volume and smoothing length
  ///
  const size_t noCells = rProf.size();
  const valueType rScale = rProf(noCells - 1);
  const valueType volScale = pow(rScale, 3.);

  const size_t noParts = PartManager.getNoLocalParts();
  for (size_t k = 0; k < noParts; k++)
    {
      pos(k, X) *= rScale;
      pos(k, Y) *= rScale;
      pos(k, Z) *= rScale;

      h(k) *= rScale;
      m(k) *= volScale;
    }


  ///
  /// set mass by looking up the density and multiplying
  /// it with the particle volume stored in the mass variable
  ///
  lut_type rhoLUT(rProf, rhoProf);
  lut_type uLUT(rProf, uProf);
  lut_type matLUT(rProf, matProf);

  for (size_t k = 0; k < noParts; k++)
    {
      const valueType curRad = sqrt(pos(k, X) * pos(k, X) +
                                    pos(k, Y) * pos(k, Y) +
                                    pos(k, Z) * pos(k, Z));

      const valueType curRho = rhoLUT(curRad);
      m(k) *= curRho;
      rho(k) = curRho;
      u(k) = uLUT(curRad);
      mat(k) = lrint( matLUT(curRad) );
    }
  
  
  sphlatch::quantsType saveQuants;

  std::cerr << " SPH density ...\n";
  density();
  saveQuants.ints += &noneigh;

  std::cerr << " calulate pressures ...\n";
  pressure();

  saveQuants.vects += &pos, &vel;
  saveQuants.scalars += &m, &rho, &h, &u, &p;
  saveQuants.ints += &id, &mat;

  PartManager.step = 0;
  PartManager.attributes["time"] = 0.;
  PartManager.attributes["gravconst"] =
    IOManager.loadAttribute("gravconst", profileFilename);
  PartManager.attributes["umin"] =
    IOManager.loadAttribute("umin", profileFilename);

  std::string outputFilename = VMap["output-file"].as<std::string>();
  std::cerr << " -> " << outputFilename << "\n";
  IOManager.saveDump(outputFilename, saveQuants);
  std::cerr << "particles saved ... \n";

  MPI::Finalize();
  return EXIT_SUCCESS;
}

