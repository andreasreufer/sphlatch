#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
//#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

// enable parallel version
#define SPHLATCH_PARALLEL

//#define SPHLATCH_RANKSPACESERIALIZE

// enable checking of bounds for the neighbour lists
#define SPHLATCH_CHECKNONEIGHBOURS

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

#include <boost/progress.hpp>

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

#include "rankspace.h"
typedef sphlatch::Rankspace neighsearch_type;

#include "eos_tillotson.h"
typedef sphlatch::Tillotson till_type;


using namespace boost::assign;
using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  po::options_description Options("Global Options");
  Options.add_options()
  ("help,h", "Produces this Help")
  ("output-file,o", po::value<std::string>(), "output  file")
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
      !VMap.count("input-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  io_type&        IOManager(io_type::instance());
  part_type&      PartManager(part_type::instance());
  comm_type&      CommManager(comm_type::instance());
  costzone_type&  CostZone(costzone_type::instance());
  till_type&      Tillotson(till_type::instance());

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType S(PartManager.S);

  valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  valvectRefType rho(PartManager.rho);
  valvectRefType u(PartManager.u);

  valvectRefType dam(PartManager.dam);
  valvectRefType epsmin(PartManager.epsmin);
  valvectRefType acoef(PartManager.acoef);
  valvectRefType mweib(PartManager.mweib);
  valvectRefType young(PartManager.young);

  idvectRefType id(PartManager.id);
  idvectRefType mat(PartManager.mat);
  idvectRefType noflaws(PartManager.noflaws);

  PartManager.useBasicSPH();
  PartManager.useEnergy();
  PartManager.useMaterials();
  PartManager.useDamage();
  PartManager.useStress();

  ///
  /// load particles
  ///
  std::string inputFilename = VMap["input-file"].as<std::string>();
  IOManager.loadDump(inputFilename);
  std::cerr << " read particle positions from: " << inputFilename << "\n";

  ///
  /// scale position, specific volume and smoothing length
  ///
  const valueType rScale = 20;
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

  /// dunite: 4
  identType surfId = 4;
  till_type::paramType surfParams = Tillotson.getMatParams(surfId);
  const valueType surfU = 1.72e9;

  ///
  /// set particle mass by multiplyin the particle volume
  /// with the desired density
  ///
  for (size_t k = 0; k < noParts; k++)
    {
      ///
      /// the flawed surface
      ///
      const valueType surfRho = Tillotson.findRho(surfU, surfId,
                                                  0., 1.e-2, 0.1, 10.);
      m(k) *= surfRho;
      rho(k) = surfRho;

      u(k) = surfU;
      mat(k) = surfId;

      acoef(k) = 0.4 * sqrt((surfParams.A + (4. / 3.) * surfParams.xmu) /
                            surfParams.rho0) / (2. * h(k));
      young(k) = 9. * surfParams.A * surfParams.xmu
                 / (3. * surfParams.A + surfParams.xmu);

      ///
      /// no initial stress
      ///
      S(k, XX) = 0.;
      S(k, XY) = 0.;
      S(k, XZ) = 0.;
      S(k, YY) = 0.;
      S(k, YZ) = 0.;
    }

  ///
  /// assign flaws to particles
  ///
  const size_t noTotFlaws = std::max(10 * noParts, 1000lu);

  std::cerr << "assigning " << noTotFlaws << " flaws to " << noParts
            << " particles\n";

  ///
  /// factor 300 from Benz setup-collision.f
  ///
  const valueType weibKV = 300. / (surfParams.cweib * volScale);
  const valueType weibExp = 1. / surfParams.pweib;

  boost::progress_display flawsProgress(noTotFlaws);
  size_t noSetFlaws = 0;
  while (noSetFlaws < noTotFlaws)
    {
      const size_t i = lrint(noParts * (rand()
                                        / static_cast<double>(RAND_MAX)));

      const valueType curEps =
        pow(weibKV * (static_cast<valueType>(noSetFlaws) + 1), weibExp);

      ///
      /// the first flaw is always the weakest,
      /// as noSetFlaws rises monotonously
      ///
      if (noflaws(i) == 0)
        epsmin(i) = curEps;

      noflaws(i)++;
      ++flawsProgress;

      /// just for temporary storage
      mweib(i) = curEps;

      noSetFlaws++;
    }

  ///
  /// set Weibull parameters
  ///
  for (size_t k = 0; k < noParts; k++)
    {
      dam(k) = 0.;

      ///
      /// no flaw was assigned, so give it one flaw with maximal strength
      ///
      if (noflaws(k) < 1)
        {
          noflaws(k) = 1;
          epsmin(k) =
            pow(weibKV * static_cast<valueType>(noTotFlaws + 1), weibExp);
        }

      if (noflaws(k) == 1)
        mweib(k) = 1.;
      else
        mweib(k) = log(static_cast<valueType>(noflaws(k)))
                   / log(mweib(k) / epsmin(k));
    }

  sphlatch::quantsType saveQuants;
  saveQuants.vects += &pos, &vel, &S;
  saveQuants.scalars += &m, &h, &rho, &u, &dam, &epsmin, &acoef, &mweib, &young;
  saveQuants.ints += &id, &mat, &noflaws;

  PartManager.step = 0;
  PartManager.attributes["time"] = 0.;

  std::string outputFilename = VMap["output-file"].as<std::string>();
  std::cerr << " -> " << outputFilename << "\n";
  IOManager.saveDump(outputFilename, saveQuants);
  std::cerr << "particles saved ... \n";

  MPI::Finalize();
  return EXIT_SUCCESS;
}

