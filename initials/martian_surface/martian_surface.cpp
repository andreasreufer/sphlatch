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
typedef sphlatch::fType fType;
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

#ifdef SPHLATCH_ANEOS
 #include "eos_aneos.h"
typedef sphlatch::ANEOS eos_type;
#else
typedef till_type eos_type;
#endif

using namespace boost::assign;
using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  po::options_description Options("Global Options");
  Options.add_options()
  ("help,h", "Produces this Help")
  ("output-file,o", po::value<std::string>(), "output  file")
  ("input-file,i", po::value<std::string>(), "input   file")
  ("surf-thick,t", po::value<fType>(), "rel. surf. thickness") 
  ("scale,s", po::value<fType>(),      "scaling constant    ");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help"))
    {
      std::cerr << Options << std::endl;
      return(EXIT_FAILURE);
    }

  if (!VMap.count("output-file") ||
      !VMap.count("input-file") ||
      !VMap.count("scale") ||
      !VMap.count("surf-thick"))
    {
      std::cerr << Options << std::endl;
      return(EXIT_FAILURE);
    }

  io_type&       IOManager(io_type::instance());
  part_type&     PartManager(part_type::instance());
  comm_type&     CommManager(comm_type::instance());
  costzone_type& CostZone(costzone_type::instance());
  till_type&     Tillotson(till_type::instance());
  eos_type&      EOS(eos_type::instance());

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
  /// keep only the lower half-sphere
  ///
  size_t noParts = PartManager.getNoLocalParts();
  for (size_t k = 0; k < noParts; k++)
    {
      if (pos(k, Z) > 0.)
        PartManager.blacklisted[k] = true;
    }

  ///
  /// exchange particle data
  ///
  CommManager.exchangeQuants.vects += &pos, &vel;
  CommManager.exchangeQuants.scalars += &m, &h;
  CommManager.exchangeQuants.ints += &id;

  CostZone.createDomainPartsIndex();
  CommManager.exchange(CostZone.domainPartsIndex,
                       CostZone.getNoGhosts());
  CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());

  ///
  /// scale position, specific volume and smoothing length
  ///
  const fType rScale = VMap["scale"].as<fType>();
  const fType volScale = pow(rScale, 3.);

  noParts = PartManager.getNoLocalParts();
  std::cerr << "keeping " << noParts << " particles \n";
  for (size_t k = 0; k < noParts; k++)
    {
      pos(k, X) *= rScale;
      pos(k, Y) *= rScale;
      pos(k, Z) *= rScale;

      h(k) *= rScale;

      m(k) *= volScale;
    }

  //const fType surfRelThickness = 0.075;
  fType surfRelThickness = VMap["surf-thick"].as<fType>();
  const fType surfThickness = surfRelThickness * rScale;

#ifdef SPHLATCH_ANEOS
  /// ice: 2
  identType baseId = 2;
#else
  /// ice: 17
  identType baseId = 17;
#endif
  till_type::paramType baseParams = Tillotson.getMatParams(17);
  const fType baseU = 5.04e9;

  /// dunite: 4
  identType surfId = 4;
  till_type::paramType surfParams = Tillotson.getMatParams(surfId);
  const fType surfU = 1.72e9;

  ///
  /// set particle mass by multiplyin the particle volume
  /// with the desired density
  ///

  const fType surfRho = EOS.findRho(surfU, surfId, 0., 1.e5, 3.00, 3.50);
  std::cout << "surface   energy: " << surfU << "   equ. density: " << surfRho << "\n";

  const fType baseRho = EOS.findRho(baseU, baseId, 0., 1.e5, 0.95, 1.20);
  std::cout << "base      energy: " << baseU << "   equ. density: " << baseRho << "\n";

  size_t noSurfParts = 0, noBaseParts = 0;
  for (size_t k = 0; k < noParts; k++)
    {
      if (pos(k, Z) < -surfThickness)
        {
          ///
          /// the unflawed base
          ///
          m(k) *= baseRho;
          rho(k) = baseRho;

          u(k) = baseU;
          mat(k) = baseId;

          acoef(k) = 0.4 * sqrt((baseParams.A + (4. / 3.) * baseParams.xmu) /
                                baseParams.rho0) / (2. * h(k));
          young(k) = 9. * baseParams.A * baseParams.xmu
                     / (3. * baseParams.A + baseParams.xmu);

          dam(k) = 0.0;
          noBaseParts++;
        }
      else
        {
          ///
          /// the flawed surface
          ///
          m(k) *= surfRho;
          rho(k) = surfRho;

          u(k) = surfU;
          mat(k) = surfId;

          acoef(k) = 0.4 * sqrt((surfParams.A + (4. / 3.) * surfParams.xmu) /
                                surfParams.rho0) / (2. * h(k));
          young(k) = 9. * surfParams.A * surfParams.xmu
                     / (3. * surfParams.A + surfParams.xmu);

          dam(k) = 0.9;
          noSurfParts++;
        }

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
  const size_t noTotSurfFlaws = std::max(10 * noSurfParts, 10000000lu);

  std::cerr << "assigning " << noTotSurfFlaws << " flaws to " << noSurfParts
            << " surface particles\n";

  ///
  /// the typical surface volume
  ///
  const fType surfVol = surfThickness * surfThickness * surfThickness;

  ///
  /// set surface flaws
  /// factor 300 from Benz setup-collision.f
  ///
  const fType surfWeibKV = 300. / (surfParams.cweib * surfVol);
  const fType surfWeibExp = 1. / surfParams.pweib;

  boost::progress_display flawsSurfProgress(noTotSurfFlaws);
  size_t noSetFlaws = 0;
  while (noSetFlaws < noTotSurfFlaws)
    {
      const size_t i = lrint(noParts * (rand()
                                        / static_cast<double>(RAND_MAX)));

      if (mat(i) == surfId)
        {
          const fType curEps =
            pow(surfWeibKV * (static_cast<fType>(noSetFlaws) + 1), surfWeibExp);

          ///
          /// the first flaw is always the weakest,
          /// as noSetFlaws rises monotonously
          ///
          if (noflaws(i) == 0)
            epsmin(i) = curEps;

          noflaws(i)++;
          ++flawsSurfProgress;

          /// just for temporary storage
          mweib(i) = curEps;

          noSetFlaws++;
        }
    }

  ///
  /// assign flaws to particles
  ///
  const size_t noTotBaseFlaws = std::max(10 * noBaseParts, 10000000lu);

  std::cerr << "assigning " << noTotBaseFlaws << " flaws to " << noBaseParts
            << " base    particles\n";

  ///
  /// the typical base volume
  ///
  const fType baseVol = pow(rScale - surfThickness, 3.);

  ///
  /// set base flaws
  /// factor 300 from Benz setup-collision.f
  ///
  const fType baseWeibKV = 300. / (baseParams.cweib * baseVol);
  const fType baseWeibExp = 1. / baseParams.pweib;

  boost::progress_display flawsBaseProgress(noTotBaseFlaws);
  noSetFlaws = 0;
  while (noSetFlaws < noTotBaseFlaws)
    {
      const size_t i = lrint(noParts * (rand()
                                        / static_cast<double>(RAND_MAX)));

      if (mat(i) == baseId)
        {
          const fType curEps =
            pow(baseWeibKV * (static_cast<fType>(noSetFlaws) + 1), baseWeibExp);

          ///
          /// the first flaw is always the weakest,
          /// as noSetFlaws rises monotonously
          ///
          if (noflaws(i) == 0)
            epsmin(i) = curEps;

          noflaws(i)++;
          ++flawsBaseProgress;

          /// just for temporary storage
          mweib(i) = curEps;

          noSetFlaws++;
        }
    }

  ///
  /// set Weibull parameters
  ///
  for (size_t k = 0; k < noParts; k++)
    {
      ///
      /// no flaw was assigned, so give it one flaw with maximal strength
      ///
      if (noflaws(k) < 1)
        {
          noflaws(k) = 1;
          if (mat(k) == surfId)
            epsmin(k) =
              pow(surfWeibKV * static_cast<fType>(noTotSurfFlaws + 1),
                  surfWeibExp);
          else
            epsmin(k) =
              pow(baseWeibKV * static_cast<fType>(noTotBaseFlaws + 1),
                  baseWeibExp);
        }

      if (noflaws(k) == 1)
        mweib(k) = 1.;
      else
        mweib(k) = log(static_cast<fType>(noflaws(k)))
        / log(mweib(k) / epsmin(k));

      ///
      /// no flaws
      ///
      //noflaws(k) = 0;
      //mweib(k) = 1.;
      //epsmin(k) = 1.e20; /// arbitrarly high strength, from Benz code
    }

  sphlatch::quantsType saveQuants;
  saveQuants.vects += &pos, &vel, &S;
  saveQuants.scalars += &m, &h, &rho, &u, &dam, &epsmin, &acoef, &mweib,
  &young;
  saveQuants.ints += &id, &mat, &noflaws;

  PartManager.step = 0;
  PartManager.attributes["time"] = 0.;

  std::string outputFilename = VMap["output-file"].as<std::string>();
  std::cerr << " -> " << outputFilename << "\n";
  IOManager.saveDump(outputFilename, saveQuants);
  std::cerr << "particles saved ... \n";

  MPI::Finalize();
  return(EXIT_SUCCESS);
}
