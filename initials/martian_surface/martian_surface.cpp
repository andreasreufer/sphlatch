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

#include <boost/progress.hpp>

namespace po = boost::program_options;

#include "typedefs.h"
typedef sphlatch::fType                  fType;
typedef sphlatch::identType              identType;
typedef sphlatch::valvectType            valvectType;

typedef sphlatch::valvectRefType         valvectRefType;
typedef sphlatch::idvectRefType          idvectRefType;
typedef sphlatch::matrixRefType          matrixRefType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager        part_type;

#include "io_manager.h"
typedef sphlatch::IOManager              io_type;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager   comm_type;

#include "costzone.h"
typedef sphlatch::CostZone               costzone_type;

#include "rankspace.h"
typedef sphlatch::Rankspace              neighsearch_type;

#include "eos_tillotson.h"
typedef sphlatch::Tillotson              till_type;

#ifdef SPHLATCH_ANEOS
 #include "eos_aneos.h"
typedef sphlatch::ANEOS                  eos_type;
#else
typedef till_type                        eos_type;
#endif

using namespace boost::assign;
using namespace sphlatch::vectindices;

struct layerT
{
   identType            mat;
   fType                top, bottom;

   // acoef has to be divided by 2h for use in the Benz&Asphaug fracture model
   fType                rho, U, acoef, young, dam;

   size_t               noParts;
   size_t               noFlaws;

   till_type::paramType matParams;

   fType                weibKV, weibExp;
};

int main(int argc, char* argv[])
{
   MPI::Init(argc, argv);

   po::options_description Options("Global Options");
   Options.add_options()
   ("help,h", "Produces this Help")
   ("output-file,o", po::value<std::string>(), "output  file")
   ("input-file,i", po::value<std::string>(), "input   file")
   ("layA-thick,a", po::value<fType>(), "surface dunite layer thickness")
   ("layB-thick,b", po::value<fType>(), "buried  ice    layer thickness")
   //("layC-thick,c", po::value<fType>(), "buried  dunite layer thickness")
   ("dun-dam,d", po::value<fType>(), "initial dunite damage")
   ("scale,s", po::value<fType>(), "scaling constant    ");

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
       !VMap.count("layA-thick") ||
       !VMap.count("layB-thick") ||
       //!VMap.count("layC-thick") ||
       !VMap.count("scale") ||
       !VMap.count("dun-dam"))
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
   CommManager.exchangeQuants.vects   += &pos, &vel;
   CommManager.exchangeQuants.scalars += &m, &h;
   CommManager.exchangeQuants.ints    += &id;

   CostZone.createDomainPartsIndex();
   CommManager.exchange(CostZone.domainPartsIndex,
                        CostZone.getNoGhosts());
   CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());

   ///
   /// scale position, specific volume and smoothing length
   ///
   const fType rScale   = VMap["scale"].as<fType>();
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

   ///
   /// define materials
   ///

#ifdef SPHLATCH_ANEOS
   /// ice: 2
   identType iceId = 2;
#else
   /// ice: 17
   identType iceId = 17;
#endif
   till_type::paramType iceParams = Tillotson.getMatParams(17);

   /// dunite: 4
   identType            dunId     = 4;
   till_type::paramType dunParams = Tillotson.getMatParams(dunId);


#ifdef SPHLATCH_ANEOS
   const fType KtoEV = 1.1604505e4;
   fType       dummy;

   const fType iceRhoGuess = 1.10;
   const fType iceT        = 240 / KtoEV;
   fType       iceU        = 0.;
   EOS.getSpecEnergy(iceRhoGuess, iceT, iceId, dummy, dummy, iceU);

   const fType dunRhoGuess = 3.32;
   const fType dunT        = 240 / KtoEV;
   fType       dunU        = 0.;
   EOS.getSpecEnergy(dunRhoGuess, dunT, dunId, dummy, dummy, dunU);
#else
   const fType iceU = 5.04e9;
   const fType dunU = 1.72e9;
#endif
   const fType iceDamage = 0.;
   const fType dunDamage = VMap["dun-dam"].as<fType>();

   std::vector<layerT> layers(3);

   std::cerr << "get equilibrium densities for ice ...";
   const fType iceRho   = EOS.findRho(iceU, iceId, 1.e2, 1.e1, 0.50, 10.0);
   const fType iceAcoef = 0.4 * sqrt(
      (iceParams.A + (4. / 3.) * iceParams.xmu) /
      iceParams.rho0);
   const fType iceYoung = 9. * iceParams.A * iceParams.xmu /
                          (3. * iceParams.A + iceParams.xmu);

   std::cerr << " done!\n";
   std::cerr << "get equilibrium densities for dunite ...";

   const fType dunRho   = EOS.findRho(dunU, dunId, 1.e2, 1.e1, 0.10, 10.0);
   const fType dunAcoef = 0.4 * sqrt(
      (dunParams.A + (4. / 3.) * dunParams.xmu) /
      dunParams.rho0);
   const fType dunYoung = 9. * dunParams.A * dunParams.xmu /
                          (3. * dunParams.A + dunParams.xmu);
   std::cerr << " done!\n";

   layers[0].mat       = dunId;
   layers[0].matParams = dunParams;
   layers[0].top       = 0.;
   layers[0].bottom    = -VMap["layA-thick"].as<fType>();
   layers[0].rho       = dunRho;
   layers[0].U         = dunU;
   layers[0].acoef     = dunAcoef;
   layers[0].young     = dunYoung;
   layers[0].dam       = dunDamage;
   layers[0].noParts   = 0;

   layers[1].mat       = iceId;
   layers[1].matParams = iceParams;
   layers[1].top       = layers[0].bottom;
   layers[1].bottom    = layers[1].top - VMap["layB-thick"].as<fType>();
   layers[1].rho       = iceRho;
   layers[1].U         = iceU;
   layers[1].acoef     = iceAcoef;
   layers[1].young     = iceYoung;
   layers[1].dam       = iceDamage;
   layers[1].noParts   = 0;

   layers[2].mat       = dunId;
   layers[2].matParams = dunParams;
   layers[2].top       = layers[1].bottom;
   layers[2].bottom    = -std::numeric_limits<fType>::infinity();
   layers[2].rho       = dunRho;
   layers[2].U         = dunU;
   layers[2].acoef     = dunAcoef;
   layers[2].young     = dunYoung;
   layers[2].dam       = dunDamage;
   layers[2].noParts   = 0;


   ///
   /// set particle mass by multiplying the particle volume
   /// with the desired density
   ///
   std::vector<size_t> layerID(noParts);
   for (size_t k = 0; k < noParts; k++)
   {
      ///
      /// no initial stress
      ///
      S(k, XX) = 0.;
      S(k, XY) = 0.;
      S(k, XZ) = 0.;
      S(k, YY) = 0.;
      S(k, YZ) = 0.;

      size_t layIdx = 0;
      layerID[k] = 0;
      for (size_t i = 0; i < 3; i++)
      {
         if ((pos(k, Z) < layers[i].top) &&
             (pos(k, Z) > layers[i].bottom))
         {
            layIdx     = i;
            layerID[k] = i;
            break;
         }
      }

      m(k)  *= layers[layIdx].rho;
      rho(k) = layers[layIdx].rho;

      u(k)   = layers[layIdx].U;
      mat(k) = layers[layIdx].mat;
      dam(k) = layers[layIdx].dam;

      acoef(k) = layers[layIdx].acoef / (2 * h(k));
      young(k) = layers[layIdx].young;

      layers[layIdx].noParts++;
   }

   for (size_t i = 0; i < 3; i++)
   {
      if (layers[i].noParts > 0)
      {
         layers[i].noFlaws = std::max(10 * layers[i].noParts, 10000000lu);

         std::cerr << "assigning " << layers[i].noFlaws << " flaws to layer " << i
                   << " with " << layers[i].noParts << " particles\n";

         const fType typVol = pow(fabs(layers[i].top - layers[i].bottom), 3);

         ///
         /// set surface flaws
         /// factor 300 from Benz setup-collision.f
         ///
         const fType weibKV  = 300. / (layers[i].matParams.cweib * typVol);
         const fType weibExp = 1. / layers[i].matParams.pweib;

         layers[i].weibKV  = weibKV;
         layers[i].weibExp = weibExp;

         boost::progress_display flawsProgress(layers[i].noFlaws);
         size_t noSetFlaws = 0;
         while (noSetFlaws < layers[i].noFlaws)
         {
            const size_t j =
               lrint(noParts * (rand() / static_cast<double>(RAND_MAX)));

            if (layerID[j] == i)
            {
               const fType curEps = pow(
                  weibKV *
                  (static_cast<fType>(noSetFlaws) + 1.),
                  weibExp);

               ///
               /// the first flaw is always the weakest,
               /// as noSetFlaws rises monotonously
               ///
               if (noflaws(j) == 0)
                  epsmin(j) = curEps;

               noflaws(j)++;
               ++flawsProgress;
               noSetFlaws++;

               ///
               /// just for temporary storage
               ///
               mweib(j) = curEps;
            }
         }
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
         const size_t curLay = layerID[k];

         epsmin(k) = pow(layers[curLay].weibKV * static_cast<fType>(
                            layers[curLay].noFlaws + 1), layers[curLay].weibExp);
      }

      if (noflaws(k) == 1)
         mweib(k) = 1.;
      else
         mweib(k) =
            log(static_cast<fType>(noflaws(k)))
            / log(mweib(k) / epsmin(k));

      ///
      /// no flaws
      ///
      //noflaws(k) = 0;
      //mweib(k) = 1.;
      //epsmin(k) = 1.e20; /// arbitrarly high strength, from Benz code
   }

   sphlatch::quantsType saveQuants;
   saveQuants.vects   += &pos, &vel, &S;
   saveQuants.scalars += &m, &h, &rho, &u,
   &dam, &epsmin, &acoef, &mweib,
   &young;
   saveQuants.ints += &id, &mat, &noflaws;

   PartManager.step = 0;
   PartManager.attributes["time"] = 0.;

   std::string outputFilename =
      VMap["output-file"].as<std::string>();
   std::cerr << " -> " << outputFilename <<
   "\n";
   IOManager.saveDump(outputFilename,
                      saveQuants);
   std::cerr << "particles saved ... \n";

   MPI::Finalize();
   return(EXIT_SUCCESS);
}
