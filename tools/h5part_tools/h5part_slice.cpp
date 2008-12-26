#define SPHLATCH_PARALLEL
#define SPHLATCH_CARTESIAN_XYZ

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

#include <boost/assign/std/vector.hpp>

namespace po = boost::program_options;

#include "typedefs.h"
typedef sphlatch::fType             fType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager   part_type;

#include "io_manager.h"
typedef sphlatch::IOManager         io_type;

#include "costzone.h"
typedef sphlatch::CostZone costzone_type;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager comm_type;

using namespace boost::assign;
using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
   MPI::Init(argc, argv);
   po::options_description Options("Global Options");

   const fType fMax = std::numeric_limits<fType>::max();

   Options.add_options() ("help,h", "Produces this Help")
   ("input-file,i", po::value<std::string>(), "input file")
   ("output-file,o", po::value<std::string>(), "output file")
   ("xmin", po::value<fType>()->default_value(-fMax),
    "minimal value for x")
   ("xmax", po::value<fType>()->default_value(fMax),
    "maximal value for x")
   ("ymin", po::value<fType>()->default_value(-fMax),
    "minimal value for y")
   ("ymax", po::value<fType>()->default_value(fMax),
    "maximal value for y")
   ("zmin", po::value<fType>()->default_value(-fMax),
    "minimal value for z")
   ("zmax", po::value<fType>()->default_value(fMax),
    "maximal value for z");

   po::variables_map VMap;
   po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
   po::notify(VMap);

   if (!VMap.count("output-file") || VMap.count("help") ||
       !VMap.count("input-file"))
   {
      std::cerr << Options << std::endl;
      return(EXIT_FAILURE);
   }

   part_type& PartManager(part_type::instance());
   io_type&   IOManager(io_type::instance());
   comm_type& CommManager(comm_type::instance());
   costzone_type& CostZone(costzone_type::instance());

   std::string inputFileName  = VMap["input-file"].as<std::string>();
   std::string outputFileName = VMap["output-file"].as<std::string>();

   const fType xmin = VMap["xmin"].as<fType>();
   const fType xmax = VMap["xmax"].as<fType>();
   const fType ymin = VMap["ymin"].as<fType>();
   const fType ymax = VMap["ymax"].as<fType>();
   const fType zmin = VMap["zmin"].as<fType>();
   const fType zmax = VMap["zmax"].as<fType>();

   using namespace sphlatch;
   using namespace boost::assign;

   quantsType saveQuants = IOManager.getQuants(inputFileName);
   PartManager.setUsedQuants(saveQuants);
   CommManager.exchangeQuants = saveQuants;

   ///
   /// load particles
   /// 
   IOManager.loadDump(inputFileName);

   ///
   /// filter particles
   ///
   matrixRefType pos(PartManager.pos);

   const size_t noParts = PartManager.getNoLocalParts();
   countsType noTotRawParts = noParts;
   CommManager.sum(noTotRawParts);


   for (size_t i = 0; i < noParts; i++)
   {
      PartManager.blacklisted[i] =
         (pos(i, X) < xmin || pos(i, X) > xmax ||
          pos(i, Y) < ymin || pos(i, Y) > ymax ||
          pos(i, Z) < zmin || pos(i, Z) > zmax);
   }

   ///
   /// exchange and save
   ///
   CostZone.createDomainPartsIndex();
   CommManager.exchange(CostZone.domainPartsIndex, 0);
   
   countsType noTotCutParts = PartManager.getNoLocalParts();
   CommManager.sum(noTotRawParts);

   std::cout << "keeping " << noTotCutParts << " of "
             << noTotRawParts << " particles\n";

   IOManager.saveDump(outputFileName, saveQuants);

   MPI::Finalize();
   return(EXIT_SUCCESS);
}
