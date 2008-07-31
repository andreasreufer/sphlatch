// some defs

#define SPHLATCH_PARALLEL

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
typedef sphlatch::valueType valueType;
typedef sphlatch::identType identType;
typedef sphlatch::matrixType matrixType;
typedef sphlatch::stringVectType stringVectType;
typedef sphlatch::valvectType valvectType;
typedef sphlatch::idvectType idvectType;

typedef sphlatch::valvectRefType valvectRefType;
typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::matrixRefType matrixRefType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

using namespace boost::assign;
using namespace sphlatch::vectindices;

///
/// quick&dirty helper functions
///
int readInt( std::fstream& _fin )
{
  union rawInt {
    char c[4];
    int i;
    } ri;
  for (size_t i = 0; i < 4; i++)
  {
    ri.c[i] = _fin.get();
  }
  return ri.i;
}

double readDbl( std::fstream& _fin )
{
  union rawDbl {
    char c[8];
    double d;
    } rd;
  for (size_t i = 0; i < 8; i++)
  {
    rd.c[i] = _fin.get();
  }
  return rd.d;
}

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);
  
  po::options_description Options("Global Options");
  Options.add_options()
  ("help,h", "Produces this Help")
  ("input-file,i" , po::value<std::string>(), "input  file")
  ("output-file,o", po::value<std::string>(), "output file");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  if (!VMap.count("output-file") && !VMap.count("input-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }
  
  part_type& PartManager(part_type::instance());
  io_type&        IOManager(io_type::instance());

  std::string outputFileName = VMap["output-file"].as<std::string>();
  std::string  inputFilename = VMap["input-file"].as<std::string>();
  
  idvectRefType id(PartManager.id);
  idvectRefType mat(PartManager.mat);

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);

  valvectRefType u(PartManager.u);
  valvectRefType h(PartManager.h);
  valvectRefType m(PartManager.m);
  valvectRefType rho(PartManager.rho);
  valvectRefType p(PartManager.p);
  
  sphlatch::quantsType saveQuants;
  saveQuants.vects += &pos, &vel;
  saveQuants.scalars += &u, &h, &m, &rho, &p;
  saveQuants.ints += &id, &mat;

  std::fstream fin;
  fin.open( inputFilename.c_str(), std::ios::in | std::ios::binary );

  ///
  /// read the "header"
  ///

  /// number of bytes in record
  readInt( fin );
  /// npart
  const size_t noParts = static_cast<size_t>( readInt( fin ) );
  /// n1
  readInt( fin );
  /// n2
  readInt( fin );
  /// t
  PartManager.attributes["time"] = static_cast<valueType>( readDbl( fin ) );
  /// trot
  readDbl( fin );
  /// tkin
  readDbl( fin );
  /// tgrav
  readDbl( fin );
  /// tterm
  readDbl( fin );

  PartManager.useBasicSPH();
  PartManager.useEnergy();
  PartManager.useMaterials();

  PartManager.setNoParts(noParts);
  PartManager.resizeAll();

  ///
  /// read the particles
  ///
  for (size_t j = 0; j < 14; j++)
  {
    for (size_t i = 0; i < noParts; i++)
      {
        switch(j)
        {
          case 0:
            id(i) = static_cast<identType>(i);
            pos(i, X) = static_cast<valueType>( readDbl( fin ) );
            break;
          case 1:
            pos(i, Y) = static_cast<valueType>( readDbl( fin ) );
            break;
          case 2:
            pos(i, Z) = static_cast<valueType>( readDbl( fin ) );
            break;
          case 3:
            vel(i, X) = static_cast<valueType>( readDbl( fin ) );
            break;
          case 4:
            vel(i, Y) = static_cast<valueType>( readDbl( fin ) );
            break;
          case 5:
            vel(i, Z) = static_cast<valueType>( readDbl( fin ) );
            break;
          case 6:
            u(i) = static_cast<valueType>( readDbl( fin ) );
            break;
          case 7:
            h(i) = static_cast<valueType>( readDbl( fin ) );
            break;
          case 8:
            m(i) = static_cast<valueType>( readDbl( fin ) );
            break;
          case 9:
            rho(i) = static_cast<valueType>( readDbl( fin ) );
            break;
          case 10:
            readDbl(fin); /// T
            break;
          case 11:
            p(i) = static_cast<valueType>( readDbl( fin ) );
            break;
          case 12:
            readDbl(fin); /// alpha?
            break;
          case 13:
            mat(i) = static_cast<identType>( lrint( readInt( fin ) ) );
            break;
        }
      } 
  }
  fin.close();

  ///
  /// conversion constants from the re,me-system to the cgs-system
  ///
  const valueType unitL = 6.348e8;   // earth radius in cm
  const valueType unitM = 5.3014e27; // earth mass in g
  const valueType unitT = 850.34;    // unit time in s

  // unit specific energy in erg/g
  const valueType unitU = unitL*unitL / (unitT*unitT);
  // unit pressure in barye
  const valueType unitP = unitM / ( unitL*unitT*unitT );
  // unit velocity in cm/s
  const valueType unitV = unitL / unitT;
  // unit density in g per cm^3
  const valueType unitRho = unitM / ( unitL*unitL*unitL );

  ///
  /// rescale the physical quantities into the cgs-system
  ///
  pos *= unitL;
  vel *= unitV;
  u   *= unitU;
  h   *= unitL;
  m   *= unitM;
  rho *= unitRho;
  p   *= unitP;
 
  ///
  /// the value for G mentioned in Michael Fuchs code in cgs
  /// (deviates from the NIST G)
  ///
  PartManager.attributes["gravconst"] = 6.6726e-8;;

  ///
  /// save everything
  ///
  IOManager.saveDump( outputFileName, saveQuants );

  MPI::Finalize();
  return EXIT_SUCCESS;
}


