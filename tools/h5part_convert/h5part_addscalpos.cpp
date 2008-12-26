// some defs

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#define SPHLATCH_TIMEDEP_ENERGY
#define SPHLATCH_TIMEDEP_SMOOTHING
#define SPHLATCH_GRAVITY
#define SPHLATCH_ANEOS
#define SPHLATCH_REMOVEESCAPING

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
typedef sphlatch::fType fType;
typedef sphlatch::valvectType valvectType;

typedef sphlatch::valvectRefType valvectRefType;
typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::matrixRefType matrixRefType;
typedef sphlatch::quantsType quantsType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

using namespace boost::assign;

int main(int argc, char* argv[])
{
  po::options_description Options("Global Options");
  Options.add_options()
  ("help,h", "Produces this Help")
  ("input-file,i" , po::value<std::string>(), "input  file");
  //("output-file,o" , po::value<std::string>(), "output file");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help") ||
      not VMap.count("input-file") )
      //not VMap.count("output-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  
  part_type& PartManager(part_type::instance());
  io_type&        IOManager(io_type::instance());

  std::string  inputFilename = VMap["input-file"].as<std::string>();
  //std::string  outputFilename = VMap["output-file"].as<std::string>();
  std::string  outputFilename = inputFilename;

  std::cout << "adding scalar position and velocity to " << inputFilename << "\n";

  PartManager.useBasicSPH();
/*  PartManager.useEnergy();
  
#ifdef SPHLATCH_TIMEDEP_ENERGY
  PartManager.useTimedepEnergy();
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
  PartManager.useTimedepH();
#endif
#ifdef SPHLATCH_GRAVITY
  PartManager.useGravity();
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
  */
  IOManager.loadDump( inputFilename );

  ///
  /// paraview doesnt like groupnames
  /// other than "Step#0"
  ///
  PartManager.step = 0;
  
  matrixRefType pos( PartManager.pos );
  matrixRefType vel( PartManager.vel );

  //IOManager.setSinglePrecOut();
  //IOManager.setDoublePrecOut();

  valvectType tmp;
  tmp.resize( pos.size1() );

  quantsType saveScalVects;
  saveScalVects.scalars += &tmp;

  PartManager.regQuantity( tmp, "pos_x" );
  tmp = sphlatch::quantColumnType( pos, 0 );
  IOManager.saveDump( outputFilename, saveScalVects );
  PartManager.unRegQuantity( tmp );
  
  PartManager.regQuantity( tmp, "pos_y" );
  tmp = sphlatch::quantColumnType( pos, 1 );
  IOManager.saveDump( outputFilename, saveScalVects );
  PartManager.unRegQuantity( tmp );

  PartManager.regQuantity( tmp, "pos_z" );
  tmp = sphlatch::quantColumnType( pos, 2 );
  IOManager.saveDump( outputFilename, saveScalVects );
  PartManager.unRegQuantity( tmp );
  
  PartManager.regQuantity( tmp, "vel_x" );
  tmp = sphlatch::quantColumnType( vel, 0 );
  IOManager.saveDump( outputFilename, saveScalVects );
  PartManager.unRegQuantity( tmp );
  
  PartManager.regQuantity( tmp, "vel_y" );
  tmp = sphlatch::quantColumnType( vel, 1 );
  IOManager.saveDump( outputFilename, saveScalVects );
  PartManager.unRegQuantity( tmp );

  PartManager.regQuantity( tmp, "vel_z" );
  tmp = sphlatch::quantColumnType( vel, 2 );
  IOManager.saveDump( outputFilename, saveScalVects );
  PartManager.unRegQuantity( tmp );
 
  ///
  /// save the rest
  ///
  /*quantsType saveQuants;
  valvectRefType h( PartManager.h );
  valvectRefType m( PartManager.m );
  valvectRefType u( PartManager.u );
  valvectRefType rho( PartManager.rho );
  saveQuants.scalars += &h, &m, &rho, &u;

  idvectRefType id( PartManager.id );
  saveQuants.ints += &id;

#ifdef SPHLATCH_TIMEDEP_ENERGY
  valvectRefType dudt( PartManager.dudt );
  saveQuants.scalars += &dudt;
#endif
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
  valvectRefType dhdt( PartManager.dhdt );
  valvectRefType divv( PartManager.divv );
  saveQuants.scalars += &dhdt, &divv;
#endif
#ifdef SPHLATCH_GRAVITY
  valvectRefType eps( PartManager.eps );
  saveQuants.scalars += &eps;
#endif
#ifdef SPHLATCH_TILLOTSON
  idvectRefType mat( PartManager.mat );
  saveQuants.ints += &mat;
#endif
#ifdef SPHLATCH_ANEOS
  valvectRefType T( PartManager.T );
  saveQuants.scalars += &T;
  
  idvectRefType mat( PartManager.mat );
  idvectRefType phase( PartManager.phase );
  saveQuants.ints += &mat, &phase;
#endif
#ifdef SPHLATCH_INTEGRATERHO
  valvectRefType drhodt( PartManager.drhodt );
  saveQuants.scalars += &drhodt;
#endif
#ifdef SPHLATCH_REMOVEESCAPING
  valvectRefType ecc( PartManager.ecc );
  saveQuants.scalars += &ecc;
#endif
  IOManager.saveDump( outputFilename, saveQuants );
  
  ///
  /// delete the /current link
  ///

  hid_t filePropList = H5Pcreate(H5P_FILE_ACCESS);
  //H5Pset_fapl_mpio(filePropList, MPI::COMM_WORLD, MPI::INFO_NULL);
  hid_t fileHandle = H5Fopen(outputFilename.c_str(),
                             H5F_ACC_RDWR, filePropList);
  H5Pclose(filePropList);

  hid_t rootGroup = H5Gopen(fileHandle, "/", H5P_DEFAULT);
  H5Ldelete(rootGroup, "/current", H5P_DEFAULT);
  
  H5Gclose(rootGroup);
  H5Fclose(fileHandle);*/

  return EXIT_SUCCESS;
}

