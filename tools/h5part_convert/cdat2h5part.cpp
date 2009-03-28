// some defs

// uncomment for single-precision calculation
#define SPHLATCH_SINGLEPREC
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
typedef sphlatch::fType             fType;
typedef sphlatch::identType         identType;
typedef sphlatch::matrixType        matrixType;
typedef sphlatch::stringVectType    stringVectType;
typedef sphlatch::valvectType       valvectType;
typedef sphlatch::idvectType        idvectType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager   part_type;

#include "io_manager.h"
typedef sphlatch::IOManager         io_type;

#include "cdat_reader.h"
typedef sphlatch::CDATreader        cdatreader_type;

using namespace boost::assign;
using namespace sphlatch::vectindices;

///
/// quick&dirty helper functions
///
void readerToPartMgr(matrixType&      _mat,
                     size_t           _j,
                     cdatreader_type& _rdr,
                     size_t           _idx)
{
   const size_t noParts = _mat.size1();

   for (size_t i = 0; i < noParts; i++)
   {
      _mat(i, _j) = _rdr.read(i, _idx);
   }
}

void readerToPartMgr(valvectType&     _vct,
                     cdatreader_type& _rdr,
                     size_t           _idx)
{
   const size_t noParts = _vct.size();

   for (size_t i = 0; i < noParts; i++)
   {
      _vct(i) = _rdr.read(i, _idx);
   }
}

void readerToPartMgr(idvectType&      _ivct,
                     cdatreader_type& _rdr,
                     size_t           _idx)
{
   const size_t noParts = _ivct.size();

   for (size_t i = 0; i < noParts; i++)
   {
      _ivct(i) = static_cast<identType>(lrint(_rdr.read(i, _idx)));
   }
}

int main(int argc, char* argv[])
{
   MPI::Init(argc, argv);

   po::options_description Options("Global Options");
   Options.add_options()
   ("help,h", "Produces this Help")
   ("input-file,i", po::value<std::string>(), "input  file")
   ("output-file,o", po::value<std::string>(), "output file");

   po::variables_map VMap;
   po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
   po::notify(VMap);

   if (VMap.count("help"))
   {
      std::cerr << Options << std::endl;
      return(EXIT_FAILURE);
   }

   if (!VMap.count("output-file") && !VMap.count("input-file"))
   {
      std::cerr << Options << std::endl;
      return(EXIT_FAILURE);
   }

   part_type&      PartManager(part_type::instance());
   io_type&        IOManager(io_type::instance());
   cdatreader_type cdatReader;

   std::string outputFileName = VMap["output-file"].as<std::string>();
   std::string inputFilename  = VMap["input-file"].as<std::string>();

   cdatReader.open(inputFilename);

   const size_t noParts = cdatReader.getNoParts();
   PartManager.setNoParts(noParts, 0);
   PartManager.step = 0;

   PartManager.attributes = cdatReader.getAttrs();

   stringVectType vars = cdatReader.getVars();

   using namespace boost::assign;
   sphlatch::quantsType saveQuants;
   for (size_t curIdx = 0; curIdx < vars.size(); curIdx++)
   {
      if (vars[curIdx] == "id")
      {
         std::cout << vars[curIdx] << " ";
         saveQuants.ints += &PartManager.id;
      }

      if (vars[curIdx] == "mat")
      {
         std::cout << vars[curIdx] << " ";
         PartManager.useMaterials();
         saveQuants.ints += &PartManager.mat;
      }

      if ((vars[curIdx] == "x") || (vars[curIdx] == "y") ||
          (vars[curIdx] == "z"))
      {
         std::cout << vars[curIdx] << " ";
         saveQuants.vects += &PartManager.pos;
      }

      if ((vars[curIdx] == "vx") || (vars[curIdx] == "vy") ||
          (vars[curIdx] == "vz"))
      {
         std::cout << vars[curIdx] << " ";
         PartManager.useGravity();
         saveQuants.vects += &PartManager.vel;
      }

      if ((vars[curIdx] == "ax") || (vars[curIdx] == "ay") ||
          (vars[curIdx] == "az"))
      {
         std::cout << vars[curIdx] << " ";
         PartManager.useGravity();
         saveQuants.vects += &PartManager.acc;
      }

      if ((vars[curIdx] == "m") || (vars[curIdx] == "mass"))
      {
         std::cout << vars[curIdx] << " ";
         PartManager.useGravity();
         saveQuants.scalars += &PartManager.m;
      }

      if ((vars[curIdx] == "graveps") || (vars[curIdx] == "eps"))
      {
         std::cout << vars[curIdx] << " ";
         PartManager.useGravity();
         saveQuants.scalars += &PartManager.eps;
      }

      if (vars[curIdx] == "h")
      {
         std::cout << vars[curIdx] << " ";
         PartManager.useBasicSPH();
         saveQuants.scalars += &PartManager.h;
      }

      if ((vars[curIdx] == "p") || (vars[curIdx] == "P"))
      {
         std::cout << vars[curIdx] << " ";
         PartManager.useBasicSPH();
         saveQuants.scalars += &PartManager.p;
      }

      if (vars[curIdx] == "rho")
      {
         std::cout << vars[curIdx] << " ";
         PartManager.useBasicSPH();
         saveQuants.scalars += &PartManager.rho;
      }

      if (vars[curIdx] == "u")
      {
         std::cout << vars[curIdx] << " ";
         PartManager.useEnergy();
         saveQuants.scalars += &PartManager.u;
      }

      if (vars[curIdx] == "T")
      {
         std::cout << vars[curIdx] << " ";
         PartManager.useTemperature();
         saveQuants.scalars += &PartManager.T;
      }

      if ((vars[curIdx] == "dam") || (vars[curIdx] == "dm"))
      {
         std::cout << vars[curIdx] << " ";
         PartManager.useDamage();
         saveQuants.scalars += &PartManager.dam;
      }

      if ((vars[curIdx] == "sxx") ||
          (vars[curIdx] == "sxy") ||
          (vars[curIdx] == "sxz") ||
          (vars[curIdx] == "syy") ||
          (vars[curIdx] == "syz"))
      {
         std::cout << vars[curIdx] << " ";
         PartManager.useStress();
         saveQuants.vects += &PartManager.S;
      }

      if (vars[curIdx] == "peakp")
      {
         std::cout << vars[curIdx] << " ";
         PartManager.usePeakPress();
         saveQuants.scalars += &PartManager.peakp;
      }
   }
   std::cout << "\n";

   PartManager.resizeAll();

   for (size_t curIdx = 0; curIdx < vars.size(); curIdx++)
   {
      std::cout << vars[curIdx] << " " << std::flush;
      if (vars[curIdx] == "id")
         readerToPartMgr(PartManager.id, cdatReader, curIdx);
      if (vars[curIdx] == "mat")
         readerToPartMgr(PartManager.mat, cdatReader, curIdx);

      if (vars[curIdx] == "x")
         readerToPartMgr(PartManager.pos, X, cdatReader, curIdx);
      if (vars[curIdx] == "y")
         readerToPartMgr(PartManager.pos, Y, cdatReader, curIdx);
      if (vars[curIdx] == "z")
         readerToPartMgr(PartManager.pos, Z, cdatReader, curIdx);

      if (vars[curIdx] == "vx")
         readerToPartMgr(PartManager.vel, X, cdatReader, curIdx);
      if (vars[curIdx] == "vy")
         readerToPartMgr(PartManager.vel, Y, cdatReader, curIdx);
      if (vars[curIdx] == "vz")
         readerToPartMgr(PartManager.vel, Z, cdatReader, curIdx);

      if (vars[curIdx] == "ax")
         readerToPartMgr(PartManager.acc, X, cdatReader, curIdx);
      if (vars[curIdx] == "ay")
         readerToPartMgr(PartManager.acc, Y, cdatReader, curIdx);
      if (vars[curIdx] == "az")
         readerToPartMgr(PartManager.acc, Z, cdatReader, curIdx);

      if ((vars[curIdx] == "m") || (vars[curIdx] == "mass"))
         readerToPartMgr(PartManager.m, cdatReader, curIdx);
      if (vars[curIdx] == "h")
         readerToPartMgr(PartManager.h, cdatReader, curIdx);
      if ((vars[curIdx] == "p") || (vars[curIdx] == "P"))
         readerToPartMgr(PartManager.p, cdatReader, curIdx);
      if (vars[curIdx] == "rho")
         readerToPartMgr(PartManager.rho, cdatReader, curIdx);
      if (vars[curIdx] == "u")
         readerToPartMgr(PartManager.u, cdatReader, curIdx);
      if (vars[curIdx] == "T")
         readerToPartMgr(PartManager.T, cdatReader, curIdx);
      if ((vars[curIdx] == "graveps") || (vars[curIdx] == "eps"))
         readerToPartMgr(PartManager.eps, cdatReader, curIdx);
      if ((vars[curIdx] == "dam") || (vars[curIdx] == "dm"))
         readerToPartMgr(PartManager.dam, cdatReader, curIdx);

      if ((vars[curIdx] == "mat"))
         readerToPartMgr(PartManager.dam, cdatReader, curIdx);

      if (vars[curIdx] == "sxx")
         readerToPartMgr(PartManager.S, XX, cdatReader, curIdx);
      if (vars[curIdx] == "sxy")
         readerToPartMgr(PartManager.S, XY, cdatReader, curIdx);
      if (vars[curIdx] == "sxz")
         readerToPartMgr(PartManager.S, XZ, cdatReader, curIdx);
      if (vars[curIdx] == "syy")
         readerToPartMgr(PartManager.S, YY, cdatReader, curIdx);
      if (vars[curIdx] == "syz")
         readerToPartMgr(PartManager.S, YZ, cdatReader, curIdx);

      if (vars[curIdx] == "peakp")
         readerToPartMgr(PartManager.peakp, cdatReader, curIdx);
   }
   std::cout << " -> " << outputFileName << "\n";
   cdatReader.close();

   IOManager.setSinglePrecOut();
   IOManager.saveDump(outputFileName, saveQuants);

   MPI::Finalize();
   return(EXIT_SUCCESS);
}
