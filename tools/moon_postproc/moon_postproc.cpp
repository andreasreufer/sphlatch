// some defs

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

// enable parallel version
#define SPHLATCH_PARALLEL

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

#include <gsl/gsl_histogram.h>

#include "typedefs.h"
typedef sphlatch::valueType fType;
typedef sphlatch::valueRefType fRefType;

typedef sphlatch::valvectType valvectType;
typedef sphlatch::valvectRefType valvectRefType;
typedef sphlatch::zerovalvectType zerovalvectType;

typedef sphlatch::matrixType matrixType;

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

#include "kernel_cubicspline.h"
typedef sphlatch::CubicSpline3D kernel_type;

#include "lookup_table1D.h"
typedef sphlatch::LookupTable1D<sphlatch::InterpolateLinear> lut_type;

using namespace sphlatch::vectindices;
using namespace boost::assign;

int main(int argc, char* argv[])
{
  ///
  /// parse program options
  ///
  po::options_description Options(
    //"<input-file> <save-time> <stop-time>\n ... or use options");
    "<input-file>\n ... or use options");

  Options.add_options()
  ("input-file,i", po::value<std::string>(),
   "input file")
  ("output-file,o", po::value<std::string>(),
   "output file")
  ("escapee-file,e", po::value<std::string>(),
   "escapee file")
  ("max-rad,r", po::value<fType>(),
   "maximal binning radius (default: auto)");

  po::positional_options_description posDesc;
  posDesc.add("input-file", 1);
  posDesc.add("output-file", 1);
  posDesc.add("escapee-file", 1);
  posDesc.add("max-rad", 1);

  po::variables_map poMap;
  po::store(po::command_line_parser(argc, argv).options(Options).positional(posDesc).run(), poMap);
  po::notify(poMap);

  if (!poMap.count("input-file") ||
      !poMap.count("output-file"))
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

  ///
  /// parse program options
  ///
  std::string loadDumpFile = poMap["input-file"].as<std::string>();
  std::string saveFile = poMap["output-file"].as<std::string>();

  std::string escapeeFile;

  if (poMap.count("escapee-file"))
    {
      escapeeFile = poMap["escapee-file"].as<std::string>();
    }
  else
    {
      escapeeFile = "escapees";
    }

  ///
  /// define what we're doing
  ///
  PartManager.useGravity();
  PartManager.useBasicSPH();
  PartManager.useEnergy();

  PartManager.useMaterials();
  PartManager.usePhase();
  PartManager.useTemperature();

  PartManager.useEccentricity();

  ///
  /// some useful references
  ///
  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);
  valvectRefType u(PartManager.u);
  valvectRefType h(PartManager.h);

  valvectRefType T(PartManager.T);
  valvectRefType rho(PartManager.rho);
  valvectRefType ecc(PartManager.ecc);

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);
  idvectRefType mat(PartManager.mat);

  fRefType dumpTime(PartManager.attributes["time"]);

  const size_t myDomain = CommManager.getMyDomain();

  ///
  /// register the quantites to be exchanged
  ///
  CommManager.exchangeQuants.vects += &pos, &vel;
  CommManager.exchangeQuants.scalars += &m, &u, &h, &T;
  CommManager.exchangeQuants.ints += &id, &mat;

  ///
  /// load particles
  ///
  IOManager.loadDump(loadDumpFile);


  ///
  /// exchange particle data
  ///
  CostZone.createDomainPartsIndex();
  CommManager.exchange(CostZone.domainPartsIndex,
                       CostZone.getNoGhosts());

  ///
  /// prepare ghost sends
  ///
  //CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());

  const size_t noParts = PartManager.getNoLocalParts();
  const size_t noTotParts = noParts + PartManager.getNoGhostParts();
  int& step(PartManager.step);
  int& substep(PartManager.substep);

  ///
  /// estimate center of protoearth by calculating
  /// the center of mass of all particles with higher than
  /// minimal core density
  ///
  const fType minCoreDensity = PartManager.attributes["mincoredens"];
  fType comX = 0., comY = 0., comZ = 0., comM = 0.;
  fType comVX = 0., comVY = 0., comVZ = 0.;
  for (size_t i = 0; i < noParts; i++)
    {
      if (rho(i) > minCoreDensity)
        {
          comX += pos(i, X) * m(i);
          comY += pos(i, Y) * m(i);
          comZ += pos(i, Z) * m(i);

          comVX += vel(i, X) * m(i);
          comVY += vel(i, Y) * m(i);
          comVZ += vel(i, Z) * m(i);

          comM += m(i);
        }
    }
  CommManager.sum(comM);

  CommManager.sum(comX);
  CommManager.sum(comY);
  CommManager.sum(comZ);
  comX /= comM;
  comY /= comM;
  comZ /= comM;

  CommManager.sum(comVX);
  CommManager.sum(comVY);
  CommManager.sum(comVZ);
  comVX /= comM;
  comVY /= comM;
  comVZ /= comM;

  std::cerr << "planet center:  pos ["
            << comX << ","
            << comY << ","
            << comZ << "]\n"
            << "                vel ["
            << comVX << ","
            << comVY << ","
            << comVZ << "]\n";
  ///
  /// determine maximal radius relative to center of mass
  ///
  fType maxRad = 0.;
  const size_t noMassBins = 10000;

  if (poMap.count("max-rad"))
    {
      maxRad = poMap["max-rad"].as<fType>();
    }
  else
    {
      for (size_t i = 0; i < noParts; i++)
        {
          const fType curRad = sqrt((pos(i, X) - comX) * (pos(i, X) - comX) +
                                    (pos(i, Y) - comY) * (pos(i, Y) - comY) +
                                    (pos(i, Z) - comZ) * (pos(i, Z) - comZ));
          maxRad = curRad > maxRad ? curRad : maxRad;
        }
      CommManager.max(maxRad);
    }

  ///
  /// bin particle masses according to radius, cumulate mass bins
  /// and sum them up globally
  ///
  valvectType massRad(noMassBins), totMassBins(noMassBins);
  for (size_t i = 0; i < noMassBins; i++)
    {
      massRad(i) = maxRad * (i / static_cast<fType>(noMassBins));
      totMassBins(i) = 0.;
    }
  
  ///
  /// store radii
  ///
  IOManager.savePrimitive(massRad, "r", saveFile);

  for (size_t i = 0; i < noParts; i++)
    {
      const fType curRad = sqrt((pos(i, X) - comX) * (pos(i, X) - comX) +
                                (pos(i, Y) - comY) * (pos(i, Y) - comY) +
                                (pos(i, Z) - comZ) * (pos(i, Z) - comZ));

      const size_t curBin = std::max(std::min(
                                       lrint((curRad / maxRad) * noMassBins
                                             - 0.5),
                                       static_cast<long int>(noMassBins - 1)),
                                     static_cast<long int>(0));

      totMassBins(curBin) += m(i);
    }
  
  ///
  /// cumulative bins
  ///
  for (size_t i = 1; i < noMassBins; i++)
    {
      totMassBins(i) += totMassBins(i - 1);
    }
  CommManager.sum(totMassBins);

  ///
  /// create lookup table for the mass profile
  ///
  lut_type massProf(massRad, totMassBins);

  ///
  /// some constants
  ///
  const fType Rearth = 6.378e8; // [cm]
  const fType Rroche = 2.9 * Rearth;

  const fType Mearth = 5.974e27; // [g]
  const fType Mmoon = 7.345e25; // [g]

  ///
  /// determine orbital elements
  ///
  const fType gravConst = PartManager.attributes["gravconst"];
  const fType maxFreeDensity = PartManager.attributes["maxfreedens"];
  const fType maxOrbitTime = PartManager.attributes["maxorbittime"];
  const fType rMax = PartManager.attributes["maxradius"];

  valvectType MdiskIron(1), MdiskDun(1), Mpe(1), Mesc(1);
  valvectType Ldisk(3), Lpe(3), Lrem(3), Lesc(3), Prem(3), Pesc(3);
  zerovalvectType zero(3), zeroBins(noMassBins);

  Lesc = zero;
  Lrem = zero;
  
  Pesc = zero;
  Prem = zero;
  
  Ldisk = zero;
  Lpe = zero;

  MdiskIron(0) = 0.;
  MdiskDun(0) = 0.;
  Mpe(0) = 0.;
  Mesc(0) = 0.;
  
  size_t noDiskIronParts = 0, noDiskDunParts = 0,
         noPlanParts = 0, noEscParts = 0;

  matrixType LBins(noMassBins, 3);

  valvectType MProfIron(noMassBins);
  MProfIron = zeroBins;
  valvectType MProfDun(noMassBins);
  MProfDun = zeroBins;

  for (size_t i = 0; i < noParts; i++)
    {
      ///
      /// calculate eccentricity vector and its magnitude
      ///
      const fType x = pos(i, X);
      const fType y = pos(i, Y);
      const fType z = pos(i, Z);

      const fType vx = vel(i, X);
      const fType vy = vel(i, Y);
      const fType vz = vel(i, Z);

      const fType rx = x - comX;
      const fType ry = y - comY;
      const fType rz = z - comZ;

      const fType rvx = vx - comVX;
      const fType rvy = vy - comVY;
      const fType rvz = vz - comVZ;

      const fType r = sqrt(rx * rx + ry * ry + rz * rz);
      const fType vivi = rvx * rvx + rvy * rvy + rvz * rvz;
      const fType rivi = rx * vx + ry * vy + rz * vz;

      const fType M = massProf(r);
      const fType mu = M * gravConst;

      const fType ex = (vivi * rx / mu) - (rivi * rvx / mu) - (rx / r);
      const fType ey = (vivi * ry / mu) - (rivi * rvy / mu) - (ry / r);
      const fType ez = (vivi * rz / mu) - (rivi * rvz / mu) - (rz / r);

      const fType e = sqrt(ex * ex + ey * ey + ez * ez);

      const fType rlx = m(i) * (ry * rvz - rz * rvy);
      const fType rly = m(i) * (rz * rvx - rx * rvz);
      const fType rlz = m(i) * (rx * rvy - ry * rvx);

      const fType lx = m(i) * (y * vz - z * vy);
      const fType ly = m(i) * (z * vx - x * vz);
      const fType lz = m(i) * (x * vy - y * vx);

      const fType px = m(i) * vx;
      const fType py = m(i) * vy;
      const fType pz = m(i) * vz;

      ///
      /// orbital energy, semi-major axis and orbital period
      ///
      const fType specE = (vivi / 2.) - (mu / r);
      const fType a = mu / (2. * specE);
      const fType Torbit = 2.*M_PI*sqrt(a * a * a / mu);

      const size_t curBin = std::max(std::min(
                                       lrint((r / maxRad) * noMassBins - 0.5),
                                       static_cast<long int>(noMassBins - 1)),
                                     static_cast<long int>(0));

      ///
      /// particle properties
      ///
      const bool isEscaping = (e > 1.) && (r > Rroche);
      const bool isDisk = (e < 1.) && (r > Rroche);
      const bool isPlanet = not (isEscaping or isDisk);
      const bool isIron = ( mat(i) == 5 );

      if (not isEscaping)
        {
          Lrem(X) += lx;
          Lrem(Y) += ly;
          Lrem(Z) += lz;
          
          Prem(X) += px;
          Prem(Y) += py;
          Prem(Z) += pz;

          LBins(curBin, X) += rlx;
          LBins(curBin, Y) += rly;
          LBins(curBin, Z) += rlz;
          
          if ( isIron )
            MProfIron(curBin) += m(i);
          else
            MProfDun(curBin) += m(i);
        }
      else
        {
          Lesc(X) += lx;
          Lesc(Y) += ly;
          Lesc(Z) += lz;

          Pesc(X) += px;
          Pesc(Y) += py;
          Pesc(Z) += pz;
        }

      if (isEscaping)
        {
          Mesc(0) += m(i);
          noEscParts++;
        }

      if (isDisk)
        {
          Ldisk(X) += rlx;
          Ldisk(Y) += rly;
          Ldisk(Z) += rlz;
  
          if ( isIron )
          {
            MdiskIron(0) += m(i);
            noDiskIronParts++;
          }
          else
          {
            MdiskDun(0) += m(i);
            noDiskDunParts++;
          }
        }

      if (isPlanet)
        {
          Lpe(X) += rlx;
          Lpe(Y) += rly;
          Lpe(Z) += rlz;
          
          Mpe(0) += m(i);
          noPlanParts++;
        }
    }

  CommManager.sum(Mpe);
  IOManager.savePrimitive(Mpe, "Mpe", saveFile);

  CommManager.sum(Lpe);
  IOManager.savePrimitive(Lpe, "Lpe", saveFile);
  
  CommManager.sum(Ldisk);
  IOManager.savePrimitive(Ldisk, "Ldisk", saveFile);

  CommManager.sum(MdiskDun);
  IOManager.savePrimitive(MdiskDun, "MdiskDun", saveFile);
  
  CommManager.sum(MdiskIron);
  IOManager.savePrimitive(MdiskIron, "MdiskIron", saveFile);
 
  valvectType Mdisk = MdiskIron + MdiskDun;
  IOManager.savePrimitive(Mdisk, "Mdisk", saveFile);

  CommManager.sum(LBins);
  IOManager.savePrimitive(LBins, "L", saveFile);

  CommManager.sum(Lrem);
  IOManager.savePrimitive(Lrem, "Lrem", saveFile);
  
  CommManager.sum(Prem);
  IOManager.savePrimitive(Prem, "Prem", saveFile);

  ///
  /// now add escapees
  ///
  std::fstream escFile;
  escFile.open(escapeeFile.c_str(), std::ios::in);

  bool lastEntry = false;
  std::string entry;
  while (escFile && !lastEntry)
    {
      escFile >> entry; // time
      const fType parTime = boost::lexical_cast<fType>(entry);

      if (parTime > dumpTime)
        {
          lastEntry = true;
          break;
        }

      //for (size_t i = 0; i < 7; i++)
      //escFile >> entry; // id, x, y, z, vx, vy, vz

      escFile >> entry; // id
      escFile >> entry; // x
      const fType x = boost::lexical_cast<fType>(entry);
      escFile >> entry; // y
      const fType y = boost::lexical_cast<fType>(entry);
      escFile >> entry; // z
      const fType z = boost::lexical_cast<fType>(entry);

      escFile >> entry; // vx
      const fType vx = boost::lexical_cast<fType>(entry);
      escFile >> entry; // vy
      const fType vy = boost::lexical_cast<fType>(entry);
      escFile >> entry; // vz
      const fType vz = boost::lexical_cast<fType>(entry);

      escFile >> entry; // lx
      //const fType lx = boost::lexical_cast<fType>(entry);
      escFile >> entry; // ly
      //const fType ly = boost::lexical_cast<fType>(entry);
      escFile >> entry; // lz
      //const fType lz = boost::lexical_cast<fType>(entry);

      escFile >> entry; // m
      const fType parM = boost::lexical_cast<fType>(entry);
      Mesc(0) += parM;

      Lesc(X) += parM * (y * vz + z * vy);
      Lesc(Y) += parM * (z * vx + x * vz);
      Lesc(Z) += parM * (x * vy + y * vx);

      Pesc(X) += parM * vx;
      Pesc(Y) += parM * vy;
      Pesc(Z) += parM * vz;
          
      noEscParts++;

      for (size_t i = 0; i < 10; i++)
        escFile >> entry; // h, u, p, noneigh, e, a, Torbit, mat, phase, T
    }
  escFile.close();

  CommManager.sum(Mesc);
  IOManager.savePrimitive(Mesc, "Mesc", saveFile);

  CommManager.sum(Lesc);
  IOManager.savePrimitive(Lesc, "Lesc", saveFile);
  
  valvectType Ltot = Lesc + Lrem;
  IOManager.savePrimitive(Ltot, "Ltot", saveFile);
  
  CommManager.sum(Pesc);
  IOManager.savePrimitive(Pesc, "Pesc", saveFile);
  
  valvectType Ptot = Pesc + Prem;
  IOManager.savePrimitive(Ptot, "Ptot", saveFile);

  std::cerr << "Mdisk:   " << (MdiskDun(0) + MdiskIron(0))/ Mmoon
            << " mlun  (" << noDiskIronParts + noDiskDunParts 
            << " particles)\n";
  std::cerr << " disk iron faction    "
            << 100.* MdiskIron(0) / ( MdiskDun(0) + MdiskIron(0) )
            << " % (" << noDiskIronParts << " particles)\n";
  std::cerr << "Mesc:    " << Mesc(0) / Mmoon
            << " mlun (" << noEscParts << " particles)\n";
  std::cerr << "Mplan:   " << Mpe(0) / Mearth
            << " me (" << noPlanParts << " particles)\n";
  std::cerr << "Total no. of particles: "
            << noPlanParts + noDiskIronParts
               + noDiskDunParts + noEscParts << "\n";

  ///
  /// store mass profile
  ///
  IOManager.savePrimitive(MProfIron, "Miron", saveFile);
  IOManager.savePrimitive(MProfDun, "Mdun", saveFile);

  valvectType MProf = MProfDun + MProfIron;
  IOManager.savePrimitive(MProf, "M", saveFile);


  ///
  /// cumulative bins
  ///
  for (size_t i = 1; i < noMassBins; i++)
    {
      MProf(i) += MProf(i - 1);
    }
  IOManager.savePrimitive(MProf, "Mcum", saveFile);

  std::cerr << "# time       "
            << "  Mpe        "
            << "  Lpe        "
            << "  Mdisk      "
            << "  Ldisk      "
            << "  Mesc       "
            << "  Lesc       "
            << "  Lrem       "
            << "  Ltot       "
            << "DiskIronFrac\n";

  const fType LpeA = sqrt(Lpe(X) * Lpe(X) +
                          Lpe(Y) * Lpe(Y) +
                          Lpe(Z) * Lpe(Z));

  const fType LdiskA = sqrt(Ldisk(X) * Ldisk(X) +
                            Ldisk(Y) * Ldisk(Y) +
                            Ldisk(Z) * Ldisk(Z));
  
  const fType LtotA = sqrt( Ltot(X) * Ltot(X) +
                            Ltot(Y) * Ltot(Y) +
                            Ltot(Z) * Ltot(Z));

  const fType LescA = sqrt( Lesc(X) * Lesc(X) +
                            Lesc(Y) * Lesc(Y) +
                            Lesc(Z) * Lesc(Z));

  const fType LremA = sqrt( Lrem(X) * Lrem(X) +
                            Lrem(Y) * Lrem(Y) +
                            Lrem(Z) * Lrem(Z));

  std::cout << std::scientific
            << dumpTime << " "
            << Mpe(0) << " "
            << LpeA << " "
            << MdiskDun(0) + MdiskIron(0) << " "
            << LdiskA << " "
            << Mesc(0) << " "
            << LescA << " "
            << LremA << " "
            << LtotA << " "
            << MdiskIron(0) / ( MdiskDun(0) + MdiskIron(0) ) << "\n";

  MPI::Finalize();
  return EXIT_SUCCESS;
}

