#ifndef SPHLATCH_EOS_TILLOTSON
#define SPHLATCH_EOS_TILLOTSON

/*
 *  eos_tillotson.h
 *
 *
 *  Created by Andreas Reufer on 26.07.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include <fstream>
#include <boost/lexical_cast.hpp>

#include "typedefs.h"
#include "eos_generic.h"

namespace sphlatch {
class Tillotson : public EOS {
public:
Tillotson()
{
  initParams("tillotson.txt");
  loadMat(1);
  Logger.stream << "init Tillotson EOS with "
                << param.size1()
                << " materials";
  Logger.flushStream();
};
~Tillotson()
{
};

static Tillotson& instance();
static Tillotson* _instance;

///
/// get the pressure for particle _i
///
valueType getPressure(const size_t& _i)
{
  const valueType curE = u(_i);
  const valueType curRho = rho(_i);
  const identType curMat = mat(_i);

  if (curMat != loadedMat)
    loadMat(curMat);

  const valueType eta = curRho / rho0;
  const valueType mu = eta - 1.;
  const valueType curER = curRho * curE;
  const valueType k1 = b / ((curE / (E0 * eta * eta)) + 1);

  ///
  /// calculate Pc
  ///
  if (curE < Ecv)
    {
      Pc = (a + k1) * curER + A * mu + B * mu * mu;
    }

  ///
  /// calculate Pe
  ///
  if (curE > Eiv)
    {
      Pe = a * curER
           + (curER * k1 + A * mu * exp(-beta * mu)) * exp(-alpha * mu * mu);
    }

  ///
  /// compressed regime
  ///
  if (curE < Eiv)
    return Pc;

  ///
  /// expanded regime
  ///
  if (curE > Ecv)
    return Pe;

  ///
  /// hybrid regime (between compressed & expanded)
  ///
  return(((curE - Eiv) * Pe + (Ecv - curE) * Pc) / (Ecv - Eiv));
};

///
/// get the speed of sound for particle _i
///
valueType getSpeedOfSound(const size_t& _i)
{
  return 0.;
};

///
/// get the temperature for particle _i
///
valueType getTemperature(const size_t& _i)
{
  return 0.;
};

private:
void loadMat(const size_t _mat)
{
  ///
  /// mat index is shifted by one to ensure
  /// compatibility with FORTRAN codes
  ///
  const size_t matIdx = _mat - 1;

  rho0 = param(matIdx, 0);
  A = param(matIdx, 1);
  B = param(matIdx, 2);
  a = param(matIdx, 3);
  b = param(matIdx, 4);
  alpha = param(matIdx, 5);
  beta = param(matIdx, 6);
  E0 = param(matIdx, 7);
  Eiv = param(matIdx, 8);
  Ecv = param(matIdx, 9);

  loadedMat = _mat;
};

///
/// load the parameter file
///
void initParams(std::string _filename)
{
  std::fstream fin;
  fin.open(_filename.c_str(), std::ios::in);

  size_t entry = 0, i = 0, noEntries = 0, noParams = 10;
  bool noEntriesRead = false;

  while (!fin.eof())
    {
      std::string str;
      fin >> str;

      if (!str.compare(0, 1, "#"))
        {
          /// ignore comments up to a size of 16384
          fin.ignore(16384, '\n');
        }
      else
        {
          if (!noEntriesRead)
            {
              noEntries = boost::lexical_cast<int>(str);
              param.resize(noEntries, noParams);
              noEntriesRead = true;
            }
          else
          if (entry < noEntries)
            {
              param(entry, i) = boost::lexical_cast<valueType>(str);
              i++;
              if (i == 10)
                {
                  i = 0;
                  entry++;
                }
            }
        }
    }
  fin.close();
};


private:
identType loadedMat;
valueType rho0, a, b, A, B, alpha, beta, E0, Eiv, Ecv;
valueType Pc, Pe;
matrixType param;
};

Tillotson* Tillotson::_instance = NULL;
Tillotson& Tillotson::instance()
{
  if (_instance == NULL)
    _instance = new Tillotson;
  return *_instance;
};
}
#endif
