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
/// the expressions for the speed of sound ( cc && ce )
/// are copied from ParaSPH. in the hybrid regime, the square
/// roots for both cc and ce are already pulled, in opposite
/// to the scheme in ParaSPH.
///
void getPressCs(const size_t& _i, valueType& _P, valueType& _cs)
{
  const valueType curE = u(_i);
  const valueType curRho = rho(_i);
  const identType curMat = mat(_i);

  if (curMat != loadedMat)
    loadMat(curMat);

  const valueType eta = curRho / rho0;
  const valueType mu = eta - 1.;
  const valueType curER = curRho * curE;
  const valueType cmin = sqrt(0.25 * A / rho0);

  const valueType k1 = 1. / ((curE / (E0 * eta * eta)) + 1);
  const valueType k3 = curE / (E0 * eta * eta);

  ///
  /// calculate Pc,cc
  ///
  if (eta > 1. || curE < Ecv)
    {
      Pc = (a + b * k1) * curER + A * mu + B * mu * mu;
      cc = sqrt((a * curE + (A + 2. * B * mu) / rho0
                 + b * curE * (3. * k3 + 1.) * k1 * k1
                 + (Pc / curRho) * (a + b * k1 * k1)));
      if (!(cc > 0.))
        cc = 0.;

      ///
      /// compressed regime or
      /// expanded, but cold regime
      ///
      if (eta > 1. || curE < Eiv)
        {
          _P = Pc;
          if (cc < cmin)
            _cs = cmin;
          else
            _cs = cc;
          return;
        }
    }

  ///
  /// if we have come so far, we
  /// are in an expanded non-cold regime
  /// -> calculate Pe,ce
  ///
  const valueType k2 = (1. / eta) - 1.;
  const valueType k4 = 1. / eta;

  const valueType exp1 = exp(-beta * k2);
  const valueType exp2 = exp(-alpha * k2 * k2);

  Pe = a * curER
       + (curER * b * k1 + A * mu * exp1) * exp2;

  ce = sqrt((b * curE * (3. * k3 + 1.) * k1 * k1
             + 2. * alpha * b * k1 * k2 * k4 * curE
             + A * exp1 * ((2. * alpha * k2 + beta) * (mu * k4 / curRho)
                           + (1. / rho0))
             ) * exp2
            + a * curE
            + (Pe / curRho) * (a + b * exp2 * k1 * k1));
  if (!(ce > 0.))
    ce = 0.;

  ///
  /// completely vaporized regime
  ///
  if (curE > Ecv)
    {
      _P = Pe;
      if (ce < cmin)
        _cs = cmin;
      else
        _cs = ce;
      return;
    }

  ///
  /// hybrid regime
  ///
  _P = ((curE - Eiv) * Pe + (Ecv - curE) * Pc) / (Ecv - Eiv);

  const valueType chy = ((curE - Eiv) * ce + (Ecv - curE) * cc) / (Ecv - Eiv);
  if (chy < cmin)
    _cs = cmin;
  else
    _cs = chy;
  return;
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

      /// ignore comments up to a size of 16384
      if (!str.compare(0, 1, "#"))
        {
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
valueType cc, ce;
matrixType param;
};

Tillotson * Tillotson::_instance = NULL;
Tillotson& Tillotson::instance()
{
  if (_instance == NULL)
    _instance = new Tillotson;
  return *_instance;
};
}
#endif

