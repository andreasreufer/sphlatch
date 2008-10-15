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
#include "err_handler.h"

namespace sphlatch {
class Tillotson : public EOS {
public:
Tillotson()
{
  initParams("tillotson.txt");
  ldMat = getMatParams(1);
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
/// struct for material constants
///
/// rho0:           initial density
/// A:              bulk modulus
/// B:              non-linear Tillotson compression coefficient
/// a,b,alpha,beta: Tillotson parameters (dimensionless)
/// E0:             initial energy
/// Eiv:            energy of incipient vaporization
/// Ecv:            energy of complete vaporization
///
/// xmu:            shear modulus
/// umelt:          melt specific energy
/// yie:            plastic yielding
/// pweib:          Weibull parameter p (called m in Benz & Asphaug 1994)
/// cweib:          Weibull parameter c (called k in Benz & Asphaug 1994)
/// J2slu:          sine of slope of J2 for undamaged material
/// coh:            cohesion
/// J2sld:          sine of slope of J2 for damaged material
///
struct paramType {
  identType id;
  valueType rho0, a, b, A, B, alpha, beta, E0, Eiv, Ecv;
  valueType xmu, umelt, yie, pweib, cweib, J2sld, coh, J2slu;
};

paramType ldMat;

///
/// get the pressure & speed of sound for particle _i
///
/// common EOS interface
///
void operator()(const size_t _i, valueType& _P, valueType& _cs)
{
  this->operator()(rho (_i), u (_i), mat (_i), _P, _cs);
}

///
/// get the pressure & speed of sound for given parameters
///
/// the expressions for the speed of sound ( cc && ce )
/// are copied from ParaSPH. in the hybrid regime, the square
/// roots for both cc and ce are already taken, in opposite
/// to the scheme in ParaSPH where the square root of the hybrid
/// formula is taken.
///
void operator()(const valueType _rho, const valueType _u,
                const identType _matId, valueType& _P, valueType& _cs)
{
  const valueType curRho = _rho;
  const valueType curE = _u;
  const identType curMatId = _matId;

  ///
  /// load the correct material, when
  /// in the cache
  ///
  if (curMatId != ldMat.id)
    ldMat = getMatParams(curMatId);

  const valueType eta = curRho / ldMat.rho0;
  const valueType mu = eta - 1.;
  const valueType curER = curRho * curE;
  const valueType cmin = sqrt(0.25 * ldMat.A / ldMat.rho0);

  const valueType k1 = 1. / ((curE / (ldMat.E0 * eta * eta)) + 1);
  const valueType k3 = curE / (ldMat.E0 * eta * eta);

  ///
  /// calculate Pc,cc
  ///
  if (eta > 1. || curE < ldMat.Ecv)
    {
      Pc = (ldMat.a + ldMat.b * k1) * curER + ldMat.A * mu + ldMat.B * mu * mu;
      cc = sqrt((ldMat.a * curE + (ldMat.A + 2. * ldMat.B * mu) / ldMat.rho0
                 + ldMat.b * curE * (3. * k3 + 1.) * k1 * k1
                 + (Pc / curRho) * (ldMat.a + ldMat.b * k1 * k1)));
      if (!(cc > 0.))
        cc = 0.;

      ///
      /// compressed regime or
      /// expanded, but cold regime
      ///
      if (eta > 1. || curE < ldMat.Eiv)
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

  const valueType exp1 = exp(-ldMat.beta * k2);
  const valueType exp2 = exp(-ldMat.alpha * k2 * k2);

  Pe = ldMat.a * curER
       + (curER * ldMat.b * k1 + ldMat.A * mu * exp1) * exp2;

  ce = sqrt((ldMat.b * curE * (3. * k3 + 1.) * k1 * k1
             + 2. * ldMat.alpha * ldMat.b * k1 * k2 * k4 * curE
             + ldMat.A * exp1 * ((2. * ldMat.alpha * k2 + ldMat.beta)
                                 * (mu * k4 / curRho)
                                 + (1. / ldMat.rho0))
             ) * exp2
            + ldMat.a * curE
            + (Pe / curRho) * (ldMat.a + ldMat.b * exp2 * k1 * k1));
  if (!(ce > 0.))
    ce = 0.;

  ///
  /// completely vaporized regime
  ///
  if (curE > ldMat.Ecv)
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
  _P = ((curE - ldMat.Eiv) * Pe + (ldMat.Ecv - curE) * Pc) / (ldMat.Ecv - ldMat.Eiv);
  const valueType chy = ((curE - ldMat.Eiv) * ce + (ldMat.Ecv - curE) * cc)
                        / (ldMat.Ecv - ldMat.Eiv);
  if (chy < cmin)
    _cs = cmin;
  else
    _cs = chy;
  return;
};

public:
paramType getMatParams(const size_t _matId)
{
  const size_t matIdx = _matId - 1;
  paramType _ret;

  _ret.rho0 = param(matIdx, 0);
  _ret.A = param(matIdx, 1);
  _ret.B = param(matIdx, 2);
  _ret.a = param(matIdx, 3);
  _ret.b = param(matIdx, 4);
  _ret.alpha = param(matIdx, 5);
  _ret.beta = param(matIdx, 6);
  _ret.E0 = param(matIdx, 7);
  _ret.Eiv = param(matIdx, 8);
  _ret.Ecv = param(matIdx, 9);

  _ret.xmu = param(matIdx, 10);
  _ret.umelt = param(matIdx, 11);
  _ret.yie = param(matIdx, 12);
  _ret.pweib = param(matIdx, 13);
  _ret.cweib = param(matIdx, 14);
  _ret.J2slu = param(matIdx, 15);
  _ret.coh = param(matIdx, 16);
  _ret.J2sld = param(matIdx, 17);

  _ret.id = _matId;

  return _ret;
}

///
/// load the parameter file
///
void initParams(std::string _filename)
{
  std::fstream fin;

  fin.open(_filename.c_str(), std::ios::in);

  if (!fin)
    throw FileNotFound(_filename);

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
              /// 18 parameters in each entry
              if (i == 18)
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
