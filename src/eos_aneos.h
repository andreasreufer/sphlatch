#ifndef SPHLATCH_EOS_ANEOS
#define SPHLATCH_EOS_ANEOS

/*
 *  eos_aneos.h
 *
 *
 *  Created by Andreas Reufer on 15.09.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include <fstream>
#include <boost/lexical_cast.hpp>

#include "typedefs.h"
#include "eos_generic.h"
#include "err_handler.h"

#ifdef SPHLATCH_ANEOS_TABLE
#include "lookup_table2D.h"
#endif

///
/// FORTRAN subroutine "ANEOSINIT" to init ANEOS
///
extern "C" void aneosinit_(const char *, /// materials file
                           int);         /// filename string length

///
/// FORTRAN subroutine "ANEOS" for p(rho,T,mat)
///
extern "C" void aneos_(const double *,  /// T (input)
                       const double *,  /// rho (input)
                       double *,    /// p
                       double *,    /// E
                       double *,    /// S
                       double *,    /// c_v
                       double *,    /// dp/dt
                       double *,    /// dp/dr
                       double *,    /// fkro (Rosseland mean opacity)
                       double *,    /// cs
                       int *,       /// phase
                       const int *, /// material (input)
                       double *,    /// fme (?)
                       double *     /// fva (?)
                       );

namespace sphlatch {
class ANEOS : public EOS {
public:
ANEOS()
{
  std::string matFilename = "aneos.input";

  aneosinit_(matFilename.c_str(), matFilename.size());

#ifdef SPHLATCH_LOGGER
  Logger.stream << "init ANEOS EOS with file "
                << matFilename;
  Logger.flushStream();
#endif

#ifdef SPHLATCH_ANEOS_TABLE
  const size_t maxMatId = 256;

  presTables.resize(maxMatId + 1);
  TempTables.resize(maxMatId + 1);
  csouTables.resize(maxMatId + 1);
  phasTables.resize(maxMatId + 1);

  tablesInit.resize(maxMatId + 1);

  for (size_t i = 0; i < maxMatId + 1; i++)
    {
      presTables[i] = NULL;
      TempTables[i] = NULL;
      csouTables[i] = NULL;
      phasTables[i] = NULL;

      tablesInit[i] = false;
    }
#endif
};

~ANEOS()
{
};

#ifdef SPHLATCH_ANEOS_TABLE
typedef LookupTable2D<InterpolateBilinear> lut_type;
std::vector<lut_type*> presTables, TempTables, csouTables, phasTables;
std::vector<bool> tablesInit;
#endif

static ANEOS& instance();
static ANEOS* _instance;

///
/// get the pressure & speed of sound for particle _i
///
/// common EOS interface for particle use
///
void operator()(const size_t _i, valueType& _P, valueType& _cs)
{
  this->operator()(rho (_i), u (_i), mat (_i), _P, _cs,
                   PartManager.T (_i), PartManager.phase (_i));
}

///
/// get the pressure & speed of sound for given parameters
///
/// common EOS interface for independent use
///
void operator()(const valueType _rho, const valueType _u, const identType _mat,
                valueType& _P, valueType& _cs)
{
  static identType tmpPhase;
  static valueType tmpT;

  this->operator()(_rho, _u, _mat, _P, _cs, tmpT, tmpPhase);
}

///
/// get the pressure & speed of sound for given parameters
///
/// specific interface
///
void operator()(const valueType _rho, const valueType _u, const identType _mat,
                valueType& _P, valueType& _cs, valueType& _T, identType& _phase)
{
#ifdef SPHLATCH_ANEOS_TABLE
  ///
  /// check whether for the desired material there are
  /// already tables available
  ///
  if (tablesInit[_mat] == false)
    {
      initTables(_mat);
    }

  const valueType curLogRho = log(_rho);
  const valueType curLogU   = log(_u);

  _P = (presTables[_mat])->operator()(curLogU, curLogRho);


#else
  iterate(_rho, _u, _mat, _P, _cs, _T, _phase);
#endif

};

void iterate(const valueType _rho, const valueType _u, const identType _mat,
                valueType& _P, valueType& _cs, valueType& _T, identType& _phase)
{
  static double curRho, curU, curP, curCs, curT;
  static int curMat, curPhase;

  curRho = static_cast<double>(_rho);
  curU = static_cast<double>(_u);
  curMat = static_cast<int>(_mat);

  rooten(curRho, curU, curMat, curT, curP, curCs, curPhase);

  _phase = static_cast<identType>(curPhase);
  _P = static_cast<valueType>(curP);
  _cs = static_cast<valueType>(curCs);
  _T = static_cast<valueType>(curT);
};

#ifdef SPHLATCH_ANEOS_TABLE
///
/// initialize tables for a material
///
void initTables(const identType _mat)
{
  ///
  /// prepare the argument vectors log(rho) and log(u)
  ///
  const valueType uMin = 1.e-5;
  const valueType uMax = 10.;
  const size_t nu = 100;
  valvectType loguVect(nu);

  const valueType dlogu = ( log(uMax) - log(uMin) )/(nu-1);
  const valueType loguMin = log(uMin);

  for (size_t i = 0; i < nu; i++)
    {
      loguVect(i) = loguMin + dlogu*static_cast<valueType>(i);
    }

  const valueType rhoMin = 1.e-3;
  const valueType rhoMax = 10.;
  const size_t nrho = 200;
  valvectType logrhoVect(nrho);

  const valueType dlogrho = ( log(rhoMax) - log(rhoMin) )/(nrho-1);
  const valueType logrhoMin = log(rhoMin);

  for (size_t i = 0; i < nrho; i++)
    {
      logrhoVect(i) = logrhoMin + dlogrho*static_cast<valueType>(i);
    }
  
  ///
  /// prepare the temporary tables
  /// (this routine may take a long time)
  ///
  matrixType presTmpTable(nu, nrho);
  matrixType TempTmpTable(nu, nrho);
  matrixType csouTmpTable(nu, nrho);
  matrixType phasTmpTable(nu, nrho);

  identType phaseInt;

  for (size_t i = 0; i < nu; i++)
  {
    for (size_t j = 0; j < nrho; j++)
    {
      const valueType curU   = exp(loguVect(i));
      const valueType curRho = exp(logrhoVect(j));

      iterate(curRho, curU, _mat,
              presTmpTable(i,j),
              csouTmpTable(i,j),
              TempTmpTable(i,j),
              phaseInt);

      //std::cout << curRho << " " << curU << " " << presTmpTable(i,j) << "\n";

      phasTmpTable(i,j) = lrint(phaseInt);
    }
  }
  
  ///
  /// instantate the lookup tables
  ///
  presTables[_mat] = new lut_type(loguVect, logrhoVect, presTmpTable);
  TempTables[_mat] = new lut_type(loguVect, logrhoVect, TempTmpTable);
  csouTables[_mat] = new lut_type(loguVect, logrhoVect, csouTmpTable);
  phasTables[_mat] = new lut_type(loguVect, logrhoVect, phasTmpTable);

  ///
  /// flag the tables for material _mat as initialized
  ///
  tablesInit[_mat] = true;
}
#endif

///
/// get p(rho,T) and u(rho,T)
///
void getSpecEnergy(const valueType _rho, const valueType _T,
                   const identType _mat,
                   valueType& _p, valueType& _cs, valueType& _u)
{
  static double T, rho, p, u, S, cv, dpdt, dpdr, fkros, cs, fme, fma;
  static int kpa, mat;

  rho = static_cast<double>(_rho);
  T = static_cast<double>(_T);
  mat = static_cast<int>(_mat);

  aneos_(&T, &rho, &p, &u, &S, &cv, &dpdt, &dpdr, &fkros,
         &cs, &kpa, &mat, &fme, &fma);

  _p = static_cast<valueType>(p);
  _cs = static_cast<valueType>(cs);
  _u = static_cast<valueType>(u);
};

private:
///
/// routine to iteratively find the temperature for a given
/// internal energy and density
///
/// copied from ParaSPH, 2004 by Bruno Nyffeler & Willy Benz
///
void rooten(const double _rhoi, const double _ui, const int _mati,
            double &_Ti, double &_pi, double &_csi, int &_kpai) const
{
  const double eps = 1.e-5;
  const double Tmin = 1.e-6; // get this as a parameter?
  const int itmax = 30;
  static double _S, _CV, _DPDT, _DPDR, _FKROS, _FME, _FMA;
  static double a, b, c, d, e, ei, fa, fb, fc, p, q, r, s, tm, tol1;

  // Initial temperature bracket (in eV)
  static double Tlb, Tub;

  Tlb = 0.001;
  Tub = 6.0;

  // Check lower boundary
  for (fa = 0.0; fa >= 0.0; Tlb *= 0.1)
    {
      // minimal temperature and enery
      if (Tlb < Tmin)
        {
          //_ui = ei; this is not the EOS job
          _Ti = a;
          return;
        }
      a = Tlb;
      aneos_(&a, &_rhoi, &_pi, &ei, &_S, &_CV, &_DPDT, &_DPDR, &_FKROS,
             &_csi, &_kpai, &_mati, &_FME, &_FMA);
      // fa = trial energy - req energy
      fa = ei - _ui;
    }
  // a: lower bound for T
  // fb: delta to lower bound energy

  // Check upper boundary
  for (fb = 0.0; fb <= 0.0; Tub *= 3.0)
    {
      if (Tub > 1.e15)
        throw TempOutOfBounds();
      b = Tub;
      aneos_(&b, &_rhoi, &_pi, &ei, &_S, &_CV, &_DPDT, &_DPDR, &_FKROS,
             &_csi, &_kpai, &_mati, &_FME, &_FMA);
      // fa = trial energy - req energy
      fb = ei - _ui;
    }
  // b: upper bound for T
  // fb: delta to upper bound energy

  // Start iteration
  fc = fb; // fc -> upper energy bound

  // just to shut up com_piler
  c = a;     // c -> lower temperature bound
  d = b - a; // d is temperature range
  e = d;
  q = 0;

  for (int i = 0; i < itmax; i++)
    {
      if (fb * fc > 0) // why should this be negative?
        {
          c = a;
          fc = fa;
          d = b - a;
          e = d;
        }
      if (fabs(fc) < fabs(fb))
        {
          a = b;
          b = c;
          c = a;
          fa = fb;
          fb = fc;
          fc = fa;
        }
      tm = 0.5 * (c - b);
      tol1 = 2. * eps * fabs(b);
      if (fabs(tm) < tol1 || fabs(fb / _ui) < eps)
        {
          _Ti = b;
          return;
        }
      if (fabs(e) > tol1 && fabs(fa) > fabs(fb))
        {
          s = fb / fa;
          if (a == c)
            {
              p = 2. * tm * s;
            }
          else
            {
              q = fa / fc; r = fb / fc;
              p = s * (2. * tm * q * (q - r) - (b - a) * (r - 1.));
              q = (q - 1.) * (r - 1.) * (s - 1.);
            }
          // there might be a problem with q here
          if (p > 0.)
            q = -q;
          p = fabs(p);
          if ((2. * p) < (3. * tm * q - fabs(tol1 * q)) && 2. * p < fabs(e * q))
            {
              e = d;
              d = p / q;
            }
          else
            {
              d = tm;
              e = d;
            }
        }
      else
        {
          d = tm;
          e = d;
        }
      a = b;
      fa = fb;
      if (fabs(d) > tol1)
        b += d;
      else
        {
          if (tm >= 0.)
            b += fabs(tol1);
          else
            b -= fabs(tol1);
        }
      aneos_(&b, &_rhoi, &_pi, &ei, &_S, &_CV, &_DPDT, &_DPDR, &_FKROS,
             &_csi, &_kpai, &_mati, &_FME, &_FMA);
      fb = ei - _ui;
    }
  ///
  /// arriving here means no convergence
  ///
  throw NoConvergence();
};

class NoConvergence : public GenericError
{
public:
NoConvergence()
{
#ifdef SPHLATCH_LOGGER
  Logger << "no convergence in ANEOS temperature iteration";
  // include particle data
#else
  std::cerr << "no convergence in ANEOS temperature iteration\n";
#endif
};

~NoConvergence()
{
};
};

class TempOutOfBounds : public GenericError
{
public:
TempOutOfBounds()
{
// include particle data
#ifdef SPHLATCH_LOGGER
  Logger << "temperature out of bounds in ANEOS";
#else
  std::cerr << "temperature out of bounds in ANEOS\n";
#endif
};

~TempOutOfBounds()
{
};
};
};

ANEOS * ANEOS::_instance = NULL;
ANEOS& ANEOS::instance()
{
  if (_instance == NULL)
    _instance = new ANEOS;
  return *_instance;
};
}
#endif
