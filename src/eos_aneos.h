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
  ///
  /// set the number of points in rho and u
  /// for the interpolation tables
  ///
  noPointsRho = 1000;
  noPointsU = 1000;

  if (PartManager.attributes.count("umin") < 1)
    throw GeneralError("umin not set for ANEOS table!");
  uMin = PartManager.attributes["umin"];

  if (PartManager.attributes.count("umax") < 1)
    throw GeneralError("umax not set for ANEOS table!");
  uMax = PartManager.attributes["umax"];

  if (PartManager.attributes.count("rhomin") < 1)
    throw GeneralError("rhomin not set for ANEOS table!");
  rhoMin = PartManager.attributes["rhomin"];

  if (PartManager.attributes.count("rhomax") < 1)
    throw GeneralError("rhomax not set for ANEOS table!");
  rhoMax = PartManager.attributes["rhomax"];

  ///
  /// the maximal mat id
  ///
  const size_t maxMatId = 32;

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
  nan = std::numeric_limits<valueType>::quiet_NaN();

  ///
  /// throw exceptions per default, keep quiet when generating tables
  ///
  quiet = false;
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

private:
valueType nan;
bool quiet;

#ifdef SPHLATCH_ANEOS_TABLE
size_t noPointsRho, noPointsU;
valueType rhoMin, rhoMax, uMin, uMax;
#endif

///
/// get the pressure & speed of sound for particle _i
///
/// common EOS interface for particle use
///
public:
void operator()(const size_t _i, valueType& _P, valueType& _cs)
{
  (*this)(rho(_i), u(_i), mat(_i), _P, _cs,
          PartManager.T(_i), PartManager.phase(_i));
}

///
/// get the pressure & speed of sound for given parameters
///
/// common EOS interface for independent use
///
public:
void operator()(const valueType _rho, const valueType _u, const identType _mat,
                valueType& _P, valueType& _cs)
{
  static identType tmpPhase;
  static valueType tmpT;

  (*this)(_rho, _u, _mat, _P, _cs, tmpT, tmpPhase);
}

///
/// get the pressure & speed of sound for given parameters
///
/// specific interface
///
public:
void operator()(const valueType _rho, const valueType _u, const identType _mat,
                valueType& _P, valueType& _cs, valueType& _T, identType& _phase)
{
#ifdef SPHLATCH_ANEOS_TABLE
  ///
  /// if rho and u are inside the talbe range, use
  /// the table. otherwise fall back to iteration.
  ///
  if (_rho > rhoMin && _rho < rhoMax &&
      _u > uMin && _u < uMax)
    {
      ///
      /// check whether for the desired material there are
      /// already tables available. if not, generate them.
      ///
      if (!tablesInit[_mat])
        initTables(_mat);

      const valueType curLogRho = log10(_rho);
      const valueType curLogU = log10(_u);

      _P = (*presTables[_mat])(curLogU, curLogRho);

      ///
      /// the first table query gives us the indices
      /// for the other table queries
      ///
      const size_t ixl = (presTables[_mat])->ixl;
      const size_t iyl = (presTables[_mat])->iyl;

      _T = (*TempTables[_mat])(curLogU, curLogRho, ixl, iyl);
      _cs = (*csouTables[_mat])(curLogU, curLogRho, ixl, iyl);

      _phase = static_cast<identType>(
        lrint((*phasTables[_mat])(curLogU, curLogRho, ixl, iyl)));
    }
  else
#endif
  ///
  /// use the slower but more accurate iteration
  ///
  iterate(_rho, _u, _mat, _P, _cs, _T, _phase);
};

private:
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
private:
void initTables(const identType _mat)
{
#ifdef SPHLATCH_LOGGER
  Logger.stream << "ANEOS generate table for mat id " << _mat;
  Logger.flushStream();
#endif
  ///
  /// prepare the argument vectors log10(rho) and log10(u)
  ///
  valvectType loguVect(noPointsU);
  const valueType dlogu = (log10(uMax) - log10(uMin)) / (noPointsU - 1);
  const valueType loguMin = log10(uMin);

  for (size_t i = 0; i < noPointsU; i++)
    {
      loguVect(i) = loguMin + dlogu * static_cast<valueType>(i);
    }

  valvectType logrhoVect(noPointsRho);
  const valueType dlogrho = (log10(rhoMax) - log10(rhoMin)) / (noPointsRho - 1);
  const valueType logrhoMin = log10(rhoMin);

  for (size_t i = 0; i < noPointsRho; i++)
    {
      logrhoVect(i) = logrhoMin + dlogrho * static_cast<valueType>(i);
    }

  ///
  /// prepare the temporary tables
  /// (this routine may take a long time)
  ///
  matrixType presTmpTable(noPointsU, noPointsRho);
  matrixType TempTmpTable(noPointsU, noPointsRho);
  matrixType csouTmpTable(noPointsU, noPointsRho);
  matrixType phasTmpTable(noPointsU, noPointsRho);

  ///
  /// inhibit NoConvergence throwing during table creation
  ///
  quiet = true;
  identType phaseInt;
  for (size_t i = 0; i < noPointsU; i++)
    {
      for (size_t j = 0; j < noPointsRho; j++)
        {
          const valueType curU = pow(10., loguVect(i));
          const valueType curRho = pow(10., logrhoVect(j));

          iterate(curRho, curU, _mat,
                  presTmpTable(i, j),
                  csouTmpTable(i, j),
                  TempTmpTable(i, j),
                  phaseInt);
          phasTmpTable(i, j) = lrint(phaseInt);
        }
    }
  quiet = false;

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
#ifdef SPHLATCH_LOGGER
  Logger.stream << "ANEOS table (" << noPointsRho << "x" << noPointsU
                << ") initialised for mat id " << _mat
                << " (range rho: " << rhoMin << "..." << rhoMax
                << ", range u: " << uMin << "..." << uMax << ")";
  Logger.flushStream();
#endif
}
#endif

///
/// get p(rho,T) and u(rho,T)
///
public:
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
        throw TempOutOfBounds(Tub, 1.e15);
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

  // just to shut up compiler
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
  /// we set all output quantities to NaN
  ///
  _Ti = nan;
  _pi = nan;
  _csi = nan;
  _kpai = -1;

  if (!quiet)
    throw NoConvergence(_rhoi, _ui, _mati);
};

private:
class NoConvergence : public GenericError
{
public:
NoConvergence(const double _rho, const double _u, const int _mat)
{
#ifdef SPHLATCH_LOGGER
  Logger.stream
#else
  std::cerr
#endif
  << "no convergence in ANEOS temperature iteration, rho = "
  << _rho << ", u = " << _u << ", mat id = " << _mat;
#ifdef SPHLATCH_LOGGER
  Logger.flushStream();
#else
  std::cerr << "\n";
#endif
};

~NoConvergence()
{
};
};

private:
class TempOutOfBounds : public GenericError
{
public:
TempOutOfBounds(const double _Tub, const double _Tmax)
{
#ifdef SPHLATCH_LOGGER
  Logger.stream
#else
  std::cerr
#endif
  << "temperature out of bounds in ANEOS, Tub = " << _Tub
  << ", Tmax = " << _Tmax;
#ifdef SPHLATCH_LOGGER
  Logger.flushStream();
#else
  std::cerr << "\n";
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
