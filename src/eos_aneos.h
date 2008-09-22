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

  Logger.stream << "init ANEOS EOS with file "
                << matFilename;
  Logger.flushStream();
};

~ANEOS()
{
};

static ANEOS& instance();
static ANEOS* _instance;

///
/// get the pressure & speed of sound for particle _i
///
/// common EOS interface
///
void operator()(const size_t _i, valueType& _P, valueType& _cs)
{
  this->operator()(rho (_i), u (_i), mat (_i), _P, _cs, PartManager.phase (_i));
}

///
/// get the pressure & speed of sound for given parameters
///
void operator()(const valueType _rho, const valueType _u, const identType _mat,
                valueType& _P, valueType& _cs, identType& _phase)
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
        {
          std::cout << "Temperature out of bounds!"; exit(1);
        }
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
  Logger << "no convergence in ANEOS temperature iteration";
  // include par_Ticle data
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
  Logger << "temperature out of bounds in ANEOS";
  // include particle data
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

