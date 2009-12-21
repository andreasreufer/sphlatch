#ifndef SPHLATCH_LAGRANGE_SPHERE1D_SOLVER
#define SPHLATCH_LAGRANGE_SPHERE1D_SOLVER

/*
 *  lagrange_sphere1D_solver.h
 *
 *
 *  Created by Andreas Reufer on 03.08.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

#ifdef SPHLATCH_TILLOTSON
#include "eos_tillotson.h"
typedef sphlatch::Tillotson eosType;
#endif

#ifdef SPHLATCH_ANEOS
#include "eos_aneos.h"
typedef sphlatch::ANEOS eosType;
#endif

#include <boost/progress.hpp>

#include <float.h>

namespace sphlatch {
class LagrangeSphere1DSolver {
public:
LagrangeSphere1DSolver(size_t _noCells) :
  EOS(eosType::instance())
{
  firstStep = true;
  time = 0.;

  gravConst = 6.67e-8;
  courantNumber = 0.3;
  friction = 1.e-9;

  noCells = _noCells;
  noEdges = noCells + 1;

  avAlpha = 1.;
  avBeta = 2.;

  rho.resize(noCells);
  q.resize(noCells);
  u.resize(noCells);
  p.resize(noCells);
  m.resize(noCells);
  cs.resize(noCells);
  pow.resize(noCells);
  mat.resize(noCells);
  rhoOld.resize(noCells);
  qOld.resize(noCells);
#ifdef SPHLATCH_ANEOS
  T.resize(noCells);
  phase.resize(noCells);
#endif

  r.resize(noEdges);
  v.resize(noEdges);
};

~LagrangeSphere1DSolver()
{
};

public:
void integrateTo(fType _t);
void densityToMass();

public:
valvectType r, v, rho, q, u, p, m, cs, pow;
idvectType mat;
#ifdef SPHLATCH_ANEOS
valvectType T;
idvectType phase;
#endif
fType time, uMin, friction;
fType gravConst, courantNumber;

private:
void integrate(fType _dtMax);
void bootstrap();
fType detTimestep();

eosType& EOS;

private:
size_t noCells, noEdges;

valvectType rhoOld, qOld;

fType avAlpha, avBeta;
fType dtp05, dtm05;

bool firstStep;
};


void LagrangeSphere1DSolver::integrateTo(fType _t)
{
  const size_t dispSteps = 1000;
  size_t curStep = 0;

  const fType totTime = _t - time;
  const fType startTime = time;

  if (firstStep)
    bootstrap();

  boost::progress_display show_progress(dispSteps);
  while (time < _t)
  {
    integrate(_t - time);

    while ( 
    ( static_cast<fType>(curStep) /
      static_cast<fType>(dispSteps) ) < ( ( time - startTime ) / totTime ) )
    {
      curStep++;
      ++show_progress;
    }
  }
}

void LagrangeSphere1DSolver::integrate(fType _dtMax)
{
  const fType pi4 = M_PI * 4.;
  
  dtp05 = std::min(detTimestep(), _dtMax);

  fType dtc00;
  if (firstStep)
    {
      dtc00 = dtp05;
      firstStep = false;
    }
  else
    dtc00 = 0.5 * (dtp05 + dtm05);

  ///
  /// update velocity
  ///
  fType mSum = 0.;
  for (size_t i = 1; i < noEdges - 1; i++)
    {
      const fType dmk = 0.5 * (m(i - 1) + m(i));
      const fType Ak = pi4 * r(i) * r(i);
      const fType Akm10 = pi4 * r(i - 1) * r(i - 1);
      const fType Akp10 = pi4 * r(i + 1) * r(i + 1);
      const fType Akm05 = 0.5 * (Akm10 + Ak);
      const fType Akp05 = 0.5 * (Akp10 + Ak);

      const fType gradP = Ak * (p(i) - p(i - 1));
      const fType gradQ = 0.5 * (q(i) * (3. * Akp05 - Ak) - q(i - 1) * (3. * Akm05 - Ak));

      mSum += m(i - 1);
      const fType accG = gravConst * mSum / (r(i) * r(i));

      v(i) -= ((gradP + gradQ) / dmk + accG + friction*v(i)) * dtc00;
      assert(!isnan(v(i)));
    }
  v(0) = 0.;
  v(noEdges - 1) = v(noEdges - 2);

  ///
  /// update position
  ///
  r(0) = 0.;
  for (size_t i = 1; i < noEdges; i++)
    {
      r(i) += v(i) * dtp05;
      assert(!isnan(r(i)));
    }

  ///
  /// calculate density
  ///
  for (size_t i = 0; i < noCells; i++)
    {
      rhoOld(i) = rho(i);
      const fType rkpow3 = r(i) * r(i) * r(i);
      const fType rkp1pow3 = r(i + 1) * r(i + 1) * r(i + 1);
      rho(i) = 3. * m(i) / (pi4 * (rkp1pow3 - rkpow3));
      assert(!isnan(rho(i)));
    }
  rho(noCells-1) = 0.; /// vacuum boundary condition

  ///
  /// update q-factor
  ///
  for (size_t i = 0; i < noCells; i++)
    {
      if (v(i + 1) < v(i))
        {
          const fType Ak = pi4 * r(i) * r(i);
          const fType Akp10 = pi4 * r(i + 1) * r(i + 1);
          const fType Akp05 = 0.5 * (Akp10 + Ak);

          const fType k1 = (v(i + 1) * (Akp05 - Akp10 / 3.)
                                - v(i) * (Akp05 - Ak / 3.)) / Akp05;
          const fType absDivV = v(i) - v(i + 1);
          const fType rhop05 = 0.5*(rho(i) + rhoOld(i));

          q(i) = -(3./2.)*rhop05*k1*( avAlpha*cs(i)+ avBeta*avBeta*absDivV );
          assert(!isinf(q(i)));
        }
        else
          q(i) = 0.;
    }

  ///
  /// update internal energy
  ///
  /// the last cell is not needed, as it is vacuum anyway
  ///
  fType pPrime = 0.;
  for (size_t i = 0; i < noCells-1; i++)
    {
      pPrime = 0.;
      const fType Ak = pi4 * r(i) * r(i);
      const fType Akp10 = pi4 * r(i + 1) * r(i + 1);
      const fType Akp05 = 0.5 * (Akp10 + Ak);

      ///
      /// try uPrime and guess pPrime
      ///
      const fType dVol = (Akp10 * v(i + 1) - Ak * v(i)) / m(i);
      const fType uPrime = u(i) - p(i) * dVol * dtp05;
#ifdef SPHLATCH_ANEOS
      EOS(rho(i), uPrime, mat(i), pPrime, cs(i), T(i), phase(i));
#else
      EOS(rho(i), uPrime, mat(i), pPrime, cs(i));
#endif
      p(i) = 0.5 * (p(i) + pPrime);
      q(i) = 0.5 * (q(i) + qOld(i));
      assert(!isnan(pPrime));
      assert(!isnan(uPrime));
      assert(!isnan(q(i)));

      ///
      /// AV heating
      ///
      const fType dQ = (v(i + 1) * (3. * Akp05 - Akp10)
                            - v(i) * (3. * Akp05 - Ak)) / m(i);
      u(i) -= (p(i) * dVol + dQ) * dtp05;
      pow(i) = -(p(i) * dVol + dQ);
      assert(!isnan(u(i)));
      assert(!isnan(pow(i)));

      if (u(i) < uMin)
        u(i) = uMin;
    }

  ///
  /// update pressure and speed of sound
  ///
  for (size_t i = 0; i < noCells-1; i++)
    {
#ifdef SPHLATCH_ANEOS
      EOS(rho(i), u(i), mat(i), p(i), cs(i), T(i), phase(i));
#else
      EOS(rho(i), u(i), mat(i), p(i), cs(i));
#endif
      assert(!isnan(p(i)));
      assert(!isnan(cs(i)));

      if (p(i) < 0.)
        p(i) = 0.;
    }
  p(noCells-1) = 0.; /// vacuum boundary condition

  dtm05 = dtp05;
  time += dtp05;
}

fType LagrangeSphere1DSolver::detTimestep()
{
  ///
  /// calculate current timestep dtp05
  ///
  fType dtCFL = std::numeric_limits<fType>::max();
  fType dtVsc = std::numeric_limits<fType>::max();
  fType dtPow = std::numeric_limits<fType>::max();

  for (size_t i = 0; i < noCells; i++)
    {
      const fType dr = r(i + 1) - r(i);
      const fType dv = v(i + 1) - v(i);
      const fType absDv = fabs(dv);

      const fType dtCFLi = dr / (cs(i) + absDv);
      if ( dtCFLi < 1.e-3 )
        std::cerr << i << ": dtCFLi " << dtCFLi  << " csi " << cs(i) << "  absDv " << absDv << "  dr " << dr << "\n";
      dtCFL = dtCFLi < dtCFL ? dtCFLi : dtCFL;

      const fType dtVsci = 1. / (4. * avBeta * fabs(v(i) + v(i + 1)) / dr);
      if ( dtVsc < 1.e-3 )
        std::cerr << i << ": dtVsc " << dtVsc << " vi " << v(i) << "  vi+1 " <<  v(i + 1) << "  dr " << dr << "\n";
      dtVsc = dtVsci < dtVsc ? dtVsci : dtVsc;

      if (pow(i) < 0.)
        {
          const fType dtPowi = -u(i) / pow(i);
          if ( dtPowi < 1.e-3 )
            std::cerr << i << ": dtPowi " << dtPowi << "  u(i) " << u(i) << "  powi " << pow(i) << "\n";
            
          dtPow = dtPowi < dtPow ? dtPowi : dtPow;
        }
    }
  dtCFL *= courantNumber;

  return std::min(std::min(dtVsc, dtCFL), dtPow);
}

void LagrangeSphere1DSolver::bootstrap()
{
  const fType pi4 = M_PI * 4.;

  ///
  /// calculate density
  ///
  for (size_t i = 0; i < noCells; i++)
    {
      const fType rkpow3 = r(i) * r(i) * r(i);
      const fType rkp1pow3 = r(i + 1) * r(i + 1) * r(i + 1);
      rho(i) = 3. * m(i) / (pi4 * (rkp1pow3 - rkpow3));
      rhoOld(i) = rho(i);
    }
  rho(noCells-1) = 0.; /// vacuum boundary condition
  rhoOld(noCells-1) = 0.; /// vacuum boundary condition
  
  ///
  /// update pressure and speed of sound
  ///
  for (size_t i = 0; i < noCells-1; i++)
    {
#ifdef SPHLATCH_ANEOS
      EOS(rho(i), u(i), mat(i), p(i), cs(i), T(i), phase(i));
#else
      EOS(rho(i), u(i), mat(i), p(i), cs(i));
#endif
    }
  p(noCells-1) = 0.; /// vacuum boundary condition
  
  ///
  /// zero stuff
  ///
  for (size_t i = 0; i < noCells; i++)
    {
      pow(i) = 0.;
      qOld(i) = 0.;
    }
}

void LagrangeSphere1DSolver::densityToMass()
{
  ///
  /// calculate mass from density
  ///
  const fType pi4 = M_PI * 4.;

  for (size_t i = 0; i < noCells; i++)
    {
      const fType rkpow3 = r(i) * r(i) * r(i);
      const fType rkp1pow3 = r(i + 1) * r(i + 1) * r(i + 1);
      m(i) = (rho(i) * pi4 * (rkp1pow3 - rkpow3)) / 3.;
    }
}
}
#endif

