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
#include "eos_tillotson.h"
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

  r.resize(noEdges);
  v.resize(noEdges);
};

~LagrangeSphere1DSolver()
{
};

public:
void integrateTo(valueType _t);
void densityToMass();

public:
valvectType r, v, rho, q, u, p, m, cs, pow;
idvectType mat;
valueType time, uMin, friction;
valueType gravConst, courantNumber;

private:
void integrate(valueType _dtMax);
void bootstrap();
valueType detTimestep();

typedef sphlatch::Tillotson eosType;
eosType& EOS;

private:
size_t noCells, noEdges;

valvectType rhoOld, qOld;

valueType avAlpha, avBeta;
valueType dtp05, dtm05;

bool firstStep;
};


void LagrangeSphere1DSolver::integrateTo(valueType _t)
{
  if (firstStep)
    bootstrap();
  while (time < _t)
    integrate(_t - time);
}

void LagrangeSphere1DSolver::integrate(valueType _dtMax)
{
  const valueType pi4 = M_PI * 4.;
  
  dtp05 = std::min(detTimestep(), _dtMax);

  valueType dtc00;
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
  valueType mSum = 0.;
  for (size_t i = 1; i < noEdges - 1; i++)
    {
      const valueType dmk = 0.5 * (m(i - 1) + m(i));
      const valueType Ak = pi4 * r(i) * r(i);
      const valueType Akm10 = pi4 * r(i - 1) * r(i - 1);
      const valueType Akp10 = pi4 * r(i + 1) * r(i + 1);
      const valueType Akm05 = 0.5 * (Akm10 + Ak);
      const valueType Akp05 = 0.5 * (Akp10 + Ak);

      const valueType gradP = Ak * (p(i) - p(i - 1));
      const valueType gradQ = 0.5 * (q(i) * (3. * Akp05 - Ak) - q(i - 1) * (3. * Akm05 - Ak));

      mSum += m(i - 1);
      const valueType accG = gravConst * mSum / (r(i) * r(i));

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
      const valueType rkpow3 = r(i) * r(i) * r(i);
      const valueType rkp1pow3 = r(i + 1) * r(i + 1) * r(i + 1);
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
          const valueType Ak = pi4 * r(i) * r(i);
          const valueType Akp10 = pi4 * r(i + 1) * r(i + 1);
          const valueType Akp05 = 0.5 * (Akp10 + Ak);

          const valueType k1 = (v(i + 1) * (Akp05 - Akp10 / 3.)
                                - v(i) * (Akp05 - Ak / 3.)) / Akp05;
          const valueType absDivV = v(i) - v(i + 1);
          const valueType rhop05 = 0.5*(rho(i) + rhoOld(i));

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
  valueType pPrime = 0.;
  for (size_t i = 0; i < noCells-1; i++)
    {
      pPrime = 0.;
      const valueType Ak = pi4 * r(i) * r(i);
      const valueType Akp10 = pi4 * r(i + 1) * r(i + 1);
      const valueType Akp05 = 0.5 * (Akp10 + Ak);

      ///
      /// try uPrime and guess pPrime
      ///
      const valueType dVol = (Akp10 * v(i + 1) - Ak * v(i)) / m(i);
      const valueType uPrime = u(i) - p(i) * dVol * dtp05;
      EOS(rho(i), uPrime, mat(i), pPrime, cs(i));
      p(i) = 0.5 * (p(i) + pPrime);
      q(i) = 0.5 * (q(i) + qOld(i));
      assert(!isnan(pPrime));
      assert(!isnan(uPrime));
      assert(!isnan(q(i)));

      ///
      /// AV heating
      ///
      const valueType dQ = (v(i + 1) * (3. * Akp05 - Akp10)
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
  for (size_t i = 0; i < noCells; i++)
    {
      EOS(rho(i), u(i), mat(i), p(i), cs(i));
      assert(!isnan(p(i)));
      assert(!isnan(cs(i)));
    }
  p(noCells-1) = 0.; /// vacuum boundary condition

  dtm05 = dtp05;
  time += dtp05;
}

valueType LagrangeSphere1DSolver::detTimestep()
{
  ///
  /// calculate current timestep dtp05
  ///
  valueType dtCFL = std::numeric_limits<valueType>::max();
  valueType dtVsc = std::numeric_limits<valueType>::max();
  valueType dtPow = std::numeric_limits<valueType>::max();

  for (size_t i = 0; i < noCells; i++)
    {
      const valueType dr = r(i + 1) - r(i);
      const valueType dv = v(i + 1) - v(i);
      const valueType absDv = fabs(dv);

      const valueType dtCFLi = dr / (cs(i) + absDv);
      dtCFL = dtCFLi < dtCFL ? dtCFLi : dtCFL;

      const valueType dtVsci = 1. / (4. * avBeta * fabs(v(i) + v(i + 1)) / dr);
      dtVsc = dtVsci < dtVsc ? dtVsci : dtVsc;

      if (pow(i) < 0.)
        {
          const valueType dtPowi = -u(i) / pow(i);
          dtPow = dtPowi < dtPow ? dtPowi : dtPow;
        }
    }
  dtCFL *= courantNumber;

  return std::min(std::min(dtVsc, dtCFL), dtPow);
}

void LagrangeSphere1DSolver::bootstrap()
{
  const valueType pi4 = M_PI * 4.;

  ///
  /// calculate density
  ///
  for (size_t i = 0; i < noCells; i++)
    {
      const valueType rkpow3 = r(i) * r(i) * r(i);
      const valueType rkp1pow3 = r(i + 1) * r(i + 1) * r(i + 1);
      rho(i) = 3. * m(i) / (pi4 * (rkp1pow3 - rkpow3));
      rhoOld(i) = rho(i);
    }
  rho(noCells-1) = 0.; /// vacuum boundary condition
  rhoOld(noCells-1) = 0.; /// vacuum boundary condition
  
  ///
  /// update pressure and speed of sound
  ///
  for (size_t i = 0; i < noCells; i++)
    {
      EOS(rho(i), u(i), mat(i), p(i), cs(i));
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
  const valueType pi4 = M_PI * 4.;

  for (size_t i = 0; i < noCells; i++)
    {
      const valueType rkpow3 = r(i) * r(i) * r(i);
      const valueType rkp1pow3 = r(i + 1) * r(i + 1) * r(i + 1);
      m(i) = (rho(i) * pi4 * (rkp1pow3 - rkpow3)) / 3.;
    }
}
}
#endif

