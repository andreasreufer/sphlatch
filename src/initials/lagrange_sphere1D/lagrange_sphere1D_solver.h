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

namespace sphlatch {
class LagrangeSphere1DSolver {
public:
LagrangeSphere1DSolver() :
EOS(eosType::instance())
{
  firstStep = true;
  time = 0.;

  gravConst = 6.67e-8;

  noCells = 1000;
  noEdges = noCells + 1;

  rho.resize(noCells);
  q.resize(noCells);
  u.resize(noCells);
  p.resize(noCells);
  m.resize(noCells);
  cs.resize(noCells);
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

public:
valvectType r, v, rho, q, u, p, m, cs;
idvectType mat;
valueType time;

private:
void integrate();
valueType detTimestep();

typedef sphlatch::Tillotson eosType;
eosType& EOS;

private:
size_t noCells, noEdges;

valvectType rhoOld, qOld;

valueType avAlpha, avBeta;
valueType dtp05, dtm05;
valueType gravConst;

bool firstStep;
};


void LagrangeSphere1DSolver::integrateTo(valueType _t)
{

  detTimestep();

}

void LagrangeSphere1DSolver::integrate()
{
  const valueType pi4 = M_PI*4.;
  const valueType dtc00 = 0.5(dtp05+dtm05);

  ///
  /// update velocity
  ///
  valueType mSum = 0.;
  for (size_t i = 1; i < noEdges - 1; i++)
  {
    const valueType dmk = 0.5*( m(i-1) + m(i) );
    const valueType Ak = pi4*r(i)*r(i);
    const valueType Akm10 = pi4*r(i-1)*r(i-1);
    const valueType Akp10 = pi4*r(i+1)*r(i+1);
    const valueType Akm05 = 0.5*(Akm10 + Ak);
    const valueType Akp05 = 0.5*(Akp10 + Ak);

    const valueType gradP = Ak*(p(i) - p(i-1));
    const valueType gradQ = 0.5*(q(i)*(3.*Akp05-Ak) - q(i-1)*(3.*Akm05-Ak));

    mSum += m(i-1);
    const valueType accG = gravConst*mSum / ( r(i)*r(i) );

    v(i) -= ( ( gradP + gradQ ) / dmk + accG )*dtc00;
  }
  v(0) = 0.;
  v(noEdges - 1) = 0.;

  ///
  /// update position
  ///
  r(0) = 0.;
  for (size_t i = 1; i < noEdges; i++)
  {
    r(i) += v(i)*dtp05;
  }
  
  ///
  /// calculate density
  ///
  for (size_t i = 0; i < noCells; i++)
  {
    const valueType rkpow3   = r(i)*r(i)*r(i);
    const valueType rkp1pow3 = r(i+1)*r(i+1)*r(i+1);
    rho(i) = 3.*m(i) / (pi4 *(rkp1pow3 - rkpow3));
  }

  ///
  /// update q-factor
  ///
  for (size_t i = 0; i < noCells; i++)
  {
    q(i) = 0.;
  }

  ///
  /// update internal energy
  ///
  valueType pPrime = 0.;
  for (size_t i = 0; i < noCells; i++)
  {
    const valueType Ak = pi4*r(i)*r(i);
    const valueType Akp10 = pi4*r(i+1)*r(i+1);
    const valueType Akp05 = 0.5*(Akp10 + Ak);

    ///
    /// try uPrime and guess pPrime
    ///
    const valueType dVol = (Akp10*v(i+1) - Ak*v(i))*(dtp05/m(i));
    const valueType uPrime = u(i) - p(i)*dVol;
    EOS(rho(i), uPrime, mat(i), pPrime, cs(i) );
    p(i) = 0.5*(p(i) + pPrime);
    q(i) = 0.5*(q(i) + qOld(i));

    ///
    /// AV heating
    ///
    const valueType dQ = (v(i+1)*(3.*Akp05-Akp10) - v(i)*(3.*Akp05-Ak))
                         *(dtp05/m(i));
    
    u(i) =- (p(i)*dVol + dQ);
  }

  ///
  /// update pressure
  ///
  for (size_t i = 0; i < noCells; i++)
  {
    EOS(rho(i), u(i), mat(i), p(i), cs(i) );
  }
}

valueType LagrangeSphere1DSolver::detTimestep()
{
  ///
  /// calculate current timestep dtp05
  ///

  dt = 0.;

  if (firstStep)
  {
    dtc00 = dt05;
    firstStep = false;
  }
  else
  {
    dtc00 = 0.5*(dtp05 + dtm05);
  }

  return 0.;
}

}
#endif

