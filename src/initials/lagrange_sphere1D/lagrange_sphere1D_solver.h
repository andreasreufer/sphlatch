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
#include "eos_generic.h"

namespace sphlatch {
class LagrangeSphere1DSolver {
public:
LagrangeSphere1DSolver()
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
valvectType r, v, rho, q, u, p, m;
idvectType mat;
valueType time;

private:
void integrate();
void detTimestep();

private:
size_t noCells, noEdges;

valvectType rhoOld, qOld;

valueType avAlpha, avBeta;
valueType dt, dtOld, dtMiddle;
valueType gravConst;

bool firstStep;
};


void LagrangeSphere1DSolver::integrateTo(valueType _t)
{


}

void LagrangeSphere1DSolver::integrate()
{
  detTimestep();

  const valueType pi4 = M_PI*4.;

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

    v(i) -= ( ( gradP + gradQ ) / dmk + accG )*dtMiddle;
  }
  v(0) = 0.;
  v(noEdges - 1) = 0.;

  ///
  /// update position
  ///
  r(0) = 0.;
  for (size_t i = 1; i < noEdges; i++)
  {
    r(i) += v(i)*dt;
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
  for (size_t i = 0; i < noCells; i++)
  {
    const valueType dmk = 0.5*( m(i-1) + m(i) );
    const valueType Ak = pi4*r(i)*r(i);
    const valueType Akp10 = pi4*r(i+1)*r(i+1);
    const valueType Akp05 = 0.5*(Akp10 + Ak);
    
  }

  ///
  /// update pressure
  ///
  for (size_t i = 0; i < noCells; i++)
  {
  }
}

void LagrangeSphere1DSolver::detTimestep()
{
  ///
  /// calculate current timestep dt
  ///

  dt = 0.;

  if (firstStep)
  {
    dtMiddle = dt;
    firstStep = false;
  }
  else
  {
    dtMiddle = 0.5*(dt + dtOld);
  }
}

}
#endif

