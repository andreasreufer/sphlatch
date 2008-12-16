#ifndef BHTREE_INTEGRATOR_VERLET_H
#define BHTREE_INTEGRATOR_VERLET_H

/*
 *  integrator_verlet.h
 *
 *
 *  Created by Andreas Reufer on 25.05.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "integrator_generic.h"

namespace sphlatch {
///
/// templates hierarchies with non-template base
/// class don't work, so we use non-templated functions here
///
/// virtual functions are a performance sin, but here we deal with
/// veeeery big loops, so that's not really a problem
///
class GenericVerlet : public GenericIntegrator {
public:
GenericVerlet()
{
};
virtual ~GenericVerlet()
{
};
public:
virtual void prepare() = 0;

///
/// drift changes the position:
///  x_1 = x_0 + v_0*dt + 0.5*a_0*dt^2
///
virtual void drift(const fType& _dt) = 0;

///
/// kick changes the velocity:
///  v_1 = v_0 + 0.5*a_0*dt + 0.5*a_1*dt;
///
virtual void kick(const fType& _dt) = 0;
};

///
/// 2nd order vectorial Verlet integrator
///
class VerletVectO2 : public GenericVerlet {
public:
VerletVectO2(matrixRefType _var,
             matrixRefType _dVar,
             matrixRefType _ddVar)
{
  var = &_var;
  dVar = &_dVar;
  ddVar = &_ddVar;

  oddVar.resize(0, var->size2());
  zero.resize(0, var->size2());

  /// take name of the second derivative variable, prefix with "o"
  /// and register it to PartManager for resizing
  PartManager.regQuantity(oddVar, "o" + PartManager.getName(*ddVar));
  /// also register to CommManager if we're gonna do parallel
#ifdef SPHLATCH_PARALLEL
  CommManager.regExchQuant(oddVar);
#endif
};

~VerletVectO2()
{
  PartManager.unRegQuantity(oddVar);
};

public:
void prepare()
{
  matrixRefType a(*ddVar);

  noParts = a.size1();
  zero.resize(noParts, a.size2());

  /// a = 0;
  a = zero;
}

void drift(const fType& _dt)
{
  matrixRangeType x(*var, rangeType(0, noParts), rangeType(0, var->size2()));
  matrixRangeType v(*dVar, rangeType(0, noParts), rangeType(0, dVar->size2()));

  matrixRefType a(*ddVar);
  matrixRefType oa(oddVar);

  /// move
  const fType dtSquare = _dt * _dt;

  x += (v * _dt) + 0.5 * (a * dtSquare);

  /// oa = a
  oa.swap(a);
}

void kick(const fType& _dt)
{
  matrixRangeType v(*dVar, rangeType(0, noParts), rangeType(0, dVar->size2()));
  matrixRefType a(*ddVar);
  matrixRefType oa(oddVar);

  /// boost
  v += 0.5 * (a * _dt + oa * _dt);
}
private:
matrixType oddVar;
matrixPtrType var, dVar, ddVar;
zeromatrixType zero;
size_t noParts;
};

///
/// 2nd order scalar Verlet integrator
///
class VerletScalO2 : public GenericVerlet {
public:
VerletScalO2(valvectRefType _var,
             valvectRefType _dVar,
             valvectRefType _ddVar)
{
  var = &_var;
  dVar = &_dVar;
  ddVar = &_ddVar;

  /// take name of the second derivative variable, prefix with "o"
  /// and register it to PartManager for resizing
  PartManager.regQuantity(oddVar, "o" + PartManager.getName(*ddVar));
  /// also register to CommManager if we're gonna do parallel
#ifdef SPHLATCH_PARALLEL
  CommManager.regExchQuant(oddVar);
#endif
};

~VerletScalO2()
{
  PartManager.unRegQuantity(oddVar);
};

public:
void prepare()
{
  valvectRefType a(*ddVar);

  noParts = a.size();
  zero.resize(noParts);

  /// a = 0;
  a = zero;
}

void drift(const fType& _dt)
{
  valvectRangeType x(*var, rangeType(0, noParts));
  valvectRangeType v(*dVar, rangeType(0, noParts));

  valvectRefType a(*ddVar);
  valvectRefType oa(oddVar);

  /// move
  const fType dtSquare = _dt * _dt;

  x += v * _dt + 0.5 * a * (dtSquare);

  /// oa = a
  oa.swap(a);
}

void kick(const fType& _dt)
{
  valvectRangeType v(*dVar, rangeType(0, noParts));
  valvectRefType a(*ddVar);
  valvectRefType oa(oddVar);

  /// boost
  v += 0.5 * (a * _dt + oa * _dt);
}
private:
valvectType oddVar;
valvectPtrType var, dVar, ddVar;
zerovalvectType zero;
size_t noParts;
};


///
/// the Verlet meta integrator
///
class VerletMetaIntegrator : public GenericMetaIntegrator
{
private:
typedef std::list<GenericVerlet*> integratorsListType;
typedef void (*voidVoidFuncPtr)(void);
typedef fType (*fTypeVoidFuncPtr)(void);

integratorsListType::iterator integItr, integEnd;
voidVoidFuncPtr derivFunc;
fTypeVoidFuncPtr timingFunc;
fType lastDt;

public:
VerletMetaIntegrator(void(*_deriv)(void), fType(*_timing)(void))
{
  derivFunc = _deriv;
  timingFunc = _timing;
};

~VerletMetaIntegrator()
{
  /// delete the integrators on meta integrator destruction
  integItr = integrators.begin();
  integEnd = integrators.end();

  while (integItr != integEnd)
    {
      delete * integItr;
      integItr++;
    }
};
public:
void bootstrap()
{
  integEnd = integrators.end();

  //PartManager.attributes["time"];

  integItr = integrators.begin();
  while (integItr != integEnd)
    {
      (*integItr)->prepare();
      integItr++;
    }

  derivFunc();
  
  const fType dt = timingFunc();
  
  integItr = integrators.begin();
  while (integItr != integEnd)
    {
      (*integItr)->drift(dt);
      integItr++;
    }
  
  PartManager.step++;
  fRefType time( PartManager.attributes["time"] );
  time += dt;
};

///
/// this is the integration mail loop
///
void integrate()
{
  integItr = integrators.begin();
  integEnd = integrators.end();

  ///
  /// resize the containers again
  ///
  integItr = integrators.begin();
  while (integItr != integEnd)
    {
      (*integItr)->prepare();
      integItr++;
    }
  PartManager.substep = 0;

  /// particles are freshly moved
  /// calculate derivative
  derivFunc();

  const fType dt = timingFunc();

  /// kick
  integItr = integrators.begin();
  while (integItr != integEnd)
    {
      (*integItr)->kick(dt);
      integItr++;
    }

  /// drift
  integItr = integrators.begin();
  while (integItr != integEnd)
    {
      (*integItr)->drift(dt);
      integItr++;
    }
  
  PartManager.step++;
  fRefType time( PartManager.attributes["time"] );
  time += dt;
};

void regIntegration(valvectRefType _var,
                    valvectRefType _dVar,
                    valvectRefType _ddVar)
{
  integrators.push_back(new VerletScalO2(_var, _dVar, _ddVar));
}

void regIntegration(matrixRefType _var,
                    matrixRefType _dVar,
                    matrixRefType _ddVar)
{
  integrators.push_back(new VerletVectO2(_var, _dVar, _ddVar));
}


private:
integratorsListType integrators;
};
};

#endif
