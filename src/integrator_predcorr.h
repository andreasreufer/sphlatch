#ifndef BHTREE_INTEGRATOR_PREDCORR_H
#define BHTREE_INTEGRATOR_PREDCORR_H

/*
 *  integrator_predcorr.h
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
class GenericPredCorr : public GenericIntegrator {
public:
GenericPredCorr()
{
};
virtual ~GenericPredCorr()
{
};
public:

///
/// bootstrap:
///  bootstrap the predictor/corrector
///
virtual void bootstrap(const valueType& _dt) = 0;

///
/// prediction step:
///  x_p = x + 0.5*( 3*v_0 + - v_-1 )*dt
///
virtual void predict(const valueType& _dt) = 0;

///
/// correction step:
///  x_1 = x_0 + 0.5(  v_p + v_0 )*dt
///
virtual void correct(const valueType& _dt) = 0;

///
/// prepare:
///  resize internal temporaries
///
virtual void prepare() = 0;
};

///
/// 1st order scalar PredCorr integrator
///
class PredCorrScalO1 : public GenericPredCorr {
public:
PredCorrScalO1(valvectRefType _var,
               valvectRefType _dVar)
{
    var =   &_var;
   dVar =  &_dVar;

  /// take name of the variables, prefix with "o"
  /// and register it to PartManager for resizing
  PartManager.regQuantity(  oVar, "o" + PartManager.getName(  *var));
  PartManager.regQuantity( odVar, "o" + PartManager.getName( *dVar));
  /// also register to CommManager if we're gonna do parallel
#ifdef SPHLATCH_PARALLEL
  CommManager.regExchQuant(  oVar);
  CommManager.regExchQuant( odVar);
#endif
};

~PredCorrScalO1()
{
  PartManager.unRegQuantity(  oVar);
  PartManager.unRegQuantity( odVar);
};

public:
void prepare(void)
{
  valvectRefType v(*dVar);

  noParts = v.size();
  if (zero.size() != noParts)
    {
      zero.resize(noParts);
    }
}

void bootstrap(const valueType& _dt)
{
  valvectRefType   v(*dVar);
  valvectRefType  ov(odVar);
  
  ///
  /// set the old variables to the current ones, with that
  /// the predictor/corretor becomes a leapfrog for the
  /// first step
  ///
  ov = v;
}

void predict(const valueType& _dt)
{
  valvectRangeType x(*var,  rangeType(0, noParts));
  valvectRefType  ox(oVar);
  
  valvectRefType   v(*dVar);
  valvectRefType  ov(odVar);

  /// ox = x
  /// x = x + 0.5*( 3v - ov )*dt
  ox = x;
  x += ( 1.5*v - 0.5*ov )*_dt;

  /// oa = a
  /// a = 0;
  ov.swap(v);
  v = zero;
}

void correct(const valueType& _dt)
{
  valvectRangeType x(*var,  rangeType(0, noParts));
  valvectRefType  ox(oVar);
  
  valvectRefType   v(*dVar);
  valvectRefType  ov(odVar);
  
  x = ox + 0.5*_dt*( v + ov );
  //v = zero;  << does that belong to the integrator?
}
private:
valvectType oVar, odVar;
valvectPtrType var, dVar;
zerovalvectType zero;
size_t noParts;
};


///
/// 2nd order vectorial PredCorr integrator
///
class PredCorrVectO2 : public GenericPredCorr {
public:
PredCorrVectO2(matrixRefType _var,
               matrixRefType _dVar,
               matrixRefType _ddVar)
{
    var =   &_var;
   dVar =  &_dVar;
  ddVar = &_ddVar;

    oVar.resize(0, var->size2());
   odVar.resize(0, var->size2());
  oddVar.resize(0, var->size2());
  
  zero.resize(0, var->size2());

  /// take name of the variables, prefix with "o"
  /// and register it to PartManager for resizing
  PartManager.regQuantity(  oVar, "o" + PartManager.getName(  *var));
  PartManager.regQuantity( odVar, "o" + PartManager.getName( *dVar));
  PartManager.regQuantity(oddVar, "o" + PartManager.getName(*ddVar));
  /// also register to CommManager if we're gonna do parallel
#ifdef SPHLATCH_PARALLEL
  CommManager.regExchQuant(  oVar);
  CommManager.regExchQuant( odVar);
  CommManager.regExchQuant(oddVar);
#endif
};

~PredCorrVectO2()
{
  PartManager.unRegQuantity(  oVar);
  PartManager.unRegQuantity( odVar);
  PartManager.unRegQuantity(oddVar);
};

public:
void prepare(void)
{
  matrixRefType a(*ddVar);

  noParts = a.size1();
  if (zero.size1() != noParts)
    {
      zero.resize(noParts, a.size2());
    }

  /// a = 0;
  //a = zero;
}

void bootstrap(const valueType& _dt)
{
  matrixRefType  a(*ddVar);
  matrixRefType oa(oddVar);
  
  matrixRangeType v(*dVar, rangeType(0, noParts), rangeType(0, dVar->size2()));
  matrixRefType  ov(odVar);

  ///
  /// set the old variables to the current ones, with that
  /// the predictor/corretor becomes a leapfrog for the
  /// first step
  ///
  oa = a;
  ov = v;
}

void predict(const valueType& _dt)
{
  matrixRangeType x(*var,  rangeType(0, noParts), rangeType(0, var->size2()));
  matrixRefType  ox(oVar);
  
  matrixRangeType v(*dVar, rangeType(0, noParts), rangeType(0, dVar->size2()));
  matrixRefType  ov(odVar);

  matrixRefType   a(*ddVar);
  matrixRefType  oa(oddVar);

  /// ox = x
  /// x = x + 0.5*( 3v - ov )*dt
  ox = x;
  x += ( 1.5*v - 0.5*ov )*_dt;

  /// ov = v
  /// v = v + 0.5*( 3a - oa )*dt
  ov = v;
  v += ( 1.5*a - 0.5*oa )*_dt;
  
  /// oa = a
  /// a = 0;
  oa.swap(a);
  a = zero;
}

void correct(const valueType& _dt)
{
  matrixRangeType x(*var,  rangeType(0, noParts), rangeType(0, var->size2()));
  matrixRefType  ox(oVar);
  
  matrixRangeType v(*dVar, rangeType(0, noParts), rangeType(0, dVar->size2()));
  matrixRefType  ov(odVar);

  matrixRefType   a(*ddVar);
  matrixRefType  oa(oddVar);
  
  x = ox + 0.5*_dt*( v + ov );
  v = ov + 0.5*_dt*( a + oa );
  a = zero;
}
private:
matrixType oVar, odVar, oddVar;
matrixPtrType var, dVar, ddVar;
zeromatrixType zero;
size_t noParts;
};

///
/// the PredCorr meta integrator
///
class PredCorrMetaIntegrator : public GenericMetaIntegrator
{
private:
typedef std::list<GenericPredCorr*> integratorsListType;
typedef void (*voidVoidFuncPtr)(void);
typedef valueType (*valueTypeVoidFuncPtr)(void);

integratorsListType::iterator integItr, integEnd;
voidVoidFuncPtr derivFunc;
valueTypeVoidFuncPtr timingFunc;

public:
///
/// instantate the metaIntegrator with a
/// function pointer to the derivation function
///
PredCorrMetaIntegrator(void(*_deriv)(void), valueType(*_timing)(void))
{
  derivFunc = _deriv;
  timingFunc = _timing;
};

~PredCorrMetaIntegrator()
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

  integItr = integrators.begin();
  while (integItr != integEnd)
    {
      (*integItr)->prepare();
      integItr++;
    }

  derivFunc();
  
  const valueType dt = timingFunc();

  integItr = integrators.begin();
  while (integItr != integEnd)
    {
      (*integItr)->bootstrap(dt);
      integItr++;
    }

  valueRefType time(PartManager.attributes["time"]);
  time += dt;
};

///
/// this is the integration mail loop
///
void integrate(void)
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

  /// particles are freshly moved
  /// derivative
  derivFunc();

  ///
  /// timestepping belongs here
  const valueType dt = timingFunc();

  /// predict
  integItr = integrators.begin();
  while (integItr != integEnd)
    {
      (*integItr)->predict(dt);
      integItr++;
    }

  valueRefType time(PartManager.attributes["time"]);
  time += dt;
  derivFunc();
  
  /// correct
  integItr = integrators.begin();
  while (integItr != integEnd)
    {
      (*integItr)->correct(dt);
      integItr++;
    }
};

void regIntegration(valvectRefType _var,
                      valvectRefType _dVar)
{
  integrators.push_back(new PredCorrScalO1(_var, _dVar));
}

void regIntegration(matrixRefType _var,
                    matrixRefType _dVar,
                    matrixRefType _ddVar)
{
  integrators.push_back(new PredCorrVectO2(_var, _dVar, _ddVar));
}


private:
integratorsListType integrators;
};
};

#endif
