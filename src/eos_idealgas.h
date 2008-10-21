#ifndef SPHLATCH_EOS_IDEALGAS
#define SPHLATCH_EOS_IDEALGAS

/*
 *  eos_idealgas.h
 *
 *
 *  Created by Andreas Reufer on 26.07.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "eos_generic.h"

namespace sphlatch {
class IdealGas : public EOS {
public:
IdealGas()
{
  loadGamma();
};

~IdealGas()
{
};

static IdealGas& instance();
static IdealGas* _instance;

///
/// get the pressure and speed of sound for particle _i
///
void operator()(const size_t _i, valueType& _P, valueType& _cs)
{
  _P = gammaone * (u(_i) * rho(_i));
  _cs = sqrt(p(_i) * gamma / rho(_i));
  //_T = uToT*u(_i);
  return;
};

public:
void loadGamma()
{
  gamma = PartManager.attributes["gamma"];
  gammaone = gamma - 1.;
#ifdef SPHLATCH_LOGGER
  Logger.stream << "ideal gas EOS with gamma " << gamma;
  Logger.flushStream();
#endif
};

private:
valueType gamma, gammaone;
};

IdealGas * IdealGas::_instance = NULL;
IdealGas& IdealGas::instance()
{
  if (_instance == NULL)
    _instance = new IdealGas;
  return *_instance;
};
}
#endif
