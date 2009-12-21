#ifndef SPHLATCH_EOS_MIEGRUENEISEN
#define SPHLATCH_EOS_MIEGRUENEISEN

/*
 *  eos_miegrueneisen.h
 *
 *
 *  Created by Andreas Reufer on 26.07.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "eos_generic.h"

namespace sphlatch {
class MieGrueneisen : public EOS {
public:
MieGrueneisen()
{
  loadGamma();
};

~MieGrueneisen()
{
};

static MieGrueneisen& instance();
static MieGrueneisen* _instance;

///
/// get the pressure and speed of sound for particle _i
///
void operator()(const size_t& _i, fType& _P, fType& _cs)
{
  _P = gammaone * ( u(_i) * rho(_i) );
  _cs = sqrt( p(_i)*gamma / rho(_i) );
  //_T = uToT*u(_i);
  return;
};

public:
void loadGamma()
{
  gamma = PartManager.attributes["gamma"];
  gammaone = gamma - 1.;
};

private:
fType gamma, gammaone;
};

MieGrueneisen* MieGrueneisen::_instance = NULL;
MieGrueneisen& MieGrueneisen::instance()
{
  if (_instance == NULL)
    _instance = new MieGrueneisen;
  return *_instance;
};
}
#endif
