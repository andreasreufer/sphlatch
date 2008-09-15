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

extern "C" void aneosinit_(const char *, /// materials file
                           int);         /// filename string length
extern "C" void aneos_(double *,  /// T
                       double *,  /// rho
                       double *,  /// p
                       double *,  /// E
                       double *,  /// S
                       double *,  /// c_v
                       double *,  /// dp/dt
                       double *,  /// dp/dr
                       double *,  /// fkro (Rosseland mean opacity)
                       double *,  /// cs
                       int *,     /// phase
                       int *,     /// material
                       double *,  /// fme (?)
                       double *  /// fva (?)
                       );

namespace sphlatch {
class ANEOS : public EOS {
public:
ANEOS()
{
  Logger.stream << "init ANEOS EOS with "
                << " materials";
  Logger.flushStream();

  std::string matFilename = "aneos.input";
  aneosinit_(matFilename.c_str(), matFilename.size());
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
  this->operator()(rho (_i), u (_i), mat (_i), _P, _cs);
}

///
/// get the pressure & speed of sound for given parameters
///
/// the expressions for the speed of sound ( cc && ce )
/// are copied from ParaSPH. in the hybrid regime, the square
/// roots for both cc and ce are already taken, in opposite
/// to the scheme in ParaSPH where the square root of the hybrid
/// formula is taken.
///
void operator()(const valueType _rho, const valueType _u, const identType _mat,
                valueType& _P, valueType& _cs)
{
  
};


private:
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

