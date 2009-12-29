#ifndef SPHLATCH_EOS_IDEALGAS
#define SPHLATCH_EOS_IDEALGAS

/*
 *  eos_idealgas.cpp
 *
 *
 *  Created by Andreas Reufer on 23.12.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include <fstream>
#include <boost/lexical_cast.hpp>

#include "typedefs.h"
#include "eos_generic.cpp"

namespace sphlatch {
template<typename _partT>
class IdealGas : public EOS {
public:
   IdealGas<_partT>() {gamma = 1.4; gammaone = 0.4;}
   ~IdealGas() {}

   static IdealGas& instance();
   static IdealGas* _instance;

private:
   fType gamma, gammaone;

///
/// get the pressure & speed of sound for particle _i
///
/// common EOS interface for particle use
///
public:
   void operator()(_partT& _part)
   {
     _part.p = gammaone * ( _part.u * _part.rho );
     _part.cs = sqrt( _part.p * gamma / _part.rho );
   }

   void setGamma(const fType _gamma)
   {
     gamma = _gamma;
     gammaone = ( gamma - 1.);
   }
};

template<typename _partT>
IdealGas<_partT>* IdealGas<_partT>::_instance = NULL;

template<typename _partT>
IdealGas<_partT>& IdealGas<_partT>::instance()
{
   if (_instance == NULL)
      _instance = new IdealGas;
   return(*_instance);
}
}
#endif
