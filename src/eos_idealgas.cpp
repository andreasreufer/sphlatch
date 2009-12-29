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
   IdealGas<_partT>()
   {}

   ~IdealGas()
   { }

   static IdealGas& instance();
   static IdealGas* _instance;

private:
   fType nan;

///
/// get the pressure & speed of sound for particle _i
///
/// common EOS interface for particle use
///
public:
   void operator()(_partT& _part)
   {
      //(*this)(_part.rho, _part.u, _part.mat, _part.p, _part.cs, _part.T,
      //        _part.phase);
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
