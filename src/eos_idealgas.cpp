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
   IdealGas<_partT>() :
      NA(6.02214179e23)
   {
      R = 9.6485e11;          // erg/eV/mol ( 1eV ^= 11604.505 K )

      k  = 1.3806504e-16;     // erg/K
      hp = 6.62606896e-27;    // erg s

      setGamma(1.4);
      setMolarmass(1.007977); // g/mol for hydrogen
   }

   ~IdealGas() { }

   static IdealGas& instance();

   static IdealGas* _instance;

private:
   fType       gamma, gammaone, M, mp, k1;
   fType       R, k, hp;
   const fType NA;

///
/// get the pressure & speed of sound for particle _i
///
/// common EOS interface for particle use
///
public:
   void operator()(_partT& _part)
   {
     // the Sakure-Tetrode equation
#ifdef SPHLATCH_TIMEDEP_ENTROPY
      _part.u = (1. / mp) *
                pow((_part.rho/mp) * exp((_part.S * mp - k1) / k), 2. / 3.);
#else
      _part.S = (1. / mp) *
                (k * log((mp / _part.rho) * pow(_part.u * mp, 3. / 2.)) + k1);
#endif
      _part.T  = (2. / 3.) * _part.u * M / R;
      _part.p  = gammaone * (_part.u * _part.rho);
      _part.cs = sqrt(_part.p * gamma / _part.rho);
   }

   void setGamma(const fType _gamma)
   {
      gamma    = _gamma;
      gammaone = (gamma - 1.);
   }

   void setMolarmass(const fType _M)
   {
      M  = _M;
      mp = M / NA;

      k1 = (3. / 2) * k * ((5. / 3.) + log(4. * M_PI * mp / (3. * hp * hp)));
   }
};

template<typename _partT>
IdealGas<_partT> * IdealGas<_partT>::_instance = NULL;

template<typename _partT>
IdealGas<_partT>& IdealGas<_partT>::instance()
{
   if (_instance == NULL)
      _instance = new IdealGas;
   return(*_instance);
}
}
#endif
