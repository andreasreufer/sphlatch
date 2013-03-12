#ifndef SPHLATCH_EOS_SUPER
#define SPHLATCH_EOS_SUPER

/*
 *  eos_super.h
 *
 *
 *  Created by Andreas Reufer on 06.09.10.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include <fstream>

#include "typedefs.h"
#include "eos_generic.cpp"

#ifdef SPHLATCH_ANEOS
#include "eos_aneos.cpp"
#endif

#include "eos_idealgas.cpp"


namespace sphlatch {

template<typename _partT>
class SuperEOS : public EOS {
public:
   SuperEOS();
   ~SuperEOS();

   static SuperEOS& instance();
   static SuperEOS* _instance;

#ifdef SPHLATCH_ANEOS
   typedef ANEOS<_partT>       aneosT;
#endif
   typedef IdealGas<_partT>    idealgasT;

public:
   void operator()(_partT& _part);

#ifdef SPHLATCH_ANEOS
   aneosT aneos;
#endif
   idealgasT idealgas;

private:
   fType nan;

};


template<typename _partT>
SuperEOS<_partT>::SuperEOS():
#ifdef SPHLATCH_ANEOS
  aneos(aneosT::instance()),
#endif
  idealgas(idealgasT::instance())
{
}

template<typename _partT>
SuperEOS<_partT>::~SuperEOS() { }

template<typename _partT>
SuperEOS<_partT> * SuperEOS<_partT>::_instance = NULL;

template<typename _partT>
SuperEOS<_partT>& SuperEOS<_partT>::instance()
{
   if (_instance == NULL)
      _instance = new SuperEOS;
  
   return(*_instance);
}

///
/// common EOS interface for particle use
///
template<typename _partT>
void SuperEOS<_partT>::operator()(_partT& _part)
{
  switch (_part.mat)
  {
    case 0:
      idealgas(_part);
      break;
#ifdef SPHLATCH_ANEOS
    default:
      aneos(_part);
      break;
#endif
  }
}

}
#endif
