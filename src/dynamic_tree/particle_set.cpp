#ifndef SPHLATCH_PARTICLE_SET_CPP
#define SPHLATCH_PARTICLE_SET_CPP

#include "particle_set.h"

namespace sphlatch {

template<typename _partT>
_partT & ParticleSet<_partT>::operator[](const size_t _i)
{
   return(parts[_i]);
}

template<typename _partT>
void ParticleSet<_partT>::resize(const size_t _i)
{
  parts.resize(_i);
}

};

#endif
