#ifndef MISC_PHYSICS_CPP
#define MISC_PHYSICS_CPP

#include "typedefs.h"
#include "particle_set.cpp"

namespace sphlatch {
template<typename _partT>
void addToCOM(ParticleSet<_partT>& _parts, vect3dT& _pos, vect3dT& _vel,
              fType& _m)
{
   const size_t nop = _parts.getNop();

   fType   oldm   = _m;
   vect3dT oldpos = _pos;
   vect3dT oldvel = _vel;

   _m   = 0.;
   _pos = 0., 0., 0.;
   _vel = 0., 0., 0.;

   for (size_t i = 0; i < nop; i++)
   {
      const fType mi = _parts[i].m;
      _m   += mi;
      _pos += mi * _parts[i].pos;
      _vel += mi * _parts[i].vel;
   }

   _m   += oldm;
   _pos += oldm * oldpos;
   _vel += oldm * oldpos;

   _pos /= _m;
   _vel /= _m;
}
};

#endif
