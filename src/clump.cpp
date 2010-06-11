#ifndef SPHLATCH_CLUMP_CPP
#define SPHLATCH_CLUMP_CPP

/*
 *  clump_finder.cpp
 *
 *
 *  Created by Andreas Reufer on 10.12.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include <list>

#include "typedefs.h"
#include "io_particle.h"

namespace sphlatch {
template<typename _partT>
class Clump : public IOPart
{
public:
   typedef std::list<_partT*>   partPtrLT;
   vect3dT pos;
   vect3dT vel;
   vect3dT L;
   vect3dT P;
   fType   m, rc, rho;

   cType id, nop;

   Clump()
   {
      clear();
   }

   ~Clump() { }

   void clear()
   {
      id  = -1;
      m   = 0.;
      rc  = 0.;
      rho = 0.;

      pos = 0., 0., 0.;
      vel = 0., 0., 0.;
      P   = 0., 0., 0.;
      L   = 0., 0., 0.;

      nop = 0;

      pList.clear();
   }

   void addParticle(_partT* _pPtr)
   {
      pList.push_back(_pPtr);
      const fType mo = m;
      const fType mi = _pPtr->m;

      nop++;

      m  += mi;
      P  += mi * (_pPtr->vel);
      pos = (mo * pos + mi * (_pPtr->pos)) / m;
      vel = P / m;
      rho = (mo * rho + mi * (_pPtr->rho)) / m;

      rc = pow(0.75 * m / (M_PI * rho), 1. / 3.);
   }

   void addParticles(const partPtrLT _pList)
   {
      pList = _pList;
      calcProperties();
   }

   bool operator<(const Clump& _rhs)
   {
      return(m < _rhs.m);
   }

   fType totEnergy(const _partT* _pPtr, const fType _G)
   {
      const vect3dT rpos = _pPtr->pos - pos;
      const vect3dT rvel = _pPtr->vel - vel;

      const fType r  = sqrt(dot(rpos, rpos));
      const fType vv = dot(rvel, rvel);

      return(0.5 * vv - (_G * m / r));
   }

   void calcProperties()
   {
      typename partPtrLT::const_iterator pItr;
      m   = 0.;
      pos = 0., 0., 0.;
      rho = 0.;
      nop = 0;

      for (pItr = pList.begin(); pItr != pList.end(); pItr++)
      {
         const fType mi = (*pItr)->m;
         m   += mi;
         P   += mi * (*pItr)->vel;
         pos += mi * (*pItr)->pos;
         rho += mi * (*pItr)->rho;
         nop++;
      }
      pos /= m;
      rho /= m;
      vel  = P / m;

      rc = pow(0.75 * m / (M_PI * rho), 1. / 3.);
   }

   void assignID(const cType _id)
   {
      typename partPtrLT::iterator pItr;
      for (pItr = pList.begin(); pItr != pList.end(); pItr++)
         (*pItr)->clumpid = _id;
      id = _id;
   }

   void calcAngularMom()
   {
      typename partPtrLT::const_iterator pItr;

      L = 0., 0., 0.;
      for (pItr = pList.begin(); pItr != pList.end(); pItr++)
      {
         const fType   mi   = (*pItr)->m;
         const vect3dT rpos = (*pItr)->pos - pos;
         const vect3dT rvel = (*pItr)->vel - vel;
         L += mi * cross(rpos, rvel);
      }
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(L, "L"));
      vars.push_back(storeVar(P, "P"));

      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(rc, "rc"));
      vars.push_back(storeVar(rho, "rho"));

      vars.push_back(storeVar(id, "id"));
      vars.push_back(storeVar(nop, "nop"));

      return(vars);
   }

private:
   partPtrLT pList;
};
}

#endif
