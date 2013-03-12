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
   fType   m, rc, rho, Ekin, Erot, Ethm, I, Erotclmp, Iclmp;
#ifdef SPHLATCH_GRAVITY
   fType Epot, Epotclmp;
#endif

   vect3dT posclmp, velclmp;
   fType   rclmp;
   fType   mclmp, mdisk, mimp, mesc;
   vect3dT Lclmp, Ldisk, Limp, Lesc;

#ifdef SPHLATCH_ANEOS
   fType mmat[33], mclmpmat[33], mdiskmat[33], mimpmat[33], mescmat[33];
#endif

#ifdef SPHLATCH_ANEOS_TRACKORIG
   fType mtarg[33], mproj[33];
#endif

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

      Ekin = 0.;
      Ethm = 0.;
      Erot = 0.;
      Erotclmp = 0.;
      I     = 0.;
      Iclmp = 0.;
#ifdef SPHLATCH_GRAVITY
      Epot = 0.;
      Epotclmp = 0.;
#endif

      mclmp = 0.;
      mdisk = 0.;
      mimp  = 0.;
      mesc  = 0.;

#ifdef SPHLATCH_ANEOS
      for (size_t i = 0; i < 33; i++)
      {
         mmat[i]     = 0.;
         mclmpmat[i] = 0.;
         mdiskmat[i] = 0.;
         mimpmat[i]  = 0.;
         mescmat[i]  = 0.;
      }
#endif

#ifdef SPHLATCH_ANEOS_TRACKORIG
      for (size_t i = 0; i < 33; i++)
      {
         mtarg[i] = 0.;
         mproj[i] = 0.;
      }
#endif


      rclmp = 0.;

      pos = 0., 0., 0.;
      vel = 0., 0., 0.;
      P   = 0., 0., 0.;
      L   = 0., 0., 0.;

      posclmp = 0., 0., 0.;
      velclmp = 0., 0., 0.;
      Lclmp   = 0., 0., 0.;
      Ldisk   = 0., 0., 0.;
      Limp    = 0., 0., 0.;
      Lesc    = 0., 0., 0.;

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
      m    = 0.;
      P    = 0., 0., 0.;
      Ekin = 0.;
      Ethm = 0.;
      pos  = 0., 0., 0.;
      rho  = 0.;
      nop  = 0;

#ifdef SPHLATCH_ANEOS
      for (size_t i = 0; i < 33; i++)
         mmat[i] = 0.;
#endif

#ifdef SPHLATCH_ANEOS_TRACKORIG
      for (size_t i = 0; i < 33; i++)
      {
         mtarg[i] = 0.;
         mproj[i] = 0.;
      }
#endif

      for (pItr = pList.begin(); pItr != pList.end(); pItr++)
      {
         const fType   mi = (*pItr)->m;
         const vect3dT vi = (*pItr)->vel;
         m    += mi;
         P    += mi * vi;
         Ekin += mi * dot(vi, vi) * 0.5;
         Ethm += mi * (*pItr)->u;
         pos  += mi * (*pItr)->pos;
         rho  += mi * (*pItr)->rho;

#ifdef SPHLATCH_ANEOS
         mmat[(*pItr)->mat] += mi;
#endif

#ifdef SPHLATCH_ANEOS_TRACKORIG
         if ((*pItr)->id > 2000000)
            mproj[(*pItr)->mat] += mi;
         else
            mtarg[(*pItr)->mat] += mi;
#endif
         nop++;
      }
      pos /= m;
      rho /= m;
      vel  = P / m;

      // FIXME: senseless large radii for gaseous objects
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

      L    = 0., 0., 0.;
      I    = 0.;
      Erot = 0.;
      for (pItr = pList.begin(); pItr != pList.end(); pItr++)
      {
         const fType   mi   = (*pItr)->m;
         const vect3dT rpos = (*pItr)->pos - pos;
         const vect3dT rvel = (*pItr)->vel - vel;

         const vect3dT rXv    = cross(rpos, rvel);
         const fType   rXvrXv = dot(rXv, rXv);
         const fType   rr     = dot(rpos, rpos);

         L    += mi * rXv;
         I    += mi * rr;
         Erot += (0.5 * mi * rXvrXv / rr);
      }
   }

   fvectT getMassFractions(const size_t _maxMat)
   {
      const size_t maxId = _maxMat + 1;
      fvectT       mfrac(maxId);

      for (size_t i = 0; i < maxId; i++)
         mfrac(i) = 0.;

      fType m = 0.;
      typename partPtrLT::const_iterator pItr;
      for (pItr = pList.begin(); pItr != pList.end(); pItr++)
      {
         const fType mi = (*pItr)->m;
         mfrac((*pItr)->mat) += mi;
         m += mi;
      }
      mfrac /= m;
      return(mfrac);
   }

   void getCentralBodyOrbits(const fType _rhoMin, const fType _G)
   {
      lowerPot potSorter;

      pList.sort(potSorter);

      // if at least 1% of all encountered particles have
      // had rho < rho_min, stop
      typename partPtrLT::const_iterator pCBend;
      const size_t lrhopmin = round(0.01 * nop);
      size_t       lrhop    = 0;
      for (pCBend = pList.begin(); pCBend != pList.end(); pCBend++)
      {
         if ((*pCBend)->rho < _rhoMin)
            lrhop++;
         if (lrhop > lrhopmin)
            break;
      }

      // pCBend now points to the first few particles with rho < rho_min, note
      // this does not necessarily mark the surface of the body
      vect3dT comcb;
      comcb = 0., 0., 0.;
      fType mcb = 0.;

      typename partPtrLT::const_iterator pItr;
      for (pItr = pList.begin(); pItr != pCBend; pItr++)
      {
         const fType mi = (*pItr)->m;
         mcb   += mi;
         comcb += mi * (*pItr)->pos;
      }
      comcb /= mcb;

      // take a first guess for the clump radius
      if (pItr == pCBend)
         pItr--;
      const vect3dT rcbguessv = (*pItr)->pos - comcb;
      const fType   rcbguess  = sqrt(dot(rcbguessv, rcbguessv));

      // now try to find the shell representing the surface of the body
      fType       rhoTol    = 0.05; // start with 5% tolerance level
      const fType shellFrac = 0.05;

      fType  rcb  = 0.;
      size_t nocb = 0;
      while (rhoTol < 1.0)
      {
         rcb  = 0.;
         nocb = 0;
         for (pItr = pList.begin(); pItr != pList.end(); pItr++)
         {
            const vect3dT rvec = (*pItr)->pos - comcb;
            const fType   r    = sqrt(dot(rvec, rvec));
            const fType   rhoi = (*pItr)->rho;

            if ((r < 3. * rcbguess) and
                (fabs((rhoi - _rhoMin) / _rhoMin) < rhoTol))
            {
               rcb += r;
               nocb++;
            }
         }

         if (nocb > shellFrac * nop)
            break;
         else
            rhoTol *= 1.5;
      }

      if (nocb == 0)
         return;

      rcb /= static_cast<fType>(nocb);

      // recalculate central body mass, center of mass and velocity
      vect3dT vcb;
      vect3dT oldcomcb = comcb;
      mcb   = 0.;
      comcb = 0., 0., 0.;
      vcb   = 0., 0., 0.;

      for (pItr = pList.begin(); pItr != pList.end(); pItr++)
      {
         const vect3dT rvec = (*pItr)->pos - oldcomcb;
         const fType   r    = sqrt(dot(rvec, rvec));
         const fType   mi   = (*pItr)->m;

         if (r < rcb)
         {
            (*pItr)->orbit = ORBITCLUMP;
            mcb           += mi;
            vcb           += mi * (*pItr)->vel;
            comcb         += mi * (*pItr)->pos;
         }
         else
            (*pItr)->orbit = ORBITNOTSET;
      }

      comcb  /= mcb;
      vcb    /= mcb;
      posclmp = comcb;
      velclmp = vcb;
      rclmp   = rcb;

      // now determine eccentricites and perihelon distances
      for (pItr = pList.begin(); pItr != pList.end(); pItr++)
      {
         const vect3dT r0v = (*pItr)->pos - comcb;
         const vect3dT v0v = (*pItr)->vel - vcb;

         const fType mi = (*pItr)->m;

         // calculate the angular momentum according to the center of mass of the clump
         const vect3dT Li = mi * cross(r0v, v0v);
         const fType Ii = mi * dot( r0v, r0v );

#ifdef SPHLATCH_ANEOS
         const iType mati = (*pItr)->mat;
#endif

         if ((*pItr)->orbit == ORBITNOTSET)
         {
            const fType r0 = sqrt(dot(r0v, r0v));
            const fType v0 = sqrt(dot(v0v, v0v));

            const fType K  = (mcb + (*pItr)->m) * _G;
            const fType k1 = r0 * v0 * v0 / K;

            const fType beta0  = asin(dot(r0v, v0v) / (r0 * v0));
            const fType theta0 = atan2(k1 * sin(beta0) * cos(beta0),
                                       (k1 * cos(beta0) * cos(beta0) - 1.));

            const vect3dT ev = (cross(v0v, cross(r0v, v0v)) / K) - (r0v / r0);
            const fType   e  = sqrt(dot(ev, ev));
            const fType   rp = r0 * ((1. + e * cos(theta0)) / (1. + e));


            (*pItr)->ecc = e;

            if (e > 1.)
            {
               (*pItr)->orbit = ORBITESCAPE;
               mesc          += mi;
               Lesc          += Li;
#ifdef SPHLATCH_ANEOS
               mescmat[mati] += mi;
#endif
            }
            else
            {
               if (rp < rcb)
               {
                  (*pItr)->orbit = ORBITREIMPACT;
                  mimp          += mi;
                  Limp          += Li;
#ifdef SPHLATCH_ANEOS
                  mimpmat[mati] += mi;
#endif
               }
               else
               {
                  (*pItr)->orbit = ORBITDISK;
                  (*pItr)->a     = rp / ( 1. - e ); 
                  mdisk         += mi;
                  Ldisk         += Li;
#ifdef SPHLATCH_ANEOS
                  mdiskmat[mati] += mi;
#endif
               }
            }
         }
         else
         {
            mclmp += mi;
            Lclmp += Li;
            Iclmp += Ii;
#ifdef SPHLATCH_ANEOS
            mclmpmat[mati] += mi;
#endif
         }
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

      vars.push_back(storeVar(Ekin, "Ekin"));
      vars.push_back(storeVar(Erot, "Erot"));
      vars.push_back(storeVar(Erotclmp, "Erotclmp"));
      vars.push_back(storeVar(Ethm, "Ethm"));
      vars.push_back(storeVar(I, "I"));
      vars.push_back(storeVar(Iclmp, "Iclmp"));
#ifdef SPHLATCH_GRAVITY
      vars.push_back(storeVar(Epot, "Epot"));
      vars.push_back(storeVar(Epotclmp, "Epotclmp"));
#endif

      vars.push_back(storeVar(id, "id"));
      vars.push_back(storeVar(nop, "nop"));

#ifdef SPHLATCH_ANEOS
      vars.push_back(storeVar(mmat, "mmat", 33));
#endif

      // only available, if getCentralBodyOrbits() was executed
      vars.push_back(storeVar(posclmp, "posclmp"));
      vars.push_back(storeVar(velclmp, "velclmp"));
      vars.push_back(storeVar(mclmp, "mclmp"));
      vars.push_back(storeVar(mdisk, "mdisk"));
      vars.push_back(storeVar(mimp, "mimp"));
      vars.push_back(storeVar(mesc, "mesc"));
      vars.push_back(storeVar(rclmp, "rclmp"));
#ifdef SPHLATCH_ANEOS
      vars.push_back(storeVar(mclmpmat, "mclmpmat", 33));
      vars.push_back(storeVar(mdiskmat, "mdiskmat", 33));
      vars.push_back(storeVar(mimpmat, "mimpmat", 33));
      vars.push_back(storeVar(mescmat, "mescmat", 33));
#endif

#ifdef SPHLATCH_ANEOS_TRACKORIG
      vars.push_back(storeVar(mtarg, "mtarg", 33));
      vars.push_back(storeVar(mproj, "mproj", 33));
#endif

      vars.push_back(storeVar(Lclmp, "Lclmp"));
      vars.push_back(storeVar(Ldisk, "Ldisk"));
      vars.push_back(storeVar(Limp, "Limp"));
      vars.push_back(storeVar(Lesc, "Lesc"));
      return(vars);
   }

public:
   class lowerPot
   {
public:
      bool operator()(_partT* _p1, _partT* _p2)
      {
         return(_p1->pot < _p2->pot);
      }
   };

private:
   partPtrLT pList;
};
}

#endif
