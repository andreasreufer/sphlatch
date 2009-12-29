#ifndef SPHLATCH_SPH_ALGORITHMS_CPP
#define SPHLATCH_SPH_ALGORITHMS_CPP

#include "typedefs.h"

namespace sphlatch {
template<typename _partT, typename _krnlT>
struct densSum
{
   _krnlT K;

   fType  rhoi;
#ifdef SPHLATCH_NONEIGH
   cType  noneighi;
#endif

   void   preSum(_partT* const _i)
   {
      rhoi = 0.;

#ifdef SPHLATCH_NONEIGH
      noneighi = 0;
#endif
   }

   void operator()(_partT* const _i,
                   const _partT* const _j,
                   const vect3dT& _rvec,
                   const fType _rr,
                   const fType _srad)
   {
      const fType r   = sqrt(_rr);
      const fType hij = 0.25 * _srad + 0.5 * (_j->h);
      const fType mj  = _j->m;

      rhoi += mj * K.value(r, hij);
#ifdef SPHLATCH_NONEIGH
      noneighi++;
#endif
   }

   void postSum(_partT* const _i)
   {
      _i->rho = rhoi;
#ifdef SPHLATCH_NONEIGH
      _i->noneigh = noneighi;
#endif
   }
};

template<typename _partT, typename _krnlT>
struct accPowSum
{
   _krnlT  K;

   void    preSum(_partT* const _i)
   {
      acci  = 0., 0., 0.;
      dudti = 0.;

      vi   = _i->vel;
      rhoi = _i->rho;
      pi   = _i->p;
      ci   = _i->cs;

      piOrhoirhoi = pi / (rhoi * rhoi);
   }

   void operator()(_partT* const _i,
                   const _partT* const _j,
                   const vect3dT& _rvec,
                   const fType _rr,
                   const fType _srad)
   {
      //FIXME: that doesn't belong here
      const fType alpha = 1.;
      const fType beta  = 2.;

      const fType r   = sqrt(_rr);
      const fType hij = 0.25 * _srad + 0.5 * (_j->h);

      const fType pj   = _j->p;
      const fType rhoj = _j->rho;
      const fType cj   = _j->cs;
      const fType mj   = _j->m;

      const vect3dT vij = vi - _j->vel;

      const fType vijrij = dot(_rvec, vij);

      fType av = 0.;

      if (vijrij < 0.)
      {
         const fType rijrij = dot(_rvec, _rvec);
         const fType rhoij  = 0.5 * (rhoi + rhoj);
         const fType cij    = 0.5 * (ci + cj);

         const fType muij = hij * vijrij / (rijrij + 0.01 * hij * hij);

         av = (-alpha * cij * muij + beta * muij * muij) / rhoij;
      }

      K.derive(r, hij, _rvec);

      const fType accTerm     = piOrhoirhoi + (pj / (rhoj * rhoj)) + av;
      const fType mjvijdivWij = mj * dot(vij, K.deriv);

      acci -= mj * accTerm * K.deriv;
      dudti += 0.5 * accTerm * mjvijdivWij;
   }

   void postSum(_partT* const _i)
   {
      _i->acc += acci;
      _i->dudt = dudti;
   }

   vect3dT vi, acci;
   fType   rhoi, pi, ci, piOrhoirhoi, dudti;
};
};
#endif
