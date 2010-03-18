#ifndef SPHLATCH_RPSPH_CPP
#define SPHLATCH_RPSPH_CPP

#include "typedefs.h"
namespace sphlatch {

template<typename _partT, typename _krnlT>
struct accRPSPH
{
   _krnlT  K;

   void    preSum(_partT* const _i)
   {
      acci  = 0., 0., 0.;
#ifdef SPHLATCH_VELDIV
      curDrhoDti = 0.;
#endif
      vi   = _i->vel;
      rhoi = _i->rho;
      pi   = _i->p;
      ci   = _i->cs;

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

      const fType accTerm = ( pj - pi ) / ( rhoj*rhoj ) + av;
#ifdef SPHLATCH_VELDIV
      const fType mjvijdivWij = mj * dot(vij, K.deriv);
      curDrhoDti += mjvijdivWij;
#endif

      acci -= mj * accTerm * K.deriv;
   }

   void postSum(_partT* const _i)
   {
      _i->acc += acci;
#ifdef SPHLATCH_VELDIV
      _i->divv = curDrhoDti / rhoi;
#endif
   }

   vect3dT vi, acci;
   fType   rhoi, pi, ci;
#ifdef SPHLATCH_VELDIV
   fType curDrhoDti;
#endif
};


};
#endif
