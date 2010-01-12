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
      const fType hij = 0.25 * _srad + 0.5 * (_j->h); // 0.5*_srad == hi
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
#ifdef SPHLATCH_TIMEDEP_ENERGY
      dudti = 0.;
#endif
#ifdef SPHLATCH_VELDIV
      curDrhoDti = 0.;
#endif
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
#ifdef SPHLATCH_VELDIV
      const fType mjvijdivWij = mj * dot(vij, K.deriv);
      curDrhoDti += mjvijdivWij;
#endif

      acci -= mj * accTerm * K.deriv;
#ifdef SPHLATCH_TIMEDEP_ENERGY
      dudti += 0.5 * accTerm * mjvijdivWij;
#endif
   }

   void postSum(_partT* const _i)
   {
      _i->acc += acci;
#ifdef SPHLATCH_VELDIV
      _i->divv = curDrhoDti / rhoi;
#endif
#ifdef SPHLATCH_TIMEDEP_ENERGY
      _i->dudt = dudti;
#endif
   }

   vect3dT vi, acci;
   fType   rhoi, pi, ci, piOrhoirhoi;
#ifdef SPHLATCH_TIMEDEP_ENERGY
   fType dudti;
#endif
#ifdef SPHLATCH_VELDIV
   fType curDrhoDti;
#endif
};

template<typename _partT>
void setDhDt(_partT& _i)
{
  const fType noNeighMin = (2./3.)*_i.noneighOpt;
  const fType noNeighMax = (5./3.)*_i.noneighOpt;

  const fType k1 = 0.5 * (1 + tanh((_i.noneigh - noNeighMin) / -5.));
  //const fType k1 = 0.;
  const fType k3 = 0.5 * (1 + tanh((_i.noneigh - noNeighMax) / 5.));
  const fType k2 = 1. - k1 - k3;

  _i.dhdt = (k1 * _i.divvmax - k3 * _i.divvmax - k2 * static_cast<fType>(1. / 3.) * _i.divv) * _i.h;
};

};
#endif
