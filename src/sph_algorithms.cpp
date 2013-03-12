#ifndef SPHLATCH_SPH_ALGORITHMS_CPP
#define SPHLATCH_SPH_ALGORITHMS_CPP

#include "typedefs.h"
namespace sphlatch {
template<typename _partT, typename _krnlT>
struct densSum
{
   _krnlT K;
#ifdef SPHLATCH_MISCIBLE
   fType  deltai;
#else
   fType  rhoi;
#endif
   void   preSum(_partT* const _i)
   {
#ifdef SPHLATCH_MISCIBLE
      deltai = 0.;
#else
      rhoi = 0.;
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

#ifdef SPHLATCH_MISCIBLE
      deltai += K.value(r, hij);
#else
      const fType mj = _j->m;
      rhoi += mj * K.value(r, hij);
#endif
   }

   void postSum(_partT* const _i)
   {
#ifdef SPHLATCH_MISCIBLE
      _i->delta = deltai;
      _i->rho   = (_i->m) * deltai;
#else
      _i->rho = rhoi;
#endif
   }
};


template<typename _partT, typename _krnlT>
struct accPowSum
{
   _krnlT  K;

   void    preSum(_partT* const _i)
   {
      acci = 0., 0., 0.;

#ifdef SPHLATCH_TIMEDEP_ENERGY
      dudti = 0.;
#endif
#ifdef SPHLATCH_TRACK_UAV
      dudtavi = 0.;
#endif
#ifdef SPHLATCH_VELDIV
      divvi = 0.;
#endif
#ifdef SPHLATCH_INTEGRATERHO
      drhodti = 0.;
#endif
#ifdef SPHLATCH_NONEIGH
      noneighi = 0;
#endif

      vi   = _i->vel;
      rhoi = _i->rho;
      pi   = _i->p;
      ci   = _i->cs;
#ifdef SPHLATCH_MISCIBLE
      dlti = _i->delta;
      mi   = _i->m;
#endif

#ifdef SPHLATCH_MISCIBLE
      piOdltidlti = pi / (dlti * dlti);
#else
      piOrhoirhoi = pi / (rhoi * rhoi);
#endif
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
      const fType cj   = _j->cs;

#ifdef SPHLATCH_MISCIBLE
      const fType dltj = _j->delta;
#else
      const fType rhoj = _j->rho;
      const fType mj = _j->m;
#endif

      const vect3dT vij    = vi - _j->vel;
      const fType   vijrij = dot(_rvec, vij);

      fType av = 0.;

      if (vijrij < 0.)
      {
         const fType rijrij = dot(_rvec, _rvec);
         const fType cij    = 0.5 * (ci + cj);

         const fType muij = hij * vijrij / (rijrij + 0.01 * hij * hij);

#ifdef SPHLATCH_MISCIBLE
         const fType dltij = 0.5 * (dlti + dltj);
         const fType mij   = 0.5 * (mi + _j->m);
         av = (-alpha * cij * muij + beta * muij * muij) * (mij / dltij);
#else
         const fType rhoij = 0.5 * (rhoi + rhoj);
         av = (-alpha * cij * muij + beta * muij * muij) / rhoij;
#endif
      }
      K.derive(r, hij, _rvec);

#ifdef SPHLATCH_MISCIBLE
      const fType accTerm = piOdltidlti + (pj / (dltj * dltj)) + av;
      //const fType accTerm = av;
 #ifdef SPHLATCH_VELDIV
      const fType vijdivWij = dot(vij, K.deriv);
      divvi -= vijdivWij;
 #endif
      acci -= accTerm * K.deriv;
 #ifdef SPHLATCH_TIMEDEP_ENERGY
      // old symmetric version
      //dudti += 0.5 * accTerm * vijdivWij;

      // dervied from cont. eq. and 1st thermodynamic law
      dudti += vijdivWij;
  #ifdef SPHLATCH_TRACK_UAV
      dudtavi += 0.5* av * vijdivWij;
  #endif
 #endif

#else // not miscible
      const fType accTerm = piOrhoirhoi + (pj / (rhoj * rhoj)) + av;
 #ifdef SPHLATCH_VELDIV
      const fType mjvijdivWij = mj * dot(vij, K.deriv);
      divvi -= mjvijdivWij / rhoj;
 #endif
      acci -= mj * accTerm * K.deriv;
 #ifdef SPHLATCH_TIMEDEP_ENERGY
      dudti += 0.5 * accTerm * mjvijdivWij;
  #ifdef SPHLATCH_TRACK_UAV
      dudtavi += 0.5 * av * mjvijdivWij;
  #endif
 #endif
#endif

#ifdef SPHLATCH_NONEIGH
      noneighi++;
#endif
   }

   void postSum(_partT* const _i)
   {
#ifdef SPHLATCH_MISCIBLE
      _i->acc += (acci / mi);
 #ifdef SPHLATCH_TRACK_ACCP
      _i->accp = (acci / mi);
 #endif
 #ifdef SPHLATCH_VELDIV
      _i->divv = (divvi / dlti);
 #endif
 #ifdef SPHLATCH_TIMEDEP_ENERGY
      // old symmetric version
      //_i->dudt = (dudti / mi);

      // dervied from cont. eq. and 1st thermodynamic law
      _i->dudt = (dudti * pi / (dlti * dlti) + dudtavi) / mi;
 #endif
 #ifdef SPHLATCH_TRACK_UAV
      _i->dudtav = (dudtavi / mi);
 #endif
 #ifdef SPHLATCH_INTEGRATERHO
      _i->drhodt = -divvi * mi * dlti;
 #endif

#else
      _i->acc += acci;
 #ifdef SPHLATCH_TRACK_ACCP
      _i->accp = acci;
 #endif
 #ifdef SPHLATCH_VELDIV
      _i->divv = divvi;
 #endif
 #ifdef SPHLATCH_TIMEDEP_ENERGY
      _i->dudt = dudti;
 #endif
 #ifdef SPHLATCH_TRACK_UAV
      _i->dudtav = dudtavi;
 #endif
 #ifdef SPHLATCH_INTEGRATERHO
      _i->drhodt = -divvi * rhoi;
 #endif
#endif
#ifdef SPHLATCH_NONEIGH
      _i->noneigh = noneighi;
#endif
   }

   vect3dT vi, acci;
   fType   rhoi, pi, ci;

#ifdef SPHLATCH_MISCIBLE
   fType   piOdltidlti, dlti, mi;
#else
   fType   piOrhoirhoi;
#endif

#ifdef SPHLATCH_TIMEDEP_ENERGY
   fType   dudti;
#endif
#ifdef SPHLATCH_TRACK_UAV
   fType   dudtavi;
#endif
#ifdef SPHLATCH_VELDIV
   fType   divvi;
#endif
#ifdef SPHLATCH_INTEGRATERHO
   fType   drhodti;
#endif
#ifdef SPHLATCH_NONEIGH
   cType   noneighi;
#endif
};



template<typename _partT>
void setDhDt(_partT& _i)
{
   const fType noNeighMin = (2. / 3.) * _i.noneighOpt;
   const fType noNeighMax = (5. / 3.) * _i.noneighOpt;

   const fType k1 = 0.5 * (1 + tanh((_i.noneigh - noNeighMin) / -5.));
   const fType k3 = 0.5 * (1 + tanh((_i.noneigh - noNeighMax) / 5.));
   const fType k2 = 1. - k1 - k3;

   _i.dhdt =
      (k1 * _i.divvmax - k3 * _i.divvmax + k2 *
       static_cast<fType>(1. / 3.) * _i.divv) * _i.h;
}
};
#endif
