#ifndef SPHLATCH_SPH_ALGORITHMS_CPP
#define SPHLATCH_SPH_ALGORITHMS_CPP

#include "typedefs.h"

namespace sphlatch {
template<typename _partT, typename _krnlT>
struct densSum
{
   _krnlT K;

   fType  rhoi;

   void   preSum(_partT* const _i)
   {
      rhoi = 0.;
   }

   void operator()(_partT* const _i,
                   const _partT* const _j,
                   const vect3dT& _rvec,
                   const fType _rr,
                   const fType _srad)
   {
      const fType r   = sqrt(_rr);
      const fType hij = 0.25*_srad + 0.5*(_j->h);
      const fType mj  = _j->m;

      rhoi += mj * K.value(r, hij);
   }

   void postSum(_partT* const _i)
   {
      _i->rho = rhoi;
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
      const fType hij = 0.25*_srad + 0.5*(_j->h);

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


      const fType accTerm     = piOrhoirhoi + (pj / (rhoj * rhoj)) + av;
      const fType mjvijdivWij = mj * dot(vij, K.deriv);

      K.derive(r, hij, _rvec);

      acci -= mj * accTerm * K.deriv;
      
      dudti += 0.5 * accTerm * mjvijdivWij;
   }

   void postSum(_partT* const _i)
   {
      _i->acc  += acci;
      _i->dudt = dudti;
   }

   vect3dT vi, acci;
   fType   rhoi, pi, ci, piOrhoirhoi, dudti;
};


///
/// 3D M4 cubic spline kernel
///
class CubicSpline3D
{
public:
   CubicSpline3D() { }
   ~CubicSpline3D() { }

   fType value(const fType& _r, const fType& _h)
   {
      const fType q = _r / _h;

      if (q > 2.)
         return(0.);
      else
      {
         const fType k1 = (1. / (_h * _h * _h * M_PI));

         if (q > 1.)
            return(0.25 * (2. - q) * (2. - q) * (2. - q) * k1);
         else
            return((1. - (3. / 2.) * q * q + 0.75 * q * q * q) * k1);
      }
   }

   void derive(const fType& _r, const fType& _h, const vect3dT& _rvec)
   {
      const fType q = _r / _h;

      // when q = 0, the kernel also has to be (0,0,0)
      if ((q > 2.) || (q == 0.))
      {
         deriv = 0., 0., 0.;
         return;
      }
      else
      {
         const fType k2 = 1. / (_h * _h * _h * _h * M_PI * _r);
         deriv = _rvec * k2;
         if (q > 1.)
         {
            const fType k3 = -0.75 * (2. - q) * (2. - q);
            deriv *= k3;
            return;
         }
         else
         {
            const fType k3 = -3. * q + 2.25 * q * q;
            deriv *= k3;
            return;
         }
      }
   }

   vect3dT deriv;
};
};

#endif

///
/// 2D M4 cubic spline kernel
///

/* class CubicSpline2D
   {

   public:
   CubicSpline2D() {};
   ~CubicSpline2D() {};

   fType value(const fType& _r, const fType& _h)
   {
   const fType q = _r / _h;

   if ( q > 2. )
   {
    return 0.;
   }
   else
   {
    const fType k1 = ( 10. / ( 7.*_h*_h*M_PI) );

    if ( q > 1. )
    {
      return 0.25*( 2. - q )*( 2. - q )*( 2. - q )*k1;
    }
    else
    {
      return ( 1. - (3./2.)*q*q + 0.75*q*q*q )*k1;
    }
   }
   }

   void derive( const fType& _r, const fType& _h,
             const fType& _rx, const fType& _ry)
   {
   const fType q = _r / _h;

   ///
   /// when q = 0, the kernel also has to be (0,0,0)
   ///
   if ( q > 2. || q == 0.)
   {
    derivX = 0.;
    derivY = 0.;
    return;
   }
   else
   {
    const fType k2 = 10. / ( 7.*_h*_h*_h*M_PI*_r );

    derivX = _rx * k2;
    derivY = _ry * k2;

    if ( q > 1. )
    {
      const fType k3 = -0.75*( 2. - q )*( 2. - q );
      derivX *= k3;
      derivY *= k3;
      return;
    }
    else
    {
      const fType k3 = - 3.*q + 2.25*q*q;
      derivX *= k3;
      derivY *= k3;
      return;
    }
   }
   }

   fType derivX, derivY;

   };*/
