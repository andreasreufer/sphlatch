#ifndef SPHLATCH_SPH_ALGORITHMS_CPP
#define SPHLATCH_SPH_ALGORITHMS_CPP

#include "typedefs.h"

namespace sphlatch {
template<typename _partT, typename _krnlT>
struct densSum
{
   _krnlT K;

   void   preSum(_partT* const _i)
   {
      _i->rho = 0.;
   }

   void postSum(_partT* const _i)
   { }

   void operator()(_partT* const _i,
                   const _partT* const _j,
                   const vect3dT& _rvec,
                   const fType _rr,
                   const fType _hi)
   {
      const fType r   = sqrt(_rr);
      const fType hij = 0.5 * (_hi + _j->h);

      _i->rho += (_j->m) * K.value(r, hij);
   }
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
