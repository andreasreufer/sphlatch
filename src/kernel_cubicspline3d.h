#ifndef KERNEL_CUBICSPLINE3D_H
#define KERNEL_CUBICSPLINE3D_H

/*
 *  kernel_cubicspline3d.h
 *
 *
 *  Created by Andreas Reufer on 17.06.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
///
/// 3D M4 cubic spline kernel
///
class CubicSpline3D
{

public:
CubicSpline3D() {};
~CubicSpline3D() {};

valueType value(const valueType& _r, const valueType& _h)
{
  const valueType q = _r / _h;

  if ( q > 2. )
  {
    return 0.;
  }
  else
  {
    const valueType k1 = ( 1. / ( _h*_h*_h*M_PI) );

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

void derive( const valueType& _r, const valueType& _h,
             const valueType& _rx, const valueType& _ry, const valueType& _rz)
{
  const valueType q = _r / _h;

  ///
  /// when q = 0, the kernel also has to be (0,0,0)
  /// 
  if ( q > 2. || q == 0.)
  {
    derivX = 0.;
    derivY = 0.;
    derivZ = 0.;
    return;
  }
  else
  {
    const valueType k2 = 1. / ( _h*_h*_h*_h*M_PI*_r );

    derivX = _rx * k2;
    derivY = _ry * k2;
    derivZ = _rz * k2;

    if ( q > 1. )
    {
      const valueType k3 = -0.75*( 2. - q )*( 2. - q );
      derivX *= k3;
      derivY *= k3;
      derivZ *= k3;
      return;
    }
    else
    {
      const valueType k3 = - 3.*q + 2.25*q*q;
      derivX *= k3;
      derivY *= k3;
      derivZ *= k3;
      return;
    }
  } 
}

valueType derivX, derivY, derivZ;

};
};
#endif

