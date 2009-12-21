#ifndef SPACEFILLINGCURVE_CARTXYZ_H
#define SPACEFILLINGCURVE_CARTXYZ_H

/*
 *  spacefillingcurve_cart.h
 *
 *
 *  Created by Andreas Reufer on 26.03.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

namespace sphlatch {
///
/// cartesian walk with X being the fastest and Z the slowest running index
///
class CartesianXYZ : public SpaceFillingCurve<CartesianXYZ>
{
public:
void curveIndexToCartIndex(const size_t &_crvIndex, size_t &_xIndex,
                           size_t &_yIndex, size_t &_zIndex)
{
  const size_t noElems1D = (1 << depth);
  const size_t noElems2D = noElems1D * noElems1D;

  _xIndex = _crvIndex % noElems1D;
  _yIndex = ((_crvIndex - _xIndex) % noElems2D) / noElems1D;
  _zIndex = (_crvIndex - _xIndex - _yIndex * noElems1D) / noElems2D;
}
};


///
/// cartesian walk with Y being the fastest and X the slowest running index
///
class CartesianYZX : public SpaceFillingCurve<CartesianYZX>
{
public:
void curveIndexToCartIndex(const size_t &_crvIndex, size_t &_xIndex,
                           size_t &_yIndex, size_t &_zIndex)
{
  const size_t noElems1D = (1 << depth);
  const size_t noElems2D = noElems1D * noElems1D;

  _yIndex = _crvIndex % noElems1D;
  _zIndex = ((_crvIndex - _yIndex) % noElems2D) / noElems1D;
  _xIndex = (_crvIndex - _yIndex - _zIndex * noElems1D) / noElems2D;
}
};

///
/// cartesian walk with Z being the fastest and Y the slowest running index
///
class CartesianZXY : public SpaceFillingCurve<CartesianZXY>
{
public:
void curveIndexToCartIndex(const size_t &_crvIndex, size_t &_xIndex,
                           size_t &_yIndex, size_t &_zIndex)
{
  const size_t noElems1D = (1 << depth);
  const size_t noElems2D = noElems1D * noElems1D;

  _zIndex = _crvIndex % noElems1D;
  _xIndex = ((_crvIndex - _zIndex) % noElems2D) / noElems1D;
  _yIndex = (_crvIndex - _zIndex - _xIndex * noElems1D) / noElems2D;
}
};
};

#endif

