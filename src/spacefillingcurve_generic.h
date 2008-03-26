#ifndef SPACEFILLINGCURVE_GENERIC_H
#define SPACEFILLINGCURVE_GENERIC_H

/*
 *  spacefillingcurve_generic.h
 *
 *
 *  Created by Andreas Reufer on 23.03.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

namespace sphlatch {
template<class T_leaftype>
class SpaceFillingCurve {

private:
T_leaftype& asLeaf()
{
  return static_cast<T_leaftype&>(*this);
}

public:
SpaceFillingCurve()
{
  depth = 0;
}

~SpaceFillingCurve()
{
}

void init(size_t _depth)
{
  depth = _depth;
  const size_t curveElems1D = (2 << (depth - 1)); // 2^(depth-1)
  const size_t curveElems2D = curveElems1D * curveElems1D;
  const size_t curveElems3D = curveElems2D * curveElems1D;

  curveCache.resize(curveElems3D);
  invCurveCache.resize(curveElems3D);

  size_t xIndex, yIndex, zIndex, cartIndex;
  for (size_t curveIndex = 0; curveIndex < curveElems3D; curveIndex++)
    {
      ///
      /// fill the caches with the space filling curve
      ///
      asLeaf().curveIndexToCartIndex(curveIndex, xIndex, yIndex, zIndex);
      cartIndex = xIndex + yIndex*curveElems1D + zIndex*curveElems2D;
      
      curveCache[curveIndex] = cartIndex;
      invCurveCache[cartIndex] = curveIndex;
    }
}

size_t getDepth()
{
  return depth;
}

size_t cartIndexToCurveIndex(const size_t &_cartIndex)
{
  return invCurveCache[_cartIndex];
}

size_t curveIndexToCartIndex(const size_t &_curveIndex)
{
  return curveCache[_curveIndex];
}

protected:
size_t depth;

private:
std::vector<size_t> curveCache, invCurveCache;

};
};

#endif

