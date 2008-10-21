#ifndef SPHLATCH_LOOKUP_TABLE1D
#define SPHLATCH_LOOKUP_TABLE1D

/*
 *  lookup_table1D.h
 *
 *  this class generates a look-up table usable like a mathe-
 *  matical function from two tables f(x) and x
 *
 *  x needs to be monotonous ascending and better
 *  be regularly spaced
 *
 *  Created by Andreas Reufer on 13.08.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
namespace sphlatch {
template<class T_leaftype>
class LookupTable1D {
public:
LookupTable1D(valvectType _x, valvectType _f)
{
  assert(_x.size() == _f.size());
  nx = _x.size();

  f.resize(nx);
  f = _f;

  x.resize(nx);
  x = _x;

  xMin = x(0);
  xMax = x(nx - 1);

  fMin = f(0);
  fMax = f(nx - 1);
};

~LookupTable1D()
{
};

public:
valueType operator()(valueType _x);

private:
T_leaftype& asLeaf()
{
  return static_cast<T_leaftype&>(*this);
};

public:
size_t ixl, ixh;

protected:
valueType xMin, xMax;
valueType fMin, fMax;
valvectType f, x;
size_t nx;
};

template<class T_leaftype>
valueType LookupTable1D<T_leaftype>::operator()(const valueType _x)
{
  ///
  /// return a constant extrapolation
  /// outside the tables range
  ///
  if (_x < xMin)
    return fMin;

  if (_x > xMax)
    return fMax;

  ///
  /// bracket the _x value
  ///
  ixl = 0;
  ixh = nx - 1;
  while ((ixh - ixl) > 1)
    {
      const size_t ixm = (ixh + ixl) / 2;

      if (x(ixm) < _x)
        ixl = ixm;
      else
        ixh = ixm;
    }

  ///
  /// interpolate
  ///
  return(asLeaf().interpolate(_x, ixl));
};

///
/// this class provides linear interpolation
/// for the look-up table
///

class InterpolateLinear : public LookupTable1D<InterpolateLinear> {
public:
valueType interpolate(const valueType _x, const size_t ixl)
{
  const valueType xl = x(ixl);
  const valueType xh = x(ixl + 1);

  const valueType t = (_x - xl) / (xh - xl);

  return((1. - t) * f(ixl) + t * f(ixl + 1));
};
};

///
/// this class provides stepwise interpolation
/// for the look-up table
///
class InterpolateStepwise : public LookupTable1D<InterpolateStepwise> {
public:
valueType interpolate(const valueType _x, const size_t ixl)
{
  const valueType xl = x(ixl);
  const valueType xh = x(ixl + 1);

  const valueType dx = xh - xl;

  if (_x < xl + 0.5 * dx)
    return f(ixl);
  else
    return f(ixl + 1);
};
};
}
#endif

