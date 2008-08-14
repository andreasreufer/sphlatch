#ifndef SPHLATCH_LOOKUP_TABLE
#define SPHLATCH_LOOKUP_TABLE

/*
 *  lookup_table.h
 *
 *  this class generates a look-up table usable like a mathe-
 *  matical function from two tables f(x) and x
 *
 *  x needs to be monotonous ascending
 *
 *  Created by Andreas Reufer on 13.08.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
namespace sphlatch {
template<class T_leaftype>
class LookupTable {
public:
LookupTable(valvectType _x, valvectType _f)
{
  assert(_x.size() == _f.size());
  noPoints = _x.size();

  f.resize(noPoints);
  f = _f;

  x.resize(noPoints);
  x = _x;

  xMin = x(0);
  xMax = x(noPoints - 1);

  yMin = f(0);
  yMax = f(noPoints - 1);
};

~LookupTable()
{
};

public:
valueType operator()(valueType _x);

private:
T_leaftype& asLeaf()
{
  return static_cast<T_leaftype&>(*this);
};

protected:
valueType xMin, xMax;
valueType yMin, yMax;
valvectType f, x;
size_t noPoints;
};

template<class T_leaftype>
valueType LookupTable<T_leaftype>::operator()(valueType _x)
{
  ///
  /// return a constant extrapolation
  /// outside the tables range
  ///
  if (_x < xMin)
    return yMin;

  if (_x > xMax)
    return yMax;

  ///
  /// bracket the _x value
  ///
  size_t iLo = 0, iHi = noPoints - 1;
  while ((iHi - iLo) > 1)
    {
      const size_t iMd = (iHi + iLo) / 2;

      if (x(iMd) < _x)
        iLo = iMd;
      else
        iHi = iMd;
    }

  ///
  /// linearly interpolate
  ///
  return ( asLeaf().interpolate( iLo, iHi, _x) );
};
}

#include "lookup_table_interpol_linear.h"
#include "lookup_table_interpol_stepwise.h"

#endif

