#ifndef SPHLATCH_LOOKUP_TABLE_INTERPOL_LINEAR
#define SPHLATCH_LOOKUP_TABLE_INTERPOL_LINEAR

/*
 *  lookup_table_interpol_linear.h
 *
 *  this class provides linear interpolation
 *  for the look-up table
 *
 *  Created by Andreas Reufer on 14.08.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
class InterpolateLinear : public LookupTable<InterpolateLinear> {
public:
valueType interpolate(size_t iLo, size_t iHi, valueType _x)
{
  const valueType xLo = x(iLo);
  const valueType xHi = x(iHi);

  const valueType kHi = (_x - xLo) / (xHi - xLo);
  const valueType kLo = 1. - kHi;

  return(kLo * f(iLo) + kHi * f(iHi));
};
};
}
#endif

