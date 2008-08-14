#ifndef SPHLATCH_LOOKUP_TABLE_INTERPOL_STEPWISE
#define SPHLATCH_LOOKUP_TABLE_INTERPOL_STEPWISE

/*
 *  lookup_table_interpol_stepwise.h
 *
 *  this class provides stepwise interpolation
 *  for the look-up table
 *
 *  Created by Andreas Reufer on 14.08.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
class InterpolateStepwise : public LookupTable<InterpolateStepwise> {
public:
valueType interpolate(size_t iLo, size_t iHi, valueType _x)
{
  const valueType xLo = x(iLo);
  const valueType xHi = x(iHi);

  const valueType dx = xHi - xLo;

  if ( _x < xLo + 0.5*dx )
    return f(iLo);
  else
    return f(iHi);
};
};
}
#endif

