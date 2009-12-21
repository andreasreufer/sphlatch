#ifndef SPHLATCH_LATTICE_HCP
#define SPHLATCH_LATTICE_HCP

/*
 *  lattice_hcp.h
 *
 *
 *  Created by Andreas Reufer on 26.07.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
class LatticeHCP {
public:
LatticeHCP(fType _spacing, fType _xMax, fType _yMax, fType _zMax)
{
  fillingFactor = sqrt(2)*M_PI / 6.;
  d = _spacing;
  d05 = 0.5 * d;
  dsq32 = d * sqrt(3.) / 2.;
  dsq34 = d * sqrt(3.) / 4.;

  xMin = -_xMax;
  yMin = -_yMax;
  zMin = -_zMax;

  xMax = _xMax;
  yMax = _yMax;
  zMax = _zMax;

  first();
};

~LatticeHCP()
{
};

public:
void first();
void next();
bool isLast;

fType xCur, yCur, zCur, rCur;
fType fillingFactor;

private:
void calcRad();
fType xOld, yOld;

fType xMin, yMin, zMin;
fType xMax, yMax, zMax;

fType d, d05, dsq32, dsq34;
};

void LatticeHCP::first()
{
  xCur = 0.; yCur = 0.; zCur = 0.;
  int ix = 0, iy = 0, iz = 0;

  while (zCur > zMin)
    {
      zCur -= dsq34;
      iz -= 1;
    }
  if (iz % 2 == 1)
    {
      xCur -= d05;
      yCur -= dsq34;
    }

  while (yCur > yMin)
    {
      yCur -= dsq32;
      iy -= 1;
    }
  if (iy % 2 == 1)
    xCur -= d05;

  while (xCur > xMin)
    {
      xCur -= d;
      ix -= 1;
    }

  xOld = xCur;
  yOld = yCur;

  isLast = false;

  calcRad();
};

void LatticeHCP::next()
{
  ///
  /// last lattice point?
  ///
  if (isLast)
    return;

  ///
  /// advance xCur
  ///
  xCur += d;

  if (xCur > xMax)
    {
      ///
      /// reset xCur and shift it
      ///
      xCur = xOld;
      if (xCur > xMin)
        xCur -= d05;
      else
        xCur += d05;
      xOld = xCur;

      ///
      /// advance yCur
      ///
      yCur += dsq32;
    }

  if (yCur > yMax)
    {
      ///
      /// reset yCur and shift it
      ///
      yCur = yOld;
      if (yCur > yMin)
        yCur -= dsq34;
      else
        yCur += dsq34;
      yOld = yCur;

      ///
      /// reset xCur and shift it
      ///
      if (xCur > xMin)
        xCur -= d05;
      else
        xCur += d05;
      xOld = xCur;

      ///
      /// advance zCur
      ///
      if (zCur > zMax)
        {
          isLast = true;
          return;
        }
      else
        zCur += dsq34;
    }
  calcRad();
};

void LatticeHCP::calcRad()
{
  ///
  /// calculate current radius
  ///
  rCur = sqrt(xCur * xCur +
              yCur * yCur +
              zCur * zCur);
};

}
#endif

