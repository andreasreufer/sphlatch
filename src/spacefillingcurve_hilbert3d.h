#ifndef SPACEFILLINGCURVE_HILBERT3D_H
#define SPACEFILLINGCURVE_HILBERT3D_H

/*
 *  spacefillingcurve_hilbert3d.h
 *
 *
 *  Created by Andreas Reufer on 23.03.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

namespace sphlatch {
class Hilbert3D : public SpaceFillingCurve<Hilbert3D>
{
public:
void curveIndexToCartIndex(const size_t &_crvIndex, size_t &_xIndex,
                           size_t &_yIndex, size_t &_zIndex)
{
  bitsetType indexBS(3 * depth, _crvIndex);
  size_t xOld = 0;
  size_t yOld = 0;
  size_t zOld = 0;

  _xIndex = 0;
  _yIndex = 0;
  _zIndex = 0;

  for (size_t curDepth = 0; curDepth < depth; curDepth++)
    {
      const size_t curDepthMult = ((2 << curDepth) >> 1);   /// thats 2^(curDepth-1)
      const size_t octant = (_crvIndex >> (3 * curDepth)) % 8;

      xOld = _xIndex;
      yOld = _yIndex;
      zOld = _zIndex;

      switch (octant)
        {
        case 0:
          _xIndex = zOld;
          _yIndex = yOld;
          _zIndex = xOld;
          break;

        case 1:
          _xIndex = yOld;
          _yIndex = xOld;
          _zIndex = zOld + curDepthMult;
          break;

        case 2:
          _xIndex = yOld;
          _yIndex = xOld + curDepthMult;
          _zIndex = zOld + curDepthMult;
          break;

        case 3:
          _xIndex = xOld;
          _yIndex = 2 * curDepthMult - 1 - yOld;
          _zIndex = curDepthMult - 1 - zOld;
          break;

        case 4:
          _xIndex = curDepthMult + xOld;
          _yIndex = 2 * curDepthMult - 1 - yOld;
          _zIndex = curDepthMult - 1 - zOld;
          break;

        case 5:
          _xIndex = 2 * curDepthMult - 1 - yOld;
          _yIndex = 2 * curDepthMult - 1 - xOld;
          _zIndex = curDepthMult + zOld;
          break;

        case 6:
          _xIndex = 2 * curDepthMult - 1 - yOld;
          _yIndex = curDepthMult - 1 - xOld;
          _zIndex = curDepthMult + zOld;
          break;

        case 7:
          _xIndex = 2 * curDepthMult - 1 - zOld;
          _yIndex = yOld;
          _zIndex = curDepthMult - 1 - xOld;
          break;
        }
    }
}
};
};

#endif

