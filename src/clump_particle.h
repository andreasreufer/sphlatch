#ifndef CLUMP_PARTICLE_H
#define CLUMP_PARTICLE_H

/*
 *  clump_particle.h
 *
 *
 *  Created by Andreas Reufer on 15.01.10
 *  Copyright 2010 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
enum clumpType
{
   CLUMPONLIST = -2,
   CLUMPNOTSET = -1,
   CLUMPNONE   = 0,
};

enum orbitType
{
  ORBITNOTSET   = 0,
  ORBITCLUMP    = 1,
  ORBITDISK     = 2,
  ORBITREIMPACT = 3,
  ORBITESCAPE   = 4
};

class clumpPart
{
public:
   iType clumpid, orbit;
   fType ecc, a;
};
};
#endif
