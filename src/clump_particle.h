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
   CLUMPNONE = -1, CLUMPNOTSET = 0
};

class clumpPart
{
public:
   iType clumpid;
};
};
#endif
