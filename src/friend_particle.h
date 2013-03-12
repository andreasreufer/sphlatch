#ifndef FRIEND_PARTICLE_H
#define FRIEND_PARTICLE_H

/*
 *  friend_particle.h
 *
 *
 *  Created by Andreas Reufer on 17.10.10
 *  Copyright 2010 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
enum friendType
{
   FRIENDONLIST = -2,
   FRIENDNOTSET = -1,
   FRIENDNONE   = 0,
};

class friendPart
{
public:
   iType friendid;
};
};
#endif
