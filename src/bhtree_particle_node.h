#ifndef BHTREE_PARTICLE_NODE_H
#define BHTREE_PARTICLE_NODE_H

/*
 *  particle_node.h
 *
 *
 *  Created by Andreas Reufer on 06.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "bhtree_generic_node.h"

namespace sphlatch {

struct particleNode : public genericNode {
  // particle position and mass
  // are cached locally
  valueType xPos, yPos, zPos;
  valueType mass;
};
};

#endif

