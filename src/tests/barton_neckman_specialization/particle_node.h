#ifndef PARTICLE_NODE_H
#define PARTICLE_NODE_H

/*
 *  particle_node.h
 *
 *
 *  Created by Andreas Reufer on 06.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
/** \brief tree node struct with an
 *           arbitrary payload
 */

struct particleNode : public genericNode {
  typedef particleProxy*      particleProxyPtrType;

  /**
   * pointers to children
   */
  particleProxyPtrType partProxy;

  /**
   * center and size of the cell
   */
  valueType xPos, yPos, zPos;
  valueType mass;
};
};

#endif

