#ifndef PARTICLE_H
#define PARTICLE_H

//#define SPHLATCH_PADTO64BYTES

/*
 *  bhtree_node_particle.h
 *
 *
 *  Created by Andreas Reufer on 20.01.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */


#include "typedefs.h"
#include "bhtree_node_particle.h"

namespace sphlatch {
class particleNode;

class particleGeneric {
public:
   typedef particleNode*   partNodePtrT;
   fType        x, y, z, m;
   partNodePtrT treeNode;
};

};
#endif
