#ifndef PARTICLE_H
#define PARTICLE_H

//#define SPHLATCH_PADTO64BYTES

/*
 *  particle.h
 *
 *
 *  Created by Andreas Reufer on 23.02.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */


#include "typedefs.h"
#include "bhtree_node_particle.h"

namespace sphlatch {
class particleNode;

class treeParticle {
public:
   typedef particleNode*   partNodePtrT;
   fType        pos[3], m;
   partNodePtrT treeNode;

   fType eps;
};

};
#endif
