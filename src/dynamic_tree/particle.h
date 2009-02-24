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

namespace sphlatch {
class particleNode;

class treeParticle {
public:
   typedef particleNode*   partNodePtrT;
   vect3dT      pos;
   fType        m, cost;
   partNodePtrT treeNode;

   fType eps;
};
};
#endif
