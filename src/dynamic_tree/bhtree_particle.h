#ifndef BHTREE_PARTICLE_H
#define BHTREE_PARTICLE_H

//#define SPHLATCH_PADTO64BYTES

/*
 *  bhtree_particle.h
 *
 *
 *  Created by Andreas Reufer on 23.02.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {

class particleNode;
///
/// basic tree particle class
///
class treeGhost {
public:
   typedef particleNode*   partNodePtrT;
   vect3dT      pos;
   fType        m, eps;
   partNodePtrT treeNode;

   idType id;
   fType  cost;

#ifdef SPHLATCH_PADTO64BYTES
private:
   char pad[8];
#endif
};

///
/// dummy resident class
/// 
class treePart : public treeGhost {
};

};
#endif
