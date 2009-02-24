#ifndef BHTREE_NODE_PARTICLE_H
#define BHTREE_NODE_PARTICLE_H

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
#include "bhtree_node_generic.h"
#include "particle.h"

namespace sphlatch {
///
/// particle node
///

class treeParticle;

class particleNode : public genericNode {
public:
   typedef treeParticle*   partPtrT;
   partPtrT partPtr;
   vect3dT  pos;
   fType    m;

   particleNode() { }
   ~particleNode() { }

   void clear();
   void update();

#ifdef SPHLATCH_PADTO64BYTES
private:
   char pad[0];
#endif
};

void particleNode::clear()
{
   genericNode::clear();

   isParticle = true;

   partPtr = NULL;
   pos     = 0., 0., 0.;
   m       = 0.;
}

void particleNode::update()
{
   assert(partPtr != NULL);
   pos = partPtr->pos;
   m   = partPtr->m;
}

};
#endif
