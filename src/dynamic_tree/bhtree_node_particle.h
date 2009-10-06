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
#include "bhtree_particle.h"

namespace sphlatch {
///
/// particle node
///

class treeGhost;

class particleNode : public genericNode {
public:
   typedef treeGhost*   partPtrT;
   
   partPtrT partPtr;
   vect3dT  pos;
   fType    m;

   particleNode() { }
   ~particleNode() { }

   void clear();
   void update();

#ifdef SPHLATCH_PADD64
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
   pos     = partPtr->pos;
   m       = partPtr->m;
}


};
#endif
