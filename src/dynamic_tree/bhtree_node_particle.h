#ifndef BHTREE_NODE_PARTICLE_H
#define BHTREE_NODE_PARTICLE_H

//#define SPHLATCH_AMD64PADDING

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

namespace sphlatch {
///
/// particle node
///
class particleNode : public genericNode {
public:
   fType xPos, yPos, zPos;
   fType mass;

   particleNode() { clear(); }
   ~particleNode() { }

   void clear();

#ifdef SPHLATCH_AMD64PADDING
private:
   char pad[0];
#endif
};

void particleNode::clear()
{
   next = NULL;
   skip = NULL;

   depth = -1;
   ident = 0;

   isParticle = true;
   isRemote   = false;
   isCZ       = false;

   xPos = 0.;
   yPos = 0.;
   zPos = 0.;
   mass = 0.;
}
};
#endif
