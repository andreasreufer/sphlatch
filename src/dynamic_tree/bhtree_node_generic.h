#ifndef BHTREE_NODE_GENERIC_H
#define BHTREE_NODE_GENERIC_H

//#define SPHLATCH_AMD64PADDING

/*
 *  bhtree_node_generic.h
 *
 *
 *  Created by Andreas Reufer on 20.01.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
///
/// base class for all nodes
///
class genericNode {
public:
   typedef genericNode*   genericNodePtrT;

   genericNodePtrT parent, next, skip;

   size_t    depth;
   identType ident;

   bool isParticle : 1;
   bool isRemote   : 1;
   bool isCZ       : 1;
   bool neighSet    : 1;

   genericNode() { clear(); }
   ~genericNode() { }

   void clear();
   genericNodePtrT operator*();
};

void genericNode::clear()
{
   next = NULL;
   skip = NULL;

   depth = -1;
   ident = 0;

   isParticle = false;
   isRemote   = false;
   isCZ       = false;
   neighSet   = false;
}

genericNode::genericNodePtrT genericNode::operator*()
{
   return(this);
}
};

#endif
