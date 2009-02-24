#ifndef BHTREE_NODE_GENERIC_H
#define BHTREE_NODE_GENERIC_H

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

   genericNodePtrT parent, next;

   idType             ident;
   short unsigned int depth;

   ///
   /// is a node a particle / a CZ cell, is it remote
   ///
   bool isParticle  : 1;
   bool isCZ        : 1;
   bool isRemote    : 1;

   ///
   /// indicate whether a CZ cell is at the bottom
   /// of the CZ tree and whether the neighbours are set
   ///
   bool atBottom    : 1;
   bool neighSet    : 1;

   ///
   /// indicates whether a particle node is settled
   /// at the correct position in the tree
   ///
   bool isSettled   : 1;

   genericNode() { }
   ~genericNode() { }

   void clear();
   genericNodePtrT operator*();

private:
};

void genericNode::clear()
{
   next = NULL;

   depth = -1;
   ident = 0;

   isParticle = false;
   isCZ       = false;
   isRemote   = false;

   atBottom = false;
   neighSet = false;

   isSettled = false;
}

genericNode::genericNodePtrT genericNode::operator*()
{
   return(this);
}
};

#endif
