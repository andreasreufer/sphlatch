#ifndef BHTREE_MONOPOLES_H
#define BHTREE_MONOPOLES_H

/*
 *  bhtree_monopoles.h
 *
 *
 *  Created by Andreas Reufer on 07.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "bhtree_monopole_node.h"

namespace sphlatch {
class Monopoles : public BHtree<Monopoles>
{

typedef genericNode* nodePtrT;
typedef particleNode* partPtrT;
typedef genericCellNode* cellPtrT;
typedef particleProxy* partProxyPtrT;

public:

///
/// allocates the root monopole cell node
///
void allocRootNode(void)
{
  rootPtr = new monopoleCellNode;
}

///
/// allocates a new monopole cell node and connects it as child _n
/// no check is performed, whether curNodePtr points to a cell node!
///
void allocNewCellChild(const size_t _n)
{
  // allocate new cell node
  monopoleCellNode* newNodePtr =
    new monopoleCellNode;
  
  // connect the new cell node to curNodePtr
  newNodePtr->parent = curNodePtr;
  static_cast<cellPtrT>(curNodePtr)->child[_n] = newNodePtr;
}

void calcGravCell()
{
}

private:
protected:

};
};

#endif
