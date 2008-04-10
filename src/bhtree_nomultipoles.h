#ifndef BHTREE_NOMULTIPOLES_H
#define BHTREE_NOMULTIPOLES_H

/*
 *  bhtree_nomultipoles.h
 *
 *
 *  Created by Andreas Reufer on 07.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "bhtree_generic_node.h"

namespace sphlatch {
class NoMultipoles : public BHtree<NoMultipoles>
{
typedef genericNode* nodePtrT;
typedef particleNode* partPtrT;
typedef genericCellNode* cellPtrT;
typedef particleProxy* partProxyPtrT;

public:
///
/// allocates the root monopole cell node
///
void allocRootNode()
{
  rootPtr = new genericCellNode;
}

///
/// report number of multipole moments
/// (includes center of mass and mass)
///
size_t noMultipoleMoments()
{
  return 0; ///  0 multipoles
}

///
/// allocates a new monopole cell node and connects it as child _n
/// no check is performed, whether curNodePtr points to a cell node!
///
void allocNewCellChild(const size_t _n)
{
// allocate new cell node
  cellPtrT newNodePtr =
    new genericCellNode;

// connect the new cell node to curNodePtr
  newNodePtr->parent = curNodePtr;
  static_cast<cellPtrT>(curNodePtr)->child[_n] = newNodePtr;
}

///
/// dummy function
///
void calcMultipole()
{
}

///
/// dummy function
///
void addCOM(valvectRefType _target, const valvectRefType _source)
{
}

///
/// dummy function
///
void addMP(valvectRefType _target, const valvectRefType _source)
{
}

///
/// dummy function
///
void cellToVect(valvectRefType _vect)
{
}

///
/// copy current particle node to a vector
/// no checks are performed, whether current node is a
/// particle nor if vector has the right size
///
void partToVect(valvectRefType _vect)
{
  _vect[X] = static_cast<partPtrT>(curNodePtr)->xPos;
  _vect[Y] = static_cast<partPtrT>(curNodePtr)->yPos;
  _vect[Z] = static_cast<partPtrT>(curNodePtr)->zPos;
  _vect[M] = static_cast<partPtrT>(curNodePtr)->mass;
}

///
/// dummy function
///
void vectToCell(valvectRefType _vect)
{
}

///
/// dummy function, abort the tree walk immediately
///
bool calcGravMAC(void)
{
  return true;
}

///
/// dummy function
///
void calcGravCell()
{
}

};
};

#endif
