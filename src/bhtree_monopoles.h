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
typedef monopoleCellNode* monoPtrT;
typedef particleProxy* partProxyPtrT;

public:
///
/// allocates the root monopole cell node
///
void allocRootNode()
{
  rootPtr = new monopoleCellNode;
}

//deprecated
///
/// resize buffers for communication
///
void prepareBuffers()
{
  localCells.resize(noToptreeCells, MSIZE);
  localIsFilled.resize(noToptreeCells);

  remoteCells.resize(noToptreeCells, MSIZE);
  remoteIsFilled.resize(noToptreeCells);
}

///
/// report number of multipole moments
/// (includes center of mass and mass)
///
size_t noMultipoleMoments()
{
  return 4; // 3 center of mass + 1 monopole = 4 multipoles
}

///
/// allocates a new monopole cell node and connects it as child _n
/// no check is performed, whether curNodePtr points to a cell node!
///
void allocNewCellChild(const size_t _n)
{
// allocate new cell node
  monoPtrT newNodePtr =
    new monopoleCellNode;

// connect the new cell node to curNodePtr
  newNodePtr->parent = curNodePtr;
  static_cast<cellPtrT>(curNodePtr)->child[_n] = newNodePtr;

// set cell vars to zero
  static_cast<monoPtrT>(newNodePtr)->mass = 0.;
  static_cast<monoPtrT>(newNodePtr)->xCom = 0.;
  static_cast<monoPtrT>(newNodePtr)->yCom = 0.;
  static_cast<monoPtrT>(newNodePtr)->zCom = 0.;
}

///
/// calculate multipole from children
/// this function does not check, whether the current node is actually a cell!
///
void calcMultipole()
{
//
// add up the contributions from the children with
// the same locality as the current node.
// all this locality business guarantees, that every particle
// contributes only ONCE to the multipole moments of the global tree
// but "ghost" cells still contain the multipoles contributed by
// their ghost children.
//
// a special case are the deepest toptree cells which are per
// definition local. if all children are ghosts, none of if contri-
// butes anything so cm gets 0 and fucks up the center of mass
// of the cell. so if nobody contributes anything, omit the addition
// to the cell.
//
  static valvectType curCell(MSIZE), childCell(MSIZE);

  for (size_t i = 0; i < MSIZE; i++)
    {
      curCell[i] = 0.;
    }

// first calculate center of mass
  for (size_t i = 0; i < 8; i++)
    {
      if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i]->isLocal
              == curNodePtr->isLocal)
            {
              goChild(i);
              if (curNodePtr->isParticle == true)
                {
                  partToVect(childCell);
                }
              else
                {
                  cellToVect(childCell);
                }
              goUp();

              addCOM(curCell, childCell);
            }
        }
    }

  /*if (childCell[MASS] > 0.)
     {
      for (size_t i = 0; i < 8; i++)
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
            {
              if (static_cast<cellPtrT>(curNodePtr)->child[i]->isLocal
                  == curNodePtr->isLocal)
                {
                  goChild(i);
                  if (curNodePtr->isParticle == true)
                    {
                      partToVect(childCell);
                    }
                  else
                    {
                      cellToVect(childCell);
                    }
                  goUp();

                  addMP(curCell, childCell);
                }
            }
        }
     }*/

// copy data to node itself ...
  if (curCell[MASS] > 0.)
    {
      vectToCell(curCell);
    }
}

///
/// add center of mass
///
void addCOM(valvectRefType _target, const valvectRefType _source)
{
  //const static does not seem to work here
  const valueType oldMass = _target[MASS];

  _target[MASS] += _source[MASS];

  if (_target[MASS] > 0.)
    {
      // const static does not seem to work here
      const valueType newMassInv = (1.0 / _target[MASS]);
      _target[CX] = newMassInv * (oldMass * _target[CX]
                                  + _source[MASS] * _source[CX]);
      _target[CY] = newMassInv * (oldMass * _target[CY]
                                  + _source[MASS] * _source[CY]);
      _target[CZ] = newMassInv * (oldMass * _target[CZ]
                                  + _source[MASS] * _source[CZ]);
    }
}

///
/// add multipole moments
/// either for a cell by its children multipole moments or to add up
/// multipoles globally in the top-tree
///
void addMP(valvectRefType _target, const valvectRefType _source)
{
  ///
  /// everything is already done in addCOM()
  ///
}

//deprecated
///
/// merge multipole moments in remote and local buffers
///
void mergeRemoteCells()
{
  static valvectType localCell, remoteCell;

  valueType oldMass, newMass;

  for (size_t i = 0; i < noToptreeCells; i++)
    {
      if (remoteIsFilled[i])
        {
          if (localIsFilled[i])
            {
              oldMass = localCells(i, MASS);
              newMass = oldMass + remoteCells(i, MASS);

              //
              // newMass may be zero for non-empty cells. this
              // happens when all children are ghosts and do not
              // contribute to the local toptree
              //
              if (newMass > 0.)
                {
                  localCells(i, MASS) = newMass;
                  localCells(i, CX) = ((oldMass / newMass) *
                                       localCells(i, CX))
                                      + ((remoteCells(i, MASS) / newMass) *
                                         remoteCells(i, CX));
                  localCells(i, CY) = ((oldMass / newMass) *
                                       localCells(i, CY))
                                      + ((remoteCells(i, MASS) / newMass) *
                                         remoteCells(i, CY));
                  localCells(i, CZ) = ((oldMass / newMass) *
                                       localCells(i, CZ))
                                      + ((remoteCells(i, MASS) / newMass) *
                                         remoteCells(i, CZ));
                }
            }
          else
            {
              localCells(i, MASS) = remoteCells(i, MASS);
              localCells(i, CX) = remoteCells(i, CX);
              localCells(i, CY) = remoteCells(i, CY);
              localCells(i, CZ) = remoteCells(i, CZ);
            }
        }
    }

  /*for (size_t i = 0; i < noToptreeCells; i++)
     {
      if (remoteIsFilled[i])
        {
          localCell = particleRowType(localCells, i);
          remoteCell = particleRowType(remoteCells, i);

          if (localIsFilled[i])
            {
              ///
              /// newMass may be zero for non-empty cells. this
              /// happens when all children are ghosts and do not
              /// contribute to the local toptree. this is handled
              /// by addCOM()
              ///
              addCOM(localCell, remoteCell);

              /// now add up multipoles relative to new COM
              addMP(localCell, remoteCell);
            }
          else
            {
              localCell = remoteCell;
            }
        }
     }*/
}

//deprecated
///
/// copy the current cell node to buffer
///
void cellToBuffer()
{
  localCells(toptreeCounter, CX) = static_cast<monoPtrT>(curNodePtr)->xCom;
  localCells(toptreeCounter, CY) = static_cast<monoPtrT>(curNodePtr)->yCom;
  localCells(toptreeCounter, CZ) = static_cast<monoPtrT>(curNodePtr)->zCom;
  localCells(toptreeCounter, MASS) = static_cast<monoPtrT>(curNodePtr)->mass;
}

//deprecated
///
/// copy buffer to current cell node
///
void bufferToCell()
{
  static_cast<monoPtrT>(curNodePtr)->xCom = localCells(toptreeCounter, CX);
  static_cast<monoPtrT>(curNodePtr)->yCom = localCells(toptreeCounter, CY);
  static_cast<monoPtrT>(curNodePtr)->zCom = localCells(toptreeCounter, CZ);
  static_cast<monoPtrT>(curNodePtr)->mass = localCells(toptreeCounter, MASS);
}

///
/// copy current cell node to a vector
/// no checks are performed, whether current node is a
/// cell nor if vector has the right size
///
void cellToVect(valvectRefType _vect)
{
  _vect[CX] = static_cast<monoPtrT>(curNodePtr)->xCom;
  _vect[CY] = static_cast<monoPtrT>(curNodePtr)->yCom;
  _vect[CZ] = static_cast<monoPtrT>(curNodePtr)->zCom;
  _vect[MASS] = static_cast<monoPtrT>(curNodePtr)->mass;
}

///
/// copy current particle node to a vector
/// no checks are performed, whether current node is a
/// particle nor if vector has the right size
///
void partToVect(valvectRefType _vect)
{
  _vect[CX] = static_cast<partPtrT>(curNodePtr)->xPos;
  _vect[CY] = static_cast<partPtrT>(curNodePtr)->yPos;
  _vect[CZ] = static_cast<partPtrT>(curNodePtr)->zPos;
  _vect[MASS] = static_cast<partPtrT>(curNodePtr)->mass;
}

///
/// copy vector to current cell node
/// no checks are performed, whether current node is a
/// cell nor if vector has the right size
///
void vectToCell(valvectRefType _vect)
{
  static_cast<monoPtrT>(curNodePtr)->xCom = _vect[CX];
  static_cast<monoPtrT>(curNodePtr)->yCom = _vect[CY];
  static_cast<monoPtrT>(curNodePtr)->zCom = _vect[CZ];
  static_cast<monoPtrT>(curNodePtr)->mass = _vect[MASS];
}


///
/// report current cell node to dumpFile stream
///
void reportMultipoles()
{
  dumpFile << static_cast<monoPtrT>(curNodePtr)->xCom << "   ";
  dumpFile << static_cast<monoPtrT>(curNodePtr)->yCom << "   ";
  dumpFile << static_cast<monoPtrT>(curNodePtr)->zCom << "   ";
  dumpFile << static_cast<monoPtrT>(curNodePtr)->mass << "   ";
}

// why is this here?
///
/// stop recursion if:
/// - current node is empty << ??
/// - MAC is fulfilled
///
bool calcGravMAC(void)
{
  // any way to speed this up?
  cellPartDist = sqrt(
    (static_cast<monoPtrT>(curNodePtr)->xCom - curGravParticleX) *
    (static_cast<monoPtrT>(curNodePtr)->xCom - curGravParticleX) +
    (static_cast<monoPtrT>(curNodePtr)->yCom - curGravParticleY) *
    (static_cast<monoPtrT>(curNodePtr)->yCom - curGravParticleY) +
    (static_cast<monoPtrT>(curNodePtr)->zCom - curGravParticleZ) *
    (static_cast<monoPtrT>(curNodePtr)->zCom - curGravParticleZ)
    );
  return(((static_cast<monoPtrT>(curNodePtr)->cellSize)
          / cellPartDist) < thetaMAC);
}


///
/// calculate acceleration due to a cell
/// no check whether current node is actually a cell!
///

void calcGravCell()
{
  // cellPartDist is already set by the MAC function

  // intermediate results for monopole term
  const valueType rInvPow3 = 1. / (cellPartDist * cellPartDist * cellPartDist);

  const valueType rx = curGravParticleX -
                       static_cast<monoPtrT>(curNodePtr)->xCom;
  const valueType ry = curGravParticleY -
                       static_cast<monoPtrT>(curNodePtr)->yCom;
  const valueType rz = curGravParticleZ -
                       static_cast<monoPtrT>(curNodePtr)->zCom;
  const valueType mass = static_cast<monoPtrT>(curNodePtr)->mass;

  // gravity due to monopole term
  curGravParticleAX -= (rInvPow3) * mass * rx;
  curGravParticleAY -= (rInvPow3) * mass * ry;
  curGravParticleAZ -= (rInvPow3) * mass * rz;
}

private:
enum monopoleIndex { CX, CY, CZ, MASS, MSIZE };

private:
};
};

#endif
