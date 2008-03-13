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
  static valueType cm, cxm, cym, czm;

  cm = 0.;
  cxm = 0.;
  cym = 0.;
  czm = 0.;

  static valueType rx, ry, rz, rr;
  rx = 0.;
  ry = 0.;
  rz = 0.;
  rr = 0.;

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
                  cm += static_cast<partPtrT>(curNodePtr)->mass;
                  cxm += (static_cast<partPtrT>(curNodePtr)->xPos) *
                         (static_cast<partPtrT>(curNodePtr)->mass);
                  cym += (static_cast<partPtrT>(curNodePtr)->yPos) *
                         (static_cast<partPtrT>(curNodePtr)->mass);
                  czm += (static_cast<partPtrT>(curNodePtr)->zPos) *
                         (static_cast<partPtrT>(curNodePtr)->mass);
                }
              else
                {
                  cm += static_cast<monoPtrT>(curNodePtr)->mass;
                  cxm += (static_cast<monoPtrT>(curNodePtr)->xCom) *
                         (static_cast<monoPtrT>(curNodePtr)->mass);
                  cym += (static_cast<monoPtrT>(curNodePtr)->yCom) *
                         (static_cast<monoPtrT>(curNodePtr)->mass);
                  czm += (static_cast<monoPtrT>(curNodePtr)->zCom) *
                         (static_cast<monoPtrT>(curNodePtr)->mass);
                }
              goUp();
            }
        }
    }

  // copy data to node itself ...
  if (cm > 0.)
    {
      static_cast<monoPtrT>(curNodePtr)->mass = cm;
      static_cast<monoPtrT>(curNodePtr)->xCom = cxm / cm;
      static_cast<monoPtrT>(curNodePtr)->yCom = cym / cm;
      static_cast<monoPtrT>(curNodePtr)->zCom = czm / cm;
    }
  else
    {
      static_cast<monoPtrT>(curNodePtr)->mass = 0.;
      static_cast<monoPtrT>(curNodePtr)->xCom = 0.;
      static_cast<monoPtrT>(curNodePtr)->yCom = 0.;
      static_cast<monoPtrT>(curNodePtr)->zCom = 0.;
    }
}

///
/// merge multipole moments in remote and local buffers
///
void mergeRemoteCells()
{
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
}

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
#ifdef SPHLATCH_TREE_PROFILE
  calcGravityCellsCounter++;
#endif
  // cellPartDist is already set by the MAC function

  // no softening for cells
  const valueType rInvPow3 = 1.0 / ( cellPartDist * cellPartDist * cellPartDist );

  const valueType rx = curGravParticleX - static_cast<monoPtrT>(curNodePtr)->xCom;
  const valueType ry = curGravParticleY - static_cast<monoPtrT>(curNodePtr)->yCom;
  const valueType rz = curGravParticleZ - static_cast<monoPtrT>(curNodePtr)->zCom;
  const valueType mass = static_cast<monoPtrT>(curNodePtr)->mass;

  // gravity due to monopole term
  curGravParticleAX -= mass * rx * rInvPow3;
  curGravParticleAY -= mass * ry * rInvPow3;
  curGravParticleAZ -= mass * rz * rInvPow3;
}

private:
enum monopoleIndex { CX, CY, CZ, MASS, MSIZE };

private:
};
};

#endif
