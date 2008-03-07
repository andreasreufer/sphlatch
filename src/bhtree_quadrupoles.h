#ifndef BHTREE_QUADRUPOLES_H
#define BHTREE_QUADRUPOLES_H

/*
 *  bhtree_quadrupoles.h
 *
 *
 *  Created by Andreas Reufer on 07.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "bhtree_quadrupole_node.h"

namespace sphlatch {
class Quadrupoles : public BHtree<Quadrupoles>
{
typedef genericNode* nodePtrT;
typedef particleNode* partPtrT;
typedef genericCellNode* cellPtrT;
typedef quadrupoleCellNode* quadPtrT;
typedef particleProxy* partProxyPtrT;

public:
///
/// allocates the root monopole cell node
///
void allocRootNode()
{
  rootPtr = new quadrupoleCellNode;
}

///
/// resize buffers for communication
///
void prepareBuffers()
{
  localCells.resize(noToptreeCells, QSIZE);
  localIsFilled.resize(noToptreeCells);

  remoteCells.resize(noToptreeCells, QSIZE);
  remoteIsFilled.resize(noToptreeCells);
}

///
/// allocates a new monopole cell node and connects it as child _n
/// no check is performed, whether curNodePtr points to a cell node!
///
void allocNewCellChild(const size_t _n)
{
  // allocate new cell node
  quadPtrT newNodePtr =
    new quadrupoleCellNode;

  // connect the new cell node to curNodePtr
  newNodePtr->parent = curNodePtr;
  static_cast<cellPtrT>(curNodePtr)->child[_n] = newNodePtr;

  // set cell vars to zero
  static_cast<quadPtrT>(newNodePtr)->mass = 0.;
  static_cast<quadPtrT>(newNodePtr)->xCom = 0.;
  static_cast<quadPtrT>(newNodePtr)->yCom = 0.;
  static_cast<quadPtrT>(newNodePtr)->zCom = 0.;

  static_cast<quadPtrT>(newNodePtr)->q11 = 0.;
  static_cast<quadPtrT>(newNodePtr)->q22 = 0.;
  static_cast<quadPtrT>(newNodePtr)->q33 = 0.;
  static_cast<quadPtrT>(newNodePtr)->q12 = 0.;
  static_cast<quadPtrT>(newNodePtr)->q13 = 0.;
  static_cast<quadPtrT>(newNodePtr)->q23 = 0.;
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
                  cm += static_cast<quadPtrT>(curNodePtr)->mass;
                  cxm += (static_cast<quadPtrT>(curNodePtr)->xCom) *
                         (static_cast<quadPtrT>(curNodePtr)->mass);
                  cym += (static_cast<quadPtrT>(curNodePtr)->yCom) *
                         (static_cast<quadPtrT>(curNodePtr)->mass);
                  czm += (static_cast<quadPtrT>(curNodePtr)->zCom) *
                         (static_cast<quadPtrT>(curNodePtr)->mass);
                }
              goUp();
            }
        }
    }

  static valueType q11, q22, q33, q12, q13, q23;
  q11 = 0.;
  q22 = 0.;
  q33 = 0.;
  q12 = 0.;
  q13 = 0.;
  q23 = 0.;

  static valueType rx, ry, rz, rr;
  rx = 0.;
  ry = 0.;
  rz = 0.;
  rr = 0.;

  if (cm > 0.)
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
                      rx = static_cast<partPtrT>(curNodePtr)->xPos
                           - (cxm / cm);
                      ry = static_cast<partPtrT>(curNodePtr)->yPos
                           - (cym / cm);
                      rz = static_cast<partPtrT>(curNodePtr)->zPos
                           - (czm / cm);
                      rr = rx * rx + ry * ry + rz * rz;

                      q11 += static_cast<partPtrT>(curNodePtr)->mass *
                             (3 * rx * rx - rr);
                      q22 += static_cast<partPtrT>(curNodePtr)->mass *
                             (3 * ry * ry - rr);
                      q33 += static_cast<partPtrT>(curNodePtr)->mass *
                             (3 * rz * rz - rr);

                      q12 += static_cast<partPtrT>(curNodePtr)->mass *
                             (3 * rx * ry);
                      q13 += static_cast<partPtrT>(curNodePtr)->mass *
                             (3 * rx * rz);
                      q23 += static_cast<partPtrT>(curNodePtr)->mass *
                             (3 * ry * rz);
                    }
                  else
                    {
                      rx = (static_cast<quadPtrT>(curNodePtr)->xCom)
                           - (cxm / cm);
                      ry = (static_cast<quadPtrT>(curNodePtr)->yCom)
                           - (cym / cm);
                      rz = (static_cast<quadPtrT>(curNodePtr)->zCom)
                           - (czm / cm);

                      rr = rx * rx + ry * ry + rz * rz;

                      q11 += (static_cast<quadPtrT>(curNodePtr)->mass *
                              (3 * rx * rx - rr)) +
                             static_cast<quadPtrT>(curNodePtr)->q11;
                      q22 += (static_cast<quadPtrT>(curNodePtr)->mass *
                              (3 * ry * ry - rr)) +
                             static_cast<quadPtrT>(curNodePtr)->q22;
                      q33 += (static_cast<quadPtrT>(curNodePtr)->mass *
                              (3 * rz * rz - rr)) +
                             static_cast<quadPtrT>(curNodePtr)->q33;

                      q12 += static_cast<quadPtrT>(curNodePtr)->mass *
                             (3 * rx * ry) +
                             static_cast<quadPtrT>(curNodePtr)->q12;
                      q13 += static_cast<quadPtrT>(curNodePtr)->mass *
                             (3 * rx * rz) +
                             static_cast<quadPtrT>(curNodePtr)->q13;
                      q23 += static_cast<quadPtrT>(curNodePtr)->mass *
                             (3 * ry * rz) +
                             static_cast<quadPtrT>(curNodePtr)->q23;
                    }
                  goUp();
                }
            }
        }
    }

  // copy data to node itself ...
  if (cm > 0.)
    {
      static_cast<quadPtrT>(curNodePtr)->mass = cm;
      static_cast<quadPtrT>(curNodePtr)->xCom = cxm / cm;
      static_cast<quadPtrT>(curNodePtr)->yCom = cym / cm;
      static_cast<quadPtrT>(curNodePtr)->zCom = czm / cm;
      static_cast<quadPtrT>(curNodePtr)->q11 = q11;
      static_cast<quadPtrT>(curNodePtr)->q22 = q22;
      static_cast<quadPtrT>(curNodePtr)->q33 = q33;
      static_cast<quadPtrT>(curNodePtr)->q12 = q12;
      static_cast<quadPtrT>(curNodePtr)->q13 = q13;
      static_cast<quadPtrT>(curNodePtr)->q23 = q23;
    }
  else
    {
      static_cast<quadPtrT>(curNodePtr)->mass = 0.;
      static_cast<quadPtrT>(curNodePtr)->xCom = 0.;
      static_cast<quadPtrT>(curNodePtr)->yCom = 0.;
      static_cast<quadPtrT>(curNodePtr)->zCom = 0.;
      static_cast<quadPtrT>(curNodePtr)->q11 = 0.;
      static_cast<quadPtrT>(curNodePtr)->q22 = 0.;
      static_cast<quadPtrT>(curNodePtr)->q33 = 0.;
      static_cast<quadPtrT>(curNodePtr)->q12 = 0.;
      static_cast<quadPtrT>(curNodePtr)->q13 = 0.;
      static_cast<quadPtrT>(curNodePtr)->q23 = 0.;
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
  localCells(toptreeCounter, CX) = static_cast<quadPtrT>(curNodePtr)->xCom;
  localCells(toptreeCounter, CY) = static_cast<quadPtrT>(curNodePtr)->yCom;
  localCells(toptreeCounter, CZ) = static_cast<quadPtrT>(curNodePtr)->zCom;
  localCells(toptreeCounter, MASS) = static_cast<quadPtrT>(curNodePtr)->mass;
  localCells(toptreeCounter, Q11) = static_cast<quadPtrT>(curNodePtr)->q11;
  localCells(toptreeCounter, Q22) = static_cast<quadPtrT>(curNodePtr)->q22;
  localCells(toptreeCounter, Q33) = static_cast<quadPtrT>(curNodePtr)->q33;
  localCells(toptreeCounter, Q12) = static_cast<quadPtrT>(curNodePtr)->q12;
  localCells(toptreeCounter, Q13) = static_cast<quadPtrT>(curNodePtr)->q13;
  localCells(toptreeCounter, Q23) = static_cast<quadPtrT>(curNodePtr)->q23;
}

///
/// copy buffer to current cell node
///
void bufferToCell()
{
  static_cast<quadPtrT>(curNodePtr)->xCom = localCells(toptreeCounter, CX);
  static_cast<quadPtrT>(curNodePtr)->yCom = localCells(toptreeCounter, CY);
  static_cast<quadPtrT>(curNodePtr)->zCom = localCells(toptreeCounter, CZ);
  static_cast<quadPtrT>(curNodePtr)->mass = localCells(toptreeCounter, MASS);
  static_cast<quadPtrT>(curNodePtr)->q11 = localCells(toptreeCounter, Q11);
  static_cast<quadPtrT>(curNodePtr)->q22 = localCells(toptreeCounter, Q22);
  static_cast<quadPtrT>(curNodePtr)->q33 = localCells(toptreeCounter, Q33);
  static_cast<quadPtrT>(curNodePtr)->q12 = localCells(toptreeCounter, Q12);
  static_cast<quadPtrT>(curNodePtr)->q13 = localCells(toptreeCounter, Q13);
  static_cast<quadPtrT>(curNodePtr)->q23 = localCells(toptreeCounter, Q23);
}

///
/// report current cell node to dumpFile stream
///
void reportMultipoles()
{
  dumpFile << static_cast<quadPtrT>(curNodePtr)->xCom << "   ";
  dumpFile << static_cast<quadPtrT>(curNodePtr)->yCom << "   ";
  dumpFile << static_cast<quadPtrT>(curNodePtr)->zCom << "   ";
  dumpFile << static_cast<quadPtrT>(curNodePtr)->mass << "   ";
  dumpFile << static_cast<quadPtrT>(curNodePtr)->q11 << "   ";
  dumpFile << static_cast<quadPtrT>(curNodePtr)->q22 << "   ";
  dumpFile << static_cast<quadPtrT>(curNodePtr)->q33 << "   ";
  dumpFile << static_cast<quadPtrT>(curNodePtr)->q12 << "   ";
  dumpFile << static_cast<quadPtrT>(curNodePtr)->q13 << "   ";
  dumpFile << static_cast<quadPtrT>(curNodePtr)->q23 << "   ";
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
    (static_cast<quadPtrT>(curNodePtr)->xCom - curGravParticleX) *
    (static_cast<quadPtrT>(curNodePtr)->xCom - curGravParticleX) +
    (static_cast<quadPtrT>(curNodePtr)->yCom - curGravParticleY) *
    (static_cast<quadPtrT>(curNodePtr)->yCom - curGravParticleY) +
    (static_cast<quadPtrT>(curNodePtr)->zCom - curGravParticleZ) *
    (static_cast<quadPtrT>(curNodePtr)->zCom - curGravParticleZ)
    );
  return(((static_cast<quadPtrT>(curNodePtr)->cellSize)
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
  const valueType cellPartDistPow3 = cellPartDist * cellPartDist * cellPartDist;

  const valueType rx = curGravParticleX - static_cast<quadPtrT>(curNodePtr)->xCom;
  const valueType ry = curGravParticleY - static_cast<quadPtrT>(curNodePtr)->yCom;
  const valueType rz = curGravParticleZ - static_cast<quadPtrT>(curNodePtr)->zCom;
  const valueType mass = static_cast<quadPtrT>(curNodePtr)->mass;

  // gravity due to monopole term
  curGravParticleAX -= mass * rx / cellPartDistPow3;
  curGravParticleAY -= mass * ry / cellPartDistPow3;
  curGravParticleAZ -= mass * rz / cellPartDistPow3;

  const valueType cellPartDistPow5 = cellPartDistPow3 * cellPartDist * cellPartDist;
  const valueType cellPartDistPow7 = cellPartDistPow5 * cellPartDist * cellPartDist;

  const valueType q11 = static_cast<quadPtrT>(curNodePtr)->q11;
  const valueType q22 = static_cast<quadPtrT>(curNodePtr)->q22;
  const valueType q33 = static_cast<quadPtrT>(curNodePtr)->q33;
  const valueType q12 = static_cast<quadPtrT>(curNodePtr)->q12;
  const valueType q13 = static_cast<quadPtrT>(curNodePtr)->q13;
  const valueType q23 = static_cast<quadPtrT>(curNodePtr)->q23;

  const valueType q1jrj = q11 * rx + q12 * ry + q13 * rz;
  const valueType q2jrj = q12 * rx + q22 * ry + q23 * rz;
  const valueType q3jrj = q13 * rx + q23 * ry + q33 * rz;
  const valueType qijrirj = q11 * rx * rx +
                            q22 * ry * ry +
                            q33 * rz * rz +
                            2. * q12 * rx * ry +
                            2. * q13 * rx * rz +
                            2. * q23 * ry * rz;

  // gravity due to quadrupole term
  curGravParticleAX -= (-1.0 * q1jrj / cellPartDistPow5)
                       + (2.5 * qijrirj * rx / cellPartDistPow7);
  curGravParticleAY -= (-1.0 * q2jrj / cellPartDistPow5)
                       + (2.5 * qijrirj * ry / cellPartDistPow7);
  curGravParticleAZ -= (-1.0 * q3jrj / cellPartDistPow5)
                       + (2.5 * qijrirj * rz / cellPartDistPow7);
}

private:
enum quadrupoleIndex { CX, CY, CZ, MASS, Q11, Q22, Q33, Q12, Q13, Q23, QSIZE };

private:
};
};

#endif
