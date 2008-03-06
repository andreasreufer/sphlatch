#ifndef BHTREE_OCTUPOLES_H
#define BHTREE_OCTUPOLES_H

/*
 *  bhtree_octupoles.h
 *
 *
 *  Created by Andreas Reufer on 07.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "bhtree_octupole_node.h"

namespace sphlatch {
class Octupoles : public BHtree<Octupoles>
{
typedef genericNode* nodePtrT;
typedef particleNode* partPtrT;
typedef genericCellNode* cellPtrT;
typedef octupoleCellNode* octuPtrT;
typedef particleProxy* partProxyPtrT;

public:
///
/// allocates the root monopole cell node
///
void allocRootNode()
{
  rootPtr = new octupoleCellNode;
}

///
/// resize buffers for communication
///
void prepareBuffers()
{
  localCells.resize(noToptreeCells, OSIZE);
  localIsFilled.resize(noToptreeCells);

  remoteCells.resize(noToptreeCells, OSIZE);
  remoteIsFilled.resize(noToptreeCells);
}

///
/// allocates a new monopole cell node and connects it as child _n
/// no check is performed, whether curNodePtr points to a cell node!
///
void allocNewCellChild(const size_t _n)
{
  // allocate new cell node
  octuPtrT newNodePtr =
    new octupoleCellNode;

  // connect the new cell node to curNodePtr
  newNodePtr->parent = curNodePtr;
  static_cast<cellPtrT>(curNodePtr)->child[_n] = newNodePtr;

  // set cell vars to zero
  static_cast<octuPtrT>(newNodePtr)->mass = 0.;
  static_cast<octuPtrT>(newNodePtr)->xCom = 0.;
  static_cast<octuPtrT>(newNodePtr)->yCom = 0.;
  static_cast<octuPtrT>(newNodePtr)->zCom = 0.;
  
  static_cast<octuPtrT>(newNodePtr)->q11 = 0.;
  static_cast<octuPtrT>(newNodePtr)->q22 = 0.;
  static_cast<octuPtrT>(newNodePtr)->q33 = 0.;
  static_cast<octuPtrT>(newNodePtr)->q12 = 0.;
  static_cast<octuPtrT>(newNodePtr)->q13 = 0.;
  static_cast<octuPtrT>(newNodePtr)->q23 = 0.;
  
  static_cast<octuPtrT>(newNodePtr)->s11 = 0.;
  static_cast<octuPtrT>(newNodePtr)->s22 = 0.;
  static_cast<octuPtrT>(newNodePtr)->s33 = 0.;
  static_cast<octuPtrT>(newNodePtr)->s12 = 0.;
  static_cast<octuPtrT>(newNodePtr)->s21 = 0.;
  static_cast<octuPtrT>(newNodePtr)->s13 = 0.;
  static_cast<octuPtrT>(newNodePtr)->s31 = 0.;
  static_cast<octuPtrT>(newNodePtr)->s23 = 0.;
  static_cast<octuPtrT>(newNodePtr)->s32 = 0.;
  static_cast<octuPtrT>(newNodePtr)->s123 = 0.;
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
                  cm += static_cast<octuPtrT>(curNodePtr)->mass;
                  cxm += (static_cast<octuPtrT>(curNodePtr)->xCom) *
                         (static_cast<octuPtrT>(curNodePtr)->mass);
                  cym += (static_cast<octuPtrT>(curNodePtr)->yCom) *
                         (static_cast<octuPtrT>(curNodePtr)->mass);
                  czm += (static_cast<octuPtrT>(curNodePtr)->zCom) *
                         (static_cast<octuPtrT>(curNodePtr)->mass);
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
                      rx = (static_cast<octuPtrT>(curNodePtr)->xCom)
                           - (cxm / cm);
                      ry = (static_cast<octuPtrT>(curNodePtr)->yCom)
                           - (cym / cm);
                      rz = (static_cast<octuPtrT>(curNodePtr)->zCom)
                           - (czm / cm);

                      rr = rx * rx + ry * ry + rz * rz;

                      q11 += (static_cast<octuPtrT>(curNodePtr)->mass *
                              (3 * rx * rx - rr)) +
                             static_cast<octuPtrT>(curNodePtr)->q11;
                      q22 += (static_cast<octuPtrT>(curNodePtr)->mass *
                              (3 * ry * ry - rr)) +
                             static_cast<octuPtrT>(curNodePtr)->q22;
                      q33 += (static_cast<octuPtrT>(curNodePtr)->mass *
                              (3 * rz * rz - rr)) +
                             static_cast<octuPtrT>(curNodePtr)->q33;

                      q12 += static_cast<octuPtrT>(curNodePtr)->mass *
                             (3 * rx * ry) +
                             static_cast<octuPtrT>(curNodePtr)->q12;
                      q13 += static_cast<octuPtrT>(curNodePtr)->mass *
                             (3 * rx * rz) +
                             static_cast<octuPtrT>(curNodePtr)->q13;
                      q23 += static_cast<octuPtrT>(curNodePtr)->mass *
                             (3 * ry * rz) +
                             static_cast<octuPtrT>(curNodePtr)->q23;
                    }
                  goUp();
                }
            }
        }
    }

  // copy data to node itself ...
  if (cm > 0.)
    {
      static_cast<octuPtrT>(curNodePtr)->mass = cm;
      static_cast<octuPtrT>(curNodePtr)->xCom = cxm / cm;
      static_cast<octuPtrT>(curNodePtr)->yCom = cym / cm;
      static_cast<octuPtrT>(curNodePtr)->zCom = czm / cm;
      static_cast<octuPtrT>(curNodePtr)->q11 = q11;
      static_cast<octuPtrT>(curNodePtr)->q22 = q22;
      static_cast<octuPtrT>(curNodePtr)->q33 = q33;
      static_cast<octuPtrT>(curNodePtr)->q12 = q12;
      static_cast<octuPtrT>(curNodePtr)->q13 = q13;
      static_cast<octuPtrT>(curNodePtr)->q23 = q23;
    }
  else
    {
      static_cast<octuPtrT>(curNodePtr)->mass = 0.;
      static_cast<octuPtrT>(curNodePtr)->xCom = 0.;
      static_cast<octuPtrT>(curNodePtr)->yCom = 0.;
      static_cast<octuPtrT>(curNodePtr)->zCom = 0.;
      static_cast<octuPtrT>(curNodePtr)->q11 = 0.;
      static_cast<octuPtrT>(curNodePtr)->q22 = 0.;
      static_cast<octuPtrT>(curNodePtr)->q33 = 0.;
      static_cast<octuPtrT>(curNodePtr)->q12 = 0.;
      static_cast<octuPtrT>(curNodePtr)->q13 = 0.;
      static_cast<octuPtrT>(curNodePtr)->q23 = 0.;
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
  localCells(toptreeCounter, CX) = static_cast<octuPtrT>(curNodePtr)->xCom;
  localCells(toptreeCounter, CY) = static_cast<octuPtrT>(curNodePtr)->yCom;
  localCells(toptreeCounter, CZ) = static_cast<octuPtrT>(curNodePtr)->zCom;
  localCells(toptreeCounter, MASS) = static_cast<octuPtrT>(curNodePtr)->mass;
  localCells(toptreeCounter, Q11) = static_cast<octuPtrT>(curNodePtr)->q11;
  localCells(toptreeCounter, Q22) = static_cast<octuPtrT>(curNodePtr)->q22;
  localCells(toptreeCounter, Q33) = static_cast<octuPtrT>(curNodePtr)->q33;
  localCells(toptreeCounter, Q12) = static_cast<octuPtrT>(curNodePtr)->q12;
  localCells(toptreeCounter, Q13) = static_cast<octuPtrT>(curNodePtr)->q13;
  localCells(toptreeCounter, Q23) = static_cast<octuPtrT>(curNodePtr)->q23;
}

///
/// copy buffer to current cell node
///
void bufferToCell()
{
  static_cast<octuPtrT>(curNodePtr)->xCom = localCells(toptreeCounter, CX);
  static_cast<octuPtrT>(curNodePtr)->yCom = localCells(toptreeCounter, CY);
  static_cast<octuPtrT>(curNodePtr)->zCom = localCells(toptreeCounter, CZ);
  static_cast<octuPtrT>(curNodePtr)->mass = localCells(toptreeCounter, MASS);
  static_cast<octuPtrT>(curNodePtr)->q11 = localCells(toptreeCounter, Q11);
  static_cast<octuPtrT>(curNodePtr)->q22 = localCells(toptreeCounter, Q22);
  static_cast<octuPtrT>(curNodePtr)->q33 = localCells(toptreeCounter, Q33);
  static_cast<octuPtrT>(curNodePtr)->q12 = localCells(toptreeCounter, Q12);
  static_cast<octuPtrT>(curNodePtr)->q13 = localCells(toptreeCounter, Q13);
  static_cast<octuPtrT>(curNodePtr)->q23 = localCells(toptreeCounter, Q23);
}

///
/// report current cell node to dumpFile stream
///
void reportMultipoles()
{
  dumpFile << static_cast<octuPtrT>(curNodePtr)->xCom << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->yCom << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->zCom << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->mass << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->q11 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->q22 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->q33 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->q12 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->q13 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->q23 << "   ";
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
    (static_cast<octuPtrT>(curNodePtr)->xCom - curGravParticleX) *
    (static_cast<octuPtrT>(curNodePtr)->xCom - curGravParticleX) +
    (static_cast<octuPtrT>(curNodePtr)->yCom - curGravParticleY) *
    (static_cast<octuPtrT>(curNodePtr)->yCom - curGravParticleY) +
    (static_cast<octuPtrT>(curNodePtr)->zCom - curGravParticleZ) *
    (static_cast<octuPtrT>(curNodePtr)->zCom - curGravParticleZ)
    );
  return(((static_cast<octuPtrT>(curNodePtr)->cellSize)
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
  valueType cellPartDistPow3;
  cellPartDistPow3 = cellPartDist * cellPartDist * cellPartDist;

  valueType rx, ry, rz, mass;
  rx = curGravParticleX - static_cast<octuPtrT>(curNodePtr)->xCom;
  ry = curGravParticleY - static_cast<octuPtrT>(curNodePtr)->yCom;
  rz = curGravParticleZ - static_cast<octuPtrT>(curNodePtr)->zCom;
  mass = static_cast<octuPtrT>(curNodePtr)->mass;

  // gravity due to monopole term
  /*curGravParticleAX -= mass * rx / cellPartDistPow3;
  curGravParticleAY -= mass * ry / cellPartDistPow3;
  curGravParticleAZ -= mass * rz / cellPartDistPow3;*/

  valueType cellPartDistPow5, cellPartDistPow7;
  cellPartDistPow5 = cellPartDistPow3 * cellPartDist * cellPartDist;
  cellPartDistPow7 = cellPartDistPow5 * cellPartDist * cellPartDist;

  valueType q11, q22, q33, q12, q13, q23;
  q11 = static_cast<octuPtrT>(curNodePtr)->q11;
  q22 = static_cast<octuPtrT>(curNodePtr)->q22;
  q33 = static_cast<octuPtrT>(curNodePtr)->q33;
  q12 = static_cast<octuPtrT>(curNodePtr)->q12;
  q13 = static_cast<octuPtrT>(curNodePtr)->q13;
  q23 = static_cast<octuPtrT>(curNodePtr)->q23;

  valueType q1jrj, q2jrj, q3jrj, qijrirj;
  q1jrj = q11 * rx + q12 * ry + q13 * rz;
  q2jrj = q12 * rx + q22 * ry + q23 * rz;
  q3jrj = q13 * rx + q23 * ry + q33 * rz;
  qijrirj = q11 * rx * rx +
            q22 * ry * ry +
            q33 * rz * rz +
            2. * q12 * rx * ry +
            2. * q13 * rx * rz +
            2. * q23 * ry * rz;

  // gravity due to qquadruupole term
  curGravParticleX -= (0.5 * q1jrj / cellPartDistPow5)
                      + (2.5 * qijrirj * rx / cellPartDistPow7);
  curGravParticleY -= (0.5 * q2jrj / cellPartDistPow5)
                      + (2.5 * qijrirj * ry / cellPartDistPow7);
  curGravParticleZ -= (0.5 * q3jrj / cellPartDistPow5)
                      + (2.5 * qijrirj * rz / cellPartDistPow7);
}

private:
enum octupoleIndex { CX, CY, CZ, MASS, Q11, Q22, Q33, Q12, Q13, Q23, 
                     S11, S22, S33, S12, S21, S13, S31, S23, S123, OSIZE };

private:
};
};

#endif
