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

///
/// calculate multipole from children
/// this function does not check, whether the current node is actually a cell!
///
valueType monopolCM, monopolCXM, monopolCYM, monopolCZM;
valueType quadQ11, quadQ22, quadQ33, quadQ12, quadQ13, quadQ23;
//valueType childDX, childDY, childDZ, childDD;
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
  // butes anything so monopolCM gets 0 and fucks up the center of mass
  // of the cell. so if nobody contributes anything, omit the addition
  // to the cell.
  //
  monopolCM = 0.;
  monopolCXM = 0.;
  monopolCYM = 0.;
  monopolCZM = 0.;

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
                  monopolCM += static_cast<partPtrT>(curNodePtr)->mass;
                  monopolCXM += (static_cast<partPtrT>(curNodePtr)->xPos) *
                                (static_cast<partPtrT>(curNodePtr)->mass);
                  monopolCYM += (static_cast<partPtrT>(curNodePtr)->yPos) *
                                (static_cast<partPtrT>(curNodePtr)->mass);
                  monopolCZM += (static_cast<partPtrT>(curNodePtr)->zPos) *
                                (static_cast<partPtrT>(curNodePtr)->mass);
                }
              else
                {
                  monopolCM += static_cast<quadPtrT>(curNodePtr)->mass;
                  monopolCXM += (static_cast<quadPtrT>(curNodePtr)->xCom) *
                                (static_cast<quadPtrT>(curNodePtr)->mass);
                  monopolCYM += (static_cast<quadPtrT>(curNodePtr)->yCom) *
                                (static_cast<quadPtrT>(curNodePtr)->mass);
                  monopolCZM += (static_cast<quadPtrT>(curNodePtr)->zCom) *
                                (static_cast<quadPtrT>(curNodePtr)->mass);
                }
              goUp();
            }
        }
    }

  /*childDX = 0.;
  childDY = 0.;
  childDZ = 0.;
  childDD = 0.;*/

  /*quadQ11 = 0.;
  quadQ22 = 0.;
  quadQ33 = 0.;
  quadQ12 = 0.;
  quadQ13 = 0.;
  quadQ23 = 0.;*/
  
  /*for (size_t i = 0; i < 8; i++)
    {
      if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i]->isLocal
              == curNodePtr->isLocal)
            {
              goChild(i);
              if (curNodePtr->isParticle == true)
                {
                  childDX = static_cast<partPtrT>(curNodePtr)->xPos
                            - (monopolCXM / monopolCM);
                  childDY = static_cast<partPtrT>(curNodePtr)->yPos
                            - (monopolCYM / monopolCM);
                  childDZ = static_cast<partPtrT>(curNodePtr)->zPos
                            - (monopolCZM / monopolCM);
                  childDD = childDX * childDX +
                    childDY * childDY + childDZ * childDZ;

                  quadQ11 += static_cast<partPtrT>(curNodePtr)->mass *
                             (3 * childDX * childDX - childDD);
                  quadQ22 += static_cast<partPtrT>(curNodePtr)->mass *
                             (3 * childDY * childDY - childDD);
                  quadQ33 += static_cast<partPtrT>(curNodePtr)->mass *
                             (3 * childDZ * childDZ - childDD);

                  quadQ12 += static_cast<partPtrT>(curNodePtr)->mass *
                             (3 * childDX * childDY);
                  quadQ13 += static_cast<partPtrT>(curNodePtr)->mass *
                             (3 * childDX * childDZ);
                  quadQ23 += static_cast<partPtrT>(curNodePtr)->mass *
                             (3 * childDY * childDZ);
                }
              else
                {
                  childDX = ( static_cast<quadPtrT>(curNodePtr)->xCom )
                            - (monopolCXM / monopolCM);
                  childDY = ( static_cast<quadPtrT>(curNodePtr)->yCom )
                            - (monopolCYM / monopolCM);
                  childDZ = ( static_cast<quadPtrT>(curNodePtr)->zCom )
                            - (monopolCZM / monopolCM);

                  childDD = childDX * childDX +
                    childDY * childDY + childDZ * childDZ;

                  quadQ11 += static_cast<quadPtrT>(curNodePtr)->mass *
                             (3 * childDX * childDX - childDD) +
                             static_cast<quadPtrT>(curNodePtr)->q11;
                  quadQ22 += static_cast<quadPtrT>(curNodePtr)->mass *
                             (3 * childDY * childDY - childDD) +
                             static_cast<quadPtrT>(curNodePtr)->q22;
                  quadQ33 += static_cast<quadPtrT>(curNodePtr)->mass *
                             (3 * childDZ * childDZ - childDD) +
                             static_cast<quadPtrT>(curNodePtr)->q33;

                  quadQ12 += static_cast<quadPtrT>(curNodePtr)->mass *
                             (3 * childDX * childDY) +
                             static_cast<quadPtrT>(curNodePtr)->q12;
                  quadQ13 += static_cast<quadPtrT>(curNodePtr)->mass *
                             (3 * childDX * childDZ) +
                             static_cast<quadPtrT>(curNodePtr)->q13;
                  quadQ23 += static_cast<quadPtrT>(curNodePtr)->mass *
                             (3 * childDY * childDZ) +
                             static_cast<quadPtrT>(curNodePtr)->q23;
                }
              goUp();
            }
        }
    }*/

  // copy data to node itself ...
  if (monopolCM > 0.)
    {
      static_cast<quadPtrT>(curNodePtr)->mass = monopolCM;
      static_cast<quadPtrT>(curNodePtr)->xCom = monopolCXM / monopolCM;
      static_cast<quadPtrT>(curNodePtr)->yCom = monopolCYM / monopolCM;
      static_cast<quadPtrT>(curNodePtr)->zCom = monopolCZM / monopolCM;
    }
  else
    {
      static_cast<quadPtrT>(curNodePtr)->mass = 0.;
      static_cast<quadPtrT>(curNodePtr)->xCom = 0.;
      static_cast<quadPtrT>(curNodePtr)->yCom = 0.;
      static_cast<quadPtrT>(curNodePtr)->zCom = 0.;
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
  /*localCells(toptreeCounter, Q11) = static_cast<quadPtrT>(curNodePtr)->q11;
  localCells(toptreeCounter, Q22) = static_cast<quadPtrT>(curNodePtr)->q22;
  localCells(toptreeCounter, Q33) = static_cast<quadPtrT>(curNodePtr)->q33;
  localCells(toptreeCounter, Q12) = static_cast<quadPtrT>(curNodePtr)->q12;
  localCells(toptreeCounter, Q13) = static_cast<quadPtrT>(curNodePtr)->q13;
  localCells(toptreeCounter, Q23) = static_cast<quadPtrT>(curNodePtr)->q23;*/
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
  /*static_cast<quadPtrT>(curNodePtr)->q11 = localCells(toptreeCounter, Q11);
  static_cast<quadPtrT>(curNodePtr)->q22 = localCells(toptreeCounter, Q22);
  static_cast<quadPtrT>(curNodePtr)->q33 = localCells(toptreeCounter, Q33);
  static_cast<quadPtrT>(curNodePtr)->q12 = localCells(toptreeCounter, Q12);
  static_cast<quadPtrT>(curNodePtr)->q13 = localCells(toptreeCounter, Q13);
  static_cast<quadPtrT>(curNodePtr)->q23 = localCells(toptreeCounter, Q23);*/
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
  cellPartDistPow3 = cellPartDist * cellPartDist * cellPartDist;

  curGravParticleAX -= (static_cast<quadPtrT>(curNodePtr)->mass) *
                       (curGravParticleX -
                        static_cast<quadPtrT>(curNodePtr)->xCom) /
                       cellPartDistPow3;
  curGravParticleAY -= (static_cast<quadPtrT>(curNodePtr)->mass) *
                       (curGravParticleY -
                        static_cast<quadPtrT>(curNodePtr)->yCom) /
                       cellPartDistPow3;
  curGravParticleAZ -= (static_cast<quadPtrT>(curNodePtr)->mass) *
                       (curGravParticleZ -
                        static_cast<quadPtrT>(curNodePtr)->zCom) /
                       cellPartDistPow3;
}

private:
enum quadrupoleIndex { CX, CY, CZ, MASS, Q11, Q22, Q33, Q12, Q13, Q23, QSIZE };

private:
};
};

#endif
