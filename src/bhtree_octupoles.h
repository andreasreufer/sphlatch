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

  static valueType rx, ry, rz, rr;
  rx = 0.;
  ry = 0.;
  rz = 0.;
  rr = 0.;

  static valueType massC;
  massC = 0.;

  static valueType q11, q22, q33, q12, q13, q23;
  q11 = 0.;
  q22 = 0.;
  q33 = 0.;
  q12 = 0.;
  q13 = 0.;
  q23 = 0.;

  static valueType rxrx, ryry, rzrz;
  rxrx = 0.;
  ryry = 0.;
  rzrz = 0.;

  static valueType q11C, q22C, q33C, q12C, q13C, q23C;
  q11C = 0.;
  q22C = 0.;
  q33C = 0.;
  q12C = 0.;
  q13C = 0.;
  q23C = 0.;

  static valueType s11, s22, s33, s12, s21, s13, s31, s23, s32, s123;
  s11 = 0.;
  s22 = 0.;
  s33 = 0.;
  s12 = 0.;
  s21 = 0.;
  s13 = 0.;
  s31 = 0.;
  s23 = 0.;
  s32 = 0.;
  s123 = 0.;

  static valueType s11C, s22C, s33C, s12C, s21C, s13C, s31C, s23C, s32C, s123C;
  s11C = 0.;
  s22C = 0.;
  s33C = 0.;
  s12C = 0.;
  s21C = 0.;
  s13C = 0.;
  s31C = 0.;
  s23C = 0.;
  s32C = 0.;
  s123C = 0.;

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
                      massC = static_cast<partPtrT>(curNodePtr)->mass;

                      rr = rx * rx + ry * ry + rz * rz;
                      rxrx = rx * rx;
                      ryry = ry * ry;
                      rzrz = rz * rz;

                      q11 += ((3. * rx * rx - rr) * massC);
                      q22 += ((3. * ry * ry - rr) * massC);
                      q33 += ((3. * rz * rz - rr) * massC);

                      q12 += (3. * rx * ry * massC);
                      q13 += (3. * rx * rz * massC);
                      q23 += (3. * ry * rz * massC);

                      s11 += (5. * rxrx - 3. * rr) * rx * massC;
                      s22 += (5. * ryry - 3. * rr) * ry * massC;
                      s33 += (5. * rzrz - 3. * rr) * rz * massC;

                      s12 += (15. * rxrx - 3. * rr) * ry * massC;
                      s21 += (15. * ryry - 3. * rr) * rx * massC;

                      s13 += (15. * rxrx - 3. * rr) * rz * massC;
                      s31 += (15. * rzrz - 3. * rr) * rx * massC;

                      s23 += (15. * ryry - 3. * rr) * rz * massC;
                      s32 += (15. * rzrz - 3. * rr) * ry * massC;

                      s123 += 15. * rx * ry * rz * massC;
                    }
                  else
                    {
                      rx = (static_cast<octuPtrT>(curNodePtr)->xCom)
                           - (cxm / cm);
                      ry = (static_cast<octuPtrT>(curNodePtr)->yCom)
                           - (cym / cm);
                      rz = (static_cast<octuPtrT>(curNodePtr)->zCom)
                           - (czm / cm);
                      massC = static_cast<octuPtrT>(curNodePtr)->mass;

                      rr = rx * rx + ry * ry + rz * rz;
                      rxrx = rx * rx;
                      ryry = ry * ry;
                      rzrz = rz * rz;

                      q11C = static_cast<octuPtrT>(curNodePtr)->q11;
                      q22C = static_cast<octuPtrT>(curNodePtr)->q22;
                      q33C = static_cast<octuPtrT>(curNodePtr)->q33;
                      q12C = static_cast<octuPtrT>(curNodePtr)->q12;
                      q13C = static_cast<octuPtrT>(curNodePtr)->q13;
                      q23C = static_cast<octuPtrT>(curNodePtr)->q23;

                      s11C = static_cast<octuPtrT>(curNodePtr)->s11;
                      s22C = static_cast<octuPtrT>(curNodePtr)->s22;
                      s33C = static_cast<octuPtrT>(curNodePtr)->s33;
                      s12C = static_cast<octuPtrT>(curNodePtr)->s12;
                      s21C = static_cast<octuPtrT>(curNodePtr)->s21;
                      s13C = static_cast<octuPtrT>(curNodePtr)->s13;
                      s31C = static_cast<octuPtrT>(curNodePtr)->s31;
                      s23C = static_cast<octuPtrT>(curNodePtr)->s23;
                      s32C = static_cast<octuPtrT>(curNodePtr)->s32;
                      s123C = static_cast<octuPtrT>(curNodePtr)->s123;

                      q11 += ((3. * rx * rx - rr) * massC) + q11C;
                      q22 += ((3. * ry * ry - rr) * massC) + q22C;
                      q33 += ((3. * rz * rz - rr) * massC) + q33C;

                      q12 += (3. * rx * ry * massC) + q12C;
                      q13 += (3. * rx * rz * massC) + q13C;
                      q23 += (3. * ry * rz * massC) + q23C;

                      s11 += (5. * rxrx - 3. * rr) * rx * massC
                             + q11C * rx * (3. / 2.)
                             - q12C * ry - q13C * rz
                             + s11C;
                      s22 += (5. * ryry - 3. * rr) * ry * massC
                             + q22C * ry * (3. / 2.)
                             - q12C * rx - q23C * rz
                             + s22C;
                      s33 += (5. * rzrz - 3. * rr) * rz * massC
                             + q33C * rz * (3. / 2.)
                             - q13C * rx - q23C * ry
                             + s33C;

                      s12 += (15. * rxrx - 3. * rr) * ry * massC
                             + 4. * rx * q12C + (5. / 2.) * ry * q11C
                             - ry * q22C - rz * q23C + s12C;
                      s21 += (15. * ryry - 3. * rr) * rx * massC
                             + 4. * ry * q12C + (5. / 2.) * rx * q22C
                             - rx * q11C - rz * q13C + s21C;

                      s13 += (15. * rxrx - 3. * rr) * rz * massC
                             + 4. * rx * q13C + (5. / 2.) * rz * q11C
                             - ry * q23C - rz * q33C + s13C;
                      s31 += (15. * rzrz - 3. * rr) * rx * massC
                             + 4. * rz * q13C + (5. / 2.) * rx * q33C
                             - rx * q11C - ry * q12C + s31C;

                      s23 += (15. * ryry - 3. * rr) * rz * massC
                             + 4. * ry * q23C + (5. / 2.) * rz * q22C
                             - rx * q13C - rz * q33C + s23C;
                      s32 += (15. * rzrz - 3. * rr) * ry * massC
                             + 4. * rz * q23C + (5. / 2.) * ry * q33C
                             - rx * q12C - ry * q22C + s32C;

                      s123 += 15. * rx * ry * rz * massC
                              + 25. * (q12C * rz + q13C * ry + q23C * rx)
                              + 15. * s123C;
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

      static_cast<octuPtrT>(curNodePtr)->s11 = s11;
      static_cast<octuPtrT>(curNodePtr)->s22 = s22;
      static_cast<octuPtrT>(curNodePtr)->s33 = s33;
      static_cast<octuPtrT>(curNodePtr)->s12 = s12;
      static_cast<octuPtrT>(curNodePtr)->s21 = s21;
      static_cast<octuPtrT>(curNodePtr)->s13 = s13;
      static_cast<octuPtrT>(curNodePtr)->s31 = s31;
      static_cast<octuPtrT>(curNodePtr)->s23 = s23;
      static_cast<octuPtrT>(curNodePtr)->s32 = s32;
      static_cast<octuPtrT>(curNodePtr)->s123 = s123;
    }
  else
    {
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

  localCells(toptreeCounter, S11) = static_cast<octuPtrT>(curNodePtr)->s11;
  localCells(toptreeCounter, S22) = static_cast<octuPtrT>(curNodePtr)->s22;
  localCells(toptreeCounter, S33) = static_cast<octuPtrT>(curNodePtr)->s33;
  localCells(toptreeCounter, S12) = static_cast<octuPtrT>(curNodePtr)->s12;
  localCells(toptreeCounter, S21) = static_cast<octuPtrT>(curNodePtr)->s21;
  localCells(toptreeCounter, S13) = static_cast<octuPtrT>(curNodePtr)->s13;
  localCells(toptreeCounter, S31) = static_cast<octuPtrT>(curNodePtr)->s31;
  localCells(toptreeCounter, S23) = static_cast<octuPtrT>(curNodePtr)->s23;
  localCells(toptreeCounter, S32) = static_cast<octuPtrT>(curNodePtr)->s32;
  localCells(toptreeCounter, S123) = static_cast<octuPtrT>(curNodePtr)->s123;
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

  static_cast<octuPtrT>(curNodePtr)->s11 = localCells(toptreeCounter, S11);
  static_cast<octuPtrT>(curNodePtr)->s22 = localCells(toptreeCounter, S22);
  static_cast<octuPtrT>(curNodePtr)->s33 = localCells(toptreeCounter, S33);
  static_cast<octuPtrT>(curNodePtr)->s12 = localCells(toptreeCounter, S12);
  static_cast<octuPtrT>(curNodePtr)->s21 = localCells(toptreeCounter, S21);
  static_cast<octuPtrT>(curNodePtr)->s13 = localCells(toptreeCounter, S13);
  static_cast<octuPtrT>(curNodePtr)->s31 = localCells(toptreeCounter, S31);
  static_cast<octuPtrT>(curNodePtr)->s23 = localCells(toptreeCounter, S23);
  static_cast<octuPtrT>(curNodePtr)->s32 = localCells(toptreeCounter, S32);
  static_cast<octuPtrT>(curNodePtr)->s123 = localCells(toptreeCounter, S123);
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

  dumpFile << static_cast<octuPtrT>(curNodePtr)->s11 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->s22 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->s33 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->s12 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->s21 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->s13 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->s31 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->s23 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->s32 << "   ";
  dumpFile << static_cast<octuPtrT>(curNodePtr)->s123 << "   ";
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
  const valueType cellPartDistPow3 = cellPartDist * cellPartDist * cellPartDist;

  const valueType rx = curGravParticleX - static_cast<octuPtrT>(curNodePtr)->xCom;
  const valueType ry = curGravParticleY - static_cast<octuPtrT>(curNodePtr)->yCom;
  const valueType rz = curGravParticleZ - static_cast<octuPtrT>(curNodePtr)->zCom;
  const valueType mass = static_cast<octuPtrT>(curNodePtr)->mass;

  // gravity due to monopole term
  curGravParticleAX -= mass * rx / cellPartDistPow3;
  curGravParticleAY -= mass * ry / cellPartDistPow3;
  curGravParticleAZ -= mass * rz / cellPartDistPow3;

  // intermediate results for quadrupole term
  const valueType cellPartDistPow5 = cellPartDistPow3 * cellPartDist * cellPartDist;
  const valueType cellPartDistPow7 = cellPartDistPow5 * cellPartDist * cellPartDist;

  const valueType q11 = static_cast<octuPtrT>(curNodePtr)->q11;
  const valueType q22 = static_cast<octuPtrT>(curNodePtr)->q22;
  const valueType q33 = static_cast<octuPtrT>(curNodePtr)->q33;
  const valueType q12 = static_cast<octuPtrT>(curNodePtr)->q12;
  const valueType q13 = static_cast<octuPtrT>(curNodePtr)->q13;
  const valueType q23 = static_cast<octuPtrT>(curNodePtr)->q23;

  // less calcs, more mem access. good idea? do some benchmarking
  const valueType rxrx = rx * rx;
  const valueType ryry = ry * ry;
  const valueType rzrz = rz * rz;

  const valueType q1jrj = q11 * rx + q12 * ry + q13 * rz;
  const valueType q2jrj = q12 * rx + q22 * ry + q23 * rz;
  const valueType q3jrj = q13 * rx + q23 * ry + q33 * rz;
  const valueType qijrirj = q11 * rxrx +
                            q22 * ryry +
                            q33 * rzrz +
                            2. * q12 * rx * ry +
                            2. * q13 * rx * rz +
                            2. * q23 * ry * rz;

  // gravity due to octupole term
  curGravParticleAX -= (-1.0 * q1jrj / cellPartDistPow5)
                       + (2.5 * qijrirj * rx / cellPartDistPow7);
  curGravParticleAY -= (-1.0 * q2jrj / cellPartDistPow5)
                       + (2.5 * qijrirj * ry / cellPartDistPow7);
  curGravParticleAZ -= (-1.0 * q3jrj / cellPartDistPow5)
                       + (2.5 * qijrirj * rz / cellPartDistPow7);

  // intermediate results for octupole term
  const valueType cellPartDistPow9 = cellPartDistPow7 * cellPartDist * cellPartDist;

  const valueType s11 = static_cast<octuPtrT>(curNodePtr)->s11;
  const valueType s22 = static_cast<octuPtrT>(curNodePtr)->s22;
  const valueType s33 = static_cast<octuPtrT>(curNodePtr)->s33;
  const valueType s12 = static_cast<octuPtrT>(curNodePtr)->s12;
  const valueType s21 = static_cast<octuPtrT>(curNodePtr)->s21;
  const valueType s13 = static_cast<octuPtrT>(curNodePtr)->s13;
  const valueType s31 = static_cast<octuPtrT>(curNodePtr)->s31;
  const valueType s23 = static_cast<octuPtrT>(curNodePtr)->s23;
  const valueType s32 = static_cast<octuPtrT>(curNodePtr)->s32;
  const valueType s123 = static_cast<octuPtrT>(curNodePtr)->s123;

  const valueType s1jrj = s11 * rx + s12 * ry + s13 * rz;
  const valueType s2jrj = s21 * rx + s22 * ry + s23 * rz;
  const valueType s3jrj = s31 * rx + s32 * ry + s33 * rz;

  const valueType si1riri = s11 * rxrx + s21 * ryry + s31 * rzrz;
  const valueType si2riri = s12 * rxrx + s22 * ryry + s32 * rzrz;
  const valueType si3riri = s13 * rxrx + s23 * ryry + s33 * rzrz;

  const valueType sijririrj = si1riri * rx + si2riri * ry + si3riri * rz;

  const valueType rxryrz = rx * ry * rz;

  // gravity due to octupole terms
  /*curGravParticleAX -= (1. / cellPartDistPow7) * (-s1jrj * rx - 0.5 * si1riri - 0.5 * s123 * ry * rz)
                       + (1. / cellPartDistPow9) * (7.0 * sijririrj * rx + 3.5 * s123 * rxryrz * rx);*/
  /*curGravParticleAY -= (1. / cellPartDistPow7) * (-s2jrj * ry - 0.5 * si2riri - 0.5 * s123 * rz * rx)
                       + (1. / cellPartDistPow9) * (7.0 * sijririrj * ry + 3.5 * s123 * rxryrz * ry);*/
  curGravParticleAZ -= (1. / cellPartDistPow7) * (-s3jrj * rz - 0.5 * si3riri - 0.5 * s123 * rx * ry)
                       + (1. / cellPartDistPow9) * (7.0 * sijririrj * rz + 3.5 * s123 * rxryrz * rz);
}

private:
enum octupoleIndex { CX, CY, CZ, MASS, Q11, Q22, Q33, Q12, Q13, Q23,
                     S11, S22, S33, S12, S21, S13, S31, S23, S32, S123, OSIZE };

private:
};
};

#endif
