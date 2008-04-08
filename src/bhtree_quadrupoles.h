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
/// report number of multipole moments
/// (includes center of mass and mass)
///
size_t noMultipoleMoments()
{
  return QSIZE; ///  3 center of mass + 1 monopole
                /// + 6 quadrupoles = 10 multipoles
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
  static valvectType curCell(QSIZE), childCell(QSIZE);

  for (size_t i = 0; i < QSIZE; i++)
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

  if (curCell[MASS] > 0.)
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
    }

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
  /// const static does not seem to work here
  const valueType oldMass = _target[MASS];

  _target[MASS] += _source[MASS];

  if (_target[MASS] > 0.)
    {
      /// const static does not seem to work here
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
  const valueType rx = _source[CX] - _target[CX];
  const valueType ry = _source[CY] - _target[CY];
  const valueType rz = _source[CZ] - _target[CZ];
  const valueType rxrx = rx * rx;
  const valueType ryry = ry * ry;
  const valueType rzrz = rz * rz;
  const valueType rr = rxrx + ryry + rzrz;

  _target[Q11] += (3. * rxrx - rr) * _source[MASS] + _source[Q11];
  _target[Q22] += (3. * ryry - rr) * _source[MASS] + _source[Q22];
  _target[Q33] += (3. * rzrz - rr) * _source[MASS] + _source[Q33];

  _target[Q12] += (3. * rx * ry * _source[MASS]) + _source[Q12];
  _target[Q13] += (3. * rx * rz * _source[MASS]) + _source[Q13];
  _target[Q23] += (3. * ry * rz * _source[MASS]) + _source[Q23];
}

///
/// copy current cell node to a vector
/// no checks are performed, whether current node is a
/// cell nor if vector has the right size
///
void cellToVect(valvectRefType _vect)
{
  _vect[CX] = static_cast<quadPtrT>(curNodePtr)->xCom;
  _vect[CY] = static_cast<quadPtrT>(curNodePtr)->yCom;
  _vect[CZ] = static_cast<quadPtrT>(curNodePtr)->zCom;
  _vect[MASS] = static_cast<quadPtrT>(curNodePtr)->mass;

  _vect[Q11] = static_cast<quadPtrT>(curNodePtr)->q11;
  _vect[Q22] = static_cast<quadPtrT>(curNodePtr)->q22;
  _vect[Q33] = static_cast<quadPtrT>(curNodePtr)->q33;
  _vect[Q12] = static_cast<quadPtrT>(curNodePtr)->q12;
  _vect[Q13] = static_cast<quadPtrT>(curNodePtr)->q13;
  _vect[Q23] = static_cast<quadPtrT>(curNodePtr)->q23;
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

  _vect[Q11] = 0.;
  _vect[Q22] = 0.;
  _vect[Q33] = 0.;
  _vect[Q12] = 0.;
  _vect[Q13] = 0.;
  _vect[Q23] = 0.;
}

///
/// copy vector to current cell node
/// no checks are performed, whether current node is a
/// cell nor if vector has the right size
///
void vectToCell(valvectRefType _vect)
{
  static_cast<quadPtrT>(curNodePtr)->xCom = _vect[CX];
  static_cast<quadPtrT>(curNodePtr)->yCom = _vect[CY];
  static_cast<quadPtrT>(curNodePtr)->zCom = _vect[CZ];
  static_cast<quadPtrT>(curNodePtr)->mass = _vect[MASS];

  static_cast<quadPtrT>(curNodePtr)->q11 = _vect[Q11];
  static_cast<quadPtrT>(curNodePtr)->q22 = _vect[Q22];
  static_cast<quadPtrT>(curNodePtr)->q33 = _vect[Q33];
  static_cast<quadPtrT>(curNodePtr)->q12 = _vect[Q12];
  static_cast<quadPtrT>(curNodePtr)->q13 = _vect[Q13];
  static_cast<quadPtrT>(curNodePtr)->q23 = _vect[Q23];
}


///
/// stop recursion if:
/// - current node is empty << ??
/// - MAC is fulfilled
///
bool calcGravMAC(void)
{
  /// any way to speed this up?
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
  /// cellPartDist is already set by the MAC function

  /// intermediate results for monopole term
  const valueType rInvPow3 = 1. / (cellPartDist * cellPartDist * cellPartDist);

  const valueType rx = curGravParticleX
                       - static_cast<quadPtrT>(curNodePtr)->xCom;
  const valueType ry = curGravParticleY
                       - static_cast<quadPtrT>(curNodePtr)->yCom;
  const valueType rz = curGravParticleZ
                       - static_cast<quadPtrT>(curNodePtr)->zCom;
  const valueType mass = static_cast<quadPtrT>(curNodePtr)->mass;

  /// gravity due to monopole term
  curGravParticleAX -= (rInvPow3) * mass * rx;
  curGravParticleAY -= (rInvPow3) * mass * ry;
  curGravParticleAZ -= (rInvPow3) * mass * rz;

  /// intermediate results for quadrupole term
  const valueType rInvPow5 = rInvPow3 / (cellPartDist * cellPartDist);
  const valueType rInvPow7 = rInvPow5 / (cellPartDist * cellPartDist);

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

  /// gravity due to quadrupole term
  curGravParticleAX += (rInvPow5) * (q1jrj) - (rInvPow7) * (2.5 * qijrirj * rx);
  curGravParticleAY += (rInvPow5) * (q2jrj) - (rInvPow7) * (2.5 * qijrirj * ry);
  curGravParticleAZ += (rInvPow5) * (q3jrj) - (rInvPow7) * (2.5 * qijrirj * rz);
}

private:
enum quadrupoleIndex { CX, CY, CZ, MASS, Q11, Q22, Q33, Q12, Q13, Q23, QSIZE };

private:
};
};

#endif
