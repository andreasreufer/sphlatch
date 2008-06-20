#ifndef BHTREE_GENERIC_H
#define BHTREE_GENERIC_H

/*
 *  bhtree_generic.h
 *
 *
 *  Created by Andreas Reufer on 08.02.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "bhtree_generic_node.h"
#include "bhtree_particle_node.h"
#include "bhtree_cell_node.h"

#ifdef SPHLATCH_PARALLEL
#include "communication_manager.h"
#endif
#include "particle_manager.h"

using namespace sphlatch::vectindices;

namespace sphlatch {
template<class T_leaftype>

class BHtree {
typedef genericNode nodeT;

typedef genericNode* nodePtrT;
typedef particleNode* partPtrT;
typedef genericCellNode* cellPtrT;

#ifdef SPHLATCH_PARALLEL
typedef sphlatch::CommunicationManager commManagerType;
#endif
typedef sphlatch::ParticleManager partManagerType;

private:
T_leaftype& asLeaf()
{
  return static_cast<T_leaftype&>(*this);
}

public:
///
/// \brief constructor:
/// - set theta, gravity constant and toptree depth
/// - check them for sanity
/// - allocate and setup RootNode and point cursor to it
/// - set counters
/// - build toptree
/// - set up buffer matrices for toptree cell nodes
///
BHtree(valueType _thetaMAC,
       valueType _gravConst,
       size_t _czDepth,
       valvectType _rootCenter,
       valueType _rootSize) :
#ifdef SPHLATCH_PARALLEL
  CommManager(commManagerType::instance()),
#endif
  PartManager(partManagerType::instance()),
  pos(PartManager.pos),
  acc(PartManager.acc),
  eps(PartManager.eps),
    m(PartManager.m)
{
  if (_thetaMAC < 0.)
    {
      std::cerr << "error: theta may not be negative! setting theta = 0.6\n";
      _thetaMAC = 0.6;
    }
  else
  if (_thetaMAC > 1.)
    {
      std::cerr << "warning: theta > 1. leads to self-acceleration!\n";
    }

  asLeaf().allocRootNode();

  rootPtr->parent = NULL;
  rootPtr->depth = 0;
  for (size_t i = 0; i < 8; i++)
    {
      static_cast<cellPtrT>(rootPtr)->child[i] = NULL;
    }

  rootPtr->ident = -1;

  static_cast<cellPtrT>(rootPtr)->xCenter = _rootCenter(0);
  static_cast<cellPtrT>(rootPtr)->yCenter = _rootCenter(1);
  static_cast<cellPtrT>(rootPtr)->zCenter = _rootCenter(2);
  static_cast<cellPtrT>(rootPtr)->cellSize = _rootSize;

  rootPtr->isParticle = false;
  rootPtr->isEmpty = true;
  rootPtr->isLocal = false;

  curNodePtr = rootPtr;

  cellCounter = 1;  // we now have the root cell
  partCounter = 0;  // ... and no particles
  ghostCounter = 0;  // ... and no ghost particles

  thetaMAC = _thetaMAC;
  gravConst = _gravConst;
  toptreeDepth = std::max(_czDepth,
                          _czDepth + lrint(ceil(-log(thetaMAC) / log(2.0))));

  if (toptreeDepth > 6)
    {
      std::cerr << "warning: toptree depth of " << toptreeDepth
                << " is pretty big!\n";
    }

  buildToptreeRecursor();
  noToptreeCells = cellCounter;
  noToptreeLeafCells = 1 << (3 * toptreeDepth); // 8^(toptreeDepth)
  noMultipoleMoments = asLeaf().noMultipoleMoments();

#ifdef SPHLATCH_PARALLEL
  localCells.resize(noToptreeLeafCells, noMultipoleMoments);
  localIsFilled.resize(noToptreeLeafCells);

  remoteCells.resize(noToptreeLeafCells, noMultipoleMoments);
  remoteIsFilled.resize(noToptreeLeafCells);

  cellVect.resize(noMultipoleMoments);
#endif
  /// ugly
  neighbourList.resize(100);
  neighDistList.resize(100);

  ///
  /// prepare the list of local particles
  ///
  noLocParts  = PartManager.getNoLocalParts();
  partProxies.resize(noLocParts);
};

///
/// destructor:
/// - postorder recurse node deletion
///
~BHtree(void)
{
  goRoot();
  empty();
  delete curNodePtr; // Seppuku!
}

protected:
#ifdef SPHLATCH_PARALLEL
commManagerType& CommManager;
#endif
partManagerType& PartManager;

protected:
nodePtrT curNodePtr, rootPtr;

///
/// variables
///
size_t cellCounter, partCounter, ghostCounter, toptreeDepth,
       noToptreeCells, noToptreeLeafCells, noMultipoleMoments;

matrixType localCells, remoteCells;
bitsetType localIsFilled, remoteIsFilled;
valvectType cellVect;

matrixRefType pos, acc;
valvectRefType eps, m;

std::vector<nodePtrT> partProxies;
size_t noLocParts;

std::fstream logFile;

/// little helpers
protected:

///
/// go up one level, works for every node
///
inline void goUp()
{
  curNodePtr = curNodePtr->parent;
}

///
/// go to child, works obviously only for cell nodes which may have childs.
/// for performance reasons there are no checks of current node type!
///
inline void goChild(const size_t _n)
{
  curNodePtr = static_cast<cellPtrT>(curNodePtr)->child[_n];
}

///
/// go to root
//
inline void goRoot()
{
  curNodePtr = rootPtr;
}

private:
///
/// create a new cell child in octant _n
///
void newCellChild(const size_t _n)
{
  if (curNodePtr->isParticle == false)
    {
      if (static_cast<cellPtrT>(curNodePtr)->child[_n] == NULL)
        {
          asLeaf().allocNewCellChild(_n);

          goChild(_n);

          curNodePtr->isParticle = false;
          curNodePtr->isEmpty = true;
          curNodePtr->depth = curNodePtr->parent->depth + 1;

          curNodePtr->ident =
            static_cast<identType>(-cellCounter - 1);
          cellCounter++;

          for (size_t i = 0; i < 8; i++)
            {
              static_cast<cellPtrT>(curNodePtr)->child[i] = NULL;
            }

          valueType curCellSize =
            static_cast<cellPtrT>(curNodePtr->parent)->cellSize / 2.;

          static_cast<cellPtrT>(curNodePtr)->cellSize = curCellSize;
          static_cast<cellPtrT>(curNodePtr)->xCenter
          = ((_n) % 2) ?
            static_cast<cellPtrT>(curNodePtr->parent)->xCenter
            + 0.5 * curCellSize :
            static_cast<cellPtrT>(curNodePtr->parent)->xCenter
            - 0.5 * curCellSize;
          static_cast<cellPtrT>(curNodePtr)->yCenter
          = ((_n >> 1) % 2) ?
            static_cast<cellPtrT>(curNodePtr->parent)->yCenter
            + 0.5 * curCellSize :
            static_cast<cellPtrT>(curNodePtr->parent)->yCenter
            - 0.5 * curCellSize;
          static_cast<cellPtrT>(curNodePtr)->zCenter
          = ((_n >> 2) % 2) ?
            static_cast<cellPtrT>(curNodePtr->parent)->zCenter
            + 0.5 * curCellSize :
            static_cast<cellPtrT>(curNodePtr->parent)->zCenter
            - 0.5 * curCellSize;

          goUp();
        }
    }
}

///
/// create a new cell child in octant _n
/// no checks are performed, whether there exist already a children of
/// the current node is actually a cell node
///
void newPartChild(const size_t _n)
{
  particleNode* newNodePtr = new particleNode;

  newNodePtr->parent = curNodePtr;
  static_cast<cellPtrT>(curNodePtr)->child[_n] = newNodePtr;

  goChild(_n);

  curNodePtr->isParticle = true;
  curNodePtr->isEmpty = true;
  curNodePtr->depth = curNodePtr->parent->depth + 1;

  goUp();
}
// end of little helpers


// top tree building stuff
private:
///
/// recursor to build toptree
///
void buildToptreeRecursor(void)
{
  if (atOrBelowToptreeStop())
    {
    }
  else
    {
      for (size_t i = 0; i < 8; i++)  // try without loop
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i] == NULL)
            {
              newCellChild(i);

              goChild(i);

              ///
              /// assume a toptree cell never to be local
              ///
              curNodePtr->isLocal = false;
              buildToptreeRecursor();
              goUp();
            }
        }
    }
}


///
/// stop recursion if:
/// - depth of current node is below toptreeDepth
///
bool atOrBelowToptreeStop(void)
{
  return(curNodePtr->depth >= toptreeDepth);
};
// end of top tree stuff


// insertParticle() stuff
public:
///
/// method to insert particle:
/// - go to root
/// - call the insertion recursor
///
void insertParticle(size_t _newPartIdx)
{
  goRoot();
  insertParticleRecursor(_newPartIdx);
      
  if ( _newPartIdx >= noLocParts ) {
    ghostCounter++;
  } else {
    partCounter++;
  }
}

private:
///
/// recursor for inserting a new particle:
/// try to insert as child of current
/// node. if child is
/// - empty, insert particle. we're done.
/// - a node, go to child and call recursor
/// - a particle, disconnect particle and
///   call recursor for existing and new
///   particle.
///
void insertParticleRecursor(size_t _newPartIdx)
{
  size_t targetOctant = 0;

  targetOctant += ( pos(_newPartIdx, X) <
                   static_cast<cellPtrT>(curNodePtr)->xCenter) ? 0 : 1;
  targetOctant += ( pos(_newPartIdx, Y) <
                   static_cast<cellPtrT>(curNodePtr)->yCenter) ? 0 : 2;
  targetOctant += ( pos(_newPartIdx, Z) <
                   static_cast<cellPtrT>(curNodePtr)->zCenter) ? 0 : 4;

///
/// If targeted child is empty, place the particle there
///
  if (static_cast<cellPtrT>(curNodePtr)->child[targetOctant] == NULL)
    {
      curNodePtr->isEmpty = false;

      newPartChild(targetOctant);
      goChild(targetOctant);

      if ( _newPartIdx >= noLocParts ) {
        curNodePtr->isLocal = false;
      } else {
        curNodePtr->isLocal = true;
        /// save the particle's proxy
        partProxies[_newPartIdx] = curNodePtr;
      }

      /// particle saves its position to node directly
      static_cast<partPtrT>(curNodePtr)->xPos = pos(_newPartIdx, X);
      static_cast<partPtrT>(curNodePtr)->yPos = pos(_newPartIdx, Y);
      static_cast<partPtrT>(curNodePtr)->zPos = pos(_newPartIdx, Z);
      static_cast<partPtrT>(curNodePtr)->mass =   m(_newPartIdx);

      /// ident saves the rowIndex of the particle
      curNodePtr->ident = _newPartIdx;

      goUp();
    }


///
/// ... or if existing child is a node, then try to place the particle
/// as a child of this node
///
  else if (not
           static_cast<cellPtrT>(curNodePtr)->child[targetOctant]->isParticle)
    {
      curNodePtr->isEmpty = false;
      goChild(targetOctant);
      insertParticleRecursor(_newPartIdx);
      goUp();
    }

///
/// ... or if existing child is a particle (ghost/nonghost), then
/// replace it by a new node and try to place the existing two
/// particles as childs of this node
  else if (static_cast<cellPtrT>(curNodePtr)->child[targetOctant]->isParticle)
    {
      /// goto child, save resident particle
      goChild(targetOctant);
      size_t residentIdx = static_cast<partPtrT>(curNodePtr)->ident;
      goUp();

      /// replace particle by cell node
      delete static_cast<cellPtrT>(curNodePtr)->child[targetOctant];
      static_cast<cellPtrT>(curNodePtr)->child[targetOctant] = NULL;
      newCellChild(targetOctant);

      /// and try to insert both particles again
      goChild(targetOctant);
      insertParticleRecursor(residentIdx);
      insertParticleRecursor(_newPartIdx);
      goUp();
    }
  else
    {
    }
}
// end of insertParticle() stuff

// calcMultipoles() stuff
public:
///
/// calculate multipoles:
///  - go to root
///  - call the multipole recursor
///  - exchange toptrees
///
void calcMultipoles(void)
{
  goRoot();
  calcMPbottomtreeRecursor();

  ///
  /// that's the part syncing the toptree leafs level
  ///
#ifdef SPHLATCH_PARALLEL
  globalCombineToptreeLeafs();
#endif

  goRoot();
  calcMPtoptreeRecursor();
}

private:
///
/// recursor for bottomtree multipole calculation
///
void calcMPbottomtreeRecursor(void)
{
  ///
  /// stop recursion if:
  /// - current node is empty
  /// - current node is a particle
  ///
  if (curNodePtr->isEmpty) /// particles are always empty,
                           /// so we can omit this check
    {
    }
  else
    {
      if (curNodePtr->isParticle == false)
        {
          for (size_t i = 0; i < 8; i++)
            {
              if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
                {
                  goChild(i);
                  calcMPbottomtreeRecursor();
                  goUp();
                }
            }

          if (curNodePtr->depth >= toptreeDepth)
            {
              /// determine locality of current node
              detLocality();
              /// calculate multipole moments (delegated to specialization)
              asLeaf().calcMultipole();
            }
        }
    }
};


///
/// determine locality of cell node
/// this function does not check, whether current node is actually a cell!
///
void detLocality(void)
{
  //
  // check locality of current cell node:
  // - if node is in the toptree, the node is local
  // - if node is below toptree, the parent node is local if
  // any child is local ( ... or ... or ... or ... )
  //
  // to check the proper working together of the costzone and
  // the toptree would be to check, that cells below toptreeDepth
  // either have only local or only non-local children, but never
  // both.
  //
  for (size_t i = 0; i < 8; i++)
    {
      if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
        {
          (curNodePtr->isLocal) |=
            static_cast<cellPtrT>(curNodePtr)->child[i]->isLocal;

          (curNodePtr->isEmpty) &=
            static_cast<cellPtrT>(curNodePtr)->child[i]->isEmpty;
        }
    }
}

void toptreeLeafsToBufferRecursor(void)
{
  ///
  /// stop recursion and copy cell to buffer, if:
  /// - toptreeDepth is reached
  ///
  if (curNodePtr->depth == toptreeDepth)
    {
      localIsFilled[toptreeLeafCounter] = (!(curNodePtr->isEmpty) && curNodePtr->isLocal);
      asLeaf().cellToVect(cellVect);
      particleRowType(localCells, toptreeLeafCounter) = cellVect;
      toptreeLeafCounter++;
    }
  else
    {
      for (size_t i = 0; i < 8; i++)
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
            {
              goChild(i);
              toptreeLeafsToBufferRecursor();
              goUp();
            }
        }
    }
};

void bufferToToptreeLeafsRecursor(void)
{
  ///
  /// stop recursion and copy buffer to cell, if:
  /// - toptreeDepth is reached
  ///
  if (curNodePtr->depth == toptreeDepth)
    {
      curNodePtr->isEmpty = !localIsFilled[toptreeLeafCounter];
      curNodePtr->isLocal = true;
      cellVect = particleRowType(localCells, toptreeLeafCounter);
      asLeaf().vectToCell(cellVect);
      toptreeLeafCounter++;
    }
  else
    {
      for (size_t i = 0; i < 8; i++)
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
            {
              goChild(i);
              bufferToToptreeLeafsRecursor();
              goUp();
            }
        }
    }
};

///
/// recursor for bottomtree multipole calculation
///
void calcMPtoptreeRecursor(void)
{
  ///
  /// stop recursion if:
  /// - toptreeDepth is reached
  ///
  if (curNodePtr->depth == toptreeDepth)
    {
    }
  else
    {
      for (size_t i = 0; i < 8; i++)
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
            {
              goChild(i);
              calcMPtoptreeRecursor();
              goUp();
            }
        }

      if (curNodePtr->depth < toptreeDepth)
        {
          /// determine locality of current node
          detLocality();
          /// calculate multipole moments (delegated to specialization)
          asLeaf().calcMultipole();
        }
    }
};



///
/// sum up the multipole cellData matrix globally
/// and distribute it again to each node. if MPI is not
/// defined, nothing is done here.
///
protected:
size_t toptreeLeafCounter;
private:
void globalCombineToptreeLeafs()
{
#ifdef SPHLATCH_PARALLEL
  const size_t myDomain = CommManager.getMyDomain();
  const size_t noDomains = CommManager.getNoDomains();

  ///
  /// magic algorithm which prepares sending and receiving queues
  /// for the summing up step and sending and receiving stacks for
  /// the distributing step.
  ///
  size_t round = 0;
  size_t remNodes = noDomains;

  std::queue<size_t> sumUpSend, sumUpRecv;
  std::stack<size_t> distrSend, distrRecv;
  while (remNodes > 1)
    {
      size_t noPairs = lrint(floor(remNodes / 2.));
      remNodes -= noPairs;
      size_t stepToNext = (1 << round);

      for (size_t i = 0; i < 2 * noPairs; i += 2)
        {
          size_t sendDomain = (noDomains - 1) - stepToNext * (i + 1);
          size_t recvDomain = (noDomains - 1) - stepToNext * i;

          if (myDomain == sendDomain)
            {
              sumUpSend.push(recvDomain);
              distrRecv.push(recvDomain);
            }
          else if (myDomain == recvDomain)
            {
              sumUpRecv.push(sendDomain);
              distrSend.push(sendDomain);
            }
          else
            {
            }
        }
      round++;
    }

  /// toptree bottom to buffers
  goRoot();
  toptreeLeafCounter = 0;
  toptreeLeafsToBufferRecursor();

  ///
  /// receive multipoles from other nodes and add
  /// them to local value
  ///
  while (!sumUpRecv.empty())
    {
      size_t recvFrom = sumUpRecv.front();
      sumUpRecv.pop();

      CommManager.recvBitset(remoteIsFilled, recvFrom);
      CommManager.recvMatrix(remoteCells, recvFrom);

      for (size_t i = 0; i < noToptreeLeafCells; i++)
        {
          if (remoteIsFilled[i])
            {
              particleRowType(localCells, i) = particleRowType(remoteCells, i);
            }
        }

      /// localIsFilled = localIsFilled or remoteIsFilled
      localIsFilled |= remoteIsFilled;
    }

  ///
  /// send local value to another node
  ///
  while (!sumUpSend.empty())
    {
      size_t sendTo = sumUpSend.front();
      sumUpSend.pop();

      CommManager.sendBitset(localIsFilled, sendTo);
      CommManager.sendMatrix(localCells, sendTo);
    }

  //
  // receive global result
  //
  while (!distrRecv.empty())
    {
      size_t recvFrom = distrRecv.top();
      distrRecv.pop();

      CommManager.recvBitset(localIsFilled, recvFrom);
      CommManager.recvMatrix(localCells, recvFrom);
    }

  //
  // ... and distribute it to other nodes
  //
  while (!distrSend.empty())
    {
      size_t sendTo = distrSend.top();
      distrSend.pop();

      CommManager.sendBitset(localIsFilled, sendTo);
      CommManager.sendMatrix(localCells, sendTo);
    }

  /// buffers to toptree bottom
  goRoot();
  toptreeLeafCounter = 0;
  bufferToToptreeLeafsRecursor();
#endif
}

// end of multipole stuff

// getParticleOrder() stuff
public:
partsIndexVectType particleOrder;
void detParticleOrder()
{
  particleOrder.resize(partCounter);
  partCounter = 0;

  goRoot();
  particleOrderRecursor();
}

private:
void particleOrderRecursor()
{
  if (curNodePtr->isParticle && curNodePtr->isLocal)
    {
      particleOrder[partCounter] = curNodePtr->ident;
      partCounter++;
    }
  else
  if (curNodePtr->isEmpty)
    {
    }
  else
    {
      for (size_t i = 0; i < 8; i++)  // try without loop
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
            {
              goChild(i);
              particleOrderRecursor();
              goUp();
            }
        }
    }
}
// end of getParticleOrder() stuff

// calcGravity() stuff
protected:
valueType curGravParticleX, curGravParticleY, curGravParticleZ;
valueType curGravParticleAX, curGravParticleAY, curGravParticleAZ;
valueType cellPartDist; //, cellPartDistPow3;
valueType thetaMAC, gravConst;
valueType epsilonSquare;
size_t calcGravityCellsCounter, calcGravityPartsCounter;

public:
///
/// calculate gravitation for a particle:
/// - load current particle data
/// - go to root
/// - call the gravity calculation recursor
/// - write back resulting acceleration
///

// what about a ref?
//void calcGravity(const partProxyPtrT _curParticle)
void calcGravity(const size_t _partIdx)
{
  curGravParticleX = pos(_partIdx, X);
  curGravParticleY = pos(_partIdx, Y);
  curGravParticleZ = pos(_partIdx, Z);

  curGravParticleAX = 0.;
  curGravParticleAY = 0.;
  curGravParticleAZ = 0.;

  epsilonSquare = eps(_partIdx)*eps(_partIdx);

#ifdef SPHLATCH_TREE_PROFILE
  calcGravityPartsCounter = 0;
  calcGravityCellsCounter = 0;
#endif

  ///
  /// trick: hide the current particle by letting it look like
  /// an empty cell node, so that it doesn't gravitate with
  /// itself. btw: a particle is always empty.
  ///
  partProxies[_partIdx]->isParticle = false;

  goRoot();
  calcGravityRecursor();

  acc(_partIdx, X) += gravConst * curGravParticleAX;
  acc(_partIdx, Y) += gravConst * curGravParticleAY;
  acc(_partIdx, Z) += gravConst * curGravParticleAZ;

  ///
  /// undo the trick above
  ///
  partProxies[_partIdx]->isParticle = true;
}

private:
void calcGravityRecursor(void)
{
  if (curNodePtr->isParticle)
    {
      calcGravParticle();
#ifdef SPHLATCH_TREE_PROFILE
      calcGravityPartsCounter++;
#endif
    }
  else
  if (curNodePtr->isEmpty)
    {
    }
  else
  if (asLeaf().calcGravMAC())
    {
      asLeaf().calcGravCell();
#ifdef SPHLATCH_TREE_PROFILE
      calcGravityCellsCounter++;
#endif
    }
  else
    {
      for (size_t i = 0; i < 8; i++)
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
            {
              goChild(i);
              calcGravityRecursor();
              goUp();
            }
        }
    }
}

///
/// calculate acceleration due to a particle
/// no check whether current node is actually a particle!
///
void calcGravParticle()
{
  const valueType partGravPartnerX = static_cast<partPtrT>(curNodePtr)->xPos;
  const valueType partGravPartnerY = static_cast<partPtrT>(curNodePtr)->yPos;
  const valueType partGravPartnerZ = static_cast<partPtrT>(curNodePtr)->zPos;
  const valueType partGravPartnerM = static_cast<partPtrT>(curNodePtr)->mass;

  const valueType partPartDistPow2 = (partGravPartnerX - curGravParticleX) *
                                     (partGravPartnerX - curGravParticleX) +
                                     (partGravPartnerY - curGravParticleY) *
                                     (partGravPartnerY - curGravParticleY) +
                                     (partGravPartnerZ - curGravParticleZ) *
                                     (partGravPartnerZ - curGravParticleZ);

  /// \todo: include spline softening
#ifdef SPHLATCH_SINGLEPREC
  const valueType cellPartDistPow3 =
    (partPartDistPow2 + epsilonSquare) *
    sqrtf(partPartDistPow2 + epsilonSquare);
#else
  const valueType cellPartDistPow3 =
    (partPartDistPow2 + epsilonSquare) *
    sqrt(partPartDistPow2 + epsilonSquare);
#endif
  //const valueType cellPartDistPow3 = partPartDistPow2 + epsilonSquare;

  // todo: use cellPartDistPow3Inv
  curGravParticleAX -= partGravPartnerM *
                       (curGravParticleX - partGravPartnerX) /
                       cellPartDistPow3;
  curGravParticleAY -= partGravPartnerM *
                       (curGravParticleY - partGravPartnerY) /
                       cellPartDistPow3;
  curGravParticleAZ -= partGravPartnerM *
                       (curGravParticleZ - partGravPartnerZ) /
                       cellPartDistPow3;
}
// end of calcGravity() stuff

// neighbour search stuff
private:
size_t noNeighbours;
// cellPartDist already definded (great, feels like FORTRAN :-)
valueType curNeighParticleX, curNeighParticleY, curNeighParticleZ,
          searchRadius, searchRadPow2, //cellCornerDist, cellPartDistPow2,
          partPartDistPow2;

public:
///
/// find neighbours:
/// - load current particle data
/// - go to current particle node
/// - go up until every neighbour is within current cell
///   ( 2h sphere around part. fully within cell )
/// - call the neighbour search recursor
/// - sort out cells which do not touch 2h sphere
/// - add particles to list
/// - brute force sort out non-neighbours
/// - return neighbours
///
partsIndexVectType neighbourList;
std::valarray<valueType> neighDistList;
void findNeighbours(size_t _curPartIdx,
                    const valueType _search_radius)
{
  ///
  /// set up some vars and go to particles node
  ///
  curNodePtr = partProxies[_curPartIdx];
  
  curNeighParticleX = static_cast<partPtrT>(curNodePtr)->xPos;
  curNeighParticleY = static_cast<partPtrT>(curNodePtr)->yPos;
  curNeighParticleZ = static_cast<partPtrT>(curNodePtr)->zPos;

  searchRadius = _search_radius;
  searchRadPow2 = _search_radius * _search_radius;
  noNeighbours = 0;

  ///
  /// go to parent cell. while the search sphere is not completely
  /// inside the cell and we still can go up, go up
  ///
  goUp();
  while (not sphereTotInsideCell() && curNodePtr->parent != NULL)
    {
      goUp();
    }
  findNeighbourRecursor();

  neighbourList[0] = noNeighbours;
}


private:
void findNeighbourRecursor()
{
  if (curNodePtr->isParticle)
    {
      partPartDistPow2 =
        (static_cast<partPtrT>(curNodePtr)->xPos - curNeighParticleX) *
        (static_cast<partPtrT>(curNodePtr)->xPos - curNeighParticleX) +
        (static_cast<partPtrT>(curNodePtr)->yPos - curNeighParticleY) *
        (static_cast<partPtrT>(curNodePtr)->yPos - curNeighParticleY) +
        (static_cast<partPtrT>(curNodePtr)->zPos - curNeighParticleZ) *
        (static_cast<partPtrT>(curNodePtr)->zPos - curNeighParticleZ);
      if (partPartDistPow2 < searchRadPow2)
        {
          ///
          /// add the particle to the neighbour list
          ///
          noNeighbours++;
          neighbourList[noNeighbours] = curNodePtr->ident;
          neighDistList[noNeighbours] = sqrt(partPartDistPow2);
        }
    }
  else
  if (curNodePtr->isEmpty)
    {
    }
  else
  ///
  /// a check whether the cell is fully inside the tree can be per-
  /// formed here. for a SPH neighbour search this almost never happens
  /// and therefore the overhead of the check is bigger than the possible
  /// gain in speed
  ///
  if (cellTotOutsideSphere())
    {
    }
  else
    {
      for (size_t i = 0; i < 8; i++)
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
            {
              goChild(i);
              findNeighbourRecursor();
              goUp();
            }
        }
    }
}

///
/// is the search sphere completely INSIDE the cell?
///
bool sphereTotInsideCell()
{
  if ((curNeighParticleX + searchRadius >
       static_cast<cellPtrT>(curNodePtr)->xCenter +
       0.5 * (static_cast<cellPtrT>(curNodePtr)->cellSize)) ||
      (curNeighParticleX - searchRadius <
       static_cast<cellPtrT>(curNodePtr)->xCenter -
       0.5 * (static_cast<cellPtrT>(curNodePtr)->cellSize)))
    {
      return false;
    }
  else
  if ((curNeighParticleY + searchRadius >
       static_cast<cellPtrT>(curNodePtr)->yCenter +
       0.5 * (static_cast<cellPtrT>(curNodePtr)->cellSize)) ||
      (curNeighParticleY - searchRadius <
       static_cast<cellPtrT>(curNodePtr)->yCenter -
       0.5 * (static_cast<cellPtrT>(curNodePtr)->cellSize)))
    {
      return false;
    }
  else
  if ((curNeighParticleZ + searchRadius >
       static_cast<cellPtrT>(curNodePtr)->zCenter +
       0.5 * (static_cast<cellPtrT>(curNodePtr)->cellSize)) ||
      (curNeighParticleZ - searchRadius <
       static_cast<cellPtrT>(curNodePtr)->zCenter -
       0.5 * (static_cast<cellPtrT>(curNodePtr)->cellSize)))
    {
      return false;
    }
  else
    {
      return true;
    }
}


///
/// is the cell completely OUTSIDE the search sphere?
///
bool cellTotOutsideSphere()
{
  const valueType cellPartDistPow2 =
    (static_cast<cellPtrT>(curNodePtr)->xCenter - curNeighParticleX) *
    (static_cast<cellPtrT>(curNodePtr)->xCenter - curNeighParticleX)
    + (static_cast<cellPtrT>(curNodePtr)->yCenter - curNeighParticleY) *
    (static_cast<cellPtrT>(curNodePtr)->yCenter - curNeighParticleY)
    + (static_cast<cellPtrT>(curNodePtr)->zCenter - curNeighParticleZ) *
    (static_cast<cellPtrT>(curNodePtr)->zCenter - curNeighParticleZ);
  const valueType cellCornerDist = static_cast<valueType>(0.866025403784439) *
                                   static_cast<cellPtrT>(curNodePtr)->cellSize;

  return(cellPartDistPow2 > (searchRadius + cellCornerDist) *
         (searchRadius + cellCornerDist));
}
// end of neighbour search stuff

// treeDOTDump() stuff
protected:
std::fstream dumpFile;
///
/// dump the tree as a
/// dot file
///
public:
void treeDOTDump(std::string _dumpFileName)
{
  dumpFile.open(_dumpFileName.c_str(), std::ios::out);

  goRoot();
  dumpFile << "digraph graphname { \n";
  treeDOTDumpRecursor();
  dumpFile << "}\n";
  dumpFile.close();
}

private:
///
/// treeDOTDump() recursor
///
void treeDOTDumpRecursor()
{
  dumpFile << curNodePtr->ident
           << " [";
  dumpFile << "label=\"\",";

  if (curNodePtr->isParticle)
    {
      //dumpFile << "label=\"\",";
      //dumpFile << "label=\"" << lrint( static_cast<partPtrT>(curNodePtr)->mass ) << "\",";
      dumpFile << "shape=circle";
    }
  else
    {
      //dumpFile << "label=\"\",";
      //dumpFile << "label=\"" << lrint( static_cast<monopoleCellNode*>(curNodePtr)->mass)
      //         << "\",";
      //dumpFile << "label=\"" << curNodePtr->ident << ": "
      //         << lrint( static_cast<monopoleCellNode*>(curNodePtr)->mass)  << "\",";
      dumpFile << "shape=box";
    }

  if ((!curNodePtr->isEmpty) or
      curNodePtr->isParticle)
    {
      dumpFile << ",style=filled";
    }

  if (curNodePtr->isLocal)
    {
      // local particle
      if (curNodePtr->isParticle)
        {
          dumpFile << ",color=green";
        }
      else
        {
          // local cell
          if (curNodePtr->depth <= toptreeDepth)
            {
              dumpFile << ",color=blue";
            }
          else
            {
              dumpFile << ",color=red";
            }
        }
    }
  else
    {
      // non-local
      dumpFile << ",color=grey";
    }

  dumpFile << ",width=0.1,height=0.1];\n";

  for (size_t i = 0; i < 8; i++)                          // try without loop
    {
      if (curNodePtr->isParticle == false)
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
            {
              dumpFile << curNodePtr->ident
                       << " -> "
                       << static_cast<cellPtrT>(curNodePtr)->child[i]->ident
                       << "\n";
              goChild(i);
              treeDOTDumpRecursor();
              goUp();
            }
        }
    }
}
// end of treeDOTDump() stuff


// treeDump() stuff
public:
void treeDump(std::string _dumpFileName)
{
  dumpFile.open(_dumpFileName.c_str(), std::ios::out);
  goRoot();
  treeDumpRecursor();
  dumpFile.close();
}

private:
void treeDumpRecursor()
{
  if (curNodePtr->isParticle)
    {
      dumpFile << "P";
    }
  else
    {
      dumpFile << "N";
    }
  if (curNodePtr->isEmpty)
    {
      dumpFile << "E";
    }
  else
    {
      dumpFile << "F";
    }
  if (curNodePtr->isLocal)
    {
      dumpFile << "L";
    }
  else
    {
      dumpFile << "R";
    }
  dumpFile << " ";

  dumpFile << std::fixed << std::right << std::setw(6);
  dumpFile << curNodePtr->ident;
  dumpFile << "   ";
  dumpFile << curNodePtr->depth;
  dumpFile << "   ";
  dumpFile << std::setprecision(9) << std::setw(15);
  dumpFile << std::scientific;

  if (curNodePtr->isParticle)
    {
      for (size_t i = 0; i < 4; i++)
        {
          dumpFile << 0. << "   ";
        }

      dumpFile << static_cast<partPtrT>(curNodePtr)->xPos << "   ";
      dumpFile << static_cast<partPtrT>(curNodePtr)->yPos << "   ";
      dumpFile << static_cast<partPtrT>(curNodePtr)->zPos << "   ";
      dumpFile << 0. << "   ";
      dumpFile << static_cast<partPtrT>(curNodePtr)->mass << "   ";

      for (size_t i = 4; i < noMultipoleMoments; i++) // no of >monopole terms
        {
          dumpFile << 0. << "   ";
        }
    }
  else
    {
      dumpFile << static_cast<cellPtrT>(curNodePtr)->xCenter << "   ";
      dumpFile << static_cast<cellPtrT>(curNodePtr)->yCenter << "   ";
      dumpFile << static_cast<cellPtrT>(curNodePtr)->zCenter << "   ";
      dumpFile << static_cast<cellPtrT>(curNodePtr)->cellSize << "   ";

      asLeaf().cellToVect(cellVect);
      for (size_t i = 0; i < noMultipoleMoments; i++)
        {
          dumpFile << cellVect[i] << "   ";
        }
    }

  dumpFile << "\n";

  for (size_t i = 0; i < 8; i++)
    {
      if (curNodePtr->isParticle == false)
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
            {
              goChild(i);
              treeDumpRecursor();
              goUp();
            }
        }
    }
}
// end of treeDump() stuff

// toptreeDump() stuff
public:
///
/// dump the toptree
///
void toptreeDump(std::string _dumpFileName)
{
  dumpFile.open(_dumpFileName.c_str(), std::ios::out);

  dumpFile << "# toptree depth is " << toptreeDepth
           << " which gives us " << noToptreeCells
           << " toptree cells\n"
           << "#ident empty depth       xCenter            yCenter   "
           << "        zCenter         sidelength            xCom    "
           << "          yCom              zCom               mass\n";
  goRoot();
  toptreeDumpRecursor();
  dumpFile.close();
}

private:
///
/// toptreeDump() recursor
///
void toptreeDumpRecursor()
{
  dumpFile << std::fixed << std::right << std::setw(6);
  dumpFile << curNodePtr->ident;
  dumpFile << "   ";
  //dumpFile << ! curNodePtr->isEmpty ;
  //dumpFile << "   ";
  dumpFile << curNodePtr->depth;
  dumpFile << "   ";
  dumpFile << std::setprecision(15) << std::setw(20);
  dumpFile << std::scientific;

  dumpFile << static_cast<cellPtrT>(curNodePtr)->xCenter << "   ";
  dumpFile << static_cast<cellPtrT>(curNodePtr)->yCenter << "   ";
  dumpFile << static_cast<cellPtrT>(curNodePtr)->zCenter << "   ";
  dumpFile << static_cast<cellPtrT>(curNodePtr)->cellSize << "   ";

  asLeaf().cellToVect(cellVect);
  for (size_t i = 0; i < noMultipoleMoments; i++)
    {
      dumpFile << cellVect[i] << "   ";
    }

  dumpFile << "\n";

  if (atOrBelowToptreeStop())
    {
    }
  else
    {
      for (size_t i = 0; i < 8; i++)             // try without loop
        {
          if (curNodePtr->isParticle == false)
            {
              if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
                {
                  goChild(i);
                  toptreeDumpRecursor();
                  goUp();
                }
            }
        }
    }
}
// end of toptreeDump() stuff

// empty() stuff
///
/// deletes the current subtree
///
private:
void empty(void)
{
  emptyRecursor();
}

///
/// recursor for emptying the
/// tree
///
void emptyRecursor()
{
  for (size_t i = 0; i < 8; i++)                          // try without loop
    {
      if (curNodePtr->isParticle == false)
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
            {
              goChild(i);
              emptyRecursor();
              goUp();

              delete static_cast<cellPtrT>(curNodePtr)->child[i];
              static_cast<cellPtrT>(curNodePtr)->child[i] = NULL;
            }
        }
    }
};
// end of empty() stuff
};
};

#endif
