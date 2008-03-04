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
#include "bhtree_particle_proxy.h"
#include "bhtree_cell_node.h"

#include "particle.h"

namespace sphlatch {
template<class T_leaftype>
class BHtree {
typedef genericNode nodeT;

typedef genericNode* nodePtrT;
typedef particleNode* partPtrT;
typedef genericCellNode* cellPtrT;
typedef particleProxy* partProxyPtrT;

private:
T_leaftype& asLeaf()
{
  return static_cast<T_leaftype&>(*this);
}

public:
/** \brief constructor:
 * - set theta, gravity constant and toptree depth
 * - check them for sanity
 * - allocate and setup RootNode and point cursor to it
 * - set counters
 * - build toptree
 * - set up buffer matrices for toptree cell nodes
 */
BHtree(valueType _thetaMAC,
       valueType _gravConst,
       size_t _czDepth,
       valvectType _rootCenter,
       valueType _rootSize)
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
  rootPtr->isLocal = true;

  curNodePtr = rootPtr;

  cellCounter = 1;  // we now have the root cell and no particles
  partCounter = 0;

  thetaMAC = _thetaMAC;
  gravConst = _gravConst;
  toptreeDepth = std::max(_czDepth,
                          _czDepth +
                          static_cast<size_t>(ceil(-log(thetaMAC) /
                                                   log(2.0))));

  if (toptreeDepth > 6)
    {
      std::cerr << "warning: toptree depth of " << toptreeDepth
                << " is pretty big!\n";
    }

  buildToptreeRecursor();
  noToptreeCells = cellCounter;

#ifdef SPHLATCH_MPI
  asLeaf().prepareBuffers();
#endif
};

///
/// destructor:
/// - postorder recurse node deletion
///
~BHtree(void)
{
  goRoot();
  empty();
  delete curNodePtr;                             // Seppuku!
}

protected:
nodePtrT curNodePtr, rootPtr;

///
/// variables
///
size_t cellCounter, partCounter, toptreeDepth, noToptreeCells;

matrixType localCells, remoteCells;
bitsetType localIsFilled, remoteIsFilled;

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
  particleNode* newNodePtr =
    new particleNode;

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
              curNodePtr->isLocal = true;
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

void insertParticle(partProxyPtrT _newPayload,
                    bool _newIsLocal)
{
  goRoot();
  insertParticleRecursor(_newPayload, _newIsLocal);
  partCounter++;
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
void insertParticleRecursor(partProxyPtrT _newPayload,
                            bool _newIsLocal)
{
  size_t targetOctant = 0;

  targetOctant += ((*_newPayload)(X) <
                   static_cast<cellPtrT>(curNodePtr)->xCenter) ? 0 : 1;
  targetOctant += ((*_newPayload)(Y) <
                   static_cast<cellPtrT>(curNodePtr)->yCenter) ? 0 : 2;
  targetOctant += ((*_newPayload)(Z) <
                   static_cast<cellPtrT>(curNodePtr)->zCenter) ? 0 : 4;

///
/// If targeted child is empty, place the particle there
///
  if (static_cast<cellPtrT>(curNodePtr)->child[targetOctant] == NULL)
    {
      curNodePtr->isEmpty = false;

      newPartChild(targetOctant);
      goChild(targetOctant);

      static_cast<partPtrT>(curNodePtr)->partProxyPtr = _newPayload;
      curNodePtr->isLocal = _newIsLocal;

/// particle saves its position to node directly
      static_cast<partPtrT>(curNodePtr)->xPos = (*_newPayload)(X);
      static_cast<partPtrT>(curNodePtr)->yPos = (*_newPayload)(Y);
      static_cast<partPtrT>(curNodePtr)->zPos = (*_newPayload)(Z);
      static_cast<partPtrT>(curNodePtr)->mass = (*_newPayload)(M);

      curNodePtr->ident =
        static_cast<identType>((*_newPayload)(PID));

///
/// don't forget to wire the nodePtr of the
/// NodeProxy back to the node
///
      static_cast<partPtrT>(curNodePtr)->partProxyPtr->nodePtr = curNodePtr;

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
      insertParticleRecursor(_newPayload, _newIsLocal);
      goUp();
    }

///
/// ... or if existing child is a particle (ghost/nonghost), then
/// replace it by a new node and try to place the existing two
/// particles as childs of this node
  else if (static_cast<cellPtrT>(curNodePtr)->child[targetOctant]->isParticle)
    {
///
/// goto child, save resident particle
/// and convert it to a node

      goChild(targetOctant);
      partProxyPtrT residentPayload =
        static_cast<partPtrT>(curNodePtr)->partProxyPtr;
      bool residentIsLocal = curNodePtr->isLocal;
      goUp();

      // replace particle by cell node
      delete static_cast<cellPtrT>(curNodePtr)->child[targetOctant];
      static_cast<cellPtrT>(curNodePtr)->child[targetOctant] = NULL;
      newCellChild(targetOctant);
      
      goChild(targetOctant);
      insertParticleRecursor(residentPayload, residentIsLocal);
      insertParticleRecursor(_newPayload, _newIsLocal);
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
  calcMultipoleRecursor();

  globalSumupMultipoles();
}

private:
///
/// recursor for multipole calculation
///
void calcMultipoleRecursor(void)
{
  if (calcMultipoleStop())
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
                  calcMultipoleRecursor();
                  goUp();
                }
            }
          detLocality();
          asLeaf().calcMultipole();
        }
    }
};


///
/// stop recursion if:
/// - current node is empty
/// - current node is a particle
///
bool calcMultipoleStop(void)
{
  return curNodePtr->isEmpty;                                                                                           // particles are per definition
  // particles are always empty, so we can omit this check
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
  if ((curNodePtr->depth) > toptreeDepth)
    {
      for (size_t i = 0; i < 8; i++)
        {
          if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
            {
              (curNodePtr->isLocal) |=
                static_cast<cellPtrT>(curNodePtr)->child[i]->isLocal;
            }
        }
    }
  else
    {
      curNodePtr->isLocal = true;
    }
}


///
/// sum up the multipole cellData matrix globally
/// and distribute it again to each node. if MPI is not
/// defined, nothing is done here.
///
protected:
size_t toptreeCounter;

private:
void globalSumupMultipoles()
{
 #ifdef SPHLATCH_MPI
  const size_t RANK = MPI::COMM_WORLD.Get_rank();
  const size_t SIZE = MPI::COMM_WORLD.Get_size();

  size_t round = 0;
  size_t remNodes = SIZE;

  const size_t noCellBytes
  = noToptreeCells * localCells.size2() * sizeof(valueType);

  std::queue<size_t> sumUpSend, sumUpRecv;
  std::stack<size_t> distrSend, distrRecv;

  //
  // magic algorithm which prepares sending and receiving queues
  // for the summing up step and sending and receiving stacks for
  // the distributing step.
  //
  while (remNodes > 1)
    {
      size_t noPairs = lrint(floor(remNodes / 2.));
      remNodes -= noPairs;
      size_t stepToNext = (1 << round);

      for (size_t i = 0; i < 2 * noPairs; i += 2)
        {
          size_t sendRank = (SIZE - 1) - stepToNext * (i + 1);
          size_t recvRank = (SIZE - 1) - stepToNext * i;

          if (RANK == sendRank)
            {
              sumUpSend.push(recvRank);
              distrRecv.push(recvRank);
            }
          else if (RANK == recvRank)
            {
              sumUpRecv.push(sendRank);
              distrSend.push(sendRank);
            }
          else
            {
            }
        }
      round++;
    }

  // toptree to buffers
  goRoot();
  toptreeCounter = 0;
  toptreeToBuffersRecursor();

  MPI::COMM_WORLD.Barrier();

  //
  // receive multipoles from other nodes and add
  // them to local value
  //
  while (!sumUpRecv.empty())
    {
      size_t recvFrom = sumUpRecv.front();
      sumUpRecv.pop();

      recvBitset(remoteIsFilled, recvFrom);

      MPI::COMM_WORLD.Recv(&remoteCells(0, 0), noCellBytes,
                           MPI_BYTE, recvFrom, RANK);

      asLeaf().mergeRemoteCells();

      // localIsFilled = localIsFilled or remoteIsFilled
      localIsFilled |= remoteIsFilled;
    }

  //
  // send local value to another node
  //
  while (!sumUpSend.empty())
    {
      size_t sendTo = sumUpSend.front();
      sumUpSend.pop();

      sendBitset(localIsFilled, sendTo);
      MPI::COMM_WORLD.Send(&localCells(0, 0), noCellBytes,
                           MPI_BYTE, sendTo, sendTo);
    }

  MPI::COMM_WORLD.Barrier();
  //
  // receive global result
  //
  while (!distrRecv.empty())
    {
      size_t recvFrom = distrRecv.top();
      distrRecv.pop();

      recvBitset(localIsFilled, recvFrom);

      MPI::COMM_WORLD.Recv(&localCells(0, 0), noCellBytes,
                           MPI_BYTE, recvFrom, RANK);
    }

  //
  // ... and distribute it to other nodes
  //
  while (!distrSend.empty())
    {
      size_t sendTo = distrSend.top();
      distrSend.pop();

      sendBitset(localIsFilled, sendTo);
      MPI::COMM_WORLD.Send(&localCells(0, 0), noCellBytes,
                           MPI_BYTE, sendTo, sendTo);
    }

  // buffers to toptree
  goRoot();
  toptreeCounter = 0;

  buffersToToptreeRecursor();
 #endif
}

#ifdef SPHLATCH_MPI
///
/// preorder recursive function to connect
/// cells to a matrix row
///
private:
void toptreeToBuffersRecursor(void)
{
  if (belowToptreeStop())
    {
    }
  else
    {
      localIsFilled[toptreeCounter] = !(curNodePtr->isEmpty);

      asLeaf().cellToBuffer();

      toptreeCounter++;

      if (curNodePtr->isParticle == false)
        {
          for (size_t i = 0; i < 8; i++)
            {
              if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
                {
                  goChild(i);
                  toptreeToBuffersRecursor();
                  goUp();
                }
            }
        }
    }
}

void buffersToToptreeRecursor(void)
{
  if (belowToptreeStop())
    {
    }
  else
    {
      curNodePtr->isEmpty = !localIsFilled[toptreeCounter];

      toptreeCounter++;

      if (curNodePtr->isParticle == false)
        {
          for (size_t i = 0; i < 8; i++)         // try without loop
            {
              if (static_cast<cellPtrT>(curNodePtr)->child[i] != NULL)
                {
                  goChild(i);
                  buffersToToptreeRecursor();
                  goUp();
                }
            }
        }
    }
}

///
/// stop recursion if:
///  - depth of current node is below toptreeDepth
///
bool belowToptreeStop(void)
{
  return(curNodePtr->depth > toptreeDepth);
}

private:
///
/// little comm helper
/// \todo move to comm_manager later on
///
void recvBitset(bitsetRefType _bitSet, size_t _sendRank)
{
  // prepare buffer
  const size_t RANK = MPI::COMM_WORLD.Get_rank();
  size_t noBlocks = lrint(ceil(static_cast<double>(noToptreeCells)
                               / (sizeof(bitsetBlockType) * 8)));
  size_t noBsBytes = noBlocks * sizeof(bitsetBlockType);

  std::vector<bitsetBlockType> recvBuff(noBlocks);

  // receive
  MPI::COMM_WORLD.Recv(&(recvBuff[0]), noBsBytes,
                       MPI_BYTE, _sendRank, RANK + 255);

  // copy back to bitset
  boost::from_block_range(recvBuff.begin(), recvBuff.end(),
                          _bitSet);
}

void sendBitset(bitsetRefType _bitSet, size_t _recvRank)
{
  // prepare buffer
  size_t noBlocks = lrint(ceil(static_cast<double>(noToptreeCells)
                               / (sizeof(bitsetBlockType) * 8)));
  size_t noBsBytes = noBlocks * sizeof(bitsetBlockType);

  std::vector<bitsetBlockType> sendBuff(noBlocks);

  // copy bitset to buffer
  boost::to_block_range(_bitSet, sendBuff.begin());

  // send
  MPI::COMM_WORLD.Send(&(sendBuff[0]), noBsBytes,
                       MPI_BYTE, _recvRank, _recvRank + 255);
}
#endif
// end of multipole stuff


// calcGravity() stuff
protected:
valueType curGravParticleX, curGravParticleY, curGravParticleZ;
valueType curGravParticleAX, curGravParticleAY, curGravParticleAZ;
valueType cellPartDist, cellPartDistPow3;
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
void calcGravity(partProxyPtrT _curParticle)
{  
  curGravParticleX = (*_curParticle)(X);
  curGravParticleY = (*_curParticle)(Y);
  curGravParticleZ = (*_curParticle)(Z);
  
  curGravParticleAX = 0.;
  curGravParticleAY = 0.;
  curGravParticleAZ = 0.;

  epsilonSquare = (*_curParticle)(GRAVEPS)*(*_curParticle)(GRAVEPS);

 #ifdef SPHLATCH_TREE_PROFILE
  calcGravityPartsCounter = 0;
  calcGravityCellsCounter = 0;
 #endif

  //
  // trick: hide the current particle by letting it look like
  // an empty cell node, so that it doesn't gravitate with
  // itself. btw: a particle is always empty.
  //
  _curParticle->nodePtr->isParticle = false;

  goRoot();
  calcGravityRecursor();

  (*_curParticle)(AX) += gravConst * curGravParticleAX;
  (*_curParticle)(AY) += gravConst * curGravParticleAY;
  (*_curParticle)(AZ) += gravConst * curGravParticleAZ;

  //
  // undo the trick above
  //
  _curParticle->nodePtr->isParticle = true;
}

private:
void calcGravityRecursor(void)
{
  if (curNodePtr->isParticle)
    {
      calcGravParticle();
    }
  else
  if (curNodePtr->isEmpty)
    {
    }
  else
  if (asLeaf().calcGravMAC())
    {
      asLeaf().calcGravCell();
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
valueType partGravPartnerX, partGravPartnerY, partGravPartnerZ,
          partGravPartnerM;
void calcGravParticle()
{
 #ifdef SPHLATCH_TREE_PROFILE
  calcGravityPartsCounter++;
 #endif
  partGravPartnerX = static_cast<partPtrT>(curNodePtr)->xPos;
  partGravPartnerY = static_cast<partPtrT>(curNodePtr)->yPos;
  partGravPartnerZ = static_cast<partPtrT>(curNodePtr)->zPos;
  partGravPartnerM = static_cast<partPtrT>(curNodePtr)->mass;

  cellPartDist = sqrt((partGravPartnerX - curGravParticleX) *
                      (partGravPartnerX - curGravParticleX) +
                      (partGravPartnerY - curGravParticleY) *
                      (partGravPartnerY - curGravParticleY) +
                      (partGravPartnerZ - curGravParticleZ) *
                      (partGravPartnerZ - curGravParticleZ)
                      );

  // todo: include spline softening
  cellPartDistPow3 = static_cast<valueType>(
    pow(cellPartDist * cellPartDist + epsilonSquare, 3. / 2.));

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
/*somecontainer_type findNeighbours(
                const nodePtrT _curParticle,
                const valueRefType _search_radius) {
        curNodePtr = _curParticle;
        size_t topDepth = <crazyformula>
        while ( curNodePtr->depth > topDepth ) {
                goUp();
        }
        neighbourFindRecursor(_curParticle);
        return neighbours;
   }*/

private:
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
           << " [label=\"\" ";

  if (curNodePtr->isParticle)
    {
      dumpFile << ",shape=circle";
    }
  else
    {
      dumpFile << ",shape=box";
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
  dumpFile << std::setprecision(9);

  if (curNodePtr->isParticle)
    {
      dumpFile << static_cast<partPtrT>(curNodePtr)->xPos << "   ";
      dumpFile << static_cast<partPtrT>(curNodePtr)->yPos << "   ";
      dumpFile << static_cast<partPtrT>(curNodePtr)->zPos << "   ";
      dumpFile << static_cast<partPtrT>(curNodePtr)->mass;
    }
  else
    {
      dumpFile << static_cast<cellPtrT>(curNodePtr)->xCenter << "   ";
      dumpFile << static_cast<cellPtrT>(curNodePtr)->yCenter << "   ";
      dumpFile << static_cast<cellPtrT>(curNodePtr)->zCenter << "   ";
      dumpFile << static_cast<cellPtrT>(curNodePtr)->cellSize << "   ";
      asLeaf().reportMultipoles();
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
/*void toptreeDump(std::string _dumpFileName)
   {
   dumpFile.open(_dumpFileName.c_str(), std::ios::out);

   dumpFile << "# toptree depth is " << toptreeDepth
           << " which gives us " << noToptreeCells
           << " toptree cells\n"
           << "#xCom            yCom             zCom"
           << "             q000           depth isEmpty\n";
   goRoot();
   toptreeDumpRecursor();
   dumpFile.close();
   }*/

private:
///
/// toptreeDump() recursor
///
/*void toptreeDumpRecursor()
   {
   dumpFile << std::scientific << std::setprecision(8);
   dumpFile << curNodePtr->xCom << "   "
           << curNodePtr->yCom << "   "
           << curNodePtr->zCom << "   "
           << curNodePtr->q000 << "   ";
   dumpFile << std::dec
           << curNodePtr->depth << "   "
           << static_cast<int>(curNodePtr->isEmpty) << "\n";

   if (atToptreeStop())
    {
    }
   else
    {
      for (size_t i = 0; i < 8; i++)             // try without loop
        {
          if (curNodePtr->child[i] != NULL)
            {
              goChild(i);
              toptreeDumpRecursor();
              goUp();
            }
        }
    }
   }*/
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
