#ifndef BHTREE_H
#define BHTREE_H

/*
 *  bhtree.h
 *
 *
 *  Created by Andreas Reufer on 12.11.07.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <vector>
#include <stack>
#include <queue>

#ifdef SPHLATCH_MPI
#include <mpi.h>
#endif

#include <fstream>

#include "typedefs.h"

#include "typedefs.h"

#include "bhtree_node_proxy.h"

#include "bhtree_node.h"

namespace sphlatch {
typedef NodeProxy NodeProxyType;
typedef NodeProxy*      NodeProxyPtrType;
typedef NodeProxy&      NodeProxyRefType;

typedef GenericOctNode<NodeProxyPtrType>* GenericOctNodePtr;

// hack
enum ParticleIndex { PID, X, Y, Z, VX, VY, VZ, AX, AY, AZ, M,
                     H, DHDT, RHO, E, P, POW, DIV_V,
                     ROTX_V, ROTY_V, ROTZ_V, Q, GRAVEPS };
enum MonopolIndex { CX, CY, CZ, Q000, MSIZE };
//enum QuadrupolIndex  { CX, CY, CZ, Q000, Q001, Q002, Q010, Q011, MSIZE};

//template <typename T>
class OctTree {
public:
/** \brief constructor:
 * - set theta, gravity constant and toptree depth
 * - check them for sanity
 * - allocate and setup RootNode and point cursor to it
 * - set counters
 * - instantiate recursors
 * - build toptree
 * - set up buffer matrices for toptree cell nodes
 */
OctTree(valueType _thetaMAC,
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

  RootPtr = new GenericOctNode<NodeProxyPtrType>;

  RootPtr->parent = NULL;
  for (size_t i = 0; i < 8; i++)
    {
      RootPtr->child[i] = NULL;
    }

  RootPtr->payload = NULL;
  RootPtr->ident = -1;

  RootPtr->xCenter = _rootCenter(0);
  RootPtr->yCenter = _rootCenter(1);
  RootPtr->zCenter = _rootCenter(2);
  RootPtr->cellSize = _rootSize;

  RootPtr->isParticle = false;
  RootPtr->isEmpty = true;
  RootPtr->isLocal = false;
  RootPtr->depth = 0;

  CurNodePtr = RootPtr;

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
  localCells.resize(noToptreeCells, MSIZE);
  localIsFilled.resize(noToptreeCells);

  remoteCells.resize(noToptreeCells, MSIZE);
  remoteIsFilled.resize(noToptreeCells);
#endif
};

/**
 * destructor:
 *  - postorder recurse node deletion
 */
~OctTree(void)
{
  goRoot();
  empty();
  delete CurNodePtr;                       // Seppuku!
}

private:
GenericOctNodePtr CurNodePtr, RootPtr;

/**
 * variables
 */
size_t cellCounter, partCounter, toptreeDepth, noToptreeCells;

matrixType localCells, remoteCells;
bitsetType localIsFilled, remoteIsFilled;

std::fstream logFile;

// little helpers
private:
/**
 * go up one level
 */
inline void goUp()
{
  CurNodePtr = CurNodePtr->parent;
}

/**
 * go to child
 */
inline void goChild(const size_t _n)
{
  CurNodePtr = CurNodePtr->child[_n];
}

/**
 * go to root
 */
inline void goRoot()
{
  CurNodePtr = RootPtr;
}

/**
 * setup coordinates and cell size as child _n
 */
void setupCoordinates(const size_t _n)
{
  if (CurNodePtr->parent != NULL)
    {
      CurNodePtr->cellSize = CurNodePtr->parent->cellSize / 2.;
      CurNodePtr->xCenter = ((_n) % 2) ?
                            CurNodePtr->parent->xCenter
                            + 0.5 * CurNodePtr->cellSize :
                            CurNodePtr->parent->xCenter
                            - 0.5 * CurNodePtr->cellSize;
      CurNodePtr->yCenter = ((_n >> 1) % 2) ?
                            CurNodePtr->parent->yCenter
                            + 0.5 * CurNodePtr->cellSize :
                            CurNodePtr->parent->yCenter
                            - 0.5 * CurNodePtr->cellSize;
      CurNodePtr->zCenter = ((_n >> 2) % 2) ?
                            CurNodePtr->parent->zCenter
                            + 0.5 * CurNodePtr->cellSize :
                            CurNodePtr->parent->zCenter
                            - 0.5 * CurNodePtr->cellSize;
    }
}

/**
 * create a new child in octant _n
 */
void newChild(const size_t _n)
{
  if (CurNodePtr->child[_n] == NULL)
    {
      //GenericOctNodePtr NewNodePtr = new GenericOctNode<T>;
      GenericOctNodePtr NewNodePtr =
        new GenericOctNode<NodeProxyPtrType>;

      NewNodePtr->parent = CurNodePtr;
      CurNodePtr->child[_n] = NewNodePtr;

      goChild(_n);

      CurNodePtr->child[0] = NULL;
      CurNodePtr->child[1] = NULL;
      CurNodePtr->child[2] = NULL;
      CurNodePtr->child[3] = NULL;
      CurNodePtr->child[4] = NULL;
      CurNodePtr->child[5] = NULL;
      CurNodePtr->child[6] = NULL;
      CurNodePtr->child[7] = NULL;

      CurNodePtr->depth = CurNodePtr->parent->depth + 1;

      setupCoordinates(_n);

      goUp();
    }
}

/**
 * convert current node to cell, if current node is a particle
 */
void converttoCell(const size_t _n)
{
  if (CurNodePtr->isParticle)
    {
      setupCoordinates(_n);

      CurNodePtr->payload = NULL;

      CurNodePtr->ident =
        static_cast<identType>(-cellCounter - 1);

      CurNodePtr->isParticle = false;
      CurNodePtr->isEmpty = true;

      CurNodePtr->isLocal = false;

      cellCounter++;
    }
};
// end of little helpers


// top tree building stuff
private:
/**
 * recursor to build toptree
 */
void buildToptreeRecursor(void)
{
  if (atToptreeStop())
    {
    }
  else
    {
      makeEmptyCells();
      for (size_t i = 0; i < 8; i++)  // try without loop
        {
          if (CurNodePtr->child[i] != NULL)
            {
              goChild(i);
              buildToptreeRecursor();
              goUp();
            }
        }
    }
}

/**
 * make empty cell nodes for the toptree
 */
void makeEmptyCells(void)
{
  for (size_t i = 0; i < 8; i++)
    {
      newChild(i);
      goChild(i);

      CurNodePtr->payload = NULL;
      CurNodePtr->ident =
        static_cast<identType>(-cellCounter - 1);

      CurNodePtr->isParticle = false;
      CurNodePtr->isEmpty = true;

      CurNodePtr->isLocal = false;

      cellCounter++;
      goUp();
    }
}

/**
 * stop recursion if:
 * - depth of current node is below toptreeDepth
 */
bool atToptreeStop(void)
{
  return(CurNodePtr->depth >= toptreeDepth);
};
// end of top tree stuff


// insertParticle() stuff
public:
/**
 * method to insert particle:
 *  - go to root
 *  - call the insertion recursor
 */
void insertParticle(NodeProxyType* _newPayload,
                    bool _newIsLocal)
{
  goRoot();
  insertParticleRecursor(_newPayload, _newIsLocal);
  partCounter++;
}

private:
/**
 * recursor for inserting a new particle:
 * try to insert as child of current
 * node. if child is
 *	- empty, insert particle. we're done.
 *	- a node, go to child and call recursor
 *	- a particle, disconnect particle and
 *	  call recursor for existing and new
 *	  particle.
 */
void insertParticleRecursor(NodeProxyType* _newPayload,
                            bool _newIsLocal)
{
  size_t targetOctant = 0;

  targetOctant += ((*_newPayload)(X) < CurNodePtr->xCenter) ? 0 : 1;
  targetOctant += ((*_newPayload)(Y) < CurNodePtr->yCenter) ? 0 : 2;
  targetOctant += ((*_newPayload)(Z) < CurNodePtr->zCenter) ? 0 : 4;

  /**
   * If targeted child is empty, place the particle there
   **/
  if (CurNodePtr->child[targetOctant] == NULL)
    {
      CurNodePtr->isEmpty = false;

      newChild(targetOctant);
      goChild(targetOctant);

      CurNodePtr->payload = _newPayload;
      CurNodePtr->isParticle = true;
      CurNodePtr->isEmpty = true;
      CurNodePtr->isLocal = _newIsLocal;

      /* particle saves its position to node directly */
      CurNodePtr->xCom = (*_newPayload)(X);
      CurNodePtr->yCom = (*_newPayload)(Y);
      CurNodePtr->zCom = (*_newPayload)(Z);
      CurNodePtr->q000 = (*_newPayload)(M);

      CurNodePtr->ident =
        static_cast<identType>((*_newPayload)(PID));

      /**
       * don't forget to wire the nodePtr of the
       * NodeProxy back to the node
       */
      CurNodePtr->payload->nodePtr = CurNodePtr;

      goUp();
    }

  /**
   * ... or if existing child is a node, then try to place the particle
   * as a child of this node
   */
  else if (!CurNodePtr->child[targetOctant]->isParticle)
    {
      CurNodePtr->isEmpty = false;
      goChild(targetOctant);
      insertParticleRecursor(_newPayload, _newIsLocal);
      goUp();
    }

  /**
   * ... or if existing child is a particle (ghost/nonghost), then
   * replace it by a new node and try to place the existing two
   * particles as childs of this node
   */
  else if (CurNodePtr->child[targetOctant]->isParticle)
    {
      /**
       * goto child, save resident particle
       * and convert it to a node
       */
      goChild(targetOctant);

      NodeProxyType*  residentPayload = CurNodePtr->payload;
      bool residentIsLocal = CurNodePtr->isLocal;

      converttoCell(targetOctant);

      insertParticleRecursor(residentPayload, residentIsLocal);
      insertParticleRecursor(_newPayload, _newIsLocal);
    }
  else
    {
      /* this never happens */
    }
}
// end of insertParticle() stuff


// calcMultipoles() stuff
public:
/**
 * calculate multipoles:
 *  - go to root
 *  - call the multipole recursor
 *  - exchange toptrees
 */
void calcMultipoles(void)
{
  goRoot();
  calcMultipoleRecursor();

  globalSumupMultipoles();
}

private:
/**
 * recursor for multipole calculation
 */
void calcMultipoleRecursor(void)
{
  if (calcMultipoleStop())
    {
    }
  else
    {
      for (size_t i = 0; i < 8; i++) // try without loop
        {
          if (CurNodePtr->child[i] != NULL)
            {
              goChild(i);
              calcMultipoleRecursor();
              goUp();
            }
        }
      calcMultipole();
    }
};

/**
 * stop recursion if:
 * - current node is empty
 * - current node is a particle
 */
bool calcMultipoleStop(void)
{
  return CurNodePtr->isEmpty; // particles are per definition
                              // empty, so we omit this check
};

/**
 * calculate multipole from children
 */
valueType monopolCM, monopolCXM, monopolCYM, monopolCZM;
void calcMultipole()
{
  monopolCM = 0.;
  monopolCXM = 0.;
  monopolCYM = 0.;
  monopolCZM = 0.;
  /**
   * check locality of current cell node:
   *  - if node is in the toptree, the node is local
   *  - if node is below toptree, the parent node is local if
   *    any child is local ( ... or ... or ... or ... )
   *
   * to check the proper working together of the costzone and
   * the toptree would be to check, that cells below toptreeDepth
   * either have only local or only non-local children, but never
   * both.
   */
  if ((CurNodePtr->depth) > toptreeDepth)
    {
      for (size_t i = 0; i < 8; i++)
        {
          if (CurNodePtr->child[i] != NULL)
            {
              (CurNodePtr->isLocal) |= CurNodePtr->child[i]->isLocal;
            }
        }
    }
  else
    {
      CurNodePtr->isLocal = true;
    }

  /**
   * add up the contributions from the children with
   * the same locality as the current node.
   * all this locality business guarantees, that every particle
   * contributes only ONCE to the multipole moments of the global tree
   * but "ghost" cells still contain the multipoles contributed by
   * their ghost children.
   *
   * a special case are the deepest toptree cells which are per
   * definition local. if all children are ghosts, none of if contri-
   * butes anything so monopolCM gets 0 and fucks up the center of mass
   * of the cell. so if nobody contributes anything, omit the addition
   * to the cell.
   */
  for (size_t i = 0; i < 8; i++)
    {
      if (CurNodePtr->child[i] != NULL)
        {
          if (CurNodePtr->child[i]->isLocal
              == CurNodePtr->isLocal)
            {
              goChild(i);
              monopolCM += CurNodePtr->q000;
              monopolCXM += (CurNodePtr->xCom) *
                            (CurNodePtr->q000);
              monopolCYM += (CurNodePtr->yCom) *
                            (CurNodePtr->q000);
              monopolCZM += (CurNodePtr->zCom) *
                            (CurNodePtr->q000);
              goUp();
            }
        }
    }

  /* copy data to node itself ... */
  if (monopolCM > 0.)
    {
      CurNodePtr->q000 = monopolCM;
      CurNodePtr->xCom = monopolCXM / monopolCM;
      CurNodePtr->yCom = monopolCYM / monopolCM;
      CurNodePtr->zCom = monopolCZM / monopolCM;
    }
  else
    {
      CurNodePtr->q000 = 0.;
      CurNodePtr->xCom = 0.;
      CurNodePtr->yCom = 0.;
      CurNodePtr->zCom = 0.;
    }
}

/**
 * sum up the multipole cellData matrix globally
 * and distribute it again to each node. if MPI is not
 * defined, nothing is done here.
 */
void globalSumupMultipoles()
{
#ifdef SPHLATCH_MPI
  const size_t RANK = MPI::COMM_WORLD.Get_rank();
  const size_t SIZE = MPI::COMM_WORLD.Get_size();

#ifdef SPHLATCH_TREE_LOGSUMUPMP
  std::string logFilename = "MPlogRank000";
  std::string rankString = boost::lexical_cast<std::string>(RANK);
  logFilename.replace(logFilename.size() - 0 - rankString.size(),
                      rankString.size(), rankString);
  double logStartTime = MPI_Wtime();
  logFile.open(logFilename.c_str(), std::ios::out);
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - logStartTime << "    start log\n";
#endif

  size_t round = 0;
  size_t remNodes = SIZE;

  const size_t noCellBytes = noToptreeCells * MSIZE * sizeof(valueType);

  std::queue<size_t> sumUpSend, sumUpRecv;
  std::stack<size_t> distrSend, distrRecv;

  /**
   * magic algorithm which prepares sending and receiving queues
   * for the summing up step and sending and receiving stacks for
   * the distributing step.
   */
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

#ifdef SPHLATCH_TREE_LOGSUMUPMP
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - logStartTime << "   "
          << " tree -> buffers \n" << std::flush;
#endif

  MPI::COMM_WORLD.Barrier();

#ifdef SPHLATCH_TREE_LOGSUMUPMP
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - logStartTime << "   "
          << " binary exchange tree up\n" << std::flush;
#endif

  /**
   * receive multipoles from other nodes and add
   * them to local value
   */
  while (!sumUpRecv.empty())
    {
      size_t recvFrom = sumUpRecv.front();
      sumUpRecv.pop();

#ifdef SPHLATCH_TREE_LOGSUMUPMP
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "   "
              << " receive bitset from " << recvFrom << "\n" << std::flush;
#endif

      recvBitset(remoteIsFilled, recvFrom);

#ifdef SPHLATCH_TREE_LOGSUMUPMP
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "   "
              << " receive cells  from " << recvFrom << "\n" << std::flush;
#endif

      MPI::COMM_WORLD.Recv(&remoteCells(0, 0), noCellBytes,
                           MPI_BYTE, recvFrom, RANK);

#ifdef SPHLATCH_TREE_LOGSUMUPMP
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "   "
              << " add up multipoles \n" << std::flush;
#endif

      valueType oldQ000, newQ000;
      for (size_t i = 0; i < noToptreeCells; i++)
        {
          if (remoteIsFilled[i])
            {
              if (localIsFilled[i])
                {
                  oldQ000 = localCells(i, Q000);
                  newQ000 = oldQ000 + remoteCells(i, Q000);

                  /**
                   * newQ000 may be zero for non-empty cells. this
                   * happens when all children are ghosts and do not
                   * contribute to the local toptree
                   */
                  if (newQ000 > 0.)
                    {
                      localCells(i, Q000) = newQ000;
                      localCells(i, CX) = ((oldQ000 / newQ000) *
                                           localCells(i, CX))
                                          + ((remoteCells(i, Q000) / newQ000) *
                                             remoteCells(i, CX));
                      localCells(i, CY) = ((oldQ000 / newQ000) *
                                           localCells(i, CY))
                                          + ((remoteCells(i, Q000) / newQ000) *
                                             remoteCells(i, CY));
                      localCells(i, CZ) = ((oldQ000 / newQ000) *
                                           localCells(i, CZ))
                                          + ((remoteCells(i, Q000) / newQ000) *
                                             remoteCells(i, CZ));
                    }
                }
              else
                {
                  localCells(i, Q000) = remoteCells(i, Q000);
                  localCells(i, CX) = remoteCells(i, CX);
                  localCells(i, CY) = remoteCells(i, CY);
                  localCells(i, CZ) = remoteCells(i, CZ);
                }
            }
        }

      // localIsFilled = localIsFilled or remoteIsFilled
      localIsFilled |= remoteIsFilled;
    }

  /**
   * send local value to another node
   */
  while (!sumUpSend.empty())
    {
      size_t sendTo = sumUpSend.front();
      sumUpSend.pop();

#ifdef SPHLATCH_TREE_LOGSUMUPMP
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "   "
              << " send bitset to " << sendTo << "\n" << std::flush;
#endif

      sendBitset(localIsFilled, sendTo);

#ifdef SPHLATCH_TREE_LOGSUMUPMP
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "   "
              << " send cells to " << sendTo << "\n" << std::flush;
#endif

      /**
       * do a synchronous send, which is non-blocking but
       * still prevents the receiving node from getting
       * bombarded by multiple sends
       * */
      MPI::COMM_WORLD.Send(&localCells(0, 0), noCellBytes,
                           MPI_BYTE, sendTo, sendTo);

#ifdef SPHLATCH_TREE_LOGSUMUPMP
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "   "
              << " send to " << sendTo << " done\n" << std::flush;
#endif
    }


  MPI::COMM_WORLD.Barrier();
#ifdef SPHLATCH_TREE_LOGSUMUPMP
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - logStartTime << "   "
          << " binary excghange tree down\n" << std::flush;
#endif

  /**
   * receive global result
   */
  while (!distrRecv.empty())
    {
      size_t recvFrom = distrRecv.top();
      distrRecv.pop();

#ifdef SPHLATCH_TREE_LOGSUMUPMP
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "   "
              << " receive bitset from " << recvFrom << "\n" << std::flush;
#endif

      recvBitset(localIsFilled, recvFrom);

#ifdef SPHLATCH_TREE_LOGSUMUPMP
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "   "
              << " receive cells  from " << recvFrom << "\n" << std::flush;
#endif

      MPI::COMM_WORLD.Recv(&localCells(0, 0), noCellBytes,
                           MPI_BYTE, recvFrom, RANK);

#ifdef SPHLATCH_TREE_LOGSUMUPMP
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "   "
              << " receive from " << recvFrom << " done\n" << std::flush;
#endif
    }

  /**
   * ... and distribute it to other nodes
   */
  while (!distrSend.empty())
    {
      size_t sendTo = distrSend.top();
      distrSend.pop();

#ifdef SPHLATCH_TREE_LOGSUMUPMP
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "   "
              << " send bitset to " << sendTo << "\n" << std::flush;
#endif

      sendBitset(localIsFilled, sendTo);

#ifdef SPHLATCH_TREE_LOGSUMUPMP
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "   "
              << " send cells to " << sendTo << "\n" << std::flush;
#endif

      /**
       * do a synchronous send, which is non-blocking but
       * still prevents the receiving node from getting
       * bombarded by multiple sends
       * */
      MPI::COMM_WORLD.Send(&localCells(0, 0), noCellBytes,
                           MPI_BYTE, sendTo, sendTo);

#ifdef SPHLATCH_TREE_LOGSUMUPMP
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "   "
              << " send to " << sendTo << " done\n" << std::flush;
#endif
    }

  // buffers to toptree
  goRoot();
  toptreeCounter = 0;

#ifdef SPHLATCH_TREE_LOGSUMUPMP
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - logStartTime << "   "
          << " buffers -> tree \n" << std::flush;
#endif

  buffersToToptreeRecursor();

#ifdef SPHLATCH_TREE_LOGSUMUPMP
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - logStartTime << "   "
          << " globalSumup() done\n" << std::flush;
  logFile.close();
#endif
#endif
}

#ifdef SPHLATCH_MPI
/**
 * preorder recursive function to connect
 * cells to a matrix row
 */
size_t toptreeCounter;
void toptreeToBuffersRecursor(void)
{
  if (belowToptreeStop())
    {
    }
  else
    {
      localIsFilled[toptreeCounter] = !(CurNodePtr->isEmpty);

      localCells(toptreeCounter, CX) = CurNodePtr->xCom;
      localCells(toptreeCounter, CY) = CurNodePtr->yCom;
      localCells(toptreeCounter, CZ) = CurNodePtr->zCom;
      localCells(toptreeCounter, Q000) = CurNodePtr->q000;

      toptreeCounter++;

      for (size_t i = 0; i < 8; i++)             // try without loop
        {
          if (CurNodePtr->child[i] != NULL)
            {
              goChild(i);
              toptreeToBuffersRecursor();
              goUp();
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
      CurNodePtr->isEmpty = !localIsFilled[toptreeCounter];

      CurNodePtr->xCom = localCells(toptreeCounter, CX);
      CurNodePtr->yCom = localCells(toptreeCounter, CY);
      CurNodePtr->zCom = localCells(toptreeCounter, CZ);
      CurNodePtr->q000 = localCells(toptreeCounter, Q000);

      toptreeCounter++;

      for (size_t i = 0; i < 8; i++)             // try without loop
        {
          if (CurNodePtr->child[i] != NULL)
            {
              goChild(i);
              buffersToToptreeRecursor();
              goUp();
            }
        }
    }
}

/**
 * stop recursion if:
 * - depth of current node is below toptreeDepth
 */
bool belowToptreeStop(void)
{
  return(CurNodePtr->depth > toptreeDepth);
}

/**
 * little comm helper
 * \todo move to comm_manager later on
 */
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
  //const size_t RANK = MPI::COMM_WORLD.Get_rank();
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
private:
NodeProxyType* curGravNodeProxy;
valueType curGravParticleX, curGravParticleY, curGravParticleZ;
valueType curGravParticleAX, curGravParticleAY, curGravParticleAZ;
valueType thetaMAC, gravConst;
valueType epsilonSquare;
size_t calcGravityCellsCounter, calcGravityPartsCounter;

public:
/**
 * calculate gravitation for a particle:
 *  - load current particle data
 *  - go to root
 *  - call the gravity calculation recursor
 *  - write back resulting acceleration
 */
void calcGravity(NodeProxyType* _curParticle)
{
  curGravNodeProxy = _curParticle;

  curGravParticleX = (*curGravNodeProxy)(X);
  curGravParticleY = (*curGravNodeProxy)(Y);
  curGravParticleZ = (*curGravNodeProxy)(Z);
  curGravParticleAX = 0.;
  curGravParticleAY = 0.;
  curGravParticleAZ = 0.;

  epsilonSquare = (*curGravNodeProxy)(GRAVEPS)*(*curGravNodeProxy)(GRAVEPS);

#ifdef SPHLATCH_TREE_PROFILE
  calcGravityPartsCounter = 0;
  calcGravityCellsCounter = 0;
#endif

  /**
   * trick: hide the current particle by letting it look like
   * an empty cell node, so that it doesn't gravitate with
   * itself. btw: a particle is always empty.
   */
  curGravNodeProxy->nodePtr->isParticle = false;

  goRoot();
  calcGravityRecursor();

  (*curGravNodeProxy)(AX) += gravConst * curGravParticleAX;
  (*curGravNodeProxy)(AY) += gravConst * curGravParticleAY;
  (*curGravNodeProxy)(AZ) += gravConst * curGravParticleAZ;

  /**
   * undo the trick above
   */
  curGravNodeProxy->nodePtr->isParticle = true;
}

private:
void calcGravityRecursor(void)
{
  if (CurNodePtr->isParticle)
    {
      calcGravParticle();
    }
  else
  if (CurNodePtr->isEmpty)
    {
    }
  else
  if (calcGravMAC())
    {
      calcGravCell();
    }
  else
    {
      for (size_t i = 0; i < 8; i++)
        {
          if (CurNodePtr->child[i] != NULL)
            {
              goChild(i);
              calcGravityRecursor();
              goUp();
            }
        }
    }
}

/**
 * stop recursion if:
 * - current node is empty << ??
 * - MAC is fulfilled
 */
valueType cellPartDist;
bool calcGravMAC(void)
{
  cellPartDist = sqrt(
    (CurNodePtr->xCom - curGravParticleX) *
    (CurNodePtr->xCom - curGravParticleX) +
    (CurNodePtr->yCom - curGravParticleY) *
    (CurNodePtr->yCom - curGravParticleY) +
    (CurNodePtr->zCom - curGravParticleZ) *
    (CurNodePtr->zCom - curGravParticleZ)
    );
  return(((CurNodePtr->cellSize) / cellPartDist) < thetaMAC);
}

/**
 * calculate acceleration due to a particle
 * todo: add grav-const
 */
valueType partGravPartnerX, partGravPartnerY, partGravPartnerZ,
          partGravPartnerM;
void calcGravParticle()
{
#ifdef SPHLATCH_TREE_PROFILE
  calcGravityPartsCounter++;
#endif
  partGravPartnerX = CurNodePtr->xCom;
  partGravPartnerY = CurNodePtr->yCom;
  partGravPartnerZ = CurNodePtr->zCom;
  partGravPartnerM = CurNodePtr->q000;

  cellPartDist = sqrt((partGravPartnerX - curGravParticleX) *
                      (partGravPartnerX - curGravParticleX) +
                      (partGravPartnerY - curGravParticleY) *
                      (partGravPartnerY - curGravParticleY) +
                      (partGravPartnerZ - curGravParticleZ) *
                      (partGravPartnerZ - curGravParticleZ)
                      );

  // todo: include spline softening
  /*cellPartDistPow3 = cellPartDist * cellPartDist * cellPartDist
   + cellPartDist * epsilonSquare;*/

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

/**
 * calculate acceleration due to a cell
 * todo: add grav-const
 */
valueType cellPartDistPow3;
void calcGravCell()
{
#ifdef SPHLATCH_TREE_PROFILE
  calcGravityCellsCounter++;
#endif
  // cellPartDist is already set by the MAC function
  
  // no softening for cells
  /*cellPartDistPow3 = cellPartDist * cellPartDist * cellPartDist
                     + cellPartDist * epsilonSquare;*/
  cellPartDistPow3 = cellPartDist * cellPartDist * cellPartDist;

  curGravParticleAX -= (CurNodePtr->q000) *
                       (curGravParticleX - CurNodePtr->xCom) /
                       cellPartDistPow3;
  curGravParticleAY -= (CurNodePtr->q000) *
                       (curGravParticleY - CurNodePtr->yCom) /
                       cellPartDistPow3;
  curGravParticleAZ -= (CurNodePtr->q000) *
                       (curGravParticleZ - CurNodePtr->zCom) /
                       cellPartDistPow3;
}
// end of calcGravity() stuff

// neighbour search stuff
public:
/**
 * find neighbours:
 *  - load current particle data
 *  - go to current particle node
 *  - go up until every neighbour is within current cell
 *  - call the neighbour search recursor
 *  - brute force sort out non-neighbours
 *  - return neighbours
 */
/*somecontainer_type findNeighbours(
                const GenericOctNodePtr _curParticle,
                const valueRefType _search_radius) {
        CurNodePtr = _curParticle;
        size_t topDepth = <crazyformula>
        while ( CurNodePtr->depth > topDepth ) {
                goUp();
        }
        neighbourFindRecursor(_curParticle);
        return neighbours;
   }*/

private:
// end of neighbour search stuff

// treeDOTDump() stuff
public:
/**
 * dump the tree as a
 * dot file
 */
std::fstream dumpFile;
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
/**
 * treeDOTDump() recursor
 */
void treeDOTDumpRecursor()
{
  dumpFile << CurNodePtr->ident
           << " [label=\"\" ";

  if (CurNodePtr->isParticle)
    {
      dumpFile << ",shape=circle";
    }
  else
    {
      dumpFile << ",shape=box";
    }

  if ((!CurNodePtr->isEmpty) or
      CurNodePtr->isParticle)
    {
      dumpFile << ",style=filled";
    }

  if (CurNodePtr->isLocal)
    {
      // local particle
      if (CurNodePtr->isParticle)
        {
          dumpFile << ",color=green";
        }
      else
        {
          // local cell
          if (CurNodePtr->depth <= toptreeDepth)
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

  for (size_t i = 0; i < 8; i++)                         // try without loop
    {
      if (CurNodePtr->child[i] != NULL)
        {
          dumpFile << CurNodePtr->ident
                   << " -> "
                   << CurNodePtr->child[i]->ident
                   << "\n";
          goChild(i);
          treeDOTDumpRecursor();
          goUp();
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
  if (CurNodePtr->isParticle)
    {
      dumpFile << "P";
    }
  else
    {
      dumpFile << "N";
    }
  if (CurNodePtr->isEmpty)
    {
      dumpFile << "E";
    }
  else
    {
      dumpFile << "F";
    }
  if (CurNodePtr->isLocal)
    {
      dumpFile << "L";
    }
  else
    {
      dumpFile << "R";
    }
  dumpFile << "   ";

  dumpFile << CurNodePtr->ident;
  dumpFile << "   ";

  if (CurNodePtr->isParticle)
    {
      dumpFile << CurNodePtr->xCom << "   ";
      dumpFile << CurNodePtr->yCom << "   ";
      dumpFile << CurNodePtr->zCom << "   ";
      dumpFile << CurNodePtr->q000;
    }
  else
    {
      dumpFile << CurNodePtr->xCenter << "   ";
      dumpFile << CurNodePtr->yCenter << "   ";
      dumpFile << CurNodePtr->zCenter << "   ";
      dumpFile << CurNodePtr->cellSize;
    }

  dumpFile << "\n";

  for (size_t i = 0; i < 8; i++)                         // try without loop
    {
      if (CurNodePtr->child[i] != NULL)
        {
          goChild(i);
          treeDumpRecursor();
          goUp();
        }
    }
}
// end of treeDump() stuff

// toptreeDump() stuff
public:
/**
 * dump the toptree
 */
void toptreeDump(std::string _dumpFileName)
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
}

private:
/**
 * toptreeDump() recursor
 */
void toptreeDumpRecursor()
{
  dumpFile << std::scientific << std::setprecision(8);
  dumpFile << CurNodePtr->xCom << "   "
           << CurNodePtr->yCom << "   "
           << CurNodePtr->zCom << "   "
           << CurNodePtr->q000 << "   ";
  dumpFile << std::dec
           << CurNodePtr->depth << "   "
           << static_cast<int>(CurNodePtr->isEmpty) << "\n";

  if (atToptreeStop())
    {
    }
  else
    {
      for (size_t i = 0; i < 8; i++)             // try without loop
        {
          if (CurNodePtr->child[i] != NULL)
            {
              goChild(i);
              toptreeDumpRecursor();
              goUp();
            }
        }
    }
}
// end of toptreeDump() stuff

// empty() stuff
/**
 * deletes the current subtree
 */
private:
void empty(void)
{
  emptyRecursor();
}

/**
 * recursor for emptying the
 * tree
 */
void emptyRecursor()
{
  for (size_t i = 0; i < 8; i++)                         // try without loop
    {
      if (CurNodePtr->child[i] != NULL)
        {
          goChild(i);
          emptyRecursor();
          goUp();

          delete CurNodePtr->child[i];
          CurNodePtr->child[i] = NULL;
        }
    }
};
// end of empty() stuff
};
};
#endif

