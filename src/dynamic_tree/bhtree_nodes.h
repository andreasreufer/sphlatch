#ifndef BHTREE_NODES_H
#define BHTREE_NODES_H

/*
 *  bhtree_nodes.h
 *
 *
 *  Created by Andreas Reufer on 07.10.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"

namespace sphlatch {
class treeGhost;
class genericNode;
class genericCellNode;
class monopoleCellNode;
class quadrupoleCellNode;
class costzoneCellNode;
class particleNode;

typedef treeGhost             treeghoT;
typedef treeGhost*            treeghoPtrT;
typedef genericNode           nodeT;
typedef genericNode*          nodePtrT;
typedef genericCellNode       gcllT;
typedef genericCellNode*      gcllPtrT;
typedef monopoleCellNode      mcllT;
typedef monopoleCellNode*     mcllPtrT;
typedef quadrupoleCellNode    qcllT;
typedef quadrupoleCellNode*   qcllPtrT;
typedef costzoneCellNode      czllT;
typedef costzoneCellNode*     czllPtrT;
typedef particleNode          pnodT;
typedef particleNode*         pnodPtrT;

///
/// a generic tree node
///
class genericNode {
public:

   nodePtrT parent, next;

   idType             ident;
   short unsigned int depth;

   ///
   /// is a node a particle / a CZ cell, is it remote
   ///
   bool isParticle  : 1;
   bool isCZ        : 1;
   bool isRemote    : 1;

   ///
   /// indicate whether a CZ cell is at the bottom
   /// of the CZ tree and whether the neighbours are set
   ///
   bool atBottom    : 1;
   bool neighSet    : 1;

   ///
   /// indicates whether a particle node is settled
   /// at the correct position in the tree
   ///
   bool isSettled    : 1;
   bool needsUpdate  : 1;

   genericNode() { }
   ~genericNode() { }

   void clear();
   nodePtrT operator*();

private:
};

///
/// generic cell node
///
class genericCellNode : public genericNode {
public:
   nodePtrT skip;
   nodePtrT child[8];

   vect3dT cen;
   fType   clSz;

   genericCellNode() { }
   ~genericCellNode() { }

   void clear();
   void inheritCellPos(size_t _n);
   bool pointInsideCell(const vect3dT& _pos);
   size_t getOctant(const vect3dT& _pos);

private:
};

///
/// cell node with monopole moments
///
class monopoleCellNode : public genericCellNode {
public:
   vect3dT com;
   fType   m;

   monopoleCellNode() { }
   ~monopoleCellNode() { }

   void clear();

#ifdef SPHLATCH_PADD64
   char pad[0];
#endif
};

///
/// cell node with quadrupole moments
///
class quadrupoleCellNode : public monopoleCellNode {
public:

   fType      q11, q22, q33, q12, q13, q23;
   countsType noParts;

   quadrupoleCellNode() { }
   ~quadrupoleCellNode() { }

   void clear();
   void initFromCZll(czllT& _czll);

   void calcMultipole();

private:
#ifdef SPHLATCH_PADD64
   char pad[38];
#endif
};

///
/// a quadrupole costzone cell
///
class costzoneCellNode : public quadrupoleCellNode {
public:

   nodePtrT neighbour[27];
   idType   domain;

   fType absCost;

   ///
   /// adopted orphans
   ///
   pnodPtrT orphFrst, orphLast;

   ///
   /// first and last nodes of CZ cell subtree
   /// (chldFrst is always the same as next)
   ///
   nodePtrT chldFrst, chldLast;

   costzoneCellNode() { }
   ~costzoneCellNode() { }

   void clear();
   void initFromCell(qcllT& _qcll);
   void adopt(pnodPtrT _pnod);
   void pushdownNeighbours();
   void delSubtree();

private:
#ifdef SPHLATCH_PADD64
   char pad[56];
#endif
};

///
/// a particle node
///
class particleNode : public genericNode {
public:
   treeghoPtrT partPtr;
   vect3dT     pos;
   fType       m;

   particleNode() { }
   ~particleNode() { }

   void clear();
   void update();

#ifdef SPHLATCH_PADD64
private:
   char pad[0];
#endif
};
};
#endif
