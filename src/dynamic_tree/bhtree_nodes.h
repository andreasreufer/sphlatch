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

typedef genericNode           nodeT;
typedef genericNode*          nodePtrT;
typedef genericCellNode       cellT;
typedef genericCellNode*      cellPtrT;
typedef monopoleCellNode      mcllT;
typedef monopoleCellNode*     mcllPtrT;
typedef quadrupoleCellNode    qcllT;
typedef quadrupoleCellNode*   qcllPtrT;
typedef costzoneCellNode      czllT;
typedef costzoneCellNode*     czllPtrT;
typedef particleNode          pnodT;
typedef particleNode*         pnodPtrT;

typedef treeGhost             treeghoT;
typedef treeGhost*            treeghoPtrT;

class particleNode;

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

   // what for should this be good?
   //countsType noParts;

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
   fType   mass;

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

   fType q11, q22, q33, q12, q13, q23;

   quadrupoleCellNode() { }
   ~quadrupoleCellNode() { }

   void clear();
   void initFromCZll(czllT& _czll);

private:
#ifdef SPHLATCH_PADD64
   char pad[46];
#endif
};


///
/// a quadrupole costzone cell
///
class costzoneCellNode : public quadrupoleCellNode {
public:

   nodePtrT neighbour[27];
   idType   domain;

   fType      absCost, relCost;
   countsType noParts;

   pnodPtrT adopFrst, adopLast;
   nodePtrT chldFrst, chldLast;

   costzoneCellNode() { }
   ~costzoneCellNode() { }

   void clear();
   void initFromCell(cellT& _cell);

   void adopt(pnodPtrT _pnod);

   void pushdownNeighbours();

private:
#ifdef SPHLATCH_PADD64
   char pad[40];
#endif
};



class particleNode : public genericNode {
public:

   partPtrT partPtr;
   vect3dT  pos;
   fType    m;

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
